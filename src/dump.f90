module dump
    !
    use atom_types
    use parameters
    use comms
    use elements
    use box
    use errors
    use forces
    use math
    use group
    use time
    use dump_logic

    implicit none

    !Dump command parameters
    integer, parameter :: max_dumps = 20
    integer, private :: dump_num 
    integer, dimension(max_dumps), private :: dump_every, dump_type, dump_group
    character(len = read_len), dimension(max_dumps), private :: base_filename, dump_name


    !Gather arrays for atoms
    integer, parameter :: max_send_count = 1024
    integer, allocatable, save :: tag_atom_gather(:), type_atom_gather(:), a_mask_gather(:)
    real(kind=wp), allocatable, save :: r_atom_gather(:,:), energy_atom_gather(:), force_atom_gather(:,:), &
                                        virial_atom_gather(:,:,:), vel_atom_gather(:,:)
    !Gather arrays for elements
    integer, allocatable, save :: tag_ele_gather(:), size_ele_gather(:), etype_gather(:), &
                                  basis_num_gather(:), basis_type_gather(:,:), e_mask_gather(:)
    real(kind=wp), allocatable, save :: r_gather(:,:,:), energy_gather(:,:), force_gather(:,:,:), virial_gather(:,:,:,:), &
                                        vel_gather(:,:,:)

    !Gather arrays for interpolated atoms
    integer, allocatable :: atomap_mask_gather(:), type_atomap_gather(:)
    real(kind = wp), allocatable, save :: r_atomap_gather(:,:)


    public 
    contains

    subroutine dump_defaults
        dump_every = huge(1)
        dump_num = 0
    end subroutine dump_defaults

    subroutine parse_dump(line)
        !Parse the dump command
        character(len = *), intent(in) :: line
        character(len = read_len) :: label, tmptxt, msg, dtype, dgroup
        integer :: iospara

        dump_num = dump_num + 1
        read(line,*, iostat=iospara) label, dump_name(dump_num), dgroup, dtype, dump_every(dump_num), base_filename(dump_num)
        if(iospara > 0) call read_error("Failure to read dump command", iospara)
        
        !Check to make sure dump file format is correct
        select case(dtype)
        case('out')
            dump_type(dump_num) = 1
        case('lmp')
            dump_type(dump_num) = 2
        case default
            write(tmptxt, *) "dump file format ", trim(adjustl(dtype)), " is not accepted"
            call command_error(tmptxt)
        end select

        write(msg, *) "Writing output to ", trim(adjustl(base_filename(dump_num))), " using ", trim(adjustl(dtype)), " file format"
        call log_msg(msg)

        !Now get rid of the extension on the filename
        base_filename(dump_num) = strip_dump_extension(base_filename(dump_num))

        !Parse the dump group
        dump_group(dump_num) = get_group_index(dgroup)
        if(dump_group(dump_num) > group_num) then 
            write(msg, *) "Group ", dgroup, " has not been created, please select an existing group to dump"
            call command_error(msg)
        end if

        !Now if dump every = 0 set it to a very large number 
        dump_every(dump_num) = normalize_dump_every(dump_every(dump_num))

        return
    end subroutine  parse_dump

    subroutine parse_undump(line)
        !This deletes a dump command from memory
        character(len=*) :: line

        integer :: i, iospara
        character(len=read_len) :: undump_name, msg, tmp

        !Parse the command to get the undump name
        read(line, *, iostat = iospara, iomsg=msg) tmp, undump_name
        if(iospara > 0) call read_error(msg, iospara)

        !Now loop over all dump commands and delete the one that matches the name
        do i = 1, dump_num
            if(dump_name(i) == undump_name) then 
                dump_name(i:dump_num-1) = dump_name(i+1:dump_num)
                dump_every(i:dump_num-1) = dump_every(i+1:dump_num)
                dump_type(i:dump_num-1) = dump_type(i+1:dump_num)
                dump_group(i:dump_num-1) = dump_group(i+1:dump_num)
                base_filename(i:dump_num-1) = base_filename(i+1:dump_num)
                dump_num = dump_num -1
                exit
            end if
        end do
        return
    end subroutine  parse_undump

    subroutine parse_write_out(line)
        !Parse the dump command
        integer :: type_for_dump, dg
        character(len = *), intent(in) :: line
        character(len = read_len) :: label, tmptxt, msg, dtype, dgroup, fname
        integer :: iospara

        fname=''
        read(line,*, iostat=iospara) label, dgroup, dtype, fname
        if(iospara > 0) call read_error("Failure to read write_out command", iospara)
        if(fname=='') call misc_error('Missing filename for write_out command' )
        
        !Check to make sure dump file format is correct
        select case(dtype)
        case('out')
            type_for_dump=1
        case('lmp')
            type_for_dump=2
        case default
            write(tmptxt, *) "dump file format ", trim(adjustl(dtype)), " is not accepted"
            call command_error(tmptxt)
        end select

        write(msg, *) "Writing output to ", trim(adjustl(fname)), " using ", trim(adjustl(dtype)), " file format"
        call log_msg(msg)

        !Parse the dump group
        dg = get_group_index(dgroup)
        if(dg > group_num) then 
            write(msg, *) "Group ", dgroup, " has not been created, please select an existing group to dump"
            call command_error(msg)
        end if

        !Now write it out
        if(atom_num > 0) call gather_at
        if((ele_num > 0).and.(type_for_dump == 1)) call gather_cg
        if((ele_num > 0).and.(type_for_dump == 2)) call gather_atomap
        !Loop over all dumps and check to see if any need to be written
        if(rank == root) then 
            select case (type_for_dump)
            case(1) 
                call write_pycac(iter, trim(adjustl(fname)), dg)
            case(2) 
                call write_lmp(trim(adjustl(fname)), dg)
            end select
        end if

        call dealloc_gather

        return
    end subroutine  parse_write_out

    subroutine write_dump(timestep, force_dump)
        !This subroutine converts the timestep to the filename and writes the correctly formatted dump file
        integer, intent(in) :: timestep
        logical, intent(in), optional :: force_dump

        integer :: i 
        character(len = read_len) :: fname
        logical :: gathered, want_dump
        
        if(present(force_dump)) then 
            want_dump = force_dump
        else 
            want_dump = .false.
        end if

        gathered = .false.
        do i = 1, dump_num
            if(should_write_dump(timestep, dump_every(i), want_dump)) then 
                !Gather arrays if needed
                if(.not.gathered) then 
                    if(atom_num > 0) call gather_at
                    if((ele_num > 0).and.any(dump_type == 1)) call gather_cg
                    if((ele_num > 0).and.any(dump_type == 2)) call gather_atomap
                    gathered = .true.
                end if
                !Loop over all dumps and check to see if any need to be written
                write(fname, *) iter
                
                if(rank == root) then 
                    select case (dump_type(i))
                    case(1) 
                        call write_pycac(iter, trim(adjustl(base_filename(i)))//trim(adjustl(fname))//".out", dump_group(i))
                    case(2) 
                        call write_lmp(trim(adjustl(base_filename(i)))//trim(adjustl(fname))//".lmp", dump_group(i))
                    end select
                end if
            end if
        end do

        if(gathered) call dealloc_gather

    end subroutine write_dump

    subroutine dealloc_gather
        !Deallocate gather arrays when no longer needed
        if(allocated(tag_ele_gather)) then 
            deallocate(tag_ele_gather, size_ele_gather, etype_gather, basis_num_gather, basis_type_gather, e_mask_gather, &
                       r_gather, energy_gather, force_gather, virial_gather, stat=allostat) 
            if(allostat > 0) call alloc_error("Failure to deallocate cg gather arrays", allostat)
            if(need_vel) deallocate(vel_gather)
        end if

        if(allocated(tag_atom_gather)) then 
            deallocate(tag_atom_gather, type_atom_gather, a_mask_gather, r_atom_gather, energy_atom_gather, &
                       force_atom_gather, virial_atom_gather, stat = allostat)
            if(allostat > 0) call alloc_error("Failure to deallocate at gather arrays", allostat)
            if(need_vel) deallocate(vel_atom_gather)

        if(allocated(type_atomap_gather)) then 
            deallocate(r_atomap_gather, type_atomap_gather, atomap_mask_gather, stat = allostat)
            if(allostat > 0) call alloc_error("Failure to deallocate atomap gather arrays", allostat)
        end if
        end if
    end subroutine dealloc_gather

    subroutine gather_cg
        integer :: n_ints, n_reals, ntot_ints, ntot_reals, packet_num, ie, ip, ireal, iint, &
                   n, je, jp, pro, inod, i, pb_in(3,max_basisnum, ng_max_node)
        real(kind=wp) :: r_nodes(3, max_basisnum, ng_max_node), force_nodes(3, max_basisnum, ng_max_node), &
                         eng_nodes(max_basisnum, ng_max_node), virial_nodes(3,3,max_basisnum,ng_max_node),&
                         vel_nodes(3, max_basisnum, ng_max_node)
        integer, allocatable :: gather_ints(:), send_ints(:)
        real(kind=wp), allocatable :: gather_reals(:), send_reals(:)

        !Allocate gather arrays if necessary
        if(rank == root) then 
            allocate(tag_ele_gather(ele_num), size_ele_gather(ele_num), etype_gather(ele_num), e_mask_gather(ele_num), & 
                     basis_num_gather(ele_num), basis_type_gather(max_basisnum, ele_num), &
                     r_gather(3, max_basisnum, node_num), energy_gather(max_basisnum, node_num), &
                     force_gather(3, max_basisnum, node_num), virial_gather(3,3,max_basisnum,node_num), stat=allostat)
            if (allostat > 0) call alloc_error("Failure to allocate gather arrays", allostat)
            if(need_vel) allocate(vel_gather(3, max_basisnum, node_num), stat=allostat)
            if (allostat > 0) call alloc_error("Failure to allocate vel gather arrays", allostat)

        end if

        !Now allocate the send_arrays
        n_ints = 5+max_basisnum
        ntot_ints = max_send_count*n_ints

        if(need_vel) then 
            n_reals = 22*max_basisnum*ng_max_node
        else
            n_reals = 19*max_basisnum*ng_max_node
        end if
        ntot_reals = max_send_count*n_reals

        allocate(send_ints(ntot_ints), send_reals(ntot_reals), stat = allostat)
        if (allostat > 0 ) call alloc_error("Failure to allocate send arrays for cg comm in gather_in_root", allostat)
        !Only allocate gather arrays on root
        if(rank == root) then 
            allocate(gather_ints(pro_num*ntot_ints), gather_reals(pro_num*ntot_reals))
            if (allostat > 0 ) call alloc_error("Failure to allocate gather arrays for cg comm in gather_cg", allostat)
        else
            allocate(gather_ints(0), gather_reals(0))
            if (allostat > 0 ) call alloc_error("Failure to allocate gather arrays for cg comm in gather_cg", allostat)
        end if

        !Now build the local send buff array
        vel_nodes=0.0_wp
        packet_num = 0
        ip = 0 
        ie = 0
        je = 0
        jp = 0
        iint = 0
        ireal = 0
        send_ints(:) = 0
        send_reals(:) = 0.0_wp
        pb_in = 0

        if (rank == root) then 
            gather_ints(:) = 0
            gather_reals(:) = 0.0_wp
        end if
        do while (ie < ele_num) 
            if( (je == ele_num_l).or.(packet_num == max_send_count) ) then 
                !If we have sent/packed all of our elements or if this packet is equal to send_counts then we communicate
                !Gather integer data
                call mpi_gather(send_ints, ntot_ints, mpi_integer, gather_ints, ntot_ints, mpi_integer, root, world, ierr)
                !Gather real data
                call mpi_gather(send_reals, ntot_reals, mpi_wp, gather_reals, ntot_reals, mpi_wp, root, world, ierr)
                !If root then loop over all received values and unpack
                if (rank == root) then 
                    outerloop:do pro= 1, pro_num
                        innerloop:do i = 1, max_send_count
                            !Index the start point for the current element in the array
                            ireal=n_reals*((pro-1)*max_send_count + (i-1))
                            iint = n_ints*((pro-1)*max_send_count + (i-1))

                            !If this buffer in the array is 0 then we move it ahead to where the next processor beings
                            if (ie == ele_num) then 
                                exit outerloop 
                            !Exit if ie == ele_num
                            else if(gather_ints(iint+1) == 0) then 
                                exit innerloop
                            end if

                            ie = ie + 1

                            !Unpack the received buffer for the current element
                            call unpack_ele_dump(gather_ints(iint+1:iint+n_ints), gather_reals(ireal+1:ireal+n_reals), &
                                                 etype_gather(ie), tag_ele_gather(ie), size_ele_gather(ie), e_mask_gather(ie), &
                                                 basis_num_gather(ie), basis_type_gather(:,ie), r_nodes, eng_nodes, &
                                                 force_nodes, virial_nodes, vel_nodes)

                            !Set nodal quantities
                            n = ng_node(etype_gather(ie))
                            r_gather(:,:,ip+1:ip+n) = r_nodes(:,:,1:n)
                            energy_gather(:,ip+1:ip+n) = eng_nodes(:,1:n)
                            force_gather(:,:,ip+1:ip+n) = force_nodes(:,:,1:n)
                            virial_gather(:,:,:,ip+1:ip+n) = virial_nodes(:,:, :,1:n)

                            if(need_vel) vel_gather(:,:,ip+1:ip+n) = vel_nodes(:,:,1:n)

                            ip = ip + n
                        end do innerloop
                    end do outerloop
                end if

                !broadcast gathered num to see if we need to continue in the loop
                call mpi_bcast(ie, 1, mpi_int, root, world, ierr)
                
                !gathered_num can't be more than ele_num
                if(ie > ele_num) then 
                    call misc_error("Gathered_num  can't be larger than ele_num ")
                end if
                !Zero out arrays to reset for next round of communications
                packet_num = 0
                send_ints(:) = 0
                send_reals(:) = 0.0_wp
                ireal = 0
                iint = 0
            else
                !We aren't ready to communicate so keep packing  
                je = je+1   
                
                if(je > ele_num_l) call misc_error("je greater than ele_num_l in gather_cg")
                if(who_has_ele(je)) then 
                    !get the node positions
                    do inod = 1, ng_node(etype(je))
                        jp = cg_node(inod, je)
                        r_nodes(:,:,inod) = r(:,:,jp)
                        eng_nodes(:,inod) = energy_eq(:,jp)
                        force_nodes(:,:,inod) = force_eq(:, :, jp)
                        virial_nodes(:,:,:,inod) = virial_eq(:,:,:,jp)
                        if(periodic) then 
                            pb_in(:,:,inod) = pb_node(:,:,jp)
                        end if
                        if(need_vel) vel_nodes(:,:,inod) = vel(:,:,jp)
                    end do
                    
                    !Get the start position for the data in the send arrays
                    ireal=n_reals*(packet_num) 
                    iint = n_ints*(packet_num)

                    call pack_ele_dump(etype(je), tag_ele(je), size_ele(je), e_mask(je), basis_num(je), basis_type(:, je), &
                                       pb_in, r_nodes, eng_nodes, force_nodes, virial_nodes, vel_nodes,&
                                       send_ints(iint+1:iint+n_ints), send_reals(ireal+1:ireal+n_reals))
                    !Increment the packet num
                    packet_num = packet_num + 1
                end if
            end if
        end do 
        return
    end subroutine gather_cg

    subroutine gather_at
        !Gather all atoms to root processor
        integer :: ia, ja, n_ints, n_reals, ntot_ints, ntot_reals, packet_num, ireal, iint, pro, i

        integer, allocatable :: gather_ints(:), send_ints(:)
        real(kind=wp) :: v(3)
        real(kind=wp), allocatable :: gather_reals(:), send_reals(:)

        !Allocate gather arrays if necessary
        if(rank == root) then 
            allocate(tag_atom_gather(atom_num), type_atom_gather(atom_num), a_mask_gather(atom_num), r_atom_gather(3, atom_num), &
                     energy_atom_gather(atom_num), force_atom_gather(3,atom_num), virial_atom_gather(3,3,atom_num), &
                     stat = allostat)
            if (allostat > 0) call alloc_error("Failure to allocate gather at arrays", allostat)

            if(need_vel) allocate(vel_atom_gather(3,atom_num), stat = allostat)
            if (allostat > 0) call alloc_error("Failure to allocate vel at arrays", allostat)

        end if

        !Now allocate the send_arrays
        n_ints = 3
        ntot_ints = max_send_count*n_ints

        if(need_vel) then 
            n_reals=22
        else
            n_reals = 19
        end if
        ntot_reals = max_send_count*n_reals

        allocate(send_ints(ntot_ints), send_reals(ntot_reals), stat = allostat)
        if (allostat > 0 ) call alloc_error("Failure to allocate send arrays for cg comm in gather_in_root", allostat)

        !Only allocate gather arrays on root
        if(rank == root) then 
            allocate(gather_ints(pro_num*ntot_ints), gather_reals(pro_num*ntot_reals))
            if (allostat > 0 ) call alloc_error("Failure to allocate gather arrays for cg comm in gather_cg", allostat)
        else
            allocate(gather_ints(0), gather_reals(0))
            if (allostat > 0 ) call alloc_error("Failure to allocate gather arrays for cg comm in gather_cg", allostat)
        end if

        !Now build the local send buff array
        packet_num = 0
        iint = 0
        ireal = 0
        ia = 0
        ja = 0
        v(:) = 0.0_wp
        send_ints(:) = 0
        send_reals(:) = 0.0_wp

        if (rank == root) then 
            gather_ints(:) = 0
            gather_reals(:) = 0.0_wp
        end if
        do while (ia < atom_num) 
            if( (ja == atom_num_l).or.(packet_num == max_send_count) ) then 
                !If we have sent/packed all of our atoms or if this packet is equal to send_counts then we communicate
                !Gather integer data
                call mpi_gather(send_ints, ntot_ints, mpi_integer, gather_ints, ntot_ints, mpi_integer, root, world, ierr)
                !Gather real data
                call mpi_gather(send_reals, ntot_reals, mpi_wp, gather_reals, ntot_reals, mpi_wp, root, world, ierr)
                !If root then loop over all received values and unpack
                if (rank == root) then 
                    outerloop:do pro= 1, pro_num
                        innerloop:do i = 1, max_send_count
                            !Index the start point for the current element in the array
                            ireal=n_reals*((pro-1)*max_send_count + (i-1))
                            iint = n_ints*((pro-1)*max_send_count + (i-1))

                            !If this buffer in the array is 0 then we move it ahead to where the next processor beings
                            !Exit if ia == atom_num
                            if (ia == atom_num) then 
                                exit outerloop 
                            else if(gather_ints(iint+1) == 0) then 
                                exit innerloop
                            end if

                            ia = ia + 1

                            !Unpack the received buffer for the current element
                            call unpack_atom_dump(gather_ints(iint+1:iint+n_ints), gather_reals(ireal+1:ireal+n_reals), &
                                                  tag_atom_gather(ia), type_atom_gather(ia), a_mask_gather(ia), &
                                                  r_atom_gather(:,ia), energy_atom_gather(ia), force_atom_gather(:,ia),&
                                                  virial_atom_gather(:,:,ia), v(1:3))

                            if(need_vel) vel_atom_gather(:,ia) = v(1:3)

                        end do innerloop
                    end do outerloop
                end if

                !broadcast gathered num to see if we need to continue in the loop
                call mpi_bcast(ia, 1, mpi_int, root, world, ierr)
                
                !ia can't be more than atom_num 
                if(ia > atom_num) then 
                    call misc_error("Gathered_num  can't be larger than atom_num ")
                end if
                !Zero out arrays to reset for next round of communications
                packet_num = 0
                send_ints(:) = 0
                send_reals(:) = 0.0_wp
                ireal = 0
                iint = 0
            else
                !We aren't ready to communicate so keep packing  
                ja = ja+1   
                
                if(ja > atom_num_l) call misc_error("ja greater than atom_num_l in gather_cg")
                !Get the start position for the data in the send arrays
                ireal=n_reals*(packet_num) 
                iint = n_ints*(packet_num)

                if(need_vel) v(:) = vel_atom(:,ja)
                call pack_atom_dump(tag_atom(ja), type_atom(ja), a_mask(ja), r_atom(:,ja), energy_atom(ja), force_atom(:,ja), &
                                    virial_atom(:,:, ja), v, send_ints(iint+1:iint+n_ints), send_reals(ireal+1:ireal+n_reals))

                !Increment the packet num
                packet_num = packet_num + 1
            end if
        end do 
        return
    end subroutine gather_at


    subroutine gather_atomap
        !This subroutine gathers atomap information. This is required to get the correct atomistic positions for
        !NVT simulations
        integer :: ia, ja, ntot_ints, ntot_reals, packet_num, ireal, iint, pro, i, &
                   atomap_mask(atomap_num_l), ie, iatom, iatomap, virt_count, ibasis

        integer, allocatable :: gather_ints(:), send_ints(:)
        real(kind=wp), allocatable :: gather_reals(:), send_reals(:)

        !Allocate gather arrays if necessary
        ntot_ints = 2*max_send_count
        ntot_reals = 3*max_send_count

        if(rank == root) then 
            allocate(r_atomap_gather(3, atomap_num), type_atomap_gather(atomap_num), atomap_mask_gather(atomap_num), &
                     stat = allostat)
            if (allostat > 0) call alloc_error("Failure to allocate gather atomap arrays", allostat)

            allocate(gather_ints(pro_num*ntot_ints), gather_reals(pro_num*ntot_reals), stat=allostat)
            if(allostat > 0) call alloc_error("Failure to allocate gather account", allostat)
        else
            allocate(gather_ints(0), gather_reals(0), stat=allostat)
            if(allostat > 0) call alloc_error("Failure to allocate gather account", allostat)

        end if

        !Now allocate the send_arrays
        allocate(send_ints(ntot_ints), send_reals(ntot_reals), stat = allostat)
        if (allostat > 0 ) call alloc_error("Failure to allocate send arrays for cg comm in gather_in_root", allostat)

        !First get the group positions for all atomaps
        do ie = 1, ele_num_l
            select case(etype(ie))
            case(1,2,3)
                virt_count=(size_ele(ie)+1)**3
            end select
            do iatom= 1, virt_count
                do ibasis = 1, basis_num(ie)
                    iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                    if(iatomap /= 0) then 
                        atomap_mask(iatomap) = e_mask(ie) 
                    end if
                end do
            end do
        end do


        packet_num = 0
        iint = 0
        ireal = 0
        ia = 0
        ja = 0
        send_ints(:) = 0
        send_reals(:) = 0.0_wp

        if(rank == root) then 
            gather_ints(:) = 0
            gather_reals(:) = 0.0_wp
        end if
        do while (ia < atomap_num) 
            if( (ja == atomap_num_l).or.(packet_num == max_send_count) ) then 
                !If we have sent/packed all of our atoms or if this packet is equal to send_counts then we communicate
                !Gather integer data
                call mpi_gather(send_ints, ntot_ints, mpi_integer, gather_ints, ntot_ints, mpi_integer, root, world, ierr)
                !Gather real data
                call mpi_gather(send_reals, ntot_reals, mpi_wp, gather_reals, ntot_reals, mpi_wp, root, world, ierr)
                !If root then loop over all received values and unpack
                if (rank == root) then 
                    outerloop:do pro= 1, pro_num
                        innerloop:do i = 1, max_send_count
                            !Index the start point for the current element in the array
                            ireal=3*((pro-1)*max_send_count + (i-1))
                            iint =2*((pro-1)*max_send_count + (i-1))

                            !If this buffer in the array is 0 then we move it ahead to where the next processor beings
                            !Exit if ia == atom_num
                            if (ia == atomap_num) then 
                                exit outerloop 
                            else if(gather_ints(iint+1) == 0) then 
                                exit innerloop
                            end if

                            ia = ia + 1

                            !Unpack the received buffer for the current element
                            type_atomap_gather(ia) = gather_ints(iint+1)
                            atomap_mask_gather(ia) = gather_ints(iint+2)
                            r_atomap_gather(:, ia) = gather_reals(ireal+1:ireal+3)
                        end do innerloop
                    end do outerloop
                end if

                !broadcast gathered num to see if we need to continue in the loop
                call mpi_bcast(ia, 1, mpi_int, root, world, ierr)
                
                !ia can't be more than atom_num 
                if(ia > atomap_num) then 
                    call misc_error("Gathered_num  can't be larger than atom_num ")
                end if
                !Zero out arrays to reset for next round of communications
                packet_num = 0
                send_ints(:) = 0
                send_reals(:) = 0.0_wp
                ireal = 0
                iint = 0
            else
                !We aren't ready to communicate so keep packing  
                ja = ja+1   
                
                if(ja > atomap_num_l) call misc_error("ja greater than atom_num_l in gather_cg")
                !Get the start position for the data in the send arrays
                ireal=3*(packet_num) 
                iint = 2*(packet_num)

                
                send_ints(iint+1) = type_atomap(ja)
                send_ints(iint+2) = atomap_mask(ja)
                send_reals(ireal+1:ireal+3) = r_atomap(:,ja)
                !Increment the packet num
                packet_num = packet_num + 1
            end if
        end do
        return
    end subroutine gather_atomap

    subroutine write_pycac(timestep, filename, g)
        integer, intent(in) :: timestep, g
        character(len=*), intent(in) :: filename
        integer :: i, ia, ie, ip, ibasis, inod, write_atoms, write_eles
        real(kind=wp) :: v_vec(6)

        !Write pycac dump file
1 format("#This is a pycac dump file created using the quasi-static CAC code"/ &
         "Can be converted using cacmb code from https://gitlab.com/aselimov/cacmb")
2 format("Atoms: ", i16, " Elements:", i16)
3 format("Atom format: id atom_type x y z pe fx fy fz v11 v22 v33 v32 v31 v21")
13 format("Atom format: id atom_type x y z pe fx fy fz v11 v22 v33 v32 v31 v21 velx vely velz")
4 format('Element format: id num_node num_basis esize')
5 format('Node format: ip ib basis_type x y z pe fx fy fz v11 v22 v33 v32 v31 v21')
15 format('Node format: ip ib basis_type x y z pe fx fy fz v11 v22 v33 v32 v31 v21 velx vely velz')
6 format("Timestep: ", i16)
7 format("Box_bd ", 3a, 6f23.15)
        
        !First get counts of atoms/elements in write group
        !If group is all then it's just atom_num and ele_num
        if (g == 1) then 
            write_atoms = atom_num
            write_eles = ele_num
        else
            write_atoms = 0
            write_eles = 0
            do i = 1, atom_num
                if(btest(a_mask_gather(i),g)) write_atoms = write_atoms + 1
            end do
            do i =1, ele_num
                if (btest(e_mask_gather(i), g)) write_eles = write_eles + 1
            end do
        end if
        !Open file if root
        open(unit=11, file=trim(adjustl(filename)), action='write', status='replace', position='rewind')
        !Write header information 
        write(11, 1)
        write(11, 6) timestep
        write(11, 2) write_atoms, write_eles
        write(11, 7) (boundary(i:i)//" ", i = 1, 3), (box_bd(i), i = 1, 6)

        !get the atom type number into the form string
        write(11, '(a, 999(f23.15,1x))') "Masses: ", (masses(i), i = 1, natom_types)

        !First write the atom section if needed
        if (atom_num > 0) then 
            !Write atom header
            if(need_vel) then 
                write(11,13)
            else
                write(11,3)
            end if
            do ia = 1, atom_num
                if(btest(a_mask_gather(ia), g)) then 
                    !Assign virial vector
                    v_vec(1) = virial_atom_gather(1,1,ia)
                    v_vec(2) = virial_atom_gather(2,2,ia)
                    v_vec(3) = virial_atom_gather(3,3,ia)
                    v_vec(4) = virial_atom_gather(3,2,ia)
                    v_vec(5) = virial_atom_gather(3,1,ia)
                    v_vec(6) = virial_atom_gather(2,1,ia)
                    
                    if (need_vel) then 
                        write(11, *) tag_atom_gather(ia), type_atom_gather(ia), r_atom_gather(:,ia), &
                                                      energy_atom_gather(ia), force_atom_gather(:,ia), v_vec(:), &
                                                      vel_atom_gather(:,ia)
                    else
                        write(11, *) tag_atom_gather(ia), type_atom_gather(ia), r_atom_gather(:,ia), &
                                                      energy_atom_gather(ia), force_atom_gather(:,ia), v_vec(:)
                    end if
                end if
            end do
        end if

        if(ele_num > 0) then 
            !Write cg header
            write(11,4)
            if(need_vel) then 
                write(11,15)
            else
                write(11,5)
            end if
            ip=0
            do ie = 1, ele_num
                if(btest(e_mask_gather(ie), g)) then 
                    !Write element information
                    write(11,'(4i16)' ) tag_ele_gather(ie), ng_node(etype_gather(ie)), basis_num_gather(ie), size_ele_gather(ie)
                    !write node information
                    do inod = 1, ng_node(etype_gather(ie))
                        ip = ip + 1    
                        do ibasis = 1, basis_num_gather(ie)
                            v_vec(1) = virial_gather(1, 1, ibasis, ip)
                            v_vec(2) = virial_gather(2, 2, ibasis, ip)
                            v_vec(3) = virial_gather(3, 3, ibasis, ip)
                            v_vec(4) = virial_gather(3, 2, ibasis, ip)
                            v_vec(5) = virial_gather(3, 1, ibasis, ip)
                            v_vec(6) = virial_gather(2, 1, ibasis, ip)
                            if(need_vel) then 
                                write(11, *) inod, ibasis, basis_type_gather(ibasis, ie), r_gather(:,ibasis,ip), &
                                            energy_gather(ibasis,ip), force_gather(:, ibasis, ip), v_vec(:), &
                                            vel_gather(:, ibasis,ip)
                            else
                                write(11, *) inod, ibasis, basis_type_gather(ibasis, ie), r_gather(:,ibasis,ip), &
                                            energy_gather(ibasis,ip), force_gather(:, ibasis, ip), v_vec(:)
                            end if
                        end do
                    end do
                end if
            end do
        end if

        close(11)
        return
    end subroutine write_pycac

    subroutine write_lmp(filename, g)
        integer, intent(in) :: g
        character(len = *) :: filename  

        integer :: i, write_num =0, max_atom_tag 

        open(unit=11, file=trim(adjustl(filename)), action='write', status='replace', position='rewind')

        !Write comment line
        write(11, '(a)') '# lmp file made using CAC code'
        write(11, '(a)') 

        !Calculate total atom number to write
        write_num = 0
        do i = 1, atom_num
            if (btest(a_mask_gather(i), g)) write_num = write_num + 1
        end do
        do i = 1, atomap_num
            if (btest(atomap_mask_gather(i), g)) write_num = write_num + 1 
        end do

        !Write total number of atoms + interpolated atoms
        write(11, '(i16, a)') write_num, ' atoms'
        !Write number of atom types
        write(11, '(i16, a)') natom_types, ' atom types'

        write(11, '(a)') 
        !Write box boundaries
        write(11, '(2f23.15, a)') box_bd(1:2), ' xlo xhi'
        write(11, '(2f23.15, a)') box_bd(3:4), ' ylo yhi'
        write(11, '(2f23.15, a)') box_bd(5:6), ' zlo zhi'

        !Masses
        write(11, '(a)') 'Masses'

        write(11, '(a)') 
        do i = 1, natom_types
            write(11, '(i16, f23.15)') i, masses(i)
        end do
        write(11, '(a)') 
        
        !Now first write atom positions
        write(11, '(a)') 'Atoms'
        write(11, '(a)') 

        do i = 1, atom_num
            if (btest(a_mask_gather(i), g)) then 
                write(11, '(2i16, 3f23.15)') tag_atom_gather(i), type_atom_gather(i), r_atom_gather(:,i)
            end if
        end do

        max_atom_tag = maxval(tag_atom_gather)
        do i = 1, atomap_num
            if (btest(atomap_mask_gather(i), g)) then 
                write(11, '(2i16, 3f23.15)') max_atom_tag+i, type_atomap_gather(i), r_atomap_gather(:,i)
            end if
        end do

        close(11)
        return
    end subroutine write_lmp

    pure function need_dump(time)
        integer, intent(in) :: time
        logical :: need_dump
        integer :: i

        need_dump = .false.
        do i = 1, dump_num
            if(mod(time, dump_every(i)) == 0) then 
                need_dump = .true.
                return
            end if
        end do
        return
    end function need_dump
end module dump
