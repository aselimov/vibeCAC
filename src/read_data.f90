module read_data
    ! This code is in charge of all the reading and writing
    use mpi
    use math
    use parameters
    use elements
    use box
    use integration
    use potential
    use dynamics
    use comms
    use logger
    use str
    use time


    public
    contains

    subroutine parse_read(line)
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, filename
        integer :: iospara

        read(line, *, iostat = iospara) tmptxt, filename
        if(iospara > 0) call read_error(" Invalid read of read_data command",iospara)

        if(rank == root) then 
            write(tmptxt,*) "Reading data from ", trim(adjustl(filename))
            call log_msg(tmptxt)
        end if
        call read_model(filename)
    end subroutine

    subroutine read_model(filename)
        !This subroutine is in charge of reading the model and distributing the model to the different processors
        character(len=read_len), intent(in) :: filename
        character(len=read_len) :: msg
        integer :: numbers(8), i, j, k
        real(kind=wp) :: box_mat_bcast(9)

        numbers(:) = 0
        if (rank == root)  then
            !First read the file if root and broadcast numbers, currently only accepts restart files
            call read_restart(filename) 
            if (ele_num > 0) atomap_max = max_basisnum*(maxval(size_ele_scatter)+1)**3

            numbers(1) = ele_num
            numbers(2) = node_num
            numbers(3) = atom_num
            numbers(4) = atomap_num
            numbers(5) = ng_max_node
            numbers(6) = max_basisnum
            numbers(7) = atomap_max
            numbers(8) = iter
            call mpi_bcast(numbers, 8, mpi_integer, root, world, ierr)
            write(msg, *) "Read ", ele_num, " elements with ", node_num, " nodes and read ", atom_num, " atoms"
            call log_msg(msg)
        else
            !If not root then wait to recieve files and then initialize numbers 
            call mpi_bcast(numbers, 8, mpi_integer, root, world, ierr)
            ele_num     = numbers(1) 
            node_num    = numbers(2)
            atom_num    = numbers(3)
            atomap_num  = numbers(4)
            ng_max_node = numbers(5) 
            max_basisnum= numbers(6) 
            atomap_max  = numbers(7)
            iter = numbers(8)

        end if

        !Communicate masses
        call mpi_bcast(masses, max_atom_types, mpi_wp, root, world, ierr)
        !Now we initialize the local variable nums        
        ele_num_l = nint((real(ele_num,wp) / pro_num) * 1.2_wp)
        node_num_l = nint((real(node_num,wp) / pro_num) * 1.2_wp)
        atom_num_l = nint((real(atom_num,wp) / pro_num) * 1.2_wp)
        atomap_num_l = nint((real(atomap_num,wp) / pro_num) * 1.2_wp)
        ele_num_lr = ele_num_l
        node_num_lr = node_num_l
        atom_num_lr = atom_num_l
        atomap_num_lr = atomap_num_l
        !Now get the largest number of atoms per element

        !Now broadcast and record box information
        call mpi_bcast(box_bd, 6, mpi_wp, root, world, ierr)
        !Calculate box length
        do i = 1, 3
            box_length(i) = box_bd(2*i) - box_bd(2*i-1)
        end do

        !Now set the original box bounds
        orig_box_bd = box_bd
        orig_box_length = box_length

        !Now broadcast box_bcs
        call mpi_bcast(boundary, 3, mpi_character, root, world, ierr)
        !Define logical periodic boundary variables
        do i = 1,3
            if (boundary(i:i) == 'p') then 
                period(i) = .true.
            else 
                period(i) = .false.
            end if
        end do
        periodic = any(period)
        
        !Now prepare box_mat buffer and broadcast
        if(rank == root) then 
            k = 1
            do i = 1,3 
                do j = 1, 3
                    box_mat_bcast(k) = box_mat(i,j)
                    k = k+1
                end do
            end do
        end if
        call mpi_bcast(box_mat_bcast, 9, mpi_wp, root, world, ierr)
        if(rank /= root) then 
            k = 1
            do i = 1,3 
                do j = 1, 3
                    box_mat(i,j) = box_mat_bcast(k) 
                    k = k+1
                end do
            end do
        end if

        !Now initialize original box_mat 
        box_mat_ori = box_mat

        !Now communicate the velocity flag
        call mpi_bcast(need_vel, 1, mpi_logical, root, world, ierr)

        !Now scatter the models
        call divide_domain

        if(ele_num > 0) then 
            call alloc_cg_arrays
            call scatter_cg
            call init_integration
        end if

        if (atom_num > 0) then 
            call alloc_at_arrays
            call scatter_at
        end if

        !Init all group
        call init_group_all
    end subroutine read_model

    subroutine read_restart(filename)
        !First we read the restart file headers
        character(len=read_len), intent(in) :: filename
        integer :: iospara, i, j, ip, ib, inod, ibasis, bnum
        character(len=2) :: paraline
        character(len=read_len) :: line, label 
        character(len=100) :: ioerror, msg
        logical :: read_vel_cg, read_vel_at

        etype_present(:) = .false.

        !open the file for reading
        open(1, file = trim(adjustl(filename)), status='old', action='read', position='rewind')

        !Read the time
        read(1, '(a)', iostat=iospara) paraline
        read(1, *, iostat=iospara, iomsg=ioerror) iter, t
        if(iospara>0) call read_error(ioerror, iospara)

        !Read element number
        read(1, '(a)', iostat=iosline) paraline 
        read(1, *, iostat=iospara) ele_num
        if(iospara>0) call read_error("Invalid read of element_num ", iospara)

        !Read max interpolated atoms per element
        read(1, '(a)', iostat=iosline) paraline 
        read(1, *, iostat=iospara) atomap_max_ele
        if(iospara>0) call read_error("Invalid read of max interpo", iospara)

        !Read node number
        read(1, '(a)', iostat=iosline) paraline
        read(1, *, iostat=iospara) node_num
        if(iospara>0) call read_error("Invalid read of node_num ", iospara)

        !Read max nodes per element and basisnum
        read(1, '(a)', iostat=iosline) paraline
        read(1, *, iostat=iospara) ng_max_node, max_basisnum
        if(iospara>0) call read_error("Invalid read of ng_max_node/max_basisnum ", iospara)

        !Read number of atoms
        read(1, '(a)', iostat=iosline) paraline
        read(1, *, iostat=iospara) atom_num
        if(iospara>0) call read_error("Invalid read of atom_num", iospara)

        !Now read the atom type masses and map them to the right potential 
        read(1, '(a)', iostat=iosline) paraline
        read(1, '(a)', iostat=iosline) line
        if(iospara>0) call read_error("Invalid read of atom_types ", iospara)
        j = tok_count(line)
        read(line, *, iostat=iospara) (masses(i), i = 1, j)

        !Now read the box boundary definition
        read(1, '(a)', iostat=iosline) paraline
        read(1, *, iostat=iospara) boundary
        if(iospara>0) then
            print *, "Error: Invalid read of boundary with error code", iospara
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        !Read box_bd
        box_bd(:) = 0.0_wp
        read(1, '(a)', iostat=iosline) paraline
        read(1, *, iostat=iospara) box_bd(1:6)
        if(iospara>0) then
            print *, "Error: Invalid read of box_bd with error code", iospara
            print *, "Last line read was: ", paraline
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if
        
        !Read box matrix 
        box_mat(:,:)=0.0_wp
        read(1, '(a)', iostat=iosline) paraline
        do i=1,3
            read( 1, *, iostat=iospara)  box_mat(1:3, i)
            if(iospara>0) then
                print *, "Error: Invalid read of box_mat with error code", iospara
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
        end do


        read_vel_cg=.false.
        read_vel_at=.false.
        ! Now read the cg part of the restart file
        if(ele_num > 0) then 
            !Initialize some coarse-grained variables and allocate scatter arrays
            call scatter_cg_array

            atomap_num = 0
            !Read info lines
            read(1, '(a)', iostat=iospara) line
            read(1, '(a)', iostat=iospara) line

            !Figure out if velocity information is included in this restart file
            read(1, '(a)', iostat=iospara) line
            read(line, *, iostat=iospara) (label, i =1 ,7)
            if(trim(adjustl(label)) == 'velx') then
                read_vel_cg = .true.
            end if

            !Now start reading the element information
            do i=1,ele_num
                read(1, '(a)', iostat=iospara) line
                read(line, *, iostat=iospara, iomsg=msg) tag_ele_scatter(i), bnum, etype_scatter(i), size_ele_scatter(i)
                if(iospara>0) call read_error(msg, iospara)

                !Set the element_type to true
                etype_present(etype_scatter(i)) = .true.

                !Now read the nodes
                do inod = 1, ng_node(etype_scatter(i))
                    do ibasis =1 , bnum
                        if(read_vel_cg) then 
                            read(1, *, iostat=iospara) ip,ib,basis_type_scatter(ib, i),r_scatter(:, ib, ip),vel_scatter(:, ib, ip)
                            if(iospara>0) then
                                print *, "Error: Invalid read of node information with error code", iospara
                                call mpi_abort(mpi_comm_world, 1, ierr)
                            end if
                        else
                            read(1, *, iostat=iospara) ip, ib, basis_type_scatter(ib,i), r_scatter(:, ib, ip)
                            if(iospara>0) then
                                print *, "Error: Invalid read of node information with error code", iospara
                                call mpi_abort(mpi_comm_world, 1, ierr)
                            end if
                        end if
                    end do
                    basis_num_scatter(i) = bnum
                end do

                !Figure out the number of atoms needed
                select case(etype_scatter(i))
                case(1,2,3)
                    atomap_num =atomap_num +  bnum*((size_ele_scatter(i)+1)**3)
                end select
            end do
     
        end if

        !Now read the atomistic information from the restart files
        if(atom_num > 0) then 
            !First allocate the data
            call scatter_at_array
            read(1, '(a)', iostat=iospara) line

            !Figure out if velocity information is included in this restart file
            read(1, '(a)', iostat=iospara) line
            read(line, *, iostat=iospara) (label, i =1 ,6)
            if(trim(adjustl(label)) == 'velx') then
                read_vel_at = .true.
            end if

            !Now read the atomistic information
            do i =1, atom_num
                if(read_vel_at) then 
                    read(1, *, iostat=iospara) tag_atom_scatter(i), type_atom_scatter(i), &
                                               r_atom_scatter(1:3,i), vel_atom_scatter(1:3,i)
                    if(iospara>0) then
                        print *, "Error: Invalid read of atom information with error code", iospara
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if
                else
                    read(1, *, iostat=iospara) tag_atom_scatter(i), type_atom_scatter(i), r_atom_scatter(1:3,i)
                    if(iospara>0) then
                        print *, "Error: Invalid read of atom information with error code", iospara
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if
            end if           
            end do

        end if

        !Figure out if we need velocity, need_vel is only true if read_vel_cg and read_vel_at are both true
        if ((ele_num > 0).and.(atom_num > 0)) then 
            need_vel = (read_vel_cg ).and. read_vel_at
        else if(ele_num > 0) then 
            need_vel = read_vel_cg
        else if(atom_num > 0) then 
            need_vel = read_vel_at
        end if

        close(1)
        return
    end subroutine read_restart


end module read_data
