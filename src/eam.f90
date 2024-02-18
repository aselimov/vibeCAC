module eam 
    !Subroutines for eam potential style
    
    use parameters
    use math
    use forces
    use elements 
    use integration
    use neighbors
    use comms
    use errors
    use atom_types

    implicit none

    integer, save :: numpts_p, numpts_d, numpts_e, pot_numf, pot_atom_types, eammpot_map(max_atom_types), &
                     eamtype_to_pair(max_atom_types, max_atom_types), eamtype_to_dens(max_atom_types, max_atom_types), &
                     pair_numf 


    integer :: nr, nrho
    real(kind = wp), save :: drho, dr, rdrho, rdr 
    real(kind = wp), save :: edens_host_atom_sum, edens_host_atom_max, edens_host_atom_min, edens_host_intpo_sum, &
                            edens_host_intpo_max, edens_host_intpo_min, edensmax
 
    real(kind = wp), save :: p_start, p_finish, d_start, d_finish, e_start, e_finish

    real(kind = wp), allocatable, save :: edens_host_atom(:), edens_host_intpo(:, :), edens_host_atomap(:), &
                                          embed_deriv_atom(:), embed_deriv_intpo(:,:), embed_deriv_atomap(:)  
    
    real(kind = wp), allocatable, save :: zr(:)
    real(kind = wp), allocatable, save :: pair_tab(:,:),  edens_tab(:,:), embed_tab(:,:) 


    real(kind=wp), allocatable, private, save :: pair_spline(:,:,:), edens_spline(:,:,:), embed_spline(:,:,:)
    real(kind=wp), allocatable, private, save :: pair_cutoff(:)

    integer, save :: type_to_pair(max_atom_types, max_atom_types), type_to_dens(max_atom_types, max_atom_types), &
                     pot_map(max_atom_types)

    character(len=2), dimension(max_atom_types) :: pair_atom_type
    !fs
    logical :: fsflag

    public
    contains

    subroutine eam_lammps(filename, alloy_flag, fs_flag)
    ! read a DYNAMO setfl formatted EAM(/ALLOY) file as used by LAMMPS
        character(len=*), intent(in) :: filename
        integer, intent(in) :: alloy_flag
        logical, intent(in), optional :: fs_flag

        integer :: n, idx, i, j, k, ncol, nline_rho, nline_r, eloop, e, ei
        real(kind =wp) :: val, mass
        character(len=5) :: lattice
        character(len=read_len) line
        real(kind = wp), dimension(10) :: test_cols ! expected tokens per line <= 10
        logical :: fs
        
        if(present(fs_flag)) then 
            fs=fs_flag
        else 
            fs=.false.
        end if
        open(21, file = filename, status = 'old', action = 'read', position = 'rewind')
        ! check if we are reading single eam or eam/alloy to determine header settings
        n = 4 
        if (alloy_flag == 1) then
            ! setfl header
            ! skip 3 comment header lines and extract atom type info
            do i=1, 3
                read(21,*)
            end do
            read(21, *, iostat = error) pot_numf, (pair_atom_type(i), i = 1, pot_numf)
            read(21,*) numpts_e, drho, numpts_d, dr, rc_off

        else
        
            read(21,*) ! skip over first line comment info

            ! standard funcfl header
            read(21,*) val, mass, val, lattice

            ! initialize fixed vals for single element
            pot_numf = 1
            pair_numf = ((pot_numf**2 + pot_numf)/2)
        end if


        read(21,*) ! skip first header section
        
        ! determine the number of tokens per line
        read(21, '(a)', iostat = error) line
        do i = 1, 10
            read(line, *, iostat=error) test_cols(1:i)
            if (error==-1) exit
        enddo
        ncol = i - 1

        ! setup loop variables, array allocations and setup
        nline_rho = numpts_e / ncol
        nline_r = numpts_d / ncol
        numpts_p = numpts_d
        pair_numf = ((pot_numf**2 + pot_numf)/2)

        allocate( pair_tab(numpts_p+n, pair_numf), stat = allostat) 
        if(allostat /= 0) call alloc_error('failure to allocate pair_tab', allostat)

        p_start = 0
        p_finish = (numpts_d-1) * dr
        pair_tab(:,:) = 0.0_wp
        pair_tab(1,:) = real(numpts_p+n, wp)
        pair_tab(4,:) = (p_finish - p_start) / real(numpts_p-1, wp)
        pair_tab(2,:) = p_start - 5.0_wp * pair_tab(4,:)
        pair_tab(3,:) = p_finish

        if(fs) then
            allocate(edens_tab(numpts_d+n,pot_numf*pot_numf),  stat = allostat)
            if(allostat /= 0) call alloc_error('failure to allocate edens_tab', allostat)
        else
            allocate(edens_tab(numpts_d+n,pot_numf), stat = allostat)
            if(allostat /= 0) call alloc_error('failure to allocate edens_tab', allostat)
        end if
        d_start = 0
        d_finish = (numpts_d-1) * dr
        edens_tab(:,:) = 0.0_wp
        edens_tab(1,:) = real(numpts_d+n, wp)
        edens_tab(4,:) = (d_finish-d_start)/real(numpts_d-1, wp)
        edens_tab(2,:) = d_start-5.0_wp*edens_tab(4,:)
        edens_tab(3,:) = d_finish

        allocate( embed_tab(numpts_e+n, pot_numf),  stat = allostat )
        if(allostat /= 0) call alloc_error("failure to allocate embed_tab", allostat)
        e_start = 0
        e_finish = (numpts_e-1) * drho
        embed_tab(:,:) = 0.0_wp
        embed_tab(1,:) = real(numpts_e+n, wp)
        embed_tab(4,:) = (e_finish-e_start)/real(numpts_e-1, wp)
        embed_tab(2,:) = e_start-5.0_wp*embed_tab(4,:)
        embed_tab(3,:) = e_finish
        edensmax=e_finish

        ! reset to start of embedding data (two lines for alloy)
        backspace 21 
        if (alloy_flag == 1) backspace 21

        if(fs) then 
            eloop = pot_numf
        else
            eloop = 1
        end if

        ei = 1
        do i = 1, pot_numf
            ! read header info
            read(21,*, iostat = error) val, val, val, lattice
            if (error == -1) call read_error("failed to read element header", allostat)
            ! read embedding data
            do j = 1, nline_rho
                idx = (j-1) * ncol + 1 + n
                read(21, '(a)', iostat = error) line
                if (ncol /= 1) then
                    read(line, *, iostat = error) (embed_tab(k, i), k = idx, idx + ncol -1)
                else
                    read(line, *, iostat = error) embed_tab(idx, i)
                endif
            end do

            ! read electron density data
            do e = 1, eloop
                do j = 1, nline_r
                    idx = (j-1) * ncol + 1 + n
                    read(21, '(a)', iostat = error) line
                    if (ncol /= 1) then
                        read(line, *, iostat = error) (edens_tab(k, ei), k = idx, idx + ncol -1)
                    else
                        read(line, *, iostat = error) edens_tab(idx,ei)
                    endif
                end do
                ei = ei + 1
            end do
        end do

        do i = 1, pair_numf
            do j = 1, nline_r
                idx = (j-1) * ncol + 1 + n
                read(21, '(a)', iostat = error) line
                if (ncol /= 1) then
                    read(line, *, iostat = error) (pair_tab(k, i), k = idx, idx + ncol -1)
                else
                    read(line, *, iostat = error) pair_tab(idx,i)
                endif
            end do

        end do
        
        !Setup permutation lists for pair potential
        eamtype_to_pair(:, :) = 0 
        k  = 1
        do i = 1, pot_numf
            do j = 1, pot_numf
                if (i >= j) then 
                    eamtype_to_pair(i,j) = k
                    eamtype_to_pair(j,i) = k
                    k = k + 1 
                else 
                    exit
                end if
            end do
        end do

        eamtype_to_dens(:,:) = 0 
        if (fs) then 
            fsflag = .true.
            ei = 1
            do i = 1, pot_numf
                do j = 1, pot_numf
                    eamtype_to_dens(i,j) = ei
                    ei = ei + 1
                end do
            end do
        else
            fsflag = .false.
            do i = 1, pot_numf
                eamtype_to_dens(i,:) = i
            end do
        end if

        close(21)
        p_finish = min(rc_off, p_finish)

    end subroutine eam_lammps
    
        subroutine alloc_eam_arrays
        !Allocate eam arrays
        if(ele_num > 0) then 
            if(allocated(edens_host_atomap)) then 
                deallocate(edens_host_atomap, edens_host_intpo, embed_deriv_atomap, embed_deriv_intpo, stat = allostat)
                if(allostat>0) call alloc_error("Failure to deallocate cg_eam arrays", allostat)
            end if
            allocate(edens_host_atomap(atomap_num_lg), edens_host_intpo(max_basisnum*max_intpo_num, ele_num_l), &
                     embed_deriv_atomap(atomap_num_lg), embed_deriv_intpo(max_basisnum*max_intpo_num, ele_num_l), stat = allostat)
            if(allostat > 0) call alloc_error("Failure to alloc cg_eam arrays", allostat)
        end if

        if(atom_num > 0 ) then 
            if(allocated(edens_host_atom)) then 
                deallocate(edens_host_atom, embed_deriv_atom, stat = allostat)
                if(allostat > 0) call alloc_error("Failure to alloc at_eam arrays", allostat)
            end if
            allocate(edens_host_atom(atom_num_lg), embed_deriv_atom(atom_num_lg), stat=allostat)
            if(allostat > 0) call alloc_error("Failure to alloc at_eam arrays", allostat)
        end if
    end subroutine alloc_eam_arrays

    subroutine comm_edens_intpo
        !Now communicate edens for all integration points to make sure every processor has the host density for 
        !all integration points in the elements that it is in charge of
        integer :: ie, je, iep, jep, ibasis
        real(kind=wp), dimension(max_intpo_num*ele_shared_num*max_basisnum) :: edens_host_array, edens_host_buff

        edens_host_array(:) = 0.0_wp
        do ie = 1, ele_num_l
            je = ele_id_shared(ie)
            
            if(je /= 0) then 
                do iep = 1, intpo_count(itype(ie))
                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1) + ibasis

                        if(who_has_intpo(jep, ie)) then 
                             
                            !First check to make sure the index fits within the element
                            if(max_intpo_num*max_basisnum*(je-1)+ jep > max_intpo_num*max_basisnum*ele_shared_num) then
                                print *, 'Error: Index of edens_host_array', &
                                         max_intpo_num*max_basisnum*(je-1)+ jep, ' is larger than', &
                                        ' array size', max_intpo_num*max_basisnum*ele_shared_num
                                call mpi_abort(mpi_comm_world, 1, ierr)
                            end if

                            edens_host_array(max_intpo_num*max_basisnum*(je-1)+jep) = edens_host_intpo(jep, ie)
                        end if
                    end do
                end do
            end if
        end do
        !Gather all edens host by summing them into edens_host_buff
        edens_host_buff(:) = 0.0_wp 
        call mpi_allreduce(edens_host_array, edens_host_buff, max_intpo_num*ele_shared_num*max_basisnum, &
                           mpi_wp, mpi_sum, world, ierr) 
        !Now set edenshost for intpo that we don't own
        do ie = 1, ele_num_l
            je = ele_id_shared(ie)
            if(je /= 0) then 
                do iep = 1, intpo_count(itype(ie))
                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1) + ibasis
                        if(.not.who_has_intpo(jep, ie)) then 
                            edens_host_intpo(jep, ie) = edens_host_buff(max_intpo_num*max_basisnum*(je-1)+jep)
                        end if
                    end do
                end do
            end if
        end do
    end subroutine comm_edens_intpo

    subroutine comm_edens_ghost_cg
        !Now communicate edenhost for all ghost atomaps

        integer :: j, k, send_n, recv_n, iatomap, jatomap, latomap, send_n_max, &
                   ireq_e, ireq_t, ireq_n, send_rank, recv_rank, atomap_num_r
        logical :: send_l, recv_l
        integer, dimension(3) :: send_coords, recv_coords
        integer, dimension(mpi_status_size) :: mstatus
        integer, allocatable :: tag_send_buff(:), tag_recv_buff(:)
        real(kind = wp), allocatable :: ed_send_buff(:), ed_recv_buff(:), em_send_buff(:), em_recv_buff(:)

        send_n_max = maxval(list_atomap_num)
        allocate(tag_send_buff(send_n_max), tag_recv_buff(send_n_max), &
                 ed_send_buff(send_n_max), ed_recv_buff(send_n_max), &
                 em_send_buff(send_n_max), em_recv_buff(send_n_max), stat = allostat)

        if(allostat > 0) call alloc_error("Failure to allocate tag/ed_send/recv_buff", allostat)
        tag_send_buff(:) = 0
        tag_recv_buff(:) = 0
        ed_send_buff(:) = 0
        ed_recv_buff(:) = 0
        em_send_buff(:) = 0
        em_recv_buff(:) = 0
        send_rank = grank
        recv_rank = grank
        atomap_num_r = send_n_max
        jatomap = atomap_num_l

        !Loop over all of the send directions
        do j = 1, 3
            do k = 1, 2
                send_coords(:) = grid_coords(:)
                recv_coords(:) = grid_coords(:)
                send_l = .true.
                recv_l = .true.

                if(num_pro(j) == 1) then
                    send_l = .false.
                    !IF period than we may need to match some ghost atoms to local atoms
                    if(period(j).eqv..true.) then
                        recv_l = .true.
                        send_n = list_atomap_num(k, j)
                        recv_n = send_n

                        !Resize if needed
                        if(recv_n > atomap_num_r) then
                            deallocate(tag_recv_buff, ed_recv_buff, stat = deallostat )
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/ed_recv_buff", deallostat) 
                            allocate(tag_recv_buff(recv_n), ed_recv_buff(recv_n), stat = allostat) 
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/ed_recv_buff", allostat)
                            tag_recv_buff(:) = 0
                            ed_recv_buff(:) = 0.0_wp
                            atomap_num_r = recv_n
                        end if

                        !Add the needed information to the receive buffs
                        do iatomap = 1, send_n
                            latomap = list_atomap(iatomap, k, j)
                            tag_recv_buff(iatomap) = tag_atomap(latomap)
                            ed_recv_buff(iatomap) = edens_host_atomap(latomap)
                            em_recv_buff(iatomap) = embed_deriv_atomap(latomap)
                        end do
                    else
                        !If not periodic then we don't need to do anything
                        recv_l = .false.
                        recv_n = 0
                    end if
                !If more than one proc per side than we have to send it
                else
                    send_coords(j) = grid_coords(j) - (-1)**k
                    recv_coords(j) = grid_coords(j) + (-1)**k
                    if(period(j).eqv..true.) then
                        call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                        call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                    else
                        !If not periodic then we have to make sure that the first/last processors don't send/recv when unneeded
                        if(k == 1) then
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if

                        else if(k == 2) then
                            if(grid_coords(j) > 0) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                              send_l = .false.
                              send_rank = -1
                            end if
                            if(grid_coords(j) < num_pro(j) - 1) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if
                        end if 
                    end if

                    !check send_l
                    if(send_l.neqv.list_atomap_logic(k, j)) then
                        print *, 'Error: send_l', send_l, ' should equal list_atomap_logic', &
                                  list_atomap_logic(k, j), ' for k', k, ' and j', j
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if

                    !prepare send_buff
                    if(send_l.eqv..true.) then
                        send_n = list_atomap_num(k, j)
                        do iatomap = 1, send_n
                            latomap = list_atomap(iatomap, k, j)
                            tag_send_buff(iatomap) = tag_atomap(latomap)
                            ed_send_buff(iatomap) = edens_host_atomap(latomap)
                            em_send_buff(iatomap) = embed_deriv_atomap(latomap)
                        end do
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, 1, grid_comm, ireq_n, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(send_n, 1, mpi_integer, send_rank, 1, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                      call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !grow recv_buff if needed
                    if(recv_l.eqv..true.) then
                        if(recv_n > atomap_num_r) then
                            deallocate(tag_recv_buff, ed_recv_buff, em_recv_buff, stat = deallostat )
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/ed_recv_buff",deallostat)

                            allocate(tag_recv_buff(recv_n), ed_recv_buff(recv_n), em_recv_buff(recv_n), stat = allostat) 
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/ed_recv_buff", allostat)

                            tag_recv_buff(:) = 0
                            ed_recv_buff(:) = 0.0_wp
                            em_recv_buff(:) = 0.0_wp
                            atomap_num_r = recv_n
                        end if
                    end if

                    !send/recv tag and ed
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(tag_recv_buff, recv_n, mpi_integer, recv_rank, 2, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff, send_n, mpi_integer, send_rank, 2, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(ed_recv_buff, recv_n, mpi_wp, recv_rank, 3, grid_comm, ireq_e, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(ed_send_buff, send_n, mpi_wp, send_rank, 3, grid_comm, ierr) 
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_e, mstatus, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(em_recv_buff, recv_n, mpi_wp, recv_rank, 3, grid_comm, ireq_e, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(em_send_buff, send_n, mpi_wp, send_rank, 3, grid_comm, ierr) 
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_e, mstatus, ierr)
                    end if
                end if

                !check to make sure the tags match for ghost atoms
                if(recv_l.eqv..true.) then
                    do iatomap = 1, recv_n
                        jatomap = jatomap + 1
                        if(tag_recv_buff(iatomap) /= tag_atomap(jatomap)) then
                            print *, 'Error: Rank', rank, ' received atomap tag', &
                                      tag_recv_buff(iatomap), ' of iatomap', iatomap, &
                                     ' from rank', recv_rank, ' which should equal tag_atomap', &
                                      tag_atomap(jatomap), ' of jatomap', jatomap, &
                                     ' when j is', j, ' and k is', k, ' in edenshost'
                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if
                        edens_host_atomap(jatomap) = ed_recv_buff(iatomap)
                        embed_deriv_atomap(jatomap) = em_recv_buff(iatomap)
                    end do
                end if
            end do
        end do

        !debug

        if(jatomap /= atomap_num_lg) then
            print *, 'Error: Wrong jatomap', jatomap, &
                     ' which should equal atomap_num_lg', atomap_num_lg
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if
    end subroutine comm_edens_ghost_cg

    subroutine comm_edens_ghost_at
        
        integer :: j, k, send_n, recv_n, ia, ja, la, send_n_max, &
                   ireq_e, ireq_t, ireq_n, send_rank, recv_rank, atom_num_r, ireq_ee
        logical :: send_l, recv_l
        integer, dimension(3) :: send_coords, recv_coords
        integer, dimension(mpi_status_size) :: mstatus
        integer, allocatable :: tag_send_buff(:), tag_recv_buff(:)
        real(kind = wp), allocatable :: ed_send_buff(:), ed_recv_buff(:), em_send_buff(:), em_recv_buff(:)

        send_n_max = maxval(list_atom_num)
        allocate( tag_send_buff(send_n_max), tag_recv_buff(send_n_max), &
                  ed_send_buff(send_n_max), ed_recv_buff(send_n_max), &
                  em_send_buff(send_n_max), em_recv_buff(send_n_max), stat = allostat )
        if(allostat /= 0) call alloc_error("Failure to allocate tag/ed_send/recv_buff", allostat)

        !Initialize some variables
        tag_send_buff(:) = 0
        tag_recv_buff(:) = 0
        ed_send_buff(:) = 0.0_wp
        ed_recv_buff(:) = 0.0_wp
        em_send_buff = 0.0_wp
        em_recv_buff = 0.0_wp
        ireq_t = 0
        ireq_e = 0
        ireq_n = 0
        ireq_ee = 0
        send_rank = grank
        recv_rank = grank
        atom_num_r = send_n_max
        ja = atom_num_l

        ! k == 1: positive
        ! k == 2: negative
        !Do same looping strategy as in ghost code, this code is almost exactly the same 
        !as the code in comm_edens_ghost_cg so check that code or the ghost code itself for 
        !more indepth description
        do j = 1, 3
            do k = 1, 2

                send_coords(:) = grid_coords(:)
                recv_coords(:) = grid_coords(:)
                send_l = .true.
                recv_l = .true.

                if(num_pro(j) == 1) then

                    send_l = .false.
                    if(period(j).eqv..true.) then
                        recv_l = .true.
                        send_n = list_atom_num(k, j)
                        recv_n = send_n

                        if(recv_n > atom_num_r) then
                            deallocate(tag_recv_buff, ed_recv_buff, em_recv_buff, stat = deallostat)

                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/ed_recv_buff", deallostat)

                            allocate(tag_recv_buff(recv_n), ed_recv_buff(recv_n), em_recv_buff(recv_n), stat = allostat) 
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/ed_recv_buff", allostat)

                            tag_recv_buff(:) = 0
                            ed_recv_buff(:) = 0.0_wp
                            em_recv_buff(:) = 0.0_wp
                            atom_num_r = recv_n

                        end if
                        !Build send_array
                        do ia = 1, send_n
                            la = list_atom(ia, k, j)
                            tag_recv_buff(ia) = tag_atom(la)
                            ed_recv_buff(ia) = edens_host_atom(la)
                            em_recv_buff(ia) = embed_deriv_atom(la)
                        end do

                    else
                        recv_l = .false.
                        recv_n = 0
                    end if
                else
                    send_coords(j) = grid_coords(j) - (-1) ** k
                    recv_coords(j) = grid_coords(j) + (-1) ** k

                    if(period(j).eqv..true.) then
                        call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                        call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                    else
                        if(k == 1) then
                            if(grid_coords(j) < num_pro(j) - 1) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if
                            if(grid_coords(j) > 0) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if

                        else if(k == 2) then
                            if(grid_coords(j) > 0) then 
                                call mpi_cart_rank(grid_comm, send_coords, send_rank, ierr)
                            else
                                send_l = .false.
                                send_rank = -1
                            end if

                            if(grid_coords(j) < num_pro(j) - 1) then
                                call mpi_cart_rank(grid_comm, recv_coords, recv_rank, ierr)
                            else
                                recv_l = .false.
                                recv_rank = -1
                            end if
                        end if 
                    end if
                    !check send_l
                    if(send_l.neqv.list_atom_logic(k, j)) then
                        print *, 'Error: Send_l', send_l, ' should equal list_atom_logic', &
                                  list_atom_logic(k, j), ' for k', k, ' and j', j
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if

                    !prepare send_buff
                    if(send_l.eqv..true.) then
                        send_n = list_atom_num(k, j)
                        do ia = 1, send_n
                            la = list_atom(ia, k, j)
                            tag_send_buff(ia) = tag_atom(la)
                            ed_send_buff(ia) = edens_host_atom(la)
                            em_send_buff(ia) = embed_deriv_atom(la)
                        end do
                    end if

                    !send/recv number
                    if(recv_l.eqv..true.) then
                        call mpi_irecv(recv_n, 1, mpi_integer, recv_rank, 0, grid_comm, ireq_n, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(send_n, 1, mpi_integer, send_rank, 0, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_n, mstatus, ierr)
                    end if

                    !update recv_buff
                    if(recv_l.eqv..true.) then
                        if(recv_n > atom_num_r) then
                            deallocate( tag_recv_buff, ed_recv_buff, em_recv_buff, stat = deallostat)
                            if(deallostat /= 0) call alloc_error("Failure to deallocate tag/ed_recv_buff", deallostat)

                            allocate(tag_recv_buff(recv_n), ed_recv_buff(recv_n), em_recv_buff(recv_n), stat = allostat) 
                            if(allostat /= 0) call alloc_error("Failure to allocate tag/ed_recv_buff", allostat)

                            tag_recv_buff(:) = 0
                            ed_recv_buff(:) = 0.0_wp
                            em_recv_buff(:) = 0.0_wp
                            atom_num_r = recv_n
                        end if
                    end if
                    !send/recv tag and ed

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(tag_recv_buff, recv_n, mpi_integer, recv_rank, 1, grid_comm, ireq_t, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(tag_send_buff, send_n, mpi_integer, send_rank, 1, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_t, mstatus, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(ed_recv_buff, recv_n, mpi_wp, recv_rank, 2, grid_comm, ireq_e, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(ed_send_buff, send_n, mpi_wp, send_rank, 2, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_e, mstatus, ierr)
                    end if

                    if(recv_l.eqv..true.) then
                        call mpi_irecv(em_recv_buff, recv_n, mpi_wp, recv_rank, 3, grid_comm, ireq_ee, ierr)
                    end if
                    if(send_l.eqv..true.) then
                        call mpi_send(em_send_buff, send_n, mpi_wp, send_rank, 3, grid_comm, ierr)
                    end if
                    if(recv_l.eqv..true.) then
                        call mpi_wait(ireq_ee, mstatus, ierr)
                    end if

                end if

                !check if the tags match
                if(recv_l.eqv..true.) then
                    do ia = 1, recv_n
                        ja = ja + 1
                        if(tag_recv_buff(ia) /= tag_atom(ja)) then
                            print *, 'Error: Rank', rank, ' received atom tag', &
                                      tag_recv_buff(ia), ' of ia', ia, &
                                     ' from rank', recv_rank, ' which should equal tag_atom', &
                                      tag_atom(ja), ' of ja', ja, &
                                     ' when j is', j, ' and k is', k, ' in edenshost'

                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if
                        edens_host_atom(ja) = ed_recv_buff(ia)
                        embed_deriv_atom(ja) = em_recv_buff(ia)
                    end do
                end if
            end do
        end do
        !debug
        if(ja /= atom_num_lg) then
            print *, 'Error: Wrong ja', ja, &
                     ' which should equal atom_num_lg', atom_num_lg

            call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        return
    end subroutine comm_edens_ghost_at

    subroutine eamspline_arrays(nd, ne, np, dn, en, pn )
        !This subroutine allocates the spline arrays                                                           
        integer, intent(in) :: nd, ne, np, dn, en, pn

        integer :: i, j, k, l, m, n
        real(kind = wp), allocatable :: ptmp(:,:,:), etmp(:,:,:), dtmp(:,:,:), ctmp(:)

        if (allocated(pair_spline)) then 
            i = size(pair_spline(1,1,:))
            j = size(embed_spline(1,1,:))
            k = size(edens_spline(1,1,:))

            l = max(np, size(pair_spline(1,:,1)))
            m = max(ne, size(embed_spline(1,:,1))) 
            n = max(nd, size(edens_spline(1,:,1))) 

            allocate(ptmp(7, l, i+pn), ctmp(i+pn), etmp(7, m, j + en), dtmp(7, n, k + dn))
            ptmp = 0.0_wp
            ctmp = 0
            etmp = 0.0_wp
            dtmp = 0.0_wp

            !Move allocation 
            ptmp(:, 1:size(pair_spline(1,:,1)), 1:size(pair_spline(1,1,:))) = pair_spline
            call move_alloc(ptmp, pair_spline)

            ctmp(1:size(pair_cutoff)) = pair_cutoff
            call move_alloc(ctmp, pair_cutoff)

            etmp(:, 1:size(embed_spline(1,:,1)), 1:size(embed_spline(1,1,:))) = embed_spline
            call move_alloc(etmp, embed_spline)

            dtmp(:, 1:size(edens_spline(1,:,1)), 1:size(edens_spline(1,1,:))) = edens_spline
            call move_alloc(dtmp, edens_spline)
        else
            allocate(pair_spline(7, np, pn), pair_cutoff(pn), embed_spline(7, ne, en), edens_spline(7, nd, dn))
            pair_spline = 0.0_wp
            embed_spline = 0.0_wp
            edens_spline = 0.0_wp
            pair_cutoff = 0.0_wp

        end if
    end subroutine eamspline_arrays


    subroutine eamarray2spline
        !This assumes that all of the potential arrays are equally spaced, get the calculation coefficients for 
        !the splines
        integer :: i, eloop, pstart, estart, dstart

        if(allocated(pair_spline)) then 
            pstart = size(pair_spline(1,1,:))
            estart = size(embed_spline(1,1,:))
            dstart = size(edens_spline(1,1,:))
        else 
            pstart = 0
            estart = 0
            dstart = 0
        end if

        if(fsflag) then 
            eloop = pot_numf*pot_numf
        else 
            eloop = pot_numf
        end if

        call eamspline_arrays(numpts_d, numpts_e, numpts_p, eloop, pot_numf, pair_numf)

        
        rdrho  = 1.0/embed_tab(4,1)
        rdr = 1.0/edens_tab(4,1)
        nr = numpts_d
        nrho = numpts_e

        do i =1, pot_numf
            call interpolate(numpts_e, embed_tab(4,i), embed_tab(5:,i), embed_spline(:,:,estart+i)) 
        end do

        do i=1, eloop
            call interpolate(numpts_d, edens_tab(4,i), edens_tab(5:,i), edens_spline(:,:,dstart+i)) 
       end do

        if (edens_tab(4,1) /= pair_tab(4,1)) then 
            call misc_error("Pair tab spacing must be the same as the electron density tab spacing")
        end if

        do i = 1, pair_numf
            call interpolate(numpts_p, pair_tab(4,i), pair_tab(5:,i), pair_spline(:,:,pstart+i)) 
        end do


        return
    end subroutine eamarray2spline
    
    subroutine set_eam_map_arrays
        character(len=read_len) :: tmptxt
        integer :: i, j
        !Now map the type arrays from the eam code
        pot_map=0
        if(atom_types_set) then 
            typeloop: do i = 1, natom_types
                do j = 1, pot_numf
                    if (atom_names(i) == trim(adjustl(pair_atom_type(j)))) then 
                        pot_map(i) = j
                        cycle typeloop
                    end if
                end do
                !If we don't get a match than log a warning
                write(tmptxt,*) "Warning: ", atom_names(i), " not defined in eam potential file"
                call log_msg(tmptxt)
            end do typeloop

            !Now map the pair types and density types to the real ones
            do i = 1, natom_types
                do j = 1, natom_types
                    if((pot_map(i) > 0).and.(pot_map(j) > 0)) then 
                        type_to_pair(i,j) = eamtype_to_pair(pot_map(i), pot_map(j))
                        types_to_pot_type(i,j) = ibset(types_to_pot_type(i,j), 1)
                        pair_cutoff(type_to_pair(i,j)) = rc_off
                    end if
                end do
            end do

            !Now do the same thing to type to the edens maps
            do i = 1, natom_types
                do j = 1, natom_types
                    if((pot_map(i) > 0).and.(pot_map(j) > 0)) then 
                        type_to_dens(i,j) = eamtype_to_dens(pot_map(i), pot_map(j))
                    end if
                end do
            end do
            
        !If we don't have atom types defined then we can set them based on the elements within the potential file  
        else
                call set_atom_types(pot_numf, pair_atom_type)
                do i = 1, natom_types
                    pot_map(i) = i
                end do
                type_to_pair = eamtype_to_pair
                type_to_dens = eamtype_to_dens
                
        end if


    end subroutine set_eam_map_arrays

    subroutine update_force_eam
        !This subroutine updates the force for all atoms and elements

        integer :: i,ie, iep, jep, ia, ja, iatom, iatomap, katomtype, latomtype, iatom_counts, ibasis, nei, &
                   ic, inod, ip, m, num_nei, en_list
        real(kind=wp) :: rl(3), rk(3), flk(3), rlk(4), rsq, edens, edens_l, edens_k, f_intpo,d_embed_l, d_embed_k, &  
                         eshape, f_atom, talliesme(3), recip, d_edens_l, &
                     d_edens_k, z2, z2p, phi, phip, psip, p, energy_sum, sum_e
        real(kind=wp), dimension(max_basisnum*max_intpo_num) ::  energy_intpo
        real(kind=wp) :: embed_atom, coef(7)
        real(kind=wp), dimension(max_basisnum*max_intpo_num, ele_num_l) :: pair_intpo, embed_intpo
        real(kind=wp), dimension(3, max_basisnum*max_intpo_num, ele_num_l) :: force_intpo
        real(kind=wp), dimension(3,3, max_basisnum*max_intpo_num, ele_num_l) :: virial_intpo
        real(kind=wp), dimension(size(pair_cutoff)) :: pair_cutoffsq
        
                        
        !eam is potential type one so figure out which neighbor list is potential type one
        en_list=1 !en_list is always one because eam is the first pair style defined
        num_nei = 0
        pair_cutoffsq = pair_cutoff*pair_cutoff
        if(.not.allocated(edens_host_atomap)) then 
            call alloc_eam_arrays
        else if((atomap_num_lg > size(edens_host_atomap))) then
            call alloc_eam_arrays
        else if((ele_num_l > size(edens_host_intpo,2))) then
            call alloc_eam_arrays
        else if(.not.allocated(edens_host_atom)) then 
            call alloc_eam_arrays
        else if((atom_num_lg > size(edens_host_atom))) then 
            call alloc_eam_arrays
        end if
        coef(:) = 0.0_wp

        if(ele_num > 0) then
            edens_host_intpo = 0.0_wp
            edens_host_atomap = 0.0_wp
        end if

        !First we update electron density for the atoms
        if(atom_num > 0) then 
            edens_host_atom = 0.0_wp
            do ia = 1, atom_num_l
                do i = 1, 3
                    rl(i) = r_atom(i, ia)
                end do
                
                !Loop over all neighbor atoms
                do nei = 1, n_at_at(ia, en_list)

                    ja = at_nei_at(nei, ia, en_list)
                    !Use symmetry of electron density
                    do i = 1, 3
                        rk(i) = r_atom(i,ja)
                    end do
                    katomtype = type_atom(ja)
                    !Update edens host for current atom
                    do i = 1, 3
                    rlk(i) = rk(i) - rl(i)
                    end do
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)

                    if(rsq < pair_cutoffsq(type_to_pair(katomtype, type_atom(ia)))) then 
                        rlk(4) = sqrt(rsq)
                        p = rlk(4)*rdr+1.0_wp
                        m = int(p)
                        m = min(m, nr-1)
                        p = p - m 
                        p = min(p, 1.0)
                        coef = edens_spline(:,m,type_to_dens(katomtype, type_atom(ia)))
                        edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                        edens_host_atom(ia) = edens_host_atom(ia) + edens

                        if(ja <= atom_num_l) then
                            coef = edens_spline(:,m,type_to_dens(type_atom(ia), katomtype))
                            edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                            edens_host_atom(ja) = edens_host_atom(ja) + edens
                        end if
                    end if
                end do
                !Loop over all neighbor atomaps
                do nei = 1, n_at_cg(ia, en_list)

                    ja = at_nei_cg(nei, ia, en_list)
                    katomtype = type_atomap(ja)
                    !Update edens host for current atom
                    rlk(1:3) = r_atomap(:,ja) - rl
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    
                    if(rsq < pair_cutoffsq(type_to_pair(katomtype, type_atom(ia)))) then 
                        rlk(4) = sqrt(rsq)
                        p = rlk(4)*rdr+1.0_wp
                        m = int(p)
                        m = min(m, nr-1)
                        p = p - m 
                        p = min(p, 1.0)
                        coef = edens_spline(:,m,type_to_dens(katomtype, type_atom(ia)))
                        edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                        edens_host_atom(ia) = edens_host_atom(ia) + edens

                        if(ja <=atomap_num_l) then
                            if(atomap_to_intpo(1,ja) /= 0) then 
                                coef = edens_spline(:,m,type_to_dens(type_atom(ia), katomtype))
                                edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                jep = atomap_to_intpo(1,ja)
                                ie = atomap_to_intpo(2,ja)
                                edens_host_intpo(jep, ie) = edens_host_intpo(jep, ie) + edens
                                edens_host_atomap(ja) = edens_host_atomap(ja) + edens
                            else if (needed_atomap(ja)) then 
                                coef = edens_spline(:,m,type_to_dens(type_atom(ia), katomtype))
                                edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                edens_host_atomap(ja) = edens_host_atomap(ja) + edens
                            end if
                        end if
                    end if
                end do

            end do
            !Now calculate all embedding energies and embedding energy derivatives
            embed_deriv_atom = 0.0_wp
            do ia = 1, atom_num_l
                if(pot_map(type_atom(ia)) > 0) then 
                    !First calculate the embedding energy contribution to the energy
                    p = edens_host_atom(ia)*rdrho+1.0_wp
                    m = int(p)
                    m = max(1, min(m, nrho-1))
                    p = p -m
                    p = min(p, 1.0)
                    coef = embed_spline(:,m,pot_map(type_atom(ia)))
                    embed_atom = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                    d_embed_l = (coef(1)*p + coef(2))*p + coef(3)
                    !Extra linear term to conserve energy when atoms are too close
                    if(edens_host_atom(ia) > edensmax) then 
                        embed_atom = embed_atom + d_embed_l*(edens_host_atom(ia) - edensmax)
                    end if
                    energy_atom(ia) = energy_atom(ia) + embed_atom
                    embed_deriv_atom(ia) = d_embed_l
                end if
            end do

            !Now communicate ghost atom edens and embed_deriv
            call comm_edens_ghost_at

        end if


        !Now we have to update the electron host densities for all integration points
        if (ele_num > 0) then 
            do ie = 1, ele_num_l
                do iep = 1, intpo_count(itype(ie))
                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1) + ibasis
                        if(who_has_intpo(jep, ie)) then

                            iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                            iatomap = cg_atomap((basis_num(ie) * (iatom-1)) + ibasis, ie)
                            do i = 1, 3
                                rl(i) = r_atomap(i, iatomap)
                            end do

                            !Loop over all atomaps
                            do nei = 1, n_cg_cg(jep, ie, en_list)
                                ja = cg_nei_cg(nei, jep, ie, en_list)
                                katomtype = type_atomap(ja)
                                rk = r_atomap(:,ja)
                                do i = 1, 3
                                    rlk(i) = rk(i) - rl(i)
                                end do
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)

                                if(rsq < pair_cutoffsq(type_to_pair(basis_type(ibasis,ie), katomtype))) then 
                                    rlk(4) = sqrt(rsq)
                                    p = rlk(4)*rdr+1.0_wp
                                    m = int(p)
                                    m = min(m, nr-1)
                                    p = p - m 
                                    p = min(p, 1.0)
                                    coef = edens_spline(:,m,type_to_dens(katomtype, basis_type(ibasis, ie)))
                                    edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                    edens_host_intpo(jep, ie) = edens_host_intpo(jep, ie) + edens
                                end if
                            end do

                            !Loop over all ghost atoms
                            do nei = 1, n_cg_at(jep,ie,en_list)
                                ja = cg_nei_at(nei, jep, ie, en_list)
                                katomtype = type_atom(ja)
                                rk = r_atom(:,ja)
                                do i = 1, 3
                                    rlk(i) = rk(i) - rl(i)
                                end do
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < pair_cutoffsq(type_to_pair(basis_type(ibasis,ie), katomtype))) then 
                                    rlk(4) = sqrt(rsq)
                                    p = rlk(4)*rdr+1.0_wp
                                    m = int(p)
                                    m = min(m, nr-1)
                                    p = p - m 
                                    p = min(p, 1.0)
                                    coef = edens_spline(:,m,type_to_dens(katomtype, basis_type(ibasis, ie)))
                                    edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                    edens_host_intpo(jep, ie) = edens_host_intpo(jep, ie) + edens
                                end if
                            end do
                        end if
                    end do
                end do
            end do

            !Now communicate the edens for the integration points
            call comm_edens_intpo

            embed_intpo = 0.0_wp
            !Now discretely calculate the densities for all the needed atomaps if desired
            if(atomap_neighbors) then 
                do iatomap = 1, atomap_num_l
                    !If it's an integration point then we already have the electron density  
                    if(atomap_to_intpo(1, iatomap) /= 0) then 
                        jep = atomap_to_intpo(1,iatomap)
                        ie = atomap_to_intpo(2,iatomap)
                        edens_host_atomap(iatomap) = edens_host_intpo(jep, ie)
                    !Otherwise if we need this density
                    else if (needed_atomap(iatomap)) then 
                        do i = 1, 3
                            rl(i) = r_atomap(i, iatomap)
                        end do

                        !Loop over all atomaps
                        do nei = 1, n_atomap_cg(iatomap, en_list)
                            ja = atomap_nei_cg(nei, iatomap, en_list)
                            katomtype = type_atomap(ja)
                            rk = r_atomap(:,ja)
                            do i = 1, 3
                                rlk(i) = rk(i) - rl(i)
                            end do
                            rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)

                            if(rsq < pair_cutoffsq(type_to_pair(type_atomap(iatomap), katomtype))) then 
                                rlk(4) = sqrt(rsq)
                                p = rlk(4)*rdr+1.0_wp
                                m = int(p)
                                m = min(m, nr-1)
                                p = p - m 
                                p = min(p, 1.0)
                                coef = edens_spline(:,m,type_to_dens(katomtype, type_atomap(iatomap)))
                                edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                edens_host_atomap(iatomap) = edens_host_atomap(iatomap) + edens
                            end if
                        end do

                        !Loop over all ghost atoms
                        do nei = 1, n_atomap_at(iatomap, en_list)
                            ja = atomap_nei_at(nei, iatomap, en_list)
                            katomtype = type_atom(ja)
                            rk = r_atom(:,ja)
                            do i = 1, 3
                                rlk(i) = rk(i) - rl(i)
                            end do
                            rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                            if(rsq < pair_cutoffsq(type_to_pair(type_atomap(iatomap), katomtype))) then 
                                rlk(4) = sqrt(rsq)
                                p = rlk(4)*rdr+1.0_wp
                                m = int(p)
                                m = min(m, nr-1)
                                p = p - m 
                                p = min(p, 1.0)
                                coef = edens_spline(:,m,type_to_dens(katomtype, type_atomap(iatomap)))
                                edens = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                edens_host_atomap(iatomap) = edens_host_atomap(iatomap) + edens
                            end if
                        end do
                    end if
                end do

                !Now loop over all atomaps and calculate the embedding energy derivative if needed
                do ia = 1, atomap_num_l
                    if(needed_atomap(ia)) then 
                        !First calculate the embedding energy contribution to the energy
                        p = edens_host_atomap(ia)*rdrho+1.0_wp
                        m = int(p)
                        m = max(1, min(m, nrho-1))
                        p = p -m
                        p = min(p, 1.0)
                        coef = embed_spline(:,m,pot_map(type_atomap(ia)))
                        embed_atom = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                        d_embed_l = (coef(1)*p + coef(2))*p + coef(3)
                        !Extra linear term to conserve energy when atoms are too close
                        if(edens_host_atomap(ia) > edensmax) then 
                            embed_atom = embed_atom + d_embed_l*(edens_host_atomap(ia) - edensmax)
                        end if
                        embed_deriv_atomap(ia) = d_embed_l

                        if(atomap_to_intpo(1,ia) /= 0) then 
                            jep = atomap_to_intpo(1,ia)
                            ie = atomap_to_intpo(2,ia)
                            embed_intpo(jep,ie) = embed_atom
                            embed_deriv_intpo(jep,ie) = d_embed_l
                        end if
                    end if
                end do
            else
                !If not doing it discretely then set the edenshost for all atomaps based on the integration point values 
                !First we add the embedding energy to integration point energy and calculate the erivative of the embedding energy

                do ie = 1, ele_num_l
                    do iep = 1, intpo_count(itype(ie))
                        do ibasis = 1, basis_num(ie)
                            if(pot_map(basis_type(ibasis, ie)) > 0) then 
                                jep = basis_num(ie)*(iep-1) + ibasis
                                !First calculate the embedding energy contribution to the energy
                                p = edens_host_intpo(jep,ie)*rdrho+1.0_wp
                                m = int(p)
                                m = max(1, min(m, nrho-1))
                                p = p -m
                                p = min(p, 1.0)
                                coef = embed_spline(:,m,pot_map(basis_type(ibasis,ie)))
                                embed_atom = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)
                                d_embed_l = (coef(1)*p + coef(2))*p + coef(3)
                                !Extra linear term to conserve energy when atoms are too close
                                if(edens_host_intpo(jep,ie) > edensmax) then 
                                    embed_atom = embed_atom + d_embed_l*(edens_host_intpo(jep,ie) - edensmax)
                                end if
                                embed_intpo(jep,ie) = embed_atom
                                embed_deriv_intpo(jep,ie) = d_embed_l
                            end if
                        end do
                    end do
                end do
                do ie = 1, ele_num_l
                    select case(etype(ie))
                    case(1,2,3)
                        iatom_counts = (size_ele(ie)+1)**3
                    end select
                    do iatom = 1, iatom_counts
                        do ibasis = 1, basis_num(ie)
                            iatomap = cg_atomap(basis_num(ie)*(iatom-1)+ibasis, ie)
                            if(iatomap /= 0) then 
                                iep = who_rep_atomap(iatom, size_to_shape(size_ele(ie)), itype(ie))
                                jep = basis_num(ie)*(iep-1) + ibasis
                                edens_host_atomap(iatomap) = edens_host_intpo(jep, ie)
                                embed_deriv_atomap(iatomap) = embed_deriv_intpo(jep,ie)
                            end if
                        end do
                    end do
                end do 
            end if
            
            !Now communicate ghost virtual atom edens
            call comm_edens_ghost_cg
        end if



        talliesme(:) = 0

        if(ele_num>0) then 
            force_intpo(:,:,:) = 0.0_wp
            pair_intpo(:,:) = 0.0_wp
            if(need_virial) then 
                virial_intpo(:,:,:,:) = 0.0_wp
            end if
        end if

        !Update atomistic region if needed 

        if(atom_num > 0) then 
            !Loop over all atoms
            do ia = 1, atom_num_l

                rl(:) = r_atom(:, ia)
                latomtype = type_atom(ia)
                edens_l = edens_host_atom(ia)
                d_embed_l = embed_deriv_atom(ia)

            
                !Now calculate the pair potential from atoms
                do nei = 1, n_at_at(ia, en_list)

                    ja = at_nei_at(nei, ia, en_list)
                    rk(:) = r_atom(:, ja)
                    katomtype = type_atom(ja)
                    edens_k = edens_host_atom(ja)

                    rlk(1:3) = rk - rl
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < pair_cutoffsq(type_to_pair(latomtype, katomtype))) then 
                        rlk(4) = sqrt(rsq)
                        !Get derivatives for edens
                        p = rlk(4)*rdr+1.0_wp
                        m = int(p)
                        m = min(m, nr-1)
                        p = p -m 
                        p = min(p, 1.0)
                        coef = edens_spline(:,m, type_to_dens(latomtype, katomtype))
                        d_edens_l  = (coef(1)*p + coef(2))*p + coef(3)
                        coef = edens_spline(:,m, type_to_dens(katomtype, latomtype))
                        d_edens_k  = (coef(1)*p + coef(2))*p + coef(3)


                        d_embed_k = embed_deriv_atom(ja)

                        !Now calculate pair energy and pair force
                        coef=pair_spline(:, m, type_to_pair(latomtype, katomtype))
                        z2p = (coef(1)*p + coef(2))*p + coef(3)
                        z2 = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)


                        !Now calculate force and energy
                        recip = 1.0/rlk(4)
                        phi = z2*recip
                        phip = z2p*recip - phi*recip
                        psip = d_embed_l*d_edens_k + d_embed_k*d_edens_l + phip
                        f_atom = psip*recip

                        !Update force
                        flk(:) = f_atom * rlk(1:3)
                        force_atom(:, ia) = force_atom(:, ia) + flk(:)

                        if(need_virial) then
                            do ic = 1, 3
                                virial_atom(:, ic, ia) = virial_atom(:, ic, ia) + flk(:) * rlk(ic) / 2.0_wp
                            end do
                        end if

                        !Update pair
                        energy_atom(ia) = energy_atom(ia) + phi/2

                        !Newtons law and symmetric pair potential
                        if(ja <= atom_num_l) then
                            force_atom(:, ja) = force_atom(:, ja) - flk(:)

                            if(need_virial) then
                                do ic = 1, 3
                                    virial_atom(:, ic, ja) = virial_atom(:, ic, ja) + flk(:) * rlk(ic) / 2.0_wp
                                end do
                            end if

                            !Update pair
                            energy_atom(ja) = energy_atom(ja) + phi/2
                        end if
                    end if
                end do
                !Now calculate the pair potential from virtual atoms
                do nei = 1, n_at_cg(ia, en_list)

                    ja = at_nei_cg(nei, ia, en_list)
                    rk(:) = r_atomap(:, ja)
                    katomtype = type_atomap(ja)
                    edens_k = edens_host_atomap(ja)

                    rlk(1:3) = rk - rl
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < pair_cutoffsq(type_to_pair(latomtype, katomtype))) then 
                        rlk(4) = sqrt(rsq)
                        !Get derivatives for edens
                        p = rlk(4)*rdr+1.0_wp
                        m = int(p)
                        m = min(m, nr-1)
                        p = p -m 
                        p = min(p, 1.0)
                        coef = edens_spline(:,m, type_to_dens(latomtype, katomtype))
                        d_edens_l  = (coef(1)*p + coef(2))*p + coef(3)
                        coef = edens_spline(:,m, type_to_dens(katomtype, latomtype))
                        d_edens_k  = (coef(1)*p + coef(2))*p + coef(3)

                        !Now calculate pair energy and pair force
                        coef=pair_spline(:, m, type_to_pair(latomtype, katomtype))
                        z2p = (coef(1)*p + coef(2))*p + coef(3)
                        z2 = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)

                        !Get embedding energy derivative for neighbor
                        d_embed_k = embed_deriv_atomap(ja)

                        !Now calculate force and energy
                        recip = 1.0/rlk(4)
                        phi = z2*recip
                        phip = z2p*recip - phi*recip
                        psip = d_embed_l*d_edens_k + d_embed_k*d_edens_l + phip
                        f_atom = psip*recip

                        !Update force
                        flk(:) = f_atom * rlk(1:3)
                        force_atom(:, ia) = force_atom(:, ia) + flk(:)

                        if(need_virial) then
                            do ic = 1, 3
                                virial_atom(:, ic, ia) = virial_atom(:, ic, ia) + flk(:) * rlk(ic) / 2.0_wp
                            end do
                        end if

                        !Update pair
                        energy_atom(ia) = energy_atom(ia) + phi/2


                        !Update for intpo neighbors
                        if(ja <= atomap_num_l) then 
                            if(atomap_to_intpo(1,ja)/=0) then 
                                jep = atomap_to_intpo(1,ja)
                                ie = atomap_to_intpo(2,ja)

                                force_intpo(:,jep,ie) = force_intpo(:,jep,ie) - flk(:)
                                if(need_virial)then 
                                    do ic = 1, 3
                                        virial_intpo(:, ic, jep,ie) = virial_intpo(:, ic, jep,ie) + flk(:) * rlk(ic) / 2.0_wp
                                    end do
                                end if
                                pair_intpo(jep,ie) = pair_intpo(jep,ie) + phi/2.0_wp
                            end if
                        end if
                    end if
                end do
            end do
        end if


        !Now actually update the force for the elements
        if(ele_num > 0) then 

            energy_intpo=0.0_wp
            !Loop over all elements
            do ie = 1, ele_num_l
                !Initialize integration point variables
                do iep = 1, intpo_count(itype(ie))

                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1)+ibasis
                        !If we have the integration point then calculate  
                        if(who_has_intpo(iep, ie).eqv..true.) then

                            iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                            iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                            latomtype = type_atomap(iatomap)
                            rl(:) = r_atomap(:, iatomap) 
                            edens_l = edens_host_intpo(jep,ie)
                            d_embed_l = embed_deriv_atomap(iatomap)

                            !Now loop over all virtual atoms
                            do nei = 1, n_cg_cg(jep, ie, en_list)

                                ja = cg_nei_cg(nei, jep, ie, en_list)

                                rk(:) = r_atomap(:, ja)
                                katomtype = type_atomap(ja)
                                edens_k = edens_host_atomap(ja)
                                  
                                rlk(1:3) = rk - rl
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < pair_cutoffsq(type_to_pair(latomtype,katomtype))) then 


                                    rlk(4) = sqrt(rsq)
                                    
                                    !Get derivatives for edens
                                    p = rlk(4)*rdr+1.0_wp
                                    m = int(p)
                                    m = min(m, nr-1)
                                    p = p -m 
                                    p = min(p, 1.0)
                                    coef = edens_spline(:,m, type_to_dens(latomtype, katomtype))
                                    d_edens_l  = (coef(1)*p + coef(2))*p + coef(3)
                                    coef = edens_spline(:,m, type_to_dens(katomtype, latomtype))
                                    d_edens_k  = (coef(1)*p + coef(2))*p + coef(3)
 
                                    !Now calculate pair energy and pair force
                                    coef=pair_spline(:, m, type_to_pair(latomtype, katomtype))
                                    z2p = (coef(1)*p + coef(2))*p + coef(3)
                                    z2 = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)

                                    !Get embedding energy derivative for neighbor
                                    d_embed_k = embed_deriv_atomap(ja)

                                    !Now calculate force and energy
                                    recip = 1.0/rlk(4)
                                    phi = z2*recip
                                    phip = z2p*recip - phi*recip
                                    psip = d_embed_l*d_edens_k + d_embed_k*d_edens_l + phip
                                    f_intpo = psip*recip
                                    
                                    !Update force
                                    flk(:) = f_intpo * rlk(1:3)

                                    force_intpo(:, jep, ie) = force_intpo(:, jep, ie) + flk(:)

                                    if(need_virial) then
                                        do ic = 1, 3
                                            virial_intpo(:, ic, jep, ie) = virial_intpo(:, ic, jep, ie) &
                                                                                + flk(:) * rlk(ic) / 2.0_wp
                                        end do
                                    end if

                                    !Update pair
                                    pair_intpo(jep, ie) = pair_intpo(jep, ie) + phi/2.0_wp
                                end if
                            end do

                            !Now calculate the pair potential from the atoms
                            do nei = 1, n_cg_at(jep, ie, en_list)

                                ja = cg_nei_at(nei, jep, ie, en_list)

                                rk(:) = r_atom(:, ja)
                                katomtype = type_atom(ja)
                                edens_k = edens_host_atom(ja)
                                  
                                rlk(1:3) = rk - rl
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < pair_cutoffsq(type_to_pair(latomtype,katomtype))) then 

                                    rlk(4) = sqrt(rsq)
                                    
                                    !Get derivatives for edens
                                    p = rlk(4)*rdr+1.0_wp
                                    m = int(p)
                                    m = min(m, nr-1)
                                    p = p -m 
                                    p = min(p, 1.0)
                                    coef = edens_spline(:,m, type_to_dens(latomtype, katomtype))
                                    d_edens_l  = (coef(1)*p + coef(2))*p + coef(3)
                                    coef = edens_spline(:,m, type_to_dens(katomtype, latomtype))
                                    d_edens_k  = (coef(1)*p + coef(2))*p + coef(3)
 
                                    !Now calculate pair energy and pair force
                                    coef=pair_spline(:, m, type_to_pair(latomtype, katomtype))
                                    z2p = (coef(1)*p + coef(2))*p + coef(3)
                                    z2 = ((coef(4)*p + coef(5))*p + coef(6))*p + coef(7)

                                    !Get embedding energy derivative for neighbor
                                    d_embed_k = embed_deriv_atom(ja)

                                    !Now calculate force and energy
                                    recip = 1.0/rlk(4)
                                    phi = z2*recip
                                    phip = z2p*recip - phi*recip
                                    psip = d_embed_l*d_edens_k + d_embed_k*d_edens_l + phip
                                    f_intpo = psip*recip
                                    
                                    !Update force
                                    flk(:) = f_intpo * rlk(1:3)
                                    force_intpo(:, jep, ie) = force_intpo(:, jep, ie) + flk(:)

                                    if(need_virial) then
                                        do ic = 1, 3
                                            virial_intpo(:, ic, jep, ie) = virial_intpo(:, ic, jep, ie) &
                                                                                + flk(:) * rlk(ic) / 2.0_wp
                                        end do
                                    end if

                                    !Update pair
                                    pair_intpo(jep, ie) = pair_intpo(jep, ie) + phi/2.0_wp
                                end if
                            end do
                        end if
                    end do
                end do

                !Now get the integration point energy
                do iep = 1, intpo_count(itype(ie))
                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1)+ibasis
                        if(who_has_intpo(iep, ie).eqv..true.) then
                                energy_intpo(jep) = pair_intpo(jep, ie) + embed_intpo(jep,ie)
                        end if
                    end do
                end do


                !update integration point values using weight

                do iep = 1, intpo_count(itype(ie))
                    do ibasis = 1, basis_num(ie)
                        jep = basis_num(ie)*(iep-1)+ibasis
                        iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                        iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)

                        if(who_has_intpo(jep, ie).eqv..true.) then
                            force_intpo(:, jep, ie) = force_intpo(:, jep, ie) &
                                                      * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))

                            if(need_virial) then
                              virial_intpo(:, :, jep, ie) = virial_intpo(:, :, jep, ie) &
                                                      * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                            end if
                            energy_intpo(jep) = energy_intpo(jep) * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                        end if
                    end do
                end do

                !calculate force, virial, and energy for nodes from integration point values
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do iep = 1, intpo_count(itype(ie))
                        do ibasis = 1, basis_num(ie)
                            jep = basis_num(ie)*(iep-1) + ibasis
                            if(who_has_intpo(iep, ie).eqv..true.) then

                                iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                                eshape = a_interpo(inod, iatom, size_to_shape(size_ele(ie)), etype(ie))
                                force(:, ibasis, ip ) = force(:, ibasis, ip) + eshape * force_intpo(:, jep, ie)

                                if(need_virial) then
                                    virial(:, :, ibasis, ip) = virial(:, :, ibasis, ip) &
                                                       + eshape * virial_intpo(:, :, jep, ie)
                                end if

                                energy(ibasis, ip) = energy(ibasis, ip) + eshape*energy_intpo(jep)

                            end if
                        end do
                    end do
                end do
            end do
            
        end if

        return
    end subroutine update_force_eam

end module eam
