module neighbors
   !This contains all of the neighbor running code
    use parameters
    use integration
    use math
    use comms
    use logger
    use forces
    use elements
    use group
    use temp
    use min_arrays

    implicit none

    logical :: reneighbored, nei_init

    !Variables associated ith the actual linked cell list data structures
    integer, private :: num_cell(3), cell_num, num_cell_max, cell_atomap_lim, cell_atom_lim, delay, last_reneighbor
    integer, allocatable :: cell_atomap(:,:), cell_atom(:,:), num_cell_atomap(:), num_cell_atom(:), v_cell_nei(:), &
                            num_celldiv_atpo(:,:), cell_neighbor(:,:), which_cell_intpo(:,:), num_neighbor_cell(:), &
                            which_cell_atom(:), which_cell_atomap(:)
    logical, allocatable :: bound_cell(:,:)
    real(kind = wp) :: rc_bin

    !Variables for the actual neighbor lists
    integer, private :: nei_lim, builds, nei_lists
    integer, allocatable :: cg_nei_at(:,:,:,:), n_cg_at(:,:,:), cg_nei_cg(:,:,:,:), n_cg_cg(:,:,:), &
                            at_nei_at(:,:,:), n_at_at(:,:), at_nei_cg(:,:,:), n_at_cg(:,:), &
                            atomap_nei_cg(:,:,:), n_atomap_cg(:,:), atomap_nei_at(:,:,:), n_atomap_at(:, :)
    
    integer, save, dimension(max_pot_types):: neilist_to_pot
    logical :: def_nei, neis_updated, atomap_neighbors, need_all_step

    !Keep position of atoms and elements last time they were checked 
    real(kind=wp), allocatable :: r_update(:,:,:), r_atom_update(:,:)

    public 
    contains


subroutine neighbors_defaults
    !set default values 
    rc_neigh = 0.0_wp
    def_nei = .false.
    rc_bin = 1.0_wp
    neis_updated = .false.
    delay = 1
    last_reneighbor = 0
    builds = 0 
    atomap_neighbors = .false.
    nei_init=.false.
    need_all_step=.true.
end subroutine neighbors_defaults

subroutine parse_neighbors(line)
    !Parse neighbor command to get neighbor rc_neigh
    character(len = *), intent(in) :: line
    character(len = read_len) :: label, args(10), msg
    integer :: iospara, i, j

    read(line, *, iostat = iospara) label, rc_bin

    if(.not.def_nei) then 
        call misc_error("Please define potential before running neighbor command")
    else if (def_nei) then 
        call set_rc_neigh(rc_neigh+rc_bin)
    end if
    
    !Read additional arguments
    j = tok_count(line)
    if (j > 2) then 
        read(line, *, iostat = iospara)  args
        i = 3
        do while(i <= j)
            select case(args(i))
            case('delay')
                i = i + 1
                read(args(i), *) delay
                if(delay < 0) then 
                    call command_error("Delay must be greater than 0 for neighbor command")
                end if
                i=i+1
                write(msg, *) "Neighbor delay set to ", delay
                call log_msg(msg)
            case('no_need_all_step')
                need_all_step=.false.
            case('need_all_step')
                need_all_step=.true.
            case default
                write(msg, *) "Keyword ", trim(adjustl(args(i))), " is not accepted in neighbor command"
                call command_error(msg)
            end select
        end do
    end if

    return
end subroutine parse_neighbors


subroutine set_rc_neigh(val)
    real(kind = wp), intent(in) :: val
    character(len = read_len) :: msg

    rc_neigh = val
    write(msg, *) "Neighbor list cutoff distance set to ", rc_neigh
    call log_msg(msg)

    !Now call the processor boundary code to set up the ghost boundaries
    call processor_bds
end subroutine set_rc_neigh

subroutine print_neighbor_info
    integer :: iep, ie
    integer :: max_nei, min_nei, max_neime, min_neime
    real(kind = wp) :: avg_neime, avg_nei
    character(len=read_len) :: msg
    !This subroutine just prints some information about the neighbor lists
    
    !Get min, max, and average number of neighbors
    min_neime = minval(n_at_at(1:atom_num_l, :)+n_at_cg(1:atom_num_l, :))
    max_neime = maxval(n_at_at(1:atom_num_l, :)+n_at_cg(1:atom_num_l, :))
    avg_neime = sum(n_at_at(1:atom_num_l, :)+n_at_cg(1:atom_num_l, :))

    do ie = 1, ele_num_l
        do iep = 1, intpo_count(itype(ie))*basis_num(ie)
            if(who_has_intpo(iep, ie) ) then  
                min_neime = min(min_neime, minval(n_cg_at(iep, ie, :)+n_cg_cg(iep,ie, :)))
                max_neime = max(max_neime, maxval(n_cg_at(iep,ie, :)+n_cg_cg(iep,ie, :)))
                avg_neime = avg_neime + sum(n_cg_at(iep, ie, :) + n_cg_cg(iep, ie, :))
            end if
        end do
    end do

    !Now reduce them 
    avg_neime = real(avg_neime,wp)/real(atom_num+intpo_num,wp)

    call mpi_reduce(min_neime, min_nei, 1, mpi_integer, mpi_min, root, world, ierr)
    call mpi_reduce(max_neime, max_nei, 1, mpi_integer, mpi_max, root, world, ierr)
    call mpi_reduce(avg_neime, avg_nei, 1, mpi_wp, mpi_sum, root, world, ierr)

    if(rank == root) then 
        write(msg, *) "Min, max, and average neighbor counts are:", min_nei, max_nei, avg_nei
        call log_msg(msg)

        write(msg, *) "Neighbor list cutoff is: ", rc_neigh
        call log_msg(msg)
    end if
end subroutine print_neighbor_info

subroutine update_arrays
    !This subroutine allocates the update arrays if necessary
    if(node_num_l > 0) then 
        !Deallocate if our arrays grew
        if(allocated(r_update)) then 
            if(node_num_l > size(r_update,3)) deallocate(r_update)
        end if

        !Allocate if necessary
        if(.not.allocated(r_update)) then 
            allocate(r_update(3,max_basisnum, node_num_lr), stat=allostat)
            if(allostat > 0 ) call alloc_error("Failure allocating r_update in update_arrays", allostat)
        end if

        !Assign r_update
        r_update(:,:,1:node_num_l) = r(:,:,1:node_num_l)
    end if

    !This subroutine allocates the update arrays if necessary
    if(atom_num_l > 0) then 
        !Deallocate if our arrays grew
        if(allocated(r_atom_update))  then 
            if (atom_num_l > size(r_atom_update,2)) deallocate(r_atom_update)
        end if

        !Allocate if necessary
        if(.not.allocated(r_atom_update)) then 
            allocate(r_atom_update(3,atom_num_lr))
            !if(allostat > 0 ) call alloc_error("Failure allocating r_atom_update in update_arrays", allostat)
        end if

        !Assign r_update
        r_atom_update(:,1:atom_num_l) = r_atom(:,1:atom_num_l)
    end if
end subroutine update_arrays

subroutine neighbor_lists
    real(kind = wp) :: t_start, t_end
    !This code initializes neighbor lists
    !If init == 0 then that means this is not the first time this was called
    ! whereas if init == 1 then it is the first time this function is called
    integer :: i, j

    !Start timer
    t_start = mpi_wtime()
    
    !First call cell initialization functions
    if(allocated(which_cell_intpo).or.allocated(which_cell_atom)) call dealloc_cell
    call cell_init
    call update_cell

    !Allocate arrays as needed
    nei_lim = 100
    if (.not. nei_init) then 
        j=1
        neilist_to_pot=0
        do i = 1, max_pot_types
            if(potential_types(i)) then 
                neilist_to_pot(j)=i
                j=j+1
            end if
        end do
    end if
    if(nei_init) call dealloc_nei_arrays
    call alloc_nei_arrays

    if(ele_num_l > 0) then
        needed_atomap=.false.
        call neighbor_list_cg
    end if

    if(atom_num_l> 0) call neighbor_list_at

    if(atomap_neighbors) then 
        call neighbor_list_atomap
    else
        !If we aren't fully calculating the cluster potential then we just set everything to false again
        if(ele_num_l > 0) needed_atomap=.false.
    end if


    !Assign update arrays
    call update_arrays

    !Get original centroid and log neighbor info
    if(.not.nei_init) then 
        call get_centroid(center_mass_ori)
        call print_neighbor_info
    end if

    t_end = mpi_wtime()

    walltime(1) = walltime(1) + (t_end - t_start)
    
    !Increment neighbor list
    builds = builds + 1
    nei_init=.true.

    return
end subroutine neighbor_lists

subroutine alloc_nei_arrays
    !Allocate neighbor arrays
    integer :: i, j
    nei_lists=0
    j=1
    do i = 1, max_pot_types
        if(potential_types(i)) then 
            nei_lists = nei_lists+1
            neilist_to_pot(j)=i
            j=j+1
        end if
    end do 

    if(ele_num_l > 0) allocate( cg_nei_cg(100, max_basisnum*max_intpo_num, ele_num_l, nei_lists), &
                                n_cg_cg(max_basisnum*max_intpo_num, ele_num_l, nei_lists), &
                                cg_nei_at(100, max_basisnum*max_intpo_num, ele_num_l, nei_lists), &
                                n_cg_at(max_basisnum*max_intpo_num, ele_num_l, nei_lists), &
                                needed_atomap(atomap_num_l), &
                                stat = allostat)

    if(atom_num_l > 0) allocate(at_nei_at(100, atom_num_l, nei_lists), &
                                n_at_at(atom_num_l, nei_lists), &
                                at_nei_cg(100, atom_num_l, nei_lists), &
                                n_at_cg(atom_num_l, nei_lists), &
                                stat = allostat)

    if(atomap_neighbors) allocate(atomap_nei_cg(100, atomap_num_l, nei_lists), &
                                  n_atomap_cg(atomap_num_l, nei_lists), &
                                  atomap_nei_at(100, atomap_num_l, nei_lists), &
                                  n_atomap_at(atomap_num_l, nei_lists), &
                                  stat = allostat)

    return 
end subroutine alloc_nei_arrays

subroutine dealloc_nei_arrays
    !deallocate neighbor arrays

    if(ele_num_l > 0) deallocate(cg_nei_cg, n_cg_cg, cg_nei_at, n_cg_at, needed_atomap, stat = allostat)

    if(atom_num_l > 0) deallocate(at_nei_at, n_at_at, at_nei_cg, n_at_cg, stat = allostat)

    if(atomap_neighbors) deallocate(atomap_nei_at, n_atomap_at, atomap_nei_cg, n_atomap_cg)
    return 
end subroutine dealloc_nei_arrays

subroutine dealloc_cell
    !Deallocate cell arrays 
    deallocate(cell_atomap, cell_atom, num_cell_atomap, num_cell_atom, num_celldiv_atpo, cell_neighbor, &
               v_cell_nei, num_neighbor_cell, bound_cell, stat = allostat)

    !Deallocate atom/atomap location arrays
    deallocate(which_cell_intpo, which_cell_atom, which_cell_atomap)

    return
end subroutine dealloc_cell

subroutine cell_init
    !This subroutine initializes various things required for the linked cell list neighbor listing code 
    integer :: i, v_cell_temp(3)

    !Calculate cell numbers along each dimension and total
    do i = 1, 3
        num_cell(i) = int(pro_length_out(i)/rc_neigh)
        if (num_cell(i) == 0) num_cell(i) = 1
    end do
    cell_num = product(num_cell)
    num_cell_max=maxval(num_cell)+1
    cell_atom_lim = 100
    cell_atomap_lim = 100
    
    !Allocate cell arrays now 
    allocate( cell_atomap(100, cell_num), &
              cell_atom(100, cell_num), &
              num_cell_atomap(cell_num), &
              num_cell_atom(cell_num), &
              num_celldiv_atpo(num_cell_max, 3), &
              cell_neighbor(27, cell_num), &
              v_cell_nei(27), &
              num_neighbor_cell(cell_num), &
              bound_cell(6,27), &
              stat = allostat)

    !Allocate atom/atomap location arrays
    allocate( which_cell_intpo(max_intpo_num*max_basisnum, ele_num_l), &
              which_cell_atom(atom_num_lg), which_cell_atomap(atomap_num_lg))

    !NOw initialize some variables for listing the neighboring cells
    v_cell_temp(1) = 1
    v_cell_temp(2) = num_cell(1)
    v_cell_temp(3) = num_cell(2)*num_cell(1)

    v_cell_nei(:) = 0.0_wp
    bound_cell(:,:) = .false.

    !Neighbors in x direction
    v_cell_nei(1) = v_cell_temp(1)
    bound_cell(2, 1) = .true.

    v_cell_nei(2) = -v_cell_temp(1)
    bound_cell(1, 2) = .true.

    !y direction

    v_cell_nei(3) = v_cell_temp(2)
    bound_cell(4, 3) = .true.

    v_cell_nei(4) = -v_cell_temp(2)
    bound_cell(3, 4) = .true.

    !x-y plane

    v_cell_nei(5) = v_cell_nei(1)+v_cell_nei(3)
    bound_cell(2, 5) = .true.
    bound_cell(4, 5) = .true.

    v_cell_nei(6) = v_cell_nei(2)+v_cell_nei(3)
    bound_cell(1, 6) = .true.
    bound_cell(4, 6) = .true.

    v_cell_nei(7) = v_cell_nei(1)+v_cell_nei(4)
    bound_cell(2, 7) = .true.
    bound_cell(3, 7) = .true.

    v_cell_nei(8) = v_cell_nei(2)+v_cell_nei(4)
    bound_cell(1, 8) = .true.
    bound_cell(3, 8) = .true.

    v_cell_nei(9) = v_cell_temp(3)
    bound_cell(6, 9) = .true.

    v_cell_nei(10) = v_cell_nei(1)+v_cell_nei(9)
    bound_cell(2, 10) = .true.
    bound_cell(6, 10) = .true.

    v_cell_nei(11) = v_cell_nei(2)+v_cell_nei(9)
    bound_cell(1, 11) = .true.
    bound_cell(6, 11) = .true.

    v_cell_nei(12) = v_cell_nei(3)+v_cell_nei(9)
    bound_cell(4, 12) = .true.
    bound_cell(6, 12) = .true.

    v_cell_nei(13) = v_cell_nei(4)+v_cell_nei(9)
    bound_cell(3, 13) = .true.
    bound_cell(6, 13) = .true.

    v_cell_nei(14) = v_cell_nei(5)+v_cell_nei(9)
    bound_cell(2, 14) = .true.
    bound_cell(4, 14) = .true.
    bound_cell(6, 14) = .true.

    v_cell_nei(15) = v_cell_nei(6)+v_cell_nei(9)
    bound_cell(1, 15) = .true.
    bound_cell(4, 15) = .true.
    bound_cell(6, 15) = .true.

    v_cell_nei(16) = v_cell_nei(7)+v_cell_nei(9)
    bound_cell(2, 16) = .true.
    bound_cell(3, 16) = .true.
    bound_cell(6, 16) = .true.

    v_cell_nei(17) = v_cell_nei(8)+v_cell_nei(9)
    bound_cell(1, 17) = .true.
    bound_cell(3, 17) = .true.
    bound_cell(6, 17) = .true.

    !lower x-y plane

    v_cell_nei(18) = -v_cell_temp(3)
    bound_cell(5, 18) = .true.

    v_cell_nei(19) = v_cell_nei(1)+v_cell_nei(18)
    bound_cell(2, 19) = .true.
    bound_cell(5, 19) = .true.

    v_cell_nei(20) = v_cell_nei(2)+v_cell_nei(18)
    bound_cell(1, 20) = .true.
    bound_cell(5, 20) = .true.

    v_cell_nei(21) = v_cell_nei(3)+v_cell_nei(18)
    bound_cell(4, 21) = .true.
    bound_cell(5, 21) = .true.

    v_cell_nei(22) = v_cell_nei(4)+v_cell_nei(18)
    bound_cell(3, 22) = .true.
    bound_cell(5, 22) = .true.

    v_cell_nei(23) = v_cell_nei(5)+v_cell_nei(18)
    bound_cell(2, 23) = .true.
    bound_cell(4, 23) = .true.
    bound_cell(5, 23) = .true.

    v_cell_nei(24) = v_cell_nei(6)+v_cell_nei(18)
    bound_cell(1, 24) = .true.
    bound_cell(4, 24) = .true.
    bound_cell(5, 24) = .true.

    v_cell_nei(25) = v_cell_nei(7)+v_cell_nei(18)
    bound_cell(2, 25) = .true.
    bound_cell(3, 25) = .true.
    bound_cell(5, 25) = .true.

    v_cell_nei(26) = v_cell_nei(8)+v_cell_nei(18)
    bound_cell(1, 26) = .true.
    bound_cell(3, 26) = .true.
    bound_cell(5, 26) = .true.

    call update_cell_neighbor
end subroutine cell_init

subroutine update_cell_neighbor
    !This code gets the neighboring cells to build the cell_neighbor list

    integer :: i, ice, jce, nei, nei_real
    integer, dimension(27) :: cell_neighbor_temp, cell_neighbor_real

    !cell_neighbor and num_neighbor_cell
    cell_neighbor(:, :) = 0
    num_neighbor_cell(:) = 0

    do ice = 1, cell_num

        nei = 0
        cell_neighbor_temp(:) = 0

        do i = 1, 27
            jce = ice + v_cell_nei(i)

            if((jce > 0).and.(jce <= cell_num)) then
              nei = nei + 1
              cell_neighbor_temp(nei) = jce
            end if
        end do

        call delete_duplicate(cell_neighbor_temp, 27, cell_neighbor_real, nei_real)

        if(nei_real == 0) then
            print *, 'Error: Cell', ice, ' does not have neighboring cell'
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        cell_neighbor(1:nei_real, ice) = cell_neighbor_real(1:nei_real)
        num_neighbor_cell(ice) = nei_real

    end do

    return
end subroutine update_cell_neighbor

subroutine update_cell
    !Sort atom and atomaps to cells

    integer :: i, j, ia, ie, iep, ice, iatom, iatomap, ibasis
    integer, dimension(3) :: num_index
    integer, allocatable :: cell_atomp_array(:, :), num_cell_intpo(:)

    !coarse-grained domain
    if(ele_num /= 0) then

        !put atomap into cells
        cell_atomap(:, :) = 0
        which_cell_atomap(:) = 0
        num_cell_atomap(:) = 0

        do iatomap = 1, atomap_num_lg
            num_index(:) = 0

            do i = 1, 3
                num_index(i) = int((r_atomap(i, iatomap) - pro_bd_out(2*i-1))/rc_neigh)
                if(num_index(i) > num_cell(i)) then

                    print *, 'Error: Atomap', iatomap, ' of rank', rank, ' with position', &
                              r_atomap(i, iatomap), ' has wrong cell index', &
                              num_index(i), ' along', i, &
                             ' direction which is larger than', num_cell(i), ' .', &
                             ' The cell boundaries are', pro_bd_out(2*i-1:2*i)
                    call mpi_abort(mpi_comm_world, 1, ierr)

                else if(num_index(i) == num_cell(i)) then
                    num_index(i) = num_cell(i) - 1

                end if
            end do

            ice = 1 + num_index(1) + num_index(2) * num_cell(1) + num_index(3) * num_cell(1) * num_cell(2)

            if((ice < 1).or. (ice > cell_num)) then
                print *, 'Error: Wrong ice', ice, ' in update_cell', &
                         ' which should be between 1 and', cell_num
                call mpi_abort(mpi_comm_world, 1, ierr)

            else
                num_cell_atomap(ice) = num_cell_atomap(ice) + 1

                if(num_cell_atomap(ice) > cell_atomap_lim) then
                    allocate(cell_atomp_array(cell_atomap_lim+20, cell_num), stat = allostat)
                    if(allostat /= 0) call alloc_error("Failure to allocate cell_atomp_array in update_cell", allostat) 

                    cell_atomp_array(1:cell_atomap_lim, :) = cell_atomap(:, :)
                    cell_atomp_array(cell_atomap_lim+1:, :) = 0
                    call move_alloc(cell_atomp_array, cell_atomap)
                    cell_atomap_lim = cell_atomap_lim+20

                end if

                cell_atomap(num_cell_atomap(ice), ice) = iatomap
                which_cell_atomap(iatomap)  = ice

            end if
        end do

        !debug

        if(maxval(num_cell_atomap) > cell_atomap_lim) then
            print *, 'Error: cell_atomap_lim, which is', cell_atomap_lim, &
                     ' should be at least', maxval(num_cell_atomap)
            call mpi_abort(mpi_comm_world, 1, ierr)

        else if(sum(num_cell_atomap) /= atomap_num_lg) then
            print *, 'Error: Wrong number of atomap in cell', sum(num_cell_atomap), &
                     ' which should be', atomap_num_lg, ' which suggests that', &
                     ' probably some atomap are not in any cell'
            call mpi_abort(mpi_comm_world, 1, ierr)

        end if

        !put intpo into cells
        allocate(num_cell_intpo(cell_num), stat = allostat)
        if(allostat /= 0) call alloc_error("Failure to allocate num_cell_intpo in update_cell", allostat)

        num_cell_intpo(:) = 0 
        which_cell_intpo(:, :) = 0

        do ie = 1, ele_num_l
            j = size_to_shape(size_ele(ie))

            do iep = 1, intpo_count(itype(ie))
                do ibasis = 1, basis_num(ie)
                    if(who_has_intpo(basis_num(ie)*(iep-1)+ibasis, ie).eqv..true.) then

                        !Get the atomap position
                        iatom = atom_intpo(iep, j, itype(ie))
                        iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                        if(iatomap == 0) then
                          print *, 'Error: Iatomap', iatomap, ' can not be zero'
                          call mpi_abort(mpi_comm_world, 1, ierr)
                        end if

                        !Get the atomap cell index
                        num_index(:) = 0
                        do i = 1, 3
                            num_index(i) = int((r_atomap(i,iatomap) - pro_bd_out(2*i-1))/rc_neigh)
                            !check to make sure that the index is correct
                            if(num_index(i) > num_cell(i)) then
                                print *, 'Error: Intpo', iep, ' of element', ie, &
                                         ' has wrong cell index', num_index(i), ' along', &
                                          i, ' direction which is larger than', num_cell(i)
                                call mpi_abort(mpi_comm_world, 1, ierr)

                            else if(num_index(i) == num_cell(i)) then
                                num_index(i) = num_cell(i) - 1

                            end if
                        end do

                        !ice is the cell index 
                        ice = 1 + num_index(1) + num_index(2) * num_cell(1) + num_index(3) * num_cell(1) * num_cell(2)

                        !Again check to make sure that the cell index is correct
                        if((ice < 1).or.(ice > cell_num)) then
                            print *, 'Error: Wrong ice', ice, ' in update_cell cg which should be between 1 and', cell_num
                            call mpi_abort(mpi_comm_world, 1, ierr)

                        else
                            which_cell_intpo(basis_num(ie)*(iep-1)+ibasis, ie) = ice
                            num_cell_intpo(ice) = num_cell_intpo(ice) + 1
                        end if
                    end if
                end do
            end do

         end do

         if(sum(num_cell_intpo) /= intpo_num_l) then
            print *, 'Error: The number of integration point contained in all', &
                     ' cells', sum(num_cell_intpo), &
                     ' does not match intpo_num_l', intpo_num_l

            call mpi_abort(mpi_comm_world, 1, ierr)
         end if

         !IF no elements then we don't do anything
    else
     cell_atomap(:, :) = 0
     num_cell_atomap(:) = 0

    end if

    !put atom into cells
    if(atom_num /= 0) then
        cell_atom(:, :) = 0
        which_cell_atom(:) = 0
        num_cell_atom(:) = 0

        do ia = 1, atom_num_lg

            num_index(:) = 0

            do i = 1, 3
                  num_index(i) = int((r_atom(i, ia) - pro_bd_out(2*i-1)) / rc_neigh)
                  if(num_index(i) > num_cell(i)) then
                      print *, 'Error: Atom', ia, ' of rank', rank, ' with position', &
                                r_atom(i, ia), &
                               ' has wrong cell index', num_index(i), ' along', i, &
                               ' direction which is larger than', num_cell(i), ' .', &
                               ' The processor outer boundaries are', pro_bd_out(2*i-1:2*i)
                      call mpi_abort(mpi_comm_world, 1, ierr)

                  else if(num_index(i) == num_cell(i)) then
                      num_index(i) = num_cell(i) - 1
                  end if
            end do

            ice = 1 + num_index(1) + num_index(2) * num_cell(1) + num_index(3) * num_cell(1) * num_cell(2)

            if((ice < 1).or.(ice > cell_num)) then
                print *, 'Error: Wrong ice', ice, ' in update_cell atomistic', &
                         ' which should be between 1 and', cell_num

                call mpi_abort(mpi_comm_world, 1, ierr)

            else
                num_cell_atom(ice) = num_cell_atom(ice) + 1

                if(num_cell_atom(ice) > cell_atom_lim) then

                    allocate( cell_atomp_array(cell_atom_lim+20, cell_num), stat = allostat)
                    if(allostat /= 0) call alloc_error("Failure to allocate cell_atom_array", allostat)

                    cell_atomp_array(1:cell_atom_lim, :) = cell_atom(:, :)
                    cell_atomp_array(cell_atom_lim+1:, :) = 0
                    call move_alloc(cell_atomp_array, cell_atom)

                    cell_atom_lim = cell_atom_lim + 20

                end if

                cell_atom(num_cell_atom(ice), ice) = ia
                which_cell_atom(ia) = ice
            end if
        end do

        !debug

        if(maxval(num_cell_atom) > cell_atom_lim) then

          print *, 'Error: cell_atom_lim, which is', cell_atom_lim, &
                   ' should be at least', maxval(num_cell_atom)

          call mpi_abort(mpi_comm_world, 1, ierr)

        else if(sum(num_cell_atom) /= atom_num_lg) then

          print *, 'Error: The number of atoms contained in all cells', &
                    sum(num_cell_atom), ' does not match atom_num_lg', atom_num_lg

          call mpi_abort(mpi_comm_world, 1, ierr)

        end if

    else
      num_cell_atom(:) = 0
      cell_atom(:, :) = 0

    end if

    return
   end subroutine update_cell

subroutine neighbor_list_cg

    integer :: ie, iep, jep, ice, nei_ice, jce, iatom, jatom, iatomap, ja, jatomap, ibasis, n_at(nei_lists), n_cg (nei_lists), i
    real(kind = wp) :: rsq, rc_neisq
    real(kind = wp), dimension(3) :: rl, rk, rlk
    integer, allocatable :: intpo_neighbor_array(:, :, :, :)

    !Initialize variables
    cg_nei_cg = 0 
    cg_nei_at = 0
    n_cg_cg = 0
    n_cg_at = 0

    rc_neisq = rc_neigh*rc_neigh
    do ie = 1, ele_num_l

        do iep = 1, intpo_count(itype(ie))
            !Get the adjusted index taking into account the basis counts
            do ibasis = 1, basis_num(ie)
                jep = basis_num(ie)*(iep-1)+ibasis

                if(who_has_intpo(jep, ie).eqv..true.) then

                    iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))

                    !The way the code is structured, atom_intpo contains only the position of integration point when 
                    !counting the lattice sites of the finite elements. cg_atomap contains the full range of interpolated 
                    !atoms, with atoms in the same basis next to each other.This code loops over both atoms at the 
                    !integration point and calculates separate neighbor lists for them from the interpolated atom list.

                    iatomap = cg_atomap(basis_num(ie) *(iatom-1) + ibasis, ie)
                    rl(:) = r_atomap(:, iatomap)
                    ice = which_cell_intpo(jep, ie)

                    n_cg = 0
                    n_at = 0
                    !Loops over all neighboring cells.
                    do nei_ice = 1, num_neighbor_cell(ice)
                        jce = cell_neighbor(nei_ice, ice)

                        if((jce < 1).or.(jce > cell_num)) then
                            print *, 'Error: Wrong neighboring cell index', jce, &
                                     ' which should be between 1 and', cell_num
                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if

                        !Loops over all the interpolated atoms within that neighboring cell
                        do jatom = 1, num_cell_atomap(jce)
                            jatomap = cell_atomap(jatom, jce)

                            if(jatomap /= iatomap) then

                                !Check to see if these are neighbors, if they are then add them to the list
                                rk(:) = r_atomap(:, jatomap)
                                rlk = rk - rl
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < rc_neisq) then
                                    !Add neighbor to list
                                    do i= 1, nei_lists
                                        if(btest(types_to_pot_type(type_atomap(iatomap), type_atomap(jatomap)), &
                                                 neilist_to_pot(i))) then
                                            n_cg(i) = n_cg(i) + 1
                                            !Check to see if we need to resize arrays
                                            if(n_cg(i) > size(cg_nei_cg,1)) then
                                                nei_lim = size(cg_nei_cg,1)
                                                allocate(intpo_neighbor_array(nei_lim+20,max_basisnum*max_intpo_num,ele_num_l, &
                                                nei_lists), stat = allostat)

                                                if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", &
                                                                                     allostat)

                                                intpo_neighbor_array(:, :, :, :) = 0
                                                intpo_neighbor_array(1:nei_lim, :, :, :) = cg_nei_cg
                                                call move_alloc(intpo_neighbor_array, cg_nei_cg)

                                            end if

                                            cg_nei_cg(n_cg(i), jep, ie, i) = jatomap
                                        end if
                                    end do
                                end if
                            end if
                        end do

                        !Now loop over all real atoms, we only need ghost atom neighbors here
                        do jatom = 1, num_cell_atom(jce)
                            ja = cell_atom(jatom, jce)
                            if(ja > atom_num_l) then 
                                rk(:) = r_atom(:, ja)

                                rlk = rk - rl
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < rc_neisq) then
                                    do i= 1, nei_lists
                                        if(btest(types_to_pot_type(type_atomap(iatomap), type_atom(ja)), &
                                                 neilist_to_pot(i))) then
                                            n_at(i) = n_at(i) + 1
                                            !Check to see if we need to resize arrays
                                            if(n_at(i) > size(cg_nei_at,1)) then
                                                nei_lim = size(cg_nei_at,1)
                                                allocate(intpo_neighbor_array(nei_lim+20,max_basisnum*max_intpo_num,ele_num_l, &
                                                nei_lists), stat = allostat)

                                                if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", &
                                                                                    allostat)

                                                intpo_neighbor_array(:, :, :, :) = 0
                                                intpo_neighbor_array(1:nei_lim, :, :, :) = cg_nei_at
                                                call move_alloc(intpo_neighbor_array, cg_nei_at)

                                            end if

                                            cg_nei_at(n_at(i), jep, ie, i) = ja
                                        end if
                                    end do
                                end if
                            end if
                        end do
                    end do

                    n_cg_at(jep, ie,:) = n_at
                    n_cg_cg(jep, ie,:) = n_cg

                    if((sum(n_at+n_cg))==0) then
                        print *, 'Warning: Integration point', iep, ' of element', ie, &
                                 ' in cell', ice, ' of rank', rank, ' does not have neighbor'
                        print *, "Integration point has position ", rl, " for atomap ", iatomap 
                        print *, "Nodes are at position:"
                        do i = 1, 8
                            print *, r(:, 1, cg_node(i, ie))
                        end do
                        print *, "Box bounds are: ", box_bd
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if
                end if
            end do
        end do
    end do

    return
end subroutine neighbor_list_cg

subroutine neighbor_list_at
    !Build neighbor list for atomistic
    integer :: ia, ja, ice, nei_ice, jce, jatom, jatomap, n_at(nei_lists), n_cg(nei_lists), i
    real(kind=wp) :: rsq, rc_neisq
    real(kind = wp), dimension(3) :: rl, rk, rlk
    integer, allocatable :: atom_neighbor_array(:, :, :)

    !Initialize arrays
    at_nei_at = 0
    n_at_at = 0
    at_nei_cg = 0 
    n_at_cg = 0

    rc_neisq = rc_neigh*rc_neigh

    do ia = 1, atom_num_l

        rl(:) = r_atom(:, ia)
        ice = which_cell_atom(ia)
        n_at = 0
        n_cg = 0

        !Loop over all neighboring cells
        do nei_ice = 1, num_neighbor_cell(ice)
            jce = cell_neighbor(nei_ice, ice)

            if((jce < 1).or.(jce > cell_num)) then
              print *, 'Error: Wrong neighboring cell index', jce, &
                       ' which should be between 1 and', cell_num
              call mpi_abort(mpi_comm_world, 1, ierr)
            end if

            !Loop over all atomaps in neighboring cells
            do jatom = 1, num_cell_atomap(jce)

                jatomap = cell_atomap(jatom, jce)
                rk(:) = r_atomap(:, jatomap)

                rlk = rk - rl
                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                if(rsq < rc_neisq) then

                    do i= 1, nei_lists
                        if(btest(types_to_pot_type(type_atom(ia), type_atomap(jatomap)), neilist_to_pot(i))) then
                            n_cg(i) = n_cg(i) + 1
                            !Check to see if we need to resize arrays
                            if(n_cg(i) > size(at_nei_cg,1)) then
                                nei_lim = size(at_nei_cg,1)
                                allocate(atom_neighbor_array(nei_lim+20, atom_num_l, nei_lists), stat = allostat)

                                if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", allostat)

                                atom_neighbor_array(:, :, :) = 0
                                atom_neighbor_array(1:nei_lim, :, :) = at_nei_cg
                                call move_alloc(atom_neighbor_array, at_nei_cg)

                            end if

                            at_nei_cg(n_cg(i), ia, i) = jatomap
                            if(jatomap<=atomap_num_l) needed_atomap(jatomap) = .true.
                        end if
                    end do
                end if
            end do

            !Loop over all atoms in neighboring cells
            do jatom = 1, num_cell_atom(jce)
                ja = cell_atom(jatom, jce)

                !Only save if ia < ja
                if(ia < ja) then

                    rk(:) = r_atom(:, ja)

                    rlk = rk - rl
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < rc_neisq) then

                        !Add neighbor to list
                        do i= 1, nei_lists
                            if(btest(types_to_pot_type(type_atom(ia), type_atom(ja)), neilist_to_pot(i))) then
                                n_at(i) = n_at(i) + 1
                                !Check to see if we need to resize arrays
                                if(n_at(i) > size(at_nei_at,1)) then
                                    nei_lim = size(at_nei_at,1)
                                    allocate(atom_neighbor_array(nei_lim+20, atom_num_l, nei_lists), stat = allostat)

                                    if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", allostat)

                                    atom_neighbor_array(:, :, :) = 0
                                    atom_neighbor_array(1:nei_lim, :, :) = at_nei_at
                                    call move_alloc(atom_neighbor_array, at_nei_at)

                                end if

                                at_nei_at(n_at(i), ia, i) = ja
                            end if
                        end do
                    end if
                end if
            end do
        end do

        n_at_at(ia,:) = n_at
        n_at_cg(ia,:) = n_cg
    end do

    return
end subroutine neighbor_list_at

subroutine neighbor_list_atomap
    !Build neighbor list for virtual atoms, this is only used in case we want to fuully calculate the cluster potential correctly
    !Instead of imposing uniform electron density on along element subdomains.
    integer :: ia, ja, ice, nei_ice, jce, jatom, jatomap, n_at(nei_lists), n_cg(nei_lists), counts, i
    real(kind=wp) :: rsq, rc_neisq
    real(kind = wp), dimension(3) :: rl, rk, rlk
    integer, allocatable :: atom_neighbor_array(:, :, :)
    character(len=read_len) :: msg

    !Initialize arrays
    atomap_nei_cg = 0 
    atomap_nei_at = 0
    n_atomap_cg = 0
    n_atomap_at = 0

    rc_neisq = rc_neigh*rc_neigh

    do ia = 1, atomap_num_l

        !Assume that an atomap that is sent as a ghost is needed, may need to change this for speed at some point
        rl(:) = r_atomap(:, ia)
        if(.not.in_block_bd(rl, pro_bd_in)) then 
            needed_atomap(ia) = .true.
        end if
        if(needed_atomap(ia)) then 
            ice = which_cell_atomap(ia)
            n_at = 0
            n_cg = 0

            !Loop over all neighboring cells
            do nei_ice = 1, num_neighbor_cell(ice)
                jce = cell_neighbor(nei_ice, ice)

                if((jce < 1).or.(jce > cell_num)) then
                  print *, 'Error: Wrong neighboring cell index', jce, &
                           ' which should be between 1 and', cell_num
                  call mpi_abort(mpi_comm_world, 1, ierr)
                end if

                !Loop over all atomaps in neighboring cells
                do jatom = 1, num_cell_atomap(jce)

                    jatomap = cell_atomap(jatom, jce)
                    rk(:) = r_atomap(:, jatomap)

                    rlk = rk - rl
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < rc_neisq) then

                        !Add neighbor to list
                        do i= 1, nei_lists
                            if(btest(types_to_pot_type(type_atomap(ia), type_atomap(jatomap)), neilist_to_pot(i))) then
                                n_cg(i) = n_cg(i) + 1
                                !Check to see if we need to resize arrays
                                if(n_cg(i) > size(atomap_nei_cg,1)) then
                                    nei_lim = size(atomap_nei_cg,1)
                                    allocate(atom_neighbor_array(nei_lim+20, atom_num_l, nei_lists), stat = allostat)

                                    if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", allostat)

                                    atom_neighbor_array(:, :, :) = 0
                                    atom_neighbor_array(1:nei_lim, :, :) = atomap_nei_cg
                                    call move_alloc(atom_neighbor_array, atomap_nei_cg)

                                end if

                                atomap_nei_cg(n_cg(i), ia, i) = jatomap
                            end if
                        end do
                    end if
                end do

                !Now loop over all real atoms, we only need ghost atom neighbors here
                do jatom = 1, num_cell_atom(jce)
                    ja = cell_atom(jatom, jce)
                    if(ja > atom_num_l) then 
                        rk(:) = r_atom(:, ja)

                        rlk = rk - rl
                        rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                        if(rsq < rc_neisq) then

                            !Add neighbor to list
                            do i= 1, nei_lists
                                if(btest(types_to_pot_type(type_atomap(ia), type_atomap(jatomap)), neilist_to_pot(i))) then
                                    n_at(i) = n_at(i) + 1
                                    !Check to see if we need to resize arrays
                                    if(n_at(i) > size(atomap_nei_at,1)) then
                                        nei_lim = size(atomap_nei_at,1)
                                        allocate(atom_neighbor_array(nei_lim+20, atom_num_l, nei_lists), stat = allostat)

                                        if(allostat /= 0) call alloc_error("Failure to allocate atom_neighbor/ac_array", allostat)

                                        atom_neighbor_array(:, :, :) = 0
                                        atom_neighbor_array(1:nei_lim, :, :) = atomap_nei_at
                                        call move_alloc(atom_neighbor_array, atomap_nei_at)

                                    end if

                                    atomap_nei_at(n_at(i), ia, i) = jatomap
                                end if
                            end do
                        end if
                    end if
                end do

                n_atomap_cg(ia, :) = n_cg
                n_atomap_at(ia, :) = n_at
            end do
        end if
    end do
    
    call mpi_reduce(count(needed_atomap), counts, 1, mpi_integer, mpi_sum, root, world, ierr)
    write(msg, *) "Calculating electron density for ", counts, " out of ", atomap_num, " virtual atoms"
    call log_msg(msg)

    return
end subroutine neighbor_list_atomap

subroutine update_neighbor(time, arg1, arg2)
    !This subroutine updates the neighbor lists when needed
    integer, intent(in) :: time
    !Arg 1 lets us know if we are calling this from the conjugate gradient code. Special code needs to be called when running cg
    logical, intent(in), optional :: arg1, arg2
    integer :: ie, i, ip, info(3), ibasis, ia, inod
    logical :: need_updateall, need_updateme, cg, force_reneighbor
    real(kind=wp) :: center_mass_disp(3), r_old(3), rc_binsq, rlk(3), rsq, t_start, t_end
                    
    character(len=read_len) :: msg

    !Start timer
    t_start = mpi_wtime()

    if(present(arg1)) then 
        cg = arg1
    else 
        cg = .false.
    end if
    
    if(present(arg2)) then 
        force_reneighbor = arg2
    else
        force_reneighbor = .false.
    end if

    !Reneighbored lets us know if we moved atoms/nodes between processors. This is important as some data
    !structures may need to be reset if this happens
    reneighbored = .false.

    !First check if we actually need to update the neighbor list by using the verlet list method
    need_updateall = .false.
    need_updateme = .false.
    
    !Get the movement of the centroid
    call get_centroid(center_mass)
    center_mass_disp = center_mass - center_mass_ori

    !Update pb_node and r when using periodic boundaries, only needed if we are exchanging

    if((ele_num > 0).and.(periodic)) then 
        do ie = 1, ele_num_l
            do inod = 1, ng_node(etype(ie))
                ip = cg_node(inod, ie)
                do ibasis = 1, basis_num(ie)
                    call cross_pb(r(:, ibasis, ip), info)
                    pb_node(:, ibasis, ip) = pb_node(:, ibasis, ip) + info(:)
                end do
            end do 
        end do
    end if

    if((atom_num > 0).and.periodic) then 
        do ia = 1, atom_num_l
            call cross_pb(r_atom(:, ia),info)
        end do
    end if


    if(force_reneighbor) then 
        need_updateall = .true.
    else if(time-last_reneighbor > delay) then 
        r_old(:) = 0.0_wp
        rc_binsq = rc_bin*rc_bin/4.0_wp
        if (ele_num > 0) then 
            ip_loop: do ip = 1, node_num_l
                do ibasis = 1, basis_num(node_cg(ip))
                    !Check to see if nodse have moved by half the bin distance
                    rlk = r(:,ibasis, ip) - r_update(:, ibasis, ip)
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq > rc_binsq) then 
                        need_updateme = .true.
                        exit ip_loop
                    end if
                end do
            end do ip_loop
        end if

        if(atom_num > 0) then 
            !Only check atoms if the elements don't trigger reneighboring code
            if(.not.need_updateme) then 
                do ia = 1, atom_num_l
                    rlk = r_atom(:,ia) - r_atom_update(:, ia)
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq > rc_binsq) then 
                        need_updateme = .true.
                        exit
                    end if
                end do
            end if
        end if

        !Gather need_update
        call mpi_allreduce(need_updateme, need_updateall, 1, mpi_logical, mpi_lor, world, ierr)
    else 
        need_updateall = .false.
    end if



    ! If we don't have to reneighbor then just  update ghost atom/virtual atom positions
    if(.not.need_updateall) then 
        !only update needed virtual atoms
        call update_virtual_atoms(need_all_step)
        call update_proc_bd(.false.)
        if(ele_num > 0) then 
            !if(tflag) call apply_perturbation
!            call ghost_cg
            call processor_atomap
        end if

        if(atom_num > 0) then 
            call processor_atomistic
        end if
        !Set the updated flag to false
        neis_updated = .false.
    else
        !update all 
        call update_virtual_atoms(.true.)
        !Now check to rescale processor boundaries if shrink-wrapped in that dimension
        call update_proc_bd

        reneighbored = .true.


        !write(msg,*) "Updating neighbor list at", time
        !call log_msg(msg)

        !Update lists
        !Check pro_length
        do i = 1, 3
            if(pro_length(i) < 2.0_wp * rc_neigh) then 
                write(msg, *) "Pro_length along the ", i, " direction ", pro_length(i), " of rank ", rank, &
                              " is smaller than", 2.0_wp*rc_neigh
                call misc_error(msg)
            end if
        end do

        if(ele_num /= 0) call update_neighbor_cg(cg)
        if(atom_num /= 0) call update_neighbor_at(cg)

        !Resize force_arrays
        call resize_force_arrays
!        if(ele_num_l > 0) then 
!            if(node_num_l > size(r_update,3)) then 
!                deallocate(r_update)
!                allocate(r_update(3, max_basisnum, node_num_l), stat = allostat)
!                call alloc_error("Failure to allocate r_update in update_neighbor", allostat)
!            end if
!            r_update(:,:,1:node_num_l) = r(:,:, 1:node_num_l)
!        end if
!
!        if(atom_num_l > 0) then 
!            if(atom_num_l > size(r_atom_update,2)) then 
!                deallocate(r_atom_update)
!                allocate(r_atom_update(3, atom_num_l), stat = allostat)
!                call alloc_error("Failure to allocate r_atom_update in update_neighbor", allostat)
!            end if
!            r_atom_update(:,1:atom_num_l) = r_atom_update(:,1:atom_num_l)
!        end if

        !Now if we are using temperature control, we apply the perturbation here
        !if(tflag) call apply_perturbation

        !Update ghosts
        if(ele_num > 0) call ghost_cg
        if(atom_num > 0) call ghost_at

        !Now update cell lists
        call neighbor_lists

        !Set the update nei flag to true
        neis_updated = .true.
    end if
    t_end = mpi_wtime()
    walltime(1) = walltime(1) + (t_end-t_start)
    return
end subroutine update_neighbor

subroutine update_neighbor_cg(arg1)

    logical, intent(in), optional :: arg1

    logical :: cg
    integer :: i, iatom, inod, ie, je, ke, ip, jp, delete_n, ibasis, send_n, send_n_sum,  &
               seg_num_real, irank, ele_num_ln, pb_in(3, max_basisnum, ng_max_node), n_ints, n_reals, &
               iint, ireal, et, tag, id, esize, bnum, btype(max_basisnum), &
               iatomap, node_num_ln, mask, n_cg, icg
    integer, dimension(pro_num) :: displs, recv_counts, counts
    real(kind=wp) :: r_interp(3,max_basisnum)
    real(kind = wp), dimension(3, max_basisnum, ng_max_node) :: r_nodes, vel_nodes, force_nodes, r0, g, h
    integer, allocatable :: delete_buff(:), delete_array(:), send_buff(:), send_array(:), &
                            recv_array(:), who_ele(:), who_ele_all(:)

    logical, allocatable :: ele_shared(:), who_has_array(:)

    real(kind = wp), allocatable :: send_reals(:), recv_reals(:), force_array(:,:,:), send_cg(:), recv_cg(:)
    integer, allocatable :: send_ints(:), recv_ints(:), integer_buff(:)
    character(len=read_len) :: msg

    if(present(arg1)) then 
        cg = arg1
    else 
        cg = .false.
    end if

    allocate(who_ele(ele_num), ele_shared(ele_num), stat = allostat) 
    if(allostat /= 0) call alloc_error("Failure to allocate who_ele and ele_shared", allostat)

    who_ele(:) = 0
    ele_shared(:) = .false.
    seg_num_real = seg_num

    !This series of loops checks to see whether any part of an element belongs to the current processor
    do ie = 1, ele_num_l

        ke = ele_glob_id(ie)
        ele_shared(ke) = .true.
    
        !Put position and pb node into continuous array
        do inod = 1, ng_node(etype(ie))
            ip = cg_node(inod, ie)
            r_nodes(:,:,inod) = r(:,:,ip)
            if(periodic) then 
                pb_in(:,:,inod) = pb_node(:,:,ip)
            else
                pb_in(:,:,inod) = 0
            end if
        end do
        !Check all atoms to get the count of virtual atoms that are within the processor boundaries
        do iatom = 1, get_virtual_count(etype(ie), size_ele(ie))
            call interp_atom(iatom, size_ele(ie), etype(ie), pb_in, basis_num(ie), r_nodes, r_interp)
            do ibasis = 1, basis_num(ie)
              if(in_block_bd(r_interp(:,ibasis), pro_bd)) then
                who_ele(ke) = who_ele(ke) + 1
              end if
            end do
        end do
    end do

!   who_ele_all, this gathers the number of atomaps in for each element in each processor and then sums them. 
!   This is to ensure that all the interpolated atoms are correctly assigned.
    allocate(who_ele_all(ele_num), stat = allostat)
    if(allostat /= 0) call alloc_error("Failure to allocate who_ele_all", allostat)


    who_ele_all(:) = 0
    call mpi_allreduce(who_ele, who_ele_all, ele_num, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    !Now check to make sure calculated virtual atom number isn't greater than the real virtual atom number
    if(rank == root) then
        if(sum(who_ele_all) > atomap_num) then
            print *, 'Error: The sum of who_ele_all', sum(who_ele_all), &
                     ' should not be larger than atomap_num', atomap_num
            call mpi_abort(mpi_comm_world, 1, ierr)
        end if
    end if

    !delete any element if none of its atomap is within pro_bd
    allocate(delete_buff(seg_num), stat = allostat)

    if(allostat /= 0) call alloc_error("Failure to allocate delete_buff", allostat)

    delete_n = 0
    delete_buff(:) = 0
    seg_num_real = seg_num

    do ie = 1, ele_num_l
        ke = ele_glob_id(ie)
        !If we don't have any virtual atoms and the element doesn't belong to use then add it to the delete list
        if((who_ele(ke) == 0).and.(who_has_ele(ie).eqv..false.)) then

            delete_n = delete_n + 1
            ele_shared(ke) = .false.
            if(delete_n > seg_num_real) then
                allocate(delete_array(seg_num_real+seg_num), stat = allostat)
                if(allostat /= 0) call alloc_error("Failure to allocate delete array in update neighbor cg", allostat)

                delete_array(1:seg_num_real) = delete_buff(:)
                delete_array(seg_num_real+1:) = 0
                call move_alloc(delete_array, delete_buff)
                seg_num_real = seg_num_real + seg_num
            end if
            delete_buff(delete_n) = ie
        end if
    end do

    !prepare send_buff when who_has_ele is true
    allocate(send_buff(seg_num), stat = allostat) 
    if(allostat /= 0) call alloc_error("Failure to allocate send buff", allostat)

    send_n = 0
    send_buff(:) = 0
    seg_num_real = seg_num
    !This loop just counts the number of elements that you information to send for. 
    do ie = 1, ele_num_l
        !If we own the element
        if(who_has_ele(ie)) then
            ke = ele_glob_id(ie)


            if(who_ele_all(ke) > basis_num(ie)*get_virtual_count(etype(ie), size_ele(ie))) then
                print *, 'Error: Who_ele_all', who_ele_all(ke), &
                         ' of global element', ke, ' of rank', rank, &
                         ' should not be larger than', get_virtual_count(etype(ie), size_ele(ie))
                call mpi_abort(mpi_comm_world, 1, ierr)

            !If the total counts is less than the real counts then that means the element is in the bounds of a 
            !processor that doesn't have it so we have to send it 
            else if(who_ele_all(ke) < basis_num(ie)*get_virtual_count(etype(ie), size_ele(ie))) then
                send_n = send_n + 1

                if(send_n > seg_num_real) then

                    allocate(send_array(seg_num_real+seg_num), stat = allostat) 
                    if(allostat /= 0) call alloc_error("Failure to allocate send_array in update_neighbor_cg", allostat)

                    send_array(1:seg_num_real) = send_buff(:)
                    send_array(seg_num_real+1:) = 0
                    call move_alloc(send_array, send_buff)
                    seg_num_real = seg_num_real + seg_num
                end if
                !Add this element to the send list
                send_buff(send_n) = ie
            end if
        end if
    end do

    !Get the total number of sends
    call mpi_allreduce(send_n, send_n_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    !Calculate counts of data to send per element
    n_ints = 6+max_basisnum
    n_reals = 3*max_basisnum*ng_max_node
    if(periodic) then 
        n_ints = n_ints+3*max_basisnum*ng_max_node
    end if
    if(need_vel) then 
        n_reals = n_reals + 3*max_basisnum*ng_max_node
    end if
    if(need_force_pre) then 
        n_reals = n_reals + 3*max_basisnum*ng_max_node
    end if

    !Allocate all send and receive buffs
    allocate(send_ints(send_n*n_ints), send_reals(send_n*n_reals), &
             recv_ints(send_n_sum*n_ints), recv_reals(send_n_sum*n_reals), &
             recv_array(send_n_sum), stat = allostat)
    if(allostat /= 0) call alloc_error("Failure to allocate send/recv_buff", allostat)
    !Initialize allocated variables
    send_ints(:) = 0
    recv_ints(:) = 0
    send_reals(:) = 0.0_wp
    recv_reals(:) = 0.0_wp
    recv_array = 0

    !If needed then initialize cg arrays
    if(cg) then 
        n_cg=3*3*max_basisnum*ng_max_node
        allocate(send_cg(send_n*n_cg), recv_cg(send_n_sum*n_cg))
    end if

    do ie = 1, send_n
        je = send_buff(ie)
        iint = n_ints*(ie-1)+1
        ireal = n_reals*(ie-1)+1

        !Get the nodal variables into sequential array
        do inod = 1, ng_node(etype(je))
            ip = cg_node(inod, je)
            r_nodes(:,:,inod) = r(:,:,ip)
            if (periodic) pb_in(:,:,inod) = pb_node(:,:,ip)
            if(need_vel) vel_nodes(:,:,inod) = vel(:,:,ip)
            if(need_force_pre) force_nodes(:,:,inod) = force_eq_pre(:, :, ip)
        end do

        
        !Pack into to the send arrays
        call pack_ele_neighbor(etype(je), tag_ele(je), ele_glob_id(je), size_ele(je), e_mask(je), basis_num(je), basis_type(:,je), &
                        r_nodes, pb_in, vel_nodes, force_nodes, send_ints(iint:iint+n_ints-1), send_reals(ireal:ireal+n_reals-1))

        if(cg) then 
            !First get the nodal cg variables into a sequential array 
            icg = n_cg*(ie-1)+1
            do inod=1, ng_node(etype(je))
                ip = cg_node(inod,je)
                r0(:,:,inod)=rzeronode(:, :, ip)
                g(:,:,inod)=gnode(:, :, ip)
                h(:,:,inod)=hnode(:, :, ip)
            end do

            !Now pack into send_cg array
            call pack_ele_cg(ng_node(etype(je)), basis_num(je), r0, g, h, send_cg(icg:icg+n_cg-1))

        end if
    end do

    recv_ints(:) = 0
    recv_reals(:) = 0.0_wp

    !If one proc then the recv_arrays are equal to the send_arrays
    if(pro_num == 1) then
        recv_ints = send_ints
        recv_reals = send_reals

    !Otherwise we have to prepare all gatherv
    else
        !Get the number of elements being send by each processor
        recv_counts(:) = 0
        call mpi_allgather(send_n, 1, mpi_integer, recv_counts, 1, mpi_integer, mpi_comm_world, ierr)

        !Get the displacements and the real data count for the integer data
        displs(:) = 0
        do irank = 1, pro_num-1
            displs(irank+1) = displs(irank) + recv_counts(irank)*n_ints
        end do
        counts = recv_counts*n_ints
        !allgatherv ints
        call mpi_allgatherv(send_ints, send_n*n_ints, mpi_integer, &
                            recv_ints, counts, displs, &
                            mpi_integer, world, ierr)

        !Get the displacements and the real data count for the real data
        displs(:) = 0
        do irank = 1, pro_num-1
            displs(irank+1) = displs(irank) + recv_counts(irank)*n_reals
        end do
        counts = recv_counts*n_reals
        !allgatherv real
        call mpi_allgatherv(send_reals, send_n*n_reals, mpi_wp, &
                            recv_reals, counts, displs, &
                            mpi_wp, world, ierr)

        if(cg) then 
            !Get the displacements and the real data count for the real data
            displs(:) = 0
            do irank = 1, pro_num-1
                displs(irank+1) = displs(irank) + recv_counts(irank)*n_cg
            end do
            counts = recv_counts*n_cg
            !allgatherv real
            call mpi_allgatherv(send_cg, send_n*n_cg, mpi_wp, &
                                recv_cg, counts, displs, &
                                mpi_wp, world, ierr)
        end if
    end if

!   Now delete delete_buff from the original array,
!   move the end of the array forward to fill the blank
!   do not delete the element when who_has_ele is true,
!   or new master proc needs to be assigned

    ele_num_ln = ele_num_l
    node_num_ln = node_num_l
    do ie = delete_n, 1, -1

        je = delete_buff(ie)
        node_num_ln = node_num_ln - ng_node(etype(je))
        if(je < ele_num_ln) then

            tag_ele(je) = tag_ele(ele_num_ln)
            size_ele(je) = size_ele(ele_num_ln)
            etype(je) = etype(ele_num_ln)
            e_mask(je) = e_mask(ele_num_ln)
            basis_num(je) = basis_num(ele_num_ln)
            basis_type(:,je) = basis_type(:, ele_num_ln)
            ele_glob_id(je) = ele_glob_id(ele_num_ln)

            !Reassign who has ele
            if(who_has_ele(je).eqv..true.) then
                print *, 'Error: Who_has_ele(je) of je', je, ' must be .false.'
                call mpi_abort(mpi_comm_world, 1, ierr)
            else 
                who_has_ele(je) = who_has_ele(ele_num_ln)
                who_has_ele(ele_num_ln) = .false.
            end if

            do inod = 1, ng_node(etype(je))

                jp = cg_node(inod, je)
                ip = cg_node(inod, ele_num_ln)
                r(:, :, jp) = r(:, :, ip)

                if(periodic.eqv..true.) then
                    pb_node(:, :, jp) = pb_node(:, :, ip)
                end if
                if(need_vel) then
                    vel(:, :, jp) = vel(:, :, ip)
                end if
                if(need_force_pre) then
                    force_eq_pre(:, :, jp) = force_eq_pre(:, :, ip)
                end if
            end do

            if(cg) then 
                do inod = 1, ng_node(etype(je))
                    jp = cg_node(inod, je)
                    ip = cg_node(inod, ele_num_ln)
                    
                    rzeronode(:, :, jp) = rzeronode(:,:, ip)
                    gnode(:, :, jp) = gnode(:,:, ip)
                    hnode(:, :, jp) = hnode(:,:, ip)
                end do
            end if

        else if(je > ele_num_ln) then
                print *, 'Error: je', je, ' should not be larger than', &
                         ' ele_num_ln', ele_num_ln
                call mpi_abort(mpi_comm_world, 1, ierr)
        end if

        ele_num_ln = ele_num_ln - 1
    end do

    if(ele_num_ln /= ele_num_l) then
        reneighbored = .true.
        !write(msg,*)'Rank', rank, ' deletes', ele_num_l - ele_num_ln, ' elements'
        !call log_msg(msg, 1, .true.)
    end if

    !append recv_array to the original array, increase array size when necessary
    je = ele_num_ln
    jp = node_num_ln
    do  ie = 1, send_n_sum
        ireal = n_reals*(ie-1)+1
        iint = n_ints*(ie-1)+1
        !Unpack the current element
        call unpack_ele_neighbor(recv_ints(iint:iint+n_ints-1), recv_reals(ireal:ireal+n_reals-1), et, tag, id, esize, mask, &
                                 bnum, btype, pb_in, r_nodes, vel_nodes, force_nodes)
        
        if(cg) then 
            icg = n_cg*(ie-1)+1
            
            call unpack_ele_cg(ng_node(et), bnum, recv_cg(icg:icg+n_cg-1), r0, g, h)
        end if
        !If ele_shared then we already have this element
        if(ele_shared(id)) then 
            recv_array(ie) = 1

        !Otherwise we need to check to see if we need to add it to our list
        else
            !Loop over all the interpolated atoms
            outloop: do iatom = 1, get_virtual_count(et, esize)
                call interp_atom(iatom, esize, et, pb_in, bnum, r_nodes, r_interp)
                do ibasis = 1, bnum
                    !Check if it's in our boundary, if it is then add it, if it isn't then go to the next one
                    if (in_block_bd(r_interp(:, ibasis), pro_bd)) then 
                        if(recv_array(ie) == 1) then 
                            write(msg, *) "Element ", ie, " has been taken by rank ", rank
                            call misc_error(msg)
                        end if

                        !Add this element 
                        je = je+1
                        !Resize element arrays if needed
                        if(je > ele_num_lr) call grow_cg_arrays(1)

                        !Add new element to arrays
                        recv_array(ie) = 1
                        tag_ele(je) = tag
                        etype(je) = et
                        ele_glob_id(je) = id
                        e_mask(je) = mask
                        size_ele(je) = esize
                        basis_num(je) = bnum
                        basis_type(:, je) = btype
                        
                        !check to see if need to resize node arrays
                        if(jp + ng_node(et) > node_num_lr) call grow_cg_arrays(2)
                        
                        if(need_force_pre) then  
                            if (jp+ng_node(et) > size(force_eq_pre,3)) then 
                                allocate(force_array(3, max_basisnum, node_num_lr), stat = allostat)
                                if(allostat > 0) call alloc_error("Failure to allocate force array in update_nei_cg", allostat)
                                force_array(:,:,1:size(force_eq_pre,3)) = force_eq_pre
                                force_array(:,:,size(force_eq_pre,3)+1:) = 0.0_wp
                                call move_alloc(force_array, force_eq_pre)
                            end if
                        end if
                        !Now assign node arrays
                        do inod = 1, ng_node(et)
                            jp = jp + 1
                            cg_node(inod, je) = jp
                            node_cg(jp) = je
                            r(:,:,jp) = r_nodes(:,:,inod)
                            if(periodic) pb_node(:, :, jp) = pb_in(:, :, inod)
                            if(need_vel) vel(:, :, jp) = vel_nodes(:, :, inod)
                            if(need_force_pre) force_eq_pre(:, :, jp) = force_nodes(:,:,inod)
                        end do
                        
                        if(cg) then 
                            if(jp+ng_node(et) > size(rzeronode,3)) call grow_ele_min_arrays
                            do inod = 1, ng_node(et)
                                jp = cg_node(inod, je)
                                rzeronode(:,:,jp) = r0(:,:,inod)
                                gnode(:,:,jp) = g(:,:,inod)
                                hnode(:,:,jp) = h(:,:,inod)
                            end do
                        end if

                        !If we added an element then we exit the loop over interpolated atoms and move to the next element
                        exit outloop
                    end if
                end do
            end do outloop
        end if
    end do

    if(je > ele_num_ln) then 
        reneighbored = .true.
        !write(msg, *) "Rank ", rank, " adds ", je-ele_num_ln, " elements" 
        !call log_msg(msg,1, .true.)
    end if

    node_num_l = jp
    ele_num_l = je

    !Rebuild itype array now that elements lists have changed
    call update_itype

    !Double check that ele_num_l summed isn't smaller than ele_num
    call mpi_reduce(ele_num_l, i, 1, mpi_integer, mpi_sum, root, mpi_comm_world, ierr)
    if(rank == root) then 
        if( i < ele_num) then 
            write(msg, *) "Total number of elements ", i, " should not be smaller than ele_num ", ele_num
            call misc_error(msg)
        end if
    end if

    !Resize tag_ele_shared, pro_shared_num, who_has_ele
    allocate(integer_buff(ele_num_l), stat=allostat)
    if (allostat>0) call alloc_error("Failure to resize tag_ele_shared in update_nei_cg", allostat)
    integer_buff(:) = 0
    call move_alloc(integer_buff, ele_id_shared)

    allocate(integer_buff(ele_num_l), stat=allostat)
    if (allostat>0) call alloc_error("Failure to resize pro_shared_num in update_nei_cg", allostat)
    integer_buff(:) = 0
    call move_alloc(integer_buff, pro_shared_num)

    if(ele_num_l > size(who_has_ele)) then 
        allocate(who_has_array(ele_num_l), stat = allostat)
        if (allostat > 0) call alloc_error("Failure to allocate who_has_array in update_nei_cg", allostat)
        who_has_array(1:size(who_has_ele)) = who_has_ele(:)
        who_has_array(size(who_has_ele)+1:) = .false.
        call move_alloc(who_has_array, who_has_ele)
    end if

    !Update the shared elements
    who_ele(:) = 0
    do ie = 1,ele_num_l
        ke = ele_glob_id(ie)
        who_ele(ke) = 1
    end do

    who_ele_all(:) = 0
    call mpi_allreduce(who_ele, who_ele_all, ele_num, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    if(rank == root) then 
        if(sum(who_ele_all) < ele_num) then 
            write(msg, *) "The sum of who_ele_all", sum(who_ele_all), " should not be smaller than ", ele_num
            call misc_error(msg)
        end if
    end if
    
    allocate(integer_buff(ele_num), stat=allostat)
    if(allostat > 0) call alloc_error("Failure to allocate integer_buff in update_nei_cg", allostat)

    integer_buff(:) = 0
    je = 0
    do ie = 1, ele_num
        !If who_ele_all(ie) > 1 that means it's a shared element 
        if(who_ele_all(ie) > 1) then 
            je = je + 1
            integer_buff(ie) = je
        end if 
    end do

    ele_shared_num = je
    ele_id_shared=0
    do ie = 1, ele_num_l
        ke = ele_glob_id(ie)
        pro_shared_num(ie) = who_ele_all(ie)
        je = integer_buff(ke)
        if(je /= 0) then 
            ele_id_shared(ie) = je
        end if
    end do

    !Now resize the atomap arrays
    if(ele_num_lr > size(cg_atomap,2)) then 
        deallocate(cg_atomap, stat = allostat)
        allocate(cg_atomap(atomap_max_ele, ele_num_lr), stat = allostat)
        if(allostat > 0) call alloc_error("Failure allocating cg_atomap in update_nei_cg", allostat)
    end if

    who_ele(:) = 0
    r_atomap(:,:) = 0.0_wp
    type_atomap(:) = 0
    cg_atomap(:,:) = 0
    atomap_num_l = 0
    !Loop over all elements 
    do ie = 1, ele_num_l
        ke = ele_glob_id(ie)

        r_nodes = 0.0_wp
        pb_in = 0
        do inod = 1, ng_node(etype(ie))
            ip = cg_node(inod, ie)
            r_nodes(:,:,inod) = r(:,:, ip)
            if(periodic) pb_in(:,:,inod) = pb_node(:,:,ip)
        end do 
        !Loop over all virtual atoms 
        do iatom = 1, get_virtual_count(etype(ie), size_ele(ie))
            call interp_atom(iatom, size_ele(ie), etype(ie), pb_in, basis_num(ie), r_nodes, r_interp)
            do ibasis = 1, basis_num(ie)
                !If it is in our boundaries then add it to the list
                if(in_block_bd(r_interp(:,ibasis),pro_bd)) then
                    atomap_num_l = atomap_num_l + 1
                    if (atomap_num_l > atomap_num_lr) call grow_cg_arrays(3)
                    
                    who_ele(ke) = who_ele(ke) + 1
                    cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie) = atomap_num_l
                    r_atomap(:, atomap_num_l) = r_interp(:, ibasis)
                    type_atomap(atomap_num_l) = basis_type(ibasis,ie)
                    
                !Debug statement to check whether it is in the box at all
                else if(.not. in_block_bd(r_interp(:,ibasis), box_bd)) then 
                    print *, "Interpolated atom ", iatom, " of element ", ie, " on rank ", rank, "is outside box boundaries"
                    print *, rank, ie, cg_node(1, ie), r(:, 1, cg_node(1, ie))
                    call mpi_abort(mpi_comm_world, 1, ierr)
                    
                end if
            end do
        end do

    end do

    !Now check to make sure all the values are correct
    if(atomap_num_l > atomap_num_lr) call misc_error("Atomap_num_lr can't be greater than atomap_num_l in update_neigh_cg")

    who_ele_all(:) = 0
    call mpi_allreduce(who_ele, who_ele_all, ele_num, mpi_integer, mpi_sum, world, ierr)
    !Check to make sure all the atomap counts are right
    do ie = 1, ele_num_l
        !If the total counts are greater than the real counts then exit with error
        if(who_has_ele(ie)) then 
            ke = ele_glob_id(ie)
            if(who_ele_all(ke) /= basis_num(ie)*get_virtual_count(etype(ie), size_ele(ie))) then 
                write(msg, *) "Who_ele_all ", who_ele_all(ke), " of global element ", ke, " should be ", &
                               basis_num(ie)*get_virtual_count(etype(ie), size_ele(ie))
                call misc_error(msg)
            end if
        end if
    end do

    if(rank == root) then 
        if(sum(who_ele_all) /= atomap_num) then 
            write(msg, *) "The sum of who_ele_all ", sum(who_ele_all), " should equal atomap_num ", atomap_num
            call misc_error(msg)
        end if
    end if

    call mpi_reduce(atomap_num_l, i, 1, mpi_integer, mpi_sum, root, world, ierr)
    if(rank == root) then 
        if(i /= atomap_num) then 
            write(msg, *) "Total atomap number ", i, " should equal ", atomap_num
            call misc_error(msg)
        end if
    end if

    !Now assign new tags
    recv_counts(:) = 0
    call mpi_allgather(atomap_num_l, 1, mpi_integer, recv_counts, 1, mpi_integer, mpi_comm_world, ierr)
    displs(:) = 0
    do irank = 1, pro_num-1
        displs(irank+1) = displs(irank) + recv_counts(irank)
    end do
    do iatomap = 1, atomap_num_l
        tag_atomap(iatomap) = displs(rank+1)+iatomap
    end do

    !Resize integration point arrays
    if(ele_num_l > size(who_has_intpo, 2)) then 
        deallocate(who_has_intpo, which_cell_intpo, stat=allostat)
        allocate(who_has_intpo(max_intpo_num, ele_num_l), which_cell_intpo(max_intpo_num, ele_num_l), stat = allostat)
        if (allostat > 0) call alloc_error("Failure to allocate/deallocate who_has_intpo in update_nei_cg", allostat)
    end if

    !update integration points
    call update_intpo
    return
end subroutine update_neighbor_cg

subroutine update_neighbor_at(arg1)
    !Transfer atoms between processors
    logical, intent(in), optional :: arg1

    logical :: cg
    integer :: i, ix, iy, iz, ia, ja, send_n, send_n_sum, &
               seg_num_real, irank, atom_num_r, atom_num_ln, n_ints, n_reals, iint, ireal, &
               tag, typ, mask, n_cg, icg
    real(kind=wp) :: rtemp(3), veltemp(3), forcetemp(3)
    integer, dimension(3) :: index_array
    integer, dimension(pro_num) :: displs, recv_counts
    real(kind = wp), dimension(3) :: r_in, vel_in, force_in, r0, g, h
    integer, allocatable :: send_buff(:), send_array(:), &
                            recv_array(:), &
                            send_ints(:), recv_ints(:)
    real(kind = wp), allocatable :: send_reals(:), recv_reals(:), force_array(:,:), send_cg(:), recv_cg(:)
    character(len=read_len) :: msg

!   update
!   r_atom, tag_atom, grain_atom
!   vel_atom (when vel_now.eqv..true.)
!   force_atom_pre (when force_pre_now.eqv..true.)
!   group_atom (when group_num /= 0)

    if(present(arg1)) then 
        cg = arg1
    else 
        cg = .false.
    end if

    allocate(send_buff(seg_num), stat = allostat)
    if(allostat /= 0) call alloc_error("Failure to allocate send_buff in update_nei_at", allostat)

    send_n = 0
    send_buff(:) = 0
    seg_num_real = seg_num

    !Figure out which atoms need to be sent
    do ia = 1, atom_num_l
        r_in(:) = r_atom(:, ia)
        if(.not.in_block_bd(r_in, pro_bd)) then

            send_n = send_n + 1
            if(send_n > seg_num_real) then

                allocate(send_array(seg_num_real+seg_num), stat = allostat) 
                if(allostat /= 0) call alloc_error("Failure to allocate send_array ", allostat)

                send_array(1:seg_num_real) = send_buff(:)
                send_array(seg_num_real+1:) = 0
                call move_alloc(send_array, send_buff)

                seg_num_real = seg_num_real + seg_num
            end if
            send_buff(send_n) = ia
        end if
    end do

    call mpi_allreduce(send_n, send_n_sum, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)

    !Calculate size of data
    n_ints = 3
    n_reals = 3
    if(need_vel) n_reals = n_reals + 3
    if(need_force_pre) n_reals = n_reals + 3

    !Allocate send and recv arrays
    allocate(send_ints(n_ints*send_n), recv_ints(n_ints*send_n_sum), &
             send_reals(n_reals*send_n), recv_reals(n_reals*send_n_sum), &
             recv_array(send_n_sum), stat = allostat)
    send_ints(:) = 0
    recv_ints(:) = 0
    send_reals(:) = 0.0_wp
    recv_reals(:) = 0.0_wp

    !If running cg then allocate cg data arrays
    if(cg) then 
        n_cg = 9
        allocate(send_cg(n_cg*send_n), recv_cg(n_cg*send_n_sum))
        send_cg = 0
        recv_cg=0
    end if
    
    !Pack atoms
    do ia = 1, send_n
        ja = send_buff(ia)
        if(need_vel) vel_in = vel_atom(:,ja)
        if(need_force_pre) force_in = force_atom_pre(:,ja)
    
        iint = n_ints*(ia-1)+1
        ireal = n_reals*(ia-1)+1
        call pack_atom_neighbor(tag_atom(ja), type_atom(ja), a_mask(ja), r_atom(:,ja), vel_in, force_in, &
                                send_ints(iint:iint+n_ints-1), send_reals(ireal:ireal+n_reals-1))
        if(cg) then
            icg = n_cg*(ia-1)+1
            call pack_atom_cg(rzeroatom(:,ja), gatom(:, ja), hatom(:,ja), send_cg(icg:icg+n_cg-1))
        end if
    end do

    !If one proc then the send arrays are equal to the recv arrays    
    if(pro_num == 1) then
        recv_ints = send_ints
        recv_reals = send_reals
        if(cg) recv_cg=send_cg
    !Otherwise we have to mpi_allgatherv the data
    else

        recv_counts(:) = 0
        call mpi_allgather(send_n, 1, mpi_integer, recv_counts, 1, mpi_integer, mpi_comm_world, ierr)

        !Get displacement and offset for integer data
        displs(:) = 0
        do irank = 1, pro_num-1
            displs(irank+1) = displs(irank) + n_ints*recv_counts(irank)
        end do
        !allgatherv ints
        call mpi_allgatherv(send_ints, n_ints*send_n, mpi_integer, &
                            recv_ints, n_ints*recv_counts, displs, &
                            mpi_integer, mpi_comm_world, ierr)

        !Get displacement and offset for real data
        displs(:) = 0
        do irank = 1, pro_num-1
            displs(irank+1) = displs(irank) + n_reals*recv_counts(irank)
        end do

        !allgatherv reals
        call mpi_allgatherv(send_reals, n_reals*send_n, mpi_wp, &
                            recv_reals, n_reals*recv_counts, displs, &
                            mpi_wp, mpi_comm_world, ierr)
        if(cg) then
            !Get displacement and offset for cg data
            displs(:) = 0
            do irank = 1, pro_num-1
                displs(irank+1) = displs(irank) + n_cg*recv_counts(irank)
            end do

            !allgatherv reals
            call mpi_allgatherv(send_cg, n_cg*send_n, mpi_wp, &
                                recv_cg, n_cg*recv_counts, displs, &
                                mpi_wp, mpi_comm_world, ierr)
        end if
            
    end if

    !delete send_buff from the original array,
    !move the end of the array to fill the blank 
    atom_num_ln = atom_num_l
    do ia = send_n, 1, -1

        ja = send_buff(ia)

        if(ja < atom_num_ln) then

            tag_atom(ja) = tag_atom(atom_num_ln)
            type_atom(ja) = type_atom(atom_num_ln)
            a_mask(ja) = a_mask(atom_num_ln)
            r_atom(:, ja) = r_atom(:, atom_num_ln)
            if(need_vel) vel_atom(:, ja) = vel_atom(:, atom_num_ln)
            if(need_force_pre) force_atom_pre(:, ja) = force_atom_pre(:, atom_num_ln)

            if(cg) then 
                rzeroatom(:,ja) = rzeroatom(:,atom_num_ln)
                gatom(:,ja) = gatom(:, atom_num_ln)
                hatom(:,ja) = hatom(:, atom_num_ln)
            end if

        else if(ja > atom_num_ln) then

          print *, 'Error: ja', ja, ' should not be larger than', &
                   ' atom_num_ln', atom_num_ln
          call mpi_abort(mpi_comm_world, 1, ierr)

        end if
        atom_num_ln = atom_num_ln - 1
    end do

    if(atom_num_ln /= atom_num_l) then
        reneighbored = .true.
    end if

    !   Now add new atoms to the end of the array
    ja = atom_num_ln
    atom_num_r = atom_num_l
    do ia = 1, send_n_sum
        ireal = n_reals*(ia-1)+1
        iint = n_ints*(ia-1)+1
        !Unpack the current atom
        call unpack_atom_neighbor(recv_ints(iint:iint+n_ints-1), recv_reals(ireal:ireal+n_reals-1), tag, typ, mask, rtemp, &
                                  veltemp, forcetemp)   
        if(cg) then
            icg = n_cg*(ia-1)+1
            call unpack_atom_cg(recv_cg(icg:icg+n_cg-1), r0, g, h)
        end if
        do ix = 0, 2
            do iy = 0, 2
                inloop: do iz = 0, 2

                    index_array(:) = [ ix, iy, iz ]
                    do i = 1, 3
                        if(index_array(i) == 0) then
                            r_in(i) = rtemp(i)

                        else if(period(i).eqv..true.) then
                            r_in(i) = rtemp(i) + (-1) ** index_array(i) * box_length(i)

                        else
                          exit inloop
                        end if
                    end do

                !Add the atom if needed
                    if(in_block_bd(r_in, pro_bd)) then
                        recv_array(ia) = 1

                        ja = ja + 1
                        !Resize arrays if needed
                        if (ja > atom_num_lr) call grow_at_arrays
                            
                        if(need_force_pre) then 
                            if(ja > size(force_atom_pre,2)) then 
                                allocate(force_array(3,atom_num_lr), stat = allostat)
                                if (allostat > 0) call alloc_error("Failure to allocate force array in update_neighbor_at", &
                                                                    allostat)
                                force_array(:,1:ja-1) = force_atom_pre(:,:)
                                force_array(:,ja:) = 0.0_wp
                                call move_alloc(force_array, force_atom_pre)
                            end if
                        end if
                        tag_atom(ja) = tag
                        type_atom(ja) = typ
                        a_mask(ja) = mask
                        r_atom(:,ja) = r_in
                        if(need_vel) vel_atom(:,ja) = veltemp
                        if(need_force_pre) then 
                            force_atom_pre(:,ja) = forcetemp
                        end if
                        
                        if(cg) then 
                            if(ja > size(gatom,2)) call grow_at_min_arrays
                            rzeroatom(:,ja) = r0
                            gatom(:,ja) = g
                            hatom(:, ja) = h
                        end if
                    end if
                end do inloop
            end do
         end do
    end do

    if(ja > atom_num_ln) then
        reneighbored = .true.
        !write(msg, *)'Rank', rank, ' adds', ja-atom_num_ln, ' atoms'
        !call log_msg(msg, 1, .true.)
    end if
    atom_num_l = ja

    !debug
    if(atom_num_l > atom_num_lr) then
        print *, 'Rank', rank, ' has atom_num_l', atom_num_l, &
                 ' which is larger than atom_num_lr', atom_num_lr
        call mpi_abort(mpi_comm_world, 1, ierr)
    end if

    !debug
    call mpi_reduce(atom_num_l, i, 1, mpi_integer, mpi_sum, root, mpi_comm_world, ierr)
    if(rank == root) then
        if(i /= atom_num) then
            write(msg, *) 'Total number of atoms', i, ' should equal atom_num', atom_num
            call misc_error(msg)
        end if
    end if

   return
end subroutine update_neighbor_at

subroutine get_centroid(centroidall)
    !This subroutine gets the average position of all virtual and real atoms
    real(kind=wp), intent(out) :: centroidall(3)
    integer :: ie, inod, ip, ia, ibasis, atomp_num, eanum
    real(kind=wp) :: centroidme(3)

    centroidme(:) = 0.0_wp
    centroidall(:) = 0.0_wp
    atomp_num = atom_num + atomap_num 
    !Get center of element positions
    if(ele_num > 0) then 
        do ie = 1, ele_num_l
            if (who_has_ele(ie)) then 
                select case(etype(ie))
                case(1,2)
                    eanum = basis_num(ie)*(size_ele(ie)+1)**3
                end select
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod,ie)
                    do ibasis = 1, basis_num(ie)
                        centroidme(:) = centroidme + r(:,ibasis,ip)*eanum/(ng_node(etype(ie))*atomp_num)
                    end do
                end do
            end if
        end do
    end if

    !Get atomistic centroid
    if(atom_num > 0) then 
        do ia = 1, atom_num_l
            centroidme(:) = centroidme(:) + r_atom(:, ia)/atomp_num
        end do
    end if

    !Sum this over all processors
    call mpi_allreduce(centroidme, centroidall, 3, mpi_wp, mpi_sum, world, ierr)

    return
end subroutine get_centroid

subroutine pack_ele_neighbor(send_etype, send_tag_ele, send_id, send_esize, mask, send_basis_num, send_basis_type, &
                    send_r_nodes, send_pb_nodes, send_vel_nodes, send_force_nodes, send_int,send_real)
    !This subroutine packs all data for one element into 2 arrays, a send_int array and a send_real array 
    !containing all information needed for elements

    integer, intent(in) :: send_etype, send_tag_ele, send_esize, send_id, mask, send_basis_num, send_basis_type(max_basisnum),&
                           send_pb_nodes(3,max_basisnum,ng_max_node)
    real(kind=wp), dimension(3, max_basisnum, ng_max_node), intent(in) :: send_r_nodes, send_vel_nodes, send_force_nodes
    integer, dimension(:), intent(out) :: send_int
    real(kind=wp), dimension(:), intent(out) :: send_real

    integer i, j, inod, ibasis

    logical :: pack_vel
    pack_vel=.false.
    
    !Initialize variables
    send_int(:) = 0  
    send_real(:) = 0.0_wp

    !Calculate send counts
    !First pack the send_int variable
    send_int(1) = send_etype
    send_int(2) = send_tag_ele
    send_int(3) = send_id
    send_int(4) = send_esize
    send_int(5) = mask
    send_int(6) = send_basis_num
    j = 6+send_basis_num
    send_int(7:j) = send_basis_type(1:send_basis_num)
    if(periodic) then 
        do inod = 1, ng_node(send_etype)
            do ibasis = 1, send_basis_num
                do i = 1, 3
                    j = j+1
                    send_int(j) = send_pb_nodes(i, ibasis, inod)
                end do
            end do
        end do
    end if
    !Now pack the send_real variable
    j = 1
    do inod = 1, ng_node(send_etype)
        do ibasis = 1, send_basis_num
            do i =1, 3
                send_real(j) = send_r_nodes(i, ibasis, inod) 
                j = j+1
                if (need_vel) then 
                    send_real(j) = send_vel_nodes(i, ibasis,inod)
                    j = j+1
                end if
                if(need_force_pre) then 
                    send_real(j) = send_force_nodes(i, ibasis, inod)
                    j=j+1
                end if
            end do
        end do 
    end do
    
    return
end subroutine pack_ele_neighbor

subroutine unpack_ele_neighbor(recv_int, recv_real, recv_etype, recv_tag_ele, recv_id, recv_esize, mask, recv_basis_num, &
                      recv_basis_type,  recv_pb_node, recv_r_nodes, recv_vel_nodes, recv_force_nodes)
    !This subroutine unpacks the arrays that are communicated 

    integer, dimension(:), intent(in) :: recv_int
    real(kind=wp), dimension(:), intent(in) :: recv_real

    integer, intent(out) :: recv_etype, recv_tag_ele, recv_id, recv_esize, mask, recv_basis_num, recv_basis_type(max_basisnum),&
                            recv_pb_node(3,max_basisnum, ng_max_node)
    real(kind=wp), dimension(3, max_basisnum, ng_max_node), intent(out) :: recv_r_nodes, recv_vel_nodes, recv_force_nodes
    

    integer i, j, inod, ibasis

    recv_basis_type(:)    = 0
    recv_pb_node(:,:,:)      = 0
    recv_r_nodes(:,:,:)   = 0.0_wp
    recv_vel_nodes(:,:,:) = 0.0_wp

    !First pack the recv_int variable
    recv_etype = recv_int(1) 
    recv_tag_ele = recv_int(2) 
    recv_id = recv_int(3)
    recv_esize = recv_int(4) 
    mask = recv_int(5)
    recv_basis_num = recv_int(6) 
    j = 6+recv_basis_num
    recv_basis_type(1:recv_basis_num) = recv_int(7:j) 
    if(periodic) then 
        do inod = 1, ng_node(recv_etype)
            do ibasis = 1, recv_basis_num
                do i = 1, 3
                    j = j+1
                    recv_pb_node(i, ibasis, inod) = recv_int(j) 
                end do
            end do
        end do
    end if
    !Now pack the recv_real variable
    j = 1
    do inod = 1,ng_node(recv_etype)
        do ibasis = 1, recv_basis_num
            do i =1, 3
                recv_r_nodes(i, ibasis, inod) = recv_real(j) 
                j = j+1
                if(need_vel) then 
                    recv_vel_nodes(i, ibasis,inod) = recv_real(j) 
                    j = j+1
                end if
                if(need_force_pre) then 
                    recv_force_nodes(i, ibasis, inod) = recv_real(j) 
                    j = j+1
                end if
            end do
        end do 
    end do
    
    return
end subroutine unpack_ele_neighbor

subroutine pack_atom_neighbor(tag_buff, type_buff, mask, r_buff, vel_buff, force_buff, send_int, send_real)
    !This subroutine packs the atom information into send_int and send_real
    integer, intent(in) :: tag_buff, type_buff, mask
    real(kind=wp), dimension(3), intent(in) :: r_buff
    real(kind=wp), intent(in), optional :: vel_buff(3)
    real(kind=wp), intent(in), optional :: force_buff(3)
    integer, dimension(:), intent(out) :: send_int
    real(kind=wp), dimension(:), intent(out) :: send_real

    integer :: i, k

    !Check to make sure vel_buff is provided if need_vel
    !First pack send_int
    send_int(1) = tag_buff
    send_int(2) = type_buff
    send_int(3) = mask

    !Now pack send_real
    k = 1 
    do i = 1,3
        send_real(k) = r_buff(i)
        k=k+1
        if(need_vel) then
            send_real(k) = vel_buff(i)
            k=k+1
        end if
        if(need_force_pre) then 
            send_real(k) = force_buff(i)
            k = k + 1
        end if
    end do

    return

end subroutine pack_atom_neighbor

subroutine unpack_atom_neighbor(recv_int, recv_real, tag_buff, type_buff, mask, r_buff, vel_buff, force_buff )
    !This subroutine packs the atom information into recv_int and recv_real
    integer, dimension(:), intent(in) :: recv_int
    real(kind=wp), dimension(:), intent(in) :: recv_real
    integer, intent(out) :: tag_buff, type_buff, mask
    real(kind=wp), dimension(3), intent(out) :: r_buff
    real(kind=wp), intent(out) :: vel_buff(3)
    real(kind=wp), intent(out) :: force_buff(3)

    integer :: i, k

    !First unpack recv_int
    tag_buff = recv_int(1) 
    type_buff= recv_int(2) 
    mask = recv_int(3)

    !Now unpack recv_real
    k = 1 
    do i = 1,3
         r_buff(i) = recv_real(k)
        k=k+1
        if(need_vel) then
            vel_buff(i) = recv_real(k) 
            k=k+1
        end if
        if(need_force_pre) then 
            force_buff(i) = recv_real(k) 
            k = k + 1
        end if
    end do
    return
end subroutine unpack_atom_neighbor


subroutine log_neighbor_info
    character(len=read_len) :: msg

    write(msg, *) "Neighbor list was built ", builds, " times in last run"
    call log_msg(msg)
    builds = 0 
end subroutine log_neighbor_info

end module neighbors
