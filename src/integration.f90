module integration
    !This module contains the code required for integration point setup and virtual atom data
    use parameters
    use comms
    use elements
    use errors
    
    implicit none

    integer :: max_intpo_num, intpo_num, intpo_num_l
    integer, save :: atomap_max_ele,  etype_count, etype_to_itype(defined_element_types), &
                     intpo_nums(4, defined_element_types), intpo_count(defined_element_types)
    integer, allocatable, save :: atom_intpo(:,:,:), who_rep_atomap(:,:,:)
    real(kind = wp), allocatable, save :: weight_intpo(:,:,:)
    real(kind = wp), allocatable, save ::  mass_mat_inv(:,:)
    character(len = 100) :: mass_mat_ty
    logical, allocatable, save :: who_has_intpo(:,:)

    !Itype is the integration array type for the elements
    integer, allocatable, save :: itype(:)
    integer, parameter :: intpo_owner_skip = 0, intpo_owner_register = 1, intpo_owner_error = 2

    public
    contains

    subroutine integration_defaults
        mass_mat_ty = 'lumped'
    end subroutine integration_defaults

    subroutine init_integration
        !Initialize the integration
        
        integer :: i, j, ie, iep, local_iep(2), global_iep(2)
        real(kind=wp) :: weight_temp
        real(kind=wp) :: v_inv_matrix(8**2)

        !First set up the mass arrays
        !Allocate inverse mass arays
        allocate(mass_mat_inv(8,8))

        !Now calculate the inverse mass matrix, this is just for FCC 8 node elements
        v_inv_matrix = [ 8.0_wp, -4.0_wp, 2.0_wp, -4.0_wp, -4.0_wp, 2.0_wp, -1.0_wp, 2.0_wp, &
                        -4.0_wp, 8.0_wp, -4.0_wp, 2.0_wp, 2.0_wp, -4.0_wp, 2.0_wp, -1.0_wp, &
                         2.0_wp, -4.0_wp, 8.0_wp, -4.0_wp, -1.0_wp, 2.0_wp, -4.0_wp, 2.0_wp, &
                        -4.0_wp, 2.0_wp, -4.0_wp, 8.0_wp, 2.0_wp, -1.0_wp, 2.0_wp, -4.0_wp, &
                        -4.0_wp, 2.0_wp, -1.0_wp, 2.0_wp, 8.0_wp, -4.0_wp, 2.0_wp, -4.0_wp, &
                         2.0_wp, -4.0_wp, 2.0_wp, -1.0_wp, -4.0_wp, 8.0_wp, -4.0_wp, 2.0_wp, &
                        -1.0_wp, 2.0_wp, -4.0_wp, 2.0_wp, 2.0_wp, -4.0_wp, 8.0_wp, -4.0_wp, &
                         2.0_wp, -1.0_wp, 2.0_wp, -4.0_wp, -4.0_wp, 2.0_wp, -4.0_wp, 8.0_wp ]

        mass_mat_inv(:,:) = reshape(v_inv_matrix, [8,8])
        
        !First communicate the etypes present, this array is only on the root processor as it's defined in read_restart
        call mpi_bcast(etype_present, defined_element_types, mpi_logical, root, world, ierr)

        !Assign the etype_to_itype array, this just maps the etype to the correct index of the itype variables
        !This is because generally all the element types won't be used in one simulation to save on memory
        j = 0
        intpo_nums(:,:) = 0
        intpo_count(:) = 0
        etype_to_itype(:) = 0
        do i = 1, defined_element_types
            if (etype_present(i)) then
                j=j+1
                etype_to_itype(i) = j 
                select case(i) 
                    !Now set up the integration point definitions for FCC elements (Will be expanded to other element types
                    !Intpo_nums has the integration point count for the integration schemes
                    !The first index si the corner integration points, the second is the number of edge integration points, 
                    !the third is the surface integration points, and the 4th is the interior integration points
                case(1)
                    !1NN element integration point count
                    intpo_nums(:,j) = (/ 8, 12, 6, 1 /) 
                    intpo_count(j) = 27
                case(2,3)
                    !2NN element integration point count
                    intpo_nums(:,j) = (/ 64, 48, 12, 1 /)
                    intpo_count(j) = 125
                end select
            end if
        end do
        
        !Calculate max number of intpo nums for elements needed for model
        max_intpo_num = maxval(intpo_count)

        !Allocate arrays
        etype_count = count(etype_present)
        allocate(who_rep_atomap((max_size+1)**3, unique_sizes, etype_count), &
                weight_intpo(max_intpo_num, unique_sizes, etype_count),  &
                atom_intpo(max_intpo_num, unique_sizes, etype_count), &
                stat = allostat)
        if(allostat > 0) call alloc_error("Failure allocating weight arrays in init_integration", allostat)

        !Now initialize the integration points for all of the elements that we need
        !At the moment this code creates the integration points for all possible combinations of esize and etype
        if(etype_present(1)) call init_rhomb(1)
        if (etype_present(2))call init_rhomb(2)
        if (etype_present(3)) call init_rhomb(2, 3)

        !Now update the itypes for all the elements
        call update_itype
        !Now double check to make sure all of the weights for all elements add up correctly and count the number of integration
        !points
        local_iep = 0
        global_iep = 0
        do ie = 1, ele_num_l
            weight_temp = 0.0_wp
            do iep = 1, intpo_count(itype(ie))
                weight_temp = weight_temp + weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
            end do
            if(.not. is_equal(weight_temp, real((size_ele(ie)+1)**3,wp))) then 
                print *, "Error: weight_temp for element ", ie, " should equal ", (size_ele(ie)+1)**3
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
            if(who_has_ele(ie)) then 
                local_iep(1)  = local_iep(1) + intpo_count(itype(ie))
                local_iep(2) = local_iep(2) + intpo_count(itype(ie))*basis_num(ie)
            end if
        end do

        !Sum to get total intpo_num
        call mpi_allreduce(local_iep, global_iep, 2, mpi_integer, mpi_sum, world, ierr)

        !Now multiply by the maximum basis_num to get the right dimensions for the arrays
        intpo_num = global_iep(2)

        !Now update integration point arrays
        call update_intpo

    end subroutine init_integration

    subroutine init_rhomb(intpo_depth, etyp)
        !Initialize one NN rhombohedral elements
        
        !intpo_depth is the integration point depth. The only difference between 1NN and 2NN elements is this integration point
        !depth
        integer, intent(in) :: intpo_depth
        integer, intent(in), optional :: etyp

        integer :: i, j, k, ide, ix, iy, iz, ie, iep, esize, et, interp_count, intpo_set_num, iatom, &
                    i_edge, i_node, i_surf, i_inner
        real(kind=wp) :: node_weight, edge_weight, surf_rect_weight, inner_weight, x, y, z

        !Variables needed for primitive unit cell definitions
        real(kind=wp) :: weight_temp, cube_temp(24), cubic_mat(3,8),  prim_cell(3,max_basisnum, ng_max_node)
        !We need a pb_in variable for the interp function even though the pb doesn't matter
        integer ::pb_in(3, max_basisnum, ng_max_node)
        !arrays needed for assigning integration points to virtual atoms
        integer :: pos_to_iatom(max_size+1, max_size+1, max_size+1), atom_rep_iatom((max_size+1)**3)
        real(kind=wp) :: all_atom_set((max_size+1)), r_in_nat(3, (max_size+1)**3), rtemp(3,max_basisnum)
        character :: msg

        !etype denotes the type of the element in the element type definition array. For the rhombohedral elements
        !it's the same as the intpo_depth
        if (present(etyp)) then 
            et=etype_to_itype(etyp)
        else
            et=etype_to_itype(intpo_depth)
        end if

        !Initialize the matrix of node positions in natural coordinates
        
        cube_temp(:) = [ -1.0_wp, -1.0_wp, -1.0_wp, &
                         1.0_wp, -1.0_wp, -1.0_wp, &
                         1.0_wp, 1.0_wp, -1.0_wp, &
                         -1.0_wp, 1.0_wp, -1.0_wp, &
                         -1.0_wp, -1.0_wp, 1.0_wp, &
                         1.0_wp, -1.0_wp, 1.0_wp, &
                         1.0_wp, 1.0_wp, 1.0_wp, &
                         -1.0_wp, 1.0_wp, 1.0_wp ]
         
        prim_cell(:,1, 1:8) = reshape(cube_temp, [ 3, 8])
        pb_in(:,:,:) = 0

        do ie = 1, unique_sizes
            esize = shape_sizes(ie) 
            node_weight = 1.0_wp
            edge_weight = esize + 1 - 2*intpo_depth
            surf_rect_weight = (esize + 1 - 2*intpo_depth)**2.0_wp
            inner_weight = (esize + 1 - 2 * intpo_depth) **3.0_wp
            
            if(rank == root) then 
                weight_temp = node_weight*intpo_nums(1, et) &
                            + edge_weight*intpo_nums(2, et) &
                            + surf_rect_weight*intpo_nums(3, et) &
                            + inner_weight*intpo_nums(4, et) 
                if(.not. is_equal(weight_temp, real((esize+1)**3, wp))) then 
                    print *, "Error: Wrong total weight of integration points: ", weight_temp, ' for element type ', et, &
                             " with esize ", esize, " should be ", (esize+1)**3
                    call mpi_abort(mpi_comm_world, 1, ierr)
                end if
            end if

            !Now figure out which iatoms are at integration point sites

            !Elements are defined from -esize/2 to  esize/2 centered around 0. Here i has the number of interpolated atoms that are either 
            ! greater than or less than 0. So for esize=4, the atoms are at -2 -1 0 1 2 and i will be 2
            all_atom_set(:) = 0
            if(mod(esize,2) == 0) then 
                i = esize/2
                do ide =1, i
                    all_atom_set(2*ide) = -1.0_wp + 2.0_wp*(ide-1)/real(esize,wp)
                    all_atom_set(2*ide+1) = 1.0_wp - 2.0_wp*(ide-1)/real(esize,wp)
                end do
            else 
                i = (esize)/2 
                all_atom_set(1)= 1/real(esize,wp)
                all_atom_set(esize+1)= -1/real(esize,wp)
                do ide =1, i
                    all_atom_set(2*ide) = -1.0_wp + 2.0_wp*(ide-1)/real(esize,wp)
                    all_atom_set(2*ide+1) = 1.0_wp - 2.0_wp*(ide-1)/real(esize,wp)
                end do
            end if            

            !all_atom set contains only the positions which will have integration points. 
            !As an example, for esize =4 and intpo_depth=1 all_atom_set is (0, -2, 2). All integration points
            !positions can be created from those 3 values in natural coordinates


            !First get the position of the interpolated atoms in the natural coordinates and find the relevant position in
            !all_atom_set
            interp_count  = (esize+1)**3
            do iatom = 1, interp_count
                call interp_atom(iatom, esize, 1, pb_in, 1, prim_cell, rtemp)
                r_in_nat(:,iatom) = rtemp(:,1)
            end do

            !Now map the all_atom_set positions to iatoms
            pos_to_iatom(:,:,:) = 0
            do k = 1, esize+1
                do j = 1, esize+1
                    do i = 1, esize +1
                        x = all_atom_set(i)
                        y = all_atom_set(j)
                        z = all_atom_set(k)
                        do iatom = 1, interp_count
                            if( is_equal(r_in_nat(1,iatom), real(x,wp)).and. &
                                is_equal(r_in_nat(2,iatom), real(y,wp)).and. &
                                is_equal(r_in_nat(3,iatom), real(z,wp))) then 

                                pos_to_iatom(i,j,k) = iatom
                                exit
                            end if
                        end do
                    end do
                end do
            end do 

            intpo_set_num = 2*intpo_depth + 1

            !Now loop over all integration point
            iep=0
            i_node = 0 
            i_edge = 0
            i_surf = 0
            i_inner= 0
            atom_rep_iatom=0
            do iz = 1, esize+1
                do iy = 1, esize+1
                    do ix = 1, esize + 1

                        !Get the position for pos_to_iatom array
                        iatom = pos_to_iatom(ix,iy,iz)

                        !These are node integration points
                        if(((ix > 1).and.(ix <= intpo_set_num)).and. &
                            ((iy >1).and.(iy <= intpo_set_num)).and. &
                            ((iz > 1).and.(iz <= intpo_set_num))) then 

                            iep = iep + 1
                            atom_intpo(iep, ie, et) = iatom
                            weight_intpo(iep, ie, et) = node_weight
                            atom_rep_iatom(iatom) = iatom
                            i_node = i_node + 1

                        !These are edge integration points 
                        else if (((ix == 1).and.((iy > 1).and.(iy<=intpo_set_num)).and.((iz > 1).and.(iz<=intpo_set_num))).or. &
                                 ((iy == 1).and.((ix > 1).and.(ix<=intpo_set_num)).and.((iz > 1).and.(iz<=intpo_set_num))).or. &
                                 ((iz == 1).and.((iy > 1).and.(iy<=intpo_set_num)).and.((ix > 1).and.(ix<=intpo_set_num)))) then 

                            iep = iep+1
                            atom_intpo(iep, ie, et) = iatom
                            weight_intpo(iep, ie, et) = edge_weight
                            atom_rep_iatom(iatom) = iatom
                            i_edge = i_edge+1

                        !These are surface integration points
                        else if(((ix == 1).and.(iy == 1).and.((iz > 1).and.(iz<=intpo_set_num))).or. &
                                (((ix > 1).and.(ix<=intpo_set_num)).and.(iy == 1).and.(iz == 1)).or. &
                                ((ix == 1).and.((iy > 1).and.(iy<=intpo_set_num)).and.(iz == 1))) then

                            iep = iep+1 
                            atom_intpo(iep, ie, et) = iatom
                            weight_intpo(iep, ie, et) = surf_rect_weight
                            atom_rep_iatom(iatom) = iatom
                            i_surf = i_surf+1

                        else if((ix == 1).and.(iy == 1).and.(iz == 1)) then
                            !Interior integration point 
                            iep = iep+1 
                            atom_intpo(iep, ie, et) = iatom
                            weight_intpo(iep, ie, et) = inner_weight
                            atom_rep_iatom(iatom) = iatom
                            i_inner = i_inner + 1

                        !Now below are points which are represented by edge integration points

                        else if(((ix == 1).or.(ix > intpo_set_num)).and. &
                                 (iy > 1).and.(iy <= intpo_set_num).and. &
                                 (iz > 1).and. (iz <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(1, iy, iz)
                            i_edge = i_edge + 1

                        else if(((iy == 1).or.(iy > intpo_set_num)).and. &
                                 (ix > 1).and. (ix <= intpo_set_num).and. &
                                 (iz > 1).and.  (iz <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(ix, 1, iz)
                            i_edge = i_edge + 1


                        else if(((iz == 1).or.(iz > intpo_set_num)).and. &
                                (ix > 1).and. (ix <= intpo_set_num).and. &
                                (iy > 1).and. (iy <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(ix, iy, 1)
                            i_edge = i_edge + 1

                        !Now below are points which are represented by surface integration points

                        else if(((ix == 1).or.(ix > intpo_set_num)).and. &
                                ((iy == 1).or.(iy > intpo_set_num)).and. &
                                (iz > 1).and.(iz <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(1,1,iz)
                            i_surf = i_surf + 1

                        else if(((iy == 1).or.(iy > intpo_set_num)).and. &
                                ((iz == 1).or.(iz > intpo_set_num)).and. &
                                (ix > 1).and. (ix <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(ix, 1, 1)
                            i_surf = i_surf + 1

                        else if(((ix == 1).or.(ix > intpo_set_num)).and. &
                                ((iz == 1).or.(iz > intpo_set_num)).and. &
                                (iy > 1).and. (iy <= intpo_set_num)) then

                            atom_rep_iatom(iatom) = pos_to_iatom(1, iy, 1)
                            i_surf = i_surf + 1


                        !Otherwise it is represented by the interior integration point
                        else 
                            atom_rep_iatom(iatom) = pos_to_iatom(1,1,1)
                            i_inner = i_inner + 1

                        end if
                    end do
                end do
            end do

            !Now verify that all the counts match up
            if ( (i_node+i_edge+i_surf+i_inner) /= interp_count) then 
                print *, "Error: total considered interpolated points should equal ", interp_count, &
                         " not ", (i_node+i_edge+i_surf+i_Inner), " in init_rhomb"
                call mpi_abort(1, mpi_comm_world, ierr)
            end if

            !We now form the mapping array which dictates where in atom_intpo our representative integration point lies
            who_rep_atomap = 0
            do iatom = 1, interp_count
                do iep = 1, intpo_count(et)
                    if (atom_rep_iatom(iatom) == atom_intpo(iep, ie, et)) then 
                        who_rep_atomap(iatom,ie, et) = iep 
                    end if
                end do
                if (atom_rep_iatom(iatom)==0) then 
                    write(msg, *) "atom_rep_iatom should not be 0 for iatom", iatom, " and esize ", esize
                    call misc_error(msg)
                end if
            end do
        end do

    end subroutine init_rhomb

    subroutine update_intpo
        !This subroutine updates the ownership arrays for the integration points
        
        integer :: ie, iep, j, iatom, iatomap, ibasis, intpo_num_sum, owner_action, owner_slot
        logical :: in_owner_block
        
        !First allocate who_has_intpo
        if(allocated(who_has_intpo)) then 
            if (ele_num_l > size(who_has_intpo,2)) deallocate(who_has_intpo)
        end if
        if(.not.allocated(who_has_intpo)) allocate(who_has_intpo(max_intpo_num*max_basisnum, ele_num_l))


        intpo_num_l = 0
        if(atomap_num_l > size(atomap_to_intpo,2)) then 
            deallocate(atomap_to_intpo)
            allocate(atomap_to_intpo(2,size(r_atomap,2)))
        end if
        atomap_to_intpo = 0
        who_has_intpo(:, :) = .false.
        do ie = 1, ele_num_l
            j = size_to_shape(size_ele(ie)) 

            do iep = 1, intpo_count(itype(ie))
                iatom = atom_intpo(iep, j, itype(ie))
                
                do ibasis = 1, basis_num(ie)
                    iatomap=cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                    owner_slot = intpo_owner_slot_index(basis_num(ie), iep, ibasis)
                    in_owner_block = .false.
                    if (iatomap /= 0) in_owner_block = in_block_bd(r_atomap(:,iatomap), pro_bd)
                    owner_action = classify_intpo_owner(iatomap, in_owner_block)

                    select case(owner_action)
                    case (intpo_owner_skip)
                        cycle
                    case (intpo_owner_register)
                        intpo_num_l = intpo_num_l + 1 
                        who_has_intpo(owner_slot, ie) = .true.
                        atomap_to_intpo(1,iatomap) = owner_slot
                        atomap_to_intpo(2,iatomap) = ie
                    case (intpo_owner_error)
                        print *, "Error: Atomap ", iatomap, " should belong to ", rank, " but doesn't"
                        print *, "Pos is ", r_atomap(:,iatomap), " and bd is ", pro_bd
                        call mpi_abort( mpi_comm_world, 1, ierr)
                    end select

                end do

            end do
        end do 

        !Now check to make sure we have the right number of intpo
        call mpi_reduce(intpo_num_l, intpo_num_sum, 1, mpi_integer, mpi_sum, root, mpi_comm_world, ierr)

        if(rank == root) then 
            if(intpo_num_sum /= intpo_num) then 
                print *, "Error: Total integration point number ", intpo_num_sum, " does not equal ", intpo_num
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
        end if

        
    end subroutine update_intpo

    pure integer function intpo_owner_slot_index(basis_count, iep, ibasis)
        integer, intent(in) :: basis_count, iep, ibasis

        intpo_owner_slot_index = basis_count*(iep-1) + ibasis
    end function intpo_owner_slot_index

    pure integer function classify_intpo_owner(iatomap, in_owner_block)
        integer, intent(in) :: iatomap
        logical, intent(in) :: in_owner_block

        if (iatomap == 0) then
            classify_intpo_owner = intpo_owner_skip
        else if (in_owner_block) then
            classify_intpo_owner = intpo_owner_register
        else
            classify_intpo_owner = intpo_owner_error
        end if
    end function classify_intpo_owner

    pure function mass_mat_coeff(esize, etype)
        integer, intent(in) :: esize, etype 
        real(kind=wp) :: mass_mat_coeff

        mass_mat_coeff = real(esize+1,wp)**3.0_wp/real(ng_node(etype),wp)
        return
    end function mass_mat_coeff

    subroutine update_itype
        !This just updates the itype array for the local element num
        integer :: i

        if (allocated(itype)) deallocate(itype)
        allocate(itype(ele_num_l))

        do i = 1, ele_num_l
            itype(i) = etype_to_itype(etype(i))
        end do
    end subroutine

end module integration
