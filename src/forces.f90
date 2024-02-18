module forces
    !This module is in charge of calculating energy, force, and virial

    use parameters
    use elements
    use integration

    implicit none

    real(kind = wp) :: energy_sum, energy_atom_sum, energy_tally, f_norm_tally, energy_equiv_tally
    !Atom arrays
    real(kind = wp), allocatable :: energy_atom(:), force_atom(:,:), force_atom_pre(:,:), virial_atom(:,:,:)

    !Element arrays
    real(kind=wp), allocatable :: energy(:,:), energy_eq(:,:), &
                                  force(:,:,:), force_eq(:,:,:), force_eq_pre(:,:,:), &
                                  virial(:,:,:,:), virial_eq(:,:,:,:)
    logical :: need_force_pre

    public 
    contains

    subroutine alloc_force_arrays
        if(atom_num > 0) then 
            if( allocated(force_atom)) then 
                deallocate(energy_atom, force_atom, virial_atom, stat=allostat)
                if (allostat > 0) call alloc_error("Failure to deallocate force_atom arrays", allostat)
            end if
            allocate(energy_atom(atom_num_l), force_atom(3, atom_num_l), &
                     virial_atom(3,3,atom_num_l), stat=allostat)                
            if (allostat > 0) call alloc_error("Failure to allocate force_atom arrays", allostat)
        end if
        if(ele_num> 0) then 
            if(allocated(force_eq)) then 
                deallocate(force_eq, energy_eq, virial_eq, stat=allostat) 
                if (allostat > 0) call alloc_error("Failure deallocating eq arrays", allostat)
            end if
            allocate(force_eq(3, max_basisnum, node_num_l), energy_eq(max_basisnum, node_num_l), &
                     virial_eq(3, 3, max_basisnum, node_num_l), stat = allostat) 
            if (allostat > 0) call alloc_error("Failure allocating init_force_arrays", allostat)
        end if
    end subroutine alloc_force_arrays

    subroutine alloc_pre_array
        if(atom_num > 0) then 
            allocate(force_atom_pre(3, atom_num_lr))
        end if
        if(ele_num>0) then 
            allocate(force_eq_pre(3, max_basisnum, node_num_lr))
        end if
    end subroutine alloc_pre_array

    subroutine resize_force_arrays
        !Resize the force arrays to match the node number
        real(kind=wp), allocatable :: force_pre_array(:,:,:), force_atom_array(:,:)
        
        call alloc_force_arrays
        !Force_eq_pre is the only array that we have to preserve data for because the other ones update prior to use
        if(need_force_pre) then 
            if(ele_num > 0) then 
                allocate(force_pre_array(3,max_basisnum, node_num_l), stat = allostat)
                if(allostat > 0) call alloc_error("Failure to allocate force_pre_array", allostat)
                force_pre_array = force_eq_pre(:, :, 1:node_num_l)
                call move_alloc(force_pre_array, force_eq_pre)
            end if

            if(atom_num > 0) then 
                allocate(force_atom_array(3,atom_num_l), stat = allostat)
                if(allostat > 0) call alloc_error("Failure to allocate force_atom_array", allostat)
                force_atom_array = force_atom_pre(:,1:atom_num_l) 
                call move_alloc(force_atom_array, force_atom_pre)
            end if
        end if

    end subroutine resize_force_arrays

    subroutine update_equiv
        !This subroutine calculates the normalized nodal quantities (referred to here as the equivalent)
        integer :: ie, ip, i, j, je, inod, jnod, ibasis, ind
        real(kind = wp), allocatable :: force_array(:), force_buff(:), virial_array(:), virial_buff(:), &
                                        energy_array(:), energy_buff(:)




        !If we have more than one processor than we need to share all this information for all elements that share data
        if(pro_num > 1) then
            allocate(force_array(3*ng_max_node*ele_shared_num*max_basisnum), &
                     energy_array(ng_max_node*ele_shared_num*max_basisnum), stat = allostat)
            if(allostat /= 0) call alloc_error("Failure to allocate force/energy_array in update_equiv", allostat)

            force_array(:) = 0.0_wp
            energy_array(:) = 0.0_wp

            if(need_virial) then
                allocate(virial_array(9*ng_max_node*ele_shared_num*max_basisnum), stat = allostat)
                if(allostat /= 0) call alloc_error("Failure to allocate virial_array in update_equiv", allostat)
                virial_array(:) = 0.0_wp
            end if

            !Assign the buffers for reducing to calculate the equivalent forces. 
            !The buffer will have empty space so the communications aren't completely 
            !efficient memory wise but it shouldn't cause very bad slowdowns.
            do ie = 1, ele_num_l
                je = ele_id_shared(ie)

                if(je /= 0) then
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis = 1, basis_num(ie)

                            !Calculate the current index  and put the force, virial, and energy into the send arrays
                            ind = 3*ng_max_node*max_basisnum*(je-1) + 3*max_basisnum*(inod-1) + 3*(ibasis-1)
                            do i = 1, 3
                                if(ind + i > 3*ng_max_node*ele_shared_num*max_basisnum) then
                                    print *, 'Error: Index of force_array', &
                                             ind+i, ' is larger than', &
                                             ' array size', 3*ng_max_node*ele_shared_num*max_basisnum
                                    call mpi_abort(mpi_comm_world, 1, ierr)
                                end if

                                force_array(ind+i) = force(i, ibasis, ip)

                                if(need_virial) then
                                    do j = 1, 3
                                        if(3*ind + 3*(i-1)+j > 9*ng_max_node*ele_shared_num*max_basisnum) then
                                            print *, 'Error: Index of virial_array', &
                                                      3*ind + 3*(i-1)+j, ' is larger than', &
                                                     ' array size', 9*ng_max_node*ele_shared_num
                                            call mpi_abort(mpi_comm_world, 1, ierr)
                                        end if
                                        virial_array(3*ind+3*(i-1)+j) = virial(j, i, ibasis, ip)
                                    end do
                                end if
                            end do

                            if(ng_max_node*max_basisnum*(je-1)+max_basisnum*(inod-1)+ibasis &
                                > ng_max_node*ele_shared_num*max_basisnum) then
                                print *, 'Error: Index of energy_array', &
                                          ng_max_node*(je-1)+max_basisnum*(inod-1)+ibasis, ' is larger than', &
                                         ' array size', ng_max_node*ele_shared_num*max_basisnum
                                call mpi_abort(mpi_comm_world, 1, ierr)

                            end if
                            energy_array(ng_max_node*max_basisnum*(je-1)+max_basisnum*(inod-1)+ibasis) = energy(ibasis, ip)
                        end do
                    end do
                end if
            end do

            !allocate force_buff, virial_buff, and energy_buff
            allocate(force_buff(3*ng_max_node*ele_shared_num*max_basisnum), &
                     energy_buff(ng_max_node*ele_shared_num*max_basisnum), stat = allostat)
            if(allostat /= 0) call alloc_error("Failure to allocate force/energy_buff in update_equiv", allostat)

            !Reduce all force_buff, energy_buff, and virial_buff
            force_buff(:) = 0.0_wp
            call mpi_allreduce(force_array, force_buff, 3*ng_max_node*ele_shared_num*max_basisnum, &
                               mpi_wp, mpi_sum, world, ierr)

            energy_buff(:) = 0.0_wp
            call mpi_allreduce(energy_array, energy_buff, ng_max_node*ele_shared_num*max_basisnum, &
                               mpi_wp, mpi_sum, world, ierr)

            if(need_virial) then
                allocate(virial_buff(9*ng_max_node*ele_shared_num*max_basisnum), stat = allostat) 
                if(allostat /= 0) call alloc_error("Failure to allocate virial_buff in update_equiv", allostat)

                virial_buff(:) = 0.0_wp
                call mpi_allreduce(virial_array, virial_buff, 9*ng_max_node*ele_shared_num*max_basisnum, &
                                   mpi_wp, mpi_sum, world, ierr)
            end if

            !update local energy_eq, force_eq, virial_eq from the buff received from all processors
            do ie = 1, ele_num_l
                je = ele_id_shared(ie)
                if(je /= 0) then
                  do inod = 1, ng_node(etype(ie))

                    ip = cg_node(inod, ie)

                    do ibasis = 1, basis_num(ie)

                      ind = 3*ng_max_node * max_basisnum * (je-1) + 3*(inod-1)*max_basisnum + 3*(ibasis-1)

                      do i = 1, 3
                        force(i, ibasis, ip) = force_buff(ind+i)

                        if(need_virial) then
                            do j = 1, 3
                                virial(j, i, ibasis, ip) = virial_buff(3*ind + 3*(i-1)+j)
                            end do
                        end if
                      end do
                      energy(ibasis, ip) = energy_buff(ng_max_node*max_basisnum*(je-1)+max_basisnum*(inod-1) + ibasis)
                    end do
                  end do
                end if
            end do
        end if

        !Reset values for eq arrays
        force_eq(:, :, :) = 0.0_wp
        if(need_virial) then
            virial_eq(:, :, :, :) = 0.0_wp
        end if
        energy_eq(:, :) = 0.0_wp

        do ie = 1, ele_num_l

            do inod = 1, ng_node(etype(ie))

                ip = cg_node(inod, ie)

                select case(mass_mat_ty) 
                case('lumped')
                    !Calculate equivalent values by dividing by the lumped mass matrix
                    do ibasis = 1, basis_num(ie)
                        force_eq(:,ibasis, ip) = force(:,ibasis, ip)/mass_mat_coeff(size_ele(ie), etype(ie))
                        if(need_virial) then
                            virial_eq(:, :, ibasis, ip) = virial(:, :, ibasis, ip) / mass_mat_coeff(size_ele(ie), etype(ie))
                        end if
                        energy_eq(ibasis, ip) = energy(ibasis, ip)/mass_mat_coeff(size_ele(ie), etype(ie))
                    end do

                case('consistent')
                    !Use a consistent mass matrix
                    do jnod = 1, ng_node(etype(ie))
                        do ibasis = 1, basis_num(ie)
                            force_eq(:, ibasis, ip) = force(:, ibasis, ip) * mass_mat_inv(inod, jnod) &
                                                      / mass_mat_coeff(size_ele(ie),etype(ie))
                          if(need_virial) then
                            virial_eq(:, :, ibasis, ip) = virial(:, :, ibasis, ip) * mass_mat_inv(inod, jnod) &
                                                           / mass_mat_coeff(size_ele(ie),etype(ie))
                          end if
                          energy_eq(ibasis, ip) = energy(ibasis, ip) * mass_mat_inv(inod, jnod) &
                                                  / mass_mat_coeff(size_ele(ie), etype(ie))
                        end do
                    end do
 
                case default
                    print *, 'Error: Mass matrix type ', mass_mat_ty, ' is not accepted'
                    call mpi_abort(mpi_comm_world, 1, ierr)
                end select
            end do
        end do
        !Deallocate arrays we no longer need
        deallocate(force, energy, stat = deallostat) 
        if(need_virial) then
            deallocate(virial, stat = deallostat) 
        end if
        if(deallostat /= 0) call alloc_error("Failure to deallocate force/energy/virial", deallostat)

        return
    end subroutine update_equiv

    subroutine communicate_force_eq
        !This subroutine communicates forces for all shared elements when needed (such as when using langevin dynamics 
        integer :: ie, ip, i, je, inod, ibasis, ind
        real(kind = wp), allocatable :: force_array(:), force_buff(:)

        !If we have more than one processor than we need to share all this information for all elements that share data
        if(pro_num > 1) then
            allocate(force_array(3*ng_max_node*ele_shared_num*max_basisnum), stat = allostat)
            if(allostat /= 0) call alloc_error("Failure to allocate force in communicate_force_eq", allostat)

            force_array(:) = 0.0_wp

            !Assign the buffers for reducing to calculate the equivalent forces. 
            !The buffer will have empty space so the communications aren't completely 
            !efficient memory wise but it shouldn't cause very bad slowdowns.
            do ie = 1, ele_num_l
                je = ele_id_shared(ie)

                if(je /= 0) then
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis = 1, basis_num(ie)
                            !Calculate the current index  and put the force, virial, and energy into the send arrays
                            ind = 3*ng_max_node*max_basisnum*(je-1) + 3*max_basisnum*(inod-1) + 3*(ibasis-1)
                            do i = 1, 3
                                if(ind + i > 3*ng_max_node*ele_shared_num*max_basisnum) then
                                    print *, 'Error: Index of force_array', &
                                             ind+i, ' is larger than', &
                                             ' array size', 3*ng_max_node*ele_shared_num*max_basisnum
                                    call mpi_abort(mpi_comm_world, 1, ierr)
                                end if
                                force_array(ind+i) = force_eq(i, ibasis, ip)
                            end do
                        end do
                    end do
                end if
            end do

            !allocate force_buff, virial_buff, and energy_buff
            allocate(force_buff(3*ng_max_node*ele_shared_num*max_basisnum), stat = allostat)
            if(allostat /= 0) call alloc_error("Failure to allocate force_buff in update_equiv", allostat)

            !Reduce all force_buff, energy_buff, and virial_buff
            force_buff(:) = 0.0_wp
            call mpi_allreduce(force_array, force_buff, 3*ng_max_node*ele_shared_num*max_basisnum, &
                               mpi_wp, mpi_sum, world, ierr)

            !update local energy_eq, force_eq, virial_eq from the buff received from all processors
            do ie = 1, ele_num_l
                je = ele_id_shared(ie)
                if(je /= 0) then
                  do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                      ind = 3*ng_max_node * max_basisnum * (je-1) + 3*(inod-1)*max_basisnum + 3*(ibasis-1)
                      do i = 1, 3
                        force_eq(i, ibasis, ip) = force_buff(ind+i)
                      end do
                    end do
                  end do
                end if
            end do
        end if
        return
    end subroutine communicate_force_eq

    subroutine reset_tallies
        !Reset energy and force tallies
        energy_tally = 0
        energy_equiv_tally = 0
        f_norm_tally = 0
    end subroutine

end module forces
