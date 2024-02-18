module quenched_dynamics

    use parameters
    use comms
    use elements
    use neighbors
    use potential 
    use time
    use vel_verlet

    implicit none
    public 
    contains

    subroutine qd(i)
        integer, intent(in) :: i 

        call update_r
        call update_neighbor(i)
        call update_force
        call quenched_vel
    end subroutine qd

    subroutine quenched_vel
        integer :: ia, ie, inod, ip, ibasis
        real(kind = wp) :: pro_vel_p, pro_force_normsq, force_normsq, force_norm, vel_p


        pro_vel_p = 0.0_wp
        pro_force_normsq = 0.0_wp

        !coarse-grained domain
        if(ele_num /= 0) then
            do ie = 1, ele_num_l
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        vel(:, ibasis, ip) = vel(:, ibasis, ip) &
                                       + 0.5_wp * (force_eq_pre(:, ibasis, ip) &
                                       + force_eq(:, ibasis, ip)) / masses(basis_type(ibasis,ie)) &
                                       * time_step / const_motion

                        if(who_has_ele(ie).eqv..true.) then
                            pro_vel_p = pro_vel_p + dot_product(vel(:, ibasis, ip), force_eq(:, ibasis, ip))
                            pro_force_normsq = pro_force_normsq + dot_product(force_eq(:, ibasis, ip), force_eq(:, ibasis, ip))
                        end if

                    end do

                end do

            end do

        end if

        !atomistic domain
        if(atom_num /= 0) then
            do ia = 1, atom_num_l
                vel_atom(:, ia) = vel_atom(:, ia) &
                                      + 0.5_wp * (force_atom_pre(:, ia) &
                                      + force_atom(:, ia)) / masses(type_atom(ia)) &
                                      * time_step / const_motion

                pro_vel_p = pro_vel_p + dot_product(vel_atom(:, ia), force_atom(:, ia))
                pro_force_normsq = pro_force_normsq + dot_product(force_atom(:, ia), force_atom(:, ia))

            end do

        end if

        !sum pro_vel_p and pro_force_normsq up to get vel_p and force_normsq
        if(pro_num == 1) then
            vel_p = pro_vel_p
            force_normsq = pro_force_normsq

        else
            call mpi_allreduce(pro_vel_p, vel_p, 1, mpi_wp, mpi_sum, world, ierr)
            call mpi_allreduce(pro_force_normsq, force_normsq, 1, mpi_wp, mpi_sum, world, ierr)

        end if

        force_norm = sqrt(force_normsq)
!       decide if the nodes/atoms should be freezed
!       coarse-grained domain

        if(ele_num /= 0) then
            do ip = 1, node_num_l
                if((vel_p < 0.0_wp).or. (abs(vel_p) < lim_zero)) then
                    vel(:, :, ip) = 0.0_wp
                else
                    !only use the component of the velocity parallel to the force
                    vel(:, :, ip) = vel_p * force_eq(:, :, ip) / force_normsq
                end if
            end do
        end if

        !atomistic domain
        if(atom_num /= 0) then
            do ia = 1, atom_num_l
                if((vel_p < 0.0_wp).or. (abs(vel_p) < lim_zero)) then
                    vel_atom(:, ia) = 0.0_wp
                else
                    !only use the component of the velocity parallel to the force
                    vel_atom(:, ia) = vel_p * force_atom(:, ia) / force_normsq
                end if
            end do
        end if
    end subroutine quenched_vel

end module quenched_dynamics
