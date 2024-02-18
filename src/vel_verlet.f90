module vel_verlet
    !This module contains all the code for running steps of the verlet algorithm
    !This is equivalent to the NVE ensemble
    use parameters
    use comms 
    use elements
    use neighbors 
    use potential
    use time

    implicit none

    public
    contains

    subroutine verlet(i)
        integer, intent(in) :: i
        !This subroutine advances the model by 1 timestep using the velocity verlet algorithm
        !Update position at time i+1
        call update_r
        !Update neighbor at time i+1
        call update_neighbor(i)
        !Update neighbor at time i+1
        call update_force
        !Update velocity at time i+1
        call update_vel(time_step)
    end subroutine verlet

    subroutine update_r
        !This subroutine updates the positions for all atoms and coarse-grained element nodes
        
        integer :: ia, ip, ibasis

        if(ele_num /= 0) then 
            do ip = 1, node_num_l
                do ibasis = 1, basis_num(node_cg(ip))
                    r(:,ibasis,ip) = r(:, ibasis, ip) + vel(:, ibasis, ip) * time_step &
                                    + 0.5_wp * force_eq(:, ibasis, ip)/ masses(basis_type(ibasis,node_cg(ip))) &
                                    * time_step**2.0_wp / const_motion
                end do
            end do
            force_eq_pre(:,:,1:node_num_l) = force_eq(:,:,1:node_num_l)
        end if

        if(atom_num /= 0) then 
            do ia = 1, atom_num_l
                r_atom(:,ia) = r_atom(:,ia) + vel_atom(:, ia)*time_step + 0.5_wp*force_atom(:,ia)/masses(type_atom(ia)) &
                              * time_step**2.0_wp / const_motion  
            end do
            force_atom_pre(:,1:atom_num_l) = force_atom(:, 1:atom_num_l)

        end if

        return
    end subroutine update_r

    subroutine update_vel(time_step)
        real(kind=wp), intent(in) :: time_step
        !update the velocity
        integer :: ie, inod, ip, ia, ibasis
        real(kind = wp) :: pro_force_normsq, avg(3) 

        pro_force_normsq = 0.0_wp
        if(ele_num /= 0) then
            do ie = 1, ele_num_l
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        vel(:, ibasis, ip) = vel(:, ibasis, ip) + 0.5_wp * (force_eq_pre(:, ibasis, ip) &
                                             + force_eq(:, ibasis, ip)) / masses(basis_type(ibasis,ie)) &
                                             * time_step / const_motion
                    end do
                end do
            end do
        end if

        !atomistic domain
        if(atom_num /= 0) then
            do ia = 1, atom_num_l
                vel_atom(:, ia) = vel_atom(:, ia) + 0.5_wp * (force_atom_pre(:, ia) &
                                + force_atom(:, ia)) / masses(type_atom(ia)) &
                                * time_step / const_motion
            end do

        end if

        !Now take out the average velocity
        avg = compute_avgvel(1)
        if(ele_num > 0) then 
            do ie = 1, ele_num_l
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        vel(:, ibasis, ip) = vel(:, ibasis, ip) - avg
                    end do
                end do
            end do
        end if

        !atomistic domain
        if(atom_num /= 0) then
            do ia = 1, atom_num_l
                vel_atom(:, ia) = vel_atom(:, ia) - avg
            end do

        end if

    end subroutine update_vel
end module vel_verlet
