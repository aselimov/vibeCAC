module dynamics

    !This code contains the code needed to run the various types of dynamics
    use parameters
    use vel_verlet
    use neighbors
    use quenched_dynamics
    use forces
    use comms
    use potential
    use dump
    use thermo
    use time
    use temp
    use berendsen
    use deform

    implicit none

    real(kind=wp), save :: itime
    character(len = 100) :: ensemble
    integer :: ensemble_option

    private :: pre_step, post_step, step
    public
    contains

    subroutine dynamics_defaults
        ensemble_option = 0 
        ensemble = "none"
    end subroutine dynamics_defaults

    subroutine parse_dynamics(line)
        character(len = *), intent(in) :: line
        
        character(len=read_len) tmptxt
        integer :: iospara

        read(line, *, iostat = iospara) tmptxt, ensemble
        if (iospara > 0) call read_error("Invalid read of dynamics command", iospara)

        select case(ensemble)
        case("none")
            ensemble_option = 0
            if(rank==root) call log_msg("Warning: no dynamics specified. Atom positions will not evolve")
        case('NVE', 'nve')
            ensemble_option = 1
        case('qd')
            ensemble_option = 2
        case('eu')
            ensemble_option = 3
        case('ld')
            ensemble_option = 4 
        case default 
            call command_error("Ensemble "//trim(adjustl(ensemble))//" is not currently accepted as an option for dynamics")
        end select

        return
    end subroutine parse_dynamics

    subroutine parse_run(line)
        !parse the run command
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt
        integer ::  iospara

        begin_step = iter
        read(line, *, iostat = iospara) tmptxt, run_steps
        if (iospara > 0) call read_error("Invalid read of run command", iospara)
        call run_dynamics(run_steps)

    end subroutine parse_run

    subroutine run_dynamics(num_steps)
        !This subroutine actually runs the dynamics steps
        integer, intent(in) :: num_steps

        integer ::i

        !Check to make sure timestep is greater than 0
        if(time_step <= 0.0_wp) then 
            call misc_error("Time_step was never set, please add command time_step val before run command")
        end if
        !Initialize the force calculation and the first dump

        if(nei_init) then 
            call update_neighbor(iter, .false., .true.)
        else
            if(ele_num > 0) call ghost_cg
            if(atom_num> 0) call ghost_at
            call neighbor_lists
        end if
        call update_force
        if(first_run) call write_dump(iter, .true.)
        call write_thermo_style
        call write_thermo_out(iter)

        do i = 1, num_steps
            iter = iter + 1 
            t = t+time_step

            call pre_step(iter)
            call step
            call post_step(iter)
        end do

        need_virial = .true.
        call update_force
        call write_dump(iter, .true.)
        call write_thermo_out(iter)

        call log_neighbor_info
        first_run = .false.

    end subroutine run_dynamics

    subroutine pre_step(i)
        !This subroutine is run before each dynamics timestep is calculated
        
        !i is the current step number
        integer, intent(in) ::i
        
        !If we are dumping this timestep then we need to calculate the virial stress
        if(need_dump(i).or.pflag) need_virial = .true.

        !Check to see if we need virial for the thermo command
        if((mod(i,thermo_every)==0).and.need_p) need_virial = .true. 

        call deform_box

        return
    end subroutine pre_step

    subroutine post_step(i)
        !This subroutine is run after the dynamics timestep is calculated
        !i is the current step number
        integer, intent(in) ::i

        !Check to see if we need to rescale velocity
        if(tflag) call rescale_v
        
        !Check to see if we need to dilate box
        if(pflag) call rescale_box

        !Dump check
        call write_dump(i)

        !Thermo if we need it 
        if(mod(i,thermo_every)==0) call write_thermo_out(i)


        need_virial = .false.

        return
    end subroutine post_step

    subroutine step
        !This subroutine calls the correct dynamics to run
        !i is the current step number
        select case(ensemble_option)
        case(0)
            continue
        case(1)
            call verlet(iter)
        case(2)
            call qd(iter)
        case(3)
            call euler(iter)        
        end select
        
    end subroutine step

    subroutine euler(iter)
        integer, intent(in) :: iter

        integer :: ip, ibasis, ia

        do ip = 1, node_num_l
            do ibasis=1, basis_num(node_cg(ip))
                vel(:, ibasis, ip) = vel(:, ibasis, ip) + time_step*ftm2v/masses(basis_type(ibasis,node_cg(ip))) &
                                                          *force_eq(:,ibasis,ip)
                r(:, ibasis, ip) = r(:, ibasis, ip) + time_step*vel(:, ibasis,ip)
            end do
        end do

        do ia = 1, atom_num_l
            vel_atom(:,ia) = vel_atom(:,ia) + time_step*ftm2v/masses(type_atom(ia))*force_atom(:,ia)
            r_atom(:,ia) = r_atom(:,ia) + time_step*vel_atom(:,ia)
        end do

        call update_neighbor(iter)
        call update_force
    end subroutine euler

end module dynamics
