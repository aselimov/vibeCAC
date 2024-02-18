module minimize

    !This code contains the code needed to run the various types of dynamics
    use parameters
    use fire
    use cg
    use neighbors
    use forces
    use comms
    use potential
    use dump
    use thermo
    use errors
    use logger
    use debug

    implicit none

    real(kind=wp), save :: force_tol, energy_tol
    integer :: min_style, max_iter, reset_num
    real(kind=wp), parameter :: eps_energy = 1d-8

    private :: pre_iter, post_iter, iterate
    public
    contains

    subroutine min_defaults
        !Set defaults
        min_style = 1
    end subroutine min_defaults

    subroutine parse_min_style(line)
        !parse the run command
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, min_string, args(20)
        integer ::iospara, j, i
        real(kind=wp) :: val

        read(line, *, iostat = iospara) tmptxt, min_string
        if (iospara > 0) call read_error("Invalid read of min_string from min_style command", iospara)

        select case(min_string)
        case('fire', 'FIRE')
            call fire_defaults
            min_style = 1
            j = tok_count(line) 
            if(j > 2) then 
                read(line, *) args(1:j)
                i = 3
                do while(i <= j)
                    select case(args(i))
                    case('alpha0', 'tmax', 'tmin')
                        i = i + 1
                        !Read param val
                        read(args(i), *) val
                        !set param val
                        call set_fire_param(args(i-1),val)
                    case default
                        write(tmptxt, *) "Cannot pass argument ", trim(adjustl(args(i))), " to min_style with style cg"
                        call misc_error(tmptxt)
                    end select
                    i = i + 1
                end do
            end if
        case('cg', 'CG')
            call command_error('Conjugate gradient not correctly working in this implementation of CAC')
            call cg_defaults
            min_style = 2
            !Now parse additional options for cg
            j = tok_count(line) 
            if(j > 2) then 
                read(line, *) args(1:j)
                i = 3
                do while(i <= j)
                    select case(args(i))
                    case('reset')
                        i = i + 1
                        !Read the number of iterations between resets
                        read(args(i), *) reset_num
                    case default
                        write(tmptxt, *) "Cannot pass argument ", trim(adjustl(args(i))), " to min_style with style cg"
                        call misc_error(tmptxt)
                    end select
                    i = i + 1
                end do
            end if
        case default 
            call command_error("Min_style "//trim(adjustl(min_string))//" is not currently accepted as an option for dynamics")
        end select

    end subroutine parse_min_style

    subroutine parse_minimize(line)
        character(len=*), intent(in) :: line
        character(len=read_len) :: tmptxt
        integer :: iospara

        read(line, *, iostat = iospara) tmptxt,  energy_tol, force_tol, max_iter
        if (iospara > 0) call read_error("Invalid read of minimize command", iospara)

        if ((force_tol < 0)) then 
            write(tmptxt, *) "Force tolerance ", force_tol, " should not be less than zero"
            call misc_error(tmptxt)
        end if
        if ((energy_tol< 0)) then 
            write(tmptxt, *) "Energy tolerance ", energy_tol, " should not be less than zero"
            call misc_error(tmptxt)
        end if
        if ((max_iter<= 0)) then 
            write(tmptxt, *) "max_iter", max_iter, " should not be less than or equal to zero"
            call misc_error(tmptxt)
        end if
        
        begin_step = iter
        call run_min
    end subroutine parse_minimize

    subroutine run_min
        !This subroutine actually runs the minimization
        integer ::i, exit_cond, delay, nevals
        real(kind=wp) :: old_energy, old_f_norm, pe, fnorm
        character(len=read_len) :: msg


        !Init neighbor list and do initial calc
        if (nei_init) then 
            call update_neighbor(iter, .false.,.true.)
        else
            if(ele_num > 0) call ghost_cg
            if(atom_num> 0) call ghost_at
            call neighbor_lists
        end if
        call update_force

        !Initialize minimizer 
        delay = 0
        select case(min_style) 
        case(1)
            !Check to make sure that the timestep has been set 
            if (is_equal(time_step, 0.0_wp)) then 
                call misc_error("Time step must be set prior to calling minimize with min_style fire")
            end if
            call fire_init
            delay = get_fire_delaystep()
            !Initialize velocities to 0
            if(node_num_l > 0) vel(:,:,:) = 0
            if(atom_num_l > 0) vel_atom(:,:) = 0
        case(2)
            call cg_init
        end select

        if(first_run) call write_dump(iter, .true.)
        call write_thermo_style
        call write_thermo_out(iter)

        exit_cond = 0

        pe = compute_pe(1, .true.)
        fnorm = compute_fnorm(1)
        do i = 1, max_iter

            iter = iter + 1
            !Save old energies
            old_energy=pe
            old_f_norm =fnorm 

            !Call pre_iterate
            call pre_iter(i)

            !Call iterate
            call iterate(i, exit_cond)

            if(min_style == 2) then 
                call update_neighbor(iter, .true.)
                call update_force

                !Now check to see if we want to reset the conjugate direction. This may slow down the convergance, but may be
                !useful for finite elements
                if(reset_num > 0) then 
                    if(mod(i, reset_num) == 0) then 
                        call reset_cg_dir
                    end if
                end if
            end if

            !Exit if the minimizer tells us to
            if(exit_cond > 0) exit 

            !Check tolerances
            if (i > delay) then 
                pe = compute_pe(1, .true.)
                fnorm = compute_fnorm(1)
                if(abs(old_energy-pe)< energy_tol*0.5_wp*(abs(pe)+abs(old_energy)+eps_energy)) then 
                    exit_cond = 1
                    exit
                else if (fnorm < force_tol) then 
                    exit_cond = 2
                    exit
                end if
            end if

            !Call post iterate
            call post_iter(i)
        end do

        !Compute virial
        need_virial = .true.
        call update_force

        !Post completion thermo and dump
        call write_thermo_out(iter)
        call write_dump(iter, .true.)

        !Post minimizer message to user
        !First print the exit conditions
        if(exit_cond < 3) then 
            select case(exit_cond)
            case(0) 
                write(msg, *) "Max Iterations"
            case(1)
                write(msg, *) "Energy Tolerance"
            case(2)
                write(msg, *) "Force Tolerance"
            end select
        else 
            select case(min_style)
            case(1)
                if(exit_cond == 3) write(msg, *) "Max P < 0 iterations"
            case(2)
                if(exit_cond == 4) write(msg, *) "Exit because search direction isn't downhill"

                if(exit_cond == 5) write(msg, *) "All search direction components are 0"

                if(exit_cond == 6) write(msg, *) "Alpha is equal to 0"

                if(exit_cond == 11) write(msg, *) "No change to force between line search iterations"
            end select
        end if
        call log_msg("Minimizer out with stopping condition: "//trim(adjustl(msg)))

        !Now print number of evaluations and iterations and clean up 
        select case(min_style)
        case(1) 
            nevals = get_fire_neval()
            call fire_clean
        case(2)
            nevals = get_cg_neval()
            call cg_clean
        end select
        write(msg,*) "Total number of iterations: ", i
        call log_msg(msg)
        write(msg, *) "Total number of force evals: ", nevals
        call log_msg(msg)

        write(msg, *) "Final and second-to-last energies: ", pe, old_energy
        call log_msg(msg)

        call log_neighbor_info

        first_run = .false.
    end subroutine run_min

    subroutine pre_iter(i)
        !This subroutine is run before each dynamics timeiter is calculated
        
        !i is the current iter number
        integer, intent(in) ::i
        
        !If we are dumping this timeiter then we need to calculate the virial stress
        if(need_dump(i)) need_virial = .true.
        if((mod(i,thermo_every)==0).and.need_p) need_virial = .true. 
        return
    end subroutine pre_iter

    subroutine post_iter(i)
        integer, intent(in) :: i

        integer :: pre_atom_num
        !This subroutine is run after the dynamics timeiter is calculated
        !i is the current iter number
        
        !Dump if we need it
        call write_dump(i)

        !Thermo if we need it 
        if(mod(i,thermo_every)==0) call write_thermo_out(iter)

        !Debug if we need it
        if(dflag) call run_debug

        need_virial = .false.
        return
    end subroutine post_iter

    subroutine iterate(i, mincode)
        !This subroutine calls the correct dynamics to run
        !i is the current iter number
        integer, intent(in) :: i
        integer, intent(out) :: mincode
        select case(min_style)
        case(1)
            call fire_iterate(i, mincode)
        case(2) 
            call cg_iterate(mincode)
        end select
        
    end subroutine iterate

end module minimize
