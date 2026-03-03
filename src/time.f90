module time
    use parameters
    use logger
    use errors
    use mpi

    implicit none
    !This subroutine contains variables and subroutines associated with timestepping
    !also contains code for timing
    
    integer :: itime_start, iter, run_steps, begin_step
    real(kind=wp) :: t, time_step, orig_time_step
    real(kind=wp) :: walltime(2), program_start, program_end

    abstract interface
        subroutine timestep_error_handler(message)
            character(len=*), intent(in) :: message
        end subroutine timestep_error_handler
    end interface
    
    public 
    contains

    pure logical function is_valid_timestep(dt)
        real(kind=wp), intent(in) :: dt

        is_valid_timestep = dt > 0.0_wp
    end function is_valid_timestep

    subroutine time_defaults
        t = 0.0_wp
        time_step = 0.001_wp
        iter = 0
        orig_time_step = 0.001_wp
        walltime=0.0_wp
    end subroutine time_defaults

    subroutine parse_timestep(line)
        character(len=*), intent(in) :: line
        real(kind=wp) :: parsed_time_step
        integer :: iospara

        call read_timestep_value(line, parsed_time_step, iospara)
        if (iospara > 0) then
            call misc_error("failed to parse timestep command")
        end if

        call validate_timestep(parsed_time_step, misc_error)

        time_step = parsed_time_step
        orig_time_step = parsed_time_step
    end subroutine parse_timestep

    subroutine read_timestep_value(line, parsed_time_step, iospara)
        character(len=*), intent(in) :: line
        real(kind=wp), intent(out) :: parsed_time_step
        integer, intent(out) :: iospara
        character(len=read_len) :: tmptxt

        read(line, *, iostat = iospara) tmptxt, parsed_time_step
    end subroutine read_timestep_value

    subroutine validate_timestep(dt, on_error)
        real(kind=wp), intent(in) :: dt
        procedure(timestep_error_handler) :: on_error

        if (.not. is_valid_timestep(dt)) then
            call on_error("time step must be greater than 0")
        end if
    end subroutine validate_timestep

    subroutine log_time
        character(len=read_len) :: msg

        program_end = mpi_wtime()
        write(msg, *) walltime(1), " spent neighboring and ", walltime(2), " spent computing forces"
        call log_msg(msg)
        write(msg, *) "Total time is:", program_end - program_start 
        call log_msg(msg)
    end subroutine log_time

end module time
