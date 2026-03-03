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
        character(len=read_len) :: tmptxt
        integer :: iospara

        read(line, *, iostat = iospara) tmptxt, time_step

        if(.not.is_valid_timestep(time_step)) then 
            call misc_error("time step must be greater than 0")
        end if

        orig_time_step = time_step
    end subroutine parse_timestep

    subroutine log_time
        character(len=read_len) :: msg

        program_end = mpi_wtime()
        write(msg, *) walltime(1), " spent neighboring and ", walltime(2), " spent computing forces"
        call log_msg(msg)
        write(msg, *) "Total time is:", program_end - program_start 
        call log_msg(msg)
    end subroutine log_time

end module time
