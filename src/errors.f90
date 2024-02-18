module errors
    !This subroutines is in charge of handling error strings and exiting the application
    use mpi
    use parameters
    use logger

    implicit none 
    public 
    contains

    subroutine alloc_error(error_string, allostat)
        integer, intent(in) :: allostat
        character(len=*), intent(in) :: error_string
        character(len=read_len) :: msg
        
        write(msg, *) "Allocate Error: ", trim(adjustl(error_string)), " with allostat ", allostat
        call exit_error(msg)
        return
    end subroutine alloc_error

   subroutine bounds_error(error_string, r, bounds)
        !Atom unexpectedly outside boundary, usually box boundary
        real(kind=wp), intent(in) :: r(3), bounds(6)
        character(len=*), intent(in) :: error_string
        character(len=read_len) :: msg

        write(msg,*) "Bounds Error: ", trim(adjustl(error_string)), " for position ", r, " and bounds", bounds
        call exit_error(msg)
        return
    end subroutine bounds_error

    subroutine read_error(error_string, stat)
        !Error when parsing string
        character(len = *), intent(in) :: error_string
        integer, intent(in) :: stat
        character(len=read_len) :: msg
        
        write(msg,*) "Read Error: ", trim(adjustl(error_string)), " with stat ", stat
        call exit_error(msg)
        return
    end subroutine read_error

    subroutine misc_error(error_string)
        !Misc_error
        character(len=*), intent(in) :: error_string
        character(len=read_len) :: msg

        write(msg,*) "Error: ", trim(adjustl(error_string))
        call exit_error(msg)
        return
    end subroutine misc_error

    subroutine command_error(error_string)
        !Error with incorrect command 
        character(len = *) :: error_string
        character(len=read_len) :: msg

        write(msg,*) "Command Error: ", trim(adjustl(error_string))
        call exit_error(msg)
        return 
    end subroutine command_error

    subroutine exit_error(msg)
        !Exit code with error
        character(len=*), intent(in) :: msg

        call log_msg(msg, 0, .true.)
        call close_log
        call mpi_abort(mpi_comm_world, 1, ierr)
        return
    end subroutine exit_error

end module errors
