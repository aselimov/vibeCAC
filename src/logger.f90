module logger
    !This is in charge of logging to stdout and to log files
    use parameters

    integer, parameter :: default_log_unit = 51
    logical,save :: log_on

    public
    contains


    subroutine log_defaults
        log_on=.true.
    end subroutine log_defaults

    subroutine init_log(log_filename, log_unit)
        !Open log file
        character(len=*), intent(in), optional :: log_filename
        integer, intent(in), optional :: log_unit
        character(len=256) :: filename
        integer :: unit_out

        filename = "cac.log"
        unit_out = default_log_unit
        if (present(log_filename)) filename = trim(adjustl(log_filename))
        if (present(log_unit)) unit_out = log_unit

        if (log_on) then 
            if (rank == root) then 
               open(unit=unit_out, file=trim(adjustl(filename)), action='write', status='replace', position='rewind')
            end if
        end if

    end subroutine init_log

    subroutine log_msg(msg, isline, all_procs, log_unit)
        !Log to stdout and log file
        character(len=*), intent(in) :: msg
        integer, intent(in), optional :: isline !if 0 then no tab in front
        logical, intent(in), optional :: all_procs !if .true. then all processors calling this write
        integer :: tab
        integer :: unit_out
        logical :: lognow
        integer, intent(in), optional :: log_unit

        unit_out = default_log_unit
        if (present(log_unit)) unit_out = log_unit

        if(log_on) then 
            lognow = .false.
            if(present(all_procs)) then 
                lognow = all_procs
            end if
            if(rank == root) lognow = .true.

            if(present(isline)) then 
                tab = isline
            else 
                tab = 1
            end if

            !write non-line information, which should be indented
            if(lognow) then 
                if(tab == 1) then 
                    print *, "    ", trim(adjustl(msg))
                    if(rank==root) write(unit_out, "(a)") "    ", trim(adjustl(msg))

                !Write line information which isn't indented 
                else if(tab == 0) then 
                    print *, trim(adjustl(msg))
                    if(rank==root) write(unit_out, "(a)") trim(adjustl(msg))
                    
                end if
            end if
        end if

        return
    end subroutine log_msg
    
    subroutine close_log(log_unit)
        !Close the log file
        integer, intent(in), optional :: log_unit
        integer :: unit_out

        unit_out = default_log_unit
        if (present(log_unit)) unit_out = log_unit

        if (log_on) then 
            if(rank == root)  close(unit_out)
        end if
    end subroutine close_log
end module logger
