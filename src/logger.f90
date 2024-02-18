module logger
    !This is in charge of logging to stdout and to log files
    use parameters

    logical,save :: log_on

    public
    contains


    subroutine log_defaults
        log_on=.true.
    end subroutine log_defaults

    subroutine init_log
        !Open log file
        if (log_on) then 
            if (rank == root) then 
               open(unit=51, file="cac.log", action='write', status='replace', position='rewind')
            end if
        end if

    end subroutine init_log

    subroutine log_msg(msg, isline, all_procs)
        !Log to stdout and log file
        character(len=*), intent(in) :: msg
        integer, intent(in), optional :: isline !if 0 then no tab in front
        logical, intent(in), optional :: all_procs !if .true. then all processors calling this write
        integer :: tab
        logical :: lognow

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
                    if(rank==root) write(51, "(a)") "    ", trim(adjustl(msg))

                !Write line information which isn't indented 
                else if(tab == 0) then 
                    print *, trim(adjustl(msg))
                    if(rank==root) write(51, "(a)") trim(adjustl(msg))
                    
                end if
            end if
        end if

        return
    end subroutine log_msg
    
    subroutine close_log
        !Close the log file
        if (log_on) then 
            if(rank == root)  close(51)
        end if
    end subroutine close_log
end module logger
