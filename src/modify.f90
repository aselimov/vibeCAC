module modify
    !This module is in charge of code which changes various aspects of current running program, mainly regarding options

    use parameters
    use neighbors 
    use elements 
    use potential

    implicit none
    
    public 
    contains

    subroutine parse_modify(line)
        character(len=*), intent(in) :: line

        integer :: i, acount, iospara
        character(len=read_len) :: args(20), msg
        logical :: l

        acount = tok_count(line)

        read(line,*) (args(i), i=1, acount)

        i = 2
        do while (i < acount)
            select case(args(i))
            case('atomap_neighbors')
                i=i+1
                read(args(i), *, iostat=iospara, iomsg=msg) l
                if(iospara > 0) call read_error(msg, iospara)
                if(l) then 
                    atomap_neighbors = .true.
                else 
                    atomap_neighbors = .false.
                end if

            case default
                write(msg, *) "Argument ", args(i), " not accepted in modify command"
            end select

            i = i + 1
        end do

        
    end subroutine parse_modify

end module modify
