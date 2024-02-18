module thermo
    !this is similar to the lammps thermo output which logs simulation information 
    use elements
    use forces
    use compute
    use str

    implicit none


    integer :: thermo_every, thermo_style_num
    character(len=read_len) :: thermo_style(100)
    logical :: need_k, need_t, need_p

    public
    contains

    subroutine thermo_defaults
        thermo_every = huge(1)
        thermo_style_num = 2 
        thermo_style(1) = "pe"
        thermo_style(2) = "fnorm"
        need_k = .false.
        need_t = .false.
        need_p = .false.
    end subroutine thermo_defaults

    subroutine parse_thermo_style(line)
        !This subroutine parses the thermo style command
        character(len=*), intent(in) :: line
        integer :: iospara, i

        character(len=read_len) :: tmptxt

        thermo_style_num = tok_count(line)-1

        need_k = .false.
        need_t = .false.
        need_p = .false.
        read(line, *, iostat=iospara) tmptxt, (thermo_style(i), i=1, thermo_style_num)

        !Now loop over all of the thermo_style items to make sure that they are acceptable
        do i= 1, thermo_style_num
            select case(thermo_style(i))
                case('pe','fnorm', 'lx', 'ly', 'lz')
                    continue
                case('ke','ke_c','ke_a')
                    need_k = .true.
                case('temp', 'temp_c', 'temp_a')
                    need_k = .true.
                    need_t = .true.
                case('px', 'py', 'pz')
                    need_p = .true.
                case default
                    write(tmptxt, *) trim(adjustl(thermo_style(i))), " is not an acceptable thermostyle command, please check " //&
                                   " the documentation for acceptable commands."    
            end select
        end do
    end subroutine parse_thermo_style

    subroutine parse_thermo(line)
        !This subroutine parses the thermo command
        character(len=*) :: line
        character(len = read_len) :: tmptxt  
        integer :: iospara

        read(line, *, iostat = iospara) tmptxt, thermo_every 
        if (iospara > 0) call read_error("Failure to read thermo command", iospara)

        !Make sure thermo_every is greater than 0
        if(thermo_every < 0) call misc_error("Thermo_every should be greater than 0 in thermo command")

        if(thermo_every == 0) then 
            thermo_every = huge(1)
        end if

    end subroutine parse_thermo

    subroutine write_thermo_style
        integer :: i
        character(len=read_len) :: style_out

        write(style_out, *) "thermo: ", (trim(adjustl(thermo_style(i)))//" ", i = 1, thermo_style_num)
        call log_msg(style_out, 0)

    end subroutine write_thermo_style

    subroutine write_thermo_out(time)
        !Preliminary write thermo out
        integer, intent(in) :: time
        integer :: i
        real(kind=wp) :: ke(2), T(3), thermo_out(100), P(3,3)
        character(len=read_len) :: msg

        thermo_out = 0.0_wp

        !Calculate temp, ke, and press if needed
         if(need_k) then 
            if(need_vel) then
                ke = compute_ke(1)
            else
                ke = 0.0_wp
            end if
        end if
        if(need_t) then 
            if(need_vel) then 
                T = compute_temp(1, ke)
            else
                T = 0.0_wp
            end if
        end if       

        if(need_p) P = compute_box_press()

        !Loop over all thermo_style commands and put the desired outputs into thermo_out
        do i= 1, thermo_style_num
            select case(thermo_style(i))
            case('pe')
                !Get pe all
                thermo_out(i) = compute_pe(1)

            case('fnorm')
                !Get fnorm all
                thermo_out(i) = compute_fnorm(1)

            case('ke')
                thermo_out(i) = ke(1) + ke(2)
            case('ke_a')
                thermo_out(i) = ke(1)
            case('ke_c')
                thermo_out(i) = ke(2)

            case('temp')
                thermo_out(i) = T(1)
            case('temp_a')
                thermo_out(i) = T(2)
            case('temp_c')
                thermo_out(i) = T(3)

            case('lx')
                thermo_out(i) = box_bd(2) - box_bd(1)
            case('ly')
                thermo_out(i) = box_bd(4) - box_bd(3)
            case('lz')
                thermo_out(i) = box_bd(6) - box_bd(5)

            case('px')
                thermo_out(i) = P(1,1)
            case('py')
                thermo_out(i) = P(2,2)
            case('pz')
                thermo_out(i) = P(3,3)

            end select

        end do
        
        if(rank == root) then 
            write(msg,*) time, thermo_out(1:thermo_style_num)
            call log_msg(msg)
        end if
        
        return
    end subroutine write_thermo_out
end module thermo
