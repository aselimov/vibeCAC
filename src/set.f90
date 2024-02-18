module set
    ! sets per-atom values as required for random alloys, temperature initializiation
    use elements
    use group
    use logger

    implicit none
    public
    contains
    
    subroutine parse_set(line)
    ! parse the 'set' command and set appropriate per-atom vals
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, g, mode, msg
        integer :: iospara, gnum, j, i
        real(kind=wp) :: args(20)


        if((ele_num == 0).and.(atom_num == 0)) call misc_error("Model must be read in before calling set command")
        
        read(line, *, iostat=iospara) tmptxt, g, mode
        if(iospara > 0) call read_error("Failure to read set command", iospara)

        gnum = get_group_index(g)
        if(gnum == 0) then 
            call misc_error("Group "//trim(adjustl(g))// " in set command has not been defined. "// &
                            "Please define group before use")
        end if

        ! read in remaining arguments
        j = tok_count(line)-3
        if (j > 0) then
            read(line, *, iostat = iospara) tmptxt, g, mode, (args(i), i = 1, j)
        end if       

        select case(mode)
        case('type')
            write(msg,*) "set ", group_counts(1,gnum), " atoms and ", group_counts(2,gnum), " elements to type: ", args(1)
            call log_msg(msg)
            call set_type(gnum, int(args(1)))
        case('fraction')
            if(group_counts(2,gnum) > 0) then
                call command_error("Set fraction cannot operate on group "//trim(adjustl(group_name(group_num))) &
                                //" containing elements")
            endif
            call set_fraction(gnum, int(args(1)), args(2))
        case('etype')
            if(group_counts(1,gnum) > 0)  call command_error("Set etype cannot operate on groups containing atoms")
            call set_etype(gnum, int(args(1)))
        case('vel')
            if(allocated(vel_atom)) then 
                write(msg,*) "set ", group_counts(1,gnum), " atoms and ", group_counts(2,gnum), " elements to vel ", args(1:3)
                call set_vel(gnum, real(args(1),wp),real(args(2),wp), real(args(3),wp))
            else
                write(msg,*) "Velocity arrays are not allocated, velocity cannot be set"
                call log_msg(msg)
            end if
            call log_msg(msg)
        end select
    end subroutine parse_set

    subroutine set_type(gnum, type)
        integer, intent(in) :: gnum, type
        integer :: ia, ip, ibasis, ie

        do ia = 1, atom_num_l
            if(btest(a_mask(ia), gnum)) then 
                type_atom(ia) = type
            end if
        end do
        do ip = 1, node_num_l
            ie = node_cg(ip)
            if(btest(e_mask(ie), gnum)) then
                do ibasis = 1, basis_num(ie)
                    basis_type(ibasis,ie) = type
                end do
            end if
        end do
    end subroutine set_type

    subroutine set_vel(gnum, vx, vy, vz)
        integer, intent(in) :: gnum
        real(kind = wp), intent(in) :: vx, vy, vz
        integer ::ia, ip, ibasis, ie
        real(kind=wp) :: v(3)

        v(1) = vx
        v(2) = vy
        v(3) = vz

        do ia = 1, atom_num_l
            if(btest(a_mask(ia), gnum)) then 
                vel_atom(:, ia) = v 
            end if
        end do
        do ip = 1, node_num_l
            ie = node_cg(ip)
            if(btest(e_mask(ie), gnum)) then
                do ibasis = 1, basis_num(ie)
                    vel(:,ibasis,ip) = v
                end do
            end if
        end do
    end subroutine set_vel

    subroutine set_fraction(gnum, type, frac)
        integer, intent(in) :: gnum, type
        real(kind=wp), intent(in) :: frac
        character(len = read_len) :: msg
        real(kind=wp) :: rand
        integer :: ia, count
        ! allocate fraction of group_counts array for randomization
        count = 0
        do ia = 1, atom_num_l
            if(btest(a_mask(ia), gnum)) then 
                call random_number(rand)
                if (rand > frac) cycle
                type_atom(ia) = type
                count = count + 1
            end if
        end do
        write(msg, *) "Set fraction swapped ", count, " atoms to type ", type
        call log_msg(msg)
        !group_counts(1,gnum)
    end subroutine set_fraction

    subroutine set_etype(gnum, et)
        integer, intent(in) :: gnum, et
        
        integer :: i
        do i = 1, ele_num_l
            if (btest(e_mask(i), gnum)) then 
                etype(i) = et
            end if
        end do

        return
    end subroutine set_etype

end module set
