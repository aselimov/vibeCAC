module force_mod 
    use forces
    use parameters
    use elements
    use group
    use str
    use time

    implicit none

    integer, parameter, private :: max_force = 20
    integer :: set_force_num, add_force_num

    integer, private :: force_group(max_force), add_force_group(max_force), add_force_every(max_force)
    logical, private :: force_mask(3, max_force), add_force_mask(3, max_force), add_force_norm(max_force), &
                        first_add(max_force)
    real(kind=wp), private :: fvec(3, max_force), add_f(3, max_force)

    public 
    contains

    subroutine init_set_force
        set_force_num = 0
        force_group = 0
        force_mask = .false.
        first_add = .false.
        fvec = 0.0_wp
    end subroutine

    subroutine init_add_force
        add_force_num = 0
        add_force_group = 0
        add_force_mask = .false.
        add_force_every(:) = 0
        add_f  = 0.0_wp
    end subroutine init_add_force

    subroutine parse_set_force(line)
        !This subroutine parses the set_force command
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, g, txtvec(3)
        integer :: iospara, i

        set_force_num = set_force_num + 1
        first_add(set_force_num) = .true.
        read(line, *, iostat=iospara) tmptxt, g, (txtvec(i), i = 1, 3)
        if(iospara > 0) call read_error("Failure to read set_force command", iospara)

        !Get group number
        force_group(set_force_num) = get_group_index(g)
        if(force_group(set_force_num) == 0) then 
            call misc_error("Group "//trim(adjustl(g))// " in set_force command has not been defined. "// &
                            "Please define group before use")
        end if

        !Read force_vec and force_mask
        do i =1, 3
            call to_lower(txtvec(i))
            if(.not.(txtvec(i) == 'null')) then 
                force_mask(i, set_force_num) = .true.
                read(txtvec(i), *, iostat = iospara) fvec(i, set_force_num)
                if(iospara > 0) call read_error("Failure to read force vector component "//trim(adjustl(txtvec(i))), iospara)
            end if
        end do

    end subroutine parse_set_force

    subroutine parse_add_force(line)
        character(len=*), intent(in) :: line
        character(len=read_len) :: tmptxt, g, txtvec(3), extra_commands(20), msg
        integer :: iospara, i, j

        add_force_num = add_force_num + 1
        read(line, *, iostat=iospara) tmptxt, g, (txtvec(i), i = 1, 3)
        if(iospara > 0) call read_error("Failure to read add_force command", iospara)

        !Get group number
        add_force_group(add_force_num) = get_group_index(g)
        if(add_force_group(add_force_num) == 0) then 
            call misc_error("Group "//trim(adjustl(g))// " in add_force command has not been defined. "// &
                            "Please define group before use")
        end if

        !Read force_vec and force_mask
        do i =1, 3
            call to_lower(txtvec(i))
            if(.not.(txtvec(i) == 'null')) then 
                add_force_mask(i, add_force_num) = .true.
                read(txtvec(i), *, iostat = iospara) add_f(i, add_force_num)
                if(iospara > 0) call read_error("Failure to read force vector component "//trim(adjustl(txtvec(i))), iospara)
            end if
        end do

        !Now get parse any extra optional commands
        j = tok_count(line)-5
        if( j > 0) then 
            read(line, *, iostat=iospara, iomsg=msg) (tmptxt, i = 1, 5), (extra_commands(i), i = 1, j)
            i = 1
            do while (i<=j)
                select case(extra_commands(i))
                case('every')
                    i = i+1
                    read(extra_commands(i), *, iostat=iospara, iomsg=msg) add_force_every(add_force_num)
                    if(iospara>0) call read_error(msg, iospara)
                case('norm')
                    add_force_norm(add_force_num) = .true.
                end select
                i=i+1
            end do
        end if

    end subroutine parse_add_force

    subroutine add_force

        integer :: i, ip, ia, ie, j
        real(kind=wp) :: f(3)
        logical :: do_now
        
        do i = 1, add_force_num
            !Check to see if we need to add force now
            do_now = .false.
            if(first_add(i)) then 
                do_now = .true.
                first_add(i) = .false.
            else if(add_force_every(i) == 0) then 
                do_now = .false.
            else if (mod(iter, add_force_every(i)) == 0) then 
                do_now = .true.
            end if

            if (do_now) then 
                f = add_f(:,i)
                !Check to see if the force is normalized, if so make sure the group contains only atoms
                if(add_force_norm(i)) then 
                    if(group_counts(2,add_force_group(i)) > 0) then 
                        call command_error("Currently cannot normalize force when addforce is used on a group containing elements")
                    end if
                    if(group_counts(1, add_force_group(i)) >0) then 
                        f(:) = f(:)/group_counts(1,add_force_group(i))
                    else
                        f=0.0_wp
                    end if
                end if
                do ia = 1, atom_num_l
                    if(btest(a_mask(ia), add_force_group(i))) then 
                        do j = 1, 3
                            if (add_force_mask(j, i)) then 
                                force_atom(j,ia) = force_atom(j,ia) + f(j)
                            end if
                        end do
                    end if
                end do 
                do ip = 1, node_num_l
                    ie = node_cg(ip)
                    if(btest(e_mask(ie),add_force_group(i))) then 
                        do j = 1, 3
                            if (add_force_mask(j,i)) then 
                                force_eq(j,:, ip) = force_eq(j,:,ip) + f(j)
                            end if
                        end do
                    end if
                end do  
            end if
        end do
        
    end subroutine

    subroutine run_set_force
        integer :: i, j, ia, ip, ie

        do i = 1, set_force_num
            do ia = 1, atom_num_l
                if(btest(a_mask(ia), force_group(i))) then 
                    do j = 1, 3
                        if (force_mask(j, i)) then 
                            force_atom(j,ia) = fvec(j,i)
                        end if
                    end do
                end if
            end do 
            do ip = 1, node_num_l
                ie = node_cg(ip)
                if(btest(e_mask(ie),force_group(i))) then 
                    do j = 1, 3
                        if (force_mask(j,i)) then 
                            force_eq(j,:, ip) = fvec(j,i)
                        end if
                    end do
                end if
            end do
        end do

    end subroutine run_set_force

end module force_mod
