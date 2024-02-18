module displace
    !This module contains subroutines which apply displacements to atoms and elements
    use forces
    use parameters 
    use elements
    use group
    use str
    use time
    use neighbors

    implicit none 

    public 
    contains
    subroutine displace_points(line)
        !This subroutine parses the displace command
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, g, txtvec(3)
        integer :: iospara, i, ia, ip, ibasis, ie, gnum
        real(kind=wp) :: disp(3) 

        if((ele_num == 0).and.(atom_num == 0)) call misc_error("Model must be read in before calling displace command")

        read(line, *, iostat=iospara) tmptxt, g, (txtvec(i), i = 1, 3)
        if(iospara > 0) call read_error("Failure to read displace command", iospara)

        !Get group number
        gnum = get_group_index(g)
        if(gnum == 0) then 
            call misc_error("Group "//trim(adjustl(g))// " in displace command has not been defined. "// &
                            "Please define group before use")
        end if

        !Read displacement vector
        do i =1, 3
            read(txtvec(i), *, iostat = iospara) disp(i)
            if(iospara > 0) call read_error("Failure to read disp vector component "//trim(adjustl(txtvec(i))), iospara)
        end do

        !Now displace points
        do ia = 1, atom_num_l
            if(btest(a_mask(ia), gnum)) then 
                r_atom(:,ia) = r_atom(:, ia) + disp
            end if
        end do
        do ip = 1, node_num_l
            ie = node_cg(ip)
            if(btest(e_mask(ie), gnum)) then
                do ibasis = 1, basis_num(ie)
                    r(:, ibasis, ip) = r(:, ibasis, ip) + disp
                end do
            end if
        end do

        call update_neighbor(iter, .false., .true.)

    end subroutine displace_points

    subroutine ramp_displace(line)
        !This subroutine parses the displace command
        character(len = *), intent(in) :: line
        character(len=read_len) :: tmptxt, g, txtvec(4), rdim, ddim
        integer :: iospara,  ia, ip, ibasis, ie, gnum, j, k
        real(kind=wp) :: rhi, rlo, disp, hi, lo, frac

        if((ele_num == 0).and.(atom_num == 0)) call misc_error("Model must be read in before calling ramp command")

        read(line, *, iostat=iospara) tmptxt, g, ddim, txtvec(1), txtvec(2), rdim, rlo, rhi
        if(iospara > 0) call read_error("Failure to read ramp command", iospara)

        !Get the ramp dimensions
        select case(ddim)
        case('x','X')
            j = 1
        case('y','Y')
            j = 2
        case('z','Z')
            j = 3
        end select

        select case(rdim)
        case('x','X')
            k = 1
        case('y','Y')
            k = 2
        case('z','Z')
            k = 3
        end select
        
        !Get the displacement magnitude
        disp = rhi-rlo

        !Parse the displacement hi and lo positions
        call parse_pos(j, txtvec(1), lo)
        call parse_pos(j, txtvec(2), hi)

        !Get group number
        gnum = get_group_index(g)
        if(gnum == 0) then 
            call misc_error("Group "//trim(adjustl(g))// " in ramp command has not been defined. "// &
                            "Please define group before use")
        end if
        
        !Now displace points
        do ia = 1, atom_num_l
            if(btest(a_mask(ia), gnum)) then 
                if(r_atom(j,ia) > hi) then 
                    frac = 1.0_wp
                else if(r_atom(j,ia) < lo) then 
                    frac = 0.0_wp
                else
                    frac = (r_atom(j,ia) - lo)/(hi-lo)
                end if
                r_atom(k,ia) = r_atom(k, ia) + frac*disp + rlo
            end if
        end do
        do ip = 1, node_num_l
            ie = node_cg(ip)
            if(btest(e_mask(ie), gnum)) then
                do ibasis = 1, basis_num(ie)
                    if(r(j,ibasis,ip) > hi) then 
                        frac = 1.0_wp
                    else if(r(j,ibasis,ip) < lo) then 
                        frac = 0.0_wp
                    else
                        frac = (r(j, ibasis, ip) - lo)/(hi-lo)
                    end if
                    r(k, ibasis, ip) = r(k, ibasis, ip) + frac*disp + rlo
                end do
            end if
        end do

        if(nei_init) then 
            call update_neighbor(iter, .false., .true.)
        end if

        return
    end subroutine ramp_displace

end module displace
