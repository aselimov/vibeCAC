module berendsen
    use elements
    use parameters
    use compute
    use forces
    use time

    implicit none

    real(kind=wp), save, private:: Ptarget(3), time_constant
    logical, save :: pflag, need_barostat(3)
    public 
    contains

    subroutine parse_berendsen(line)
        !This subroutine parses the press command
        character(len=*), intent(in) :: line
        character(len = read_len) :: tmptxt, msg, press_targets(3)
        integer :: iospara, i

        read(line, *, iomsg = msg, iostat = iospara) tmptxt, time_constant, (press_targets(i), i = 1, 3)
        if(iospara > 0) call read_error(msg, iospara)

        need_barostat=.false.
        do i = 1,  3
            if ((trim(adjustl(press_targets(i))) == 'null').or.(trim(adjustl(press_targets(i))) == 'NULL')) then 
                Ptarget(i)=0
            else
                read(press_targets(i), *) Ptarget(i)
                need_barostat(i)=.true.
            end if
        end do
        pflag = .true.

        write(msg, *) "Using Berendsen barostat to target pressure ", Ptarget, " for dimensions", need_barostat
        call log_msg(msg)
        return
    end subroutine parse_berendsen

    subroutine rescale_box
        !This function actually rescales the box for barostatting.
        integer :: i, j
        real(kind = wp) :: pressure(3,3), identity(3,3), mu(3,3), cntr(3), cdamp

        !First compute the box pressure
        pressure = compute_box_press()

        !Initialize some needed variables 
        identity(:,:) = 0.0_wp
        forall(i=1:3) identity(i,i) = 1.0_wp
        mu(:,:) = 0.0_wp

        !Now calculate the scaling factor
        cdamp = time_step/time_constant
        forall(i=1:3)mu(i,i) = (identity(i,i) - cdamp*(Ptarget(i) - pressure(i,i)))**(1.0_wp/3.0_wp)
        !We set the scaling factor  to 1 in non-periodic dimensions
        do i = 1, 3
            if ((.not.period(i)).or.(.not.need_barostat(i))) mu(i,i) = 1.0_wp
        end do
        !Now get the box center
        cntr = 0.0_wp
        do i = 1,3
            cntr(i) = box_length(i)/2 + box_bd(2*i-1)
        end do
        !Now scale the box boundaries
        do i = 1, 6
            j = (i+1)/2
            box_bd(i) = cntr(j) + mu(j,j)*(box_bd(i)-cntr(j)) 
        end do

        !update box length
        do i = 1, 3
            box_length(i) = box_bd(2*i) - box_bd(2*i-1)
        end do

        !Now scale the atomic and nodal positions
        if(atom_num > 0) then
            do i = 1, atom_num_l
                do j = 1, 3
                    r_atom(j,i) = cntr(j) + mu(j,j)*(r_atom(j,i)-cntr(j))
                end do
            end do
        end if
        if(ele_num > 0) then 
            do i = 1, node_num_l
                do j=1, 3
                    r(j,:,i) = cntr(j) + mu(j,j)*(r(j,:,i)-cntr(j))
                end do
            end do
        end if

        !Now scale the box matrix
        box_mat = matmul(box_mat,mu)

        !Now update the processor boundaries
        pro_length(:) = 0.0_wp
        pro_bd(:) = 0.0_wp

        do i = 1, 3
            !Define processor boundaries
            pro_length(i) = box_length(i) / num_pro(i)
            pro_bd(2*i-1) = box_bd(2*i-1) + grid_coords(i) * pro_length(i)
            pro_bd(2*i) = box_bd(2*i-1) + (grid_coords(i) + 1) * pro_length(i)

            if(grid_coords(i) == 0) then
                !If the processor has coordinate 0 it's the bottom one so set the boundary equal to the box boundary
                pro_bd(2*i-1) = box_bd(2*i-1)
            end if

            if(grid_coords(i) == num_pro(i) - 1) then
                !If the processor has coordinate num_pro(i) - 1 it's the top one so set the boundary equal to top boundary
                pro_bd(2*i) = box_bd(2*i)
            end if

            pro_length(i) = pro_bd(2*i) - pro_bd(2*i-1)
        end do
        call processor_bds

        return
    end subroutine rescale_box

end module berendsen
