
module box
    !This module contains all information related to simulation box definitions
   
    use parameters
    use errors

    implicit none
    
    !Box definition variables
    real(kind=wp) :: box_bd(6), orig_box_bd(6), box_length(3), orig_box_length(3)
    real(kind=wp) :: box_mat_ori(3,3), box_mat(3,3)
    !Boundary definitions
    character(len=3) :: boundary
    logical :: periodic
    logical, dimension(3) :: period

    real(kind=wp), dimension(3) :: center_mass, center_mass_ori, center_mass_disp

    public 
    contains

    subroutine parse_boundary(line)
        !This subroutine parses the optional boundary command
        character(len=*), intent(in) :: line
        integer :: iospara, i
        character(len=read_len) :: msg, tmptxt
        read(line, *, iostat = iospara, iomsg = msg) tmptxt, boundary
        if(iospara > 0) call read_error(msg, iospara)

        do i = 1, 3
            select case(boundary(i:i))
            case('s', 'S')
                period(i) = .false.
            case('p','P')
                period(i) = .true.
            case default
                call command_error("Boundary must be either p or s")
            end select
        end do

        if(any(period)) then
            periodic = .true.
            if(.not.allocated(pb_node))then
                allocate(pb_node(3,max_basisnum, node_num_lr))
            end if
        else 
            periodic = .false.
        end if


    end subroutine parse_boundary

    subroutine cross_pb(r, info)
        !This subroutine sets which periodic boundaries the nodes have crossed
        real(kind = wp), dimension(3), intent(inout) :: r
        integer, dimension(3), intent(out) :: info

        integer :: i

        info(:) = 0
        do i = 1, 3
            if(period(i).eqv..true.) then
                if(r(i) < box_bd(i*2-1)) then
                    r(i) = r(i) + box_length(i)
                    info(i) = 1

                else if(r(i) > box_bd(i*2)) then
                    r(i) = r(i) - box_length(i)
                    info(i) = -1

                end if
            end if
        end do
        return 
    end subroutine cross_pb

    subroutine restore_pb(r_in, pb_move, ipb)

        real(kind = wp), dimension(3), intent(inout) :: r_in
        integer, dimension(3), intent(in) :: pb_move
        logical, optional, intent(out) :: ipb

        integer :: i

        if(present(ipb)) ipb = .false.
        do i = 1, 3
            if((period(i).eqv..true.).and. (pb_move(i) /= 0)) then
                r_in(i) = r_in(i) - pb_move(i) * box_length(i)
                if(present(ipb)) ipb = .true.
            end if
       end do

       return
    end subroutine restore_pb

    function in_box_bd(r)
        !this function just checks to make sure position r is within the overall simulation box
        real(kind=wp), intent(in) :: r(3)
        logical :: in_box_bd

        integer :: i
        in_box_bd = .true.
        do i = 1,3
            if (r(i) < box_bd(2*i-1)) then 
                !Lower bound
                in_box_bd = .false.
                return
            else if (r(i) > box_bd(2*i)) then 
                !Upper bound
                in_box_bd = .false.
                return
            end if
        end do

        return
    end function in_box_bd

    recursive subroutine parse_pos(i, pos_string, pos)
        !This subroutine parses the pos command allowing for command which include inf
        integer, intent(in) ::  i !The dimension of the position
        character(len=*), intent(in) :: pos_string !The position string
        real(kind=dp), intent(out) :: pos !The output parsed position value

        integer :: iospara
        real(kind=dp) :: rand, rone, rtwo
        character(len=100) :: cone, ctwo

        iospara = 0
        
        if(trim(adjustl(pos_string)) == 'inf') then 
            pos=box_bd(2*i)

        else if(trim(adjustl(pos_string)) == '-inf') then
            pos=box_bd(2*i-1)

        else if (trim(adjustl(pos_string)) == 'rand') then 
            call random_number(rand)
            pos = (box_bd(2*i)-box_bd(2*i-1))*rand + box_bd(2*i-1)

        else if (index(pos_string,'rand')>0) then 
            call random_number(rand)
            cone = pos_string(index(pos_string, '[')+1:index(pos_string,':')-1)
            call parse_pos(i, cone, rone)
            ctwo = pos_string(index(pos_string, ':')+1:index(pos_string,']')-1)
            call parse_pos(i, ctwo, rtwo)
            pos = (rtwo - rone)*rand + rone            

        else if ((index(pos_string,'-') > 0).and.(index(pos_string,'inf')>0)) then 
            !Now extract the number we are reducing from infinity
            if(index(pos_string,'inf') < index(pos_string,'-')) then 
                read(pos_string(index(pos_string,'-')+1:), *, iostat=iospara) pos
            else
                read(pos_string(1:index(pos_string,'-')-1), *, iostat=iospara) pos
            end if
            pos = box_bd(2*i) - pos

        else if ((index(pos_string,'+') > 0).and.(index(pos_string,'inf')>0)) then 
            !Now extract the number we are reducing from infinity
           if(index(pos_string,'inf') < index(pos_string,'+')) then 
                read(pos_string(index(pos_string,'+')+1:), *, iostat=iospara) pos
            else
                read(pos_string(1:index(pos_string,'+')-1), *, iostat=iospara) pos
            end if
            pos = box_bd(2*i-1) + pos

        else if ((index(pos_string,'*') > 0).and.(index(pos_string,'inf')>0)) then 
            !Now extract the number we are reducing from infinity
           if(index(pos_string,'inf') < index(pos_string,'*')) then 
                read(pos_string(index(pos_string,'*')+1:), *, iostat=iospara) pos
            else
                read(pos_string(1:index(pos_string,'*')-1), *, iostat=iospara) pos
            end if
            pos = (box_bd(2*i)-box_bd(2*i-1))*pos + box_bd(2*i-1)

        else
            read(pos_string, *, iostat=iospara) pos

        end if

        if (iospara > 0) then 
            print *, "Error reading position argument ", trim(adjustl(pos_string)), ". Please reformat and try again."
        end if
    end subroutine parse_pos

    subroutine update_box_mat
        !This subroutine updates the box matrix 
        integer :: i
        
        box_mat = 0.0_wp
        do i = 1, 3
            box_mat(i,i) = box_length(i)
        end do

        return
    end subroutine update_box_mat
end module box
