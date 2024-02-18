module deform
    !This module applies deformations to the box which is then applied to all of the atoms

    use parameters
    use box
    use elements 
    use neighbors
    use time
    use comms

    implicit none
    
    integer :: deform_num
    integer, private, parameter :: max_deform_num = 10
    real(kind=wp),private :: rate(max_deform_num)
    integer, private :: ddim(max_deform_num), nnevery(max_deform_num), exclude_gnum(max_deform_num)
    
    public 
    contains

    subroutine parse_deform(line)
        !Parse the deform command
        character(len=*), intent(in) :: line
        character(len=read_len) :: tmptxt, deform_dim, msg, excludeg, args(20)
        integer :: iospara, i, d, j, now_excludeg
        real(kind = dp) :: strain, def
        real(kind=dp), allocatable :: r_zero_atom(:, :), r_zero_node(:,:,:)
        
        
        read(line, *, iostat = iospara, iomsg = msg) tmptxt, tmptxt
        if(iospara > 0) call read_error(msg, iospara)
        select case(tmptxt)
        case('now')
            now_excludeg=0
            read(line, *, iostat= iospara, iomsg = msg) tmptxt, tmptxt, deform_dim, strain
            if(iospara > 0) call read_error(msg, iospara)
            j=tok_count(line)
            if (j>4) then
                read(line, *) args(1:j)
                i=5
                do while(i<=j)
                    select case(args(i))
                        case('exclude')
                            i=i+1
                            read(args(i), *) excludeg
                            i=i+1
                            now_excludeg = get_group_index(excludeg)
                            allocate(r_zero_atom(3,atom_num_l), r_zero_node(3, max_basisnum, node_num_l))
                            if(atom_num_l > 0) r_zero_atom = r_atom
                            if(ele_num_l > 0) r_zero_node = r
                            write(msg, *) "Excluding group ", trim(adjustl(excludeg)), ' from deform command' 
                            call log_msg(msg)
                        case default
                            write(msg, *) "Keyword ", trim(adjustl(args(i))), " is not accepted in neighbor command"
                            call command_error(msg)
                    end select
                end do
            end if
            !Now convert the dimension to an integer for use later
            select case(deform_dim)
            case('x','X')
                d = 1
            case('y','Y')
                d = 2
            case('z','Z')
                d = 3
            end select
            !First we convert atoms and nodes to fractional coordinates if necessary
            call x2frac

            def = strain*orig_box_length(d)
            box_bd(2*d-1) = box_bd(2*d-1) - 0.5_wp*def
            box_bd(2*d) = box_bd(2*d) + 0.5_wp*def
            box_length(d) = box_bd(2*d) - box_bd(2*d-1)
            !Now update proc bds
            do i = 1, 3
                pro_length(i) = box_length(i) / num_pro(i)
                pro_bd(2*i-1) = box_bd(2*i-1) + grid_coords(i) * pro_length(i)
                pro_bd(2*i) = box_bd(2*i-1) + (grid_coords(i) + 1) * pro_length(i)
            end do
            call processor_bds        
            !Convert back to real coordinates if we applied a deformation
            call frac2x

            if(now_excludeg > 0) then 
                do i = 1, atom_num_l
                    if(btest(a_mask(i),now_excludeg)) then 
                        r_atom(:, i)=r_zero_atom(:,i)
                    end if
                end do
                do i = 1, ele_num_l
                    if(btest(e_mask(i), now_excludeg)) then 
                        r(:,:,i) = r_zero_node(:,:,i) 
                    end if
                end do
            end if

            write(msg, *) "Increasing box bounds in ", trim(adjustl(deform_dim)), " by ", def
            call log_msg(msg)

            call update_neighbor(iter, .false.,.true.)

        case default
            deform_num=deform_num + 1
            read(line, *, iostat= iospara, iomsg = msg)tmptxt, nnevery(deform_num), deform_dim, rate(deform_num)
            if(iospara > 0) call read_error(msg, iospara)

            !Now convert the dimension to an integer for use later
            select case(deform_dim)
            case('x','X')
                ddim(deform_num) = 1
            case('y','Y')
                ddim(deform_num) = 2
            case('z','Z')
                ddim(deform_num) = 3
            end select
        end select

        return

    end subroutine parse_deform

    subroutine deform_box
        !This subroutine updates the box boundaries and then rescales the atomic positions
        
        integer :: i
        real(kind=wp) :: def
        logical :: convert 
        character(len=1) :: def_dim(3)
        character(len=read_len) :: msg

        def_dim = ['x','y','z']


        !First update the box length
        !Now update the box boundaries
        convert = .false.
        do i = 1, deform_num

            if(mod(iter, nnevery(i)) == 0) then 
                !First c
                if(.not.convert) then 
                    !First we convert atoms and nodes to factional coordinates if necessary
                    call x2frac
                    convert = .true.
                end if

                def = time_step*rate(i)*orig_box_length(ddim(i))
                box_bd(2*ddim(i)-1) = box_bd(2*ddim(i)-1) - 0.5_wp*def
                box_bd(2*(ddim(i))) = box_bd(2*(ddim(i))) + 0.5_wp*def
                box_length(ddim(i)) = box_bd(2*ddim(i)) - box_bd(2*ddim(i)-1)

                write(msg, *) "Increasing box bounds in ", trim(adjustl(def_dim(ddim(1)))), " by ", def, " for deform ", i
                call log_msg(msg)
            end if
        end do


        !Now update proc bds
        if(convert) then
            do i = 1, 3
                pro_length(i) = box_length(i) / num_pro(i)
                pro_bd(2*i-1) = box_bd(2*i-1) + grid_coords(i) * pro_length(i)
                pro_bd(2*i) = box_bd(2*i-1) + (grid_coords(i) + 1) * pro_length(i)
            end do
            call processor_bds        

            !Convert back to real coordinates if we applied a deformation
             call frac2x

            call update_neighbor(iter, .false.,.true.)
        end if
        return
    end subroutine deform_box

end module deform
