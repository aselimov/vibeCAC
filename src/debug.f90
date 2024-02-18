module debug
    !This module is used to run various debugging commands which can be used during simulations for additional correctness checks
    use parameters
    use forces
    use elements
    use str

    implicit none

    logical, private :: check_ele_flag
    logical :: dflag
    public 
    contains
    

    subroutine parse_debug(line)
        !This subroutine parses the debug command
        character(len=*), intent(in) :: line 
        integer :: i, j, iospara
        character(len=read_len) :: tmptxt, commands(10), msg

        j = tok_count(line)-1
        dflag=.true.
        
        !Initialize flags
        check_ele_flag=.false.

        read(line, *, iostat=iospara, iomsg = msg) tmptxt, (commands(i), i =1, j)

        !Now set the debug flags
        do i = 1, j
            select case(commands(i))
            case('check_ele')
                check_ele_flag= .true.
            end select
        end do

    end subroutine parse_debug

    subroutine run_debug
        !This is the main running function which calls all of the sub debug subroutines
        if(check_ele_flag) call check_ele
    end subroutine run_debug

    subroutine check_ele
        !This code goes over every single element and compares the positions among all elements. 
        !This code is very slow
        integer :: ie, ke, inod, i, j, k
        real(kind=wp) :: rbuff(3*ng_max_node), rbuffall(3*ng_max_node*pro_num)
        real(kind=wp) :: velbuff(3*ng_max_node), velbuffall(3*ng_max_node*pro_num)
        logical :: have_ele, have_ele_all(pro_num)
        character(len=read_len) :: msg

        rbuff=0.0_wp
        velbuff=0.0_wp
        !Loop over all global tags
        do ke = 1, ele_num
            !Loop over all local tags to find the matching element
            do ie = 1, ele_num_l
                if(ele_glob_id(ie) == ke) exit
            end do

            if (ie > ele_num_l) then 
                have_ele=.false.
                rbuff = 0.0_wp
            else
                have_ele=.true.
                do inod =1, ng_node(etype(ie))
                    do i = 1, 3
                        rbuff(3*(inod-1)+i) = r(i, 1, cg_node(inod,ie))
                        if(need_vel) velbuff(3*(inod-1)+i) = vel(i, 1, cg_node(inod,ie))
                    end do
                end do
            end if

            rbuffall=0.0_wp
            velbuffall=0.0_wp
            !Now communicate this element 
            call mpi_gather(have_ele, 1, mpi_logical, have_ele_all, 1, mpi_logical, root, world, ierr)
            call mpi_gather(rbuff, 3*ng_max_node, mpi_wp, rbuffall, 3*ng_max_node, mpi_wp, root, world, ierr)
            if(need_vel) call mpi_gather(velbuff, 3*ng_max_node, mpi_wp, velbuffall, 3*ng_max_node, mpi_wp, root, world, ierr)

            !Now check the elements
            if(rank==root) then 
                if(count(have_ele_all) > 1) then 
                    do i = 1, pro_num
                        if(have_ele_all(i)) then 
                            do j=1, pro_num
                                if((j>i) .and. (have_ele_all(j))) then 
                                    do k = 1, 3*ng_max_node
                                        if(.not.is_equal(rbuffall((i-1)*ng_max_node*3 + k), &
                                                         rbuffall((j-1)*ng_max_node*3 + k))) then
                                            write(msg, *) "Element with global tag ", ke, " does not have matching nodal ", &
                                                          "positions for processors ", i-1, " and ", j-1, " with values ", &
                                                          rbuffall((i-1)*ng_max_node*3 + k), " and ", &
                                                          rbuffall((j-1)*ng_max_node*3 + k)
                                            call misc_error(msg)

                                        end if
                                        if(need_vel) then 
                                            if(.not.is_equal(velbuffall((i-1)*ng_max_node*3 + k), &
                                                             velbuffall((j-1)*ng_max_node*3 + k))) then
                                                write(msg, *) "Element with global tag ", ke, " does not have matching nodal ", &
                                                              "velocities for processors ", i-1, " and ", j-1, " with values ", &
                                                              velbuffall((i-1)*ng_max_node*3 + k), " and ", &
                                                              velbuffall((j-1)*ng_max_node*3 + k)
                                                call misc_error(msg)
                                            end if
                                        end if
                                    end do
                                end if
                            end do
                        end if
                    end do
                end if
            end if
        end do
    end subroutine check_ele
end module debug
