module group

    use parameters
    use elements
    use box
    use str
    use comms

    implicit none

    integer, parameter :: max_group_num = 20
    integer, save :: group_num, group_counts(4,max_group_num)
    integer, dimension(10, max_group_num), save :: group_parameters_int
    real(kind=wp), dimension(10, max_group_num) :: group_parameters_real
    character(len=50), dimension(10, max_group_num) :: group_parameters_string
    character(len=read_len), dimension(max_group_num) :: group_name, group_string, group_shape, which_in_group
    

    public 
    contains 


    subroutine group_defaults
        group_parameters_int = 0
        group_parameters_real=0.0_wp
        group_counts = 0
        group_name(1) = 'all'
        group_num = 1
    end subroutine group_defaults

    subroutine init_group_all
        !This initializes the all group
        integer :: i, group_countme(4)
        group_countme = 0

        if(ele_num > 0) e_mask=0
        if(atom_num > 0) a_mask=0
        group_countme(1) = atom_num_l
        group_countme(2) = ele_num_l
        do i = 1, atom_num_l
            a_mask(i)= ibset(a_mask(i), 1)
        end do
        do i = 1, ele_num_l
            e_mask(i)= ibset(e_mask(i), 1)
            group_countme(3) = group_countme(3) + ng_node(etype(i))*basis_num(i)
            group_countme(4) = group_countme(4) + basis_num(i)*(size_ele(i)+1)**3
        end do
        call mpi_allreduce(group_countme, group_counts(:,1), 4, mpi_integer, mpi_sum, world, ierr)
    end subroutine init_group_all

    function get_group_index(gname)
        !This subroutine returns the index of a group
        character(len=*), intent(in) :: gname
        integer :: get_group_index
        integer :: i

        get_group_index = 0
        do i = 1, group_num
            if (trim(group_name(i)) == trim(gname)) then 
                get_group_index=i  
                exit
            end if
        end do
        return
    end function get_group_index

    subroutine parse_group(line)
        !This subroutine parses the group command
        character(len=*), intent(in) :: line 
        integer :: iospara, i, j, shape_args
        character(len=read_len) :: err, tmptxt
        !First figure out what group_shape
        group_num = group_num + 1
        
        j = tok_count(line)
        
        read(line, *, iomsg=err, iostat=iospara) tmptxt, group_name(group_num), which_in_group(group_num), group_shape(group_num),&
                                                 (group_parameters_string(i,group_num), i = 1, j)
        if(iospara > 0) call read_error(err, iospara)

        !Check to make sure which_in_group is correct
        select case(which_in_group(group_num))
        case('elements', 'atoms', 'all')
            continue
        case default
            call command_error("Group command doesn't accept selection based on "//trim(adjustl(which_in_group(group_num))) &
                             //" should be elements, atoms, or all")
        end select

        !Now make sure that the group_shape is one of the acceptable options and check the inputted parameters
        select case(group_shape(group_num)) 
        case('block', 'sphere', 'type', 'file')
            continue
        case default
            call command_error("Group command doesn't accept shape "//trim(adjustl(group_shape(group_num))) &
                             //" should be block, sphere, or type")
        end select


        !Now parse the group parameters
        select case(group_shape(group_num))
        case('block')
            do i =1, 6
                call parse_pos(int((i+1)/2), group_parameters_string(i, group_num), group_parameters_real(i, group_num))
            end do
            shape_args = 6
        
        case('sphere')
            do i =1, 3
                call parse_pos(i, group_parameters_string(i, group_num), group_parameters_real(i, group_num))
            end do
            !Now get the radius
            read(group_parameters_string(4, group_num), *, iostat = iospara, iomsg = err) group_parameters_real(i, group_num)
            if(iospara > 0) call read_error(err, iospara)
            group_parameters_real(4, group_num) = group_parameters_real(4, group_num) * group_parameters_real(4, group_num)
            shape_args=4

        case('type')
            read(group_parameters_string(1,group_num), *, iostat = iospara, iomsg = err) group_parameters_int(1, group_num)
            if(iospara > 0) call read_error(err, iospara)
            shape_args = 1

        case('file')
            continue
        end select

        call assign_group(group_num)
    end subroutine parse_group
    
    subroutine assign_group(i)
        !This subroutine finds all elements and/or atoms within the group boundaries
        !specified by the user.
        integer, intent(in) :: i
        integer ::  inod, ibasis, ia, ie, ip, group_countme(4), j, reada, reade, intag
        integer, allocatable :: amap(:), emap(:)
        real(kind=dp) :: rc(3)
        character(len = read_len) :: msg

        group_countme = 0
        select case(trim(adjustl(group_shape(i))))
        case('block')
            write(msg, *) "Group ", trim(group_name(i)), " has block shape with boundaries: ", group_parameters_real(1:6, i)
        case('sphere')
            write(msg, *) "Group ", trim(group_name(i)), " has sphere shape with centroid ", group_parameters_real(1:3, i), & 
                          " and radius ", group_parameters_real(4, i)
        case('type')
            write(msg, *) "Group ", trim(group_name(i)), " selects all atoms of type ", group_parameters_int(1,i)

        case('file')
            write(msg, *) "Group ", trim(group_name(i)), " is read in from file ", trim(adjustl(group_parameters_string(1,i)))
        end select
        call log_msg(msg)


        !Assign atom groups
        if(trim(adjustl(group_shape(i))) == "file") then 
            !First make a map for the tags to position in the array
            
            open(unit=11, file=trim(adjustl(group_parameters_string(1,i))), status='old', action='read', position='rewind')
            !First read in number of atoms/elements to read
            read(11,*) reada, reade

            !Now first read the atoms
            if(reada > 0) then 
                !If we are actually interested
                if( (trim(adjustl(which_in_group(i))) == 'atoms').or. (trim(adjustl(which_in_group(i))) == 'all'))then 
                    allocate(amap(maxval(tag_atom)))
                    amap=0
                    !First create the tag to index map
                    do j = 1, atom_num_l
                        amap(tag_atom(j)) = j 
                    end do

                    !Now read in the data and set the group mask
                    do j = 1, reada
                        read(11,*) intag
                        if(intag <= size(amap)) then 
                            if(amap(intag) > 0) then 
                                a_mask(amap(intag)) = ibset(a_mask(amap(intag)),i)
                                group_countme(1) = group_countme(1) + 1
                            end if
                        end if
                    end do
                else
                    !Otherwise just read until you are done with the atom reading
                    do j = 1, reada
                        read(11,*) intag
                    end do
                end if
            end if

            if(reade > 0) then 
                if( (trim(adjustl(which_in_group(i))) == 'elements').or. (trim(adjustl(which_in_group(i))) == 'all'))then 
                    !First create the tag to index map
                    allocate(emap(maxval(tag_ele)))
                    emap =0
                    do j = 1, ele_num_l
                        emap(tag_ele(j)) = j
                    end do

                    !Now read in the data and set the group mask
                    do j=1, reade
                        read(11,*) intag
                        if(intag <= size(emap)) then 
                            if(emap(intag) > 0) then 
                                e_mask(emap(intag)) = ibset(e_mask(emap(intag)), i)
                                ie=emap(intag)
                                if(who_has_ele(ie)) then 
                                    group_countme(2) = group_countme(2) + 1
                                    group_countme(3) = group_countme(3) + ng_node(etype(ie))*basis_num(ie)
                                    group_countme(4) = group_countme(4) + basis_num(ie) * (size_ele(ie)+1)**3
                                end if
                            end if
                        end if
                    end do
                end if
            end if
            close(11)
        else
            if((trim(adjustl(which_in_group(i))) == 'atoms').or. (trim(adjustl(which_in_group(i))) == 'all')) then 
                do ia = 1, atom_num_l
                    if( trim(adjustl(group_shape(i))) == 'type') then 
                        if(type_atom(ia) == group_parameters_int(1,i)) then 
                            a_mask(ia) = ibset(a_mask(ia), i)
                            group_countme(1) = group_countme(1) + 1
                        end if

                    else if (in_group(i, r_atom(:,ia))) then 
                        a_mask(ia) = ibset(a_mask(ia), i)
                        group_countme(1) = group_countme(1) + 1
                    end if
                end do
            end if

            !Assign element groups
            if((trim(adjustl(which_in_group(i))) == 'elements').or. (trim(adjustl(which_in_group(i))) == 'all')) then 
                j=0
                do ie = 1, ele_num_l
                    if( trim(adjustl(group_shape(i))) == 'type') then 
                        if(any(basis_type(:,ie) == group_parameters_int(1,i))) then 
                            e_mask(ie) = ibset(e_mask(ie), i)
                            if(who_has_ele(ie)) then 
                                group_countme(2) = group_countme(2) + 1
                                group_countme(3) = group_countme(3) + ng_node(etype(ie))*basis_num(ie)
                                group_countme(4) = group_countme(4) + basis_num(ie)*(size_ele(ie)+1)**3
                            end if
                        end if

                    else 
                        !Otherwise get the element centroid 
                        rc = 0.0_wp
                        do inod = 1, ng_node(etype(ie))
                            ip = cg_node(inod, ie)
                            do ibasis = 1, basis_num(ie)
                                rc = rc + r(:,ibasis, ip)
                            end do
                        end do
                        rc = rc/(ng_node(etype(ie))*basis_num(ie))

                        !If centroid is in group then we set the mask
                        if (in_group(i, rc))  then 
                            e_mask(ie) =  ibset(e_mask(ie), i)
                            if(who_has_ele(ie)) then 
                                group_countme(2) = group_countme(2) + 1
                                group_countme(3) = group_countme(3) + ng_node(etype(ie))*basis_num(ie)
                                group_countme(4) = group_countme(4) + basis_num(ie)*(size_ele(ie)+1)**3
                            end if
                        end if
                    end if
                end do
            end if
        end if

        call mpi_allreduce(group_countme, group_counts(:,i), 4, mpi_integer, mpi_sum, world, ierr)
        write(msg,*) "Group ", trim(group_name(i)), " with num ", i, " has ", group_counts(1,i), " atoms and ", &
                      group_counts(2,i), " elements"
        call log_msg(msg)


    end subroutine assign_group

    pure function in_group(g, r)
        !Return true if in group
        integer, intent(in) :: g
        real(kind=wp), intent(in) :: r(3)

        real(kind=wp) :: rsq, rdiff(3)
        logical :: in_group

        in_group = .false.

        select case(trim(adjustl(group_shape(g))))
        case('block')
            in_group = in_block_bd(r(:), group_parameters_real(1:6, g))

        case('sphere')
            rdiff = r-group_parameters_real(1:3, g)
            rsq = rdiff(1)*rdiff(1) + rdiff(2)*rdiff(2) + rdiff(3)*rdiff(3)
            if (rsq < group_parameters_real(4, g)) in_group = .true.
        end select

    end function in_group

    subroutine write_group(line)
        !This subroutine writes out the group information to a file named after the group 
        character(len=*),intent(in) :: line
        
        character(len=read_len) :: group_name, txtholder
        integer :: g, at, el, a, e, i, counts_atom(pro_num), displs_atom(pro_num), counts_ele(pro_num), displs_ele(pro_num)
        integer, allocatable :: gatoms(:), gelements(:), gatoms_gather(:), gelements_gather(:)

        read(line, *) txtholder, group_name

        g=get_group_index(group_name) 

        !Now allocate the variables which will contain the tags
        at=group_counts(1,g) 
        el=group_counts(2,g)
        allocate(gatoms(at), gelements(el)) 
        
        !Allocate gather variables on root
        if(rank==root) then 
            allocate(gatoms_gather(at), gelements_gather(el))
        else
            allocate(gatoms_gather(1), gelements_gather(1))
        end if
        
        !Now add the atoms and element tags within the group to the list 
        if(atom_num > 0) then
            a=0
            do i = 1, atom_num_l
                if (btest(a_mask(i), g)) then 
                    a=a+1
                    gatoms(a) = tag_atom(i)
                end if
            end do
        end if
        if(ele_num > 0) then 
            e=0
            do i = 1, ele_num_l
                if(who_has_ele(i)) then 
                    if(btest(e_mask(i), g)) then 
                        e=e+1
                        gelements(e) = tag_ele(i)
                    end if
                end if
            end do
        end if

        !First gather counts and displacements for atoms
        if(atom_num > 0 ) then 
            displs_atom=0
            counts_atom=0
            call mpi_gather(a, 1, mpi_integer, counts_atom, 1, mpi_integer, root, world, ierr) 
            if(rank==root) then 
                displs_atom=0
                do i =2, pro_num
                    displs_atom(i)=displs_atom(i-1)+counts_atom(i-1)
                end do
                if(at /= sum(counts_atom)) call misc_error("Sum of counts doesn't match atoms in group in write_group")
            end if
        end if

        !Now on root gather counts and displacements for elements
        if(ele_num > 0) then 
            counts_ele=0
            displs_ele=0
            call mpi_gather(e, 1, mpi_integer, counts_ele, 1, mpi_integer, root, world, ierr) 
            if(rank==root) then 
                displs_ele=0
                do i =2, pro_num
                    displs_ele(i)=displs_ele(i-1)+counts_ele(i-1)
                end do
                if(el /= sum(counts_ele)) call misc_error("Sum of counts doesn't match eles in group in write_group")
            end if
        end if

        !Now gatherv the tags to root
        call mpi_gatherv(gatoms(1:a), a, mpi_integer, gatoms_gather, counts_atom, displs_atom, mpi_integer, root, world, ierr) 
        call mpi_gatherv(gelements(1:e), e, mpi_integer, gelements_gather, &
                         counts_ele, displs_ele, mpi_integer, root, world, ierr) 
        
        !Now write out the data if you are root 
        if(rank==root) then 
            open(unit=11, file=trim(adjustl(group_name))//"_tags.dat", action="write", status='replace', position='rewind')
            write(11,*) at, el

            if(atom_num > 0) then
                do i=1, at
                    write(11, *) gatoms_gather(i)          
                end do
            end if
            if(ele_num > 0) then 
                do i=1, el
                    write(11, *) gelements_gather(i)          
                end do
            end if
            close(11)
        end if
    end subroutine write_group

end module group
