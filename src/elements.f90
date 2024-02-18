module elements

    use parameters
    use errors
    use box
    use math

    implicit none
    !This module contains all arrays and subroutines for processing atom and element data

    !Variables to track the various global values
    integer, save :: ele_num, node_num, atom_num, ng_max_node, max_basisnum, atomap_num, atomap_max

    logical :: need_vel !this variable denotes if we need velocity arrays

    !Variables to track various local values
    integer, save :: ele_num_l, node_num_l, atom_num_l, atomap_num_l, ele_num_lr, node_num_lr, atom_num_lr, atomap_num_lr, &
                     atomap_num_lg, atom_num_lg

    !Current accepted element type definitions
    integer, parameter :: defined_element_types=3 
    character(len=5), parameter, dimension(defined_element_types) :: element_types = (/ '1NN', '2NN', '20C' /) 
    logical, dimension(defined_element_types) :: etype_present 
    integer, parameter, dimension(defined_element_types) :: ng_node = (/ 8, 8, 20 /) !Number of nodes per element
    integer :: max_size !Max of size_ele

    !Coarse-grained elements
    integer, allocatable, save :: tag_ele(:), size_ele(:), etype(:), basis_num(:), basis_type(:,:), &
                                  cg_node(:,:), node_cg(:), pb_node(:,:,:), cg_atomap(:,:), ele_glob_id(:), e_mask(:)
    real(kind=wp), allocatable, save :: r(:,:,:), vel(:,:,:) !allocate to 3, max_basis_num, node_num_l, so it's nodal not
                                                             ! by element

    !Interpolated atoms
    integer, allocatable, save :: tag_atomap(:), type_atomap(:), atomap_to_intpo(:,:)
    real(kind=wp), allocatable, save :: r_atomap(:,:)
     
    !Interpolation variable
    real(kind=wp), allocatable, save :: a_interpo(:,:,:,:)
    integer :: unique_sizes 
    integer, allocatable :: size_to_shape(:), shape_sizes(:)
    logical, allocatable :: needed_atomap(:)

    !Real atoms
    integer, allocatable :: tag_atom(:), type_atom(:), a_mask(:)
    real(kind=wp), allocatable :: r_atom(:,:), vel_atom(:,:)
    public 
    contains

    subroutine alloc_at_arrays
        !Allocate atom arrays
        integer :: allostat

        if(allocated(tag_atom)) deallocate(tag_atom, type_atom, a_mask, r_atom)

        allocate( tag_atom(atom_num_l), &
                  type_atom(atom_num_l), &
                  a_mask(atom_num_l), &
                  r_atom(3,atom_num_l), &
                  stat=allostat)
        if(allostat > 0) call alloc_error("Failed to allocate at_arrays", allostat)


        tag_atom(:) = 0
        type_atom(:) = 0
        a_mask(:) = 0
        r_atom = 0.0_wp

        if(need_vel) then 
            if(allocated(vel_atom)) deallocate(vel_atom)
            allocate(vel_atom(3,atom_num_l), stat=allostat)
            if(allostat > 0) call alloc_error("Failed to allocate vel_atom", allostat)
            vel_atom = 0.0_wp
        end if
        return

    end subroutine alloc_at_arrays

    subroutine dealloc_at_arrays
        !deallocate atom arrays
        integer :: allostat
        
        deallocate( tag_atom, type_atom, a_mask, r_atom, stat=allostat)
        if(allostat > 0) call alloc_error("Failed to deallocate at_arrays", allostat)

        if(need_vel) then 
            deallocate(vel_atom, stat=allostat)
            if(allostat > 0) call alloc_error("Failed to deallocate vel_atom", allostat)
        end if
        return

    end subroutine dealloc_at_arrays

    subroutine alloc_velocity_arrays
        if((ele_num > 0) .and..not.allocated(vel)) then 
            allocate(vel(3, max_basisnum, node_num_lr), stat = allostat)
            if (allostat > 0) call alloc_error("Failure to allocate vel in alloc_velocity array", allostat)
            vel(:,:,:) = 0.0_wp
        end if
        if((atom_num > 0).and..not.allocated(vel_atom)) then 
            allocate(vel_atom(3, atom_num_lr), stat = allostat)
            if (allostat > 0) call alloc_error("Failure to allocate vel_atom in alloc_velocity array", allostat)
            vel_atom(:,:) = 0.0_wp
        end if
    end subroutine alloc_velocity_arrays
        
    subroutine dealloc_velocity_arrays
        if((ele_num > 0) .and..not.allocated(vel)) then 
            deallocate(vel, stat = allostat)
            if (allostat > 0) call alloc_error("Failure to deallocate vel in alloc_velocity array", allostat)
        end if
        if((atom_num > 0).and..not.allocated(vel_atom)) then 
            deallocate(vel_atom, stat = allostat)
            if (allostat > 0) call alloc_error("Failure to deallocate vel_atom in alloc_velocity array", allostat)
        end if
    end subroutine dealloc_velocity_arrays

    subroutine grow_at_arrays(seg_real)
        !Grow atom arrays
        integer, intent(in), optional :: seg_real

        integer :: grow_by
        integer, allocatable :: tag_array(:), type_array(:), mask_array(:)
        real(kind=wp), allocatable :: r_array(:,:), vel_array(:,:)

        if (present(seg_real)) then
            grow_by = seg_real
        else 
            grow_by = seg_num
        end if

        allocate( tag_array(atom_num_lr+grow_by), &
                  type_array(atom_num_lr+grow_by),& 
                  mask_array(atom_num_lr+grow_by),&
                  r_array(3,atom_num_lr+grow_by), &
                  stat=allostat)  
        if(allostat > 0) call alloc_error("Failure allocating atom arrays in grow_at_arrays", allostat) 

        tag_array(1:atom_num_lr) = tag_atom
        tag_array(atom_num_lr+1:) = 0
        call move_alloc(tag_array, tag_atom)

        type_array(1:atom_num_lr) = type_atom
        type_array(atom_num_lr+1:) = 0
        call move_alloc(type_array, type_atom)

        mask_array(1:atom_num_lr) = a_mask
        mask_array(atom_num_lr+1:) = 0
        call move_alloc(mask_array, a_mask)

        r_array(:,1:atom_num_lr) = r_atom
        r_array(:,atom_num_lr+1:) = 0.0_wp
        call move_alloc(r_array, r_atom)

        if(need_vel) then 
            allocate(vel_array(3,atom_num_lr+grow_by), stat=allostat)
            if(allostat > 0) call alloc_error("Failed allocating vel_array in grow_at_arrays", allostat)
            vel_array(:,1:atom_num_lr) = vel_atom
            vel_array(:,atom_num_lr+1:) = 0.0_wp
            call move_alloc(vel_array, vel_atom)
        end if
        
        atom_num_lr = atom_num_lr +grow_by
        return
    end subroutine grow_at_arrays

    subroutine alloc_cg_arrays
        !Allocate arrays needed for coarse-grained elements 
        integer :: allostat
        
        if(allocated(tag_ele)) deallocate(tag_ele, size_ele, ele_glob_id, etype, cg_node, cg_atomap, basis_num, basis_type, &
                                          node_cg, r)  
        allocate( &
            tag_ele(ele_num_l), &
            size_ele(ele_num_l), &
            ele_glob_id(ele_num_l), &
            etype(ele_num_l), &
            e_mask(ele_num_l), &
            cg_node(ng_max_node, ele_num_l), &
            cg_atomap(atomap_max, ele_num_l),   &
            basis_num(ele_num_l), &
            basis_type(max_basisnum, ele_num_l), &
            node_cg(node_num_l), &
            r(3,max_basisnum, node_num_l), &
            stat=allostat )
        if (allostat > 0) call alloc_error("Failed allocating cg_arrays", allostat)
        tag_ele(:) = 0
        size_ele(:) = 0
        etype(:) = 0
        e_mask(:) = 0
        ele_glob_id(:) = 0
        cg_node(:,:) =0
        cg_atomap(:,:) = 0
        basis_num(:) = 0
        basis_type(:,:) = 0
        node_cg(:) = 0
        r(:,:,:) = 0.0_wp


        if(periodic) then 
            if (allocated(pb_node)) deallocate(pb_node)
            allocate(pb_node(3, max_basisnum, node_num_l), stat=allostat) 
            if (allostat > 0) call alloc_error("Failed allcoating pb_node in alloc_cg_arrays", allostat)
            pb_node(:,:,:) = 0
        end if

        if(need_vel) then 
            if (allocated(vel)) deallocate(vel)
            allocate(vel(3, max_basisnum, node_num_l))
            if (allostat > 0) call alloc_error("Failed allocating cg vel array", allostat)
            vel(:,:,:) = 0.0_wp
        end if

        !Interpolation arrays
        if(allocated(tag_atomap)) deallocate(tag_atomap, type_atomap, atomap_to_intpo, r_atomap)
        allocate (&
            tag_atomap(atomap_num_l),           &
            type_atomap(atomap_num_l),          &
            atomap_to_intpo(2, atomap_num_l),   &
            r_atomap(3, atomap_num_l),          &
            stat=allostat &
        )

        if (allostat > 0) call alloc_error("Failed allocating interp_arrays", allostat) 
        tag_atomap(:) = 0
        type_atomap = 0
        atomap_to_intpo = 0
        r_atomap(:,:) = 0.0_wp

    end subroutine alloc_cg_arrays

    subroutine dealloc_cg_arrays
        !Deallocate arrays needed for coarse-grained elements 
        integer :: allostat
        deallocate( tag_ele, size_ele, ele_glob_id, etype, e_mask, cg_node, cg_atomap, basis_num, basis_type, &
                  node_cg, r, stat=allostat )
        if (allostat > 0) call alloc_error("Failed deallocating cg_arrays", allostat)
        if(periodic) then 
            deallocate(pb_node, stat=allostat) 
            if (allostat > 0) call alloc_error("Failed deallcoating pb_node in alloc_cg_arrays", allostat)
        end if

        if(need_vel) then 
            deallocate(vel)
            if (allostat > 0) call alloc_error("Failed deallocating cg vel array", allostat)
        end if

        !Interpolation arrays
        deallocate(tag_atomap, type_atomap, r_atomap, atomap_to_intpo, stat=allostat)

        if (allostat > 0) call alloc_error("Failed allocating interp_arrays", allostat) 

        return
    end subroutine dealloc_cg_arrays


    subroutine grow_cg_arrays(opt, seg_real)
        !This subroutine just grows the element arrays if needed

        !opt is 1 for growing element arrays, 2 for growing node arrays, 3 for growing atomap arrays
        !seg_real is the growth size, usually compared to 
        integer, intent(in) :: opt
        integer, intent(in), optional :: seg_real

        !temp element arrays
        integer, allocatable ::tag_ele_array(:), size_ele_array(:), etype_array(:), cg_node_array(:,:), basis_num_array(:), &
        basis_type_array(:,:), cg_atomap_array(:,:), mask_array(:)

        !Temp node arrays
        integer, allocatable :: node_cg_array(:),  pb_node_array(:,:,:), ele_glob_id_array(:)
        real(kind=wp), allocatable :: r_array(:,:,:), vel_array(:,:,:)

        !temp atomap arrays 
        integer, allocatable :: tag_atomap_array(:), type_atomap_array(:), atomap_to_intpo_array(:,:)
        real(kind=wp), allocatable :: r_atomap_array(:,:)

        integer :: grow_by
                                 
        if (present(seg_real)) then
            grow_by = seg_real
        else 
            grow_by = seg_num
        end if

        if (opt == 1) then 
            allocate( &
               tag_ele_array(ele_num_lr+grow_by), &
               ele_glob_id_array(ele_num_lr+grow_by), &
               size_ele_array(ele_num_lr+grow_by), &
               etype_array(ele_num_lr+grow_by), &
               mask_array(ele_num_lr+grow_by), &
               cg_node_array(ng_max_node, ele_num_lr+grow_by), &
               cg_atomap_array(atomap_max, ele_num_lr+grow_by), &
               basis_num_array(ele_num_lr+grow_by), &
               basis_type_array(max_basisnum, ele_num_lr+grow_by), &
               stat=allostat &
            )
            if(allostat > 0) call alloc_error("Failed allocation of temp ele arrays in grow_cg_arrays", allostat)

            tag_ele_array(1:ele_num_lr) = tag_ele(:)
            tag_ele_array(ele_num_lr+1:) = 0
            call move_alloc(tag_ele_array, tag_ele)

            ele_glob_id_array(1:ele_num_lr) = ele_glob_id(:)
            ele_glob_id_array(ele_num_lr+1:) = 0
            call move_alloc(ele_glob_id_array, ele_glob_id)

            size_ele_array(1:ele_num_lr) = size_ele(:)
            size_ele_array(ele_num_lr+1:) = 0
            call move_alloc(size_ele_array, size_ele)

            etype_array(1:ele_num_lr) = etype(:)
            etype_array(ele_num_lr+1:) = 0
            call move_alloc(etype_array, etype)

            mask_array(1:ele_num_lr) = e_mask
            mask_array(ele_num_lr+1:)= 0
            call move_alloc(mask_array, e_mask)

            cg_node_array(:,1:ele_num_lr) = cg_node(:,:)
            cg_node_array(:,ele_num_lr+1:) = 0
            call move_alloc(cg_node_array, cg_node)

            cg_atomap_array(:,1:ele_num_lr) = cg_atomap(:,:)
            cg_atomap_array(:,ele_num_lr+1:) = 0
            call move_alloc(cg_atomap_array, cg_atomap)

            basis_num_array(1:ele_num_lr) = basis_num(:)
            basis_num_array(ele_num_lr+1:) = 0
            call move_alloc(basis_num_array, basis_num)

            basis_type_array(:,1:ele_num_lr) = basis_type(:,:)
            basis_type_array(:,ele_num_lr+1:) = 0
            call move_alloc(basis_type_array, basis_type)

            ele_num_lr = ele_num_lr + grow_by
        else if (opt==2) then 
            allocate( node_cg_array(node_num_lr + grow_by),    &
                     r_array(3, max_basisnum, node_num_lr + grow_by), &
                     stat=allostat)

            node_cg_array(1:node_num_lr) = node_cg
            node_cg_array(node_num_lr+1:) = 0
            call move_alloc(node_cg_array, node_cg)

            r_array(:,:,1:node_num_lr) = r(:,:,1:node_num_lr)
            r_array(:,:,node_num_lr+1:) = 0.0_wp
            call move_alloc(r_array, r)

            if(need_vel) then 
                allocate( vel_array(3, max_basisnum, node_num_lr + grow_by), &
                          stat=allostat)

                if(allostat > 0) call alloc_error("Failure allocating vel in grow_cg_arrays",allostat)
                vel_array(:,:,1:node_num_lr) = vel(:,:,1:node_num_lr)
                vel_array(:,:,node_num_lr+1:) = 0.0_wp
                call move_alloc(vel_array, vel)
            end if

            if(periodic) then 
                allocate( pb_node_array(3, max_basisnum, node_num_lr + grow_by), &
                          stat=allostat)
                if(allostat > 0) call alloc_error("Failure allocating pb_node in grow_cg_arrays",allostat)
                pb_node_array(:,:,1:node_num_lr) = pb_node
                pb_node_array(:,:,node_num_lr+1:) = 0.0_wp
                call move_alloc(pb_node_array, pb_node)
            end if

            node_num_lr = node_num_lr + grow_by

        else if(opt==3) then 
            !Resize atomap arrays
            allocate(tag_atomap_array(atomap_num_lr + grow_by), &
                     type_atomap_array(atomap_num_lr + grow_by), &
                     atomap_to_intpo_array(2, size(atomap_to_intpo,2)+grow_by), &
                     r_atomap_array(3,atomap_num_lr + grow_by), &
                     stat = allostat)
            if(allostat > 0) call alloc_error("Failure allocating temp atomap arrays in grow_cg_arrays", allostat)

            tag_atomap_array(1:atomap_num_lr) = tag_atomap(:)
            tag_atomap_array(atomap_num_lr+1:) = 0
            call move_alloc(tag_atomap_array, tag_atomap)

            type_atomap_array(1:atomap_num_lr) = type_atomap(:)
            type_atomap_array(atomap_num_lr+1:) = 0
            call move_alloc(type_atomap_array, type_atomap)

            atomap_to_intpo_array = 0
            atomap_to_intpo_array(:, 1:size(atomap_to_intpo,2)) = atomap_to_intpo
            call move_alloc(atomap_to_intpo_array, atomap_to_intpo)

            r_atomap_array(:,1:atomap_num_lr) = r_atomap(:,:)
            r_atomap_array(:,atomap_num_lr+1:) = 0.0_wp
            call move_alloc(r_atomap_array, r_atomap)

            atomap_num_lr = atomap_num_lr + grow_by
        end if

        return
    end subroutine grow_cg_arrays

    subroutine setup_interpolation(max_esize, unique_enum, unique_esize)
        !This subroutine sets up the interpolation arrays
        integer, intent(in) :: max_esize, unique_enum
        integer, dimension(unique_enum), intent(in) :: unique_esize

        integer :: i, ix, iy, iz, esize, iatom
        real(kind = wp) :: dr, ds, dt

        !First allocate important arrays
        allocate( a_interpo(maxval(ng_node), atomap_max, unique_enum, defined_element_types), &
                  size_to_shape(max_esize),  &
                  shape_sizes(max_esize), &
                  stat = allostat)
        if (allostat > 0) call alloc_error("Failed to allocate shape arrays in setup_interpolation", allostat)

        !Save the maximum esize as we need it later
        max_size = max_esize

        !Now loop over all unique elements. We duplicate the interpolation shape function for all element types, this leads to some
        !extra memory usage, but it shoudl be small enough to where it doesn't matter
        size_to_shape(:) = 0
        unique_sizes = unique_enum
        a_interpo(:,:,:,:) = 0.0_wp
        do i = 1, unique_enum
            esize = unique_esize(i)
            shape_sizes(i) =unique_esize(i)
            !size_to_shape maps the element size to the shape function index
            size_to_shape(esize) = i
            iatom=0
            
            !Cube element shape functions
            do iz = 1, esize+1
                dt=-1.0_dp +(iz-1)*(2.0_dp/(esize))   
                do iy = 1, esize+1
                    ds=-1.0_dp +(iy-1)*(2.0_dp/(esize))   
                    do ix = 1, esize+1
                        dr=-1.0_dp +(ix-1)*(2.0_dp/(esize))   
                        iatom = iatom + 1
                        a_interpo(1,iatom,i,1:2) = (1.0_wp-dr)*(1.0_wp-ds)*(1.0_wp-dt)/8.0_wp
                        a_interpo(2,iatom,i,1:2) = (1.0_wp+dr)*(1.0_wp-ds)*(1.0_wp-dt)/8.0_wp
                        a_interpo(3,iatom,i,1:2) = (1.0_wp+dr)*(1.0_wp+ds)*(1.0_wp-dt)/8.0_wp
                        a_interpo(4,iatom,i,1:2) = (1.0_wp-dr)*(1.0_wp+ds)*(1.0_wp-dt)/8.0_wp
                        a_interpo(5,iatom,i,1:2) = (1.0_wp-dr)*(1.0_wp-ds)*(1.0_wp+dt)/8.0_wp
                        a_interpo(6,iatom,i,1:2) = (1.0_wp+dr)*(1.0_wp-ds)*(1.0_wp+dt)/8.0_wp
                        a_interpo(7,iatom,i,1:2) = (1.0_wp+dr)*(1.0_wp+ds)*(1.0_wp+dt)/8.0_wp
                        a_interpo(8,iatom,i,1:2) = (1.0_wp-dr)*(1.0_wp+ds)*(1.0_wp+dt)/8.0_wp
                    end do
                end do
            end do

            !serendipity 20 node element shape function
            iatom = 0
            do iz = 1, esize+1
                dt = (iz - ((real(esize,wp)/2.0_wp) + 1.0_wp)) / (real(esize,wp) / 2.0_wp)
                do iy = 1, esize+1
                    ds = (iy - ((real(esize,wp)/2.0_wp) + 1.0_wp)) / (real(esize,wp) / 2.0_wp)
                    do ix = 1, esize+1
                        dr = (ix - ((real(esize,wp) / 2.0_wp) + 1.0_wp)) / (real(esize,wp) / 2.0_wp)
                        iatom = iatom + 1
                        
                        !Corner nodes
                        a_interpo(1,iatom,i,3) = (1.0_dp-dr)*(1.0_dp-ds)*(1.0_dp-dt)*(-dr-ds-dt-2)/8.0_dp
                        a_interpo(2, iatom,i,3) = (1.0_dp+dr)*(1.0_dp-ds)*(1.0_dp-dt)*(dr-ds-dt-2)/8.0_dp
                        a_interpo(3, iatom,i,3) = (1.0_dp+dr)*(1.0_dp+ds)*(1.0_dp-dt)*(dr+ds-dt-2)/8.0_dp
                        a_interpo(4, iatom,i,3) = (1.0_dp-dr)*(1.0_dp+ds)*(1.0_dp-dt)*(-dr+ds-dt-2)/8.0_dp
                        a_interpo(5, iatom,i,3) = (1.0_dp-dr)*(1.0_dp-ds)*(1.0_dp+dt)*(-dr-ds+dt-2)/8.0_dp
                        a_interpo(6, iatom,i,3) = (1.0_dp+dr)*(1.0_dp-ds)*(1.0_dp+dt)*(dr-ds+dt-2)/8.0_dp
                        a_interpo(7, iatom,i,3) = (1.0_dp+dr)*(1.0_dp+ds)*(1.0_dp+dt)*(dr+ds+dt-2)/8.0_dp
                        a_interpo(8, iatom,i,3) = (1.0_dp-dr)*(1.0_dp+ds)*(1.0_dp+dt)*(-dr+ds+dt-2)/8.0_dp

                        !Side nodes, first node r is zero
                        a_interpo(9, iatom,i,3) =  (1-dr*dr)*(1-ds)*(1-dt)/4.0_dp
                        a_interpo(11, iatom,i,3) = (1-dr*dr)*(1+ds)*(1-dt)/4.0_dp
                        a_interpo(17, iatom,i,3) = (1-dr*dr)*(1-ds)*(1+dt)/4.0_dp
                        a_interpo(19, iatom,i,3) = (1-dr*dr)*(1+ds)*(1+dt)/4.0_dp

                        !node s is zero
                        a_interpo(10, iatom,i,3) = (1+dr)*(1-ds*ds)*(1-dt)/4.0_dp
                        a_interpo(12, iatom,i,3) = (1-dr)*(1-ds*ds)*(1-dt)/4.0_dp
                        a_interpo(18, iatom,i,3) = (1+dr)*(1-ds*ds)*(1+dt)/4.0_dp
                        a_interpo(20, iatom,i,3) = (1-dr)*(1-ds*ds)*(1+dt)/4.0_dp

                        !node t is zero
                        a_interpo(13, iatom,i,3) = (1-dr)*(1-ds)*(1-dt*dt)/4.0_dp
                        a_interpo(14, iatom,i,3) = (1+dr)*(1-ds)*(1-dt*dt)/4.0_dp
                        a_interpo(15, iatom,i,3) = (1+dr)*(1+ds)*(1-dt*dt)/4.0_dp
                        a_interpo(16, iatom,i,3) = (1-dr)*(1+ds)*(1-dt*dt)/4.0_dp

                    end do
                end do
            end do           
        end do 

        return
    end subroutine setup_interpolation
    
    subroutine dealloc_shape_arrays
        !this subroutine deallocates shape arrays. Will be needed if element sizes change for some reason.
        if(allocated(a_interpo)) then 
            deallocate( a_interpo, size_to_shape, shape_sizes, stat = allostat)
            if (allostat > 0) call alloc_error("Failed deallocating shape arrays in dealloc_shape_arrays", allostat)
        end if

    end subroutine

    subroutine interp_atom(iatom, esize, ele_type, pb_in, basisnum, r_nodes, r_iatom)
        !This subroutine interpolates atom positions based on the iatom number
        
        !ARGUMENTS
        integer, intent(in) :: iatom !Atom to interpolate
        integer, intent(in) :: esize !Size of element, needed for getting the correct interpolation array
        integer, intent(in) :: ele_type
        integer, dimension(3, max_basisnum, ng_max_node), intent(in) :: pb_in !Periodic boundary info
        integer, intent(in) :: basisnum !Number of basis atoms at nodes
        real(kind=wp), dimension(3, max_basisnum, ng_max_node), intent(in) :: r_nodes !Contiguous array of nodal positions
        real(kind=wp), dimension(3,max_basisnum), intent(out) :: r_iatom !Position of all interpolated atoms at lat point
        
        !INTERNAL VARIABLES 
        integer :: inod, ibasis, info(3), node_num
        real(kind=wp) :: rout(3), shape_fun(ng_node(ele_type))
        logical :: ipb, ip_change(max_basisnum,ng_node(ele_type))

        node_num = ng_node(ele_type)
        !Initialization
        r_iatom(:,:) = 0.0_wp    
        ipb =.false.
        ip_change(:,:) = .false.
        shape_fun = a_interpo(1:node_num, iatom, size_to_shape(esize), ele_type)

        do inod = 1, node_num
            do ibasis = 1, basisnum
                !Figure out periodic boundaries for the nodes
                rout=r_nodes(:,ibasis, inod)
                if(periodic) then 
                    !Restore periodic boundaries if needed
                    call restore_pb(rout, pb_in(:,ibasis,inod), ipb)
                    ip_change(ibasis,inod) = ipb
                end if

                r_iatom(:,ibasis) = r_iatom(:,ibasis) + rout*shape_fun(inod)
            end do 
        end do

        !Now apply periodic boundaries if needed
        
        if(any(ip_change)) then 
            do ibasis = 1, basisnum
                call cross_pb(r_iatom(:,ibasis), info)
                if(.not. in_box_bd(r_iatom(:, ibasis))) then 
                    call bounds_error("Atom incorrectly interpolated", r_iatom, box_bd)
                end if
            end do
        end if

        return
    end subroutine interp_atom



    subroutine update_virtual_atoms(update_all)
        !This subroutine updates the virtual atoms from nodal positions
        logical, intent(in) :: update_all

        integer :: ie, virt_count, iatom, ibasis, jatomap, iatomap, inod, ip, pb_in(3, max_basisnum, ng_max_node)
        real(kind = wp) :: r_interp(3, max_basisnum), r_nodes(3, max_basisnum, ng_max_node)

        !Loop over all atoms
        jatomap = 0
        do ie = 1, ele_num_l
            select case(etype(ie))
            case(1,2,3)
                virt_count = (size_ele(ie)+1)**3    
            end select

            !Get the nodal values into a connected array
            r_nodes(:,:,:) = 0.0_wp
            pb_in(:,:,:) = 0
            do inod = 1, ng_node(etype(ie))
                ip = cg_node(inod, ie)
                r_nodes(:,:,inod) = r(:,:,ip)
                if(periodic) pb_in(:,:,inod) = pb_node(:,:,ip)
            end do
            
            do iatom = 1, virt_count

                !interpolate the atoms
                call interp_atom(iatom, size_ele(ie), etype(ie), pb_in, basis_num(ie), r_nodes, r_interp)
                do ibasis = 1, basis_num(ie)
                    !Now get the real atom position
                    iatomap=cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie) 

                    if(iatomap /= 0) then 
                        jatomap = jatomap+1
                        !This just checks to make sure that the cg_atomap is sorted, ie the first virtual atom we have 
                        ! is first in r_atomap and the second is second etc.
                        if(iatomap /= jatomap) then
                            print *, "Error: Iatomap", iatomap, " should equal", jatomap
                            call mpi_abort(mpi_comm_world, 1, ierr)
                        end if
                        
                        r_atomap(:,iatomap) = r_interp(:,ibasis)
                    end if
                end do
            end do
        end do
        return
    end subroutine update_virtual_atoms

    pure function get_virtual_count(et, esize)
        !Return the number of virtual lattice points 
        integer, intent(in) :: et, esize
        integer :: get_virtual_count

        select case(et)
        case(1,2,3)
            get_virtual_count = (esize+1)**3
        end select
        return
    end function get_virtual_count

    subroutine x2frac
        !This subroutine changes r arrays to be fractional coordinates

        integer :: i, j, ibasis

        do i = 1, atom_num_l
            do j = 1,3 
                r_atom(j, i) = (r_atom(j,i)-box_bd(2*j-1))/(box_length(j))
            end do
        end do

        do i = 1, node_num_l
            do ibasis= 1, basis_num(node_cg(i))
                do j = 1, 3
                    r(j, ibasis,  i) = (r(j,ibasis,i)-box_bd(2*j-1))/(box_length(j))
                end do
            end do
        end do
        return
    end subroutine x2frac

    subroutine frac2x
        !This subroutine changes the fractional coordinates back to real coordinates
        integer :: i, j, ibasis

        do i = 1, atom_num_l
            do j = 1,3 
                r_atom(j,i) =  r_atom(j,i)*box_length(j) + box_bd(2*j-1)
            end do
        end do

        do i = 1, node_num_l
            do ibasis= 1, basis_num(node_cg(i))
                do j = 1, 3
                    r(j,ibasis,i) =  r(j, ibasis, i)*box_length(j) + box_bd(2*j-1)
                end do
            end do
        end do
        return
    end subroutine frac2x

end module elements
