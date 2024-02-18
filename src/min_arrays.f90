module min_arrays
    !This module contains all of the important min_arrays

    use parameters
    use elements
    use errors

    implicit none

    !Arrrays needed for conjugate gradient minimization
    real(kind=wp), allocatable, save :: hatom(:,:), hnode(:,:,:), gatom(:,:), gnode(:,:,:), rzeroatom(:,:), rzeronode(:,:,:)

    public 
    contains

    subroutine alloc_min_arrays
        !Initialize cg variables
        if(ele_num > 0) then 
            if(allocated(rzeronode)) deallocate(rzeronode, hnode, gnode)
            allocate(hnode(3, max_basisnum, node_num_l), gnode(3, max_basisnum, node_num_l), &
                     rzeronode(3, max_basisnum, node_num_l))
        end if 

        if(atom_num > 0) then 
            if(allocated(rzeroatom)) deallocate(rzeroatom, hatom, gatom)
            allocate(hatom(3,atom_num_l), gatom(3,atom_num_l), rzeroatom(3, atom_num_l))
        end if
    
    end subroutine alloc_min_arrays

    subroutine grow_at_min_arrays(seg_real)
        !Grow atom min arrays
        integer, intent(in), optional :: seg_real

        integer :: seg_num
        real(kind=wp), allocatable :: r_array(:,:), g_array(:,:), h_array(:,:)

        
        if(present(seg_real)) then 
            seg_num = seg_real
        else
            seg_num = 1024
        end if

        allocate(r_array(3, size(rzeroatom,2)+seg_num), g_array(3, size(gatom,2)+seg_num), h_array(3, size(hatom,2)+seg_num))

        r_array=0.0_wp
        g_array=0.0_wp
        h_array=0.0_wp

        r_array(:,1:size(rzeroatom,2)) = rzeroatom
        g_array(:,1:size(gatom,2)) = gatom
        h_array(:,1:size(hatom,2)) = hatom

        call move_alloc(r_array, rzeroatom)
        call move_alloc(g_array, gatom)
        call move_alloc(h_array, hatom)

        return
    end subroutine grow_at_min_arrays

    subroutine grow_ele_min_arrays(seg_real)
        !Grow atom min arrays
        integer, intent(in), optional :: seg_real

        integer :: seg_num
        real(kind=wp), allocatable :: r_array(:,:,:), g_array(:,:,:), h_array(:,:,:)

        if(present(seg_real)) then 
            seg_num = seg_real
        else
            seg_num = 1024
        end if

        allocate(r_array(3, max_basisnum, size(rzeronode,3)+seg_num), &
                 g_array(3, max_basisnum, size(gnode,3)+seg_num), h_array(3, max_basisnum, size(hnode,3)+seg_num))

        r_array=0.0_wp
        g_array=0.0_wp
        h_array=0.0_wp

        r_array(:,:,1:size(rzeronode,3)) = rzeronode
        g_array(:,:,1:size(gnode,3)) = gnode
        h_array(:,:,1:size(hnode,3)) = hnode

        call move_alloc(r_array, rzeronode)
        call move_alloc(g_array, gnode)
        call move_alloc(h_array, hnode)

        return
    end subroutine grow_ele_min_arrays



    subroutine pack_atom_cg(r0,g,h, send_cg)
        !This subroutine packs the atom information for conjugate gradient minimization
        real(kind=wp), dimension(3), intent(in) :: r0, g, h
        real(kind=wp), dimension(:), intent(out) :: send_cg
    
        integer :: i, k
    
        !Now pack send_cg
        k = 1 
        do i = 1,3
            send_cg(k) = r0(i)
            k=k+1
            send_cg(k) = g(i)
            k=k+1
            send_cg(k) = h(i)
            k=k+1
        end do
    
        return
    
    end subroutine pack_atom_cg

    subroutine unpack_atom_cg(recv_cg, r0, g, h)
        !This subroutine unpacks the atom information for cg minimization
        real(kind=wp), dimension(:), intent(in) :: recv_cg
        real(kind=wp), dimension(3), intent(out) :: r0, g, h
    
        integer :: i, k
    
        !Now unpack recv_cg
        k = 1 
        do i = 1,3
            r0(i) = recv_cg(k)
            k=k+1
            g(i) = recv_cg(k)
            k=k+1
            h(i) = recv_cg(k)
            k=k+1
        end do
    
        return
    
    end subroutine unpack_atom_cg

    subroutine pack_ele_cg(n, b, r0,g,h, send_cg)
        !This subroutine packs the atom information for conjugate gradient minimization
        integer, intent(in) :: n, b
        real(kind=wp), dimension(3, max_basisnum, ng_max_node), intent(in) :: r0, g, h
        real(kind=wp), dimension(:), intent(out) :: send_cg
    
        integer :: i, k, inod, ibasis
    
        !Now pack send_cg
        k = 1 
        do inod = 1, n
            do ibasis = 1, b
                do i = 1,3
                    send_cg(k) = r0(i,ibasis,inod)
                    k=k+1
                    send_cg(k) = g(i,ibasis,inod)
                    k=k+1
                    send_cg(k) = h(i,ibasis,inod)
                    k=k+1
                end do
            end do
        end do 

        return
    end subroutine pack_ele_cg

    subroutine unpack_ele_cg(n, b, recv_cg, r0, g, h)
        !This subroutine packs the atom information for conjugate gradient minimization
        integer, intent(in) :: n, b
        real(kind=wp), dimension(:), intent(in) :: recv_cg
        real(kind=wp), dimension(3, max_basisnum, ng_max_node), intent(out) :: r0, g, h
    
        integer :: i, k, inod, ibasis
    
        !Now pack send_cg
        k = 1 
        do inod = 1, n
            do ibasis = 1, b
                do i = 1,3
                    r0(i,ibasis,inod)=recv_cg(k)
                    k=k+1
                    g(i,ibasis,inod)=recv_cg(k)
                    k=k+1
                    h(i,ibasis,inod)=recv_cg(k)
                    k=k+1
                end do
            end do
        end do 

        return
    end subroutine unpack_ele_cg
end module min_arrays
