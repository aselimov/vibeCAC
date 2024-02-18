module cg

    use parameters
    use comms
    use elements 
    use forces  
    use time
    use atom_types
    use neighbors
    use potential
    use integration
    use min_arrays
    
    implicit none

    integer :: neval, nlimit
    real(kind=wp), private, save :: alpha_max, alpha_reduce, backtrack_slope, quadratic_tol, emach, eps_quad, dmax, gg


    public 
    contains

    subroutine cg_defaults
        alpha_max = 1.0_wp
        alpha_reduce = 0.25_wp
        backtrack_slope = 0.4_wp
        quadratic_tol = 0.1_wp
        emach = 1.0d-8
        eps_quad = 1.0d-28
        dmax = 0.1
        nlimit = min(huge(1),atom_num+node_num)
    end subroutine cg_defaults

    subroutine cg_init
        integer :: ia, ibasis, ip, ie, inod
        real(kind = wp) :: ggme 

        !Initialize cg variables
        call alloc_min_arrays
        if(atom_num > 0) then 
            hatom = force_atom(:,1:atom_num_l)
            gatom = hatom
        end if
        if(ele_num > 0) then 
            hnode = force_eq(:,:,1:node_num_l)
            gnode = hnode
        end if
        
        !Initialize number of evaluations to 0
        neval = 0

        !Initialize force norm_squared
        ggme = 0
        do ia = 1, atom_num_l
            ggme = ggme + force_atom(1, ia)*force_atom(1,ia) + force_atom(2,ia)*force_atom(2,ia) &
                        + force_atom(3,ia)*force_atom(3,ia)
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod,ie)
                    do ibasis = 1, basis_num(node_cg(ip))
                        ggme = ggme + force_eq(1, ibasis,ip)*force_eq(1,ibasis,ip) &
                                    + force_eq(2,ibasis,ip)*force_eq(2,ibasis,ip) & 
                                    + force_eq(3,ibasis,ip)*force_eq(3,ibasis,ip)
                    end do
                end do
            end if
        end do
       
        call mpi_allreduce(ggme, gg, 1, mpi_wp, mpi_sum, world, ierr)
    end subroutine cg_init

    subroutine cg_clean
        if(atom_num > 0) deallocate(hatom, gatom, rzeroatom)
        if(ele_num > 0) deallocate(hnode, gnode, rzeronode)
    end subroutine cg_clean

    subroutine reset_cg_dir
        !This simple subroutine just resets the cg search direction
        if(atom_num > 0) then 
            hatom = force_atom(:,1:atom_num_l)
            gatom = hatom
        end if
        if(ele_num > 0) then 
            hnode = force_eq(:,:,1:node_num_l)
            gnode = hnode
        end if
    end subroutine reset_cg_dir

    subroutine resize_min_array
        !This subroutine resizes the min arrays if needed
        real(kind=wp), allocatable :: gatom_array(:,:), hatom_array(:,:), gnode_array(:,:,:), hnode_array(:,:,:)
        if (atom_num > 0) then 
            allocate(gatom_array(3, atom_num_l), hatom_array(3, atom_num_l))
            gatom_array=0.0_wp
            call move_alloc(gatom_array, gatom)
            hatom_array=0.0_wp
            call move_alloc(hatom_array, hatom)

            hatom = force_atom(:,1:atom_num_l)
            gatom = hatom
        end if

        if (ele_num > 0) then 
            allocate(gnode_array(3, max_basisnum, node_num_l), hnode_array(3, max_basisnum, node_num_l))
            gnode_array=0.0_wp
            call move_alloc(gnode_array, gnode)
            hnode_array=0.0_wp
            call move_alloc(hnode_array, hnode)

            hnode = force_eq(:,:,1:node_num_l)
            gnode = hnode
        end if

    end subroutine resize_min_array


    subroutine cg_iterate(code)
        integer, intent(out) :: code

        integer :: ie, inod, ibasis, ip, ia
        real(kind= wp) :: beta, dotme(2), dotall(2), fdotf, eprevious

        !Line minimization along direction h from current position
        eprevious = compute_pe(1)
        call linesearch_backtrack(code)

        if(code > 0) return

        !Calculate force dot products
        dotme(:) = 0
        do ia = 1, atom_num_l
            dotme(1) = dotme(1) + force_atom(1, ia)*force_atom(1,ia) + force_atom(2,ia)*force_atom(2,ia) &
                                + force_atom(3,ia)*force_atom(3,ia)

            dotme(2) = dotme(2) + force_atom(1, ia)*gatom(1,ia) + force_atom(2,ia)*gatom(2,ia) &
                                + force_atom(3,ia)*gatom(3,ia)
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(node_cg(ip))
                        dotme(1) = dotme(1) + force_eq(1, ibasis,ip)*force_eq(1,ibasis,ip) &
                                            + force_eq(2,ibasis,ip)*force_eq(2,ibasis,ip) &
                                            + force_eq(3,ibasis,ip)*force_eq(3,ibasis,ip)

                        dotme(2) = dotme(2) + force_eq(1, ibasis,ip)*gnode(1,ibasis,ip) &
                                            + force_eq(2,ibasis,ip)*gnode(2,ibasis,ip) &
                                            + force_eq(3,ibasis,ip)*gnode(3,ibasis,ip)
                    end do
                end do
            end if
        end do       
        call mpi_allreduce(dotme, dotall, 2, mpi_wp, mpi_sum, world, ierr)
        fdotf = dotall(1)

        ! update new search direction h from new f = -Grad(x) and old g
        ! this is Polak-Ribieri formulation
        ! beta = dotall[0]/gg would be Fletcher-Reeves
        ! reinitialize CG every ndof iterations by setting beta = 0.0
        beta = max(0.0_wp,(dotall(1)-dotall(2))/gg)
        if( mod((iter+1),nlimit) == 0) beta = 0.0_wp
        gg = dotall(1)

        !Update g and h
        do ia = 1, atom_num_l
            gatom(:,ia) = force_atom(:,ia)
            hatom(:,ia) = gatom(:, ia) + beta*hatom(:,ia)
        end do
        do ip = 1, node_num_l
            do ibasis = 1, basis_num(node_cg(ip))
                gnode(:,ibasis,ip) = force_eq(:,ibasis,ip)
                hnode(:,ibasis,ip) = gnode(:, ibasis,ip) + beta*hnode(:,ibasis,ip)
            end do
        end do

        !Reinitialize CG if new search direction is not downhill
        dotme(1) = 0.0_wp
        do ia = 1, atom_num_l
            dotme(1) = dotme(1) + gatom(1,ia)*hatom(1,ia) + gatom(2,ia)*hatom(2,ia) + gatom(3,ia)*hatom(3,ia)
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        dotme(1) = dotme(1) + gnode(1,ibasis,ip)*hnode(1,ibasis,ip) + gnode(2,ibasis,ip)*hnode(2,ibasis,ip) &
                                        + gnode(3,ibasis,ip)*hnode(3,ibasis,ip)
                    end do
                end do
            end if
        end do
        call mpi_allreduce(dotme(1), dotall(1), 1, mpi_wp, mpi_sum, world, ierr)
        if(dotall(1) <= 0.0) then 
            if(atom_num > 0) hatom=gatom
            if(ele_num > 0) hnode=gnode
        end if
        
    end subroutine cg_iterate

    subroutine linesearch_backtrack(code)
        !linemin: backtracking line search (Proc 3.1, p 41 in Nocedal and Wright)
        !uses no gradient info, but should be very robust
        !start at maxdist, backtrack until energy decrease is sufficient
        integer, intent(out) :: code

        integer ::  ia, ip, ibasis, inod, ie
        real(kind=wp) :: fdothall, fdothme, hme, hmaxall, de_ideal, de, ecurrent, eoriginal, alpha

        !First calculate fdothall which is the projection of search dir along downhill gradient
        fdothme = 0.0_wp
        do ia = 1, atom_num_l
            fdothme = fdothme + force_atom(1,ia)*hatom(1,ia) + force_atom(2,ia)*hatom(2,ia) + force_atom(3, ia)*hatom(3,ia)
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        fdothme = fdothme + force_eq(1,ibasis,ip)*hnode(1,ibasis,ip) + force_eq(2,ibasis,ip)*hnode(2,ibasis,ip) &
                                          + force_eq(3,ibasis,ip)*hnode(3,ibasis,ip)
                    end do
                end do
            end if
        end do
        call mpi_allreduce(fdothme, fdothall, 1, mpi_wp, mpi_sum, world, ierr)
        !If search direction isn't downhill return with error
        if(fdothall <= 0.0_wp) then 
            code = 4
            return
        end if

        ! set alpha so no dof is changed by more than max allowed amount
        ! for atom coords, max amount = dmax
        ! for extra per-atom dof, max amount = extra_max[]
        ! for extra global dof, max amount is set by fix
        ! also insure alpha <= ALPHA_MAX
        ! else will have to backtrack from huge value when forces are tiny
        ! if all search dir components are already 0.0, exit with error
        hme = 0.0_wp
        do ia = 1, atom_num_l
            hme = max(hme, hatom(1,ia), hatom(2,ia), hatom(3,ia))
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(node_cg(ip))
                        hme = max(hme, hnode(1,ibasis,ip), hnode(2,ibasis,ip), hnode(3,ibasis,ip))
                    end do
                end do
            end if
        end do       
        call mpi_allreduce(hme, hmaxall, 1, mpi_wp, mpi_max, world, ierr)
        if(is_equal(hmaxall, 0.0_wp)) then 
            code = 5
            return
        end if
        !Set starting alpha
        alpha = min(alpha_max, dmax/hmaxall)

        !Save initial coordinates
        if(atom_num_l > 0) rzeroatom(:,1:atom_num_l) = r_atom(:,1:atom_num_l)
        if(node_num_l > 0) rzeronode(:,:,1:node_num_l) = r(:,:,1:node_num_l)
        eoriginal = compute_pe(1)

        !Backtrack until decrease in energy is sufficient
        do while(.true.)
            !Step forward along h
            if (alpha > 0.0_wp) then 
                do ia = 1, atom_num_l
                    r_atom(:,ia) = rzeroatom(:,ia) + alpha*hatom(:,ia) 
                end do

                do ip = 1, node_num_l
                    do ibasis = 1, basis_num(node_cg(ip))
                        r(:,ibasis,ip) = rzeronode(:,ibasis, ip) + alpha*hnode(:,ibasis, ip)
                    end do
                end do
                neval = neval + 1
                call update_neighbor(iter, .true.)
                call update_force
                ecurrent=compute_pe(1)

            end if

            !Quit with success if energy change is better than ideal
            de_ideal = -backtrack_slope*alpha*fdothall
            de = ecurrent-eoriginal
            if(de <= de_ideal) then 
                code = 0
                return
            end if

            !Reduce alpha
            alpha = alpha * alpha_reduce
            !Backtracked too much, reset to starting point
            ! if de is positive exit with error'
            ! if de is negative than we have reached energy tol
            if((alpha <= lim_zero) .or. (de_ideal >= EMACH)) then 
                if(atom_num_l > 0) r_atom(:, 1:atom_num_l) = rzeroatom(:, 1:atom_num_l)
                if(ele_num_l > 0) r(:,:,1:node_num_l) = rzeronode(:,:,1:node_num_l)
                call update_neighbor(iter, .true.)
                call update_force
                if(de < 0.0) then 
                    code = 1
                    return
                else
                    code = 6 
                    return
                end if
            end if
        end do
    end subroutine linesearch_backtrack



    subroutine linesearch_quadratic(code)
        !----------------------------------------------------------------------
        ! linemin: quadratic line search (adapted from Dennis and Schnabel)
        ! The objective function is approximated by a quadratic
        ! function in alpha, for sufficiently small alpha.
        ! This idea is the same as that used in the well-known secant
        ! method. However, since the change in the objective function
        ! (difference of two finite numbers) is not known as accurately
        ! as the gradient (which is close to zero), all the expressions
        ! are written in terms of gradients. In this way, we can converge
        ! the LAMMPS forces much closer to zero.
        !
        ! We know E,Eprev,fh,fhprev. The Taylor series about alpha_prev
        ! truncated at the quadratic term is:
        !
        !     E = Eprev - del_alpha*fhprev + (1/2)del_alpha^2*Hprev
        !
        ! and
        !
        !     fh = fhprev - del_alpha*Hprev
        !
        ! where del_alpha = alpha-alpha_prev
        !
        ! We solve these two equations for Hprev and E=Esolve, giving:
        !
        !     Esolve = Eprev - del_alpha*(f+fprev)/2
        !
        ! We define relerr to be:
        !
        !      relerr = |(Esolve-E)/Eprev|
        !             = |1.0 - (0.5*del_alpha*(f+fprev)+E)/Eprev|
        !
        ! If this is accurate to within a reasonable tolerance, then
        ! we go ahead and use a secant step to fh = 0:
        !
        !      alpha0 = alpha - (alpha-alphaprev)*fh/delfh;
        !
        !------------------------------------------------------------------------- 
        integer, intent(out) :: code

        integer ::  ia, ip, ibasis, ie, inod
        real(kind=wp) :: fdothall, fdothme, hme, hmaxall, de_ideal, de, ecurrent, eoriginal, alpha, alphaprev, fhprev, &
                         engprev, dotme(2), dotall(2), fh, ff, delfh, relerr, alpha0, alphamax

        !First calculate fdothall which is the projection of search dir along downhill gradient
        fdothme = 0.0_wp
        do ia = 1, atom_num_l
            fdothme = fdothme + force_atom(1,ia)*hatom(1,ia) + force_atom(2,ia)*hatom(2,ia) + force_atom(3, ia)*hatom(3,ia)
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie)) 
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        fdothme = fdothme + force_eq(1,ibasis,ip)*hnode(1,ibasis,ip) &
                                          + force_eq(2,ibasis,ip)*hnode(2,ibasis,ip) &
                                          + force_eq(3,ibasis,ip)*hnode(3,ibasis,ip)
                    end do
                end do
            end if
        end do
        call mpi_allreduce(fdothme, fdothall, 1, mpi_wp, mpi_sum, world, ierr)
        !If search direction isn't downhill return with error
        if(fdothall <= 0.0_wp) then 
            code = 4
            return
        end if


        ! set alpha so no dof is changed by more than max allowed amount
        ! for atom coords, max amount = dmax
        ! for extra per-atom dof, max amount = extra_max[]
        ! for extra global dof, max amount is set by fix
        ! also insure alpha <= ALPHA_MAX
        ! else will have to backtrack from huge value when forces are tiny
        ! if all search dir components are already 0.0, exit with error
        hme = 0.0_wp
        do ia = 1, atom_num_l
            hme = max(hme, hatom(1,ia), hatom(2,ia), hatom(3,ia))
        end do
        do ie = 1, ele_num_l
            if(who_has_ele(ie)) then 
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        hme = max(hme, hnode(1,ibasis,ip), hnode(2,ibasis,ip), hnode(3,ibasis,ip))
                    end do
                end do
            end if
        end do       
        call mpi_allreduce(hme, hmaxall, 1, mpi_wp, mpi_max, world, ierr)
        if(is_equal(hmaxall, 0.0_wp)) then 
            code = 5
            return
        end if
        !Get alphamax
        alphamax = min(alpha_max, dmax/hmaxall)

        !Save initial coordinates
        if(atom_num_l > 0) rzeroatom(:,1:atom_num_l) = r_atom(:,1:atom_num_l)
        if(node_num_l > 0) rzeronode(:,:,1:node_num_l) = r(:,:,1:node_num_l)
        eoriginal = compute_pe(1)
        alpha = alphamax
        fhprev = fdothall
        engprev = eoriginal
        alphaprev = 0.0_wp

        !Backtrack until decrease in energy is sufficient
        do while(.true.)
            
            
            !Step forward along h
            if (alpha > 0.0_wp) then 
                do ia = 1, atom_num_l
                    r_atom(:,ia) = rzeroatom(:,ia) + alpha*hatom(:,ia) 
                end do

                do ip = 1, node_num_l
                    do ibasis = 1, basis_num(node_cg(ip))
                        r(:,ibasis,ip) = rzeronode(:,ibasis, ip) + alpha*hnode(:,ibasis, ip)
                    end do
                end do
                neval = neval + 1
                call update_neighbor(iter, .true.)
                call update_force
                ecurrent=compute_pe(1)
            end if

            ! Compute new fh, alpha, delfh
            dotme(:) = 0
            do ia = 1, atom_num_l
                dotme(1) = dotme(1) + force_atom(1, ia)*force_atom(1,ia) + force_atom(2,ia)*force_atom(2,ia) &
                                    + force_atom(3,ia)*force_atom(3,ia)

                dotme(2) = dotme(2) + force_atom(1, ia)*hatom(1,ia) + force_atom(2,ia)*hatom(2,ia) &
                                    + force_atom(3,ia)*hatom(3,ia)
            end do
            do ie = 1, ele_num_l
                if(who_has_ele(ie)) then 
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis = 1, basis_num(ie)
                            dotme(1) = dotme(1) + force_eq(1, ibasis,ip)*force_eq(1,ibasis,ip) &
                                                + force_eq(2,ibasis,ip)*force_eq(2,ibasis,ip) &
                                                + force_eq(3,ibasis,ip)*force_eq(3,ibasis,ip)

                            dotme(2) = dotme(2) + force_eq(1, ibasis,ip)*hnode(1,ibasis,ip) &
                                                + force_eq(2,ibasis,ip)*hnode(2,ibasis,ip) &
                                                + force_eq(3,ibasis,ip)*hnode(3,ibasis,ip)
                        end do
                    end do
                end if
            end do       
            call mpi_allreduce(dotme, dotall, 2, mpi_wp, mpi_sum, world, ierr)

            ff = dotall(1)
            fh = dotall(2)
            delfh = fh-fhprev

            !If fh or delfh is epsilon, reset to starting point and exit with error
            if((abs(fh) < eps_quad) .or. (abs(delfh) < eps_quad)) then 
                if(atom_num_l > 0) r_atom(:, 1:atom_num_l) = rzeroatom(:, 1:atom_num_l)
                if(node_num_l > 0) r(:,:,1:node_num_l) = rzeronode(:, :, 1:node_num_l)
                call update_neighbor(iter, .true.)
                call update_force

                code = 11
                return
            end if

            !!Check if ready for quadratic projection, quivalent to secant method
            !alpha0 = projected alpha
            relerr = abs(1.0_wp-(0.5_wp*(alpha-alphaprev)*(fh+fhprev)+ecurrent)/engprev)
            alpha0 = alpha - (alpha-alphaprev)*fh/delfh
            if((relerr <= quadratic_tol).and.(alpha0 > 0.0).and.(alpha0<alphamax)) then 

                !Step forward along h
                do ia = 1, atom_num_l
                    r_atom(:,ia) = rzeroatom(:,ia) + alpha0*hatom(:,ia) 
                end do

                do ip = 1, node_num_l
                    do ibasis = 1, basis_num(node_cg(ip))
                        r(:,ibasis,ip) = rzeronode(:,ibasis, ip) + alpha0*hnode(:,ibasis, ip)
                    end do
                end do
                neval = neval + 1
                call update_neighbor(iter, .true.)
                call update_force
                ecurrent=compute_pe(1)
                if((ecurrent-eoriginal) < emach) then 
                    code = 0
                    return 
                end if
            end if

            !Quit with success if energy change is better than ideal
            de_ideal = -backtrack_slope*alpha*fdothall
            de = ecurrent-eoriginal
            if(de <= de_ideal) then 
                code = 0
                return
            end if

            !Save previous state
            fhprev = fh
            engprev = ecurrent
            alphaprev = alpha

            !Reduce alpha
            alpha = alpha * alpha_reduce

            !Backtracked too much, reset to starting point
            ! if de is positive exit with error'
            ! if de is negative than we have reached energy tol
            if((alpha <= 0.0) .or. (de_ideal >= EMACH)) then 
                if(atom_num_l > 0) r_atom(:, 1:atom_num_l) = rzeroatom
                if(ele_num_l > 0) r(:,:,1:node_num_l) = rzeronode

                call update_neighbor(iter, .true.)
                call update_force

                if(de < 0.0) then 
                    code = 1
                    return
                else
                    code = 6 
                    return
                end if
            end if
        end do
    end subroutine linesearch_quadratic

    pure function get_cg_neval()
        integer :: get_cg_neval

        get_cg_neval = neval
        return
    end function get_cg_neval

end module cg
