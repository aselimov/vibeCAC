module langevin
    !This module contains all the code for running steps of the verlet algorithm
    !This is equivalent to the NVE ensemble
    use parameters
    use comms 
    use elements
    use forces
    use time
    use atom_types
    use group

    implicit none

    !Arguments passed in by users
    real(kind=wp), private, save :: Tb, damp, T0 
    real(kind=wp), private, save :: gfactor1(max_atom_types), gfactor2(max_atom_types)
    integer, save :: langevin_group(20), lnum 
    public
    contains

    subroutine langevin_defaults
        lnum = 0 
        langevin_group = 0
    end subroutine langevin_defaults

    subroutine parse_langevin(line)
        !Parse the langevin thermostat command 
        character(len=*), intent(in) :: line

        character(len=read_len) :: tmptxt, g
        integer :: i

        if(tok_count(line) < 5) call command_error("Missing arguments for dynamics langevin command")
        read(line, *) tmptxt, g, T0, Tb, damp

        !Get the group for this langevin call
        lnum = lnum + 1
        langevin_group(lnum) = get_group_index(g)
        if(langevin_group(lnum) == 0) then 
            write(tmptxt,*) "Group name ", g, " not defined in langevin command"
            call command_error(tmptxt)
        end if

        !Now get force prefactors for all atom types
        do i = 1, natom_types
            gfactor1(i) = -masses(i)/damp/ftm2v
            gfactor2(i) = sqrt(masses(i)) * sqrt(24.0*boltzmann/damp/time_step/const_motion)/ftm2v
        end do

        return
    end subroutine parse_langevin

    subroutine langevin_post_force
        integer :: i, ia, ie, ibasis, inod, ip
        real(kind = wp) :: T_target, delta, tsqrt, fdrag(3), fran(3), rand, gam1, gam2

        delta = (iter - begin_step)/(run_steps)

        T_target = T0 + delta*(Tb-T0)
        tsqrt=sqrt(T_target)

        !update force with langevin for all atoms
        !Loop over all langevin calls
        do i = 1, lnum
            do ia = 1, atom_num_l
                if(btest(a_mask(ia), langevin_group(lnum))) then 
                    gam1 = gfactor1(type_atom(i))
                    gam2 = gfactor2(type_atom(i))*tsqrt

                    call random_number(rand)
                    fran(1) = gam2*(rand-0.5)
                    call random_number(rand)
                    fran(2) = gam2*(rand-0.5)
                    call random_number(rand)
                    fran(3) = gam2*(rand-0.5)

                    fdrag(1) = gam1*vel_atom(1,ia)
                    fdrag(2) = gam1*vel_atom(2,ia)
                    fdrag(3) = gam1*vel_atom(3,ia)

                    force_atom(1,ia) = force_atom(1,ia) + fdrag(1) + fran(1)
                    force_atom(2,ia) = force_atom(2,ia) + fdrag(2) + fran(2)
                    force_atom(3,ia) = force_atom(3,ia) + fdrag(3) + fran(3)
                end if
            end do
            
            do ie = 1, ele_num_l
                if(btest(e_mask(ie), langevin_group(lnum)).and.(who_has_ele(ie))) then
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod,ie)
                        do ibasis = 1, basis_num(ie)
                            gam1 = gfactor1(basis_type(ibasis,ie))
                            gam2 = gfactor2(basis_type(ibasis,ie))*tsqrt

                            call random_number(rand)
                            fran(1) = gam2*(rand-0.5)
                            call random_number(rand)
                            fran(2) = gam2*(rand-0.5)
                            call random_number(rand)
                            fran(3) = gam2*(rand-0.5)

                            fdrag(1) = gam1*vel_atom(1,ia)
                            fdrag(2) = gam1*vel_atom(2,ia)
                            fdrag(3) = gam1*vel_atom(3,ia)

                            force_eq(1,ibasis,ip) = force_eq(1,ibasis,ip) + fdrag(1) + fran(1)
                            force_eq(2,ibasis,ip) = force_eq(2,ibasis,ip) + fdrag(2) + fran(2)
                            force_eq(3,ibasis,ip) = force_eq(3,ibasis,ip) + fdrag(3) + fran(3)                     
                        end do
                    end do
                end if
            end do
        end do

        if(ele_num > 0) call communicate_force_eq

    end subroutine langevin_post_force

end module langevin
