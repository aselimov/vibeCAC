module fire 

    use parameters
    use comms
    use elements 
    use forces  
    use time
    use atom_types
    use neighbors
    use potential
    use debug
    
    implicit none

    real(kind=wp), private, save :: tmax, tmin, dtgrow, dtshrink, alpha0, alphashrink, alpha, last_negative, dmax, &
                                    dtmax, dtmin
    integer, private, save :: delaystep, vdfmax, vdotf_negatif, neval
    logical, private, save :: initialdelay, flagv0

   real(kind = wp), parameter :: &
     f_inc_fire = 1.1_wp, &
     f_dec_fire = 0.5_wp, &
     f_alpha_fire = 0.99_wp, &
     alpha_start_fire = 0.1_wp

    real(kind=wp), save :: alpha_fire


    public 
    contains

    subroutine fire_defaults
        !This sets the default parameters for fire. These are taken from the default parameters used in lammps
        !Presented in https://doi.org/10.1016%2Fj.commatsci.2020.109584
        !Max timstep is tmax x delta t_start
        tmax = 10.0_wp
        !Minimum timestep is tmin x delta t_start
        tmin = 0.02_wp
        !Number of steps to wait after P<0 before increasing deltat
        delaystep = 20
        !Factor by which del t increases/decrases
        dtgrow = 1.1_wp
        dtshrink = 0.5_wp
        !Coefficient for mixing velocity and force vectors and the factor by which it decreases
        alpha0 = 0.25_wp
        alpha = alpha0
        alphashrink = 0.99_wp
        !Exit after vdfmax consecutive iterations with P(t) < 0
        vdfmax = 2000 !Inertia correction halfstepback = .true.  activates initial delay in modifying delta t and alpha
        initialdelay = .true.
        dmax = 0.1

    end subroutine fire_defaults

    subroutine set_fire_param(param, val)
        character(len=*), intent(in) :: param
        real(kind=wp), intent(in) :: val

        select case(param)
        case('alpha0')
            alpha0=val
        case('tmax')
            tmax=val
        case('tmin')
            tmin=val
        end select

    end subroutine set_fire_param

    subroutine fire_init
        !Initialize fire variables
        last_negative = 0
        vdotf_negatif = 0
        flagv0=.true.
        neval = 0
        alpha = alpha0
        dtmax = tmax*time_step
        dtmin = tmin*time_step
    end subroutine fire_init

    subroutine fire_iterate(i, code)
        integer, intent(in) :: i
        integer, intent(out) :: code

        integer :: ie, ibasis, ip, ia
        real(kind= wp) :: vdotfme, vdotfall, vdotvme, vdotvall, fdotfme, fdotfall, scale1, scale2, dtf, dtfm, &
                          dtvone, vmax, dtv
        logical :: delayflag

        !First calculate vdotfall 
        vdotfme = 0.0_wp
        code=0
        do ip = 1, node_num_l
            ie = node_cg(ip)
            if(who_has_ele(ie)) then 
                do ibasis = 1, basis_num(ie)
                    !We multiply by the energy of all virtual atoms for the energy, which is force_eq times the nodal mass
                    vdotfme = vdotfme + vel(1,ibasis,ip)*force_eq(1,ibasis,ip)+vel(2,ibasis,ip)*force_eq(2,ibasis,ip)  &
                                                 + vel(3,ibasis,ip)*force_eq(3,ibasis,ip)
                end do
            end if
        end do
        do ia = 1, atom_num_l
            vdotfme = vdotfme + force_atom(1,ia)*vel_atom(1,ia) + force_atom(2,ia)*vel_atom(2,ia) &
                              + force_atom(3,ia)*vel_atom(3,ia)
        end do
        call mpi_allreduce(vdotfme, vdotfall, 1, mpi_wp, mpi_sum, world, ierr)


        !if (v dot f) > 0:
        !v = (1-alpha) v + alpha |v| Fhat
        !|v| = length of v, Fhat = unit f
        !Only: (1-alpha) and alpha |v| Fhat is calculated here
        !the modificatin of v is made within the integration, after v update
        !if more than delaystep since v dot f was negative:
        !increase timestep, update global timestep and decrease alpha
        if (vdotfall > 0) then 
            vdotvme = 0.0_wp
            fdotfme = 0.0_wp
            vdotf_negatif = 0
            do ip = 1, node_num_l
                ie = node_cg(ip)
                if (who_has_ele(ie)) then 
                    do ibasis = 1, basis_num(ie)
                        !We multiply by the force/velocity of all virtual atoms for the energy, which is force_eq times the nodal mass
                        vdotvme = vdotvme + vel(1,ibasis,ip)*vel(1,ibasis,ip)+vel(2,ibasis,ip)*vel(2,ibasis,ip) &
                                                     + vel(3,ibasis,ip)*vel(3,ibasis,ip)

                        fdotfme = fdotfme + force_eq(1,ibasis,ip)*force_eq(1,ibasis,ip) &
                                          + force_eq(2,ibasis,ip)*force_eq(2,ibasis,ip) & 
                                          + force_eq(3,ibasis,ip)*force_eq(3,ibasis,ip)
                        
                    end do
                end if
            end do
            do ia = 1, atom_num_l
                vdotvme = vdotvme + vel_atom(1,ia)*vel_atom(1,ia) + vel_atom(2,ia)*vel_atom(2,ia) &
                                  + vel_atom(3,ia)*vel_atom(3,ia)
                fdotfme = fdotfme + force_atom(1,ia)*force_atom(1,ia) + force_atom(2,ia)*force_atom(2,ia) &
                                  + force_atom(3,ia)*force_atom(3,ia)
            end do
            call mpi_allreduce(vdotvme, vdotvall, 1, mpi_wp, mpi_sum, world, ierr)       
            call mpi_allreduce(fdotfme, fdotfall, 1, mpi_wp, mpi_sum, world, ierr)       
            
            !Calculate scaling factors
            scale1 = 1.0 - alpha
            if(fdotfall <= 1d-20) then 
                scale2 = 0.0_wp
            else 
                scale2 = alpha * sqrt(vdotvall/fdotfall)
            end if

            if(i - last_negative > delaystep) then 
                time_step = min(time_step*dtgrow, dtmax)
                alpha = alpha * alphashrink
            end if
        else
            last_negative = i
            delayflag = .true.
            if ((i - 1 < delaystep).and.initialdelay) delayflag = .false.
            if(delayflag) then 
                alpha = alpha0
                if(time_step*dtshrink >= dtmin) time_step = time_step * dtshrink
            end if
            
            !Check stopping criterion
            vdotf_negatif = vdotf_negatif + 1
            if ((vdfmax > 0).and.(vdotf_negatif > vdfmax)) then 
                code = 3
                return
            end if

            !Apply intertia correcection
            do ip = 1, node_num_l
                do ibasis = 1, basis_num(ie)
                    r(:,ibasis, ip) = r(:, ibasis, ip) - 0.5_wp*time_step*vel(:,ibasis,ip)
                end do
            end do
            do ia = 1, atom_num_l
                r_atom(:,ia) = r_atom(:, ia) - 0.5_wp*time_step*vel_atom(:,ia)
            end do

            !Zero velocities
            if(node_num_l > 0) vel(:,:,:) = 0.0_wp
            if(atom_num_l > 0 ) vel_atom(:,:) = 0.0_wp
            flagv0 = .true.
        end if

        !Evaluate velocity to determine whether dtv has to be limited. Required when v is reset
        if(flagv0) then 
            dtf = time_step*ftm2v
            call update_neighbor(iter)
            call update_force
            neval = neval + 1 
            do ip = 1, node_num_l
                ie = node_cg(ip)
                do ibasis = 1, basis_num(ie)
                    dtfm = dtf/masses(basis_type(ibasis,ie))
                    vel(:,ibasis, ip)  = dtfm * force_eq(:, ibasis, ip)
                end do
            end do
            do ia = 1, atom_num_l
                dtfm = dtf/masses(type_atom(ia))
                vel_atom(:,ia) = dtfm*force_atom(:,ia)
            end do           
        end if

        !Limit timestep so no particle moves further than dmax
        dtvone = time_step
        do ip = 1, node_num_l
            if(who_has_ele(node_cg(ip))) then 
                vmax = maxval(vel(:,:,ip))
                if(dtvone*vmax > dmax) dtvone = dmax/vmax
            end if
        end do
        do ia = 1, atom_num_l
            vmax = maxval(vel_atom(:,ia))
            if(dtvone*vmax > dmax) dtvone = dmax/vmax
        end do           
        call mpi_allreduce(dtvone, dtv, 1, mpi_wp, mpi_min, world, ierr)

        !Reset velocities if necessary
        if(flagv0) then 
            if(node_num_l > 0) vel(:,:,:) = 0.0_wp
            if(atom_num_l > 0) vel_atom(:,:) = 0.0_wp
        end if
        
        !Now do the integration step, this is the lammps default semi-implicit Euler scheme

        dtf =  dtv * ftm2v 
         do ip = 1, node_num_l
            do ibasis = 1, basis_num(node_cg(ip))
                dtfm = dtf/masses(basis_type(ibasis,node_cg(ip)))
                vel(:,ibasis, ip)  = vel(:, ibasis, ip) + dtfm * force_eq(:, ibasis, ip)
                if (vdotfall > 0.0_wp) vel(:,ibasis,ip) = scale1*vel(:,ibasis,ip) + scale2*force_eq(:,ibasis,ip)
                r(:,ibasis,ip) = r(:, ibasis, ip) + dtv*vel(:, ibasis,ip)
            end do
        end do

        do ia = 1, atom_num_l
            dtfm = dtf/masses(type_atom(ia))
            vel_atom(:,ia) = vel_atom(:,ia) + dtfm*force_atom(:,ia)
            if(vdotfall > 0.0_wp) vel_atom(:,ia) = scale1*vel_atom(:,ia) + scale2*force_atom(:,ia)
            r_atom(:,ia) = r_atom(:,ia) + dtv*vel_atom(:,ia)
        end do                  

        call update_neighbor(iter)
        call update_force

        neval = neval + 1
!        do ip = 1, node_num_l
!            do ibasis = 1, basis_num(ie)
!                dtfm = dtf/masses(basis_type(ibasis,ie))
!                vel(:,ibasis, ip)  = vel(:, ibasis, ip) + dtfm * force_eq(:, ibasis, ip)
!            end do
!        end do
!        do ia = 1, atom_num_l
!            dtfm = dtf/masses(type_atom(ia))
!            vel_atom(:,ia) = vel_atom(:,ia) + dtfm*force_atom(:,ia)
!        end do                  

        !Set the velocity evaluation flog 
        flagv0 = .false.

    end subroutine fire_iterate

    pure function get_fire_delaystep()
        integer :: get_fire_delaystep

        get_fire_delaystep = delaystep
        return
    end function

    pure function get_fire_neval()
        integer :: get_fire_neval

        get_fire_neval = neval
        return
    end function

    subroutine fire_clean
        time_step = orig_time_step
    end subroutine fire_clean
end module
