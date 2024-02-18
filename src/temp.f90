module temp 

    use parameters
    use mpi 
    use comms
    use elements 
    use group
    use compute 
    use time
    use atom_types

    implicit none

    integer :: Nf(2), gnum
    real(kind=wp) :: Tinit=0.0_wp !This is the target temp for the velocity initialization subroutine

    !Target temp and kinetic energy for rescaling subroutines, we have three Ktarget values currently for testing
    real(kind=wp) :: Ttarget = 0.0_wp, Ktarget(3), urms(max_atom_types)

    real(kind=wp) :: time_constant !The time constant used to calculate the rescaling coefficient
    real(kind=wp) :: omega !This is the average frequency used for applying an oscillation to the interpolated atoms

    logical :: tflag
    integer, private :: delay

    public 
    contains

    subroutine parse_temp(line)
        !This subroutine parses the temperature command
        character(len = *), intent(in) :: line
        character(len = read_len) :: tmptxt, group_id, command, args(10)
        integer :: iospara,  i, g
        real(kind = wp) :: K(3), T 

        delay = 1.0_wp
        read(line, *, iostat = iospara) tmptxt, command

        !Initialize random number generator
        call random_seed()

        select case(command)
        case('create')
            !Read the temperature command 
            read(line, *, iostat = iospara) tmptxt, command, group_id, T

        case('control')
            omega = 0.0_wp
            !Read the temperature command 
            read(line, *, iostat = iospara) tmptxt, command, group_id, T, time_constant

            !Check for optional arguments 
            read(line, *, iostat = iospara) args
            i = 5
            do while(i < tok_count(line))
                i = i + 1
                select case(args(i))
                case('delay')
                    i = i + 1
                    read(args(i), *) delay
                end select
            end do
            
        end select

        g = get_group_index(group_id)
        if(g == 0) then 
            call misc_error("Missing group name in temp command")
        end if
        !Nf(1) is atomic dof and Nf(2) is elemental dof
        Nf(1) = 3*group_counts(1, g)
        Nf(2) = 3*group_counts(4, g)

        !1st component is target energy for all degrees of freedom
        K(1) = (Nf(1)+Nf(2))*boltzmann*T/2.0_wp
        !2nd is for atomic dof
        K(2) = (Nf(1))*boltzmann*T/2.0_wp
        K(3) = (Nf(2))*boltzmann*T/2.0_wp

        select case(command)
        case('create')
            write(tmptxt, *) "Initializing velocities to match temperature", T, " with kinetic energy",  K(1)
            call log_msg(tmptxt)
            call init_vel(g,K)
        case('control')
            Ktarget = K
            Ttarget = T
            gnum = g

            !if (Nf(2) > 0) then
            !    !Check to see if we have elements, if we do then we have to pass a frequency suboption dictating the frequency of
            !    !oscillation to apply to the elements
            !    read(line, *, iostat = iospara, iomsg = msg) tmptxt, command, group_id, Ttarget, time_constant, tmptxt, om
            !    if(iospara > 0) call read_error(msg, iospara)
            !    !Inputted frequency must be converted to angular frequency, the input frequency should be given in THz
            !    omega = omega*2*pi*10.0_wp**12

            !    !Now calculate the root mean-squared displacement magnitude for each different atom type

            !    do i = 1, natom_types
            !        urms(i) = sqrt((3*10d20*boltzmann*electron_volt*Ttarget)/(masses(i)*amu*omega**2))
            !    end do
            !end if

            write(tmptxt, *) "Running temperature control: Target Ke is ", Ktarget(1), " target temp is ", &
                              Ttarget, " with a time", "constant of ", time_constant
            call log_msg(tmptxt)

            if(omega>lim_small) then 
                write(tmptxt,*) "Applying a mean square displacement of ", (urms(i), i =1, natom_types), " for atom types ", &
                                (i, i =1, natom_types)
                call log_msg(tmptxt)
            end if
            tflag = .true.
        end select

        return
    end subroutine parse_temp

    subroutine init_vel(g, Ktar)
        !This subroutine initializes velocities to a random gaussian distribution around the desired temperature
        real(kind=wp), intent(in) :: Ktar(3)
        integer, intent(in) :: g
        integer :: i, j, ibasis, id, inod, ip, ie
        real(kind=wp) :: alpha, Kme, Kall, vel_buff(3*ele_num*max_basisnum*ng_max_node)
    
        Kme = 0.0_wp
        
        !Allocate velocity arrays if not already allocated
        if(.not. need_vel) call alloc_velocity_arrays

        if(atom_num > 0) then 
            do i=1, atom_num_l
                if (btest(a_mask(i), g)) then 
                    do j = 1, 3
                        vel_atom(j,i) = gasdev()
                        Kme = Kme + 0.5_wp * masses(type_atom(i))*const_motion*(vel_atom(j,i)*vel_atom(j,i))
                    end do
                end if
            end do
        end if

        if(ele_num > 0) then 
            do ie = 1, ele_num_l
                if(who_has_ele(ie).and.btest(e_mask(ie), g)) then 
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis = 1, basis_num(ie) 
                            do j = 1, 3
                                vel(j,ibasis,ip) = gasdev()        
                                Kme = Kme + 0.5_wp * masses(basis_type(ibasis,ie))*const_motion &
                                                   * mass_mat_coeff(size_ele(ie), etype(ie))*(vel(j,ibasis,ip)*vel(j,ibasis,ip)) 
                            end do
                        end do
                    end do
                end if
            end do
        end if

        !Now sum all of the kinetic energies from all processors, calculate the scaling factor and apply it to the velocities
        call mpi_allreduce(Kme, Kall, 1, mpi_wp, mpi_sum, world, ierr)
        alpha=sqrt(Ktar(1)/Kall)
        if(atom_num > 0) vel_atom = alpha*vel_atom
        
        if(ele_num > 0) then 
            vel=alpha*vel
            !Now communicate the element node velocities for shared elements
            vel_buff = 0
            do ie = 1, ele_num_l
                if(who_has_ele(ie)) then 
                    id = ele_glob_id(ie)
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis = 1, basis_num(ie)
                            do j = 1, 3  
                                vel_buff(3*(id-1)*ng_max_node*max_basisnum + 3*(inod-1)*max_basisnum + 3*(ibasis-1)+j) &
                                = vel(j, ibasis, ip)
                            end do
                        end do
                    end do
                end if
            end do

            call mpi_allreduce(mpi_in_place, vel_buff, 3*ng_max_node*ele_num*max_basisnum, mpi_wp, mpi_sum, world, ierr)
            do ie = 1, ele_num_l
                id = ele_glob_id(ie)
                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis = 1, basis_num(ie)
                        do j = 1, 3  
                            vel(j, ibasis, ip) &
                            =vel_buff(3*(id-1)*ng_max_node*max_basisnum + 3*(inod-1)*max_basisnum + 3*(ibasis-1)+j) 
                        end do
                    end do
                end do
            end do
        end if

    end subroutine init_vel


    subroutine rescale_v
    !rescale velocity
    integer :: i, j, ibasis, ie
    real(kind=wp) :: alpha, avg_vel(3), ke(2), rr, rsum, c1, c2, max_vel(3)


    if (mod(iter, delay) == 0) then
        !First calculate the kinetic energy 
        avg_vel = compute_avgvel(gnum)
        ke = compute_ke(gnum, avg_vel)

        if(Nf(1) > 0) then 
            if(rank == root) then 
                !Sample the target kinetic energy from the canonical ensemble
                if (is_equal(time_constant, 0.0_wp)) then 
                    c1=0
                else
                    c1 = exp(-1/(time_constant))
                end if
                c2 = (1.0-c1)*Ktarget(2)/(ke(1))/(Nf(1))
                rr = gasdev()
                rsum = sumnoises(Nf(1)-1)
                !This sampling formulation comes from:
                !Bussi, Giovanni et al. "Canonical sampling through velocity rescaling" The Journal of Chemical Physics (2007)

                alpha = sqrt(c1 + c2*(rr*rr + rsum) + 2.0_wp*rr*sqrt(c1*c2))
            end if

            call mpi_bcast(alpha, 1, mpi_wp, root, world, ierr)
            

            !Now rescale the velocities 

            max_vel=0.0_wp
            do i = 1, atom_num_l
                do j = 1, 3
                    vel_atom(j,i) = avg_vel(j) + alpha*(vel_atom(j,i)-avg_vel(j))
                end do 
                if(norm2(vel_atom(:, i)) > norm2(max_vel)) max_vel = vel_atom(:,i)
            end do

        end if

        if (Nf(2) > 0) then 
            if(rank == root) then 
                !Sample the target kinetic energy from the canonical ensemble
                if (is_equal(time_constant, 0.0_wp)) then 
                    c1=0
                else
                    c1 = exp(-1/(time_constant))
                end if
                c2 = (1.0-c1)*Ktarget(3)/(ke(2))/(Nf(2))
                rr = gasdev()
                rsum = sumnoises(Nf(2)-1)
                !This sampling formulation comes from:
                !Bussi, Giovanni et al. "Canonical sampling through velocity rescaling" The Journal of Chemical Physics (2007)

                alpha = sqrt(c1 + c2*(rr*rr + rsum) + 2.0_wp*rr*sqrt(c1*c2))
            end if

            call mpi_bcast(alpha, 1, mpi_wp, root, world, ierr)
            do i = 1, node_num_l
                ie = node_cg(i)
                do ibasis=1, basis_num(ie)
                    do j = 1, 3
                        vel(j, ibasis, i) = avg_vel(j)+alpha*(vel(j,ibasis,i)-avg_vel(j))
                    end do  
                end do
            end do
        end if
    end if
    return
    end subroutine rescale_v

    

    function sumnoises(nn)
        !Returns the sum of n independent gaussian noises squared. This is equivalent to a chi-squared dist
        implicit none
        integer, intent(in) :: nn
        real(kind=wp) :: sumnoises
        
        if (nn==0) then 
            sumnoises=0.0_wp
        else if (nn==1) then 
            sumnoises=gasdev()**2.0_wp
        else if (modulo(nn,2)==0) then 
            sumnoises=2.0_wp*gamdev(nn/2)
        else
            sumnoises=2.0_wp*gamdev((nn-1)/2) + gasdev()**2.0_wp
        end if
    end function sumnoises

    function gamdev(ia)
        !gamma-distributed random number, implemented as described in numerical recipes

        implicit none
        integer, intent(in) :: ia
        integer :: j
        real(kind=wp) :: am, e, s, v1, v2, x, y, rand
        real(kind=wp) :: gamdev

        if(ia < 1) then 
            print *, "Bad argument to function gamdev"
            call mpi_abort(mpi_comm_world, 1, ierr)
        else if (ia < 6) then 
            x = 1
            do j = 1, ia
                call random_number(rand)
                x = x*rand
            end do

            if( x < lim_small) then
                x = 708.4
            else
                x = -log(x)
            end if
        else
            do while(.true.)
                call random_number(rand)
                v1=rand
                call random_number(rand)
                v2=2.0_wp*rand-1.0_wp
                if((v1*v1+v2*v2)>1.0_wp) cycle
                y = v2/v1
                am=ia-1
                s=sqrt(2.0_wp*am+1.0_wp)
                x=s*y+am
                
                if ( x <= 0) cycle 

                if (((am*log(x/am)-s*y) < -700) .or.(v1<0.00001)) cycle
                e = (1.0_wp + y**2.0_wp)*exp(am*log(x/am)-s*y)
                call random_number(rand)
                if(rand > e) cycle
                exit 
            end do
        end if
        gamdev=x
    end function gamdev

    function gasdev()                                                 
        ! gaussian-distributed random number, implemented as described in numerical recipes
        
        implicit none
        integer, save :: iset=0
        real(kind=wp), save :: gset
        real(kind=wp) :: rand, fac, rsq, v1, v2
        real(kind=wp) gasdev

        if(iset ==0) then 
            do while(.true.)
                call random_number(rand)
                v1=2.0_wp*rand-1.0_wp
                call random_number(rand)
                v2=2.0_wp*rand-1.0_wp
                rsq = v1**2.0_wp + v2**2.0_wp
                if((rsq>=1.0_wp).or.((rsq>-lim_small).and.(rsq < lim_small))) then 
                    cycle
                else 
                    exit
                end if
            end do
            fac=sqrt(-2.0_wp*log(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
        else
            gasdev=gset
            iset=0
        end if
    end function gasdev

    subroutine apply_perturbation
        !This subroutine applies the appropriate perturbation to the interpolated atoms of an element
            
        integer :: ie, iatom, virt_count, iatomap, ibasis
        real(kind = wp) :: bvec(3,3), r1, r2, r3, v(3), max_disp

        max_disp = 0
        do ie = 1, ele_num_l
            select case(etype(ie))
                case(1,2)
                    virt_count = (size_ele(ie)+1)**3
            end select

            !First calculate the lattice basis vectors from the element nodal positions
            bvec(:,1) = r(:,1,cg_node(2,ie)) - r(:,1, cg_node(1,ie))
            bvec(:,2) = r(:,1,cg_node(4,ie)) - r(:,1, cg_node(1,ie))
            bvec(:,3) = r(:,1,cg_node(5,ie)) - r(:,1, cg_node(1,ie))

            !Now normalize the lattice vectors
            bvec(:,1) = bvec(:,1)/norm2(bvec(:,1))
            bvec(:,2) = bvec(:,2)/norm2(bvec(:,2))
            bvec(:,3) = bvec(:,3)/norm2(bvec(:,3))

            !Loop over all interpolated atoms 
            do iatom = 1, virt_count
                do ibasis = 1, basis_num(ie)
                    iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                    if (iatomap /= 0) then 
                        !Get 3 random numbers 
                        r1 = gasdev()
                        r2 = gasdev()
                        r3 = gasdev()

                        !Now scale the lattice vectors by the gaussian number 
                        v = 0.0_wp
                        v = v+r1*bvec(:,1)
                        v = v+r2*bvec(:,2)
                        v = v+r3*bvec(:,3)

                        !Add the vectors together to get the perturbation vector and then scale the magnitude to the mean sqaured
                        !oscillation magnitude
                        
                        v = v*(urms(type_atomap(iatomap))/norm2(v))
                        max_disp = norm2(v)

                        r_atomap(:, iatomap) = r_atomap(:, iatomap) + v
                    end if
                end do
            end do  
        end do 

    end subroutine apply_perturbation 
end module temp
