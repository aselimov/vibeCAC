module compute

    !This module contains all of the compute subroutines
    use mpi
    use parameters
    use comms
    use atom_types
    use forces
    use elements
    use group
    use box

    implicit none 

    public 
    contains

    function compute_pe(group, arg1)
        !This function returns the summed energy for all atoms and elements, for elements it adds the contribution for all 
        !interpolated atoms
        integer, intent(in) :: group
        logical, intent(in), optional :: arg1 
        real(kind=wp) :: compute_pe
        integer :: ia, ie, ip, ibasis, inod
        real(kind=wp) :: pe, peall
        logical :: weigh_ele
       
        if(present(arg1)) then 
            weigh_ele=arg1
        else 
            weigh_ele=.true.
        end if

        pe = 0
        do ia = 1, atom_num_l
            if (btest(a_mask(ia), group)) then 
                pe = pe + energy_atom(ia)
            end if
        end do 

        if(weigh_ele) then 
            do ie = 1, ele_num_l
                if(btest(e_mask(ie), group).and.who_has_ele(ie)) then 
                    do inod =1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis=1, basis_num(ie)
                            pe = pe + energy_eq(ibasis, ip)*mass_mat_coeff(size_ele(ie), etype(ie))
                        end do
                    end do
                end if
            end do
        else
            do ie = 1, ele_num_l
                if(btest(e_mask(ie), group).and.who_has_ele(ie)) then 
                    do inod =1, ng_node(etype(ie))
                        ip = cg_node(inod, ie)
                        do ibasis=1, basis_num(ie)
                            pe = pe + energy_eq(ibasis, ip)
                        end do
                    end do
                end if
            end do
        end if
        
        call mpi_allreduce(pe, peall, 1, mpi_wp, mpi_sum, world, ierr)
        compute_pe = peall
        return
    end function compute_pe
   
    function compute_fnorm(group)
        !This function returns the global norm of the force vector for all atoms and elements, 
        !for elements it adds the contribution for all interpolated atoms
        integer, intent(in) :: group
        real(kind=wp) :: compute_fnorm
        integer :: ia, ie, ip, ibasis, inod
        real(kind=wp) ::  fme, f(3), fall
        
        fme = 0.0_wp
        do ia = 1, atom_num_l
            if (btest(a_mask(ia), group)) then 
                fme =  fme + force_atom(1,ia)*force_atom(1,ia) + force_atom(2, ia)*force_atom(2,ia) &
                      +force_atom(3,ia)*force_atom(3,ia)                                  
            end if
        end do 

        do ie = 1, ele_num_l
            if(btest(e_mask(ie), group).and.who_has_ele(ie)) then 
                do inod =1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis=1, basis_num(ie)
                        f = force_eq(:,ibasis, ip)
                        fme = fme + (f(1)*f(1) + f(2)*f(2) + f(3)*f(3))
                    end do
                end do
            end if
        end do
        
        call mpi_allreduce(fme, fall, 1, mpi_wp, mpi_sum, world, ierr)
        compute_fnorm = sqrt(fall)
        return
    end function compute_fnorm

    function compute_avgvel(group)
        !This function returns the average velocity of the system
        integer, intent(in) :: group
        real(kind=wp) :: compute_avgvel(3)
        integer :: ia, ie, ip, ibasis, inod
        real(kind=wp) :: vme(3), vall(3)

        vme=0.0_wp
        do ia =1, atom_num_l
            if (btest(a_mask(ia), group)) then 
                vme = vme + vel_atom(:,ia)
            end if
        end do

        do ie = 1, ele_num_l
            if(btest(e_mask(ie), group).and.who_has_ele(ie)) then 
                do inod =1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis=1, basis_num(ie)
                        vme = vme + vel(:, ibasis,ip)
                    end do
                end do
            end if
        end do

        call mpi_allreduce(vme(:), vall(:), 3, mpi_wp, mpi_sum, world, ierr)

        compute_avgvel = vall/(group_counts(1,group)+group_counts(3,group))

        return
    end function compute_avgvel

    function compute_ke(group, in_avg_vel)
        !This subroutine calculates the kinetic energy of the group
        integer, intent(in) :: group
        real(kind=wp), intent(in), optional  :: in_avg_vel(3)

        real(kind=wp) :: compute_ke(2)
        integer :: ia, j, ie, ip, ibasis, inod
        real(kind=wp) :: kme(2), avgvel(3), kall(2)


        if(present(in_avg_vel)) then 
            avgvel = in_avg_vel
        else
            avgvel = compute_avgvel(group)
        end if
        
        kme = 0.0_wp
        do ia = 1, atom_num_l
            if (btest(a_mask(ia), group)) then 
                do j = 1,3
                    kme(1) = kme(1) + 0.5_wp*masses(type_atom(ia))*const_motion*(vel_atom(j,ia) - avgvel(j))**2.0_wp
                end do
            end if
        end do 

        do ie = 1, ele_num_l
            if(btest(e_mask(ie), group).and.who_has_ele(ie)) then 
                do inod =1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do ibasis=1, basis_num(ie)
                        do j = 1,3
                            kme(2) = kme(2) + mass_mat_coeff(size_ele(ie), etype(ie))*0.5_wp*masses(basis_type(ibasis,ie)) &
                                              *const_motion*(vel(j,ibasis,ip) - avgvel(j))**2.0_wp
                        end do
                    end do 
                end do
            end if
        end do

        call mpi_allreduce(kme, kall, 2, mpi_wp, mpi_sum, world, ierr)

        compute_ke = kall
        return
    end function compute_ke

    function compute_temp(group, in_ke)
        !This function calculates the temperature for the current model
        integer, intent(in) :: group
        real(kind=wp), intent(in), optional :: in_ke(2)

        real(kind=wp) ::  compute_temp(3)
        real(kind=wp) :: ke(2)

        if(present(in_ke)) then 
            ke = in_ke
        else
            ke = compute_ke(group)
        end if

        compute_temp = 0

        compute_temp(1) = (2.0_wp*(ke(1)+ke(2)))/(3*(group_counts(1,group)+group_counts(4,group))*boltzmann)
        if (atom_num > 0) compute_temp(2) = (2.0_wp*(ke(1)))/(3*(group_counts(1,group))*boltzmann)
        if (ele_num > 0) compute_temp(3) = (2.0_wp*(ke(2)))/(3*group_counts(4,group)*boltzmann)
        return
    end function compute_temp

    function compute_box_press()
        !This calculates the box pressure. Can only be done on the full box which is why it takes no inputs
        real(kind=wp) :: compute_box_press(3,3)
            
        integer :: ia, ibasis, ip, inod, ie, i, j, k
        real(kind=wp) :: vir_buffme(9), vir_buffall(9), box_volume

        vir_buffme(:) = 0.0_wp
        !Sum virial stress in the atom region
        if(atom_num_l > 0) then 
            do ia = 1, atom_num_l
                k = 1
                do j = 1,3
                    do i = 1,3
                        vir_buffme(k) = vir_buffme(k) - virial_atom(i,j,ia)
                        k = k+1
                    end do
                end do
            end do
        end if

        if(ele_num_l > 0) then 
            do ie = 1, ele_num_l
                if(who_has_ele(ie)) then 
                    do inod = 1, ng_node(etype(ie))
                        ip = cg_node(inod,ie)
                        do ibasis = 1, basis_num(ie)
                            k = 1
                            do j = 1,3
                                do i = 1,3
                                    vir_buffme(k) = vir_buffme(k) - virial_eq(i,j,ibasis,ip)*mass_mat_coeff(size_ele(ie), itype(ie))
                                    k = k + 1
                                end do
                            end do
                        end do
                    end do  
                end if
            end do
        end if

        !Now sum the virial over all processors
        call mpi_allreduce(vir_buffme, vir_buffall, 9, mpi_wp, mpi_sum, world, ierr)
    
        box_volume = box_length(1)*box_length(2)*box_length(3)
        k = 1
        do j = 1,3
            do i = 1,3
                compute_box_press(i,j) = nktv2p*vir_buffall(k)/box_volume
            end do
        end do
        return
    end function compute_box_press

end module compute
