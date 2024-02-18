module morse
    
    !subroutines for pairtab potential style
    
    use parameters
    use math 
    use forces 
    use elements 
    use integration 
    use neighbors 
    use comms 
    use errors
    use atom_types

    implicit none

    real(kind=wp), allocatable, private, save :: pair_spline(:,:,:)
    real(kind=wp), private, dimension(max_atom_types, max_atom_types), save :: pair_cutoffsq, d0, alpha, r0, morse1

    public 
    contains

    subroutine add_morse_potential(type1, type2, ind0, inalpha, inr0, incutoff)
        integer, intent(in) :: type1, type2
        real(kind=wp), intent(in) :: ind0, inalpha, inr0, incutoff

        integer :: i,j, mn

        d0(type1, type2) = ind0
        d0(type2, type1) = ind0
        alpha(type1, type2) = inalpha
        alpha(type2, type1) = inalpha
        r0(type1, type2) = inr0
        r0(type2, type1) = inr0
        pair_cutoffsq(type1, type2)=incutoff*incutoff
        pair_cutoffsq(type2, type1)=incutoff*incutoff
        morse1(type1, type2) = 2*inalpha*ind0
        morse1(type2, type1) = 2*inalpha*ind0

        if((types_to_pot_type(type1, type2) > 0).or.(types_to_pot_type(type2,type1) > 0)) then 
            call misc_error("Can't define multiple potentials for single pair interaction")
        end if
        types_to_pot_type(type1,type2)=ibset(types_to_pot_type(type1,type2),2)
        types_to_pot_type(type2,type1)=ibset(types_to_pot_type(type2,type1),2)

        if(.not.def_nei) then 
            def_nei=.true.
            rc_neigh=incutoff
        else
            rc_neigh=max(rc_neigh, incutoff)
        end if

    end subroutine add_morse_potential
    
    subroutine update_force_morse
        integer :: i, ia, ie, iep, jep, ja, iatomap, jatomap, iatom, latomtype, katomtype, nei, mn_list, ibasis, ic, inod, ip
        real(kind=wp) :: rl(3), rk(3), rlk(4), dr, d_exp, pair, flk(3), eshape, rsq

        real(kind=wp), dimension(max_basisnum*max_intpo_num, ele_num_l) :: energy_intpo
        real(kind=wp), dimension(3, max_basisnum*max_intpo_num, ele_num_l) :: force_intpo
        real(kind=wp), dimension(3,3, max_basisnum*max_intpo_num, ele_num_l) :: virial_intpo

        mn_list=0
        do i = 1,2
            if(potential_types(i)) mn_list=mn_list+1
        end do 
        energy_intpo=0.0_wp
        force_intpo=0.0_wp
        virial_intpo=0.0_wp

        !First calculate energy/force for atoms
        if (atom_num > 0) then 
            do ia=1, atom_num_l
                rl(:)=r_atom(:, ia)
                latomtype=type_atom(ia)

                do nei = 1, n_at_at(ia, mn_list)
                    ja= at_nei_at(nei,ia,mn_list)
                    rk=r_atom(:, ja)
                    katomtype=type_atom(ja)

                    rlk(1:3)=rl-rk
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < pair_cutoffsq(latomtype, katomtype)) then 
                        rlk(4)=sqrt(rsq)
                        dr=rlk(4)-r0(latomtype, katomtype)
                        d_exp=exp(-alpha(latomtype, katomtype)*dr)
                        flk=rlk(1:3)*(morse1(latomtype, katomtype)*(d_exp*d_exp-d_exp)/rlk(4))
                        pair=d0(latomtype, katomtype)*(d_exp*d_exp-2.0_wp*d_exp) 

                        energy_atom(ia)=energy_atom(ia)+pair/2.0_wp

                        force_atom(:, ia) = force_atom(:, ia) + flk

                        if(need_virial) then
                            do ic = 1, 3
                                virial_atom(:, ic, ia) = virial_atom(:, ic, ia) + flk(:) * rlk(ic) / 2.0_wp
                            end do
                        end if
                        if (ja <= atom_num_l) then 
                            energy_atom(ja)=energy_atom(ja) + pair/2.0_wp
                            force_atom(:, ja) = force_atom(:, ja) - flk
                            if(need_virial) then
                                do ic = 1, 3
                                    virial_atom(:, ic, ja) = virial_atom(:, ic, ja) + flk(:) * rlk(ic) / 2.0_wp
                                end do
                            end if
                        end if
                    end if
                end do
                do nei=1, n_at_cg(ia, mn_list)
                    ja = at_nei_cg(nei, ia, mn_list)
                    rk(:) = r_atomap(:, ja)
                    katomtype = type_atomap(ja)

                    rlk(1:3) = rl - rk
                    rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                    if(rsq < pair_cutoffsq(latomtype, katomtype)) then 
                        rlk(4)=sqrt(rsq)
                        dr=rlk(4)-r0(latomtype, katomtype)
                        d_exp=exp(-alpha(latomtype, katomtype)*dr)
                        flk=rlk(1:3)*(morse1(latomtype, katomtype)*(d_exp*d_exp-d_exp)/rlk(4))
                        pair=d0(latomtype, katomtype)*(d_exp*d_exp-2.0_wp*d_exp) 
                        energy_atom(ia)=energy_atom(ia)+pair/2.0_wp
                        force_atom(:, ia) = force_atom(:, ia) + flk

                        if(need_virial) then
                            do ic = 1, 3
                                virial_atom(:, ic, ia) = virial_atom(:, ic, ia) + flk(:) * rlk(ic) / 2.0_wp
                            end do
                        end if

                        if(ja <= atomap_num_l) then 
                            if(atomap_to_intpo(1,ja)/=0) then 
                                jep = atomap_to_intpo(1,ja)
                                ie = atomap_to_intpo(2,ja)

                                force_intpo(:,jep,ie) = force_intpo(:,jep,ie) - flk(:)
                                if(need_virial)then 
                                    do ic = 1, 3
                                        virial_intpo(:, ic, jep,ie) = virial_intpo(:, ic, jep,ie) + flk(:) * rlk(ic) / 2.0_wp
                                    end do
                                end if
                                energy_intpo(jep,ie) = energy_intpo(jep,ie) + pair/2.0_wp
                            end if

                        end if
                    end if
                end do
            end do 
        end if

        if(ele_num_l > 0) then 
            do ie=1, ele_num_l
                do iep=1, intpo_count(itype(ie))
                    do ibasis=1, basis_num(ie)
                        jep=basis_num(ie)*(iep-1)+ibasis
                        if(who_has_intpo(iep, ie)) then 
                            iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                            iatomap = cg_atomap(basis_num(ie)*(iatom-1) + ibasis, ie)
                            latomtype = type_atomap(iatomap)
                            rl(:) = r_atomap(:, iatomap) 
                            
                            do nei=1, n_cg_cg(jep, ie, mn_list)
                                ja=cg_nei_cg(nei, jep, ie, mn_list)

                                rk(:) = r_atomap(:, ja)
                                katomtype = type_atomap(ja)

                                rlk(1:3) = rl - rk
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < pair_cutoffsq(latomtype,katomtype)) then 

                                    rlk(4) = sqrt(rsq)

                                    dr=rlk(4)-r0(latomtype, katomtype)
                                    d_exp=exp(-alpha(latomtype, katomtype)*dr)
                                    flk=rlk(1:3)*(morse1(latomtype, katomtype)*(d_exp*d_exp-d_exp)/rlk(4))
                                    pair=d0(latomtype, katomtype)*(d_exp*d_exp-2.0_wp*d_exp) 

                                    force_intpo(:, jep, ie) = force_intpo(:, jep, ie) + flk(:)
                                    if(need_virial) then
                                        do ic = 1, 3
                                            virial_intpo(:, ic, jep, ie) = virial_intpo(:, ic, jep, ie) + flk(:) * rlk(ic) / 2.0_wp
                                        end do
                                    end if
                                    energy_intpo(jep, ie) = energy_intpo(jep, ie)+ pair/2.0_wp
                                end if
                            end do

                            do nei=1, n_cg_at(jep, ie, mn_list)

                                ja = cg_nei_at(nei, jep, ie, mn_list)
                                rk(:) = r_atom(:, ja)
                                katomtype = type_atom(ja)

                                rlk(1:3) = rl - rk
                                rsq = rlk(1)*rlk(1) + rlk(2)*rlk(2) + rlk(3)*rlk(3)
                                if(rsq < pair_cutoffsq(latomtype,katomtype)) then 
                                    rlk(4) = sqrt(rsq)

                                    dr=rlk(4)-r0(latomtype, katomtype)
                                    d_exp=exp(-alpha(latomtype, katomtype)*dr)
                                    flk=rlk(1:3)*(morse1(latomtype, katomtype)*(d_exp*d_exp-d_exp)/rlk(4))
                                    pair=d0(latomtype, katomtype)*(d_exp*d_exp-2.0_wp*d_exp) 

                                    force_intpo(:, jep, ie) = force_intpo(:, jep, ie) + flk(:)
                                    if(need_virial) then
                                        do ic = 1, 3
                                            virial_intpo(:, ic, jep, ie) = virial_intpo(:, ic, jep, ie) + flk(:) * rlk(ic) / 2.0_wp
                                        end do
                                    end if
                                    energy_intpo(jep, ie) = energy_intpo(jep, ie)+ pair/2.0_wp

                                end if
                            end do
                        end if
                    end do
                end do

                do iep= 1, intpo_count(itype(ie))
                    do ibasis=1, basis_num(ie)
                        jep=basis_num(ie)*(iep-1)+ibasis
                        iatom=atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                        iatomap=cg_atomap(basis_num(ie)*(iatom-1)+ ibasis, ie)

                        if(who_has_intpo(jep, ie).eqv..true.) then
                            force_intpo(:, jep, ie) = force_intpo(:, jep, ie) &
                                                      * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))

                            if(need_virial) then
                              virial_intpo(:, :, jep, ie) = virial_intpo(:, :, jep, ie) &
                                                      * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                            end if
                            energy_intpo(jep, ie) = energy_intpo(jep, ie) &
                                                    * weight_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                        end if
                    end do
                end do

                do inod = 1, ng_node(etype(ie))
                    ip = cg_node(inod, ie)
                    do iep = 1, intpo_count(itype(ie))
                        do ibasis = 1, basis_num(ie)
                            jep = basis_num(ie)*(iep-1) + ibasis
                            if(who_has_intpo(iep, ie).eqv..true.) then

                                iatom = atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie))
                                eshape = a_interpo(inod, iatom, size_to_shape(size_ele(ie)), etype(ie))
                                force(:, ibasis, ip ) = force(:, ibasis, ip) + eshape * force_intpo(:, jep, ie)

                                if(need_virial) then
                                    virial(:, :, ibasis, ip) = virial(:, :, ibasis, ip) &
                                                       + eshape * virial_intpo(:, :, jep, ie)
                                end if

                                energy(ibasis, ip) = energy(ibasis, ip) + eshape*energy_intpo(jep, ie)

                            end if
                        end do
                    end do
                end do

            end do
        end if
    end subroutine update_force_morse

end module morse
