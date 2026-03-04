module potential

    use parameters
    use eam
    use morse
    use forces
    use comms
    use str
    use atom_types
    use force_mod

    implicit none

    !This code handles controlling which potentials are called and processes potential input
    integer :: pnum, enum, fnum 

    abstract interface
        subroutine force_adapter_proc
        end subroutine force_adapter_proc

        real(kind=wp) function walltime_adapter_proc()
            import :: wp
        end function walltime_adapter_proc
    end interface

    procedure(force_adapter_proc), pointer :: pre_force_adapter => pre_force
    procedure(force_adapter_proc), pointer :: post_force_adapter => post_force
    procedure(force_adapter_proc), pointer :: update_force_eam_adapter => update_force_eam
    procedure(force_adapter_proc), pointer :: update_force_morse_adapter => update_force_morse
    procedure(force_adapter_proc), pointer :: update_equiv_adapter => update_equiv
    procedure(walltime_adapter_proc), pointer :: walltime_reader => default_walltime

    !Spline arrays
    public 
    contains
    
    subroutine potential_defaults
        potential_types=.false.
        rdr = 0.0_wp
        rdrho = 0.0_wp
        nr = 0.0_wp
        nrho = 0.0_wp
        pnum = 0 
        enum = 0 
        fnum = 0
        type_to_pair = 0
        type_to_dens = 0 
        pot_map = 0
        call reset_update_force_adapters
    end subroutine potential_defaults

    subroutine parse_potential(line)
        character(len = *), intent(in) :: line
        character(len = 100) :: tmptxt, tmptxt2
        character(len = 100) potential_type
        character(len = read_len) :: filename
        character(len=2), dimension(max_atom_types) :: names
        integer :: iospara, i, j, mn, itype, jtype, lim
        real(kind=wp) :: ind0, inalpha, inr0,cutoff

        real(kind=wp), allocatable :: mpair(:)


        read(line, *, iostat = iospara) tmptxt, potential_type
        select case(potential_type)
        case('eam','fs')

            potential_types(1)=.true.

            read(line, *, iostat = iospara) tmptxt, potential_type, filename
            if(iospara> 0) call read_error("Error parsing potential command", iospara)
            !Now if we have more than 2 tokens we read them in as atom types unless we have already set the atom types
            j = tok_count(line)-3
            if ((j > 0).and.(.not.atom_types_set)) then
                read(line, *, iostat = iospara) tmptxt, potential_type, filename, (names(i), i = 1, j)

                !and pass it to atom_types 
                call set_atom_types(j, names)
            end if                   
            
            if(potential_type == "eam") then 
                call eam_lammps(filename, 1)
            else if(potential_type == "fs") then 
                call eam_lammps(filename, 1, .true.)
            end if

            call eamarray2spline

            !Set cutoff radiuses
            rc_off = min(d_finish, p_finish)
            rc_sq = rc_off**2
            rc_min = max(d_start, p_start, lim_small)
            if (.not.def_nei) then 
                def_nei = .true.
                rc_neigh = rc_off
            end if

            call set_eam_map_arrays

        case('morse')
            potential_types(2)=.true.
            !atom types have to be set before morse potential is called, it usually needs to be set with the types command

            if(.not.atom_types_set) call misc_error("Atom types need to be set using types command before morse is used")

            read(line, *) tmptxt, tmptxt2, itype, jtype, ind0, inalpha, inr0, cutoff
            call add_morse_potential(itype, jtype, ind0, inalpha, inr0, cutoff)

        case default
            call command_error(trim(adjustl(potential_type))//" is not an acceptable potential type")
        end select

        !Log the atom types
        !call log_types


    end subroutine parse_potential

    subroutine alloc_potential_arrays
        !Initialize potential arrays
        if (potential_types(1)) then
            call alloc_eam_arrays
        end if
    end subroutine alloc_potential_arrays

    subroutine update_force

        real(kind=wp) :: t_start, t_end

        call pre_force_adapter
        t_start = walltime_reader()

        !First zero arrays
        if(ele_num > 0) then 
            if(allocated(force)) then 
                if(node_num_l > size(force, 3)) then
                    deallocate(force, energy)
                    allocate(force(3, max_basisnum, node_num_l), energy(max_basisnum, node_num_l), stat = allostat)
                    if(allostat /= 0) call alloc_error("Failure to allocate force/energy", allostat)
                    if(need_virial) then 
                        deallocate(virial)
                        allocate(virial(3, 3, max_basisnum, node_num_l), stat = allostat)
                        if(allostat /= 0) call alloc_error("Failure to allocate virial", allostat)
                    end if
                end if 
            else
                allocate(force(3, max_basisnum, node_num_l), energy(max_basisnum, node_num_l), stat = allostat)
                if(allostat /= 0) call alloc_error("Failure to allocate force/energy", allostat)
                if(need_virial) then 
                    allocate(virial(3, 3, max_basisnum, node_num_l), stat = allostat)
                    if(allostat /= 0) call alloc_error("Failure to allocate virial", allostat)
                end if
            end if

            !Initialize force and energy
            force(:, :, :) = 0.0_wp
            energy(:, :) = 0.0_wp
            !Allocate virial if needed
            if(need_virial) then
                virial(:, :, :, :) = 0.0_wp
            end if

        end if

        if(atom_num > 0) then 
            force_atom=0.0_wp
            energy_atom=0.0_wp
            if(need_virial) then 
                virial_atom=0.0_wp
            end if
        end if

        if(potential_types(1)) call update_force_eam_adapter
        
        if(potential_types(2)) call update_force_morse_adapter

        !Now call update_equiv to update equivalent node value arrays if there are finite elements
        if(ele_num > 0) call update_equiv_adapter
        t_end = walltime_reader()

        walltime(2) = walltime(2) + (t_end-t_start)
        call post_force_adapter
    end subroutine update_force

    subroutine set_update_force_adapters(pre_adapter, post_adapter, eam_adapter, morse_adapter, equiv_adapter, walltime_adapter)
        procedure(force_adapter_proc), optional :: pre_adapter, post_adapter, eam_adapter, morse_adapter, equiv_adapter
        procedure(walltime_adapter_proc), optional :: walltime_adapter

        if (present(pre_adapter)) pre_force_adapter => pre_adapter
        if (present(post_adapter)) post_force_adapter => post_adapter
        if (present(eam_adapter)) update_force_eam_adapter => eam_adapter
        if (present(morse_adapter)) update_force_morse_adapter => morse_adapter
        if (present(equiv_adapter)) update_equiv_adapter => equiv_adapter
        if (present(walltime_adapter)) walltime_reader => walltime_adapter
    end subroutine set_update_force_adapters

    subroutine reset_update_force_adapters
        pre_force_adapter => pre_force
        post_force_adapter => post_force
        update_force_eam_adapter => update_force_eam
        update_force_morse_adapter => update_force_morse
        update_equiv_adapter => update_equiv
        walltime_reader => default_walltime
    end subroutine reset_update_force_adapters

    real(kind=wp) function default_walltime()
        default_walltime = mpi_wtime()
    end function default_walltime

    subroutine pre_force
        return
    end subroutine pre_force

    subroutine post_force
        if(set_force_num > 0) call run_set_force
        if((add_force_num > 0)) call add_force
    end subroutine post_force

!    subroutine pairtab(line)
!        !This subroutine reads in tabulated morse potential files
!        character(len=*), intent(in) :: line
!        
!        integer :: i, itype, jtype, mn
!        real(kind=wp) :: mdr
!        real(kind=wp), allocatable :: mpair(:)
!        character(len = read_len) :: tmptxt, filename
!
!        read(line, *) tmptxt, itype, jtype, filename
!        open(unit=21, file=trim(adjustl(filename)), status='old', action='read', position='rewind')
!        !First read the step size and number of units, these start from the first step
!
!        read(21, *) mn, mdr
!
!        !Now check to make sure the morse step size is equivalent to the eam step size if the eam has already been defined
!        if(rdr > 0.0_wp) then 
!            if(.not.is_equal(1/mdr, rdr)) call misc_error("Step size for morse potential must match step size for eam potential")
!        else 
!            rdr = 0.0
!        end if
!
!        !Now read in the potential table
!        allocate(mpair(mn))
!        do i = 1, mn
!            read(21,*) mpair(i) 
!        end do
!
!        !Now add the pair potential to the spline arrays
!        call spline_arrays(mn,mn,mn,1,1,1)
!        !
!        embed_spline(:,:, size(embed_spline,3)) = 0.0_wp
!        edens_spline(:,:, size(edens_spline,3)) = 0.0_wp
!
!        call interpolate(mn, mdr, mpair, pair_spline(:,:,size(pair_spline,3)))
!        
!        return
!
!    end subroutine pairtab

end module potential
