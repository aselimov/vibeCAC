program at_main
    use mpi
    use comms
    use input_parser 
    use integration
    use neighbors
    use logger
    use minimize 
    use time
    use group
    use dump
    use thermo
    use temp
    use berendsen
    use potential 
    use debug
    use deform
    use langevin
    use read_data
    use compute
    use read_dump
    use parameters

    implicit none 

    real(kind=dp) :: pe_real, fnorm_real, pe_check_sum, pe_real_sum
    real(kind=dp) :: pe_check, fnorm_check
    integer(kind=dp), allocatable :: tags_map_real(:), tags_map_check(:)
    integer :: i,j, real_ind
    logical :: failed, test_results(6)

    character(len=100) :: command

    !initialize mpi
    call comm_init
    !Initialize time
    call time_defaults

    !Start timer
    program_start = mpi_wtime()

    !Initialize log
    log_on=.true.
    call init_log

    !Initialize some default values
    deform_num = 0

    !Initialize some flags 
    tflag = .false.
    pflag = .false.
    dflag = .false.
    first_run = .true.

    !Now call all necessary default setting subroutines
    call potential_defaults
    call integration_defaults 
    call neighbors_defaults
    call min_defaults
    call group_defaults
    call thermo_defaults
    call dump_defaults
    call langevin_defaults
    call dynamics_defaults

!   ------------------------------------------------------------
!   Execute tests 
!   ------------------------------------------------------------
    !Results are calculated from lammps simulations. These calculate the stacking fault for Cu at a displacement of 7 angstroms 
    !along the 211 direction only allowing displacement along the z direction
    
    !Now read in the data, we pass commands using strings to the necessary functions
    !Start with the pure atomistic case 
    test_results=.false.
    command='read_data post_swap.restart'
    call parse_read(command)
    command='types Fe S Ni Cr He H'
    call parse_types(command)
    command='potential fs FeNiCrH19_Aug12_V.t'
    call parse_potential(command)
    boundary='ppp'
    period=(/ .true., .true., .true. /)

    command='potential morse 1 5 0.000529628 1.410836 4.657221 3.0'
    call parse_potential(command)
    command='potential morse 3 5 0.001324599 1.412238 4.347442 3.0'
    call parse_potential(command)
    command='potential morse 4 5 0.000811794 1.399679 4.517206 3.0'
    call parse_potential(command)
    command='potential morse 5 5 0.000204048 4.655408 2.490410 1.9'
    call parse_potential(command)
    command='potential morse 5 6 0.038285972 1.493264 2.398987 2.0'
    call parse_potential(command)

    command='potential morse 1  2 0.00000000 0.000000 0.000000 0.000 '
    call parse_potential(command)
    command='potential morse 2  2 0.00000000 0.000000 0.000000 0.000'
    call parse_potential(command)
    command='potential morse 3  2 0.00000000 0.000000 0.000000 0.000 '
    call parse_potential(command)
    command='potential morse 4  2 0.00000000 0.000000 0.000000 0.000 '
    call parse_potential(command)
    command='potential morse 5  2 0.00000000 0.000000 0.000000 0.000 '
    call parse_potential(command)
    command='potential morse 6  2 0.00000000 0.000000 0.000000 0.000 '
    call parse_potential(command)

    command='neighbor 1.5'
    call parse_neighbors(command)
    command='run 0'
    call pre_calc(1)
    call parse_run(command)

    pe_check=compute_pe(1)
    fnorm_check=compute_fnorm(1)
    pe_real=-5986558.6831048_wp
    fnorm_real=1239.54116121321_wp
    failed=.false.
    if (rank==root) then 
        if (abs(pe_check-pe_real) < 1d-10) then
            print *, '[OK] - Potential energy calculation failed for atomistic model'
        else
            print *, '[FAILED] - Potential energy calculation failed for atomistic model with real=', pe_real, &
                     ' and calculated=', pe_check
            failed=.true.
        end if
        if (abs(fnorm_check-fnorm_real) < 1d-10) then
            print *, '[OK] - Force norm calculation failed for atomistic model'
        else
            print *, '[FAILED] - Force norm calculation failed for atomistic model with real=', fnorm_real, &
                     ' and calculated=', fnorm_check
            failed=.true.
        end if
    end if
    
    if(failed .and.pro_num==1) then 
        allocate(tags_map_check(atom_num), tags_map_real(atom_num))
        tags_map_check=0
        tags_map_real=0
        call read_in_dump('0.dump')

        do i=1,atom_num
            tags_map_check(tag_atom(i)) = i
            tags_map_real(read_atom_id(i)) = i
        end do

        pe_check_sum=0.0_dp
        pe_real_sum =0.0_dp

        do i=1, atom_num
            !Check position
            real_ind=tags_map_real(tag_atom(i)) 
            if(tag_atom(i) /= read_atom_id(real_ind)) then 
                print *, "Check and real tags don't match with check tags", tag_atom(i), ' and real tags ', read_atom_id(real_ind)
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
            do j=1,3
                if(abs(r_atom(j, i) - read_atom_pos(j, real_ind))> 1d-10) then 
                    print *, "[FAILED] - Positions don't match for atom ", tag_atom(i), ' for dimension ', j, ' with check ', &
                             r_atom(j,i), ' and real ', read_atom_pos(j,real_ind)
                    call mpi_abort(mpi_comm_world, 1, ierr)
                end if
            end do
            do j=1,3
                if(abs(force_atom(j, i) - read_atom_force(j, real_ind))> 1d-12) then 
                    print *, "[FAILED] - Forces don't match for atom ", tag_atom(i), ' for dimension ', j, ' with check ', &
                             force_atom(j,i), ' and real ', read_atom_force(j,real_ind)
                    call mpi_abort(mpi_comm_world, 1, ierr)
                end if
            end do
            
            pe_check_sum = pe_check_sum + energy_atom(i)
            pe_real_sum = pe_real_sum + read_atom_energy(real_ind)
            if(abs(energy_atom(i) - read_atom_energy(real_ind))> 1d-12) then 
                print *, "[FAILED] - Energies don't match for atom ", tag_atom(i), ' for dimension ', j, ' with check ', &
                         energy_atom(i), ' and real ', read_atom_energy(real_ind)
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if
        end do
    end if
    
    call log_time
    call close_log

    call mpi_finalize(ierr)


end program at_main


