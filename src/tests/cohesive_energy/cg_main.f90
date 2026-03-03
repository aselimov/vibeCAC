program cg_main
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

    implicit none 

    real(kind=dp) :: pe_real, fnorm_real
    real(kind=dp) :: pe_check, fnorm_check
    logical :: test_results(2)

    character(len=100) :: command

    !initialize mpi
    call comm_init
    !Initialize time
    call time_defaults

    !Start timer
    program_start = mpi_wtime()

    !Initialize log
    log_on=.false.
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
    !Results are calculated from lammps simulations. These are cohesive energy calculations meaning that a pure block of Cu is
    !generated with a lattice constant of 3.615 \AA and the resulting energy and force_norm are copied here 
    
    !Now read in the data, we pass commands using strings to the necessary functions
    !Start with the pure atomistic case 
    test_results=.false.


    command='read_data Cu.restart'
    call parse_read(command)
    command='potential eam Cu_mishin1.eam.alloy Cu'
    call parse_potential(command)
    command='neighbor 1.5'
    call parse_neighbors(command)
    command='run 0'
    call pre_calc(1)
    call parse_run(command)

    pe_check=compute_pe(1)
    fnorm_check = compute_fnorm(1)

    pe_real = -278728.47
    fnorm_real = 1.6437592d-11

    if (rank==root) then 
        if (abs(pe_check - pe_real) < 0.01) then
            print *, '[OK] - Potential energy calculation success for coarse-grained model'
            test_results(1) = .true.
        else
            print *, '[FAILED] - Potential energy calculation failure for coarse-grained model. CAC = ', pe_check, &
                     ' reference = ',pe_real
        end if

        if (abs(fnorm_check - fnorm_real) < 1d-5) then
            print *, '[OK] - Force norm calculation success for coarse-grained model'
            test_results(2) = .true.
        else
            print *, '[FAILED] - Force norm calculation failure for coarse-grained model. CAC = ', fnorm_check, &
                     ' reference = ', fnorm_real
        end if

    end if

    !close log
    call log_time
    call close_log

    call mpi_bcast(test_results, int(size(test_results)), mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_finalize(ierr)

    if (all(test_results)) then 
        stop
    else 
        stop 2
    end if

end program cg_main
