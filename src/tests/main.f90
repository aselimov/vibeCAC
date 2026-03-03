program main
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
    use refinement
    use langevin

    implicit none 
    !initialize mpi
    call comm_init
    !Initialize time
    call time_defaults

    !Start timer
    program_start = mpi_wtime()

    !Initialize log
    call init_log

    !Initialize some default values
    deform_num = 0
    refinement_check = 0

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

    ! Execute test

    !close log
    call log_time
    call close_log

    call mpi_finalize(ierr)

   stop
end program main
