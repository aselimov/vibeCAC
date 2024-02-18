! Copyright (c) 2017-2018 Georgia Institute of Technology. All Rights Reserved
! Redistributing this source code is prohibited. This is a testing version of CAC.

! This source code is provided as is, with no warranties or representations of accuracy 
! or suitability for any application, and with no expectation of user support.

! Please alert Alex Selimov (aselimov3@gatech.edu) with any bugs or changes made to source code.
! Written by Shuozhi Xu (shuozhixu@ucsb.edu)
program main

    use mpi
    use comms
    use parameters
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

    implicit none

    !initialize mpi
    call comm_init
    !Initialize time
    call time_defaults

    !Start timer
    program_start = mpi_wtime()

    !Initialize log
    call log_defaults
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

    !Run input parsing code
    call read_input

    !close log
    call log_time
    call close_log

    call mpi_finalize(ierr)

   stop
   end program main
