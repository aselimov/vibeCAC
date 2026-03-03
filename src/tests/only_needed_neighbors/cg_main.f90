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
    use refinement
    use langevin
    use read_data
    use compute

    implicit none 

    real(kind=dp) :: pre_eng_real, post_eng_real, pre_fnorm_real, post_fnorm_real, &
                     pre_eng_check, post_eng_check, pre_fnorm_check, post_fnorm_check, time1, time2, times(3)
    logical :: test_results(4)

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
    refinement_check = 0
    times= (/ 50.212986147999999_dp, 15.410064985000000_dp, 10.219318316999999_dp /)
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
    pre_eng_real=-1163009.8481146714_dp
    pre_fnorm_real=1.6997644167762397E-002_dp
    post_eng_real=-1163007.6250000000
    post_fnorm_real= 9.5590002449639979E-004


    time1= mpi_wtime()
    command='read_data Cu.restart'
    call parse_read(command)

    command='potential eam Cu_mishin1.eam.alloy Cu'
    call parse_potential(command)

    command='neighbor 1.5 no_need_all_step'
    call parse_neighbors(command)

    boundary='pps'
    period=(/ .true., .true., .false. /)

    command='group top all block -inf inf -inf inf 112 inf'
    call parse_group(command)

    command='setforce all 0 0 NULL'
    call init_set_force 
    call parse_set_force(command)

    command='thermo 1'
    call parse_thermo(command)

    command='run 0'
    call pre_calc(1)
    call parse_run(command)
    pre_eng_check = compute_pe(1)
    pre_fnorm_check=compute_fnorm(1)
    
    if (abs(pre_eng_check-pre_eng_real) > 1d-8) then
        print *, '[FAILED] - Failed initial energy calculation with reference =',pre_eng_real, ' and calculated =', pre_eng_check
    else
        print *, '[OK] - Passed initial energy calculation'
        test_results(1)=.true.
    end if

    if (abs(pre_fnorm_check-pre_fnorm_real) > 1d-8) then
        print *, '[FAILED] - Failed initial fnorm calculation with reference =',pre_fnorm_real,' and calculated =',pre_fnorm_check
    else
        print *, '[OK] - Passed initial fnorm calculation'
        test_results(2)=.true.
    end if

    command='min_style fire'
    call parse_min_style(command)
    command='minimize 0 1e-3 10000'
    call pre_calc(2)
    call parse_minimize(command)

    post_eng_check = compute_pe(1)
    post_fnorm_check=compute_fnorm(1)

    if (abs(post_eng_check-post_eng_real) > 1d-8) then
        print *, '[FAILED] - Failed post min energy calculation with reference =', post_eng_check, &
                 ' and calculated =', post_eng_real
    else
        print *, '[OK] - Passed post min energy calculation'
        test_results(3)=.true.
    end if

    if (abs(post_fnorm_check-post_fnorm_real) > 1d-8) then
        print *, '[FAILED] - Failed post min fnorm calculation with reference =', post_fnorm_check, &
                 ' and calculated =', post_fnorm_real
    else
        print *, '[OK] - Passed post min fnorm calculation'
        test_results(4)=.true.
    end if

    time2= mpi_wtime()

    print *, time2-time1
    if (pro_num==1) then
        if (time2-time1 < times(1)) then
            print *, '[OK] - Ran faster than default'
        else 
            print *, '[FAILED] - Ran slower than default with calculated =', time2-time1, ' and reference =', times(1)
        end if
    else if(pro_num==4) then 
        if (time2-time1 < times(2)) then
            print *, '[OK] - Ran faster than default'
        else 
            print *, '[FAILED] - Ran slower than default with calculated =', time2-time1, ' and reference =', times(2)
        end if
    else if(pro_num==8) then 
        if (time2-time1 < times(3)) then
            print *, '[OK] - Ran faster than default'
        else 
            print *, '[FAILED] - Ran slower than default with calculated =', time2-time1, ' and reference =', times(3)
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
