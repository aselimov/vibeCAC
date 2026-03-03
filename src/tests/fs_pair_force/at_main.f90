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

    implicit none 

    real(kind=dp) :: pe_real(41), fnorm_real(41), disp(41)
    real(kind=dp) :: pe_check, fnorm_check, tol
    integer :: i 

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
    !Results are calculated from lammps simulations. This just calculates the pair potential and force for a variety of distances
    pe_real=(/ 9999999999999.36_dp, 4999999999994.27_dp, 3333333333324.69_dp, 2499999999989.81_dp, 148563932270.089_dp, &
               1011.08997194053_dp, 596.038415317792_dp, 363.991600783288_dp, 228.016243702275_dp, 145.468794123426_dp, & 
               91.3444829326932_dp, 54.8476920689363_dp, 31.3780764754464_dp, 16.8229149338112_dp, 8.05110069207755_dp, &
               2.88312947957085_dp, -0.0700081594788218_dp, -1.68738220816374_dp, -2.60982936368697_dp, -3.16924095123063_dp, &
               -3.4803996851366_dp, -3.62425195185428_dp, -3.65562190775608_dp, -3.61068864697084_dp, -3.51286811646782_dp, &
               -3.37727142856666_dp, -3.2139597162888_dp, -3.03020065585815_dp, -2.83191620979554_dp, -2.62448984588025_dp, &
               -2.41309678346414_dp, -2.20271002002876_dp, -1.99790660926228_dp, -1.80257265434805_dp, -1.61961468076475_dp, &
               -1.45078133440576_dp, -1.29665778818291_dp, -1.15689080270187_dp, -1.03069016385057_dp, -0.917639038816428_dp,&
               -0.818342471722933_dp /)

    fnorm_real=(/ 141421356237402.0_dp, 35355339059381.7_dp, 15713484026397.7_dp, 8838834764846.81_dp, 3.04764252645022e15_dp, &
                  7871.00817142953_dp, 4289.4479966714_dp, 2466.11632252638_dp, 1475.35590761412_dp, 910.148585089464_dp, &
                  630.101061629911_dp, 413.34959866945_dp, 260.266808620324_dp, 158.777274773767_dp, 94.4632351559387_dp, &
                  54.8823819422511_dp, 30.6868752165717_dp, 16.2946598921435_dp, 10.1635149969605_dp, 5.92831633941791_dp,&
                  3.06130912073685_dp, 1.13489763162427_dp, 0.163729290392537_dp, 1.05303021909981_dp, 1.67889058367453_dp, &
                  2.13316884751354_dp, 2.46918733987136_dp, 2.71452780729784_dp, 2.88118653809808_dp, 2.97351396226005_dp, &
                  2.99381142485739_dp, 2.94596135211483_dp, 2.8375196719926_dp, 2.6804653839289_dp, 2.49024243781296_dp, &
                  2.28376032165085_dp, 2.07645802246115_dp, 1.87880785158624_dp, 1.6919163099495_dp, 1.50331195438987_dp,&
                  1.30921316226896_dp /)

    command='read_data Cu.restart'
    call parse_read(command)
    boundary='sss'
    periodic=.false.
    period=(/.false., .false., .false. /)
    command='potential fs CuZr_mm.eam.fs Cu Zr'
    call parse_potential(command)
    command='neighbor 1.5'
    call parse_neighbors(command)
    call pre_calc(1)
    call update_neighbor(0,.false., .true.)
    r_atom(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)

    do i=1,41
        r_atom(:,2) = (/ 0.0_dp, i*0.1_dp, 0.0_dp /)
        call update_force
        if (i < 6) then
            tol = 0.1
        else 
            tol = 1d-8
        end if
        if (abs(pe_real(i) - sum(energy_atom(1:2)))>tol) then 
            print *, "[FAILED] - Potential energies don't match for reference=",pe_real(i), ' and calculated=', sum(energy_atom)
            call mpi_abort( mpi_comm_world, 1, ierr)
        end if
        fnorm_check=sqrt(force_atom(1,1)*force_atom(1,1) + force_atom(2,1)*force_atom(2,1) + force_atom(3,1)*force_atom(3,1) + &
              force_atom(1,2)*force_atom(1,2) + force_atom(2,2)*force_atom(2,2) + force_atom(3,2)*force_atom(3,2))

        if (abs(fnorm_real(i) - fnorm_check)>tol*100) then 
            print *, "[FAILED] - Fnorm don't match for reference=",fnorm_real(i), ' and calculated=', fnorm_check 
            call mpi_abort( mpi_comm_world, 1, ierr)
        end if
    end do
    
    print *, '[OK] - Passed single element energy and force calculations'

    !close log
    call log_time
    call close_log

    call mpi_finalize(ierr)

end program at_main
