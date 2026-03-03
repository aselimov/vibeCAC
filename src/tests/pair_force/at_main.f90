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
    real(kind=dp) :: pe_check, fnorm_check
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
    pe_real=(/ 15388.0454684558_dp, 8374.6283633242_dp, 4569.8422870799_dp, 2512.12811789192_dp, 1400.00150466955_dp, &
               796.12509483028_dp, 463.731883455318_dp, 276.130567901411_dp, 166.588487816659_dp, 100.546070859581_dp, &
               60.1843215976997_dp, 35.8340691418639_dp, 21.2951562604956_dp, 12.4145839028601_dp, 6.83956861996021_dp, &
               3.27722804969108_dp, 1.00769113647554_dp, -0.395080444281487_dp, -1.22070289920207_dp, -1.67257364105928_dp, &
               -1.87943987981488_dp, -1.92457707963509_dp, -1.86582300804121_dp, -1.745568600395_dp, -1.59280137044289_dp, &
               -1.42655726420376_dp, -1.25912837159664_dp, -1.09820422662717_dp, -0.948324617700717_dp, -0.811881855102902_dp, &
               -0.689821070616971_dp, -0.582133404766365_dp, -0.48820420492193_dp, -0.407058046788199_dp, -0.337529489742935_dp, &
               -0.278380045487459_dp, -0.228376154447178_dp, -0.186339007877918_dp, -0.151174220067881_dp, -0.121887278177993_dp, & 
               -0.0975891508733144_dp /)

    fnorm_real=(/ 132443.86782818_dp, 71971.4371927161_dp, 38969.0108217772_dp, 21054.7114572366_dp, 11395.1195173019_dp, &
                  6222.86002535959_dp, 3466.01492027412_dp, 1991.20972943528_dp, 1185.10342520004_dp, 722.929909830406_dp,&
                  440.386389723179_dp, 262.624106684267_dp, 158.318324494481_dp, 98.2153211012881_dp, 62.4076280271208_dp,&
                  39.9843914532297_dp, 25.1755541378172_dp, 15.1715850347069_dp, 8.65680579604669_dp, 4.41996586841663_dp,&
                  1.62281434035806_dp, 0.212506305689457_dp, 1.34876677066294_dp, 1.98532227370447_dp, 2.29145713390841_dp,&
                  2.3821557200658_dp, 2.33565882981283_dp, 2.20549437663258_dp, 2.02826438741258_dp, 1.82880624253313_dp, &
                  1.62369567931723_dp, 1.42367702433597_dp, 1.23538061181006_dp, 1.06255338636312_dp, 0.906948345840069_dp, &
                  0.768969585211864_dp, 0.648139308372852_dp, 0.543433735681048_dp, 0.453521965407942_dp, 0.37693299373884_dp,&
                  0.312169755436428_dp /)

    command='read_data Cu.restart'
    call parse_read(command)
    boundary='sss'
    periodic=.false.
    period=(/.false., .false., .false. /)
    command='potential eam Cu_mishin1.eam.alloy Cu'
    call parse_potential(command)
    command='neighbor 1.5'
    call parse_neighbors(command)
    call pre_calc(1)
    call update_neighbor(0,.false., .true.)
    r_atom(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)

    do i=1,41
        r_atom(:,2) = (/ 0.0_dp, i*0.1_dp, 0.0_dp /)
        call update_force
        if (abs(pe_real(i) - sum(energy_atom(1:2)))>1d-9) then 
            print *, "[FAILED] - Potential energies don't match for reference=",pe_real(i), ' and calculated=', sum(energy_atom)
            call mpi_abort( mpi_comm_world, 1, ierr)
        end if
        fnorm_check=sqrt(force_atom(1,1)*force_atom(1,1) + force_atom(2,1)*force_atom(2,1) + force_atom(3,1)*force_atom(3,1) + &
              force_atom(1,2)*force_atom(1,2) + force_atom(2,2)*force_atom(2,2) + force_atom(3,2)*force_atom(3,2))

        if (abs(fnorm_real(i) - fnorm_check)>1d-9) then 
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
