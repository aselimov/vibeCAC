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
    !Results are calculated from lammps simulations. This just calculates the pair potential and force for a variety of distances

    pe_real=(/ 202.839586065967_dp, 152.895955213894_dp, 115.240959256127_dp, 86.8521085273478_dp, 65.4501706918583_dp, &
               49.3164166389302_dp, 37.1547998765048_dp, 27.9880220245331_dp, 21.079153090186_dp, 15.8725235262972_dp, &
               11.9491498034614_dp, 8.99312017480803_dp, 6.76624585122577_dp, 5.08894534884884_dp, 3.82582942327178_dp, &
               2.87483081390638_dp, 2.15900718894046_dp, 1.62035998306567_dp, 1.21517343336307_dp, 0.910499997309645_dp, &
               0.681510249872969_dp, 0.509494671390344_dp, 0.38035701038665_dp, 0.283478325697607_dp, 0.21086054011041_dp, &
               0.156480756210957_dp, 0.115804491354539_dp, 0.0854187379355531_dp, 0.062755369492367_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
               0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)

    fnorm_real=(/ 810.730621055971_dp, 611.260838995818_dp, 460.850883135648_dp, 347.436684620931_dp, 261.920496162452_dp, &
                  197.441569826202_dp, 148.826146358199_dp, 112.172671092505_dp, 84.5389899850447_dp, 63.7064531265885_dp, &
                  48.0020173361934_dp, 36.1640881789499_dp, 27.2413475809743_dp, 20.5164571194941_dp, 15.4485209537088_dp, &
                  11.6296960420211_dp, 8.75247128279655_dp, 6.58499240845068_dp, 4.95245440033986_dp, 3.72306956684669_dp, &
                  2.79748622265997_dp, 2.10080952279828_dp, 1.57658461172363_dp, 1.18225956633655_dp, 0.885764251347544_dp, &
                  0.662930675958203_dp, 0.495547913477021_dp, 0.369895528800911_dp, 0.275637831293647_dp, 0.0_dp, 0.0_dp, 0.0_dp,&
                  0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /)

    command='read_data Cu.restart'
    call parse_read(command)
    boundary='sss'
    periodic=.false.
    period=(/.false., .false., .false. /)
    command='types Av He'
    call parse_types(command)
    command='potential morse 1 2 0.000529628 1.410836 4.657221 3.0'
    call parse_potential(command)
    command='potential morse 1 1 0.000529628 1.410836 4.657221 3.0'
    call parse_potential(command)
    command='potential morse 2 2 0.000529628 1.410836 4.657221 3.0'
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
