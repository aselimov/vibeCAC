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
    integer :: iep, ibasis, jep, iatom, iatomap, intpo_neighbor(125, 14), nei, ia, ja, ie
    
    logical :: test_results(2)

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
    !Results are calculated from lammps simulations. These are cohesive energy calculations meaning that a pure block of Cu is
    !generated with a lattice constant of 3.615 \AA and the resulting energy and force_norm are copied here 
    
    !Now read in the data, we pass commands using strings to the necessary functions
    !Start with the pure atomistic case 
    test_results=.false.


    command='read_data Cu.restart'
    call parse_read(command)
    command='potential eam Cu_mishin1.eam.alloy Cu'
    call parse_potential(command)
    command='neighbor 0.0'
    call parse_neighbors(command)
    command='run 0'
    call pre_calc(1)
    call parse_run(command)

    !First check the neighbor lists to make sure that the neighbor list 
    intpo_neighbor=0
    do ia=1, atom_num_l
        do nei=1, n_at_cg(ia, 1)
            ja= at_nei_cg(nei, ia, 1)
            if( ja <= atomap_num_l) then
                if (atomap_to_intpo(1,ja)/=0) then
                        iep=atomap_to_intpo(1,ja)
                        ie=atomap_to_intpo(2,ja)
                        intpo_neighbor(iep, ie) = intpo_neighbor(iep, ie) + 1
                end if
            end if
        end do
    end do


    do ie= 1, ele_num
        do iep =1, 125
            iatom=atom_intpo(iep, size_to_shape(size_ele(ie)), itype(ie)) 
            iatomap=cg_atomap(iatom, ie)
            intpo_neighbor(iep, ie)= intpo_neighbor(iep, ie) + n_cg_cg( iep, ie, 1) +  n_cg_at(iep, ie, 1)
            if (intpo_neighbor(iep, ie) /= 54) then
                    print *, '[FAILED] - Intpo should have 54 neighbors not ', intpo_neighbor(iep, ie), ' for iep=',iep
                    call mpi_abort(1, mpi_comm_world, ierr)
            end if
        end do
    end do

    print *, '[OK] - Passed neighbor list check for integration points'

    !Now check node energyes

    !close log
    call log_time
    call close_log

    call mpi_bcast(test_results, int(size(test_results)), mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_finalize(ierr)

end program cg_main
