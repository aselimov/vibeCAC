module input_parser
    !This code is in charge of reading the code and running functions that mattter

    use mpi 
    use comms
    use parameters
    use potential
    use eam
    use read_data
    use neighbors
    use logger
    use dump
    use group
    use minimize
    use force_mod
    use displace
    use set
    use deform
    use modify
    use langevin
    use input_validation

    implicit none

    logical, private :: initialized

    public 
    contains

    subroutine read_input
        !Actually read the input, every processor should have access to the input arguments
        !so we should be able to loop over these arguments
        
        character(len=read_len) :: line, label, req_commands(2), msg, loop_commands(100)
        integer :: iosline, iospara, i, loop_times, nloopcom, nloop, j
        logical :: req_flags(2), first_flag(2), read_loop_flag, run_loop_flag

        initialized = .false.
        req_flags(:) = .false.
        read_loop_flag = .false.
        run_loop_flag = .false.
        nloopcom = 0
        first_flag = .true.
        req_commands(1) =  'read'
        req_commands(2) = 'potential'
        iosline = 0
        do while (iosline == 0)
            if(run_loop_flag) then 
                if(j > nloopcom) then 
                    j = 1
                    nloop = nloop + 1  
                    if(nloop > loop_times) then 
                        run_loop_flag = .false.
                        nloopcom = 0
                        cycle
                    end if
                end if

                line = loop_commands(j)
                j = j+1
            else
                if(rank == root) then 
                    line=''
                    read(*, '(a)', iostat = iosline) line
                    if(iosline > 0) then
                        print *, 'Error: Wrong reading input file line ', line, &
                                 ' because of', iosline
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if
            
                    !scan(line, '#' /= 1) is used to skip all lines starting with #
                    if((scan(line, '#') /= 1).and.(line /= '')) then

                       !This communication needs to be done over mpi_comm_world to make sure every processor
                        !is running the commands
                        call mpi_bcast(line, read_len, mpi_character, root, mpi_comm_world, ierr)
                    else 
                        cycle
                    end if
                else
                    call mpi_bcast(line, read_len, mpi_character, root, mpi_comm_world, ierr) 
                    !Check for the exit input loop code 
                    if (line == "exit input loop") exit
                end if
            end if

            !Get the command 
            read(line, *, iostat = iospara) label
            if(iospara > 0) then
                print *, 'Error: Wrong reading input file label ', label, &
                         ' because of', iospara
                call mpi_abort(mpi_comm_world, 1, ierr)
            end if

            !Log the command 
            call log_msg(line,0)

            !run the command
            if(read_loop_flag) then 
                if(label == 'endloop') then 
                    read_loop_flag = .false.    
                    run_loop_flag= .true.
                    nloop = 1
                    j = 1
                else 
                    nloopcom = nloopcom+1
                    loop_commands(nloopcom) = trim(adjustl(line))
                end if
            else
                select case(label)
                case('loop')
                    !Get the number of times to loop
                    read(line, *) label, loop_times
                    read_loop_flag = .true.
                    nloopcom = 0

                case('read_data')
                    call parse_read(line)
                    initialized=.false.
                    need_vel=.false.
                    need_force_pre=.false.
                    req_flags(1) = .true.
                case('potential')
                    call parse_potential(line)
                    req_flags(2) = .true.
                case('neighbor')
                    call parse_neighbors(line)
                case('write_out')
                    call parse_write_out(line)
                case('dump')
                    call parse_dump(line)
                case('undump')
                    call parse_undump(line)
                case('timestep')
                    call parse_timestep(line)
                case('modify')
                    call parse_modify(line)
                case('group')
                    call parse_group(line)
                case('mass')
                    call parse_mass(line)
                case('debug')
                    call parse_debug(line)
                case('langevin')
                    call parse_langevin(line)
                case('setforce')
                    if(first_flag(1)) then 
                        call init_set_force
                        first_flag(1) = .false.
                    end if
                    call parse_set_force(line)
                case('addforce')
                    if(first_flag(2)) then 
                        call init_add_force
                        first_flag(2) = .false.
                    end if
                    call parse_add_force(line)
                case('displace')
                    call displace_points(line)
                case('ramp')
                    call ramp_displace(line)
                case('boundary')
                    call parse_boundary(line)

                    !Now regenerate the grid comm
                    call mpi_cart_create(world, 3, num_pro, period, .false., grid_comm, ierr)
                    call mpi_comm_rank(grid_comm, grank, ierr)
                    call mpi_cart_coords(grid_comm, grank, 3, grid_coords, ierr) 

                case('thermo')
                    call parse_thermo(line)
                case('types')
                    call parse_types(line)
                case('min_style')
                    call parse_min_style(line)
                case('minimize')
                    select case(min_style)
                    case(1)
                        !If we are using fire than we need all of the info
                        call pre_calc(2)
                    case(2)
                        !If we are using cg than we only need the virial
                        call pre_calc(3)
                    end select
                    call parse_minimize(line)
                case('run')
                    i = first_missing_required_command(req_flags)
                    if (i > 0) then
                        print *, "Error: must call ", req_commands(i), " before run command"
                        call mpi_abort(mpi_comm_world, 1, ierr)
                    end if
                    !Prepare for calculation
                    call pre_calc(1)
                    call parse_run(line)
                case('dynamics')
                    call parse_dynamics(line)
                case('set')
                    call parse_set(line)
                case('temp')
                    call parse_temp(line)
                case('press')
                    call parse_berendsen(line)
                case('unpress')
                    pflag=.false.
                case('deform')
                    call parse_deform(line)
                case('thermo_style')
                    call parse_thermo_style(line)
                case('write_group')
                    call write_group(line)
                case default
                    write(msg, *) 'Input parameter label ', trim(adjustl(label)), &
                                  ' in line ', trim(adjustl(line)), ' is not accepted'
                    call misc_error(msg)
                end select
            end if
        end do


        !Because root will be the only processor that has access to stdin we have to broadcast a code to all the 
        !other processors for them to exit the input loop
        if(rank == root) then 
            line = "exit input loop"
            call mpi_bcast(line, read_len, mpi_character, root, mpi_comm_world, ierr) 
        end if

        return
    end subroutine read_input

    subroutine pre_calc(runtype)
        !This is the code which distributes the model and builds neighbor lists before calculation

        !Runtype dictates what the command is being called for calculation. The option for runtype changes how flags are set
        integer, intent(in) :: runtype

        !Set flags and allocate variables
        need_virial = .true.
        select case(runtype)
        case(1)
            if(.not.need_vel) then 
                need_vel = .true.
                call alloc_velocity_arrays 
            end if

            if(.not.need_force_pre) then 
                need_force_pre = .true.
                call alloc_pre_array
            end if
        case(2)
            if(.not.need_vel) then 
                need_vel = .true.
                call alloc_velocity_arrays 
            end if
            if(need_force_pre) then 
                need_force_pre=.false.
                if(atom_num > 0) deallocate(force_atom_pre)
                if(ele_num>0) deallocate(force_eq_pre)
            end if
        case(3)
            if(need_vel) then 
                need_vel = .false.
                call dealloc_velocity_arrays 
            end if
        end select
        
        if (.not. initialized) then 
            !Set up the potential map

            !Allocate necessary force and potential arrays
            call alloc_force_arrays
            call alloc_potential_arrays

            initialized = .true.

        end if

    end subroutine pre_calc
    
end module input_parser

