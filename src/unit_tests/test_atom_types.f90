program test_atom_types
    use atom_types, only: init_atom_types, set_atom_types, parse_types, natom_types, atom_names, atom_types_set, &
                          types_to_pot_type
    use assertions, only: assert_true, assert_false, assert_equal_int, assert_equal_str
    implicit none

    integer :: failures

    failures = 0

    call test_init_atom_types_resets_flags_and_maps(failures)
    call test_set_atom_types_assigns_count_and_names(failures)
    call test_parse_types_reads_names_from_command(failures)

    if (failures > 0) then
        print '(A,1X,I0)', 'Unit tests failed:', failures
        stop 1
    end if

    print '(A)', 'All atom_types unit tests passed.'

contains

    subroutine test_init_atom_types_resets_flags_and_maps(failures)
        integer, intent(inout) :: failures

        ! Arrange
        atom_types_set = .true.
        types_to_pot_type = 9

        ! Act
        call init_atom_types

        ! Assert
        call assert_false(atom_types_set, 'init_atom_types clears atom_types_set flag', failures)
        call assert_true(all(types_to_pot_type == 0), 'init_atom_types zeroes type-to-potential mapping', failures)
    end subroutine test_init_atom_types_resets_flags_and_maps

    subroutine test_set_atom_types_assigns_count_and_names(failures)
        integer, intent(inout) :: failures
        character(len=2) :: names(3)

        ! Arrange
        names = ['Cu', 'Ni', 'Al']
        call init_atom_types

        ! Act
        call set_atom_types(3, names)

        ! Assert
        call assert_equal_int(natom_types, 3, 'set_atom_types stores number of atom types', failures)
        call assert_equal_str(atom_names(1), 'Cu', 'set_atom_types stores first name', failures)
        call assert_equal_str(atom_names(2), 'Ni', 'set_atom_types stores second name', failures)
        call assert_equal_str(atom_names(3), 'Al', 'set_atom_types stores third name', failures)
        call assert_true(atom_types_set, 'set_atom_types marks atom types initialized', failures)
    end subroutine test_set_atom_types_assigns_count_and_names

    subroutine test_parse_types_reads_names_from_command(failures)
        integer, intent(inout) :: failures

        ! Arrange
        call init_atom_types

        ! Act
        call parse_types('types Fe C O')

        ! Assert
        call assert_equal_int(natom_types, 3, 'parse_types infers type count from tokens', failures)
        call assert_equal_str(atom_names(1), 'Fe', 'parse_types reads first atom name', failures)
        call assert_equal_str(atom_names(2), 'C', 'parse_types reads second atom name', failures)
        call assert_equal_str(atom_names(3), 'O', 'parse_types reads third atom name', failures)
        call assert_true(atom_types_set, 'parse_types marks atom types initialized', failures)
    end subroutine test_parse_types_reads_names_from_command

end program test_atom_types
