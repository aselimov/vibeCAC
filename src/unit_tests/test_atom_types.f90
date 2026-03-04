program test_atom_types
    use atom_types, only: init_atom_types, set_atom_types, parse_types, natom_types, atom_names, atom_types_set, &
                          types_to_pot_type, max_atom_types
    use assertions, only: assert_true, assert_false, assert_equal_int, assert_equal_str
    use fixtures, only: begin_suite, end_suite
    implicit none

    integer :: failures

    call begin_suite(failures)

    call test_init_atom_types_resets_flags_and_maps(failures)
    call test_set_atom_types_assigns_count_and_names(failures)
    call test_set_atom_types_accepts_max_atom_types(failures)
    call test_parse_types_reads_names_from_command(failures)
    call test_parse_types_accepts_missing_type_names(failures)

    call end_suite(failures, 'test_atom_types')

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

    subroutine test_set_atom_types_accepts_max_atom_types(failures)
        integer, intent(inout) :: failures
        character(len=2) :: names(max_atom_types)

        ! Arrange
        names = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'B0']
        call init_atom_types

        ! Act
        call set_atom_types(max_atom_types, names)

        ! Assert
        call assert_equal_int(natom_types, max_atom_types, 'set_atom_types supports max_atom_types boundary', failures)
        call assert_equal_str(atom_names(1), 'A1', 'set_atom_types keeps first boundary name', failures)
        call assert_equal_str(atom_names(max_atom_types), 'B0', 'set_atom_types keeps last boundary name', failures)
        call assert_true(atom_types_set, 'set_atom_types marks initialized at max boundary', failures)
    end subroutine test_set_atom_types_accepts_max_atom_types

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

    subroutine test_parse_types_accepts_missing_type_names(failures)
        integer, intent(inout) :: failures

        ! Arrange
        call init_atom_types

        ! Act
        call parse_types('types')

        ! Assert
        call assert_equal_int(natom_types, 0, 'parse_types keeps zero type count for missing names', failures)
        call assert_true(atom_types_set, 'parse_types still marks initialized for zero-type input', failures)
    end subroutine test_parse_types_accepts_missing_type_names

end program test_atom_types
