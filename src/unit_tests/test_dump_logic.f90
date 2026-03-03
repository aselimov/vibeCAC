program test_dump_logic
    use assertions
    use dump_logic

    implicit none

    integer :: failures

    failures = 0

    call test_normalize_dump_every_zero_maps_to_huge(failures)
    call test_normalize_dump_every_positive_unchanged(failures)
    call test_should_write_dump_on_dump_interval(failures)
    call test_should_write_dump_on_force_dump(failures)
    call test_strip_dump_extension_removes_dot_suffix(failures)
    call test_strip_dump_extension_keeps_name_without_dot(failures)

    if (failures > 0) then
        print '(A,I0)', 'test_dump_logic failed with ', failures
        stop 1
    end if

    print '(A)', 'test_dump_logic passed'

contains

    subroutine test_normalize_dump_every_zero_maps_to_huge(failures)
        integer, intent(inout) :: failures
        integer :: normalized

        ! Arrange / Act
        normalized = normalize_dump_every(0)

        ! Assert
        call assert_equal_int(normalized, huge(1), 'dump_every=0 maps to huge(1)', failures)
    end subroutine test_normalize_dump_every_zero_maps_to_huge

    subroutine test_normalize_dump_every_positive_unchanged(failures)
        integer, intent(inout) :: failures
        integer :: normalized

        ! Arrange / Act
        normalized = normalize_dump_every(25)

        ! Assert
        call assert_equal_int(normalized, 25, 'positive dump_every remains unchanged', failures)
    end subroutine test_normalize_dump_every_positive_unchanged

    subroutine test_should_write_dump_on_dump_interval(failures)
        integer, intent(inout) :: failures
        logical :: should_dump

        ! Arrange / Act
        should_dump = should_write_dump(100, 20, .false.)

        ! Assert
        call assert_true(should_dump, 'writes dump when timestep is divisible by dump_every', failures)
    end subroutine test_should_write_dump_on_dump_interval

    subroutine test_should_write_dump_on_force_dump(failures)
        integer, intent(inout) :: failures
        logical :: should_dump

        ! Arrange / Act
        should_dump = should_write_dump(101, 20, .true.)

        ! Assert
        call assert_true(should_dump, 'writes dump when force_dump is true', failures)
    end subroutine test_should_write_dump_on_force_dump

    subroutine test_strip_dump_extension_removes_dot_suffix(failures)
        integer, intent(inout) :: failures
        character(len=32) :: stripped

        ! Arrange / Act
        stripped = strip_dump_extension('traj.out')

        ! Assert
        call assert_equal_str(trim(adjustl(stripped)), 'traj', 'strip_dump_extension removes extension', failures)
    end subroutine test_strip_dump_extension_removes_dot_suffix

    subroutine test_strip_dump_extension_keeps_name_without_dot(failures)
        integer, intent(inout) :: failures
        character(len=32) :: stripped

        ! Arrange / Act
        stripped = strip_dump_extension('trajectory')

        ! Assert
        call assert_equal_str(trim(adjustl(stripped)), 'trajectory', 'strip_dump_extension keeps name without dot', failures)
    end subroutine test_strip_dump_extension_keeps_name_without_dot

end program test_dump_logic
