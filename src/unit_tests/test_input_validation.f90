program test_input_validation
    use assertions
    use input_validation
    use fixtures, only: begin_suite, end_suite, make_required_flags

    implicit none

    integer :: failures

    call begin_suite(failures)

    call test_first_missing_required_command_none_missing(failures)
    call test_first_missing_required_command_first_missing(failures)
    call test_first_missing_required_command_second_missing(failures)

    call end_suite(failures, 'test_input_validation')

contains

    subroutine test_first_missing_required_command_none_missing(failures)
        integer, intent(inout) :: failures
        logical :: req_flags(2)
        integer :: missing_index

        ! Arrange
        req_flags = make_required_flags(.true., .true.)

        ! Act
        missing_index = first_missing_required_command(req_flags)

        ! Assert
        call assert_equal_int(missing_index, 0, 'no missing required commands', failures)
    end subroutine test_first_missing_required_command_none_missing

    subroutine test_first_missing_required_command_first_missing(failures)
        integer, intent(inout) :: failures
        logical :: req_flags(2)
        integer :: missing_index

        ! Arrange
        req_flags = make_required_flags(.false., .true.)

        ! Act
        missing_index = first_missing_required_command(req_flags)

        ! Assert
        call assert_equal_int(missing_index, 1, 'first required command detected as missing', failures)
    end subroutine test_first_missing_required_command_first_missing

    subroutine test_first_missing_required_command_second_missing(failures)
        integer, intent(inout) :: failures
        logical :: req_flags(2)
        integer :: missing_index

        ! Arrange
        req_flags = make_required_flags(.true., .false.)

        ! Act
        missing_index = first_missing_required_command(req_flags)

        ! Assert
        call assert_equal_int(missing_index, 2, 'second required command detected as missing', failures)
    end subroutine test_first_missing_required_command_second_missing

end program test_input_validation
