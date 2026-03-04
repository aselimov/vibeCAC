program test_time
    use time, only: is_valid_timestep, read_timestep_value, validate_timestep
    use parameters, only: wp
    use assertions, only: assert_true, assert_false, assert_equal_int, assert_close_real
    use fixtures, only: begin_suite, end_suite
    implicit none

    integer :: failures
    logical :: error_called

    call begin_suite(failures)
    error_called = .false.

    call test_is_valid_timestep_accepts_positive(failures)
    call test_is_valid_timestep_rejects_zero(failures)
    call test_is_valid_timestep_rejects_negative(failures)
    call test_read_timestep_value_parses_numeric_argument(failures)
    call test_read_timestep_value_rejects_nonnumeric_argument(failures)
    call test_validate_timestep_calls_handler_on_invalid_value(failures)
    call test_validate_timestep_does_not_call_handler_on_valid_value(failures)

    call end_suite(failures, 'test_time')

contains

    subroutine capture_error(message)
        character(len=*), intent(in) :: message
        if (len_trim(message) < 0) stop 1
        error_called = .true.
    end subroutine capture_error

    subroutine test_is_valid_timestep_accepts_positive(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_true(is_valid_timestep(1.0e-6_wp), 'positive timestep should be accepted', failures)
    end subroutine test_is_valid_timestep_accepts_positive

    subroutine test_is_valid_timestep_rejects_zero(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_false(is_valid_timestep(0.0_wp), 'zero timestep should be rejected', failures)
    end subroutine test_is_valid_timestep_rejects_zero

    subroutine test_is_valid_timestep_rejects_negative(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_false(is_valid_timestep(-1.0e-6_wp), 'negative timestep should be rejected', failures)
    end subroutine test_is_valid_timestep_rejects_negative

    subroutine test_read_timestep_value_parses_numeric_argument(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: parsed_dt
        integer :: iospara

        ! Arrange / Act
        call read_timestep_value('timestep 0.0125', parsed_dt, iospara)

        ! Assert
        call assert_equal_int(iospara, 0, 'timestep parser should read valid command', failures)
        call assert_close_real(parsed_dt, 0.0125_wp, 1.0e-12_wp, 'timestep parser should return parsed value', failures)
    end subroutine test_read_timestep_value_parses_numeric_argument

    subroutine test_read_timestep_value_rejects_nonnumeric_argument(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: parsed_dt
        integer :: iospara

        ! Arrange / Act
        call read_timestep_value('timestep abc', parsed_dt, iospara)

        ! Assert
        call assert_true(iospara > 0, 'timestep parser should fail when numeric argument is non-numeric', failures)
    end subroutine test_read_timestep_value_rejects_nonnumeric_argument

    subroutine test_validate_timestep_calls_handler_on_invalid_value(failures)
        integer, intent(inout) :: failures

        ! Arrange
        error_called = .false.

        ! Act
        call validate_timestep(0.0_wp, capture_error)

        ! Assert
        call assert_true(error_called, 'invalid timestep should trigger injected error handler', failures)
    end subroutine test_validate_timestep_calls_handler_on_invalid_value

    subroutine test_validate_timestep_does_not_call_handler_on_valid_value(failures)
        integer, intent(inout) :: failures

        ! Arrange
        error_called = .false.

        ! Act
        call validate_timestep(1.0e-6_wp, capture_error)

        ! Assert
        call assert_false(error_called, 'valid timestep should not trigger injected error handler', failures)
    end subroutine test_validate_timestep_does_not_call_handler_on_valid_value

end program test_time
