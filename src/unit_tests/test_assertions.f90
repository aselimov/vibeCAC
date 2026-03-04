program test_assertions
    use parameters, only: wp
    use assertions, only: assert_equal_int_array, assert_close_real_array, assert_close_real_tol, assert_str_contains
    implicit none

    integer :: failures

    failures = 0

    call test_int_array_assertion_accepts_equal_arrays(failures)
    call test_real_array_assertion_honors_abs_and_rel_tolerance(failures)
    call test_real_tolerance_assertion_accepts_relative_match(failures)
    call test_expected_failure_message_pattern_assertion(failures)

    if (failures > 0) then
        print '(A,1X,I0)', 'Unit tests failed:', failures
        stop 1
    end if

    print '(A)', 'All assertion helper tests passed.'

contains

    subroutine test_int_array_assertion_accepts_equal_arrays(failures)
        integer, intent(inout) :: failures
        integer, dimension(4) :: actual, expected

        actual = [1, 2, 3, 4]
        expected = [1, 2, 3, 4]

        call assert_equal_int_array(actual, expected, 'integer arrays with same values should match', failures)
    end subroutine test_int_array_assertion_accepts_equal_arrays

    subroutine test_real_array_assertion_honors_abs_and_rel_tolerance(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: actual, expected

        expected = [1.0_wp, 2.0_wp, 1.0e6_wp]
        actual = [1.0_wp + 1.0e-9_wp, 2.0_wp - 1.0e-9_wp, 1.0e6_wp + 5.0e-4_wp]

        call assert_close_real_array(actual, expected, 1.0e-8_wp, 1.0e-9_wp, &
            'real array comparison should support combined abs/rel tolerance', failures)
    end subroutine test_real_array_assertion_honors_abs_and_rel_tolerance

    subroutine test_real_tolerance_assertion_accepts_relative_match(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: actual, expected

        expected = 1.0e8_wp
        actual = expected + 5.0e-2_wp

        call assert_close_real_tol(actual, expected, 1.0e-6_wp, 1.0e-9_wp, &
            'large magnitude values should use relative tolerance', failures)
    end subroutine test_real_tolerance_assertion_accepts_relative_match

    subroutine test_expected_failure_message_pattern_assertion(failures)
        integer, intent(inout) :: failures
        character(len=128) :: message_text

        message_text = 'ERROR: timestep must be > 0'
        call assert_str_contains(trim(message_text), 'must be > 0', &
            'expected failure message pattern should be matched', failures)
    end subroutine test_expected_failure_message_pattern_assertion

end program test_assertions
