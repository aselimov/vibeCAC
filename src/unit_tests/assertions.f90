module assertions
    use parameters, only: wp
    implicit none

    private
    public :: assert_true, assert_false, assert_equal_int, assert_equal_str, assert_close_real
    public :: assert_equal_int_array, assert_close_real_array, assert_close_real_tol
    public :: assert_str_contains

contains

    subroutine assert_true(condition, message, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            print '(A)', 'FAIL: ' // trim(message) // ' expected=.true. actual=.false.'
        end if
    end subroutine assert_true

    subroutine assert_false(condition, message, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        call assert_true(.not. condition, message, failures)
    end subroutine assert_false

    subroutine assert_equal_int(actual, expected, message, failures)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (actual /= expected) then
            failures = failures + 1
            print '(A,1X,I0,1X,A,1X,I0)', 'FAIL: ' // trim(message) // ' expected=', expected, 'actual=', actual
        end if
    end subroutine assert_equal_int

    subroutine assert_equal_str(actual, expected, message, failures)
        character(len=*), intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (trim(actual) /= trim(expected)) then
            failures = failures + 1
            print '(A)', 'FAIL: ' // trim(message) // ' expected="' // trim(expected) // '" actual="' // trim(actual) // '"'
        end if
    end subroutine assert_equal_str

    subroutine assert_close_real(actual, expected, eps, message, failures)
        real(kind=wp), intent(in) :: actual, expected, eps
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        call assert_close_real_tol(actual, expected, eps, 0.0_wp, message, failures)
    end subroutine assert_close_real

    subroutine assert_close_real_tol(actual, expected, abs_tol, rel_tol, message, failures)
        real(kind=wp), intent(in) :: actual, expected, abs_tol, rel_tol
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures
        real(kind=wp) :: diff, limit

        diff = abs(actual - expected)
        limit = max(abs_tol, rel_tol * max(abs(actual), abs(expected)))

        if (diff > limit) then
            failures = failures + 1
            print '(A,1X,ES12.5,1X,A,1X,ES12.5,1X,A,1X,ES12.5,1X,A,1X,ES12.5)', &
                'FAIL: ' // trim(message) // ' expected=', expected, 'actual=', actual, 'abs_diff=', diff, 'allowed=', limit
        end if
    end subroutine assert_close_real_tol

    subroutine assert_equal_int_array(actual, expected, message, failures)
        integer, dimension(:), intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures
        integer :: i

        if (size(actual) /= size(expected)) then
            failures = failures + 1
            print '(A,1X,I0,1X,A,1X,I0)', 'FAIL: ' // trim(message) // ' expected_size=', size(expected), 'actual_size=', size(actual)
            return
        end if

        do i = 1, size(actual)
            if (actual(i) /= expected(i)) then
                failures = failures + 1
                print '(A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)', &
                    'FAIL: ' // trim(message) // ' idx=', i, 'expected=', expected(i), 'actual=', actual(i)
                return
            end if
        end do
    end subroutine assert_equal_int_array

    subroutine assert_close_real_array(actual, expected, abs_tol, rel_tol, message, failures)
        real(kind=wp), dimension(:), intent(in) :: actual, expected
        real(kind=wp), intent(in) :: abs_tol, rel_tol
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures
        integer :: i
        real(kind=wp) :: diff, limit

        if (size(actual) /= size(expected)) then
            failures = failures + 1
            print '(A,1X,I0,1X,A,1X,I0)', 'FAIL: ' // trim(message) // ' expected_size=', size(expected), 'actual_size=', size(actual)
            return
        end if

        do i = 1, size(actual)
            diff = abs(actual(i) - expected(i))
            limit = max(abs_tol, rel_tol * max(abs(actual(i)), abs(expected(i))))
            if (diff > limit) then
                failures = failures + 1
                print '(A,1X,I0,1X,A,1X,ES12.5,1X,A,1X,ES12.5,1X,A,1X,ES12.5,1X,A,1X,ES12.5)', &
                    'FAIL: ' // trim(message) // ' idx=', i, 'expected=', expected(i), 'actual=', actual(i), 'abs_diff=', diff, 'allowed=', limit
                return
            end if
        end do
    end subroutine assert_close_real_array

    subroutine assert_str_contains(actual, expected_substring, message, failures)
        character(len=*), intent(in) :: actual, expected_substring
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (index(trim(actual), trim(expected_substring)) == 0) then
            failures = failures + 1
            print '(A)', 'FAIL: ' // trim(message) // ' expected_pattern="' // trim(expected_substring) // '" actual="' // trim(actual) // '"'
        end if
    end subroutine assert_str_contains

end module assertions
