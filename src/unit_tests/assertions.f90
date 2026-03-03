module assertions
    use parameters, only: wp
    implicit none

    private
    public :: assert_true, assert_false, assert_equal_int, assert_equal_str, assert_close_real

contains

    subroutine assert_true(condition, message, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            print '(A)', 'FAIL: ' // trim(message)
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
            print '(A,1X,I0,1X,A,1X,I0)', 'FAIL: ' // trim(message) // ' actual=', actual, 'expected=', expected
        end if
    end subroutine assert_equal_int

    subroutine assert_equal_str(actual, expected, message, failures)
        character(len=*), intent(in) :: actual, expected
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (trim(actual) /= trim(expected)) then
            failures = failures + 1
            print '(A)', 'FAIL: ' // trim(message) // ' actual="' // trim(actual) // '" expected="' // trim(expected) // '"'
        end if
    end subroutine assert_equal_str

    subroutine assert_close_real(actual, expected, eps, message, failures)
        real(kind=wp), intent(in) :: actual, expected, eps
        character(len=*), intent(in) :: message
        integer, intent(inout) :: failures

        if (abs(actual - expected) > eps) then
            failures = failures + 1
            print '(A,1X,ES12.5,1X,A,1X,ES12.5)', 'FAIL: ' // trim(message) // ' actual=', actual, 'expected=', expected
        end if
    end subroutine assert_close_real

end module assertions
