module fixtures
    implicit none

    public :: begin_suite, end_suite, make_required_flags

contains

    subroutine begin_suite(failures)
        integer, intent(out) :: failures

        failures = 0
    end subroutine begin_suite

    subroutine end_suite(failures, suite_name)
        integer, intent(in) :: failures
        character(len=*), intent(in) :: suite_name

        if (failures > 0) then
            print '(A,1X,I0)', trim(suite_name)//' failed with', failures
            stop 1
        end if

        print '(A)', trim(suite_name)//' passed'
    end subroutine end_suite

    function make_required_flags(first_present, second_present) result(req_flags)
        logical, intent(in) :: first_present, second_present
        logical :: req_flags(2)

        req_flags = [first_present, second_present]
    end function make_required_flags

end module fixtures
