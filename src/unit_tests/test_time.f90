program test_time
    use time, only: is_valid_timestep
    use parameters, only: wp
    use assertions, only: assert_true, assert_false
    implicit none

    integer :: failures

    failures = 0

    call test_is_valid_timestep_accepts_positive(failures)
    call test_is_valid_timestep_rejects_zero(failures)
    call test_is_valid_timestep_rejects_negative(failures)

    if (failures > 0) then
        print '(A,1X,I0)', 'Unit tests failed:', failures
        stop 1
    end if

    print '(A)', 'All time unit tests passed.'

contains

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

end program test_time
