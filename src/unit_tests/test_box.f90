program test_box
    use parameters, only: wp
    use box, only: box_bd, box_length, period, cross_pb, restore_pb
    use assertions, only: assert_true, assert_false, assert_equal_int_array, assert_close_real_array
    use fixtures, only: begin_suite, end_suite
    implicit none

    integer :: failures

    call begin_suite(failures)

    call test_cross_pb_wraps_lower_and_upper_periodic(failures)
    call test_cross_pb_keeps_exact_boundaries(failures)
    call test_cross_pb_regression_does_not_wrap_exact_upper_bound(failures)
    call test_cross_pb_ignores_nonperiodic_out_of_bounds(failures)
    call test_restore_pb_applies_only_periodic_moves(failures)
    call test_restore_pb_nonperiodic_moves_do_not_set_flag(failures)

    call end_suite(failures, 'test_box')

contains

    subroutine set_test_box(lower, upper)
        real(kind=wp), intent(in) :: lower, upper

        box_bd = [lower, upper, lower, upper, lower, upper]
        box_length = [upper - lower, upper - lower, upper - lower]
    end subroutine set_test_box

    subroutine test_cross_pb_wraps_lower_and_upper_periodic(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: info(3)

        ! Arrange
        call set_test_box(0.0_wp, 10.0_wp)
        period = [.true., .true., .true.]
        r = [-0.5_wp, 10.5_wp, 5.0_wp]

        ! Act
        call cross_pb(r, info)

        ! Assert
        call assert_close_real_array(r, [9.5_wp, 0.5_wp, 5.0_wp], 1.0e-12_wp, 0.0_wp, &
            'cross_pb should wrap coordinates that cross periodic boundaries', failures)
        call assert_equal_int_array(info, [1, -1, 0], 'cross_pb should report boundary crossing direction', failures)
    end subroutine test_cross_pb_wraps_lower_and_upper_periodic

    subroutine test_cross_pb_keeps_exact_boundaries(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: info(3)

        ! Arrange
        call set_test_box(0.0_wp, 10.0_wp)
        period = [.true., .true., .true.]
        r = [0.0_wp, 10.0_wp, 7.5_wp]

        ! Act
        call cross_pb(r, info)

        ! Assert
        call assert_close_real_array(r, [0.0_wp, 10.0_wp, 7.5_wp], 1.0e-12_wp, 0.0_wp, &
            'cross_pb should not move coordinates exactly on box bounds', failures)
        call assert_equal_int_array(info, [0, 0, 0], 'cross_pb should not flag exact-boundary coordinates', failures)
    end subroutine test_cross_pb_keeps_exact_boundaries

    subroutine test_cross_pb_regression_does_not_wrap_exact_upper_bound(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: info(3)

        ! Regression: BUGREF(local:20260303-cross-pb-upper-bound) exact upper-bound coordinates previously got wrapped due to inclusive bound checks.
        call set_test_box(0.0_wp, 10.0_wp)
        period = [.false., .true., .false.]
        r = [4.0_wp, 10.0_wp, 6.0_wp]

        call cross_pb(r, info)

        call assert_close_real_array(r, [4.0_wp, 10.0_wp, 6.0_wp], 1.0e-12_wp, 0.0_wp, &
            'cross_pb regression: exact upper bound remains unchanged', failures)
        call assert_equal_int_array(info, [0, 0, 0], 'cross_pb regression: exact upper bound is not flagged', failures)
    end subroutine test_cross_pb_regression_does_not_wrap_exact_upper_bound

    subroutine test_cross_pb_ignores_nonperiodic_out_of_bounds(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: info(3)

        ! Arrange
        call set_test_box(0.0_wp, 10.0_wp)
        period = [.false., .true., .false.]
        r = [-1.0_wp, 11.0_wp, 12.0_wp]

        ! Act
        call cross_pb(r, info)

        ! Assert
        call assert_close_real_array(r, [-1.0_wp, 1.0_wp, 12.0_wp], 1.0e-12_wp, 0.0_wp, &
            'cross_pb should only wrap periodic dimensions', failures)
        call assert_equal_int_array(info, [0, -1, 0], 'cross_pb should only report periodic-dimension crossings', failures)
    end subroutine test_cross_pb_ignores_nonperiodic_out_of_bounds

    subroutine test_restore_pb_applies_only_periodic_moves(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: pb_move(3)
        logical :: ipb

        ! Arrange
        box_length = [10.0_wp, 20.0_wp, 30.0_wp]
        period = [.true., .false., .true.]
        r = [5.0_wp, 7.0_wp, 9.0_wp]
        pb_move = [1, 1, -2]

        ! Act
        call restore_pb(r, pb_move, ipb)

        ! Assert
        call assert_close_real_array(r, [-5.0_wp, 7.0_wp, 69.0_wp], 1.0e-12_wp, 0.0_wp, &
            'restore_pb should apply periodic boundary displacements only where periodic', failures)
        call assert_true(ipb, 'restore_pb should mark ipb when periodic displacement is applied', failures)
    end subroutine test_restore_pb_applies_only_periodic_moves

    subroutine test_restore_pb_nonperiodic_moves_do_not_set_flag(failures)
        integer, intent(inout) :: failures
        real(kind=wp) :: r(3)
        integer :: pb_move(3)
        logical :: ipb

        ! Arrange
        box_length = [10.0_wp, 20.0_wp, 30.0_wp]
        period = [.false., .true., .false.]
        r = [1.0_wp, 2.0_wp, 3.0_wp]
        pb_move = [2, 0, -1]

        ! Act
        call restore_pb(r, pb_move, ipb)

        ! Assert
        call assert_close_real_array(r, [1.0_wp, 2.0_wp, 3.0_wp], 1.0e-12_wp, 0.0_wp, &
            'restore_pb should ignore nonperiodic dimensions even with nonzero displacement', failures)
        call assert_false(ipb, 'restore_pb should keep ipb false when no periodic displacement is applied', failures)
    end subroutine test_restore_pb_nonperiodic_moves_do_not_set_flag

end program test_box
