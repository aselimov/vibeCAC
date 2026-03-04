program test_math
    use parameters, only: wp, dp
    use math, only: identity_mat, cross_product, triple_product, in_block_bd, in_cutoff, neighbor_dis_free, dumb_interp
    use assertions, only: assert_true, assert_false, assert_equal_int, assert_close_real
    use fixtures, only: begin_suite, end_suite
    implicit none

    integer :: failures

    call begin_suite(failures)

    call test_identity_matrix(failures)
    call test_cross_product_basis(failures)
    call test_triple_product_right_handed(failures)
    call test_triple_product_zero_for_coplanar(failures)
    call test_in_block_bd_interior_and_upper_boundary(failures)
    call test_in_block_bd_optional_dimension_mask(failures)
    call test_in_cutoff_threshold_behavior(failures)
    call test_neighbor_dis_free_components_and_norm(failures)
    call test_dumb_interp_midpoint(failures)

    call end_suite(failures, 'test_math')

contains

    subroutine test_identity_matrix(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3,3) :: id

        ! Arrange
        id = identity_mat(3)

        ! Act / Assert
        call assert_equal_int(count(abs(id - reshape([ &
            1.0_wp, 0.0_wp, 0.0_wp, &
            0.0_wp, 1.0_wp, 0.0_wp, &
            0.0_wp, 0.0_wp, 1.0_wp], [3,3])) < 1.0e-12_wp), 9, &
            'test_identity_matrix should return a 3x3 identity matrix', failures)
    end subroutine test_identity_matrix

    subroutine test_cross_product_basis(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: ex, ey, cross

        ! Arrange
        ex = [1.0_wp, 0.0_wp, 0.0_wp]
        ey = [0.0_wp, 1.0_wp, 0.0_wp]

        ! Act
        cross = cross_product(ex, ey)

        ! Assert
        call assert_close_real(cross(1), 0.0_wp, 1.0e-12_wp, 'cross basis x component', failures)
        call assert_close_real(cross(2), 0.0_wp, 1.0e-12_wp, 'cross basis y component', failures)
        call assert_close_real(cross(3), 1.0_wp, 1.0e-12_wp, 'cross basis z component', failures)
    end subroutine test_cross_product_basis

    subroutine test_triple_product_right_handed(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: ex, ey, ez
        real(kind=wp) :: tp

        ! Arrange
        ex = [1.0_wp, 0.0_wp, 0.0_wp]
        ey = [0.0_wp, 1.0_wp, 0.0_wp]
        ez = [0.0_wp, 0.0_wp, 1.0_wp]

        ! Act
        tp = triple_product(ex, ey, ez)

        ! Assert
        call assert_close_real(tp, 1.0_wp, 1.0e-12_wp, 'triple product unit basis', failures)
    end subroutine test_triple_product_right_handed

    subroutine test_triple_product_zero_for_coplanar(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: a, b, c
        real(kind=wp) :: tp

        ! Arrange
        a = [1.0_wp, 2.0_wp, 0.0_wp]
        b = [2.0_wp, 4.0_wp, 0.0_wp]
        c = [3.0_wp, -1.0_wp, 0.0_wp]

        ! Act
        tp = triple_product(a, b, c)

        ! Assert
        call assert_close_real(tp, 0.0_wp, 1.0e-12_wp, 'triple product coplanar vectors', failures)
        call assert_false(abs(tp) > 1.0e-12_wp, 'triple product coplanar should be near zero', failures)
        call assert_true(abs(tp) <= 1.0e-12_wp, 'triple product coplanar in tolerance', failures)
    end subroutine test_triple_product_zero_for_coplanar

    subroutine test_in_block_bd_interior_and_upper_boundary(failures)
        integer, intent(inout) :: failures
        real(kind=dp), dimension(3) :: v_inside, v_at_upper
        real(kind=dp), dimension(6) :: bounds

        ! Arrange
        bounds = [0.0_dp, 10.0_dp, -5.0_dp, 5.0_dp, 1.0_dp, 2.0_dp]
        v_inside = [5.0_dp, 0.0_dp, 1.5_dp]
        v_at_upper = [5.0_dp, 0.0_dp, 2.0_dp]

        ! Act / Assert
        call assert_true(in_block_bd(v_inside, bounds), 'interior point should be inside bounds', failures)
        call assert_false(in_block_bd(v_at_upper, bounds), 'upper boundary should be excluded', failures)
    end subroutine test_in_block_bd_interior_and_upper_boundary

    subroutine test_in_block_bd_optional_dimension_mask(failures)
        integer, intent(inout) :: failures
        real(kind=dp), dimension(3) :: v
        real(kind=dp), dimension(6) :: bounds
        logical, dimension(3) :: check_bd

        ! Arrange
        bounds = [0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 2.0_dp]
        v = [1.0_dp, 5.0_dp, 1.0_dp]
        check_bd = [.true., .false., .true.]

        ! Act / Assert
        call assert_true(in_block_bd(v, bounds, check_bd), 'disabled dimension should be ignored by in_block_bd', failures)
    end subroutine test_in_block_bd_optional_dimension_mask

    subroutine test_in_cutoff_threshold_behavior(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: rl, rk_at_rcmin, rk_at_rcoff, rk_outside

        ! Arrange
        rl = [0.0_wp, 0.0_wp, 0.0_wp]
        rk_at_rcmin = [1.0_wp, 0.0_wp, 0.0_wp]
        rk_at_rcoff = [2.0_wp, 0.0_wp, 0.0_wp]
        rk_outside = [2.1_wp, 0.0_wp, 0.0_wp]

        ! Act / Assert
        call assert_true(in_cutoff(rl, rk_at_rcmin, 1.0_wp, 2.0_wp), 'distance at rcmin should be included', failures)
        call assert_true(in_cutoff(rl, rk_at_rcoff, 1.0_wp, 2.0_wp), 'distance at rcoff should be included', failures)
        call assert_false(in_cutoff(rl, rk_outside, 1.0_wp, 2.0_wp), 'distance beyond rcoff should be excluded', failures)
    end subroutine test_in_cutoff_threshold_behavior

    subroutine test_neighbor_dis_free_components_and_norm(failures)
        integer, intent(inout) :: failures
        real(kind=wp), dimension(3) :: rl, rk
        real(kind=wp), dimension(4) :: dis

        ! Arrange
        rl = [1.0_wp, 2.0_wp, 3.0_wp]
        rk = [4.0_wp, 6.0_wp, 3.0_wp]

        ! Act
        dis = neighbor_dis_free(rl, rk)

        ! Assert
        call assert_close_real(dis(1), 3.0_wp, 1.0e-12_wp, 'neighbor_dis_free dx', failures)
        call assert_close_real(dis(2), 4.0_wp, 1.0e-12_wp, 'neighbor_dis_free dy', failures)
        call assert_close_real(dis(3), 0.0_wp, 1.0e-12_wp, 'neighbor_dis_free dz', failures)
        call assert_close_real(dis(4), 5.0_wp, 1.0e-12_wp, 'neighbor_dis_free norm', failures)
    end subroutine test_neighbor_dis_free_components_and_norm

    subroutine test_dumb_interp_midpoint(failures)
        integer, intent(inout) :: failures
        real(kind=dp) :: y

        ! Arrange / Act
        y = dumb_interp(1.5_dp, 1.0_dp, 10.0_dp, 2.0_dp, 14.0_dp)

        ! Assert
        call assert_close_real(y, 12.0_wp, 1.0e-12_wp, 'dumb_interp should linearly interpolate midpoint value', failures)
    end subroutine test_dumb_interp_midpoint

end program test_math
