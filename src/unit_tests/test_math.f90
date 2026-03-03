program test_math
    use parameters, only: wp
    use math, only: identity_mat, cross_product, triple_product
    use assertions, only: assert_true, assert_false, assert_equal_int, assert_close_real
    implicit none

    integer :: failures

    failures = 0

    call test_identity_matrix(failures)
    call test_cross_product_basis(failures)
    call test_triple_product_right_handed(failures)
    call test_triple_product_zero_for_coplanar(failures)

    if (failures > 0) then
        print '(A,1X,I0)', 'Unit tests failed:', failures
        stop 1
    end if

    print '(A)', 'All math unit tests passed.'

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

end program test_math
