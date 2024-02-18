module math
    !This contains math helper functions 
    use mpi
    use parameters

    implicit none

    public

    contains 
    subroutine sort_array(array_in, n, array_out, array_order)
        !sort 1d real array from minimum to maximum

        integer, intent(in) :: n
        integer, dimension(n), intent(out) :: array_order
        real(kind = wp), dimension(n), intent(in) :: array_in
        real(kind = wp), dimension(n), intent(out) :: array_out


        integer :: i, j, temp_int
        real(kind = wp) :: temp

        array_out(:) = array_in(:)
        do i = 1, n
            array_order(i) = i
        end do

        do i = 1, n
            do j = i+1, n
                if(array_out(i) > array_out(j)) then
                    temp = array_out(j)
                    array_out(j) = array_out(i)
                    array_out(i) = temp
                    temp_int = array_order(j)
                    array_order(j) = array_order(i)
                    array_order(i) = temp_int
                end if
            end do
        end do
        return
    end subroutine sort_array

    pure function identity_mat(n)
        !real identity matrix of n by n

        implicit none

        integer, intent(in) :: n
        real(kind = wp), dimension(n, n) :: identity_mat

        integer :: i

        identity_mat(:, :) = 0.0_wp

        do i = 1, n
            identity_mat(i, i) = 1.0_wp
        end do

        return
    end function identity_mat
 
    pure function cross_product(a, b)
        !cross product between two 3*1 vectors
        implicit none
        real(kind = wp), dimension(3), intent(in) :: a, b
        real(kind = wp), dimension(3) :: cross_product

        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - a(2) * b(1)
        return
     end function cross_product

    pure function triple_product(a, b, c)
        !triple product between three 3*1 vectors
        implicit none
        real(kind = wp), dimension(3), intent(in) :: a, b, c
        real(kind = wp) :: triple_product

        triple_product = dot_product(a, cross_product(b, c))

        return
    end function triple_product

    subroutine matrix_normal(a, n, a_nor)
        !normalize an n by n matrix
        implicit none
        integer, intent(in) :: n
        real(kind = wp), dimension(n, n), intent(in) :: a
        real(kind = wp), dimension(n, n), intent(out) :: a_nor

        real(kind = wp), dimension(n) :: v
        integer :: i

        a_nor(:, :) = a(:, :)
        do i = 1, n
            v(:) = a(:, i)
            a_nor(:, i) = v(:) / norm2(v)
        end do
        return
    end subroutine matrix_normal

     subroutine matrix_lattice(a, n, iorh)
        !check if a matrix is orthogonal and obey right hand rule
        implicit none
        integer, intent(in) :: n
        real(kind = wp), dimension(n, n), intent(in) :: a
        logical, dimension(2), intent(out) :: iorh

        real(kind = wp), dimension(n) :: v, v_k
        integer :: i, j

        iorh(:) = .true.
        i_loop: do i = 1, n
            do j = i + 1, n

                if(abs(dot_product(a(:, i), a(:, j))) > lim_small) then
                    iorh(1) = .false.
                end if

                if(j == i + 1) then
                    v(:) = cross_product(a(:, i), a(:, j))
                    v_k(:) = v(:) - a(:, mod(j, n)+1)

                else if((i == 1).and.(j == n)) then
                    v(:) = cross_product(a(:, j), a(:, i))
                    v_k(:) = v(:) - a(:, i+1)
                
                end if

                if(norm2(v_k(:)) > lim_small) then
                    iorh(2) = .false.
                end if

                if(all(iorh).eqv..false.) then
                    exit i_loop
                end if
            end do
        end do i_loop
        return

    end subroutine matrix_lattice

    subroutine matrix_inverse(a, n, a_inv)
        !inverse an n by n matrix
        implicit none
        integer, intent(in) :: n
        real(kind = wp), dimension(n, n), intent(in) :: a
        real(kind = wp), dimension(n, n), intent(out) :: a_inv
        
        integer :: i, j, k, piv_loc
        real(kind = wp) :: coeff, sum_l, sum_u
        real(kind = wp), dimension(n) :: b, x, y, b_piv
        real(kind = wp), dimension(n, n) :: l, u, p
        real(kind = wp), allocatable :: v(:), u_temp(:), l_temp(:), p_temp(:)

        l(:, :) = identity_mat(n)
        u(:, :) = a(:, :)
        p(:, :) = identity_mat(n)

        !LU decomposition with partial pivoting
        do j = 1, n-1

            allocate( v(n-j+1), stat = allostat )
            if(allostat /=0 ) then
                print *, 'Fail to allocate v in matrix_inverse'
                call mpi_abort(1,mpi_comm_world,ierr)
            end if

            v(:) = u(j:n, j)
            if(maxval(abs(v)) < lim_zero) then
                print *, 'Fail to inverse matrix', a
                call mpi_abort(1,mpi_comm_world,ierr)
            end if

            piv_loc = maxloc(abs(v), 1)

            deallocate( v,  stat = deallostat )
            if(deallostat /=0 ) then
                print *, 'Fail to deallocate v in matrix_inverse'
                call mpi_abort(1,mpi_comm_world,ierr)
            end if

            !partial pivoting

            if(piv_loc /= 1) then

                allocate( u_temp(n-j+1), p_temp(n), stat = allostat) 
                if(allostat /=0 ) then
                    print *, 'Fail to allocate p_temp and/or u_temp in matrix_inverse'
                    call mpi_abort(1,mpi_comm_world,ierr)
                end if

                u_temp(:) = u(j, j:n)
                u(j, j:n) = u(piv_loc+j-1, j:n)
                u(piv_loc+j-1, j:n) = u_temp(:)
                p_temp(:) = p(j, :)
                p(j, :) = p(piv_loc+j-1, :)
                p(piv_loc+j-1, :) = p_temp(:)

                deallocate( u_temp, p_temp, stat = deallostat )
                if(deallostat /=0 ) then
                    print *, 'Fail to deallocate p_temp and/or u_temp in matrix_inverse'
                    call mpi_abort(1,mpi_comm_world,ierr)
                end if

                if(j > 1) then
                    allocate( l_temp(j-1), stat = allostat )
                    if(allostat /= 0) then
                        print *, 'Fail to allocate l_temp in matrix_inverse'
                        call mpi_abort(1,mpi_comm_world,ierr)
                    end if

                    l_temp(:) = l(j, 1:j-1)
                    l(j, 1:j-1) = l(piv_loc+j-1, 1:j-1)
                    l(piv_loc+j-1, 1:j-1) = l_temp(:)
                    
                    deallocate(l_temp, stat = deallostat) 
                    if(deallostat /=0 ) then
                        print *, 'Fail to deallocate l_temp in matrix_inverse'
                        call mpi_abort(1,mpi_comm_world,ierr)
                    end if
                end if
            end if

            !LU decomposition
            do i = j+1, n
                coeff = u(i, j)/u(j, j)
                l(i, j) = coeff
                u(i, j:n) = u(i, j:n)-coeff*u(j, j:n)
            end do
        end do

        a_inv(:, :) = 0.0_wp

        do j = 1, n
            b(:) = 0.0_wp
            b(j) = 1.0_wp
            b_piv(:) = matmul(p, b)

            !Now we have LUx = b_piv
            !the first step is to solve y from Ly = b_piv
            !forward substitution
            do i = 1, n
                if(i == 1) then
                    y(i) = b_piv(i)/l(i, i)
                else
                    sum_l = 0
                    do k = 1, i-1
                        sum_l = sum_l+l(i, k)*y(k)
                    end do
                    y(i) = (b_piv(i)-sum_l)/l(i, i)
                end if
            end do

        !then we solve x from ux = y
        !backward subsitution
            do i = n, 1, -1
                if(i == n) then
                    x(i) = y(i)/u(i, i)
                else
                    sum_u = 0
                    do k = i+1, n
                        sum_u = sum_u+u(i, k)*x(k)
                    end do
                    x(i) = (y(i)-sum_u)/u(i, i)
                end if
            end do

!      put x into j column of a_inv
            a_inv(:, j) = x(:)

        end do

        return
     end subroutine matrix_inverse
 
!    three point interpolation
 
     function interp3(x, tab) 
        
        !this comes from dl poly
        implicit none

        real(kind = wp), intent(in) :: x
        real(kind = wp), dimension(:), intent(in) :: tab
        real(kind = wp) :: interp3

        integer :: l
        real(kind = wp) :: rdr, rrr, ppp, gk0, gk1, gk2, t1, t2
 
 
        interp3 = 0.0_wp

        if(x > tab(3)) then
            interp3 = 0

        else

            rdr = 1.0_wp / tab(4)
            rrr = x - tab(2)

            l = min(nint(rrr*rdr), nint(tab(1))-1)
            if(l < 5) then
                interp3 = tab(5)
            end if

            ppp = rrr * rdr - real(l, wp)
            gk0 = tab(l-1)
            gk1 = tab(l)
            gk2 = tab(l+1)
            t1 = gk1 + (gk1 - gk0) * ppp
            t2 = gk1 + (gk2 - gk1) * ppp

            if(ppp < 0.0_wp) then
                interp3 = t1 + 0.5_wp * (t2 - t1) * (ppp + 1.0_wp)

            else if(l == 5) then
                interp3 = t2

            else
                interp3 = t2 + 0.5_wp * (t2 - t1) * (ppp - 1.0_wp)

            end if
        end if
 
        return
     end function interp3
 
 
     subroutine delete_duplicate(array_in, n, array_out, m)
        !delete duplicates in an 1-d integer array, the others are 0
        implicit none

        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: array_in
        integer, intent(out) :: m
        integer, dimension(n), intent(out) :: array_out

        integer :: i, j
        integer, dimension(n) :: array_temp
 
 
        m = 0
        array_out(:) = 0
        array_temp(:) = array_in(:)
 
        do i = 1, n
            if(array_temp(i) > 0) then
                do j = i+1, n
                    if(array_temp(i) == array_temp(j)) then
                        array_temp(j) = 0
                    end if
                end do

                m = m+1
                array_out(m) = array_temp(i)
 
            end if
        end do

        return
     end subroutine delete_duplicate
 

     subroutine normal_distribution(mean, dev, rand)
        !return a normally distributed random number
        implicit none
        
        real(kind = wp), intent(in) :: mean, dev
        real(kind = wp), intent(out) :: rand

        real(kind = wp) :: rand_uf, rand_us

        rand_uf = -1.0_wp

        do while(rand_uf < lim_zero)
            call random_number(rand_uf)
        end do

        if(rand_uf <= 0.5_wp) then
            call random_number(rand_us)
            rand = sqrt(-2.0_wp * log(rand_uf)) * cos(2.0_wp * pi * rand_us)

        else
            call random_number(rand_us)
            rand = sqrt(-2.0_wp * log(rand_uf)) * sin(2.0_wp * pi * rand_us)

        end if

        rand = mean + dev * rand

        return
    end subroutine normal_distribution

    function is_equal(A,B)
        real(kind=dp), intent(in) :: A, B
        logical :: is_equal

        if( (abs(A-B)) <= zerotol) then 
            is_equal =  .true.
        else
            is_equal = .false.
        end if
        return
    end function

    function dumb_interp(test_x, xlo, ylo, xhi, yhi)
        !Simple interpolation used just for validation in unit tests
        real(kind=dp), intent(in) :: test_x, xlo, ylo, xhi, yhi
        real(kind=dp) :: dumb_interp  
        
        dumb_interp = ((test_x - xlo)/(xhi-xlo))*(yhi-ylo) + ylo
        return
    end function

    pure function in_block_bd(v, block_bd, check_bd)
        !This function determines whether a point is within a block in 3d
        !Input/output
        real(kind=dp), dimension(3), intent(in) :: v
        real(kind=dp), dimension(6), intent(in) :: block_bd
        logical, dimension(3), optional, intent(in) :: check_bd
        logical :: in_block_bd

        !Other variables
        integer :: i 
        logical :: bool_bd(3)

        in_block_bd = .true.
        if(present(check_bd)) then 
            bool_bd = check_bd
        else
            bool_bd(:)=.true.
        end if


        do i =1 ,3
            if(bool_bd(i)) then 
                !Check upper bound
                if(v(i) >= (block_bd(2*i)) ) then 
                    in_block_bd =.false.
                    exit
                !Check lower bound
                else if (v(i) < (block_bd(2*i-1))) then 
                    in_block_bd = .false.
                    exit
                end if
            end if
        end do

    end function in_block_bd

    pure function in_cutoff(rl, rk, rcmin, rcoff)
        !This function returns true when rl and rk are separated by a distance between rcmin and rcoff
        real(kind=wp), dimension(3), intent(in) :: rl, rk
        real(kind=wp), intent(in) :: rcmin, rcoff
        logical :: in_cutoff

        real(kind=wp) :: rlk
        rlk=norm2(rk-rl)

        if (( rlk > rcoff).or.(rlk<rcmin)) then 
            in_cutoff=.false.
        else
            in_cutoff=.true.
        end if

        return
    end function

    pure function neighbor_dis_free(rl, rk)
        real(kind=wp), dimension(3), intent(in) :: rl, rk
        real(kind=wp), dimension(4) :: neighbor_dis_free

        neighbor_dis_free(1:3) = rk-rl
        neighbor_dis_free(4) = norm2(neighbor_dis_free(1:3))

        return
    end function

    subroutine interpolate(npoints, delta, tab, spline)
        !Interpolate a tab array using cubic splines
        integer, intent(in) :: npoints
        real(kind=wp), intent(in) :: delta
        real(kind = wp), dimension(npoints) :: tab
        real(kind=wp), dimension(7,npoints), intent(out) :: spline

        integer :: i

        spline(7,:) = tab(:)
        spline(6,1) = spline(7,2) - spline(7,1)
        spline(6,2) = 0.5_wp*(spline(7,3) - spline(7,1))
        spline(6, npoints-1) = 0.5_wp*(spline(7, npoints) - spline(7, npoints-2))
        spline(6, npoints) = spline(7,npoints) - spline(7,npoints-1)

        do i = 3, npoints-2
            spline(6,i) = ((spline(7,i-2)-spline(7,i+2)) + 8.0_wp*(spline(7,i+1)-spline(7,i-1)))/12.0_wp
        end do

        do i = 1, npoints-1
            spline(5,i) = 3.0_wp*(spline(7,i+1) - spline(7,i)) - 2.0_wp*spline(6,i) - spline(6,i+1)
            spline(4,i) = spline(6,i) + spline(6, i+1) - 2.0_wp*(spline(7,i+1) - spline(7,i))
        end do
        spline(5,npoints) = 0.0_wp
        spline(4,npoints) = 0.0_wp

        do i =1, npoints
            spline(3,i) = spline(6,i)/delta
            spline(2,i) = 2.0_wp*spline(5,i)/delta
            spline(1,i) = 3.0_wp*spline(4,i)/delta
        end do

    end subroutine interpolate

end module math
