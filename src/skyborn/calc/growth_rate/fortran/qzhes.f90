subroutine qzhes(nm, n, a, b)
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none

    integer :: i, j, k, l, n, lb, l1, nm, nk1, nm1, nm2
    real(real64) :: a(nm, n), b(nm, n)
    real(real64) :: r, s, t, u1, u2, v1, v2, rho

    !     this subroutine is the first step of the qz algorithm
    !     for solving generalized matrix eigenvalue problems,
    !     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
    !
    !     this subroutine accepts a pair of real general matrices and
    !     reduces one of them to upper hessenberg form and the other
    !     to upper triangular form using orthogonal transformations.
    !     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
    !
    !     on input
    !
    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.
    !
    !        n is the order of the matrices.
    !
    !        a contains a real general matrix.
    !
    !        b contains a real general matrix.
    !
    !     on output
    !
    !        a has been reduced to upper hessenberg form.  the elements
    !          below the first subdiagonal have been set to zero.
    !
    !        b has been reduced to upper triangular form.  the elements
    !          below the main diagonal have been set to zero.
    !
    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory
    !
    !     this version dated august 1983.
    !
    !     ------------------------------------------------------------------

    !     .......... reduce b to upper triangular form ..........
    if (n <= 1) return

    nm1 = n - 1
    do l = 1, nm1
        l1 = l + 1
        s = 0.0_real64

        do i = l1, n
            s = s + abs(b(i, l))
        end do

        if (s == 0.0_real64) cycle

        s = s + abs(b(l, l))
        r = 0.0_real64
        do i = l, n
            b(i, l) = b(i, l) / s
            r = r + b(i, l) ** 2
        end do

        r = sign(sqrt(r), b(l, l))
        b(l, l) = b(l, l) + r
        rho = r * b(l, l)

        do j = l1, n
            t = 0.0_real64
            do i = l, n
                t = t + b(i, l) * b(i, j)
            end do

            t = -t / rho
            do i = l, n
                b(i, j) = b(i, j) + t * b(i, l)
            end do
        end do

        do j = 1, n
            t = 0.0_real64
            do i = l, n
                t = t + b(i, l) * a(i, j)
            end do

            t = -t / rho
            do i = l, n
                a(i, j) = a(i, j) + t * b(i, l)
            end do
        end do

        b(l, l) = -s * r
        do i = l1, n
            b(i, l) = 0.0_real64
        end do
    end do

    !     .......... reduce a to upper hessenberg form, while
    !                keeping b triangular ..........
    if (n == 2) return

    nm2 = n - 2
    do k = 1, nm2
        nk1 = nm1 - k

        !     .......... for l=n-1 step -1 until k+1 do -- ..........
        do lb = 1, nk1
            l = n - lb
            l1 = l + 1

            !     .......... zero a(l+1,k) ..........
            s = abs(a(l, k)) + abs(a(l1, k))
            if (s == 0.0_real64) cycle

            u1 = a(l, k) / s
            u2 = a(l1, k) / s
            r = sign(sqrt(u1 * u1 + u2 * u2), u1)
            v1 = -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1

            do j = k, n
                t = a(l, j) + u2 * a(l1, j)
                a(l, j) = a(l, j) + t * v1
                a(l1, j) = a(l1, j) + t * v2
            end do

            a(l1, k) = 0.0_real64
            do j = l, n
                t = b(l, j) + u2 * b(l1, j)
                b(l, j) = b(l, j) + t * v1
                b(l1, j) = b(l1, j) + t * v2
            end do

            !     .......... zero b(l+1,l) ..........
            s = abs(b(l1, l1)) + abs(b(l1, l))
            if (s == 0.0_real64) cycle

            u1 = b(l1, l1) / s
            u2 = b(l1, l) / s
            r = sign(sqrt(u1 * u1 + u2 * u2), u1)
            v1 = -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1

            do i = 1, l1
                t = b(i, l1) + u2 * b(i, l)
                b(i, l1) = b(i, l1) + t * v1
                b(i, l) = b(i, l) + t * v2
            end do

            b(l1, l) = 0.0_real64
            do i = 1, n
                t = a(i, l1) + u2 * a(i, l)
                a(i, l1) = a(i, l1) + t * v1
                a(i, l) = a(i, l) + t * v2
            end do

        end do
    end do
end subroutine qzhes
