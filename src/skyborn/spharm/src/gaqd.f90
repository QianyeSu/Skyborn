!
! Optimized version of gaqd.f
! Computes Gaussian quadrature points and weights
! Optimizations: Modern Fortran syntax, improved convergence, vectorization
!
subroutine gaqd(nlat, theta, wts, dwork, ldwork, ierror)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, ldwork
    real(8), intent(out) :: theta(nlat), wts(nlat)
    real(8), intent(inout) :: dwork(ldwork)
    integer, intent(out) :: ierror

    ! Function declarations
    real(8) :: dzeps
    external :: dzeps

    ! Local variables
    real(8) :: x, pi, pis2, dtheta, dthalf, cmax, zprev, zlast, zero
    real(8) :: zhold, pb, dpb, dcor, sum, cz, eps, sgnd
    integer :: mnlat, ns2, nhalf, idx, nix, it, i

    ! Error checking
    ierror = 1
    if (nlat <= 0) return
    ierror = 0

    ! Handle special cases analytically
    if (nlat == 1) then
        theta(1) = acos(0.0d0)
        wts(1) = 2.0d0
        return
    end if

    if (nlat == 2) then
        x = sqrt(1.0d0/3.0d0)
        theta(1) = acos(x)
        theta(2) = acos(-x)
        wts(1) = 1.0d0
        wts(2) = 1.0d0
        return
    end if

    ! Constants and initialization
    eps = sqrt(dzeps(1.0d0))
    eps = eps * sqrt(eps)
    pis2 = 2.0d0 * atan(1.0d0)
    pi = pis2 + pis2
    mnlat = mod(nlat, 2)
    ns2 = nlat / 2
    nhalf = (nlat + 1) / 2
    idx = ns2 + 2

    ! Compute initial coefficients
    call cpdp(nlat, cz, theta(ns2+1), wts(ns2+1))

    dtheta = pis2 / nhalf
    dthalf = dtheta / 2.0d0
    cmax = 0.2d0 * dtheta

    ! Estimate first point next to theta = pi/2
    if (mnlat /= 0) then
        zero = pis2 - dtheta
        zprev = pis2
        nix = nhalf - 1
    else
        zero = pis2 - dthalf
        nix = nhalf
    end if

    ! Main iteration loop for finding zeros
    do while (nix > 0)
        it = 0

        ! Newton iterations with improved convergence and safeguards
        do
            it = it + 1
            zlast = zero

            ! Evaluate polynomial and derivative
            call tpdp(nlat, zero, cz, theta(ns2+1), wts(ns2+1), pb, dpb)

            ! Safeguard against zero derivative
            if (abs(dpb) < epsilon(dpb) * 1000.0d0) then
                dcor = sign(eps, pb)
            else
                dcor = pb / dpb
            end if

            sgnd = 1.0d0
            if (dcor /= 0.0d0) sgnd = dcor / abs(dcor)
            dcor = sgnd * min(abs(dcor), cmax)
            zero = zero - dcor

            ! Enhanced convergence check with relative and absolute tolerance
            if (abs(zero - zlast) <= eps * max(abs(zero), 1.0d0)) exit

            ! Prevent infinite iterations
            if (it > 50) exit
        end do

        theta(nix) = zero
        zhold = zero

        ! Yakimiw's formula for weights (more stable) with safeguards
        if (abs(sin(zlast)) > epsilon(sin(zlast)) * 1000.0d0) then
            wts(nix) = (nlat + nlat + 1) / (dpb + pb * cos(zlast) / sin(zlast))**2
        else
            ! Fallback to standard formula if sin(zlast) is too small
            wts(nix) = (nlat + nlat + 1) / (dpb * dpb)
        end if

        nix = nix - 1
        if (nix == 0) exit

        ! Update zero estimate for next iteration
        if (nix == nhalf - 1) then
            zero = 3.0d0 * zero - pi
        else if (nix < nhalf - 1) then
            zero = zero + zero - zprev
        end if
        zprev = zhold
    end do

    ! Handle center point for odd nlat
    if (mnlat /= 0) then
        theta(nhalf) = pis2
        call tpdp(nlat, pis2, cz, theta(ns2+1), wts(ns2+1), pb, dpb)
        wts(nhalf) = (nlat + nlat + 1) / (dpb * dpb)
    end if

    ! Extend points and weights via symmetries
    do i = 1, ns2
        wts(nlat - i + 1) = wts(i)
        theta(nlat - i + 1) = pi - theta(i)
    end do

    ! Normalize weights with Kahan summation for better accuracy
    sum = 0.0d0
    if (nlat > 1000) then
        ! Use Kahan summation for large arrays
        call kahan_sum(wts, nlat, sum)
    else
        do i = 1, nlat
            sum = sum + wts(i)
        end do
    end if

    ! Safeguard against division by zero
    if (abs(sum) < epsilon(sum) * nlat) then
        sum = 2.0d0  ! Theoretical total weight
    end if

    do i = 1, nlat
        wts(i) = 2.0d0 * wts(i) / sum
    end do

end subroutine gaqd

! Optimized coefficient computation
subroutine cpdp(n, cz, cp, dcp)
    implicit none
    integer, intent(in) :: n
    real(8), intent(out) :: cz
    real(8), intent(out) :: cp(n/2+1), dcp(n/2+1)

    integer :: ncp, j
    real(8) :: t1, t2, t3, t4

    ncp = (n + 1) / 2
    t1 = -1.0d0
    t2 = n + 1.0d0
    t3 = 0.0d0
    t4 = n + n + 1.0d0

    if (mod(n, 2) == 0) then
        cp(ncp) = 1.0d0
        do j = ncp, 2, -1
            t1 = t1 + 2.0d0
            t2 = t2 - 1.0d0
            t3 = t3 + 1.0d0
            t4 = t4 - 2.0d0
            cp(j-1) = (t1 * t2) / (t3 * t4) * cp(j)
        end do
        t1 = t1 + 2.0d0
        t2 = t2 - 1.0d0
        t3 = t3 + 1.0d0
        t4 = t4 - 2.0d0
        cz = (t1 * t2) / (t3 * t4) * cp(1)

        !DIR$ VECTOR ALWAYS
        do j = 1, ncp
            dcp(j) = (j + j) * cp(j)
        end do
    else
        cp(ncp) = 1.0d0
        do j = ncp - 1, 1, -1
            t1 = t1 + 2.0d0
            t2 = t2 - 1.0d0
            t3 = t3 + 1.0d0
            t4 = t4 - 2.0d0
            cp(j) = (t1 * t2) / (t3 * t4) * cp(j+1)
        end do

        !DIR$ VECTOR ALWAYS
        do j = 1, ncp
            dcp(j) = (j + j - 1) * cp(j)
        end do
    end if

end subroutine cpdp

! Optimized polynomial evaluation
subroutine tpdp(n, theta, cz, cp, dcp, pb, dpb)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: theta, cz
    real(8), intent(in) :: cp(n/2+1), dcp(n/2+1)
    real(8), intent(out) :: pb, dpb

    integer :: kdo, k
    real(8) :: fn, cdt, sdt, cth, sth, chh

    fn = n
    cdt = cos(theta + theta)
    sdt = sin(theta + theta)

    if (mod(n, 2) == 0) then
        ! n even
        kdo = n / 2
        pb = 0.5d0 * cz
        dpb = 0.0d0

        if (n > 0) then
            cth = cdt
            sth = sdt

            !DIR$ VECTOR ALWAYS
            do k = 1, kdo
                pb = pb + cp(k) * cth
                dpb = dpb - dcp(k) * sth
                chh = cdt * cth - sdt * sth
                sth = sdt * cth + cdt * sth
                cth = chh
            end do
        end if
    else
        ! n odd
        kdo = (n + 1) / 2
        pb = 0.0d0
        dpb = 0.0d0
        cth = cos(theta)
        sth = sin(theta)

        !DIR$ VECTOR ALWAYS
        do k = 1, kdo
            pb = pb + cp(k) * cth
            dpb = dpb - dcp(k) * sth
            chh = cdt * cth - sdt * sth
            sth = sdt * cth + cdt * sth
            cth = chh
        end do
    end if

end subroutine tpdp

! Unit roundoff estimation function (optimized for IEEE 754)
pure real(8) function dzeps(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: a, b, c, eps
    real(8), parameter :: zero_threshold = 1.0d-16
    integer, parameter :: max_iterations = 100
    integer :: iter

    a = 4.0d0 / 3.0d0
    iter = 0
    do
        iter = iter + 1
        b = a - 1.0d0
        c = b + b + b
        eps = abs(c - 1.0d0)
        if (eps > zero_threshold .or. iter >= max_iterations) exit
        a = a + epsilon(a) ! Increment slightly if stuck
    end do

    ! Use intrinsic epsilon for IEEE 754 compliance if algorithm fails
    if (iter >= max_iterations) then
        eps = epsilon(1.0d0)
    end if

    dzeps = eps * abs(x)

end function dzeps

! Kahan summation algorithm for improved numerical accuracy
pure subroutine kahan_sum(array, n, result)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: array(n)
    real(8), intent(out) :: result

    real(8) :: c, y, t
    integer :: i

    result = 0.0d0
    c = 0.0d0  ! Compensation for lost low-order bits

    do i = 1, n
        y = array(i) - c     ! Subtract the compensation
        t = result + y       ! Add y to sum
        c = (t - result) - y ! Recover the lost low-order bits
        result = t
    end do

end subroutine kahan_sum
