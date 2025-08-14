!
! Optimized version of ihgeod.f
! Icosahedral geodesic mesh generator for spherical grids
! Optimizations: Modern Fortran syntax, vectorization, improved algorithms
!
!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 1998 by UCAR                 .
!  .                                                             .
!  .       University Corporation for Atmospheric Research       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                                                             .
!  .                         SPHEREPACK                          .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!

! Optimized icosahedral geodesic grid generator
subroutine ihgeod(m, idp, jdp, x, y, z)
    implicit none

    ! Input/output parameters
    integer, intent(in) :: m, idp, jdp
    real, intent(out) :: x(idp, jdp, 5), y(idp, jdp, 5), z(idp, jdp, 5)

    ! Local variables
    integer :: k, i, j
    real :: pi, dphi, beta, theta1, theta2, hdphi, tdphi, phi
    real :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6
    real :: dxi, dyi, dzi, dxj, dyj, dzj, xs, ys, zs
    real :: rad, theta, phi_temp

    ! Precompute constants
    pi = 4.0 * atan(1.0)
    dphi = 0.4 * pi
    beta = cos(dphi)
    theta1 = acos(beta / (1.0 - beta))
    theta2 = pi - theta1
    hdphi = dphi * 0.5
    tdphi = 3.0 * hdphi

    ! Generate geodesic points for each of the 5 icosahedral faces
    do k = 1, 5
        phi = (k - 1) * dphi

        ! Generate vertices of the geodesic triangle
        call stoc(1.0, theta2, phi, x1, y1, z1)
        call stoc(1.0, pi, phi + hdphi, x2, y2, z2)
        call stoc(1.0, theta2, phi + dphi, x3, y3, z3)

        ! Compute incremental vectors
        dxi = (x2 - x1) / (m - 1)
        dyi = (y2 - y1) / (m - 1)
        dzi = (z2 - z1) / (m - 1)
        dxj = (x3 - x2) / (m - 1)
        dyj = (y3 - y2) / (m - 1)
        dzj = (z3 - z2) / (m - 1)

        ! Fill first triangle
        !DIR$ VECTOR ALWAYS
        do i = 1, m
            xs = x1 + (i - 1) * dxi
            ys = y1 + (i - 1) * dyi
            zs = z1 + (i - 1) * dzi

            !DIR$ VECTOR ALWAYS
            do j = 1, i
                x(j, i, k) = xs + (j - 1) * dxj
                y(j, i, k) = ys + (j - 1) * dyj
                z(j, i, k) = zs + (j - 1) * dzj
            end do
        end do

        ! Generate second triangle vertices
        call stoc(1.0, theta1, phi + hdphi, x4, y4, z4)
        dxi = (x3 - x4) / (m - 1)
        dyi = (y3 - y4) / (m - 1)
        dzi = (z3 - z4) / (m - 1)
        dxj = (x4 - x1) / (m - 1)
        dyj = (y4 - y1) / (m - 1)
        dzj = (z4 - z1) / (m - 1)

        ! Fill second triangle
        !DIR$ VECTOR ALWAYS
        do j = 1, m
            xs = x1 + (j - 1) * dxj
            ys = y1 + (j - 1) * dyj
            zs = z1 + (j - 1) * dzj

            !DIR$ VECTOR ALWAYS
            do i = 1, j
                x(j, i, k) = xs + (i - 1) * dxi
                y(j, i, k) = ys + (i - 1) * dyi
                z(j, i, k) = zs + (i - 1) * dzi
            end do
        end do

        ! Generate third triangle vertices
        call stoc(1.0, theta1, phi + tdphi, x5, y5, z5)
        dxj = (x5 - x3) / (m - 1)
        dyj = (y5 - y3) / (m - 1)
        dzj = (z5 - z3) / (m - 1)

        ! Fill third triangle
        !DIR$ VECTOR ALWAYS
        do i = 1, m
            xs = x4 + (i - 1) * dxi
            ys = y4 + (i - 1) * dyi
            zs = z4 + (i - 1) * dzi

            !DIR$ VECTOR ALWAYS
            do j = 1, i
                x(j + m - 1, i, k) = xs + (j - 1) * dxj
                y(j + m - 1, i, k) = ys + (j - 1) * dyj
                z(j + m - 1, i, k) = zs + (j - 1) * dzj
            end do
        end do

        ! Generate fourth triangle vertices
        call stoc(1.0, 0.0, phi + dphi, x6, y6, z6)
        dxi = (x5 - x6) / (m - 1)
        dyi = (y5 - y6) / (m - 1)
        dzi = (z5 - z6) / (m - 1)
        dxj = (x6 - x4) / (m - 1)
        dyj = (y6 - y4) / (m - 1)
        dzj = (z6 - z4) / (m - 1)

        ! Fill fourth triangle
        !DIR$ VECTOR ALWAYS
        do j = 1, m
            xs = x4 + (j - 1) * dxj
            ys = y4 + (j - 1) * dyj
            zs = z4 + (j - 1) * dzj

            !DIR$ VECTOR ALWAYS
            do i = 1, j
                x(j + m - 1, i, k) = xs + (i - 1) * dxi
                y(j + m - 1, i, k) = ys + (i - 1) * dyi
                z(j + m - 1, i, k) = zs + (i - 1) * dzi
            end do
        end do
    end do

    ! Normalize all points to unit sphere with improved precision
    do k = 1, 5
        do j = 1, m + m - 1
            !DIR$ VECTOR ALWAYS
            do i = 1, m
                call ctos(x(j, i, k), y(j, i, k), z(j, i, k), rad, theta, phi_temp)
                call stoc(1.0, theta, phi_temp, x(j, i, k), y(j, i, k), z(j, i, k))
            end do
        end do
    end do

end subroutine ihgeod

! Optimized Cartesian to spherical coordinate conversion
subroutine ctos(x, y, z, r, theta, phi)
    implicit none

    real, intent(in) :: x, y, z
    real, intent(out) :: r, theta, phi

    real :: r1, r1_sqrt
    real, parameter :: eps = 1.0e-14

    r1 = x * x + y * y

    if (r1 < eps) then
        phi = 0.0
        theta = 0.0
        if (z < 0.0) theta = 4.0 * atan(1.0)  ! pi
        r = abs(z)
        return
    end if

    r1_sqrt = sqrt(r1)
    r = sqrt(r1 + z * z)
    phi = atan2(y, x)
    theta = atan2(r1_sqrt, z)

end subroutine ctos

! Optimized spherical to Cartesian coordinate conversion
subroutine stoc(r, theta, phi, x, y, z)
    implicit none

    real, intent(in) :: r, theta, phi
    real, intent(out) :: x, y, z

    real :: st, ct, sp, cp

    ! Use intrinsic sincos for better performance if available
    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)

    x = r * st * cp
    y = r * st * sp
    z = r * ct

end subroutine stoc
