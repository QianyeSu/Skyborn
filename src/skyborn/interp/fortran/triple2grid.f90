! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-17
! File    : triple2grid.f90
! Purpose : Map scattered triple-list observations onto a rectilinear target
!           grid through the legacy `triple2grid` entry points.
! Notes   : This free-form source preserves the historical `trip2grd2` and
!           `trip2grd3` behaviors, including their exact-match prepass and
!           nearest-point assignment semantics.
! =============================================================================
!
! QUICK REFERENCE
! PURPOSE
!    MAP AN IRREGULAR (X, Y, Z) TRIPLE LIST ONTO A TARGET RECTILINEAR
!    GRID USING THE HISTORICAL NCL `triple2grid` ENTRY POINTS.
!    THIS IS THE MODERN FREE-FORM REPLACEMENT FOR
!    `archive/triple2grid.f`, WHILE PRESERVING THE LEGACY PUBLIC API.
!
! EXPECTED INPUT SHAPES
!    XI(KZ), YI(KZ), ZI(KZ)
!                 INPUT TRIPLE COORDINATES AND VALUES
!    GX(MX), GY(NY)
!                 TARGET GRID X / Y AXES
!
! FLAGS
!    ZMSG       - MISSING VALUE SENTINEL
!    DOMAIN     - OUTLIER PADDING SCALE FOR THE OVERSIZED "BIG GRID"
!    LOOP       - 0: PREFER `trip2grd2`
!                 1: USE `trip2grd3`
!    METHOD     - 0: PLANE DISTANCE
!                 1: GREAT-CIRCLE DISTANCE
!    DISTMX     - MAXIMUM ALLOWED ASSIGNMENT DISTANCE
!
! OUTPUT
!    GRID(MX, NY) HOLDS THE NEAREST-POINT RESULT.
!    IER=0 ON SUCCESS.
!    IER=-1 / -11 WHEN ONLY THE EXACT-MATCH PREPASS WAS NEEDED.

subroutine triple2grid1(kz, xi, yi, zi, zmsg, mx, ny, gx, gy, grid, domain, loop, method, distmx, mx2, ny2, &
                        x, y, z, gbigx, gbigy, gbigxy, ier)
    implicit none

    integer, intent(in) :: mx, ny, kz, mx2, ny2, loop, method
    integer, intent(out) :: ier
    integer :: m, n, k, kout, kpts, mflag, nflag
    double precision, intent(in) :: gx(mx), gy(ny), distmx, domain
    double precision, intent(in) :: xi(kz), yi(kz), zi(kz), zmsg
    double precision, intent(inout) :: x(kz), y(kz), z(kz)
    double precision, intent(inout) :: gbigx(mx2), gbigy(ny2), gbigxy(mx2, ny2)
    double precision, intent(out) :: grid(mx, ny)
    double precision :: dd, ddeps, ddcrit

    ier = 0
    ddeps = 1.0d-3
    if (distmx <= 0.0d0) then
        ddcrit = 1.0d20
    else
        ddcrit = distmx
    end if

    ! Strip missing inputs while counting values that lie outside the target
    ! grid bounds.
    kpts = 0
    kout = 0
    do k = 1, kz
        if (zi(k) /= zmsg) then
            kpts = kpts + 1
            x(kpts) = xi(k)
            y(kpts) = yi(k)
            z(kpts) = zi(k)
            if (xi(k) < gx(1) .or. xi(k) > gx(mx) .or. yi(k) < gy(1) .or. yi(k) > gy(ny)) then
                kout = kout + 1
            end if
        end if
    end do

    if (kpts == 0) then
        do n = 1, ny
            do m = 1, mx
                grid(m, n) = zmsg
            end do
        end do
        return
    end if

    ! Detect approximately equal target-grid spacing to enable the faster
    ! index arithmetic used by TRIP2GRD2 / TRIP2GRD3.
    mflag = 1
    dd = abs(gx(2) - gx(1))
    do m = 2, mx - 1
        if (dd < (gx(m + 1) - gx(m)) - ddeps .or. dd > (gx(m + 1) - gx(m)) + ddeps) then
            mflag = 0
            exit
        end if
    end do

    nflag = 1
    dd = abs(gy(2) - gy(1))
    do n = 2, ny - 1
        if (dd < (gy(n + 1) - gy(n)) - ddeps .or. dd > (gy(n + 1) - gy(n)) + ddeps) then
            nflag = 0
            exit
        end if
    end do

    ! NCL retained LOOP=0 because TRIP2GRD3 has a long-standing correctness
    ! issue on some platforms. Keep the historical branch anyway because the
    ! direct F2PY interface exposes both paths.
    if (kout == 0) then
        if (loop == 0) then
            call trip2grd2(kpts, x, y, z, zmsg, mx, ny, gx, gy, grid, &
                           mflag, nflag, method, ddcrit, ier)
        else
            call trip2grd3(kpts, x, y, z, zmsg, mx, ny, gx, gy, grid, &
                           mflag, nflag, method, ddcrit, ier)
        end if
    else
        ! Expand the target domain by one extra ring so out-of-bounds input
        ! points can still influence the edge cells, matching the legacy path.
        do n = 1, ny
            gbigy(n + 1) = gy(n)
        end do

        do m = 1, mx
            gbigx(m + 1) = gx(m)
        end do

        gbigy(1) = gy(1) - domain * (gy(2) - gy(1))
        gbigy(ny2) = gy(ny) + domain * (gy(ny) - gy(ny - 1))
        gbigx(1) = gx(1) - domain * (gx(2) - gx(1))
        gbigx(mx2) = gx(mx) + domain * (gx(mx) - gx(mx - 1))

        if (loop == 0) then
            call trip2grd2(kpts, x, y, z, zmsg, mx2, ny2, gbigx, gbigy, gbigxy, &
                           mflag, nflag, method, ddcrit, ier)
        else
            call trip2grd3(kpts, x, y, z, zmsg, mx2, ny2, gbigx, gbigy, gbigxy, &
                           mflag, nflag, method, ddcrit, ier)
        end if

        ! Copy the interior of the oversized grid back to the requested domain.
        do n = 1, ny
            do m = 1, mx
                grid(m, n) = gbigxy(m + 1, n + 1)
            end do
        end do
    end if
end subroutine triple2grid1


! QUICK REFERENCE
! PURPOSE
!    ASSIGN EACH OBSERVATION TO THE NEAREST ELIGIBLE TARGET GRID POINT.
!    THIS IS THE HISTORICAL `TRIP2GRD2` PATH THAT NCL PREFERS OVER
!    `TRIP2GRD3` BECAUSE OF A LONG-STANDING CORRECTNESS ISSUE IN THAT
!    ALTERNATIVE KERNEL.
!
! FLAGS
!    MFLAG/NFLAG - APPROXIMATE EQUAL-SPACING FLAGS FOR GXOUT / GYOUT
!    METHOD      - 0: PLANE DISTANCE
!                  1: GREAT-CIRCLE DISTANCE
!    DDCRIT      - MAXIMUM ALLOWED ASSIGNMENT DISTANCE
!
! OUTPUT
!    GOUT(MX, NY) HOLDS THE NEAREST-POINT RESULT.
!    IER=-1 WHEN THE EXACT-MATCH PREPASS ASSIGNED EVERY INPUT POINT.
subroutine trip2grd2(kz, x, y, z, zmsg, mx, ny, gxout, gyout, gout, mflag, nflag, method, ddcrit, ier)
    implicit none

    integer, intent(in) :: mx, ny, kz, mflag, nflag, method
    integer, intent(out) :: ier
    integer :: m, n, k, mm, nn, ksum, kpts
    double precision, intent(in) :: gxout(mx), gyout(ny), x(kz), y(kz), z(kz), zmsg, ddcrit
    double precision, intent(out) :: gout(mx, ny)
    double precision :: dout(mx, ny)
    double precision :: dd, xx, yy, slpy, slpx, dx, dy, atmp, ylat, re, rad
    logical :: found_exact

    ier = 0

    re = 6371.2200d0
    rad = 4.0d0 * atan(1.0d0) / 180.0d0

    kpts = kz
    dx = abs(gxout(3) - gxout(2))
    dy = abs(gyout(3) - gyout(2))

    dout = 1.0d20
    gout = zmsg

    ! First preserve exact observation-to-grid-point matches.
    ksum = 0
    do k = 1, kpts
        if (mflag == 1 .and. nflag == 1) then
            mm = nint((x(k) - gxout(1)) / dx) + 1
            nn = nint((y(k) - gyout(1)) / dy) + 1
            if (mm < 1 .or. mm > mx .or. x(k) /= gxout(mm)) mm = -1
            if (nn < 1 .or. nn > ny .or. y(k) /= gyout(nn)) nn = -1
            if (mm >= 1 .and. nn >= 1) then
                gout(mm, nn) = z(k)
                dout(mm, nn) = 0.0d0
                ksum = ksum + 1
                cycle
            end if
        else
            found_exact = .false.
            do n = 1, ny
                if (y(k) == gyout(n)) then
                    do m = 1, mx
                        if (x(k) == gxout(m)) then
                            gout(m, n) = z(k)
                            dout(m, n) = 0.0d0
                            ksum = ksum + 1
                            found_exact = .true.
                            exit
                        end if
                    end do
                    if (found_exact) exit
                end if
            end do
        end if
    end do

    if (ksum == kpts) then
        ier = -1
        return
    end if

    ! Otherwise assign each observation to the nearest eligible grid point.
    ! When multiple observations compete for one cell, keep the closest one.
    ksum = 0
    do k = 1, kpts
        nn = -1
        if (nflag == 1) then
            nn = int(((y(k) - gyout(1)) / dy)) + 2
        else
            do n = 1, ny - 1
                if (y(k) >= gyout(n) .and. y(k) < gyout(n + 1)) then
                    dy = y(k) - gyout(n)
                    slpy = dy / (gyout(n + 1) - gyout(n))
                    yy = n + slpy
                    nn = nint(yy)
                    exit
                end if
            end do
        end if

        mm = -1
        if (mflag == 1) then
            mm = int(((x(k) - gxout(1)) / dx)) + 2
        else
            do m = 1, mx - 1
                if (x(k) >= gxout(m) .and. x(k) < gxout(m + 1)) then
                    dx = x(k) - gxout(m)
                    slpx = dx / (gxout(m + 1) - gxout(m))
                    xx = m + slpx
                    mm = nint(xx)
                    exit
                end if
            end do
        end if

        if (mm >= 1 .and. mm <= mx .and. nn >= 1 .and. nn <= ny) then
            if (method == 0) then
                dd = sqrt(slpx ** 2 + slpy ** 2)
            else
                ylat = y(k) * rad
                atmp = sin(ylat) * sin(gxout(mm) * rad) + cos(ylat) * cos(gyout(nn) * rad) * cos((x(k) - gxout(mm)) * rad)
                atmp = min(1.0d0, max(-1.0d0, atmp))
                dd = acos(atmp) * re
            end if

            if (dd < dout(mm, nn) .and. dd < ddcrit) then
                ksum = ksum + 1
                gout(mm, nn) = z(k)
                dout(mm, nn) = dd
            end if
        end if
    end do
end subroutine trip2grd2


! QUICK REFERENCE
! PURPOSE
!    HISTORICAL ALTERNATIVE TO `TRIP2GRD2` THAT LOOPS OVER TARGET GRID
!    CELLS AND SEARCHES ALL OBSERVATIONS FOR THE CLOSEST MATCH.
!    THIS ENTRY POINT IS RETAINED FOR COMPATIBILITY WITH THE LEGACY
!    INTERFACE EVEN THOUGH NCL USUALLY FORCES `LOOP=0`.
!
! FLAGS
!    METHOD  - 0: PLANE DISTANCE
!              1: GREAT-CIRCLE DISTANCE
!    DDCRIT  - MAXIMUM ALLOWED ASSIGNMENT DISTANCE
!
! OUTPUT
!    GOUT(MX, NY) HOLDS THE NEAREST-POINT RESULT.
!    IER=-11 WHEN THE EXACT-MATCH PREPASS ASSIGNED EVERY INPUT POINT.
subroutine trip2grd3(kz, x, y, z, zmsg, mx, ny, gxout, gyout, gout, mflag, nflag, method, ddcrit, ier)
    implicit none

    integer, intent(in) :: mx, ny, kz, mflag, nflag, method
    integer, intent(out) :: ier
    integer :: m, n, k, ksum, kpts, mm, nn
    double precision, intent(in) :: gxout(mx), gyout(ny), x(kz), y(kz), z(kz), zmsg, ddcrit
    double precision, intent(out) :: gout(mx, ny)
    double precision :: dout(mx, ny), dist(kz), re, rad, rlat, atmp, dx, dy
    logical :: found_exact

    ier = 0
    kpts = kz
    dx = abs(gxout(3) - gxout(2))
    dy = abs(gyout(3) - gyout(2))

    dout = 1.0d20
    gout = zmsg

    ksum = 0
    do k = 1, kpts
        if (mflag == 1 .and. nflag == 1) then
            mm = nint((x(k) - gxout(1)) / dx) + 1
            nn = nint((y(k) - gyout(1)) / dy) + 1
            if (mm < 1 .or. mm > mx .or. x(k) /= gxout(mm)) mm = -1
            if (nn < 1 .or. nn > ny .or. y(k) /= gyout(nn)) nn = -1
            if (mm >= 1 .and. nn >= 1) then
                gout(mm, nn) = z(k)
                dout(mm, nn) = 0.0d0
                ksum = ksum + 1
                cycle
            end if
        else
            found_exact = .false.
            do n = 1, ny
                if (y(k) == gyout(n)) then
                    do m = 1, mx
                        if (x(k) == gxout(m)) then
                            gout(m, n) = z(k)
                            dout(m, n) = 0.0d0
                            ksum = ksum + 1
                            found_exact = .true.
                            exit
                        end if
                    end do
                    if (found_exact) exit
                end if
            end do
        end if
    end do

    if (ksum == kpts) then
        ier = -11
        return
    end if

    re = 6371.2200d0
    rad = 4.0d0 * atan(1.0d0) / 180.0d0

    ! For every target cell that did not already receive an exact match,
    ! search the full observation list for the closest eligible input point.
    do n = 1, ny
        do m = 1, mx
            if (dout(m, n) /= 0.0d0) then
                if (method == 0) then
                    do k = 1, kpts
                        dist(k) = sqrt((x(k) - gxout(m)) ** 2 + (y(k) - gyout(n)) ** 2)
                    end do
                else
                    rlat = gyout(n) * rad
                    do k = 1, kpts
                        ! Clamp the cosine argument so round-off does not
                        ! push ACOS outside [-1, 1].
                        atmp = sin(rlat) * sin(y(k) * rad) + &
                               cos(rlat) * cos(y(k) * rad) * &
                               cos((x(k) - gxout(m)) * rad)
                        atmp = min(1.0d0, max(-1.0d0, atmp))
                        dist(k) = acos(atmp) * re
                    end do
                end if

                do k = 1, kpts
                    if (z(k) /= zmsg) then
                        if (dist(k) < dout(m, n) .and. dist(k) < ddcrit) then
                            dout(m, n) = dist(k)
                            gout(m, n) = z(k)
                        end if
                    end if
                end do
            end if
        end do
    end do
end subroutine trip2grd3
