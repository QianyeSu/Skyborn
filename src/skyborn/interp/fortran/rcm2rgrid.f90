! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-17
! File    : rcm2rgrid.f90
! Purpose : Remap between curvilinear and rectilinear grids while preserving
!           the legacy Skyborn/NCL interpolation entry points.
! Notes   : The helper module centralizes exact-hit, cell-search, and
!           short-gap cleanup logic so both public directions share the same
!           compatibility-sensitive arithmetic.
! =============================================================================
!
module rcm2rgrid_kernels_core
    use rcm_geodesy_core, only : dgcdist_core
    use linmsg_core, only : dlinmsg_interior_short_gaps
    implicit none

    integer, parameter :: real64 = selected_real_kind(15, 307)
contains

    ! Check whether a 1-D coordinate is strictly increasing so the caller can
    ! safely enable the monotonic-grid anchor fast path.
    pure logical function is_strictly_increasing(x) result(ok)
        real(real64), intent(in) :: x(:)
        integer :: idx

        ok = .true.
        if (size(x) <= 1) return

        do idx = 1, size(x) - 1
            if (x(idx + 1) <= x(idx)) then
                ok = .false.
                return
            end if
        end do
    end function is_strictly_increasing

    ! Return the lower bracketing index for a value on a strictly increasing
    ! coordinate using a binary search.
    pure integer function lower_bracket_increasing(x, value) result(idx)
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: value
        integer :: lo, hi, mid, n

        n = size(x)
        if (n <= 1) then
            idx = 1
            return
        end if

        if (value <= x(1)) then
            idx = 1
            return
        end if

        if (value >= x(n)) then
            idx = n - 1
            return
        end if

        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x(mid) <= value) then
                lo = mid
            else
                hi = mid
            end if
        end do

        idx = lo
    end function lower_bracket_increasing

    ! Shift the lower-bracket guess back to the earliest cell corner that can
    ! still span a stride-k curvilinear interpolation box.
    pure integer function earliest_stride_candidate(x, value, stride) result(idx)
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: value
        integer, intent(in) :: stride
        integer :: base, n

        n = size(x)
        if (n <= stride) then
            idx = 1
            return
        end if

        base = lower_bracket_increasing(x, value)
        idx = max(1, min(n - stride, base - stride + 1))
    end function earliest_stride_candidate

    ! Build earliest stride-k anchors for a monotonic batch of target
    ! coordinates by walking the source coordinate once instead of binary
    ! searching every target independently.
    pure subroutine build_stride_candidates_increasing(x, values, stride, idx)
        real(real64), intent(in) :: x(:), values(:)
        integer, intent(in) :: stride
        integer, intent(out) :: idx(size(values))
        integer :: base, cursor, j, n

        n = size(x)
        if (size(values) == 0) return

        if (n <= stride) then
            idx = 1
            return
        end if

        cursor = 1
        do j = 1, size(values)
            if (values(j) <= x(1)) then
                base = 1
            else if (values(j) >= x(n)) then
                base = n - 1
                cursor = n - 1
            else
                do while (cursor < n - 1 .and. x(cursor + 1) <= values(j))
                    cursor = cursor + 1
                end do
                base = cursor
            end if

            idx(j) = max(1, min(n - stride, base - stride + 1))
        end do
    end subroutine build_stride_candidates_increasing

    ! Search a small local curvilinear neighborhood for an exact source-grid
    ! point match before falling back to the global scan.
    subroutine find_exact_curv_local(xi, yi, xo, yo, eps, ix_lo, ix_hi, iy_lo, iy_hi, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:, :)
        real(real64), intent(in) :: yi(:, :)
        real(real64), intent(in) :: xo, yo, eps
        integer, intent(in) :: ix_lo, ix_hi, iy_lo, iy_hi
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = iy_lo, iy_hi
            do ix = ix_lo, ix_hi
                if (xo >= (xi(ix, iy) - eps) .and. xo <= (xi(ix, iy) + eps) .and. &
                    yo >= (yi(ix, iy) - eps) .and. yo <= (yi(ix, iy) + eps)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_exact_curv_local

    ! Full-grid exact-match scan for curvilinear source coordinates. This is
    ! only used when the fast-path local search does not succeed.
    subroutine find_exact_curv_full(xi, yi, xo, yo, eps, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:, :)
        real(real64), intent(in) :: yi(:, :)
        real(real64), intent(in) :: xo, yo, eps
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = 1, size(yi, 2)
            do ix = 1, size(xi, 1)
                if (xo >= (xi(ix, iy) - eps) .and. xo <= (xi(ix, iy) + eps) .and. &
                    yo >= (yi(ix, iy) - eps) .and. yo <= (yi(ix, iy) + eps)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_exact_curv_full

    ! Search a small local window for the curvilinear source cell whose
    ! corners bracket the requested output point.
    subroutine find_curv_cell_local(xi, yi, xo, yo, ix_lo, ix_hi, iy_lo, iy_hi, k, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:, :)
        real(real64), intent(in) :: yi(:, :)
        real(real64), intent(in) :: xo, yo
        integer, intent(in) :: ix_lo, ix_hi, iy_lo, iy_hi, k
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = iy_lo, iy_hi
            do ix = ix_lo, ix_hi
                if (xo >= xi(ix, iy) .and. xo <= xi(ix + k, iy) .and. &
                    yo >= yi(ix, iy) .and. yo <= yi(ix, iy + k)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_curv_cell_local

    ! Full-grid curvilinear cell search used as the robust fallback when the
    ! monotonic anchor heuristic cannot identify a local candidate.
    subroutine find_curv_cell_full(xi, yi, xo, yo, k, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:, :)
        real(real64), intent(in) :: yi(:, :)
        real(real64), intent(in) :: xo, yo
        integer, intent(in) :: k
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = 1, size(yi, 2) - k
            do ix = 1, size(xi, 1) - k
                if (xo >= xi(ix, iy) .and. xo <= xi(ix + k, iy) .and. &
                    yo >= yi(ix, iy) .and. yo <= yi(ix, iy + k)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_curv_cell_full

    ! Interpolate one output point from the four corners of a curvilinear cell.
    ! The output vector is updated only where it still holds xmsg so an earlier
    ! exact-match pass can preserve non-missing fields while letting missing
    ! fields fall through to interpolation, matching the historical kernel.
    subroutine interpolate_curv_cell(fi, xi, yi, xo, yo, ix, iy, k, xmsg, ncrit, out)
        real(real64), intent(in) :: fi(:, :, :)
        real(real64), intent(in) :: xi(:, :)
        real(real64), intent(in) :: yi(:, :)
        real(real64), intent(in) :: xo, yo, xmsg
        integer, intent(in) :: ix, iy, k, ncrit
        real(real64), intent(inout) :: out(size(fi, 3))
        real(real64) :: fw(2, 2), weights(2, 2), sumf, sumw
        real(real64) :: dist11, dist21, dist12, dist22
        integer :: ng, m, n, nw

        dist11 = dgcdist_core(yo, xo, yi(ix, iy), xi(ix, iy), 2)
        dist21 = dgcdist_core(yo, xo, yi(ix + k, iy), xi(ix + k, iy), 2)
        dist12 = dgcdist_core(yo, xo, yi(ix, iy + k), xi(ix, iy + k), 2)
        dist22 = dgcdist_core(yo, xo, yi(ix + k, iy + k), xi(ix + k, iy + k), 2)

        weights(1, 1) = merge(huge(1.0_real64), (1.0_real64 / dist11) ** 2, dist11 == 0.0_real64)
        weights(2, 1) = merge(huge(1.0_real64), (1.0_real64 / dist21) ** 2, dist21 == 0.0_real64)
        weights(1, 2) = merge(huge(1.0_real64), (1.0_real64 / dist12) ** 2, dist12 == 0.0_real64)
        weights(2, 2) = merge(huge(1.0_real64), (1.0_real64 / dist22) ** 2, dist22 == 0.0_real64)

        do ng = 1, size(fi, 3)
            if (out(ng) /= xmsg) cycle

            fw(1, 1) = fi(ix, iy, ng)
            fw(2, 1) = fi(ix + k, iy, ng)
            fw(1, 2) = fi(ix, iy + k, ng)
            fw(2, 2) = fi(ix + k, iy + k, ng)

            if (dist11 == 0.0_real64 .and. fw(1, 1) /= xmsg) then
                out(ng) = fw(1, 1)
                cycle
            end if
            if (dist21 == 0.0_real64 .and. fw(2, 1) /= xmsg) then
                out(ng) = fw(2, 1)
                cycle
            end if
            if (dist12 == 0.0_real64 .and. fw(1, 2) /= xmsg) then
                out(ng) = fw(1, 2)
                cycle
            end if
            if (dist22 == 0.0_real64 .and. fw(2, 2) /= xmsg) then
                out(ng) = fw(2, 2)
                cycle
            end if

            nw = 0
            sumf = 0.0_real64
            sumw = 0.0_real64
            do n = 1, 2
                do m = 1, 2
                    if (fw(m, n) /= xmsg) then
                        sumf = sumf + fw(m, n) * weights(m, n)
                        sumw = sumw + weights(m, n)
                        nw = nw + 1
                    end if
                end do
            end do

            if (nw >= ncrit .and. sumw > 0.0_real64) then
                out(ng) = sumf / sumw
            end if
        end do
    end subroutine interpolate_curv_cell

    ! Fill unresolved holes along each output row using the same semantics as
    ! DLINMSG(..., MFLAG=0, MPTCRT=2), but with a dedicated short-gap kernel
    ! because this cleanup pass only needs to patch interior runs of length
    ! one or two.
    subroutine fill_missing_rows(fo, xmsg)
        real(real64), intent(inout) :: fo(:, :, :)
        real(real64), intent(in) :: xmsg
        integer :: ng, ny

        do ng = 1, size(fo, 3)
            do ny = 1, size(fo, 2)
                call dlinmsg_interior_short_gaps(fo(:, ny, ng), size(fo, 1), xmsg)
            end do
        end do
    end subroutine fill_missing_rows

    ! Search the local regular-grid neighborhood for an exact coordinate hit.
    subroutine find_exact_regular_local(xi, yi, xo, yo, eps, ix_anchor, iy_anchor, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:)
        real(real64), intent(in) :: yi(:)
        real(real64), intent(in) :: xo, yo, eps
        integer, intent(in) :: ix_anchor, iy_anchor
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix_lo, ix_hi, iy_lo, iy_hi, ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        ix_lo = max(1, ix_anchor)
        ix_hi = min(size(xi), ix_anchor + 1)
        iy_lo = max(1, iy_anchor)
        iy_hi = min(size(yi), iy_anchor + 1)

        do iy = iy_lo, iy_hi
            do ix = ix_lo, ix_hi
                if (xo >= (xi(ix) - eps) .and. xo <= (xi(ix) + eps) .and. &
                    yo >= (yi(iy) - eps) .and. yo <= (yi(iy) + eps)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_exact_regular_local

    ! Full-grid exact-match search for the regular-grid source path.
    subroutine find_exact_regular_full(xi, yi, xo, yo, eps, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:)
        real(real64), intent(in) :: yi(:)
        real(real64), intent(in) :: xo, yo, eps
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = 1, size(yi)
            do ix = 1, size(xi)
                if (xo >= (xi(ix) - eps) .and. xo <= (xi(ix) + eps) .and. &
                    yo >= (yi(iy) - eps) .and. yo <= (yi(iy) + eps)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_exact_regular_full

    ! Brute-force search for the regular-grid cell that contains one output
    ! point. This only runs when the monotonic fast path cannot resolve it.
    subroutine find_regular_cell_full(xi, yi, xo, yo, found, ix_hit, iy_hit)
        real(real64), intent(in) :: xi(:)
        real(real64), intent(in) :: yi(:)
        real(real64), intent(in) :: xo, yo
        logical, intent(out) :: found
        integer, intent(out) :: ix_hit, iy_hit
        integer :: ix, iy

        found = .false.
        ix_hit = 1
        iy_hit = 1

        do iy = 1, size(yi) - 1
            do ix = 1, size(xi) - 1
                if (xo >= xi(ix) .and. xo < xi(ix + 1) .and. &
                    yo >= yi(iy) .and. yo < yi(iy + 1)) then
                    ix_hit = ix
                    iy_hit = iy
                    found = .true.
                    return
                end if
            end do
        end do
    end subroutine find_regular_cell_full

    ! Interpolate one regular-grid cell. The complete-data case uses bilinear
    ! interpolation; missing-corner cases fall back to inverse-distance
    ! weighting so the historical missing-value behavior is preserved.
    subroutine interpolate_regular_cell(fi, xi, yi, xo, yo, ix, iy, xmsg, ncrit, out)
        real(real64), intent(in) :: fi(:, :, :)
        real(real64), intent(in) :: xi(:)
        real(real64), intent(in) :: yi(:)
        real(real64), intent(in) :: xo, yo, xmsg
        integer, intent(in) :: ix, iy, ncrit
        real(real64), intent(inout) :: out(size(fi, 3))
        real(real64) :: fw(2, 2), weights(2, 2), sumf, sumw
        real(real64) :: slope_x, slope_y
        real(real64) :: dist11, dist21, dist12, dist22
        integer :: ng, m, n, nw

        do ng = 1, size(fi, 3)
            if (out(ng) /= xmsg) cycle

            fw(1, 1) = fi(ix, iy, ng)
            fw(2, 1) = fi(ix + 1, iy, ng)
            fw(1, 2) = fi(ix, iy + 1, ng)
            fw(2, 2) = fi(ix + 1, iy + 1, ng)

            if (fw(1, 1) /= xmsg .and. fw(2, 1) /= xmsg .and. &
                fw(1, 2) /= xmsg .and. fw(2, 2) /= xmsg) then
                slope_x = (xo - xi(ix)) / (xi(ix + 1) - xi(ix))
                slope_y = (yo - yi(iy)) / (yi(iy + 1) - yi(iy))
                out(ng) = (fw(1, 1) + slope_x * (fw(2, 1) - fw(1, 1))) + &
                          slope_y * ((fw(1, 2) + slope_x * (fw(2, 2) - fw(1, 2))) - &
                                     (fw(1, 1) + slope_x * (fw(2, 1) - fw(1, 1))))
                cycle
            end if

            dist11 = dgcdist_core(yo, xo, yi(iy), xi(ix), 2)
            dist21 = dgcdist_core(yo, xo, yi(iy), xi(ix + 1), 2)
            dist12 = dgcdist_core(yo, xo, yi(iy + 1), xi(ix), 2)
            dist22 = dgcdist_core(yo, xo, yi(iy + 1), xi(ix + 1), 2)

            weights(1, 1) = merge(huge(1.0_real64), (1.0_real64 / dist11) ** 2, dist11 == 0.0_real64)
            weights(2, 1) = merge(huge(1.0_real64), (1.0_real64 / dist21) ** 2, dist21 == 0.0_real64)
            weights(1, 2) = merge(huge(1.0_real64), (1.0_real64 / dist12) ** 2, dist12 == 0.0_real64)
            weights(2, 2) = merge(huge(1.0_real64), (1.0_real64 / dist22) ** 2, dist22 == 0.0_real64)

            nw = 0
            sumf = 0.0_real64
            sumw = 0.0_real64
            do n = 1, 2
                do m = 1, 2
                    if (fw(m, n) /= xmsg) then
                        sumf = sumf + fw(m, n) * weights(m, n)
                        sumw = sumw + weights(m, n)
                        nw = nw + 1
                    end if
                end do
            end do

            if (ncrit >= 1 .and. sumw > 0.0_real64) then
                out(ng) = sumf / sumw
            end if
        end do
    end subroutine interpolate_regular_cell

end module rcm2rgrid_kernels_core


! QUICK REFERENCE
! PURPOSE
!    INTERPOLATE A CURVILINEAR SOURCE GRID TO A RECTILINEAR TARGET GRID.
!    THIS IS THE MODERN FREE-FORM REPLACEMENT FOR THE HISTORICAL
!    `archive/rcm2rgrid.f` FIXED-FORM KERNEL, WHILE PRESERVING ITS
!    PUBLIC ENTRY POINT NAME.
!
! EXPECTED INPUT SHAPES
!    YI(NXI,NYI)      - SOURCE LATITUDES ON THE CURVILINEAR GRID
!    XI(NXI,NYI)      - SOURCE LONGITUDES ON THE CURVILINEAR GRID
!    FI(NXI,NYI,NGRD) - SOURCE DATA FIELDS
!    YO(NYO)          - TARGET LATITUDES ON THE RECTILINEAR GRID
!    XO(NXO)          - TARGET LONGITUDES ON THE RECTILINEAR GRID
!
! FLAGS
!    NCRIT  - MINIMUM NUMBER OF VALID CELL CORNERS REQUIRED FOR A VALUE
!             1 THROUGH 4, CLIPPED INTERNALLY TO THAT RANGE
!    OPT    - RESERVED FOR LEGACY API COMPATIBILITY
!    XMSG   - MISSING VALUE SENTINEL
!
! OUTPUT
!    FO(NXO,NYO,NGRD) HOLDS THE INTERPOLATED RECTILINEAR RESULT.
!    IER=0 ON SUCCESS, IER=1 WHEN AN INPUT GRID IS TOO SMALL.
subroutine drcm2rgrid(ngrd, nyi, nxi, yi, xi, fi, nyo, yo, nxo, xo, fo, xmsg, ncrit, opt, ier)
    use rcm2rgrid_kernels_core, only : real64, is_strictly_increasing, build_stride_candidates_increasing, &
        find_exact_curv_local, find_exact_curv_full, find_curv_cell_local, find_curv_cell_full, &
        interpolate_curv_cell, fill_missing_rows
    implicit none

    integer, intent(in) :: ngrd, nyi, nxi, nyo, nxo, ncrit, opt
    integer, intent(out) :: ier
    real(real64), intent(in) :: yi(nxi, nyi), xi(nxi, nyi), fi(nxi, nyi, ngrd)
    real(real64), intent(in) :: yo(nyo), xo(nxo), xmsg
    real(real64), intent(out) :: fo(nxo, nyo, ngrd)

    integer :: ix_hit, iy_hit, k, ncrit_use, nx, ny
    integer :: x_anchor(nxo), x_cell_hi(nxo), x_cell_lo(nxo), x_exact_hi(nxo), x_exact_lo(nxo)
    integer :: y_anchor(nyo), y_cell_hi(nyo), y_cell_lo(nyo), y_exact_hi(nyo), y_exact_lo(nyo)
    real(real64) :: x_ref(nxi), y_ref(nyi)
    real(real64), parameter :: eps = 1.0e-4_real64
    logical :: exact_found, found, use_fast_path

    ier = 0
    if (nxi <= 1 .or. nyi <= 1 .or. nxo <= 1 .or. nyo <= 1) then
        ier = 1
        return
    end if

    k = 2
    ncrit_use = max(1, min(4, ncrit))
    fo = xmsg

    ! The fast path uses monotonic reference slices from the first row/column
    ! to seed a narrow local search window around each output point.
    x_ref = xi(:, 1)
    y_ref = yi(1, :)
    use_fast_path = is_strictly_increasing(x_ref) .and. is_strictly_increasing(y_ref) .and. &
                    is_strictly_increasing(xo) .and. is_strictly_increasing(yo)

    if (use_fast_path) then
        call build_stride_candidates_increasing(x_ref, xo, k, x_anchor)
        call build_stride_candidates_increasing(y_ref, yo, k, y_anchor)
        do nx = 1, nxo
            x_exact_lo(nx) = max(1, x_anchor(nx) - k)
            x_exact_hi(nx) = min(nxi, x_anchor(nx) + 2 * k)
            x_cell_lo(nx) = max(1, x_anchor(nx) - 1)
            x_cell_hi(nx) = min(nxi - k, x_anchor(nx) + 2)
        end do
        do ny = 1, nyo
            y_exact_lo(ny) = max(1, y_anchor(ny) - k)
            y_exact_hi(ny) = min(nyi, y_anchor(ny) + 2 * k)
            y_cell_lo(ny) = max(1, y_anchor(ny) - 1)
            y_cell_hi(ny) = min(nyi - k, y_anchor(ny) + 2)
        end do
    else
        x_anchor = 1
        y_anchor = 1
        x_exact_lo = 1
        x_exact_hi = 1
        x_cell_lo = 1
        x_cell_hi = 1
        y_exact_lo = 1
        y_exact_hi = 1
        y_cell_lo = 1
        y_cell_hi = 1
    end if

    do ny = 1, nyo
        do nx = 1, nxo
            exact_found = .false.
            if (use_fast_path) then
                call find_exact_curv_local( &
                    xi, yi, xo(nx), yo(ny), eps, &
                    x_exact_lo(nx), x_exact_hi(nx), y_exact_lo(ny), y_exact_hi(ny), &
                    exact_found, ix_hit, iy_hit &
                )
            end if

            if (exact_found) then
                fo(nx, ny, :) = fi(ix_hit, iy_hit, :)
                if (all(fo(nx, ny, :) /= xmsg)) cycle
            end if

            found = .false.
            if (use_fast_path) then
                call find_curv_cell_local( &
                    xi, yi, xo(nx), yo(ny), &
                    x_cell_lo(nx), x_cell_hi(nx), y_cell_lo(ny), y_cell_hi(ny), k, &
                    found, ix_hit, iy_hit &
                )
            end if

            if (.not. found) then
                if (.not. exact_found) then
                    call find_exact_curv_full(xi, yi, xo(nx), yo(ny), eps, exact_found, ix_hit, iy_hit)
                    if (exact_found) then
                        fo(nx, ny, :) = fi(ix_hit, iy_hit, :)
                        if (all(fo(nx, ny, :) /= xmsg)) cycle
                    end if
                end if

                call find_curv_cell_full(xi, yi, xo(nx), yo(ny), k, found, ix_hit, iy_hit)
            end if

            if (found) then
                call interpolate_curv_cell(fi, xi, yi, xo(nx), yo(ny), ix_hit, iy_hit, k, xmsg, ncrit_use, &
                                           fo(nx, ny, :))
            end if
        end do
    end do

    call fill_missing_rows(fo, xmsg)
end subroutine drcm2rgrid


! QUICK REFERENCE
! PURPOSE
!    INTERPOLATE A RECTILINEAR SOURCE GRID TO A CURVILINEAR TARGET GRID.
!    THIS KEEPS THE LEGACY `drgrid2rcm` ENTRY POINT WHILE USING MODERN
!    FREE-FORM FORTRAN AND A MONOTONIC-GRID ANCHOR FAST PATH.
!
! EXPECTED INPUT SHAPES
!    YI(NYI)          - SOURCE LATITUDES ON THE RECTILINEAR GRID
!    XI(NXI)          - SOURCE LONGITUDES ON THE RECTILINEAR GRID
!    FI(NXI,NYI,NGRD) - SOURCE DATA FIELDS
!    YO(NXO,NYO)      - TARGET LATITUDES ON THE CURVILINEAR GRID
!    XO(NXO,NYO)      - TARGET LONGITUDES ON THE CURVILINEAR GRID
!
! FLAGS
!    NCRIT  - MINIMUM NUMBER OF VALID CELL CORNERS REQUIRED FOR A VALUE
!    OPT    - RESERVED FOR LEGACY API COMPATIBILITY
!    XMSG   - MISSING VALUE SENTINEL
!
! OUTPUT
!    FO(NXO,NYO,NGRD) HOLDS THE INTERPOLATED CURVILINEAR RESULT.
!    IER=0 ON SUCCESS, IER=1 WHEN AN INPUT GRID IS TOO SMALL.
subroutine drgrid2rcm(ngrd, nyi, nxi, yi, xi, fi, nyo, nxo, yo, xo, fo, xmsg, ncrit, opt, ier)
    use rcm2rgrid_kernels_core, only : real64, is_strictly_increasing, lower_bracket_increasing, &
        find_exact_regular_local, find_exact_regular_full, find_regular_cell_full, interpolate_regular_cell
    implicit none

    integer, intent(in) :: ngrd, nyi, nxi, nyo, nxo, ncrit, opt
    integer, intent(out) :: ier
    real(real64), intent(in) :: yi(nyi), xi(nxi), fi(nxi, nyi, ngrd)
    real(real64), intent(in) :: yo(nxo, nyo), xo(nxo, nyo), xmsg
    real(real64), intent(out) :: fo(nxo, nyo, ngrd)

    integer :: nx, ny, ix_hit, iy_hit, ix_anchor, iy_anchor
    real(real64), parameter :: eps = 1.0e-3_real64
    logical :: exact_found, found, use_fast_path

    ier = 0
    if (nxi <= 1 .or. nyi <= 1 .or. nxo <= 1 .or. nyo <= 1) then
        ier = 1
        return
    end if

    fo = xmsg
    use_fast_path = is_strictly_increasing(xi) .and. is_strictly_increasing(yi)

    do ny = 1, nyo
        do nx = 1, nxo
            exact_found = .false.
            found = .false.

            if (use_fast_path) then
                ix_anchor = lower_bracket_increasing(xi, xo(nx, ny))
                iy_anchor = lower_bracket_increasing(yi, yo(nx, ny))
                call find_exact_regular_local(xi, yi, xo(nx, ny), yo(nx, ny), eps, ix_anchor, iy_anchor, &
                                              exact_found, ix_hit, iy_hit)
                if (exact_found) then
                    fo(nx, ny, :) = fi(ix_hit, iy_hit, :)
                    if (all(fo(nx, ny, :) /= xmsg)) cycle
                end if

                if (xo(nx, ny) >= xi(ix_anchor) .and. xo(nx, ny) < xi(ix_anchor + 1) .and. &
                    yo(nx, ny) >= yi(iy_anchor) .and. yo(nx, ny) < yi(iy_anchor + 1)) then
                    found = .true.
                    ix_hit = ix_anchor
                    iy_hit = iy_anchor
                end if
            end if

            if (.not. found) then
                if (.not. exact_found) then
                    call find_exact_regular_full(xi, yi, xo(nx, ny), yo(nx, ny), eps, exact_found, ix_hit, iy_hit)
                    if (exact_found) then
                        fo(nx, ny, :) = fi(ix_hit, iy_hit, :)
                        if (all(fo(nx, ny, :) /= xmsg)) cycle
                    end if
                end if

                call find_regular_cell_full(xi, yi, xo(nx, ny), yo(nx, ny), found, ix_hit, iy_hit)
            end if

            if (found) then
                call interpolate_regular_cell(fi, xi, yi, xo(nx, ny), yo(nx, ny), ix_hit, iy_hit, xmsg, ncrit, &
                                              fo(nx, ny, :))
            end if
        end do
    end do
end subroutine drgrid2rcm
