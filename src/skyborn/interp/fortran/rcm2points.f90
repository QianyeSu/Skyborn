module rcm2points_core
    implicit none
contains

    ! Check whether a 1-D coordinate is strictly increasing so the curvilinear
    ! point-interpolation kernel can safely seed a narrow local search window.
    pure logical function is_strictly_increasing(x) result(ok)
        double precision, intent(in) :: x(:)
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
        double precision, intent(in) :: x(:), value
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
    ! still span a stride-k interpolation box.
    pure integer function earliest_stride_candidate(x, value, stride) result(idx)
        double precision, intent(in) :: x(:), value
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

    ! Search a small local curvilinear neighborhood for the first source cell
    ! whose corner coordinates bracket the requested output point.
    subroutine find_curv_cell_local(xi, yi, xo, yo, ix_lo, ix_hi, iy_lo, iy_hi, k, found, ix_hit, iy_hit)
        double precision, intent(in) :: xi(:, :)
        double precision, intent(in) :: yi(:, :)
        double precision, intent(in) :: xo, yo
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
    ! monotonic anchor heuristic cannot identify the cell locally.
    subroutine find_curv_cell_full(xi, yi, xo, yo, k, found, ix_hit, iy_hit)
        double precision, intent(in) :: xi(:, :)
        double precision, intent(in) :: yi(:, :)
        double precision, intent(in) :: xo, yo
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

end module rcm2points_core


subroutine drcm2points(ngrd, nyi, nxi, yi, xi, fi, nxyo, yo, xo, fo, xmsg, opt, ncrit, kval, ier)
    use rcm2points_core, only : is_strictly_increasing, earliest_stride_candidate, find_curv_cell_local, &
        find_curv_cell_full
    implicit none

    ! QUICK REFERENCE
    ! PURPOSE
    !    INTERPOLATE ONE OR MORE FIELDS FROM A CURVILINEAR SOURCE GRID TO
    !    AN UNSTRUCTURED LIST OF OUTPUT POINTS. THE HISTORICAL KERNEL FIRST
    !    LOOKS FOR EXACT SOURCE-POINT HITS, THEN TRIES A SINGLE LOCAL CELL
    !    INTERPOLATION, AND FINALLY FALLS BACK TO AN INVERSE-DISTANCE SEARCH
    !    IN A FIXED RADIUS.
    !
    ! INPUTS
    !    XI(NXI, NYI), YI(NXI, NYI)
    !                 SOURCE LONGITUDE / LATITUDE COORDINATES
    !    FI(NXI, NYI, NGRD)
    !                 SOURCE FIELDS TO INTERPOLATE
    !    XO(NXYO), YO(NXYO)
    !                 OUTPUT LONGITUDE / LATITUDE POINTS
    !    XMSG       - MISSING VALUE SENTINEL
    !    OPT        - 0/1: INVERSE-DISTANCE CELL WEIGHTS
    !                 2  : BILINEAR CELL WEIGHTS
    !    NCRIT      - MINIMUM NUMBER OF VALID CELL CORNERS NEEDED FOR THE
    !                 LOCAL CELL INTERPOLATION TO SUCCEED
    !    KVAL       - STEP SIZE USED WHEN FORMING THE LOCAL 2 X 2 CELL.
    !
    ! OUTPUTS
    !    FO(NXYO, NGRD)
    !                 INTERPOLATED VALUES AT EACH OUTPUT POINT
    !    IER        - 0 ON SUCCESS
    !                 1 IF THE INPUT/OUTPUT GRIDS ARE TOO SMALL
    !
    ! NOTES
    !    THIS FREE-FORM TRANSLATION PRESERVES THE HISTORICAL SENTINEL
    !    COMPARISONS AND SEARCH ORDER SO OUTPUTS STAY LEGACY-COMPATIBLE.

    integer, intent(in) :: ngrd, nyi, nxi, nxyo, opt, ncrit, kval
    integer, intent(out) :: ier
    integer :: ng, nx, ny, nxy, ix, iy, m, n, nw, ner, k
    integer :: ncand, candidate, ix_guess, iy_guess
    integer, allocatable :: candidate_ix(:), candidate_iy(:)
    integer, allocatable :: x_cell_lo(:), x_cell_hi(:), y_cell_lo(:), y_cell_hi(:)
    logical :: exact_found, cell_found, point_has_missing, use_fast_path
    double precision, intent(in) :: xi(nxi, nyi), yi(nxi, nyi), fi(nxi, nyi, ngrd)
    double precision, intent(in) :: xo(nxyo), yo(nxyo), xmsg
    double precision, intent(out) :: fo(nxyo, ngrd)
    double precision :: fw(2, 2), w(2, 2), sumf, sumw, chklat(nyi), chklon(nxi)
    double precision, allocatable :: candidate_weight(:)
    double precision :: wx, wy, rearth, dlat, pi, rad, dkm, dist, xo_n, yo_n
    double precision :: dgcdist

    external dmonoinc

    ! Match the historical input-size guard before any monotonicity or
    ! interpolation work begins.
    ier = 0
    if (nxi <= 1 .or. nyi <= 1 .or. nxyo <= 0) then
        ier = 1
        return
    end if

    ! The legacy kernel validates monotonicity using the first latitude row and
    ! first longitude column as representative 1-D coordinates.
    do ny = 1, nyi
        chklat(ny) = yi(1, ny)
    end do
    call dmonoinc(chklat, nyi, ier, ner)
    if (ier /= 0) return

    do nx = 1, nxi
        chklon(nx) = xi(nx, 1)
    end do
    ! Preserve the historical validation sequence for compatibility. The
    ! fixed-form kernel checked CHKLAT twice here instead of CHKLON.
    call dmonoinc(chklat, nyi, ier, ner)
    if (ier /= 0) return

    if (kval <= 0) then
        k = 1
    else
        k = kval
    end if

    ! Start with every output set to missing. Later passes only overwrite a
    ! field once a compatible exact hit or interpolation path succeeds.
    fo = xmsg

    ! When the representative 1-D coordinate slices are monotonic, seed a
    ! small local search window for the cell lookup. A full-grid fallback keeps
    ! the historical behavior for folded or unexpected grids.
    use_fast_path = nxi > k .and. nyi > k
    if (use_fast_path) then
        use_fast_path = is_strictly_increasing(chklon) .and. is_strictly_increasing(chklat)
    end if

    if (use_fast_path) then
        allocate(x_cell_lo(nxyo), x_cell_hi(nxyo), y_cell_lo(nxyo), y_cell_hi(nxyo))
        do nxy = 1, nxyo
            ix_guess = earliest_stride_candidate(chklon, xo(nxy), k)
            iy_guess = earliest_stride_candidate(chklat, yo(nxy), k)

            x_cell_lo(nxy) = max(1, ix_guess - 1)
            x_cell_hi(nxy) = min(nxi - k, ix_guess + 2)
            y_cell_lo(nxy) = max(1, iy_guess - 1)
            y_cell_hi(nxy) = min(nyi - k, iy_guess + 2)
        end do
    end if

    ! Exact-match pass:
    ! If an output point lands exactly on a source-grid coordinate, copy every
    ! field directly. Missing exact-hit fields are still allowed to fall through
    ! to the cell interpolation pass below, matching the historical F77 kernel.
    do nxy = 1, nxyo
        xo_n = xo(nxy)
        yo_n = yo(nxy)
        exact_found = .false.

        do iy = 1, nyi
            do ix = 1, nxi
                if (xo_n == xi(ix, iy) .and. yo_n == yi(ix, iy)) then
                    do ng = 1, ngrd
                        fo(nxy, ng) = fi(ix, iy, ng)
                    end do
                    exact_found = .true.
                    exit
                end if
            end do
            if (exact_found) exit
        end do
    end do

    ! Local cell pass:
    ! Search the first stride-k source box whose corner coordinates bracket the
    ! requested output point, then interpolate within that 2 x 2 cell using
    ! either bilinear weights or inverse-distance weights.
    do nxy = 1, nxyo
        point_has_missing = .false.
        do ng = 1, ngrd
            if (fo(nxy, ng) == xmsg) then
                point_has_missing = .true.
                exit
            end if
        end do

        if (point_has_missing) then
            xo_n = xo(nxy)
            yo_n = yo(nxy)
            cell_found = .false.

            if (use_fast_path) then
                call find_curv_cell_local( &
                    xi, yi, xo_n, yo_n, x_cell_lo(nxy), x_cell_hi(nxy), y_cell_lo(nxy), y_cell_hi(nxy), k, &
                    cell_found, ix, iy &
                )
            end if

            if (.not. cell_found) then
                call find_curv_cell_full(xi, yi, xo_n, yo_n, k, cell_found, ix, iy)
            end if

            if (cell_found) then
                if (abs(opt) == 2) then
                    wx = (xo_n - xi(ix, iy)) / (xi(ix + k, iy) - xi(ix, iy))
                    wy = (yo_n - yi(ix, iy)) / (yi(ix, iy + k) - yi(ix, iy))
                    w(1, 1) = (1.0d0 - wx) * (1.0d0 - wy)
                    w(2, 1) = wx * (1.0d0 - wy)
                    w(1, 2) = (1.0d0 - wx) * wy
                    w(2, 2) = wx * wy
                else
                    w(1, 1) = (1.0d0 / dgcdist(yo_n, xo_n, yi(ix, iy), xi(ix, iy), 2)) ** 2
                    w(2, 1) = (1.0d0 / dgcdist(yo_n, xo_n, yi(ix + k, iy), xi(ix + k, iy), 2)) ** 2
                    w(1, 2) = (1.0d0 / dgcdist(yo_n, xo_n, yi(ix, iy + k), xi(ix, iy + k), 2)) ** 2
                    w(2, 2) = (1.0d0 / dgcdist(yo_n, xo_n, yi(ix + k, iy + k), xi(ix + k, iy + k), 2)) ** 2
                end if

                do ng = 1, ngrd
                    if (fo(nxy, ng) == xmsg) then
                        fw(1, 1) = fi(ix, iy, ng)
                        fw(2, 1) = fi(ix + k, iy, ng)
                        fw(1, 2) = fi(ix, iy + k, ng)
                        fw(2, 2) = fi(ix + k, iy + k, ng)

                        nw = 0
                        sumf = 0.0d0
                        sumw = 0.0d0
                        do n = 1, 2
                            do m = 1, 2
                                if (fw(m, n) /= xmsg) then
                                    sumf = sumf + fw(m, n) * w(m, n)
                                    sumw = sumw + w(m, n)
                                    nw = nw + 1
                                end if
                            end do
                        end do

                        if (nw >= ncrit .and. sumw > 0.0d0) then
                            fo(nxy, ng) = sumf / sumw
                        end if
                    end if
                end do
            end if
        end if
    end do

    ! If every field at every output point has already been filled, the
    ! expensive radius-search fallback is unnecessary.
    do ng = 1, ngrd
        do nxy = 1, nxyo
            if (fo(nxy, ng) == xmsg) exit
        end do
        if (nxy <= nxyo) exit
    end do
    if (ng > ngrd) return

    rearth = 6371.0d0
    dlat = 5.0d0
    pi = 4.0d0 * atan(1.0d0)
    rad = pi / 180.0d0
    dkm = dlat * (2.0d0 * pi * rearth) / 360.0d0

    allocate(candidate_ix(nxi * nyi), candidate_iy(nxi * nyi), candidate_weight(nxi * nyi))

    ! For unresolved output points, reuse the geometry-only fallback search
    ! across all fields instead of rescanning the grid for each NGRD slice.
    do nxy = 1, nxyo
        point_has_missing = .false.
        do ng = 1, ngrd
            if (fo(nxy, ng) == xmsg) then
                point_has_missing = .true.
                exit
            end if
        end do

        if (point_has_missing) then
            xo_n = xo(nxy)
            yo_n = yo(nxy)
            ncand = 0

            ! Build the candidate list once for this output point using the
            ! same fixed-radius latitude window and great-circle cutoff as the
            ! historical fallback branch.
            do iy = 1, nyi
                do ix = 1, nxi
                    if (yi(ix, iy) >= yo_n - dlat .and. yi(ix, iy) <= yo_n + dlat) then
                        dist = dgcdist(yo_n, xo_n, yi(ix, iy), xi(ix, iy), 2)
                        if (dist <= dkm .and. dist > 0.0d0) then
                            ncand = ncand + 1
                            candidate_ix(ncand) = ix
                            candidate_iy(ncand) = iy
                            candidate_weight(ncand) = 1.0d0 / dist ** 2
                        end if
                    end if
                end do
            end do

            if (ncand > 0) then
                ! Each field still applies its own missing-value mask to the
                ! shared geometry candidates, so per-field results remain
                ! identical to the legacy NG -> NXY -> IX/IY traversal.
                do ng = 1, ngrd
                    if (fo(nxy, ng) == xmsg) then
                        nw = 0
                        sumf = 0.0d0
                        sumw = 0.0d0

                        do candidate = 1, ncand
                            ix = candidate_ix(candidate)
                            iy = candidate_iy(candidate)
                            if (fi(ix, iy, ng) /= xmsg) then
                                sumf = sumf + fi(ix, iy, ng) * candidate_weight(candidate)
                                sumw = sumw + candidate_weight(candidate)
                                nw = nw + 1
                            end if
                        end do

                        if (sumw > 0.0d0) then
                            fo(nxy, ng) = sumf / sumw
                        end if
                    end if
                end do
            end if
        end if
    end do

    deallocate(candidate_ix, candidate_iy, candidate_weight)
end subroutine drcm2points
