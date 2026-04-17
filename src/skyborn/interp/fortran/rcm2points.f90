subroutine drcm2points(ngrd, nyi, nxi, yi, xi, fi, nxyo, yo, xo, fo, xmsg, opt, ncrit, kval, ier)
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
    integer :: ncand, candidate
    integer, allocatable :: candidate_ix(:), candidate_iy(:)
    logical :: exact_hit(nxyo), exact_found, cell_found, point_has_missing
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
    ! point once a compatible exact hit or interpolation path succeeds.
    fo = xmsg
    exact_hit = .false.

    ! Exact-match pass:
    ! If an output point lands exactly on a source-grid coordinate, copy every
    ! field directly and skip later interpolation for that output point.
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
                    exact_hit(nxy) = .true.
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
        if (.not. exact_hit(nxy)) then
            xo_n = xo(nxy)
            yo_n = yo(nxy)
            cell_found = .false.

            do iy = 1, nyi - k
                do ix = 1, nxi - k
                    if (xo_n >= xi(ix, iy) .and. xo_n <= xi(ix + k, iy) .and. &
                        yo_n >= yi(ix, iy) .and. yo_n <= yi(ix, iy + k)) then

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

                        cell_found = .true.
                        exit
                    end if
                end do
                if (cell_found) exit
            end do
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
