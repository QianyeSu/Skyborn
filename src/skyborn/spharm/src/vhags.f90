! =============================================================================
!
!                  copyright (c) 1998 by UCAR
!
!       University Corporation for Atmospheric Research
!
!                      all rights reserved
!
!                         SPHEREPACK
!
! ... file vhags.f
!
!     Modernized to Fortran 2008 standard with SIMD optimizations.
!     - Converted to free-form source.
!     - All subroutines are external (global), without a containing module.
!     - Replaced all archaic control flow (GOTO, computed GOTO) with modern
!       constructs (SELECT CASE, DO loops).
!     - Added 'implicit none', intent attributes, and portable precision kinds.
!     - Used assumed-size arrays dimension(*) to maintain F77 compatibility.
!
!     CRITICAL BUG FIX & PERFORMANCE OPTIMIZATION:
!     - Implemented conditional initialization of output spectral coefficient
!       arrays (br, bi, cr, ci) in vhags1, mirroring the exact F77 logic to
!       prevent numerical instability from uninitialized memory.
!     - Applied !$OMP SIMD directives to inner loops for explicit vectorization.
!     - Utilized temporary scalar variables to improve cache efficiency.
!     - Corrected all subroutine names from the previous erroneous version.
!
! =============================================================================


!> @brief Main vector spherical harmonic analysis routine on a Gaussian grid.
subroutine vhags(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhags, lvhags, work, lwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw, mdab, ndab, lvhags, lwork
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(in) :: wvhags(*)
    real, intent(inout) :: work(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, mmax, idz, lzimn, lmn, idv, lnl, ist
    integer :: jw1, jw2, jw3, iw1, iw2, iw3, iw4

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if
    if (ityp < 0 .or. ityp > 8) then; ierror = 3; return; end if
    if (nt < 0) then; ierror = 4; return; end if

    imid = (nlat + 1) / 2
    if ((ityp <= 2 .and. idvw < nlat) .or. (ityp > 2 .and. idvw < imid)) then
        ierror = 5; return
    end if

    if (jdvw < nlon) then; ierror = 6; return; end if

    mmax = min(nlat, (nlon + 1) / 2)
    if (mdab < mmax) then; ierror = 7; return; end if
    if (ndab < nlat) then; ierror = 8; return; end if

    idz = (mmax * (nlat + nlat - mmax + 1)) / 2
    lzimn = idz * imid
    if (lvhags < lzimn + lzimn + nlon + 15) then; ierror = 9; return; end if

    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt * idv * nlon
    if (lwork < lnl + lnl + idv * nlon) then; ierror = 10; return; end if

    ierror = 0

    ! --- Workspace pointer setup ---
    lmn = nlat * (nlat + 1) / 2
    jw1 = 1
    jw2 = jw1 + imid * lmn
    jw3 = jw2 + imid * lmn

    ist = 0
    if (ityp <= 2) ist = imid
    iw1 = ist + 1
    iw2 = lnl + 1
    iw3 = iw2 + ist
    iw4 = iw2 + lnl

    ! --- Call the core computational routine ---
    call vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhags(jw1), wvhags(jw2), wvhags(jw3))

end subroutine vhags


!> @brief Core vector harmonic analysis computation routine on a Gaussian grid with SIMD.
subroutine vhags1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
                  ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, imid, idvw, jdvw, mdab, ndab, idv, idz
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(inout) :: ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), wo(idv, nlon, *)
    real, intent(inout) :: work(*)
    real, intent(in) :: vb(imid, *), wb(imid, *), wrfft(*)

    ! Local variables
    integer :: nlp1, mlat, mmax, imm1, ndo1, ndo2, itypp
    integer :: k, i, j, mp1, np1, m, mb, mp2
    real :: tsn, fsn
    real :: br_val_0, cr_val_0
    real :: br_val_m, bi_val_m, cr_val_m, ci_val_m

    ! --- Precompute constants ---
    nlp1 = nlat + 1
    tsn = 2.0 / real(nlon)
    fsn = 4.0 / real(nlon)
    mlat = mod(nlat, 2)
    mmax = min(nlat, (nlon + 1) / 2)
    imm1 = imid
    if (mlat /= 0) imm1 = imid - 1

    ! --- Decompose grid into even and odd components ---
    if (ityp > 2) then
        do k = 1, nt
            do j = 1, nlon
                do i = 1, imm1
                    ve(i, j, k) = fsn * v(i, j, k)
                    vo(i, j, k) = fsn * v(i, j, k)
                    we(i, j, k) = fsn * w(i, j, k)
                    wo(i, j, k) = fsn * w(i, j, k)
                end do
            end do
        end do
    else
        do k = 1, nt
            do j = 1, nlon
                do i = 1, imm1
                    ve(i, j, k) = tsn * (v(i, j, k) + v(nlp1 - i, j, k))
                    vo(i, j, k) = tsn * (v(i, j, k) - v(nlp1 - i, j, k))
                    we(i, j, k) = tsn * (w(i, j, k) + w(nlp1 - i, j, k))
                    wo(i, j, k) = tsn * (w(i, j, k) - w(nlp1 - i, j, k))
                end do
            end do
        end do
    end if

    if (mlat /= 0) then
        do k = 1, nt
            do j = 1, nlon
                ve(imid, j, k) = tsn * v(imid, j, k)
                we(imid, j, k) = tsn * w(imid, j, k)
            end do
        end do
    end if

    ! --- Perform Forward Fourier Transform ---
    do k = 1, nt
        call hrfftf(idv, nlon, ve(1, 1, k), idv, wrfft, work)
        call hrfftf(idv, nlon, we(1, 1, k), idv, wrfft, work)
    end do

    ndo1 = nlat
    ndo2 = nlat
    if (mlat /= 0) ndo1 = nlat - 1
    if (mlat == 0) ndo2 = nlat - 1

    ! --- CRITICAL FIX: Conditionally initialize spectral coefficient arrays ---
    if (ityp /= 2 .and. ityp /= 5 .and. ityp /= 8) then
        do k = 1, nt
            do mp1 = 1, mmax
                do np1 = mp1, nlat
                    br(mp1, np1, k) = 0.0
                    bi(mp1, np1, k) = 0.0
                end do
            end do
        end do
    end if
    if (ityp /= 1 .and. ityp /= 4 .and. ityp /= 7) then
        do k = 1, nt
            do mp1 = 1, mmax
                do np1 = mp1, nlat
                    cr(mp1, np1, k) = 0.0
                    ci(mp1, np1, k) = 0.0
                end do
            end do
        end do
    end if

    itypp = ityp + 1

    ! --- Perform Legendre transform via SELECT CASE on symmetry type ---
    select case (itypp)

    case (1) ! ityp=0: no symmetries
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                br_val_0 = br(1, np1, k); cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0, cr_val_0)
                do i = 1, imid
                    br_val_0 = br_val_0 + vb(i, np1) * ve(i, 1, k)
                    cr_val_0 = cr_val_0 - vb(i, np1) * we(i, 1, k)
                end do
                br(1, np1, k) = br_val_0; cr(1, np1, k) = cr_val_0
            end do
            do np1 = 3, ndo1, 2
                br_val_0 = br(1, np1, k); cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0, cr_val_0)
                do i = 1, imm1
                    br_val_0 = br_val_0 + vb(i, np1) * vo(i, 1, k)
                    cr_val_0 = cr_val_0 - vb(i, np1) * wo(i, 1, k)
                end do
                br(1, np1, k) = br_val_0; cr(1, np1, k) = cr_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m, cr_val_m, ci_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * vo(i, 2*mp1-2, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * vo(i, 2*mp1-1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                                cr_val_m = cr_val_m - vb(i, np1+mb) * wo(i, 2*mp1-2, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * wo(i, 2*mp1-1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m, cr_val_m, ci_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * ve(i, 2*mp1-2, k) + wb(i, np1+mb) * wo(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * ve(i, 2*mp1-1, k) - wb(i, np1+mb) * wo(i, 2*mp1-2, k)
                                cr_val_m = cr_val_m - vb(i, np1+mb) * we(i, 2*mp1-2, k) + wb(i, np1+mb) * vo(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * we(i, 2*mp1-1, k) - wb(i, np1+mb) * vo(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (2) ! ityp=1: no symmetries, cr=ci=0
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imid
                    br_val_0 = br_val_0 + vb(i, np1) * ve(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
            do np1 = 3, ndo1, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imm1
                    br_val_0 = br_val_0 + vb(i, np1) * vo(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * vo(i, 2*mp1-2, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * vo(i, 2*mp1-1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * ve(i, 2*mp1-2, k) + wb(i, np1+mb) * wo(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * ve(i, 2*mp1-1, k) - wb(i, np1+mb) * wo(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (3) ! ityp=2: no symmetries, br=bi=0
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imid
                    cr_val_0 = cr_val_0 - vb(i, np1) * we(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
            do np1 = 3, ndo1, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imm1
                    cr_val_0 = cr_val_0 - vb(i, np1) * wo(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * wo(i, 2*mp1-2, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * wo(i, 2*mp1-1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * we(i, 2*mp1-2, k) + wb(i, np1+mb) * vo(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * we(i, 2*mp1-1, k) - wb(i, np1+mb) * vo(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (4) ! ityp=3: v even, w odd
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imid
                    br_val_0 = br_val_0 + vb(i, np1) * ve(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
            do np1 = 3, ndo1, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imm1
                    cr_val_0 = cr_val_0 - vb(i, np1) * wo(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * wo(i, 2*mp1-2, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * wo(i, 2*mp1-1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * ve(i, 2*mp1-2, k) + wb(i, np1+mb) * wo(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * ve(i, 2*mp1-1, k) - wb(i, np1+mb) * wo(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (5) ! ityp=4: v even, w odd, cr=ci=0
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imid
                    br_val_0 = br_val_0 + vb(i, np1) * ve(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * ve(i, 2*mp1-2, k) + wb(i, np1+mb) * wo(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * ve(i, 2*mp1-1, k) - wb(i, np1+mb) * wo(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1+mb) * ve(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (6) ! ityp=5: v even, w odd, br=bi=0
        ! m=0
        do k = 1, nt
            do np1 = 3, ndo1, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imm1
                    cr_val_0 = cr_val_0 - vb(i, np1) * wo(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * wo(i, 2*mp1-2, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * wo(i, 2*mp1-1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1+mb) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1+mb) * ve(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (7) ! ityp=6: v odd, w even
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imid
                    cr_val_0 = cr_val_0 - vb(i, np1) * we(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
            do np1 = 3, ndo1, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imm1
                    br_val_0 = br_val_0 + vb(i, np1) * vo(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * vo(i, 2*mp1-2, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * vo(i, 2*mp1-1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * we(i, 2*mp1-2, k) + wb(i, np1+mb) * vo(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * we(i, 2*mp1-1, k) - wb(i, np1+mb) * vo(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (8) ! ityp=7: v odd, w even, cr=ci=0
        ! m=0
        do k = 1, nt
            do np1 = 3, ndo1, 2
                br_val_0 = br(1, np1, k)
                !$OMP SIMD REDUCTION(+:br_val_0)
                do i = 1, imm1
                    br_val_0 = br_val_0 + vb(i, np1) * vo(i, 1, k)
                end do
                br(1, np1, k) = br_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            br_val_m = br(mp1, np1, k); bi_val_m = bi(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:br_val_m, bi_val_m)
                            do i = 1, imm1
                                br_val_m = br_val_m + vb(i, np1+mb) * vo(i, 2*mp1-2, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi_val_m = bi_val_m + vb(i, np1+mb) * vo(i, 2*mp1-1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end do
                            br(mp1, np1, k) = br_val_m; bi(mp1, np1, k) = bi_val_m
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1+mb) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1+mb) * we(i, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (9) ! ityp=8: v odd, w even, br=bi=0
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                cr_val_0 = cr(1, np1, k)
                !$OMP SIMD REDUCTION(+:cr_val_0)
                do i = 1, imid
                    cr_val_0 = cr_val_0 - vb(i, np1) * we(i, 1, k)
                end do
                cr(1, np1, k) = cr_val_0
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                mb = m * nlat - (m * (m + 1)) / 2
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            cr_val_m = cr(mp1, np1, k); ci_val_m = ci(mp1, np1, k)
                            !$OMP SIMD REDUCTION(+:cr_val_m, ci_val_m)
                            do i = 1, imm1
                                cr_val_m = cr_val_m - vb(i, np1+mb) * we(i, 2*mp1-2, k) + wb(i, np1+mb) * vo(i, 2*mp1-1, k)
                                ci_val_m = ci_val_m - vb(i, np1+mb) * we(i, 2*mp1-1, k) - wb(i, np1+mb) * vo(i, 2*mp1-2, k)
                            end do
                            cr(mp1, np1, k) = cr_val_m; ci(mp1, np1, k) = ci_val_m
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1+mb) * we(i, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    end select

end subroutine vhags1


!> @brief Initialization routine for vector harmonic analysis on a Gaussian grid.
subroutine vhagsi(nlat, nlon, wvhags, lvhags, dwork, ldwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, lvhags, ldwork
    real, intent(out) :: wvhags(*)
    double precision, intent(inout) :: dwork(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, lwk, i
    integer :: jw1, jw2, jw3, iw1, iw2, iw3, iw4

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if

    imid = (nlat + 1) / 2
    if (lvhags < 2 * (imid * (nlat * (nlat+1)/2)) + nlon + 15) then; ierror = 3; return; end if

    if (ldwork < (nlat * (3 * nlat + 9) + 2) / 2) then; ierror = 4; return; end if
    ierror = 0

    ! --- Workspace pointer setup ---
    jw1 = 1
    jw2 = jw1 + imid * (nlat * (nlat + 1) / 2)
    jw3 = jw2 + imid * (nlat * (nlat + 1) / 2)

    iw1 = 1
    iw2 = iw1 + nlat
    iw3 = iw2 + nlat
    iw4 = iw3 + 3 * imid * nlat

    ! --- Call initialization core routine and FFT setup ---
    call vhgai1(nlat, imid, wvhags(jw1), wvhags(jw2), &
                dwork(iw1), dwork(iw2), dwork(iw3), dwork(iw4))

    call hrffti(nlon, wvhags(jw3))

end subroutine vhagsi


!> @brief Core initialization routine for vector basis functions with Gaussian weights.
subroutine vhgai1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, imid
    real, intent(out) :: vb(imid, *), wb(imid, *)
    double precision, intent(inout) :: dthet(*), dwts(*), dpbar(imid, nlat, 3), work(*)

    ! Local variables
    integer, parameter :: real64 = kind(1.0d0)
    integer :: lwk, ierror, i, n, nm, nz, np, m, ix, iy
    real(kind=real64) :: abel, bbel, cbel, ssqr2, dcf

    ! --- Compute Gauss points and weights ---
    lwk = nlat * (nlat + 2)
    call gaqd(nlat, dthet, dwts, dpbar, lwk, ierror)

    ! --- Compute associated Legendre functions and basis functions ---
    ssqr2 = 1.0_real64 / sqrt(2.0_real64)
    do i = 1, imid
        dpbar(i, 1, 1) = ssqr2
    end do
    vb(:, 1) = 0.0
    wb(:, 1) = 0.0

    do n = 1, nlat - 1
        nm = mod(n - 2, 3) + 1
        nz = mod(n - 1, 3) + 1
        np = mod(n, 3) + 1

        ! m=0
        call dnlfk(0, n, work)
        do i = 1, imid
            call dnlft(0, n, dthet(i), work, dpbar(i, 1, np))
        end do

        ! m=1
        call dnlfk(1, n, work)
        do i = 1, imid
            call dnlft(1, n, dthet(i), work, dpbar(i, 2, np))
        end do

        ! m=2 to n (recursive computation)
        if (n >= 2) then
            do m = 2, n
                abel = sqrt(real((2*n+1)*(m+n-2)*(m+n-3), kind=real64) / real((2*n-3)*(m+n-1)*(m+n), kind=real64))
                bbel = sqrt(real((2*n+1)*(n-m-1)*(n-m), kind=real64) / real((2*n-3)*(m+n-1)*(m+n), kind=real64))
                cbel = sqrt(real((n-m+1)*(n-m+2), kind=real64) / real((m+n-1)*(m+n), kind=real64))

                if (m < n - 1) then
                    do i = 1, imid
                        dpbar(i, m+1, np) = abel * dpbar(i, m-1, nm) + bbel * dpbar(i, m+1, nm) - cbel * dpbar(i, m-1, np)
                    end do
                else
                    do i = 1, imid
                        dpbar(i, m+1, np) = abel * dpbar(i, m-1, nm) - cbel * dpbar(i, m-1, np)
                    end do
                end if
            end do
        end if

        ! --- Compute and store vb and wb, including Gaussian weights ---
        ix = indx(0, n)
        iy = indx(n, n)
        do i = 1, imid
            vb(i, ix) = -dpbar(i, 2, np) * dwts(i)
            vb(i, iy) = dpbar(i, n, np) / sqrt(real(2*(n+1), kind=real64)) * dwts(i)
        end do

        if (n >= 2) then
            dcf = sqrt(real(4*n*(n+1), kind=real64))
            do m = 1, n - 1
                ix = indx(m, n)
                abel = sqrt(real((n+m)*(n-m+1), kind=real64)) / dcf
                bbel = sqrt(real((n-m)*(n+m+1), kind=real64)) / dcf
                do i = 1, imid
                    vb(i, ix) = (abel * dpbar(i, m, np) - bbel * dpbar(i, m+2, np)) * dwts(i)
                end do
            end do
        end if

        ix = indx(0, n)
        do i = 1, imid
            wb(i, ix) = 0.0_real64
        end do

        dcf = sqrt(real(2*n+1, kind=real64) / real(4*n*(n+1)*(2*n-1), kind=real64))
        do m = 1, n
            ix = indx(m, n)
            abel = dcf * sqrt(real((n+m)*(n+m-1), kind=real64))
            bbel = dcf * sqrt(real((n-m)*(n-m-1), kind=real64))
            if (m < n - 1) then
                do i = 1, imid
                    wb(i, ix) = (abel * dpbar(i, m, nz) + bbel * dpbar(i, m+2, nz)) * dwts(i)
                end do
            else
                do i = 1, imid
                    wb(i, ix) = abel * dpbar(i, m, nz) * dwts(i)
                end do
            end if
        end do
    end do

contains

    !> @brief Internal function to compute 1D index from (m,n).
    !> This replaces the F77 statement function for clarity and safety.
    function indx(m, n)
        integer, intent(in) :: m, n
        integer :: indx
        indx = m * nlat - (m * (m + 1)) / 2 + n + 1
    end function indx

end subroutine vhgai1
