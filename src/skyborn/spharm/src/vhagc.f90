! =============================================================================
!
!                  copyright (c) 2025 by Qianye Su
!
!       University Corporation for Atmospheric Research
!
!                      all rights reserved
!
!                         SPHEREPACK
!
! ... file vhagc.f
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
!     - Replicated the exact conditional initialization logic for the output
!       spectral coefficient arrays (br, bi, cr, ci) from the original
!       F77 code. This prevents catastrophic numerical errors from uninitialized memory.
!     - Applied !$OMP SIMD directives to all computationally intensive inner
!       loops to explicitly enable vectorization for significant speedup.
!     - Used temporary scalar variables within loops to promote register usage
!       and reduce memory access, enhancing SIMD efficiency.
!
! =============================================================================


!> @brief Main vector spherical harmonic analysis routine on a Gaussian grid.
subroutine vhagc(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhagc, lvhagc, work, lwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw, mdab, ndab, lvhagc, lwork
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(in) :: wvhagc(*)
    real, intent(inout) :: work(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, mmax, lzz1, labc, idv, lnl, ist
    integer :: iw1, iw2, iw3, iw4, iw5, lwzvin, jw1, jw2, jw3

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

    lzz1 = 2 * nlat * imid
    labc = 3 * (max(mmax - 2, 0) * (nlat + nlat - mmax - 1)) / 2
    if (lvhagc < 2 * (lzz1 + labc) + nlon + imid + 15) then; ierror = 9; return; end if

    if (ityp <= 2 .and. lwork < nlat * (4 * nlon * nt + 6 * imid)) then; ierror = 10; return; end if
    if (ityp > 2 .and. lwork < imid * (4 * nlon * nt + 6 * nlat)) then; ierror = 10; return; end if

    ierror = 0

    ! --- Workspace pointer setup ---
    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt * idv * nlon
    ist = 0
    if (ityp <= 2) ist = imid

    iw1 = ist + 1
    iw2 = lnl + 1
    iw3 = iw2 + ist
    iw4 = iw2 + lnl
    iw5 = iw4 + 3 * imid * nlat

    lwzvin = lzz1 + labc
    jw1 = (nlat + 1) / 2 + 1
    jw2 = jw1 + lwzvin
    jw3 = jw2 + lwzvin

    ! --- Call the core computational routine ---
    call vhagc1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), work(iw5), wvhagc, wvhagc(jw1), wvhagc(jw2), wvhagc(jw3))

end subroutine vhagc


!> @brief Core vector harmonic analysis computation routine on a Gaussian grid with SIMD
subroutine vhagc1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
                  ndab, br, bi, cr, ci, idv, ve, vo, we, wo, vb, wb, wts, wvbin, wwbin, wrfft)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, imid, idvw, jdvw, mdab, ndab, idv
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(inout) :: ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), wo(idv, nlon, *)
    real, intent(inout) :: vb(imid, nlat, 3), wb(imid, nlat, 3)
    real, intent(in) :: wts(*), wvbin(*), wwbin(*), wrfft(*)

    ! Local variables
    integer :: nlp1, mlat, mmax, imm1, ndo1, ndo2, itypp
    integer :: k, i, j, mp1, np1, m, mp2, iv, iw
    real :: tsn, fsn, tv, tw
    real :: tvo1, tvo2, tve1, tve2, two1, two2, twe1, twe2

    ! --- Precompute constants ---
    nlp1 = nlat + 1
    tsn = 2.0 / real(nlon)
    fsn = 4.0 / real(nlon)
    mlat = mod(nlat, 2)
    mmax = min(nlat, (nlon + 1) / 2)
    imm1 = imid
    if (mlat /= 0) imm1 = imid - 1

    ! --- Decompose grid into even and odd components about the equator ---
    if (ityp <= 2) then
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
    else
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
        call hrfftf(idv, nlon, ve(1, 1, k), idv, wrfft, vb)
        call hrfftf(idv, nlon, we(1, 1, k), idv, wrfft, vb)
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
        call vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tv = ve(i, 1, k) * wts(i)
                    tw = we(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tv = vo(i, 1, k) * wts(i)
                    tw = wo(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mp2 = mp1 + 1
                call vbin(0, nlat, nlon, m, vb, iv, wvbin)
                call wbin(0, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tvo2 + wb(i, np1, iw) * twe1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tvo1 - wb(i, np1, iw) * twe2
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * two2 + wb(i, np1, iw) * tve1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * two1 - wb(i, np1, iw) * tve2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1, iw) * we(i, 2*mp1-1, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1, iw) * we(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1, iw) * ve(i, 2*mp1-1, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1, iw) * ve(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tve2 + wb(i, np1, iw) * two1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tve1 - wb(i, np1, iw) * two2
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * twe2 + wb(i, np1, iw) * tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * twe1 - wb(i, np1, iw) * tvo2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-2, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-1, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-2, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (2) ! ityp=1: no symmetries, cr=ci=0
        call vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tv = ve(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tv = vo(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(0, nlat, nlon, m, vb, iv, wvbin)
                call wbin(0, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tvo2 + wb(i, np1, iw) * twe1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tvo1 - wb(i, np1, iw) * twe2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1, iw) * we(i, 2*mp1-1, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1, iw) * we(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tve2 + wb(i, np1, iw) * two1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tve1 - wb(i, np1, iw) * two2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-2, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (3) ! ityp=2: no symmetries, br=bi=0
        call vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tw = we(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tw = wo(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(0, nlat, nlon, m, vb, iv, wvbin)
                call wbin(0, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * two2 + wb(i, np1, iw) * tve1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * two1 - wb(i, np1, iw) * tve2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1, iw) * ve(i, 2*mp1-1, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1, iw) * ve(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * twe2 + wb(i, np1, iw) * tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * twe1 - wb(i, np1, iw) * tvo2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-2, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (4) ! ityp=3: v even, w odd
        call vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tv = ve(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tw = wo(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(0, nlat, nlon, m, vb, iv, wvbin)
                call wbin(0, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * two2 + wb(i, np1, iw) * tve1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * two1 - wb(i, np1, iw) * tve2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1, iw) * ve(i, 2*mp1-1, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1, iw) * ve(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tve2 + wb(i, np1, iw) * two1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tve1 - wb(i, np1, iw) * two2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-2, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (5) ! ityp=4: v even, w odd, cr=ci=0
        call vbin(1, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tv = ve(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(1, nlat, nlon, m, vb, iv, wvbin)
                call wbin(1, nlat, nlon, m, wb, iw, wwbin)
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tve2 + wb(i, np1, iw) * two1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tve1 - wb(i, np1, iw) * two2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-2, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * ve(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (6) ! ityp=5: v even, w odd, br=bi=0
        call vbin(2, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tw = wo(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(2, nlat, nlon, m, vb, iv, wvbin)
                call wbin(2, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                two1 = wo(i, 2*mp1-1, k) * wts(i); two2 = wo(i, 2*mp1-2, k) * wts(i)
                                tve1 = ve(i, 2*mp1-1, k) * wts(i); tve2 = ve(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * two2 + wb(i, np1, iw) * tve1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * two1 - wb(i, np1, iw) * tve2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) + wb(i, np1, iw) * ve(i, 2*mp1-1, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - wb(i, np1, iw) * ve(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (7) ! ityp=6: v odd, w even
        call vbin(0, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tw = we(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tv = vo(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(0, nlat, nlon, m, vb, iv, wvbin)
                call wbin(0, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tvo2 + wb(i, np1, iw) * twe1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tvo1 - wb(i, np1, iw) * twe2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1, iw) * we(i, 2*mp1-1, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1, iw) * we(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * twe2 + wb(i, np1, iw) * tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * twe1 - wb(i, np1, iw) * tvo2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-2, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (8) ! ityp=7: v odd, w even, cr=ci=0
        call vbin(2, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    tv = vo(i, 1, k) * wts(i)
                    br(1, np1, k) = br(1, np1, k) + vb(i, np1, iv) * tv
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(2, nlat, nlon, m, vb, iv, wvbin)
                call wbin(2, nlat, nlon, m, wb, iw, wwbin)
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                br(mp1, np1, k) = br(mp1, np1, k) + vb(i, np1, iv) * tvo2 + wb(i, np1, iw) * twe1
                                bi(mp1, np1, k) = bi(mp1, np1, k) + vb(i, np1, iv) * tvo1 - wb(i, np1, iw) * twe2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                br(mp1, np1, k) = br(mp1, np1, k) + wb(i, np1, iw) * we(i, 2*mp1-1, k) * wts(i)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - wb(i, np1, iw) * we(i, 2*mp1-2, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (9) ! ityp=8: v odd, w even, br=bi=0
        call vbin(1, nlat, nlon, 0, vb, iv, wvbin)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    tw = we(i, 1, k) * wts(i)
                    cr(1, np1, k) = cr(1, np1, k) - vb(i, np1, iv) * tw
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1; mp2 = mp1 + 1
                call vbin(1, nlat, nlon, m, vb, iv, wvbin)
                call wbin(1, nlat, nlon, m, wb, iw, wwbin)
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                twe1 = we(i, 2*mp1-1, k) * wts(i); twe2 = we(i, 2*mp1-2, k) * wts(i)
                                tvo1 = vo(i, 2*mp1-1, k) * wts(i); tvo2 = vo(i, 2*mp1-2, k) * wts(i)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * twe2 + wb(i, np1, iw) * tvo1
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * twe1 - wb(i, np1, iw) * tvo2
                            end do
                            if (mlat /= 0) then
                                i = imid
                                cr(mp1, np1, k) = cr(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-2, k) * wts(i)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - vb(i, np1, iv) * we(i, 2*mp1-1, k) * wts(i)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    end select
end subroutine vhagc1


!> @brief Initialization routine for vector harmonic analysis on a Gaussian grid.
subroutine vhagci(nlat, nlon, wvhagc, lvhagc, dwork, ldwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, lvhagc, ldwork
    real, intent(out) :: wvhagc(*)
    double precision, intent(inout) :: dwork(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, lzz1, mmax, labc, lwk, jw1, jw2, jw3, iwrk, iw1, lwvbin, iw2, iw3, i

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if

    imid = (nlat + 1) / 2
    lzz1 = 2 * nlat * imid
    mmax = min(nlat, (nlon + 1) / 2)
    labc = 3 * (max(mmax - 2, 0) * (nlat + nlat - mmax - 1)) / 2
    if (lvhagc < 2 * (lzz1 + labc) + nlon + imid + 15) then; ierror = 3; return; end if

    if (ldwork < 2 * nlat * (nlat + 1) + 1) then; ierror = 4; return; end if
    ierror = 0

    ! --- Compute Gaussian points and weights (double precision) ---
    lwk = 2 * nlat * (nlat + 2)
    jw1 = 1
    jw2 = jw1 + nlat
    jw3 = jw2 + nlat
    call gaqd(nlat, dwork(jw1), dwork(jw2), dwork(jw3), lwk, ierror)

    ! --- Set single precision weights in the first part of wvhagc ---
    do i = 1, imid
        wvhagc(i) = dwork(nlat + i)
    end do

    ! --- Initialize vector basis functions and FFT workspace ---
    iwrk = (nlat + 1) / 2 + 1
    iw1 = imid + 1
    call vbgint(nlat, nlon, dwork, wvhagc(iw1), dwork(iwrk))

    lwvbin = lzz1 + labc
    iw2 = iw1 + lwvbin
    call wbgint(nlat, nlon, dwork, wvhagc(iw2), dwork(iwrk))

    iw3 = iw2 + lwvbin
    call hrffti(nlon, wvhagc(iw3))

end subroutine vhagci
