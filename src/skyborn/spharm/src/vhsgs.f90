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
! ... file vhsgs.f
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
!     - Added explicit initialization of workspace arrays (ve, vo, we, wo)
!       to zero at the start of the core computation, preventing errors
!       from uninitialized memory.
!     - Applied !$OMP SIMD directives to all computationally intensive inner
!       loops for explicit vectorization and significant speedup.
!     - Used temporary scalar variables to cache coefficients within loops,
!       reducing memory bandwidth and improving SIMD performance.
!
! =============================================================================


!> @brief Main vector spherical harmonic synthesis routine.
subroutine vhsgs(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhsgs, lvhsgs, work, lwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw, mdab, ndab, lvhsgs, lwork
    real, intent(out) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(in) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(in) :: wvhsgs(*)
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
    if (lvhsgs < lzimn + lzimn + nlon + 15) then; ierror = 9; return; end if

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
    call vhsgs1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhsgs(jw1), wvhsgs(jw2), wvhsgs(jw3))

end subroutine vhsgs


!> @brief Core vector harmonic synthesis computation routine with SIMD.
subroutine vhsgs1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
                  ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, vb, wb, wrfft)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, imid, idvw, jdvw, mdab, ndab, idv, idz
    real, intent(out) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(in) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(inout) :: ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), wo(idv, nlon, *)
    real, intent(inout) :: work(*)
    real, intent(in) :: vb(imid, *), wb(imid, *), wrfft(*)

    ! Local variables
    integer :: nlp1, mlat, mmax, imm1, ndo1, ndo2, itypp
    integer :: k, i, j, mp1, np1, m, mb, mp2, mn
    real :: br_val, bi_val, cr_val, ci_val

    ! --- Precompute constants ---
    nlp1 = nlat + 1
    mlat = mod(nlat, 2)
    mmax = min(nlat, (nlon + 1) / 2)
    imm1 = imid
    if (mlat /= 0) imm1 = imid - 1

    ! --- CRITICAL: Initialize workspace arrays to zero ---
    do k = 1, nt
        do j = 1, nlon
            do i = 1, idv
                ve(i, j, k) = 0.0
                we(i, j, k) = 0.0
            end do
        end do
    end do
    if (ityp <= 2) then
        do k = 1, nt
            do j = 1, nlon
                do i = 1, idv
                    vo(i, j, k) = 0.0
                    wo(i, j, k) = 0.0
                end do
            end do
        end do
    end if

    ndo1 = nlat
    ndo2 = nlat
    if (mlat /= 0) ndo1 = nlat - 1
    if (mlat == 0) ndo2 = nlat - 1

    itypp = ityp + 1

    ! --- Legendre transform via SELECT CASE on symmetry type ---
    select case (itypp)

    case (1) ! ityp=0: no symmetries
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                br_val = br(1, np1, k); cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    ve(i, 1, k) = ve(i, 1, k) + br_val * vb(i, np1)
                    we(i, 1, k) = we(i, 1, k) - cr_val * vb(i, np1)
                end do
            end do
            do np1 = 3, ndo1, 2
                br_val = br(1, np1, k); cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    vo(i, 1, k) = vo(i, 1, k) + br_val * vb(i, np1)
                    wo(i, 1, k) = wo(i, 1, k) - cr_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
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
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    ve(i, 1, k) = ve(i, 1, k) + br_val * vb(i, np1)
                end do
            end do
            do np1 = 3, ndo1, 2
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    vo(i, 1, k) = vo(i, 1, k) + br_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
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
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    we(i, 1, k) = we(i, 1, k) - cr_val * vb(i, np1)
                end do
            end do
            do np1 = 3, ndo1, 2
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    wo(i, 1, k) = wo(i, 1, k) - cr_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
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
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    ve(i, 1, k) = ve(i, 1, k) + br_val * vb(i, np1)
                end do
            end do
            do np1 = 3, ndo1, 2
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    wo(i, 1, k) = wo(i, 1, k) - cr_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
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
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    ve(i, 1, k) = ve(i, 1, k) + br_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb(i, mn)
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
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    wo(i, 1, k) = wo(i, 1, k) - cr_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb(i, mn)
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
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    we(i, 1, k) = we(i, 1, k) - cr_val * vb(i, np1)
                end do
            end do
            do np1 = 3, ndo1, 2
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    vo(i, 1, k) = vo(i, 1, k) + br_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
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
                br_val = br(1, np1, k)
                !$OMP SIMD
                do i = 1, imm1
                    vo(i, 1, k) = vo(i, 1, k) + br_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            br_val = br(mp1, np1, k); bi_val = bi(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) + br_val * vb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + bi_val * vb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb(i, mn)
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
                cr_val = cr(1, np1, k)
                !$OMP SIMD
                do i = 1, imid
                    we(i, 1, k) = we(i, 1, k) - cr_val * vb(i, np1)
                end do
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
                            mn = mb + np1
                            cr_val = cr(mp1, np1, k); ci_val = ci(mp1, np1, k)
                            !$OMP SIMD
                            do i = 1, imm1
                                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) - ci_val * wb(i, mn)
                                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + cr_val * wb(i, mn)
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end do
                            if (mlat /= 0) then
                                i = imid
                                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb(i, mn)
                                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb(i, mn)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    end select

    ! --- Perform Inverse Fourier Transform and combine components ---
    do k = 1, nt
        call hrfftb(idv, nlon, ve(1, 1, k), idv, wrfft, work)
        call hrfftb(idv, nlon, we(1, 1, k), idv, wrfft, work)
    end do

    if (ityp > 2) then
        do k = 1, nt
            do j = 1, nlon
                !$OMP SIMD
                do i = 1, idv
                    v(i, j, k) = 0.5 * ve(i, j, k)
                    w(i, j, k) = 0.5 * we(i, j, k)
                end do
            end do
        end do
    else
        do k = 1, nt
            do j = 1, nlon
                !$OMP SIMD
                do i = 1, imm1
                    v(i, j, k) = 0.5 * (ve(i, j, k) + vo(i, j, k))
                    w(i, j, k) = 0.5 * (we(i, j, k) + wo(i, j, k))
                    v(nlp1 - i, j, k) = 0.5 * (ve(i, j, k) - vo(i, j, k))
                    w(nlp1 - i, j, k) = 0.5 * (we(i, j, k) - wo(i, j, k))
                end do
            end do
        end do
    end if

    if (mlat /= 0) then
        do k = 1, nt
            do j = 1, nlon
                v(imid, j, k) = 0.5 * ve(imid, j, k)
                w(imid, j, k) = 0.5 * we(imid, j, k)
            end do
        end do
    end if

end subroutine vhsgs1


!> @brief Initialization routine for vector harmonic synthesis.
subroutine vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, lvhsgs, ldwork
    real, intent(out) :: wvhsgs(*)
    double precision, intent(inout) :: dwork(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, lmn, jw1, jw2, jw3, iw1, iw2, iw3, iw4

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if

    imid = (nlat + 1) / 2
    lmn = nlat * (nlat + 1) / 2
    if (lvhsgs < 2 * (imid * lmn) + nlon + 15) then; ierror = 3; return; end if

    if (ldwork < (nlat * 3 * (nlat + 3) + 2) / 2) then; ierror = 4; return; end if
    ierror = 0

    ! --- Workspace pointer setup ---
    jw1 = 1
    jw2 = jw1 + imid * lmn
    jw3 = jw2 + imid * lmn

    iw1 = 1
    iw2 = iw1 + nlat
    iw3 = iw2 + nlat
    iw4 = iw3 + 3 * imid * nlat

    ! --- Call initialization core routine and FFT setup ---
    call vhgsi1(nlat, imid, wvhsgs(jw1), wvhsgs(jw2), &
                dwork(iw1), dwork(iw2), dwork(iw3), dwork(iw4))

    call hrffti(nlon, wvhsgs(jw3))

end subroutine vhsgsi


!> @brief Core initialization routine for vector basis functions.
subroutine vhgsi1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, imid
    real, intent(out) :: vb(imid, *), wb(imid, *)
    double precision, intent(inout) :: dthet(*), dwts(*), dpbar(imid, nlat, 3), work(*)

    ! Local variables
    integer, parameter :: real64 = kind(1.0d0)
    integer :: lwk, ierror, i, n, nm, nz, np, m, ix, iy, ldim
    real(kind=real64) :: abel, bbel, cbel, ssqr2, dcf

    ! Statement function for index calculation (as in original F77 code)
    integer :: indx
    indx(m, n, ldim) = (m-1) * ldim - (m-1)*(m-2)/2 + (n-m+1)

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

        ! --- Compute vb (derivative of pbar) ---
        ix = n * (n + 1) / 2
        iy = ix + n
        do i = 1, imid
            vb(i, ix) = -dpbar(i, 2, np)
            vb(i, iy) = dpbar(i, n, np) / sqrt(real(2*(n+1), kind=real64))
        end do

        if (n >= 2) then
            dcf = sqrt(real(4*n*(n+1), kind=real64))
            do m = 1, n - 1
                ix = (m-1)*nlat - (m-1)*(m-2)/2 + (n-m+1)
                abel = sqrt(real((n+m)*(n-m+1), kind=real64)) / dcf
                bbel = sqrt(real((n-m)*(n+m+1), kind=real64)) / dcf
                do i = 1, imid
                    vb(i, ix) = abel * dpbar(i, m, np) - bbel * dpbar(i, m+2, np)
                end do
            end do
        end if

        ! --- Compute wb (m*pbar/sin(theta)) ---
        ix = n * (n + 1) / 2
        do i = 1, imid
            wb(i, ix) = 0.0_real64
        end do

        dcf = sqrt(real(2*n+1, kind=real64) / real(4*n*(n+1)*(2*n-1), kind=real64))
        do m = 1, n
            ix = (m-1)*nlat - (m-1)*(m-2)/2 + (n-m+1)
            abel = dcf * sqrt(real((n+m)*(n+m-1), kind=real64))
            bbel = dcf * sqrt(real((n-m)*(n-m-1), kind=real64))
            if (m < n - 1) then
                do i = 1, imid
                    wb(i, ix) = abel * dpbar(i, m, nz) + bbel * dpbar(i, m+2, nz)
                end do
            else
                do i = 1, imid
                    wb(i, ix) = abel * dpbar(i, m, nz)
                end do
            end if
        end do
    end do

end subroutine vhgsi1
