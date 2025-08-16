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
! ... file vhaes.f
!
!     Modernized to Fortran 2008 standard.
!     - Converted to free-form source.
!     - All subroutines are external (global), without a containing module.
!     - Added 'implicit none' and explicit variable declarations.
!     - Replaced archaic control flow (GOTO, computed GOTO) with modern
!       constructs (IF/SELECT CASE, DO loops).
!     - Used intent attributes for arguments.
!     - Used assumed-size arrays dimension(*) to maintain F77 compatibility
!       and avoid the need for explicit interfaces.
!
!     CRITICAL BUG FIX:
!     - Replicated the exact conditional initialization logic for the output
!       spectral coefficient arrays (br, bi, cr, ci) from the original
!       F77 code inside vhaes1. The lack of this precise initialization
!       in previous modernization attempts was the source of major numerical
!       errors (170%+ RMSE) and NaN values.
!
! =============================================================================
! This file contains optimized versions of:
!   - vhaes: Main vector spherical harmonic analysis routine
!   - vhaes1: Core computation routine
!   - vhaesi: Initialization routine
!   - vea1: Core initialization routine

!> @brief Optimized main vector spherical harmonic analysis routine
!> Modernized from FORTRAN 77 to Fortran 2008+ standards
subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)
    implicit none

    ! Input parameters
    integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw
    integer, intent(in) :: mdab, ndab, lvhaes, lwork
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(in) :: wvhaes(*)

    ! Output parameters
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *)
    real, intent(out) :: cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(inout) :: work(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, mmax, idz, lzimn, idv, lnl, ist
    integer :: iw1, iw2, iw3, iw4, jw1, jw2

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
    if (lvhaes < lzimn + lzimn + nlon + 15) then; ierror = 9; return; end if

    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt * idv * nlon
    if (lwork < lnl + lnl + idv * nlon) then; ierror = 10; return; end if

    ierror = 0

    ! --- Workspace pointer setup ---
    ist = 0
    if (ityp <= 2) ist = imid
    iw1 = ist + 1
    iw2 = lnl + 1
    iw3 = iw2 + ist
    iw4 = iw2 + lnl
    jw1 = lzimn + 1
    jw2 = jw1 + lzimn

    ! --- Call the core computational routine ---
    call vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhaes, wvhaes(jw1), wvhaes(jw2))

end subroutine vhaes


!> @brief Core vector harmonic analysis computation routine
subroutine vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
                  ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, zv, zw, wrfft)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, imid, idvw, jdvw, mdab, ndab, idv, idz
    real, intent(in) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(out) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(inout) :: ve(idv, nlon, *), vo(idv, nlon, *), we(idv, nlon, *), wo(idv, nlon, *)
    real, intent(inout) :: work(*)
    real, intent(in) :: zv(idz, *), zw(idz, *), wrfft(*)

    ! Local variables
    integer :: nlp1, mlat, mmax, imm1, ndo1, ndo2, itypp
    integer :: k, i, j, mp1, np1, m, mb, mp2
    real :: tsn, fsn

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
                !$OMP SIMD
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
                !$OMP SIMD
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

    ! --- Perform Forward Fourier Transform along longitude lines ---
    do k = 1, nt
        call hrfftf(idv, nlon, ve(1, 1, k), idv, wrfft, work)
        call hrfftf(idv, nlon, we(1, 1, k), idv, wrfft, work)
    end do

    ndo1 = nlat
    ndo2 = nlat
    if (mlat /= 0) ndo1 = nlat - 1
    if (mlat == 0) ndo2 = nlat - 1

    ! =======================================================================
    ! CRITICAL FIX: Conditionally initialize spectral coefficient arrays to zero.
    ! This logic exactly matches the original F77 code and prevents the use
    ! of uninitialized memory in the accumulation loops, which was the cause
    ! of the 170% error.
    ! =======================================================================
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
                !$OMP SIMD
                do i = 1, imid
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (2) ! ityp=1: no symmetries, cr=ci=0 (curl-free)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    case (3) ! ityp=2: no symmetries, br=bi=0 (divergence-free)
        ! m=0
        do k = 1, nt
            do np1 = 2, ndo2, 2
                !$OMP SIMD
                do i = 1, imid
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
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
                !$OMP SIMD
                do i = 1, imid
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
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
                !$OMP SIMD
                do i = 1, imid
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
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
                !$OMP SIMD
                do i = 1, imm1
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
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
                !$OMP SIMD
                do i = 1, imid
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
                end do
            end do
            do np1 = 3, ndo1, 2
                !$OMP SIMD
                do i = 1, imm1
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                            end if
                        end do
                    end do
                end if
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
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
                !$OMP SIMD
                do i = 1, imm1
                    br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                if (mp1 <= ndo1) then
                    do k = 1, nt
                        do np1 = mp1, ndo1, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                                bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
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
                !$OMP SIMD
                do i = 1, imid
                    cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
                end do
            end do
        end do
        ! m>0
        if (mmax >= 2) then
            do mp1 = 2, mmax
                m = mp1 - 1
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
                mp2 = mp1 + 1
                if (mp2 <= ndo2) then
                    do k = 1, nt
                        do np1 = mp2, ndo2, 2
                            !$OMP SIMD
                            do i = 1, imm1
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                            end do
                            if (mlat /= 0) then
                                cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                                ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                            end if
                        end do
                    end do
                end if
            end do
        end if

    end select
end subroutine vhaes1

!> @brief Initialization routine for vector harmonic analysis
!> Initializes the workspace array wvhaes for repeated use by vhaes
!> @param[in] nlat Number of colatitudes on full sphere
!> @param[in] nlon Number of longitude points
!> @param[out] wvhaes Workspace array to be initialized
!> @param[in] lvhaes Dimension of wvhaes array
!> @param[inout] work Temporary workspace array
!> @param[in] lwork Dimension of work array
!> @param[inout] dwork Double precision workspace
!> @param[in] ldwork Dimension of dwork array
!> @param[out] ierror Error flag (0=success, >0=error code)
subroutine vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, ldwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, lvhaes, lwork, ldwork
    real, intent(out) :: wvhaes(*)
    real, intent(inout) :: work(*)
    double precision, intent(inout) :: dwork(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: mmax, imid, lzimn, labc, iw1, idz

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if

    mmax = min(nlat, (nlon + 1) / 2)
    imid = (nlat + 1) / 2
    lzimn = (imid * mmax * (nlat + nlat - mmax + 1)) / 2
    if (lvhaes < lzimn + lzimn + nlon + 15) then; ierror = 3; return; end if

    labc = 3 * (max(mmax - 2, 0) * (nlat + nlat - mmax - 1)) / 2
    if (lwork < 5 * nlat * imid + labc) then; ierror = 4; return; end if

    if (ldwork < 2 * (nlat + 1)) then; ierror = 5; return; end if

    ierror = 0

    ! --- Workspace pointer setup and initialization calls ---
    iw1 = 3 * nlat * imid + 1
    idz = (mmax * (nlat + nlat - mmax + 1)) / 2

    call vea1(nlat, nlon, imid, wvhaes, wvhaes(lzimn + 1), idz, &
              work, work(iw1), dwork)

    call hrffti(nlon, wvhaes(2 * lzimn + 1))

end subroutine vhaesi


!> @brief Core initialization routine for vector analysis basis functions
subroutine vea1(nlat, nlon, imid, zv, zw, idz, zin, wzvin, dwork)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, imid, idz
    real, intent(out) :: zv(idz, *), zw(idz, *)
    real, intent(inout) :: zin(imid, nlat, 3), wzvin(*)
    double precision, intent(inout) :: dwork(*)

    ! Local variables
    integer :: mmax, mp1, m, np1, mn, i, i3

    mmax = min(nlat, (nlon + 1) / 2)

    ! --- Initialize V-component analysis functions (zv) ---
    call zvinit(nlat, nlon, wzvin, dwork)
    do mp1 = 1, mmax
        m = mp1 - 1
        call zvin(0, nlat, nlon, m, zin, i3, wzvin)
        do np1 = mp1, nlat
            mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
            do i = 1, imid
                zv(mn, i) = zin(i, np1, i3)
            end do
        end do
    end do

    ! --- Initialize W-component analysis functions (zw) ---
    call zwinit(nlat, nlon, wzvin, dwork)
    do mp1 = 1, mmax
        m = mp1 - 1
        call zwin(0, nlat, nlon, m, zin, i3, wzvin)
        do np1 = mp1, nlat
            mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
            do i = 1, imid
                zw(mn, i) = zin(i, np1, i3)
            end do
        end do
    end do

end subroutine vea1
