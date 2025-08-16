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
! ... file vhses.f
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


!
!     this file contains code and documentation for subroutines
!     vhses and vhsesi
!
! ... files which must be loaded with vhses.f
!
!     sphcom.f, hrfft.f
!
!
!     subroutine vhses(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci,
!    +                 mdab,ndab,wvhses,lvhses,work,lwork,ierror)
!
!     subroutine vhses performs the vector spherical harmonic synthesis
!     of the arrays br, bi, cr, and ci and stores the result in the
!     arrays v and w. v(i,j) and w(i,j) are the colatitudinal
!     (measured from the north pole) and east longitudinal components
!     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1)
!     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral
!     representation of (v,w) is given below at output parameters v,w.
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     ityp   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i,j),w(i,j) for i=1,...,nlat and
!                 j=1,...,nlon.
!
!            = 1  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i,j),w(i,j) for i=1,...,nlat and
!                 j=1,...,nlon. the curl of (v,w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
!                 the coefficients cr and ci are zero.
!
!            = 2  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 arrays v(i,j),w(i,j) for i=1,...,nlat and
!                 j=1,...,nlon. the divergence of (v,w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
!                 the coefficients br and bi are zero.
!
!            = 3  v is symmetric and w is antisymmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!            = 4  v is symmetric and w is antisymmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!                 the curl of (v,w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
!                 the coefficients cr and ci are zero.
!
!            = 5  v is symmetric and w is antisymmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!                 the divergence of (v,w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
!                 the coefficients br and bi are zero.
!
!            = 6  v is antisymmetric and w is symmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!            = 7  v is antisymmetric and w is symmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!                 the curl of (v,w) is zero. that is,
!                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0.
!                 the coefficients cr and ci are zero.
!
!            = 8  v is antisymmetric and w is symmetric about the
!                 equator. the synthesis is performed on the northern
!                 hemisphere only.  i.e., if nlat is odd the synthesis
!                 is performed on the arrays v(i,j),w(i,j) for
!                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is
!                 even the synthesis is performed on the the arrays
!                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!                 the divergence of (v,w) is zero. i.e.,
!                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0.
!                 the coefficients br and bi are zero.
!
!
!     nt     the number of syntheses.  in the program that calls vhses,
!            the arrays v,w,br,bi,cr, and ci can be three dimensional
!            in which case multiple syntheses will be performed.
!            the third index is the synthesis index which assumes the
!            values k=1,...,nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that all the arrays are two
!            dimensional.
!
!     idvw   the first dimension of the arrays v,w as it appears in
!            the program that calls vhaes. if ityp .le. 2 then idvw
!            must be at least nlat.  if ityp .gt. 2 and nlat is
!            even then idvw must be at least nlat/2. if ityp .gt. 2
!            and nlat is odd then idvw must be at least (nlat+1)/2.
!
!     jdvw   the second dimension of the arrays v,w as it appears in
!            the program that calls vhses. jdvw must be at least nlon.
!
!     br,bi  two or three dimensional arrays (see input parameter nt)
!     cr,ci  that contain the vector spherical harmonic coefficients
!            in the spectral representation of v(i,j) and w(i,j) given
!            below at the discription of output parameters v and w.
!
!     mdab   the first dimension of the arrays br,bi,cr, and ci as it
!            appears in the program that calls vhses. mdab must be at
!            least min0(nlat,nlon/2) if nlon is even or at least
!            min0(nlat,(nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays br,bi,cr, and ci as it
!            appears in the program that calls vhses. ndab must be at
!            least nlat.
!
!     wvhses an array which must be initialized by subroutine vhsesi.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlon and nlat remain unchanged.  wvhses must
!            not be altered between calls of vhses.
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min0(nlat,nlon/2) if nlon is even or
!               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                 l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if ityp .le. 2 then lwork must be at least
!
!                       (2*nt+1)*nlat*nlon
!
!            if ityp .gt. 2 then lwork must be at least
!
!                        (2*nt+1)*l2*nlon
!
!     **************************************************************
!
!     output parameters
!
!     v,w    two or three dimensional arrays (see input parameter nt)
!            in which the synthesis is stored. v is the colatitudinal
!            component and w is the east longitudinal component.
!            v(i,j),w(i,j) contain the components at colatitude
!            theta(i) = (i-1)*pi/(nlat-1) and longitude phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            the input parameter ityp. v and w are computed from the
!            formulas given below
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of ityp
!            = 4  error in the specification of nt
!            = 5  error in the specification of idvw
!            = 6  error in the specification of jdvw
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lvhses
!            = 10 error in the specification of lwork
!
! ************************************************************
!
!> @brief Main vector spherical harmonic synthesis routine on an equally spaced grid.
subroutine vhses(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhses, lvhses, work, lwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw, mdab, ndab, lvhses, lwork
    real, intent(out) :: v(idvw, jdvw, *), w(idvw, jdvw, *)
    real, intent(in) :: br(mdab, ndab, *), bi(mdab, ndab, *), cr(mdab, ndab, *), ci(mdab, ndab, *)
    real, intent(in) :: wvhses(*)
    real, intent(inout) :: work(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, mmax, idz, lzimn, idv, lnl, ist
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
    if (lvhses < lzimn + lzimn + nlon + 15) then; ierror = 9; return; end if

    idv = nlat
    if (ityp > 2) idv = imid
    lnl = nt * idv * nlon
    if (lwork < lnl + lnl + idv * nlon) then; ierror = 10; return; end if

    ierror = 0

    ! --- Workspace pointer setup ---
    jw1 = 1
    jw2 = lzimn + 1
    jw3 = jw2 + lzimn

    ist = 0
    if (ityp <= 2) ist = imid
    iw1 = ist + 1
    iw2 = lnl + 1
    iw3 = iw2 + ist
    iw4 = iw2 + lnl

    ! --- Call the core computational routine ---
    call vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
                br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
                work(iw4), idz, wvhses(jw1), wvhses(jw2), wvhses(jw3))

end subroutine vhses


!> @brief Core vector harmonic synthesis computation routine with SIMD.
subroutine vhses1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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
                mb = m * (nlat - 1) - (m * (m - 1)) / 2
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

end subroutine vhses1


!
!     subroutine vhsesi(nlat,nlon,wvhses,lvhses,work,lwork,dwork,
!    +                  ldwork,ierror)
!
!     subroutine vhsesi initializes the array wvhses which can then be
!     used repeatedly by subroutine vhses until nlat or nlon is changed.
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than zero. the axisymmetric case corresponds to nlon=1.
!            the efficiency of the computation is improved when nlon
!            is a product of small prime numbers.
!
!     lvhses the dimension of the array wvhses as it appears in the
!            program that calls vhses. define
!
!               l1 = min0(nlat,nlon/2) if nlon is even or
!               l1 = min0(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lvhses must be at least
!
!                  l1*l2*(nlat+nlat-l1+1)+nlon+15
!
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls vhses. lwork must be at least
!
!              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
!
!     dwork  an unsaved double precision work space
!
!     ldwork the length of the array dwork as it appears in the
!            program that calls vhsesi.  ldwork must be at least
!            2*(nlat+1)
!
!
!     **************************************************************
!
!     output parameters
!
!     wvhses an array which is initialized for use by subroutine vhses.
!            once initialized, wvhses can be used repeatedly by vhses
!            as long as nlat or nlon remain unchanged.  wvhses must not
!            be altered between calls of vhses.
!
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lvhses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
! *****************************************
!> @brief Initialization routine for vector harmonic synthesis.
subroutine vhsesi(nlat, nlon, wvhses, lvhses, work, lwork, dwork, ldwork, ierror)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, lvhses, lwork, ldwork
    real, intent(out) :: wvhses(*)
    real, intent(inout) :: work(*)
    double precision, intent(inout) :: dwork(*)
    integer, intent(out) :: ierror

    ! Local variables
    integer :: imid, mmax, lzimn, labc, idz
    integer :: jw1, jw2, jw3, iw1

    ! --- Input validation ---
    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if

    imid = (nlat + 1) / 2
    mmax = min(nlat, (nlon + 1) / 2)
    lzimn = (imid * mmax * (nlat + nlat - mmax + 1)) / 2
    if (lvhses < lzimn + lzimn + nlon + 15) then; ierror = 3; return; end if

    labc = 3 * (max(mmax - 2, 0) * (nlat + nlat - mmax - 1)) / 2
    if (lwork < 5 * nlat * imid + labc) then; ierror = 4; return; end if

    if (ldwork < 2 * (nlat + 1)) then; ierror = 5; return; end if
    ierror = 0

    ! --- Workspace pointer setup and initialization calls ---
    idz = (mmax * (nlat + nlat - mmax + 1)) / 2
    jw1 = 1
    jw2 = lzimn + 1
    jw3 = jw2 + lzimn

    iw1 = 3 * nlat * imid + 1

    call ves1(nlat, nlon, imid, wvhses, wvhses(jw2), idz, work, work(iw1), dwork)

    call hrffti(nlon, wvhses(jw3))

end subroutine vhsesi


!> @brief Core initialization routine for vector basis functions on an equally spaced grid.
subroutine ves1(nlat, nlon, imid, vb, wb, idz, vin, wzvin, dwork)
    implicit none

    ! Argument definitions
    integer, intent(in) :: nlat, nlon, imid, idz
    real, intent(out) :: vb(imid, *), wb(imid, *)
    real, intent(inout) :: vin(imid, nlat, 3), wzvin(*)
    double precision, intent(inout) :: dwork(*)

    ! Local variables
    integer :: mmax, mp1, m, np1, mn, i, i3

    mmax = min(nlat, (nlon + 1) / 2)

    ! --- Initialize V-component basis functions ---
    call vbinit(nlat, nlon, wzvin, dwork)
    do mp1 = 1, mmax
        m = mp1 - 1
        call vbin(0, nlat, nlon, m, vin, i3, wzvin)
        do np1 = mp1, nlat
            mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
            do i = 1, imid
                vb(i, mn) = vin(i, np1, i3)
            end do
        end do
    end do

    ! --- Initialize W-component basis functions ---
    call wbinit(nlat, nlon, wzvin, dwork)
    do mp1 = 1, mmax
        m = mp1 - 1
        call wbin(0, nlat, nlon, m, vin, i3, wzvin)
        do np1 = mp1, nlat
            mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
            do i = 1, imid
                wb(i, mn) = vin(i, np1, i3)
            end do
        end do
    end do

end subroutine ves1
