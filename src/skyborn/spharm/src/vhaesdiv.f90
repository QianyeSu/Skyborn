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
! ... file vhaesdiv.f90
!
!     Internal divergence-only blocked analysis wrapper for the regular/stored
!     backend. This keeps the mathematical kernel in `vhaes1` unchanged while
!     reducing the live worker size for large `nt` workloads.
!
! =============================================================================

subroutine vhaesdiv(nlat, nlon, nt, v, w, idvw, jdvw, br, bi, mdab, ndab, &
                    wvhaes, lvhaes, nbatch, work, lwork, ierror)
    implicit none

    integer, intent(in) :: nlat, nlon, nt, idvw, jdvw, mdab, ndab, lvhaes
    integer, intent(in) :: nbatch, lwork
    real, intent(in) :: v(idvw, jdvw, nt), w(idvw, jdvw, nt)
    real, intent(out) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
    real, intent(in) :: wvhaes(lvhaes)
    real, intent(inout) :: work(lwork)
    integer, intent(out) :: ierror

    integer :: imid, mmax, idz, lzimn, block_nt, active_lnl
    integer :: jw1, jw2, iw1, iw2, iw3, iw4, kstart, kstop, nblk
    real :: dummy_cr(1, 1, max(1, min(nt, max(1, nbatch))))
    real :: dummy_ci(1, 1, max(1, min(nt, max(1, nbatch))))

    if (nlat < 3) then; ierror = 1; return; end if
    if (nlon < 1) then; ierror = 2; return; end if
    if (nt < 0) then; ierror = 4; return; end if

    imid = (nlat + 1) / 2
    if (idvw < nlat) then; ierror = 5; return; end if
    if (jdvw < nlon) then; ierror = 6; return; end if

    mmax = min(nlat, (nlon + 1) / 2)
    if (mdab < mmax) then; ierror = 7; return; end if
    if (ndab < nlat) then; ierror = 8; return; end if

    idz = (mmax * (nlat + nlat - mmax + 1)) / 2
    lzimn = idz * imid
    if (lvhaes < lzimn + lzimn + nlon + 15) then; ierror = 9; return; end if

    block_nt = max(1, min(nt, max(1, nbatch)))
    active_lnl = block_nt * nlat * nlon
    if (lwork < 2 * active_lnl + nlat * nlon) then
        ierror = 10
        return
    end if

    ierror = 0

    iw1 = imid + 1
    iw2 = active_lnl + 1
    iw3 = iw2 + imid
    iw4 = iw2 + active_lnl
    jw1 = lzimn + 1
    jw2 = jw1 + lzimn

    if (nt == 0) return

    do kstart = 1, nt, block_nt
        kstop = min(kstart + block_nt - 1, nt)
        nblk = kstop - kstart + 1
        call vhaes1( &
            nlat, nlon, 1, nblk, imid, idvw, jdvw, &
            v(:, :, kstart:kstop), w(:, :, kstart:kstop), &
            mdab, ndab, br(:, :, kstart:kstop), bi(:, :, kstart:kstop), &
            dummy_cr(:, :, 1:nblk), dummy_ci(:, :, 1:nblk), &
            nlat, work, work(iw1), work(iw2), work(iw3), work(iw4), idz, &
            wvhaes, wvhaes(jw1), wvhaes(jw2))
    end do
end subroutine vhaesdiv
