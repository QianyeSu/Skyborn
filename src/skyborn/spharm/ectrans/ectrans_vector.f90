module vector_backend
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  implicit none
  private

  public :: vrtdiv_analysis_stub
  public :: uv_synthesis_stub
  public :: gradient_synthesis_stub

  interface
    subroutine vhgsi1(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
      integer, intent(in) :: nlat, imid
      real, intent(out) :: vb(imid, *)
      real, intent(out) :: wb(imid, *)
      double precision, intent(inout) :: dthet(*), dwts(*), dpbar(imid, nlat, 3), work(*)
    end subroutine vhgsi1

    subroutine hrffti(n, wsave)
      integer, intent(in) :: n
      real, intent(out) :: wsave(2 * n + 15)
    end subroutine hrffti

    subroutine hrfftb(m, n, r, mdimr, wsave, work)
      integer, intent(in) :: m, n, mdimr
      real, intent(inout) :: r(mdimr, n)
      real, intent(in) :: wsave(2 * n + 15)
      real, intent(inout) :: work(m, n)
    end subroutine hrfftb
  end interface

contains

  subroutine vrtdiv_analysis_stub( &
    ndgl, nloen, ngptot, ntrunc, nt, ugrid, vgrid, vrtspec_r, vrtspec_i, divspec_r, divspec_i, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: ugrid(ngptot * nt)
    real(c_double), intent(in) :: vgrid(ngptot * nt)
    real(c_double), intent(out) :: vrtspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: vrtspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    integer(c_int), intent(out) :: ierror

    vrtspec_r = 0.0_c_double
    vrtspec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double
    ierror = -1_c_int
  end subroutine vrtdiv_analysis_stub

  subroutine uv_synthesis_stub( &
    ndgl, nloen, ngptot, rsphere, ntrunc, nt, vrtspec_r, vrtspec_i, divspec_r, divspec_i, ugrid, vgrid, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    real(c_double), value, intent(in) :: rsphere
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: vrtspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: vrtspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: ugrid(ngptot * nt)
    real(c_double), intent(out) :: vgrid(ngptot * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: imid, lmn, ldwork, mmax, ndo1, ndo2, mlat, imm1
    integer(c_int) :: idz, ilat, it, j, nlon_lat, offset, nlp1, source_ilat, current_nlon
    real, allocatable :: vb(:,:), wb(:,:)
    real, allocatable :: br(:,:,:), bi(:,:,:), cr(:,:,:), ci(:,:,:)
    real, allocatable :: ve(:,:,:), vo(:,:,:), we(:,:,:), wo(:,:,:)
    real, allocatable :: vpack(:,:), wpack(:,:), fft_work(:,:), wsave(:)
    double precision, allocatable :: dtheta(:), dwts(:), dpbar(:,:,:), dwork(:)

    if (ndgl < 3_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int .or. ntrunc < 0_c_int) then
      ugrid = 0.0_c_double
      vgrid = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    if (ntrunc > ndgl - 1_c_int) then
      ugrid = 0.0_c_double
      vgrid = 0.0_c_double
      ierror = 2_c_int
      return
    end if

    imid = (ndgl + 1_c_int) / 2_c_int
    lmn = ndgl * (ndgl + 1_c_int) / 2_c_int
    ldwork = (ndgl * 3_c_int * (ndgl + 3_c_int) + 2_c_int) / 2_c_int
    idz = imid * lmn
    mmax = min(ndgl, (maxval(nloen) + 1_c_int) / 2_c_int)
    mlat = mod(ndgl, 2_c_int)
    imm1 = imid
    if (mlat /= 0_c_int) imm1 = imid - 1_c_int
    nlp1 = ndgl + 1_c_int
    ndo1 = ndgl
    ndo2 = ndgl
    if (mlat /= 0_c_int) ndo1 = ndgl - 1_c_int
    if (mlat == 0_c_int) ndo2 = ndgl - 1_c_int

    allocate(vb(imid, idz))
    allocate(wb(imid, idz))
    allocate(dtheta(ndgl), dwts(ndgl), dpbar(imid, ndgl, 3), dwork(ldwork))
    allocate(br(ndgl, ndgl, nt), bi(ndgl, ndgl, nt), cr(ndgl, ndgl, nt), ci(ndgl, ndgl, nt))
    allocate(ve(imid, 2 * mmax, nt), vo(imid, 2 * mmax, nt), we(imid, 2 * mmax, nt), wo(imid, 2 * mmax, nt))

    call vhgsi1(ndgl, imid, vb, wb, dtheta, dwts, dpbar, dwork)

    br = 0.0
    bi = 0.0
    cr = 0.0
    ci = 0.0
    call unpack_vrtdiv_to_vector_harmonics(ndgl, ntrunc, nt, rsphere, vrtspec_r, vrtspec_i, divspec_r, divspec_i, br, bi, cr, ci)

    ve = 0.0
    vo = 0.0
    we = 0.0
    wo = 0.0
    call accumulate_gaussian_vector_modes( &
      ndgl, nt, imid, imm1, mlat, mmax, ndo1, ndo2, vb, wb, br, bi, cr, ci, ve, vo, we, wo)

    ugrid = 0.0_c_double
    vgrid = 0.0_c_double
    offset = 0_c_int
    current_nlon = -1_c_int

    do ilat = 1_c_int, ndgl
      nlon_lat = nloen(ilat)
      if (nlon_lat /= current_nlon) then
        if (allocated(vpack)) then
          deallocate(vpack, wpack, fft_work, wsave)
        end if
        allocate(vpack(2, nlon_lat), wpack(2, nlon_lat), fft_work(2, nlon_lat), wsave(2 * nlon_lat + 15))
        call hrffti(nlon_lat, wsave)
        current_nlon = nlon_lat
      end if

      source_ilat = ilat
      if (ilat > imid) source_ilat = nlp1 - ilat

      do it = 1_c_int, nt
        vpack = 0.0
        wpack = 0.0
        vpack(1, 1:min(size(ve, 2), nlon_lat)) = ve(source_ilat, 1:min(size(ve, 2), nlon_lat), it)
        wpack(1, 1:min(size(we, 2), nlon_lat)) = we(source_ilat, 1:min(size(we, 2), nlon_lat), it)

        if (ilat /= imid .or. mlat == 0_c_int) then
          vpack(2, 1:min(size(vo, 2), nlon_lat)) = vo(source_ilat, 1:min(size(vo, 2), nlon_lat), it)
          wpack(2, 1:min(size(wo, 2), nlon_lat)) = wo(source_ilat, 1:min(size(wo, 2), nlon_lat), it)
          call hrfftb(2, nlon_lat, vpack, 2, wsave, fft_work)
          call hrfftb(2, nlon_lat, wpack, 2, wsave, fft_work)
        else
          call hrfftb(1, nlon_lat, vpack(1:1, 1:nlon_lat), 1, wsave, fft_work(1:1, 1:nlon_lat))
          call hrfftb(1, nlon_lat, wpack(1:1, 1:nlon_lat), 1, wsave, fft_work(1:1, 1:nlon_lat))
        end if

        do j = 1_c_int, nlon_lat
          if (ilat < imid .or. (mlat == 0_c_int .and. ilat <= imid)) then
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(wpack(1, j) + wpack(2, j), c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = -0.5_c_double * real(vpack(1, j) + vpack(2, j), c_double)
          else if (ilat == imid .and. mlat /= 0_c_int) then
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(wpack(1, j), c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = -0.5_c_double * real(vpack(1, j), c_double)
          else
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(wpack(1, j) - wpack(2, j), c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = -0.5_c_double * real(vpack(1, j) - vpack(2, j), c_double)
          end if
        end do
      end do
      offset = offset + nlon_lat
    end do

    if (allocated(vpack)) deallocate(vpack, wpack, fft_work, wsave)
    deallocate(vb, wb, dtheta, dwts, dpbar, dwork, br, bi, cr, ci, ve, vo, we, wo)
    ierror = 0_c_int
  end subroutine uv_synthesis_stub

  subroutine gradient_synthesis_stub( &
    ndgl, nloen, ngptot, rsphere, ntrunc, nt, chispec_r, chispec_i, ugrad, vgrad, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    real(c_double), value, intent(in) :: rsphere
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: chispec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: chispec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: ugrad(ngptot * nt)
    real(c_double), intent(out) :: vgrad(ngptot * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ncoeff, nm, it, ierr_local
    real(c_double) :: rsphere_inv_sq, nreal, lap_factor
    real(c_double), allocatable :: zerospec_r(:), zerospec_i(:), divspec_r(:), divspec_i(:)

    if (ndgl < 3_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int .or. ntrunc < 0_c_int) then
      ugrad = 0.0_c_double
      vgrad = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    allocate(zerospec_r(ncoeff * nt), zerospec_i(ncoeff * nt), divspec_r(ncoeff * nt), divspec_i(ncoeff * nt))
    zerospec_r = 0.0_c_double
    zerospec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double

    rsphere_inv_sq = 1.0_c_double / (rsphere * rsphere)
    call apply_lap_to_scalar_spec(ntrunc, nt, chispec_r, chispec_i, rsphere_inv_sq, divspec_r, divspec_i)

    call uv_synthesis_stub( &
      ndgl, nloen, ngptot, rsphere, ntrunc, nt, &
      zerospec_r, zerospec_i, divspec_r, divspec_i, &
      ugrad, vgrad, ierr_local)

    deallocate(zerospec_r, zerospec_i, divspec_r, divspec_i)
    ierror = ierr_local
  end subroutine gradient_synthesis_stub

  subroutine apply_lap_to_scalar_spec( &
    ntrunc, nt, spec_r, spec_i, rsphere_inv_sq, lap_r, lap_i)
    integer(c_int), intent(in) :: ntrunc, nt
    real(c_double), intent(in) :: spec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: spec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: rsphere_inv_sq
    real(c_double), intent(out) :: lap_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: lap_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)

    integer(c_int) :: m, n, nm, nmstrt, it
    real(c_double) :: nreal, lap_factor

    lap_r = 0.0_c_double
    lap_i = 0.0_c_double

    do it = 1_c_int, nt
      nmstrt = 0_c_int
      do m = 1_c_int, ntrunc + 1_c_int
        do n = m, ntrunc + 1_c_int
          nm = nmstrt + n - m + 1_c_int
          nreal = real(n, c_double)
          lap_factor = -nreal * real(n - 1_c_int, c_double) * rsphere_inv_sq
          lap_r((nm - 1_c_int) * nt + it) = lap_factor * spec_r((nm - 1_c_int) * nt + it)
          lap_i((nm - 1_c_int) * nt + it) = lap_factor * spec_i((nm - 1_c_int) * nt + it)
        end do
        nmstrt = nmstrt + ntrunc - m + 2_c_int
      end do
    end do
  end subroutine apply_lap_to_scalar_spec

  subroutine unpack_vrtdiv_to_vector_harmonics( &
    ndgl, ntrunc, nt, rsphere, vrtspec_r, vrtspec_i, divspec_r, divspec_i, br, bi, cr, ci)
    integer(c_int), intent(in) :: ndgl, ntrunc, nt
    real(c_double), intent(in) :: rsphere
    real(c_double), intent(in) :: vrtspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: vrtspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real, intent(out) :: br(ndgl, ndgl, nt), bi(ndgl, ndgl, nt), cr(ndgl, ndgl, nt), ci(ndgl, ndgl, nt)

    integer(c_int) :: i, m, n, nm, nmstrt
    real(c_double) :: scale, nreal, factor

    br = 0.0
    bi = 0.0
    cr = 0.0
    ci = 0.0
    scale = 0.5_c_double

    do i = 1_c_int, nt
      nmstrt = 0_c_int
      do m = 1_c_int, ntrunc + 1_c_int
        do n = m, ntrunc + 1_c_int
          nm = nmstrt + n - m + 1_c_int
          nreal = real(n - 1_c_int, c_double)
          factor = 1.0_c_double
          if (nreal > 0.0_c_double) then
            factor = rsphere / sqrt(nreal * (nreal + 1.0_c_double))
          end if
          br(m, n, i) = real(-divspec_r((nm - 1_c_int) * nt + i) * factor / scale, kind=kind(br))
          bi(m, n, i) = real(-divspec_i((nm - 1_c_int) * nt + i) * factor / scale, kind=kind(bi))
          cr(m, n, i) = real(vrtspec_r((nm - 1_c_int) * nt + i) * factor / scale, kind=kind(cr))
          ci(m, n, i) = real(vrtspec_i((nm - 1_c_int) * nt + i) * factor / scale, kind=kind(ci))
        end do
        nmstrt = nmstrt + ntrunc - m + 2_c_int
      end do
    end do
  end subroutine unpack_vrtdiv_to_vector_harmonics

  subroutine accumulate_gaussian_vector_modes( &
    ndgl, nt, imid, imm1, mlat, mmax, ndo1, ndo2, vb, wb, br, bi, cr, ci, ve, vo, we, wo)
    integer(c_int), intent(in) :: ndgl, nt, imid, imm1, mlat, mmax, ndo1, ndo2
    real, intent(in) :: vb(imid, *), wb(imid, *)
    real, intent(in) :: br(ndgl, ndgl, nt), bi(ndgl, ndgl, nt), cr(ndgl, ndgl, nt), ci(ndgl, ndgl, nt)
    real, intent(inout) :: ve(imid, 2 * mmax, nt), vo(imid, 2 * mmax, nt), we(imid, 2 * mmax, nt), wo(imid, 2 * mmax, nt)

    integer(c_int) :: k, i, mp1, np1, m, mb, mp2, mn
    real :: br_val, bi_val, cr_val, ci_val, vb_val, wb_val

    do k = 1_c_int, nt
      do np1 = 2_c_int, ndo2, 2_c_int
        br_val = br(1, np1, k)
        cr_val = cr(1, np1, k)
        do i = 1_c_int, imid
          vb_val = vb(i, np1)
          ve(i, 1, k) = ve(i, 1, k) + br_val * vb_val
          we(i, 1, k) = we(i, 1, k) - cr_val * vb_val
        end do
      end do
      do np1 = 3_c_int, ndo1, 2_c_int
        br_val = br(1, np1, k)
        cr_val = cr(1, np1, k)
        do i = 1_c_int, imm1
          vb_val = vb(i, np1)
          vo(i, 1, k) = vo(i, 1, k) + br_val * vb_val
          wo(i, 1, k) = wo(i, 1, k) - cr_val * vb_val
        end do
      end do
    end do

    if (mmax >= 2_c_int) then
      do mp1 = 2_c_int, mmax
        m = mp1 - 1_c_int
        mp2 = mp1 + 1_c_int
        mb = m * ndgl - (m * (m + 1_c_int)) / 2_c_int
        if (mp1 <= ndo1) then
          do k = 1_c_int, nt
            do np1 = mp1, ndo1, 2_c_int
              mn = mb + np1
              br_val = br(mp1, np1, k)
              bi_val = bi(mp1, np1, k)
              cr_val = cr(mp1, np1, k)
              ci_val = ci(mp1, np1, k)
              do i = 1_c_int, imm1
                vb_val = vb(i, mn)
                wb_val = wb(i, mn)
                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) + br_val * vb_val
                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb_val
                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + bi_val * vb_val
                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb_val
                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - cr_val * vb_val
                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb_val
                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) - ci_val * vb_val
                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb_val
              end do
              if (mlat /= 0_c_int) then
                i = imid
                wb_val = wb(i, mn)
                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) - ci_val * wb_val
                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + cr_val * wb_val
                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - bi_val * wb_val
                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) + br_val * wb_val
              end if
            end do
          end do
        end if
        if (mp2 <= ndo2) then
          do k = 1_c_int, nt
            do np1 = mp2, ndo2, 2_c_int
              mn = mb + np1
              br_val = br(mp1, np1, k)
              bi_val = bi(mp1, np1, k)
              cr_val = cr(mp1, np1, k)
              ci_val = ci(mp1, np1, k)
              do i = 1_c_int, imm1
                vb_val = vb(i, mn)
                wb_val = wb(i, mn)
                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb_val
                vo(i, 2*mp1-2, k) = vo(i, 2*mp1-2, k) - ci_val * wb_val
                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb_val
                vo(i, 2*mp1-1, k) = vo(i, 2*mp1-1, k) + cr_val * wb_val
                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb_val
                wo(i, 2*mp1-2, k) = wo(i, 2*mp1-2, k) - bi_val * wb_val
                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb_val
                wo(i, 2*mp1-1, k) = wo(i, 2*mp1-1, k) + br_val * wb_val
              end do
              if (mlat /= 0_c_int) then
                i = imid
                vb_val = vb(i, mn)
                ve(i, 2*mp1-2, k) = ve(i, 2*mp1-2, k) + br_val * vb_val
                ve(i, 2*mp1-1, k) = ve(i, 2*mp1-1, k) + bi_val * vb_val
                we(i, 2*mp1-2, k) = we(i, 2*mp1-2, k) - cr_val * vb_val
                we(i, 2*mp1-1, k) = we(i, 2*mp1-1, k) - ci_val * vb_val
              end if
            end do
          end do
        end if
      end do
    end if
  end subroutine accumulate_gaussian_vector_modes

end module vector_backend
