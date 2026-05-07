module vector_backend
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  implicit none
  private

  public :: vrtdiv_analysis_stub
  public :: uv_synthesis_stub
  public :: gradient_synthesis_stub

  interface
    subroutine vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierror)
      integer, intent(in) :: nlat, nlon, lvhsgs, ldwork
      real, intent(out) :: wvhsgs(*)
      double precision, intent(inout) :: dwork(*)
      integer, intent(out) :: ierror
    end subroutine vhsgsi
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
    ndgl, nloen, ngptot, ntrunc, nt, vrtspec_r, vrtspec_i, divspec_r, divspec_i, ugrid, vgrid, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: vrtspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: vrtspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: ugrid(ngptot * nt)
    real(c_double), intent(out) :: vgrid(ngptot * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: imid, lmn, lvhsgs, ldwork, mmax, ndo1, ndo2, mlat, imm1
    integer(c_int) :: idz, ilat, it, j, m, mp1, np1, mb, mn, nlon_lat, offset
    integer(c_int) :: ncoeff, ierr_local
    real :: br_val, bi_val, cr_val, ci_val, vb_val, wb_val
    real :: ve_val, vo_val, we_val, wo_val
    real(c_double) :: twopi, lon_step, lon
    real(c_double), allocatable :: cos_table(:,:), sin_table(:,:)
    real, allocatable :: wvhsgs(:), vb(:,:), wb(:,:)
    real, allocatable :: br(:,:,:), bi(:,:,:), cr(:,:,:), ci(:,:,:)
    real, allocatable :: ve(:,:,:), vo(:,:,:), we(:,:,:), wo(:,:,:)
    double precision, allocatable :: dwork(:)

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

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    imid = (ndgl + 1_c_int) / 2_c_int
    lmn = ndgl * (ndgl + 1_c_int) / 2_c_int
    lvhsgs = 2_c_int * (imid * lmn) + maxval(nloen) + 15_c_int
    ldwork = (ndgl * 3_c_int * (ndgl + 3_c_int) + 2_c_int) / 2_c_int
    idz = imid * lmn
    mmax = min(ndgl, (maxval(nloen) + 1_c_int) / 2_c_int)
    mlat = mod(ndgl, 2_c_int)
    imm1 = imid
    if (mlat /= 0_c_int) imm1 = imid - 1_c_int
    ndo1 = ndgl
    ndo2 = ndgl
    if (mlat /= 0_c_int) ndo1 = ndgl - 1_c_int
    if (mlat == 0_c_int) ndo2 = ndgl - 1_c_int

    allocate(wvhsgs(lvhsgs))
    allocate(dwork(ldwork))
    allocate(vb(imid, idz))
    allocate(wb(imid, idz))
    allocate(br(ndgl, ndgl, nt), bi(ndgl, ndgl, nt), cr(ndgl, ndgl, nt), ci(ndgl, ndgl, nt))
    allocate(ve(imid, 2 * mmax, nt), vo(imid, 2 * mmax, nt), we(imid, 2 * mmax, nt), wo(imid, 2 * mmax, nt))

    call vhsgsi(ndgl, maxval(nloen), wvhsgs, lvhsgs, dwork, ldwork, ierr_local)
    if (ierr_local /= 0_c_int) then
      ugrid = 0.0_c_double
      vgrid = 0.0_c_double
      ierror = 100_c_int + ierr_local
      deallocate(wvhsgs, dwork, vb, wb, br, bi, cr, ci, ve, vo, we, wo)
      return
    end if

    vb(:, :) = reshape(wvhsgs(1:idz), shape(vb))
    wb(:, :) = reshape(wvhsgs(idz + 1:2 * idz), shape(wb))

    br = 0.0
    bi = 0.0
    cr = 0.0
    ci = 0.0
    call unpack_vrtdiv_to_vector_harmonics(ndgl, ntrunc, nt, vrtspec_r, vrtspec_i, divspec_r, divspec_i, br, bi, cr, ci)

    ve = 0.0
    vo = 0.0
    we = 0.0
    wo = 0.0
    call accumulate_gaussian_vector_modes( &
      ndgl, nt, imid, imm1, mlat, mmax, ndo1, ndo2, vb, wb, br, bi, cr, ci, ve, vo, we, wo)

    ugrid = 0.0_c_double
    vgrid = 0.0_c_double
    offset = 0_c_int
    twopi = 2.0_c_double * acos(-1.0_c_double)

    do ilat = 1_c_int, ndgl
      nlon_lat = nloen(ilat)
      allocate(cos_table(0:mmax-1, nlon_lat), sin_table(0:mmax-1, nlon_lat))
      lon_step = twopi / real(nlon_lat, c_double)
      do j = 1_c_int, nlon_lat
        lon = real(j - 1_c_int, c_double) * lon_step
        cos_table(0, j) = 1.0_c_double
        sin_table(0, j) = 0.0_c_double
        do m = 1_c_int, mmax - 1_c_int
          cos_table(m, j) = cos(real(m, c_double) * lon)
          sin_table(m, j) = sin(real(m, c_double) * lon)
        end do
      end do

      do it = 1_c_int, nt
        do j = 1_c_int, nlon_lat
          ve_val = 0.0
          vo_val = 0.0
          we_val = 0.0
          wo_val = 0.0
          do mp1 = 1_c_int, mmax
            ve_val = ve_val + ve(min(ilat, imid), 2 * mp1 - 1_c_int, it) * real(cos_table(mp1 - 1_c_int, j), kind=kind(ve_val))
            we_val = we_val + we(min(ilat, imid), 2 * mp1 - 1_c_int, it) * real(cos_table(mp1 - 1_c_int, j), kind=kind(we_val))
            if (mp1 > 1_c_int) then
              ve_val = ve_val + ve(min(ilat, imid), 2 * mp1 - 2_c_int, it) * real(sin_table(mp1 - 1_c_int, j), kind=kind(ve_val))
              we_val = we_val + we(min(ilat, imid), 2 * mp1 - 2_c_int, it) * real(sin_table(mp1 - 1_c_int, j), kind=kind(we_val))
              vo_val = vo_val + vo(min(ilat, imid), 2 * mp1 - 1_c_int, it) * real(cos_table(mp1 - 1_c_int, j), kind=kind(vo_val))
              vo_val = vo_val + vo(min(ilat, imid), 2 * mp1 - 2_c_int, it) * real(sin_table(mp1 - 1_c_int, j), kind=kind(vo_val))
              wo_val = wo_val + wo(min(ilat, imid), 2 * mp1 - 1_c_int, it) * real(cos_table(mp1 - 1_c_int, j), kind=kind(wo_val))
              wo_val = wo_val + wo(min(ilat, imid), 2 * mp1 - 2_c_int, it) * real(sin_table(mp1 - 1_c_int, j), kind=kind(wo_val))
            end if
          end do

          if (ilat < imid .or. (mlat == 0_c_int .and. ilat <= imid)) then
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(ve_val + vo_val, c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(we_val + wo_val, c_double)
          else if (ilat == imid .and. mlat /= 0_c_int) then
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(ve_val, c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(we_val, c_double)
          else
            ugrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(ve_val - vo_val, c_double)
            vgrid((offset + j - 1_c_int) * nt + it) = 0.5_c_double * real(we_val - wo_val, c_double)
          end if
        end do
      end do
      offset = offset + nlon_lat
      deallocate(cos_table, sin_table)
    end do

    deallocate(wvhsgs, dwork, vb, wb, br, bi, cr, ci, ve, vo, we, wo)
    ierror = 0_c_int
  end subroutine uv_synthesis_stub

  subroutine gradient_synthesis_stub( &
    ndgl, nloen, ngptot, ntrunc, nt, chispec_r, chispec_i, ugrad, vgrad, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: chispec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(in) :: chispec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: ugrad(ngptot * nt)
    real(c_double), intent(out) :: vgrad(ngptot * nt)
    integer(c_int), intent(out) :: ierror

    ugrad = 0.0_c_double
    vgrad = 0.0_c_double
    ierror = -1_c_int
  end subroutine gradient_synthesis_stub

  subroutine unpack_vrtdiv_to_vector_harmonics( &
    ndgl, ntrunc, nt, vrtspec_r, vrtspec_i, divspec_r, divspec_i, br, bi, cr, ci)
    integer(c_int), intent(in) :: ndgl, ntrunc, nt
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
            factor = 1.0_c_double / sqrt(nreal * (nreal + 1.0_c_double))
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
