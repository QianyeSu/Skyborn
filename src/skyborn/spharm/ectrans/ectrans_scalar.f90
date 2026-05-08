module scalar_backend
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex
  implicit none
  private

  public :: scalar_analysis_stub
  public :: scalar_synthesis_stub
  public :: scalar_fourier_stub
  public :: scalar_block_solve_stub
  public :: weighted_block_solve_stub

  integer(c_int), save :: basis_cache_ready = 0_c_int
  integer(c_int), save :: cache_ndgl = 0_c_int
  integer(c_int), save :: cache_ngptot = 0_c_int
  integer(c_int), save :: cache_ntrunc = -1_c_int
  integer(c_int), allocatable, save :: cache_nloen(:)
  real(c_double), allocatable, save :: cache_mu(:)
  real(c_double), allocatable, save :: cache_basis_real(:)
  real(c_double), allocatable, save :: cache_basis_imag(:)

contains

  subroutine scalar_block_solve_stub( &
    ndgl, nblock, nt, weights, basis_real, basis_imag, observed_real, observed_imag, &
    solution_real, solution_imag, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl, nblock, nt
    real(c_double), intent(in) :: weights(ndgl)
    real(c_double), intent(in) :: basis_real(ndgl * nblock)
    real(c_double), intent(in) :: basis_imag(ndgl * nblock)
    real(c_double), intent(in) :: observed_real(ndgl * nt)
    real(c_double), intent(in) :: observed_imag(ndgl * nt)
    real(c_double), intent(out) :: solution_real(nblock * nt)
    real(c_double), intent(out) :: solution_imag(nblock * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ilat, iblock, jblock, it, info
    integer(c_int) :: basis_idx, obs_idx, sol_idx
    complex(c_double_complex), allocatable :: gram(:,:), rhs(:,:), basis_row(:), observed_row(:)
    complex(c_double_complex) :: coeff_i
    real(c_double) :: weight

    if (ndgl < 1_c_int .or. nblock < 1_c_int .or. nt < 1_c_int) then
      solution_real = 0.0_c_double
      solution_imag = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    allocate(gram(nblock, nblock))
    allocate(rhs(nblock, nt))
    allocate(basis_row(nblock))
    allocate(observed_row(nt))

    gram = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)
    rhs = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)

    do ilat = 1_c_int, ndgl
      weight = weights(ilat)

      do iblock = 1_c_int, nblock
        basis_idx = (ilat - 1_c_int) * nblock + iblock
        basis_row(iblock) = cmplx(basis_real(basis_idx), basis_imag(basis_idx), kind=c_double_complex)
      end do

      do it = 1_c_int, nt
        obs_idx = (ilat - 1_c_int) * nt + it
        observed_row(it) = cmplx(observed_real(obs_idx), observed_imag(obs_idx), kind=c_double_complex)
      end do

      do iblock = 1_c_int, nblock
        coeff_i = conjg(basis_row(iblock))
        do jblock = 1_c_int, nblock
          gram(iblock, jblock) = gram(iblock, jblock) + weight * coeff_i * basis_row(jblock)
        end do
        do it = 1_c_int, nt
          rhs(iblock, it) = rhs(iblock, it) + weight * coeff_i * observed_row(it)
        end do
      end do
    end do

    call solve_complex_system(gram, rhs, nblock, nt, info)
    if (info /= 0_c_int) then
      solution_real = 0.0_c_double
      solution_imag = 0.0_c_double
      ierror = info
      deallocate(gram, rhs, basis_row, observed_row)
      return
    end if

    do iblock = 1_c_int, nblock
      do it = 1_c_int, nt
        sol_idx = (iblock - 1_c_int) * nt + it
        solution_real(sol_idx) = real(rhs(iblock, it), c_double)
        solution_imag(sol_idx) = aimag(rhs(iblock, it))
      end do
    end do

    deallocate(gram, rhs, basis_row, observed_row)
    ierror = 0_c_int
  end subroutine scalar_block_solve_stub

  subroutine weighted_block_solve_stub( &
    nrow, nblock, nt, weights, basis_real, basis_imag, observed_real, observed_imag, &
    solution_real, solution_imag, ierror) bind(C)
    integer(c_int), value, intent(in) :: nrow, nblock, nt
    real(c_double), intent(in) :: weights(nrow)
    real(c_double), intent(in) :: basis_real(nrow * nblock)
    real(c_double), intent(in) :: basis_imag(nrow * nblock)
    real(c_double), intent(in) :: observed_real(nrow * nt)
    real(c_double), intent(in) :: observed_imag(nrow * nt)
    real(c_double), intent(out) :: solution_real(nblock * nt)
    real(c_double), intent(out) :: solution_imag(nblock * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: irow, iblock, jblock, it, info
    integer(c_int) :: basis_idx, obs_idx, sol_idx
    complex(c_double_complex), allocatable :: gram(:,:), rhs(:,:), basis_row(:), observed_row(:)
    complex(c_double_complex) :: coeff_i
    real(c_double) :: weight

    if (nrow < 1_c_int .or. nblock < 1_c_int .or. nt < 1_c_int) then
      solution_real = 0.0_c_double
      solution_imag = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    allocate(gram(nblock, nblock))
    allocate(rhs(nblock, nt))
    allocate(basis_row(nblock))
    allocate(observed_row(nt))

    gram = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)
    rhs = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)

    do irow = 1_c_int, nrow
      weight = weights(irow)

      do iblock = 1_c_int, nblock
        basis_idx = (irow - 1_c_int) * nblock + iblock
        basis_row(iblock) = cmplx(basis_real(basis_idx), basis_imag(basis_idx), kind=c_double_complex)
      end do

      do it = 1_c_int, nt
        obs_idx = (irow - 1_c_int) * nt + it
        observed_row(it) = cmplx(observed_real(obs_idx), observed_imag(obs_idx), kind=c_double_complex)
      end do

      do iblock = 1_c_int, nblock
        coeff_i = conjg(basis_row(iblock))
        do jblock = 1_c_int, nblock
          gram(iblock, jblock) = gram(iblock, jblock) + weight * coeff_i * basis_row(jblock)
        end do
        do it = 1_c_int, nt
          rhs(iblock, it) = rhs(iblock, it) + weight * coeff_i * observed_row(it)
        end do
      end do
    end do

    call solve_complex_system(gram, rhs, nblock, nt, info)
    if (info /= 0_c_int) then
      solution_real = 0.0_c_double
      solution_imag = 0.0_c_double
      ierror = info
      deallocate(gram, rhs, basis_row, observed_row)
      return
    end if

    do iblock = 1_c_int, nblock
      do it = 1_c_int, nt
        sol_idx = (iblock - 1_c_int) * nt + it
        solution_real(sol_idx) = real(rhs(iblock, it), c_double)
        solution_imag(sol_idx) = aimag(rhs(iblock, it))
      end do
    end do

    deallocate(gram, rhs, basis_row, observed_row)
    ierror = 0_c_int
  end subroutine weighted_block_solve_stub

  subroutine scalar_fourier_stub( &
    ndgl, nloen, ngptot, mmax, nt, datagrid, fourier_real, fourier_imag, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot, mmax, nt
    real(c_double), intent(in) :: datagrid(ngptot * nt)
    real(c_double), intent(out) :: fourier_real(ndgl * (mmax + 1_c_int) * nt)
    real(c_double), intent(out) :: fourier_imag(ndgl * (mmax + 1_c_int) * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ilat, it, m, j, nlon_lat, offset, idx
    real(c_double) :: lon_step, angle, value
    complex(c_double_complex) :: accum, phase

    if (ndgl < 1_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int .or. mmax < 0_c_int) then
      fourier_real = 0.0_c_double
      fourier_imag = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    fourier_real = 0.0_c_double
    fourier_imag = 0.0_c_double
    offset = 0_c_int

    do ilat = 1_c_int, ndgl
      nlon_lat = nloen(ilat)
      lon_step = 2.0_c_double * acos(-1.0_c_double) / real(nlon_lat, c_double)

      do it = 1_c_int, nt
        do m = 0_c_int, mmax
          accum = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)
          if (m < nlon_lat) then
            do j = 0_c_int, nlon_lat - 1_c_int
              angle = -real(m * j, c_double) * lon_step
              phase = cmplx(cos(angle), sin(angle), kind=c_double_complex)
              value = datagrid((offset + j) * nt + it)
              accum = accum + cmplx(value, 0.0_c_double, kind=c_double_complex) * phase
            end do
            accum = accum / real(nlon_lat, c_double)
          end if
          idx = ((ilat - 1_c_int) * (mmax + 1_c_int) + m) * nt + it
          fourier_real(idx) = real(accum, c_double)
          fourier_imag(idx) = aimag(accum)
        end do
      end do

      offset = offset + nlon_lat
    end do

    ierror = 0_c_int
  end subroutine scalar_fourier_stub

  subroutine scalar_analysis_stub( &
    ndgl, nloen, mu, weights, ngptot, ntrunc, nt, datagrid, dataspec_packed, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    real(c_double), intent(in) :: mu(ndgl)
    real(c_double), intent(in) :: weights(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: datagrid(ngptot * nt)
    real(c_double), intent(out) :: dataspec_packed(((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ncoeff, nspec2, m, nblock, start_idx
    integer(c_int) :: ilat, icol, it, coeff_index, ierr_local
    real(c_double), allocatable :: observed_real(:), observed_imag(:)
    real(c_double), allocatable :: basis_block_real(:), basis_block_imag(:)
    real(c_double), allocatable :: observed_block_real(:), observed_block_imag(:)
    real(c_double), allocatable :: solution_real(:), solution_imag(:)

    if (ndgl < 1_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int .or. ntrunc < 0_c_int) then
      dataspec_packed = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    nspec2 = 2_c_int * ncoeff

    call ensure_scalar_basis_cache(ndgl, nloen, mu, ngptot, ntrunc, ierr_local)
    if (ierr_local /= 0_c_int) then
      dataspec_packed = 0.0_c_double
      ierror = ierr_local
      return
    end if

    allocate(observed_real(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(observed_imag(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(basis_block_real(ndgl * (ntrunc + 1_c_int)))
    allocate(basis_block_imag(ndgl * (ntrunc + 1_c_int)))
    allocate(observed_block_real(ndgl * nt))
    allocate(observed_block_imag(ndgl * nt))
    allocate(solution_real((ntrunc + 1_c_int) * nt))
    allocate(solution_imag((ntrunc + 1_c_int) * nt))

    dataspec_packed = 0.0_c_double

    call scalar_fourier_stub( &
      ndgl, nloen, ngptot, ntrunc, nt, datagrid, observed_real, observed_imag, ierr_local)
    if (ierr_local /= 0_c_int) then
      dataspec_packed = 0.0_c_double
      ierror = ierr_local
      deallocate( &
        observed_real, observed_imag, basis_block_real, basis_block_imag, &
        observed_block_real, observed_block_imag, solution_real, solution_imag)
      return
    end if

    start_idx = 1_c_int
    do m = 0_c_int, ntrunc
      nblock = ntrunc - m + 1_c_int

      do ilat = 1_c_int, ndgl
        do icol = 1_c_int, nblock
          coeff_index = start_idx + icol - 1_c_int
          basis_block_real((ilat - 1_c_int) * nblock + icol) = &
            cache_basis_real((ilat - 1_c_int) * ncoeff + coeff_index)
          basis_block_imag((ilat - 1_c_int) * nblock + icol) = &
            cache_basis_imag((ilat - 1_c_int) * ncoeff + coeff_index)
        end do
        do it = 1_c_int, nt
          observed_block_real((ilat - 1_c_int) * nt + it) = &
            observed_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
          observed_block_imag((ilat - 1_c_int) * nt + it) = &
            observed_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
        end do
      end do

      call scalar_block_solve_stub( &
        ndgl, nblock, nt, weights, &
        basis_block_real(1:ndgl * nblock), &
        basis_block_imag(1:ndgl * nblock), &
        observed_block_real, observed_block_imag, &
        solution_real(1:nblock * nt), &
        solution_imag(1:nblock * nt), ierr_local)
      if (ierr_local /= 0_c_int) then
        dataspec_packed = 0.0_c_double
        ierror = ierr_local
        deallocate( &
          observed_real, observed_imag, basis_block_real, basis_block_imag, &
          observed_block_real, observed_block_imag, solution_real, solution_imag)
        return
      end if

      do icol = 1_c_int, nblock
        coeff_index = start_idx + icol - 1_c_int
        do it = 1_c_int, nt
          dataspec_packed((2_c_int * coeff_index - 2_c_int) * nt + it) = &
            solution_real((icol - 1_c_int) * nt + it)
          dataspec_packed((2_c_int * coeff_index - 1_c_int) * nt + it) = &
            solution_imag((icol - 1_c_int) * nt + it)
        end do
      end do

      start_idx = start_idx + nblock
    end do

    deallocate( &
      observed_real, observed_imag, basis_block_real, basis_block_imag, &
      observed_block_real, observed_block_imag, solution_real, solution_imag)
    ierror = 0_c_int
  end subroutine scalar_analysis_stub

  subroutine scalar_synthesis_stub( &
    ndgl, nloen, mu, ngptot, ntrunc, nt, dataspec_packed, datagrid, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    real(c_double), intent(in) :: mu(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: dataspec_packed(((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) * nt)
    real(c_double), intent(out) :: datagrid(ngptot * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ilat, ilon, it, m, nblock, coeff_index, ncoeff
    integer(c_int) :: nlon_lat, offset, start_idx, icol
    real(c_double) :: lon_rad, lon_step, basis_val, value
    real(c_double) :: cos_angle, sin_angle, two_pi
    real(c_double), allocatable :: fourier_real(:), fourier_imag(:)
    integer(c_int) :: ierr_local

    if (ndgl < 3_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int) then
      datagrid = 0.0_c_double
      ierror = 1_c_int
      return
    end if

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    two_pi = 2.0_c_double * acos(-1.0_c_double)

    call ensure_scalar_basis_cache(ndgl, nloen, mu, ngptot, ntrunc, ierr_local)
    if (ierr_local /= 0_c_int) then
      datagrid = 0.0_c_double
      ierror = ierr_local
      return
    end if

    allocate(fourier_real(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(fourier_imag(ndgl * (ntrunc + 1_c_int) * nt))

    datagrid = 0.0_c_double
    fourier_real = 0.0_c_double
    fourier_imag = 0.0_c_double

    start_idx = 1_c_int
    do m = 0_c_int, ntrunc
      nblock = ntrunc - m + 1_c_int
      do ilat = 1_c_int, ndgl
        do icol = 1_c_int, nblock
          coeff_index = start_idx + icol - 1_c_int
          basis_val = cache_basis_real((ilat - 1_c_int) * ncoeff + coeff_index)
          do it = 1_c_int, nt
            fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) = &
              fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) + &
              basis_val * dataspec_packed((2_c_int * coeff_index - 2_c_int) * nt + it)
            fourier_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) = &
              fourier_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) + &
              basis_val * dataspec_packed((2_c_int * coeff_index - 1_c_int) * nt + it)
          end do
        end do
      end do
      start_idx = start_idx + nblock
    end do

    do it = 1_c_int, nt
      offset = 0_c_int
      do ilat = 1_c_int, ndgl
        nlon_lat = nloen(ilat)
        lon_step = two_pi / real(nlon_lat, c_double)
        do ilon = 0_c_int, nlon_lat - 1_c_int
          lon_rad = real(ilon, c_double) * lon_step
          value = fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int)) * nt + it)
          do m = 1_c_int, ntrunc
            cos_angle = cos(real(m, c_double) * lon_rad)
            sin_angle = sin(real(m, c_double) * lon_rad)
            value = value + 2.0_c_double * ( &
              fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) * cos_angle - &
              fourier_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it) * sin_angle &
            )
          end do
          datagrid((offset + ilon) * nt + it) = value
        end do
        offset = offset + nlon_lat
      end do
    end do

    deallocate(fourier_real, fourier_imag)
    ierror = 0_c_int
  end subroutine scalar_synthesis_stub

  pure real(c_double) function mu_to_lat(mu) result(lat_deg)
    real(c_double), intent(in) :: mu
    lat_deg = asin(max(-1.0_c_double, min(1.0_c_double, mu))) * 180.0_c_double / acos(-1.0_c_double)
  end function mu_to_lat

  subroutine build_legfunc_triangle(legfunc, lat_deg, ntrunc)
    integer(c_int), intent(in) :: ntrunc
    real(c_double), intent(in) :: lat_deg
    real(c_double), intent(out) :: legfunc((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int)

    integer(c_int) :: m, n, nm, nmstrt
    real(c_double) :: theta, pi
    real(c_double), allocatable :: cp(:)

    pi = 4.0_c_double * atan(1.0_c_double)
    theta = 0.5_c_double * pi - (pi / 180.0_c_double) * lat_deg

    legfunc = 0.0_c_double
    nmstrt = 0_c_int
    allocate(cp((ntrunc / 2_c_int) + 1_c_int))

    do m = 1_c_int, ntrunc + 1_c_int
      do n = m, ntrunc + 1_c_int
        nm = nmstrt + n - m + 1_c_int
        call alfk_local(n - 1_c_int, m - 1_c_int, cp)
        call lfpt_local(n - 1_c_int, m - 1_c_int, theta, cp, legfunc(nm))
      end do
      nmstrt = nmstrt + ntrunc - m + 2_c_int
    end do

    deallocate(cp)
  end subroutine build_legfunc_triangle

  subroutine alfk_local(n, m, cp)
    integer(c_int), intent(in) :: n, m
    real(c_double), intent(out) :: cp(n / 2_c_int + 1_c_int)

    integer(c_int) :: ma, nmms2, i, l
    real(c_double) :: fnum, fnmh, pm1, t1, t2, fden, nex
    real(c_double) :: cp2, fnnp1, fnmsq, fk, a1, b1, c1
    real(c_double), parameter :: sc10 = 1024.0_c_double
    real(c_double), parameter :: sc20 = sc10 * sc10
    real(c_double), parameter :: sc40 = sc20 * sc20

    cp(1) = 0.0_c_double
    ma = abs(m)
    if (ma > n) then
      return
    end if

    select case (n)
    case (0_c_int)
      cp(1) = sqrt(2.0_c_double)
      return
    case (1_c_int)
      if (ma == 0_c_int) then
        cp(1) = sqrt(1.5_c_double)
      else
        cp(1) = sqrt(0.75_c_double)
        if (m == -1_c_int) then
          cp(1) = -cp(1)
        end if
      end if
      return
    end select

    if (mod(n + ma, 2_c_int) == 0_c_int) then
      nmms2 = (n - ma) / 2_c_int
      fnum = real(n + ma + 1_c_int, c_double)
      fnmh = real(n - ma + 1_c_int, c_double)
      pm1 = 1.0_c_double
    else
      nmms2 = (n - ma - 1_c_int) / 2_c_int
      fnum = real(n + ma + 2_c_int, c_double)
      fnmh = real(n - ma + 2_c_int, c_double)
      pm1 = -1.0_c_double
    end if

    t1 = 1.0_c_double / sc20
    nex = 20.0_c_double
    fden = 2.0_c_double

    if (nmms2 >= 1_c_int) then
      do i = 1_c_int, nmms2
        t1 = fnum * t1 / fden
        if (t1 > sc20) then
          t1 = t1 / sc40
          nex = nex + 40.0_c_double
        end if
        fnum = fnum + 2.0_c_double
        fden = fden + 2.0_c_double
      end do
    end if

    t1 = t1 / (2.0_c_double ** real(n - 1_c_int - int(nex, c_int), c_double))
    if (mod(ma / 2_c_int, 2_c_int) /= 0_c_int) then
      t1 = -t1
    end if

    t2 = 1.0_c_double
    if (ma > 0_c_int) then
      do i = 1_c_int, ma
        t2 = fnmh * t2 / (fnmh + pm1)
        fnmh = fnmh + 2.0_c_double
      end do
    end if

    cp2 = t1 * sqrt((real(n, c_double) + 0.5_c_double) * t2)
    fnnp1 = real(n * (n + 1_c_int), c_double)
    fnmsq = fnnp1 - 2.0_c_double * real(ma * ma, c_double)
    l = (n + 1_c_int) / 2_c_int
    if (mod(n, 2_c_int) == 0_c_int .and. mod(ma, 2_c_int) == 0_c_int) then
      l = l + 1_c_int
    end if
    cp(l) = cp2

    if (m < 0_c_int .and. mod(ma, 2_c_int) /= 0_c_int) then
      cp(l) = -cp(l)
    end if

    if (l <= 1_c_int) then
      return
    end if

    fk = real(n, c_double)
    a1 = (fk - 2.0_c_double) * (fk - 1.0_c_double) - fnnp1
    b1 = 2.0_c_double * (fk * fk - fnmsq)
    cp(l - 1_c_int) = b1 * cp(l) / a1

    do while (l > 2_c_int)
      l = l - 1_c_int
      fk = fk - 2.0_c_double
      a1 = (fk - 2.0_c_double) * (fk - 1.0_c_double) - fnnp1
      b1 = -2.0_c_double * (fk * fk - fnmsq)
      c1 = (fk + 1.0_c_double) * (fk + 2.0_c_double) - fnnp1
      cp(l - 1_c_int) = -(b1 * cp(l) + c1 * cp(l + 1_c_int)) / a1
    end do
  end subroutine alfk_local

  subroutine lfpt_local(n, m, theta, cp, pb)
    integer(c_int), intent(in) :: n, m
    real(c_double), intent(in) :: theta
    real(c_double), intent(in) :: cp(*)
    real(c_double), intent(out) :: pb

    integer(c_int) :: ma, np1, nmod, mmod, kdo, k
    real(c_double) :: cdt, sdt, ct, st, sum, cth

    pb = 0.0_c_double
    ma = abs(m)
    if (ma > n) then
      return
    end if

    if (n <= 0_c_int) then
      if (ma <= 0_c_int) then
        pb = sqrt(0.5_c_double)
      end if
      return
    end if

    np1 = n + 1_c_int
    nmod = mod(n, 2_c_int)
    mmod = mod(ma, 2_c_int)

    if (nmod == 0_c_int) then
      if (mmod == 0_c_int) then
        kdo = n / 2_c_int + 1_c_int
        cdt = cos(theta + theta)
        sdt = sin(theta + theta)
        ct = 1.0_c_double
        st = 0.0_c_double
        sum = 0.5_c_double * cp(1)
        do k = 2_c_int, kdo
          cth = cdt * ct - sdt * st
          st = sdt * ct + cdt * st
          ct = cth
          sum = sum + cp(k) * ct
        end do
        pb = sum
      else
        kdo = n / 2_c_int
        cdt = cos(theta + theta)
        sdt = sin(theta + theta)
        ct = 1.0_c_double
        st = 0.0_c_double
        sum = 0.0_c_double
        do k = 1_c_int, kdo
          cth = cdt * ct - sdt * st
          st = sdt * ct + cdt * st
          ct = cth
          sum = sum + cp(k) * st
        end do
        pb = sum
      end if
    else
      kdo = (n + 1_c_int) / 2_c_int
      cdt = cos(theta + theta)
      sdt = sin(theta + theta)
      ct = cos(theta)
      st = -sin(theta)
      sum = 0.0_c_double
      if (mmod == 0_c_int) then
        do k = 1_c_int, kdo
          cth = cdt * ct - sdt * st
          st = sdt * ct + cdt * st
          ct = cth
          sum = sum + cp(k) * ct
        end do
      else
        do k = 1_c_int, kdo
          cth = cdt * ct - sdt * st
          st = sdt * ct + cdt * st
          ct = cth
          sum = sum + cp(k) * st
        end do
      end if
      pb = sum
    end if
  end subroutine lfpt_local

  subroutine solve_complex_system(matrix, rhs, nsize, nrhs, info)
    integer(c_int), intent(in) :: nsize, nrhs
    complex(c_double_complex), intent(inout) :: matrix(nsize, nsize)
    complex(c_double_complex), intent(inout) :: rhs(nsize, nrhs)
    integer(c_int), intent(out) :: info

    integer(c_int) :: ipivot, jpivot, irhs, pivot_row
    complex(c_double_complex) :: factor, temp
    real(c_double) :: pivot_abs

    info = 0_c_int

    do ipivot = 1_c_int, nsize
      pivot_row = ipivot
      pivot_abs = abs(matrix(ipivot, ipivot))
      do jpivot = ipivot + 1_c_int, nsize
        if (abs(matrix(jpivot, ipivot)) > pivot_abs) then
          pivot_abs = abs(matrix(jpivot, ipivot))
          pivot_row = jpivot
        end if
      end do

      if (pivot_abs <= tiny(1.0_c_double)) then
        info = 10_c_int + ipivot
        return
      end if

      if (pivot_row /= ipivot) then
        do jpivot = ipivot, nsize
          temp = matrix(ipivot, jpivot)
          matrix(ipivot, jpivot) = matrix(pivot_row, jpivot)
          matrix(pivot_row, jpivot) = temp
        end do
        do irhs = 1_c_int, nrhs
          temp = rhs(ipivot, irhs)
          rhs(ipivot, irhs) = rhs(pivot_row, irhs)
          rhs(pivot_row, irhs) = temp
        end do
      end if

      do jpivot = ipivot + 1_c_int, nsize
        factor = matrix(jpivot, ipivot) / matrix(ipivot, ipivot)
        matrix(jpivot, ipivot) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double_complex)
        do irhs = ipivot + 1_c_int, nsize
          matrix(jpivot, irhs) = matrix(jpivot, irhs) - factor * matrix(ipivot, irhs)
        end do
        do irhs = 1_c_int, nrhs
          rhs(jpivot, irhs) = rhs(jpivot, irhs) - factor * rhs(ipivot, irhs)
        end do
      end do
    end do

    do ipivot = nsize, 1_c_int, -1_c_int
      if (abs(matrix(ipivot, ipivot)) <= tiny(1.0_c_double)) then
        info = 100_c_int + ipivot
        return
      end if
      do irhs = 1_c_int, nrhs
        temp = rhs(ipivot, irhs)
        do jpivot = ipivot + 1_c_int, nsize
          temp = temp - matrix(ipivot, jpivot) * rhs(jpivot, irhs)
        end do
        rhs(ipivot, irhs) = temp / matrix(ipivot, ipivot)
      end do
    end do
  end subroutine solve_complex_system

  logical function scalar_basis_cache_matches(ndgl, nloen, mu, ngptot, ntrunc) result(matches)
    integer(c_int), intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    real(c_double), intent(in) :: mu(ndgl)
    integer(c_int), intent(in) :: ngptot, ntrunc

    matches = .false.
    if (basis_cache_ready == 0_c_int) then
      return
    end if
    if (.not. allocated(cache_nloen) .or. .not. allocated(cache_mu)) then
      return
    end if
    if (cache_ndgl /= ndgl .or. cache_ngptot /= ngptot .or. cache_ntrunc /= ntrunc) then
      return
    end if
    if (.not. all(cache_nloen == nloen)) then
      return
    end if
    if (.not. all(abs(cache_mu - mu) <= 1.0e-14_c_double)) then
      return
    end if
    matches = .true.
  end function scalar_basis_cache_matches

  subroutine reset_scalar_basis_cache()
    if (allocated(cache_nloen)) then
      deallocate(cache_nloen)
    end if
    if (allocated(cache_mu)) then
      deallocate(cache_mu)
    end if
    if (allocated(cache_basis_real)) then
      deallocate(cache_basis_real)
    end if
    if (allocated(cache_basis_imag)) then
      deallocate(cache_basis_imag)
    end if
    basis_cache_ready = 0_c_int
    cache_ndgl = 0_c_int
    cache_ngptot = 0_c_int
    cache_ntrunc = -1_c_int
  end subroutine reset_scalar_basis_cache

  subroutine ensure_scalar_basis_cache(ndgl, nloen, mu, ngptot, ntrunc, ierror)
    integer(c_int), intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    real(c_double), intent(in) :: mu(ndgl)
    integer(c_int), intent(in) :: ngptot, ntrunc
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ncoeff, m, nblock, start_idx
    integer(c_int) :: icol, ilat, coeff_index
    real(c_double) :: lat_deg
    real(c_double), allocatable :: legfunc_all(:,:)

    if (scalar_basis_cache_matches(ndgl, nloen, mu, ngptot, ntrunc)) then
      ierror = 0_c_int
      return
    end if

    call reset_scalar_basis_cache()

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int

    allocate(cache_nloen(ndgl))
    allocate(cache_mu(ndgl))
    allocate(cache_basis_real(ndgl * ncoeff))
    allocate(cache_basis_imag(ndgl * ncoeff))
    cache_nloen = nloen
    cache_mu = mu
    cache_basis_real = 0.0_c_double
    cache_basis_imag = 0.0_c_double

    allocate(legfunc_all(ncoeff, ndgl))
    do ilat = 1_c_int, ndgl
      lat_deg = mu_to_lat(mu(ilat))
      call build_legfunc_triangle(legfunc_all(:, ilat), lat_deg, ntrunc)
    end do

    start_idx = 1_c_int
    do m = 0_c_int, ntrunc
      nblock = ntrunc - m + 1_c_int
      do ilat = 1_c_int, ndgl
        do icol = 1_c_int, nblock
          coeff_index = start_idx + icol - 1_c_int
          cache_basis_real((ilat - 1_c_int) * ncoeff + coeff_index) = &
            real(legfunc_all(coeff_index, ilat), c_double)
          cache_basis_imag((ilat - 1_c_int) * ncoeff + coeff_index) = 0.0_c_double
        end do
      end do

      start_idx = start_idx + nblock
    end do

    deallocate(legfunc_all)

    basis_cache_ready = 1_c_int
    cache_ndgl = ndgl
    cache_ngptot = ngptot
    cache_ntrunc = ntrunc
    ierror = 0_c_int
  end subroutine ensure_scalar_basis_cache

end module scalar_backend
