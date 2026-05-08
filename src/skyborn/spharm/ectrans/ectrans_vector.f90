module vector_backend
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_double_complex
  use scalar_backend, only : scalar_fourier_stub
  implicit none
  private

  public :: vrtdiv_analysis_stub
  public :: uv_synthesis_stub
  public :: gradient_synthesis_stub

  integer(c_int), save :: vector_basis_cache_ready = 0_c_int
  integer(c_int), save :: vector_cache_ndgl = 0_c_int
  integer(c_int), save :: vector_cache_ngptot = 0_c_int
  integer(c_int), save :: vector_cache_ntrunc = -1_c_int
  integer(c_int), allocatable, save :: vector_cache_nloen(:)
  real(c_double), save :: vector_cache_rsphere = 0.0_c_double
  real(c_double), allocatable, save :: vector_cache_u_basis_real(:)
  real(c_double), allocatable, save :: vector_cache_u_basis_imag(:)
  real(c_double), allocatable, save :: vector_cache_v_basis_real(:)
  real(c_double), allocatable, save :: vector_cache_v_basis_imag(:)

  interface
    subroutine gaqd(nlat, theta, wts, dwork, ldwork, ierror)
      integer, intent(in) :: nlat, ldwork
      double precision, intent(out) :: theta(nlat), wts(nlat)
      double precision, intent(inout) :: dwork(ldwork)
      integer, intent(out) :: ierror
    end subroutine gaqd

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
    ndgl, nloen, weights, ngptot, rsphere, ntrunc, nt, ugrid, vgrid, vrtspec_r, vrtspec_i, divspec_r, divspec_i, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    real(c_double), intent(in) :: weights(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    real(c_double), value, intent(in) :: rsphere
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), intent(in) :: ugrid(ngptot * nt)
    real(c_double), intent(in) :: vgrid(ngptot * nt)
    real(c_double), intent(out) :: vrtspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: vrtspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_r((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_i((((ntrunc + 1_c_int) * (ntrunc + 2_c_int)) / 2_c_int) * nt)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ncoeff, nt_basis, m, nblock, ncol, start_idx, coeff_index, ierr_local
    integer(c_int) :: ilat, icol, it, base_idx, u_row_offset, v_row_offset
    integer(c_int) :: active_count, active_idx
    real(c_double) :: column_norm
    real(c_double), allocatable :: u_fourier_real(:), u_fourier_imag(:)
    real(c_double), allocatable :: v_fourier_real(:), v_fourier_imag(:)
    real(c_double), allocatable :: basis_block_real(:), basis_block_imag(:)
    real(c_double), allocatable :: observed_block_real(:), observed_block_imag(:)
    real(c_double), allocatable :: solution_real(:), solution_imag(:)
    real(c_double), allocatable :: vector_weights(:)
    real(c_double), allocatable :: reduced_basis_real(:), reduced_basis_imag(:)
    real(c_double), allocatable :: reduced_solution_real(:), reduced_solution_imag(:)
    integer(c_int), allocatable :: active_columns(:)

    vrtspec_r = 0.0_c_double
    vrtspec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double

    if (ndgl < 3_c_int .or. ngptot < 1_c_int .or. nt < 1_c_int .or. ntrunc < 0_c_int) then
      ierror = 1_c_int
      return
    end if
    if (ntrunc > ndgl - 1_c_int) then
      ierror = 2_c_int
      return
    end if

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    nt_basis = 2_c_int * ncoeff

    call ensure_vector_basis_cache(ndgl, nloen, ngptot, ntrunc, rsphere, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      return
    end if

    allocate(u_fourier_real(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(u_fourier_imag(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(v_fourier_real(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(v_fourier_imag(ndgl * (ntrunc + 1_c_int) * nt))
    allocate(basis_block_real((2_c_int * ndgl) * (2_c_int * (ntrunc + 1_c_int))))
    allocate(basis_block_imag((2_c_int * ndgl) * (2_c_int * (ntrunc + 1_c_int))))
    allocate(observed_block_real((2_c_int * ndgl) * nt))
    allocate(observed_block_imag((2_c_int * ndgl) * nt))
    allocate(solution_real((2_c_int * (ntrunc + 1_c_int)) * nt))
    allocate(solution_imag((2_c_int * (ntrunc + 1_c_int)) * nt))
    allocate(vector_weights(2_c_int * ndgl))
    allocate(reduced_basis_real((2_c_int * ndgl) * (2_c_int * (ntrunc + 1_c_int))))
    allocate(reduced_basis_imag((2_c_int * ndgl) * (2_c_int * (ntrunc + 1_c_int))))
    allocate(reduced_solution_real((2_c_int * (ntrunc + 1_c_int)) * nt))
    allocate(reduced_solution_imag((2_c_int * (ntrunc + 1_c_int)) * nt))
    allocate(active_columns(2_c_int * (ntrunc + 1_c_int)))

    call scalar_fourier_stub( &
      ndgl, nloen, ngptot, ntrunc, nt, ugrid, u_fourier_real, u_fourier_imag, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      deallocate( &
        u_fourier_real, u_fourier_imag, v_fourier_real, v_fourier_imag, &
        basis_block_real, basis_block_imag, observed_block_real, observed_block_imag, &
        solution_real, solution_imag, vector_weights, &
        reduced_basis_real, reduced_basis_imag, reduced_solution_real, reduced_solution_imag, &
        active_columns)
      return
    end if
    call scalar_fourier_stub( &
      ndgl, nloen, ngptot, ntrunc, nt, vgrid, v_fourier_real, v_fourier_imag, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      deallocate( &
        u_fourier_real, u_fourier_imag, v_fourier_real, v_fourier_imag, &
        basis_block_real, basis_block_imag, observed_block_real, observed_block_imag, &
        solution_real, solution_imag, vector_weights, &
        reduced_basis_real, reduced_basis_imag, reduced_solution_real, reduced_solution_imag, &
        active_columns)
      return
    end if

    vector_weights(1:ndgl) = weights
    vector_weights(ndgl + 1_c_int:2_c_int * ndgl) = weights

    start_idx = 1_c_int
    do m = 0_c_int, ntrunc
      nblock = ntrunc - m + 1_c_int
      ncol = 2_c_int * nblock
      basis_block_real = 0.0_c_double
      basis_block_imag = 0.0_c_double
      observed_block_real = 0.0_c_double
      observed_block_imag = 0.0_c_double
      solution_real = 0.0_c_double
      solution_imag = 0.0_c_double

      do ilat = 1_c_int, ndgl
        base_idx = ((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt_basis
        u_row_offset = (ilat - 1_c_int) * ncol
        v_row_offset = (ndgl + ilat - 1_c_int) * ncol

        do icol = 1_c_int, nblock
          coeff_index = start_idx + icol - 1_c_int
          basis_block_real(u_row_offset + icol) = &
            vector_cache_u_basis_real(base_idx + coeff_index)
          basis_block_imag(u_row_offset + icol) = &
            vector_cache_u_basis_imag(base_idx + coeff_index)
          basis_block_real(u_row_offset + nblock + icol) = &
            vector_cache_u_basis_real(base_idx + ncoeff + coeff_index)
          basis_block_imag(u_row_offset + nblock + icol) = &
            vector_cache_u_basis_imag(base_idx + ncoeff + coeff_index)

          basis_block_real(v_row_offset + icol) = &
            vector_cache_v_basis_real(base_idx + coeff_index)
          basis_block_imag(v_row_offset + icol) = &
            vector_cache_v_basis_imag(base_idx + coeff_index)
          basis_block_real(v_row_offset + nblock + icol) = &
            vector_cache_v_basis_real(base_idx + ncoeff + coeff_index)
          basis_block_imag(v_row_offset + nblock + icol) = &
            vector_cache_v_basis_imag(base_idx + ncoeff + coeff_index)
        end do

        do it = 1_c_int, nt
          observed_block_real((ilat - 1_c_int) * nt + it) = &
            u_fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
          observed_block_imag((ilat - 1_c_int) * nt + it) = &
            u_fourier_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
          observed_block_real((ndgl + ilat - 1_c_int) * nt + it) = &
            v_fourier_real(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
          observed_block_imag((ndgl + ilat - 1_c_int) * nt + it) = &
            v_fourier_imag(((ilat - 1_c_int) * (ntrunc + 1_c_int) + m) * nt + it)
        end do
      end do

      active_count = 0_c_int
      do icol = 1_c_int, ncol
        column_norm = 0.0_c_double
        do ilat = 1_c_int, 2_c_int * ndgl
          base_idx = (ilat - 1_c_int) * ncol + icol
          column_norm = max( &
            column_norm, &
            abs(basis_block_real(base_idx)) + abs(basis_block_imag(base_idx)))
        end do
        if (column_norm > 0.0_c_double) then
          active_count = active_count + 1_c_int
          active_columns(active_count) = icol
        end if
      end do

      if (active_count > 0_c_int) then
        reduced_basis_real = 0.0_c_double
        reduced_basis_imag = 0.0_c_double
        reduced_solution_real = 0.0_c_double
        reduced_solution_imag = 0.0_c_double

        do ilat = 1_c_int, 2_c_int * ndgl
          do active_idx = 1_c_int, active_count
            icol = active_columns(active_idx)
            reduced_basis_real((ilat - 1_c_int) * active_count + active_idx) = &
              basis_block_real((ilat - 1_c_int) * ncol + icol)
            reduced_basis_imag((ilat - 1_c_int) * active_count + active_idx) = &
              basis_block_imag((ilat - 1_c_int) * ncol + icol)
          end do
        end do

        call solve_weighted_block_local( &
          2_c_int * ndgl, active_count, nt, vector_weights, &
          reduced_basis_real(1:2_c_int * ndgl * active_count), &
          reduced_basis_imag(1:2_c_int * ndgl * active_count), &
          observed_block_real, observed_block_imag, &
          reduced_solution_real(1:active_count * nt), &
          reduced_solution_imag(1:active_count * nt), ierr_local)
        if (ierr_local /= 0_c_int) then
          ierror = ierr_local
          deallocate( &
            u_fourier_real, u_fourier_imag, v_fourier_real, v_fourier_imag, &
            basis_block_real, basis_block_imag, observed_block_real, observed_block_imag, &
            solution_real, solution_imag, vector_weights, &
            reduced_basis_real, reduced_basis_imag, reduced_solution_real, reduced_solution_imag, &
            active_columns)
          return
        end if

        do active_idx = 1_c_int, active_count
          icol = active_columns(active_idx)
          do it = 1_c_int, nt
            solution_real((icol - 1_c_int) * nt + it) = &
              reduced_solution_real((active_idx - 1_c_int) * nt + it)
            solution_imag((icol - 1_c_int) * nt + it) = &
              reduced_solution_imag((active_idx - 1_c_int) * nt + it)
          end do
        end do
      end if

      do icol = 1_c_int, nblock
        coeff_index = start_idx + icol - 1_c_int
        do it = 1_c_int, nt
          vrtspec_r((coeff_index - 1_c_int) * nt + it) = &
            solution_real((icol - 1_c_int) * nt + it)
          vrtspec_i((coeff_index - 1_c_int) * nt + it) = &
            solution_imag((icol - 1_c_int) * nt + it)
          divspec_r((coeff_index - 1_c_int) * nt + it) = &
            solution_real((nblock + icol - 1_c_int) * nt + it)
          divspec_i((coeff_index - 1_c_int) * nt + it) = &
            solution_imag((nblock + icol - 1_c_int) * nt + it)
        end do
      end do

      start_idx = start_idx + nblock
    end do

    deallocate( &
      u_fourier_real, u_fourier_imag, v_fourier_real, v_fourier_imag, &
      basis_block_real, basis_block_imag, observed_block_real, observed_block_imag, &
      solution_real, solution_imag, vector_weights, &
      reduced_basis_real, reduced_basis_imag, reduced_solution_real, reduced_solution_imag, &
      active_columns)
    ierror = 0_c_int
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

    call vhgsi1_local(ndgl, imid, vb, wb, dtheta, dwts, dpbar, dwork)

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

  subroutine vhgsi1_local(nlat, imid, vb, wb, dthet, dwts, dpbar, work)
    integer(c_int), intent(in) :: nlat, imid
    real, intent(out) :: vb(imid, *), wb(imid, *)
    double precision, intent(inout) :: dthet(*), dwts(*), dpbar(imid, nlat, 3), work(*)

    integer(c_int) :: lwk, ierror, i, n, nm, nz, np, m, ix, iy
    real(c_double) :: abel, bbel, cbel, ssqr2, dcf

    lwk = nlat * (nlat + 2_c_int)
    call gaqd(nlat, dthet, dwts, dpbar, lwk, ierror)

    ssqr2 = 1.0_c_double / sqrt(2.0_c_double)
    do i = 1_c_int, imid
      dpbar(i, 1, 1) = ssqr2
    end do
    vb(:, 1) = 0.0
    wb(:, 1) = 0.0

    do n = 1_c_int, nlat - 1_c_int
      nm = mod(n - 2_c_int, 3_c_int) + 1_c_int
      nz = mod(n - 1_c_int, 3_c_int) + 1_c_int
      np = mod(n, 3_c_int) + 1_c_int

      call dnlfk_local(0_c_int, n, work)
      do i = 1_c_int, imid
        call dnlft_local(0_c_int, n, dthet(i), work, dpbar(i, 1, np))
      end do

      call dnlfk_local(1_c_int, n, work)
      do i = 1_c_int, imid
        call dnlft_local(1_c_int, n, dthet(i), work, dpbar(i, 2, np))
      end do

      if (n >= 2_c_int) then
        do m = 2_c_int, n
          abel = sqrt( &
            real((2_c_int * n + 1_c_int) * (m + n - 2_c_int) * (m + n - 3_c_int), c_double) / &
            real((2_c_int * n - 3_c_int) * (m + n - 1_c_int) * (m + n), c_double))
          bbel = sqrt( &
            real((2_c_int * n + 1_c_int) * (n - m - 1_c_int) * (n - m), c_double) / &
            real((2_c_int * n - 3_c_int) * (m + n - 1_c_int) * (m + n), c_double))
          cbel = sqrt( &
            real((n - m + 1_c_int) * (n - m + 2_c_int), c_double) / &
            real((m + n - 1_c_int) * (m + n), c_double))

          if (m < n - 1_c_int) then
            do i = 1_c_int, imid
              dpbar(i, m + 1_c_int, np) = &
                abel * dpbar(i, m - 1_c_int, nm) + &
                bbel * dpbar(i, m + 1_c_int, nm) - &
                cbel * dpbar(i, m - 1_c_int, np)
            end do
          else
            do i = 1_c_int, imid
              dpbar(i, m + 1_c_int, np) = &
                abel * dpbar(i, m - 1_c_int, nm) - &
                cbel * dpbar(i, m - 1_c_int, np)
            end do
          end if
        end do
      end if

      ix = vector_indx_local(nlat, 0_c_int, n)
      iy = vector_indx_local(nlat, n, n)
      do i = 1_c_int, imid
        vb(i, ix) = -real(dpbar(i, 2, np), kind(vb))
        vb(i, iy) = real(dpbar(i, n, np) / sqrt(real(2_c_int * (n + 1_c_int), c_double)), kind(vb))
      end do

      if (n >= 2_c_int) then
        dcf = sqrt(real(4_c_int * n * (n + 1_c_int), c_double))
        do m = 1_c_int, n - 1_c_int
          ix = vector_indx_local(nlat, m, n)
          abel = sqrt(real((n + m) * (n - m + 1_c_int), c_double)) / dcf
          bbel = sqrt(real((n - m) * (n + m + 1_c_int), c_double)) / dcf
          do i = 1_c_int, imid
            vb(i, ix) = real(abel * dpbar(i, m, np) - bbel * dpbar(i, m + 2_c_int, np), kind(vb))
          end do
        end do
      end if

      ix = vector_indx_local(nlat, 0_c_int, n)
      do i = 1_c_int, imid
        wb(i, ix) = 0.0
      end do

      dcf = sqrt(real(2_c_int * n + 1_c_int, c_double) / real(4_c_int * n * (n + 1_c_int) * (2_c_int * n - 1_c_int), c_double))
      do m = 1_c_int, n
        ix = vector_indx_local(nlat, m, n)
        abel = dcf * sqrt(real((n + m) * (n + m - 1_c_int), c_double))
        bbel = dcf * sqrt(real((n - m) * (n - m - 1_c_int), c_double))
        if (m < n - 1_c_int) then
          do i = 1_c_int, imid
            wb(i, ix) = real(abel * dpbar(i, m, nz) + bbel * dpbar(i, m + 2_c_int, nz), kind(wb))
          end do
        else
          do i = 1_c_int, imid
            wb(i, ix) = real(abel * dpbar(i, m, nz), kind(wb))
          end do
        end if
      end do
    end do
  end subroutine vhgsi1_local

  integer(c_int) function vector_indx_local(nlat, m, n) result(indx)
    integer(c_int), intent(in) :: nlat, m, n
    indx = m * nlat - (m * (m + 1_c_int)) / 2_c_int + n + 1_c_int
  end function vector_indx_local

  subroutine dnlfk_local(m, n, cp)
    integer(c_int), intent(in) :: m, n
    double precision, intent(out) :: cp(*)

    integer(c_int) :: ma, nmms2, l, nex, i
    real(c_double) :: fnum, fden, fnmh, a1, b1, c1, cp2, fnnp1, fnmsq, fk
    real(c_double) :: t1, t2, pm1
    real(c_double), parameter :: sc10 = 1024.0_c_double
    real(c_double), parameter :: sc20 = sc10 * sc10
    real(c_double), parameter :: sc40 = sc20 * sc20

    cp(1) = 0.0d0
    ma = abs(m)
    if (ma > n) return

    select case (n)
    case (: -1_c_int)
      return
    case (0_c_int)
      cp(1) = sqrt(2.0d0)
      return
    case (1_c_int)
      if (ma == 0_c_int) then
        cp(1) = sqrt(1.5d0)
      else
        cp(1) = sqrt(0.75d0)
        if (m == -1_c_int) cp(1) = -cp(1)
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
    nex = 20_c_int
    fden = 2.0_c_double
    if (nmms2 >= 1_c_int) then
      do i = 1_c_int, nmms2
        t1 = fnum * t1 / fden
        if (t1 > sc20) then
          t1 = t1 / sc40
          nex = nex + 40_c_int
        end if
        fnum = fnum + 2.0_c_double
        fden = fden + 2.0_c_double
      end do
    end if

    t1 = t1 / (2.0_c_double ** real(n - 1_c_int - nex, c_double))
    if (mod(ma / 2_c_int, 2_c_int) /= 0_c_int) t1 = -t1

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
    if (mod(n, 2_c_int) == 0_c_int .and. mod(ma, 2_c_int) == 0_c_int) l = l + 1_c_int
    cp(l) = cp2

    if (m < 0_c_int) then
      if (mod(ma, 2_c_int) /= 0_c_int) cp(l) = -cp(l)
    end if
    if (l <= 1_c_int) return

    fk = real(n, c_double)
    a1 = (fk - 2.0_c_double) * (fk - 1.0_c_double) - fnnp1
    b1 = 2.0_c_double * (fk * fk - fnmsq)
    cp(l - 1_c_int) = b1 * cp(l) / a1

    do l = l - 1_c_int, 2_c_int, -1_c_int
      fk = fk - 2.0_c_double
      a1 = (fk - 2.0_c_double) * (fk - 1.0_c_double) - fnnp1
      b1 = -2.0_c_double * (fk * fk - fnmsq)
      c1 = (fk + 1.0_c_double) * (fk + 2.0_c_double) - fnnp1
      cp(l - 1_c_int) = -(b1 * cp(l) + c1 * cp(l + 1_c_int)) / a1
    end do
  end subroutine dnlfk_local

  subroutine dnlft_local(m, n, theta, cp, pb)
    integer(c_int), intent(in) :: m, n
    double precision, intent(in) :: theta
    double precision, intent(in) :: cp(*)
    double precision, intent(out) :: pb

    integer(c_int) :: kdo, k
    real(c_double) :: cdt, sdt, cth, sth, chh

    cdt = cos(theta + theta)
    sdt = sin(theta + theta)

    if (mod(n, 2_c_int) == 0_c_int) then
      kdo = n / 2_c_int
      if (mod(m, 2_c_int) == 0_c_int) then
        pb = 0.5d0 * cp(1)
        if (n == 0_c_int) return
        cth = cdt
        sth = sdt
        do k = 1_c_int, kdo
          pb = pb + cp(k + 1_c_int) * cth
          chh = cdt * cth - sdt * sth
          sth = sdt * cth + cdt * sth
          cth = chh
        end do
      else
        pb = 0.0d0
        if (kdo == 0_c_int) return
        cth = cdt
        sth = sdt
        do k = 1_c_int, kdo
          pb = pb + cp(k) * sth
          chh = cdt * cth - sdt * sth
          sth = sdt * cth + cdt * sth
          cth = chh
        end do
      end if
    else
      kdo = (n + 1_c_int) / 2_c_int
      cth = cos(theta)
      sth = sin(theta)
      pb = 0.0d0
      if (mod(m, 2_c_int) == 0_c_int) then
        do k = 1_c_int, kdo
          pb = pb + cp(k) * cth
          chh = cdt * cth - sdt * sth
          sth = sdt * cth + cdt * sth
          cth = chh
        end do
      else
        do k = 1_c_int, kdo
          pb = pb + cp(k) * sth
          chh = cdt * cth - sdt * sth
          sth = sdt * cth + cdt * sth
          cth = chh
        end do
      end if
    end if
  end subroutine dnlft_local

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

  subroutine solve_weighted_block_local( &
    nrow, nblock, nt, weights, basis_real, basis_imag, observed_real, observed_imag, &
    solution_real, solution_imag, ierror)
    integer(c_int), intent(in) :: nrow, nblock, nt
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

    call solve_complex_system_local(gram, rhs, nblock, nt, info)
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
  end subroutine solve_weighted_block_local

  subroutine solve_complex_system_local(matrix, rhs, nsize, nrhs, info)
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
  end subroutine solve_complex_system_local

  logical function vector_basis_cache_matches(ndgl, nloen, ngptot, ntrunc, rsphere) result(matches)
    integer(c_int), intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), intent(in) :: ngptot, ntrunc
    real(c_double), intent(in) :: rsphere

    matches = .false.
    if (vector_basis_cache_ready == 0_c_int) then
      return
    end if
    if (.not. allocated(vector_cache_nloen)) then
      return
    end if
    if (.not. allocated(vector_cache_u_basis_real) .or. .not. allocated(vector_cache_u_basis_imag)) then
      return
    end if
    if (.not. allocated(vector_cache_v_basis_real) .or. .not. allocated(vector_cache_v_basis_imag)) then
      return
    end if
    if (vector_cache_ndgl /= ndgl .or. vector_cache_ngptot /= ngptot .or. vector_cache_ntrunc /= ntrunc) then
      return
    end if
    if (.not. all(vector_cache_nloen == nloen)) then
      return
    end if
    if (abs(vector_cache_rsphere - rsphere) > 1.0e-12_c_double) then
      return
    end if
    matches = .true.
  end function vector_basis_cache_matches

  subroutine reset_vector_basis_cache()
    if (allocated(vector_cache_nloen)) then
      deallocate(vector_cache_nloen)
    end if
    if (allocated(vector_cache_u_basis_real)) then
      deallocate(vector_cache_u_basis_real)
    end if
    if (allocated(vector_cache_u_basis_imag)) then
      deallocate(vector_cache_u_basis_imag)
    end if
    if (allocated(vector_cache_v_basis_real)) then
      deallocate(vector_cache_v_basis_real)
    end if
    if (allocated(vector_cache_v_basis_imag)) then
      deallocate(vector_cache_v_basis_imag)
    end if
    vector_basis_cache_ready = 0_c_int
    vector_cache_ndgl = 0_c_int
    vector_cache_ngptot = 0_c_int
    vector_cache_ntrunc = -1_c_int
    vector_cache_rsphere = 0.0_c_double
  end subroutine reset_vector_basis_cache

  subroutine ensure_vector_basis_cache(ndgl, nloen, ngptot, ntrunc, rsphere, ierror)
    integer(c_int), intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), intent(in) :: ngptot, ntrunc
    real(c_double), intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ncoeff, nt_basis, coeff_index, ierr_local
    real(c_double), allocatable :: vrtspec_r(:), vrtspec_i(:), divspec_r(:), divspec_i(:)
    real(c_double), allocatable :: basis_ugrid(:), basis_vgrid(:)

    if (vector_basis_cache_matches(ndgl, nloen, ngptot, ntrunc, rsphere)) then
      ierror = 0_c_int
      return
    end if

    call reset_vector_basis_cache()

    ncoeff = (ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int
    nt_basis = 2_c_int * ncoeff

    allocate(vector_cache_nloen(ndgl))
    vector_cache_nloen = nloen
    allocate(vrtspec_r(ncoeff * nt_basis), vrtspec_i(ncoeff * nt_basis))
    allocate(divspec_r(ncoeff * nt_basis), divspec_i(ncoeff * nt_basis))
    allocate(basis_ugrid(ngptot * nt_basis), basis_vgrid(ngptot * nt_basis))
    allocate(vector_cache_u_basis_real(ndgl * (ntrunc + 1_c_int) * nt_basis))
    allocate(vector_cache_u_basis_imag(ndgl * (ntrunc + 1_c_int) * nt_basis))
    allocate(vector_cache_v_basis_real(ndgl * (ntrunc + 1_c_int) * nt_basis))
    allocate(vector_cache_v_basis_imag(ndgl * (ntrunc + 1_c_int) * nt_basis))

    vrtspec_r = 0.0_c_double
    vrtspec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double
    do coeff_index = 1_c_int, ncoeff
      vrtspec_r((coeff_index - 1_c_int) * nt_basis + coeff_index) = 1.0_c_double
      divspec_r((coeff_index - 1_c_int) * nt_basis + ncoeff + coeff_index) = 1.0_c_double
    end do

    call uv_synthesis_stub( &
      ndgl, nloen, ngptot, rsphere, ntrunc, nt_basis, &
      vrtspec_r, vrtspec_i, divspec_r, divspec_i, &
      basis_ugrid, basis_vgrid, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      deallocate(vrtspec_r, vrtspec_i, divspec_r, divspec_i, basis_ugrid, basis_vgrid)
      call reset_vector_basis_cache()
      return
    end if

    call scalar_fourier_stub( &
      ndgl, nloen, ngptot, ntrunc, nt_basis, basis_ugrid, &
      vector_cache_u_basis_real, vector_cache_u_basis_imag, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      deallocate(vrtspec_r, vrtspec_i, divspec_r, divspec_i, basis_ugrid, basis_vgrid)
      call reset_vector_basis_cache()
      return
    end if
    call scalar_fourier_stub( &
      ndgl, nloen, ngptot, ntrunc, nt_basis, basis_vgrid, &
      vector_cache_v_basis_real, vector_cache_v_basis_imag, ierr_local)
    if (ierr_local /= 0_c_int) then
      ierror = ierr_local
      deallocate(vrtspec_r, vrtspec_i, divspec_r, divspec_i, basis_ugrid, basis_vgrid)
      call reset_vector_basis_cache()
      return
    end if

    deallocate(vrtspec_r, vrtspec_i, divspec_r, divspec_i, basis_ugrid, basis_vgrid)
    vector_basis_cache_ready = 1_c_int
    vector_cache_ndgl = ndgl
    vector_cache_ngptot = ngptot
    vector_cache_ntrunc = ntrunc
    vector_cache_rsphere = rsphere
    ierror = 0_c_int
  end subroutine ensure_vector_basis_cache

end module vector_backend
