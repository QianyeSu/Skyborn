!
! Reduced-Gaussian scalar transforms for the public spharm spectral layout.
!
! This file adds a minimal packed reduced-Gaussian scalar analysis/synthesis
! path directly under the spharm backend. The packed grid is ordered
! north-to-south by latitude circle, with all longitudes from one circle stored
! contiguously before the next circle.
!
! The public spectral layout matches the existing spharm complex triangular
! coefficient array:
!   dataspec(nm, it), nm = 1..(ntrunc+1)*(ntrunc+2)/2
!
module reduced_gaussian_scalar_mod
    implicit none
    private

    interface
        subroutine getlegfunc(legfunc, lat, ntrunc)
            integer, intent(in) :: ntrunc
            real, intent(in) :: lat
            real, intent(out) :: legfunc((ntrunc + 1) * (ntrunc + 2) / 2)
        end subroutine getlegfunc
    end interface

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)

    public :: reduced_gaussian_grdtospec_impl
    public :: reduced_gaussian_spectogrd_impl

    integer, save :: cache_ready = 0
    integer, save :: cache_ndgl = 0
    integer, save :: cache_ntrunc = -1
    real(sp), allocatable, save :: cache_lats(:)
    real(sp), allocatable, save :: cache_legfunc(:,:)
    real(sp), allocatable, save :: cache_public_scale(:)

contains

    pure integer function ncoeff_from_ntrunc(ntrunc) result(ncoeff)
        integer, intent(in) :: ntrunc
        ncoeff = (ntrunc + 1) * (ntrunc + 2) / 2
    end function ncoeff_from_ntrunc

    pure integer function infer_ntrunc_from_ncoeff(ncoeff) result(ntrunc)
        integer, intent(in) :: ncoeff
        real(dp) :: disc
        disc = sqrt(1.0_dp + 8.0_dp * real(ncoeff, dp))
        ntrunc = nint((-3.0_dp + disc) * 0.5_dp)
    end function infer_ntrunc_from_ncoeff

    pure integer function row_mmax(nlon_lat, ntrunc) result(mmax)
        integer, intent(in) :: nlon_lat, ntrunc
        mmax = min(ntrunc, (nlon_lat - 1) / 2)
    end function row_mmax

    logical function reduced_gaussian_cache_matches(ndgl, lats, ntrunc) result(matches)
        integer, intent(in) :: ndgl, ntrunc
        real(sp), intent(in) :: lats(ndgl)

        matches = .false.
        if (cache_ready == 0) return
        if (.not. allocated(cache_lats) .or. .not. allocated(cache_legfunc)) return
        if (cache_ndgl /= ndgl .or. cache_ntrunc /= ntrunc) return
        if (.not. all(abs(cache_lats - lats) <= 1.0e-5_sp)) return
        matches = .true.
    end function reduced_gaussian_cache_matches

    subroutine reset_reduced_gaussian_cache()
        if (allocated(cache_lats)) deallocate(cache_lats)
        if (allocated(cache_legfunc)) deallocate(cache_legfunc)
        if (allocated(cache_public_scale)) deallocate(cache_public_scale)
        cache_ready = 0
        cache_ndgl = 0
        cache_ntrunc = -1
    end subroutine reset_reduced_gaussian_cache

    subroutine ensure_reduced_gaussian_cache(ndgl, lats, ntrunc, ierror)
        integer, intent(in) :: ndgl, ntrunc
        real(sp), intent(in) :: lats(ndgl)
        integer, intent(out) :: ierror

        integer :: ilat, icoeff, ncoeff
        real(sp), allocatable :: public_legfunc(:)
        real(dp), allocatable :: scale_sum(:)
        integer, allocatable :: scale_count(:)
        real(dp) :: public_val

        if (reduced_gaussian_cache_matches(ndgl, lats, ntrunc)) then
            ierror = 0
            return
        end if

        call reset_reduced_gaussian_cache()

        ncoeff = ncoeff_from_ntrunc(ntrunc)
        allocate(cache_lats(ndgl))
        allocate(cache_legfunc(ncoeff, ndgl))
        allocate(cache_public_scale(ncoeff))
        allocate(public_legfunc(ncoeff))
        allocate(scale_sum(ncoeff))
        allocate(scale_count(ncoeff))

        cache_lats = lats
        scale_sum = 0.0_dp
        scale_count = 0
        do ilat = 1, ndgl
            call build_reduced_legendre_triangle(cache_legfunc(:, ilat), lats(ilat), ntrunc)
            call getlegfunc(public_legfunc, lats(ilat), ntrunc)
            do icoeff = 1, ncoeff
                public_val = real(public_legfunc(icoeff), dp)
                if (abs(public_val) > 1.0e-6_dp) then
                    scale_sum(icoeff) = scale_sum(icoeff) + &
                        real(cache_legfunc(icoeff, ilat), dp) / public_val
                    scale_count(icoeff) = scale_count(icoeff) + 1
                end if
            end do
        end do

        do icoeff = 1, ncoeff
            if (scale_count(icoeff) > 0) then
                cache_public_scale(icoeff) = real( &
                    scale_sum(icoeff) / real(scale_count(icoeff), dp), &
                    sp &
                )
            else
                cache_public_scale(icoeff) = 1.0_sp
            end if
        end do

        deallocate(public_legfunc, scale_sum, scale_count)
        cache_ready = 1
        cache_ndgl = ndgl
        cache_ntrunc = ntrunc
        ierror = 0
    end subroutine ensure_reduced_gaussian_cache

    subroutine build_reduced_legendre_triangle(legfunc, lat_deg, ntrunc)
        integer, intent(in) :: ntrunc
        real(sp), intent(in) :: lat_deg
        real(sp), intent(out) :: legfunc((ntrunc + 1) * (ntrunc + 2) / 2)

        integer :: waves, twowaves, jn, jk, jm, n, idx
        real(dp) :: mu, zcos2, lat, zan, zsqp, zcospar, zsinpar
        real(dp) :: zcosfak, zsinfak, jnmjk
        real(dp) :: zq, zwm2, zw, zwq, zq2m1, zwm2q2, z2q2
        real(dp) :: zcnm, zdnm, zenm
        real(dp), allocatable :: ztemp1(:), ztemp2(:), prev2(:), prev1(:), row(:)

        waves = ntrunc + 1
        twowaves = 2 * waves
        allocate(ztemp1(twowaves), ztemp2(twowaves), prev2(twowaves), prev1(twowaves), row(twowaves))

        mu = sin(real(lat_deg, dp) * acos(-1.0_dp) / 180.0_dp)
        zcos2 = sqrt(max(0.0_dp, 1.0_dp - mu * mu))
        lat = acos(max(-1.0_dp, min(1.0_dp, mu)))
        zan = 1.0_dp

        ztemp1 = 0.0_dp
        ztemp2 = 0.0_dp
        ztemp1(1) = 0.5_dp

        do jn = 2, twowaves
            zsqp = 1.0_dp / sqrt(real((jn - 1) + (jn - 1) * (jn - 1), dp))
            zan = zan * sqrt(1.0_dp - 1.0_dp / (4.0_dp * real((jn - 1) * (jn - 1), dp)))
            zcospar = cos(lat * real(jn - 1, dp))
            zsinpar = sin(lat * real(jn - 1, dp)) * real(jn - 1, dp) * zsqp
            zcosfak = 1.0_dp

            do jk = 2, jn - 2, 2
                jnmjk = real((jn - 1) - jk, dp)
                zcosfak = zcosfak * real(jk - 1, dp) * real((jn - 1) + (jn - 1 - jk) + 2, dp) / &
                    (real(jk, dp) * real((jn - 1) + (jn - 1 - jk) + 1, dp))
                zsinfak = zcosfak * zsqp * jnmjk
                zcospar = zcospar + zcosfak * cos(lat * jnmjk)
                zsinpar = zsinpar + zsinfak * sin(lat * jnmjk)
            end do

            if (mod(jn - 1, 2) == 0) then
                zcosfak = zcosfak * real(((jn - 1) - 1) * ((jn - 1) + 2), dp) / real((jn - 1) * jn, dp)
                zcospar = zcospar + zcosfak * 0.5_dp
            end if

            ztemp1(jn) = zan * zcospar
            ztemp2(jn - 1) = zan * zsinpar
        end do

        legfunc = 0.0_sp
        idx = 1
        do n = 0, ntrunc
            legfunc(idx) = real(ztemp1(n + 1), sp)
            idx = idx + 1
        end do

        prev2 = ztemp1
        prev1 = ztemp2
        idx = ntrunc + 2
        do n = 1, ntrunc
            legfunc(idx) = real(prev1(n), sp)
            idx = idx + 1
        end do

        do jm = 2, waves - 1
            row = 0.0_dp
            row(1) = sqrt(1.0_dp + 1.0_dp / real(jm + jm, dp)) * zcos2 * prev1(1)
            do jn = 2, twowaves - jm
                zq = real(jm + jm + jn - 3, dp)
                zwm2 = zq + real(jn - 1, dp)
                zw = zwm2 + 2.0_dp
                zwq = zw * zq
                zq2m1 = zq * zq - 1.0_dp
                zwm2q2 = zwm2 * zq2m1
                z2q2 = zq2m1 * 2.0_dp
                zcnm = sqrt((zwq * (zq - 2.0_dp)) / (zwm2q2 - z2q2))
                zdnm = sqrt((zwq * real(jn, dp)) / zwm2q2)
                zenm = sqrt(zw * real(jn - 1, dp) / ((zq + 1.0_dp) * zwm2))
                row(jn) = zcnm * prev2(jn) - mu * (zdnm * prev2(jn + 1) - zenm * row(jn - 1))
            end do

            do n = jm, ntrunc
                legfunc(idx) = real(row(n - jm + 1), sp)
                idx = idx + 1
            end do
            prev2 = prev1
            prev1 = row
        end do

        deallocate(ztemp1, ztemp2, prev2, prev1, row)
    end subroutine build_reduced_legendre_triangle

    subroutine prepare_fft_cache(ndgl, nloen, unique_nlon, row_fft_map, fft_offsets, fft_work, nuniq)
        integer, intent(in) :: ndgl
        integer, intent(in) :: nloen(ndgl)
        integer, intent(out) :: unique_nlon(ndgl)
        integer, intent(out) :: row_fft_map(ndgl)
        integer, intent(out) :: fft_offsets(ndgl)
        real(sp), allocatable, intent(out) :: fft_work(:)
        integer, intent(out) :: nuniq

        integer :: ilat, iu, total_len
        logical :: found

        nuniq = 0
        unique_nlon = 0
        row_fft_map = 0
        fft_offsets = 0

        do ilat = 1, ndgl
            found = .false.
            do iu = 1, nuniq
                if (unique_nlon(iu) == nloen(ilat)) then
                    row_fft_map(ilat) = iu
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                nuniq = nuniq + 1
                unique_nlon(nuniq) = nloen(ilat)
                row_fft_map(ilat) = nuniq
            end if
        end do

        total_len = 0
        do iu = 1, nuniq
            fft_offsets(iu) = total_len + 1
            total_len = total_len + unique_nlon(iu) + 15
        end do

        allocate(fft_work(total_len))
        do iu = 1, nuniq
            call hrffti( &
                unique_nlon(iu), &
                fft_work(fft_offsets(iu) : fft_offsets(iu) + unique_nlon(iu) + 14) &
            )
        end do
    end subroutine prepare_fft_cache

    subroutine solve_complex_system(matrix, rhs, nsize, nrhs, info)
        integer, intent(in) :: nsize, nrhs
        complex(dp), intent(inout) :: matrix(nsize, nsize)
        complex(dp), intent(inout) :: rhs(nsize, nrhs)
        integer, intent(out) :: info

        integer :: ipivot, jpivot, irhs, pivot_row
        complex(dp) :: factor, temp
        real(dp) :: pivot_abs

        info = 0

        do ipivot = 1, nsize
            pivot_row = ipivot
            pivot_abs = abs(matrix(ipivot, ipivot))
            do jpivot = ipivot + 1, nsize
                if (abs(matrix(jpivot, ipivot)) > pivot_abs) then
                    pivot_abs = abs(matrix(jpivot, ipivot))
                    pivot_row = jpivot
                end if
            end do

            if (pivot_abs <= tiny(1.0_dp)) then
                info = 10 + ipivot
                return
            end if

            if (pivot_row /= ipivot) then
                do jpivot = ipivot, nsize
                    temp = matrix(ipivot, jpivot)
                    matrix(ipivot, jpivot) = matrix(pivot_row, jpivot)
                    matrix(pivot_row, jpivot) = temp
                end do
                do irhs = 1, nrhs
                    temp = rhs(ipivot, irhs)
                    rhs(ipivot, irhs) = rhs(pivot_row, irhs)
                    rhs(pivot_row, irhs) = temp
                end do
            end if

            do jpivot = ipivot + 1, nsize
                factor = matrix(jpivot, ipivot) / matrix(ipivot, ipivot)
                matrix(jpivot, ipivot) = cmplx(0.0_dp, 0.0_dp, kind=dp)
                do irhs = ipivot + 1, nsize
                    matrix(jpivot, irhs) = matrix(jpivot, irhs) - factor * matrix(ipivot, irhs)
                end do
                do irhs = 1, nrhs
                    rhs(jpivot, irhs) = rhs(jpivot, irhs) - factor * rhs(ipivot, irhs)
                end do
            end do
        end do

        do ipivot = nsize, 1, -1
            if (abs(matrix(ipivot, ipivot)) <= tiny(1.0_dp)) then
                info = 100 + ipivot
                return
            end if
            do irhs = 1, nrhs
                temp = rhs(ipivot, irhs)
                do jpivot = ipivot + 1, nsize
                    temp = temp - matrix(ipivot, jpivot) * rhs(jpivot, irhs)
                end do
                rhs(ipivot, irhs) = temp / matrix(ipivot, ipivot)
            end do
        end do
    end subroutine solve_complex_system

    subroutine weighted_complex_block_solve(nrow, nblock, nt, weights, basis, observed_real, observed_imag, &
                                            solution_real, solution_imag, ierror)
        integer, intent(in) :: nrow, nblock, nt
        real(sp), intent(in) :: weights(nrow)
        real(sp), intent(in) :: basis(nrow, nblock)
        real(sp), intent(in) :: observed_real(nrow, nt)
        real(sp), intent(in) :: observed_imag(nrow, nt)
        real(sp), intent(out) :: solution_real(nblock, nt)
        real(sp), intent(out) :: solution_imag(nblock, nt)
        integer, intent(out) :: ierror

        integer :: irow, iblock, jblock, it, info
        real(dp) :: weight_dp, basis_i, basis_j
        complex(dp), allocatable :: gram(:,:), rhs(:,:)

        if (nrow < 1 .or. nblock < 1 .or. nt < 1) then
            solution_real = 0.0_sp
            solution_imag = 0.0_sp
            ierror = 1
            return
        end if

        allocate(gram(nblock, nblock))
        allocate(rhs(nblock, nt))
        gram = cmplx(0.0_dp, 0.0_dp, kind=dp)
        rhs = cmplx(0.0_dp, 0.0_dp, kind=dp)

        do irow = 1, nrow
            weight_dp = real(weights(irow), dp)
            do iblock = 1, nblock
                basis_i = real(basis(irow, iblock), dp)
                do jblock = 1, nblock
                    basis_j = real(basis(irow, jblock), dp)
                    gram(iblock, jblock) = gram(iblock, jblock) + weight_dp * basis_i * basis_j
                end do
                do it = 1, nt
                    rhs(iblock, it) = rhs(iblock, it) + weight_dp * basis_i * &
                        cmplx(real(observed_real(irow, it), dp), real(observed_imag(irow, it), dp), kind=dp)
                end do
            end do
        end do

        call solve_complex_system(gram, rhs, nblock, nt, info)
        if (info /= 0) then
            solution_real = 0.0_sp
            solution_imag = 0.0_sp
            ierror = info
            deallocate(gram, rhs)
            return
        end if

        do iblock = 1, nblock
            do it = 1, nt
                solution_real(iblock, it) = real(rhs(iblock, it), sp)
                solution_imag(iblock, it) = real(aimag(rhs(iblock, it)), sp)
            end do
        end do

        deallocate(gram, rhs)
        ierror = 0
    end subroutine weighted_complex_block_solve

    subroutine reduced_gaussian_grdtospec_impl(datagrid, nloen, weights, lats, ndgl, ngptot, ntrunc, nt, dataspec, ierror)
        integer, intent(in) :: ndgl, ngptot, ntrunc, nt
        real(sp), intent(in) :: datagrid(ngptot, nt)
        integer, intent(in) :: nloen(ndgl)
        real(sp), intent(in) :: weights(ndgl), lats(ndgl)
        complex, intent(out) :: dataspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
        integer, intent(out) :: ierror

        integer :: ilat, it, m, nblock, coeff_index, start_idx
        integer :: nactive, active_idx, ierr_local, nuniq, max_nblock
        integer :: nlon_lat, mmax_lat, iu
        integer, allocatable :: lat_offsets(:), active_rows(:)
        integer :: unique_nlon(ndgl), row_fft_map(ndgl), fft_offsets(ndgl)
        real(sp), allocatable :: fft_work(:)
        real(sp), allocatable :: row_buffer(:,:), row_work(:,:)
        real(sp), allocatable :: fourier_real(:,:,:), fourier_imag(:,:,:)
        real(sp), allocatable :: basis_block(:,:), observed_block_real(:,:), observed_block_imag(:,:)
        real(sp), allocatable :: solution_real(:,:), solution_imag(:,:), weights_active(:)

        if (ndgl < 3 .or. ngptot < 1 .or. nt < 1 .or. ntrunc < 0) then
            dataspec = cmplx(0.0_sp, 0.0_sp)
            ierror = 1
            return
        end if
        if (ntrunc > ndgl - 1) then
            dataspec = cmplx(0.0_sp, 0.0_sp)
            ierror = 2
            return
        end if
        if (sum(nloen) /= ngptot) then
            dataspec = cmplx(0.0_sp, 0.0_sp)
            ierror = 3
            return
        end if
        if (any(nloen <= 0)) then
            dataspec = cmplx(0.0_sp, 0.0_sp)
            ierror = 4
            return
        end if

        call ensure_reduced_gaussian_cache(ndgl, lats, ntrunc, ierr_local)
        if (ierr_local /= 0) then
            dataspec = cmplx(0.0_sp, 0.0_sp)
            ierror = 5
            return
        end if

        allocate(lat_offsets(ndgl))
        lat_offsets(1) = 0
        do ilat = 2, ndgl
            lat_offsets(ilat) = lat_offsets(ilat - 1) + nloen(ilat - 1)
        end do

        call prepare_fft_cache(ndgl, nloen, unique_nlon, row_fft_map, fft_offsets, fft_work, nuniq)

        allocate(fourier_real(ndgl, ntrunc + 1, nt))
        allocate(fourier_imag(ndgl, ntrunc + 1, nt))
        fourier_real = 0.0_sp
        fourier_imag = 0.0_sp

        allocate(row_buffer(1, maxval(nloen)))
        allocate(row_work(1, maxval(nloen)))

        do ilat = 1, ndgl
            nlon_lat = nloen(ilat)
            iu = row_fft_map(ilat)
            mmax_lat = row_mmax(nlon_lat, ntrunc)
            do it = 1, nt
                row_buffer = 0.0_sp
                row_work = 0.0_sp
                row_buffer(1, 1:nlon_lat) = datagrid(lat_offsets(ilat) + 1 : lat_offsets(ilat) + nlon_lat, it)
                call hrfftf(1, nlon_lat, row_buffer, 1, &
                            fft_work(fft_offsets(iu) : fft_offsets(iu) + nlon_lat + 14), row_work)
                row_buffer(1, 1:nlon_lat) = row_buffer(1, 1:nlon_lat) * (2.0_sp / real(nlon_lat, sp))
                fourier_real(ilat, 1, it) = row_buffer(1, 1)
                do m = 1, mmax_lat
                    fourier_real(ilat, m + 1, it) = row_buffer(1, 2 * m)
                    fourier_imag(ilat, m + 1, it) = row_buffer(1, 2 * m + 1)
                end do
            end do
        end do

        max_nblock = ntrunc + 1
        allocate(active_rows(ndgl))
        allocate(basis_block(ndgl, max_nblock))
        allocate(observed_block_real(ndgl, nt))
        allocate(observed_block_imag(ndgl, nt))
        allocate(solution_real(max_nblock, nt))
        allocate(solution_imag(max_nblock, nt))
        allocate(weights_active(ndgl))

        dataspec = cmplx(0.0_sp, 0.0_sp)
        start_idx = 1
        do m = 0, ntrunc
            nblock = ntrunc - m + 1
            nactive = 0
            do ilat = 1, ndgl
                if (row_mmax(nloen(ilat), ntrunc) >= m) then
                    nactive = nactive + 1
                    active_rows(nactive) = ilat
                end if
            end do
            if (nactive < nblock) then
                dataspec = cmplx(0.0_sp, 0.0_sp)
                ierror = 20 + m
                deallocate(lat_offsets, fft_work, fourier_real, fourier_imag, row_buffer, row_work, &
                           active_rows, basis_block, observed_block_real, observed_block_imag, &
                           solution_real, solution_imag, weights_active)
                return
            end if

            do active_idx = 1, nactive
                ilat = active_rows(active_idx)
                weights_active(active_idx) = weights(ilat)
                do coeff_index = 1, nblock
                    basis_block(active_idx, coeff_index) = cache_legfunc(start_idx + coeff_index - 1, ilat)
                end do
                do it = 1, nt
                    observed_block_real(active_idx, it) = fourier_real(ilat, m + 1, it)
                    observed_block_imag(active_idx, it) = fourier_imag(ilat, m + 1, it)
                end do
            end do

            call weighted_complex_block_solve(nactive, nblock, nt, weights_active(1:nactive), &
                                              basis_block(1:nactive, 1:nblock), &
                                              observed_block_real(1:nactive, 1:nt), &
                                              observed_block_imag(1:nactive, 1:nt), &
                                              solution_real(1:nblock, 1:nt), solution_imag(1:nblock, 1:nt), ierr_local)
            if (ierr_local /= 0) then
                dataspec = cmplx(0.0_sp, 0.0_sp)
                ierror = 100 + ierr_local
                deallocate(lat_offsets, fft_work, fourier_real, fourier_imag, row_buffer, row_work, &
                           active_rows, basis_block, observed_block_real, observed_block_imag, &
                           solution_real, solution_imag, weights_active)
                return
            end if

            do coeff_index = 1, nblock
                do it = 1, nt
                    dataspec(start_idx + coeff_index - 1, it) = &
                        (0.5_sp * cache_public_scale(start_idx + coeff_index - 1)) * cmplx( &
                            solution_real(coeff_index, it), solution_imag(coeff_index, it))
                end do
            end do
            start_idx = start_idx + nblock
        end do

        deallocate(lat_offsets, fft_work, fourier_real, fourier_imag, row_buffer, row_work, &
                   active_rows, basis_block, observed_block_real, observed_block_imag, &
                   solution_real, solution_imag, weights_active)
        ierror = 0
    end subroutine reduced_gaussian_grdtospec_impl

    subroutine reduced_gaussian_spectogrd_impl(dataspec, nloen, lats, ndgl, ngptot, nmdim, nt, datagrid, ierror)
        integer, intent(in) :: ndgl, ngptot, nmdim, nt
        complex, intent(in) :: dataspec(nmdim, nt)
        integer, intent(in) :: nloen(ndgl)
        real(sp), intent(in) :: lats(ndgl)
        real(sp), intent(out) :: datagrid(ngptot, nt)
        integer, intent(out) :: ierror

        integer :: ilat, it, m, ntrunc, nblock, coeff_index, start_idx
        integer :: nlon_lat, mmax_lat, iu, nuniq
        integer, allocatable :: lat_offsets(:)
        integer :: unique_nlon(ndgl), row_fft_map(ndgl), fft_offsets(ndgl)
        real(sp), allocatable :: fft_work(:)
        real(sp), allocatable :: row_buffer(:,:), row_work(:,:)
        real(sp), allocatable :: a(:,:,:), b(:,:,:)
        real(sp), allocatable :: fourier_real(:,:,:), fourier_imag(:,:,:)
        complex, allocatable :: local_spec(:,:)
        integer :: expected_ncoeff, ierr_local

        if (ndgl < 3 .or. ngptot < 1 .or. nt < 1 .or. nmdim < 1) then
            datagrid = 0.0_sp
            ierror = 1
            return
        end if
        if (sum(nloen) /= ngptot) then
            datagrid = 0.0_sp
            ierror = 2
            return
        end if

        ntrunc = infer_ntrunc_from_ncoeff(nmdim)
        expected_ncoeff = ncoeff_from_ntrunc(ntrunc)
        if (expected_ncoeff /= nmdim) then
            datagrid = 0.0_sp
            ierror = 3
            return
        end if
        if (ntrunc > ndgl - 1) then
            datagrid = 0.0_sp
            ierror = 4
            return
        end if

        call ensure_reduced_gaussian_cache(ndgl, lats, ntrunc, ierr_local)
        if (ierr_local /= 0) then
            datagrid = 0.0_sp
            ierror = 5
            return
        end if

        allocate(lat_offsets(ndgl))
        lat_offsets(1) = 0
        do ilat = 2, ndgl
            lat_offsets(ilat) = lat_offsets(ilat - 1) + nloen(ilat - 1)
        end do

        call prepare_fft_cache(ndgl, nloen, unique_nlon, row_fft_map, fft_offsets, fft_work, nuniq)

        allocate(local_spec(nmdim, nt))
        do coeff_index = 1, nmdim
            local_spec(coeff_index, :) = dataspec(coeff_index, :) / cache_public_scale(coeff_index)
        end do

        allocate(a(ndgl, ndgl, nt))
        allocate(b(ndgl, ndgl, nt))
        a = 0.0_sp
        b = 0.0_sp
        call onedtotwod(local_spec, a, b, ndgl, nmdim, nt)

        allocate(fourier_real(ndgl, ntrunc + 1, nt))
        allocate(fourier_imag(ndgl, ntrunc + 1, nt))
        fourier_real = 0.0_sp
        fourier_imag = 0.0_sp

        start_idx = 1
        do m = 0, ntrunc
            nblock = ntrunc - m + 1
            do ilat = 1, ndgl
                do coeff_index = 1, nblock
                    do it = 1, nt
                        fourier_real(ilat, m + 1, it) = fourier_real(ilat, m + 1, it) + &
                            cache_legfunc(start_idx + coeff_index - 1, ilat) * a(m + 1, m + coeff_index, it)
                        fourier_imag(ilat, m + 1, it) = fourier_imag(ilat, m + 1, it) + &
                            cache_legfunc(start_idx + coeff_index - 1, ilat) * b(m + 1, m + coeff_index, it)
                    end do
                end do
            end do
            start_idx = start_idx + nblock
        end do

        allocate(row_buffer(1, maxval(nloen)))
        allocate(row_work(1, maxval(nloen)))
        datagrid = 0.0_sp

        do ilat = 1, ndgl
            nlon_lat = nloen(ilat)
            iu = row_fft_map(ilat)
            mmax_lat = row_mmax(nlon_lat, ntrunc)
            do it = 1, nt
                row_buffer = 0.0_sp
                row_work = 0.0_sp
                row_buffer(1, 1) = fourier_real(ilat, 1, it)
                do m = 1, mmax_lat
                    row_buffer(1, 2 * m) = fourier_real(ilat, m + 1, it)
                    row_buffer(1, 2 * m + 1) = fourier_imag(ilat, m + 1, it)
                end do
                call hrfftb(1, nlon_lat, row_buffer, 1, &
                            fft_work(fft_offsets(iu) : fft_offsets(iu) + nlon_lat + 14), row_work)
                datagrid(lat_offsets(ilat) + 1 : lat_offsets(ilat) + nlon_lat, it) = 0.5_sp * row_buffer(1, 1:nlon_lat)
            end do
        end do

        deallocate(lat_offsets, fft_work, local_spec, a, b, fourier_real, fourier_imag, row_buffer, row_work)
        ierror = 0
    end subroutine reduced_gaussian_spectogrd_impl

end module reduced_gaussian_scalar_mod

subroutine reduced_gaussian_grdtospec(datagrid, nloen, weights, lats, dataspec, ndgl, ngptot, ntrunc, nt, ierror)
    use reduced_gaussian_scalar_mod, only : reduced_gaussian_grdtospec_impl
    implicit none

    integer, intent(in) :: ndgl, ngptot, ntrunc, nt
    real, intent(in) :: datagrid(ngptot, nt)
    integer, intent(in) :: nloen(ndgl)
    real, intent(in) :: weights(ndgl), lats(ndgl)
    complex, intent(out) :: dataspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    call reduced_gaussian_grdtospec_impl(datagrid, nloen, weights, lats, ndgl, ngptot, ntrunc, nt, dataspec, ierror)
end subroutine reduced_gaussian_grdtospec

subroutine reduced_gaussian_spectogrd(dataspec, nloen, lats, datagrid, ndgl, ngptot, nmdim, nt, ierror)
    use reduced_gaussian_scalar_mod, only : reduced_gaussian_spectogrd_impl
    implicit none

    integer, intent(in) :: ndgl, ngptot, nmdim, nt
    complex, intent(in) :: dataspec(nmdim, nt)
    integer, intent(in) :: nloen(ndgl)
    real, intent(in) :: lats(ndgl)
    real, intent(out) :: datagrid(ngptot, nt)
    integer, intent(out) :: ierror

    call reduced_gaussian_spectogrd_impl(dataspec, nloen, lats, ndgl, ngptot, nmdim, nt, datagrid, ierror)
end subroutine reduced_gaussian_spectogrd
