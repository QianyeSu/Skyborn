!
! Direct scalar transforms for packed reduced Gaussian grids.
!
! These helpers use the same public spectral coefficient convention as the
! existing Spharmt Gaussian backend, but replace the rectangular longitude
! FFT with row-wise Fourier sums over each reduced latitude circle.
!

module reduced_gaussian_fft_cache_mod
    implicit none

    integer, save :: reduced_gaussian_fft_cache_nlat = 0
    integer, save :: reduced_gaussian_fft_cache_max_wsave = 0
    integer, allocatable, save :: reduced_gaussian_fft_cache_pl(:)
    real, allocatable, save :: reduced_gaussian_fft_wsave(:, :)

contains

    subroutine reduced_gaussian_prepare_fft_cache(pl, nlat, ierror)
        implicit none

        integer, intent(in) :: nlat
        integer, intent(in) :: pl(nlat)
        integer, intent(out) :: ierror

        integer :: j, max_nlon, max_wsave, alloc_status
        logical :: cache_matches

        external :: hrffti

        ierror = 1
        if (nlat < 1) return

        max_nlon = 0
        do j = 1, nlat
            if (pl(j) < 1) return
            if (pl(j) > max_nlon) max_nlon = pl(j)
        end do
        max_wsave = 2 * max_nlon + 15

        cache_matches = allocated(reduced_gaussian_fft_cache_pl) .and. &
            allocated(reduced_gaussian_fft_wsave) .and. &
            reduced_gaussian_fft_cache_nlat == nlat .and. &
            reduced_gaussian_fft_cache_max_wsave >= max_wsave
        if (cache_matches) then
            if (all(reduced_gaussian_fft_cache_pl(1:nlat) == pl(1:nlat))) then
                ierror = 0
                return
            end if
        end if

!$omp critical
        cache_matches = allocated(reduced_gaussian_fft_cache_pl) .and. &
            allocated(reduced_gaussian_fft_wsave) .and. &
            reduced_gaussian_fft_cache_nlat == nlat .and. &
            reduced_gaussian_fft_cache_max_wsave >= max_wsave
        if (cache_matches) then
            if (all(reduced_gaussian_fft_cache_pl(1:nlat) == pl(1:nlat))) then
                ierror = 0
            else
                cache_matches = .false.
            end if
        end if

        if (.not. cache_matches) then
            if (allocated(reduced_gaussian_fft_cache_pl)) then
                deallocate(reduced_gaussian_fft_cache_pl)
            end if
            if (allocated(reduced_gaussian_fft_wsave)) then
                deallocate(reduced_gaussian_fft_wsave)
            end if

            allocate(reduced_gaussian_fft_cache_pl(nlat), &
                     reduced_gaussian_fft_wsave(max_wsave, nlat), &
                     stat=alloc_status)
            if (alloc_status /= 0) then
                reduced_gaussian_fft_cache_nlat = 0
                reduced_gaussian_fft_cache_max_wsave = 0
                ierror = 2
            else
                reduced_gaussian_fft_cache_pl(:) = pl(:)
                reduced_gaussian_fft_wsave(:, :) = 0.0
                do j = 1, nlat
                    call hrffti(pl(j), reduced_gaussian_fft_wsave(1, j))
                end do
                reduced_gaussian_fft_cache_nlat = nlat
                reduced_gaussian_fft_cache_max_wsave = max_wsave
                ierror = 0
            end if
        end if
!$omp end critical

    end subroutine reduced_gaussian_prepare_fft_cache

end module reduced_gaussian_fft_cache_mod


subroutine reduced_gaussian_legendre_basis(basis, nlat, ntrunc, ierror)
    implicit none

    integer, intent(in) :: nlat, ntrunc
    real, intent(out) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    integer, intent(out) :: ierror

    integer :: i, js, m, n, nm, nmdim, l, late, nlon_basis
    integer :: lwork, ldwork, lshags, alloc_status, mode, km
    real :: parity, value
    real, allocatable :: wshags(:), work(:), pmn(:, :, :)
    double precision, allocatable :: dwork(:)

    external :: shagsp, legin

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2
    basis(:, :) = 0.0

    l = ntrunc + 1
    late = (nlat + 1) / 2
    nlon_basis = max(4, 2 * ntrunc)
    lshags = nlat * (2 * late + 3 * l - 2) + &
        3 * l * (1 - l) / 2 + nlon_basis + 15
    lwork = 4 * nlat * (nlat + 2) + 2
    ldwork = nlat * (nlat + 4)

    allocate(wshags(lshags), work(lwork), dwork(ldwork), &
             pmn(nlat, late, 3), stat=alloc_status)
    ierror = 3
    if (alloc_status /= 0) return

    call shagsp(nlat, nlon_basis, wshags, lshags, dwork, ldwork, ierror)
    if (ierror /= 0) then
        ierror = 4
        deallocate(wshags, work, dwork, pmn)
        return
    end if

    pmn(:, :, :) = 0.0
    mode = 0
    do m = 0, ntrunc
        call legin(mode, l, nlat, m, wshags, pmn, km)
        parity = 1.0
        do n = m, ntrunc
            nm = m * (2 * ntrunc - m + 3) / 2 + n - m + 1
            do i = 1, late
                value = pmn(n + 1, i, km)
                basis(i, nm) = value
                js = nlat - i + 1
                if (js /= i) basis(js, nm) = parity * value
            end do
            parity = -parity
        end do
    end do

    deallocate(wshags, work, dwork, pmn)
    ierror = 0

end subroutine reduced_gaussian_legendre_basis


subroutine reduced_gaussian_legendre_derivative_from_basis(basis, dbasis, &
                                                           nlat, ntrunc, &
                                                           ierror)
    implicit none

    integer, intent(in) :: nlat, ntrunc
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(out) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    integer, intent(out) :: ierror

    integer :: j, m, n, nm, prev_nm, ldwork, lw
    double precision :: theta(nlat), weights(nlat)
    double precision :: dwork(nlat * (nlat + 2))
    double precision :: sin_theta, cos_theta, coeff, deriv

    external :: gaqd

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    dbasis(:, :) = 0.0

    ierror = 0
    ldwork = nlat * (nlat + 2)
    lw = ldwork
    call gaqd(nlat, theta, weights, dwork, lw, ierror)
    if (ierror /= 0) return

!$omp parallel do schedule(dynamic) private(m, n, nm, prev_nm, j, &
!$omp& sin_theta, cos_theta, coeff, deriv)
    do m = 0, ntrunc
        do n = max(1, m), ntrunc
            nm = m * (2 * ntrunc - m + 3) / 2 + n - m + 1
            coeff = 0.0d0
            prev_nm = 0
            if (n > m) then
                prev_nm = m * (2 * ntrunc - m + 3) / 2 + (n - 1) - m + 1
                coeff = sqrt(dble((2 * n + 1) * (n * n - m * m)) / &
                             dble(2 * n - 1))
            end if

            do j = 1, nlat
                sin_theta = sin(theta(j))
                cos_theta = cos(theta(j))
                deriv = dble(n) * cos_theta * dble(basis(j, nm))
                if (prev_nm > 0) then
                    deriv = deriv - coeff * dble(basis(j, prev_nm))
                end if
                dbasis(j, nm) = real(deriv / sin_theta)
            end do
        end do
    end do
!$omp end parallel do

end subroutine reduced_gaussian_legendre_derivative_from_basis


subroutine reduced_gaussian_legendre_derivative_basis(dbasis, nlat, ntrunc, ierror)
    implicit none

    integer, intent(in) :: nlat, ntrunc
    real, intent(out) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    integer, intent(out) :: ierror

    integer :: nmdim, alloc_status
    real, allocatable :: basis(:, :)

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2
    dbasis(:, :) = 0.0

    allocate(basis(nlat, nmdim), stat=alloc_status)
    ierror = 3
    if (alloc_status /= 0) return

    call reduced_gaussian_legendre_basis(basis, nlat, ntrunc, ierror)
    if (ierror /= 0) then
        deallocate(basis)
        return
    end if

    call reduced_gaussian_legendre_derivative_from_basis( &
        basis, dbasis, nlat, ntrunc, ierror)
    deallocate(basis)

end subroutine reduced_gaussian_legendre_derivative_basis


subroutine reduced_gaussian_grdtospec(datagrid, pl, weights, basis, dataspec, &
                                      ngptot, nlat, ntrunc, nt, ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: datagrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    complex, intent(out) :: dataspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total, twom
    integer :: max_nlon, mlimit, alloc_status
    integer :: offsets(nlat + 1)
    double precision :: inv_nlon
    double precision :: coeff_real, coeff_imag, weighted_basis
    double precision, allocatable :: cos_fft(:, :), sin_fft(:, :)
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftf

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nt < 1) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    ierror = 4
    if (total /= ngptot) return

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 5
    if (alloc_status /= 0) return

    allocate(cos_fft(0:ntrunc, nlat), sin_fft(0:ntrunc, nlat), &
             stat=alloc_status)
    ierror = 6
    if (alloc_status /= 0) return

    dataspec(:, :) = (0.0, 0.0)

    do k = 1, nt
!$omp parallel private(j, i, m, offset, row_nlon, mlimit, inv_nlon, twom, &
!$omp& fft_rows, fft_work)
        allocate(fft_rows(1, max_nlon), fft_work(1, max_nlon))
!$omp do schedule(dynamic)
        do j = 1, nlat
            offset = offsets(j)
            row_nlon = pl(j)
            mlimit = min(ntrunc, row_nlon / 2)
            inv_nlon = 1.0d0 / dble(row_nlon)
            fft_rows(:, 1:row_nlon) = 0.0
            do i = 1, row_nlon
                fft_rows(1, i) = datagrid(offset + i, k)
            end do

            call hrfftf(1, row_nlon, fft_rows, 1, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            cos_fft(0, j) = dble(fft_rows(1, 1)) * inv_nlon
            sin_fft(0, j) = 0.0d0
            do m = 1, mlimit
                twom = 2 * m
                cos_fft(m, j) = dble(fft_rows(1, twom)) * inv_nlon
                if (twom < row_nlon) then
                    sin_fft(m, j) = dble(fft_rows(1, twom + 1)) * inv_nlon
                else
                    sin_fft(m, j) = 0.0d0
                end if
            end do
        end do
!$omp end do
        deallocate(fft_rows, fft_work)
!$omp end parallel

!$omp parallel do schedule(dynamic) private(m, j, n, nm, row_nlon, &
!$omp& coeff_real, coeff_imag, weighted_basis)
        do m = 0, ntrunc
            do j = 1, nlat
                row_nlon = pl(j)
                if (2 * m <= row_nlon) then
                    coeff_real = cos_fft(m, j)
                    coeff_imag = sin_fft(m, j)
                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        weighted_basis = dble(weights(j)) * dble(basis(j, nm))
                        dataspec(nm, k) = dataspec(nm, k) + &
                            cmplx(real(coeff_real * weighted_basis), &
                                  real(coeff_imag * weighted_basis))
                    end do
                end if
            end do
        end do
!$omp end parallel do
    end do

    deallocate(cos_fft, sin_fft)
    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_grdtospec


subroutine reduced_gaussian_spectogrd(dataspec, pl, basis, datagrid, &
                                      nmdim, nlat, ntrunc, nt, ngptot, ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    complex, intent(in) :: dataspec(nmdim, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(out) :: datagrid(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, js, i, k, m, n, nm, offset, offset_s, row_nlon, total
    integer :: max_nlon, mlimit, alloc_status, pair_count, half_nlat
    integer :: offsets(nlat + 1)
    logical :: symmetric_pl
    double precision :: scrm_real, scrm_imag, scrm_real_s, scrm_imag_s
    double precision :: term_real, term_imag, pval, parity
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftb

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    ierror = 5
    if (total /= ngptot) return

    pair_count = nlat / 2
    half_nlat = (nlat + 1) / 2
    symmetric_pl = .true.
    do j = 1, pair_count
        if (pl(j) /= pl(nlat - j + 1)) symmetric_pl = .false.
    end do

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 6
    if (alloc_status /= 0) return

    datagrid(:, :) = 0.0

!$omp parallel private(j, js, i, k, m, n, nm, offset, offset_s, &
!$omp& row_nlon, mlimit, scrm_real, scrm_imag, scrm_real_s, &
!$omp& scrm_imag_s, term_real, term_imag, pval, parity, &
!$omp& fft_rows, fft_work)
    allocate(fft_rows(2 * nt, max_nlon), fft_work(2 * nt, max_nlon))
    if (symmetric_pl) then
!$omp do schedule(dynamic)
    do j = 1, half_nlat
        offset = offsets(j)
        if (j <= pair_count) then
            js = nlat - j + 1
            offset_s = offsets(js)
        else
            js = j
            offset_s = offset
        end if
        row_nlon = pl(j)
        mlimit = min(ntrunc, row_nlon / 2)
        fft_rows(:, 1:row_nlon) = 0.0

        do k = 1, nt
            do m = 0, mlimit
                scrm_real = 0.0d0
                scrm_imag = 0.0d0
                scrm_real_s = 0.0d0
                scrm_imag_s = 0.0d0
                nm = spectral_offset(m, ntrunc)
                parity = 1.0d0
                do n = m, ntrunc
                    nm = nm + 1
                    pval = dble(basis(j, nm))
                    term_real = dble(real(dataspec(nm, k))) * pval
                    term_imag = dble(aimag(dataspec(nm, k))) * pval
                    scrm_real = scrm_real + term_real
                    scrm_imag = scrm_imag + term_imag
                    scrm_real_s = scrm_real_s + parity * term_real
                    scrm_imag_s = scrm_imag_s + parity * term_imag
                    parity = -parity
                end do

                if (m == 0) then
                    fft_rows(k, 1) = real(2.0d0 * scrm_real)
                    if (j <= pair_count) then
                        fft_rows(nt + k, 1) = real(2.0d0 * scrm_real_s)
                    end if
                else
                    fft_rows(k, 2 * m) = real(2.0d0 * scrm_real)
                    if (j <= pair_count) then
                        fft_rows(nt + k, 2 * m) = real(2.0d0 * scrm_real_s)
                    end if
                    if (2 * m < row_nlon) then
                        fft_rows(k, 2 * m + 1) = real(2.0d0 * scrm_imag)
                        if (j <= pair_count) then
                            fft_rows(nt + k, 2 * m + 1) = &
                                real(2.0d0 * scrm_imag_s)
                        end if
                    end if
                end if
            end do
        end do

        if (j <= pair_count) then
            call hrfftb(2 * nt, row_nlon, fft_rows, 2 * nt, &
                reduced_gaussian_fft_wsave(1, j), fft_work)
        else
            call hrfftb(nt, row_nlon, fft_rows, 2 * nt, &
                reduced_gaussian_fft_wsave(1, j), fft_work)
        end if

        do k = 1, nt
            do i = 1, row_nlon
                datagrid(offset + i, k) = 0.5 * fft_rows(k, i)
                if (j <= pair_count) then
                    datagrid(offset_s + i, k) = 0.5 * fft_rows(nt + k, i)
                end if
            end do
        end do
    end do
!$omp end do
    else
!$omp do schedule(dynamic)
    do j = 1, nlat
        offset = offsets(j)
        row_nlon = pl(j)
        mlimit = min(ntrunc, row_nlon / 2)
        fft_rows(:, 1:row_nlon) = 0.0

        do k = 1, nt
            do m = 0, mlimit
                scrm_real = 0.0d0
                scrm_imag = 0.0d0
                nm = spectral_offset(m, ntrunc)
                do n = m, ntrunc
                    nm = nm + 1
                    pval = dble(basis(j, nm))
                    scrm_real = scrm_real + dble(real(dataspec(nm, k))) * pval
                    scrm_imag = scrm_imag + dble(aimag(dataspec(nm, k))) * pval
                end do

                if (m == 0) then
                    fft_rows(k, 1) = real(2.0d0 * scrm_real)
                else
                    fft_rows(k, 2 * m) = real(2.0d0 * scrm_real)
                    if (2 * m < row_nlon) then
                        fft_rows(k, 2 * m + 1) = real(2.0d0 * scrm_imag)
                    end if
                end if
            end do
        end do

        call hrfftb(nt, row_nlon, fft_rows, 2 * nt, &
            reduced_gaussian_fft_wsave(1, j), fft_work)

        do k = 1, nt
            do i = 1, row_nlon
                datagrid(offset + i, k) = 0.5 * fft_rows(k, i)
            end do
        end do
    end do
!$omp end do
    end if
    deallocate(fft_rows, fft_work)
!$omp end parallel

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_spectogrd


subroutine reduced_gaussian_spectogrd_pair(speca, specb, pl, basis, &
                                           grida, gridb, nmdim, nlat, &
                                           ntrunc, nt, ngptot, ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    complex, intent(in) :: speca(nmdim, nt)
    complex, intent(in) :: specb(nmdim, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(out) :: grida(ngptot, nt)
    real, intent(out) :: gridb(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, js, i, k, m, n, nm, offset, offset_s, row_nlon, total
    integer :: max_nlon, mlimit, alloc_status, pair_count, half_nlat
    integer :: offsets(nlat + 1)
    integer :: row_a, row_b, row_a_s, row_b_s
    logical :: symmetric_pl
    double precision :: a_real, a_imag, b_real, b_imag
    double precision :: a_real_s, a_imag_s, b_real_s, b_imag_s
    double precision :: a_term_real, a_term_imag, b_term_real, b_term_imag
    double precision :: pval, parity
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftb

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    pair_count = nlat / 2
    half_nlat = (nlat + 1) / 2
    symmetric_pl = .true.
    do j = 1, pair_count
        if (pl(j) /= pl(nlat - j + 1)) symmetric_pl = .false.
    end do

    ierror = 5
    if (total /= ngptot) return

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 6
    if (alloc_status /= 0) return

    grida(:, :) = 0.0
    gridb(:, :) = 0.0

!$omp parallel private(j, js, i, k, m, n, nm, offset, offset_s, &
!$omp& row_nlon, mlimit, alloc_status, row_a, row_b, row_a_s, row_b_s, &
!$omp& a_real, a_imag, b_real, b_imag, a_real_s, a_imag_s, &
!$omp& b_real_s, b_imag_s, a_term_real, a_term_imag, b_term_real, &
!$omp& b_term_imag, pval, parity, fft_rows, fft_work)
    allocate(fft_rows(4 * nt, max_nlon), fft_work(4 * nt, max_nlon), &
             stat=alloc_status)
    if (alloc_status == 0) then
    if (symmetric_pl) then
!$omp do schedule(dynamic)
    do j = 1, half_nlat
        if (j <= pair_count) then
            js = nlat - j + 1
            offset = offsets(j)
            offset_s = offsets(js)
            row_nlon = pl(j)
            mlimit = min(ntrunc, row_nlon / 2)
            fft_rows(:, 1:row_nlon) = 0.0

            do k = 1, nt
                row_a = 4 * k - 3
                row_b = 4 * k - 2
                row_a_s = 4 * k - 1
                row_b_s = 4 * k
                do m = 0, mlimit
                    a_real = 0.0d0
                    a_imag = 0.0d0
                    b_real = 0.0d0
                    b_imag = 0.0d0
                    a_real_s = 0.0d0
                    a_imag_s = 0.0d0
                    b_real_s = 0.0d0
                    b_imag_s = 0.0d0
                    nm = spectral_offset(m, ntrunc)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(j, nm))
                        a_term_real = dble(real(speca(nm, k))) * pval
                        a_term_imag = dble(aimag(speca(nm, k))) * pval
                        b_term_real = dble(real(specb(nm, k))) * pval
                        b_term_imag = dble(aimag(specb(nm, k))) * pval
                        a_real = a_real + a_term_real
                        a_imag = a_imag + a_term_imag
                        b_real = b_real + b_term_real
                        b_imag = b_imag + b_term_imag
                        a_real_s = a_real_s + parity * a_term_real
                        a_imag_s = a_imag_s + parity * a_term_imag
                        b_real_s = b_real_s + parity * b_term_real
                        b_imag_s = b_imag_s + parity * b_term_imag
                        parity = -parity
                    end do

                    if (m == 0) then
                        fft_rows(row_a, 1) = real(2.0d0 * a_real)
                        fft_rows(row_b, 1) = real(2.0d0 * b_real)
                        fft_rows(row_a_s, 1) = real(2.0d0 * a_real_s)
                        fft_rows(row_b_s, 1) = real(2.0d0 * b_real_s)
                    else
                        fft_rows(row_a, 2 * m) = real(2.0d0 * a_real)
                        fft_rows(row_b, 2 * m) = real(2.0d0 * b_real)
                        fft_rows(row_a_s, 2 * m) = real(2.0d0 * a_real_s)
                        fft_rows(row_b_s, 2 * m) = real(2.0d0 * b_real_s)
                        if (2 * m < row_nlon) then
                            fft_rows(row_a, 2 * m + 1) = real(2.0d0 * a_imag)
                            fft_rows(row_b, 2 * m + 1) = real(2.0d0 * b_imag)
                            fft_rows(row_a_s, 2 * m + 1) = &
                                real(2.0d0 * a_imag_s)
                            fft_rows(row_b_s, 2 * m + 1) = &
                                real(2.0d0 * b_imag_s)
                        end if
                    end if
                end do
            end do

            call hrfftb(4 * nt, row_nlon, fft_rows, 4 * nt, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            do k = 1, nt
                row_a = 4 * k - 3
                row_b = 4 * k - 2
                row_a_s = 4 * k - 1
                row_b_s = 4 * k
                do i = 1, row_nlon
                    grida(offset + i, k) = 0.5 * fft_rows(row_a, i)
                    gridb(offset + i, k) = 0.5 * fft_rows(row_b, i)
                    grida(offset_s + i, k) = 0.5 * fft_rows(row_a_s, i)
                    gridb(offset_s + i, k) = 0.5 * fft_rows(row_b_s, i)
                end do
            end do
        else
            offset = offsets(j)
            row_nlon = pl(j)
            mlimit = min(ntrunc, row_nlon / 2)
            fft_rows(:, 1:row_nlon) = 0.0

            do k = 1, nt
                row_a = 4 * k - 3
                row_b = 4 * k - 2
                do m = 0, mlimit
                    a_real = 0.0d0
                    a_imag = 0.0d0
                    b_real = 0.0d0
                    b_imag = 0.0d0
                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(j, nm))
                        a_real = a_real + dble(real(speca(nm, k))) * pval
                        a_imag = a_imag + dble(aimag(speca(nm, k))) * pval
                        b_real = b_real + dble(real(specb(nm, k))) * pval
                        b_imag = b_imag + dble(aimag(specb(nm, k))) * pval
                    end do

                    if (m == 0) then
                        fft_rows(row_a, 1) = real(2.0d0 * a_real)
                        fft_rows(row_b, 1) = real(2.0d0 * b_real)
                    else
                        fft_rows(row_a, 2 * m) = real(2.0d0 * a_real)
                        fft_rows(row_b, 2 * m) = real(2.0d0 * b_real)
                        if (2 * m < row_nlon) then
                            fft_rows(row_a, 2 * m + 1) = real(2.0d0 * a_imag)
                            fft_rows(row_b, 2 * m + 1) = real(2.0d0 * b_imag)
                        end if
                    end if
                end do
            end do

            call hrfftb(4 * nt, row_nlon, fft_rows, 4 * nt, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            do k = 1, nt
                row_a = 4 * k - 3
                row_b = 4 * k - 2
                do i = 1, row_nlon
                    grida(offset + i, k) = 0.5 * fft_rows(row_a, i)
                    gridb(offset + i, k) = 0.5 * fft_rows(row_b, i)
                end do
            end do
        end if
    end do
!$omp end do
    else
!$omp do schedule(dynamic)
    do j = 1, nlat
        offset = offsets(j)
        row_nlon = pl(j)
        mlimit = min(ntrunc, row_nlon / 2)
        fft_rows(:, 1:row_nlon) = 0.0

        do k = 1, nt
            row_a = 4 * k - 3
            row_b = 4 * k - 2
            do m = 0, mlimit
                a_real = 0.0d0
                a_imag = 0.0d0
                b_real = 0.0d0
                b_imag = 0.0d0
                nm = spectral_offset(m, ntrunc)
                do n = m, ntrunc
                    nm = nm + 1
                    pval = dble(basis(j, nm))
                    a_real = a_real + dble(real(speca(nm, k))) * pval
                    a_imag = a_imag + dble(aimag(speca(nm, k))) * pval
                    b_real = b_real + dble(real(specb(nm, k))) * pval
                    b_imag = b_imag + dble(aimag(specb(nm, k))) * pval
                end do

                if (m == 0) then
                    fft_rows(row_a, 1) = real(2.0d0 * a_real)
                    fft_rows(row_b, 1) = real(2.0d0 * b_real)
                else
                    fft_rows(row_a, 2 * m) = real(2.0d0 * a_real)
                    fft_rows(row_b, 2 * m) = real(2.0d0 * b_real)
                    if (2 * m < row_nlon) then
                        fft_rows(row_a, 2 * m + 1) = real(2.0d0 * a_imag)
                        fft_rows(row_b, 2 * m + 1) = real(2.0d0 * b_imag)
                    end if
                end if
            end do
        end do

        call hrfftb(4 * nt, row_nlon, fft_rows, 4 * nt, &
            reduced_gaussian_fft_wsave(1, j), fft_work)

        do k = 1, nt
            row_a = 4 * k - 3
            row_b = 4 * k - 2
            do i = 1, row_nlon
                grida(offset + i, k) = 0.5 * fft_rows(row_a, i)
                gridb(offset + i, k) = 0.5 * fft_rows(row_b, i)
            end do
        end do
    end do
!$omp end do
    end if
    end if
    if (allocated(fft_rows)) deallocate(fft_rows, fft_work)
!$omp end parallel

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_spectogrd_pair


subroutine reduced_gaussian_getgrad(dataspec, pl, basis, dbasis, sin_theta, &
                                    ugrad, vgrad, nmdim, nlat, ntrunc, nt, &
                                    ngptot, rsphere, ierror)
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(in) :: dbasis(nlat, nmdim)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(in) :: dataspec(nmdim, nt)
    real, intent(out) :: ugrad(ngptot, nt)
    real, intent(out) :: vgrad(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total, mlimit
    integer :: offsets(nlat + 1)
    double precision :: angle_step, inv_r, two_inv_r, two_inv_r_sin
    double precision :: cos_angle, sin_angle, cos_step, sin_step, cos_next
    double precision :: spec_real, spec_imag, pval, dpval
    double precision :: uvalue, vvalue
    double precision :: sum_p_real, sum_p_imag, sum_dp_real, sum_dp_imag
    double precision :: m_grad_scale, dp_grad_scale
    double precision :: u_cos(0:ntrunc, nt), u_sin(0:ntrunc, nt)
    double precision :: v_cos(0:ntrunc, nt), v_sin(0:ntrunc, nt)
    double precision, parameter :: twopi = 6.283185307179586476925286766559d0

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    ierror = 5
    if (rsphere <= 0.0) return

    total = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
        offsets(j + 1) = total
    end do

    ierror = 6
    if (total /= ngptot) return

    inv_r = 1.0d0 / dble(rsphere)
    two_inv_r = 2.0d0 * inv_r
    ugrad(:, :) = 0.0
    vgrad(:, :) = 0.0

!$omp parallel do schedule(dynamic) private(j, i, k, m, n, nm, offset, &
!$omp& row_nlon, angle_step, two_inv_r_sin, cos_angle, &
!$omp& sin_angle, cos_step, sin_step, cos_next, spec_real, spec_imag, &
!$omp& pval, dpval, uvalue, vvalue, sum_p_real, sum_p_imag, &
!$omp& sum_dp_real, sum_dp_imag, mlimit, m_grad_scale, dp_grad_scale, &
!$omp& u_cos, u_sin, v_cos, v_sin)
    do j = 1, nlat
        offset = offsets(j)
        row_nlon = pl(j)
        mlimit = min(ntrunc, row_nlon / 2)
        two_inv_r_sin = two_inv_r / dble(sin_theta(j))

        u_cos(0:mlimit, :) = 0.0d0
        u_sin(0:mlimit, :) = 0.0d0
        v_cos(0:mlimit, :) = 0.0d0
        v_sin(0:mlimit, :) = 0.0d0

        do k = 1, nt
            do m = 0, mlimit
                sum_p_real = 0.0d0
                sum_p_imag = 0.0d0
                sum_dp_real = 0.0d0
                sum_dp_imag = 0.0d0
                nm = spectral_offset(m, ntrunc)
                do n = m, ntrunc
                    nm = nm + 1
                    spec_real = dble(real(dataspec(nm, k)))
                    spec_imag = dble(aimag(dataspec(nm, k)))
                    pval = dble(basis(j, nm))
                    dpval = dble(dbasis(j, nm))
                    sum_p_real = sum_p_real + spec_real * pval
                    sum_p_imag = sum_p_imag + spec_imag * pval
                    sum_dp_real = sum_dp_real + spec_real * dpval
                    sum_dp_imag = sum_dp_imag + spec_imag * dpval
                end do

                if (m == 0) then
                    v_cos(m, k) = -sum_dp_real * inv_r
                else
                    m_grad_scale = dble(m) * two_inv_r_sin
                    dp_grad_scale = two_inv_r
                    u_cos(m, k) = -m_grad_scale * sum_p_imag
                    u_sin(m, k) = -m_grad_scale * sum_p_real
                    v_cos(m, k) = -dp_grad_scale * sum_dp_real
                    v_sin(m, k) = dp_grad_scale * sum_dp_imag
                end if
            end do
        end do

        do k = 1, nt
            do i = 1, row_nlon
                ugrad(offset + i, k) = 0.0
                vgrad(offset + i, k) = real(v_cos(0, k))
            end do

            do m = 1, mlimit
                angle_step = twopi * dble(m) / dble(row_nlon)
                cos_step = cos(angle_step)
                sin_step = sin(angle_step)
                cos_angle = 1.0d0
                sin_angle = 0.0d0
                do i = 1, row_nlon
                    uvalue = u_cos(m, k) * cos_angle + &
                        u_sin(m, k) * sin_angle
                    vvalue = v_cos(m, k) * cos_angle + &
                        v_sin(m, k) * sin_angle
                    ugrad(offset + i, k) = ugrad(offset + i, k) + real(uvalue)
                    vgrad(offset + i, k) = vgrad(offset + i, k) + real(vvalue)
                    cos_next = cos_angle * cos_step - sin_angle * sin_step
                    sin_angle = sin_angle * cos_step + cos_angle * sin_step
                    cos_angle = cos_next
                end do
            end do
        end do
    end do
!$omp end parallel do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getgrad


subroutine reduced_gaussian_getgrad_pair(speca, specb, pl, basis, dbasis, &
                                         sin_theta, a_ugrad, a_vgrad, &
                                         b_ugrad, b_vgrad, nmdim, nlat, &
                                         ntrunc, nt, ngptot, rsphere, ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(in) :: dbasis(nlat, nmdim)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(in) :: speca(nmdim, nt)
    complex, intent(in) :: specb(nmdim, nt)
    real, intent(out) :: a_ugrad(ngptot, nt)
    real, intent(out) :: a_vgrad(ngptot, nt)
    real, intent(out) :: b_ugrad(ngptot, nt)
    real, intent(out) :: b_vgrad(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, js, i, k, m, n, nm, offset, offset_s, row_nlon, total
    integer :: max_nlon, mlimit, alloc_status, pair_count, half_nlat
    integer :: offsets(nlat + 1)
    logical :: symmetric_pl
    double precision :: inv_r, two_inv_r, two_inv_r_sin
    double precision :: a_spec_real, a_spec_imag, b_spec_real, b_spec_imag
    double precision :: a_p_real, a_p_imag, a_dp_real, a_dp_imag
    double precision :: b_p_real, b_p_imag, b_dp_real, b_dp_imag
    double precision :: pval, dpval
    double precision :: a_sum_p_real, a_sum_p_imag
    double precision :: a_sum_dp_real, a_sum_dp_imag
    double precision :: b_sum_p_real, b_sum_p_imag
    double precision :: b_sum_dp_real, b_sum_dp_imag
    double precision :: a_sum_p_real_s, a_sum_p_imag_s
    double precision :: a_sum_dp_real_s, a_sum_dp_imag_s
    double precision :: b_sum_p_real_s, b_sum_p_imag_s
    double precision :: b_sum_dp_real_s, b_sum_dp_imag_s
    double precision :: m_grad_scale, dp_grad_scale
    double precision :: parity
    double precision :: a_u_cos(0:ntrunc, nt), a_u_sin(0:ntrunc, nt)
    double precision :: a_v_cos(0:ntrunc, nt), a_v_sin(0:ntrunc, nt)
    double precision :: b_u_cos(0:ntrunc, nt), b_u_sin(0:ntrunc, nt)
    double precision :: b_v_cos(0:ntrunc, nt), b_v_sin(0:ntrunc, nt)
    double precision :: a_u_cos_s(0:ntrunc, nt), a_u_sin_s(0:ntrunc, nt)
    double precision :: a_v_cos_s(0:ntrunc, nt), a_v_sin_s(0:ntrunc, nt)
    double precision :: b_u_cos_s(0:ntrunc, nt), b_u_sin_s(0:ntrunc, nt)
    double precision :: b_v_cos_s(0:ntrunc, nt), b_v_sin_s(0:ntrunc, nt)
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftb

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    ierror = 5
    if (rsphere <= 0.0) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    pair_count = nlat / 2
    half_nlat = (nlat + 1) / 2
    symmetric_pl = .true.
    do j = 1, pair_count
        if (pl(j) /= pl(nlat - j + 1)) symmetric_pl = .false.
    end do

    ierror = 6
    if (total /= ngptot) return

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 7
    if (alloc_status /= 0) return

    inv_r = 1.0d0 / dble(rsphere)
    two_inv_r = 2.0d0 * inv_r
    a_ugrad(:, :) = 0.0
    a_vgrad(:, :) = 0.0
    b_ugrad(:, :) = 0.0
    b_vgrad(:, :) = 0.0

!$omp parallel private(j, js, i, k, m, n, nm, offset, offset_s, row_nlon, mlimit, &
!$omp& alloc_status, two_inv_r_sin, a_spec_real, a_spec_imag, &
!$omp& b_spec_real, b_spec_imag, a_p_real, a_p_imag, a_dp_real, &
!$omp& a_dp_imag, b_p_real, b_p_imag, b_dp_real, b_dp_imag, pval, &
!$omp& dpval, a_sum_p_real, a_sum_p_imag, a_sum_dp_real, a_sum_dp_imag, &
!$omp& b_sum_p_real, b_sum_p_imag, b_sum_dp_real, b_sum_dp_imag, &
!$omp& a_sum_p_real_s, a_sum_p_imag_s, a_sum_dp_real_s, &
!$omp& a_sum_dp_imag_s, b_sum_p_real_s, b_sum_p_imag_s, &
!$omp& b_sum_dp_real_s, b_sum_dp_imag_s, m_grad_scale, dp_grad_scale, &
!$omp& a_u_cos, a_u_sin, a_v_cos, a_v_sin, b_u_cos, b_u_sin, &
!$omp& b_v_cos, b_v_sin, a_u_cos_s, a_u_sin_s, a_v_cos_s, a_v_sin_s, &
!$omp& b_u_cos_s, b_u_sin_s, b_v_cos_s, b_v_sin_s, parity, &
!$omp& fft_rows, fft_work)
    allocate(fft_rows(4, max_nlon), fft_work(4, max_nlon), &
             stat=alloc_status)
    if (alloc_status == 0) then
    if (symmetric_pl) then
!$omp do schedule(dynamic)
    do j = 1, half_nlat
        if (j <= pair_count) then
            js = nlat - j + 1
            offset = offsets(j)
            offset_s = offsets(js)
            row_nlon = pl(j)
            two_inv_r_sin = two_inv_r / dble(sin_theta(j))
            mlimit = min(ntrunc, row_nlon / 2)

            a_u_cos(0:mlimit, :) = 0.0d0
            a_u_sin(0:mlimit, :) = 0.0d0
            a_v_cos(0:mlimit, :) = 0.0d0
            a_v_sin(0:mlimit, :) = 0.0d0
            b_u_cos(0:mlimit, :) = 0.0d0
            b_u_sin(0:mlimit, :) = 0.0d0
            b_v_cos(0:mlimit, :) = 0.0d0
            b_v_sin(0:mlimit, :) = 0.0d0
            a_u_cos_s(0:mlimit, :) = 0.0d0
            a_u_sin_s(0:mlimit, :) = 0.0d0
            a_v_cos_s(0:mlimit, :) = 0.0d0
            a_v_sin_s(0:mlimit, :) = 0.0d0
            b_u_cos_s(0:mlimit, :) = 0.0d0
            b_u_sin_s(0:mlimit, :) = 0.0d0
            b_v_cos_s(0:mlimit, :) = 0.0d0
            b_v_sin_s(0:mlimit, :) = 0.0d0

            do k = 1, nt
                do m = 0, mlimit
                    a_sum_p_real = 0.0d0
                    a_sum_p_imag = 0.0d0
                    a_sum_dp_real = 0.0d0
                    a_sum_dp_imag = 0.0d0
                    b_sum_p_real = 0.0d0
                    b_sum_p_imag = 0.0d0
                    b_sum_dp_real = 0.0d0
                    b_sum_dp_imag = 0.0d0
                    a_sum_p_real_s = 0.0d0
                    a_sum_p_imag_s = 0.0d0
                    a_sum_dp_real_s = 0.0d0
                    a_sum_dp_imag_s = 0.0d0
                    b_sum_p_real_s = 0.0d0
                    b_sum_p_imag_s = 0.0d0
                    b_sum_dp_real_s = 0.0d0
                    b_sum_dp_imag_s = 0.0d0
                    nm = spectral_offset(m, ntrunc)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        a_spec_real = dble(real(speca(nm, k)))
                        a_spec_imag = dble(aimag(speca(nm, k)))
                        b_spec_real = dble(real(specb(nm, k)))
                        b_spec_imag = dble(aimag(specb(nm, k)))
                        pval = dble(basis(j, nm))
                        dpval = dble(dbasis(j, nm))
                        a_p_real = a_spec_real * pval
                        a_p_imag = a_spec_imag * pval
                        a_dp_real = a_spec_real * dpval
                        a_dp_imag = a_spec_imag * dpval
                        b_p_real = b_spec_real * pval
                        b_p_imag = b_spec_imag * pval
                        b_dp_real = b_spec_real * dpval
                        b_dp_imag = b_spec_imag * dpval
                        a_sum_p_real = a_sum_p_real + a_p_real
                        a_sum_p_imag = a_sum_p_imag + a_p_imag
                        a_sum_dp_real = a_sum_dp_real + a_dp_real
                        a_sum_dp_imag = a_sum_dp_imag + a_dp_imag
                        b_sum_p_real = b_sum_p_real + b_p_real
                        b_sum_p_imag = b_sum_p_imag + b_p_imag
                        b_sum_dp_real = b_sum_dp_real + b_dp_real
                        b_sum_dp_imag = b_sum_dp_imag + b_dp_imag
                        a_sum_p_real_s = a_sum_p_real_s + parity * a_p_real
                        a_sum_p_imag_s = a_sum_p_imag_s + parity * a_p_imag
                        a_sum_dp_real_s = a_sum_dp_real_s - parity * a_dp_real
                        a_sum_dp_imag_s = a_sum_dp_imag_s - parity * a_dp_imag
                        b_sum_p_real_s = b_sum_p_real_s + parity * b_p_real
                        b_sum_p_imag_s = b_sum_p_imag_s + parity * b_p_imag
                        b_sum_dp_real_s = b_sum_dp_real_s - parity * b_dp_real
                        b_sum_dp_imag_s = b_sum_dp_imag_s - parity * b_dp_imag
                        parity = -parity
                    end do

                    if (m == 0) then
                        a_v_cos(m, k) = -a_sum_dp_real * inv_r
                        b_v_cos(m, k) = -b_sum_dp_real * inv_r
                        a_v_cos_s(m, k) = -a_sum_dp_real_s * inv_r
                        b_v_cos_s(m, k) = -b_sum_dp_real_s * inv_r
                    else
                        m_grad_scale = dble(m) * two_inv_r_sin
                        dp_grad_scale = two_inv_r
                        a_u_cos(m, k) = -m_grad_scale * a_sum_p_imag
                        a_u_sin(m, k) = -m_grad_scale * a_sum_p_real
                        a_v_cos(m, k) = -dp_grad_scale * a_sum_dp_real
                        a_v_sin(m, k) = dp_grad_scale * a_sum_dp_imag
                        b_u_cos(m, k) = -m_grad_scale * b_sum_p_imag
                        b_u_sin(m, k) = -m_grad_scale * b_sum_p_real
                        b_v_cos(m, k) = -dp_grad_scale * b_sum_dp_real
                        b_v_sin(m, k) = dp_grad_scale * b_sum_dp_imag
                        a_u_cos_s(m, k) = -m_grad_scale * a_sum_p_imag_s
                        a_u_sin_s(m, k) = -m_grad_scale * a_sum_p_real_s
                        a_v_cos_s(m, k) = -dp_grad_scale * a_sum_dp_real_s
                        a_v_sin_s(m, k) = dp_grad_scale * a_sum_dp_imag_s
                        b_u_cos_s(m, k) = -m_grad_scale * b_sum_p_imag_s
                        b_u_sin_s(m, k) = -m_grad_scale * b_sum_p_real_s
                        b_v_cos_s(m, k) = -dp_grad_scale * b_sum_dp_real_s
                        b_v_sin_s(m, k) = dp_grad_scale * b_sum_dp_imag_s
                    end if
                end do
            end do

            do k = 1, nt
                fft_rows(:, 1:row_nlon) = 0.0
                fft_rows(2, 1) = 2.0 * real(a_v_cos(0, k))
                fft_rows(4, 1) = 2.0 * real(b_v_cos(0, k))

                do m = 1, mlimit
                    fft_rows(1, 2 * m) = real(a_u_cos(m, k))
                    fft_rows(2, 2 * m) = real(a_v_cos(m, k))
                    fft_rows(3, 2 * m) = real(b_u_cos(m, k))
                    fft_rows(4, 2 * m) = real(b_v_cos(m, k))
                    if (2 * m < row_nlon) then
                        fft_rows(1, 2 * m + 1) = -real(a_u_sin(m, k))
                        fft_rows(2, 2 * m + 1) = -real(a_v_sin(m, k))
                        fft_rows(3, 2 * m + 1) = -real(b_u_sin(m, k))
                        fft_rows(4, 2 * m + 1) = -real(b_v_sin(m, k))
                    end if
                end do

                call hrfftb(4, row_nlon, fft_rows, 4, &
                    reduced_gaussian_fft_wsave(1, j), fft_work)

                do i = 1, row_nlon
                    a_ugrad(offset + i, k) = 0.5 * fft_rows(1, i)
                    a_vgrad(offset + i, k) = 0.5 * fft_rows(2, i)
                    b_ugrad(offset + i, k) = 0.5 * fft_rows(3, i)
                    b_vgrad(offset + i, k) = 0.5 * fft_rows(4, i)
                end do

                fft_rows(:, 1:row_nlon) = 0.0
                fft_rows(2, 1) = 2.0 * real(a_v_cos_s(0, k))
                fft_rows(4, 1) = 2.0 * real(b_v_cos_s(0, k))

                do m = 1, mlimit
                    fft_rows(1, 2 * m) = real(a_u_cos_s(m, k))
                    fft_rows(2, 2 * m) = real(a_v_cos_s(m, k))
                    fft_rows(3, 2 * m) = real(b_u_cos_s(m, k))
                    fft_rows(4, 2 * m) = real(b_v_cos_s(m, k))
                    if (2 * m < row_nlon) then
                        fft_rows(1, 2 * m + 1) = -real(a_u_sin_s(m, k))
                        fft_rows(2, 2 * m + 1) = -real(a_v_sin_s(m, k))
                        fft_rows(3, 2 * m + 1) = -real(b_u_sin_s(m, k))
                        fft_rows(4, 2 * m + 1) = -real(b_v_sin_s(m, k))
                    end if
                end do

                call hrfftb(4, row_nlon, fft_rows, 4, &
                    reduced_gaussian_fft_wsave(1, js), fft_work)

                do i = 1, row_nlon
                    a_ugrad(offset_s + i, k) = 0.5 * fft_rows(1, i)
                    a_vgrad(offset_s + i, k) = 0.5 * fft_rows(2, i)
                    b_ugrad(offset_s + i, k) = 0.5 * fft_rows(3, i)
                    b_vgrad(offset_s + i, k) = 0.5 * fft_rows(4, i)
                end do
            end do
        else
            offset = offsets(j)
            row_nlon = pl(j)
            two_inv_r_sin = two_inv_r / dble(sin_theta(j))
            mlimit = min(ntrunc, row_nlon / 2)

            a_u_cos(0:mlimit, :) = 0.0d0
            a_u_sin(0:mlimit, :) = 0.0d0
            a_v_cos(0:mlimit, :) = 0.0d0
            a_v_sin(0:mlimit, :) = 0.0d0
            b_u_cos(0:mlimit, :) = 0.0d0
            b_u_sin(0:mlimit, :) = 0.0d0
            b_v_cos(0:mlimit, :) = 0.0d0
            b_v_sin(0:mlimit, :) = 0.0d0

            do k = 1, nt
                do m = 0, mlimit
                    a_sum_p_real = 0.0d0
                    a_sum_p_imag = 0.0d0
                    a_sum_dp_real = 0.0d0
                    a_sum_dp_imag = 0.0d0
                    b_sum_p_real = 0.0d0
                    b_sum_p_imag = 0.0d0
                    b_sum_dp_real = 0.0d0
                    b_sum_dp_imag = 0.0d0
                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        a_spec_real = dble(real(speca(nm, k)))
                        a_spec_imag = dble(aimag(speca(nm, k)))
                        b_spec_real = dble(real(specb(nm, k)))
                        b_spec_imag = dble(aimag(specb(nm, k)))
                        pval = dble(basis(j, nm))
                        dpval = dble(dbasis(j, nm))
                        a_sum_p_real = a_sum_p_real + a_spec_real * pval
                        a_sum_p_imag = a_sum_p_imag + a_spec_imag * pval
                        a_sum_dp_real = a_sum_dp_real + a_spec_real * dpval
                        a_sum_dp_imag = a_sum_dp_imag + a_spec_imag * dpval
                        b_sum_p_real = b_sum_p_real + b_spec_real * pval
                        b_sum_p_imag = b_sum_p_imag + b_spec_imag * pval
                        b_sum_dp_real = b_sum_dp_real + b_spec_real * dpval
                        b_sum_dp_imag = b_sum_dp_imag + b_spec_imag * dpval
                    end do

                    if (m == 0) then
                        a_v_cos(m, k) = -a_sum_dp_real * inv_r
                        b_v_cos(m, k) = -b_sum_dp_real * inv_r
                    else
                        m_grad_scale = dble(m) * two_inv_r_sin
                        dp_grad_scale = two_inv_r
                        a_u_cos(m, k) = -m_grad_scale * a_sum_p_imag
                        a_u_sin(m, k) = -m_grad_scale * a_sum_p_real
                        a_v_cos(m, k) = -dp_grad_scale * a_sum_dp_real
                        a_v_sin(m, k) = dp_grad_scale * a_sum_dp_imag
                        b_u_cos(m, k) = -m_grad_scale * b_sum_p_imag
                        b_u_sin(m, k) = -m_grad_scale * b_sum_p_real
                        b_v_cos(m, k) = -dp_grad_scale * b_sum_dp_real
                        b_v_sin(m, k) = dp_grad_scale * b_sum_dp_imag
                    end if
                end do
            end do

            do k = 1, nt
                fft_rows(:, 1:row_nlon) = 0.0
                fft_rows(2, 1) = 2.0 * real(a_v_cos(0, k))
                fft_rows(4, 1) = 2.0 * real(b_v_cos(0, k))

                do m = 1, mlimit
                    fft_rows(1, 2 * m) = real(a_u_cos(m, k))
                    fft_rows(2, 2 * m) = real(a_v_cos(m, k))
                    fft_rows(3, 2 * m) = real(b_u_cos(m, k))
                    fft_rows(4, 2 * m) = real(b_v_cos(m, k))
                    if (2 * m < row_nlon) then
                        fft_rows(1, 2 * m + 1) = -real(a_u_sin(m, k))
                        fft_rows(2, 2 * m + 1) = -real(a_v_sin(m, k))
                        fft_rows(3, 2 * m + 1) = -real(b_u_sin(m, k))
                        fft_rows(4, 2 * m + 1) = -real(b_v_sin(m, k))
                    end if
                end do

                call hrfftb(4, row_nlon, fft_rows, 4, &
                    reduced_gaussian_fft_wsave(1, j), fft_work)

                do i = 1, row_nlon
                    a_ugrad(offset + i, k) = 0.5 * fft_rows(1, i)
                    a_vgrad(offset + i, k) = 0.5 * fft_rows(2, i)
                    b_ugrad(offset + i, k) = 0.5 * fft_rows(3, i)
                    b_vgrad(offset + i, k) = 0.5 * fft_rows(4, i)
                end do
            end do
        end if
    end do
!$omp end do
    else
!$omp do schedule(dynamic)
    do j = 1, nlat
        offset = offsets(j)
        row_nlon = pl(j)
        two_inv_r_sin = two_inv_r / dble(sin_theta(j))
        mlimit = min(ntrunc, row_nlon / 2)

        a_u_cos(0:mlimit, :) = 0.0d0
        a_u_sin(0:mlimit, :) = 0.0d0
        a_v_cos(0:mlimit, :) = 0.0d0
        a_v_sin(0:mlimit, :) = 0.0d0
        b_u_cos(0:mlimit, :) = 0.0d0
        b_u_sin(0:mlimit, :) = 0.0d0
        b_v_cos(0:mlimit, :) = 0.0d0
        b_v_sin(0:mlimit, :) = 0.0d0

        do k = 1, nt
            do m = 0, mlimit
                a_sum_p_real = 0.0d0
                a_sum_p_imag = 0.0d0
                a_sum_dp_real = 0.0d0
                a_sum_dp_imag = 0.0d0
                b_sum_p_real = 0.0d0
                b_sum_p_imag = 0.0d0
                b_sum_dp_real = 0.0d0
                b_sum_dp_imag = 0.0d0
                nm = spectral_offset(m, ntrunc)
                do n = m, ntrunc
                    nm = nm + 1
                    a_spec_real = dble(real(speca(nm, k)))
                    a_spec_imag = dble(aimag(speca(nm, k)))
                    b_spec_real = dble(real(specb(nm, k)))
                    b_spec_imag = dble(aimag(specb(nm, k)))
                    pval = dble(basis(j, nm))
                    dpval = dble(dbasis(j, nm))
                    a_sum_p_real = a_sum_p_real + a_spec_real * pval
                    a_sum_p_imag = a_sum_p_imag + a_spec_imag * pval
                    a_sum_dp_real = a_sum_dp_real + a_spec_real * dpval
                    a_sum_dp_imag = a_sum_dp_imag + a_spec_imag * dpval
                    b_sum_p_real = b_sum_p_real + b_spec_real * pval
                    b_sum_p_imag = b_sum_p_imag + b_spec_imag * pval
                    b_sum_dp_real = b_sum_dp_real + b_spec_real * dpval
                    b_sum_dp_imag = b_sum_dp_imag + b_spec_imag * dpval
                end do

                if (m == 0) then
                    a_v_cos(m, k) = -a_sum_dp_real * inv_r
                    b_v_cos(m, k) = -b_sum_dp_real * inv_r
                else
                    m_grad_scale = dble(m) * two_inv_r_sin
                    dp_grad_scale = two_inv_r
                    a_u_cos(m, k) = -m_grad_scale * a_sum_p_imag
                    a_u_sin(m, k) = -m_grad_scale * a_sum_p_real
                    a_v_cos(m, k) = -dp_grad_scale * a_sum_dp_real
                    a_v_sin(m, k) = dp_grad_scale * a_sum_dp_imag
                    b_u_cos(m, k) = -m_grad_scale * b_sum_p_imag
                    b_u_sin(m, k) = -m_grad_scale * b_sum_p_real
                    b_v_cos(m, k) = -dp_grad_scale * b_sum_dp_real
                    b_v_sin(m, k) = dp_grad_scale * b_sum_dp_imag
                end if
            end do
        end do

        do k = 1, nt
            fft_rows(:, 1:row_nlon) = 0.0
            fft_rows(2, 1) = 2.0 * real(a_v_cos(0, k))
            fft_rows(4, 1) = 2.0 * real(b_v_cos(0, k))

            do m = 1, mlimit
                fft_rows(1, 2 * m) = real(a_u_cos(m, k))
                fft_rows(2, 2 * m) = real(a_v_cos(m, k))
                fft_rows(3, 2 * m) = real(b_u_cos(m, k))
                fft_rows(4, 2 * m) = real(b_v_cos(m, k))
                if (2 * m < row_nlon) then
                    fft_rows(1, 2 * m + 1) = -real(a_u_sin(m, k))
                    fft_rows(2, 2 * m + 1) = -real(a_v_sin(m, k))
                    fft_rows(3, 2 * m + 1) = -real(b_u_sin(m, k))
                    fft_rows(4, 2 * m + 1) = -real(b_v_sin(m, k))
                end if
            end do

            call hrfftb(4, row_nlon, fft_rows, 4, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            do i = 1, row_nlon
                a_ugrad(offset + i, k) = 0.5 * fft_rows(1, i)
                a_vgrad(offset + i, k) = 0.5 * fft_rows(2, i)
                b_ugrad(offset + i, k) = 0.5 * fft_rows(3, i)
                b_vgrad(offset + i, k) = 0.5 * fft_rows(4, i)
            end do
        end do
    end do
!$omp end do
    end if
    end if
    if (allocated(fft_rows)) deallocate(fft_rows, fft_work)
!$omp end parallel

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getgrad_pair


subroutine reduced_gaussian_getvrtdivspec(ugrid, vgrid, pl, weights, basis, &
                                          dbasis, sin_theta, vrtspec, divspec, &
                                          ngptot, nlat, ntrunc, nt, rsphere, &
                                          ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: ugrid(ngptot, nt)
    real, intent(in) :: vgrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis((ntrunc + 1) * (ntrunc + 2) / 2, nlat)
    real, intent(in) :: dbasis((ntrunc + 1) * (ntrunc + 2) / 2, nlat)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(out) :: vrtspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    complex, intent(out) :: divspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    integer :: j, js, i, k, m, n, nm, offset, row_nlon, total, first_j, twom
    integer :: max_nlon, mlimit, alloc_status, pair_count, half_nlat
    integer :: offsets(nlat + 1), first_pair_for_m(0:ntrunc)
    logical :: symmetric_pl, center_active(0:ntrunc)
    double precision :: inv_nlon, inv_r
    double precision :: u_mean, v_mean, u_cos, u_sin, v_cos, v_sin
    double precision :: u_mean_s, v_mean_s, u_cos_s, u_sin_s, v_cos_s, v_sin_s
    double precision :: pval, dpval, wval, div_real, div_imag, vrt_real, vrt_imag
    double precision :: p_scale, dp_scale, p_weighted, dp_weighted
    double precision :: parity, p_u_cos, p_u_sin, p_v_cos, p_v_sin
    double precision :: dp_u_cos, dp_u_sin, dp_v_cos, dp_v_sin
    double precision :: div_real_acc(0:ntrunc), div_imag_acc(0:ntrunc)
    double precision :: vrt_real_acc(0:ntrunc), vrt_imag_acc(0:ntrunc)
    double precision :: weight_inv_r(nlat), weight_inv_r_sin(nlat)
    double precision, allocatable :: u_cos_fft(:, :), u_sin_fft(:, :)
    double precision, allocatable :: v_cos_fft(:, :), v_sin_fft(:, :)
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftf

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nt < 1) return

    ierror = 4
    if (rsphere <= 0.0) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    pair_count = nlat / 2
    half_nlat = (nlat + 1) / 2
    symmetric_pl = .true.
    do j = 1, pair_count
        if (pl(j) /= pl(nlat - j + 1)) symmetric_pl = .false.
    end do

    ierror = 5
    if (total /= ngptot) return

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 6
    if (alloc_status /= 0) return

    allocate(u_cos_fft(nlat, 0:ntrunc), u_sin_fft(nlat, 0:ntrunc), &
             v_cos_fft(nlat, 0:ntrunc), v_sin_fft(nlat, 0:ntrunc), &
             stat=alloc_status)
    ierror = 7
    if (alloc_status /= 0) return

    inv_r = 1.0d0 / dble(rsphere)
    do j = 1, nlat
        weight_inv_r(j) = dble(weights(j)) * inv_r
        weight_inv_r_sin(j) = weight_inv_r(j) / dble(sin_theta(j))
    end do
    do m = 0, ntrunc
        first_j = 1
        twom = 2 * m
        do while (first_j <= pair_count .and. twom > pl(first_j))
            first_j = first_j + 1
        end do
        first_pair_for_m(m) = first_j
        center_active(m) = half_nlat > pair_count .and. twom <= pl(half_nlat)
    end do
    vrtspec(:, :) = (0.0, 0.0)
    divspec(:, :) = (0.0, 0.0)

    do k = 1, nt
!$omp parallel private(j, i, m, offset, row_nlon, mlimit, inv_nlon, &
!$omp& twom, fft_rows, fft_work)
        allocate(fft_rows(2, max_nlon), fft_work(2, max_nlon))
!$omp do schedule(dynamic)
        do j = 1, nlat
            offset = offsets(j)
            row_nlon = pl(j)
            mlimit = min(ntrunc, row_nlon / 2)
            inv_nlon = 1.0d0 / dble(row_nlon)
            fft_rows(:, 1:row_nlon) = 0.0
            do i = 1, row_nlon
                fft_rows(1, i) = ugrid(offset + i, k)
                fft_rows(2, i) = vgrid(offset + i, k)
            end do

            call hrfftf(2, row_nlon, fft_rows, 2, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            u_cos_fft(j, 0) = dble(fft_rows(1, 1)) * inv_nlon
            v_cos_fft(j, 0) = dble(fft_rows(2, 1)) * inv_nlon
            do m = 1, mlimit
                twom = 2 * m
                u_cos_fft(j, m) = dble(fft_rows(1, twom)) * inv_nlon
                v_cos_fft(j, m) = dble(fft_rows(2, twom)) * inv_nlon
                if (twom < row_nlon) then
                    u_sin_fft(j, m) = -dble(fft_rows(1, twom + 1)) * inv_nlon
                    v_sin_fft(j, m) = -dble(fft_rows(2, twom + 1)) * inv_nlon
                else
                    u_sin_fft(j, m) = 0.0d0
                    v_sin_fft(j, m) = 0.0d0
                end if
            end do
        end do
!$omp end do
        deallocate(fft_rows, fft_work)
!$omp end parallel

        if (symmetric_pl) then
!$omp parallel do schedule(dynamic) private(m, j, js, n, nm, row_nlon, first_j, &
!$omp& u_mean, v_mean, u_cos, u_sin, v_cos, v_sin, &
!$omp& u_mean_s, v_mean_s, u_cos_s, u_sin_s, v_cos_s, v_sin_s, &
!$omp& pval, dpval, wval, div_real, div_imag, vrt_real, vrt_imag, &
!$omp& p_scale, dp_scale, p_weighted, dp_weighted, parity, &
!$omp& p_u_cos, p_u_sin, p_v_cos, p_v_sin, dp_u_cos, dp_u_sin, &
!$omp& dp_v_cos, dp_v_sin, div_real_acc, div_imag_acc, &
!$omp& vrt_real_acc, vrt_imag_acc)
        do m = 0, ntrunc
            div_real_acc(m:ntrunc) = 0.0d0
            div_imag_acc(m:ntrunc) = 0.0d0
            vrt_real_acc(m:ntrunc) = 0.0d0
            vrt_imag_acc(m:ntrunc) = 0.0d0
            if (m == 0) then
                do j = 1, pair_count
                    js = nlat - j + 1
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)
                    u_mean_s = u_cos_fft(js, 0)
                    v_mean_s = v_cos_fft(js, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(nm, j))
                        dp_weighted = dp_scale * dpval
                        div_real = dp_weighted * (v_mean - parity * v_mean_s)
                        vrt_real = -dp_weighted * (u_mean - parity * u_mean_s)
                        div_real_acc(n) = div_real_acc(n) + div_real
                        vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                        parity = -parity
                    end do
                end do

                if (half_nlat > pair_count) then
                    j = half_nlat
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(nm, j))
                        dp_weighted = dp_scale * dpval
                        div_real = dp_weighted * v_mean
                        vrt_real = -dp_weighted * u_mean
                        div_real_acc(n) = div_real_acc(n) + div_real
                        vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                    end do
                end if
            else
                first_j = first_pair_for_m(m)
                do j = first_j, pair_count
                    row_nlon = pl(j)
                    js = nlat - j + 1
                    u_cos = u_cos_fft(j, m)
                    u_sin = u_sin_fft(j, m)
                    v_cos = v_cos_fft(j, m)
                    v_sin = v_sin_fft(j, m)
                    u_cos_s = u_cos_fft(js, m)
                    u_sin_s = u_sin_fft(js, m)
                    v_cos_s = v_cos_fft(js, m)
                    v_sin_s = v_sin_fft(js, m)

                    nm = spectral_offset(m, ntrunc)
                    p_scale = dble(m) * weight_inv_r_sin(j)
                    dp_scale = weight_inv_r(j)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(nm, j))
                        dpval = dble(dbasis(nm, j))
                        p_weighted = p_scale * pval
                        dp_weighted = dp_scale * dpval
                        p_u_sin = u_sin + parity * u_sin_s
                        p_u_cos = u_cos + parity * u_cos_s
                        p_v_sin = v_sin + parity * v_sin_s
                        p_v_cos = v_cos + parity * v_cos_s
                        dp_u_sin = u_sin - parity * u_sin_s
                        dp_u_cos = u_cos - parity * u_cos_s
                        dp_v_sin = v_sin - parity * v_sin_s
                        dp_v_cos = v_cos - parity * v_cos_s
                        div_real = p_weighted * p_u_sin + dp_weighted * dp_v_cos
                        div_imag = p_weighted * p_u_cos - dp_weighted * dp_v_sin
                        vrt_real = -dp_weighted * dp_u_cos + p_weighted * p_v_sin
                        vrt_imag = dp_weighted * dp_u_sin + p_weighted * p_v_cos
                        div_real_acc(n) = div_real_acc(n) + div_real
                        div_imag_acc(n) = div_imag_acc(n) + div_imag
                        vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                        vrt_imag_acc(n) = vrt_imag_acc(n) + vrt_imag
                        parity = -parity
                    end do
                end do

                if (center_active(m)) then
                    j = half_nlat
                    u_cos = u_cos_fft(j, m)
                    u_sin = u_sin_fft(j, m)
                    v_cos = v_cos_fft(j, m)
                    v_sin = v_sin_fft(j, m)

                    nm = spectral_offset(m, ntrunc)
                    p_scale = dble(m) * weight_inv_r_sin(j)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(nm, j))
                        dpval = dble(dbasis(nm, j))
                        p_weighted = p_scale * pval
                        dp_weighted = dp_scale * dpval
                        div_real = p_weighted * u_sin + dp_weighted * v_cos
                        div_imag = p_weighted * u_cos - dp_weighted * v_sin
                        vrt_real = -dp_weighted * u_cos + p_weighted * v_sin
                        vrt_imag = dp_weighted * u_sin + p_weighted * v_cos
                        div_real_acc(n) = div_real_acc(n) + div_real
                        div_imag_acc(n) = div_imag_acc(n) + div_imag
                        vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                        vrt_imag_acc(n) = vrt_imag_acc(n) + vrt_imag
                    end do
                end if
            end if
            nm = spectral_offset(m, ntrunc)
            do n = m, ntrunc
                nm = nm + 1
                divspec(nm, k) = cmplx(real(div_real_acc(n)), &
                    real(div_imag_acc(n)))
                vrtspec(nm, k) = cmplx(real(vrt_real_acc(n)), &
                    real(vrt_imag_acc(n)))
            end do
        end do
!$omp end parallel do
        else
!$omp parallel do schedule(dynamic) private(m, j, n, nm, row_nlon, &
!$omp& u_mean, v_mean, u_cos, u_sin, v_cos, v_sin, pval, &
!$omp& dpval, wval, div_real, div_imag, vrt_real, vrt_imag, p_scale, &
!$omp& dp_scale, p_weighted, dp_weighted, div_real_acc, div_imag_acc, &
!$omp& vrt_real_acc, vrt_imag_acc)
        do m = 0, ntrunc
            div_real_acc(m:ntrunc) = 0.0d0
            div_imag_acc(m:ntrunc) = 0.0d0
            vrt_real_acc(m:ntrunc) = 0.0d0
            vrt_imag_acc(m:ntrunc) = 0.0d0
            if (m == 0) then
                do j = 1, nlat
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(nm, j))
                        dp_weighted = dp_scale * dpval
                        div_real = dp_weighted * v_mean
                        vrt_real = -dp_weighted * u_mean
                        div_real_acc(n) = div_real_acc(n) + div_real
                        vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                    end do
                end do
            else
                do j = 1, nlat
                    row_nlon = pl(j)
                    if (2 * m <= row_nlon) then
                        u_cos = u_cos_fft(j, m)
                        u_sin = u_sin_fft(j, m)
                        v_cos = v_cos_fft(j, m)
                        v_sin = v_sin_fft(j, m)

                        nm = spectral_offset(m, ntrunc)
                        p_scale = dble(m) * weight_inv_r_sin(j)
                        dp_scale = weight_inv_r(j)
                        do n = m, ntrunc
                            nm = nm + 1
                            pval = dble(basis(nm, j))
                            dpval = dble(dbasis(nm, j))
                            p_weighted = p_scale * pval
                            dp_weighted = dp_scale * dpval
                            div_real = p_weighted * u_sin + dp_weighted * v_cos
                            div_imag = p_weighted * u_cos - dp_weighted * v_sin
                            vrt_real = -dp_weighted * u_cos + p_weighted * v_sin
                            vrt_imag = dp_weighted * u_sin + p_weighted * v_cos
                            div_real_acc(n) = div_real_acc(n) + div_real
                            div_imag_acc(n) = div_imag_acc(n) + div_imag
                            vrt_real_acc(n) = vrt_real_acc(n) + vrt_real
                            vrt_imag_acc(n) = vrt_imag_acc(n) + vrt_imag
                        end do
                    end if
                end do
            end if
            nm = spectral_offset(m, ntrunc)
            do n = m, ntrunc
                nm = nm + 1
                divspec(nm, k) = cmplx(real(div_real_acc(n)), &
                    real(div_imag_acc(n)))
                vrtspec(nm, k) = cmplx(real(vrt_real_acc(n)), &
                    real(vrt_imag_acc(n)))
            end do
        end do
!$omp end parallel do
        end if
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getvrtdivspec


subroutine reduced_gaussian_getvrtspec(ugrid, vgrid, pl, weights, basis, &
                                       dbasis, sin_theta, vrtspec, &
                                       ngptot, nlat, ntrunc, nt, rsphere, &
                                       ierror)
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: ugrid(ngptot, nt)
    real, intent(in) :: vgrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(out) :: vrtspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    call reduced_gaussian_getvecspec_component(ugrid, vgrid, pl, weights, &
        basis, dbasis, sin_theta, vrtspec, ngptot, nlat, ntrunc, nt, &
        rsphere, 1, ierror)

end subroutine reduced_gaussian_getvrtspec


subroutine reduced_gaussian_getdivspec(ugrid, vgrid, pl, weights, basis, &
                                       dbasis, sin_theta, divspec, &
                                       ngptot, nlat, ntrunc, nt, rsphere, &
                                       ierror)
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: ugrid(ngptot, nt)
    real, intent(in) :: vgrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(out) :: divspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    call reduced_gaussian_getvecspec_component(ugrid, vgrid, pl, weights, &
        basis, dbasis, sin_theta, divspec, ngptot, nlat, ntrunc, nt, &
        rsphere, 2, ierror)

end subroutine reduced_gaussian_getdivspec


subroutine reduced_gaussian_getvecspec_component(ugrid, vgrid, pl, weights, &
                                                 basis, dbasis, sin_theta, &
                                                 outspec, ngptot, nlat, &
                                                 ntrunc, nt, rsphere, &
                                                 component, ierror)
    use reduced_gaussian_fft_cache_mod, only: reduced_gaussian_prepare_fft_cache, &
        reduced_gaussian_fft_wsave
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt, component
    real, intent(in) :: ugrid(ngptot, nt)
    real, intent(in) :: vgrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(out) :: outspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    integer :: j, js, i, k, m, n, nm, offset, row_nlon, total, first_j, twom
    integer :: max_nlon, mlimit, alloc_status, pair_count, half_nlat
    integer :: offsets(nlat + 1), first_pair_for_m(0:ntrunc)
    logical :: symmetric_pl, center_active(0:ntrunc)
    double precision :: inv_nlon, inv_r
    double precision :: u_mean, v_mean, u_cos, u_sin, v_cos, v_sin
    double precision :: u_mean_s, v_mean_s, u_cos_s, u_sin_s
    double precision :: v_cos_s, v_sin_s
    double precision :: pval, dpval, wval, value_real, value_imag
    double precision :: p_scale, dp_scale, p_weighted, dp_weighted
    double precision :: parity, p_u_cos, p_u_sin, p_v_cos, p_v_sin
    double precision :: dp_u_cos, dp_u_sin, dp_v_cos, dp_v_sin
    double precision :: value_real_acc(0:ntrunc), value_imag_acc(0:ntrunc)
    double precision :: weight_inv_r(nlat), weight_inv_r_sin(nlat)
    double precision, allocatable :: u_cos_fft(:, :), u_sin_fft(:, :)
    double precision, allocatable :: v_cos_fft(:, :), v_sin_fft(:, :)
    real, allocatable :: fft_rows(:, :), fft_work(:, :)

    external :: hrfftf

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nt < 1) return

    ierror = 4
    if (rsphere <= 0.0) return

    ierror = 5
    if (component < 1 .or. component > 2) return

    total = 0
    max_nlon = 0
    offsets(1) = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
        if (pl(j) > max_nlon) max_nlon = pl(j)
        offsets(j + 1) = total
    end do

    pair_count = nlat / 2
    half_nlat = (nlat + 1) / 2
    symmetric_pl = .true.
    do j = 1, pair_count
        if (pl(j) /= pl(nlat - j + 1)) symmetric_pl = .false.
    end do

    ierror = 6
    if (total /= ngptot) return

    call reduced_gaussian_prepare_fft_cache(pl, nlat, alloc_status)
    ierror = 7
    if (alloc_status /= 0) return

    allocate(u_cos_fft(nlat, 0:ntrunc), u_sin_fft(nlat, 0:ntrunc), &
             v_cos_fft(nlat, 0:ntrunc), v_sin_fft(nlat, 0:ntrunc), &
             stat=alloc_status)
    ierror = 8
    if (alloc_status /= 0) return

    inv_r = 1.0d0 / dble(rsphere)
    do j = 1, nlat
        weight_inv_r(j) = dble(weights(j)) * inv_r
        weight_inv_r_sin(j) = weight_inv_r(j) / dble(sin_theta(j))
    end do
    do m = 0, ntrunc
        first_j = 1
        twom = 2 * m
        do while (first_j <= pair_count .and. twom > pl(first_j))
            first_j = first_j + 1
        end do
        first_pair_for_m(m) = first_j
        center_active(m) = half_nlat > pair_count .and. twom <= pl(half_nlat)
    end do
    outspec(:, :) = (0.0, 0.0)

    do k = 1, nt
!$omp parallel private(j, i, m, offset, row_nlon, mlimit, inv_nlon, &
!$omp& twom, fft_rows, fft_work)
        allocate(fft_rows(2, max_nlon), fft_work(2, max_nlon))
!$omp do schedule(dynamic)
        do j = 1, nlat
            offset = offsets(j)
            row_nlon = pl(j)
            mlimit = min(ntrunc, row_nlon / 2)
            inv_nlon = 1.0d0 / dble(row_nlon)
            fft_rows(:, 1:row_nlon) = 0.0
            do i = 1, row_nlon
                fft_rows(1, i) = ugrid(offset + i, k)
                fft_rows(2, i) = vgrid(offset + i, k)
            end do

            call hrfftf(2, row_nlon, fft_rows, 2, &
                reduced_gaussian_fft_wsave(1, j), fft_work)

            u_cos_fft(j, 0) = dble(fft_rows(1, 1)) * inv_nlon
            v_cos_fft(j, 0) = dble(fft_rows(2, 1)) * inv_nlon
            do m = 1, mlimit
                twom = 2 * m
                u_cos_fft(j, m) = dble(fft_rows(1, twom)) * inv_nlon
                v_cos_fft(j, m) = dble(fft_rows(2, twom)) * inv_nlon
                if (twom < row_nlon) then
                    u_sin_fft(j, m) = -dble(fft_rows(1, twom + 1)) * inv_nlon
                    v_sin_fft(j, m) = -dble(fft_rows(2, twom + 1)) * inv_nlon
                else
                    u_sin_fft(j, m) = 0.0d0
                    v_sin_fft(j, m) = 0.0d0
                end if
            end do
        end do
!$omp end do
        deallocate(fft_rows, fft_work)
!$omp end parallel

        if (symmetric_pl) then
!$omp parallel do schedule(dynamic) private(m, j, js, n, nm, row_nlon, first_j, &
!$omp& u_mean, v_mean, u_cos, u_sin, v_cos, v_sin, &
!$omp& u_mean_s, v_mean_s, u_cos_s, u_sin_s, v_cos_s, v_sin_s, &
!$omp& pval, dpval, wval, value_real, value_imag, p_scale, dp_scale, &
!$omp& p_weighted, dp_weighted, parity, p_u_cos, p_u_sin, p_v_cos, &
!$omp& p_v_sin, dp_u_cos, dp_u_sin, dp_v_cos, dp_v_sin, &
!$omp& value_real_acc, value_imag_acc)
        do m = 0, ntrunc
            value_real_acc(m:ntrunc) = 0.0d0
            value_imag_acc(m:ntrunc) = 0.0d0
            if (m == 0) then
                do j = 1, pair_count
                    js = nlat - j + 1
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)
                    u_mean_s = u_cos_fft(js, 0)
                    v_mean_s = v_cos_fft(js, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(j, nm))
                        dp_weighted = dp_scale * dpval
                        if (component == 1) then
                            value_real = -dp_weighted * &
                                (u_mean - parity * u_mean_s)
                        else
                            value_real = dp_weighted * &
                                (v_mean - parity * v_mean_s)
                        end if
                        value_real_acc(n) = value_real_acc(n) + value_real
                        parity = -parity
                    end do
                end do

                if (half_nlat > pair_count) then
                    j = half_nlat
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(j, nm))
                        dp_weighted = dp_scale * dpval
                        if (component == 1) then
                            value_real = -dp_weighted * u_mean
                        else
                            value_real = dp_weighted * v_mean
                        end if
                        value_real_acc(n) = value_real_acc(n) + value_real
                    end do
                end if
            else
                first_j = first_pair_for_m(m)
                do j = first_j, pair_count
                    row_nlon = pl(j)
                    js = nlat - j + 1
                    u_cos = u_cos_fft(j, m)
                    u_sin = u_sin_fft(j, m)
                    v_cos = v_cos_fft(j, m)
                    v_sin = v_sin_fft(j, m)
                    u_cos_s = u_cos_fft(js, m)
                    u_sin_s = u_sin_fft(js, m)
                    v_cos_s = v_cos_fft(js, m)
                    v_sin_s = v_sin_fft(js, m)

                    nm = spectral_offset(m, ntrunc)
                    p_scale = dble(m) * weight_inv_r_sin(j)
                    dp_scale = weight_inv_r(j)
                    parity = 1.0d0
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(j, nm))
                        dpval = dble(dbasis(j, nm))
                        p_weighted = p_scale * pval
                        dp_weighted = dp_scale * dpval
                        if (component == 1) then
                            p_v_sin = v_sin + parity * v_sin_s
                            p_v_cos = v_cos + parity * v_cos_s
                            dp_u_sin = u_sin - parity * u_sin_s
                            dp_u_cos = u_cos - parity * u_cos_s
                            value_real = -dp_weighted * dp_u_cos + &
                                p_weighted * p_v_sin
                            value_imag = dp_weighted * dp_u_sin + &
                                p_weighted * p_v_cos
                        else
                            p_u_sin = u_sin + parity * u_sin_s
                            p_u_cos = u_cos + parity * u_cos_s
                            dp_v_sin = v_sin - parity * v_sin_s
                            dp_v_cos = v_cos - parity * v_cos_s
                            value_real = p_weighted * p_u_sin + &
                                dp_weighted * dp_v_cos
                            value_imag = p_weighted * p_u_cos - &
                                dp_weighted * dp_v_sin
                        end if
                        value_real_acc(n) = value_real_acc(n) + value_real
                        value_imag_acc(n) = value_imag_acc(n) + value_imag
                        parity = -parity
                    end do
                end do

                if (center_active(m)) then
                    j = half_nlat
                    u_cos = u_cos_fft(j, m)
                    u_sin = u_sin_fft(j, m)
                    v_cos = v_cos_fft(j, m)
                    v_sin = v_sin_fft(j, m)

                    nm = spectral_offset(m, ntrunc)
                    p_scale = dble(m) * weight_inv_r_sin(j)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(j, nm))
                        dpval = dble(dbasis(j, nm))
                        p_weighted = p_scale * pval
                        dp_weighted = dp_scale * dpval
                        if (component == 1) then
                            value_real = -dp_weighted * u_cos + &
                                p_weighted * v_sin
                            value_imag = dp_weighted * u_sin + &
                                p_weighted * v_cos
                        else
                            value_real = p_weighted * u_sin + &
                                dp_weighted * v_cos
                            value_imag = p_weighted * u_cos - &
                                dp_weighted * v_sin
                        end if
                        value_real_acc(n) = value_real_acc(n) + value_real
                        value_imag_acc(n) = value_imag_acc(n) + value_imag
                    end do
                end if
            end if
            nm = spectral_offset(m, ntrunc)
            do n = m, ntrunc
                nm = nm + 1
                outspec(nm, k) = cmplx(real(value_real_acc(n)), &
                    real(value_imag_acc(n)))
            end do
        end do
!$omp end parallel do
        else
!$omp parallel do schedule(dynamic) private(m, j, n, nm, row_nlon, &
!$omp& u_mean, v_mean, u_cos, u_sin, v_cos, v_sin, pval, &
!$omp& dpval, wval, value_real, value_imag, p_scale, dp_scale, &
!$omp& p_weighted, dp_weighted, value_real_acc, value_imag_acc)
        do m = 0, ntrunc
            value_real_acc(m:ntrunc) = 0.0d0
            value_imag_acc(m:ntrunc) = 0.0d0
            if (m == 0) then
                do j = 1, nlat
                    u_mean = u_cos_fft(j, 0)
                    v_mean = v_cos_fft(j, 0)

                    nm = spectral_offset(m, ntrunc)
                    dp_scale = weight_inv_r(j)
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(j, nm))
                        dp_weighted = dp_scale * dpval
                        if (component == 1) then
                            value_real = -dp_weighted * u_mean
                        else
                            value_real = dp_weighted * v_mean
                        end if
                        value_real_acc(n) = value_real_acc(n) + value_real
                    end do
                end do
            else
                do j = 1, nlat
                    row_nlon = pl(j)
                    if (2 * m <= row_nlon) then
                        u_cos = u_cos_fft(j, m)
                        u_sin = u_sin_fft(j, m)
                        v_cos = v_cos_fft(j, m)
                        v_sin = v_sin_fft(j, m)

                        nm = spectral_offset(m, ntrunc)
                        p_scale = dble(m) * weight_inv_r_sin(j)
                        dp_scale = weight_inv_r(j)
                        do n = m, ntrunc
                            nm = nm + 1
                            pval = dble(basis(j, nm))
                            dpval = dble(dbasis(j, nm))
                            p_weighted = p_scale * pval
                            dp_weighted = dp_scale * dpval
                            if (component == 1) then
                                value_real = -dp_weighted * u_cos + &
                                    p_weighted * v_sin
                                value_imag = dp_weighted * u_sin + &
                                    p_weighted * v_cos
                            else
                                value_real = p_weighted * u_sin + &
                                    dp_weighted * v_cos
                                value_imag = p_weighted * u_cos - &
                                    dp_weighted * v_sin
                            end if
                            value_real_acc(n) = value_real_acc(n) + value_real
                            value_imag_acc(n) = value_imag_acc(n) + value_imag
                        end do
                    end if
                end do
            end if
            nm = spectral_offset(m, ntrunc)
            do n = m, ntrunc
                nm = nm + 1
                outspec(nm, k) = cmplx(real(value_real_acc(n)), &
                    real(value_imag_acc(n)))
            end do
        end do
!$omp end parallel do
        end if
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getvecspec_component
