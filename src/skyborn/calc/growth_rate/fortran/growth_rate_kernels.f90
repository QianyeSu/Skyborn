! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-05-03
! File    : growth_rate_kernels.f90
! Purpose : Compiled Chemke-style baroclinic and barotropic growth-rate cores.
! =============================================================================

module growth_rate_kernels_core
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
    implicit none
    private

    ! shared physical constants and numerical controls for the growth-rate
    ! kernels so the barotropic and baroclinic paths use the same sampling.
    real(real64), parameter :: omega = 7.292e-5_real64
    real(real64), parameter :: radius = 6.371e6_real64
    real(real64), parameter :: gas_constant_dry = 287.04_real64
    real(real64), parameter :: heat_capacity = 1004.7_real64
    real(real64), parameter :: kappa = gas_constant_dry / heat_capacity
    real(real64), parameter :: reference_pressure_pa = 1.0e5_real64
    real(real64), parameter :: wavenumber_step = 1.0e-7_real64
    integer, parameter :: nwavenumbers = 200

    public :: real64
    public :: gas_constant_dry
    public :: kappa
    public :: reference_pressure_pa
    public :: nan_value
    public :: nwavenumbers
    public :: omega
    public :: radius
    public :: dgeev_lwork
    public :: interp_pressure_profile_no_extrap
    public :: solve_tridiagonal_system_eig_max_imag
    public :: wavenumber_step
    public :: finalize_baroc_growth
    public :: finalize_barot_growth
    public :: apply_centered_running_mean

    ! lapack interfaces used by the dense eigenvalue and tridiagonal solve
    ! helpers below.
    interface
        subroutine dgtsv(n, nrhs, dl, d, du, b, ldb, info)
            import :: real64
            integer, intent(in) :: n, nrhs, ldb
            real(real64), intent(inout) :: dl(*), d(*), du(*), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgtsv

        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
                         work, lwork, info)
            import :: real64
            character(len=1), intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(real64), intent(inout) :: a(lda, *)
            real(real64), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine dgeev
    end interface

contains

    ! quick reference
    ! purpose
    !    return a quiet nan sentinel for failed growth-rate calculations.
    ! expected input shapes
    !    none.
    ! output
    !    value - quiet nan sentinel.
    pure real(real64) function nan_value() result(value)
        value = ieee_value(0.0_real64, ieee_quiet_nan)
    end function nan_value


    ! quick reference
    ! purpose
    !    query lapack dgeev for a safe workspace size for an n x n matrix.
    ! expected input shapes
    !    n - matrix order used to size the temporary query arrays.
    ! output
    !    lwork - recommended workspace length for dgeev.
    integer function dgeev_lwork(n) result(lwork)
        integer, intent(in) :: n

        real(real64), allocatable :: a(:,:), wr(:), wi(:)
        real(real64) :: vl(1, 1), vr(1, 1), work_query(1)
        integer :: info

        allocate(a(n, n), wr(n), wi(n))
        a = 0.0_real64
        work_query = 0.0_real64
        call dgeev('N', 'N', n, a, n, wr, wi, vl, 1, vr, 1, work_query, -1, info)
        if (info == 0) then
            lwork = max(1, int(work_query(1)))
        else
            lwork = max(1, 8 * n)
        end if
        deallocate(a, wr, wi)
    end function dgeev_lwork


    ! quick reference
    ! purpose
    !    solve the tridiagonal system, convert it to a dense eigenproblem,
    !    and return the largest imaginary eigenvalue component.
    ! expected input shapes
    !    eig_matrix(n,n) - dense matrix overwritten by dgtsv and passed to dgeev.
    !    dl(n-1)         - lower diagonal band of the tridiagonal matrix.
    !    d(n)            - main diagonal band of the tridiagonal matrix.
    !    du(n-1)         - upper diagonal band of the tridiagonal matrix.
    !    wr(n)           - real parts of the lapack eigenvalues.
    !    wi(n)           - imaginary parts of the lapack eigenvalues.
    !    work(lwork)     - lapack workspace for dgeev.
    ! output
    !    max_imag - largest imaginary eigenvalue component.
    !    info      - lapack status flag from dgtsv/dgeev.
    subroutine solve_tridiagonal_system_eig_max_imag( &
        eig_matrix, dl, d, du, wr, wi, work, lwork, max_imag, info &
    )
        real(real64), intent(inout) :: eig_matrix(:, :), dl(:), d(:), du(:), wr(:), wi(:), work(:)
        integer, intent(in) :: lwork
        real(real64), intent(out) :: max_imag
        integer, intent(out) :: info

        real(real64) :: vl(1, 1), vr(1, 1)
        integer :: i, n
        logical :: has_valid

        n = size(eig_matrix, 1)
        call dgtsv(n, n, dl, d, du, eig_matrix, n, info)
        if (info /= 0) then
            max_imag = nan_value()
            return
        end if

        call dgeev('N', 'N', n, eig_matrix, n, wr, wi, vl, 1, vr, 1, work, lwork, info)
        if (info /= 0) then
            max_imag = nan_value()
            return
        end if

        max_imag = nan_value()
        has_valid = .false.
        do i = 1, n
            if (.not. has_valid) then
                max_imag = wi(i)
                has_valid = .true.
            else
                max_imag = max(max_imag, wi(i))
            end if
        end do

        if (.not. has_valid) then
            info = 1
            max_imag = nan_value()
        end if
    end subroutine solve_tridiagonal_system_eig_max_imag


    ! quick reference
    ! purpose
    !    linearly or log-linearly interpolate one monotonic pressure profile
    !    to a new monotonic pressure grid without extrapolation.
    ! expected input shapes
    !    source_pressure(nin) - source pressure levels.
    !    source_values(nin) - field values on source_pressure.
    !    target_pressure(nout) - target pressure levels.
    ! output
    !    target_values(nout) - interpolated output values.
    !    ier - status flag: 0 on success, 1 when the source coordinate is not
    !          strictly monotonic, 2 when the target coordinate is not strictly
    !          increasing, 3 when log interpolation sees nonpositive pressure,
    !          4 when a target level falls outside the source range.
    pure subroutine interp_pressure_profile_no_extrap( &
        source_pressure, source_values, target_pressure, interp_kind, target_values, ier &
    )
        real(real64), intent(in) :: source_pressure(:), source_values(:), target_pressure(:)
        integer, intent(in) :: interp_kind
        real(real64), intent(out) :: target_values(size(target_pressure))
        integer, intent(out) :: ier

        real(real64) :: pressure_work(size(source_pressure)), values_work(size(source_values))
        real(real64) :: delta_pressure, source_hi, source_lo, target_level, value_hi, value_lo
        integer :: idx, jdx, nin, nout

        nin = size(source_pressure)
        nout = size(target_pressure)
        ier = 0
        target_values = 0.0_real64

        if (nin < 2) then
            ier = 1
            return
        end if

        if (source_pressure(1) > source_pressure(nin)) then
            do idx = 1, nin
                pressure_work(idx) = source_pressure(nin - idx + 1)
                values_work(idx) = source_values(nin - idx + 1)
            end do
        else
            pressure_work = source_pressure
            values_work = source_values
        end if

        do idx = 1, nin - 1
            if (pressure_work(idx) >= pressure_work(idx + 1)) then
                ier = 1
                return
            end if
        end do

        if (nout > 1) then
            do idx = 1, nout - 1
                if (target_pressure(idx) >= target_pressure(idx + 1)) then
                    ier = 2
                    return
                end if
            end do
        end if

        if (interp_kind == 2) then
            if (any(pressure_work <= 0.0_real64) .or. any(target_pressure <= 0.0_real64)) then
                ier = 3
                return
            end if
        end if

        jdx = 1
        do idx = 1, nout
            target_level = target_pressure(idx)
            if (target_level < pressure_work(1) .or. target_level > pressure_work(nin)) then
                ier = 4
                return
            end if

            do while (jdx < nin - 1 .and. target_level > pressure_work(jdx + 1))
                jdx = jdx + 1
            end do

            if (target_level == pressure_work(jdx)) then
                target_values(idx) = values_work(jdx)
                cycle
            end if

            if (target_level == pressure_work(jdx + 1)) then
                target_values(idx) = values_work(jdx + 1)
                cycle
            end if

            source_lo = pressure_work(jdx)
            source_hi = pressure_work(jdx + 1)
            value_lo = values_work(jdx)
            value_hi = values_work(jdx + 1)
            if (interp_kind == 1) then
                delta_pressure = source_hi - source_lo
                target_values(idx) = value_lo + &
                    ((value_hi - value_lo) / delta_pressure) * (target_level - source_lo)
            else
                delta_pressure = log(source_hi) - log(source_lo)
                target_values(idx) = value_lo + &
                    ((value_hi - value_lo) / delta_pressure) * (log(target_level) - log(source_lo))
            end if
        end do
    end subroutine interp_pressure_profile_no_extrap


    ! quick reference
    ! purpose
    !    collapse the barotropic growth-rate spectrum to its maximum finite
    !    value.
    ! expected input shapes
    !    growth(nwavenumbers) - sampled growth-rate spectrum across wavenumbers.
    ! output
    !    max_growth - maximum finite growth-rate sample.
    !    ier - status flag: 0 on success, 200 when no finite sample exists.
    subroutine finalize_barot_growth(growth, max_growth, ier)
        real(real64), intent(in) :: growth(:)
        real(real64), intent(out) :: max_growth
        integer, intent(out) :: ier

        integer :: i
        logical :: has_valid

        has_valid = .false.
        max_growth = nan_value()
        do i = 1, size(growth)
            if (ieee_is_nan(growth(i))) cycle
            if (.not. has_valid) then
                max_growth = growth(i)
                has_valid = .true.
            else
                max_growth = max(max_growth, growth(i))
            end if
        end do

        if (has_valid) then
            ier = 0
        else
            ier = 200
            max_growth = nan_value()
        end if
    end subroutine finalize_barot_growth


    ! quick reference
    ! purpose
    !    smooth a one-dimensional profile with a centered running mean that
    !    ignores nan samples.
    ! expected input shapes
    !    values(n) - raw growth-rate or diagnostic samples.
    !    window - centered window width; values <= 1 leave the series unchanged.
    ! output
    !    smoothed(n) - averaged profile.
    subroutine apply_centered_running_mean(values, window, smoothed)
        real(real64), intent(in) :: values(:)
        integer, intent(in) :: window
        real(real64), intent(out) :: smoothed(size(values))

        integer :: count_valid, half_window, i, j, j_end, j_start
        real(real64) :: value_sum

        if (window <= 1) then
            smoothed = values
            return
        end if

        half_window = window / 2
        do i = 1, size(values)
            j_start = max(1, i - half_window)
            j_end = min(size(values), i + half_window)
            value_sum = 0.0_real64
            count_valid = 0
            do j = j_start, j_end
                if (ieee_is_nan(values(j))) cycle
                value_sum = value_sum + values(j)
                count_valid = count_valid + 1
            end do

            if (count_valid > 0) then
                smoothed(i) = value_sum / real(count_valid, real64)
            else
                smoothed(i) = nan_value()
            end if
        end do
    end subroutine apply_centered_running_mean


    ! quick reference
    ! purpose
    !    smooth the baroclinic growth-rate spectrum, find the last turning
    !    point, and return the peak value up to that point.
    ! expected input shapes
    !    growth(nwavenumbers) - sampled growth-rate spectrum across wavenumbers.
    !    smooth_window - centered smoothing width applied before peak detection.
    ! output
    !    max_growth - maximum growth up to the last turning point.
    !    ier - status flag: 0 on success, 300 when no turning point is found,
    !         301 when no valid sample remains after smoothing.
    subroutine finalize_baroc_growth(growth, smooth_window, max_growth, ier)
        real(real64), intent(in) :: growth(:)
        integer, intent(in) :: smooth_window
        real(real64), intent(out) :: max_growth
        integer, intent(out) :: ier

        integer :: i, last_turn
        logical :: has_valid
        real(real64) :: growth_used(size(growth))

        call apply_centered_running_mean(growth, smooth_window, growth_used)

        last_turn = 0
        do i = 1, size(growth_used) - 1
            if (ieee_is_nan(growth_used(i)) .or. ieee_is_nan(growth_used(i + 1))) cycle
            if (growth_used(i) - growth_used(i + 1) > 0.0_real64) last_turn = i
        end do

        if (last_turn == 0) then
            ier = 300
            max_growth = nan_value()
            return
        end if

        has_valid = .false.
        max_growth = nan_value()
        do i = 1, last_turn + 1
            if (ieee_is_nan(growth_used(i))) cycle
            if (.not. has_valid) then
                max_growth = growth_used(i)
                has_valid = .true.
            else
                max_growth = max(max_growth, growth_used(i))
            end if
        end do

        if (has_valid) then
            ier = 0
        else
            ier = 301
            max_growth = nan_value()
        end if
    end subroutine finalize_baroc_growth

end module growth_rate_kernels_core


! quick reference
! purpose
!    compute the chemke-style barotropic growth-rate spectrum for a latitude
!    profile and return the maximum growth.
! expected input shapes
!    lat(nlat) - latitude coordinates in degrees for the zonal-mean profile.
!    u(nlat) - zonal wind profile on the same latitude grid.
! units
!    lat - degrees
!    u - m s-1
! output
!    max_growth - maximum finite barotropic growth rate over the sampled
!                 wavenumbers.
!    ier - status flag: 0 on success, 1 when nlat < 3, 200 when no finite
!          sample exists.
subroutine dbarot_growth_rate_1d(lat, u, max_growth, ier, nlat)
    use growth_rate_kernels_core, only : &
        dgeev_lwork, real64, nan_value, nwavenumbers, omega, radius, &
        solve_tridiagonal_system_eig_max_imag, wavenumber_step, finalize_barot_growth
    implicit none

    integer, intent(in) :: nlat
    real(real64), intent(in) :: lat(nlat), u(nlat)
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a_diag_linear(nlat), a_lower_coeff(nlat - 1), a_upper_coeff(nlat - 1)
    real(real64) :: b_diag_base(nlat), b_lower_coeff(nlat - 1), b_upper_coeff(nlat - 1), beta_arr(nlat)
    real(real64) :: dy(nlat - 1), dy1(nlat - 1), dy2(nlat - 2), growth(nwavenumbers)
    real(real64) :: kval, kval2, kval3, y(nlat), y1, y2, growth_k
    real(real64), allocatable :: eig_matrix(:, :), b_diag(:), b_lower(:), b_upper(:), wr(:), wi(:), eig_work(:)
    integer :: info, k, lwork, n

    if (nlat < 3) then
        ier = 1
        max_growth = nan_value()
        return
    end if

    y = (acos(-1.0_real64) / 180.0_real64) * lat * radius
    dy1 = (y(2:) + y(:nlat - 1)) / 2.0_real64
    dy2 = dy1(2:) - dy1(:nlat - 2)
    dy = (y(2:) - y(:nlat - 1)) * cos((acos(-1.0_real64) / 180.0_real64) * ((lat(2:) + lat(:nlat - 1)) / 2.0_real64))
    y1 = (y(2) - y(1)) * cos((acos(-1.0_real64) / 180.0_real64) * ((lat(2) + lat(1)) / 2.0_real64))
    y2 = (y(nlat) - y(nlat - 1)) * cos((acos(-1.0_real64) / 180.0_real64) * ((lat(nlat) + lat(nlat - 1)) / 2.0_real64))
    beta_arr = 2.0_real64 * omega * cos((acos(-1.0_real64) / 180.0_real64) * lat) / radius

    a_diag_linear(1) = beta_arr(1) + u(1) / y1 * (-1.0_real64 / y1 - 1.0_real64 / dy(1)) - &
        (1.0_real64 / y1) * (((-u(1)) / y1) - ((u(1) - u(2)) / dy(1)))
    b_diag_base(1) = (1.0_real64 / y1) * (-1.0_real64 / y1 - 1.0_real64 / dy(1))
    a_upper_coeff(1) = u(1) / (y1 * dy(1))
    b_upper_coeff(1) = 1.0_real64 / (y1 * dy(1))

    do n = 2, nlat - 1
        a_lower_coeff(n - 1) = u(n) / (dy2(n - 1) * dy(n - 1))
        a_upper_coeff(n) = u(n) / (dy2(n - 1) * dy(n))
        a_diag_linear(n) = beta_arr(n) + &
            u(n) / dy2(n - 1) * (-1.0_real64 / dy(n - 1) - 1.0_real64 / dy(n)) - &
            (1.0_real64 / dy2(n - 1)) * (((u(n - 1) - u(n)) / dy(n - 1)) - ((u(n) - u(n + 1)) / dy(n)))
        b_lower_coeff(n - 1) = 1.0_real64 / (dy2(n - 1) * dy(n - 1))
        b_upper_coeff(n) = 1.0_real64 / (dy2(n - 1) * dy(n))
        b_diag_base(n) = (1.0_real64 / dy2(n - 1)) * (-1.0_real64 / dy(n - 1) - 1.0_real64 / dy(n))
    end do

    a_lower_coeff(nlat - 1) = u(nlat) / (y2 * dy(nlat - 1))
    b_lower_coeff(nlat - 1) = 1.0_real64 / (y2 * dy(nlat - 1))
    a_diag_linear(nlat) = beta_arr(nlat) + u(nlat) / y2 * (-1.0_real64 / dy(nlat - 1) - 1.0_real64 / y2) - &
        (1.0_real64 / y2) * (((u(nlat - 1) - u(nlat)) / dy(nlat - 1)) - u(nlat) / y2)
    b_diag_base(nlat) = (1.0_real64 / y2) * (-1.0_real64 / dy(nlat - 1) - 1.0_real64 / y2)

    lwork = dgeev_lwork(nlat)
    allocate(eig_matrix(nlat, nlat), b_diag(nlat), b_lower(nlat - 1), b_upper(nlat - 1), wr(nlat), wi(nlat), eig_work(lwork))

    growth = nan_value()
    do k = 1, nwavenumbers
        kval = real(k - 1, real64) * wavenumber_step
        kval2 = kval * kval
        kval3 = kval2 * kval
        eig_matrix = 0.0_real64
        b_lower = b_lower_coeff
        b_diag = b_diag_base - kval2
        b_upper = b_upper_coeff

        do n = 1, nlat
            eig_matrix(n, n) = kval * a_diag_linear(n) - kval3 * u(n)
        end do
        do n = 1, nlat - 1
            eig_matrix(n + 1, n) = kval * a_lower_coeff(n)
            eig_matrix(n, n + 1) = kval * a_upper_coeff(n)
        end do

        call solve_tridiagonal_system_eig_max_imag( &
            eig_matrix, b_lower, b_diag, b_upper, wr, wi, eig_work, lwork, growth_k, info &
        )
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(eig_matrix, b_diag, b_lower, b_upper, wr, wi, eig_work)
    call finalize_barot_growth(growth, max_growth, ier)
end subroutine dbarot_growth_rate_1d


! quick reference
! purpose
!    compute the chemke-style baroclinic growth-rate spectrum for one vertical
!    profile and return the smoothed peak growth.
! expected input shapes
!    u(nlev) - zonal wind profile with height.
!    theta(nlev) - potential temperature profile on the same vertical grid.
!    pressure(nlev) - pressure profile used to build the vertical spacing.
!    temperature(nlev) - temperature profile used to compute density.
!    lat - latitude in degrees for the coriolis terms.
!    smooth_window - centered smoothing width applied before peak detection.
! units
!    lat - degrees
!    pressure - pa
!    temperature - k
!    u - m s-1
!    theta - k
! output
!    max_growth - maximum baroclinic growth rate up to the last turning point.
!    ier - status flag: 0 on success, 1 when nlev < 3, 2 when smooth_window < 1,
!          300 when no turning point is found, 301 when no valid sample remains
!          after smoothing.
subroutine dbaroc_growth_rate_1d(u, theta, pressure, temperature, lat, smooth_window, max_growth, ier, nlev)
    use growth_rate_kernels_core, only : &
        dgeev_lwork, real64, gas_constant_dry, nan_value, nwavenumbers, omega, &
        radius, solve_tridiagonal_system_eig_max_imag, wavenumber_step, finalize_baroc_growth
    implicit none

    integer, intent(in) :: nlev
    integer, intent(in) :: smooth_window
    real(real64), intent(in) :: u(nlev), theta(nlev), pressure(nlev), temperature(nlev), lat
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a_diag_linear(nlev), a_lower_coeff(nlev - 1), a_upper_coeff(nlev - 1)
    real(real64) :: b_diag_base(nlev), b_lower_coeff(nlev - 1), b_upper_coeff(nlev - 1), dz1(nlev - 1), dz2(nlev)
    real(real64) :: growth(nwavenumbers), growth_k, kval
    real(real64) :: f, f2, beta, kval2, kval3, nz(nlev - 1), rho(nlev)
    real(real64), allocatable :: eig_matrix(:, :), b_diag(:), b_lower(:), b_upper(:), wr(:), wi(:), eig_work(:)
    integer :: info, k, lwork, n

    if (nlev < 3) then
        ier = 1
        max_growth = nan_value()
        return
    end if

    if (smooth_window < 1) then
        ier = 2
        max_growth = nan_value()
        return
    end if

    dz1 = (pressure(2:) + pressure(:nlev - 1)) / 2.0_real64
    dz2(1) = pressure(2) - pressure(1)
    dz2(2:nlev - 1) = dz1(2:) - dz1(:nlev - 2)
    dz2(nlev) = pressure(nlev) - pressure(nlev - 1)
    f = 2.0_real64 * omega * sin((acos(-1.0_real64) / 180.0_real64) * lat)
    f2 = f * f
    beta = 2.0_real64 * omega * cos((acos(-1.0_real64) / 180.0_real64) * lat) / radius
    rho = pressure / (gas_constant_dry * temperature)
    nz = -((theta(2:) - theta(:nlev - 1)) * (1.0_real64 / (rho(:nlev - 1) * theta(:nlev - 1))))

    a_diag_linear(1) = beta + u(1) * f2 / dz2(1) * (-1.0_real64 / nz(1)) - &
        f2 / dz2(1) * ((u(2) - u(1)) / nz(1))
    b_diag_base(1) = f2 / dz2(1) * (-1.0_real64 / nz(1))
    a_upper_coeff(1) = u(1) * f2 / (dz2(1) * nz(1))
    b_upper_coeff(1) = f2 / (dz2(1) * nz(1))

    do n = 2, nlev - 1
        a_lower_coeff(n - 1) = u(n) * f2 / (dz2(n) * nz(n - 1))
        a_upper_coeff(n) = u(n) * f2 / (dz2(n) * nz(n))
        a_diag_linear(n) = beta + u(n) * f2 / dz2(n) * (-1.0_real64 / nz(n - 1) - 1.0_real64 / nz(n)) - &
            f2 / dz2(n) * (((u(n - 1) - u(n)) / nz(n - 1)) - ((u(n) - u(n + 1)) / nz(n)))
        b_lower_coeff(n - 1) = f2 / (dz2(n) * nz(n - 1))
        b_upper_coeff(n) = f2 / (dz2(n) * nz(n))
        b_diag_base(n) = f2 / dz2(n) * (-1.0_real64 / nz(n - 1) - 1.0_real64 / nz(n))
    end do

    a_lower_coeff(nlev - 1) = u(nlev) * f2 / (dz2(nlev) * nz(nlev - 1))
    b_lower_coeff(nlev - 1) = f2 / (dz2(nlev) * nz(nlev - 1))
    a_diag_linear(nlev) = beta + u(nlev) * f2 / dz2(nlev) * (-1.0_real64 / nz(nlev - 1)) - &
        f2 / dz2(nlev) * ((u(nlev - 1) - u(nlev)) / nz(nlev - 1))
    b_diag_base(nlev) = f2 / dz2(nlev) * (-1.0_real64 / nz(nlev - 1))

    lwork = dgeev_lwork(nlev)
    allocate(eig_matrix(nlev, nlev), b_diag(nlev), b_lower(nlev - 1), b_upper(nlev - 1), wr(nlev), wi(nlev), eig_work(lwork))

    growth = nan_value()
    do k = 1, nwavenumbers
        kval = real(k - 1, real64) * wavenumber_step
        kval2 = kval * kval
        kval3 = kval2 * kval
        eig_matrix = 0.0_real64
        b_lower = b_lower_coeff
        b_diag = b_diag_base - kval2
        b_upper = b_upper_coeff

        do n = 1, nlev
            eig_matrix(n, n) = kval * a_diag_linear(n) - kval3 * u(n)
        end do
        do n = 1, nlev - 1
            eig_matrix(n + 1, n) = kval * a_lower_coeff(n)
            eig_matrix(n, n + 1) = kval * a_upper_coeff(n)
        end do

        call solve_tridiagonal_system_eig_max_imag( &
            eig_matrix, b_lower, b_diag, b_upper, wr, wi, eig_work, lwork, growth_k, info &
        )
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(eig_matrix, b_diag, b_lower, b_upper, wr, wi, eig_work)
    call finalize_baroc_growth(growth, smooth_window, max_growth, ier)
end subroutine dbaroc_growth_rate_1d


! quick reference
! purpose
!    compute chemke-style baroclinic growth rates for multiple independent
!    atmospheric profiles without looping in Python.
! expected input shapes
!    u_input(nlev_in,nprofile) - zonal-wind profiles.
!    temperature_input(nlev_in,nprofile) - temperature profiles.
!    source_pressure(nlev_in) - shared source pressure coordinate.
!    target_pressure(nlev_out,nprofile) - per-profile solver pressure grids.
!    lat(nprofile) - latitude per profile.
! units
!    lat - degrees
!    source_pressure - pa
!    target_pressure - pa
!    temperature_input - k
!    u_input - m s-1
! output
!    growth(nprofile) - maximum baroclinic growth rate for each profile.
!    ier(nprofile) - profile status flag: 0 on success, 100 when interpolation
!                    cannot be performed without extrapolation, otherwise the
!                    dbaroc_growth_rate_1d status code.
subroutine dbaroc_growth_rate_profiles( &
    u_input, temperature_input, source_pressure, target_pressure, lat, &
    interp_kind, smooth_window, missing_value, growth, ier, nlev_in, nprofile, nlev_out &
)
    use growth_rate_kernels_core, only : &
        real64, kappa, reference_pressure_pa, interp_pressure_profile_no_extrap
    implicit none

    integer, intent(in) :: nlev_in, nprofile, nlev_out
    integer, intent(in) :: interp_kind, smooth_window
    real(real64), intent(in) :: u_input(nlev_in, nprofile), temperature_input(nlev_in, nprofile)
    real(real64), intent(in) :: source_pressure(nlev_in), target_pressure(nlev_out, nprofile), lat(nprofile)
    real(real64), intent(in) :: missing_value
    real(real64), intent(out) :: growth(nprofile)
    integer, intent(out) :: ier(nprofile)

    real(real64) :: temperature_solver(nlev_out), theta_solver(nlev_out), u_solver(nlev_out)
    integer :: col, interp_ier, profile_ier

    growth = missing_value
    ier = 0
    do col = 1, nprofile
        call interp_pressure_profile_no_extrap( &
            source_pressure, u_input(:, col), target_pressure(:, col), interp_kind, u_solver, interp_ier &
        )
        if (interp_ier /= 0) then
            ier(col) = 100
            cycle
        end if

        call interp_pressure_profile_no_extrap( &
            source_pressure, temperature_input(:, col), target_pressure(:, col), interp_kind, temperature_solver, interp_ier &
        )
        if (interp_ier /= 0) then
            ier(col) = 100
            cycle
        end if

        theta_solver = temperature_solver * (reference_pressure_pa / target_pressure(:, col)) ** kappa
        call dbaroc_growth_rate_1d( &
            u_solver, theta_solver, target_pressure(:, col), temperature_solver, &
            lat(col), smooth_window, growth(col), profile_ier, nlev_out &
        )
        ier(col) = profile_ier
    end do
end subroutine dbaroc_growth_rate_profiles
