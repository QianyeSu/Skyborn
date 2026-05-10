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
    public :: interp_pressure_profile_no_extrap
    public :: generalized_eig_max_imag
    public :: wavenumber_step
    public :: finalize_baroc_growth
    public :: finalize_barot_growth
    public :: apply_centered_running_mean

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
    !    run a small internal QZ generalized-eigenvalue solve and return the
    !    largest finite imaginary component of lambda = (alfr + i alfi) / beta.
    ! expected input shapes
    !    a_matrix(n,n) - dense generalized-eigenvalue left matrix.
    !    b_matrix(n,n) - dense generalized-eigenvalue right matrix.
    !    alfr(n), alfi(n), beta(n) - real-valued eigenvalue work arrays.
    ! output
    !    max_imag - largest finite imaginary eigenvalue component.
    !    info     - zero on success, otherwise the QZ iteration status.
    subroutine generalized_eig_max_imag(a_matrix, b_matrix, alfr, alfi, beta, max_imag, info)
        real(real64), intent(inout) :: a_matrix(:, :), b_matrix(:, :)
        real(real64), intent(out) :: alfr(:), alfi(:), beta(:), max_imag
        integer, intent(out) :: info

        integer :: i, n
        logical :: has_valid
        real(real64) :: imag_part
        external :: qzhes, qzit, qzval

        n = size(a_matrix, 1)

        call qzhes(n, n, a_matrix, b_matrix)
        call qzit(n, n, a_matrix, b_matrix, 0.0_real64, info)
        if (info /= 0) then
            max_imag = nan_value()
            return
        end if

        call qzval(n, n, a_matrix, b_matrix, alfr, alfi, beta)

        max_imag = nan_value()
        has_valid = .false.
        do i = 1, n
            if (abs(beta(i)) <= tiny(1.0_real64)) cycle
            imag_part = alfi(i) / beta(i)
            if (.not. has_valid) then
                max_imag = imag_part
                has_valid = .true.
            else
                max_imag = max(max_imag, imag_part)
            end if
        end do

        if (.not. has_valid) then
            info = 1
            max_imag = nan_value()
        else
            info = 0
        end if
    end subroutine generalized_eig_max_imag


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
subroutine dbarot_growth_rate_1d_impl(lat, u, max_growth, ier, nlat)
    use growth_rate_kernels_core, only : &
        real64, nan_value, nwavenumbers, omega, radius, generalized_eig_max_imag, &
        wavenumber_step, finalize_barot_growth
    implicit none

    integer, intent(in) :: nlat
    real(real64), intent(in) :: lat(nlat), u(nlat)
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a_diag_linear(nlat), a_lower_coeff(nlat - 1), a_upper_coeff(nlat - 1)
    real(real64) :: b_diag_base(nlat), b_lower_coeff(nlat - 1), b_upper_coeff(nlat - 1), beta_arr(nlat)
    real(real64) :: dy(nlat - 1), dy1(nlat - 1), dy2(nlat - 2), growth(nwavenumbers)
    real(real64) :: kval, kval2, kval3, y(nlat), y1, y2, growth_k
    real(real64), allocatable :: a_matrix(:, :), b_matrix(:, :)
    real(real64), allocatable :: alfr(:), alfi(:), beta_eig(:)
    integer :: info, k, n

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

    allocate( &
        a_matrix(nlat, nlat), b_matrix(nlat, nlat), &
        alfr(nlat), alfi(nlat), beta_eig(nlat) &
    )

    growth = nan_value()
    do k = 1, nwavenumbers
        kval = real(k - 1, real64) * wavenumber_step
        kval2 = kval * kval
        kval3 = kval2 * kval
        a_matrix = 0.0_real64
        b_matrix = 0.0_real64

        do n = 1, nlat
            a_matrix(n, n) = kval * a_diag_linear(n) - kval3 * u(n)
            b_matrix(n, n) = b_diag_base(n) - kval2
        end do
        do n = 1, nlat - 1
            a_matrix(n + 1, n) = kval * a_lower_coeff(n)
            a_matrix(n, n + 1) = kval * a_upper_coeff(n)
            b_matrix(n + 1, n) = b_lower_coeff(n)
            b_matrix(n, n + 1) = b_upper_coeff(n)
        end do

        call generalized_eig_max_imag(a_matrix, b_matrix, alfr, alfi, beta_eig, growth_k, info)
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(a_matrix, b_matrix, alfr, alfi, beta_eig)
    call finalize_barot_growth(growth, max_growth, ier)
end subroutine dbarot_growth_rate_1d_impl


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
!    wavenumber_mode - 1 for the fixed high-resolution grid, 2 for the
!                      lower-resolution longitude-derived grid.
!    wavenumber_count - number of zonal wavenumber samples.
!    zonal_length - zonal domain length in meters when wavenumber_mode = 2.
! units
!    f_cor - s-1
!    beta - m-1 s-1
!    pressure - pa
!    temperature - k
!    u - m s-1
!    theta - k
! output
!    max_growth - maximum baroclinic growth rate up to the last turning point.
!    ier - status flag: 0 on success, 1 when nlev < 3, 2 when smooth_window < 1,
!          3 when wavenumber_count < 2, 4 when wavenumber_mode is invalid,
!          5 when a low-resolution zonal length is not positive, 300 when no
!          turning point is found, 301 when no valid sample remains after smoothing.
subroutine dbaroc_growth_rate_1d_impl( &
    u, theta, pressure, temperature, f_cor, beta, smooth_window, &
    wavenumber_mode, wavenumber_count, zonal_length, max_growth, ier, nlev &
)
    use growth_rate_kernels_core, only : &
        real64, gas_constant_dry, nan_value, generalized_eig_max_imag, &
        wavenumber_step, finalize_baroc_growth
    implicit none

    integer, intent(in) :: nlev
    integer, intent(in) :: smooth_window, wavenumber_mode, wavenumber_count
    real(real64), intent(in) :: u(nlev), theta(nlev), pressure(nlev), temperature(nlev), f_cor, beta
    real(real64), intent(in) :: zonal_length
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a_diag_linear(nlev), a_lower_coeff(nlev - 1), a_upper_coeff(nlev - 1)
    real(real64) :: b_diag_base(nlev), b_lower_coeff(nlev - 1), b_upper_coeff(nlev - 1), dz1(nlev - 1), dz2(nlev)
    real(real64) :: growth_k, kval, wavenumber_factor
    real(real64) :: f, f2, kval2, kval3, nz(nlev - 1), rho(nlev)
    real(real64), allocatable :: a_matrix(:, :), b_matrix(:, :), growth(:)
    real(real64), allocatable :: alfr(:), alfi(:), beta_eig(:)
    integer :: info, k, n

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

    if (wavenumber_count < 2) then
        ier = 3
        max_growth = nan_value()
        return
    end if

    if (wavenumber_mode /= 1 .and. wavenumber_mode /= 2) then
        ier = 4
        max_growth = nan_value()
        return
    end if

    if (wavenumber_mode == 2 .and. zonal_length <= 0.0_real64) then
        ier = 5
        max_growth = nan_value()
        return
    end if

    dz1 = (pressure(2:) + pressure(:nlev - 1)) / 2.0_real64
    dz2(1) = pressure(2) - pressure(1)
    dz2(2:nlev - 1) = dz1(2:) - dz1(:nlev - 2)
    dz2(nlev) = pressure(nlev) - pressure(nlev - 1)
    f = f_cor
    f2 = f * f
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

    allocate( &
        a_matrix(nlev, nlev), b_matrix(nlev, nlev), &
        alfr(nlev), alfi(nlev), beta_eig(nlev), growth(wavenumber_count) &
    )

    if (wavenumber_mode == 1) then
        wavenumber_factor = wavenumber_step
    else
        wavenumber_factor = 2.0_real64 * acos(-1.0_real64) / zonal_length
    end if

    growth = nan_value()
    do k = 1, wavenumber_count
        kval = real(k - 1, real64) * wavenumber_factor
        kval2 = kval * kval
        kval3 = kval2 * kval
        a_matrix = 0.0_real64
        b_matrix = 0.0_real64

        do n = 1, nlev
            a_matrix(n, n) = kval * a_diag_linear(n) - kval3 * u(n)
            b_matrix(n, n) = b_diag_base(n) - kval2
        end do
        do n = 1, nlev - 1
            a_matrix(n + 1, n) = kval * a_lower_coeff(n)
            a_matrix(n, n + 1) = kval * a_upper_coeff(n)
            b_matrix(n + 1, n) = b_lower_coeff(n)
            b_matrix(n, n + 1) = b_upper_coeff(n)
        end do

        call generalized_eig_max_imag(a_matrix, b_matrix, alfr, alfi, beta_eig, growth_k, info)
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(a_matrix, b_matrix, alfr, alfi, beta_eig)
    call finalize_baroc_growth(growth, smooth_window, max_growth, ier)
    deallocate(growth)
end subroutine dbaroc_growth_rate_1d_impl


! quick reference
! purpose
!    compute chemke-style baroclinic growth rates for multiple independent
!    atmospheric profiles without looping in Python.
! expected input shapes
!    u_input(nlev_in,nprofile) - zonal-wind profiles.
!    temperature_input(nlev_in,nprofile) - temperature profiles.
!    source_pressure(nlev_in) - shared source pressure coordinate.
!    target_pressure(nlev_out,nprofile) - per-profile solver pressure grids.
!    f_cor(nprofile) - coriolis parameter per profile.
!    beta(nprofile) - planetary-vorticity gradient per profile.
! units
!    f_cor - s-1
!    beta - m-1 s-1
!    source_pressure - pa
!    target_pressure - pa
!    temperature_input - k
!    u_input - m s-1
! output
!    growth(nprofile) - maximum baroclinic growth rate for each profile.
!    ier(nprofile) - profile status flag: 0 on success, 100 when interpolation
!                    cannot be performed without extrapolation, otherwise the
!                    dbaroc_growth_rate_1d status code.
subroutine dbaroc_growth_rate_profiles_impl( &
    u_input, temperature_input, source_pressure, target_pressure, f_cor, beta, &
    interp_kind, smooth_window, wavenumber_mode, wavenumber_count, zonal_length, &
    growth, ier, nlev_in, nprofile, nlev_out &
)
    use growth_rate_kernels_core, only : &
        real64, kappa, reference_pressure_pa, interp_pressure_profile_no_extrap, nan_value
    implicit none

    integer, intent(in) :: nlev_in, nprofile, nlev_out
    integer, intent(in) :: interp_kind, smooth_window, wavenumber_mode, wavenumber_count
    real(real64), intent(in) :: u_input(nlev_in, nprofile), temperature_input(nlev_in, nprofile)
    real(real64), intent(in) :: source_pressure(nlev_in), target_pressure(nlev_out, nprofile)
    real(real64), intent(in) :: f_cor(nprofile), beta(nprofile), zonal_length(nprofile)
    real(real64), intent(out) :: growth(nprofile)
    integer, intent(out) :: ier(nprofile)

    real(real64) :: temperature_solver(nlev_out), theta_solver(nlev_out), u_solver(nlev_out)
    integer :: col, interp_ier, profile_ier

    growth = nan_value()
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
        call dbaroc_growth_rate_1d_impl( &
            u_solver, theta_solver, target_pressure(:, col), temperature_solver, &
            f_cor(col), beta(col), smooth_window, wavenumber_mode, wavenumber_count, &
            zonal_length(col), growth(col), profile_ier, nlev_out &
        )
        ier(col) = profile_ier
    end do
end subroutine dbaroc_growth_rate_profiles_impl


! quick reference
! purpose
!    expose the barotropic growth-rate kernel through a C-interoperable
!    wrapper while keeping the scientific implementation in Fortran.
subroutine dbarot_growth_rate_1d(lat, u, max_growth, ier, nlat) bind(C, name="dbarot_growth_rate_1d")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: lat, u
    integer(c_int), value, intent(in) :: nlat
    real(c_double), intent(out) :: max_growth
    integer(c_int), intent(out) :: ier

    integer :: nlat_f
    real(real64), pointer :: lat_view(:), u_view(:)
    real(real64) :: max_growth_impl
    integer :: ier_impl

    nlat_f = int(nlat)
    call c_f_pointer(lat, lat_view, [nlat_f])
    call c_f_pointer(u, u_view, [nlat_f])

    call dbarot_growth_rate_1d_impl(lat_view, u_view, max_growth_impl, ier_impl, nlat_f)
    max_growth = max_growth_impl
    ier = ier_impl
end subroutine dbarot_growth_rate_1d


! quick reference
! purpose
!    expose the baroclinic single-profile kernel through a C-interoperable
!    wrapper while keeping the scientific implementation in Fortran.
subroutine dbaroc_growth_rate_1d( &
    u, theta, pressure, temperature, f_cor, beta, smooth_window, &
    wavenumber_mode, wavenumber_count, zonal_length, max_growth, ier, nlev &
) bind(C, name="dbaroc_growth_rate_1d")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: u, theta, pressure, temperature
    real(c_double), value, intent(in) :: f_cor, beta, zonal_length
    integer(c_int), value, intent(in) :: smooth_window, wavenumber_mode, wavenumber_count, nlev
    real(c_double), intent(out) :: max_growth
    integer(c_int), intent(out) :: ier

    integer :: nlev_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f
    real(real64), pointer :: u_view(:), theta_view(:), pressure_view(:), temperature_view(:)
    real(real64) :: max_growth_impl
    integer :: ier_impl

    nlev_f = int(nlev)
    smooth_window_f = int(smooth_window)
    wavenumber_mode_f = int(wavenumber_mode)
    wavenumber_count_f = int(wavenumber_count)

    call c_f_pointer(u, u_view, [nlev_f])
    call c_f_pointer(theta, theta_view, [nlev_f])
    call c_f_pointer(pressure, pressure_view, [nlev_f])
    call c_f_pointer(temperature, temperature_view, [nlev_f])

    call dbaroc_growth_rate_1d_impl( &
        u_view, theta_view, pressure_view, temperature_view, &
        f_cor, beta, smooth_window_f, wavenumber_mode_f, wavenumber_count_f, zonal_length, &
        max_growth_impl, ier_impl, nlev_f &
    )
    max_growth = max_growth_impl
    ier = ier_impl
end subroutine dbaroc_growth_rate_1d


! quick reference
! purpose
!    expose the batched baroclinic kernel through a C-interoperable wrapper
!    while keeping the scientific implementation in Fortran.
subroutine dbaroc_growth_rate_profiles( &
    u_input, temperature_input, source_pressure, target_pressure, f_cor, beta, &
    interp_kind, smooth_window, wavenumber_mode, wavenumber_count, zonal_length, &
    growth, ier, nlev_in, nprofile, nlev_out &
) bind(C, name="dbaroc_growth_rate_profiles")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: u_input, temperature_input, source_pressure, target_pressure
    type(c_ptr), value, intent(in) :: f_cor, beta, zonal_length, growth, ier
    integer(c_int), value, intent(in) :: interp_kind, smooth_window, wavenumber_mode, wavenumber_count
    integer(c_int), value, intent(in) :: nlev_in, nprofile, nlev_out

    integer :: nlev_in_f, nprofile_f, nlev_out_f
    integer :: interp_kind_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f
    real(real64), pointer :: u_view(:, :), temperature_view(:, :), source_pressure_view(:)
    real(real64), pointer :: target_pressure_view(:, :), f_cor_view(:), beta_view(:), zonal_length_view(:)
    real(real64), pointer :: growth_view(:)
    integer(c_int), pointer :: ier_view(:)

    nlev_in_f = int(nlev_in)
    nprofile_f = int(nprofile)
    nlev_out_f = int(nlev_out)
    interp_kind_f = int(interp_kind)
    smooth_window_f = int(smooth_window)
    wavenumber_mode_f = int(wavenumber_mode)
    wavenumber_count_f = int(wavenumber_count)

    call c_f_pointer(u_input, u_view, [nlev_in_f, nprofile_f])
    call c_f_pointer(temperature_input, temperature_view, [nlev_in_f, nprofile_f])
    call c_f_pointer(source_pressure, source_pressure_view, [nlev_in_f])
    call c_f_pointer(target_pressure, target_pressure_view, [nlev_out_f, nprofile_f])
    call c_f_pointer(f_cor, f_cor_view, [nprofile_f])
    call c_f_pointer(beta, beta_view, [nprofile_f])
    call c_f_pointer(zonal_length, zonal_length_view, [nprofile_f])
    call c_f_pointer(growth, growth_view, [nprofile_f])
    call c_f_pointer(ier, ier_view, [nprofile_f])
    call dbaroc_growth_rate_profiles_impl( &
        u_view, temperature_view, source_pressure_view, target_pressure_view, f_cor_view, beta_view, &
        interp_kind_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f, zonal_length_view, &
        growth_view, ier_view, nlev_in_f, nprofile_f, nlev_out_f &
    )
end subroutine dbaroc_growth_rate_profiles
