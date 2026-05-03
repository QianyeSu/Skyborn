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

    real(real64), parameter :: omega = 7.292e-5_real64
    real(real64), parameter :: radius = 6.371e6_real64
    real(real64), parameter :: gas_constant_dry = 287.04_real64
    real(real64), parameter :: wavenumber_step = 1.0e-7_real64
    integer, parameter :: nwavenumbers = 200

    public :: real64
    public :: gas_constant_dry
    public :: nan_value
    public :: nwavenumbers
    public :: omega
    public :: radius
    public :: dgeev_lwork
    public :: solve_tridiagonal_system_eig_max_imag
    public :: wavenumber_step
    public :: finalize_baroc_growth
    public :: finalize_barot_growth

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

    pure real(real64) function nan_value() result(value)
        value = ieee_value(0.0_real64, ieee_quiet_nan)
    end function nan_value


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


    subroutine finalize_baroc_growth(growth, max_growth, ier)
        real(real64), intent(in) :: growth(:)
        real(real64), intent(out) :: max_growth
        integer, intent(out) :: ier

        integer :: i, last_turn
        logical :: has_valid

        last_turn = 0
        do i = 1, size(growth) - 1
            if (ieee_is_nan(growth(i)) .or. ieee_is_nan(growth(i + 1))) cycle
            if (growth(i) - growth(i + 1) > 0.0_real64) last_turn = i
        end do

        if (last_turn == 0) then
            ier = 300
            max_growth = nan_value()
            return
        end if

        has_valid = .false.
        max_growth = nan_value()
        do i = 1, last_turn + 1
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
            ier = 301
            max_growth = nan_value()
        end if
    end subroutine finalize_baroc_growth

end module growth_rate_kernels_core


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


subroutine dbaroc_growth_rate_1d(u, theta, pressure, temperature, lat, max_growth, ier, nlev)
    use growth_rate_kernels_core, only : &
        dgeev_lwork, real64, gas_constant_dry, nan_value, nwavenumbers, omega, &
        radius, solve_tridiagonal_system_eig_max_imag, wavenumber_step, finalize_baroc_growth
    implicit none

    integer, intent(in) :: nlev
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
    call finalize_baroc_growth(growth, max_growth, ier)
end subroutine dbaroc_growth_rate_1d
