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
    public :: solve_linear_system_eig_max_imag
    public :: wavenumber_step
    public :: finalize_baroc_growth
    public :: finalize_barot_growth

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: real64
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            real(real64), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv

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


    subroutine solve_linear_system_eig_max_imag( &
        a_in, b_in, eig_matrix, b_work, ipiv, wr, wi, work, lwork, max_imag, info &
    )
        real(real64), intent(in) :: a_in(:, :), b_in(:, :)
        real(real64), intent(inout) :: eig_matrix(:, :), b_work(:, :), wr(:), wi(:), work(:)
        integer, intent(inout) :: ipiv(:)
        integer, intent(in) :: lwork
        real(real64), intent(out) :: max_imag
        integer, intent(out) :: info

        real(real64) :: vl(1, 1), vr(1, 1)
        integer :: i, n
        logical :: has_valid

        n = size(a_in, 1)
        eig_matrix = a_in
        b_work = b_in

        call dgesv(n, n, b_work, n, ipiv, eig_matrix, n, info)
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
    end subroutine solve_linear_system_eig_max_imag


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
        solve_linear_system_eig_max_imag, wavenumber_step, finalize_barot_growth
    implicit none

    integer, intent(in) :: nlat
    real(real64), intent(in) :: lat(nlat), u(nlat)
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a(nlat, nlat), b(nlat, nlat), beta_arr(nlat)
    real(real64) :: dy(nlat - 1), dy1(nlat - 1), dy2(nlat - 2), growth(nwavenumbers)
    real(real64) :: kval, y(nlat), y1, y2, growth_k
    real(real64), allocatable :: eig_matrix(:, :), b_work(:, :), wr(:), wi(:), eig_work(:)
    integer, allocatable :: ipiv(:)
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

    lwork = dgeev_lwork(nlat)
    allocate(eig_matrix(nlat, nlat), b_work(nlat, nlat), wr(nlat), wi(nlat), eig_work(lwork), ipiv(nlat))

    growth = nan_value()
    do k = 1, nwavenumbers
        kval = real(k - 1, real64) * wavenumber_step
        a = 0.0_real64
        b = 0.0_real64

        do n = 1, nlat
            if (n == nlat) then
                a(n, n) = kval * beta_arr(n) - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) / y2 * (-1.0_real64 / dy(n - 1) - 1.0_real64 / y2) - &
                    kval / y2 * (((u(n - 1) - u(n)) / dy(n - 1)) - u(n) / y2)
                a(n, n - 1) = (kval * u(n)) / (y2 * dy(n - 1))
                b(n, n) = -(kval**2) + (1.0_real64 / y2) * (-1.0_real64 / dy(n - 1) - 1.0_real64 / y2)
                b(n, n - 1) = 1.0_real64 / (y2 * dy(n - 1))
            else if (n == 1) then
                a(n, n) = kval * beta_arr(n) - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) / y1 * (-1.0_real64 / y1 - 1.0_real64 / dy(n)) - &
                    kval / y1 * (((-u(n)) / y1) - ((u(n) - u(n + 1)) / dy(n)))
                a(n, n + 1) = (kval * u(n)) / (y1 * dy(n))
                b(n, n) = -(kval**2) + (1.0_real64 / y1) * (-1.0_real64 / y1 - 1.0_real64 / dy(n))
                b(n, n + 1) = 1.0_real64 / (y1 * dy(n))
            else
                a(n, n - 1) = (kval * u(n)) / (dy2(n - 1) * dy(n - 1))
                a(n, n + 1) = (kval * u(n)) / (dy2(n - 1) * dy(n))
                a(n, n) = kval * beta_arr(n) - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) / dy2(n - 1) * (-1.0_real64 / dy(n - 1) - 1.0_real64 / dy(n)) - &
                    kval / dy2(n - 1) * (((u(n - 1) - u(n)) / dy(n - 1)) - ((u(n) - u(n + 1)) / dy(n)))
                b(n, n - 1) = 1.0_real64 / (dy2(n - 1) * dy(n - 1))
                b(n, n + 1) = 1.0_real64 / (dy2(n - 1) * dy(n))
                b(n, n) = -(kval**2) + (1.0_real64 / dy2(n - 1)) * (-1.0_real64 / dy(n - 1) - 1.0_real64 / dy(n))
            end if
        end do

        call solve_linear_system_eig_max_imag( &
            a, b, eig_matrix, b_work, ipiv, wr, wi, eig_work, lwork, growth_k, info &
        )
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(eig_matrix, b_work, wr, wi, eig_work, ipiv)
    call finalize_barot_growth(growth, max_growth, ier)
end subroutine dbarot_growth_rate_1d


subroutine dbaroc_growth_rate_1d(u, theta, pressure, temperature, lat, max_growth, ier, nlev)
    use growth_rate_kernels_core, only : &
        dgeev_lwork, real64, gas_constant_dry, nan_value, nwavenumbers, omega, &
        radius, solve_linear_system_eig_max_imag, wavenumber_step, finalize_baroc_growth
    implicit none

    integer, intent(in) :: nlev
    real(real64), intent(in) :: u(nlev), theta(nlev), pressure(nlev), temperature(nlev), lat
    real(real64), intent(out) :: max_growth
    integer, intent(out) :: ier

    real(real64) :: a(nlev, nlev), b(nlev, nlev), dz1(nlev - 1), dz2(nlev)
    real(real64) :: growth(nwavenumbers), growth_k, kval
    real(real64) :: f, beta, nz(nlev - 1), rho(nlev)
    real(real64), allocatable :: eig_matrix(:, :), b_work(:, :), wr(:), wi(:), eig_work(:)
    integer, allocatable :: ipiv(:)
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
    beta = 2.0_real64 * omega * cos((acos(-1.0_real64) / 180.0_real64) * lat) / radius
    rho = pressure / (gas_constant_dry * temperature)
    nz = -((theta(2:) - theta(:nlev - 1)) * (1.0_real64 / (rho(:nlev - 1) * theta(:nlev - 1))))

    lwork = dgeev_lwork(nlev)
    allocate(eig_matrix(nlev, nlev), b_work(nlev, nlev), wr(nlev), wi(nlev), eig_work(lwork), ipiv(nlev))

    growth = nan_value()
    do k = 1, nwavenumbers
        kval = real(k - 1, real64) * wavenumber_step
        a = 0.0_real64
        b = 0.0_real64

        do n = 1, nlev
            if (n == nlev) then
                a(n, n) = kval * beta - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) * (f**2) / dz2(n) * (-1.0_real64 / nz(n - 1)) - &
                    kval * (f**2) / dz2(n) * ((u(n - 1) - u(n)) / nz(n - 1))
                a(n, n - 1) = (kval * u(n)) * (f**2) / (dz2(n) * nz(n - 1))
                b(n, n) = -(kval**2) + (f**2) / dz2(n) * (-1.0_real64 / nz(n - 1))
                b(n, n - 1) = (f**2) / (dz2(n) * nz(n - 1))
            else if (n == 1) then
                a(n, n) = kval * beta - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) * (f**2) / dz2(n) * (-1.0_real64 / nz(n)) - &
                    kval * (f**2) / dz2(n) * ((u(n + 1) - u(n)) / nz(n))
                a(n, n + 1) = (kval * u(n)) * (f**2) / (dz2(n) * nz(n))
                b(n, n) = -(kval**2) + (f**2) / dz2(n) * (-1.0_real64 / nz(n))
                b(n, n + 1) = (f**2) / (dz2(n) * nz(n))
            else
                a(n, n - 1) = (kval * u(n)) * (f**2) / (dz2(n) * nz(n - 1))
                a(n, n + 1) = (kval * u(n)) * (f**2) / (dz2(n) * nz(n))
                a(n, n) = kval * beta - (kval**2) * (kval * u(n)) + &
                    (kval * u(n)) * (f**2) / dz2(n) * (-1.0_real64 / nz(n - 1) - 1.0_real64 / nz(n)) - &
                    kval * (f**2) / dz2(n) * (((u(n - 1) - u(n)) / nz(n - 1)) - ((u(n) - u(n + 1)) / nz(n)))
                b(n, n - 1) = (f**2) / (dz2(n) * nz(n - 1))
                b(n, n + 1) = (f**2) / (dz2(n) * nz(n))
                b(n, n) = -(kval**2) + (f**2) / dz2(n) * (-1.0_real64 / nz(n - 1) - 1.0_real64 / nz(n))
            end if
        end do

        call solve_linear_system_eig_max_imag( &
            a, b, eig_matrix, b_work, ipiv, wr, wi, eig_work, lwork, growth_k, info &
        )
        if (info == 0) growth(k) = growth_k
    end do

    deallocate(eig_matrix, b_work, wr, wi, eig_work, ipiv)
    call finalize_baroc_growth(growth, max_growth, ier)
end subroutine dbaroc_growth_rate_1d
