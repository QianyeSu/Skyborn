! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-05-08
! File    : isentropic_kernels.f90
! Purpose : C-order kernel for interpolation to isentropic surfaces. The public
!           Python wrapper handles xarray metadata while this file keeps the hot
!           vertical profile loops in Fortran.
! =============================================================================
!
module isentropic_kernels_core
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: real64
    public :: interp_isentropic_flat_profile

contains

    pure logical function valid_value(value, spvl) result(is_valid)
        real(real64), intent(in) :: value, spvl

        is_valid = value == value
        if (spvl == spvl) then
            is_valid = is_valid .and. value /= spvl
        end if
    end function valid_value


    pure subroutine sort_pairs(x, y)
        real(real64), intent(inout) :: x(:), y(:)

        integer :: i, j, n
        real(real64) :: tx, ty

        n = size(x)
        do i = 2, n
            tx = x(i)
            ty = y(i)
            j = i - 1
            do while (j >= 1 .and. x(j) > tx)
                x(j + 1) = x(j)
                y(j + 1) = y(j)
                j = j - 1
            end do
            x(j + 1) = tx
            y(j + 1) = ty
        end do
    end subroutine sort_pairs


    pure real(real64) function interp_sorted_profile(x, y, nvalid, target, kxtrp, spvl) result(value)
        real(real64), intent(in) :: x(:), y(:), target, spvl
        integer, intent(in) :: nvalid, kxtrp

        integer :: k
        real(real64) :: x0, x1, y0, y1

        value = spvl
        if (nvalid < 2) return

        if (target < x(1)) then
            if (kxtrp /= 0) value = y(1)
            return
        end if
        if (target > x(nvalid)) then
            if (kxtrp /= 0) value = y(nvalid)
            return
        end if

        do k = 1, nvalid
            if (target == x(k)) then
                value = y(k)
                return
            end if
        end do

        do k = 1, nvalid - 1
            if (target >= x(k) .and. target <= x(k + 1)) then
                x0 = x(k)
                x1 = x(k + 1)
                y0 = y(k)
                y1 = y(k + 1)
                if (x1 == x0) then
                    value = y0
                else
                    value = y0 + (y1 - y0) * (target - x0) / (x1 - x0)
                end if
                return
            end if
        end do
    end function interp_sorted_profile


    pure subroutine interp_isentropic_flat_profile( &
        data_flat, temp_flat, pressure_flat, theta_levels, output_flat, &
        p0, kappa, spvl, kxtrp, base_in, base_out, inner, ninner, nlev, ntheta)
        real(real64), intent(in) :: data_flat(:), temp_flat(:), pressure_flat(:), theta_levels(:)
        real(real64), intent(inout) :: output_flat(:)
        real(real64), intent(in) :: p0, kappa, spvl
        integer, intent(in) :: kxtrp, base_in, base_out, inner, ninner, nlev, ntheta

        integer :: k, out, idx_in, idx_out, nvalid
        real(real64) :: theta(nlev), values(nlev)

        nvalid = 0
        idx_in = base_in + inner
        do k = 1, nlev
            if (valid_value(data_flat(idx_in), spvl) .and. &
                valid_value(temp_flat(idx_in), spvl) .and. &
                valid_value(pressure_flat(idx_in), spvl) .and. pressure_flat(idx_in) > 0.0_real64) then
                nvalid = nvalid + 1
                theta(nvalid) = temp_flat(idx_in) * (p0 / pressure_flat(idx_in)) ** kappa
                values(nvalid) = data_flat(idx_in)
            end if
            idx_in = idx_in + ninner
        end do

        idx_out = base_out + inner
        do out = 1, ntheta
            output_flat(idx_out) = spvl
            idx_out = idx_out + ninner
        end do
        if (nvalid < 2) return

        call sort_pairs(theta(:nvalid), values(:nvalid))
        idx_out = base_out + inner
        do out = 1, ntheta
            output_flat(idx_out) = interp_sorted_profile( &
                theta(:nvalid), values(:nvalid), nvalid, theta_levels(out), kxtrp, spvl)
            idx_out = idx_out + ninner
        end do
    end subroutine interp_isentropic_flat_profile

end module isentropic_kernels_core


subroutine dinterp_to_isentropic_corder_into_impl( &
    data_flat, temp_flat, pressure_flat, theta_levels, output_flat, &
    p0, kappa, spvl, kxtrp, nouter, nlev, ninner, ntheta)
    use isentropic_kernels_core, only : real64, interp_isentropic_flat_profile
    implicit none

    integer, intent(in) :: nouter, nlev, ninner, ntheta, kxtrp
    real(real64), intent(in) :: data_flat(nouter * nlev * ninner)
    real(real64), intent(in) :: temp_flat(nouter * nlev * ninner)
    real(real64), intent(in) :: pressure_flat(nouter * nlev * ninner)
    real(real64), intent(in) :: theta_levels(ntheta), p0, kappa, spvl
    real(real64), intent(inout) :: output_flat(nouter * ntheta * ninner)

    integer :: outer, inner, base_in, base_out

    do outer = 0, nouter - 1
        base_in = outer * nlev * ninner + 1
        base_out = outer * ntheta * ninner + 1
        do inner = 0, ninner - 1
            call interp_isentropic_flat_profile( &
                data_flat, temp_flat, pressure_flat, theta_levels, output_flat, &
                p0, kappa, spvl, kxtrp, base_in, base_out, inner, ninner, nlev, ntheta)
        end do
    end do
end subroutine dinterp_to_isentropic_corder_into_impl


subroutine dinterp_to_isentropic_corder_into( &
    data_flat, temp_flat, pressure_flat, theta_levels, output_flat, &
    p0, kappa, spvl, kxtrp, nouter, nlev, ninner, ntheta) &
    bind(C, name="dinterp_to_isentropic_corder_into")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use isentropic_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data_flat, temp_flat, pressure_flat, theta_levels, output_flat
    real(c_double), value, intent(in) :: p0, kappa, spvl
    integer(c_int), value, intent(in) :: kxtrp, nouter, nlev, ninner, ntheta

    real(real64), pointer :: data_flat_view(:), temp_flat_view(:), pressure_flat_view(:), theta_levels_view(:)
    real(real64), pointer :: output_flat_view(:)
    integer :: nouter_f, nlev_f, ninner_f, ntheta_f, kxtrp_f

    nouter_f = int(nouter)
    nlev_f = int(nlev)
    ninner_f = int(ninner)
    ntheta_f = int(ntheta)
    kxtrp_f = int(kxtrp)

    call c_f_pointer(data_flat, data_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(temp_flat, temp_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(pressure_flat, pressure_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(theta_levels, theta_levels_view, [ntheta_f])
    call c_f_pointer(output_flat, output_flat_view, [nouter_f * ntheta_f * ninner_f])

    call dinterp_to_isentropic_corder_into_impl( &
        data_flat_view, temp_flat_view, pressure_flat_view, theta_levels_view, output_flat_view, &
        p0, kappa, spvl, kxtrp_f, nouter_f, nlev_f, ninner_f, ntheta_f &
    )
end subroutine dinterp_to_isentropic_corder_into
