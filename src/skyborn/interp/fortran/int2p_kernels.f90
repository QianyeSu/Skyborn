! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-05-03
! File    : int2p_kernels.f90
! Purpose : Modern one-dimensional pressure-coordinate interpolation kernels.
! Notes   : This is a clean-room reimplementation of the archived
!           `archive/int2p_dp.f` behavior for production use. The archived
!           F77 source is retained only as provenance and validation reference.
! =============================================================================

module int2p_kernels_core
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    public :: real64
    public :: compact_valid_levels
    public :: copy_descending
    public :: extrapolate_log_low_compat
    public :: interpolate_value
    public :: is_increasing
    public :: locate_exact_match_desc
    public :: locate_interval_desc

contains

    pure logical function is_increasing(values) result(flag)
        real(real64), intent(in) :: values(:)

        if (size(values) < 2) then
            flag = .false.
        else
            flag = values(1) < values(2)
        end if
    end function is_increasing


    pure subroutine copy_descending(source, target, reverse_input)
        real(real64), intent(in) :: source(:)
        real(real64), intent(out) :: target(size(source))
        logical, intent(in) :: reverse_input

        integer :: idx, n

        n = size(source)
        if (reverse_input) then
            do idx = 1, n
                target(idx) = source(n - idx + 1)
            end do
        else
            target = source
        end if
    end subroutine copy_descending


    pure subroutine compact_valid_levels(p_in, x_in, xmsg, p_valid, x_valid, nvalid)
        real(real64), intent(in) :: p_in(:), x_in(:), xmsg
        real(real64), intent(out) :: p_valid(size(p_in)), x_valid(size(x_in))
        integer, intent(out) :: nvalid

        integer :: idx

        nvalid = 0
        do idx = 1, size(p_in)
            if (p_in(idx) /= xmsg .and. x_in(idx) /= xmsg) then
                nvalid = nvalid + 1
                p_valid(nvalid) = p_in(idx)
                x_valid(nvalid) = x_in(idx)
            end if
        end do
    end subroutine compact_valid_levels


    pure integer function locate_exact_match_desc(target, levels) result(idx)
        real(real64), intent(in) :: target
        real(real64), intent(in) :: levels(:)

        integer :: high, low, mid

        idx = 0
        low = 1
        high = size(levels)
        do while (low <= high)
            mid = low + (high - low) / 2
            if (target == levels(mid)) then
                idx = mid
                return
            end if

            if (target > levels(mid)) then
                high = mid - 1
            else
                low = mid + 1
            end if
        end do
    end function locate_exact_match_desc


    pure integer function locate_interval_desc(target, levels) result(idx)
        real(real64), intent(in) :: target
        real(real64), intent(in) :: levels(:)

        integer :: high, low, mid, nlev
        real(real64) :: p_hi, p_lo

        idx = 0
        nlev = size(levels)
        if (nlev < 2) return

        if (target >= levels(1) .or. target <= levels(nlev)) return

        low = 1
        high = nlev - 1
        do while (low <= high)
            mid = low + (high - low) / 2
            p_hi = levels(mid)
            p_lo = levels(mid + 1)

            if (target < p_hi .and. target > p_lo) then
                idx = mid
                return
            end if

            if (target > p_hi) then
                high = mid - 1
            else
                low = mid + 1
            end if
        end do
    end function locate_interval_desc


    pure real(real64) function interpolate_linear( &
        x_hi, x_lo, p_target, p_hi, p_lo &
    ) result(value)
        real(real64), intent(in) :: x_hi, x_lo, p_target, p_hi, p_lo

        value = x_lo + ((x_hi - x_lo) / (p_hi - p_lo)) * (p_target - p_lo)
    end function interpolate_linear


    pure real(real64) function interpolate_log( &
        x_hi, x_lo, p_target, p_hi, p_lo &
    ) result(value)
        real(real64), intent(in) :: x_hi, x_lo, p_target, p_hi, p_lo

        value = x_lo + ((x_hi - x_lo) / (log(p_hi) - log(p_lo))) * &
            (log(p_target) - log(p_lo))
    end function interpolate_log


    pure real(real64) function extrapolate_log_low_compat( &
        x_hi, x_lo, p_target, p_hi, p_lo &
    ) result(value)
        real(real64), intent(in) :: x_hi, x_lo, p_target, p_hi, p_lo

        ! Preserve the archived DINT2P low-end log extrapolation convention
        ! exactly so the modern kernel stays numerically compatible with the
        ! historical reference implementation.
        value = x_lo + ((x_lo - x_hi) / (log(p_lo) - log(p_hi))) * &
            (log(p_target) - log(p_hi))
    end function extrapolate_log_low_compat


    pure real(real64) function interpolate_value( &
        x_hi, x_lo, p_target, p_hi, p_lo, loglin &
    ) result(value)
        real(real64), intent(in) :: x_hi, x_lo, p_target, p_hi, p_lo
        integer, intent(in) :: loglin

        if (loglin == 1) then
            value = interpolate_linear(x_hi, x_lo, p_target, p_hi, p_lo)
        else
            value = interpolate_log(x_hi, x_lo, p_target, p_hi, p_lo)
        end if
    end function interpolate_value

end module int2p_kernels_core


subroutine dinterp_pressure_1d_impl(ppin, xxin, ppout, xxout, linlog, xmsg, ier, npin, npout)
    use int2p_kernels_core, only : &
        real64, compact_valid_levels, copy_descending, extrapolate_log_low_compat, &
        interpolate_value, is_increasing, locate_exact_match_desc, &
        locate_interval_desc
    implicit none

    integer, intent(in) :: npin, npout, linlog
    real(real64), intent(in) :: ppin(npin), xxin(npin), ppout(npout), xmsg
    real(real64), intent(out) :: xxout(npout)
    integer, intent(out) :: ier

    real(real64) :: pin_desc(npin), xin_desc(npin), pout_desc(npout)
    real(real64) :: p_valid(npin), x_valid(npin), xout_desc(npout)
    integer :: exact_idx, idx, interval_idx, loglin, nvalid
    logical :: reverse_input, reverse_output

    xxout = xmsg
    ier = 0
    if (npin < 2 .or. npout < 1) then
        ier = 1
        return
    end if

    reverse_input = is_increasing(ppin)
    reverse_output = .false.
    if (npout > 1) reverse_output = is_increasing(ppout)

    call copy_descending(ppin, pin_desc, reverse_input)
    call copy_descending(xxin, xin_desc, reverse_input)
    call copy_descending(ppout, pout_desc, reverse_output)
    call compact_valid_levels(pin_desc, xin_desc, xmsg, p_valid, x_valid, nvalid)

    if (nvalid < 2) then
        ier = 1000
        return
    end if

    loglin = abs(linlog)
    xout_desc = xmsg
    do idx = 1, npout
        if (pout_desc(idx) == xmsg) cycle

        exact_idx = locate_exact_match_desc(pout_desc(idx), p_valid(:nvalid))
        if (exact_idx > 0) then
            xout_desc(idx) = x_valid(exact_idx)
            cycle
        end if

        if (pout_desc(idx) > p_valid(1)) then
            if (linlog < 0) then
                xout_desc(idx) = interpolate_value( &
                    x_valid(1), x_valid(2), pout_desc(idx), &
                    p_valid(1), p_valid(2), loglin &
                )
            end if
            cycle
        end if

        if (pout_desc(idx) < p_valid(nvalid)) then
            if (linlog < 0) then
                if (loglin == 1) then
                    xout_desc(idx) = interpolate_value( &
                        x_valid(nvalid - 1), x_valid(nvalid), pout_desc(idx), &
                        p_valid(nvalid - 1), p_valid(nvalid), loglin &
                    )
                else
                    xout_desc(idx) = extrapolate_log_low_compat( &
                        x_valid(nvalid - 1), x_valid(nvalid), pout_desc(idx), &
                        p_valid(nvalid - 1), p_valid(nvalid) &
                    )
                end if
            end if
            cycle
        end if

        interval_idx = locate_interval_desc(pout_desc(idx), p_valid(:nvalid))
        if (interval_idx > 0) then
            xout_desc(idx) = interpolate_value( &
                x_valid(interval_idx), x_valid(interval_idx + 1), pout_desc(idx), &
                p_valid(interval_idx), p_valid(interval_idx + 1), loglin &
            )
        end if
    end do

    if (reverse_output) then
        do idx = 1, npout
            xxout(idx) = xout_desc(npout - idx + 1)
        end do
    else
        xxout = xout_desc
    end if
end subroutine dinterp_pressure_1d_impl
