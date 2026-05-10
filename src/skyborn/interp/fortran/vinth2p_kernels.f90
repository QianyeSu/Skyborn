! =============================================================================
! Author  : Qianye Su
! Copyright (C) 2026 by Qianye Su
! Created : 2026-04-17
! File    : vinth2p_kernels.f90
! Purpose : Collect the modern free-form vertical interpolation kernels used
!           for hybrid-to-pressure and sigma-to-hybrid remapping.
! Notes   : This file keeps the legacy interpolation formulas, below-ground
!           extrapolation behavior, and F2PY-facing entry points while sharing
!           arithmetic between column-major and flat C-order call paths.
! =============================================================================
!
module vinth2p_kernels_core
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    ! Keep the scalar interpolation and extrapolation formulas in one module so
    ! the column-major and flat C-order entry points share identical arithmetic.
    real(real64), parameter :: rd = 287.04_real64
    real(real64), parameter :: ginv = 1.0_real64 / 9.80616_real64
    real(real64), parameter :: alpha = 0.0065_real64 * rd * ginv

    public :: real64
    public :: build_input_pressures
    public :: build_double_log_levels
    public :: compute_delta_pressure_columns
    public :: compute_geopotential_height_columns
    public :: compute_geopotential_height_flat_column
    public :: compute_pressure_at_hybrid_levels_columns
    public :: convert_levels_to_mb
    public :: extrapolate_geopotential
    public :: extrapolate_temperature
    public :: interpolate_column_ecmwf
    public :: interpolate_column_sigma
    public :: interpolate_flat_column_sigma
    public :: interpolate_flat_column_ecmwf
    public :: interpolate_flat_column_nodes
    public :: interpolate_column_nodes
    public :: interpolate_value
    public :: is_below_lowest_model_level
    public :: levels_ascending
    public :: legacy_ecmwf_bracketing_floor
    public :: locate_bracketing_level_ordered
    public :: lowest_bracketing_level_index
    public :: lowest_model_level_index
    public :: locate_bracketing_level

contains

    ! Convert public pressure levels from Pa to the legacy mb arithmetic used
    ! by the original vinth2p kernels.
    pure subroutine convert_levels_to_mb(plevo, plevo_mb)
        real(real64), intent(in) :: plevo(:)
        real(real64), intent(out) :: plevo_mb(size(plevo))

        plevo_mb = plevo * 0.01_real64
    end subroutine convert_levels_to_mb


    ! Build the hybrid pressure profile for one column from the A/B
    ! coefficients, reference pressure, and surface pressure.
    pure subroutine build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        real(real64), intent(in) :: hbcofa(:), hbcofb(:)
        real(real64), intent(in) :: p0_mb, psfc_mb
        real(real64), intent(out) :: plevi(size(hbcofa))

        plevi = hbcofa * p0_mb + hbcofb * psfc_mb
    end subroutine build_input_pressures


    ! Compute hybrid-layer pressure thickness directly without materializing
    ! the full adjacent pressure fields in temporary arrays.
    pure subroutine compute_delta_pressure_columns(psfc, hbcofa, hbcofb, p0, dph)
        real(real64), intent(in) :: psfc(:), hbcofa(:), hbcofb(:), p0
        real(real64), intent(out) :: dph(size(hbcofa) - 1, size(psfc))

        integer :: col, k
        real(real64) :: delta_a(size(hbcofa) - 1), delta_b(size(hbcofa) - 1)

        do k = 1, size(hbcofa) - 1
            delta_a(k) = hbcofa(k) - hbcofa(k + 1)
            delta_b(k) = hbcofb(k) - hbcofb(k + 1)
        end do

        do col = 1, size(psfc)
            do k = 1, size(hbcofa) - 1
                dph(k, col) = abs(p0 * delta_a(k) + delta_b(k) * psfc(col))
            end do
        end do
    end subroutine compute_delta_pressure_columns


    ! Compute hybrid pressures directly for every output level and column.
    pure subroutine compute_pressure_at_hybrid_levels_columns(psfc, hbcofa, hbcofb, p0, pressure)
        real(real64), intent(in) :: psfc(:), hbcofa(:), hbcofb(:), p0
        real(real64), intent(out) :: pressure(size(hbcofa), size(psfc))

        integer :: col
        real(real64) :: scaled_a(size(hbcofa))

        scaled_a = p0 * hbcofa
        do col = 1, size(psfc)
            pressure(:, col) = scaled_a + hbcofb * psfc(col)
        end do
    end subroutine compute_pressure_at_hybrid_levels_columns


    ! Compute geopotential height on full model levels from virtual temperature,
    ! hybrid coefficients, surface pressure, and surface geopotential.
    !
    ! This implementation follows the ECMWF ERA5 model-level geopotential
    ! recipe for:
    !   - virtual temperature: Tv = T * (1 + 0.609133 * q)
    !   - bottom-up integration from surface geopotential
    !   - the special top-layer shortcut only when the top interface pressure
    !     is exactly zero, as in the ERA5 half-level definition
    !
    ! Reference:
    !   ERA5: compute pressure and geopotential on model levels,
    !   geopotential height and geometric height
    !   https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    !
    ! Use the general hypsometric full-layer integral whenever the top
    ! interface pressure is finite. Keep the legacy ECMWF/IFS top-layer
    ! shortcut only for coordinate systems whose top interface is exactly
    ! zero-pressure.
    pure subroutine compute_geopotential_height_columns( &
        temp, q, psfc, phis, hyai, hybi, hyam, hybm, p0, z3)
        real(real64), intent(in) :: temp(:, :), q(:, :)
        real(real64), intent(in) :: psfc(:), phis(:), hyai(:), hybi(:), hyam(:), hybm(:), p0
        real(real64), intent(out) :: z3(size(hyam), size(psfc))

        integer :: col, k, nlev, nlevi
        real(real64), parameter :: rd = 287.06_real64
        real(real64), parameter :: g = 9.80665_real64
        real(real64) :: p_half(size(hyai))
        real(real64) :: p_upper, p_lower, dlog_p, alpha_full
        real(real64) :: phi_half, tv_level, phi_full

        nlev = size(hyam)
        nlevi = size(hyai)
        do col = 1, size(psfc)
            do k = 1, nlevi
                p_half(k) = hyai(k) * p0 + hybi(k) * psfc(col)
            end do

            phi_half = phis(col)
            do k = nlev, 1, -1
                p_upper = p_half(k)
                p_lower = p_half(k + 1)

                if (k == 1 .and. p_upper <= 0.0_real64) then
                    dlog_p = log(p_lower / 0.1_real64)
                    alpha_full = log(2.0_real64)
                else
                    dlog_p = log(p_lower / p_upper)
                    alpha_full = 1.0_real64 - (p_upper / (p_lower - p_upper)) * dlog_p
                end if

                tv_level = rd * temp(k, col) * (1.0_real64 + 0.609133_real64 * q(k, col))
                phi_full = phi_half + tv_level * alpha_full
                z3(k, col) = phi_full / g
                phi_half = phi_half + tv_level * dlog_p
            end do
        end do
    end subroutine compute_geopotential_height_columns


    ! Compute one flat C-order column of geopotential height in place.
    pure subroutine compute_geopotential_height_flat_column( &
        temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi, p0, &
        base_in, base_out, inner, ninner, nlev)
        real(real64), intent(in) :: temp_flat(:), q_flat(:)
        real(real64), intent(inout) :: z3_flat(:)
        real(real64), intent(in) :: psfc, phis, hyai(:), hybi(:), p0
        integer, intent(in) :: base_in, base_out, inner, ninner, nlev

        integer :: idx_in, idx_out, k
        real(real64), parameter :: rd = 287.06_real64
        real(real64), parameter :: g = 9.80665_real64
        real(real64) :: p_half(size(hyai))
        real(real64) :: p_upper, p_lower, dlog_p, alpha_full
        real(real64) :: phi_half, tv_level, phi_full

        do k = 1, size(hyai)
            p_half(k) = hyai(k) * p0 + hybi(k) * psfc
        end do

        phi_half = phis
        idx_in = base_in + inner + (nlev - 1) * ninner
        idx_out = base_out + inner + (nlev - 1) * ninner
        do k = nlev, 1, -1
            p_upper = p_half(k)
            p_lower = p_half(k + 1)

            if (k == 1 .and. p_upper <= 0.0_real64) then
                dlog_p = log(p_lower / 0.1_real64)
                alpha_full = log(2.0_real64)
            else
                dlog_p = log(p_lower / p_upper)
                alpha_full = 1.0_real64 - (p_upper / (p_lower - p_upper)) * dlog_p
            end if

            tv_level = rd * temp_flat(idx_in) * (1.0_real64 + 0.609133_real64 * q_flat(idx_in))
            phi_full = phi_half + tv_level * alpha_full
            z3_flat(idx_out) = phi_full / g
            phi_half = phi_half + tv_level * dlog_p
            idx_in = idx_in - ninner
            idx_out = idx_out - ninner
        end do
    end subroutine compute_geopotential_height_flat_column


    ! Build the flat-buffer input offsets for one physical column so the hot
    ! loops can index levels by lookup instead of repeating stride multiplies.
    pure subroutine build_flat_level_offsets(base_in, inner, ninner, offsets)
        integer, intent(in) :: base_in, inner, ninner
        integer, intent(out) :: offsets(:)

        integer :: idx, k

        idx = base_in + inner
        do k = 1, size(offsets)
            offsets(k) = idx
            idx = idx + ninner
        end do
    end subroutine build_flat_level_offsets


    ! Precompute the legacy double-log transform for one coordinate vector.
    pure subroutine build_double_log_levels(levels_mb, levels_double_log)
        real(real64), intent(in) :: levels_mb(:)
        real(real64), intent(out) :: levels_double_log(size(levels_mb))

        integer :: k

        do k = 1, size(levels_mb)
            levels_double_log(k) = double_log_pressure(levels_mb(k))
        end do
    end subroutine build_double_log_levels


    ! Build the exact legacy double-log denominator for each adjacent source
    ! layer pair from precomputed transformed coordinates.
    pure subroutine build_double_log_interval_denominators( &
        levels_double_log, double_log_denominators &
    )
        real(real64), intent(in) :: levels_double_log(:)
        real(real64), intent(out) :: double_log_denominators(size(levels_double_log) - 1)

        double_log_denominators = &
            levels_double_log(2:) - levels_double_log(:size(levels_double_log) - 1)
    end subroutine build_double_log_interval_denominators


    ! Return whether the source coordinate increases from top to bottom.
    pure logical function levels_ascending(plevi) result(is_ascending)
        real(real64), intent(in) :: plevi(:)

        ! Some tests and callers hand us the vertical axis in descending order.
        is_ascending = plevi(1) <= plevi(size(plevi))
    end function levels_ascending


    ! Return the index of the physically lowest model level for either
    ! ascending or descending source-coordinate order.
    pure integer function lowest_model_level_index(plevi) result(idx)
        real(real64), intent(in) :: plevi(:)

        if (levels_ascending(plevi)) then
            idx = size(plevi)
        else
            idx = 1
        end if
    end function lowest_model_level_index


    ! Return the upper index of the lowest usable bracketing pair.
    pure integer function lowest_bracketing_level_index(plevi) result(idx)
        real(real64), intent(in) :: plevi(:)

        if (levels_ascending(plevi)) then
            idx = size(plevi) - 1
        else
            idx = 1
        end if
    end function lowest_bracketing_level_index


    ! Check whether a requested pressure lies below the lowest model level.
    pure logical function is_below_lowest_model_level(plev_out, plevi) result(is_below)
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)

        is_below = plev_out > plevi(lowest_model_level_index(plevi))
    end function is_below_lowest_model_level


    ! Convenience wrapper that detects source-coordinate order before calling
    ! the ordered binary-search helper.
    pure integer function locate_bracketing_level(plev_out, plevi) result(kp)
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)

        kp = locate_bracketing_level_ordered(plev_out, plevi, levels_ascending(plevi))
    end function locate_bracketing_level


    ! Locate the upper index of the bracketing source interval for one target
    ! coordinate value, assuming the caller already knows the source order.
    pure integer function locate_bracketing_level_ordered( &
        plev_out, plevi, is_ascending &
    ) result(kp)
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)
        logical, intent(in) :: is_ascending

        integer :: high, low, mid, nlev

        nlev = size(plevi)
        if (nlev < 2) then
            kp = 0
            return
        end if

        ! Guard the top and bottom boundary pairs explicitly before the binary
        ! search so exact legacy edge behavior is preserved.
        if (is_ascending) then
            if (plev_out <= plevi(1)) then
                kp = 1
                return
            end if

            if (plev_out >= plevi(nlev - 1)) then
                kp = nlev - 1
                return
            end if
        else
            if (plev_out >= plevi(1)) then
                kp = 1
                return
            end if

            if (plev_out <= plevi(nlev - 1)) then
                kp = nlev - 1
                return
            end if
        end if

        low = 1
        high = nlev - 1
        do while (low < high)
            mid = low + (high - low) / 2
            if (is_ascending) then
                if (plev_out <= plevi(mid + 1)) then
                    high = mid
                else
                    low = mid + 1
                end if
            else
                if (plev_out >= plevi(mid + 1)) then
                    high = mid
                else
                    low = mid + 1
                end if
            end if
        end do

        kp = low
    end function locate_bracketing_level_ordered


    ! Match the archived ECMWF node kernel near the model top. The legacy F77
    ! code starts the interior search with KP = 1 and increments before the
    ! first comparison, so the first interior bracket it can return is the
    ! second model-layer pair. Keep this helper only for archived F77 parity:
    ! the public GeoCAT-compatible ECMWF path must not use it for plev <= p_sfc,
    ! because GeoCAT keeps the standard ordered bracket for in-column remapping
    ! and only applies the ECMWF formulation to below-ground extrapolation.
    pure integer function legacy_ecmwf_bracketing_floor( &
        kp, plev_out, plevi, is_ascending &
    ) result(adjusted_kp)
        integer, intent(in) :: kp
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)
        logical, intent(in) :: is_ascending

        adjusted_kp = kp
        if (.not. is_ascending) return
        if (size(plevi) < 3) return
        if (plev_out <= plevi(1)) return

        if (adjusted_kp < 2) adjusted_kp = 2
    end function legacy_ecmwf_bracketing_floor


    ! Legacy double-log transform used by vinth2p "log-log" interpolation.
    pure real(real64) function double_log_pressure(pressure_mb) result(value)
        real(real64), intent(in) :: pressure_mb

        value = log(log(pressure_mb + 2.72_real64))
    end function double_log_pressure


    ! Legacy linear interpolation on pressure levels.
    pure real(real64) function interpolate_value_linear( &
        lower_value, upper_value, plev_out, plev_lower, plev_upper &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out, plev_lower, plev_upper

        value = lower_value + &
            (upper_value - lower_value) * &
            (plev_out - plev_lower) / &
            (plev_upper - plev_lower)
    end function interpolate_value_linear


    ! Legacy log-pressure interpolation on pressure levels.
    pure real(real64) function interpolate_value_log( &
        lower_value, upper_value, plev_out, plev_lower, plev_upper &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out, plev_lower, plev_upper

        value = lower_value + &
            (upper_value - lower_value) * &
            log(plev_out / plev_lower) / &
            log(plev_upper / plev_lower)
    end function interpolate_value_log


    ! Legacy double-log interpolation on pressure levels.
    pure real(real64) function interpolate_value_loglog( &
        lower_value, upper_value, plev_out, plev_lower, plev_upper &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out, plev_lower, plev_upper

        value = lower_value + &
            (upper_value - lower_value) * &
            (double_log_pressure(plev_out) - double_log_pressure(plev_lower)) / &
            (double_log_pressure(plev_upper) - double_log_pressure(plev_lower))
    end function interpolate_value_loglog


    ! Legacy double-log interpolation that reuses precomputed transformed
    ! coordinates for the target and bracketing source pair.
    pure real(real64) function interpolate_value_loglog_precomputed( &
        lower_value, upper_value, plev_out_double_log, plev_lower_double_log, &
        double_log_denominator &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out_double_log, plev_lower_double_log
        real(real64), intent(in) :: double_log_denominator

        value = lower_value + &
            (upper_value - lower_value) * &
            (plev_out_double_log - plev_lower_double_log) / &
            double_log_denominator
    end function interpolate_value_loglog_precomputed


    ! Apply the legacy vinth2p interpolation formula for one bracketing pair.
    pure real(real64) function interpolate_value( &
        lower_value, upper_value, plev_out, plev_lower, plev_upper, intyp, spvl &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out, plev_lower, plev_upper, spvl
        integer, intent(in) :: intyp

        select case (intyp)
        case (1)
            value = interpolate_value_linear( &
                lower_value, upper_value, plev_out, plev_lower, plev_upper &
            )
        case (2)
            value = interpolate_value_log( &
                lower_value, upper_value, plev_out, plev_lower, plev_upper &
            )
        case (3)
            value = interpolate_value_loglog( &
                lower_value, upper_value, plev_out, plev_lower, plev_upper &
            )
        case default
            value = spvl
        end select
    end function interpolate_value


    ! Interpolate one vertical column from hybrid levels to requested pressure
    ! levels without ECMWF below-ground extrapolation.
    subroutine interpolate_column_nodes( &
        dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
        intyp, psfc, spvl, kxtrp &
    )
        real(real64), intent(in) :: dati_col(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: plevo_loglog(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl
        real(real64), intent(out) :: dato_col(size(plevo_mb))
        integer, intent(in) :: intyp, kxtrp

        real(real64) :: bottom_pressure, data_lower, data_upper
        real(real64) :: plevi(size(dati_col)), plev_lower, plev_upper
        real(real64) :: plev_lower_loglog, plev_upper_loglog, psfc_mb
        integer :: bottom_idx, bottom_pair_idx, k, kp, kp_hint
        logical :: is_ascending, use_descending_walk

        ! The Python bridge writes the exact sentinel value, so direct equality
        ! is intentional here.
        if (psfc == spvl) then
            dato_col = spvl
            return
        end if

        ! Build hybrid pressures once per column, then reuse them for every
        ! requested output pressure level.
        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pair_idx = lowest_bracketing_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)
        use_descending_walk = is_ascending
        if (use_descending_walk) then
            do k = 2, size(plevo_mb)
                if (plevo_mb(k) > plevo_mb(k - 1)) then
                    use_descending_walk = .false.
                    exit
                end if
            end do
        end if
        kp_hint = bottom_pair_idx

        select case (intyp)
        case (1)
            do k = 1, size(plevo_mb)
                ! Below-ground requests either return the sentinel or use the same
                ! lowest bracketing pair as the legacy node kernel.
                if (plevo_mb(k) > bottom_pressure) then
                    if (kxtrp == 0) then
                        dato_col(k) = spvl
                    else
                        kp = bottom_pair_idx
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower = plevi(kp)
                        plev_upper = plevi(kp + 1)
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            (plevo_mb(k) - plev_lower) / &
                            (plev_upper - plev_lower)
                    end if
                else
                    if (use_descending_walk) then
                        do while (kp_hint > 1 .and. plevo_mb(k) < plevi(kp_hint))
                            kp_hint = kp_hint - 1
                        end do
                        kp = kp_hint
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                    end if
                    if (kp < 1) then
                        dato_col(k) = spvl
                    else
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower = plevi(kp)
                        plev_upper = plevi(kp + 1)
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            (plevo_mb(k) - plev_lower) / &
                            (plev_upper - plev_lower)
                    end if
                end if
            end do
        case (2)
            do k = 1, size(plevo_mb)
                if (plevo_mb(k) > bottom_pressure) then
                    if (kxtrp == 0) then
                        dato_col(k) = spvl
                    else
                        kp = bottom_pair_idx
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower = plevi(kp)
                        plev_upper = plevi(kp + 1)
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            log(plevo_mb(k) / plev_lower) / &
                            log(plev_upper / plev_lower)
                    end if
                else
                    if (use_descending_walk) then
                        do while (kp_hint > 1 .and. plevo_mb(k) < plevi(kp_hint))
                            kp_hint = kp_hint - 1
                        end do
                        kp = kp_hint
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                    end if
                    if (kp < 1) then
                        dato_col(k) = spvl
                    else
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower = plevi(kp)
                        plev_upper = plevi(kp + 1)
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            log(plevo_mb(k) / plev_lower) / &
                            log(plev_upper / plev_lower)
                    end if
                end if
            end do
        case (3)
            do k = 1, size(plevo_mb)
                if (plevo_mb(k) > bottom_pressure) then
                    if (kxtrp == 0) then
                        dato_col(k) = spvl
                    else
                        kp = bottom_pair_idx
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower_loglog = double_log_pressure(plevi(kp))
                        plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            (plevo_loglog(k) - plev_lower_loglog) / &
                            (plev_upper_loglog - plev_lower_loglog)
                    end if
                else
                    if (use_descending_walk) then
                        do while (kp_hint > 1 .and. plevo_mb(k) < plevi(kp_hint))
                            kp_hint = kp_hint - 1
                        end do
                        kp = kp_hint
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                    end if
                    if (kp < 1) then
                        dato_col(k) = spvl
                    else
                        data_lower = dati_col(kp)
                        data_upper = dati_col(kp + 1)
                        plev_lower_loglog = double_log_pressure(plevi(kp))
                        plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                        dato_col(k) = data_lower + &
                            (data_upper - data_lower) * &
                            (plevo_loglog(k) - plev_lower_loglog) / &
                            (plev_upper_loglog - plev_lower_loglog)
                    end if
                end if
            end do
        case default
            dato_col = spvl
        end select
    end subroutine interpolate_column_nodes


    ! Interpolate one vertical column from hybrid levels to requested pressure
    ! levels with ECMWF below-ground extrapolation when requested.
    subroutine interpolate_column_ecmwf( &
        dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
        intyp, psfc, spvl, kxtrp, varflg, tbot, phis &
    )
        real(real64), intent(in) :: dati_col(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: plevo_loglog(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl, tbot, phis
        real(real64), intent(out) :: dato_col(size(plevo_mb))
        integer, intent(in) :: intyp, kxtrp, varflg

        real(real64) :: bottom_pressure, data_lower, data_upper
        real(real64) :: plevi(size(dati_col)), plev_lower, plev_upper
        real(real64) :: plev_lower_loglog, plev_upper_loglog, psfc_mb
        integer :: bottom_idx, k, kp
        logical :: is_ascending

        ! The Python bridge writes the exact sentinel value, so direct equality
        ! is intentional here.
        if (psfc == spvl) then
            dato_col = spvl
            return
        end if

        ! GeoCAT-compatible ECMWF behavior keeps the standard ordered bracket
        ! for plev <= p_sfc and reserves the ECMWF formulation for the
        ! below-ground branch only. Do not apply legacy_ecmwf_bracketing_floor()
        ! here: that F77-only quirk shifts 5/7 hPa requests off the topmost
        ! bracket pair and produces large high-top mismatches versus GeoCAT.
        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)

        select case (intyp)
        case (1)
            if (kxtrp == 0) then
                do k = 1, size(plevo_mb)
                    if (plevo_mb(k) > bottom_pressure) then
                        dato_col(k) = spvl
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_col(k) = spvl
                        else
                            data_lower = dati_col(kp)
                            data_upper = dati_col(kp + 1)
                            plev_lower = plevi(kp)
                            plev_upper = plevi(kp + 1)
                            dato_col(k) = data_lower + &
                                (data_upper - data_lower) * &
                                (plevo_mb(k) - plev_lower) / &
                                (plev_upper - plev_lower)
                        end if
                    end if
                end do
            else
                select case (varflg)
                case (1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_temperature( &
                                dati_col(bottom_idx), plevi(bottom_idx), &
                                plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_mb(k) - plev_lower) / &
                                    (plev_upper - plev_lower)
                            end if
                        end if
                    end do
                case (-1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_geopotential( &
                                tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_mb(k) - plev_lower) / &
                                    (plev_upper - plev_lower)
                            end if
                        end if
                    end do
                case default
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = dati_col(bottom_idx)
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_mb(k) - plev_lower) / &
                                    (plev_upper - plev_lower)
                            end if
                        end if
                    end do
                end select
            end if
        case (2)
            if (kxtrp == 0) then
                do k = 1, size(plevo_mb)
                    if (plevo_mb(k) > bottom_pressure) then
                        dato_col(k) = spvl
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_col(k) = spvl
                        else
                            data_lower = dati_col(kp)
                            data_upper = dati_col(kp + 1)
                            plev_lower = plevi(kp)
                            plev_upper = plevi(kp + 1)
                            dato_col(k) = data_lower + &
                                (data_upper - data_lower) * &
                                log(plevo_mb(k) / plev_lower) / &
                                log(plev_upper / plev_lower)
                        end if
                    end if
                end do
            else
                select case (varflg)
                case (1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_temperature( &
                                dati_col(bottom_idx), plevi(bottom_idx), &
                                plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    log(plevo_mb(k) / plev_lower) / &
                                    log(plev_upper / plev_lower)
                            end if
                        end if
                    end do
                case (-1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_geopotential( &
                                tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    log(plevo_mb(k) / plev_lower) / &
                                    log(plev_upper / plev_lower)
                            end if
                        end if
                    end do
                case default
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = dati_col(bottom_idx)
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower = plevi(kp)
                                plev_upper = plevi(kp + 1)
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    log(plevo_mb(k) / plev_lower) / &
                                    log(plev_upper / plev_lower)
                            end if
                        end if
                    end do
                end select
            end if
        case (3)
            if (kxtrp == 0) then
                do k = 1, size(plevo_mb)
                    if (plevo_mb(k) > bottom_pressure) then
                        dato_col(k) = spvl
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_col(k) = spvl
                        else
                            data_lower = dati_col(kp)
                            data_upper = dati_col(kp + 1)
                            plev_lower_loglog = double_log_pressure(plevi(kp))
                            plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                            dato_col(k) = data_lower + &
                                (data_upper - data_lower) * &
                                (plevo_loglog(k) - plev_lower_loglog) / &
                                (plev_upper_loglog - plev_lower_loglog)
                        end if
                    end if
                end do
            else
                select case (varflg)
                case (1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_temperature( &
                                dati_col(bottom_idx), plevi(bottom_idx), &
                                plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower_loglog = double_log_pressure(plevi(kp))
                                plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_loglog(k) - plev_lower_loglog) / &
                                    (plev_upper_loglog - plev_lower_loglog)
                            end if
                        end if
                    end do
                case (-1)
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = extrapolate_geopotential( &
                                tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower_loglog = double_log_pressure(plevi(kp))
                                plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_loglog(k) - plev_lower_loglog) / &
                                    (plev_upper_loglog - plev_lower_loglog)
                            end if
                        end if
                    end do
                case default
                    do k = 1, size(plevo_mb)
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_col(k) = dati_col(bottom_idx)
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_col(k) = spvl
                            else
                                data_lower = dati_col(kp)
                                data_upper = dati_col(kp + 1)
                                plev_lower_loglog = double_log_pressure(plevi(kp))
                                plev_upper_loglog = double_log_pressure(plevi(kp + 1))
                                dato_col(k) = data_lower + &
                                    (data_upper - data_lower) * &
                                    (plevo_loglog(k) - plev_lower_loglog) / &
                                    (plev_upper_loglog - plev_lower_loglog)
                            end if
                        end if
                    end do
                end select
            end if
        case default
            dato_col = spvl
        end select
    end subroutine interpolate_column_ecmwf


    ! Interpolate one vertical column from sigma levels to target sigma values.
    ! This keeps the same linear/log formulas as the Python/MetPy path but
    ! executes the per-column loop in Fortran.
    subroutine interpolate_column_sigma( &
        dati_col, dato_col, sigmai, sigmao_col, intyp, spvl &
    )
        real(real64), intent(in) :: dati_col(:), sigmai(:), sigmao_col(:)
        real(real64), intent(out) :: dato_col(size(sigmao_col))
        real(real64), intent(in) :: spvl
        integer, intent(in) :: intyp

        logical :: is_ascending
        integer :: k, kp, nlevi

        nlevi = size(sigmai)
        if (nlevi < 2) then
            dato_col = spvl
            return
        end if

        is_ascending = levels_ascending(sigmai)
        do k = 1, size(sigmao_col)
            if (is_ascending) then
                if (sigmao_col(k) < sigmai(1) .or. sigmao_col(k) > sigmai(nlevi)) then
                    dato_col(k) = spvl
                    cycle
                end if
            else
                if (sigmao_col(k) > sigmai(1) .or. sigmao_col(k) < sigmai(nlevi)) then
                    dato_col(k) = spvl
                    cycle
                end if
            end if

            kp = locate_bracketing_level_ordered(sigmao_col(k), sigmai, is_ascending)
            if (kp < 1) then
                dato_col(k) = spvl
            else if (intyp == 2 .and. &
                (sigmao_col(k) <= 0.0_real64 .or. sigmai(kp) <= 0.0_real64 .or. sigmai(kp + 1) <= 0.0_real64)) then
                dato_col(k) = spvl
            else
                dato_col(k) = interpolate_value( &
                    dati_col(kp), dati_col(kp + 1), &
                    sigmao_col(k), sigmai(kp), sigmai(kp + 1), intyp, spvl &
                )
            end if
        end do
    end subroutine interpolate_column_sigma


    ! Flat-buffer helper for the C-order sigma->hybrid path. The target sigma
    ! values are computed on the fly from ps, hyam, and hybm so Python does not
    ! need to materialize a full sigma-target volume.
    subroutine interpolate_flat_column_sigma( &
        dati_flat, dato_flat, sigmai, hyam, hybm, p0, psfc, intyp, spvl, &
        base_in, base_out, inner, ninner, nlevi, nlevo &
    )
        real(real64), intent(in) :: dati_flat(:), sigmai(:), hyam(:), hybm(:)
        real(real64), intent(inout) :: dato_flat(:)
        real(real64), intent(in) :: p0, psfc, spvl
        integer, intent(in) :: intyp, base_in, base_out, inner
        integer, intent(in) :: ninner, nlevi, nlevo

        logical :: is_ascending
        integer :: input_idx, k, kp, output_idx
        real(real64) :: sigma_out

        is_ascending = levels_ascending(sigmai)
        do k = 1, nlevo
            output_idx = base_out + (k - 1) * ninner + inner
            sigma_out = hyam(k) * p0 / psfc + hybm(k)

            if (is_ascending) then
                if (sigma_out < sigmai(1) .or. sigma_out > sigmai(nlevi)) then
                    dato_flat(output_idx) = spvl
                    cycle
                end if
            else
                if (sigma_out > sigmai(1) .or. sigma_out < sigmai(nlevi)) then
                    dato_flat(output_idx) = spvl
                    cycle
                end if
            end if

            kp = locate_bracketing_level_ordered(sigma_out, sigmai, is_ascending)
            if (kp < 1) then
                dato_flat(output_idx) = spvl
            else if (intyp == 2 .and. &
                (sigma_out <= 0.0_real64 .or. sigmai(kp) <= 0.0_real64 .or. sigmai(kp + 1) <= 0.0_real64)) then
                dato_flat(output_idx) = spvl
            else
                input_idx = base_in + (kp - 1) * ninner + inner
                dato_flat(output_idx) = interpolate_value( &
                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                    sigma_out, sigmai(kp), sigmai(kp + 1), intyp, spvl &
                )
            end if
        end do
    end subroutine interpolate_flat_column_sigma


    ! Flat-buffer helper for the C-order hybrid->pressure path.
    subroutine interpolate_flat_column_nodes( &
        dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
        intyp, psfc, spvl, kxtrp, base_in, base_out, inner, ninner, nlevi, nlevo &
    )
        real(real64), intent(in) :: dati_flat(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: plevo_loglog(:)
        real(real64), intent(inout) :: dato_flat(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl
        integer, intent(in) :: intyp, kxtrp, base_in, base_out, inner
        integer, intent(in) :: ninner, nlevi, nlevo

        real(real64) :: bottom_pressure, plevi(nlevi), psfc_mb
        integer :: bottom_idx, bottom_pair_idx, bottom_pair_input_idx
        integer :: input_idx, input_offsets(nlevi), k, kp, output_idx
        logical :: is_ascending

        if (psfc == spvl) then
            output_idx = base_out + inner
            do k = 1, nlevo
                dato_flat(output_idx) = spvl
                output_idx = output_idx + ninner
            end do
            return
        end if

        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        call build_flat_level_offsets(base_in, inner, ninner, input_offsets)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pair_idx = lowest_bracketing_level_index(plevi)
        bottom_pair_input_idx = input_offsets(bottom_pair_idx)
        bottom_pressure = plevi(bottom_idx)

        select case (intyp)
        case (1)
            output_idx = base_out + inner
            do k = 1, nlevo
                if (plevo_mb(k) > bottom_pressure) then
                    if (kxtrp == 0) then
                        dato_flat(output_idx) = spvl
                    else
                        dato_flat(output_idx) = interpolate_value_linear( &
                            dati_flat(bottom_pair_input_idx), dati_flat(bottom_pair_input_idx + ninner), &
                            plevo_mb(k), plevi(bottom_pair_idx), plevi(bottom_pair_idx + 1) &
                        )
                    end if
                else
                    kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                    if (kp < 1) then
                        dato_flat(output_idx) = spvl
                    else
                        input_idx = input_offsets(kp)
                        dato_flat(output_idx) = interpolate_value_linear( &
                            dati_flat(input_idx), dati_flat(input_idx + ninner), &
                            plevo_mb(k), plevi(kp), plevi(kp + 1) &
                        )
                    end if
                end if
                output_idx = output_idx + ninner
            end do
        case (2)
            output_idx = base_out + inner
            do k = 1, nlevo
                if (plevo_mb(k) > bottom_pressure) then
                    if (kxtrp == 0) then
                        dato_flat(output_idx) = spvl
                    else
                        dato_flat(output_idx) = interpolate_value_log( &
                            dati_flat(bottom_pair_input_idx), dati_flat(bottom_pair_input_idx + ninner), &
                            plevo_mb(k), plevi(bottom_pair_idx), plevi(bottom_pair_idx + 1) &
                        )
                    end if
                else
                    kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                    if (kp < 1) then
                        dato_flat(output_idx) = spvl
                    else
                        input_idx = input_offsets(kp)
                        dato_flat(output_idx) = interpolate_value_log( &
                            dati_flat(input_idx), dati_flat(input_idx + ninner), &
                            plevo_mb(k), plevi(kp), plevi(kp + 1) &
                        )
                    end if
                end if
                output_idx = output_idx + ninner
            end do
        case (3)
            block
                real(real64) :: plevi_loglog(nlevi), plevi_loglog_den(nlevi - 1)
                call build_double_log_levels(plevi, plevi_loglog)
                call build_double_log_interval_denominators(plevi_loglog, plevi_loglog_den)
                output_idx = base_out + inner
                do k = 1, nlevo
                    if (plevo_mb(k) > bottom_pressure) then
                        if (kxtrp == 0) then
                            dato_flat(output_idx) = spvl
                        else
                            dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                dati_flat(bottom_pair_input_idx), dati_flat(bottom_pair_input_idx + ninner), &
                                plevo_loglog(k), plevi_loglog(bottom_pair_idx), &
                                plevi_loglog_den(bottom_pair_idx) &
                            )
                        end if
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_flat(output_idx) = spvl
                        else
                            input_idx = input_offsets(kp)
                            dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                dati_flat(input_idx), dati_flat(input_idx + ninner), plevo_loglog(k), &
                                plevi_loglog(kp), plevi_loglog_den(kp) &
                            )
                        end if
                    end if
                    output_idx = output_idx + ninner
                end do
            end block
        case default
            output_idx = base_out + inner
            do k = 1, nlevo
                dato_flat(output_idx) = spvl
                output_idx = output_idx + ninner
            end do
        end select
    end subroutine interpolate_flat_column_nodes


    ! Flat-buffer helper for the C-order ECMWF hybrid->pressure path.
    subroutine interpolate_flat_column_ecmwf( &
        dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
        intyp, psfc, spvl, kxtrp, varflg, tbot, phis, base_in, base_out, inner, ninner, &
        nlevi, nlevo &
    )
        real(real64), intent(in) :: dati_flat(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: plevo_loglog(:)
        real(real64), intent(inout) :: dato_flat(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl, tbot, phis
        integer, intent(in) :: intyp, kxtrp, varflg, base_in, base_out, inner
        integer, intent(in) :: ninner, nlevi, nlevo

        real(real64) :: bottom_pressure, plevi(nlevi), psfc_mb
        integer :: bottom_idx, bottom_input_idx, input_idx, input_offsets(nlevi)
        integer :: k, kp, output_idx
        logical :: is_ascending

        if (psfc == spvl) then
            output_idx = base_out + inner
            do k = 1, nlevo
                dato_flat(output_idx) = spvl
                output_idx = output_idx + ninner
            end do
            return
        end if

        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        call build_flat_level_offsets(base_in, inner, ninner, input_offsets)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_input_idx = input_offsets(bottom_idx)
        bottom_pressure = plevi(bottom_idx)

        ! Keep the same ordered bracket as the plain node kernel while the
        ! target level stays inside the model column. Only plev > p_sfc should
        ! switch to the ECMWF below-ground formulation, matching GeoCAT's
        ! public Python semantics.

        select case (intyp)
        case (1)
            if (kxtrp == 0) then
                output_idx = base_out + inner
                do k = 1, nlevo
                    if (plevo_mb(k) > bottom_pressure) then
                        dato_flat(output_idx) = spvl
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_flat(output_idx) = spvl
                        else
                            input_idx = input_offsets(kp)
                            dato_flat(output_idx) = interpolate_value_linear( &
                                dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                plevo_mb(k), plevi(kp), plevi(kp + 1) &
                            )
                        end if
                    end if
                    output_idx = output_idx + ninner
                end do
            else
                select case (varflg)
                case (1)
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = extrapolate_temperature( &
                                dati_flat(bottom_input_idx), plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_linear( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                case (-1)
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = extrapolate_geopotential( &
                                tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_linear( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                case default
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = dati_flat(bottom_input_idx)
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_linear( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                end select
            end if
        case (2)
            if (kxtrp == 0) then
                output_idx = base_out + inner
                do k = 1, nlevo
                    if (plevo_mb(k) > bottom_pressure) then
                        dato_flat(output_idx) = spvl
                    else
                        kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                        if (kp < 1) then
                            dato_flat(output_idx) = spvl
                        else
                            input_idx = input_offsets(kp)
                            dato_flat(output_idx) = interpolate_value_log( &
                                dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                plevo_mb(k), plevi(kp), plevi(kp + 1) &
                            )
                        end if
                    end if
                    output_idx = output_idx + ninner
                end do
            else
                select case (varflg)
                case (1)
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = extrapolate_temperature( &
                                dati_flat(bottom_input_idx), plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_log( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                case (-1)
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = extrapolate_geopotential( &
                                tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                            )
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_log( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                case default
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = dati_flat(bottom_input_idx)
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_log( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), &
                                    plevo_mb(k), plevi(kp), plevi(kp + 1) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                end select
            end if
        case (3)
            block
                real(real64) :: plevi_loglog(nlevi), plevi_loglog_den(nlevi - 1)
                call build_double_log_levels(plevi, plevi_loglog)
                call build_double_log_interval_denominators(plevi_loglog, plevi_loglog_den)
                if (kxtrp == 0) then
                    output_idx = base_out + inner
                    do k = 1, nlevo
                        if (plevo_mb(k) > bottom_pressure) then
                            dato_flat(output_idx) = spvl
                        else
                            kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                            if (kp < 1) then
                                dato_flat(output_idx) = spvl
                            else
                                input_idx = input_offsets(kp)
                                dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                    dati_flat(input_idx), dati_flat(input_idx + ninner), plevo_loglog(k), &
                                    plevi_loglog(kp), plevi_loglog_den(kp) &
                                )
                            end if
                        end if
                        output_idx = output_idx + ninner
                    end do
                else
                    select case (varflg)
                    case (1)
                        output_idx = base_out + inner
                        do k = 1, nlevo
                            if (plevo_mb(k) > bottom_pressure) then
                                dato_flat(output_idx) = extrapolate_temperature( &
                                    dati_flat(bottom_input_idx), plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                                )
                            else
                                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                                if (kp < 1) then
                                    dato_flat(output_idx) = spvl
                                else
                                    input_idx = input_offsets(kp)
                                    dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                        dati_flat(input_idx), dati_flat(input_idx + ninner), plevo_loglog(k), &
                                        plevi_loglog(kp), plevi_loglog_den(kp) &
                                    )
                                end if
                            end if
                            output_idx = output_idx + ninner
                        end do
                    case (-1)
                        output_idx = base_out + inner
                        do k = 1, nlevo
                            if (plevo_mb(k) > bottom_pressure) then
                                dato_flat(output_idx) = extrapolate_geopotential( &
                                    tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                                )
                            else
                                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                                if (kp < 1) then
                                    dato_flat(output_idx) = spvl
                                else
                                    input_idx = input_offsets(kp)
                                    dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                        dati_flat(input_idx), dati_flat(input_idx + ninner), plevo_loglog(k), &
                                        plevi_loglog(kp), plevi_loglog_den(kp) &
                                    )
                                end if
                            end if
                            output_idx = output_idx + ninner
                        end do
                    case default
                        output_idx = base_out + inner
                        do k = 1, nlevo
                            if (plevo_mb(k) > bottom_pressure) then
                                dato_flat(output_idx) = dati_flat(bottom_input_idx)
                            else
                                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                                if (kp < 1) then
                                    dato_flat(output_idx) = spvl
                                else
                                    input_idx = input_offsets(kp)
                                    dato_flat(output_idx) = interpolate_value_loglog_precomputed( &
                                        dati_flat(input_idx), dati_flat(input_idx + ninner), plevo_loglog(k), &
                                        plevi_loglog(kp), plevi_loglog_den(kp) &
                                    )
                                end if
                            end if
                            output_idx = output_idx + ninner
                        end do
                    end select
                end if
            end block
        case default
            output_idx = base_out + inner
            do k = 1, nlevo
                dato_flat(output_idx) = spvl
                output_idx = output_idx + ninner
            end do
        end select
    end subroutine interpolate_flat_column_ecmwf


    ! ECMWF below-ground temperature extrapolation for one target pressure.
    pure real(real64) function extrapolate_temperature( &
        tbot, plev_bottom, plev_out, psfc_mb, phi_sfc &
    ) result(value)
        real(real64), intent(in) :: tbot, plev_bottom, plev_out, psfc_mb, phi_sfc

        real(real64) :: alnp, hgt, t0, tplat, tprime0, tstar

        tstar = tbot * (1.0_real64 + alpha * ((psfc_mb / plev_bottom) - 1.0_real64))
        hgt = phi_sfc * ginv

        if (hgt < 2000.0_real64) then
            alnp = alpha * log(plev_out / psfc_mb)
        else
            t0 = tstar + 0.0065_real64 * hgt
            tplat = min(t0, 298.0_real64)

            if (hgt <= 2500.0_real64) then
                tprime0 = 0.002_real64 * &
                    ((2500.0_real64 - hgt) * t0 + (hgt - 2000.0_real64) * tplat)
            else
                tprime0 = tplat
            end if

            if (tprime0 < tstar) then
                alnp = 0.0_real64
            else
                alnp = rd * (tprime0 - tstar) / phi_sfc * log(plev_out / psfc_mb)
            end if
        end if

        value = tstar * &
            (1.0_real64 + alnp + 0.5_real64 * alnp**2 + (1.0_real64 / 6.0_real64) * alnp**3)
    end function extrapolate_temperature


    ! ECMWF below-ground geopotential extrapolation for one target pressure.
    pure real(real64) function extrapolate_geopotential( &
        tbot, plev_bottom, plev_out, psfc_mb, phi_sfc &
    ) result(value)
        real(real64), intent(in) :: tbot, plev_bottom, plev_out, psfc_mb, phi_sfc

        real(real64) :: alnp, alph, hgt, t0, tstar

        hgt = phi_sfc * ginv
        tstar = tbot * (1.0_real64 + alpha * ((psfc_mb / plev_bottom) - 1.0_real64))
        t0 = tstar + 0.0065_real64 * hgt

        if (tstar <= 290.5_real64 .and. t0 > 290.5_real64) then
            alph = rd / phi_sfc * (290.5_real64 - tstar)
        else if (tstar > 290.5_real64 .and. t0 > 290.5_real64) then
            alph = 0.0_real64
            tstar = 0.5_real64 * (290.5_real64 + tstar)
        else
            alph = alpha
        end if

        if (tstar < 255.0_real64) then
            tstar = 0.5_real64 * (tstar + 255.0_real64)
        end if

        alnp = alph * log(plev_out / psfc_mb)
        value = hgt - rd * tstar * ginv * log(plev_out / psfc_mb) * &
            (1.0_real64 + 0.5_real64 * alnp + (1.0_real64 / 6.0_real64) * alnp**2)
    end function extrapolate_geopotential

end module vinth2p_kernels_core


! QUICK REFERENCE
! PURPOSE
!    GENERIC COLUMN-MAJOR ENTRY POINT FOR HYBRID-TO-PRESSURE
!    INTERPOLATION WITHOUT ECMWF BELOW-GROUND EXTRAPOLATION.
!
! EXPECTED INPUT SHAPES
!    DATI(NLEVI,NCOL)   - INPUT FIELD ON HYBRID LEVELS
!    DATO(NLEVO,NCOL)   - OUTPUT FIELD ON PRESSURE LEVELS
!    HBCOFA(NLEVI)      - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)      - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)       - TARGET PRESSURE LEVELS
!    PSFC(NCOL)         - SURFACE PRESSURE FOR EACH COLUMN
!
! UNITS
!    P0     - PA ON ENTRY, CONVERTED TO MB INTERNALLY
!    PLEVO  - PA ON ENTRY, CONVERTED TO MB INTERNALLY
!    PSFC   - PA
!    DATI   - CALLER'S FIELD UNITS, PRESERVED IN DATO
!
! FLAGS
!    INTYP  - 1: LINEAR, 2: LOG, 3: LOG-LOG
!    KXTRP  - 0: RETURN SPVL BELOW GROUND
!             1: EXTRAPOLATE USING THE LOWEST BRACKETING PAIR
!    SPVL   - MISSING VALUE SENTINEL
!
! OUTPUT
!    DATO HOLDS THE INTERPOLATED PRESSURE-LEVEL FIELD FOR EACH COLUMN.
subroutine dvinth2p_nodes_pa( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo)
    use vinth2p_kernels_core, only : &
        real64, build_double_log_levels, convert_levels_to_mb, interpolate_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl

    real(real64) :: plevo_mb(nlevo), plevo_loglog(nlevo), p0_mb
    integer :: j

    ! This is the generic column-major entry point used by the fallback Python
    ! bridge and by the explicit output-buffer wrapper below.
    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)
    if (intyp == 3) call build_double_log_levels(plevo_mb, plevo_loglog)

    do j = 1, ncol
        call interpolate_column_nodes( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
            intyp, psfc(j), spvl, kxtrp &
        )
    end do
end subroutine dvinth2p_nodes_pa


! QUICK REFERENCE
! PURPOSE
!    GENERIC COLUMN-MAJOR ENTRY POINT FOR HYBRID-TO-PRESSURE
!    INTERPOLATION WITH ECMWF BELOW-GROUND EXTRAPOLATION.
!
! EXPECTED INPUT SHAPES
!    DATI(NLEVI,NCOL)   - INPUT FIELD ON HYBRID LEVELS
!    DATO(NLEVO,NCOL)   - OUTPUT FIELD ON PRESSURE LEVELS
!    HBCOFA(NLEVI)      - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)      - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)       - TARGET PRESSURE LEVELS
!    PSFC(NCOL)         - SURFACE PRESSURE FOR EACH COLUMN
!    TBOT(NCOL)         - TEMPERATURE AT THE LOWEST MODEL LEVEL
!    PHIS(NCOL)         - SURFACE GEOPOTENTIAL
!
! UNITS
!    P0     - PA ON ENTRY, CONVERTED TO MB INTERNALLY
!    PLEVO  - PA ON ENTRY, CONVERTED TO MB INTERNALLY
!    PSFC   - PA
!    TBOT   - K
!    PHIS   - M2 S-2
!
! FLAGS
!    INTYP  - 1: LINEAR, 2: LOG, 3: LOG-LOG
!    KXTRP  - 0: RETURN SPVL BELOW GROUND
!             1: APPLY ECMWF BELOW-GROUND EXTRAPOLATION
!    VARFLG - +1: TEMPERATURE, -1: GEOPOTENTIAL, 0: OTHER
!             VARIABLES USE THE LOWEST MODEL LEVEL BELOW GROUND
!    SPVL   - MISSING VALUE SENTINEL
!
! OUTPUT
!    DATO HOLDS THE INTERPOLATED OR EXTRAPOLATED PRESSURE-LEVEL FIELD
!    FOR EACH COLUMN.
subroutine dvinth2p_ecmwf_nodes_pa( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo, varflg, tbot, phis)
    use vinth2p_kernels_core, only : &
        real64, build_double_log_levels, convert_levels_to_mb, interpolate_column_ecmwf
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl
    real(real64), intent(in) :: tbot(ncol), phis(ncol)

    real(real64) :: plevo_mb(nlevo), plevo_loglog(nlevo), p0_mb
    integer :: j

    ! This is the generic column-major ECMWF entry point used by the fallback
    ! Python bridge and by the explicit output-buffer wrapper below.
    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)
    if (intyp == 3) call build_double_log_levels(plevo_mb, plevo_loglog)

    do j = 1, ncol
        call interpolate_column_ecmwf( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
            intyp, psfc(j), spvl, kxtrp, varflg, tbot(j), phis(j) &
        )
    end do
end subroutine dvinth2p_ecmwf_nodes_pa


! QUICK REFERENCE
! PURPOSE
!    THIN F2PY-FRIENDLY WRAPPER THAT WRITES INTO A CALLER-PROVIDED
!    COLUMN-MAJOR OUTPUT BUFFER TO AVOID RETURN-ARRAY ALLOCATION.
!
! EXPECTED INPUT SHAPES
!    DATI(NLEVI,NCOL)   - INPUT FIELD ON HYBRID LEVELS
!    DATO(NLEVO,NCOL)   - PREALLOCATED OUTPUT BUFFER
!    HBCOFA(NLEVI)      - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)      - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)       - TARGET PRESSURE LEVELS
!    PSFC(NCOL)         - SURFACE PRESSURE FOR EACH COLUMN
subroutine dvinth2p_nodes_pa_into( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo)
    use vinth2p_kernels_core, only : real64
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(inout) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl

    call dvinth2p_nodes_pa( &
        dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
        nlevi, ncol, nlevo &
    )
end subroutine dvinth2p_nodes_pa_into


! QUICK REFERENCE
! PURPOSE
!    THIN F2PY-FRIENDLY WRAPPER THAT WRITES INTO A CALLER-PROVIDED
!    COLUMN-MAJOR OUTPUT BUFFER FOR THE ECMWF EXTRAPOLATION PATH.
!
! EXPECTED INPUT SHAPES
!    DATI(NLEVI,NCOL)   - INPUT FIELD ON HYBRID LEVELS
!    DATO(NLEVO,NCOL)   - PREALLOCATED OUTPUT BUFFER
!    HBCOFA(NLEVI)      - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)      - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)       - TARGET PRESSURE LEVELS
!    PSFC(NCOL)         - SURFACE PRESSURE FOR EACH COLUMN
!    TBOT(NCOL)         - TEMPERATURE AT THE LOWEST MODEL LEVEL
!    PHIS(NCOL)         - SURFACE GEOPOTENTIAL
subroutine dvinth2p_ecmwf_nodes_pa_into( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo, varflg, tbot, phis)
    use vinth2p_kernels_core, only : real64
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(inout) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl
    real(real64), intent(in) :: tbot(ncol), phis(ncol)

    call dvinth2p_ecmwf_nodes_pa( &
        dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
        nlevi, ncol, nlevo, varflg, tbot, phis &
    )
end subroutine dvinth2p_ecmwf_nodes_pa_into


! QUICK REFERENCE
! PURPOSE
!    GENERIC ENTRY POINT FOR HYBRID LAYER-THICKNESS CALCULATION USING A
!    COLUMN-MAJOR [NLEV-1, NCOL] OUTPUT BUFFER.
!
! EXPECTED INPUT SHAPES
!    PSFC(NCOL)                 - SURFACE PRESSURE FOR EACH COLUMN
!    DPH(NLEVO,NCOL)            - OUTPUT PRESSURE THICKNESS
!    HBCOFA(NLEV)               - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEV)               - HYBRID "B" COEFFICIENTS
!
! UNITS
!    P0     - PA
!    PSFC   - PA
!    DPH    - PA
subroutine ddelta_pressure_hybrid_pa(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo)
    use vinth2p_kernels_core, only : real64, compute_delta_pressure_columns
    implicit none

    integer, intent(in) :: ncol, nlev, nlevo
    real(real64), intent(in) :: psfc(ncol), hbcofa(nlev), hbcofb(nlev), p0
    real(real64), intent(out) :: dph(nlevo, ncol)

    call compute_delta_pressure_columns(psfc, hbcofa, hbcofb, p0, dph)
end subroutine ddelta_pressure_hybrid_pa


! QUICK REFERENCE
! PURPOSE
!    THIN F2PY-FRIENDLY WRAPPER THAT WRITES HYBRID DELTA-PRESSURE RESULTS
!    INTO A CALLER-PROVIDED COLUMN-MAJOR OUTPUT BUFFER.
subroutine ddelta_pressure_hybrid_pa_into(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo)
    use vinth2p_kernels_core, only : real64
    implicit none

    integer, intent(in) :: ncol, nlev, nlevo
    real(real64), intent(in) :: psfc(ncol), hbcofa(nlev), hbcofb(nlev), p0
    real(real64), intent(inout) :: dph(nlevo, ncol)

    call ddelta_pressure_hybrid_pa(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo)
end subroutine ddelta_pressure_hybrid_pa_into


! QUICK REFERENCE
! PURPOSE
!    GENERIC ENTRY POINT FOR HYBRID PRESSURE CALCULATION USING A
!    COLUMN-MAJOR [NLEV, NCOL] OUTPUT BUFFER.
!
! EXPECTED INPUT SHAPES
!    PSFC(NCOL)                 - SURFACE PRESSURE FOR EACH COLUMN
!    PRESSURE(NLEV,NCOL)        - OUTPUT PRESSURE AT EACH HYBRID LEVEL
!    HBCOFA(NLEV)               - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEV)               - HYBRID "B" COEFFICIENTS
!
! UNITS
!    P0        - PA
!    PSFC      - PA
!    PRESSURE  - PA
subroutine dpressure_at_hybrid_levels_pa(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev)
    use vinth2p_kernels_core, only : real64, compute_pressure_at_hybrid_levels_columns
    implicit none

    integer, intent(in) :: ncol, nlev
    real(real64), intent(in) :: psfc(ncol), hbcofa(nlev), hbcofb(nlev), p0
    real(real64), intent(out) :: pressure(nlev, ncol)

    call compute_pressure_at_hybrid_levels_columns(psfc, hbcofa, hbcofb, p0, pressure)
end subroutine dpressure_at_hybrid_levels_pa


! QUICK REFERENCE
! PURPOSE
!    THIN F2PY-FRIENDLY WRAPPER THAT WRITES HYBRID PRESSURES INTO A
!    CALLER-PROVIDED COLUMN-MAJOR OUTPUT BUFFER.
subroutine dpressure_at_hybrid_levels_pa_into(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev)
    use vinth2p_kernels_core, only : real64
    implicit none

    integer, intent(in) :: ncol, nlev
    real(real64), intent(in) :: psfc(ncol), hbcofa(nlev), hbcofb(nlev), p0
    real(real64), intent(inout) :: pressure(nlev, ncol)

    call dpressure_at_hybrid_levels_pa(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev)
end subroutine dpressure_at_hybrid_levels_pa_into


! QUICK REFERENCE
! PURPOSE
!    FAST C-ORDER ENTRY POINT FOR HYBRID-LEVEL GEOPOTENTIAL HEIGHT.
!    INPUT AND OUTPUT ARRAYS STAY IN FLAT NUMPY C-ORDER BUFFERS SO PYTHON
!    DOES NOT NEED TO PACK THEM INTO COLUMN-MAJOR STORAGE.
subroutine dgeopotential_height_hybrid_corder_pa_into( &
    temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi, p0, &
    nouter, nlev, ninner)
    use vinth2p_kernels_core, only : real64, compute_geopotential_height_flat_column
    implicit none

    integer, intent(in) :: nouter, nlev, ninner
    real(real64), intent(in) :: temp_flat(nouter * nlev * ninner)
    real(real64), intent(in) :: q_flat(nouter * nlev * ninner)
    real(real64), intent(inout) :: z3_flat(nouter * nlev * ninner)
    real(real64), intent(in) :: psfc(nouter * ninner), phis(nouter * ninner)
    real(real64), intent(in) :: hyai(nlev + 1), hybi(nlev + 1)
    real(real64), intent(in) :: p0

    integer :: base_in, base_out, col_idx, inner, outer

    do outer = 0, nouter - 1
        base_in = outer * nlev * ninner
        base_out = outer * nlev * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call compute_geopotential_height_flat_column( &
                temp_flat, q_flat, z3_flat, psfc(col_idx), phis(col_idx), hyai, hybi, p0, &
                base_in, base_out, inner, ninner, nlev &
            )
        end do
    end do
end subroutine dgeopotential_height_hybrid_corder_pa_into


! QUICK REFERENCE
! PURPOSE
!    GENERIC COLUMN-MAJOR ENTRY POINT FOR INTERPOLATION FROM SOURCE
!    SIGMA LEVELS TO TARGET HYBRID-IMPLIED SIGMA LEVELS.
!
! EXPECTED INPUT SHAPES
!    DATI(NLEVI,NCOL)   - INPUT FIELD ON SIGMA LEVELS
!    DATO(NLEVO,NCOL)   - OUTPUT FIELD ON TARGET HYBRID LEVELS
!    SIGMAI(NLEVI)      - SOURCE SIGMA COORDINATES
!    SIGMAO(NLEVO,NCOL) - TARGET SIGMA COORDINATES FOR EACH COLUMN
!
! FLAGS
!    INTYP  - 1: LINEAR, 2: LOG
!    SPVL   - MISSING VALUE SENTINEL
subroutine dsigma2hybrid_nodes( &
    dati, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo)
    use vinth2p_kernels_core, only : real64, interpolate_column_sigma
    implicit none

    integer, intent(in) :: intyp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: sigmai(nlevi), sigmao(nlevo, ncol), spvl

    integer :: j

    do j = 1, ncol
        call interpolate_column_sigma( &
            dati(:, j), dato(:, j), sigmai, sigmao(:, j), intyp, spvl &
        )
    end do
end subroutine dsigma2hybrid_nodes


! QUICK REFERENCE
! PURPOSE
!    THIN F2PY-FRIENDLY WRAPPER THAT WRITES SIGMA->HYBRID RESULTS INTO
!    A CALLER-PROVIDED COLUMN-MAJOR OUTPUT BUFFER.
subroutine dsigma2hybrid_nodes_into( &
    dati, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo)
    use vinth2p_kernels_core, only : real64
    implicit none

    integer, intent(in) :: intyp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(inout) :: dato(nlevo, ncol)
    real(real64), intent(in) :: sigmai(nlevi), sigmao(nlevo, ncol), spvl

    call dsigma2hybrid_nodes( &
        dati, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo &
    )
end subroutine dsigma2hybrid_nodes_into


! QUICK REFERENCE
! PURPOSE
!    C-ORDER FAST PATH FOR SIGMA->HYBRID INTERPOLATION. INPUT DATA IS
!    PASSED AS A FLAT [NOUTER, NLEVI, NINNER] BUFFER, AND THE TARGET
!    HYBRID-IMPLIED SIGMA VALUES ARE COMPUTED INSIDE FORTRAN FROM
!    PS, HYAM, HYBM, AND P0.
subroutine dsigma2hybrid_nodes_corder_into( &
    dati_flat, dato_flat, sigmai, hyam, hybm, p0, psfc, intyp, spvl, &
    nouter, nlevi, ninner, nlevo)
    use vinth2p_kernels_core, only : &
        real64, interpolate_flat_column_sigma
    implicit none

    integer, intent(in) :: intyp, nouter, nlevi, ninner, nlevo
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: sigmai(nlevi), hyam(nlevo), hybm(nlevo)
    real(real64), intent(in) :: p0, psfc(nouter * ninner), spvl

    integer :: base_in, base_out, col_idx, inner, outer

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call interpolate_flat_column_sigma( &
                dati_flat, dato_flat, sigmai, hyam, hybm, p0, psfc(col_idx), intyp, spvl, &
                base_in, base_out, inner, ninner, nlevi, nlevo &
            )
        end do
    end do
end subroutine dsigma2hybrid_nodes_corder_into


! QUICK REFERENCE
! PURPOSE
!    FAST PATH FOR NUMPY C-ORDER INPUTS. DATA IS PASSED AS A FLAT 1-D
!    BUFFER REPRESENTING [NOUTER, NLEVI, NINNER] IN C ORDER SO PYTHON
!    DOES NOT NEED TO TRANSPOSE INTO COLUMN-MAJOR STORAGE.
!
! EXPECTED INPUT SHAPES
!    DATI_FLAT(NOUTER*NLEVI*NINNER) - FLAT INPUT BUFFER
!    DATO_FLAT(NOUTER*NLEVO*NINNER) - FLAT OUTPUT BUFFER
!    HBCOFA(NLEVI)                  - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)                  - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)                   - TARGET PRESSURE LEVELS
!    PSFC(NOUTER*NINNER)            - SURFACE PRESSURE FOR EACH COLUMN
!
! INDEXING NOTE
!    FOR A GIVEN OUTER INDEX AND INNER INDEX, ONE PHYSICAL COLUMN IS
!    FORMED BY WALKING THROUGH THE NLEVI LEVELS IN THE FLAT BUFFER.
subroutine dvinth2p_nodes_corder_pa_into( &
    dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nouter, nlevi, ninner, nlevo)
    use vinth2p_kernels_core, only : &
        real64, build_double_log_levels, convert_levels_to_mb, interpolate_flat_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl

    ! This flat-array entry point exists purely to avoid Python-side
    ! transpose/packing when the caller already has NumPy C-order data with the
    ! interpolation axis in place.
    real(real64) :: plevo_mb(nlevo), plevo_loglog(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)
    if (intyp == 3) call build_double_log_levels(plevo_mb, plevo_loglog)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call interpolate_flat_column_nodes( &
                dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
                intyp, psfc(col_idx), spvl, kxtrp, &
                base_in, base_out, inner, ninner, nlevi, nlevo &
            )
        end do
    end do
end subroutine dvinth2p_nodes_corder_pa_into


! QUICK REFERENCE
! PURPOSE
!    FAST C-ORDER ECMWF ENTRY POINT. THIS MATCHES THE GENERIC ECMWF
!    KERNEL BUT KEEPS INPUT AND OUTPUT IN FLAT NUMPY C-ORDER BUFFERS.
!
! EXPECTED INPUT SHAPES
!    DATI_FLAT(NOUTER*NLEVI*NINNER) - FLAT INPUT BUFFER
!    DATO_FLAT(NOUTER*NLEVO*NINNER) - FLAT OUTPUT BUFFER
!    HBCOFA(NLEVI)                  - HYBRID "A" COEFFICIENTS
!    HBCOFB(NLEVI)                  - HYBRID "B" COEFFICIENTS
!    PLEVO(NLEVO)                   - TARGET PRESSURE LEVELS
!    PSFC(NOUTER*NINNER)            - SURFACE PRESSURE FOR EACH COLUMN
!    TBOT(NOUTER*NINNER)            - TEMPERATURE AT THE LOWEST MODEL
!                                     LEVEL
!    PHIS(NOUTER*NINNER)            - SURFACE GEOPOTENTIAL
!
! FLAGS
!    VARFLG - +1: TEMPERATURE, -1: GEOPOTENTIAL, 0: OTHER
!             VARIABLES USE THE LOWEST MODEL LEVEL BELOW GROUND
subroutine dvinth2p_ecmwf_nodes_corder_pa_into( &
    dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nouter, nlevi, ninner, nlevo, varflg, tbot, phis)
    use vinth2p_kernels_core, only : &
        real64, build_double_log_levels, convert_levels_to_mb, interpolate_flat_column_ecmwf
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo, varflg
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl
    real(real64), intent(in) :: tbot(nouter * ninner), phis(nouter * ninner)

    ! This flat-array ECMWF entry point mirrors the generic column-major kernel
    ! but keeps the data in its original NumPy C-order layout.
    real(real64) :: plevo_mb(nlevo), plevo_loglog(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)
    if (intyp == 3) call build_double_log_levels(plevo_mb, plevo_loglog)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call interpolate_flat_column_ecmwf( &
                dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, plevo_loglog, &
                intyp, psfc(col_idx), spvl, kxtrp, varflg, tbot(col_idx), phis(col_idx), &
                base_in, base_out, inner, ninner, nlevi, nlevo &
            )
        end do
    end do
end subroutine dvinth2p_ecmwf_nodes_corder_pa_into


! QUICK REFERENCE
! PURPOSE
!    DIRECT C-ABI WRAPPERS FOR THE PUBLIC `vinth2p_kernels` ENTRY POINTS.
!    THESE KEEP THE SCIENTIFIC KERNELS ABOVE UNCHANGED WHILE EXPOSING THEM
!    THROUGH THE DIRECT `Fortran/C + Meson` BUILD PATH USED BY SKYBORN.
subroutine dvinth2p_nodes_pa_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                               nlevi, ncol, nlevo) bind(C, name="dvinth2p_nodes_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])

    call dvinth2p_nodes_pa(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                           nlevi, ncol, nlevo)
end subroutine dvinth2p_nodes_pa_c


subroutine dvinth2p_nodes_pa_into_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                    nlevi, ncol, nlevo) bind(C, name="dvinth2p_nodes_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])

    call dvinth2p_nodes_pa_into(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                                nlevi, ncol, nlevo)
end subroutine dvinth2p_nodes_pa_into_c


subroutine dvinth2p_ecmwf_nodes_pa_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                     nlevi, ncol, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(tbot, tbot_f, [ncol])
    call c_f_pointer(phis, phis_f, [ncol])

    call dvinth2p_ecmwf_nodes_pa(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, kxtrp, &
                                 nlevi, ncol, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_pa_c


subroutine dvinth2p_ecmwf_nodes_pa_into_c(dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
                                          nlevi, ncol, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati, dato, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(tbot, tbot_f, [ncol])
    call c_f_pointer(phis, phis_f, [ncol])

    call dvinth2p_ecmwf_nodes_pa_into(dati_f, dato_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, &
                                      kxtrp, nlevi, ncol, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_pa_into_c


subroutine ddelta_pressure_hybrid_pa_c(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo) &
        bind(C, name="ddelta_pressure_hybrid_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, dph, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev, nlevo
    real(real64), pointer :: psfc_f(:), dph_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(dph, dph_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call ddelta_pressure_hybrid_pa(psfc_f, dph_f, hbcofa_f, hbcofb_f, p0, ncol, nlev, nlevo)
end subroutine ddelta_pressure_hybrid_pa_c


subroutine ddelta_pressure_hybrid_pa_into_c(psfc, dph, hbcofa, hbcofb, p0, ncol, nlev, nlevo) &
        bind(C, name="ddelta_pressure_hybrid_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, dph, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev, nlevo
    real(real64), pointer :: psfc_f(:), dph_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(dph, dph_f, [nlevo, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call ddelta_pressure_hybrid_pa_into(psfc_f, dph_f, hbcofa_f, hbcofb_f, p0, ncol, nlev, nlevo)
end subroutine ddelta_pressure_hybrid_pa_into_c


subroutine dpressure_at_hybrid_levels_pa_c(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev) &
        bind(C, name="dpressure_at_hybrid_levels_pa_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, pressure, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev
    real(real64), pointer :: psfc_f(:), pressure_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(pressure, pressure_f, [nlev, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call dpressure_at_hybrid_levels_pa(psfc_f, pressure_f, hbcofa_f, hbcofb_f, p0, ncol, nlev)
end subroutine dpressure_at_hybrid_levels_pa_c


subroutine dpressure_at_hybrid_levels_pa_into_c(psfc, pressure, hbcofa, hbcofb, p0, ncol, nlev) &
        bind(C, name="dpressure_at_hybrid_levels_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: psfc, pressure, hbcofa, hbcofb
    real(c_double), value :: p0
    integer(c_int), value :: ncol, nlev
    real(real64), pointer :: psfc_f(:), pressure_f(:, :), hbcofa_f(:), hbcofb_f(:)

    call c_f_pointer(psfc, psfc_f, [ncol])
    call c_f_pointer(pressure, pressure_f, [nlev, ncol])
    call c_f_pointer(hbcofa, hbcofa_f, [nlev])
    call c_f_pointer(hbcofb, hbcofb_f, [nlev])

    call dpressure_at_hybrid_levels_pa_into(psfc_f, pressure_f, hbcofa_f, hbcofb_f, p0, ncol, nlev)
end subroutine dpressure_at_hybrid_levels_pa_into_c


subroutine dgeopotential_height_hybrid_corder_pa_into_c(temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi, p0, &
                                                        nouter, nlev, ninner) &
        bind(C, name="dgeopotential_height_hybrid_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: temp_flat, q_flat, z3_flat, psfc, phis, hyai, hybi
    real(c_double), value :: p0
    integer(c_int), value :: nouter, nlev, ninner
    real(real64), pointer :: temp_flat_f(:), q_flat_f(:), z3_flat_f(:), psfc_f(:), phis_f(:), hyai_f(:), hybi_f(:)

    call c_f_pointer(temp_flat, temp_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(q_flat, q_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(z3_flat, z3_flat_f, [nouter * nlev * ninner])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])
    call c_f_pointer(phis, phis_f, [nouter * ninner])
    call c_f_pointer(hyai, hyai_f, [nlev + 1])
    call c_f_pointer(hybi, hybi_f, [nlev + 1])

    call dgeopotential_height_hybrid_corder_pa_into(temp_flat_f, q_flat_f, z3_flat_f, psfc_f, phis_f, hyai_f, hybi_f, &
                                                    p0, nouter, nlev, ninner)
end subroutine dgeopotential_height_hybrid_corder_pa_into_c


subroutine dsigma2hybrid_nodes_c(dsigmai, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dsigmai, dato, sigmai, sigmao
    real(c_double), value :: spvl
    integer(c_int), value :: intyp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), sigmai_f(:), sigmao_f(:, :)

    call c_f_pointer(dsigmai, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(sigmao, sigmao_f, [nlevo, ncol])

    call dsigma2hybrid_nodes(dati_f, dato_f, sigmai_f, sigmao_f, intyp, spvl, nlevi, ncol, nlevo)
end subroutine dsigma2hybrid_nodes_c


subroutine dsigma2hybrid_nodes_into_c(dsigmai, dato, sigmai, sigmao, intyp, spvl, nlevi, ncol, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dsigmai, dato, sigmai, sigmao
    real(c_double), value :: spvl
    integer(c_int), value :: intyp, nlevi, ncol, nlevo
    real(real64), pointer :: dati_f(:, :), dato_f(:, :), sigmai_f(:), sigmao_f(:, :)

    call c_f_pointer(dsigmai, dati_f, [nlevi, ncol])
    call c_f_pointer(dato, dato_f, [nlevo, ncol])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(sigmao, sigmao_f, [nlevo, ncol])

    call dsigma2hybrid_nodes_into(dati_f, dato_f, sigmai_f, sigmao_f, intyp, spvl, nlevi, ncol, nlevo)
end subroutine dsigma2hybrid_nodes_into_c


subroutine dsigma2hybrid_nodes_corder_into_c(dati_flat, dato_flat, sigmai, hyam, hybm, p0, psfc, intyp, spvl, &
                                             nouter, nlevi, ninner, nlevo) &
        bind(C, name="dsigma2hybrid_nodes_corder_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, sigmai, hyam, hybm, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, nouter, nlevi, ninner, nlevo
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), sigmai_f(:), hyam_f(:), hybm_f(:), psfc_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(sigmai, sigmai_f, [nlevi])
    call c_f_pointer(hyam, hyam_f, [nlevo])
    call c_f_pointer(hybm, hybm_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])

    call dsigma2hybrid_nodes_corder_into(dati_flat_f, dato_flat_f, sigmai_f, hyam_f, hybm_f, p0, psfc_f, intyp, spvl, &
                                         nouter, nlevi, ninner, nlevo)
end subroutine dsigma2hybrid_nodes_corder_into_c


subroutine dvinth2p_nodes_corder_pa_into_c(dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, &
                                           kxtrp, nouter, nlevi, ninner, nlevo) &
        bind(C, name="dvinth2p_nodes_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, hbcofa, hbcofb, plevo, psfc
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nouter, nlevi, ninner, nlevo
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])

    call dvinth2p_nodes_corder_pa_into(dati_flat_f, dato_flat_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, spvl, &
                                       kxtrp, nouter, nlevi, ninner, nlevo)
end subroutine dvinth2p_nodes_corder_pa_into_c


subroutine dvinth2p_ecmwf_nodes_corder_pa_into_c(dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, &
                                                 kxtrp, nouter, nlevi, ninner, nlevo, varflg, tbot, phis) &
        bind(C, name="dvinth2p_ecmwf_nodes_corder_pa_into_c")
    use iso_c_binding
    use vinth2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value :: dati_flat, dato_flat, hbcofa, hbcofb, plevo, psfc, tbot, phis
    real(c_double), value :: p0, spvl
    integer(c_int), value :: intyp, kxtrp, nouter, nlevi, ninner, nlevo, varflg
    real(real64), pointer :: dati_flat_f(:), dato_flat_f(:), hbcofa_f(:), hbcofb_f(:), plevo_f(:), psfc_f(:)
    real(real64), pointer :: tbot_f(:), phis_f(:)

    call c_f_pointer(dati_flat, dati_flat_f, [nouter * nlevi * ninner])
    call c_f_pointer(dato_flat, dato_flat_f, [nouter * nlevo * ninner])
    call c_f_pointer(hbcofa, hbcofa_f, [nlevi])
    call c_f_pointer(hbcofb, hbcofb_f, [nlevi])
    call c_f_pointer(plevo, plevo_f, [nlevo])
    call c_f_pointer(psfc, psfc_f, [nouter * ninner])
    call c_f_pointer(tbot, tbot_f, [nouter * ninner])
    call c_f_pointer(phis, phis_f, [nouter * ninner])

    call dvinth2p_ecmwf_nodes_corder_pa_into(dati_flat_f, dato_flat_f, hbcofa_f, hbcofb_f, p0, plevo_f, intyp, psfc_f, &
                                             spvl, kxtrp, nouter, nlevi, ninner, nlevo, varflg, tbot_f, phis_f)
end subroutine dvinth2p_ecmwf_nodes_corder_pa_into_c
