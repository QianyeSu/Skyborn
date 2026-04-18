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


    ! Legacy double-log transform used by vinth2p "log-log" interpolation.
    pure real(real64) function double_log_pressure(pressure_mb) result(value)
        real(real64), intent(in) :: pressure_mb

        value = log(log(pressure_mb + 2.72_real64))
    end function double_log_pressure


    ! Apply the legacy vinth2p interpolation formula for one bracketing pair.
    pure real(real64) function interpolate_value( &
        lower_value, upper_value, plev_out, plev_lower, plev_upper, intyp, spvl &
    ) result(value)
        real(real64), intent(in) :: lower_value, upper_value
        real(real64), intent(in) :: plev_out, plev_lower, plev_upper, spvl
        integer, intent(in) :: intyp

        select case (intyp)
        case (1)
            value = lower_value + &
                (upper_value - lower_value) * &
                (plev_out - plev_lower) / &
                (plev_upper - plev_lower)
        case (2)
            value = lower_value + &
                (upper_value - lower_value) * &
                log(plev_out / plev_lower) / &
                log(plev_upper / plev_lower)
        case (3)
            value = lower_value + &
                (upper_value - lower_value) * &
                (double_log_pressure(plev_out) - double_log_pressure(plev_lower)) / &
                (double_log_pressure(plev_upper) - double_log_pressure(plev_lower))
        case default
            value = spvl
        end select
    end function interpolate_value


    ! Interpolate one vertical column from hybrid levels to requested pressure
    ! levels without ECMWF below-ground extrapolation.
    subroutine interpolate_column_nodes( &
        dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, intyp, psfc, spvl, kxtrp &
    )
        real(real64), intent(in) :: dati_col(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl
        real(real64), intent(out) :: dato_col(size(plevo_mb))
        integer, intent(in) :: intyp, kxtrp

        real(real64) :: bottom_pressure, plevi(size(dati_col)), psfc_mb
        integer :: bottom_idx, bottom_pair_idx, k, kp
        logical :: is_ascending

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

        do k = 1, size(plevo_mb)
            ! Below-ground requests either return the sentinel or use the same
            ! lowest bracketing pair as the legacy node kernel.
            if (plevo_mb(k) > bottom_pressure) then
                if (kxtrp == 0) then
                    dato_col(k) = spvl
                else
                    kp = bottom_pair_idx
                    dato_col(k) = interpolate_value( &
                        dati_col(kp), dati_col(kp + 1), &
                        plevo_mb(k), plevi(kp), plevi(kp + 1), intyp, spvl &
                    )
                end if
            else
                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                if (kp < 1) then
                    dato_col(k) = spvl
                else
                    dato_col(k) = interpolate_value( &
                        dati_col(kp), dati_col(kp + 1), &
                        plevo_mb(k), plevi(kp), plevi(kp + 1), intyp, spvl &
                    )
                end if
            end if
        end do
    end subroutine interpolate_column_nodes


    ! Interpolate one vertical column from hybrid levels to requested pressure
    ! levels with ECMWF below-ground extrapolation when requested.
    subroutine interpolate_column_ecmwf( &
        dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, intyp, psfc, spvl, kxtrp, &
        varflg, tbot, phis &
    )
        real(real64), intent(in) :: dati_col(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl, tbot, phis
        real(real64), intent(out) :: dato_col(size(plevo_mb))
        integer, intent(in) :: intyp, kxtrp, varflg

        real(real64) :: bottom_pressure, plevi(size(dati_col)), psfc_mb
        integer :: bottom_idx, k, kp
        logical :: is_ascending

        ! The Python bridge writes the exact sentinel value, so direct equality
        ! is intentional here.
        if (psfc == spvl) then
            dato_col = spvl
            return
        end if

        ! ECMWF extrapolation uses the lowest physical model level of the
        ! current column, but interior interpolation still follows the same
        ! bracketing logic as the plain node kernel.
        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)

        do k = 1, size(plevo_mb)
            if (plevo_mb(k) > bottom_pressure) then
                if (kxtrp == 0) then
                    dato_col(k) = spvl
                else
                    select case (varflg)
                    case (1)
                        dato_col(k) = extrapolate_temperature( &
                            dati_col(bottom_idx), plevi(bottom_idx), &
                            plevo_mb(k), psfc_mb, phis &
                        )
                    case (-1)
                        dato_col(k) = extrapolate_geopotential( &
                            tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                        )
                    case default
                        dato_col(k) = dati_col(bottom_idx)
                    end select
                end if
            else
                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                if (kp < 1) then
                    dato_col(k) = spvl
                else
                    dato_col(k) = interpolate_value( &
                        dati_col(kp), dati_col(kp + 1), &
                        plevo_mb(k), plevi(kp), plevi(kp + 1), intyp, spvl &
                    )
                end if
            end if
        end do
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
        dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, intyp, psfc, spvl, kxtrp, &
        base_in, base_out, inner, ninner, nlevi, nlevo &
    )
        real(real64), intent(in) :: dati_flat(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(inout) :: dato_flat(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl
        integer, intent(in) :: intyp, kxtrp, base_in, base_out, inner
        integer, intent(in) :: ninner, nlevi, nlevo

        real(real64) :: bottom_pressure, plevi(nlevi), psfc_mb
        integer :: bottom_idx, bottom_pair_idx, input_idx, k, kp, output_idx
        logical :: is_ascending

        if (psfc == spvl) then
            do k = 1, nlevo
                output_idx = base_out + (k - 1) * ninner + inner
                dato_flat(output_idx) = spvl
            end do
            return
        end if

        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pair_idx = lowest_bracketing_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)

        do k = 1, nlevo
            output_idx = base_out + (k - 1) * ninner + inner
            if (plevo_mb(k) > bottom_pressure) then
                if (kxtrp == 0) then
                    dato_flat(output_idx) = spvl
                else
                    input_idx = base_in + (bottom_pair_idx - 1) * ninner + inner
                    dato_flat(output_idx) = interpolate_value( &
                        dati_flat(input_idx), dati_flat(input_idx + ninner), &
                        plevo_mb(k), plevi(bottom_pair_idx), plevi(bottom_pair_idx + 1), &
                        intyp, spvl &
                    )
                end if
            else
                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                if (kp < 1) then
                    dato_flat(output_idx) = spvl
                else
                    input_idx = base_in + (kp - 1) * ninner + inner
                    dato_flat(output_idx) = interpolate_value( &
                        dati_flat(input_idx), dati_flat(input_idx + ninner), &
                        plevo_mb(k), plevi(kp), plevi(kp + 1), intyp, spvl &
                    )
                end if
            end if
        end do
    end subroutine interpolate_flat_column_nodes


    ! Flat-buffer helper for the C-order ECMWF hybrid->pressure path.
    subroutine interpolate_flat_column_ecmwf( &
        dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, intyp, psfc, spvl, kxtrp, &
        varflg, tbot, phis, base_in, base_out, inner, ninner, nlevi, nlevo &
    )
        real(real64), intent(in) :: dati_flat(:), hbcofa(:), hbcofb(:), plevo_mb(:)
        real(real64), intent(inout) :: dato_flat(:)
        real(real64), intent(in) :: p0_mb, psfc, spvl, tbot, phis
        integer, intent(in) :: intyp, kxtrp, varflg, base_in, base_out, inner
        integer, intent(in) :: ninner, nlevi, nlevo

        real(real64) :: bottom_pressure, plevi(nlevi), psfc_mb
        integer :: bottom_idx, input_idx, k, kp, output_idx
        logical :: is_ascending

        if (psfc == spvl) then
            do k = 1, nlevo
                output_idx = base_out + (k - 1) * ninner + inner
                dato_flat(output_idx) = spvl
            end do
            return
        end if

        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)

        do k = 1, nlevo
            output_idx = base_out + (k - 1) * ninner + inner
            if (plevo_mb(k) > bottom_pressure) then
                if (kxtrp == 0) then
                    dato_flat(output_idx) = spvl
                else
                    select case (varflg)
                    case (1)
                        input_idx = base_in + (bottom_idx - 1) * ninner + inner
                        dato_flat(output_idx) = extrapolate_temperature( &
                            dati_flat(input_idx), plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                        )
                    case (-1)
                        dato_flat(output_idx) = extrapolate_geopotential( &
                            tbot, plevi(bottom_idx), plevo_mb(k), psfc_mb, phis &
                        )
                    case default
                        input_idx = base_in + (bottom_idx - 1) * ninner + inner
                        dato_flat(output_idx) = dati_flat(input_idx)
                    end select
                end if
            else
                kp = locate_bracketing_level_ordered(plevo_mb(k), plevi, is_ascending)
                if (kp < 1) then
                    dato_flat(output_idx) = spvl
                else
                    input_idx = base_in + (kp - 1) * ninner + inner
                    dato_flat(output_idx) = interpolate_value( &
                        dati_flat(input_idx), dati_flat(input_idx + ninner), &
                        plevo_mb(k), plevi(kp), plevi(kp + 1), intyp, spvl &
                    )
                end if
            end if
        end do
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
        real64, convert_levels_to_mb, interpolate_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl

    real(real64) :: plevo_mb(nlevo), p0_mb
    integer :: j

    ! This is the generic column-major entry point used by the fallback Python
    ! bridge and by the explicit output-buffer wrapper below.
    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do j = 1, ncol
        call interpolate_column_nodes( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, &
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
        real64, convert_levels_to_mb, interpolate_column_ecmwf
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo, varflg
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl
    real(real64), intent(in) :: tbot(ncol), phis(ncol)

    real(real64) :: plevo_mb(nlevo), p0_mb
    integer :: j

    ! This is the generic column-major ECMWF entry point used by the fallback
    ! Python bridge and by the explicit output-buffer wrapper below.
    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do j = 1, ncol
        call interpolate_column_ecmwf( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, &
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
        real64, convert_levels_to_mb, interpolate_flat_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl

    ! This flat-array entry point exists purely to avoid Python-side
    ! transpose/packing when the caller already has NumPy C-order data with the
    ! interpolation axis in place.
    real(real64) :: plevo_mb(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call interpolate_flat_column_nodes( &
                dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, &
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
        real64, convert_levels_to_mb, interpolate_flat_column_ecmwf
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo, varflg
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl
    real(real64), intent(in) :: tbot(nouter * ninner), phis(nouter * ninner)

    ! This flat-array ECMWF entry point mirrors the generic column-major kernel
    ! but keeps the data in its original NumPy C-order layout.
    real(real64) :: plevo_mb(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            call interpolate_flat_column_ecmwf( &
                dati_flat, dato_flat, hbcofa, hbcofb, p0_mb, plevo_mb, &
                intyp, psfc(col_idx), spvl, kxtrp, varflg, tbot(col_idx), phis(col_idx), &
                base_in, base_out, inner, ninner, nlevi, nlevo &
            )
        end do
    end do
end subroutine dvinth2p_ecmwf_nodes_corder_pa_into
