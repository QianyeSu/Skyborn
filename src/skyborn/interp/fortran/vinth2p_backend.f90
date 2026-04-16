module vinth2p_backend_core
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none
    private

    real(real64), parameter :: rd = 287.04_real64
    real(real64), parameter :: ginv = 1.0_real64 / 9.80616_real64
    real(real64), parameter :: alpha = 0.0065_real64 * rd * ginv

    public :: real64
    public :: build_input_pressures
    public :: convert_levels_to_mb
    public :: extrapolate_geopotential
    public :: extrapolate_temperature
    public :: interpolate_column_ecmwf
    public :: interpolate_column_nodes
    public :: interpolate_value
    public :: is_below_lowest_model_level
    public :: levels_ascending
    public :: locate_bracketing_level_ordered
    public :: lowest_bracketing_level_index
    public :: lowest_model_level_index
    public :: locate_bracketing_level

contains

    pure subroutine convert_levels_to_mb(plevo, plevo_mb)
        real(real64), intent(in) :: plevo(:)
        real(real64), intent(out) :: plevo_mb(size(plevo))

        plevo_mb = plevo * 0.01_real64
    end subroutine convert_levels_to_mb


    pure subroutine build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        real(real64), intent(in) :: hbcofa(:), hbcofb(:)
        real(real64), intent(in) :: p0_mb, psfc_mb
        real(real64), intent(out) :: plevi(size(hbcofa))

        plevi = hbcofa * p0_mb + hbcofb * psfc_mb
    end subroutine build_input_pressures


    pure logical function levels_ascending(plevi) result(is_ascending)
        real(real64), intent(in) :: plevi(:)

        is_ascending = plevi(1) <= plevi(size(plevi))
    end function levels_ascending


    pure integer function lowest_model_level_index(plevi) result(idx)
        real(real64), intent(in) :: plevi(:)

        if (levels_ascending(plevi)) then
            idx = size(plevi)
        else
            idx = 1
        end if
    end function lowest_model_level_index


    pure integer function lowest_bracketing_level_index(plevi) result(idx)
        real(real64), intent(in) :: plevi(:)

        if (levels_ascending(plevi)) then
            idx = size(plevi) - 1
        else
            idx = 1
        end if
    end function lowest_bracketing_level_index


    pure logical function is_below_lowest_model_level(plev_out, plevi) result(is_below)
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)

        is_below = plev_out > plevi(lowest_model_level_index(plevi))
    end function is_below_lowest_model_level


    pure integer function locate_bracketing_level(plev_out, plevi) result(kp)
        real(real64), intent(in) :: plev_out
        real(real64), intent(in) :: plevi(:)

        kp = locate_bracketing_level_ordered(plev_out, plevi, levels_ascending(plevi))
    end function locate_bracketing_level


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


    pure real(real64) function double_log_pressure(pressure_mb) result(value)
        real(real64), intent(in) :: pressure_mb

        value = log(log(pressure_mb + 2.72_real64))
    end function double_log_pressure


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

        if (psfc == spvl) then
            dato_col = spvl
            return
        end if

        psfc_mb = psfc * 0.01_real64
        call build_input_pressures(hbcofa, hbcofb, p0_mb, psfc_mb, plevi)
        is_ascending = levels_ascending(plevi)
        bottom_idx = lowest_model_level_index(plevi)
        bottom_pair_idx = lowest_bracketing_level_index(plevi)
        bottom_pressure = plevi(bottom_idx)

        do k = 1, size(plevo_mb)
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

        if (psfc == spvl) then
            dato_col = spvl
            return
        end if

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

end module vinth2p_backend_core


subroutine dvinth2p_nodes_pa( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo)
    use vinth2p_backend_core, only : &
        real64, convert_levels_to_mb, interpolate_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nlevi, ncol, nlevo
    real(real64), intent(in) :: dati(nlevi, ncol)
    real(real64), intent(out) :: dato(nlevo, ncol)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(ncol), spvl

    real(real64) :: plevo_mb(nlevo), p0_mb
    integer :: j

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do j = 1, ncol
        call interpolate_column_nodes( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, &
            intyp, psfc(j), spvl, kxtrp &
        )
    end do
end subroutine dvinth2p_nodes_pa


subroutine dvinth2p_ecmwf_nodes_pa( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo, varflg, tbot, phis)
    use vinth2p_backend_core, only : &
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

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do j = 1, ncol
        call interpolate_column_ecmwf( &
            dati(:, j), dato(:, j), hbcofa, hbcofb, p0_mb, plevo_mb, &
            intyp, psfc(j), spvl, kxtrp, varflg, tbot(j), phis(j) &
        )
    end do
end subroutine dvinth2p_ecmwf_nodes_pa


subroutine dvinth2p_nodes_pa_into( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo)
    use vinth2p_backend_core, only : real64
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


subroutine dvinth2p_ecmwf_nodes_pa_into( &
    dati, dato, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nlevi, ncol, nlevo, varflg, tbot, phis)
    use vinth2p_backend_core, only : real64
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


subroutine dvinth2p_nodes_corder_pa_into( &
    dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nouter, nlevi, ninner, nlevo)
    use vinth2p_backend_core, only : real64, convert_levels_to_mb, interpolate_column_nodes
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl

    real(real64) :: dati_col(nlevi), dato_col(nlevo), plevo_mb(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, k, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            do k = 1, nlevi
                dati_col(k) = dati_flat(base_in + (k - 1) * ninner + inner)
            end do
            call interpolate_column_nodes( &
                dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, &
                intyp, psfc(col_idx), spvl, kxtrp &
            )
            do k = 1, nlevo
                dato_flat(base_out + (k - 1) * ninner + inner) = dato_col(k)
            end do
        end do
    end do
end subroutine dvinth2p_nodes_corder_pa_into


subroutine dvinth2p_ecmwf_nodes_corder_pa_into( &
    dati_flat, dato_flat, hbcofa, hbcofb, p0, plevo, intyp, psfc, spvl, kxtrp, &
    nouter, nlevi, ninner, nlevo, varflg, tbot, phis)
    use vinth2p_backend_core, only : real64, convert_levels_to_mb, interpolate_column_ecmwf
    implicit none

    integer, intent(in) :: intyp, kxtrp, nouter, nlevi, ninner, nlevo, varflg
    real(real64), intent(in) :: dati_flat(nouter * nlevi * ninner)
    real(real64), intent(inout) :: dato_flat(nouter * nlevo * ninner)
    real(real64), intent(in) :: hbcofa(nlevi), hbcofb(nlevi)
    real(real64), intent(in) :: p0, plevo(nlevo), psfc(nouter * ninner), spvl
    real(real64), intent(in) :: tbot(nouter * ninner), phis(nouter * ninner)

    real(real64) :: dati_col(nlevi), dato_col(nlevo), plevo_mb(nlevo), p0_mb
    integer :: base_in, base_out, col_idx, inner, k, outer

    p0_mb = p0 * 0.01_real64
    call convert_levels_to_mb(plevo, plevo_mb)

    do outer = 0, nouter - 1
        base_in = outer * nlevi * ninner
        base_out = outer * nlevo * ninner
        do inner = 1, ninner
            col_idx = outer * ninner + inner
            do k = 1, nlevi
                dati_col(k) = dati_flat(base_in + (k - 1) * ninner + inner)
            end do
            call interpolate_column_ecmwf( &
                dati_col, dato_col, hbcofa, hbcofb, p0_mb, plevo_mb, &
                intyp, psfc(col_idx), spvl, kxtrp, varflg, tbot(col_idx), phis(col_idx) &
            )
            do k = 1, nlevo
                dato_flat(base_out + (k - 1) * ninner + inner) = dato_col(k)
            end do
        end do
    end do
end subroutine dvinth2p_ecmwf_nodes_corder_pa_into
