!
! Fixed version - Exact match with pcmin_2013.f algorithm
! Based on Kerry Emanuel's PCMIN subroutine
!

subroutine calculate_pi_single_profile_fixed(sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
                                             num_levels, actual_levels, min_pressure, max_wind, error_flag)
    implicit none

    ! Input parameters
    integer, intent(in) :: num_levels, actual_levels
    real, intent(in) :: sst_in                       ! Sea surface temperature (K)
    real, intent(in) :: psl_in                       ! Sea level pressure (Pa)
    real, intent(in) :: pressure_levels(num_levels)  ! Pressure levels (mb)
    real, intent(in) :: temp_in(num_levels)          ! Temperature (K)
    real, intent(in) :: mixing_ratio_in(num_levels)  ! Mixing ratio (kg/kg)

    ! Output parameters
    real, intent(out) :: min_pressure                ! Minimum central pressure (mb)
    real, intent(out) :: max_wind                    ! Maximum surface wind speed (m/s)
    integer, intent(out) :: error_flag

    ! Constants matching pcmin_2013.f EXACTLY
    real, parameter :: CKCD = 0.9                    ! Ratio of Ck to CD
    real, parameter :: SIG = 0.0                     ! Buoyancy parameter
    real, parameter :: B_EXPONENT = 2.0              ! Eye velocity profile exponent
    real, parameter :: WIND_REDUCTION = 0.8          ! Gradient to 10m wind factor
    real, parameter :: RD = 287.04                   ! Gas constant for dry air
    real, parameter :: UNDEF = -9.99e33

    ! Local variables
    real :: sst_celsius, psl_mb, sst_kelvin, es0
    real :: temp_kelvin(num_levels), mixing_ratio_decimal(num_levels)
    real :: cape_env, cape_rm, cape_sat
    real :: tob_env, tob_rm, tob_sat
    real :: parcel_temp, parcel_mixing_ratio, parcel_pressure
    real :: pressure_estimate, pressure_new
    real :: dissipation_ratio, rs0
    real :: tv1, tvav, cat, catfac
    real :: available_energy
    integer :: surface_level, iteration_count, cape_error_flag
    integer :: k
    logical :: converged

    ! Convert units
    psl_mb = psl_in * 0.01                          ! Pa to mb
    sst_celsius = sst_in - 273.15                   ! K to C
    sst_kelvin = sst_in

    ! Check SST threshold
    if (sst_celsius <= 5.0) then
        min_pressure = UNDEF
        max_wind = UNDEF
        error_flag = 0
        return
    end if

    ! Normalize quantities (matching PCMIN)
    es0 = 6.112 * exp(17.67 * sst_celsius / (243.5 + sst_celsius))

    do k = 1, num_levels
        mixing_ratio_decimal(k) = mixing_ratio_in(k)  ! Already in decimal (kg/kg)
        temp_kelvin(k) = temp_in(k)                   ! Already in K
    end do

    ! Set defaults
    max_wind = 0.0
    min_pressure = psl_mb
    error_flag = 1
    surface_level = 1

    ! Calculate environmental CAPE
    parcel_temp = temp_kelvin(surface_level)
    parcel_mixing_ratio = mixing_ratio_decimal(surface_level)
    parcel_pressure = pressure_levels(surface_level)

    call CAPE_FIXED(parcel_temp, parcel_mixing_ratio, parcel_pressure, &
                    temp_kelvin, mixing_ratio_decimal, pressure_levels, &
                    num_levels, actual_levels, SIG, &
                    cape_env, tob_env, cape_error_flag)

    if (cape_error_flag /= 1) then
        error_flag = 2
        return
    end if

    ! Iterative solution for minimum pressure
    iteration_count = 0
    pressure_estimate = 950.0  ! Initial guess (matching PCMIN)
    converged = .false.

    do while (.not. converged .and. iteration_count < 1000)
        iteration_count = iteration_count + 1

        ! CAPE at radius of maximum winds
        parcel_temp = temp_kelvin(surface_level)
        parcel_pressure = min(pressure_estimate, 1000.0)

        ! EXACT formula from pcmin_2013.f line 103
        parcel_mixing_ratio = 0.622 * mixing_ratio_decimal(surface_level) * psl_mb / &
                             (parcel_pressure * (0.622 + mixing_ratio_decimal(surface_level)) - &
                              mixing_ratio_decimal(surface_level) * psl_mb)

        call CAPE_FIXED(parcel_temp, parcel_mixing_ratio, parcel_pressure, &
                        temp_kelvin, mixing_ratio_decimal, pressure_levels, &
                        num_levels, actual_levels, SIG, &
                        cape_rm, tob_rm, cape_error_flag)

        if (cape_error_flag /= 1) then
            error_flag = 2
            return
        end if

        ! Saturated CAPE at radius of maximum winds
        parcel_temp = sst_kelvin
        parcel_pressure = min(pressure_estimate, 1000.0)

        ! EXACT formula from pcmin_2013.f line 111
        parcel_mixing_ratio = 0.622 * es0 / (parcel_pressure - es0)

        call CAPE_FIXED(parcel_temp, parcel_mixing_ratio, parcel_pressure, &
                        temp_kelvin, mixing_ratio_decimal, pressure_levels, &
                        num_levels, actual_levels, SIG, &
                        cape_sat, tob_sat, cape_error_flag)

        if (cape_error_flag /= 1) then
            error_flag = 2
            return
        end if

        ! Dissipation ratio
        dissipation_ratio = sst_kelvin / tob_sat

        ! Initial estimate of minimum pressure (EXACT from pcmin_2013.f)
        rs0 = parcel_mixing_ratio
        tv1 = temp_kelvin(1) * (1.0 + mixing_ratio_decimal(1)/0.622) / (1.0 + mixing_ratio_decimal(1))
        tvav = 0.5 * (tv1 + sst_kelvin * (1.0 + rs0/0.622) / (1.0 + rs0))

        ! EXACT formula from pcmin_2013.f line 123
        cat = cape_rm - cape_env + 0.5 * CKCD * dissipation_ratio * (cape_sat - cape_rm)
        cat = max(cat, 0.0)
        pressure_new = psl_mb * exp(-cat / (RD * tvav))

        ! Test for convergence (EXACT threshold from pcmin_2013.f)
        if (abs(pressure_new - pressure_estimate) <= 0.2) then
            converged = .true.

            ! Final calculation with corrected factor (EXACT from pcmin_2013.f)
            catfac = 0.5 * (1.0 + 1.0/B_EXPONENT)
            cat = cape_rm - cape_env + CKCD * dissipation_ratio * catfac * (cape_sat - cape_rm)
            cat = max(cat, 0.0)
            min_pressure = psl_mb * exp(-cat / (RD * tvav))
        else
            pressure_estimate = pressure_new

            ! Check for unrealistic values (EXACT from pcmin_2013.f)
            if (iteration_count > 1000 .or. pressure_estimate < 400.0) then
                min_pressure = psl_mb
                error_flag = 0
                return
            end if
        end if
    end do

    ! Check if iteration limit was reached
    if (.not. converged) then
        min_pressure = psl_mb
        error_flag = 0
        return
    end if

    ! Calculate maximum wind speed (EXACT from pcmin_2013.f)
    available_energy = max(0.0, cape_sat - cape_rm)
    max_wind = WIND_REDUCTION * sqrt(CKCD * dissipation_ratio * available_energy)

end subroutine calculate_pi_single_profile_fixed


! CAPE subroutine matching pcmin_2013.f logic exactly
subroutine CAPE_FIXED(parcel_temp, parcel_mixing_ratio, parcel_pressure, &
                      temp_profile, mixing_ratio_profile, pressure_profile, &
                      array_size, num_points, buoyancy_param, &
                      cape_value, outflow_temp, error_flag)
    implicit none

    ! Input parameters
    real, intent(in) :: parcel_temp                  ! Initial parcel temperature (K)
    real, intent(in) :: parcel_mixing_ratio          ! Initial parcel mixing ratio (decimal)
    real, intent(in) :: parcel_pressure              ! Initial parcel pressure (mb)
    real, intent(in) :: temp_profile(array_size)     ! Environmental temperature profile (K)
    real, intent(in) :: mixing_ratio_profile(array_size) ! Environmental mixing ratio profile (decimal)
    real, intent(in) :: pressure_profile(array_size) ! Pressure profile (mb)
    integer, intent(in) :: array_size, num_points
    real, intent(in) :: buoyancy_param               ! Buoyancy parameter (0-1)

    ! Output parameters
    real, intent(out) :: cape_value                  ! Calculated CAPE (J/kg)
    real, intent(out) :: outflow_temp                ! Temperature at level of neutral buoyancy (K)
    integer, intent(out) :: error_flag

    ! Constants (EXACT from pcmin_2013.f)
    real, parameter :: CPD = 1005.7
    real, parameter :: CPV = 1870.0
    real, parameter :: CL = 2500.0    ! Note: pcmin_2013.f uses 2500, not 4190
    real, parameter :: CPVMCL = CPV - CL
    real, parameter :: RV = 461.5
    real, parameter :: RD = 287.04
    real, parameter :: EPS = RD/RV
    real, parameter :: ALV0 = 2.501e6

    ! Local variables
    real :: virtual_temp_diff(100)
    real :: parcel_temp_celsius, saturation_pressure, vapor_pressure
    real :: relative_humidity, latent_heat, reversible_entropy
    real :: chi, lifted_condensation_pressure
    real :: parcel_temp_lifted, parcel_mixing_ratio_lifted
    real :: temp_celsius, saturation_pressure_env
    real :: entropy_rate_temp, vapor_pressure_parcel
    real :: entropy_parcel, adaptation_param, temp_new
    real :: mean_mixing_ratio, virtual_temp_parcel
    real :: positive_area, negative_area
    real :: pressure_factor, area_above_lnb, pressure_lnb, temp_lnb
    integer :: level, min_level, level_neutral_buoyancy
    integer :: iteration_count, ncmax
    logical :: converged

    ! Initialize
    cape_value = 0.0
    outflow_temp = temp_profile(1)
    error_flag = 1

    ! Validate input
    if (parcel_mixing_ratio < 1.0e-6 .or. parcel_temp < 200.0) then
        error_flag = 0
        return
    end if

    ! Calculate parcel properties
    parcel_temp_celsius = parcel_temp - 273.15
    saturation_pressure = 6.112 * exp(17.67 * parcel_temp_celsius / (243.5 + parcel_temp_celsius))
    vapor_pressure = parcel_mixing_ratio * parcel_pressure / (EPS + parcel_mixing_ratio)
    relative_humidity = min(vapor_pressure / saturation_pressure, 1.0)
    latent_heat = ALV0 + CPVMCL * parcel_temp_celsius

    ! Reversible entropy
    reversible_entropy = (CPD + parcel_mixing_ratio * CL) * log(parcel_temp) - &
                        RD * log(parcel_pressure - vapor_pressure) + &
                        latent_heat * parcel_mixing_ratio / parcel_temp - &
                        parcel_mixing_ratio * RV * log(relative_humidity)

    ! Lifted condensation level
    chi = parcel_temp / (1669.0 - 122.0 * relative_humidity - parcel_temp)
    lifted_condensation_pressure = parcel_pressure * (relative_humidity ** chi)

    ! Initialize
    ncmax = 0
    do level = 1, num_points
        virtual_temp_diff(level) = 0.0
    end do

    ! Find minimum level
    min_level = 1000000
    do level = 1, num_points
        if (pressure_profile(level) < 59.0 .or. &
            pressure_profile(level) >= parcel_pressure) cycle
        min_level = min(min_level, level)
    end do

    if (min_level == 1000000) return

    ! Calculate parcel ascent
    do level = min_level, num_points
        if (pressure_profile(level) < 59.0 .or. &
            pressure_profile(level) >= parcel_pressure) cycle

        if (pressure_profile(level) >= lifted_condensation_pressure) then
            ! Below LCL - dry adiabatic
            parcel_temp_lifted = parcel_temp * (pressure_profile(level) / parcel_pressure) ** (RD/CPD)
            parcel_mixing_ratio_lifted = parcel_mixing_ratio

            virtual_temp_parcel = parcel_temp_lifted * (1.0 + parcel_mixing_ratio_lifted/EPS) / &
                                 (1.0 + parcel_mixing_ratio_lifted)
            virtual_temp_diff(level) = virtual_temp_parcel - &
                                      temp_profile(level) * (1.0 + mixing_ratio_profile(level)/EPS) / &
                                      (1.0 + mixing_ratio_profile(level))
        else
            ! Above LCL - moist adiabatic
            parcel_temp_lifted = temp_profile(level)
            temp_celsius = temp_profile(level) - 273.15
            saturation_pressure_env = 6.112 * exp(17.67 * temp_celsius / (243.5 + temp_celsius))
            parcel_mixing_ratio_lifted = EPS * saturation_pressure_env / (pressure_profile(level) - saturation_pressure_env)

            ! Iterative calculation
            iteration_count = 0
            converged = .false.

            do while (.not. converged .and. iteration_count < 500)
                iteration_count = iteration_count + 1

                latent_heat = ALV0 + CPVMCL * (parcel_temp_lifted - 273.15)
                entropy_rate_temp = (CPD + parcel_mixing_ratio * CL + &
                                    latent_heat * latent_heat * parcel_mixing_ratio_lifted / (RV * parcel_temp_lifted * parcel_temp_lifted)) / &
                                    parcel_temp_lifted
                vapor_pressure_parcel = parcel_mixing_ratio_lifted * pressure_profile(level) / (EPS + parcel_mixing_ratio_lifted)
                entropy_parcel = (CPD + parcel_mixing_ratio * CL) * log(parcel_temp_lifted) - &
                                RD * log(pressure_profile(level) - vapor_pressure_parcel) + &
                                latent_heat * parcel_mixing_ratio_lifted / parcel_temp_lifted

                if (iteration_count < 3) then
                    adaptation_param = 0.3
                else
                    adaptation_param = 1.0
                end if

                temp_new = parcel_temp_lifted + adaptation_param * (reversible_entropy - entropy_parcel) / entropy_rate_temp

                if (abs(temp_new - parcel_temp_lifted) <= 0.001) then
                    converged = .true.
                else
                    parcel_temp_lifted = temp_new
                    temp_celsius = parcel_temp_lifted - 273.15
                    saturation_pressure_env = 6.112 * exp(17.67 * temp_celsius / (243.5 + temp_celsius))

                    if (iteration_count > 500 .or. saturation_pressure_env > (pressure_profile(level) - 1.0)) then
                        error_flag = 2
                        return
                    end if

                    parcel_mixing_ratio_lifted = EPS * saturation_pressure_env / (pressure_profile(level) - saturation_pressure_env)
                end if
            end do

            if (.not. converged) then
                error_flag = 2
                return
            end if

            ncmax = max(iteration_count, ncmax)

            ! Buoyancy
            mean_mixing_ratio = buoyancy_param * parcel_mixing_ratio_lifted + (1.0 - buoyancy_param) * parcel_mixing_ratio
            virtual_temp_parcel = parcel_temp_lifted * (1.0 + parcel_mixing_ratio_lifted/EPS) / (1.0 + mean_mixing_ratio)
            virtual_temp_diff(level) = virtual_temp_parcel - &
                                      temp_profile(level) * (1.0 + mixing_ratio_profile(level)/EPS) / &
                                      (1.0 + mixing_ratio_profile(level))
        end if
    end do

    ! Find level of neutral buoyancy
    positive_area = 0.0
    negative_area = 0.0

    level_neutral_buoyancy = 1
    do level = num_points, min_level, -1
        if (virtual_temp_diff(level) > 0.0) then
            level_neutral_buoyancy = max(level_neutral_buoyancy, level)
        end if
    end do

    if (level_neutral_buoyancy == 1) return

    ! Calculate CAPE
    if (level_neutral_buoyancy > 1) then
        do level = min_level + 1, level_neutral_buoyancy
            pressure_factor = RD * (virtual_temp_diff(level) + virtual_temp_diff(level-1)) * &
                             (pressure_profile(level-1) - pressure_profile(level)) / &
                             (pressure_profile(level) + pressure_profile(level-1))
            positive_area = positive_area + max(pressure_factor, 0.0)
            negative_area = negative_area - min(pressure_factor, 0.0)
        end do

        ! Area between parcel level and first level
        pressure_factor = RD * (parcel_pressure - pressure_profile(min_level)) / &
                         (parcel_pressure + pressure_profile(min_level))
        positive_area = positive_area + pressure_factor * max(virtual_temp_diff(min_level), 0.0)
        negative_area = negative_area - pressure_factor * min(virtual_temp_diff(min_level), 0.0)

        ! Residual area above LNB
        area_above_lnb = 0.0
        outflow_temp = temp_profile(level_neutral_buoyancy)

        if (level_neutral_buoyancy < num_points) then
            pressure_lnb = (pressure_profile(level_neutral_buoyancy+1) * virtual_temp_diff(level_neutral_buoyancy) - &
                           pressure_profile(level_neutral_buoyancy) * virtual_temp_diff(level_neutral_buoyancy+1)) / &
                           (virtual_temp_diff(level_neutral_buoyancy) - virtual_temp_diff(level_neutral_buoyancy+1))
            area_above_lnb = RD * virtual_temp_diff(level_neutral_buoyancy) * &
                            (pressure_profile(level_neutral_buoyancy) - pressure_lnb) / &
                            (pressure_profile(level_neutral_buoyancy) + pressure_lnb)
            temp_lnb = (temp_profile(level_neutral_buoyancy) * (pressure_lnb - pressure_profile(level_neutral_buoyancy+1)) + &
                       temp_profile(level_neutral_buoyancy+1) * (pressure_profile(level_neutral_buoyancy) - pressure_lnb)) / &
                       (pressure_profile(level_neutral_buoyancy) - pressure_profile(level_neutral_buoyancy+1))
            outflow_temp = temp_lnb
        end if

        cape_value = positive_area + area_above_lnb - negative_area
        cape_value = max(cape_value, 0.0)
    end if

end subroutine CAPE_FIXED
