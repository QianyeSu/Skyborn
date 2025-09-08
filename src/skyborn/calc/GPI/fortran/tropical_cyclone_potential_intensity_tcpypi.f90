!
! Modern Fortran implementation of tcpyPI (Tropical Cyclone Potential Intensity)
! Translated from Kerry Emanuel's Python tcpyPI implementation
! Compatible with tropical_cyclone_potential_intensity.f90 unit conventions
!
! INPUT UNITS (matching tropical_cyclone_potential_intensity.f90):
!   - Sea surface temperature: [K]
!   - Pressure levels: [hPa]
!   - Temperature profiles: [K]
!   - Mixing ratio: [kg/kg]
!
! Author: Translated by Qianye Su
! Original: Kerry Emanuel et al.
!

!--------------------------------------------------------------------------
! Physical and Thermodynamic Constants (from tcpyPI constants.py)
!--------------------------------------------------------------------------

subroutine get_physical_constants(cpd, cpv, cl, rv, rd, eps, alv0)
    ! Return physical constants used in tcpyPI calculations
    implicit none

    real, intent(out) :: cpd    ! Specific heat of dry air [J/kg·K]
    real, intent(out) :: cpv    ! Specific heat of water vapor [J/kg·K]
    real, intent(out) :: cl     ! Specific heat of liquid water [J/kg·K]
    real, intent(out) :: rv     ! Gas constant of water vapor [J/kg·K]
    real, intent(out) :: rd     ! Gas constant of dry air [J/kg·K]
    real, intent(out) :: eps    ! Ratio of gas constants [unitless]
    real, intent(out) :: alv0   ! Latent heat of vaporization at 0°C [J/kg]

    cpd = 1005.7
    cpv = 1870.0
    cl = 2500.0      ! Modified from original 4190.0
    rv = 461.5
    rd = 287.04
    eps = rd / rv
    alv0 = 2.501e6

end subroutine get_physical_constants

!--------------------------------------------------------------------------
! Temperature Conversion Utilities
!--------------------------------------------------------------------------

subroutine temperature_kelvin_to_celsius(temp_kelvin, temp_celsius)
    !---------------------------------------------------------------------------
    ! PURPOSE: Convert temperature from Kelvin to Celsius
    ! FORMULA: T(°C) = T(K) - 273.15
    !---------------------------------------------------------------------------
    implicit none

    ! INPUT:
    real, intent(in) :: temp_kelvin    ! Temperature [K]
    ! OUTPUT:
    real, intent(out) :: temp_celsius   ! Temperature [°C]

    temp_celsius = temp_kelvin - 273.15
end subroutine temperature_kelvin_to_celsius

subroutine temperature_celsius_to_kelvin(temp_celsius, temp_kelvin)
    ! Convert temperature from Celsius to Kelvin
    implicit none

    real, intent(in) :: temp_celsius
    real, intent(out) :: temp_kelvin

    temp_kelvin = temp_celsius + 273.15
end subroutine temperature_celsius_to_kelvin

!--------------------------------------------------------------------------
! Saturation Vapor Pressure Calculations
!--------------------------------------------------------------------------

subroutine saturation_vapor_pressure_magnus(temp_celsius, es)
    !---------------------------------------------------------------------------
    ! PURPOSE: Calculate saturated water vapor pressure using Magnus formula
    ! REFERENCE: August-Roche-Magnus formula (Clausius-Clapeyron approximation)
    ! FORMULA: es = 6.112 * exp(17.67*T/(243.5+T)) where T is in °C
    !---------------------------------------------------------------------------
    implicit none

    ! INPUT:
    real, intent(in) :: temp_celsius    ! Temperature [°C]
    ! OUTPUT:
    real, intent(out) :: es             ! Saturation vapor pressure [hPa]

    ! Magnus formula constants
    real, parameter :: es_base = 6.112
    real, parameter :: magnus_a = 17.67
    real, parameter :: magnus_b = 243.5

    es = es_base * exp(magnus_a * temp_celsius / (magnus_b + temp_celsius))

end subroutine saturation_vapor_pressure_magnus

!--------------------------------------------------------------------------
! Moisture and Vapor Pressure Calculations
!--------------------------------------------------------------------------

subroutine vapor_pressure_from_mixing_ratio(mixing_ratio, pressure, vapor_pressure)
    ! Calculate vapor pressure from mixing ratio and total pressure
    implicit none

    real, intent(in) :: mixing_ratio     ! Mixing ratio [kg/kg]
    real, intent(in) :: pressure         ! Total pressure [hPa]
    real, intent(out) :: vapor_pressure  ! Vapor pressure [hPa]

    real :: eps

    eps = 287.04 / 461.5  ! RD/RV
    vapor_pressure = mixing_ratio * pressure / (eps + mixing_ratio)

end subroutine vapor_pressure_from_mixing_ratio

subroutine mixing_ratio_from_vapor_pressure(vapor_pressure, pressure, mixing_ratio)
    ! Calculate mixing ratio from vapor pressure and total pressure
    implicit none

    real, intent(in) :: vapor_pressure   ! Vapor pressure [hPa]
    real, intent(in) :: pressure         ! Total pressure [hPa]
    real, intent(out) :: mixing_ratio    ! Mixing ratio [kg/kg]

    real :: eps

    eps = 287.04 / 461.5  ! RD/RV
    mixing_ratio = eps * vapor_pressure / (pressure - vapor_pressure)

end subroutine mixing_ratio_from_vapor_pressure

!--------------------------------------------------------------------------
! Latent Heat and Thermodynamic Properties
!--------------------------------------------------------------------------

subroutine latent_heat_vaporization(temp_celsius, latent_heat)
    ! Calculate latent heat of vaporization as a function of temperature
    implicit none

    real, intent(in) :: temp_celsius     ! Temperature [°C]
    real, intent(out) :: latent_heat     ! Latent heat [J/kg]

    real, parameter :: alv0 = 2.501e6  ! Latent heat at 0°C
    real, parameter :: cpv = 1870.0    ! Specific heat of vapor
    real, parameter :: cl = 2500.0     ! Modified specific heat of liquid

    latent_heat = alv0 + (cpv - cl) * temp_celsius

end subroutine latent_heat_vaporization

!--------------------------------------------------------------------------
! Density Temperature Calculations
!--------------------------------------------------------------------------

subroutine density_temperature(temperature, total_mixing_ratio, vapor_mixing_ratio, &
                              density_temp)
    ! Calculate density (virtual) temperature in K
    implicit none

    real, intent(in) :: temperature            ! Temperature [K]
    real, intent(in) :: total_mixing_ratio     ! Total water mixing ratio [kg/kg]
    real, intent(in) :: vapor_mixing_ratio     ! Water vapor mixing ratio [kg/kg]
    real, intent(out) :: density_temp          ! Density temperature [K]

    real :: eps

    eps = 287.04 / 461.5  ! RD/RV
    density_temp = temperature * (1.0 + vapor_mixing_ratio/eps) / &
                   (1.0 + total_mixing_ratio)

end subroutine density_temperature

!--------------------------------------------------------------------------
! Entropy Calculations
!--------------------------------------------------------------------------

subroutine specific_entropy(temperature, mixing_ratio, pressure, entropy)
    ! Calculate total specific entropy per unit mass of dry air (E94, Eqn. 4.5.9)
    implicit none

    real, intent(in) :: temperature      ! Temperature [K]
    real, intent(in) :: mixing_ratio     ! Mixing ratio [kg/kg]
    real, intent(in) :: pressure         ! Pressure [hPa]
    real, intent(out) :: entropy         ! Specific entropy [J/kg·K]

    real :: temp_celsius, vapor_pressure, sat_vapor_pressure
    real :: relative_humidity, latent_heat
    real, parameter :: cpd = 1005.7
    real, parameter :: cl = 2500.0
    real, parameter :: rd = 287.04
    real, parameter :: rv = 461.5
    call temperature_kelvin_to_celsius(temperature, temp_celsius)
    call vapor_pressure_from_mixing_ratio(mixing_ratio, pressure, vapor_pressure)
    call saturation_vapor_pressure_magnus(temp_celsius, sat_vapor_pressure)

    relative_humidity = min(vapor_pressure / sat_vapor_pressure, 1.0)
    call latent_heat_vaporization(temp_celsius, latent_heat)

    entropy = (cpd + mixing_ratio * cl) * log(temperature) - &
              rd * log(pressure - vapor_pressure) + &
              latent_heat * mixing_ratio / temperature - &
              mixing_ratio * rv * log(relative_humidity)

end subroutine specific_entropy

!--------------------------------------------------------------------------
! Lifting Condensation Level (LCL) Calculation
!--------------------------------------------------------------------------

subroutine lifting_condensation_level_pressure(parcel_temp, relative_humidity, &
                                              parcel_pressure, lcl_pressure)
    ! Calculate empirical lifting condensation level pressure
    implicit none

    real, intent(in) :: parcel_temp        ! Parcel temperature [K]
    real, intent(in) :: relative_humidity  ! Relative humidity [0-1]
    real, intent(in) :: parcel_pressure    ! Parcel pressure [hPa]
    real, intent(out) :: lcl_pressure      ! LCL pressure [hPa]

    real :: chi
    real, parameter :: lcl_a = 1669.0
    real, parameter :: lcl_b = 122.0

    chi = parcel_temp / (lcl_a - lcl_b * relative_humidity - parcel_temp)
    lcl_pressure = parcel_pressure * (relative_humidity ** chi)

end subroutine lifting_condensation_level_pressure

!--------------------------------------------------------------------------
! CORE tcpyPI Thermodynamic Functions (Exactly Matching Python Implementation)
!--------------------------------------------------------------------------

subroutine entropy_s(temperature, mixing_ratio, pressure, entropy)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate total specific entropy per unit mass of dry air
    !   Following Emanuel (1994) Equation 4.5.9
    !
    ! FORMULA:
    !   S = (Cpd + R*Cl)*ln(T) - Rd*ln(P-e) + Lv*R/T - R*Rv*ln(RH)
    !   where e is vapor pressure and RH is relative humidity
    !---------------------------------------------------------------------------
    implicit none

    ! INPUT:
    real, intent(in) :: temperature
        ! Temperature [K]
    real, intent(in) :: mixing_ratio
        ! Water vapor mixing ratio [kg/kg]
    real, intent(in) :: pressure
        ! Pressure [hPa]

    ! OUTPUT:
    real, intent(out) :: entropy
        ! Specific entropy [J/kg·K]

    real :: ev_vapor_pressure, es_sat_pressure, relative_humidity, alv_latent_heat
    real :: temp_celsius

    ! Physical constants exactly from tcpyPI constants.py
    real, parameter :: cpd = 1005.7       ! CPD
    real, parameter :: cl = 2500.0        ! CL (modified value)
    real, parameter :: rd = 287.04        ! RD
    real, parameter :: rv = 461.5         ! RV

    ! Line 189: EV=ev(R,P)
    call vapor_pressure_from_mixing_ratio(mixing_ratio, pressure, ev_vapor_pressure)

    ! Line 190: ES=es_cc(T-273.15)
    call temperature_kelvin_to_celsius(temperature, temp_celsius)
    call saturation_vapor_pressure_magnus(temp_celsius, es_sat_pressure)

    ! Line 191: RH=min([EV/ES,1.0])
    relative_humidity = min(ev_vapor_pressure / es_sat_pressure, 1.0)

    ! Line 192: ALV=Lv(T-273.15)
    call latent_heat_vaporization(temp_celsius, alv_latent_heat)

    ! Line 193: S=(CPD+R*CL)*log(T)-RD*log(P-EV)+ALV*R/T-R*RV*log(RH)
    entropy = (cpd + mixing_ratio * cl) * log(temperature) - &
              rd * log(pressure - ev_vapor_pressure) + &
              alv_latent_heat * mixing_ratio / temperature - &
              mixing_ratio * rv * log(relative_humidity)

end subroutine entropy_s

subroutine trho(temperature, total_mixing_ratio, vapor_mixing_ratio, density_temp)
    ! Calculate density temperature in K - exactly matching tcpyPI utilities.Trho()
    ! Line 214: return T*(1.+R/constants.EPS)/(1.+RT)
    implicit none

    real, intent(in) :: temperature           ! Temperature [K]
    real, intent(in) :: total_mixing_ratio    ! Total water mixing ratio [kg/kg]
    real, intent(in) :: vapor_mixing_ratio    ! Water vapor mixing ratio [kg/kg]
    real, intent(out) :: density_temp         ! Density temperature [K]

    real, parameter :: eps = 0.622            ! RD/RV from constants.py

    density_temp = temperature * (1.0 + vapor_mixing_ratio / eps) / (1.0 + total_mixing_ratio)

end subroutine trho

subroutine e_plcl(parcel_temp, relative_humidity, parcel_pressure, lcl_pressure)
    ! Calculate empirical lifting condensation level pressure
    ! Exactly matching tcpyPI utilities.e_pLCL() - line 234
    implicit none

    real, intent(in) :: parcel_temp        ! Parcel temperature [K]
    real, intent(in) :: relative_humidity  ! Relative humidity [0-1]
    real, intent(in) :: parcel_pressure    ! Parcel pressure [hPa]
    real, intent(out) :: lcl_pressure      ! LCL pressure [hPa]

    ! Constants from tcpyPI constants.py (lines 22-23)
    real, parameter :: lcl_a = 1669.0
    real, parameter :: lcl_b = 122.0

    ! Line 234: return PP*(RH**(TP/(constants.A-constants.B*RH-TP)))
    lcl_pressure = parcel_pressure * (relative_humidity ** &
                   (parcel_temp / (lcl_a - lcl_b * relative_humidity - parcel_temp)))

end subroutine e_plcl

subroutine solve_temperature_from_entropy(target_entropy, pressure, dry_mixing_ratio, &
                                                 initial_temp, final_temp, final_mixing_ratio, &
                                                 converged_flag)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Solve for temperature that gives a specified saturated entropy value
    !   using Newton-Raphson iteration
    !
    ! METHOD:
    !   Iteratively solves SG(T) = S where SG is saturated entropy
    !   Uses Newton-Raphson method with adaptive step size
    !
    ! CONVERGENCE:
    !   Tolerance: 0.001 K
    !   Maximum iterations: 500
    !---------------------------------------------------------------------------
    implicit none

    ! INPUT:
    real, intent(in) :: target_entropy
        ! Target saturated entropy [J/kg/K] to solve for
    real, intent(in) :: pressure
        ! Ambient pressure [hPa]
    real, intent(in) :: dry_mixing_ratio
        ! Parcel mixing ratio for dry component [kg/kg]
    real, intent(in) :: initial_temp
        ! Initial guess for temperature [K]

    ! OUTPUT:
    real, intent(out) :: final_temp
        ! Temperature [K] that satisfies SG(T) = target_entropy
    real, intent(out) :: final_mixing_ratio
        ! Computed saturated mixing ratio [kg/kg] at final_temp
    integer, intent(out) :: converged_flag
        ! Convergence status: 1=success, 2=did not converge

    ! Local variables matching tcpyPI variable names
    real :: tgnew, tg, tjc, es, rg, enew, alv, sl, em, sg, ap
    integer :: nc

    ! Physical constants from tcpyPI constants.py
    real, parameter :: cpd = 1005.7
    real, parameter :: cl = 2500.0
    real, parameter :: rd = 287.04
    real, parameter :: rv = 461.5
    real, parameter :: convergence_tol = 0.001
    integer, parameter :: max_iterations = 500

    ! Initial default values before loop (lines 333-337)
    tgnew = initial_temp
    call temperature_kelvin_to_celsius(initial_temp, tjc)
    call saturation_vapor_pressure_magnus(tjc, es)
    call mixing_ratio_from_vapor_pressure(es, pressure, rg)

    ! Set loop counter and initial condition (lines 345-346)
    nc = 0
    tg = 0.0

    ! Loop until converges or bails out (line 349)
    do while (abs(tgnew - tg) > convergence_tol)

        ! Parcel temperature and mixing ratio during this iteration (lines 351-355)
        tg = tgnew
        call temperature_kelvin_to_celsius(tg, tjc)
        call saturation_vapor_pressure_magnus(tjc, enew)
        call mixing_ratio_from_vapor_pressure(enew, pressure, rg)

        ! Increase iteration count (line 358)
        nc = nc + 1

        ! Calculate estimates of rates of change of entropy with temperature (lines 365-371)
        call latent_heat_vaporization(tjc, alv)

        ! Calculate rate of change of entropy with temperature, s_ell (line 367)
        sl = (cpd + dry_mixing_ratio * cl + alv * alv * rg / (rv * tg * tg)) / tg

        ! Calculate vapor pressure from mixing ratio (line 368)
        call vapor_pressure_from_mixing_ratio(rg, pressure, em)

        ! Calculate saturated entropy, s_k (line 371)
        ! Note: last term vanishes with saturation (RH=1)
        sg = (cpd + dry_mixing_ratio * cl) * log(tg) - &
             rd * log(pressure - em) + alv * rg / tg

        ! Convergence speed varies with iteration count (lines 374-379)
        if (nc < 3) then
            ap = 0.3  ! converge slowly with smaller step
        else
            ap = 1.0  ! speed up when nearing convergence
        end if

        ! Find new temperature (line 381)
        tgnew = tg + ap * (target_entropy - sg) / sl

        ! Check for convergence failure (lines 386-388)
        if (nc > max_iterations .or. enew > (pressure - 1.0)) then
            converged_flag = 2  ! Did not converge
            final_temp = tg
            final_mixing_ratio = rg
            return
        end if

    end do

    ! Success (lines 393-394)
    converged_flag = 1
    final_temp = tg
    final_mixing_ratio = rg

end subroutine solve_temperature_from_entropy

!--------------------------------------------------------------------------
! Main tcpyPI Potential Intensity Calculation
!--------------------------------------------------------------------------

subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                               pressure_levels_hpa, temp_celsius_profile, &
                                               mixing_ratio_gkg_profile, num_levels, &
                                               max_wind_speed, min_pressure_hpa, &
                                               outflow_temp, outflow_level, &
                                               status_flag, ckcd, ascent_flag, diss_flag, &
                                               v_reduc, miss_handle)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate maximum potential intensity of tropical cyclones given
    !   environmental conditions (Bister & Emanuel 2002 algorithm)
    !
    ! REFERENCE:
    !   - Emanuel, K. (1994): Atmospheric Convection, Oxford University Press
    !   - Bister, M. and K. Emanuel (1998): Dissipative heating and hurricane intensity
    !   - Bister, M. and K. Emanuel (2002): Low frequency variability of tropical
    !     cyclone potential intensity
    !
    ! ALGORITHM:
    !   Iteratively finds the pressure at radius of maximum winds that maximizes
    !   the difference between CAPE at the eyewall and environmental CAPE
    !---------------------------------------------------------------------------
    implicit none

    !---------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    !---------------------------------------------------------------------------
    integer, intent(in) :: num_levels
        ! Number of vertical levels in atmospheric profile

    real, intent(in) :: sst_celsius
        ! Sea surface temperature [°C]
        ! Valid range: 5.0 to 100.0°C (checked internally)

    real, intent(in) :: msl_pressure_hpa
        ! Mean sea level pressure [hPa]
        ! Typical range: 900-1050 hPa

    real, intent(in) :: pressure_levels_hpa(num_levels)
        ! Pressure levels [hPa], ordered from surface to top
        ! Must be monotonically decreasing with index

    real, intent(in) :: temp_celsius_profile(num_levels)
        ! Temperature profile [°C] at each pressure level
        ! Should extend to at least tropopause

    real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
        ! Water vapor mixing ratio [g/kg] at each pressure level
        ! Missing values can be set to 0.0 above boundary layer

    real, intent(in), optional :: ckcd
        ! Ratio of exchange coefficients Ck/Cd [dimensionless]
        ! Default: 0.9 (based on Wing et al. 2015)
        ! Typical range: 0.5-2.0

    integer, intent(in), optional :: ascent_flag
        ! Type of parcel ascent:
        !   0 = Reversible ascent (default)
        !   1 = Pseudo-adiabatic ascent

    integer, intent(in), optional :: diss_flag
        ! Include dissipative heating:
        !   0 = No dissipative heating
        !   1 = Include dissipative heating (default)

    real, intent(in), optional :: v_reduc
        ! Reduction factor from gradient to 10m winds [dimensionless]
        ! Default: 0.8 (Emanuel 2000, Powell 1980)
        ! Typical range: 0.7-0.9

    integer, intent(in), optional :: miss_handle
        ! How to handle missing values:
        !   0 = Ignore NaN, compute PI anyway
        !   1 = Return missing if NaN found (default)

    !---------------------------------------------------------------------------
    ! OUTPUT PARAMETERS:
    !---------------------------------------------------------------------------
    real, intent(out) :: max_wind_speed
        ! Maximum potential wind speed [m/s]
        ! At 10m height, reduced from gradient wind

    real, intent(out) :: min_pressure_hpa
        ! Minimum central pressure [hPa]
        ! At center of tropical cyclone

    real, intent(out) :: outflow_temp
        ! Temperature at outflow level [K]
        ! Temperature at level of neutral buoyancy

    real, intent(out) :: outflow_level
        ! Pressure of outflow level [hPa]
        ! Level of neutral buoyancy for eyewall parcel

    integer, intent(out) :: status_flag
        ! Computation status:
        !   1 = Success
        !   0 = No convergence or invalid input
        !   2 = CAPE routine failed to converge
        !   3 = Missing data in profile

    ! Local variables - exactly following tcpyPI variable naming and logic
    real :: sst_kelvin                          ! SSTK in tcpyPI
    real, dimension(:), allocatable :: temp_kelvin_profile     ! T in tcpyPI
    real, dimension(:), allocatable :: mixing_ratio_kgkg_profile ! R in tcpyPI (converted from g/kg to kg/kg)

    ! tcpyPI algorithm variables (matching Python names)
    real :: capea, capem, capems                 ! Environmental, max wind, saturated CAPE
    real :: toms, lnbs                          ! Outflow temp and level from CAPE calls
    real :: rat                                 ! Thermodynamic efficiency ratio
    real :: pm, pmold, pnew                     ! Pressure iteration variables
    real :: rs0                                 ! Surface saturation mixing ratio
    real :: tv0, tvsst, tvav                    ! Density temperatures
    real :: cat, catfac                         ! CAT factor variables
    real :: fac                                 ! Wind calculation factor
    real :: es0                                 ! Saturated vapor pressure at SST
    integer :: np                               ! Iteration counter
    integer :: nk                               ! Surface level index
    integer :: ifl_cape                         ! CAPE function flag

    ! Local parameter values (with defaults matching tcpyPI)
    real :: ckcd_local
    integer :: ascent_flag_local, diss_flag_local, miss_handle_local
    real :: v_reduc_local
    real, parameter :: ptop_local = 50.0  ! Fixed top pressure [hPa]

    ! Physical constants (exactly from tcpyPI constants.py)
    real, parameter :: missing_value = -9999.0
    real, parameter :: rd_constant = 287.04      ! Gas constant for dry air
    real, parameter :: b_constant = 2.0          ! Exponent b from constants.py
    real, parameter :: pressure_tolerance = 0.5  ! Convergence tolerance [hPa]
    real, parameter :: max_iterations = 200      ! Maximum iterations
    real, parameter :: min_pressure_limit = 400.0 ! Minimum allowed pressure [hPa]
    real, parameter :: eps_const = 0.622         ! Rd/Rv constant for optimization
    real, parameter :: inv_eps = 1.607717        ! 1.0/0.622 pre-computed for speed
    real, parameter :: inv_rd = 0.0034835        ! 1.0/287.04 pre-computed for speed
    real, parameter :: inv_b = 0.5               ! 1.0/b_constant pre-computed

    ! Set default parameter values (matching tcpyPI defaults)
    ckcd_local = 0.9
    if (present(ckcd)) ckcd_local = ckcd

    ascent_flag_local = 0
    if (present(ascent_flag)) ascent_flag_local = ascent_flag

    diss_flag_local = 1
    if (present(diss_flag)) diss_flag_local = diss_flag

    v_reduc_local = 0.8
    if (present(v_reduc)) v_reduc_local = v_reduc

    ! ptop_local is now a fixed parameter = 50.0 hPa

    miss_handle_local = 1
    if (present(miss_handle)) miss_handle_local = miss_handle

    ! Initialize outputs
    max_wind_speed = missing_value
    min_pressure_hpa = missing_value
    outflow_temp = missing_value
    outflow_level = missing_value
    status_flag = 0

    ! Allocate local arrays
    allocate(temp_kelvin_profile(num_levels))
    allocate(mixing_ratio_kgkg_profile(num_levels))

    ! Unit conversions exactly following tcpyPI (lines 532-534)
    call temperature_celsius_to_kelvin(sst_celsius, sst_kelvin)  ! SSTK=utilities.T_Ctok(SSTC)
    do nk = 1, num_levels
        call temperature_celsius_to_kelvin(temp_celsius_profile(nk), temp_kelvin_profile(nk)) ! T=utilities.T_Ctok(TC)
        mixing_ratio_kgkg_profile(nk) = mixing_ratio_gkg_profile(nk) * 0.001  ! R=R*0.001
    end do

    ! Validation checks exactly following tcpyPI (lines 536-564)
    ! CHECK 1a: do SSTs exceed 5C?
    ! CHECK 1b: are SSTs less than 100C (if not, indicative of input in kelvin)
    if (sst_celsius <= 5.0 .or. sst_celsius > 100.0) then
        status_flag = 0  ! IFL=0
        if (allocated(temp_kelvin_profile)) deallocate(temp_kelvin_profile)
        if (allocated(mixing_ratio_kgkg_profile)) deallocate(mixing_ratio_kgkg_profile)
        return
    end if

    ! CHECK 2a: do Temperature profiles exceed 100K?
    ! CHECK 2b: are Temperatures in Celsius less than 100C (if not, indicative of input in kelvin)
    if (minval(temp_kelvin_profile) <= 100.0 .or. maxval(temp_celsius_profile) > 100.0) then
        status_flag = 0  ! IFL=0
        return
    end if

    ! Set Missing mixing ratios to zero g/g, following Kerry's BE02 algorithm (line 567)
    where (mixing_ratio_kgkg_profile /= missing_value)
        mixing_ratio_kgkg_profile = max(mixing_ratio_kgkg_profile, 0.0)
    elsewhere
        mixing_ratio_kgkg_profile = 0.0
    end where

    ! Saturated water vapor pressure from Clausius-Clapeyron relation/August-Roche-Magnus formula (line 571)
    call saturation_vapor_pressure_magnus(sst_celsius, es0)

    ! Define the level from which parcels lifted (first pressure level) (line 574)
    nk = 1

    ! *** Find environmental CAPE *** (lines 577-584)
    call calculate_cape(temp_kelvin_profile(nk), mixing_ratio_kgkg_profile(nk), pressure_levels_hpa(nk), &
                              temp_kelvin_profile, mixing_ratio_kgkg_profile, pressure_levels_hpa, &
                              num_levels, ascent_flag_local, miss_handle_local, &
                              capea, toms, lnbs, ifl_cape)

    ! if the CAPE function tripped a flag, set the output IFL to it (lines 585-587)
    if (ifl_cape /= 1) then
        status_flag = ifl_cape
        return
    end if

    ! *** Begin iteration to find minimum pressure *** (lines 590-598)
    np = 0              ! loop counter
    pm = 970.0          ! initial pressure estimate
    pmold = pm         ! initial condition from minimum pressure
    pnew = 0.0         ! initial condition from minimum pressure
    status_flag = 1    ! Default flag for CAPE calculation

    ! loop until convergence or bail out (line 601)
    do while (abs(pnew - pmold) > pressure_tolerance)

        ! *** Find CAPE at radius of maximum winds *** (lines 604-616)
        pm = min(pm, 1000.0)
        ! Find the mixing ratio with the average of the lowest level pressure and MSL (line 609)
        rs0 = eps_const * mixing_ratio_kgkg_profile(nk) * msl_pressure_hpa / &
              (pm * (eps_const + mixing_ratio_kgkg_profile(nk)) - mixing_ratio_kgkg_profile(nk) * msl_pressure_hpa)

        call calculate_cape(temp_kelvin_profile(nk), rs0, pm, &
                                  temp_kelvin_profile, mixing_ratio_kgkg_profile, pressure_levels_hpa, &
                                  num_levels, ascent_flag_local, miss_handle_local, &
                                  capem, toms, lnbs, ifl_cape)

        if (ifl_cape /= 1) then
            status_flag = ifl_cape
            return
        end if

        ! *** Find saturation CAPE at radius of maximum winds *** (lines 620-631)
        call mixing_ratio_from_vapor_pressure(es0, pm, rs0)
        call calculate_cape(sst_kelvin, rs0, pm, &
                                  temp_kelvin_profile, mixing_ratio_kgkg_profile, pressure_levels_hpa, &
                                  num_levels, ascent_flag_local, miss_handle_local, &
                                  capems, toms, lnbs, ifl_cape)

        if (ifl_cape /= 1) then
            status_flag = ifl_cape
            return
        end if

        ! Store outflow temperature and level (lines 630-631)
        outflow_temp = toms
        outflow_level = lnbs

        ! Calculate the proxy for TC efficiency (BE02, EQN. 1-3) (lines 632-636)
        ! Use multiplication instead of division
        rat = sst_kelvin * (1.0 / toms)  ! Compiler may optimize this
        if (diss_flag_local == 0) then
            rat = 1.0
        end if

        ! Initial estimate of pressure at radius of maximum winds (lines 641-651)
        ! Density temperatures - inline calculation for performance
        ! tv0 = T*(1+R/eps)/(1+RT) where RT=R for saturated parcel
        ! Use multiplication by reciprocal instead of division
        tv0 = temp_kelvin_profile(nk) * (1.0 + mixing_ratio_kgkg_profile(nk)*inv_eps) * &
              (1.0 / (1.0 + mixing_ratio_kgkg_profile(nk)))
        tvsst = sst_kelvin * (1.0 + rs0*inv_eps) * (1.0 / (1.0 + rs0))
        tvav = 0.5 * (tv0 + tvsst)

        ! Converge toward CAPE*-CAPEM (BE02, EQN 3-4)
        cat = (capem - capea) + 0.5 * ckcd_local * rat * (capems - capem)
        cat = max(cat, 0.0)

        ! Iterate on pressure - use pre-computed inverse
        pnew = msl_pressure_hpa * exp(-cat * inv_rd / tvav)

        ! Update iteration variables (lines 655-661)
        pmold = pm        ! store the previous step's pressure
        pm = pnew        ! store the current step's pressure
        np = np + 1      ! increase iteration count

        ! *** If the routine does not converge, set IFL=0 and return missing PI *** (lines 664-672)
        if (np > max_iterations .or. pm < min_pressure_limit) then
            max_wind_speed = missing_value
            min_pressure_hpa = missing_value
            status_flag = 0  ! IFL=0
            outflow_temp = missing_value
            outflow_level = missing_value
            return
        end if

    end do

    ! Once converged, calculate final values following tcpyPI (lines 674-687)
    catfac = 0.5 * (1.0 + inv_b)  ! Use pre-computed inverse
    cat = (capem - capea) + ckcd_local * rat * catfac * (capems - capem)
    cat = max(cat, 0.0)

    ! Calculate the minimum pressure at the eye of the storm (BE02 EQN. 4)
    min_pressure_hpa = msl_pressure_hpa * exp(-cat * inv_rd / tvav)

    ! Calculate the potential intensity (BE02 EQN. 3)
    fac = max(0.0, (capems - capem))
    max_wind_speed = v_reduc_local * sqrt(ckcd_local * rat * fac)

    ! Success
    status_flag = 1

    ! Deallocate arrays
    if (allocated(temp_kelvin_profile)) deallocate(temp_kelvin_profile)
    if (allocated(mixing_ratio_kgkg_profile)) deallocate(mixing_ratio_kgkg_profile)

end subroutine calculate_potential_intensity

!--------------------------------------------------------------------------
! CAPE Calculation Subroutine (following tcpyPI cape() function exactly)
!--------------------------------------------------------------------------

subroutine calculate_cape(tp, rp, pp, t, r, p, nlvl, ascent_flag, miss_handle, &
                                 caped, tob, lnb, iflag)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate Convective Available Potential Energy (CAPE) of a parcel
    !   given parcel and environmental conditions following Emanuel (1994)
    !
    ! REFERENCE:
    !   Emanuel, K. (1994): Atmospheric Convection, Equation 6.3.6
    !
    ! ALGORITHM:
    !   1. Lift parcel dry adiabatically to LCL
    !   2. Continue lifting moist adiabatically above LCL
    !   3. Integrate positive buoyancy to find CAPE
    !   4. Find level of neutral buoyancy (LNB)
    !---------------------------------------------------------------------------
    implicit none

    !---------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    !---------------------------------------------------------------------------
    real, intent(in) :: tp
        ! Parcel temperature [K] - Initial temperature of lifted parcel

    real, intent(in) :: rp
        ! Parcel mixing ratio [kg/kg] - Initial water vapor mixing ratio

    real, intent(in) :: pp
        ! Parcel pressure [hPa] - Initial pressure level of parcel

    integer, intent(in) :: nlvl
        ! Number of vertical levels in arrays

    real, intent(in) :: t(nlvl)
        ! Environmental temperature profile [K] - Ordered from surface to top

    real, intent(in) :: r(nlvl)
        ! Environmental mixing ratio profile [kg/kg] - Water vapor at each level

    real, intent(in) :: p(nlvl)
        ! Environmental pressure levels [hPa] - Monotonically decreasing

    integer, intent(in) :: ascent_flag
        ! Type of parcel ascent:
        !   0 = Reversible (retains all water)
        !   1 = Pseudo-adiabatic (immediate precipitation)

    integer, intent(in) :: miss_handle
        ! Missing value handling: 0=Skip NaN, 1=Return missing if NaN

    !---------------------------------------------------------------------------
    ! OUTPUT PARAMETERS:
    !---------------------------------------------------------------------------
    real, intent(out) :: caped
        ! CAPE [J/kg] - Integrated positive buoyancy of parcel

    real, intent(out) :: tob
        ! Temperature at LNB [K] - Outflow temperature of parcel

    real, intent(out) :: lnb
        ! LNB pressure [hPa] - Pressure where buoyancy becomes zero

    integer, intent(out) :: iflag
        ! Status: 1=Success, 0=Invalid, 2=No convergence, 3=Missing values

    ! Local variables exactly matching tcpyPI variable names (lines 94-299)
    real :: tvrdif(nlvl)                         ! Virtual temperature difference array
    real :: tpc, esp, evp, rh, s, plcl           ! Parcel properties and entropy
    real :: pa, na, pat                          ! Positive area, negative area, residual area
    real :: tg, rg, tlvr, tvenv, tenv            ! Parcel and environmental properties
    real :: rmean                                 ! Mixed parcel mixing ratio
    real :: pfac, pma, pinb                      ! CAPE integration factors
    integer :: first_valid, n_levels, j, jmin, inb
    logical :: valid_i(nlvl)
    integer :: nc_temp, convergence_status

    ! Physical constants
    real, parameter :: rd_const = 287.04
    real, parameter :: cpd_const = 1005.7
    real, parameter :: missing_val = -9999.0
    real, parameter :: small_mixing = 1.0e-6
    real, parameter :: min_temp = 200.0
    real, parameter :: eps_const = 0.622         ! Rd/Rv for optimization
    real, parameter :: inv_eps = 1.607717        ! 1.0/0.622 pre-computed for speed
    real, parameter :: rd_cpd_ratio = 0.28571    ! rd_const/cpd_const pre-computed
    real, parameter :: ptop = 50.0               ! Fixed top pressure [hPa]

    ! Initialize outputs
    caped = 0.0
    tob = t(1)
    lnb = p(1)
    iflag = 1

    ! Initialize loop variables
    jmin = 1000000
    inb = 0
    n_levels = nlvl
    first_valid = 1

    ! Check for missing values (lines 95-120 in Python)
    valid_i(:) = .true.
    do j = 1, nlvl
        ! Check for NaN using standard Fortran approach
        if (t(j) /= t(j) .or. r(j) /= r(j)) then
            valid_i(j) = .false.
        end if
    end do

    ! Find first valid level
    do j = 1, nlvl
        if (valid_i(j)) then
            first_valid = j
            exit
        end if
    end do

    ! Handle missing values according to flag
    if (count(valid_i) /= nlvl) then
        if (miss_handle /= 0) then
            caped = missing_val
            tob = missing_val
            lnb = missing_val
            iflag = 3
            return
        end if
    end if

    ! Find level closest to ptop (line 124)
    n_levels = nlvl
    do j = 1, nlvl
        if (abs(p(j) - ptop) < abs(p(n_levels) - ptop)) then
            n_levels = j
        end if
    end do

    ! Check pressure ordering (lines 136-143)
    if (n_levels >= 3) then
        if (p(2) - p(1) > 0.0) then
            caped = 0.0
            tob = missing_val
            lnb = missing_val
            iflag = 0
            return
        end if
    end if

    ! Check parcel suitability (lines 145-152)
    if (rp < small_mixing .or. tp < min_temp) then
        caped = 0.0
        tob = missing_val
        lnb = missing_val
        iflag = 0
        return
    end if

    ! Calculate parcel properties (lines 158-176)
    call temperature_kelvin_to_celsius(tp, tpc)
    call saturation_vapor_pressure_magnus(tpc, esp)
    call vapor_pressure_from_mixing_ratio(rp, pp, evp)
    rh = min(evp / esp, 1.0)
    call entropy_s(tp, rp, pp, s)
    call e_plcl(tp, rh, pp, plcl)

    ! Initialize CAPE calculation variables
    pa = 0.0
    na = 0.0
    tvrdif(:) = 0.0

    ! Main updraft loop (lines 190-232)
    do j = first_valid, n_levels
        jmin = min(jmin, j)

        if (p(j) >= plcl) then
            ! Below LCL (lines 199-208)
            tg = tp * (p(j) / pp) ** rd_cpd_ratio  ! Use pre-computed ratio
            rg = rp
            ! Inline density temperature calculation for performance
            ! Pre-compute reciprocals to avoid division
            tlvr = tg * (1.0 + rg*inv_eps) * (1.0 / (1.0 + rg))
            tvenv = t(j) * (1.0 + r(j)*inv_eps) * (1.0 / (1.0 + r(j)))
            tvrdif(j) = tlvr - tvenv
        else
            ! Above LCL (lines 213-231)
            call solve_temperature_from_entropy(s, p(j), rp, t(j), &
                                                       tg, rg, convergence_status)
            if (convergence_status == 2) then
                caped = 0.0
                tob = t(1)
                lnb = p(1)
                iflag = 2
                return
            end if

            ! Calculate buoyancy
            rmean = real(ascent_flag) * rg + (1.0 - real(ascent_flag)) * rp
            ! Inline density temperature calculation for performance
            ! Pre-compute reciprocals to avoid division
            tlvr = tg * (1.0 + rg*inv_eps) * (1.0 / (1.0 + rmean))
            tenv = t(j) * (1.0 + r(j)*inv_eps) * (1.0 / (1.0 + r(j)))
            tvrdif(j) = tlvr - tenv
        end if
    end do

    ! Find maximum level of positive buoyancy (lines 243-246)
    inb = 0
    do j = n_levels, jmin, -1
        if (tvrdif(j) > 0.0) then
            inb = max(inb, j)
        end if
    end do

    ! Check if LNB is above surface (lines 248-256)
    if (inb == 0) then
        caped = 0.0
        tob = t(1)
        lnb = 0.0
        iflag = 1
        return
    end if

    ! Calculate positive and negative areas (lines 265-268)
    ! OpenMP SIMD optimization for vectorized accumulation
    !$OMP SIMD REDUCTION(+:pa) REDUCTION(+:na) PRIVATE(pfac)
    do j = jmin + 1, inb
        ! Pre-compute reciprocal to avoid division
        pfac = rd_const * (tvrdif(j) + tvrdif(j-1)) * (p(j-1) - p(j)) * (1.0 / (p(j) + p(j-1)))
        pa = pa + max(pfac, 0.0)
        na = na - min(pfac, 0.0)
    end do
    !$OMP END SIMD

    ! Area between parcel pressure and first level (lines 273-276)
    pma = pp + p(jmin)
    pfac = rd_const * (pp - p(jmin)) * (1.0 / pma)
    pa = pa + pfac * max(tvrdif(jmin), 0.0)
    na = na - pfac * min(tvrdif(jmin), 0.0)

    ! Find residual positive area and LNB (lines 282-289)
    pat = 0.0
    tob = t(inb)
    lnb = p(inb)

    if (inb < n_levels) then
        pinb = (p(inb+1) * tvrdif(inb) - p(inb) * tvrdif(inb+1)) / &
               (tvrdif(inb) - tvrdif(inb+1))
        lnb = pinb
        pat = rd_const * tvrdif(inb) * (p(inb) - pinb) / (p(inb) + pinb)
        tob = (t(inb) * (pinb - p(inb+1)) + t(inb+1) * (p(inb) - pinb)) / &
              (p(inb) - p(inb+1))
    end if

    ! Final CAPE calculation (lines 294-297)
    caped = pa + pat - na
    caped = max(caped, 0.0)
    iflag = 1

end subroutine calculate_cape

!--------------------------------------------------------------------------
! 3D Gridded Data Processing (spatial grid at each time)
!--------------------------------------------------------------------------

subroutine calculate_pi_3d_grid(sst_grid, msl_grid, pressure_levels, &
                                temp_grid, mixing_ratio_grid, &
                                nlat, nlon, nlevels, &
                                vmax_grid, pmin_grid, &
                                missing_value, ckcd, ascent_flag, &
                                diss_flag, v_reduc)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate potential intensity for 3D gridded atmospheric data
    !   Processing spatial grid (lat x lon) with vertical profiles
    !
    ! NOTE ON UNITS - DIFFERENT FROM tropical_cyclone_potential_intensity.f90:
    !   This routine uses meteorological units matching tcpyPI:
    !   - Temperature: °C (not K)
    !   - Pressure: hPa (not Pa)
    !   - Mixing ratio: g/kg (not kg/kg)
    !
    ! DIMENSIONS:
    !   Input: (nlevels, nlat, nlon) for atmospheric fields
    !   Output: (nlat, nlon) for PI fields
    !---------------------------------------------------------------------------
    implicit none

    ! Interface for subroutine with optional arguments
    interface
        subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                                pressure_levels_hpa, temp_celsius_profile, &
                                                mixing_ratio_gkg_profile, num_levels, &
                                                max_wind_speed, min_pressure_hpa, &
                                                outflow_temp, outflow_level, &
                                                status_flag, ckcd, ascent_flag, diss_flag, &
                                                v_reduc, miss_handle)
            integer, intent(in) :: num_levels
            real, intent(in) :: sst_celsius
            real, intent(in) :: msl_pressure_hpa
            real, intent(in) :: pressure_levels_hpa(num_levels)
            real, intent(in) :: temp_celsius_profile(num_levels)
            real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
            real, intent(out) :: max_wind_speed
            real, intent(out) :: min_pressure_hpa
            real, intent(out) :: outflow_temp
            real, intent(out) :: outflow_level
            integer, intent(out) :: status_flag
            real, intent(in), optional :: ckcd
            integer, intent(in), optional :: ascent_flag
            integer, intent(in), optional :: diss_flag
            real, intent(in), optional :: v_reduc
            integer, intent(in), optional :: miss_handle
        end subroutine calculate_potential_intensity
    end interface

    !---------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    !---------------------------------------------------------------------------
    integer, intent(in) :: nlat, nlon, nlevels
        ! Grid dimensions: latitude, longitude, vertical levels

    real, intent(in) :: sst_grid(nlat, nlon)
        ! Sea surface temperature [°C] - NOT Kelvin!
        ! Missing values should be set to missing_value

    real, intent(in) :: msl_grid(nlat, nlon)
        ! Mean sea level pressure [hPa] - NOT Pa!

    real, intent(in) :: pressure_levels(nlevels)
        ! Pressure levels [hPa], same for all grid points
        ! Ordered from surface to top (decreasing)

    real, intent(in) :: temp_grid(nlevels, nlat, nlon)
        ! Temperature [°C] at each level and grid point - NOT Kelvin!

    real, intent(in) :: mixing_ratio_grid(nlevels, nlat, nlon)
        ! Water vapor mixing ratio [g/kg] - NOT kg/kg!

    real, intent(in) :: missing_value
        ! Value indicating missing data (e.g., -9999.0)

    real, intent(in), optional :: ckcd
        ! Ck/Cd ratio (default: 0.9)

    integer, intent(in), optional :: ascent_flag
        ! 0=reversible, 1=pseudo-adiabatic (default: 0)

    integer, intent(in), optional :: diss_flag
        ! 0=no dissipative heating, 1=include (default: 1)

    real, intent(in), optional :: v_reduc
        ! Wind reduction factor (default: 0.8)

    !---------------------------------------------------------------------------
    ! OUTPUT PARAMETERS:
    !---------------------------------------------------------------------------
    real, intent(out) :: vmax_grid(nlat, nlon)
        ! Maximum potential wind speed [m/s] at each grid point

    real, intent(out) :: pmin_grid(nlat, nlon)
        ! Minimum central pressure [hPa] at each grid point

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES:
    !---------------------------------------------------------------------------
    integer :: i, j, k, status_flag
    real :: temp_profile(nlevels), mixing_profile(nlevels)
    real :: outflow_temp, outflow_level
    logical :: valid_point

    ! Process optional parameters
    real :: ckcd_local, v_reduc_local
    integer :: ascent_local, diss_local
    real, parameter :: ptop_local = 50.0  ! Fixed top pressure [hPa]

    ! Set defaults for optional parameters
    ckcd_local = 0.9
    if (present(ckcd)) ckcd_local = ckcd

    ascent_local = 0
    if (present(ascent_flag)) ascent_local = ascent_flag

    diss_local = 1
    if (present(diss_flag)) diss_local = diss_flag

    v_reduc_local = 0.8
    if (present(v_reduc)) v_reduc_local = v_reduc

    ! ptop_local is now a fixed parameter = 50.0 hPa

    ! Initialize output arrays with missing values
    vmax_grid(:,:) = missing_value
    pmin_grid(:,:) = missing_value

    ! Loop over all grid points
    ! Note: OpenMP parallelization removed - parallelize at higher level for better efficiency
    do j = 1, nlon
        do i = 1, nlat

            ! Check for missing SST or MSL
            if (abs(sst_grid(i,j) - missing_value) < 1.0 .or. &
                abs(msl_grid(i,j) - missing_value) < 1.0) then
                cycle  ! Skip this point
            end if

            ! Check if SST is reasonable for TC formation
            if (sst_grid(i,j) < 5.0 .or. sst_grid(i,j) > 100.0) then
                cycle  ! Skip unrealistic SST values (changed from 5.0 to 0.0 to match tcpyPI)
            end if

            ! Extract vertical profile at this grid point
            valid_point = .true.
            do k = 1, nlevels
                temp_profile(k) = temp_grid(k, i, j)
                mixing_profile(k) = mixing_ratio_grid(k, i, j)

                ! Check for missing values in profile
                if (abs(temp_profile(k) - missing_value) < 1.0 .or. &
                    abs(mixing_profile(k) - missing_value) < 1.0) then
                    ! Set to zero for upper levels (common practice)
                    if (pressure_levels(k) < 500.0) then
                        mixing_profile(k) = 0.0
                    else
                        valid_point = .false.
                        exit
                    end if
                end if
            end do

            if (.not. valid_point) cycle

            ! Calculate PI for this grid point
            call calculate_potential_intensity( &
                sst_grid(i,j), msl_grid(i,j), &
                pressure_levels, temp_profile, mixing_profile, &
                nlevels, &
                vmax_grid(i,j), pmin_grid(i,j), &
                outflow_temp, outflow_level, status_flag, &
                ckcd_local, ascent_local, diss_local, &
                v_reduc_local, 1)

            ! Check for calculation failure
            if (status_flag /= 1) then
                vmax_grid(i,j) = missing_value
                pmin_grid(i,j) = missing_value
            end if

        end do
    end do

end subroutine calculate_pi_3d_grid

!--------------------------------------------------------------------------
! 4D Gridded Data Processing (time series of spatial grids)
!--------------------------------------------------------------------------

subroutine calculate_pi_4d_grid(sst_grid, msl_grid, pressure_levels, &
                                temp_grid, mixing_ratio_grid, &
                                ntime, nlat, nlon, nlevels, &
                                vmax_grid, pmin_grid, &
                                missing_value, ckcd, ascent_flag, &
                                diss_flag, v_reduc)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate potential intensity for 4D gridded atmospheric data
    !   Processing time series of spatial grids with vertical profiles
    !
    ! NOTE ON UNITS - DIFFERENT FROM tropical_cyclone_potential_intensity.f90:
    !   This routine uses meteorological units matching tcpyPI:
    !   - Temperature: °C (not K)
    !   - Pressure: hPa (not Pa)
    !   - Mixing ratio: g/kg (not kg/kg)
    !
    ! DIMENSIONS:
    !   Input: (ntime, nlevels, nlat, nlon) for atmospheric fields
    !   Output: (ntime, nlat, nlon) for PI fields
    !---------------------------------------------------------------------------
    implicit none

    ! Interface for subroutine with optional arguments
    interface
        subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                                pressure_levels_hpa, temp_celsius_profile, &
                                                mixing_ratio_gkg_profile, num_levels, &
                                                max_wind_speed, min_pressure_hpa, &
                                                outflow_temp, outflow_level, &
                                                status_flag, ckcd, ascent_flag, diss_flag, &
                                                v_reduc, miss_handle)
            integer, intent(in) :: num_levels
            real, intent(in) :: sst_celsius
            real, intent(in) :: msl_pressure_hpa
            real, intent(in) :: pressure_levels_hpa(num_levels)
            real, intent(in) :: temp_celsius_profile(num_levels)
            real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
            real, intent(out) :: max_wind_speed
            real, intent(out) :: min_pressure_hpa
            real, intent(out) :: outflow_temp
            real, intent(out) :: outflow_level
            integer, intent(out) :: status_flag
            real, intent(in), optional :: ckcd
            integer, intent(in), optional :: ascent_flag
            integer, intent(in), optional :: diss_flag
            real, intent(in), optional :: v_reduc
            integer, intent(in), optional :: miss_handle
        end subroutine calculate_potential_intensity
    end interface

    !---------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    !---------------------------------------------------------------------------
    integer, intent(in) :: ntime, nlat, nlon, nlevels
        ! Grid dimensions: time, latitude, longitude, vertical levels

    real, intent(in) :: sst_grid(ntime, nlat, nlon)
        ! Sea surface temperature [°C] time series - NOT Kelvin!

    real, intent(in) :: msl_grid(ntime, nlat, nlon)
        ! Mean sea level pressure [hPa] time series - NOT Pa!

    real, intent(in) :: pressure_levels(nlevels)
        ! Pressure levels [hPa], constant in time

    real, intent(in) :: temp_grid(ntime, nlevels, nlat, nlon)
        ! Temperature [°C] 4D field - NOT Kelvin!

    real, intent(in) :: mixing_ratio_grid(ntime, nlevels, nlat, nlon)
        ! Water vapor mixing ratio [g/kg] 4D field - NOT kg/kg!

    real, intent(in) :: missing_value
        ! Value indicating missing data

    real, intent(in), optional :: ckcd
    integer, intent(in), optional :: ascent_flag, diss_flag
    real, intent(in), optional :: v_reduc

    !---------------------------------------------------------------------------
    ! OUTPUT PARAMETERS:
    !---------------------------------------------------------------------------
    real, intent(out) :: vmax_grid(ntime, nlat, nlon)
        ! Maximum potential wind speed [m/s] time series

    real, intent(out) :: pmin_grid(ntime, nlat, nlon)
        ! Minimum central pressure [hPa] time series

    !---------------------------------------------------------------------------
    ! LOCAL VARIABLES:
    !---------------------------------------------------------------------------
    integer :: t, i, j, k, status_flag
    real :: temp_profile(nlevels), mixing_profile(nlevels)
    real :: outflow_temp, outflow_level
    logical :: valid_point

    ! Process optional parameters
    real :: ckcd_local, v_reduc_local
    integer :: ascent_local, diss_local
    real, parameter :: ptop_local = 50.0  ! Fixed top pressure [hPa]

    ! Set defaults
    ckcd_local = 0.9
    if (present(ckcd)) ckcd_local = ckcd

    ascent_local = 0
    if (present(ascent_flag)) ascent_local = ascent_flag

    diss_local = 1
    if (present(diss_flag)) diss_local = diss_flag

    v_reduc_local = 0.8
    if (present(v_reduc)) v_reduc_local = v_reduc

    ! ptop_local is now a fixed parameter = 50.0 hPa

    ! Initialize output arrays
    vmax_grid(:,:,:) = missing_value
    pmin_grid(:,:,:) = missing_value

    ! Loop over all time steps and grid points
    ! MEMORY-OPTIMAL LOOP ORDER: j->i->t (following tropical_cyclone_potential_intensity.f90)
    ! Based on Fortran column-major storage and cache optimization
    ! Note: OpenMP parallelization removed - parallelize at higher level for better efficiency
    do j = 1, nlon          ! longitude outermost
        do i = 1, nlat      ! latitude middle
            do t = 1, ntime ! time innermost

                ! Check for missing SST or MSL
                if (abs(sst_grid(t,i,j) - missing_value) < 1.0 .or. &
                    abs(msl_grid(t,i,j) - missing_value) < 1.0) then
                    cycle
                end if

                ! Check SST validity
                if (sst_grid(t,i,j) < 5.0 .or. sst_grid(t,i,j) > 100.0) then
                    cycle
                end if

                ! Extract vertical profile
                valid_point = .true.
                do k = 1, nlevels
                    temp_profile(k) = temp_grid(t, k, i, j)
                    mixing_profile(k) = mixing_ratio_grid(t, k, i, j)

                    ! Handle missing values
                    if (abs(temp_profile(k) - missing_value) < 1.0 .or. &
                        abs(mixing_profile(k) - missing_value) < 1.0) then
                        if (pressure_levels(k) < 500.0) then
                            mixing_profile(k) = 0.0
                        else
                            valid_point = .false.
                            exit
                        end if
                    end if
                end do

                if (.not. valid_point) cycle

                ! Calculate PI
                call calculate_potential_intensity( &
                    sst_grid(t,i,j), msl_grid(t,i,j), &
                    pressure_levels, temp_profile, mixing_profile, &
                    nlevels, &
                    vmax_grid(t,i,j), pmin_grid(t,i,j), &
                    outflow_temp, outflow_level, status_flag, &
                    ckcd_local, ascent_local, diss_local, &
                    v_reduc_local, 1)

                ! Handle calculation failure
                if (status_flag /= 1) then
                    vmax_grid(t,i,j) = missing_value
                    pmin_grid(t,i,j) = missing_value
                end if

            end do
        end do
    end do

end subroutine calculate_pi_4d_grid

!--------------------------------------------------------------------------
! 3D Gridded Data Processing WITH MISSING VALUE HANDLING
!--------------------------------------------------------------------------

subroutine calculate_pi_3d_grid_with_missing(sst_grid, msl_grid, pressure_levels, &
                                            temp_grid, mixing_ratio_grid, &
                                            nlat, nlon, nlevels, &
                                            vmax_grid, pmin_grid, missing_value)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate potential intensity for 3D gridded data with optimized
    !   missing value handling (e.g., land points, missing data)
    !
    ! OPTIMIZATION:
    !   - Skips computation for land points (missing SST)
    !   - Handles partial vertical profiles
    !   - Pre-checks for physical constraints
    !---------------------------------------------------------------------------
    implicit none

    ! Interface for subroutine with optional arguments
    interface
        subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                                pressure_levels_hpa, temp_celsius_profile, &
                                                mixing_ratio_gkg_profile, num_levels, &
                                                max_wind_speed, min_pressure_hpa, &
                                                outflow_temp, outflow_level, &
                                                status_flag, ckcd, ascent_flag, diss_flag, &
                                                v_reduc, miss_handle)
            integer, intent(in) :: num_levels
            real, intent(in) :: sst_celsius
            real, intent(in) :: msl_pressure_hpa
            real, intent(in) :: pressure_levels_hpa(num_levels)
            real, intent(in) :: temp_celsius_profile(num_levels)
            real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
            real, intent(out) :: max_wind_speed
            real, intent(out) :: min_pressure_hpa
            real, intent(out) :: outflow_temp
            real, intent(out) :: outflow_level
            integer, intent(out) :: status_flag
            real, intent(in), optional :: ckcd
            integer, intent(in), optional :: ascent_flag
            integer, intent(in), optional :: diss_flag
            real, intent(in), optional :: v_reduc
            integer, intent(in), optional :: miss_handle
        end subroutine calculate_potential_intensity
    end interface

    integer, intent(in) :: nlat, nlon, nlevels
    real, intent(in) :: sst_grid(nlat, nlon)
    real, intent(in) :: msl_grid(nlat, nlon)
    real, intent(in) :: pressure_levels(nlevels)
    real, intent(in) :: temp_grid(nlevels, nlat, nlon)
    real, intent(in) :: mixing_ratio_grid(nlevels, nlat, nlon)
    real, intent(in) :: missing_value

    real, intent(out) :: vmax_grid(nlat, nlon)
    real, intent(out) :: pmin_grid(nlat, nlon)

    ! Local variables
    integer :: i, j, k, status_flag, valid_start
    real :: temp_profile(nlevels), mixing_profile(nlevels)
    real :: outflow_temp, outflow_level

    ! Initialize outputs
    vmax_grid(:,:) = missing_value
    pmin_grid(:,:) = missing_value

    ! Process each grid point
    do j = 1, nlon
        do i = 1, nlat
            ! OPTIMIZATION: Check SST first (skip land points)
            if (abs(sst_grid(i,j) - missing_value) < 1.0) then
                cycle  ! Skip land point
            end if

            ! Check physical constraint
            if (sst_grid(i,j) < 5.0 .or. sst_grid(i,j) > 100.0) then
                cycle
            end if

            ! Check MSL validity
            if (abs(msl_grid(i,j) - missing_value) < 1.0) then
                cycle
            end if

            ! Find valid vertical data range
            valid_start = nlevels + 1
            do k = 1, nlevels
                if (abs(temp_grid(k,i,j) - missing_value) > 1.0 .and. &
                    abs(mixing_ratio_grid(k,i,j) - missing_value) > 1.0) then
                    valid_start = min(valid_start, k)
                    temp_profile(k) = temp_grid(k, i, j)
                    mixing_profile(k) = mixing_ratio_grid(k, i, j)
                else if (pressure_levels(k) < 500.0) then
                    ! Set upper levels to zero mixing ratio
                    temp_profile(k) = temp_grid(k, i, j)
                    mixing_profile(k) = 0.0
                else
                    temp_profile(k) = missing_value
                    mixing_profile(k) = missing_value
                end if
            end do

            ! Skip if no valid data
            if (valid_start > nlevels) cycle

            ! Calculate PI with valid subset
            call calculate_potential_intensity( &
                sst_grid(i,j), msl_grid(i,j), &
                pressure_levels(valid_start:), &
                temp_profile(valid_start:), &
                mixing_profile(valid_start:), &
                nlevels - valid_start + 1, &
                vmax_grid(i,j), pmin_grid(i,j), &
                outflow_temp, outflow_level, status_flag, &
                0.9, 0, 1, 0.8, 1)

            ! Check status
            if (status_flag /= 1) then
                vmax_grid(i,j) = missing_value
                pmin_grid(i,j) = missing_value
            end if
        end do
    end do

end subroutine calculate_pi_3d_grid_with_missing

!--------------------------------------------------------------------------
! 4D Gridded Data Processing WITH MISSING VALUE HANDLING
!--------------------------------------------------------------------------

subroutine calculate_pi_4d_grid_with_missing(sst_grid, msl_grid, pressure_levels, &
                                            temp_grid, mixing_ratio_grid, &
                                            ntime, nlat, nlon, nlevels, &
                                            vmax_grid, pmin_grid, missing_value)
    !---------------------------------------------------------------------------
    ! PURPOSE:
    !   Calculate potential intensity for 4D gridded data with optimized
    !   missing value handling for time series
    !---------------------------------------------------------------------------
    implicit none

    ! Interface for subroutine with optional arguments
    interface
        subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                                pressure_levels_hpa, temp_celsius_profile, &
                                                mixing_ratio_gkg_profile, num_levels, &
                                                max_wind_speed, min_pressure_hpa, &
                                                outflow_temp, outflow_level, &
                                                status_flag, ckcd, ascent_flag, diss_flag, &
                                                v_reduc, miss_handle)
            integer, intent(in) :: num_levels
            real, intent(in) :: sst_celsius
            real, intent(in) :: msl_pressure_hpa
            real, intent(in) :: pressure_levels_hpa(num_levels)
            real, intent(in) :: temp_celsius_profile(num_levels)
            real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
            real, intent(out) :: max_wind_speed
            real, intent(out) :: min_pressure_hpa
            real, intent(out) :: outflow_temp
            real, intent(out) :: outflow_level
            integer, intent(out) :: status_flag
            real, intent(in), optional :: ckcd
            integer, intent(in), optional :: ascent_flag
            integer, intent(in), optional :: diss_flag
            real, intent(in), optional :: v_reduc
            integer, intent(in), optional :: miss_handle
        end subroutine calculate_potential_intensity
    end interface

    integer, intent(in) :: ntime, nlat, nlon, nlevels
    real, intent(in) :: sst_grid(ntime, nlat, nlon)
    real, intent(in) :: msl_grid(ntime, nlat, nlon)
    real, intent(in) :: pressure_levels(nlevels)
    real, intent(in) :: temp_grid(ntime, nlevels, nlat, nlon)
    real, intent(in) :: mixing_ratio_grid(ntime, nlevels, nlat, nlon)
    real, intent(in) :: missing_value

    real, intent(out) :: vmax_grid(ntime, nlat, nlon)
    real, intent(out) :: pmin_grid(ntime, nlat, nlon)

    ! Local variables
    integer :: t, i, j, k, status_flag, valid_start
    real :: temp_profile(nlevels), mixing_profile(nlevels)
    real :: outflow_temp, outflow_level

    ! Initialize outputs
    vmax_grid(:,:,:) = missing_value
    pmin_grid(:,:,:) = missing_value

    ! Memory-optimal loop order
    do j = 1, nlon
        do i = 1, nlat
            do t = 1, ntime
                ! OPTIMIZATION: Check SST first
                if (abs(sst_grid(t,i,j) - missing_value) < 1.0) then
                    cycle
                end if

                ! Check physical constraint
                if (sst_grid(t,i,j) < 5.0 .or. sst_grid(t,i,j) > 100.0) then
                    cycle
                end if

                ! Check MSL validity
                if (abs(msl_grid(t,i,j) - missing_value) < 1.0) then
                    cycle
                end if

                ! Find valid vertical data range
                valid_start = nlevels + 1
                do k = 1, nlevels
                    if (abs(temp_grid(t,k,i,j) - missing_value) > 1.0 .and. &
                        abs(mixing_ratio_grid(t,k,i,j) - missing_value) > 1.0) then
                        valid_start = min(valid_start, k)
                        temp_profile(k) = temp_grid(t, k, i, j)
                        mixing_profile(k) = mixing_ratio_grid(t, k, i, j)
                    else if (pressure_levels(k) < 500.0) then
                        temp_profile(k) = temp_grid(t, k, i, j)
                        mixing_profile(k) = 0.0
                    else
                        temp_profile(k) = missing_value
                        mixing_profile(k) = missing_value
                    end if
                end do

                ! Skip if no valid data
                if (valid_start > nlevels) cycle

                ! Calculate PI with valid subset
                call calculate_potential_intensity( &
                    sst_grid(t,i,j), msl_grid(t,i,j), &
                    pressure_levels(valid_start:), &
                    temp_profile(valid_start:), &
                    mixing_profile(valid_start:), &
                    nlevels - valid_start + 1, &
                    vmax_grid(t,i,j), pmin_grid(t,i,j), &
                    outflow_temp, outflow_level, status_flag, &
                    0.9, 0, 1, 0.8, 1)

                ! Check status
                if (status_flag /= 1) then
                    vmax_grid(t,i,j) = missing_value
                    pmin_grid(t,i,j) = missing_value
                end if
            end do
        end do
    end do

end subroutine calculate_pi_4d_grid_with_missing

!--------------------------------------------------------------------------
! Test Program
!--------------------------------------------------------------------------

program test_tcpypi_modern_fortran
    implicit none

    ! Interface for subroutine with optional arguments
    interface
        subroutine calculate_potential_intensity(sst_celsius, msl_pressure_hpa, &
                                                       pressure_levels_hpa, temp_celsius_profile, &
                                                       mixing_ratio_gkg_profile, num_levels, &
                                                       max_wind_speed, min_pressure_hpa, &
                                                       outflow_temp, outflow_level, &
                                                       status_flag, ckcd, ascent_flag, diss_flag, &
                                                       v_reduc, miss_handle)
            integer, intent(in) :: num_levels
            real, intent(in) :: sst_celsius
            real, intent(in) :: msl_pressure_hpa
            real, intent(in) :: pressure_levels_hpa(num_levels)
            real, intent(in) :: temp_celsius_profile(num_levels)
            real, intent(in) :: mixing_ratio_gkg_profile(num_levels)
            real, intent(in), optional :: ckcd
            integer, intent(in), optional :: ascent_flag
            integer, intent(in), optional :: diss_flag
            real, intent(in), optional :: v_reduc
            integer, intent(in), optional :: miss_handle
            real, intent(out) :: max_wind_speed
            real, intent(out) :: min_pressure_hpa
            real, intent(out) :: outflow_temp
            real, intent(out) :: outflow_level
            integer, intent(out) :: status_flag
        end subroutine calculate_potential_intensity
    end interface

    ! Test variables

    ! Test PI calculation variables
    integer, parameter :: test_levels = 31  ! Use same number of levels as pcmin test
    real :: sst_test, msl_pressure_test
    real :: pressure_profile_test(test_levels)
    real :: temp_profile_test(test_levels)
    real :: mixing_profile_test(test_levels)
    real :: max_wind_out, min_pressure_out, outflow_temp_out, outflow_level_out
    integer :: status_out, i

    ! Expected values from test_pcmin.f90
    real, parameter :: EXPECTED_VMAX = 42.15   ! m/s
    real, parameter :: EXPECTED_PMIN = 985.57  ! hPa

    write(*, '(A)') '================================================='
    write(*, '(A)') 'tcpyPI FORTRAN IMPLEMENTATION VALIDATION TEST'
    write(*, '(A)') 'Location: 0°N, 150°W (Equatorial Pacific)'
    write(*, '(A)') 'Using same data as test_pcmin.f90'
    write(*, '(A)') '================================================='
    write(*, *)

    ! Use exact same test data as test_pcmin.f90
    sst_test = 25.73                      ! Sea surface temperature (°C)
    msl_pressure_test = 1009.75            ! Sea level pressure (hPa)

    ! Initialize pressure levels (hPa) - same as test_pcmin.f90
    pressure_profile_test = (/1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
                              750., 725., 700., 650., 600., 550., 500., 450., 400., 350., &
                              300., 250., 200., 150., 100.,  70.,  50.,  40.,  30.,  20., 10./)

    ! Initialize temperature profile (°C) - same as test_pcmin.f90
    temp_profile_test = (/24.73, 22.66, 20.73, 19.19, 18.26, 18.26, 18.25, 17.49, 16.43, 15.14, &
                          13.58, 11.75,  9.82,  6.66,  3.66,  0.11, -3.85, -8.49, -13.96, -20.58, &
                         -29.12, -39.51, -52.10, -67.19, -84.29, -82.89, -71.98, -65.16, -56.62, &
                         -51.35, -44.32/)

    ! Initialize mixing ratio profile (g/kg) - same as test_pcmin.f90
    mixing_profile_test = (/15.573, 15.043, 14.485, 13.549, 12.102,  9.528,  6.931,  5.712,  4.992, &
                            4.770,  5.050,  5.534,  5.743,  4.403,  2.715,  1.946,  1.560,  0.969, &
                            0.507,  0.330,  0.277,  0.159,  0.054,  0.010,  0.001,  0.002,  0.003, &
                            0.003,  0.003,  0.003,  0.003/)

    write(*, '(A)') 'INPUT DATA:'
    write(*, '(A,F8.2,A)') ' SST = ', sst_test, ' °C'
    write(*, '(A,F8.2,A)') ' SLP = ', msl_pressure_test, ' hPa'
    write(*, '(A,I3)')     ' Pressure levels = ', test_levels
    write(*, *)

    write(*, '(A)') 'Calling PI calculation...'

    ! Test PI calculation with tcpyPI units - call with all parameters
    call calculate_potential_intensity(sst_test, msl_pressure_test, &  ! SST in °C, MSL in hPa
                                             pressure_profile_test, temp_profile_test, &
                                             mixing_profile_test, test_levels, &
                                             max_wind_out, min_pressure_out, &
                                             outflow_temp_out, outflow_level_out, &
                                             status_out, 0.9, 0, 1, 0.8, 1)

    write(*, '(A)') 'PI calculation completed.'
    write(*, *)

    write(*, '(A)') 'RESULTS:'
    write(*, '(A)') '--------'
    write(*, '(A, I0)') 'Status Flag (IFL): ', status_out

    if (status_out == 1) then
        write(*, '(A, F8.2, A)') 'VMAX (Maximum Wind Speed): ', max_wind_out, ' m/s'
        write(*, '(A, F8.2, A)') 'PMIN (Minimum Pressure):   ', min_pressure_out, ' hPa'
        write(*, '(A, F8.2, A)') 'TO (Outflow Temperature):  ', outflow_temp_out, ' K'
        write(*, '(A, F8.2, A)') 'OTL (Outflow Level):       ', outflow_level_out, ' hPa'
        write(*, *)

        write(*, '(A)') 'VALIDATION AGAINST EXPECTED VALUES:'
        write(*, '(A)') '-----------------------------------'
        write(*, '(A, F8.2, A)') 'Expected VMAX: ', EXPECTED_VMAX, ' m/s'
        write(*, '(A, F8.2, A)') 'Expected PMIN: ', EXPECTED_PMIN, ' hPa'
        write(*, *)

        write(*, '(A, F8.2, A)') 'Error in VMAX: ', abs(max_wind_out - EXPECTED_VMAX), ' m/s'
        write(*, '(A, F8.2, A)') 'Error in PMIN: ', abs(min_pressure_out - EXPECTED_PMIN), ' hPa'
        write(*, *)

        ! Assessment
        if (abs(max_wind_out - EXPECTED_VMAX) < 5.0 .and. abs(min_pressure_out - EXPECTED_PMIN) < 5.0) then
            write(*, '(A)') 'SUCCESS: Results match expected values within tolerance (< 5 units)'
        else
            write(*, '(A)') 'WARNING: Results deviate from expected values'
        end if
    else
        write(*, '(A)') 'ERROR: PI Calculation failed!'
        write(*, '(A)') 'Debug information:'
        if (status_out == 0) then
            write(*, '(A)') '  Status 0: No convergence or invalid input'
        else if (status_out == 2) then
            write(*, '(A)') '  Status 2: CAPE routine failed to converge'
        else if (status_out == 3) then
            write(*, '(A)') '  Status 3: Missing data in profile'
        end if
    end if

    write(*, *)
    write(*, '(A)') '================================================'
    write(*, '(A)') 'tcpyPI test completed!'

end program test_tcpypi_modern_fortran
