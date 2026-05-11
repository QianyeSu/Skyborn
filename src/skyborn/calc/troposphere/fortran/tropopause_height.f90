! Multi-dimensional tropopause calculation for gridded atmospheric data
! Global subroutines for f2py compatibility
! Processes (lat, lon, level, [time]) arrays to compute tropopause properties
!
! This module contains high-performance Fortran implementations for calculating
! WMO tropopause properties from atmospheric data using the thermal tropopause
! definition. It includes optimized routines for 1D profiles, 2D cross-sections,
! 3D spatial data, and 4D time series analysis.
!
! *** DATA REQUIREMENTS ***
! - CRITICAL: This module requires ISOBARIC (constant pressure level) data
! - Temperature and pressure data must be provided on the same pressure levels
! - Pressure levels MUST be sorted in ASCENDING order:
!   (from low pressure/high altitude to high pressure/low altitude)
!   Example: [10, 20, 50, 100, 200, 300, 500, 700, 850, 1000] hPa
! - For model level data, first interpolate to pressure levels
!
! *** WMO TROPOPAUSE ALGORITHM (STATTROP) ***
!
! Mathematical Formulation:
! ========================
!
! The tropopause calculation follows the WMO (1957) definition:
! "The first tropopause is defined as the lowest level at which
! the lapse rate decreases to 2 deg C per kilometer or less,
! provided also the average lapse rate between this level and
! all higher levels within 2 kilometers does not exceed 2 deg C"
!
! Lapse Rate Calculation:
! ----------------------
! The lapse rate gamma = -dT/dz [K/km] is calculated using the hydrostatic approximation:
!
! From hydrostatic equation: dz = -dp/(g*rho) = -dp/p * R_d/g * T = -d(ln p) * R_d*T/g
!
! Therefore: gamma = -dT/dz = dT/T * g/R_d * 1/d(ln p) = d(ln T)/d(ln p) * CONST
!
! Where: CONST = 1000 * g / R_d = 1000 * 9.80665 / 287.04 ≈ 34.16 K
!        g = gravitational acceleration = 9.80665 m/s²
!        R_d = specific gas constant for dry air = 287.04 J/(kg·K)
!
! For discrete levels:
!   LAPSE(i) = [ln(T(i)) - ln(T(i+1))] / [ln(p(i)) - ln(p(i+1))] * CONST
!   PHALF(i) = [p(i) + p(i+1)] / 2  (mid-level pressure for interpolation)
!
! Tropopause Detection Algorithm:
! ===============================
!
! 1. PRIMARY CRITERION (Level-by-level check):
!    Find lowest level where: LAPSE(i) < LAPSE_CRITERION (default: 2.0 K/km)
!    AND pressure < P_MAX (450 hPa, to exclude boundary layer effects)
!
! 2. SECONDARY CRITERION (2 km layer average):
!    Calculate 2 km pressure difference: P_2KM = P_TROP * exp(-DELTAZ*CONST/T)
!    Where DELTAZ = 2.0 km
!    Average lapse rate in 2 km layer above tropopause must also be < 2.0 K/km
!
! 3. PRESSURE INTERPOLATION:
!    Linear interpolation in log(pressure) space:
!    P1 = ln(PHALF(i)), P2 = ln(PHALF(i+1))
!    WEIGHT = (LAPSE_CRITERION - LAPSE(i)) / (LAPSE(i+1) - LAPSE(i))
!    P_TROP = exp(P1 + WEIGHT * (P2 - P1))
!
! 4. HEIGHT CALCULATION:
!    Using US Standard Atmosphere (1976) barometric formula:
!    For troposphere (h < 11 km): T = T₀ - L*h, where L = 0.0065 K/m
!    h = (T₀ - T) / L = (288.15 - T(P_trop)) / 0.0065
!    More complex calculation for stratosphere using temperature profile
!
! 5. CONSTRAINTS:
!    - Minimum tropopause pressure: 85 hPa (≈ 16 km altitude)
!    - Maximum search pressure: 450 hPa (≈ 6 km altitude)
!    - If no tropopause found, return missing value (-999.0)
!
! Physical Constants Used:
! =======================
! g    = 9.80665 m/s²           (standard gravity)
! R_d  = 287.04 J/(kg·K)        (specific gas constant for dry air)
! CONST = 1000 * g/R_d = 34.16 K (lapse rate conversion factor)
! T₀   = 288.15 K               (sea level standard temperature)
! L    = 0.0065 K/m             (standard temperature lapse rate)
!
! References:
! ==========
! - WMO (1992): International meteorological vocabulary, Genf, 784pp.
! - Randel WJ, Wu F, Gaffen DJ, JGR, 105, 15509-15523, 2000.
! - US Standard Atmosphere, 1976 (for height calculations)
! - Original algorithm: Dominik Brunner, V1.0 Aug 2000
!   Built upon routine stattrop by Peter van Velthoven, KNMI
!   Modified for multi-dimensional processing and optimized with OpenMP
!
! IMPORTANT: Pressure levels must be in ascending order
! (from low pressure/high altitude to high pressure/low altitude)

module tropopause_height_mod
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    implicit none

contains

subroutine tropopause_grid_3d( &
    nlat, nlon, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
)
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: nlat, nlon, nlev, nlevm, punit

    ! INPUT ARRAYS
    real(c_double), intent(in) :: pfull(nlev)                ! 1D pressure levels [hPa or Pa]
    real(c_double), intent(in) :: tfull(nlat, nlon, nlev)    ! Temperature [K]
    real(c_double), intent(in) :: tmsg, lapsec               ! Missing value, lapse criterion

    ! OUTPUT ARRAYS
    real(c_double), intent(out) :: ptrop_hpa(nlat, nlon)     ! Tropopause pressure [hPa]
    real(c_double), intent(out) :: htrop_m(nlat, nlon)       ! Tropopause height [m]
    integer(c_int), intent(out) :: itrop(nlat, nlon)         ! Tropopause level index
    real(c_double), intent(out) :: lapse_rate(nlat, nlon)    ! Tropopause lapse rate [K/km]
    integer(c_int), intent(out) :: success(nlat, nlon)       ! Success flag

    ! LOCAL VARIABLES
    integer :: i, j, k
    integer :: itrop_temp
    real(c_double) :: profile_t(nlev)
    real(c_double) :: work_lapse(nlevm), work_phalf(nlevm), work_pmb(nlev)
    real(c_double) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    !
    ! Calculate WMO tropopause properties for 3D gridded atmospheric data
    !
    ! *** CRITICAL: This routine ONLY accepts ISOBARIC (constant pressure level) data ***
    ! If you have model level data, you MUST interpolate to pressure levels first!
    !
    ! This routine processes spatial atmospheric data (lat, lon, level) to identify
    ! the thermal tropopause at each grid point using the WMO definition.
    !
    ! Input Parameters:
    ! ================
    ! NLAT          : Number of latitude points
    ! NLON          : Number of longitude points
    ! NLEV          : Number of vertical pressure levels
    ! NLEVM         : NLEV + 1 (for work arrays)
    ! PFULL(NLEV)   : 1D array of pressure levels [hPa or Pa]
    !                 *** MUST be in ASCENDING order (high altitude → surface) ***
    !                 Example: [10, 20, 50, 100, 200, 300, 500, 700, 850, 1000] hPa
    !                 *** MUST be ISOBARIC levels (constant pressure) ***
    ! TFULL(NLAT,NLON,NLEV) : Temperature array [K] on isobaric levels
    !                         Same pressure levels as PFULL for all grid points
    ! TMSG          : Missing value for invalid data (typically -999.0)
    ! LAPSEC        : WMO lapse rate criterion [K/km] (typically 2.0)
    ! PUNIT         : Pressure unit flag (0=hPa, 1=Pa)
    !
    ! Output Parameters:
    ! =================
    ! PTROP_HPA(NLAT,NLON)    : Tropopause pressure [hPa]
    ! HTROP_M(NLAT,NLON)      : Tropopause height [m above sea level]
    ! ITROP(NLAT,NLON)        : Tropopause level index (0-based, -999 if not found)
    ! LAPSE_RATE(NLAT,NLON)   : Lapse rate at tropopause [K/km]
    ! SUCCESS(NLAT,NLON)      : .TRUE. if tropopause found, .FALSE. otherwise
    !
    ! Algorithm Implementation (from STATTROP):
    ! ========================================
    ! 1. For each grid point, extract vertical profile from TFULL
    ! 2. Since PFULL is 1D and constant for all grid points (isobaric data),
    !    we use the same pressure levels for all profiles
    ! 3. Calculate lapse rate: LAPSE(i) = CONST * [ln(T(i)) - ln(T(i+1))] / [ln(p(i)) - ln(p(i+1))]
    ! 4. Apply WMO criteria to find tropopause
    ! 5. Interpolate tropopause pressure and calculate height
    !
    ! OpenMP parallelization over grid points for performance
    !
    ! Initialize output arrays
    ptrop_hpa = tmsg
    htrop_m = tmsg
    itrop = -999_c_int
    lapse_rate = tmsg
    success = 0_c_int

    !$OMP PARALLEL DO PRIVATE(i, j, k, itrop_temp, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(nlat, nlon, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    !$OMP ptrop_hpa, htrop_m, itrop, lapse_rate, success) COLLAPSE(2)
    do j = 1, nlon
        do i = 1, nlat
            ! Extract vertical temperature profile for this grid point
            ! No need to copy PFULL since it's the same 1D array for all grid points (isobaric data)
            valid_profile = .true.
            do k = 1, nlev
                profile_t(k) = tfull(i, j, k)
                ! Check for missing values (pressure check done once outside loop for efficiency)
                if (abs(profile_t(k) - tmsg) < 1.0d-10) then
                    valid_profile = .false.
                    exit
                end if
            end do

            ! Check pressure levels for missing values (done once per grid point)
            if (valid_profile) then
                do k = 1, nlev
                    if (abs(pfull(k) - tmsg) < 1.0d-10) then
                        valid_profile = .false.
                        exit
                    end if
                end do
            end if

            ! Process valid profiles only
            if (valid_profile) then
                ! Calculate tropopause pressure and level
                ! Pass PFULL directly (1D isobaric levels, same for all grid points)
                itrop_temp = -999
                call STATTROPX( &
                    nlev, nlevm, pfull, profile_t, tmsg, lapsec, punit, &
                    ptrop_temp, itrop_temp, work_lapse, work_phalf, work_pmb &
                )

                ! Check if tropopause was found
                if (abs(ptrop_temp - tmsg) > 1.0d-10) then
                    ! Convert pressure units if needed
                    if (punit == 0) then
                        ptrop_hpa(i, j) = ptrop_temp
                    else
                        ptrop_hpa(i, j) = ptrop_temp * 0.01d0
                    end if

                    ! Calculate height using US Standard Atmosphere 1976
                    ! DSTDATMP expects array inputs, so we use temporary variables
                    temp_p(1) = ptrop_hpa(i, j)
                    call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
                    htrop_m(i, j) = temp_h(1)
                    itrop(i, j) = int(itrop_temp, c_int)

                    ! Calculate lapse rate at tropopause level
                    if (itrop_temp >= 1 .and. itrop_temp < nlev - 1) then
                        lapse_rate(i, j) = work_lapse(itrop_temp + 1)
                    end if

                    success(i, j) = 1_c_int
                end if
            end if
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine tropopause_grid_3d

subroutine tropopause_grid_4d( &
    nlat, nlon, nlev, ntime, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
)
    !
    ! Calculate WMO tropopause properties for 4D gridded atmospheric time series
    !
    ! *** CRITICAL: This routine ONLY accepts ISOBARIC (constant pressure level) data ***
    ! If you have model level data, you MUST interpolate to pressure levels first!
    !
    ! Input Parameters:
    ! ================
    ! NLAT          : Number of latitude points
    ! NLON          : Number of longitude points
    ! NLEV          : Number of vertical pressure levels
    ! NTIME         : Number of time steps
    ! NLEVM         : NLEV + 1 (for work arrays)
    ! PFULL(NLEV)   : 1D array of pressure levels [hPa or Pa]
    !                 *** MUST be in ASCENDING order (high altitude → surface) ***
    !                 *** MUST be ISOBARIC levels (constant pressure) ***
    ! TFULL(NLAT,NLON,NLEV,NTIME) : Temperature array [K] on isobaric levels
    ! TMSG          : Missing value for invalid data (typically -999.0)
    ! LAPSEC        : WMO lapse rate criterion [K/km] (typically 2.0)
    ! PUNIT         : Pressure unit flag (0=hPa, 1=Pa)
    !
    ! Output Parameters:
    ! =================
    ! PTROP_HPA(NLAT,NLON,NTIME)  : Tropopause pressure [hPa]
    ! HTROP_M(NLAT,NLON,NTIME)    : Tropopause height [m above sea level]
    ! ITROP(NLAT,NLON,NTIME)      : Tropopause level index (0-based, -999 if not found)
    ! LAPSE_RATE(NLAT,NLON,NTIME) : Lapse rate at tropopause [K/km]
    ! SUCCESS(NLAT,NLON,NTIME)    : .TRUE. if tropopause found, .FALSE. otherwise
    !
    implicit none

    integer, intent(in) :: nlat, nlon, nlev, ntime, nlevm, punit
    real(c_double), intent(in) :: pfull(nlev)
    real(c_double), intent(in) :: tfull(nlat, nlon, nlev, ntime)
    real(c_double), intent(in) :: tmsg, lapsec
    real(c_double), intent(out) :: ptrop_hpa(nlat, nlon, ntime)
    real(c_double), intent(out) :: htrop_m(nlat, nlon, ntime)
    integer(c_int), intent(out) :: itrop(nlat, nlon, ntime)
    real(c_double), intent(out) :: lapse_rate(nlat, nlon, ntime)
    integer(c_int), intent(out) :: success(nlat, nlon, ntime)

    integer :: t

    do t = 1, ntime
        call tropopause_grid_3d( &
            nlat, nlon, nlev, nlevm, pfull, tfull(:, :, :, t), tmsg, lapsec, punit, &
            ptrop_hpa(:, :, t), htrop_m(:, :, t), itrop(:, :, t), lapse_rate(:, :, t), success(:, :, t) &
        )
    end do
end subroutine tropopause_grid_4d

subroutine tropopause_grid_2d( &
    nspatial, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
)
    !
    ! Calculate WMO tropopause properties for 2D atmospheric cross-sections
    !
    ! *** CRITICAL: This routine ONLY accepts ISOBARIC (constant pressure level) data ***
    ! If you have model level data, you MUST interpolate to pressure levels first!
    !
    ! This routine handles 2D cross-sections like (level, lat) or (level, lon)
    !
    ! Input Parameters:
    ! ================
    ! NSPATIAL      : Number of spatial points (latitude or longitude)
    ! NLEV          : Number of vertical pressure levels
    ! NLEVM         : NLEV + 1 (for work arrays)
    ! PFULL(NLEV)   : 1D array of pressure levels [hPa or Pa]
    !                 *** MUST be in ASCENDING order (high altitude → surface) ***
    !                 *** MUST be ISOBARIC levels (constant pressure) ***
    ! TFULL(NSPATIAL,NLEV) : Temperature array [K] on isobaric levels
    ! TMSG          : Missing value for invalid data (typically -999.0)
    ! LAPSEC        : WMO lapse rate criterion [K/km] (typically 2.0)
    ! PUNIT         : Pressure unit flag (0=hPa, 1=Pa)
    !
    ! Output Parameters:
    ! =================
    ! PTROP_HPA(NSPATIAL)  : Tropopause pressure [hPa]
    ! HTROP_M(NSPATIAL)    : Tropopause height [m above sea level]
    ! ITROP(NSPATIAL)      : Tropopause level index (0-based, -999 if not found)
    ! LAPSE_RATE(NSPATIAL) : Lapse rate at tropopause [K/km]
    ! SUCCESS(NSPATIAL)    : .TRUE. if tropopause found, .FALSE. otherwise
    !
    implicit none

    integer, intent(in) :: nspatial, nlev, nlevm, punit
    real(c_double), intent(in) :: pfull(nlev)
    real(c_double), intent(in) :: tfull(nspatial, nlev)
    real(c_double), intent(in) :: tmsg, lapsec
    real(c_double), intent(out) :: ptrop_hpa(nspatial)
    real(c_double), intent(out) :: htrop_m(nspatial)
    integer(c_int), intent(out) :: itrop(nspatial)
    real(c_double), intent(out) :: lapse_rate(nspatial)
    integer(c_int), intent(out) :: success(nspatial)

    integer :: i, k, itrop_temp
    real(c_double) :: profile_t(nlev)
    real(c_double) :: work_lapse(nlevm), work_phalf(nlevm), work_pmb(nlev)
    real(c_double) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ptrop_hpa = tmsg
    htrop_m = tmsg
    itrop = -999_c_int
    lapse_rate = tmsg
    success = 0_c_int

    !$OMP PARALLEL DO PRIVATE(i, k, itrop_temp, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(nspatial, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    !$OMP ptrop_hpa, htrop_m, itrop, lapse_rate, success)
    do i = 1, nspatial
        valid_profile = .true.
        do k = 1, nlev
            profile_t(k) = tfull(i, k)
            if (abs(profile_t(k) - tmsg) < 1.0d-10) then
                valid_profile = .false.
                exit
            end if
        end do

        if (valid_profile) then
            do k = 1, nlev
                if (abs(pfull(k) - tmsg) < 1.0d-10) then
                    valid_profile = .false.
                    exit
                end if
            end do
        end if

        if (valid_profile) then
            itrop_temp = -999
            call STATTROPX( &
                nlev, nlevm, pfull, profile_t, tmsg, lapsec, punit, &
                ptrop_temp, itrop_temp, work_lapse, work_phalf, work_pmb &
            )

            if (abs(ptrop_temp - tmsg) > 1.0d-10) then
                if (punit == 0) then
                    ptrop_hpa(i) = ptrop_temp
                else
                    ptrop_hpa(i) = ptrop_temp * 0.01d0
                end if

                temp_p(1) = ptrop_hpa(i)
                call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
                htrop_m(i) = temp_h(1)
                itrop(i) = int(itrop_temp, c_int)

                if (itrop_temp >= 1 .and. itrop_temp < nlev - 1) then
                    lapse_rate(i) = work_lapse(itrop_temp + 1)
                end if

                success(i) = 1_c_int
            end if
        end if
    end do
    !$OMP END PARALLEL DO
end subroutine tropopause_grid_2d

subroutine tropopause_profile_1d( &
    nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
)
    !
    ! Calculate WMO tropopause properties for a single 1D vertical profile.
    !
    ! Input Parameters:
    ! ================
    ! NLEV          : Number of vertical pressure levels
    ! NLEVM         : NLEV + 1 (for work arrays)
    ! PFULL(NLEV)   : 1D array of pressure levels [hPa or Pa]
    ! TFULL(NLEV)    : Temperature profile [K]
    ! TMSG          : Missing value for invalid data (typically -999.0)
    ! LAPSEC        : WMO lapse rate criterion [K/km] (typically 2.0)
    ! PUNIT         : Pressure unit flag (0=hPa, 1=Pa)
    !
    ! Output Parameters:
    ! =================
    ! PTROP_HPA     : Tropopause pressure [hPa]
    ! HTROP_M       : Tropopause height [m]
    ! ITROP         : Tropopause level index
    ! LAPSE_RATE    : Tropopause lapse rate [K/km]
    ! SUCCESS       : Success flag
    !
    ! Algorithm:
    ! ==========
    ! 1. Check the profile for missing values
    ! 2. Call STATTROPX to compute tropopause pressure and level
    ! 3. Convert pressure to hPa if needed
    ! 4. Convert pressure to height using DSTDATMP
    ! 5. Extract the lapse rate at the detected tropopause level
    !
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: nlev, nlevm, punit

    ! INPUT ARRAYS - vertical profile
    real(c_double), intent(in) :: pfull(nlev)      ! Pressure [hPa or Pa]
    real(c_double), intent(in) :: tfull(nlev)      ! Temperature [K]
    real(c_double), intent(in) :: tmsg, lapsec     ! Missing value, lapse criterion

    ! OUTPUT SCALARS
    real(c_double), intent(out) :: ptrop_hpa       ! Tropopause pressure [hPa]
    real(c_double), intent(out) :: htrop_m         ! Tropopause height [m]
    integer(c_int), intent(out) :: itrop           ! Tropopause level index
    real(c_double), intent(out) :: lapse_rate      ! Tropopause lapse rate [K/km]
    integer(c_int), intent(out) :: success         ! Success flag

    ! LOCAL VARIABLES
    integer :: k
    integer :: itrop_temp
    real(c_double) :: work_lapse(nlevm), work_phalf(nlevm), work_pmb(nlev)
    real(c_double) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize outputs
    ptrop_hpa = tmsg
    htrop_m = tmsg
    itrop = -999_c_int
    lapse_rate = tmsg
    success = 0_c_int

    ! Check for missing values in profile
    valid_profile = .true.
    do k = 1, nlev
        if (abs(pfull(k) - tmsg) < 1e-10_dp .or. &
            abs(tfull(k) - tmsg) < 1e-10_dp) then
            valid_profile = .false.
            exit
        end if
    end do

    ! Process valid profile only
    if (valid_profile) then
        ! Calculate tropopause pressure and level
        call STATTROPX(nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
                      ptrop_temp, itrop_temp, work_lapse, work_phalf, work_pmb)

        ! Check if tropopause was found
        if (abs(ptrop_temp - tmsg) > 1e-10_dp) then
            ! Convert pressure units if needed
            if (punit == 0) then
                ptrop_hpa = ptrop_temp  ! Already in hPa
            else
                ptrop_hpa = ptrop_temp * 0.01_dp  ! Pa to hPa
            end if

            ! Calculate height using US Standard Atmosphere 1976
            temp_p(1) = ptrop_hpa
            call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
            htrop_m = temp_h(1)
            itrop = int(itrop_temp, c_int)

            ! Calculate lapse rate at tropopause level
            if (itrop_temp >= 1 .and. itrop_temp < nlev-1) then
                lapse_rate = work_lapse(itrop_temp+1)  ! Convert from 0-based to 1-based
            end if

            success = 1_c_int
        end if
    end if

end subroutine tropopause_profile_1d

end module tropopause_height_mod
