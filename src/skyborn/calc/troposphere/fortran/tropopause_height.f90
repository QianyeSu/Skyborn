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

subroutine TROPOPAUSE_GRID_3D(NLAT, NLON, NLEV, NLEVM, &
                             PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
                             PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
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
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: NLAT, NLON, NLEV, NLEVM, PUNIT

    ! INPUT ARRAYS
    real(dp), intent(in) :: PFULL(NLEV)                ! 1D pressure levels [hPa or Pa]
    real(dp), intent(in) :: TFULL(NLAT, NLON, NLEV)    ! Temperature [K]
    real(dp), intent(in) :: TMSG, LAPSEC               ! Missing value, lapse criterion

    ! OUTPUT ARRAYS
    real(dp), intent(out) :: PTROP_HPA(NLAT, NLON)     ! Tropopause pressure [hPa]
    real(dp), intent(out) :: HTROP_M(NLAT, NLON)       ! Tropopause height [m]
    integer, intent(out) :: ITROP(NLAT, NLON)          ! Tropopause level index
    real(dp), intent(out) :: LAPSE_RATE(NLAT, NLON)    ! Tropopause lapse rate [K/km]
    logical, intent(out) :: SUCCESS(NLAT, NLON)        ! Success flag

    ! LOCAL VARIABLES
    integer :: i, j, k
    real(dp) :: profile_t(NLEV)
    real(dp) :: work_lapse(NLEVM), work_phalf(NLEVM), work_pmb(NLEV)
    real(dp) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize output arrays
    PTROP_HPA = TMSG
    HTROP_M = TMSG
    ITROP = -999
    LAPSE_RATE = TMSG
    SUCCESS = .false.

    !$OMP PARALLEL DO PRIVATE(i, j, k, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(NLAT, NLON, NLEV, NLEVM, PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
    !$OMP PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS) COLLAPSE(2)
    do j = 1, NLON
        do i = 1, NLAT
            ! Extract vertical temperature profile for this grid point
            ! No need to copy PFULL since it's the same 1D array for all grid points (isobaric data)
            valid_profile = .true.
            do k = 1, NLEV
                profile_t(k) = TFULL(i, j, k)

                ! Check for missing values (pressure check done once outside loop for efficiency)
                if (abs(profile_t(k) - TMSG) < 1e-10_dp) then
                    valid_profile = .false.
                    exit
                end if
            end do

            ! Check pressure levels for missing values (done once per grid point)
            if (valid_profile) then
                do k = 1, NLEV
                    if (abs(PFULL(k) - TMSG) < 1e-10_dp) then
                        valid_profile = .false.
                        exit
                    end if
                end do
            end if

            ! Process valid profiles only
            if (valid_profile) then
                ! Calculate tropopause pressure and level
                ! Pass PFULL directly (1D isobaric levels, same for all grid points)
                call STATTROPX(NLEV, NLEVM, PFULL, profile_t, TMSG, LAPSEC, PUNIT, &
                              ptrop_temp, ITROP(i,j), work_lapse, work_phalf, work_pmb)

                ! Check if tropopause was found
                if (abs(ptrop_temp - TMSG) > 1e-10_dp) then
                    ! Convert pressure units if needed
                    if (PUNIT == 0) then
                        PTROP_HPA(i,j) = ptrop_temp  ! Already in hPa
                    else
                        PTROP_HPA(i,j) = ptrop_temp * 0.01_dp  ! Pa to hPa
                    end if

                    ! Calculate height using US Standard Atmosphere 1976
                    ! DSTDATMP expects array inputs, so we use temporary variables
                    temp_p(1) = PTROP_HPA(i,j)
                    call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
                    HTROP_M(i,j) = temp_h(1)

                    ! Calculate lapse rate at tropopause level
                    if (ITROP(i,j) >= 1 .and. ITROP(i,j) < NLEV-1) then
                        LAPSE_RATE(i,j) = work_lapse(ITROP(i,j)+1)  ! Convert from 0-based to 1-based
                    end if

                    SUCCESS(i,j) = .true.
                end if
            end if
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine TROPOPAUSE_GRID_3D

subroutine TROPOPAUSE_GRID_4D(NLAT, NLON, NLEV, NTIME, NLEVM, &
                             PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
                             PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
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

    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: NLAT, NLON, NLEV, NTIME, NLEVM, PUNIT

    ! INPUT ARRAYS
    real(dp), intent(in) :: PFULL(NLEV)                     ! 1D pressure levels [hPa or Pa]
    real(dp), intent(in) :: TFULL(NLAT, NLON, NLEV, NTIME)  ! Temperature
    real(dp), intent(in) :: TMSG, LAPSEC

    ! OUTPUT ARRAYS
    real(dp), intent(out) :: PTROP_HPA(NLAT, NLON, NTIME)
    real(dp), intent(out) :: HTROP_M(NLAT, NLON, NTIME)
    integer, intent(out) :: ITROP(NLAT, NLON, NTIME)
    real(dp), intent(out) :: LAPSE_RATE(NLAT, NLON, NTIME)
    logical, intent(out) :: SUCCESS(NLAT, NLON, NTIME)

    ! LOCAL VARIABLES
    integer :: t

    ! Process each time step
    do t = 1, NTIME
        call TROPOPAUSE_GRID_3D(NLAT, NLON, NLEV, NLEVM, &
                               PFULL, TFULL(:,:,:,t), &
                               TMSG, LAPSEC, PUNIT, &
                               PTROP_HPA(:,:,t), HTROP_M(:,:,t), ITROP(:,:,t), &
                               LAPSE_RATE(:,:,t), SUCCESS(:,:,t))
    end do

end subroutine TROPOPAUSE_GRID_4D

subroutine TROPOPAUSE_GRID_2D(NSPATIAL, NLEV, NLEVM, &
                             PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
                             PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
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

    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: NSPATIAL, NLEV, NLEVM, PUNIT

    ! INPUT ARRAYS
    real(dp), intent(in) :: PFULL(NLEV)              ! 1D pressure levels [hPa or Pa]
    real(dp), intent(in) :: TFULL(NSPATIAL, NLEV)    ! Temperature [K]
    real(dp), intent(in) :: TMSG, LAPSEC             ! Missing value, lapse criterion

    ! OUTPUT ARRAYS
    real(dp), intent(out) :: PTROP_HPA(NSPATIAL)     ! Tropopause pressure [hPa]
    real(dp), intent(out) :: HTROP_M(NSPATIAL)       ! Tropopause height [m]
    integer, intent(out) :: ITROP(NSPATIAL)          ! Tropopause level index
    real(dp), intent(out) :: LAPSE_RATE(NSPATIAL)    ! Tropopause lapse rate [K/km]
    logical, intent(out) :: SUCCESS(NSPATIAL)        ! Success flag

    ! LOCAL VARIABLES
    integer :: i, k
    real(dp) :: profile_t(NLEV)
    real(dp) :: work_lapse(NLEVM), work_phalf(NLEVM), work_pmb(NLEV)
    real(dp) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize output arrays
    PTROP_HPA = TMSG
    HTROP_M = TMSG
    ITROP = -999
    LAPSE_RATE = TMSG
    SUCCESS = .false.

    !$OMP PARALLEL DO PRIVATE(i, k, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(NSPATIAL, NLEV, NLEVM, PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
    !$OMP PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
    do i = 1, NSPATIAL
        ! Extract vertical temperature profile for this spatial point
        ! No need to copy PFULL since it's the same 1D array for all spatial points (isobaric data)
        valid_profile = .true.
        do k = 1, NLEV
            profile_t(k) = TFULL(i, k)

            ! Check for missing values (pressure check done once outside loop for efficiency)
            if (abs(profile_t(k) - TMSG) < 1e-10_dp) then
                valid_profile = .false.
                exit
            end if
        end do

        ! Check pressure levels for missing values (done once per spatial point)
        if (valid_profile) then
            do k = 1, NLEV
                if (abs(PFULL(k) - TMSG) < 1e-10_dp) then
                    valid_profile = .false.
                    exit
                end if
            end do
        end if

        ! Process valid profiles only
        if (valid_profile) then
            ! Calculate tropopause pressure and level
            ! Pass PFULL directly (1D isobaric levels, same for all spatial points)
            call STATTROPX(NLEV, NLEVM, PFULL, profile_t, TMSG, LAPSEC, PUNIT, &
                          ptrop_temp, ITROP(i), work_lapse, work_phalf, work_pmb)

            ! Check if tropopause was found
            if (abs(ptrop_temp - TMSG) > 1e-10_dp) then
                ! Convert pressure units if needed
                if (PUNIT == 0) then
                    PTROP_HPA(i) = ptrop_temp  ! Already in hPa
                else
                    PTROP_HPA(i) = ptrop_temp * 0.01_dp  ! Pa to hPa
                end if

                ! Calculate height using US Standard Atmosphere 1976
                temp_p(1) = PTROP_HPA(i)
                call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
                HTROP_M(i) = temp_h(1)

                ! Calculate lapse rate at tropopause level
                if (ITROP(i) >= 1 .and. ITROP(i) < NLEV-1) then
                    LAPSE_RATE(i) = work_lapse(ITROP(i)+1)  ! Convert from 0-based to 1-based
                end if

                SUCCESS(i) = .true.
            end if
        end if
    end do
    !$OMP END PARALLEL DO

end subroutine TROPOPAUSE_GRID_2D

subroutine TROPOPAUSE_PROFILE_1D(NLEV, NLEVM, &
                                PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
                                PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)

    ! INPUT DIMENSIONS
    integer, intent(in) :: NLEV, NLEVM, PUNIT

    ! INPUT ARRAYS - vertical profile
    real(dp), intent(in) :: PFULL(NLEV)      ! Pressure [hPa or Pa]
    real(dp), intent(in) :: TFULL(NLEV)      ! Temperature [K]
    real(dp), intent(in) :: TMSG, LAPSEC     ! Missing value, lapse criterion

    ! OUTPUT SCALARS
    real(dp), intent(out) :: PTROP_HPA       ! Tropopause pressure [hPa]
    real(dp), intent(out) :: HTROP_M         ! Tropopause height [m]
    integer, intent(out) :: ITROP            ! Tropopause level index
    real(dp), intent(out) :: LAPSE_RATE      ! Tropopause lapse rate [K/km]
    logical, intent(out) :: SUCCESS          ! Success flag

    ! LOCAL VARIABLES
    integer :: k
    real(dp) :: work_lapse(NLEVM), work_phalf(NLEVM), work_pmb(NLEV)
    real(dp) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize outputs
    PTROP_HPA = TMSG
    HTROP_M = TMSG
    ITROP = -999
    LAPSE_RATE = TMSG
    SUCCESS = .false.

    ! Check for missing values in profile
    valid_profile = .true.
    do k = 1, NLEV
        if (abs(PFULL(k) - TMSG) < 1e-10_dp .or. &
            abs(TFULL(k) - TMSG) < 1e-10_dp) then
            valid_profile = .false.
            exit
        end if
    end do

    ! Process valid profile only
    if (valid_profile) then
        ! Calculate tropopause pressure and level
        call STATTROPX(NLEV, NLEVM, PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
                      ptrop_temp, ITROP, work_lapse, work_phalf, work_pmb)

        ! Check if tropopause was found
        if (abs(ptrop_temp - TMSG) > 1e-10_dp) then
            ! Convert pressure units if needed
            if (PUNIT == 0) then
                PTROP_HPA = ptrop_temp  ! Already in hPa
            else
                PTROP_HPA = ptrop_temp * 0.01_dp  ! Pa to hPa
            end if

            ! Calculate height using US Standard Atmosphere 1976
            temp_p(1) = PTROP_HPA
            call DSTDATMP(1, temp_p, temp_t, temp_d, temp_h)
            HTROP_M = temp_h(1)

            ! Calculate lapse rate at tropopause level
            if (ITROP >= 1 .and. ITROP < NLEV-1) then
                LAPSE_RATE = work_lapse(ITROP+1)  ! Convert from 0-based to 1-based
            end if

            SUCCESS = .true.
        end if
    end if

end subroutine TROPOPAUSE_PROFILE_1D
