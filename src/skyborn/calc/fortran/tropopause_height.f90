! Multi-dimensional tropopause calculation for gridded atmospheric data
! Global subroutines for f2py compatibility
! Processes (lat, lon, level, [time]) arrays to compute tropopause properties
!
! This module contains high-performance Fortran implementations for calculating
! WMO tropopause properties from atmospheric data using the thermal tropopause
! definition. It includes optimized routines for 1D profiles, 2D cross-sections,
! 3D spatial data, and 4D time series analysis.
!
! The tropopause calculation follows the WMO (1957) definition:
! "The first tropopause is defined as the lowest level at which
! the lapse rate decreases to 2 deg C per kilometer or less,
! provided also the average lapse rate between this level and
! all higher levels within 2 kilometers does not exceed 2 deg C"
!
! References:
! - WMO (1992): International meteorological vocabulary, Genf, 784pp.
! - Randel WJ, Wu F, Gaffen DJ, JGR, 105, 15509-15523, 2000.
! - US Standard Atmosphere, 1976 (for height calculations)
! - Original tropopause algorithm: Dominik Brunner, V1.0 Aug 2000
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
    ! This routine processes spatial atmospheric data (lat, lon, level) to identify
    ! the thermal tropopause at each grid point using the WMO definition.
    !
    ! Algorithm:
    ! 1. For each grid point, calculate lapse rate gamma = -dT/dz (K/km)
    !    Using hydrostatic approximation: gamma = CONST * d(ln T)/d(ln p)
    !    where CONST = 1000 * g / R_gas = 1000 * 9.80665 / 287.04
    !
    ! 2. Search from surface to top for lowest level where:
    !    - Lapse rate < threshold (default 2.0 K/km)
    !    - Pressure < 450 hPa (exclude boundary layer effects)
    !    - Average lapse rate in 2 km layer above also < threshold
    !
    ! 3. Interpolate tropopause pressure linearly in log(p) space
    !
    ! 4. Calculate height using US Standard Atmosphere 1976
    !
    ! 5. Enforce minimum tropopause pressure of 85 hPa
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
    real(dp) :: profile_p(NLEV), profile_t(NLEV)
    real(dp) :: work_lapse(NLEVM), work_phalf(NLEVM), work_pmb(NLEV)
    real(dp) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize output arrays
    PTROP_HPA = TMSG
    HTROP_M = TMSG
    ITROP = -999
    LAPSE_RATE = TMSG
    SUCCESS = .false.

    !$OMP PARALLEL DO PRIVATE(i, j, k, profile_p, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(NLAT, NLON, NLEV, NLEVM, PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
    !$OMP PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS) COLLAPSE(2)
    do j = 1, NLON
        do i = 1, NLAT
            ! Extract vertical profile for this grid point
            valid_profile = .true.
            do k = 1, NLEV
                profile_p(k) = PFULL(k)         ! 1D pressure array
                profile_t(k) = TFULL(i, j, k)

                ! Check for missing values
                if (abs(profile_p(k) - TMSG) < 1e-10_dp .or. &
                    abs(profile_t(k) - TMSG) < 1e-10_dp) then
                    valid_profile = .false.
                    exit
                end if
            end do

            ! Process valid profiles only
            if (valid_profile) then
                ! Calculate tropopause pressure and level
                call STATTROPX(NLEV, NLEVM, profile_p, profile_t, TMSG, LAPSEC, PUNIT, &
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
    real(dp) :: profile_p(NLEV), profile_t(NLEV)
    real(dp) :: work_lapse(NLEVM), work_phalf(NLEVM), work_pmb(NLEV)
    real(dp) :: ptrop_temp, temp_p(1), temp_t(1), temp_d(1), temp_h(1)
    logical :: valid_profile

    ! Initialize output arrays
    PTROP_HPA = TMSG
    HTROP_M = TMSG
    ITROP = -999
    LAPSE_RATE = TMSG
    SUCCESS = .false.

    !$OMP PARALLEL DO PRIVATE(i, k, profile_p, profile_t, work_lapse, &
    !$OMP work_phalf, work_pmb, ptrop_temp, temp_p, temp_t, temp_d, temp_h, valid_profile) &
    !$OMP SHARED(NSPATIAL, NLEV, NLEVM, PFULL, TFULL, TMSG, LAPSEC, PUNIT, &
    !$OMP PTROP_HPA, HTROP_M, ITROP, LAPSE_RATE, SUCCESS)
    do i = 1, NSPATIAL
        ! Extract vertical profile for this spatial point
        valid_profile = .true.
        do k = 1, NLEV
            profile_p(k) = PFULL(k)         ! 1D pressure array
            profile_t(k) = TFULL(i, k)

            ! Check for missing values
            if (abs(profile_p(k) - TMSG) < 1e-10_dp .or. &
                abs(profile_t(k) - TMSG) < 1e-10_dp) then
                valid_profile = .false.
                exit
            end if
        end do

        ! Process valid profiles only
        if (valid_profile) then
            ! Calculate tropopause pressure and level
            call STATTROPX(NLEV, NLEVM, profile_p, profile_t, TMSG, LAPSEC, PUNIT, &
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
