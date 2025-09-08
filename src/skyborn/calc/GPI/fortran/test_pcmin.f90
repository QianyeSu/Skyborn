program test_pcmin_comparison
    implicit none

    ! Comprehensive comparison test for PCMIN implementations:
    ! 1. Original Emanuel pcmin_2013.f (PCMIN subroutine)
    ! 2. Modernized tropical_cyclone_potential_intensity.f90 (calculate_pi_single_profile)

    integer, parameter :: NL = 31  ! Number of pressure levels
    integer, parameter :: NTEST = 5  ! Number of test cases

    ! Arrays for atmospheric profiles
    real(4), dimension(NL) :: P_LEVELS, T_PROFILE, R_PROFILE
    real(4), dimension(NL) :: T_WORK, R_WORK  ! Working arrays for PCMIN (gets modified)

    ! Test case parameters
    real(4), dimension(NTEST) :: SST_TESTS, SLP_TESTS
    character(len=50), dimension(NTEST) :: TEST_NAMES

    ! Output variables for original PCMIN
    real(4) :: PMIN_ORIG, VMAX_ORIG
    integer :: IFL_ORIG

    ! Output variables for modernized calculate_pi_single_profile
    real(4) :: PMIN_NEW, VMAX_NEW
    integer :: IFL_NEW

    ! Statistics variables
    real(4) :: pmin_diff, vmax_diff, pmin_percent, vmax_percent
    real(4) :: max_pmin_diff, max_vmax_diff
    real(4) :: avg_pmin_diff, avg_vmax_diff
    integer :: i, j, nvalid

    write(*,*) '========================================================================='
    write(*,*) '           PCMIN IMPLEMENTATION COMPARISON TEST'
    write(*,*) '========================================================================='
    write(*,*) 'Comparing:'
    write(*,*) '  1. Original Emanuel pcmin_2013.f (PCMIN)'
    write(*,*) '  2. Modernized tropical_cyclone_potential_intensity.f90'
    write(*,*) '========================================================================='
    write(*,*)

    ! Initialize pressure levels (hPa)
    P_LEVELS = [1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
                 750., 725., 700., 650., 600., 550., 500., 450., 400., 350., &
                 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30.,  20., 10.]

    ! Define test cases
    TEST_NAMES(1) = "Weak tropical conditions (SST=26°C)"
    SST_TESTS(1) = 26.0
    SLP_TESTS(1) = 1010.0

    TEST_NAMES(2) = "Moderate tropical conditions (SST=28°C)"
    SST_TESTS(2) = 28.0
    SLP_TESTS(2) = 1008.0

    TEST_NAMES(3) = "Strong tropical conditions (SST=30°C)"
    SST_TESTS(3) = 30.0
    SLP_TESTS(3) = 1005.0

    TEST_NAMES(4) = "Extreme conditions (SST=32°C)"
    SST_TESTS(4) = 32.0
    SLP_TESTS(4) = 1003.0

    TEST_NAMES(5) = "Marginal conditions (SST=25°C)"
    SST_TESTS(5) = 25.0
    SLP_TESTS(5) = 1012.0

    ! Initialize statistics
    max_pmin_diff = 0.0
    max_vmax_diff = 0.0
    avg_pmin_diff = 0.0
    avg_vmax_diff = 0.0
    nvalid = 0

    ! Header for results table
    write(*,'(A)') 'TEST CASE                                  | SST(°C) | SLP(hPa) |   PCMIN RESULTS   |  MODERN RESULTS   | DIFFERENCES'
    write(*,'(A)') '                                           |         |          | PMIN    VMAX      | PMIN    VMAX      | ΔPMIN   ΔVMAX'
    write(*,'(A)') '-------------------------------------------|---------|----------|-------------------|-------------------|-------------'

    ! Run all test cases
    do i = 1, NTEST
        ! Generate realistic temperature and mixing ratio profiles based on SST
        call generate_atmospheric_profile(SST_TESTS(i), T_PROFILE, R_PROFILE)

        ! Make working copies for PCMIN (it modifies the arrays)
        T_WORK = T_PROFILE
        R_WORK = R_PROFILE

        ! Call original PCMIN
        call PCMIN(SST_TESTS(i), SLP_TESTS(i), P_LEVELS, T_WORK, R_WORK, &
                   NL, NL, PMIN_ORIG, VMAX_ORIG, IFL_ORIG)

        ! Call modernized calculate_pi_single_profile (fixed version)
        ! Note: Need to convert units properly
        call calculate_pi_single_profile_fixed( &
            SST_TESTS(i) + 273.15,      & ! SST in K
            SLP_TESTS(i) * 100.0,        & ! SLP in Pa
            P_LEVELS,                    & ! Pressure in mb (same as PCMIN)
            T_PROFILE + 273.15,          & ! Temperature in K
            R_PROFILE * 0.001,           & ! Mixing ratio in kg/kg
            NL, NL, PMIN_NEW, VMAX_NEW, IFL_NEW)

        ! Calculate differences
        if (IFL_ORIG == 1 .and. IFL_NEW == 1) then
            pmin_diff = abs(PMIN_NEW - PMIN_ORIG)
            vmax_diff = abs(VMAX_NEW - VMAX_ORIG)

            if (PMIN_ORIG > 0) pmin_percent = 100.0 * pmin_diff / PMIN_ORIG
            if (VMAX_ORIG > 0) vmax_percent = 100.0 * vmax_diff / VMAX_ORIG

            max_pmin_diff = max(max_pmin_diff, pmin_diff)
            max_vmax_diff = max(max_vmax_diff, vmax_diff)
            avg_pmin_diff = avg_pmin_diff + pmin_diff
            avg_vmax_diff = avg_vmax_diff + vmax_diff
            nvalid = nvalid + 1

            ! Print results
            write(*,'(A43,A1,F8.2,A1,F9.2,A1,F7.2,F9.2,A1,F7.2,F9.2,A1,F7.2,F8.2)') &
                TEST_NAMES(i), '|', SST_TESTS(i), '|', SLP_TESTS(i), '|', &
                PMIN_ORIG, VMAX_ORIG, '|', PMIN_NEW, VMAX_NEW, '|', &
                pmin_diff, vmax_diff
        else
            write(*,'(A43,A1,F8.2,A1,F9.2,A1,A18,A1,A18,A1,A13)') &
                TEST_NAMES(i), '|', SST_TESTS(i), '|', SLP_TESTS(i), '|', &
                '  FAILED         ', '|', '  FAILED         ', '|', '     N/A     '
        end if
    end do

    write(*,'(A)') '========================================================================='
    write(*,*)

    ! Calculate averages
    if (nvalid > 0) then
        avg_pmin_diff = avg_pmin_diff / real(nvalid)
        avg_vmax_diff = avg_vmax_diff / real(nvalid)
    end if

    ! Print detailed statistics
    write(*,*) 'DETAILED STATISTICS:'
    write(*,*) '-------------------'
    write(*,'(A,I3,A,I3)') ' Valid comparisons: ', nvalid, ' out of ', NTEST
    write(*,*)

    if (nvalid > 0) then
        write(*,*) 'Pressure differences (PMIN):'
        write(*,'(A,F10.4,A)') '  Maximum difference: ', max_pmin_diff, ' hPa'
        write(*,'(A,F10.4,A)') '  Average difference: ', avg_pmin_diff, ' hPa'
        write(*,*)

        write(*,*) 'Wind speed differences (VMAX):'
        write(*,'(A,F10.4,A)') '  Maximum difference: ', max_vmax_diff, ' m/s'
        write(*,'(A,F10.4,A)') '  Average difference: ', avg_vmax_diff, ' m/s'
        write(*,*)

        ! Assessment
        write(*,*) 'ASSESSMENT:'
        write(*,*) '-----------'

        if (max_pmin_diff < 0.5 .and. max_vmax_diff < 0.5) then
            write(*,*) '✓ EXCELLENT: Both implementations produce nearly identical results'
            write(*,*) '  Maximum differences < 0.5 units'
        else if (max_pmin_diff < 2.0 .and. max_vmax_diff < 2.0) then
            write(*,*) '✓ GOOD: Implementations are consistent within acceptable tolerances'
            write(*,*) '  Maximum differences < 2.0 units'
        else if (max_pmin_diff < 5.0 .and. max_vmax_diff < 5.0) then
            write(*,*) '⚠ WARNING: Some differences detected between implementations'
            write(*,*) '  Maximum differences < 5.0 units'
        else
            write(*,*) '✗ ALERT: Significant differences detected!'
            write(*,*) '  Maximum differences > 5.0 units'
            write(*,*) '  Further investigation recommended'
        end if
    else
        write(*,*) '✗ ERROR: No valid comparisons could be made'
    end if

    write(*,*)
    write(*,*) '========================================================================='
    write(*,*) 'Test completed!'
    write(*,*) '========================================================================='

contains

    subroutine generate_atmospheric_profile(sst, temp, mixr)
        ! Generate realistic atmospheric profiles based on SST
        real(4), intent(in) :: sst
        real(4), intent(out) :: temp(NL), mixr(NL)
        integer :: k
        real(4) :: t_surface, q_surface, lapse_rate

        ! Surface temperature slightly cooler than SST
        t_surface = sst - 1.0

        ! Surface mixing ratio based on SST (approximate saturation)
        q_surface = 4.5 * exp(0.0687 * (sst - 273.15))  ! Simplified C-C relation
        q_surface = min(q_surface, 20.0)  ! Cap at 20 g/kg

        ! Generate temperature profile with realistic lapse rate
        do k = 1, NL
            if (P_LEVELS(k) >= 850.0) then
                ! Boundary layer - small lapse rate
                lapse_rate = 6.0  ! K/km
                temp(k) = t_surface - lapse_rate * (1000.0 - P_LEVELS(k)) / 100.0
            else if (P_LEVELS(k) >= 200.0) then
                ! Troposphere - standard lapse rate
                lapse_rate = 6.5
                temp(k) = t_surface - 10.0 - lapse_rate * (850.0 - P_LEVELS(k)) / 70.0
            else
                ! Stratosphere - isothermal or slight warming
                temp(k) = -55.0 - (200.0 - P_LEVELS(k)) / 10.0
            end if
        end do

        ! Generate mixing ratio profile (exponential decay with height)
        do k = 1, NL
            if (P_LEVELS(k) >= 700.0) then
                mixr(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 1.5
            else
                mixr(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 3.0
            end if
            mixr(k) = max(mixr(k), 0.001)  ! Minimum value
        end do

    end subroutine generate_atmospheric_profile

end program test_pcmin_comparison
