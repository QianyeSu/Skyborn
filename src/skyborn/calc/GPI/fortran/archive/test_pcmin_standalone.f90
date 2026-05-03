program test_pcmin_standalone
    ! Standalone test - compile separately with each implementation
    implicit none

    integer, parameter :: NL = 31  ! Number of pressure levels
    integer, parameter :: NTEST = 5  ! Number of test cases

    ! Arrays for atmospheric profiles
    real(4), dimension(NL) :: P_LEVELS, T_PROFILE, R_PROFILE
    real(4), dimension(NL) :: T_WORK, R_WORK  ! Working arrays

    ! Test case parameters
    real(4), dimension(NTEST) :: SST_TESTS, SLP_TESTS
    character(len=50), dimension(NTEST) :: TEST_NAMES

    ! Output variables
    real(4) :: PMIN, VMAX
    integer :: IFL

    integer :: i, k
    real(4) :: t_surface, q_surface, lapse_rate

    write(*,*) '========================================================================='
    write(*,*) '           PCMIN STANDALONE TEST'
    write(*,*) '========================================================================='
    write(*,*) 'Testing implementation with 5 different conditions'
    write(*,*) '========================================================================='
    write(*,*)

    ! Initialize pressure levels (hPa)
    P_LEVELS = [1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
                 750., 725., 700., 650., 600., 550., 500., 450., 400., 350., &
                 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30.,  20., 10.]

    ! Define test cases
    TEST_NAMES(1) = "Weak tropical (SST=26°C)"
    SST_TESTS(1) = 26.0
    SLP_TESTS(1) = 1010.0

    TEST_NAMES(2) = "Moderate tropical (SST=28°C)"
    SST_TESTS(2) = 28.0
    SLP_TESTS(2) = 1008.0

    TEST_NAMES(3) = "Strong tropical (SST=30°C)"
    SST_TESTS(3) = 30.0
    SLP_TESTS(3) = 1005.0

    TEST_NAMES(4) = "Extreme (SST=32°C)"
    SST_TESTS(4) = 32.0
    SLP_TESTS(4) = 1003.0

    TEST_NAMES(5) = "Marginal (SST=25°C)"
    SST_TESTS(5) = 25.0
    SLP_TESTS(5) = 1012.0

    ! Header
    write(*,'(A)') 'TEST CASE                      | SST(°C) | SLP(hPa) | PMIN(hPa) | VMAX(m/s) | Status'
    write(*,'(A)') '-------------------------------|---------|----------|-----------|-----------|-------'

    ! Run all test cases
    do i = 1, NTEST
        ! Generate realistic atmospheric profiles based on SST
        t_surface = SST_TESTS(i) - 1.0
        q_surface = 4.5 * exp(0.0687 * SST_TESTS(i))
        q_surface = min(q_surface, 20.0)

        ! Temperature profile
        do k = 1, NL
            if (P_LEVELS(k) >= 850.0) then
                lapse_rate = 6.0
                T_PROFILE(k) = t_surface - lapse_rate * (1000.0 - P_LEVELS(k)) / 100.0
            else if (P_LEVELS(k) >= 200.0) then
                lapse_rate = 6.5
                T_PROFILE(k) = t_surface - 10.0 - lapse_rate * (850.0 - P_LEVELS(k)) / 70.0
            else
                T_PROFILE(k) = -55.0 - (200.0 - P_LEVELS(k)) / 10.0
            end if
        end do

        ! Mixing ratio profile
        do k = 1, NL
            if (P_LEVELS(k) >= 700.0) then
                R_PROFILE(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 1.5
            else
                R_PROFILE(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 3.0
            end if
            R_PROFILE(k) = max(R_PROFILE(k), 0.001)
        end do

        ! Make working copies
        T_WORK = T_PROFILE
        R_WORK = R_PROFILE

#ifdef USE_ORIGINAL
        ! Call original PCMIN
        call PCMIN(SST_TESTS(i), SLP_TESTS(i), P_LEVELS, T_WORK, R_WORK, &
                   NL, NL, PMIN, VMAX, IFL)
#else
        ! Call modernized version (fixed)
        call calculate_pi_single_profile_fixed( &
            SST_TESTS(i) + 273.15,      & ! SST in K
            SLP_TESTS(i) * 100.0,        & ! SLP in Pa
            P_LEVELS,                    & ! Pressure in mb
            T_PROFILE + 273.15,          & ! Temperature in K
            R_PROFILE * 0.001,           & ! Mixing ratio in kg/kg
            NL, NL, PMIN, VMAX, IFL)
#endif

        ! Print results
        if (IFL == 1) then
            write(*,'(A30,A1,F8.2,A1,F9.2,A1,F10.2,A1,F10.2,A1,A7)') &
                TEST_NAMES(i)(1:30), '|', SST_TESTS(i), '|', SLP_TESTS(i), '|', &
                PMIN, '|', VMAX, '|', '  OK   '
        else
            write(*,'(A30,A1,F8.2,A1,F9.2,A1,A11,A1,A11,A1,A7)') &
                TEST_NAMES(i)(1:30), '|', SST_TESTS(i), '|', SLP_TESTS(i), '|', &
                '   FAILED  ', '|', '   FAILED  ', '|', 'FAILED '
        end if
    end do

    write(*,'(A)') '========================================================================='
    write(*,*) 'Test completed!'
    write(*,*) '========================================================================='

end program test_pcmin_standalone
