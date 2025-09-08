program debug_cape_comparison
    ! Debug program to trace CAPE calculations in both implementations
    implicit none

    integer, parameter :: NL = 31
    real(4), dimension(NL) :: P_LEVELS, T_PROFILE, R_PROFILE
    real(4), dimension(NL) :: T_WORK, R_WORK

    ! Test with a single case
    real(4) :: SST_TEST = 28.0  ! °C
    real(4) :: SLP_TEST = 1008.0  ! hPa

    ! Output variables for PCMIN
    real(4) :: PMIN_ORIG, VMAX_ORIG
    integer :: IFL_ORIG

    ! Output variables for modernized
    real(4) :: PMIN_NEW, VMAX_NEW
    integer :: IFL_NEW

    ! CAPE tracking variables
    real(4) :: CAPE_ENV_ORIG, CAPE_ENV_NEW
    real(4) :: CAPE_RM_ORIG, CAPE_RM_NEW
    real(4) :: CAPE_SAT_ORIG, CAPE_SAT_NEW
    real(4) :: TOB_ENV_ORIG, TOB_ENV_NEW

    integer :: k
    real(4) :: t_surface, q_surface, lapse_rate

    write(*,*) '========================================================================='
    write(*,*) '                  CAPE CALCULATION DEBUG COMPARISON'
    write(*,*) '========================================================================='
    write(*,'(A,F6.2,A)') ' Test Case: SST = ', SST_TEST, ' °C'
    write(*,'(A,F7.2,A)') '           SLP = ', SLP_TEST, ' hPa'
    write(*,*) '========================================================================='
    write(*,*)

    ! Initialize pressure levels (hPa)
    P_LEVELS = [1000., 975., 950., 925., 900., 875., 850., 825., 800., 775., &
                 750., 725., 700., 650., 600., 550., 500., 450., 400., 350., &
                 300., 250., 200., 150., 100.,  70.,  50.,  40.,  30.,  20., 10.]

    ! Generate atmospheric profile
    t_surface = SST_TEST - 1.0
    q_surface = 4.5 * exp(0.0687 * SST_TEST)
    q_surface = min(q_surface, 20.0)

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

    do k = 1, NL
        if (P_LEVELS(k) >= 700.0) then
            R_PROFILE(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 1.5
        else
            R_PROFILE(k) = q_surface * (P_LEVELS(k) / 1000.0) ** 3.0
        end if
        R_PROFILE(k) = max(R_PROFILE(k), 0.001)
    end do

    write(*,*) 'ATMOSPHERIC PROFILE (first 10 levels):'
    write(*,*) 'Level   P(hPa)    T(°C)     R(g/kg)'
    write(*,*) '-----   -------   -------   --------'
    do k = 1, 10
        write(*,'(I5, 3F10.2)') k, P_LEVELS(k), T_PROFILE(k), R_PROFILE(k)
    end do
    write(*,*)

    ! Test CAPE calculations directly
    write(*,*) '========================================================================='
    write(*,*) '                        CAPE CALCULATIONS'
    write(*,*) '========================================================================='

    ! Calculate environmental CAPE with original PCMIN
    T_WORK = T_PROFILE
    R_WORK = R_PROFILE

    ! Convert units for PCMIN (already in °C and g/kg)
    call debug_pcmin_cape(SST_TEST, SLP_TEST, P_LEVELS, T_WORK, R_WORK, NL, &
                         CAPE_ENV_ORIG, CAPE_RM_ORIG, CAPE_SAT_ORIG, &
                         TOB_ENV_ORIG, PMIN_ORIG, VMAX_ORIG, IFL_ORIG)

    write(*,*) 'ORIGINAL PCMIN RESULTS:'
    write(*,*) '-----------------------'
    write(*,'(A,F10.2,A)') ' Environmental CAPE: ', CAPE_ENV_ORIG, ' J/kg'
    write(*,'(A,F10.2,A)') ' CAPE at Rmax:       ', CAPE_RM_ORIG, ' J/kg'
    write(*,'(A,F10.2,A)') ' Saturated CAPE:     ', CAPE_SAT_ORIG, ' J/kg'
    write(*,'(A,F10.2,A)') ' Outflow Temperature:', TOB_ENV_ORIG, ' K'
    write(*,'(A,F10.2,A)') ' PMIN:               ', PMIN_ORIG, ' hPa'
    write(*,'(A,F10.2,A)') ' VMAX:               ', VMAX_ORIG, ' m/s'
    write(*,*)

    ! Calculate with modernized version
    call debug_modern_cape(SST_TEST + 273.15, SLP_TEST * 100.0, P_LEVELS, &
                          T_PROFILE + 273.15, R_PROFILE * 0.001, NL, &
                          CAPE_ENV_NEW, CAPE_RM_NEW, CAPE_SAT_NEW, &
                          TOB_ENV_NEW, PMIN_NEW, VMAX_NEW, IFL_NEW)

    write(*,*) 'MODERNIZED VERSION RESULTS:'
    write(*,*) '---------------------------'
    write(*,'(A,F10.2,A)') ' Environmental CAPE: ', CAPE_ENV_NEW, ' J/kg'
    write(*,'(A,F10.2,A)') ' CAPE at Rmax:       ', CAPE_RM_NEW, ' J/kg'
    write(*,'(A,F10.2,A)') ' Saturated CAPE:     ', CAPE_SAT_NEW, ' J/kg'
    write(*,'(A,F10.2,A)') ' Outflow Temperature:', TOB_ENV_NEW, ' K'
    write(*,'(A,F10.2,A)') ' PMIN:               ', PMIN_NEW, ' hPa'
    write(*,'(A,F10.2,A)') ' VMAX:               ', VMAX_NEW, ' m/s'
    write(*,*)

    write(*,*) '========================================================================='
    write(*,*) '                           DIFFERENCES'
    write(*,*) '========================================================================='
    write(*,'(A,F10.2,A,F6.1,A)') ' ΔCAPE_env:  ', CAPE_ENV_NEW - CAPE_ENV_ORIG, &
                                   ' J/kg (', 100.0*(CAPE_ENV_NEW - CAPE_ENV_ORIG)/CAPE_ENV_ORIG, '%)'
    write(*,'(A,F10.2,A,F6.1,A)') ' ΔCAPE_rm:   ', CAPE_RM_NEW - CAPE_RM_ORIG, &
                                   ' J/kg (', 100.0*(CAPE_RM_NEW - CAPE_RM_ORIG)/CAPE_RM_ORIG, '%)'
    write(*,'(A,F10.2,A,F6.1,A)') ' ΔCAPE_sat:  ', CAPE_SAT_NEW - CAPE_SAT_ORIG, &
                                   ' J/kg (', 100.0*(CAPE_SAT_NEW - CAPE_SAT_ORIG)/CAPE_SAT_ORIG, '%)'
    write(*,'(A,F10.2,A)')         ' ΔT_outflow: ', TOB_ENV_NEW - TOB_ENV_ORIG, ' K'
    write(*,'(A,F10.2,A)')         ' ΔPMIN:      ', PMIN_NEW - PMIN_ORIG, ' hPa'
    write(*,'(A,F10.2,A)')         ' ΔVMAX:      ', VMAX_NEW - VMAX_ORIG, ' m/s'
    write(*,*)

    ! Analysis
    write(*,*) 'KEY OBSERVATIONS:'
    write(*,*) '-----------------'
    if (abs(CAPE_ENV_NEW - CAPE_ENV_ORIG) > 100.0) then
        write(*,*) '⚠ Large difference in Environmental CAPE calculation'
    end if
    if (abs(CAPE_SAT_NEW - CAPE_SAT_ORIG) > 100.0) then
        write(*,*) '⚠ Large difference in Saturated CAPE calculation'
    end if
    if (abs(TOB_ENV_NEW - TOB_ENV_ORIG) > 5.0) then
        write(*,*) '⚠ Significant difference in Outflow Temperature'
    end if

    write(*,*) '========================================================================='

contains

    subroutine debug_pcmin_cape(sst, slp, p, t, r, n, cape_env, cape_rm, cape_sat, tob, pmin, vmax, ifl)
        ! Wrapper to extract CAPE values from PCMIN calculation
        real(4), intent(in) :: sst, slp
        real(4), intent(in) :: p(n), t(n), r(n)
        integer, intent(in) :: n
        real(4), intent(out) :: cape_env, cape_rm, cape_sat, tob
        real(4), intent(out) :: pmin, vmax
        integer, intent(out) :: ifl

        real(4) :: t_work(n), r_work(n)
        real(4) :: tp, rp, pp, sstk, es0
        real(4) :: toa, tom, toms
        integer :: iflag

        ! Make working copies
        t_work = t
        r_work = r

        ! Call PCMIN to get overall results
        call PCMIN(sst, slp, p, t_work, r_work, n, n, pmin, vmax, ifl)

        ! Now calculate individual CAPE values
        ! Note: PCMIN modifies arrays, so use fresh copies
        t_work = t
        r_work = r

        ! Convert for CAPE calculation (PCMIN expects °C, converts internally)
        sstk = sst + 273.15
        es0 = 6.112 * exp(17.67 * sst / (243.5 + sst))

        ! Convert to internal units (same as PCMIN does)
        r_work = r_work * 0.001  ! g/kg to kg/kg
        t_work = t_work + 273.15  ! °C to K

        ! Environmental CAPE
        tp = t_work(1)
        rp = r_work(1)
        pp = p(1)
        call CAPE_PCMIN(tp, rp, pp, t_work, r_work, p, n, n, 0.0, cape_env, toa, iflag)
        tob = toa

        ! CAPE at radius of maximum winds (using pmin from calculation)
        tp = t_work(1)
        pp = min(pmin, 1000.0)
        rp = 0.622 * r_work(1) * slp / (pp * (0.622 + r_work(1)) - r_work(1) * slp)
        call CAPE_PCMIN(tp, rp, pp, t_work, r_work, p, n, n, 0.0, cape_rm, tom, iflag)

        ! Saturated CAPE
        tp = sstk
        pp = min(pmin, 1000.0)
        rp = 0.622 * es0 / (pp - es0)
        call CAPE_PCMIN(tp, rp, pp, t_work, r_work, p, n, n, 0.0, cape_sat, toms, iflag)

    end subroutine debug_pcmin_cape

    subroutine debug_modern_cape(sst_k, slp_pa, p, t_k, r_kgkg, n, cape_env, cape_rm, cape_sat, tob, pmin, vmax, ifl)
        ! Wrapper to extract CAPE values from modernized calculation
        real(4), intent(in) :: sst_k, slp_pa
        real(4), intent(in) :: p(n), t_k(n), r_kgkg(n)
        integer, intent(in) :: n
        real(4), intent(out) :: cape_env, cape_rm, cape_sat, tob
        real(4), intent(out) :: pmin, vmax
        integer, intent(out) :: ifl

        real(4) :: tp, rp, pp
        real(4) :: sst_c, slp_mb, es0
        real(4) :: t_c(n), r_gkg(n)
        real(4) :: toa, tom, toms
        integer :: iflag

        ! Convert units for internal use
        sst_c = sst_k - 273.15
        slp_mb = slp_pa * 0.01
        t_c = t_k - 273.15
        r_gkg = r_kgkg * 1000.0

        ! Call modernized version for overall results
        call calculate_pi_single_profile(sst_k, slp_pa, p, t_k, r_kgkg, &
                                        n, n, pmin, vmax, ifl)

        ! Calculate saturation vapor pressure
        es0 = 6.112 * exp(17.67 * sst_c / (243.5 + sst_c))

        ! Environmental CAPE
        tp = t_k(1)
        rp = r_kgkg(1)
        pp = p(1)
        call CAPE(tp, rp, pp, t_k, r_kgkg, p, n, n, 0.0, cape_env, toa, iflag)
        tob = toa

        ! CAPE at radius of maximum winds
        tp = t_k(1)
        pp = min(pmin, 1000.0)
        rp = 0.622 * r_kgkg(1) * slp_mb / (pp * (0.622 + r_kgkg(1)) - r_kgkg(1) * slp_mb)
        call CAPE(tp, rp, pp, t_k, r_kgkg, p, n, n, 0.0, cape_rm, tom, iflag)

        ! Saturated CAPE
        tp = sst_k
        pp = min(pmin, 1000.0)
        rp = 0.622 * es0 / (pp - es0)
        call CAPE(tp, rp, pp, t_k, r_kgkg, p, n, n, 0.0, cape_sat, toms, iflag)

    end subroutine debug_modern_cape

end program debug_cape_comparison
