program diagnose_cape_detail
    ! Detailed diagnosis of CAPE calculation differences
    implicit none

    ! Test a single parcel ascent
    real(4) :: SST = 28.0  ! °C
    real(4) :: SLP = 1008.0  ! hPa
    real(4) :: PM_TEST = 978.42  ! Test pressure from PCMIN result

    ! Thermodynamic variables
    real(4) :: SSTK, ES0, RP_PCMIN, RP_MODERN
    real(4) :: R1, PSL
    real(4) :: RP_SAT_PCMIN, RP_SAT_MODERN
    real(4) :: PP_SAT
    real(4) :: T1, TV1_PCMIN, TV1_MODERN
    real(4) :: CAPEA, CAPEM, CAPEMS, CAT_PCMIN, CAT_MODERN
    real(4) :: CKCD, RAT
    real(4) :: TVAV, PMIN_PCMIN, PMIN_MODERN

    write(*,*) '========================================================================='
    write(*,*) '              DETAILED CAPE CALCULATION DIAGNOSIS'
    write(*,*) '========================================================================='
    write(*,*)

    ! Convert to internal units
    SSTK = SST + 273.15
    ES0 = 6.112 * exp(17.67 * SST / (243.5 + SST))
    PSL = SLP  ! Already in mb

    write(*,*) 'INPUT CONDITIONS:'
    write(*,*) '-----------------'
    write(*,'(A,F8.2,A)') ' SST:        ', SST, ' °C'
    write(*,'(A,F8.2,A)') ' SSTK:       ', SSTK, ' K'
    write(*,'(A,F8.2,A)') ' SLP:        ', SLP, ' hPa'
    write(*,'(A,F8.4,A)') ' ES0:        ', ES0, ' hPa'
    write(*,'(A,F8.2,A)') ' PM (test):  ', PM_TEST, ' hPa'
    write(*,*)

    ! Test case: Surface mixing ratio
    R1 = 0.020  ! 20 g/kg in decimal (0.020 kg/kg)

    write(*,*) 'MIXING RATIO CALCULATIONS AT RADIUS OF MAX WINDS:'
    write(*,*) '------------------------------------------------'
    write(*,'(A,F10.6)') ' Surface R (decimal): ', R1
    write(*,*)

    ! PCMIN formula (line 103 of pcmin_2013.f)
    ! RP=0.622*R(NK)*PSL/(PP*(0.622+R(NK))-R(NK)*PSL)
    RP_PCMIN = 0.622 * R1 * PSL / (PM_TEST * (0.622 + R1) - R1 * PSL)

    write(*,*) 'PCMIN CALCULATION:'
    write(*,'(A,F10.6)') ' Numerator:   0.622 * R1 * PSL = ', 0.622 * R1 * PSL
    write(*,'(A,F10.6)') ' Denom term1: PP * (0.622 + R1) = ', PM_TEST * (0.622 + R1)
    write(*,'(A,F10.6)') ' Denom term2: R1 * PSL = ', R1 * PSL
    write(*,'(A,F10.6)') ' Denominator: ', PM_TEST * (0.622 + R1) - R1 * PSL
    write(*,'(A,F10.6)') ' RP_PCMIN = ', RP_PCMIN
    write(*,*)

    ! Modern formula (should be identical)
    RP_MODERN = 0.622 * R1 * PSL / (PM_TEST * (0.622 + R1) - R1 * PSL)

    write(*,*) 'MODERN CALCULATION:'
    write(*,'(A,F10.6)') ' RP_MODERN = ', RP_MODERN
    write(*,*)

    write(*,*) 'SATURATED MIXING RATIO CALCULATIONS:'
    write(*,*) '------------------------------------'

    ! Saturated mixing ratio at PM
    ! RP=0.622*ES0/(PP-ES0)

    PP_SAT = min(PM_TEST, 1000.0)

    RP_SAT_PCMIN = 0.622 * ES0 / (PP_SAT - ES0)

    write(*,'(A,F8.4)') ' PP for saturation: ', PP_SAT
    write(*,'(A,F10.6)') ' Numerator:   0.622 * ES0 = ', 0.622 * ES0
    write(*,'(A,F10.6)') ' Denominator: PP - ES0 = ', PP_SAT - ES0
    write(*,'(A,F10.6)') ' RP_SAT = ', RP_SAT_PCMIN
    write(*,*)

    ! Check virtual temperature calculations
    write(*,*) 'VIRTUAL TEMPERATURE CALCULATIONS:'
    write(*,*) '--------------------------------'

    T1 = 300.15  ! Example temperature in K

    ! PCMIN formula: TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
    TV1_PCMIN = T1 * (1.0 + R1/0.622) / (1.0 + R1)

    ! Modern equivalent
    TV1_MODERN = T1 * (1.0 + R1/0.622) / (1.0 + R1)

    write(*,'(A,F8.2,A)') ' Test temperature: ', T1, ' K'
    write(*,'(A,F8.4)') ' TV (PCMIN):  ', TV1_PCMIN
    write(*,'(A,F8.4)') ' TV (Modern): ', TV1_MODERN
    write(*,'(A,F8.4)') ' Difference:  ', TV1_MODERN - TV1_PCMIN
    write(*,*)

    ! Check the key energy calculation
    write(*,*) 'ENERGY CALCULATION (CAT):'
    write(*,*) '------------------------'

    ! Example values (from debug output)
    CAPEA = 3365.16
    CAPEM = 4005.66
    CAPEMS = 6373.77
    CKCD = 0.9
    RAT = SSTK / 233.35  ! Using outflow temperature

    ! PCMIN formula (line 123): CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
    CAT_PCMIN = CAPEM - CAPEA + 0.5 * CKCD * RAT * (CAPEMS - CAPEM)

    write(*,'(A,F10.2)') ' CAPEA:  ', CAPEA
    write(*,'(A,F10.2)') ' CAPEM:  ', CAPEM
    write(*,'(A,F10.2)') ' CAPEMS: ', CAPEMS
    write(*,'(A,F10.4)') ' CKCD:   ', CKCD
    write(*,'(A,F10.4)') ' RAT:    ', RAT
    write(*,*)
    write(*,'(A,F10.2)') ' CAPEM - CAPEA = ', CAPEM - CAPEA
    write(*,'(A,F10.2)') ' CAPEMS - CAPEM = ', CAPEMS - CAPEM
    write(*,'(A,F10.2)') ' 0.5*CKCD*RAT*(CAPEMS-CAPEM) = ', 0.5 * CKCD * RAT * (CAPEMS - CAPEM)
    write(*,'(A,F10.2)') ' CAT_PCMIN = ', CAT_PCMIN
    write(*,*)

    ! Now with the values from modern version
    CAPEM = 3919.47
    CAPEMS = 6250.79

    CAT_MODERN = CAPEM - CAPEA + 0.5 * CKCD * RAT * (CAPEMS - CAPEM)

    write(*,*) 'With Modern CAPE values:'
    write(*,'(A,F10.2)') ' CAPEM (modern):  ', CAPEM
    write(*,'(A,F10.2)') ' CAPEMS (modern): ', CAPEMS
    write(*,'(A,F10.2)') ' CAT_MODERN = ', CAT_MODERN
    write(*,*)
    write(*,'(A,F10.2)') ' Difference in CAT: ', CAT_MODERN - CAT_PCMIN
    write(*,*)

    ! Calculate resulting pressure
    TVAV = 287.0  ! Example average virtual temperature

    PMIN_PCMIN = PSL * exp(-CAT_PCMIN / (287.04 * TVAV))
    PMIN_MODERN = PSL * exp(-CAT_MODERN / (287.04 * TVAV))

    write(*,*) 'PRESSURE CALCULATION:'
    write(*,*) '--------------------'
    write(*,'(A,F10.2)') ' PMIN (PCMIN):  ', PMIN_PCMIN
    write(*,'(A,F10.2)') ' PMIN (Modern): ', PMIN_MODERN
    write(*,'(A,F10.2)') ' Difference:    ', PMIN_MODERN - PMIN_PCMIN

    write(*,*)
    write(*,*) '========================================================================='
    write(*,*) 'KEY FINDING: The difference in CAPE values at Rmax and saturated'
    write(*,*) 'conditions leads to different CAT values, which then produces'
    write(*,*) 'different minimum pressures.'
    write(*,*) '========================================================================='

end program diagnose_cape_detail
