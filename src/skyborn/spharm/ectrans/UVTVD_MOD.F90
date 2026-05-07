module UVTVD_MOD
contains
  subroutine UVTVD(KM, KFIELD, PEPSNM, PU, PV, PVOR, PDIV)
    use PARKIND1, only : JPIM, JPRB
    use TPM_DIM, only : R
    use TPM_FIELDS, only : F
    implicit none

    integer(kind=JPIM), intent(in) :: KFIELD, KM
    real(kind=JPRB), intent(in) :: PEPSNM(0:R%NTMAX + 2)
    real(kind=JPRB), intent(out) :: PVOR(:,:), PDIV(:,:)
    real(kind=JPRB), intent(inout) :: PU(:,:), PV(:,:)

    integer(kind=JPIM) :: II, IN, IR, J, JN, ITMAX
    real(kind=JPRB) :: ZKM
    real(kind=JPRB) :: ZN(-1:R%NTMAX + 3)

    ZKM = real(KM, JPRB)
    ITMAX = R%NTMAX
    ZN(KM - 1:ITMAX + 3) = F%RN(KM - 1:ITMAX + 3)

    IN = F%NLTN(KM - 1)
    do J = 1, 2 * KFIELD
      PU(IN, J) = 0.0_JPRB
      PV(IN, J) = 0.0_JPRB
    end do

    if (KM /= 0) then
      do JN = KM, ITMAX
        IN = ITMAX + 2 - JN
        do J = 1, KFIELD
          IR = 2 * J - 1
          II = IR + 1
          PVOR(IN, IR) = -ZKM * PV(IN, II) - &
                          ZN(JN) * PEPSNM(JN + 1) * PU(IN - 1, IR) + &
                          ZN(JN + 1) * PEPSNM(JN) * PU(IN + 1, IR)
          PVOR(IN, II) = ZKM * PV(IN, IR) - &
                          ZN(JN) * PEPSNM(JN + 1) * PU(IN - 1, II) + &
                          ZN(JN + 1) * PEPSNM(JN) * PU(IN + 1, II)
          PDIV(IN, IR) = -ZKM * PU(IN, II) + &
                          ZN(JN) * PEPSNM(JN + 1) * PV(IN - 1, IR) - &
                          ZN(JN + 1) * PEPSNM(JN) * PV(IN + 1, IR)
          PDIV(IN, II) = ZKM * PU(IN, IR) + &
                          ZN(JN) * PEPSNM(JN + 1) * PV(IN - 1, II) - &
                          ZN(JN + 1) * PEPSNM(JN) * PV(IN + 1, II)
        end do
      end do
    else
      do JN = KM, ITMAX
        IN = ITMAX + 2 - JN
        do J = 1, KFIELD
          IR = 2 * J - 1
          PVOR(IN, IR) = -ZN(JN) * PEPSNM(JN + 1) * PU(IN - 1, IR) + &
                           ZN(JN + 1) * PEPSNM(JN) * PU(IN + 1, IR)
          PDIV(IN, IR) = ZN(JN) * PEPSNM(JN + 1) * PV(IN - 1, IR) - &
                           ZN(JN + 1) * PEPSNM(JN) * PV(IN + 1, IR)
        end do
      end do
    end if
  end subroutine UVTVD
end module UVTVD_MOD
