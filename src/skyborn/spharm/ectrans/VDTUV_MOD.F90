module VDTUV_MOD
contains
  subroutine VDTUV(KM, KFIELD, PEPSNM, PVOR, PDIV, PU, PV)
    use PARKIND1, only : JPIM, JPRB
    use TPM_DIM, only : R
    use TPM_FIELDS, only : F
    implicit none

    integer(kind=JPIM), intent(in) :: KM, KFIELD
    real(kind=JPRB), intent(in) :: PEPSNM(0:R%NTMAX + 2)
    real(kind=JPRB), intent(in) :: PVOR(:,:), PDIV(:,:)
    real(kind=JPRB), intent(out) :: PU(:,:), PV(:,:)

    integer(kind=JPIM) :: II, IJ, IR, J, JN, ISMAX, JI
    real(kind=JPRB) :: ZKM
    real(kind=JPRB) :: ZN(-1:R%NTMAX + 4)
    real(kind=JPRB) :: ZLAPIN(-1:R%NSMAX + 4)
    real(kind=JPRB) :: ZEPSNM(-1:R%NSMAX + 4)

    ZKM = real(KM, JPRB)
    ISMAX = R%NSMAX
    do JN = KM - 1, ISMAX + 2
      IJ = ISMAX + 3 - JN
      ZN(IJ) = F%RN(JN)
      ZLAPIN(IJ) = F%RLAPIN(JN)
      if (JN >= 0) ZEPSNM(IJ) = PEPSNM(JN)
    end do
    ZN(0) = F%RN(ISMAX + 3)

    if (KM == 0) then
      do J = 1, KFIELD
        IR = 2 * J - 1
        do JI = 2, ISMAX + 3 - KM
          PU(JI, IR) = ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PVOR(JI + 1, IR) - &
                        ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PVOR(JI - 1, IR)
          PV(JI, IR) = -ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PDIV(JI + 1, IR) + &
                         ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PDIV(JI - 1, IR)
        end do
      end do
    else
      do J = 1, KFIELD
        IR = 2 * J - 1
        II = IR + 1
        do JI = 2, ISMAX + 3 - KM
          PU(JI, IR) = -ZKM * ZLAPIN(JI) * PDIV(JI, II) + &
                        ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PVOR(JI + 1, IR) - &
                        ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PVOR(JI - 1, IR)
          PU(JI, II) = ZKM * ZLAPIN(JI) * PDIV(JI, IR) + &
                        ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PVOR(JI + 1, II) - &
                        ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PVOR(JI - 1, II)
          PV(JI, IR) = -ZKM * ZLAPIN(JI) * PVOR(JI, II) - &
                        ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PDIV(JI + 1, IR) + &
                        ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PDIV(JI - 1, IR)
          PV(JI, II) = ZKM * ZLAPIN(JI) * PVOR(JI, IR) - &
                        ZN(JI + 1) * ZEPSNM(JI) * ZLAPIN(JI + 1) * PDIV(JI + 1, II) + &
                        ZN(JI - 2) * ZEPSNM(JI - 1) * ZLAPIN(JI - 1) * PDIV(JI - 1, II)
        end do
      end do
    end if
  end subroutine VDTUV
end module VDTUV_MOD
