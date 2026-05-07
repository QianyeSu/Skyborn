module VD2UV_MOD
contains
  subroutine VD2UV(KM, KMLOC, KF_UV, KLEI2, PSPVOR, PSPDIV, PU, PV)
    use PARKIND1, only : JPIM, JPRB, JPHOOK
    use YOMHOOK, only : LHOOK, DR_HOOK
    use TPM_CONSTANTS
    use TPM_DIM, only : R
    use TPM_DISTR, only : D
    use PREPSNM_MOD, only : PREPSNM
    use PRFI1B_MOD, only : PRFI1B
    use VDTUV_MOD, only : VDTUV
    implicit none

    integer(kind=JPIM), intent(in) :: KM, KMLOC, KF_UV, KLEI2
    real(kind=JPRB), intent(in) :: PSPVOR(:,:), PSPDIV(:,:)
    real(kind=JPRB), intent(out) :: PU(:,:), PV(:,:)

    real(kind=JPRB) :: ZIA(R%NLEI1, KLEI2)
    real(kind=JPRB) :: ZEPSNM(0:R%NTMAX + 2), ZA_R
    integer(kind=JPIM) :: JFLD, IVORL, IVORU, IDIVL, IDIVU, IUL, IUU, IVL, IVU
    integer(kind=JPIM) :: ILCM, IOFF, J, INM, IR, II
    real(kind=JPHOOK) :: ZHOOK_HANDLE

    if (LHOOK) call DR_HOOK('VD2UV_MOD', 0, ZHOOK_HANDLE)

    call PREPSNM(KM, KMLOC, ZEPSNM)

    if (KF_UV > 0) then
      ZIA = 0.0_JPRB
      IVORL = 1
      IVORU = 2 * KF_UV
      IDIVL = 2 * KF_UV + 1
      IDIVU = 4 * KF_UV
      IUL = 4 * KF_UV + 1
      IUU = 6 * KF_UV
      IVL = 6 * KF_UV + 1
      IVU = 8 * KF_UV

      call PRFI1B(KM, ZIA(:, IVORL:IVORU), PSPVOR, KF_UV)
      call PRFI1B(KM, ZIA(:, IDIVL:IDIVU), PSPDIV, KF_UV)
      call VDTUV(KM, KF_UV, ZEPSNM, ZIA(:, IVORL:IVORU), ZIA(:, IDIVL:IDIVU), &
                 ZIA(:, IUL:IUU), ZIA(:, IVL:IVU))

      ILCM = R%NSMAX + 1 - KM
      IOFF = D%NASM0(KM)
      ZA_R = 1.0_JPRB / RA
      do J = 1, ILCM
        INM = IOFF + (ILCM - J) * 2
        do JFLD = 1, KF_UV
          IR = 2 * (JFLD - 1) + 1
          II = IR + 1
          PU(JFLD, INM) = ZIA(J + 2, IR + IUL - 1) * ZA_R
          PU(JFLD, INM + 1) = ZIA(J + 2, II + IUL - 1) * ZA_R
          PV(JFLD, INM) = ZIA(J + 2, IR + IVL - 1) * ZA_R
          PV(JFLD, INM + 1) = ZIA(J + 2, II + IVL - 1) * ZA_R
        end do
      end do
    end if

    if (LHOOK) call DR_HOOK('VD2UV_MOD', 1, ZHOOK_HANDLE)
  end subroutine VD2UV
end module VD2UV_MOD
