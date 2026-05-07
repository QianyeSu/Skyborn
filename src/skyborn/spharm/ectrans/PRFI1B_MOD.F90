module PRFI1B_MOD
contains
  subroutine PRFI1B(KM, PIA, PSPEC, KFIELDS, KFLDPTR)
    use PARKIND1, only : JPIM, JPRB
    use TPM_DIM, only : R
    use TPM_DISTR, only : D
    implicit none

    integer(kind=JPIM), intent(in) :: KM, KFIELDS
    real(kind=JPRB), intent(in) :: PSPEC(:,:)
    real(kind=JPRB), intent(out) :: PIA(:,:)
    integer(kind=JPIM), intent(in), optional :: KFLDPTR(:)

    integer(kind=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF, IFLD

    ILCM = R%NSMAX + 1 - KM
    IOFF = D%NASM0(KM)

    if (present(KFLDPTR)) then
      do JFLD = 1, KFIELDS
        IR = 2 * (JFLD - 1) + 1
        II = IR + 1
        IFLD = KFLDPTR(JFLD)
        do J = 1, ILCM
          INM = IOFF + (ILCM - J) * 2
          PIA(J + 2, IR) = PSPEC(IFLD, INM)
          PIA(J + 2, II) = PSPEC(IFLD, INM + 1)
        end do
      end do
    else
      do J = 1, ILCM
        INM = IOFF + (ILCM - J) * 2
        do JFLD = 1, KFIELDS
          IR = 2 * (JFLD - 1) + 1
          II = IR + 1
          PIA(J + 2, IR) = PSPEC(JFLD, INM)
          PIA(J + 2, II) = PSPEC(JFLD, INM + 1)
        end do
      end do
    end if

    do JFLD = 1, 2 * KFIELDS
      PIA(1, JFLD) = 0.0_JPRB
      PIA(2, JFLD) = 0.0_JPRB
      PIA(ILCM + 3, JFLD) = 0.0_JPRB
    end do
  end subroutine PRFI1B
end module PRFI1B_MOD
