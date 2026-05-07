module PREPSNM_MOD
contains
  subroutine PREPSNM(KM, KMLOC, PEPSNM)
    use PARKIND1, only : JPIM, JPRB
    use TPM_DIM, only : R
    use TPM_FIELDS, only : F
    use TPM_DISTR, only : D
    implicit none

    integer(kind=JPIM), intent(in) :: KM, KMLOC
    real(kind=JPRB), intent(out) :: PEPSNM(0:R%NTMAX + 2)
    integer(kind=JPIM) :: JN

    if (KM > 0) PEPSNM(0:KM - 1) = 0.0_JPRB
    do JN = KM, R%NTMAX + 2
      PEPSNM(JN) = F%REPSNM(D%NPMT(KM) + KMLOC - KM + JN)
    end do
  end subroutine PREPSNM
end module PREPSNM_MOD
