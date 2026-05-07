module VD2UV_CTL_MOD
contains
  subroutine VD2UV_CTL(KF_UV, PSPVOR, PSPDIV, PU, PV)
    use PARKIND1, only : JPIM, JPRB
    use TPM_DISTR, only : D
    use VD2UV_MOD, only : VD2UV
    implicit none

    integer(kind=JPIM), intent(in) :: KF_UV
    real(kind=JPRB), intent(in) :: PSPVOR(:,:), PSPDIV(:,:)
    real(kind=JPRB), intent(out) :: PU(:,:), PV(:,:)

    integer(kind=JPIM) :: JM, IM, ILEI2

    call GSTATS(102, 0)
    ILEI2 = 8 * KF_UV
    call GSTATS(1647, 0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM, IM)
    do JM = 1, D%NUMP
      IM = D%MYMS(JM)
      call VD2UV(IM, JM, KF_UV, ILEI2, PSPVOR, PSPDIV, PU, PV)
    end do
!$OMP END PARALLEL DO
    call GSTATS(1647, 1)
    call GSTATS(102, 1)
  end subroutine VD2UV_CTL
end module VD2UV_CTL_MOD
