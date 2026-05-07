module SUWAVEDI_MOD
contains
  subroutine SUWAVEDI(KSMAX, KTMAX, KPRTRW, KMYSETW, KASM0, KSPOLEGL, KPROCM, &
                      KUMPP, KSPEC, KSPEC2, KSPEC2MX, KPOSSP, KMYMS, &
                      KPTRMS, KALLMS, KDIM0G)
    use PARKIND1, only : JPIM
    implicit none

    integer(kind=JPIM), intent(in) :: KSMAX, KTMAX, KPRTRW, KMYSETW
    integer(kind=JPIM), optional, intent(out) :: KSPEC, KSPEC2, KSPEC2MX, KSPOLEGL
    integer(kind=JPIM), optional, intent(out) :: KASM0(0:KSMAX), KPROCM(0:KSMAX)
    integer(kind=JPIM), optional, intent(out) :: KUMPP(KPRTRW), KMYMS(KSMAX + 1)
    integer(kind=JPIM), optional, intent(out) :: KPOSSP(KPRTRW + 1), KPTRMS(KPRTRW)
    integer(kind=JPIM), optional, intent(out) :: KALLMS(KSMAX + 1), KDIM0G(0:KSMAX)

    integer(kind=JPIM) :: IK, IL, IND, IPOS, ISPEC2P, JA, JM, JMLOC, IM
    integer(kind=JPIM) :: ISPOLEGL, ISPEC2MX_LOCAL
    integer(kind=JPIM) :: IASM0(0:KSMAX), IPROCM(0:KSMAX)
    integer(kind=JPIM) :: IUMPP(KPRTRW), IMYMS(KSMAX + 1), IPOSSP(KPRTRW + 1)
    integer(kind=JPIM) :: IPTRMS(KPRTRW), IALLMS(KSMAX + 1), IDIM0G(0:KSMAX)
    integer(kind=JPIM) :: ISPEC_LOCAL(KPRTRW), IC(KPRTRW)

    ISPEC_LOCAL(:) = 0
    IUMPP(:) = 0
    IASM0(:) = -99
    ISPOLEGL = 0

    IL = 1
    IND = 1
    IK = 0
    IPOS = 1
    do JM = 0, KSMAX
      IK = IK + IND
      if (IK > KPRTRW) then
        IK = KPRTRW
        IND = -1
      else if (IK < 1) then
        IK = 1
        IND = 1
      end if
      IPROCM(JM) = IK
      ISPEC_LOCAL(IK) = ISPEC_LOCAL(IK) + KSMAX - JM + 1
      IUMPP(IK) = IUMPP(IK) + 1
      if (IK == KMYSETW) then
        ISPOLEGL = ISPOLEGL + KTMAX + 1 - JM + 1
        IMYMS(IL) = JM
        IASM0(JM) = IPOS
        IPOS = IPOS + (KSMAX - JM + 1) * 2
        IL = IL + 1
      end if
    end do

    IPOSSP(1) = 1
    ISPEC2P = 2 * ISPEC_LOCAL(1)
    ISPEC2MX_LOCAL = ISPEC2P
    IPTRMS(1) = 1
    do JA = 2, KPRTRW
      IPOSSP(JA) = IPOSSP(JA - 1) + ISPEC2P
      ISPEC2P = 2 * ISPEC_LOCAL(JA)
      ISPEC2MX_LOCAL = max(ISPEC2MX_LOCAL, ISPEC2P)
      IPTRMS(JA) = IPTRMS(JA - 1) + IUMPP(JA - 1)
    end do
    IPOSSP(KPRTRW + 1) = IPOSSP(KPRTRW) + ISPEC2P

    IC(:) = 0
    do JM = 0, KSMAX
      IALLMS(IC(IPROCM(JM)) + IPTRMS(IPROCM(JM))) = JM
      IC(IPROCM(JM)) = IC(IPROCM(JM)) + 1
    end do

    IPOS = 1
    do JA = 1, KPRTRW
      do JMLOC = 1, IUMPP(JA)
        IM = IALLMS(IPTRMS(JA) + JMLOC - 1)
        IDIM0G(IM) = IPOS
        IPOS = IPOS + (KSMAX + 1 - IM) * 2
      end do
    end do

    if (present(KSPEC)) KSPEC = ISPEC_LOCAL(KMYSETW)
    if (present(KSPEC2)) KSPEC2 = 2 * ISPEC_LOCAL(KMYSETW)
    if (present(KSPEC2MX)) KSPEC2MX = ISPEC2MX_LOCAL
    if (present(KSPOLEGL)) KSPOLEGL = ISPOLEGL
    if (present(KASM0)) KASM0(:) = IASM0(:)
    if (present(KPROCM)) KPROCM(:) = IPROCM(:)
    if (present(KUMPP)) KUMPP(:) = IUMPP(:)
    if (present(KMYMS)) KMYMS(:) = IMYMS(:)
    if (present(KPOSSP)) KPOSSP(:) = IPOSSP(:)
    if (present(KPTRMS)) KPTRMS(:) = IPTRMS(:)
    if (present(KALLMS)) KALLMS(:) = IALLMS(:)
    if (present(KDIM0G)) KDIM0G(:) = IDIM0G(:)
  end subroutine SUWAVEDI
end module SUWAVEDI_MOD
