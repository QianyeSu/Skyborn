module vd2uv_pilot_state_mod
  use PARKIND1, only : JPIM, JPRB
  use TPM_CONSTANTS, only : RA
  use TPM_DIM, only : DIM_RESOL, R
  use TPM_DISTR, only : DISTR_RESOL, D, NPRTRV, MYSETV
  use TPM_FIELDS, only : FIELDS_RESOL, F
  use SUWAVEDI_MOD, only : SUWAVEDI
  implicit none
  private

  public :: ensure_vd2uv_pilot_state

  integer(kind=JPIM), save :: cached_ntrunc = -1
  real(kind=JPRB), save :: cached_rsphere = -1.0_JPRB

contains

  subroutine ensure_vd2uv_pilot_state(ntrunc, rsphere, ierror)
    integer(kind=JPIM), intent(in) :: ntrunc
    real(kind=JPRB), intent(in) :: rsphere
    integer(kind=JPIM), intent(out) :: ierror

    integer(kind=JPIM) :: numpp(1), myms_all(ntrunc + 1)
    integer(kind=JPIM) :: jm, jn, im, icount, inm

    ierror = 0
    if (associated(R) .and. associated(D) .and. associated(F)) then
      if (cached_ntrunc == ntrunc .and. cached_rsphere == rsphere) return
    end if

    call clear_vd2uv_pilot_state()

    allocate(DIM_RESOL(1))
    allocate(DISTR_RESOL(1))
    allocate(FIELDS_RESOL(1))
    R => DIM_RESOL(1)
    D => DISTR_RESOL(1)
    F => FIELDS_RESOL(1)

    RA = rsphere
    NPRTRV = 1
    MYSETV = 1

    R%NSMAX = ntrunc
    R%NTMAX = ntrunc
    R%NLEI1 = R%NSMAX + 4 + mod(R%NSMAX + 5, 2)

    allocate(D%NASM0(0:ntrunc))
    call SUWAVEDI(ntrunc, ntrunc, 1, 1, KASM0=D%NASM0, KUMPP=numpp, KMYMS=myms_all)
    D%NUMP = numpp(1)
    allocate(D%MYMS(D%NUMP))
    D%MYMS(:) = myms_all(1:D%NUMP)
    allocate(D%NPMT(0:ntrunc))

    inm = 0
    do jm = 1, D%NUMP
      im = D%MYMS(jm)
      D%NPMT(im) = inm
      inm = inm + R%NTMAX + 2 - im
    end do

    icount = 0
    do jm = 1, D%NUMP
      im = D%MYMS(jm)
      do jn = im, R%NTMAX + 2
        icount = icount + 1
      end do
    end do

    allocate(F%REPSNM(icount))
    allocate(F%RN(-1:R%NTMAX + 3))
    allocate(F%RLAPIN(-1:R%NSMAX + 2))

    icount = 0
    do jm = 1, D%NUMP
      im = D%MYMS(jm)
      do jn = im, R%NTMAX + 2
        icount = icount + 1
        F%REPSNM(icount) = sqrt(real(jn * jn - im * im, JPRB) / real(4 * jn * jn - 1, JPRB))
      end do
    end do

    do jn = -1, R%NTMAX + 3
      F%RN(jn) = real(jn, JPRB)
    end do

    F%RLAPIN(:) = 0.0_JPRB
    F%RLAPIN(-1) = 0.0_JPRB
    F%RLAPIN(0) = 0.0_JPRB
    do jn = 1, R%NSMAX + 2
      F%RLAPIN(jn) = -(RA * RA / real(jn * (jn + 1), JPRB))
    end do

    cached_ntrunc = ntrunc
    cached_rsphere = rsphere
  end subroutine ensure_vd2uv_pilot_state

  subroutine clear_vd2uv_pilot_state()
    if (associated(F)) then
      if (allocated(F%REPSNM)) deallocate(F%REPSNM)
      if (allocated(F%RN)) deallocate(F%RN)
      if (allocated(F%RLAPIN)) deallocate(F%RLAPIN)
    end if
    if (associated(D)) then
      if (allocated(D%MYMS)) deallocate(D%MYMS)
      if (allocated(D%NASM0)) deallocate(D%NASM0)
      if (allocated(D%NPMT)) deallocate(D%NPMT)
    end if
    if (allocated(FIELDS_RESOL)) deallocate(FIELDS_RESOL)
    if (allocated(DISTR_RESOL)) deallocate(DISTR_RESOL)
    if (allocated(DIM_RESOL)) deallocate(DIM_RESOL)
    nullify(F)
    nullify(D)
    nullify(R)
  end subroutine clear_vd2uv_pilot_state
end module vd2uv_pilot_state_mod
