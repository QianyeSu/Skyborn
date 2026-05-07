module runtime_state_mod
  use PARKIND1, only : JPIM, JPRB
  use TPM_CONSTANTS, only : RA
  use TPM_DIM, only : DIM_RESOL, R
  use TPM_DISTR, only : DISTR_RESOL, D, NPRTRV, MYSETV
  use TPM_FIELDS, only : FIELDS_RESOL, F
  use TPM_GEOMETRY, only : GEOM_RESOL, G
  use SUWAVEDI_MOD, only : SUWAVEDI
  implicit none
  private

  public :: ensure_runtime_state

  interface
    subroutine gaqd(nlat, theta, wts, dwork, ldwork, ierror)
      integer, intent(in) :: nlat, ldwork
      real(8), intent(out) :: theta(nlat), wts(nlat)
      real(8), intent(inout) :: dwork(ldwork)
      integer, intent(out) :: ierror
    end subroutine gaqd
  end interface

  integer(kind=JPIM), save :: cached_ntrunc = -1
  real(kind=JPRB), save :: cached_rsphere = -1.0_JPRB

contains

  subroutine ensure_runtime_state(ntrunc, rsphere, ierror)
    integer(kind=JPIM), intent(in) :: ntrunc
    real(kind=JPRB), intent(in) :: rsphere
    integer(kind=JPIM), intent(out) :: ierror

    integer(kind=JPIM) :: numpp(1), myms_all(ntrunc + 1)
    integer(kind=JPIM) :: jm, jn, im, icount, inm
    integer(kind=JPIM) :: ndgl, gaqd_ierror
    real(kind=JPRB) :: cache_tol, sin_theta
    real(kind=8), allocatable :: theta(:), wts(:), dwork(:)

    ierror = 0
    cache_tol = 10.0_JPRB * epsilon(1.0_JPRB) * max(1.0_JPRB, abs(rsphere))
    if (associated(R) .and. associated(D) .and. associated(F) .and. associated(G)) then
      if (cached_ntrunc == ntrunc .and. abs(cached_rsphere - rsphere) <= cache_tol) return
    end if

    call clear_runtime_state()

    allocate(DIM_RESOL(1))
    allocate(DISTR_RESOL(1))
    allocate(FIELDS_RESOL(1))
    allocate(GEOM_RESOL(1))
    R => DIM_RESOL(1)
    D => DISTR_RESOL(1)
    F => FIELDS_RESOL(1)
    G => GEOM_RESOL(1)

    RA = rsphere
    NPRTRV = 1
    MYSETV = 1

    ndgl = ntrunc + 1_JPIM
    R%NSMAX = ntrunc
    R%NTMAX = ntrunc
    R%NDGL = ndgl
    R%NDGNH = (ndgl + 1_JPIM) / 2_JPIM
    R%NLEI1 = R%NSMAX + 4 + mod(R%NSMAX + 5, 2)

    allocate(G%NDGLU(0:ntrunc))
    G%NDGLU = R%NDGNH

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
    allocate(F%NLTN(-1:R%NTMAX + 3))
    allocate(F%RW(R%NDGL))
    allocate(F%RACTHE(R%NDGL))

    allocate(theta(R%NDGL), wts(R%NDGL), dwork(R%NDGL * (R%NDGL + 2_JPIM)))
    call gaqd(int(R%NDGL), theta, wts, dwork, int(size(dwork)), gaqd_ierror)
    if (gaqd_ierror /= 0_JPIM) then
      deallocate(theta, wts, dwork)
      call clear_runtime_state()
      ierror = 10_JPIM + gaqd_ierror
      return
    end if
    do jn = 1, R%NDGL
      F%RW(jn) = real(wts(jn), JPRB)
      sin_theta = sin(real(theta(jn), JPRB))
      F%RACTHE(jn) = 1.0_JPRB / (RA * sin_theta)
    end do
    deallocate(theta, wts, dwork)

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
      F%NLTN(jn) = R%NTMAX + 2 - jn
    end do

    F%RLAPIN(:) = 0.0_JPRB
    F%RLAPIN(-1) = 0.0_JPRB
    F%RLAPIN(0) = 0.0_JPRB
    do jn = 1, R%NSMAX + 2
      F%RLAPIN(jn) = -(RA * RA / real(jn * (jn + 1), JPRB))
    end do

    cached_ntrunc = ntrunc
    cached_rsphere = rsphere
  end subroutine ensure_runtime_state

  subroutine clear_runtime_state()
    if (associated(F)) then
      if (allocated(F%RW)) deallocate(F%RW)
      if (allocated(F%RACTHE)) deallocate(F%RACTHE)
      if (allocated(F%REPSNM)) deallocate(F%REPSNM)
      if (allocated(F%RN)) deallocate(F%RN)
      if (allocated(F%RLAPIN)) deallocate(F%RLAPIN)
      if (allocated(F%NLTN)) deallocate(F%NLTN)
    end if
    if (associated(G)) then
      if (allocated(G%NDGLU)) deallocate(G%NDGLU)
    end if
    if (associated(D)) then
      if (allocated(D%MYMS)) deallocate(D%MYMS)
      if (allocated(D%NASM0)) deallocate(D%NASM0)
      if (allocated(D%NPMT)) deallocate(D%NPMT)
    end if
    if (allocated(FIELDS_RESOL)) deallocate(FIELDS_RESOL)
    if (allocated(GEOM_RESOL)) deallocate(GEOM_RESOL)
    if (allocated(DISTR_RESOL)) deallocate(DISTR_RESOL)
    if (allocated(DIM_RESOL)) deallocate(DIM_RESOL)
    nullify(F)
    nullify(G)
    nullify(D)
    nullify(R)
  end subroutine clear_runtime_state
end module runtime_state_mod
