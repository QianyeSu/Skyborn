module vordiv_to_uv_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use PARKIND1, only : JPIM, JPRB
  use TPM_DISTR, only : D
  use VD2UV_CTL_MOD, only : VD2UV_CTL
  use runtime_state_mod, only : ensure_runtime_state
  implicit none
  private

  public :: vordiv_to_uv

contains

  subroutine vordiv_to_uv( &
    ntrunc, nt, rsphere, &
    vrtspec_r, vrtspec_i, divspec_r, divspec_i, &
    uspec_r, uspec_i, vspec_r, vspec_i, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), value, intent(in) :: rsphere
    real(c_double), intent(in) :: vrtspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: vrtspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: uspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: uspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: vspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: vspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    integer(c_int), intent(out) :: ierror

    integer(kind=JPIM) :: ncoeff, nspec2, ierr_local
    integer(kind=JPIM) :: m, n, it, idx, inm
    real(kind=JPRB), allocatable :: pspvor(:,:), pspdiv(:,:), pspu(:,:), pspv(:,:)

    uspec_r = 0.0_c_double
    uspec_i = 0.0_c_double
    vspec_r = 0.0_c_double
    vspec_i = 0.0_c_double

    if (ntrunc < 0_c_int .or. nt < 1_c_int) then
      ierror = 1_c_int
      return
    end if

    call ensure_runtime_state(int(ntrunc, JPIM), real(rsphere, JPRB), ierr_local)
    if (ierr_local /= 0_JPIM) then
      ierror = int(ierr_local, c_int)
      return
    end if

    ncoeff = (int(ntrunc, JPIM) + 1_JPIM) * (int(ntrunc, JPIM) + 2_JPIM) / 2_JPIM
    nspec2 = 2_JPIM * ncoeff
    allocate(pspvor(nt, nspec2), pspdiv(nt, nspec2), pspu(nt, nspec2), pspv(nt, nspec2))
    pspvor = 0.0_JPRB
    pspdiv = 0.0_JPRB
    pspu = 0.0_JPRB
    pspv = 0.0_JPRB

    idx = 0_JPIM
    do m = 0_JPIM, int(ntrunc, JPIM)
      do n = m, int(ntrunc, JPIM)
        idx = idx + 1_JPIM
        inm = D%NASM0(m) + 2_JPIM * (n - m)
        do it = 1_JPIM, int(nt, JPIM)
          pspvor(it, inm) = vrtspec_r((idx - 1_JPIM) * int(nt, JPIM) + it)
          pspvor(it, inm + 1_JPIM) = vrtspec_i((idx - 1_JPIM) * int(nt, JPIM) + it)
          pspdiv(it, inm) = divspec_r((idx - 1_JPIM) * int(nt, JPIM) + it)
          pspdiv(it, inm + 1_JPIM) = divspec_i((idx - 1_JPIM) * int(nt, JPIM) + it)
        end do
      end do
    end do

    call VD2UV_CTL(int(nt, JPIM), pspvor, pspdiv, pspu, pspv)

    idx = 0_JPIM
    do m = 0_JPIM, int(ntrunc, JPIM)
      do n = m, int(ntrunc, JPIM)
        idx = idx + 1_JPIM
        inm = D%NASM0(m) + 2_JPIM * (n - m)
        do it = 1_JPIM, int(nt, JPIM)
          uspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) = pspu(it, inm)
          uspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) = pspu(it, inm + 1_JPIM)
          vspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) = pspv(it, inm)
          vspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) = pspv(it, inm + 1_JPIM)
        end do
      end do
    end do

    deallocate(pspvor, pspdiv, pspu, pspv)
    ierror = 0_c_int
  end subroutine vordiv_to_uv
end module vordiv_to_uv_mod
