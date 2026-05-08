module uv_to_vordiv_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use PARKIND1, only : JPIM, JPRB
  use TPM_CONSTANTS, only : RA
  use TPM_DISTR, only : D
  use PREPSNM_MOD, only : PREPSNM
  use PRFI1B_MOD, only : PRFI1B
  use UPDSPB_MOD, only : UPDSPB
  use UVTVD_MOD, only : UVTVD
  use runtime_state_mod, only : ensure_runtime_state
  implicit none
  private

  public :: uv_to_vordiv

contains

  subroutine uv_to_vordiv( &
    ntrunc, nt, rsphere, &
    uspec_r, uspec_i, vspec_r, vspec_i, &
    vrtspec_r, vrtspec_i, divspec_r, divspec_i, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_double), value, intent(in) :: rsphere
    real(c_double), intent(in) :: uspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: uspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: vspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: vspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: vrtspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: vrtspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: divspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    integer(c_int), intent(out) :: ierror

    integer(kind=JPIM) :: ncoeff, nspec2, ierr_local
    integer(kind=JPIM) :: m, n, it, idx, inm
    real(kind=JPRB) :: za_inv
    real(kind=JPRB), allocatable :: pspu(:,:), pspv(:,:), pspvor(:,:), pspdiv(:,:), zia(:,:), zep(:)

    vrtspec_r = 0.0_c_double
    vrtspec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double

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
    allocate(pspu(nt, nspec2), pspv(nt, nspec2), pspvor(nt, nspec2), pspdiv(nt, nspec2))
    pspu = 0.0_JPRB
    pspv = 0.0_JPRB
    pspvor = 0.0_JPRB
    pspdiv = 0.0_JPRB
    za_inv = 1.0_JPRB / RA

    idx = 0_JPIM
    do m = 0_JPIM, int(ntrunc, JPIM)
      do n = m, int(ntrunc, JPIM)
        idx = idx + 1_JPIM
        inm = D%NASM0(m) + 2_JPIM * (n - m)
        do it = 1_JPIM, int(nt, JPIM)
          pspu(it, inm) = uspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) * za_inv
          pspu(it, inm + 1_JPIM) = uspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) * za_inv
          pspv(it, inm) = vspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) * za_inv
          pspv(it, inm + 1_JPIM) = vspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) * za_inv
        end do
      end do
    end do

    allocate(zia(int(ntrunc, JPIM) + 4_JPIM + mod(int(ntrunc, JPIM) + 5_JPIM, 2_JPIM), 8_JPIM * int(nt, JPIM)))
    allocate(zep(0:int(ntrunc, JPIM) + 2_JPIM))
    do m = 0_JPIM, int(ntrunc, JPIM)
      zia = 0.0_JPRB
      call PREPSNM(m, m + 1_JPIM, zep)
      call PRFI1B(m, zia(:, 1:2 * int(nt, JPIM)), pspu, int(nt, JPIM))
      call PRFI1B(m, zia(:, 2 * int(nt, JPIM) + 1:4 * int(nt, JPIM)), pspv, int(nt, JPIM))
      call UVTVD(m, int(nt, JPIM), zep, zia(:, 1:2 * int(nt, JPIM)), zia(:, 2 * int(nt, JPIM) + 1:4 * int(nt, JPIM)), &
                 zia(:, 4 * int(nt, JPIM) + 1:6 * int(nt, JPIM)), &
                 zia(:, 6 * int(nt, JPIM) + 1:8 * int(nt, JPIM)))
      call UPDSPB(m, int(nt, JPIM), zia(:, 4 * int(nt, JPIM) + 1:6 * int(nt, JPIM)), pspvor)
      call UPDSPB(m, int(nt, JPIM), zia(:, 6 * int(nt, JPIM) + 1:8 * int(nt, JPIM)), pspdiv)
    end do

    idx = 0_JPIM
    do m = 0_JPIM, int(ntrunc, JPIM)
      do n = m, int(ntrunc, JPIM)
        idx = idx + 1_JPIM
        inm = D%NASM0(m) + 2_JPIM * (n - m)
        do it = 1_JPIM, int(nt, JPIM)
          vrtspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) = pspvor(it, inm)
          vrtspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) = pspvor(it, inm + 1_JPIM)
          divspec_r((idx - 1_JPIM) * int(nt, JPIM) + it) = pspdiv(it, inm)
          divspec_i((idx - 1_JPIM) * int(nt, JPIM) + it) = pspdiv(it, inm + 1_JPIM)
        end do
      end do
    end do

    deallocate(zia, zep, pspu, pspv, pspvor, pspdiv)
    ierror = 0_c_int
  end subroutine uv_to_vordiv
end module uv_to_vordiv_mod
