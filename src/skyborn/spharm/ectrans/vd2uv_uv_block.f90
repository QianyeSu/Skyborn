module vd2uv_uv_block_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use PARKIND1, only : JPIM, JPRB
  use TPM_DISTR, only : D
  use PREPSNM_MOD, only : PREPSNM
  use PRFI1B_MOD, only : PRFI1B
  use VDTUV_MOD, only : VDTUV
  use runtime_state_mod, only : ensure_runtime_state
  implicit none
  private

  public :: vd2uv_uv_block

contains

  subroutine vd2uv_uv_block( &
    ntrunc, km, nt, rsphere, &
    vrtspec_r, vrtspec_i, divspec_r, divspec_i, &
    poa1_out, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: km
    integer(c_int), value, intent(in) :: nt
    real(c_double), value, intent(in) :: rsphere
    real(c_double), intent(in) :: vrtspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: vrtspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(in) :: divspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * nt)
    real(c_double), intent(out) :: poa1_out( &
      (ntrunc + 4_c_int + mod(ntrunc + 5_c_int, 2_c_int)) * (4_c_int * nt))
    integer(c_int), intent(out) :: ierror

    integer(kind=JPIM) :: ncoeff, nspec2, nlei1, ierr_local
    integer(kind=JPIM) :: m, n, it, idx, inm, row, col
    real(kind=JPRB), allocatable :: pspvor(:,:), pspdiv(:,:), zia(:,:), zep(:)

    poa1_out = 0.0_c_double

    if (ntrunc < 0_c_int .or. km < 0_c_int .or. km > ntrunc .or. nt < 1_c_int) then
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
    nlei1 = int(ntrunc, JPIM) + 4_JPIM + mod(int(ntrunc, JPIM) + 5_JPIM, 2_JPIM)

    allocate(pspvor(int(nt, JPIM), nspec2), pspdiv(int(nt, JPIM), nspec2))
    allocate(zia(nlei1, 8_JPIM * int(nt, JPIM)))
    allocate(zep(0:int(ntrunc, JPIM) + 2_JPIM))
    pspvor = 0.0_JPRB
    pspdiv = 0.0_JPRB
    zia = 0.0_JPRB

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

    call PREPSNM(int(km, JPIM), int(km, JPIM) + 1_JPIM, zep)
    call PRFI1B(int(km, JPIM), zia(:, 1:2 * int(nt, JPIM)), pspvor, int(nt, JPIM))
    call PRFI1B(int(km, JPIM), zia(:, 2 * int(nt, JPIM) + 1:4 * int(nt, JPIM)), pspdiv, int(nt, JPIM))
    call VDTUV( &
      int(km, JPIM), int(nt, JPIM), zep, &
      zia(:, 1:2 * int(nt, JPIM)), zia(:, 2 * int(nt, JPIM) + 1:4 * int(nt, JPIM)), &
      zia(:, 4 * int(nt, JPIM) + 1:6 * int(nt, JPIM)), zia(:, 6 * int(nt, JPIM) + 1:8 * int(nt, JPIM)))

    do row = 1_JPIM, nlei1
      do col = 1_JPIM, 4_JPIM * int(nt, JPIM)
        poa1_out((row - 1_JPIM) * (4_JPIM * int(nt, JPIM)) + col) = &
          zia(row, 4_JPIM * int(nt, JPIM) + col)
      end do
    end do

    deallocate(pspvor, pspdiv, zia, zep)
    ierror = 0_c_int
  end subroutine vd2uv_uv_block

end module vd2uv_uv_block_mod
