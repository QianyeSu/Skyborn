module poa1_to_vordiv_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use PARKIND1, only : JPIM, JPRB
  use TPM_DISTR, only : D
  use PREPSNM_MOD, only : PREPSNM
  use UPDSPB_MOD, only : UPDSPB
  use UVTVD_MOD, only : UVTVD
  use runtime_state_mod, only : ensure_runtime_state
  implicit none
  private

  public :: poa1_to_vordiv

contains

  subroutine poa1_to_vordiv( &
    ntrunc, km, kf_uv, rsphere, poa1_in, &
    vrtspec_r, vrtspec_i, divspec_r, divspec_i, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: km
    integer(c_int), value, intent(in) :: kf_uv
    real(c_double), value, intent(in) :: rsphere
    real(c_double), intent(in) :: poa1_in((ntrunc + 4_c_int + mod(ntrunc + 5_c_int, 2_c_int)) * (4_c_int * kf_uv))
    real(c_double), intent(out) :: vrtspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * kf_uv)
    real(c_double), intent(out) :: vrtspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * kf_uv)
    real(c_double), intent(out) :: divspec_r(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * kf_uv)
    real(c_double), intent(out) :: divspec_i(((ntrunc + 1_c_int) * (ntrunc + 2_c_int) / 2_c_int) * kf_uv)
    integer(c_int), intent(out) :: ierror

    integer(kind=JPIM) :: ierr_local, ncoeff, nspec2, nlei1
    integer(kind=JPIM) :: row, col, m, n, it, idx, inm
    real(kind=JPRB), allocatable :: zep(:)
    real(kind=JPRB), allocatable :: pu(:,:), pv(:,:), pvor(:,:), pdiv(:,:)
    real(kind=JPRB), allocatable :: pspvor(:,:), pspdiv(:,:)

    vrtspec_r = 0.0_c_double
    vrtspec_i = 0.0_c_double
    divspec_r = 0.0_c_double
    divspec_i = 0.0_c_double

    if (ntrunc < 0_c_int .or. km < 0_c_int .or. km > ntrunc .or. kf_uv < 1_c_int) then
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

    allocate(zep(0:int(ntrunc, JPIM) + 2_JPIM))
    allocate(pu(nlei1, 2 * int(kf_uv, JPIM)))
    allocate(pv(nlei1, 2 * int(kf_uv, JPIM)))
    allocate(pvor(nlei1, 2 * int(kf_uv, JPIM)))
    allocate(pdiv(nlei1, 2 * int(kf_uv, JPIM)))
    allocate(pspvor(int(kf_uv, JPIM), nspec2))
    allocate(pspdiv(int(kf_uv, JPIM), nspec2))

    pu = 0.0_JPRB
    pv = 0.0_JPRB
    pvor = 0.0_JPRB
    pdiv = 0.0_JPRB
    pspvor = 0.0_JPRB
    pspdiv = 0.0_JPRB

    do row = 1_JPIM, nlei1
      do col = 1_JPIM, 2_JPIM * int(kf_uv, JPIM)
        pu(row, col) = poa1_in((row - 1_JPIM) * (4_JPIM * int(kf_uv, JPIM)) + col)
        pv(row, col) = poa1_in((row - 1_JPIM) * (4_JPIM * int(kf_uv, JPIM)) + 2_JPIM * int(kf_uv, JPIM) + col)
      end do
    end do

    call PREPSNM(int(km, JPIM), int(km, JPIM) + 1_JPIM, zep)
    call UVTVD(int(km, JPIM), int(kf_uv, JPIM), zep, pu, pv, pvor, pdiv)
    call UPDSPB(int(km, JPIM), int(kf_uv, JPIM), pvor, pspvor)
    call UPDSPB(int(km, JPIM), int(kf_uv, JPIM), pdiv, pspdiv)

    if (km == 0_c_int) then
      pspvor(:, D%NASM0(0)) = 0.0_JPRB
      pspdiv(:, D%NASM0(0)) = 0.0_JPRB
    end if

    idx = 0_JPIM
    do m = 0_JPIM, int(ntrunc, JPIM)
      do n = m, int(ntrunc, JPIM)
        idx = idx + 1_JPIM
        inm = D%NASM0(m) + 2_JPIM * (n - m)
        do it = 1_JPIM, int(kf_uv, JPIM)
          vrtspec_r((idx - 1_JPIM) * int(kf_uv, JPIM) + it) = pspvor(it, inm)
          vrtspec_i((idx - 1_JPIM) * int(kf_uv, JPIM) + it) = pspvor(it, inm + 1_JPIM)
          divspec_r((idx - 1_JPIM) * int(kf_uv, JPIM) + it) = pspdiv(it, inm)
          divspec_i((idx - 1_JPIM) * int(kf_uv, JPIM) + it) = pspdiv(it, inm + 1_JPIM)
        end do
      end do
    end do

    deallocate(zep, pu, pv, pvor, pdiv, pspvor, pspdiv)
    ierror = 0_c_int
  end subroutine poa1_to_vordiv

end module poa1_to_vordiv_mod
