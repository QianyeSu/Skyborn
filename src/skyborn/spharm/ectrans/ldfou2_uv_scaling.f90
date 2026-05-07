module ldfou2_uv_scaling_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  use PARKIND1, only : JPIM, JPRB
  use TPM_DIM, only : R
  use LDFOU2_MOD, only : LDFOU2
  use runtime_state_mod, only : ensure_runtime_state
  implicit none
  private

  public :: ldfou2_uv_scaling

contains

  subroutine ldfou2_uv_scaling( &
    ntrunc, km, kf_uv, rsphere, &
    paia_in, psia_in, paia_out, psia_out, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: km
    integer(c_int), value, intent(in) :: kf_uv
    real(c_double), value, intent(in) :: rsphere
    real(c_double), intent(in) :: paia_in((4_c_int * kf_uv) * ((ntrunc + 2_c_int) / 2_c_int))
    real(c_double), intent(in) :: psia_in((4_c_int * kf_uv) * ((ntrunc + 2_c_int) / 2_c_int))
    real(c_double), intent(out) :: paia_out((4_c_int * kf_uv) * ((ntrunc + 2_c_int) / 2_c_int))
    real(c_double), intent(out) :: psia_out((4_c_int * kf_uv) * ((ntrunc + 2_c_int) / 2_c_int))
    integer(c_int), intent(out) :: ierror

    integer(kind=JPIM) :: ierr_local, ndgnh, kled2
    integer(kind=JPIM) :: jgl, j, idx
    real(kind=JPRB), allocatable :: paia(:,:), psia(:,:)

    paia_out = 0.0_c_double
    psia_out = 0.0_c_double

    if (ntrunc < 0_c_int .or. km < 0_c_int .or. kf_uv < 1_c_int) then
      ierror = 1_c_int
      return
    end if

    call ensure_runtime_state(int(ntrunc, JPIM), real(rsphere, JPRB), ierr_local)
    if (ierr_local /= 0_JPIM) then
      ierror = int(ierr_local, c_int)
      return
    end if

    if (km > int(ntrunc, JPIM)) then
      ierror = 2_c_int
      return
    end if

    ndgnh = R%NDGNH
    kled2 = 4_JPIM * int(kf_uv, JPIM)
    allocate(paia(kled2, ndgnh), psia(kled2, ndgnh))

    do j = 1_JPIM, kled2
      do jgl = 1_JPIM, ndgnh
        idx = (j - 1_JPIM) * ndgnh + jgl
        paia(j, jgl) = paia_in(idx)
        psia(j, jgl) = psia_in(idx)
      end do
    end do

    call LDFOU2(int(km, JPIM), int(kf_uv, JPIM), paia, psia)

    do j = 1_JPIM, kled2
      do jgl = 1_JPIM, ndgnh
        idx = (j - 1_JPIM) * ndgnh + jgl
        paia_out(idx) = paia(j, jgl)
        psia_out(idx) = psia(j, jgl)
      end do
    end do

    deallocate(paia, psia)
    ierror = 0_c_int
  end subroutine ldfou2_uv_scaling
end module ldfou2_uv_scaling_mod
