module ledir_dgemm_mod
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  implicit none
  private

  public :: ledir_dgemm

contains

  subroutine ledir_dgemm( &
    ntrunc, km, kfc, kdglu, &
    paia, psia, rpnma, rpnms, pw, poa1, ierror) bind(C)
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: km
    integer(c_int), value, intent(in) :: kfc
    integer(c_int), value, intent(in) :: kdglu
    real(c_double), intent(in) :: paia(kfc * kdglu)
    real(c_double), intent(in) :: psia(kfc * kdglu)
    real(c_double), intent(in) :: rpnma(kdglu * ((ntrunc - km + 2_c_int) / 2_c_int))
    real(c_double), intent(in) :: rpnms(kdglu * ((ntrunc - km + 3_c_int) / 2_c_int))
    real(c_double), intent(in) :: pw(kdglu)
    real(c_double), intent(out) :: poa1( &
      (ntrunc + 4_c_int + mod(ntrunc + 5_c_int, 2_c_int)) * kfc)
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ia, is, ila, ils, iskip
    integer(c_int) :: jk, j, col, row
    real(c_double) :: accum

    poa1 = 0.0_c_double

    if (ntrunc < 0_c_int .or. km < 0_c_int .or. km > ntrunc) then
      ierror = 1_c_int
      return
    end if
    if (kfc < 1_c_int .or. kdglu < 1_c_int) then
      ierror = 2_c_int
      return
    end if

    ia = 1_c_int + mod(ntrunc - km + 2_c_int, 2_c_int)
    is = 1_c_int + mod(ntrunc - km + 1_c_int, 2_c_int)
    ila = (ntrunc - km + 2_c_int) / 2_c_int
    ils = (ntrunc - km + 3_c_int) / 2_c_int
    iskip = 1_c_int
    if (km == 0_c_int) iskip = 2_c_int

    do jk = 1_c_int, kfc, iskip
      do col = 1_c_int, ila
        accum = 0.0_c_double
        do j = 1_c_int, kdglu
          accum = accum + rpnma((j - 1_c_int) * ila + col) * &
            paia((jk - 1_c_int) * kdglu + j) * pw(j)
        end do
        row = ia + (col - 1_c_int) * 2_c_int
        poa1((row - 1_c_int) * kfc + jk) = accum
      end do

      do col = 1_c_int, ils
        accum = 0.0_c_double
        do j = 1_c_int, kdglu
          accum = accum + rpnms((j - 1_c_int) * ils + col) * &
            psia((jk - 1_c_int) * kdglu + j) * pw(j)
        end do
        row = is + (col - 1_c_int) * 2_c_int
        poa1((row - 1_c_int) * kfc + jk) = accum
      end do
    end do

    ierror = 0_c_int
  end subroutine ledir_dgemm

end module ledir_dgemm_mod
