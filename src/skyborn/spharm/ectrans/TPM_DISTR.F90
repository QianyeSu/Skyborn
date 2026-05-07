module TPM_DISTR
  use PARKIND1, only : JPIM
  implicit none

  integer(kind=JPIM) :: NPRTRV = 1
  integer(kind=JPIM) :: MYSETV = 1

  type DISTR_TYPE
    integer(kind=JPIM) :: NUMP = 0
    integer(kind=JPIM), allocatable :: MYMS(:)
    integer(kind=JPIM), allocatable :: NASM0(:)
    integer(kind=JPIM), allocatable :: NPMT(:)
  end type DISTR_TYPE

  type(DISTR_TYPE), allocatable, target :: DISTR_RESOL(:)
  type(DISTR_TYPE), pointer :: D => null()
end module TPM_DISTR
