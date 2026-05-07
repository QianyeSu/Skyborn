module TPM_DIM
  use PARKIND1, only : JPIM
  implicit none

  type DIM_TYPE
    integer(kind=JPIM) :: NSMAX = 0
    integer(kind=JPIM) :: NTMAX = 0
    integer(kind=JPIM) :: NLEI1 = 0
    integer(kind=JPIM) :: NDGL = 0
    integer(kind=JPIM) :: NDGNH = 0
  end type DIM_TYPE

  type(DIM_TYPE), allocatable, target :: DIM_RESOL(:)
  type(DIM_TYPE), pointer :: R => null()
end module TPM_DIM
