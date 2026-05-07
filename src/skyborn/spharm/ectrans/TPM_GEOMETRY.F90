module TPM_GEOMETRY
  use PARKIND1, only : JPIM
  implicit none

  type GEOM_TYPE
    integer(kind=JPIM), allocatable :: NDGLU(:)
  end type GEOM_TYPE

  type(GEOM_TYPE), allocatable, target :: GEOM_RESOL(:)
  type(GEOM_TYPE), pointer :: G => null()
end module TPM_GEOMETRY
