module TPM_FIELDS
  use PARKIND1, only : JPIM, JPRB
  implicit none

  type FIELDS_TYPE
    real(kind=JPRB), allocatable :: RW(:)
    real(kind=JPRB), allocatable :: RACTHE(:)
    real(kind=JPRB), allocatable :: REPSNM(:)
    real(kind=JPRB), allocatable :: RN(:)
    real(kind=JPRB), allocatable :: RLAPIN(:)
    integer(kind=JPIM), allocatable :: NLTN(:)
  end type FIELDS_TYPE

  type(FIELDS_TYPE), allocatable, target :: FIELDS_RESOL(:)
  type(FIELDS_TYPE), pointer :: F => null()
end module TPM_FIELDS
