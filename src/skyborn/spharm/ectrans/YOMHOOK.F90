module YOMHOOK
  use PARKIND1, only : JPHOOK
  implicit none

  logical :: LHOOK = .false.

contains

  subroutine DR_HOOK(CDNAME, KFLAG, PHANDLE)
    character(len=*), intent(in) :: CDNAME
    integer, intent(in) :: KFLAG
    real(kind=JPHOOK), intent(inout) :: PHANDLE

    PHANDLE = PHANDLE
  end subroutine DR_HOOK
end module YOMHOOK
