module gridmeta
  use, intrinsic :: iso_c_binding, only : c_int
  implicit none
  private

  public :: validate_nloen
  public :: build_grid_layout

contains

  subroutine validate_nloen(ndgl, nloen, ngptot, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), intent(out) :: ngptot
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ilat

    ngptot = 0_c_int
    ierror = 0_c_int

    if (ndgl < 3_c_int) then
      ierror = 1_c_int
      return
    end if

    do ilat = 1_c_int, ndgl
      if (nloen(ilat) <= 0_c_int) then
        ierror = 2_c_int
        ngptot = 0_c_int
        return
      end if
      ngptot = ngptot + nloen(ilat)
    end do
  end subroutine validate_nloen

  subroutine build_grid_layout( &
    ndgl, nloen, ngptot, lat_offsets, max_nlon, ierror) bind(C)
    integer(c_int), value, intent(in) :: ndgl
    integer(c_int), intent(in) :: nloen(ndgl)
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), intent(out) :: lat_offsets(ndgl)
    integer(c_int), intent(out) :: max_nlon
    integer(c_int), intent(out) :: ierror

    integer(c_int) :: ilat
    integer(c_int) :: running_total

    ierror = 0_c_int
    max_nlon = 0_c_int
    running_total = 0_c_int

    if (ndgl < 1_c_int) then
      ierror = 3_c_int
      return
    end if

    do ilat = 1_c_int, ndgl
      if (nloen(ilat) <= 0_c_int) then
        ierror = 4_c_int
        lat_offsets = 0_c_int
        max_nlon = 0_c_int
        return
      end if
      lat_offsets(ilat) = running_total
      running_total = running_total + nloen(ilat)
      if (nloen(ilat) > max_nlon) then
        max_nlon = nloen(ilat)
      end if
    end do

    if (running_total /= ngptot) then
      ierror = 5_c_int
      lat_offsets = 0_c_int
      max_nlon = 0_c_int
      return
    end if
  end subroutine build_grid_layout

end module gridmeta
