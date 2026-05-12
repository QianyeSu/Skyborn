subroutine dcape_profile_c(nlev, pressure_hpa, temperature_c, dewpoint_c, dcape) bind(C, name="dcape_profile_c")
  use iso_c_binding, only: c_double, c_int
  use dcape_mod, only: dcape_profile_impl
  implicit none

  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: pressure_hpa(nlev)
  real(c_double), intent(in) :: temperature_c(nlev)
  real(c_double), intent(in) :: dewpoint_c(nlev)
  real(c_double), intent(out) :: dcape

  call dcape_profile_impl(int(nlev), pressure_hpa, temperature_c, dewpoint_c, dcape)
end subroutine dcape_profile_c

subroutine dcape_grid_c(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, out) bind(C, name="dcape_grid_c")
  use iso_c_binding, only: c_double, c_int
  use dcape_mod, only: dcape_grid_impl
  implicit none

  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
  real(c_double), intent(out) :: out(nlat, nlon)

  call dcape_grid_impl(int(nlev), int(nlat), int(nlon), pressure_3d, t_3d, td_3d, out)
end subroutine dcape_grid_c
