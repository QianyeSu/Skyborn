subroutine cape_profile_c(nlev, pressure_hpa, temperature_c, dewpoint_c, cape, cin, lfc_p, el_p) bind(C, name="cape_profile_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: cape_profile_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: pressure_hpa(nlev)
  real(c_double), intent(in) :: temperature_c(nlev)
  real(c_double), intent(in) :: dewpoint_c(nlev)
  real(c_double), intent(out) :: cape, cin, lfc_p, el_p
  call cape_profile_impl(int(nlev), pressure_hpa, temperature_c, dewpoint_c, cape, cin, lfc_p, el_p)
end subroutine cape_profile_c

subroutine cape_grid_c(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, cape2, cin2, lfc2, el2) bind(C, name="cape_grid_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: cape_grid_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
  real(c_double), intent(out) :: cape2(nlat, nlon)
  real(c_double), intent(out) :: cin2(nlat, nlon)
  real(c_double), intent(out) :: lfc2(nlat, nlon)
  real(c_double), intent(out) :: el2(nlat, nlon)
  call cape_grid_impl(int(nlev), int(nlat), int(nlon), pressure_3d, t_3d, td_3d, cape2, cin2, lfc2, el2)
end subroutine cape_grid_c

subroutine parcel_profile_c(nlev, pressure_hpa, temperature_c, dewpoint_c, t_par_c) bind(C, name="parcel_profile_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: parcel_profile_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: pressure_hpa(nlev)
  real(c_double), intent(in) :: temperature_c(nlev)
  real(c_double), intent(in) :: dewpoint_c(nlev)
  real(c_double), intent(out) :: t_par_c(nlev)
  call parcel_profile_impl(int(nlev), pressure_hpa, temperature_c, dewpoint_c, t_par_c)
end subroutine parcel_profile_c

subroutine parcel_profile_grid_c(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, out_3d) bind(C, name="parcel_profile_grid_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: parcel_profile_grid_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
  real(c_double), intent(out) :: out_3d(nlev, nlat, nlon)
  call parcel_profile_grid_impl(int(nlev), int(nlat), int(nlon), pressure_3d, t_3d, td_3d, out_3d)
end subroutine parcel_profile_grid_c

subroutine most_unstable_parcel_c(nlev, pressure_hpa, temperature_c, dewpoint_c, depth_hpa, mup_p, mup_t, mup_td, mup_idx) bind(C, name="most_unstable_parcel_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: most_unstable_parcel_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: pressure_hpa(nlev)
  real(c_double), intent(in) :: temperature_c(nlev)
  real(c_double), intent(in) :: dewpoint_c(nlev)
  real(c_double), value, intent(in) :: depth_hpa
  real(c_double), intent(out) :: mup_p, mup_t, mup_td
  integer(c_int), intent(out) :: mup_idx
  call most_unstable_parcel_impl(int(nlev), pressure_hpa, temperature_c, dewpoint_c, depth_hpa, mup_p, mup_t, mup_td, mup_idx)
end subroutine most_unstable_parcel_c

subroutine most_unstable_parcel_grid_c(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, depth_hpa, out_p3, out_t3, out_td3, out_idx3) bind(C, name="most_unstable_parcel_grid_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: most_unstable_parcel_grid_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
  real(c_double), value, intent(in) :: depth_hpa
  real(c_double), intent(out) :: out_p3(nlev, nlat, nlon)
  real(c_double), intent(out) :: out_t3(nlev, nlat, nlon)
  real(c_double), intent(out) :: out_td3(nlev, nlat, nlon)
  integer(c_int), intent(out) :: out_idx3(nlev, nlat, nlon)
  call most_unstable_parcel_grid_impl(int(nlev), int(nlat), int(nlon), pressure_3d, t_3d, td_3d, depth_hpa, out_p3, out_t3, out_td3, out_idx3)
end subroutine most_unstable_parcel_grid_c

subroutine mucape_c(nlev, pressure_hpa, temperature_c, dewpoint_c, depth_hpa, cape, cin, lfc_p, el_p) bind(C, name="mucape_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: mucape_profile_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: pressure_hpa(nlev)
  real(c_double), intent(in) :: temperature_c(nlev)
  real(c_double), intent(in) :: dewpoint_c(nlev)
  real(c_double), value, intent(in) :: depth_hpa
  real(c_double), intent(out) :: cape, cin, lfc_p, el_p
  call mucape_profile_impl(int(nlev), pressure_hpa, temperature_c, dewpoint_c, depth_hpa, cape, cin, lfc_p, el_p)
end subroutine mucape_c

subroutine mucape_grid_c(nlev, nlat, nlon, pressure_3d, t_3d, td_3d, depth_hpa, cape2, cin2, lfc2, el2) bind(C, name="mucape_grid_c")
  use iso_c_binding, only: c_double, c_int
  use cape_mod, only: mucape_grid_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: pressure_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: t_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: td_3d(nlev, nlat, nlon)
  real(c_double), value, intent(in) :: depth_hpa
  real(c_double), intent(out) :: cape2(nlat, nlon)
  real(c_double), intent(out) :: cin2(nlat, nlon)
  real(c_double), intent(out) :: lfc2(nlat, nlon)
  real(c_double), intent(out) :: el2(nlat, nlon)
  call mucape_grid_impl(int(nlev), int(nlat), int(nlon), pressure_3d, t_3d, td_3d, depth_hpa, cape2, cin2, lfc2, el2)
end subroutine mucape_grid_c
