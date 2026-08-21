subroutine srh_profile_c(nlev, height_m, u_ms, v_ms, depth_m, bottom_m, storm_u_ms, storm_v_ms, srh_pos, srh_neg, srh_total) bind(C, name="srh_profile_c")
  use iso_c_binding, only: c_double, c_int
  use srh_mod, only: srh_profile_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev
  real(c_double), intent(in) :: height_m(nlev)
  real(c_double), intent(in) :: u_ms(nlev)
  real(c_double), intent(in) :: v_ms(nlev)
  real(c_double), value, intent(in) :: depth_m, bottom_m, storm_u_ms, storm_v_ms
  real(c_double), intent(out) :: srh_pos, srh_neg, srh_total
  call srh_profile_impl(int(nlev), height_m, u_ms, v_ms, depth_m, bottom_m, storm_u_ms, storm_v_ms, srh_pos, srh_neg, srh_total)
end subroutine srh_profile_c

subroutine srh_grid_c(nlev, nlat, nlon, height_3d, u_3d, v_3d, depth_m, bottom_m, storm_u_ms, storm_v_ms, pos2, neg2, tot2) bind(C, name="srh_grid_c")
  use iso_c_binding, only: c_double, c_int
  use srh_mod, only: srh_grid_impl
  implicit none
  integer(c_int), value, intent(in) :: nlev, nlat, nlon
  real(c_double), intent(in) :: height_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: u_3d(nlev, nlat, nlon)
  real(c_double), intent(in) :: v_3d(nlev, nlat, nlon)
  real(c_double), value, intent(in) :: depth_m, bottom_m, storm_u_ms, storm_v_ms
  real(c_double), intent(out) :: pos2(nlat, nlon)
  real(c_double), intent(out) :: neg2(nlat, nlon)
  real(c_double), intent(out) :: tot2(nlat, nlon)
  call srh_grid_impl(int(nlev), int(nlat), int(nlon), height_3d, u_3d, v_3d, depth_m, bottom_m, storm_u_ms, storm_v_ms, pos2, neg2, tot2)
end subroutine srh_grid_c
