! =============================================================================
! File    : tropopause_height_bindings.f90
! Purpose : C interoperability bindings split from tropopause_height.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine tropopause_grid_3d_c( &
    nlat, nlon, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
) bind(C, name="tropopause_grid_3d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use tropopause_height_mod, only : tropopause_grid_3d
    implicit none

    integer(c_int), value, intent(in) :: nlat, nlon, nlev, nlevm, punit
    type(c_ptr), value, intent(in) :: pfull, tfull, ptrop_hpa, htrop_m, itrop, lapse_rate, success
    real(c_double), value, intent(in) :: tmsg, lapsec

    integer :: nlat_f, nlon_f, nlev_f, nlevm_f
    real(c_double), pointer :: pfull_view(:)
    real(c_double), pointer :: tfull_view(:, :, :)
    real(c_double), pointer :: ptrop_hpa_view(:, :)
    real(c_double), pointer :: htrop_m_view(:, :)
    integer(c_int), pointer :: itrop_view(:, :)
    real(c_double), pointer :: lapse_rate_view(:, :)
    integer(c_int), pointer :: success_view(:, :)

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    nlev_f = int(nlev)
    nlevm_f = int(nlevm)

    call c_f_pointer(pfull, pfull_view, [nlev_f])
    call c_f_pointer(tfull, tfull_view, [nlat_f, nlon_f, nlev_f])
    call c_f_pointer(ptrop_hpa, ptrop_hpa_view, [nlat_f, nlon_f])
    call c_f_pointer(htrop_m, htrop_m_view, [nlat_f, nlon_f])
    call c_f_pointer(itrop, itrop_view, [nlat_f, nlon_f])
    call c_f_pointer(lapse_rate, lapse_rate_view, [nlat_f, nlon_f])
    call c_f_pointer(success, success_view, [nlat_f, nlon_f])

    call tropopause_grid_3d( &
        nlat_f, nlon_f, nlev_f, nlevm_f, pfull_view, tfull_view, tmsg, lapsec, int(punit), &
        ptrop_hpa_view, htrop_m_view, itrop_view, lapse_rate_view, success_view &
    )
end subroutine tropopause_grid_3d_c


subroutine tropopause_grid_4d_c( &
    nlat, nlon, nlev, ntime, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
) bind(C, name="tropopause_grid_4d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use tropopause_height_mod, only : tropopause_grid_4d
    implicit none

    integer(c_int), value, intent(in) :: nlat, nlon, nlev, ntime, nlevm, punit
    type(c_ptr), value, intent(in) :: pfull, tfull, ptrop_hpa, htrop_m, itrop, lapse_rate, success
    real(c_double), value, intent(in) :: tmsg, lapsec

    integer :: nlat_f, nlon_f, nlev_f, ntime_f, nlevm_f
    real(c_double), pointer :: pfull_view(:)
    real(c_double), pointer :: tfull_view(:, :, :, :)
    real(c_double), pointer :: ptrop_hpa_view(:, :, :)
    real(c_double), pointer :: htrop_m_view(:, :, :)
    integer(c_int), pointer :: itrop_view(:, :, :)
    real(c_double), pointer :: lapse_rate_view(:, :, :)
    integer(c_int), pointer :: success_view(:, :, :)

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    nlev_f = int(nlev)
    ntime_f = int(ntime)
    nlevm_f = int(nlevm)

    call c_f_pointer(pfull, pfull_view, [nlev_f])
    call c_f_pointer(tfull, tfull_view, [nlat_f, nlon_f, nlev_f, ntime_f])
    call c_f_pointer(ptrop_hpa, ptrop_hpa_view, [nlat_f, nlon_f, ntime_f])
    call c_f_pointer(htrop_m, htrop_m_view, [nlat_f, nlon_f, ntime_f])
    call c_f_pointer(itrop, itrop_view, [nlat_f, nlon_f, ntime_f])
    call c_f_pointer(lapse_rate, lapse_rate_view, [nlat_f, nlon_f, ntime_f])
    call c_f_pointer(success, success_view, [nlat_f, nlon_f, ntime_f])

    call tropopause_grid_4d( &
        nlat_f, nlon_f, nlev_f, ntime_f, nlevm_f, pfull_view, tfull_view, tmsg, lapsec, int(punit), &
        ptrop_hpa_view, htrop_m_view, itrop_view, lapse_rate_view, success_view &
    )
end subroutine tropopause_grid_4d_c


subroutine tropopause_grid_2d_c( &
    nspatial, nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
) bind(C, name="tropopause_grid_2d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use tropopause_height_mod, only : tropopause_grid_2d
    implicit none

    integer(c_int), value, intent(in) :: nspatial, nlev, nlevm, punit
    type(c_ptr), value, intent(in) :: pfull, tfull, ptrop_hpa, htrop_m, itrop, lapse_rate, success
    real(c_double), value, intent(in) :: tmsg, lapsec

    integer :: nspatial_f, nlev_f, nlevm_f
    real(c_double), pointer :: pfull_view(:)
    real(c_double), pointer :: tfull_view(:, :)
    real(c_double), pointer :: ptrop_hpa_view(:)
    real(c_double), pointer :: htrop_m_view(:)
    integer(c_int), pointer :: itrop_view(:)
    real(c_double), pointer :: lapse_rate_view(:)
    integer(c_int), pointer :: success_view(:)

    nspatial_f = int(nspatial)
    nlev_f = int(nlev)
    nlevm_f = int(nlevm)

    call c_f_pointer(pfull, pfull_view, [nlev_f])
    call c_f_pointer(tfull, tfull_view, [nspatial_f, nlev_f])
    call c_f_pointer(ptrop_hpa, ptrop_hpa_view, [nspatial_f])
    call c_f_pointer(htrop_m, htrop_m_view, [nspatial_f])
    call c_f_pointer(itrop, itrop_view, [nspatial_f])
    call c_f_pointer(lapse_rate, lapse_rate_view, [nspatial_f])
    call c_f_pointer(success, success_view, [nspatial_f])

    call tropopause_grid_2d( &
        nspatial_f, nlev_f, nlevm_f, pfull_view, tfull_view, tmsg, lapsec, int(punit), &
        ptrop_hpa_view, htrop_m_view, itrop_view, lapse_rate_view, success_view &
    )
end subroutine tropopause_grid_2d_c


subroutine tropopause_profile_1d_c( &
    nlev, nlevm, pfull, tfull, tmsg, lapsec, punit, &
    ptrop_hpa, htrop_m, itrop, lapse_rate, success &
) bind(C, name="tropopause_profile_1d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use tropopause_height_mod, only : tropopause_profile_1d
    implicit none

    integer(c_int), value, intent(in) :: nlev, nlevm, punit
    type(c_ptr), value, intent(in) :: pfull, tfull, ptrop_hpa, htrop_m, itrop, lapse_rate, success
    real(c_double), value, intent(in) :: tmsg, lapsec

    integer :: nlev_f, nlevm_f
    real(c_double), pointer :: pfull_view(:)
    real(c_double), pointer :: tfull_view(:)
    real(c_double), pointer :: ptrop_hpa_view
    real(c_double), pointer :: htrop_m_view
    integer(c_int), pointer :: itrop_view
    real(c_double), pointer :: lapse_rate_view
    integer(c_int), pointer :: success_view

    nlev_f = int(nlev)
    nlevm_f = int(nlevm)

    call c_f_pointer(pfull, pfull_view, [nlev_f])
    call c_f_pointer(tfull, tfull_view, [nlev_f])
    call c_f_pointer(ptrop_hpa, ptrop_hpa_view)
    call c_f_pointer(htrop_m, htrop_m_view)
    call c_f_pointer(itrop, itrop_view)
    call c_f_pointer(lapse_rate, lapse_rate_view)
    call c_f_pointer(success, success_view)

    call tropopause_profile_1d( &
        nlev_f, nlevm_f, pfull_view, tfull_view, tmsg, lapsec, int(punit), &
        ptrop_hpa_view, htrop_m_view, itrop_view, lapse_rate_view, success_view &
    )
end subroutine tropopause_profile_1d_c
