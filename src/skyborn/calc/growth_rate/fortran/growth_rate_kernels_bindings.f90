! =============================================================================
! File    : growth_rate_kernels_bindings.f90
! Purpose : C interoperability bindings split from growth_rate_kernels.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine dbarot_growth_rate_1d_c(lat, u, max_growth, ier, nlat) bind(C, name="dbarot_growth_rate_1d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: lat, u
    integer(c_int), value, intent(in) :: nlat
    real(c_double), intent(out) :: max_growth
    integer(c_int), intent(out) :: ier

    integer :: nlat_f
    real(real64), pointer :: lat_view(:), u_view(:)
    real(real64) :: max_growth_impl
    integer :: ier_impl

    nlat_f = int(nlat)
    call c_f_pointer(lat, lat_view, [nlat_f])
    call c_f_pointer(u, u_view, [nlat_f])

    call dbarot_growth_rate_1d(lat_view, u_view, max_growth_impl, ier_impl, nlat_f)
    max_growth = max_growth_impl
    ier = ier_impl
end subroutine dbarot_growth_rate_1d_c


subroutine dbaroc_growth_rate_1d_c( &
    u, theta, pressure, temperature, f_cor, beta, smooth_window, &
    wavenumber_mode, wavenumber_count, zonal_length, max_growth, ier, nlev &
) bind(C, name="dbaroc_growth_rate_1d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: u, theta, pressure, temperature
    real(c_double), value, intent(in) :: f_cor, beta, zonal_length
    integer(c_int), value, intent(in) :: smooth_window, wavenumber_mode, wavenumber_count, nlev
    real(c_double), intent(out) :: max_growth
    integer(c_int), intent(out) :: ier

    integer :: nlev_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f
    real(real64), pointer :: u_view(:), theta_view(:), pressure_view(:), temperature_view(:)
    real(real64) :: max_growth_impl
    integer :: ier_impl

    nlev_f = int(nlev)
    smooth_window_f = int(smooth_window)
    wavenumber_mode_f = int(wavenumber_mode)
    wavenumber_count_f = int(wavenumber_count)

    call c_f_pointer(u, u_view, [nlev_f])
    call c_f_pointer(theta, theta_view, [nlev_f])
    call c_f_pointer(pressure, pressure_view, [nlev_f])
    call c_f_pointer(temperature, temperature_view, [nlev_f])

    call dbaroc_growth_rate_1d( &
        u_view, theta_view, pressure_view, temperature_view, &
        f_cor, beta, smooth_window_f, wavenumber_mode_f, wavenumber_count_f, zonal_length, &
        max_growth_impl, ier_impl, nlev_f &
    )
    max_growth = max_growth_impl
    ier = ier_impl
end subroutine dbaroc_growth_rate_1d_c


subroutine dbaroc_growth_rate_profiles_c( &
    u_input, temperature_input, source_pressure, target_pressure, f_cor, beta, &
    interp_kind, smooth_window, wavenumber_mode, wavenumber_count, zonal_length, &
    growth, ier, nlev_in, nprofile, nlev_out &
) bind(C, name="dbaroc_growth_rate_profiles_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use growth_rate_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: u_input, temperature_input, source_pressure, target_pressure
    type(c_ptr), value, intent(in) :: f_cor, beta, zonal_length, growth, ier
    integer(c_int), value, intent(in) :: interp_kind, smooth_window, wavenumber_mode, wavenumber_count
    integer(c_int), value, intent(in) :: nlev_in, nprofile, nlev_out

    integer :: nlev_in_f, nprofile_f, nlev_out_f
    integer :: interp_kind_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f
    real(real64), pointer :: u_view(:, :), temperature_view(:, :), source_pressure_view(:)
    real(real64), pointer :: target_pressure_view(:, :), f_cor_view(:), beta_view(:), zonal_length_view(:)
    real(real64), pointer :: growth_view(:)
    integer(c_int), pointer :: ier_view(:)

    nlev_in_f = int(nlev_in)
    nprofile_f = int(nprofile)
    nlev_out_f = int(nlev_out)
    interp_kind_f = int(interp_kind)
    smooth_window_f = int(smooth_window)
    wavenumber_mode_f = int(wavenumber_mode)
    wavenumber_count_f = int(wavenumber_count)

    call c_f_pointer(u_input, u_view, [nlev_in_f, nprofile_f])
    call c_f_pointer(temperature_input, temperature_view, [nlev_in_f, nprofile_f])
    call c_f_pointer(source_pressure, source_pressure_view, [nlev_in_f])
    call c_f_pointer(target_pressure, target_pressure_view, [nlev_out_f, nprofile_f])
    call c_f_pointer(f_cor, f_cor_view, [nprofile_f])
    call c_f_pointer(beta, beta_view, [nprofile_f])
    call c_f_pointer(zonal_length, zonal_length_view, [nprofile_f])
    call c_f_pointer(growth, growth_view, [nprofile_f])
    call c_f_pointer(ier, ier_view, [nprofile_f])

    call dbaroc_growth_rate_profiles( &
        u_view, temperature_view, source_pressure_view, target_pressure_view, f_cor_view, beta_view, &
        interp_kind_f, smooth_window_f, wavenumber_mode_f, wavenumber_count_f, zonal_length_view, &
        growth_view, ier_view, nlev_in_f, nprofile_f, nlev_out_f &
    )
end subroutine dbaroc_growth_rate_profiles_c
