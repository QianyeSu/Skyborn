! =============================================================================
! File    : tropical_cyclone_potential_intensity_bindings.f90
! Purpose : C interoperability bindings split from tropical_cyclone_potential_intensity.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine calculate_pi_gridded_data_c(sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
                                       min_pressure, max_wind, error_flag, nlat, nlon, num_levels) &
                                       bind(C, name="calculate_pi_gridded_data_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    integer(c_int), value, intent(in) :: nlat, nlon, num_levels

    real(c_float), pointer :: sst_view(:, :), psl_view(:, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :), mixing_ratio_view(:, :, :)
    real(c_float), pointer :: min_pressure_view(:, :), max_wind_view(:, :)
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)

    call c_f_pointer(sst_in, sst_view, [nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)

    call calculate_pi_gridded_data( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, min_pressure_view, max_wind_view, error_flag_view &
    )
end subroutine calculate_pi_gridded_data_c

subroutine calculate_pi_gridded_with_missing_c(sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
                                               min_pressure, max_wind, error_flag, nlat, nlon, num_levels) &
                                               bind(C, name="calculate_pi_gridded_with_missing_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    integer(c_int), value, intent(in) :: nlat, nlon, num_levels

    real(c_float), pointer :: sst_view(:, :), psl_view(:, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :), mixing_ratio_view(:, :, :)
    real(c_float), pointer :: min_pressure_view(:, :), max_wind_view(:, :)
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)

    call c_f_pointer(sst_in, sst_view, [nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)

    call calculate_pi_gridded_with_missing( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, min_pressure_view, max_wind_view, error_flag_view &
    )
end subroutine calculate_pi_gridded_with_missing_c

subroutine calculate_pi_gridded_diagnostics_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd, &
    outflow_source_flag, ckcd_in, nlat, nlon, num_levels &
) bind(C, name="calculate_pi_gridded_diagnostics_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    type(c_ptr), value, intent(in) :: outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd
    integer(c_int), value, intent(in) :: outflow_source_flag, nlat, nlon, num_levels
    real(c_float), value, intent(in) :: ckcd_in

    real(c_float), pointer :: sst_view(:, :), psl_view(:, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :), mixing_ratio_view(:, :, :)
    real(c_float), pointer :: min_pressure_view(:, :), max_wind_view(:, :)
    real(c_float), pointer :: outflow_temp_view(:, :), outflow_level_view(:, :)
    real(c_float), pointer :: lnpi_view(:, :), lneff_view(:, :), lndiseq_view(:, :), lnckcd_view
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, outflow_flag_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    outflow_flag_f = int(outflow_source_flag)

    call c_f_pointer(sst_in, sst_view, [nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)
    call c_f_pointer(outflow_temp, outflow_temp_view, [nlat_f, nlon_f])
    call c_f_pointer(outflow_level, outflow_level_view, [nlat_f, nlon_f])
    call c_f_pointer(lnpi, lnpi_view, [nlat_f, nlon_f])
    call c_f_pointer(lneff, lneff_view, [nlat_f, nlon_f])
    call c_f_pointer(lndiseq, lndiseq_view, [nlat_f, nlon_f])
    call c_f_pointer(lnckcd, lnckcd_view)

    call calculate_pi_gridded_diagnostics( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, min_pressure_view, max_wind_view, error_flag_view, &
        outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view, &
        outflow_flag_f, ckcd_in &
    )
end subroutine calculate_pi_gridded_diagnostics_c

subroutine calculate_pi_gridded_diagnostics_with_missing_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd, &
    outflow_source_flag, ckcd_in, nlat, nlon, num_levels &
) bind(C, name="calculate_pi_gridded_diagnostics_with_missing_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    type(c_ptr), value, intent(in) :: outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd
    integer(c_int), value, intent(in) :: outflow_source_flag, nlat, nlon, num_levels
    real(c_float), value, intent(in) :: ckcd_in

    real(c_float), pointer :: sst_view(:, :), psl_view(:, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :), mixing_ratio_view(:, :, :)
    real(c_float), pointer :: min_pressure_view(:, :), max_wind_view(:, :)
    real(c_float), pointer :: outflow_temp_view(:, :), outflow_level_view(:, :)
    real(c_float), pointer :: lnpi_view(:, :), lneff_view(:, :), lndiseq_view(:, :), lnckcd_view
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, outflow_flag_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    outflow_flag_f = int(outflow_source_flag)

    call c_f_pointer(sst_in, sst_view, [nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)
    call c_f_pointer(outflow_temp, outflow_temp_view, [nlat_f, nlon_f])
    call c_f_pointer(outflow_level, outflow_level_view, [nlat_f, nlon_f])
    call c_f_pointer(lnpi, lnpi_view, [nlat_f, nlon_f])
    call c_f_pointer(lneff, lneff_view, [nlat_f, nlon_f])
    call c_f_pointer(lndiseq, lndiseq_view, [nlat_f, nlon_f])
    call c_f_pointer(lnckcd, lnckcd_view)

    call calculate_pi_gridded_diagnostics_with_missing( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, min_pressure_view, max_wind_view, error_flag_view, &
        outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view, &
        outflow_flag_f, ckcd_in &
    )
end subroutine calculate_pi_gridded_diagnostics_with_missing_c

subroutine calculate_pi_single_profile_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, actual_levels, &
    min_pressure, max_wind, error_flag, num_levels &
) bind(C, name="calculate_pi_single_profile_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    real(c_float), value, intent(in) :: sst_in, psl_in
    type(c_ptr), value, intent(in) :: pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    integer(c_int), value, intent(in) :: actual_levels, num_levels

    real(c_float), pointer :: pressure_view(:), temp_view(:), mixing_ratio_view(:)
    real(c_float), pointer :: min_pressure_view, max_wind_view
    integer(c_int), pointer :: error_flag_view
    integer :: actual_levels_f, num_levels_f

    actual_levels_f = int(actual_levels)
    num_levels_f = int(num_levels)

    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f])
    call c_f_pointer(min_pressure, min_pressure_view)
    call c_f_pointer(max_wind, max_wind_view)
    call c_f_pointer(error_flag, error_flag_view)

    call calculate_pi_single_profile( &
        sst_in, psl_in, pressure_view, temp_view, mixing_ratio_view, &
        num_levels_f, actual_levels_f, min_pressure_view, max_wind_view, error_flag_view &
    )
end subroutine calculate_pi_single_profile_c

subroutine calculate_pi_profile_diagnostics_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, actual_levels, &
    min_pressure, max_wind, error_flag, outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd, &
    outflow_source_flag, ckcd_in, num_levels &
) bind(C, name="calculate_pi_profile_diagnostics_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    real(c_float), value, intent(in) :: sst_in, psl_in, ckcd_in
    type(c_ptr), value, intent(in) :: pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    type(c_ptr), value, intent(in) :: outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd
    integer(c_int), value, intent(in) :: actual_levels, outflow_source_flag, num_levels

    real(c_float), pointer :: pressure_view(:), temp_view(:), mixing_ratio_view(:)
    real(c_float), pointer :: min_pressure_view, max_wind_view
    real(c_float), pointer :: outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view
    integer(c_int), pointer :: error_flag_view
    integer :: actual_levels_f, outflow_flag_f, num_levels_f

    actual_levels_f = int(actual_levels)
    outflow_flag_f = int(outflow_source_flag)
    num_levels_f = int(num_levels)

    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_levels_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_levels_f])
    call c_f_pointer(min_pressure, min_pressure_view)
    call c_f_pointer(max_wind, max_wind_view)
    call c_f_pointer(error_flag, error_flag_view)
    call c_f_pointer(outflow_temp, outflow_temp_view)
    call c_f_pointer(outflow_level, outflow_level_view)
    call c_f_pointer(lnpi, lnpi_view)
    call c_f_pointer(lneff, lneff_view)
    call c_f_pointer(lndiseq, lndiseq_view)
    call c_f_pointer(lnckcd, lnckcd_view)

    call calculate_pi_profile_diagnostics( &
        sst_in, psl_in, pressure_view, temp_view, mixing_ratio_view, &
        num_levels_f, actual_levels_f, min_pressure_view, max_wind_view, error_flag_view, &
        outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view, &
        outflow_flag_f, ckcd_in &
    )
end subroutine calculate_pi_profile_diagnostics_c

subroutine calculate_pi_4d_data_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, nlat, nlon, num_levels, num_times &
) bind(C, name="calculate_pi_4d_data_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    integer(c_int), value, intent(in) :: nlat, nlon, num_levels, num_times

    real(c_float), pointer :: sst_view(:, :, :), psl_view(:, :, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :, :), mixing_ratio_view(:, :, :, :)
    real(c_float), pointer :: min_pressure_view(:, :, :), max_wind_view(:, :, :)
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, num_times_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    num_times_f = int(num_times)

    call c_f_pointer(sst_in, sst_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)

    call calculate_pi_4d_data( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, num_times_f, min_pressure_view, max_wind_view, error_flag_view &
    )
end subroutine calculate_pi_4d_data_c

subroutine calculate_pi_4d_with_missing_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, nlat, nlon, num_levels, num_times &
) bind(C, name="calculate_pi_4d_with_missing_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    integer(c_int), value, intent(in) :: nlat, nlon, num_levels, num_times

    real(c_float), pointer :: sst_view(:, :, :), psl_view(:, :, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :, :), mixing_ratio_view(:, :, :, :)
    real(c_float), pointer :: min_pressure_view(:, :, :), max_wind_view(:, :, :)
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, num_times_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    num_times_f = int(num_times)

    call c_f_pointer(sst_in, sst_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)

    call calculate_pi_4d_with_missing( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, num_times_f, min_pressure_view, max_wind_view, error_flag_view &
    )
end subroutine calculate_pi_4d_with_missing_c

subroutine calculate_pi_4d_diagnostics_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd, &
    outflow_source_flag, ckcd_in, nlat, nlon, num_levels, num_times &
) bind(C, name="calculate_pi_4d_diagnostics_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    type(c_ptr), value, intent(in) :: outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd
    integer(c_int), value, intent(in) :: outflow_source_flag, nlat, nlon, num_levels, num_times
    real(c_float), value, intent(in) :: ckcd_in

    real(c_float), pointer :: sst_view(:, :, :), psl_view(:, :, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :, :), mixing_ratio_view(:, :, :, :)
    real(c_float), pointer :: min_pressure_view(:, :, :), max_wind_view(:, :, :)
    real(c_float), pointer :: outflow_temp_view(:, :, :), outflow_level_view(:, :, :)
    real(c_float), pointer :: lnpi_view(:, :, :), lneff_view(:, :, :), lndiseq_view(:, :, :), lnckcd_view
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, num_times_f, outflow_flag_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    num_times_f = int(num_times)
    outflow_flag_f = int(outflow_source_flag)

    call c_f_pointer(sst_in, sst_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)
    call c_f_pointer(outflow_temp, outflow_temp_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(outflow_level, outflow_level_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lnpi, lnpi_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lneff, lneff_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lndiseq, lndiseq_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lnckcd, lnckcd_view)

    call calculate_pi_4d_diagnostics( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, num_times_f, min_pressure_view, max_wind_view, error_flag_view, &
        outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view, &
        outflow_flag_f, ckcd_in &
    )
end subroutine calculate_pi_4d_diagnostics_c

subroutine calculate_pi_4d_diagnostics_with_missing_c( &
    sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in, &
    min_pressure, max_wind, error_flag, outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd, &
    outflow_source_flag, ckcd_in, nlat, nlon, num_levels, num_times &
) bind(C, name="calculate_pi_4d_diagnostics_with_missing_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: sst_in, psl_in, pressure_levels, temp_in, mixing_ratio_in
    type(c_ptr), value, intent(in) :: min_pressure, max_wind, error_flag
    type(c_ptr), value, intent(in) :: outflow_temp, outflow_level, lnpi, lneff, lndiseq, lnckcd
    integer(c_int), value, intent(in) :: outflow_source_flag, nlat, nlon, num_levels, num_times
    real(c_float), value, intent(in) :: ckcd_in

    real(c_float), pointer :: sst_view(:, :, :), psl_view(:, :, :), pressure_view(:)
    real(c_float), pointer :: temp_view(:, :, :, :), mixing_ratio_view(:, :, :, :)
    real(c_float), pointer :: min_pressure_view(:, :, :), max_wind_view(:, :, :)
    real(c_float), pointer :: outflow_temp_view(:, :, :), outflow_level_view(:, :, :)
    real(c_float), pointer :: lnpi_view(:, :, :), lneff_view(:, :, :), lndiseq_view(:, :, :), lnckcd_view
    integer(c_int), pointer :: error_flag_view
    integer :: nlat_f, nlon_f, num_levels_f, num_times_f, outflow_flag_f

    nlat_f = int(nlat)
    nlon_f = int(nlon)
    num_levels_f = int(num_levels)
    num_times_f = int(num_times)
    outflow_flag_f = int(outflow_source_flag)

    call c_f_pointer(sst_in, sst_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(psl_in, psl_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(pressure_levels, pressure_view, [num_levels_f])
    call c_f_pointer(temp_in, temp_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(mixing_ratio_in, mixing_ratio_view, [num_times_f, num_levels_f, nlat_f, nlon_f])
    call c_f_pointer(min_pressure, min_pressure_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(max_wind, max_wind_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(error_flag, error_flag_view)
    call c_f_pointer(outflow_temp, outflow_temp_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(outflow_level, outflow_level_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lnpi, lnpi_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lneff, lneff_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lndiseq, lndiseq_view, [num_times_f, nlat_f, nlon_f])
    call c_f_pointer(lnckcd, lnckcd_view)

    call calculate_pi_4d_diagnostics_with_missing( &
        sst_view, psl_view, pressure_view, temp_view, mixing_ratio_view, &
        nlat_f, nlon_f, num_levels_f, num_times_f, min_pressure_view, max_wind_view, error_flag_view, &
        outflow_temp_view, outflow_level_view, lnpi_view, lneff_view, lndiseq_view, lnckcd_view, &
        outflow_flag_f, ckcd_in &
    )
end subroutine calculate_pi_4d_diagnostics_with_missing_c

subroutine cape_c(parcel_temp, parcel_mixing_ratio, parcel_pressure, temp_profile, mixing_ratio_profile, &
                  pressure_profile, num_points, buoyancy_param, cape_value, outflow_temp, outflow_level, &
                  error_flag, array_size) bind(C, name="cape_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_float, c_int, c_ptr
    implicit none

    real(c_float), value, intent(in) :: parcel_temp, parcel_mixing_ratio, parcel_pressure, buoyancy_param
    type(c_ptr), value, intent(in) :: temp_profile, mixing_ratio_profile, pressure_profile
    type(c_ptr), value, intent(in) :: cape_value, outflow_temp, outflow_level, error_flag
    integer(c_int), value, intent(in) :: array_size, num_points

    real(c_float), pointer :: temp_profile_view(:), mixing_ratio_profile_view(:), pressure_profile_view(:)
    real(c_float), pointer :: cape_value_view, outflow_temp_view, outflow_level_view
    integer(c_int), pointer :: error_flag_view
    integer :: array_size_f, num_points_f

    array_size_f = int(array_size)
    num_points_f = int(num_points)

    call c_f_pointer(temp_profile, temp_profile_view, [array_size_f])
    call c_f_pointer(mixing_ratio_profile, mixing_ratio_profile_view, [array_size_f])
    call c_f_pointer(pressure_profile, pressure_profile_view, [array_size_f])
    call c_f_pointer(cape_value, cape_value_view)
    call c_f_pointer(outflow_temp, outflow_temp_view)
    call c_f_pointer(outflow_level, outflow_level_view)
    call c_f_pointer(error_flag, error_flag_view)

    call CAPE( &
        parcel_temp, parcel_mixing_ratio, parcel_pressure, &
        temp_profile_view, mixing_ratio_profile_view, pressure_profile_view, &
        array_size_f, num_points_f, buoyancy_param, &
        cape_value_view, outflow_temp_view, outflow_level_view, error_flag_view &
    )
end subroutine cape_c
