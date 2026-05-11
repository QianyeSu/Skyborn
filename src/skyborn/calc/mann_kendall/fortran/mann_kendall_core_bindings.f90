! =============================================================================
! File    : mann_kendall_core_bindings.f90
! Purpose : C interoperability bindings split from mann_kendall_core.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine mk_score_var_batch_c(data, s_values, var_values, modified, ntime, nseries) bind(C, name="mk_score_var_batch_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, s_values, var_values
    integer(c_int), value, intent(in) :: modified, ntime, nseries

    integer :: ntime_f, nseries_f, modified_f
    real(real64), pointer :: data_view(:, :), s_values_view(:), var_values_view(:)

    ntime_f = int(ntime)
    nseries_f = int(nseries)
    modified_f = int(modified)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])

    call mk_score_var_batch(data_view, s_values_view, var_values_view, modified_f, ntime_f, nseries_f)
end subroutine mk_score_var_batch_c


subroutine sen_slope_batch_c(data, slopes, ntime, nseries) bind(C, name="sen_slope_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, slopes
    integer(c_int), value, intent(in) :: ntime, nseries

    integer :: ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), slopes_view(:)

    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call sen_slope_batch(data_view, slopes_view, ntime_f, nseries_f)
end subroutine sen_slope_batch_c


subroutine grouped_sen_slope_batch_c(data, slopes, period, ntime, nseries) bind(C, name="grouped_sen_slope_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, slopes
    integer(c_int), value, intent(in) :: period, ntime, nseries

    integer :: period_f, ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), slopes_view(:)

    period_f = int(period)
    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call grouped_sen_slope_batch(data_view, slopes_view, period_f, ntime_f, nseries_f)
end subroutine grouped_sen_slope_batch_c


subroutine grouped_correlated_stats_batch_c( &
    data, s_values, var_values, denom, period, ntime, nseries &
) bind(C, name="grouped_correlated_stats_batch_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, s_values, var_values, denom
    integer(c_int), value, intent(in) :: period, ntime, nseries

    integer :: period_f, ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), s_values_view(:), var_values_view(:), denom_view

    period_f = int(period)
    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(denom, denom_view)

    call grouped_correlated_stats_batch( &
        data_view, s_values_view, var_values_view, denom_view, period_f, ntime_f, nseries_f &
    )
end subroutine grouped_correlated_stats_batch_c


subroutine partial_stats_batch_c( &
    response, covariate, s_values, var_values, tau_values, ntime, nseries &
) bind(C, name="partial_stats_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: response, covariate, s_values, var_values, tau_values
    integer(c_int), value, intent(in) :: ntime, nseries

    integer :: ntime_f, nseries_f
    real(real64), pointer :: response_view(:, :), covariate_view(:, :)
    real(real64), pointer :: s_values_view(:), var_values_view(:), tau_values_view(:)

    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(response, response_view, [ntime_f, nseries_f])
    call c_f_pointer(covariate, covariate_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(tau_values, tau_values_view, [nseries_f])

    call partial_stats_batch( &
        response_view, covariate_view, s_values_view, var_values_view, tau_values_view, ntime_f, nseries_f &
    )
end subroutine partial_stats_batch_c


subroutine partial_stats_sen_batch_c( &
    response, covariate, s_values, var_values, tau_values, slopes, ntime, nseries &
) bind(C, name="partial_stats_sen_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: response, covariate, s_values, var_values, tau_values, slopes
    integer(c_int), value, intent(in) :: ntime, nseries

    integer :: ntime_f, nseries_f
    real(real64), pointer :: response_view(:, :), covariate_view(:, :)
    real(real64), pointer :: s_values_view(:), var_values_view(:), tau_values_view(:), slopes_view(:)

    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(response, response_view, [ntime_f, nseries_f])
    call c_f_pointer(covariate, covariate_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(tau_values, tau_values_view, [nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call partial_stats_sen_batch( &
        response_view, covariate_view, s_values_view, var_values_view, tau_values_view, slopes_view, ntime_f, nseries_f &
    )
end subroutine partial_stats_sen_batch_c


subroutine mk_score_var_sen_batch_c( &
    data, s_values, var_values, slopes, modified, ntime, nseries &
) bind(C, name="mk_score_var_sen_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, s_values, var_values, slopes
    integer(c_int), value, intent(in) :: modified, ntime, nseries

    integer :: modified_f, ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), s_values_view(:), var_values_view(:), slopes_view(:)

    modified_f = int(modified)
    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call mk_score_var_sen_batch( &
        data_view, s_values_view, var_values_view, slopes_view, modified_f, ntime_f, nseries_f &
    )
end subroutine mk_score_var_sen_batch_c


subroutine mk_yue_wang_score_var_sen_batch_c( &
    data, s_values, var_values, slopes, lag, ntime, nseries &
) bind(C, name="mk_yue_wang_score_var_sen_batch_c")
    use, intrinsic :: iso_c_binding, only : c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, s_values, var_values, slopes
    integer(c_int), value, intent(in) :: lag, ntime, nseries

    integer :: lag_f, ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), s_values_view(:), var_values_view(:), slopes_view(:)

    lag_f = int(lag)
    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call mk_yue_wang_score_var_sen_batch( &
        data_view, s_values_view, var_values_view, slopes_view, lag_f, ntime_f, nseries_f &
    )
end subroutine mk_yue_wang_score_var_sen_batch_c


subroutine mk_hamed_rao_score_var_sen_batch_c( &
    data, s_values, var_values, slopes, interval, lag, ntime, nseries &
) bind(C, name="mk_hamed_rao_score_var_sen_batch_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use mann_kendall_core_mod, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data, s_values, var_values, slopes
    real(c_double), value, intent(in) :: interval
    integer(c_int), value, intent(in) :: lag, ntime, nseries

    integer :: lag_f, ntime_f, nseries_f
    real(real64), pointer :: data_view(:, :), s_values_view(:), var_values_view(:), slopes_view(:)

    lag_f = int(lag)
    ntime_f = int(ntime)
    nseries_f = int(nseries)

    call c_f_pointer(data, data_view, [ntime_f, nseries_f])
    call c_f_pointer(s_values, s_values_view, [nseries_f])
    call c_f_pointer(var_values, var_values_view, [nseries_f])
    call c_f_pointer(slopes, slopes_view, [nseries_f])

    call mk_hamed_rao_score_var_sen_batch( &
        data_view, s_values_view, var_values_view, slopes_view, interval, lag_f, ntime_f, nseries_f &
    )
end subroutine mk_hamed_rao_score_var_sen_batch_c
