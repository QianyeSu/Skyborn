subroutine liang_single_c( &
    y1, y2, t21, tau21, z, dh1_star, dh1_noise, nm, npt, ierr &
) bind(C, name="liang_single_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use causality_core_mod, only : liang_single
    implicit none

    type(c_ptr), value, intent(in) :: y1, y2
    real(c_double), intent(out) :: t21, tau21, z, dh1_star, dh1_noise
    integer(c_int), value, intent(in) :: nm, npt
    integer(c_int), intent(out) :: ierr
    integer :: ierr_f

    real(c_double), pointer :: y1_view(:), y2_view(:)

    call c_f_pointer(y1, y1_view, [int(nm)])
    call c_f_pointer(y2, y2_view, [int(nm)])
    call liang_single( &
        y1_view, y2_view, int(nm), int(npt), &
        t21, tau21, z, dh1_star, dh1_noise, ierr_f &
    )
    ierr = int(ierr_f, kind=c_int)
end subroutine liang_single_c


subroutine liang_batch_c(y1, y2, t21, tau21, nm, nsim, npt, ierr) &
    bind(C, name="liang_batch_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use causality_core_mod, only : liang_batch
    implicit none

    type(c_ptr), value, intent(in) :: y1, y2, t21, tau21
    integer(c_int), value, intent(in) :: nm, nsim, npt
    integer(c_int), intent(out) :: ierr
    integer :: ierr_f

    real(c_double), pointer :: y1_view(:, :), y2_view(:, :)
    real(c_double), pointer :: t21_view(:), tau21_view(:)

    call c_f_pointer(y1, y1_view, [int(nm), int(nsim)])
    call c_f_pointer(y2, y2_view, [int(nm), int(nsim)])
    call c_f_pointer(t21, t21_view, [int(nsim)])
    call c_f_pointer(tau21, tau21_view, [int(nsim)])
    call liang_batch( &
        y1_view, y2_view, int(nm), int(nsim), int(npt), &
        t21_view, tau21_view, ierr_f &
    )
    ierr = int(ierr_f, kind=c_int)
end subroutine liang_batch_c
