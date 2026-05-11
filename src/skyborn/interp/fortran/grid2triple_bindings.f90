! =============================================================================
! File    : grid2triple_bindings.f90
! Purpose : C interoperability bindings split from grid2triple.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine grid2triple(x, y, z, d, ld, zmsg, ier, mx, ny, ldmax) bind(C, name="grid2triple")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    implicit none

    type(c_ptr), value, intent(in) :: x, y, z, d, ld, ier
    integer(c_int), value, intent(in) :: mx, ny, ldmax
    real(c_double), value, intent(in) :: zmsg

    real(c_double), pointer :: x_view(:), y_view(:), z_view(:, :)
    real(c_double), pointer :: d_view(:, :)
    integer(c_int), pointer :: ld_view, ier_view
    integer :: mx_f, ny_f, ldmax_f, ld_f, ier_f

    mx_f = int(mx)
    ny_f = int(ny)
    ldmax_f = int(ldmax)

    call c_f_pointer(x, x_view, [mx_f])
    call c_f_pointer(y, y_view, [ny_f])
    call c_f_pointer(z, z_view, [mx_f, ny_f])
    call c_f_pointer(d, d_view, [ldmax_f, 3])
    call c_f_pointer(ld, ld_view)
    call c_f_pointer(ier, ier_view)

    call grid2triple_impl(x_view, y_view, z_view, mx_f, ny_f, d_view, ldmax_f, ld_f, zmsg, ier_f)
    ld_view = ld_f
    ier_view = ier_f
end subroutine grid2triple
