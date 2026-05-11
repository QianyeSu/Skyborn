! =============================================================================
! File    : triple2grid_bindings.f90
! Purpose : C interoperability bindings split from triple2grid.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine triple2grid1(kz, xi, yi, zi, zmsg, mx, ny, gx, gy, grid, domain, loop, method, distmx, mx2, ny2, &
                        x, y, z, gbigx, gbigy, gbigxy, ier) bind(C, name="triple2grid1")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    implicit none

    integer(c_int), value, intent(in) :: kz, mx, ny, loop, method, mx2, ny2
    type(c_ptr), value, intent(in) :: xi, yi, zi, gx, gy, grid, x, y, z, gbigx, gbigy, gbigxy, ier
    real(c_double), value, intent(in) :: zmsg, domain, distmx

    real(c_double), pointer :: xi_view(:), yi_view(:), zi_view(:), gx_view(:), gy_view(:)
    real(c_double), pointer :: grid_view(:, :)
    real(c_double), pointer :: x_view(:), y_view(:), z_view(:), gbigx_view(:), gbigy_view(:), gbigxy_view(:, :)
    integer(c_int), pointer :: ier_view
    integer :: kz_f, mx_f, ny_f, loop_f, method_f, mx2_f, ny2_f, ier_f

    kz_f = int(kz)
    mx_f = int(mx)
    ny_f = int(ny)
    loop_f = int(loop)
    method_f = int(method)
    mx2_f = int(mx2)
    ny2_f = int(ny2)

    call c_f_pointer(xi, xi_view, [kz_f])
    call c_f_pointer(yi, yi_view, [kz_f])
    call c_f_pointer(zi, zi_view, [kz_f])
    call c_f_pointer(gx, gx_view, [mx_f])
    call c_f_pointer(gy, gy_view, [ny_f])
    call c_f_pointer(grid, grid_view, [mx_f, ny_f])
    call c_f_pointer(x, x_view, [kz_f])
    call c_f_pointer(y, y_view, [kz_f])
    call c_f_pointer(z, z_view, [kz_f])
    call c_f_pointer(gbigx, gbigx_view, [mx2_f])
    call c_f_pointer(gbigy, gbigy_view, [ny2_f])
    call c_f_pointer(gbigxy, gbigxy_view, [mx2_f, ny2_f])
    call c_f_pointer(ier, ier_view)

    call triple2grid1_impl( &
        kz_f, xi_view, yi_view, zi_view, zmsg, mx_f, ny_f, gx_view, gy_view, grid_view, &
        domain, loop_f, method_f, distmx, mx2_f, ny2_f, x_view, y_view, z_view, gbigx_view, gbigy_view, &
        gbigxy_view, ier_f &
    )
    ier_view = ier_f
end subroutine triple2grid1


subroutine trip2grd2(kz, x, y, z, zmsg, mx, ny, gxout, gyout, gout, mflag, nflag, method, ddcrit, ier) &
    bind(C, name="trip2grd2")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    implicit none

    integer(c_int), value, intent(in) :: kz, mx, ny, mflag, nflag, method
    type(c_ptr), value, intent(in) :: x, y, z, gxout, gyout, gout, ier
    real(c_double), value, intent(in) :: zmsg, ddcrit

    real(c_double), pointer :: x_view(:), y_view(:), z_view(:), gxout_view(:), gyout_view(:), gout_view(:, :)
    integer(c_int), pointer :: ier_view
    integer :: kz_f, mx_f, ny_f, mflag_f, nflag_f, method_f, ier_f

    kz_f = int(kz)
    mx_f = int(mx)
    ny_f = int(ny)
    mflag_f = int(mflag)
    nflag_f = int(nflag)
    method_f = int(method)

    call c_f_pointer(x, x_view, [kz_f])
    call c_f_pointer(y, y_view, [kz_f])
    call c_f_pointer(z, z_view, [kz_f])
    call c_f_pointer(gxout, gxout_view, [mx_f])
    call c_f_pointer(gyout, gyout_view, [ny_f])
    call c_f_pointer(gout, gout_view, [mx_f, ny_f])
    call c_f_pointer(ier, ier_view)

    call trip2grd2_impl(kz_f, x_view, y_view, z_view, zmsg, mx_f, ny_f, gxout_view, gyout_view, gout_view, &
                        mflag_f, nflag_f, method_f, ddcrit, ier_f)
    ier_view = ier_f
end subroutine trip2grd2


subroutine trip2grd3(kz, x, y, z, zmsg, mx, ny, gxout, gyout, gout, mflag, nflag, method, ddcrit, ier) &
    bind(C, name="trip2grd3")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    implicit none

    integer(c_int), value, intent(in) :: kz, mx, ny, mflag, nflag, method
    type(c_ptr), value, intent(in) :: x, y, z, gxout, gyout, gout, ier
    real(c_double), value, intent(in) :: zmsg, ddcrit

    real(c_double), pointer :: x_view(:), y_view(:), z_view(:), gxout_view(:), gyout_view(:), gout_view(:, :)
    integer(c_int), pointer :: ier_view
    integer :: kz_f, mx_f, ny_f, mflag_f, nflag_f, method_f, ier_f

    kz_f = int(kz)
    mx_f = int(mx)
    ny_f = int(ny)
    mflag_f = int(mflag)
    nflag_f = int(nflag)
    method_f = int(method)

    call c_f_pointer(x, x_view, [kz_f])
    call c_f_pointer(y, y_view, [kz_f])
    call c_f_pointer(z, z_view, [kz_f])
    call c_f_pointer(gxout, gxout_view, [mx_f])
    call c_f_pointer(gyout, gyout_view, [ny_f])
    call c_f_pointer(gout, gout_view, [mx_f, ny_f])
    call c_f_pointer(ier, ier_view)

    call trip2grd3_impl(kz_f, x_view, y_view, z_view, zmsg, mx_f, ny_f, gxout_view, gyout_view, gout_view, &
                        mflag_f, nflag_f, method_f, ddcrit, ier_f)
    ier_view = ier_f
end subroutine trip2grd3
