! =============================================================================
! File    : rcm2points_bindings.f90
! Purpose : C interoperability bindings split from rcm2points.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine drcm2points(ngrd, nyi, nxi, yi, xi, fi, nxyo, yo, xo, fo, xmsg, opt, ncrit, kval, ier) bind(C, name="drcm2points")
    use iso_c_binding
    implicit none

    integer(c_int), value :: ngrd, nyi, nxi, nxyo, opt, ncrit, kval
    integer(c_int) :: ier
    real(c_double), value :: xmsg
    type(c_ptr), value :: yi, xi, fi, yo, xo, fo
    real(c_double), pointer :: yi_f(:,:), xi_f(:,:), fi_f(:,:,:)
    real(c_double), pointer :: yo_f(:), xo_f(:), fo_f(:,:)

    call c_f_pointer(yi, yi_f, [nxi, nyi])
    call c_f_pointer(xi, xi_f, [nxi, nyi])
    call c_f_pointer(fi, fi_f, [nxi, nyi, ngrd])
    call c_f_pointer(yo, yo_f, [nxyo])
    call c_f_pointer(xo, xo_f, [nxyo])
    call c_f_pointer(fo, fo_f, [nxyo, ngrd])

    call drcm2points_impl(ngrd, nyi, nxi, yi_f, xi_f, fi_f, nxyo, yo_f, xo_f, fo_f, xmsg, opt, ncrit, kval, ier)
end subroutine drcm2points
