! =============================================================================
! File    : rcm2rgrid_bindings.f90
! Purpose : C interoperability bindings split from rcm2rgrid.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine drcm2rgrid(ngrd, nyi, nxi, yi, xi, fi, nyo, yo, nxo, xo, fo, xmsg, ncrit, opt, ier) &
        bind(C, name="drcm2rgrid")
    use iso_c_binding
    implicit none

    integer(c_int), value :: ngrd, nyi, nxi, nyo, nxo, ncrit, opt
    integer(c_int) :: ier
    real(c_double), value :: xmsg
    type(c_ptr), value :: yi, xi, fi, yo, xo, fo
    real(c_double), pointer :: yi_f(:,:), xi_f(:,:), fi_f(:,:,:)
    real(c_double), pointer :: yo_f(:), xo_f(:), fo_f(:,:,:)

    call c_f_pointer(yi, yi_f, [nxi, nyi])
    call c_f_pointer(xi, xi_f, [nxi, nyi])
    call c_f_pointer(fi, fi_f, [nxi, nyi, ngrd])
    call c_f_pointer(yo, yo_f, [nyo])
    call c_f_pointer(xo, xo_f, [nxo])
    call c_f_pointer(fo, fo_f, [nxo, nyo, ngrd])

    call drcm2rgrid_impl(ngrd, nyi, nxi, yi_f, xi_f, fi_f, nyo, yo_f, nxo, xo_f, fo_f, xmsg, ncrit, opt, ier)
end subroutine drcm2rgrid


subroutine drgrid2rcm(ngrd, nyi, nxi, yi, xi, fi, nyo, nxo, yo, xo, fo, xmsg, ncrit, opt, ier) &
        bind(C, name="drgrid2rcm")
    use iso_c_binding
    implicit none

    integer(c_int), value :: ngrd, nyi, nxi, nyo, nxo, ncrit, opt
    integer(c_int) :: ier
    real(c_double), value :: xmsg
    type(c_ptr), value :: yi, xi, fi, yo, xo, fo
    real(c_double), pointer :: yi_f(:), xi_f(:), fi_f(:,:,:)
    real(c_double), pointer :: yo_f(:,:), xo_f(:,:), fo_f(:,:,:)

    call c_f_pointer(yi, yi_f, [nyi])
    call c_f_pointer(xi, xi_f, [nxi])
    call c_f_pointer(fi, fi_f, [nxi, nyi, ngrd])
    call c_f_pointer(yo, yo_f, [nxo, nyo])
    call c_f_pointer(xo, xo_f, [nxo, nyo])
    call c_f_pointer(fo, fo_f, [nxo, nyo, ngrd])

    call drgrid2rcm_impl(ngrd, nyi, nxi, yi_f, xi_f, fi_f, nyo, nxo, yo_f, xo_f, fo_f, xmsg, ncrit, opt, ier)
end subroutine drgrid2rcm
