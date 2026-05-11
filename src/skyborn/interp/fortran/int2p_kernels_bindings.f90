! =============================================================================
! File    : int2p_kernels_bindings.f90
! Purpose : C interoperability bindings split from int2p_kernels.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine dinterp_pressure_1d(ppin, xxin, ppout, xxout, linlog, xmsg, ier, npin, npout) &
    bind(C, name="dinterp_pressure_1d")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use int2p_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: ppin, xxin, ppout, xxout, ier
    integer(c_int), value, intent(in) :: linlog, npin, npout
    real(c_double), value, intent(in) :: xmsg

    real(real64), pointer :: ppin_view(:), xxin_view(:), ppout_view(:), xxout_view(:)
    integer(c_int), pointer :: ier_view
    integer :: npin_f, npout_f, linlog_f, ier_f

    npin_f = int(npin)
    npout_f = int(npout)
    linlog_f = int(linlog)

    call c_f_pointer(ppin, ppin_view, [npin_f])
    call c_f_pointer(xxin, xxin_view, [npin_f])
    call c_f_pointer(ppout, ppout_view, [npout_f])
    call c_f_pointer(xxout, xxout_view, [npout_f])
    call c_f_pointer(ier, ier_view)

    call dinterp_pressure_1d_impl( &
        ppin_view, xxin_view, ppout_view, xxout_view, linlog_f, xmsg, ier_f, npin_f, npout_f &
    )
    ier_view = ier_f
end subroutine dinterp_pressure_1d
