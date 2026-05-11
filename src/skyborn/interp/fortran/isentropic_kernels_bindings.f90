! =============================================================================
! File    : isentropic_kernels_bindings.f90
! Purpose : C interoperability bindings split from isentropic_kernels.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine dinterp_to_isentropic_corder_into( &
    data_flat, temp_flat, pressure_flat, theta_levels, output_flat, &
    p0, kappa, spvl, kxtrp, nouter, nlev, ninner, ntheta) &
    bind(C, name="dinterp_to_isentropic_corder_into")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use isentropic_kernels_core, only : real64
    implicit none

    type(c_ptr), value, intent(in) :: data_flat, temp_flat, pressure_flat, theta_levels, output_flat
    real(c_double), value, intent(in) :: p0, kappa, spvl
    integer(c_int), value, intent(in) :: kxtrp, nouter, nlev, ninner, ntheta

    real(real64), pointer :: data_flat_view(:), temp_flat_view(:), pressure_flat_view(:), theta_levels_view(:)
    real(real64), pointer :: output_flat_view(:)
    integer :: nouter_f, nlev_f, ninner_f, ntheta_f, kxtrp_f

    nouter_f = int(nouter)
    nlev_f = int(nlev)
    ninner_f = int(ninner)
    ntheta_f = int(ntheta)
    kxtrp_f = int(kxtrp)

    call c_f_pointer(data_flat, data_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(temp_flat, temp_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(pressure_flat, pressure_flat_view, [nouter_f * nlev_f * ninner_f])
    call c_f_pointer(theta_levels, theta_levels_view, [ntheta_f])
    call c_f_pointer(output_flat, output_flat_view, [nouter_f * ntheta_f * ninner_f])

    call dinterp_to_isentropic_corder_into_impl( &
        data_flat_view, temp_flat_view, pressure_flat_view, theta_levels_view, output_flat_view, &
        p0, kappa, spvl, kxtrp_f, nouter_f, nlev_f, ninner_f, ntheta_f &
    )
end subroutine dinterp_to_isentropic_corder_into
