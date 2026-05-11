module spherepack_reduced_bindings_mod
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer, c_int, c_float, c_double, c_float_complex
    implicit none

contains

subroutine reduced_gaqd_c(nlat, theta, wts, ierror) bind(C, name="reduced_gaqd_c")
    integer(c_int), value, intent(in) :: nlat
    type(c_ptr), value, intent(in) :: theta
    type(c_ptr), value, intent(in) :: wts
    integer(c_int), intent(out) :: ierror

    integer :: nlat_f
    integer :: ldwork
    integer :: ierror_f
    real(c_double), pointer :: theta_view(:)
    real(c_double), pointer :: wts_view(:)
    real(c_double), allocatable :: dwork(:)

    nlat_f = int(nlat)
    ldwork = nlat_f * (nlat_f + 2)

    call c_f_pointer(theta, theta_view, [nlat_f])
    call c_f_pointer(wts, wts_view, [nlat_f])
    allocate(dwork(ldwork))

    call gaqd(nlat_f, theta_view, wts_view, dwork, ldwork, ierror_f)
    ierror = int(ierror_f, kind=c_int)

    deallocate(dwork)
end subroutine reduced_gaqd_c

subroutine reduced_lap_c(dataspec, dataspec_lap, nmdim, nt, rsphere) bind(C, name="reduced_lap_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: dataspec_lap
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    complex(c_float_complex), pointer :: dataspec_view(:, :)
    complex(c_float_complex), pointer :: dataspec_lap_view(:, :)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(dataspec_lap, dataspec_lap_view, [int(nmdim), int(nt)])

    call lap(dataspec_view, dataspec_lap_view, int(nmdim), int(nt), rsphere)
end subroutine reduced_lap_c

subroutine reduced_invlap_c(dataspec, dataspec_ilap, nmdim, nt, rsphere) bind(C, name="reduced_invlap_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: dataspec_ilap
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    complex(c_float_complex), pointer :: dataspec_view(:, :)
    complex(c_float_complex), pointer :: dataspec_ilap_view(:, :)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(dataspec_ilap, dataspec_ilap_view, [int(nmdim), int(nt)])

    call invlap(dataspec_view, dataspec_ilap_view, int(nmdim), int(nt), rsphere)
end subroutine reduced_invlap_c

subroutine reduced_multsmoothfact_c(dataspec, dataspec_smooth, smooth, nlat, nmdim, nt) bind(C, name="reduced_multsmoothfact_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: dataspec_smooth
    type(c_ptr), value, intent(in) :: smooth
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt

    complex(c_float_complex), pointer :: dataspec_view(:, :)
    complex(c_float_complex), pointer :: dataspec_smooth_view(:, :)
    real(c_float), pointer :: smooth_view(:)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(dataspec_smooth, dataspec_smooth_view, [int(nmdim), int(nt)])
    call c_f_pointer(smooth, smooth_view, [int(nlat)])

    call multsmoothfact(dataspec_view, dataspec_smooth_view, smooth_view, int(nlat), int(nmdim), int(nt))
end subroutine reduced_multsmoothfact_c

subroutine reduced_gaussian_legendre_basis_c(basis, nlat, ntrunc, ierror) bind(C, name="reduced_gaussian_legendre_basis_c")
    type(c_ptr), value, intent(in) :: basis
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), intent(out) :: ierror

    integer :: nmdim
    integer :: ierror_f
    real(c_float), pointer :: basis_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(basis, basis_view, [int(nlat), nmdim])
    call reduced_gaussian_legendre_basis(basis_view, int(nlat), int(ntrunc), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_legendre_basis_c

subroutine reduced_gaussian_legendre_derivative_from_basis_c(basis, dbasis, nlat, ntrunc, ierror) bind(C, name="reduced_gaussian_legendre_derivative_from_basis_c")
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), intent(out) :: ierror

    integer :: nmdim
    integer :: ierror_f
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(basis, basis_view, [int(nlat), nmdim])
    call c_f_pointer(dbasis, dbasis_view, [int(nlat), nmdim])
    call reduced_gaussian_legendre_derivative_from_basis( &
        basis_view, dbasis_view, int(nlat), int(ntrunc), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_legendre_derivative_from_basis_c

subroutine reduced_gaussian_legendre_derivative_basis_c(dbasis, nlat, ntrunc, ierror) bind(C, name="reduced_gaussian_legendre_derivative_basis_c")
    type(c_ptr), value, intent(in) :: dbasis
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), intent(out) :: ierror

    integer :: nmdim
    integer :: ierror_f
    real(c_float), pointer :: dbasis_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(dbasis, dbasis_view, [int(nlat), nmdim])
    call reduced_gaussian_legendre_derivative_basis( &
        dbasis_view, int(nlat), int(ntrunc), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_legendre_derivative_basis_c

subroutine reduced_gaussian_grdtospec_c(datagrid, pl, weights, basis, dataspec, ngptot, nlat, ntrunc, nt, ierror) bind(C, name="reduced_gaussian_grdtospec_c")
    type(c_ptr), value, intent(in) :: datagrid
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: weights
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dataspec
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    integer :: nmdim
    real(c_float), pointer :: datagrid_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: weights_view(:)
    real(c_float), pointer :: basis_view(:, :)
    complex(c_float_complex), pointer :: dataspec_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(datagrid, datagrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(weights, weights_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [int(nlat), nmdim])
    call c_f_pointer(dataspec, dataspec_view, [nmdim, int(nt)])

    call reduced_gaussian_grdtospec( &
        datagrid_view, pl_view, weights_view, basis_view, dataspec_view, &
        int(ngptot), int(nlat), int(ntrunc), int(nt), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_grdtospec_c

subroutine reduced_gaussian_spectogrd_c(dataspec, pl, basis, datagrid, nmdim, nlat, ntrunc, nt, ngptot, ierror) bind(C, name="reduced_gaussian_spectogrd_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: datagrid
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    complex(c_float_complex), pointer :: dataspec_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: datagrid_view(:, :)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(datagrid, datagrid_view, [int(ngptot), int(nt)])

    call reduced_gaussian_spectogrd( &
        dataspec_view, pl_view, basis_view, datagrid_view, int(nmdim), int(nlat), &
        int(ntrunc), int(nt), int(ngptot), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_spectogrd_c

subroutine reduced_gaussian_spectogrd_pair_c(speca, specb, pl, basis, grida, gridb, nmdim, nlat, ntrunc, nt, ngptot, ierror) bind(C, name="reduced_gaussian_spectogrd_pair_c")
    type(c_ptr), value, intent(in) :: speca
    type(c_ptr), value, intent(in) :: specb
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: grida
    type(c_ptr), value, intent(in) :: gridb
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    complex(c_float_complex), pointer :: speca_view(:, :)
    complex(c_float_complex), pointer :: specb_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: grida_view(:, :)
    real(c_float), pointer :: gridb_view(:, :)

    call c_f_pointer(speca, speca_view, [int(nmdim), int(nt)])
    call c_f_pointer(specb, specb_view, [int(nmdim), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(grida, grida_view, [int(ngptot), int(nt)])
    call c_f_pointer(gridb, gridb_view, [int(ngptot), int(nt)])

    call reduced_gaussian_spectogrd_pair( &
        speca_view, specb_view, pl_view, basis_view, grida_view, gridb_view, &
        int(nmdim), int(nlat), int(ntrunc), int(nt), int(ngptot), ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_spectogrd_pair_c

subroutine reduced_gaussian_getgrad_c(dataspec, pl, basis, dbasis, sin_theta, ugrad, vgrad, nmdim, nlat, ntrunc, nt, ngptot, rsphere, ierror) bind(C, name="reduced_gaussian_getgrad_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    type(c_ptr), value, intent(in) :: sin_theta
    type(c_ptr), value, intent(in) :: ugrad
    type(c_ptr), value, intent(in) :: vgrad
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    integer(c_int), value, intent(in) :: ngptot
    real(c_float), value, intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    complex(c_float_complex), pointer :: dataspec_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)
    real(c_float), pointer :: sin_theta_view(:)
    real(c_float), pointer :: ugrad_view(:, :)
    real(c_float), pointer :: vgrad_view(:, :)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(dbasis, dbasis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(sin_theta, sin_theta_view, [int(nlat)])
    call c_f_pointer(ugrad, ugrad_view, [int(ngptot), int(nt)])
    call c_f_pointer(vgrad, vgrad_view, [int(ngptot), int(nt)])

    call reduced_gaussian_getgrad( &
        dataspec_view, pl_view, basis_view, dbasis_view, sin_theta_view, &
        ugrad_view, vgrad_view, int(nmdim), int(nlat), int(ntrunc), int(nt), &
        int(ngptot), rsphere, ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_getgrad_c

subroutine reduced_gaussian_getgrad_pair_c(speca, specb, pl, basis, dbasis, sin_theta, a_ugrad, a_vgrad, b_ugrad, b_vgrad, nmdim, nlat, ntrunc, nt, ngptot, rsphere, ierror) bind(C, name="reduced_gaussian_getgrad_pair_c")
    type(c_ptr), value, intent(in) :: speca
    type(c_ptr), value, intent(in) :: specb
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    type(c_ptr), value, intent(in) :: sin_theta
    type(c_ptr), value, intent(in) :: a_ugrad
    type(c_ptr), value, intent(in) :: a_vgrad
    type(c_ptr), value, intent(in) :: b_ugrad
    type(c_ptr), value, intent(in) :: b_vgrad
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    integer(c_int), value, intent(in) :: ngptot
    real(c_float), value, intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    complex(c_float_complex), pointer :: speca_view(:, :)
    complex(c_float_complex), pointer :: specb_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)
    real(c_float), pointer :: sin_theta_view(:)
    real(c_float), pointer :: a_ugrad_view(:, :)
    real(c_float), pointer :: a_vgrad_view(:, :)
    real(c_float), pointer :: b_ugrad_view(:, :)
    real(c_float), pointer :: b_vgrad_view(:, :)

    call c_f_pointer(speca, speca_view, [int(nmdim), int(nt)])
    call c_f_pointer(specb, specb_view, [int(nmdim), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(dbasis, dbasis_view, [int(nmdim), int(nlat)])
    call c_f_pointer(sin_theta, sin_theta_view, [int(nlat)])
    call c_f_pointer(a_ugrad, a_ugrad_view, [int(ngptot), int(nt)])
    call c_f_pointer(a_vgrad, a_vgrad_view, [int(ngptot), int(nt)])
    call c_f_pointer(b_ugrad, b_ugrad_view, [int(ngptot), int(nt)])
    call c_f_pointer(b_vgrad, b_vgrad_view, [int(ngptot), int(nt)])

    call reduced_gaussian_getgrad_pair( &
        speca_view, specb_view, pl_view, basis_view, dbasis_view, sin_theta_view, &
        a_ugrad_view, a_vgrad_view, b_ugrad_view, b_vgrad_view, int(nmdim), &
        int(nlat), int(ntrunc), int(nt), int(ngptot), rsphere, ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_getgrad_pair_c

subroutine reduced_gaussian_getvrtdivspec_c(ugrid, vgrid, pl, weights, basis, dbasis, sin_theta, vrtspec, divspec, ngptot, nlat, ntrunc, nt, rsphere, ierror) bind(C, name="reduced_gaussian_getvrtdivspec_c")
    type(c_ptr), value, intent(in) :: ugrid
    type(c_ptr), value, intent(in) :: vgrid
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: weights
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    type(c_ptr), value, intent(in) :: sin_theta
    type(c_ptr), value, intent(in) :: vrtspec
    type(c_ptr), value, intent(in) :: divspec
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    integer :: nmdim
    real(c_float), pointer :: ugrid_view(:, :)
    real(c_float), pointer :: vgrid_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: weights_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)
    real(c_float), pointer :: sin_theta_view(:)
    complex(c_float_complex), pointer :: vrtspec_view(:, :)
    complex(c_float_complex), pointer :: divspec_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(ugrid, ugrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(vgrid, vgrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(weights, weights_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [nmdim, int(nlat)])
    call c_f_pointer(dbasis, dbasis_view, [nmdim, int(nlat)])
    call c_f_pointer(sin_theta, sin_theta_view, [int(nlat)])
    call c_f_pointer(vrtspec, vrtspec_view, [nmdim, int(nt)])
    call c_f_pointer(divspec, divspec_view, [nmdim, int(nt)])

    call reduced_gaussian_getvrtdivspec( &
        ugrid_view, vgrid_view, pl_view, weights_view, basis_view, dbasis_view, &
        sin_theta_view, vrtspec_view, divspec_view, int(ngptot), int(nlat), &
        int(ntrunc), int(nt), rsphere, ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_getvrtdivspec_c

subroutine reduced_gaussian_getvrtspec_c(ugrid, vgrid, pl, weights, basis, dbasis, sin_theta, vrtspec, ngptot, nlat, ntrunc, nt, rsphere, ierror) bind(C, name="reduced_gaussian_getvrtspec_c")
    type(c_ptr), value, intent(in) :: ugrid
    type(c_ptr), value, intent(in) :: vgrid
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: weights
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    type(c_ptr), value, intent(in) :: sin_theta
    type(c_ptr), value, intent(in) :: vrtspec
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    integer :: nmdim
    real(c_float), pointer :: ugrid_view(:, :)
    real(c_float), pointer :: vgrid_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: weights_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)
    real(c_float), pointer :: sin_theta_view(:)
    complex(c_float_complex), pointer :: vrtspec_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(ugrid, ugrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(vgrid, vgrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(weights, weights_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [nmdim, int(nlat)])
    call c_f_pointer(dbasis, dbasis_view, [nmdim, int(nlat)])
    call c_f_pointer(sin_theta, sin_theta_view, [int(nlat)])
    call c_f_pointer(vrtspec, vrtspec_view, [nmdim, int(nt)])

    call reduced_gaussian_getvrtspec( &
        ugrid_view, vgrid_view, pl_view, weights_view, basis_view, dbasis_view, &
        sin_theta_view, vrtspec_view, int(ngptot), int(nlat), int(ntrunc), &
        int(nt), rsphere, ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_getvrtspec_c

subroutine reduced_gaussian_getdivspec_c(ugrid, vgrid, pl, weights, basis, dbasis, sin_theta, divspec, ngptot, nlat, ntrunc, nt, rsphere, ierror) bind(C, name="reduced_gaussian_getdivspec_c")
    type(c_ptr), value, intent(in) :: ugrid
    type(c_ptr), value, intent(in) :: vgrid
    type(c_ptr), value, intent(in) :: pl
    type(c_ptr), value, intent(in) :: weights
    type(c_ptr), value, intent(in) :: basis
    type(c_ptr), value, intent(in) :: dbasis
    type(c_ptr), value, intent(in) :: sin_theta
    type(c_ptr), value, intent(in) :: divspec
    integer(c_int), value, intent(in) :: ngptot
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere
    integer(c_int), intent(out) :: ierror

    integer :: ierror_f
    integer :: nmdim
    real(c_float), pointer :: ugrid_view(:, :)
    real(c_float), pointer :: vgrid_view(:, :)
    integer, pointer :: pl_view(:)
    real(c_float), pointer :: weights_view(:)
    real(c_float), pointer :: basis_view(:, :)
    real(c_float), pointer :: dbasis_view(:, :)
    real(c_float), pointer :: sin_theta_view(:)
    complex(c_float_complex), pointer :: divspec_view(:, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2
    call c_f_pointer(ugrid, ugrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(vgrid, vgrid_view, [int(ngptot), int(nt)])
    call c_f_pointer(pl, pl_view, [int(nlat)])
    call c_f_pointer(weights, weights_view, [int(nlat)])
    call c_f_pointer(basis, basis_view, [nmdim, int(nlat)])
    call c_f_pointer(dbasis, dbasis_view, [nmdim, int(nlat)])
    call c_f_pointer(sin_theta, sin_theta_view, [int(nlat)])
    call c_f_pointer(divspec, divspec_view, [nmdim, int(nt)])

    call reduced_gaussian_getdivspec( &
        ugrid_view, vgrid_view, pl_view, weights_view, basis_view, dbasis_view, &
        sin_theta_view, divspec_view, int(ngptot), int(nlat), int(ntrunc), &
        int(nt), rsphere, ierror_f)
    ierror = int(ierror_f, kind=c_int)
end subroutine reduced_gaussian_getdivspec_c

end module spherepack_reduced_bindings_mod
