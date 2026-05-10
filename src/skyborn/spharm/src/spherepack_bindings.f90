module spherepack_bindings_mod
    use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer, c_int, c_float, c_double, c_float_complex
    implicit none

contains

subroutine gaqd_c(nlat, theta, wts, ierror) bind(C, name="gaqd_c")
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
end subroutine gaqd_c

subroutine getlegfunc_c(legfunc, lat, ntrunc) bind(C, name="getlegfunc_c")
    type(c_ptr), value, intent(in) :: legfunc
    real(c_float), value, intent(in) :: lat
    integer(c_int), value, intent(in) :: ntrunc

    integer :: ntrunc_f
    integer :: nmdim
    real(c_float), pointer :: legfunc_view(:)

    ntrunc_f = int(ntrunc)
    nmdim = (ntrunc_f + 1) * (ntrunc_f + 2) / 2

    call c_f_pointer(legfunc, legfunc_view, [nmdim])
    call getlegfunc(legfunc_view, lat, ntrunc_f)
end subroutine getlegfunc_c

subroutine specintrp_c(rlon, ntrunc, datnm, pnm, ob) bind(C, name="specintrp_c")
    real(c_float), value, intent(in) :: rlon
    integer(c_int), value, intent(in) :: ntrunc
    type(c_ptr), value, intent(in) :: datnm
    type(c_ptr), value, intent(in) :: pnm
    real(c_float), intent(out) :: ob

    integer :: ntrunc_f
    integer :: nmdim
    complex(c_float_complex), pointer :: datnm_view(:)
    real(c_float), pointer :: pnm_view(:)
    complex(c_float_complex), allocatable :: scrm(:)

    ntrunc_f = int(ntrunc)
    nmdim = (ntrunc_f + 1) * (ntrunc_f + 2) / 2

    call c_f_pointer(datnm, datnm_view, [nmdim])
    call c_f_pointer(pnm, pnm_view, [nmdim])
    allocate(scrm(ntrunc_f + 1))

    call specintrp(rlon, ntrunc_f, datnm_view, scrm, pnm_view, ob)

    deallocate(scrm)
end subroutine specintrp_c

subroutine lap_c(dataspec, dataspec_lap, nmdim, nt, rsphere) bind(C, name="lap_c")
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
end subroutine lap_c

subroutine invlap_c(dataspec, dataspec_ilap, nmdim, nt, rsphere) bind(C, name="invlap_c")
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
end subroutine invlap_c

subroutine multsmoothfact_c(dataspec, dataspec_smooth, smooth, nlat, nmdim, nt) bind(C, name="multsmoothfact_c")
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
end subroutine multsmoothfact_c

subroutine onedtotwod_c(dataspec, a, b, nlat, nmdim, nt) bind(C, name="onedtotwod_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt

    complex(c_float_complex), pointer :: dataspec_view(:, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)

    call c_f_pointer(dataspec, dataspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])

    call onedtotwod(dataspec_view, a_view, b_view, int(nlat), int(nmdim), int(nt))
end subroutine onedtotwod_c

subroutine onedtotwod_vrtdiv_c(vrtspec, divspec, br, bi, cr, ci, nlat, nmdim, nt, rsphere) bind(C, name="onedtotwod_vrtdiv_c")
    type(c_ptr), value, intent(in) :: vrtspec
    type(c_ptr), value, intent(in) :: divspec
    type(c_ptr), value, intent(in) :: br
    type(c_ptr), value, intent(in) :: bi
    type(c_ptr), value, intent(in) :: cr
    type(c_ptr), value, intent(in) :: ci
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    complex(c_float_complex), pointer :: vrtspec_view(:, :)
    complex(c_float_complex), pointer :: divspec_view(:, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)

    call c_f_pointer(vrtspec, vrtspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(divspec, divspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])

    call onedtotwod_vrtdiv(vrtspec_view, divspec_view, br_view, bi_view, cr_view, ci_view, int(nlat), int(nmdim), int(nt), rsphere)
end subroutine onedtotwod_vrtdiv_c

subroutine onedtotwod_vrt_c(vrtspec, cr, ci, nlat, nmdim, nt, rsphere) bind(C, name="onedtotwod_vrt_c")
    type(c_ptr), value, intent(in) :: vrtspec
    type(c_ptr), value, intent(in) :: cr
    type(c_ptr), value, intent(in) :: ci
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    complex(c_float_complex), pointer :: vrtspec_view(:, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)

    call c_f_pointer(vrtspec, vrtspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])

    call onedtotwod_vrt(vrtspec_view, cr_view, ci_view, int(nlat), int(nmdim), int(nt), rsphere)
end subroutine onedtotwod_vrt_c

subroutine onedtotwod_div_c(divspec, br, bi, nlat, nmdim, nt, rsphere) bind(C, name="onedtotwod_div_c")
    type(c_ptr), value, intent(in) :: divspec
    type(c_ptr), value, intent(in) :: br
    type(c_ptr), value, intent(in) :: bi
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: nmdim
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    complex(c_float_complex), pointer :: divspec_view(:, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)

    call c_f_pointer(divspec, divspec_view, [int(nmdim), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])

    call onedtotwod_div(divspec_view, br_view, bi_view, int(nlat), int(nmdim), int(nt), rsphere)
end subroutine onedtotwod_div_c

subroutine twodtooned_c(dataspec, a, b, nlat, ntrunc, nt) bind(C, name="twodtooned_c")
    type(c_ptr), value, intent(in) :: dataspec
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt

    integer :: nmdim
    complex(c_float_complex), pointer :: dataspec_view(:, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2

    call c_f_pointer(dataspec, dataspec_view, [nmdim, int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])

    call twodtooned(dataspec_view, a_view, b_view, int(nlat), int(ntrunc), int(nt))
end subroutine twodtooned_c

subroutine twodtooned_vrtdiv_c(vrtspec, divspec, br, bi, cr, ci, nlat, ntrunc, nt, rsphere) bind(C, name="twodtooned_vrtdiv_c")
    type(c_ptr), value, intent(in) :: vrtspec
    type(c_ptr), value, intent(in) :: divspec
    type(c_ptr), value, intent(in) :: br
    type(c_ptr), value, intent(in) :: bi
    type(c_ptr), value, intent(in) :: cr
    type(c_ptr), value, intent(in) :: ci
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    integer :: nmdim
    complex(c_float_complex), pointer :: vrtspec_view(:, :)
    complex(c_float_complex), pointer :: divspec_view(:, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2

    call c_f_pointer(vrtspec, vrtspec_view, [nmdim, int(nt)])
    call c_f_pointer(divspec, divspec_view, [nmdim, int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])

    call twodtooned_vrtdiv(vrtspec_view, divspec_view, br_view, bi_view, cr_view, ci_view, int(nlat), int(ntrunc), int(nt), rsphere)
end subroutine twodtooned_vrtdiv_c

subroutine twodtooned_vrt_c(vrtspec, cr, ci, nlat, ntrunc, nt, rsphere) bind(C, name="twodtooned_vrt_c")
    type(c_ptr), value, intent(in) :: vrtspec
    type(c_ptr), value, intent(in) :: cr
    type(c_ptr), value, intent(in) :: ci
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    integer :: nmdim
    complex(c_float_complex), pointer :: vrtspec_view(:, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2

    call c_f_pointer(vrtspec, vrtspec_view, [nmdim, int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])

    call twodtooned_vrt(vrtspec_view, cr_view, ci_view, int(nlat), int(ntrunc), int(nt), rsphere)
end subroutine twodtooned_vrt_c

subroutine twodtooned_div_c(divspec, br, bi, nlat, ntrunc, nt, rsphere) bind(C, name="twodtooned_div_c")
    type(c_ptr), value, intent(in) :: divspec
    type(c_ptr), value, intent(in) :: br
    type(c_ptr), value, intent(in) :: bi
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: ntrunc
    integer(c_int), value, intent(in) :: nt
    real(c_float), value, intent(in) :: rsphere

    integer :: nmdim
    complex(c_float_complex), pointer :: divspec_view(:, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)

    nmdim = (int(ntrunc) + 1) * (int(ntrunc) + 2) / 2

    call c_f_pointer(divspec, divspec_view, [nmdim, int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])

    call twodtooned_div(divspec_view, br_view, bi_view, int(nlat), int(ntrunc), int(nt), rsphere)
end subroutine twodtooned_div_c

subroutine ihgeod_c(m, x, y, z) bind(C, name="ihgeod_c")
    integer(c_int), value, intent(in) :: m
    type(c_ptr), value, intent(in) :: x
    type(c_ptr), value, intent(in) :: y
    type(c_ptr), value, intent(in) :: z

    integer :: m_f
    integer :: idp
    integer :: jdp
    real(c_float), pointer :: x_view(:, :, :)
    real(c_float), pointer :: y_view(:, :, :)
    real(c_float), pointer :: z_view(:, :, :)

    m_f = int(m)
    idp = 2 * m_f - 1
    jdp = m_f

    call c_f_pointer(x, x_view, [idp, jdp, 5])
    call c_f_pointer(y, y_view, [idp, jdp, 5])
    call c_f_pointer(z, z_view, [idp, jdp, 5])

    call ihgeod(m_f, idp, jdp, x_view, y_view, z_view)
end subroutine ihgeod_c

subroutine shaesi_c(nlat, nlon, lshaes, lwork, ldwork, wshaes, ierror) bind(C, name="shaesi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshaes, lwork, ldwork
    type(c_ptr), value, intent(in) :: wshaes
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshaes_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshaes, wshaes_view, [int(lshaes)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call shaesi(int(nlat), int(nlon), wshaes_view, int(lshaes), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine shaesi_c

subroutine shsesi_c(nlat, nlon, lshses, lwork, ldwork, wshses, ierror) bind(C, name="shsesi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshses, lwork, ldwork
    type(c_ptr), value, intent(in) :: wshses
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshses_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshses, wshses_view, [int(lshses)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call shsesi(int(nlat), int(nlon), wshses_view, int(lshses), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine shsesi_c

subroutine shagsi_c(nlat, nlon, lshags, lwork, ldwork, wshags, ierror) bind(C, name="shagsi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshags, lwork, ldwork
    type(c_ptr), value, intent(in) :: wshags
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshags_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshags, wshags_view, [int(lshags)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call shagsi(int(nlat), int(nlon), wshags_view, int(lshags), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine shagsi_c

subroutine shsgsi_c(nlat, nlon, lshsgs, lwork, ldwork, wshsgs, ierror) bind(C, name="shsgsi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshsgs, lwork, ldwork
    type(c_ptr), value, intent(in) :: wshsgs
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshsgs_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshsgs, wshsgs_view, [int(lshsgs)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call shsgsi(int(nlat), int(nlon), wshsgs_view, int(lshsgs), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine shsgsi_c

subroutine shaeci_c(nlat, nlon, lshaec, ldwork, wshaec, ierror) bind(C, name="shaeci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshaec, ldwork
    type(c_ptr), value, intent(in) :: wshaec
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshaec_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshaec, wshaec_view, [int(lshaec)])
    allocate(dwork(int(ldwork)))
    call shaeci(int(nlat), int(nlon), wshaec_view, int(lshaec), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine shaeci_c

subroutine shseci_c(nlat, nlon, lshsec, ldwork, wshsec, ierror) bind(C, name="shseci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshsec, ldwork
    type(c_ptr), value, intent(in) :: wshsec
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshsec_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshsec, wshsec_view, [int(lshsec)])
    allocate(dwork(int(ldwork)))
    call shseci(int(nlat), int(nlon), wshsec_view, int(lshsec), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine shseci_c

subroutine shagci_c(nlat, nlon, lshagc, ldwork, wshagc, ierror) bind(C, name="shagci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshagc, ldwork
    type(c_ptr), value, intent(in) :: wshagc
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshagc_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshagc, wshagc_view, [int(lshagc)])
    allocate(dwork(int(ldwork)))
    call shagci(int(nlat), int(nlon), wshagc_view, int(lshagc), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine shagci_c

subroutine shsgci_c(nlat, nlon, lshsgc, ldwork, wshsgc, ierror) bind(C, name="shsgci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lshsgc, ldwork
    type(c_ptr), value, intent(in) :: wshsgc
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wshsgc_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wshsgc, wshsgc_view, [int(lshsgc)])
    allocate(dwork(int(ldwork)))
    call shsgci(int(nlat), int(nlon), wshsgc_view, int(lshsgc), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine shsgci_c

subroutine shaes_c(g, nlat, nlon, nt, wshaes, lshaes, lwork, a, b, ierror) bind(C, name="shaes_c")
    type(c_ptr), value, intent(in) :: g, wshaes, a, b
    integer(c_int), value, intent(in) :: nlat, nlon, nt, lshaes, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: wshaes_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wshaes, wshaes_view, [int(lshaes)])
    allocate(work(int(lwork)))
    call shaes(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshaes_view, int(lshaes), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shaes_c

subroutine shags_c(g, nlat, nlon, nt, wshags, lshags, lwork, a, b, ierror) bind(C, name="shags_c")
    type(c_ptr), value, intent(in) :: g, wshags, a, b
    integer(c_int), value, intent(in) :: nlat, nlon, nt, lshags, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: wshags_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wshags, wshags_view, [int(lshags)])
    allocate(work(int(lwork)))
    call shags(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshags_view, int(lshags), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shags_c

subroutine shaec_c(g, nlat, nlon, nt, wshaec, lshaec, lwork, a, b, ierror) bind(C, name="shaec_c")
    type(c_ptr), value, intent(in) :: g, wshaec, a, b
    integer(c_int), value, intent(in) :: nlat, nlon, nt, lshaec, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: wshaec_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wshaec, wshaec_view, [int(lshaec)])
    allocate(work(int(lwork)))
    call shaec(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshaec_view, int(lshaec), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shaec_c

subroutine shagc_c(g, nlat, nlon, nt, wshagc, lshagc, lwork, a, b, ierror) bind(C, name="shagc_c")
    type(c_ptr), value, intent(in) :: g, wshagc, a, b
    integer(c_int), value, intent(in) :: nlat, nlon, nt, lshagc, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: wshagc_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wshagc, wshagc_view, [int(lshagc)])
    allocate(work(int(lwork)))
    call shagc(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshagc_view, int(lshagc), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shagc_c

subroutine shses_c(nlon, a, b, nlat, nt, wshses, lshses, lwork, g, ierror) bind(C, name="shses_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, lshses, lwork
    type(c_ptr), value, intent(in) :: a, b, wshses, g
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: wshses_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wshses, wshses_view, [int(lshses)])
    allocate(work(int(lwork)))
    call shses(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshses_view, int(lshses), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shses_c

subroutine shsgs_c(nlon, a, b, nlat, nt, wshsgs, lshsgs, lwork, g, ierror) bind(C, name="shsgs_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, lshsgs, lwork
    type(c_ptr), value, intent(in) :: a, b, wshsgs, g
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: wshsgs_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wshsgs, wshsgs_view, [int(lshsgs)])
    allocate(work(int(lwork)))
    call shsgs(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshsgs_view, int(lshsgs), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shsgs_c

subroutine shsec_c(nlon, a, b, nlat, nt, wshsec, lshsec, lwork, g, ierror) bind(C, name="shsec_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, lshsec, lwork
    type(c_ptr), value, intent(in) :: a, b, wshsec, g
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: wshsec_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wshsec, wshsec_view, [int(lshsec)])
    allocate(work(int(lwork)))
    call shsec(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshsec_view, int(lshsec), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shsec_c

subroutine shsgc_c(nlon, a, b, nlat, nt, wshsgc, lshsgc, lwork, g, ierror) bind(C, name="shsgc_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, lshsgc, lwork
    type(c_ptr), value, intent(in) :: a, b, wshsgc, g
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: a_view(:, :, :)
    real(c_float), pointer :: b_view(:, :, :)
    real(c_float), pointer :: g_view(:, :, :)
    real(c_float), pointer :: wshsgc_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(a, a_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(b, b_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(g, g_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wshsgc, wshsgc_view, [int(lshsgc)])
    allocate(work(int(lwork)))
    call shsgc(int(nlat), int(nlon), 0, int(nt), g_view, int(nlat), int(nlon), a_view, b_view, int(nlat), int(nlat), wshsgc_view, int(lshsgc), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine shsgc_c

subroutine vhaesi_c(nlat, nlon, lvhaes, lwork, ldwork, wvhaes, ierror) bind(C, name="vhaesi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhaes, lwork, ldwork
    type(c_ptr), value, intent(in) :: wvhaes
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhaes_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhaes, wvhaes_view, [int(lvhaes)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call vhaesi(int(nlat), int(nlon), wvhaes_view, int(lvhaes), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine vhaesi_c

subroutine vhagsi_c(nlat, nlon, lvhags, ldwork, wvhags, ierror) bind(C, name="vhagsi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhags, ldwork
    type(c_ptr), value, intent(in) :: wvhags
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhags_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhags, wvhags_view, [int(lvhags)])
    allocate(dwork(int(ldwork)))
    call vhagsi(int(nlat), int(nlon), wvhags_view, int(lvhags), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhagsi_c

subroutine vhsesi_c(nlat, nlon, lvhses, lwork, ldwork, wvhses, ierror) bind(C, name="vhsesi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhses, lwork, ldwork
    type(c_ptr), value, intent(in) :: wvhses
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhses_view(:)
    real(c_float), allocatable :: work(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhses, wvhses_view, [int(lvhses)])
    allocate(work(int(lwork)), dwork(int(ldwork)))
    call vhsesi(int(nlat), int(nlon), wvhses_view, int(lvhses), work, int(lwork), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work, dwork)
end subroutine vhsesi_c

subroutine vhsgsi_c(nlat, nlon, lvhsgs, ldwork, wvhsgs, ierror) bind(C, name="vhsgsi_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhsgs, ldwork
    type(c_ptr), value, intent(in) :: wvhsgs
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhsgs_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhsgs, wvhsgs_view, [int(lvhsgs)])
    allocate(dwork(int(ldwork)))
    call vhsgsi(int(nlat), int(nlon), wvhsgs_view, int(lvhsgs), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhsgsi_c

subroutine vhaeci_c(nlat, nlon, lvhaec, ldwork, wvhaec, ierror) bind(C, name="vhaeci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhaec, ldwork
    type(c_ptr), value, intent(in) :: wvhaec
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhaec_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhaec, wvhaec_view, [int(lvhaec)])
    allocate(dwork(int(ldwork)))
    call vhaeci(int(nlat), int(nlon), wvhaec_view, int(lvhaec), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhaeci_c

subroutine vhagci_c(nlat, nlon, lvhagc, ldwork, wvhagc, ierror) bind(C, name="vhagci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhagc, ldwork
    type(c_ptr), value, intent(in) :: wvhagc
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhagc_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhagc, wvhagc_view, [int(lvhagc)])
    allocate(dwork(int(ldwork)))
    call vhagci(int(nlat), int(nlon), wvhagc_view, int(lvhagc), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhagci_c

subroutine vhseci_c(nlat, nlon, lvhsec, ldwork, wvhsec, ierror) bind(C, name="vhseci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhsec, ldwork
    type(c_ptr), value, intent(in) :: wvhsec
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhsec_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhsec, wvhsec_view, [int(lvhsec)])
    allocate(dwork(int(ldwork)))
    call vhseci(int(nlat), int(nlon), wvhsec_view, int(lvhsec), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhseci_c

subroutine vhsgci_c(nlat, nlon, lvhsgc, ldwork, wvhsgc, ierror) bind(C, name="vhsgci_c")
    integer(c_int), value, intent(in) :: nlat, nlon, lvhsgc, ldwork
    type(c_ptr), value, intent(in) :: wvhsgc
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: wvhsgc_view(:)
    real(c_double), allocatable :: dwork(:)
    integer :: ierror_f

    call c_f_pointer(wvhsgc, wvhsgc_view, [int(lvhsgc)])
    allocate(dwork(int(ldwork)))
    call vhsgci(int(nlat), int(nlon), wvhsgc_view, int(lvhsgc), dwork, int(ldwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(dwork)
end subroutine vhsgci_c

subroutine vhaes_c(v, w, nlat, nlon, nt, ityp, wvhaes, lvhaes, lwork, br, bi, cr, ci, ierror) bind(C, name="vhaes_c")
    type(c_ptr), value, intent(in) :: v, w, wvhaes, br, bi, cr, ci
    integer(c_int), value, intent(in) :: nlat, nlon, nt, ityp, lvhaes, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: wvhaes_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wvhaes, wvhaes_view, [int(lvhaes)])
    allocate(work(int(lwork)))
    call vhaes(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhaes_view, int(lvhaes), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhaes_c

subroutine vhags_c(v, w, nlat, nlon, nt, ityp, wvhags, lvhags, lwork, br, bi, cr, ci, ierror) bind(C, name="vhags_c")
    type(c_ptr), value, intent(in) :: v, w, wvhags, br, bi, cr, ci
    integer(c_int), value, intent(in) :: nlat, nlon, nt, ityp, lvhags, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: wvhags_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wvhags, wvhags_view, [int(lvhags)])
    allocate(work(int(lwork)))
    call vhags(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhags_view, int(lvhags), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhags_c

subroutine vhaec_c(v, w, nlat, nlon, nt, ityp, wvhaec, lvhaec, lwork, br, bi, cr, ci, ierror) bind(C, name="vhaec_c")
    type(c_ptr), value, intent(in) :: v, w, wvhaec, br, bi, cr, ci
    integer(c_int), value, intent(in) :: nlat, nlon, nt, ityp, lvhaec, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: wvhaec_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wvhaec, wvhaec_view, [int(lvhaec)])
    allocate(work(int(lwork)))
    call vhaec(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhaec_view, int(lvhaec), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhaec_c

subroutine vhagc_c(v, w, nlat, nlon, nt, ityp, wvhagc, lvhagc, lwork, br, bi, cr, ci, ierror) bind(C, name="vhagc_c")
    type(c_ptr), value, intent(in) :: v, w, wvhagc, br, bi, cr, ci
    integer(c_int), value, intent(in) :: nlat, nlon, nt, ityp, lvhagc, lwork
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: wvhagc_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(wvhagc, wvhagc_view, [int(lvhagc)])
    allocate(work(int(lwork)))
    call vhagc(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhagc_view, int(lvhagc), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhagc_c

subroutine vhses_c(nlon, br, bi, cr, ci, nlat, nt, ityp, wvhses, lvhses, lwork, v, w, ierror) bind(C, name="vhses_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, ityp, lvhses, lwork
    type(c_ptr), value, intent(in) :: br, bi, cr, ci, wvhses, v, w
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: wvhses_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wvhses, wvhses_view, [int(lvhses)])
    allocate(work(int(lwork)))
    call vhses(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhses_view, int(lvhses), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhses_c

subroutine vhsgs_c(nlon, br, bi, cr, ci, nlat, nt, ityp, wvhsgs, lvhsgs, lwork, v, w, ierror) bind(C, name="vhsgs_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, ityp, lvhsgs, lwork
    type(c_ptr), value, intent(in) :: br, bi, cr, ci, wvhsgs, v, w
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: wvhsgs_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wvhsgs, wvhsgs_view, [int(lvhsgs)])
    allocate(work(int(lwork)))
    call vhsgs(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhsgs_view, int(lvhsgs), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhsgs_c

subroutine vhsec_c(nlon, br, bi, cr, ci, nlat, nt, ityp, wvhsec, lvhsec, lwork, v, w, ierror) bind(C, name="vhsec_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, ityp, lvhsec, lwork
    type(c_ptr), value, intent(in) :: br, bi, cr, ci, wvhsec, v, w
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: wvhsec_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wvhsec, wvhsec_view, [int(lvhsec)])
    allocate(work(int(lwork)))
    call vhsec(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhsec_view, int(lvhsec), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhsec_c

subroutine vhsgc_c(nlon, br, bi, cr, ci, nlat, nt, ityp, wvhsgc, lvhsgc, lwork, v, w, ierror) bind(C, name="vhsgc_c")
    integer(c_int), value, intent(in) :: nlon, nlat, nt, ityp, lvhsgc, lwork
    type(c_ptr), value, intent(in) :: br, bi, cr, ci, wvhsgc, v, w
    integer(c_int), intent(out) :: ierror

    real(c_float), pointer :: br_view(:, :, :)
    real(c_float), pointer :: bi_view(:, :, :)
    real(c_float), pointer :: cr_view(:, :, :)
    real(c_float), pointer :: ci_view(:, :, :)
    real(c_float), pointer :: v_view(:, :, :)
    real(c_float), pointer :: w_view(:, :, :)
    real(c_float), pointer :: wvhsgc_view(:)
    real(c_float), allocatable :: work(:)
    integer :: ierror_f

    call c_f_pointer(br, br_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(bi, bi_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(cr, cr_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(ci, ci_view, [int(nlat), int(nlat), int(nt)])
    call c_f_pointer(v, v_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(w, w_view, [int(nlat), int(nlon), int(nt)])
    call c_f_pointer(wvhsgc, wvhsgc_view, [int(lvhsgc)])
    allocate(work(int(lwork)))
    call vhsgc(int(nlat), int(nlon), int(ityp), int(nt), v_view, w_view, int(nlat), int(nlon), br_view, bi_view, cr_view, ci_view, int(nlat), int(nlat), wvhsgc_view, int(lvhsgc), work, int(lwork), ierror_f)
    ierror = int(ierror_f, kind=c_int)
    deallocate(work)
end subroutine vhsgc_c

end module spherepack_bindings_mod
