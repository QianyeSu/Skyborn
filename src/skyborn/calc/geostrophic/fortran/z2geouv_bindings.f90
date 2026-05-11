! =============================================================================
! File    : z2geouv_bindings.f90
! Purpose : C interoperability bindings split from z2geouv.f90.
! Notes   : Keep bind(C) entry points separate from the scientific kernels.
! =============================================================================

subroutine z2geouv_c(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C, name="z2geouv_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use geostrophicwind_mod, only : z2geouv
    implicit none

    type(c_ptr), value, intent(in) :: z, glon, glat, ug, vg
    integer(c_int), value, intent(in) :: nlat, mlon, iopt
    real(c_double), value, intent(in) :: zmsg

    integer :: nlat_f, mlon_f
    real(c_double), pointer :: z_view(:, :)
    real(c_double), pointer :: glon_view(:)
    real(c_double), pointer :: glat_view(:)
    real(c_double), pointer :: ug_view(:, :)
    real(c_double), pointer :: vg_view(:, :)

    nlat_f = int(nlat)
    mlon_f = int(mlon)

    call c_f_pointer(z, z_view, [nlat_f, mlon_f])
    call c_f_pointer(glon, glon_view, [mlon_f])
    call c_f_pointer(glat, glat_view, [nlat_f])
    call c_f_pointer(ug, ug_view, [nlat_f, mlon_f])
    call c_f_pointer(vg, vg_view, [nlat_f, mlon_f])

    call z2geouv(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine z2geouv_c


subroutine z2geouv_3d_c(z, nlat, mlon, n3rd, zmsg, glon, glat, ug, vg, iopt) bind(C, name="z2geouv_3d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use geostrophicwind_mod, only : z2geouv_3d
    implicit none

    type(c_ptr), value, intent(in) :: z, glon, glat, ug, vg
    integer(c_int), value, intent(in) :: nlat, mlon, n3rd, iopt
    real(c_double), value, intent(in) :: zmsg

    integer :: nlat_f, mlon_f, n3rd_f
    real(c_double), pointer :: z_view(:, :, :)
    real(c_double), pointer :: glon_view(:)
    real(c_double), pointer :: glat_view(:)
    real(c_double), pointer :: ug_view(:, :, :)
    real(c_double), pointer :: vg_view(:, :, :)

    nlat_f = int(nlat)
    mlon_f = int(mlon)
    n3rd_f = int(n3rd)

    call c_f_pointer(z, z_view, [nlat_f, mlon_f, n3rd_f])
    call c_f_pointer(glon, glon_view, [mlon_f])
    call c_f_pointer(glat, glat_view, [nlat_f])
    call c_f_pointer(ug, ug_view, [nlat_f, mlon_f, n3rd_f])
    call c_f_pointer(vg, vg_view, [nlat_f, mlon_f, n3rd_f])

    call z2geouv_3d(z_view, nlat_f, mlon_f, n3rd_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine z2geouv_3d_c


subroutine z2geouv_4d_c(z, nlat, mlon, n3rd, n4th, zmsg, glon, glat, ug, vg, iopt) bind(C, name="z2geouv_4d_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use geostrophicwind_mod, only : z2geouv_4d
    implicit none

    type(c_ptr), value, intent(in) :: z, glon, glat, ug, vg
    integer(c_int), value, intent(in) :: nlat, mlon, n3rd, n4th, iopt
    real(c_double), value, intent(in) :: zmsg

    integer :: nlat_f, mlon_f, n3rd_f, n4th_f
    real(c_double), pointer :: z_view(:, :, :, :)
    real(c_double), pointer :: glon_view(:)
    real(c_double), pointer :: glat_view(:)
    real(c_double), pointer :: ug_view(:, :, :, :)
    real(c_double), pointer :: vg_view(:, :, :, :)

    nlat_f = int(nlat)
    mlon_f = int(mlon)
    n3rd_f = int(n3rd)
    n4th_f = int(n4th)

    call c_f_pointer(z, z_view, [nlat_f, mlon_f, n3rd_f, n4th_f])
    call c_f_pointer(glon, glon_view, [mlon_f])
    call c_f_pointer(glat, glat_view, [nlat_f])
    call c_f_pointer(ug, ug_view, [nlat_f, mlon_f, n3rd_f, n4th_f])
    call c_f_pointer(vg, vg_view, [nlat_f, mlon_f, n3rd_f, n4th_f])

    call z2geouv_4d( &
        z_view, nlat_f, mlon_f, n3rd_f, n4th_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt) &
    )
end subroutine z2geouv_4d_c


subroutine zuvnew_c(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C, name="zuvnew_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use geostrophicwind_mod, only : zuvnew
    implicit none

    type(c_ptr), value, intent(in) :: z, glon, glat, ug, vg
    integer(c_int), value, intent(in) :: nlat, mlon, iopt
    real(c_double), value, intent(in) :: zmsg

    integer :: nlat_f, mlon_f
    real(c_double), pointer :: z_view(:, :)
    real(c_double), pointer :: glon_view(:)
    real(c_double), pointer :: glat_view(:)
    real(c_double), pointer :: ug_view(:, :)
    real(c_double), pointer :: vg_view(:, :)

    nlat_f = int(nlat)
    mlon_f = int(mlon)

    call c_f_pointer(z, z_view, [nlat_f, mlon_f])
    call c_f_pointer(glon, glon_view, [mlon_f])
    call c_f_pointer(glat, glat_view, [nlat_f])
    call c_f_pointer(ug, ug_view, [nlat_f, mlon_f])
    call c_f_pointer(vg, vg_view, [nlat_f, mlon_f])

    call zuvnew(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine zuvnew_c


subroutine z2guv_c(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C, name="z2guv_c")
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    use geostrophicwind_mod, only : z2guv
    implicit none

    type(c_ptr), value, intent(in) :: z, glon, glat, ug, vg
    integer(c_int), value, intent(in) :: nlat, mlon, iopt
    real(c_double), value, intent(in) :: zmsg

    integer :: nlat_f, mlon_f
    real(c_double), pointer :: z_view(:, :)
    real(c_double), pointer :: glon_view(:)
    real(c_double), pointer :: glat_view(:)
    real(c_double), pointer :: ug_view(:, :)
    real(c_double), pointer :: vg_view(:, :)

    nlat_f = int(nlat)
    mlon_f = int(mlon)

    call c_f_pointer(z, z_view, [nlat_f, mlon_f])
    call c_f_pointer(glon, glon_view, [mlon_f])
    call c_f_pointer(glat, glat_view, [nlat_f])
    call c_f_pointer(ug, ug_view, [nlat_f, mlon_f])
    call c_f_pointer(vg, vg_view, [nlat_f, mlon_f])

    call z2guv(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine z2guv_c
