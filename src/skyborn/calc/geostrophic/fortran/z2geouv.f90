module geostrophicwind_mod
    use, intrinsic :: iso_c_binding, only : c_double, c_f_pointer, c_int, c_ptr
    implicit none

contains

subroutine z2geouv(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C)
    type(c_ptr), value, intent(in) :: z
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: mlon
    real(c_double), value, intent(in) :: zmsg
    type(c_ptr), value, intent(in) :: glon
    type(c_ptr), value, intent(in) :: glat
    type(c_ptr), value, intent(in) :: ug
    type(c_ptr), value, intent(in) :: vg
    integer(c_int), value, intent(in) :: iopt

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

    call z2geouv_impl(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine z2geouv

subroutine z2geouv_3d(z, nlat, mlon, n3rd, zmsg, glon, glat, ug, vg, iopt) bind(C)
    type(c_ptr), value, intent(in) :: z
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: mlon
    integer(c_int), value, intent(in) :: n3rd
    real(c_double), value, intent(in) :: zmsg
    type(c_ptr), value, intent(in) :: glon
    type(c_ptr), value, intent(in) :: glat
    type(c_ptr), value, intent(in) :: ug
    type(c_ptr), value, intent(in) :: vg
    integer(c_int), value, intent(in) :: iopt

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

    call z2geouv_3d_impl( &
        z_view, nlat_f, mlon_f, n3rd_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt) &
    )
end subroutine z2geouv_3d

subroutine z2geouv_4d(z, nlat, mlon, n3rd, n4th, zmsg, glon, glat, ug, vg, iopt) bind(C)
    type(c_ptr), value, intent(in) :: z
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: mlon
    integer(c_int), value, intent(in) :: n3rd
    integer(c_int), value, intent(in) :: n4th
    real(c_double), value, intent(in) :: zmsg
    type(c_ptr), value, intent(in) :: glon
    type(c_ptr), value, intent(in) :: glat
    type(c_ptr), value, intent(in) :: ug
    type(c_ptr), value, intent(in) :: vg
    integer(c_int), value, intent(in) :: iopt

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

    call z2geouv_4d_impl( &
        z_view, nlat_f, mlon_f, n3rd_f, n4th_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt) &
    )
end subroutine z2geouv_4d

subroutine zuvnew(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C)
    type(c_ptr), value, intent(in) :: z
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: mlon
    real(c_double), value, intent(in) :: zmsg
    type(c_ptr), value, intent(in) :: glon
    type(c_ptr), value, intent(in) :: glat
    type(c_ptr), value, intent(in) :: ug
    type(c_ptr), value, intent(in) :: vg
    integer(c_int), value, intent(in) :: iopt

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

    call zuvnew_impl(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine zuvnew

subroutine z2guv(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt) bind(C)
    type(c_ptr), value, intent(in) :: z
    integer(c_int), value, intent(in) :: nlat
    integer(c_int), value, intent(in) :: mlon
    real(c_double), value, intent(in) :: zmsg
    type(c_ptr), value, intent(in) :: glon
    type(c_ptr), value, intent(in) :: glat
    type(c_ptr), value, intent(in) :: ug
    type(c_ptr), value, intent(in) :: vg
    integer(c_int), value, intent(in) :: iopt

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

    call z2guv_impl(z_view, nlat_f, mlon_f, zmsg, glon_view, glat_view, ug_view, vg_view, int(iopt))
end subroutine z2guv

subroutine z2geouv_impl(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    integer, intent(in) :: nlat
    integer, intent(in) :: mlon
    integer, intent(in) :: iopt
    real(c_double), intent(in) :: z(nlat, mlon)
    real(c_double), intent(in) :: zmsg
    real(c_double), intent(in) :: glat(nlat)
    real(c_double), intent(in) :: glon(mlon)
    real(c_double), intent(out) :: ug(nlat, mlon)
    real(c_double), intent(out) :: vg(nlat, mlon)

    if (nlat < 2 .or. mlon < 2) then
        ug = zmsg
        vg = zmsg
        return
    end if

    if (glat(2) > glat(1)) then
        call z2guv_impl(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    else
        call zuvnew_impl(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    end if
end subroutine z2geouv_impl

subroutine z2geouv_3d_impl(z, nlat, mlon, n3rd, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    integer, intent(in) :: nlat
    integer, intent(in) :: mlon
    integer, intent(in) :: n3rd
    integer, intent(in) :: iopt
    real(c_double), intent(in) :: z(nlat, mlon, n3rd)
    real(c_double), intent(in) :: zmsg
    real(c_double), intent(in) :: glat(nlat)
    real(c_double), intent(in) :: glon(mlon)
    real(c_double), intent(out) :: ug(nlat, mlon, n3rd)
    real(c_double), intent(out) :: vg(nlat, mlon, n3rd)

    integer :: k

    !$omp parallel do private(k) if(n3rd > 4)
    do k = 1, n3rd
        call z2geouv_impl(z(:, :, k), nlat, mlon, zmsg, glon, glat, ug(:, :, k), vg(:, :, k), iopt)
    end do
    !$omp end parallel do
end subroutine z2geouv_3d_impl

subroutine z2geouv_4d_impl(z, nlat, mlon, n3rd, n4th, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    integer, intent(in) :: nlat
    integer, intent(in) :: mlon
    integer, intent(in) :: n3rd
    integer, intent(in) :: n4th
    integer, intent(in) :: iopt
    real(c_double), intent(in) :: z(nlat, mlon, n3rd, n4th)
    real(c_double), intent(in) :: zmsg
    real(c_double), intent(in) :: glat(nlat)
    real(c_double), intent(in) :: glon(mlon)
    real(c_double), intent(out) :: ug(nlat, mlon, n3rd, n4th)
    real(c_double), intent(out) :: vg(nlat, mlon, n3rd, n4th)

    integer :: k, t

    !$omp parallel do private(k, t) collapse(2) if(n3rd * n4th > 16)
    do t = 1, n4th
        do k = 1, n3rd
            call z2geouv_impl( &
                z(:, :, k, t), nlat, mlon, zmsg, glon, glat, ug(:, :, k, t), vg(:, :, k, t), iopt &
            )
        end do
    end do
    !$omp end parallel do
end subroutine z2geouv_4d_impl

subroutine zuvnew_impl(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    integer, intent(in) :: nlat
    integer, intent(in) :: mlon
    integer, intent(in) :: iopt
    real(c_double), intent(in) :: z(nlat, mlon)
    real(c_double), intent(in) :: zmsg
    real(c_double), intent(in) :: glat(nlat)
    real(c_double), intent(in) :: glon(mlon)
    real(c_double), intent(out) :: ug(nlat, mlon)
    real(c_double), intent(out) :: vg(nlat, mlon)

    real(c_double) :: ztmp(nlat, mlon), glattmp(nlat)
    real(c_double) :: utmp(nlat), vtmp(nlat)
    integer :: nl, ml, nlat1

    nlat1 = nlat + 1

    do nl = 1, nlat
        glattmp(nlat1 - nl) = glat(nl)
        !$omp simd
        do ml = 1, mlon
            ztmp(nlat1 - nl, ml) = z(nl, ml)
        end do
        !$omp end simd
    end do

    call z2guv_impl(ztmp, nlat, mlon, zmsg, glon, glattmp, ug, vg, iopt)

    do ml = 1, mlon
        !$omp simd
        do nl = 1, nlat
            utmp(nlat1 - nl) = ug(nl, ml)
            vtmp(nlat1 - nl) = vg(nl, ml)
        end do
        !$omp end simd
        !$omp simd
        do nl = 1, nlat
            ug(nl, ml) = utmp(nl)
            vg(nl, ml) = vtmp(nl)
        end do
        !$omp end simd
    end do
end subroutine zuvnew_impl

subroutine z2guv_impl(z, nlat, mlon, zmsg, glon, glat, ug, vg, iopt)
    implicit none

    integer, intent(in) :: nlat
    integer, intent(in) :: mlon
    integer, intent(in) :: iopt
    real(c_double), intent(in) :: z(nlat, mlon)
    real(c_double), intent(in) :: zmsg
    real(c_double), intent(in) :: glat(nlat)
    real(c_double), intent(in) :: glon(mlon)
    real(c_double), intent(out) :: ug(nlat, mlon)
    real(c_double), intent(out) :: vg(nlat, mlon)

    real(c_double) :: zz(0:nlat + 1, 0:mlon + 1)
    real(c_double) :: gf(nlat), dx(nlat), dy(nlat)
    real(c_double) :: inv_dx(nlat), inv_dy(nlat)
    real(c_double), parameter :: gr = 9.80616_c_double
    real(c_double), parameter :: omega = 7.292e-5_c_double
    real(c_double), parameter :: re = 6371220.0_c_double
    real(c_double), parameter :: pi = 4.0_c_double * atan(1.0_c_double)
    real(c_double), parameter :: rad = pi / 180.0_c_double
    real(c_double), parameter :: gr_over_two_omega = gr / (2.0_c_double * omega)
    real(c_double) :: rr, dlon, dlon2
    real(c_double) :: sin_lat, cos_lat, abs_lat
    integer :: nl, ml
    logical :: valid_lat(nlat)

    ug = zmsg
    vg = zmsg

    if (nlat < 2 .or. mlon < 2) then
        return
    end if

    rr = re * rad
    dlon = rr * (glon(2) - glon(1))
    dlon2 = 2.0_c_double * dlon

    !$omp simd
    do nl = 1, nlat
        abs_lat = abs(glat(nl))
        valid_lat(nl) = (glat(nl) /= 0.0_c_double) .and. (abs_lat <= 89.999_c_double)

        if (valid_lat(nl)) then
            sin_lat = sin(glat(nl) * rad)
            cos_lat = cos(glat(nl) * rad)
            gf(nl) = gr_over_two_omega / sin_lat
            dx(nl) = dlon2 * cos_lat
            if (dx(nl) /= 0.0_c_double) then
                inv_dx(nl) = 1.0_c_double / dx(nl)
            else
                inv_dx(nl) = 0.0_c_double
            end if
        else
            gf(nl) = zmsg
            dx(nl) = 0.0_c_double
            inv_dx(nl) = 0.0_c_double
        end if
    end do
    !$omp end simd

    if (glat(1) < -89.999_c_double) then
        dx(1) = 0.0_c_double
        inv_dx(1) = 0.0_c_double
    end if
    if (glat(nlat) > 89.999_c_double) then
        dx(nlat) = 0.0_c_double
        inv_dx(nlat) = 0.0_c_double
    end if

    dy(1) = (glat(2) - glat(1)) * rr
    inv_dy(1) = 1.0_c_double / dy(1)

    !$omp simd
    do nl = 2, nlat - 1
        dy(nl) = (glat(nl + 1) - glat(nl - 1)) * rr
        inv_dy(nl) = 1.0_c_double / dy(nl)
    end do
    !$omp end simd

    dy(nlat) = (glat(nlat) - glat(nlat - 1)) * rr
    inv_dy(nlat) = 1.0_c_double / dy(nlat)

    do nl = 1, nlat
        !$omp simd
        do ml = 1, mlon
            zz(nl, ml) = z(nl, ml)
        end do
        !$omp end simd
    end do

    !$omp simd
    do ml = 1, mlon
        zz(0, ml) = zz(1, ml)
        zz(nlat + 1, ml) = zz(nlat, ml)
    end do
    !$omp end simd

    if (iopt == 1) then
        !$omp simd
        do nl = 0, nlat + 1
            zz(nl, 0) = zz(nl, mlon)
            zz(nl, mlon + 1) = zz(nl, 1)
        end do
        !$omp end simd
    else
        !$omp simd
        do nl = 0, nlat + 1
            zz(nl, 0) = zz(nl, 1)
            zz(nl, mlon + 1) = zz(nl, mlon)
        end do
        !$omp end simd
    end if

    do nl = 1, nlat
        if (valid_lat(nl)) then
            !$omp simd
            do ml = 1, mlon
                if (zz(nl + 1, ml) /= zmsg .and. zz(nl - 1, ml) /= zmsg .and. &
                    zz(nl, ml + 1) /= zmsg .and. zz(nl, ml - 1) /= zmsg) then
                    ug(nl, ml) = -gf(nl) * (zz(nl + 1, ml) - zz(nl - 1, ml)) * inv_dy(nl)
                    vg(nl, ml) = gf(nl) * (zz(nl, ml + 1) - zz(nl, ml - 1)) * inv_dx(nl)
                else
                    if (zz(nl + 1, ml) /= zmsg .and. zz(nl - 1, ml) /= zmsg) then
                        ug(nl, ml) = -gf(nl) * (zz(nl + 1, ml) - zz(nl - 1, ml)) * inv_dy(nl)
                    end if
                    if (zz(nl, ml + 1) /= zmsg .and. zz(nl, ml - 1) /= zmsg) then
                        vg(nl, ml) = gf(nl) * (zz(nl, ml + 1) - zz(nl, ml - 1)) * inv_dx(nl)
                    end if
                end if
            end do
            !$omp end simd
        end if
    end do

    if (iopt /= 1) then
        !$omp simd
        do nl = 1, nlat
            if (vg(nl, 1) /= zmsg) vg(nl, 1) = vg(nl, 1) * 2.0_c_double
            if (vg(nl, mlon) /= zmsg) vg(nl, mlon) = vg(nl, mlon) * 2.0_c_double
        end do
        !$omp end simd
    end if
end subroutine z2guv_impl

end module geostrophicwind_mod
