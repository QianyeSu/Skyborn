!
! Direct scalar transforms for packed reduced Gaussian grids.
!
! These helpers use the same public spectral coefficient convention as the
! existing Spharmt Gaussian backend, but replace the rectangular longitude
! FFT with row-wise Fourier sums over each reduced latitude circle.
!

subroutine reduced_gaussian_legendre_basis(basis, nlat, ntrunc, ierror)
    implicit none

    integer, intent(in) :: nlat, ntrunc
    real, intent(out) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    integer, intent(out) :: ierror

    integer :: j, m, n, nm, nmdim, ldwork, lw
    double precision :: theta(nlat), weights(nlat)
    double precision :: dwork(nlat * (nlat + 2))
    double precision :: cp(nlat + 1), pb

    external :: gaqd, dnlfk, dnlft

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2
    basis(:, :) = 0.0

    ierror = 0
    ldwork = nlat * (nlat + 2)
    lw = ldwork
    call gaqd(nlat, theta, weights, dwork, lw, ierror)
    if (ierror /= 0) return

    nm = 0
    do m = 0, ntrunc
        do n = m, ntrunc
            nm = nm + 1
            call dnlfk(m, n, cp)
            do j = 1, nlat
                call dnlft(m, n, theta(j), cp, pb)
                basis(j, nm) = real(pb)
            end do
        end do
    end do

end subroutine reduced_gaussian_legendre_basis


subroutine reduced_gaussian_legendre_derivative_basis(dbasis, nlat, ntrunc, ierror)
    implicit none

    integer, intent(in) :: nlat, ntrunc
    real, intent(out) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    integer, intent(out) :: ierror

    integer :: j, m, n, nm, nmdim, ldwork, lw
    double precision :: theta(nlat), weights(nlat)
    double precision :: dwork(nlat * (nlat + 2))
    double precision :: cp(nlat + 1), dpb

    external :: gaqd, dnlfk, dnlftd

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2
    dbasis(:, :) = 0.0

    ierror = 0
    ldwork = nlat * (nlat + 2)
    lw = ldwork
    call gaqd(nlat, theta, weights, dwork, lw, ierror)
    if (ierror /= 0) return

    nm = 0
    do m = 0, ntrunc
        do n = m, ntrunc
            nm = nm + 1
            call dnlfk(m, n, cp)
            do j = 1, nlat
                call dnlftd(m, n, theta(j), cp, dpb)
                dbasis(j, nm) = real(dpb)
            end do
        end do
    end do

end subroutine reduced_gaussian_legendre_derivative_basis


subroutine reduced_gaussian_grdtospec(datagrid, pl, weights, basis, dataspec, &
                                      ngptot, nlat, ntrunc, nt, ierror)
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: datagrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    complex, intent(out) :: dataspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total
    double precision :: angle, angle_step, inv_nlon
    double precision :: cos_angle, sin_angle
    double precision :: coeff_real, coeff_imag, value, weighted_basis
    double precision, parameter :: twopi = 6.283185307179586476925286766559d0

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nt < 1) return

    total = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        total = total + pl(j)
    end do

    ierror = 4
    if (total /= ngptot) return

    dataspec(:, :) = (0.0, 0.0)
    offset = 0

    do j = 1, nlat
        row_nlon = pl(j)
        inv_nlon = 1.0d0 / dble(row_nlon)

        do m = 0, ntrunc
            if (m == 0) then
                do k = 1, nt
                    coeff_real = 0.0d0
                    do i = 1, row_nlon
                        coeff_real = coeff_real + dble(datagrid(offset + i, k))
                    end do
                    coeff_real = coeff_real * inv_nlon

                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        weighted_basis = dble(weights(j)) * dble(basis(j, nm))
                        dataspec(nm, k) = dataspec(nm, k) + &
                            cmplx(real(coeff_real * weighted_basis), 0.0)
                    end do
                end do
            else if (2 * m <= row_nlon) then
                angle_step = twopi * dble(m) * inv_nlon
                do k = 1, nt
                    coeff_real = 0.0d0
                    coeff_imag = 0.0d0
                    do i = 1, row_nlon
                        angle = angle_step * dble(i - 1)
                        cos_angle = cos(angle)
                        sin_angle = sin(angle)
                        value = dble(datagrid(offset + i, k))
                        coeff_real = coeff_real + value * cos_angle
                        coeff_imag = coeff_imag - value * sin_angle
                    end do
                    coeff_real = coeff_real * inv_nlon
                    coeff_imag = coeff_imag * inv_nlon

                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        weighted_basis = dble(weights(j)) * dble(basis(j, nm))
                        dataspec(nm, k) = dataspec(nm, k) + &
                            cmplx(real(coeff_real * weighted_basis), &
                                  real(coeff_imag * weighted_basis))
                    end do
                end do
            end if
        end do

        offset = offset + row_nlon
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_grdtospec


subroutine reduced_gaussian_spectogrd(dataspec, pl, basis, datagrid, &
                                      nmdim, nlat, ntrunc, nt, ngptot, ierror)
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    complex, intent(in) :: dataspec(nmdim, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(out) :: datagrid(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total
    double precision :: angle, angle_step, value
    double precision :: scrm_real, scrm_imag
    double precision, parameter :: twopi = 6.283185307179586476925286766559d0

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    total = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        total = total + pl(j)
    end do

    ierror = 5
    if (total /= ngptot) return

    offset = 0
    do j = 1, nlat
        row_nlon = pl(j)
        do k = 1, nt
            do i = 1, row_nlon
                value = 0.0d0
                do m = 0, ntrunc
                    scrm_real = 0.0d0
                    scrm_imag = 0.0d0
                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        scrm_real = scrm_real + dble(real(dataspec(nm, k))) * &
                            dble(basis(j, nm))
                        scrm_imag = scrm_imag + dble(aimag(dataspec(nm, k))) * &
                            dble(basis(j, nm))
                    end do

                    if (m == 0) then
                        value = value + scrm_real
                    else
                        angle_step = twopi * dble(m) / dble(row_nlon)
                        angle = angle_step * dble(i - 1)
                        value = value + 2.0d0 * &
                            (scrm_real * cos(angle) - scrm_imag * sin(angle))
                    end if
                end do
                datagrid(offset + i, k) = real(value)
            end do
        end do
        offset = offset + row_nlon
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_spectogrd


subroutine reduced_gaussian_getgrad(dataspec, pl, basis, dbasis, sin_theta, &
                                    ugrad, vgrad, nmdim, nlat, ntrunc, nt, &
                                    ngptot, rsphere, ierror)
    implicit none

    integer, intent(in) :: nmdim, nlat, ntrunc, nt, ngptot
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: basis(nlat, nmdim)
    real, intent(in) :: dbasis(nlat, nmdim)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(in) :: dataspec(nmdim, nt)
    real, intent(out) :: ugrad(ngptot, nt)
    real, intent(out) :: vgrad(ngptot, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total
    double precision :: angle, angle_step, inv_r, inv_r_sin
    double precision :: cos_angle, sin_angle, cos_step, sin_step, cos_next
    double precision :: spec_real, spec_imag, pval, dpval
    double precision :: uvalue, vvalue
    double precision :: sum_p_real, sum_p_imag, sum_dp_real, sum_dp_imag
    double precision :: u_cos(0:ntrunc, nt), u_sin(0:ntrunc, nt)
    double precision :: v_cos(0:ntrunc, nt), v_sin(0:ntrunc, nt)
    double precision, parameter :: twopi = 6.283185307179586476925286766559d0

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nmdim /= (ntrunc + 1) * (ntrunc + 2) / 2) return

    ierror = 4
    if (nt < 1) return

    ierror = 5
    if (rsphere <= 0.0) return

    total = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
    end do

    ierror = 6
    if (total /= ngptot) return

    inv_r = 1.0d0 / dble(rsphere)
    ugrad(:, :) = 0.0
    vgrad(:, :) = 0.0

    offset = 0
    do j = 1, nlat
        row_nlon = pl(j)
        inv_r_sin = inv_r / dble(sin_theta(j))

        u_cos(:, :) = 0.0d0
        u_sin(:, :) = 0.0d0
        v_cos(:, :) = 0.0d0
        v_sin(:, :) = 0.0d0

        do k = 1, nt
            do m = 0, ntrunc
                sum_p_real = 0.0d0
                sum_p_imag = 0.0d0
                sum_dp_real = 0.0d0
                sum_dp_imag = 0.0d0
                nm = spectral_offset(m, ntrunc)
                do n = m, ntrunc
                    nm = nm + 1
                    spec_real = dble(real(dataspec(nm, k)))
                    spec_imag = dble(aimag(dataspec(nm, k)))
                    pval = dble(basis(j, nm))
                    dpval = dble(dbasis(j, nm))
                    sum_p_real = sum_p_real + spec_real * pval
                    sum_p_imag = sum_p_imag + spec_imag * pval
                    sum_dp_real = sum_dp_real + spec_real * dpval
                    sum_dp_imag = sum_dp_imag + spec_imag * dpval
                end do

                if (m == 0) then
                    v_cos(m, k) = -sum_dp_real * inv_r
                else
                    u_cos(m, k) = -2.0d0 * dble(m) * sum_p_imag * inv_r_sin
                    u_sin(m, k) = -2.0d0 * dble(m) * sum_p_real * inv_r_sin
                    v_cos(m, k) = -2.0d0 * sum_dp_real * inv_r
                    v_sin(m, k) = 2.0d0 * sum_dp_imag * inv_r
                end if
            end do
        end do

        do k = 1, nt
            do i = 1, row_nlon
                ugrad(offset + i, k) = 0.0
                vgrad(offset + i, k) = real(v_cos(0, k))
            end do

            do m = 1, ntrunc
                angle_step = twopi * dble(m) / dble(row_nlon)
                cos_step = cos(angle_step)
                sin_step = sin(angle_step)
                cos_angle = 1.0d0
                sin_angle = 0.0d0
                do i = 1, row_nlon
                    uvalue = u_cos(m, k) * cos_angle + &
                        u_sin(m, k) * sin_angle
                    vvalue = v_cos(m, k) * cos_angle + &
                        v_sin(m, k) * sin_angle
                    ugrad(offset + i, k) = ugrad(offset + i, k) + real(uvalue)
                    vgrad(offset + i, k) = vgrad(offset + i, k) + real(vvalue)
                    cos_next = cos_angle * cos_step - sin_angle * sin_step
                    sin_angle = sin_angle * cos_step + cos_angle * sin_step
                    cos_angle = cos_next
                end do
            end do
        end do
        offset = offset + row_nlon
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getgrad


subroutine reduced_gaussian_getvrtdivspec(ugrid, vgrid, pl, weights, basis, &
                                          dbasis, sin_theta, vrtspec, divspec, &
                                          ngptot, nlat, ntrunc, nt, rsphere, &
                                          ierror)
    implicit none

    integer, intent(in) :: ngptot, nlat, ntrunc, nt
    real, intent(in) :: ugrid(ngptot, nt)
    real, intent(in) :: vgrid(ngptot, nt)
    integer, intent(in) :: pl(nlat)
    real, intent(in) :: weights(nlat)
    real, intent(in) :: basis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: dbasis(nlat, (ntrunc + 1) * (ntrunc + 2) / 2)
    real, intent(in) :: sin_theta(nlat)
    real, intent(in) :: rsphere
    complex, intent(out) :: vrtspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    complex, intent(out) :: divspec((ntrunc + 1) * (ntrunc + 2) / 2, nt)
    integer, intent(out) :: ierror

    integer :: j, i, k, m, n, nm, offset, row_nlon, total
    double precision :: angle, angle_step, inv_nlon, inv_r, inv_r_sin
    double precision :: cos_angle, sin_angle, cos_step, sin_step, cos_next
    double precision :: uvalue, vvalue
    double precision :: u_mean, v_mean, u_cos, u_sin, v_cos, v_sin
    double precision :: pval, dpval, wval, div_real, div_imag, vrt_real, vrt_imag
    double precision, parameter :: twopi = 6.283185307179586476925286766559d0

    ierror = 1
    if (nlat < 3) return

    ierror = 2
    if (ntrunc < 0 .or. ntrunc > nlat - 1) return

    ierror = 3
    if (nt < 1) return

    ierror = 4
    if (rsphere <= 0.0) return

    total = 0
    do j = 1, nlat
        if (pl(j) < 1) return
        if (sin_theta(j) <= 0.0) return
        total = total + pl(j)
    end do

    ierror = 5
    if (total /= ngptot) return

    inv_r = 1.0d0 / dble(rsphere)
    vrtspec(:, :) = (0.0, 0.0)
    divspec(:, :) = (0.0, 0.0)
    offset = 0

    do j = 1, nlat
        row_nlon = pl(j)
        inv_nlon = 1.0d0 / dble(row_nlon)
        inv_r_sin = inv_r / dble(sin_theta(j))
        wval = dble(weights(j))

        do m = 0, ntrunc
            if (m == 0) then
                do k = 1, nt
                    u_mean = 0.0d0
                    v_mean = 0.0d0
                    do i = 1, row_nlon
                        u_mean = u_mean + dble(ugrid(offset + i, k))
                        v_mean = v_mean + dble(vgrid(offset + i, k))
                    end do
                    u_mean = u_mean * inv_nlon
                    v_mean = v_mean * inv_nlon

                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        dpval = dble(dbasis(j, nm))
                        div_real = wval * dpval * inv_r * v_mean
                        vrt_real = -wval * dpval * inv_r * u_mean
                        divspec(nm, k) = divspec(nm, k) + cmplx(real(div_real), 0.0)
                        vrtspec(nm, k) = vrtspec(nm, k) + cmplx(real(vrt_real), 0.0)
                    end do
                end do
            else if (2 * m <= row_nlon) then
                angle_step = twopi * dble(m) * inv_nlon
                cos_step = cos(angle_step)
                sin_step = sin(angle_step)
                do k = 1, nt
                    u_cos = 0.0d0
                    u_sin = 0.0d0
                    v_cos = 0.0d0
                    v_sin = 0.0d0
                    cos_angle = 1.0d0
                    sin_angle = 0.0d0
                    do i = 1, row_nlon
                        uvalue = dble(ugrid(offset + i, k))
                        vvalue = dble(vgrid(offset + i, k))
                        u_cos = u_cos + uvalue * cos_angle
                        u_sin = u_sin + uvalue * sin_angle
                        v_cos = v_cos + vvalue * cos_angle
                        v_sin = v_sin + vvalue * sin_angle
                        cos_next = cos_angle * cos_step - sin_angle * sin_step
                        sin_angle = sin_angle * cos_step + cos_angle * sin_step
                        cos_angle = cos_next
                    end do
                    u_cos = u_cos * inv_nlon
                    u_sin = u_sin * inv_nlon
                    v_cos = v_cos * inv_nlon
                    v_sin = v_sin * inv_nlon

                    nm = spectral_offset(m, ntrunc)
                    do n = m, ntrunc
                        nm = nm + 1
                        pval = dble(basis(j, nm))
                        dpval = dble(dbasis(j, nm))
                        div_real = wval * (dble(m) * pval * inv_r_sin * u_sin + &
                                           dpval * inv_r * v_cos)
                        div_imag = wval * (dble(m) * pval * inv_r_sin * u_cos - &
                                           dpval * inv_r * v_sin)
                        vrt_real = wval * (-dpval * inv_r * u_cos + &
                                           dble(m) * pval * inv_r_sin * v_sin)
                        vrt_imag = wval * (dpval * inv_r * u_sin + &
                                           dble(m) * pval * inv_r_sin * v_cos)
                        divspec(nm, k) = divspec(nm, k) + &
                            cmplx(real(div_real), real(div_imag))
                        vrtspec(nm, k) = vrtspec(nm, k) + &
                            cmplx(real(vrt_real), real(vrt_imag))
                    end do
                end do
            end if
        end do

        offset = offset + row_nlon
    end do

    ierror = 0

contains

    integer function spectral_offset(m_value, trunc)
        integer, intent(in) :: m_value, trunc
        spectral_offset = m_value * (2 * trunc - m_value + 3) / 2
    end function spectral_offset

end subroutine reduced_gaussian_getvrtdivspec
