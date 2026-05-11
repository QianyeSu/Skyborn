!
! Optimized version of onedtotwod.f
! Converts spectral coefficients to grid coefficients
! Optimizations: Modern Fortran syntax, vectorization, reduced operations
!
subroutine onedtotwod(dataspec, a, b, nlat, nmdim, nt)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, nmdim, nt
    complex, intent(in) :: dataspec(nmdim, nt)
    real, intent(out) :: a(nlat, nlat, nt), b(nlat, nlat, nt)

    ! Local variables
    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale_inv
    complex :: spec_val

    ! Calculate truncation from nmdim
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    ! The direct C wrapper allocates these arrays without zero-initialization.
    ! SPHEREPACK synthesis expects the unused triangular slots to be zero.
    a(:, :, :) = 0.0
    b(:, :, :) = 0.0

    ! Pre-compute inverse scale factor
    scale_inv = 2.0  ! 1.0 / 0.5

    ! Main computation loops - optimized for cache efficiency
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, spec_val)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1

                ! Extract complex value once and scale
                spec_val = dataspec(nm, i) * scale_inv

                ! Split into real and imaginary parts
                a(m, n, i) = real(spec_val)
                b(m, n, i) = aimag(spec_val)
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod

!
! Converts vorticity and divergence spectral coefficients to vector-harmonic
! coefficients in one pass.
!
subroutine onedtotwod_vrtdiv(vrtspec, divspec, br, bi, cr, ci, &
                             nlat, nmdim, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: vrtspec(nmdim, nt), divspec(nmdim, nt)
    real, intent(out) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)
    real, intent(out) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)

    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: vrt_real, vrt_imag, div_real, div_imag

    scale = 0.5
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    br(:, :, :) = 0.0
    bi(:, :, :) = 0.0
    cr(:, :, :) = 0.0
    ci(:, :, :) = 0.0

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, div_real, div_imag, vrt_real, vrt_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))
                div_real = real(divspec(nm, i))
                div_imag = aimag(divspec(nm, i))
                vrt_real = real(vrtspec(nm, i))
                vrt_imag = aimag(vrtspec(nm, i))
                br(m, n, i) = -n_factor * div_real / scale
                bi(m, n, i) = -n_factor * div_imag / scale
                cr(m, n, i) = n_factor * vrt_real / scale
                ci(m, n, i) = n_factor * vrt_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_vrtdiv

!
! Converts vorticity spectral coefficients to vector-harmonic coefficients only.
!
subroutine onedtotwod_vrt(vrtspec, cr, ci, nlat, nmdim, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: vrtspec(nmdim, nt)
    real, intent(out) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)

    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: vrt_real, vrt_imag

    scale = 0.5
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    cr(:, :, :) = 0.0
    ci(:, :, :) = 0.0

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, vrt_real, vrt_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))
                vrt_real = real(vrtspec(nm, i))
                vrt_imag = aimag(vrtspec(nm, i))
                cr(m, n, i) = n_factor * vrt_real / scale
                ci(m, n, i) = n_factor * vrt_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_vrt

!
! Converts divergence spectral coefficients to vector-harmonic coefficients only.
!
subroutine onedtotwod_div(divspec, br, bi, nlat, nmdim, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: divspec(nmdim, nt)
    real, intent(out) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)

    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: div_real, div_imag

    scale = 0.5
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    br(:, :, :) = 0.0
    bi(:, :, :) = 0.0

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, div_real, div_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))
                div_real = real(divspec(nm, i))
                div_imag = aimag(divspec(nm, i))
                br(m, n, i) = -n_factor * div_real / scale
                bi(m, n, i) = -n_factor * div_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_div
