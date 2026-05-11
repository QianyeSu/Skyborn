!
! Optimized version of twodtooned.f
! Converts grid coefficients to spectral coefficients
! Optimizations: Modern Fortran syntax, vectorization, reduced operations
!
subroutine twodtooned(dataspec, a, b, nlat, ntrunc, nt)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: a(nlat, nlat, nt), b(nlat, nlat, nt)
    complex, intent(out) :: dataspec((ntrunc+1)*(ntrunc+2)/2, nt)

    ! Local variables
    integer :: m, n, nm, nmstrt, i
    real, parameter :: scale = 0.5

    ! Main computation loops - optimized for cache efficiency
    dataspec(:, :) = (0.0, 0.0)
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1

                ! Combine real and imaginary parts with scaling
                dataspec(nm, i) = scale * cmplx(a(m, n, i), b(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned

!
! Converts vector harmonic coefficients to vorticity and divergence spectral
! coefficients in one pass.
!
subroutine twodtooned_vrtdiv(vrtspec, divspec, br, bi, cr, ci, &
                             nlat, ntrunc, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: rsphere
    real, intent(in) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)
    real, intent(in) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)
    complex, intent(out) :: vrtspec((ntrunc+1)*(ntrunc+2)/2, nt)
    complex, intent(out) :: divspec((ntrunc+1)*(ntrunc+2)/2, nt)

    integer :: m, n, nm, nmstrt, i
    real :: scale, rsphere_inv, n_real, n_factor

    scale = 0.5
    rsphere_inv = 1.0 / rsphere

    vrtspec(:, :) = (0.0, 0.0)
    divspec(:, :) = (0.0, 0.0)

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = sqrt(n_real * (n_real - 1.0)) * rsphere_inv * scale
                divspec(nm, i) = -n_factor * cmplx(br(m, n, i), bi(m, n, i))
                vrtspec(nm, i) = n_factor * cmplx(cr(m, n, i), ci(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned_vrtdiv

!
! Converts vector harmonic coefficients to vorticity spectral coefficients only.
!
subroutine twodtooned_vrt(vrtspec, cr, ci, nlat, ntrunc, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: rsphere
    real, intent(in) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)
    complex, intent(out) :: vrtspec((ntrunc+1)*(ntrunc+2)/2, nt)

    integer :: m, n, nm, nmstrt, i
    real :: scale, rsphere_inv, n_real, n_factor

    scale = 0.5
    rsphere_inv = 1.0 / rsphere

    vrtspec(:, :) = (0.0, 0.0)

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = sqrt(n_real * (n_real - 1.0)) * rsphere_inv * scale
                vrtspec(nm, i) = n_factor * cmplx(cr(m, n, i), ci(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned_vrt

!
! Converts vector harmonic coefficients to divergence spectral coefficients only.
!
subroutine twodtooned_div(divspec, br, bi, nlat, ntrunc, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: rsphere
    real, intent(in) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)
    complex, intent(out) :: divspec((ntrunc+1)*(ntrunc+2)/2, nt)

    integer :: m, n, nm, nmstrt, i
    real :: scale, rsphere_inv, n_real, n_factor

    scale = 0.5
    rsphere_inv = 1.0 / rsphere

    divspec(:, :) = (0.0, 0.0)

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = sqrt(n_real * (n_real - 1.0)) * rsphere_inv * scale
                divspec(nm, i) = -n_factor * cmplx(br(m, n, i), bi(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned_div
