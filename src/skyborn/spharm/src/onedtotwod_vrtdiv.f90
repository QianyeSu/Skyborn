!
! Optimized version of onedtotwod_vrtdiv.f
! Converts spectral coefficients to grid coefficients for vorticity and divergence
! Optimizations: Modern Fortran syntax, vectorization, reduced function calls
!
subroutine onedtotwod_vrtdiv(vrtspec, divspec, br, bi, cr, ci, &
                            nlat, nmdim, nt, rsphere)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: vrtspec(nmdim, nt), divspec(nmdim, nt)
    real, intent(out) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)
    real, intent(out) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)

    ! Local variables
    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: vrt_real, vrt_imag, div_real, div_imag

    ! Constants
    scale = 0.5

    ! Calculate truncation from nmdim
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    ! Main computation loops - optimized order for cache efficiency
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, div_real, div_imag, vrt_real, vrt_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)

                ! Pre-compute common factor - CORRECTED to match .f version
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))

                ! Extract real and imaginary parts once
                div_real = real(divspec(nm, i))
                div_imag = aimag(divspec(nm, i))
                vrt_real = real(vrtspec(nm, i))
                vrt_imag = aimag(vrtspec(nm, i))

                ! Compute coefficients with corrected operations
                br(m, n, i) = -n_factor * div_real / scale
                bi(m, n, i) = -n_factor * div_imag / scale
                cr(m, n, i) = n_factor * vrt_real / scale
                ci(m, n, i) = n_factor * vrt_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_vrtdiv
