!
! Optimized version of twodtooned_vrtdiv.f
! Converts grid coefficients to spectral coefficients for vorticity and divergence
! Optimizations: Modern Fortran syntax, vectorization, reduced function calls
!
subroutine twodtooned_vrtdiv(vrtspec, divspec, br, bi, cr, ci, &
                            nlat, ntrunc, nt, rsphere)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: rsphere
    real, intent(in) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)
    real, intent(in) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)
    complex, intent(out) :: vrtspec((ntrunc+1)*(ntrunc+2)/2, nt)
    complex, intent(out) :: divspec((ntrunc+1)*(ntrunc+2)/2, nt)

    ! Local variables
    integer :: m, n, nm, nmstrt, i
    real :: scale, rsphere_inv, n_real, n_factor

    ! Constants
    scale = 0.5
    rsphere_inv = 1.0 / rsphere

    ! Main computation loops - optimized for cache efficiency
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)

                ! Pre-compute common factor
                n_factor = sqrt(n_real * (n_real - 1.0)) * rsphere_inv * scale

                ! Compute spectral coefficients with optimized operations
                divspec(nm, i) = -n_factor * cmplx(br(m, n, i), bi(m, n, i))
                vrtspec(nm, i) = n_factor * cmplx(cr(m, n, i), ci(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned_vrtdiv
