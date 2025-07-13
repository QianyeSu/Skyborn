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
    !$OMP PARALLEL DO PRIVATE(nmstrt, m, n, nm) SHARED(ntrunc)
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !DIR$ VECTOR ALWAYS
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1

                ! Combine real and imaginary parts with scaling
                dataspec(nm, i) = scale * cmplx(a(m, n, i), b(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine twodtooned
