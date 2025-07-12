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
    ntrunc = int(-1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim))))

    ! Pre-compute inverse scale factor
    scale_inv = 2.0  ! 1.0 / 0.5

    ! Main computation loops - optimized for cache efficiency
    !$OMP PARALLEL DO PRIVATE(nmstrt, m, n, nm, spec_val) &
    !$OMP             SHARED(ntrunc)
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !DIR$ VECTOR ALWAYS
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
    !$OMP END PARALLEL DO

end subroutine onedtotwod
