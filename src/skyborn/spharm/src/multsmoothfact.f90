!
! Optimized version of multsmoothfact.f
! Applies smoothing factors to spectral coefficients
! Optimizations: Modern Fortran syntax, vectorization, cache-friendly access
!
subroutine multsmoothfact(dataspec, dataspec_smooth, smooth, nlat, nmdim, nt)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nlat, nmdim, nt
    complex, intent(in) :: dataspec(nmdim, nt)
    real, intent(in) :: smooth(nlat)
    complex, intent(out) :: dataspec_smooth(nmdim, nt)

    ! Local variables
    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: smooth_factor

    ! Calculate truncation from nmdim
    ntrunc = int(-1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim))))

    ! Main computation loops - optimized for vectorization and cache efficiency
    !$OMP PARALLEL DO PRIVATE(nmstrt, m, n, nm, smooth_factor) &
    !$OMP             SHARED(ntrunc, smooth)
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !DIR$ VECTOR ALWAYS
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                smooth_factor = smooth(n)

                ! Apply smoothing factor
                dataspec_smooth(nm, i) = dataspec(nm, i) * smooth_factor
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do
    !$OMP END PARALLEL DO

end subroutine multsmoothfact
