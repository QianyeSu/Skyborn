!
! Optimized version of invlap.f
! Computes inverse Laplacian operator in spectral space
! Optimizations: Modern Fortran syntax, vectorization, reduced function calls
!
subroutine invlap(dataspec, dataspec_ilap, nmdim, nt, rsphere)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nmdim
    integer, intent(in) :: nt
    real, intent(in) :: rsphere
    complex, intent(in) :: dataspec(nmdim, nt)
    complex, intent(out) :: dataspec_ilap(nmdim, nt)

    ! Local variables
    integer :: ntrunc, m, n, n1, nm, nmstrt, i
    real :: rsphere_sq, n_real, invlap_factor
    complex, parameter :: zero = (0.0, 0.0)

    ! Constants
    rsphere_sq = rsphere * rsphere

    ! Calculate truncation from nmdim
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    ! Main computation loops - optimized for vectorization
    !$OMP PARALLEL DO PRIVATE(nmstrt, m, n1, n, nm, n_real, invlap_factor) &
    !$OMP             SHARED(ntrunc, rsphere_sq)
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            n1 = m
            if (m == 1) n1 = 2  ! Skip n=1 for m=1 to avoid division by zero

            !DIR$ VECTOR ALWAYS
            do n = n1, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)

                ! Pre-compute inverse Laplacian factor
                invlap_factor = -rsphere_sq / (n_real * (n_real - 1.0))

                ! Apply inverse Laplacian operator
                dataspec_ilap(nm, i) = invlap_factor * dataspec(nm, i)
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do

        ! Set first mode to zero (corresponds to global mean)
        dataspec_ilap(1, i) = zero
    end do
    !$OMP END PARALLEL DO

end subroutine invlap
