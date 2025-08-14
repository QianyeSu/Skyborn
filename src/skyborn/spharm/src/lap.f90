!
! Optimized version of lap.f
! Computes Laplacian operator in spectral space
! Optimizations: Modern Fortran syntax, vectorization, reduced function calls
!
subroutine lap(dataspec, dataspec_lap, nmdim, nt, rsphere)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: dataspec(nmdim, nt)
    complex, intent(out) :: dataspec_lap(nmdim, nt)

    ! Local variables
    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: rsphere_inv_sq, n_real, lap_factor

    ! Constants
    rsphere_inv_sq = 1.0 / (rsphere * rsphere)

    ! Calculate truncation from nmdim
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    ! Main computation loops - optimized for vectorization
    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !DIR$ VECTOR ALWAYS
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)

                ! Pre-compute Laplacian factor
                lap_factor = -n_real * (n_real - 1.0) * rsphere_inv_sq

                ! Apply Laplacian operator
                dataspec_lap(nm, i) = lap_factor * dataspec(nm, i)
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine lap
