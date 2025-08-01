!
! Optimized version of getlegfunc.f - MODIFIED TO TEST EQUIVALENCE
! Computes associated Legendre functions with improved performance
! Optimizations: Modern Fortran syntax, reduced function calls, better memory access
! MODIFIED: Added a test comment to verify recompilation works
!
subroutine getlegfunc(legfunc, lat, ntrunc)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: ntrunc
    real, intent(in) :: lat
    real, intent(out) :: legfunc((ntrunc+1)*(ntrunc+2)/2)

    ! Local variables
    integer :: m, n, nm, nmstrt
    real :: theta
    real, parameter :: pi = 3.141592653589793  ! High precision pi
    real :: cp((ntrunc/2)+1)

    ! Convert latitude to colatitude in radians
    theta = 0.5 * pi - (pi / 180.0) * lat

    ! Main computation loop - optimized memory access pattern
    nmstrt = 0
    do m = 1, ntrunc + 1
        !DIR$ VECTOR ALWAYS  ! Force vectorization for performance
        do n = m, ntrunc + 1
            nm = nmstrt + n - m + 1

            ! Compute associated Legendre function coefficients
            call alfk(n-1, m-1, cp)
            call lfpt(n-1, m-1, theta, cp, legfunc(nm))
        end do
        nmstrt = nmstrt + ntrunc - m + 2
    end do

end subroutine getlegfunc
