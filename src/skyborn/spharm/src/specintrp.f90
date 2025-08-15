!
! Optimized version of specintrp.f
! Spectral interpolation function
! Optimizations: Modern Fortran syntax, vectorization, trigonometric optimization
!
subroutine specintrp(rlon, ntrunc, datnm, scrm, pnm, ob)
    implicit none

    ! Input/Output parameters
    integer, intent(in) :: ntrunc
    real, intent(in) :: rlon
    complex, intent(in) :: datnm((ntrunc+1)*(ntrunc+2)/2)
    real, intent(in) :: pnm((ntrunc+1)*(ntrunc+2)/2)
    complex, intent(out) :: scrm(ntrunc+1)
    real, intent(out) :: ob

    ! Local variables
    integer :: m, n, nm, nmstrt, mwaves
    real :: m_rlon, cos_m_rlon, sin_m_rlon
    complex, parameter :: zero = (0.0, 0.0)

    mwaves = ntrunc + 1

    ! First pass: compute Fourier coefficients in latitude direction
    nmstrt = 0
    do m = 1, mwaves
        scrm(m) = zero

        !$OMP SIMD PRIVATE(nm)
        do n = 1, mwaves - m + 1
            nm = nmstrt + n
            scrm(m) = scrm(m) + datnm(nm) * pnm(nm)
        end do
        nmstrt = nmstrt + mwaves - m + 1
    end do

    ! Second pass: sum Fourier series in longitude direction
    ! Start with m=1 term (no longitude dependence)
    ob = real(scrm(1))

    ! Add m>=2 terms with optimized trigonometric evaluation
    !$OMP SIMD PRIVATE(m_rlon, cos_m_rlon, sin_m_rlon)
    do m = 2, mwaves
        m_rlon = real(m-1) * rlon
        cos_m_rlon = cos(m_rlon)
        sin_m_rlon = sin(m_rlon)

        ob = ob + 2.0 * (real(scrm(m)) * cos_m_rlon - aimag(scrm(m)) * sin_m_rlon)
    end do

end subroutine specintrp
