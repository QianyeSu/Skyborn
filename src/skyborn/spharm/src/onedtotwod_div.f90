!
! Converts divergence spectral coefficients to vector-harmonic coefficients only.
!
subroutine onedtotwod_div(divspec, br, bi, nlat, nmdim, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: divspec(nmdim, nt)
    real, intent(out) :: br(nlat, nlat, nt), bi(nlat, nlat, nt)

    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: div_real, div_imag

    scale = 0.5
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, div_real, div_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))
                div_real = real(divspec(nm, i))
                div_imag = aimag(divspec(nm, i))
                br(m, n, i) = -n_factor * div_real / scale
                bi(m, n, i) = -n_factor * div_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_div
