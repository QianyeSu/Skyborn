!
! Converts vorticity spectral coefficients to vector-harmonic coefficients only.
!
subroutine onedtotwod_vrt(vrtspec, cr, ci, nlat, nmdim, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, nmdim, nt
    real, intent(in) :: rsphere
    complex, intent(in) :: vrtspec(nmdim, nt)
    real, intent(out) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)

    integer :: ntrunc, m, n, nm, nmstrt, i
    real :: scale, n_real, n_factor
    real :: vrt_real, vrt_imag

    scale = 0.5
    ntrunc = -1.5 + 0.5 * sqrt(9.0 - 8.0 * (1.0 - real(nmdim)))

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor, vrt_real, vrt_imag)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = rsphere / sqrt(n_real * (n_real - 1.0))
                vrt_real = real(vrtspec(nm, i))
                vrt_imag = aimag(vrtspec(nm, i))
                cr(m, n, i) = n_factor * vrt_real / scale
                ci(m, n, i) = n_factor * vrt_imag / scale
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine onedtotwod_vrt
