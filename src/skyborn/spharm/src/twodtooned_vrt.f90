!
! Converts vector harmonic coefficients to vorticity spectral coefficients only.
!
subroutine twodtooned_vrt(vrtspec, cr, ci, nlat, ntrunc, nt, rsphere)
    implicit none

    integer, intent(in) :: nlat, ntrunc, nt
    real, intent(in) :: rsphere
    real, intent(in) :: cr(nlat, nlat, nt), ci(nlat, nlat, nt)
    complex, intent(out) :: vrtspec((ntrunc+1)*(ntrunc+2)/2, nt)

    integer :: m, n, nm, nmstrt, i
    real :: scale, rsphere_inv, n_real, n_factor

    scale = 0.5
    rsphere_inv = 1.0 / rsphere

    do i = 1, nt
        nmstrt = 0
        do m = 1, ntrunc + 1
            !$OMP SIMD PRIVATE(nm, n_real, n_factor)
            do n = m, ntrunc + 1
                nm = nmstrt + n - m + 1
                n_real = real(n)
                n_factor = sqrt(n_real * (n_real - 1.0)) * rsphere_inv * scale
                vrtspec(nm, i) = n_factor * cmplx(cr(m, n, i), ci(m, n, i))
            end do
            nmstrt = nmstrt + ntrunc - m + 2
        end do
    end do

end subroutine twodtooned_vrt
