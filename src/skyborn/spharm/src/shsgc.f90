!> @file shsgc.f90
!> @brief SPHEREPACK Spherical harmonic synthesis (Gaussian cosine) - OPTIMIZED for modern Fortran
!> @author SPHEREPACK team, optimized by Qianye Su
!> @date 2025
!
!> OPTIMIZATION NOTES:
!> - Modernized from FORTRAN 77 to Fortran 2008+ standards
!> - Added explicit variable declarations with intent specifications
!> - Replaced all GOTO statements with structured control flow
!> - Optimized memory access patterns for better cache efficiency
!> - Added vectorization hints for compiler optimization
!> - Pure serial implementation with SIMD directives for auto-vectorization
!> - Precomputed constants and eliminated redundant calculations
!> - Maintained 100% mathematical accuracy with original algorithms
!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 1998 by UCAR                 .
!  .                                                             .
!  .       University Corporation for Atmospheric Research       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                         SPHEREPACK                          .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

module shsgc_mod
   use sphcom_mod
   use hrfft_mod
   implicit none
   private

   ! F2PY compatible parameters
   integer, parameter :: wp = kind(1.0d0)  ! double precision
   integer, parameter :: ip = kind(1)      ! default integer
   real(wp), parameter :: ZERO = 0.0_wp
   real(wp), parameter :: HALF = 0.5_wp
   real(wp), parameter :: ONE = 1.0_wp
   real(wp), parameter :: TWO = 2.0_wp

   ! Public interfaces
   public :: shsgc, shsgci

contains

   !> @brief Spherical harmonic synthesis on Gaussian grid with computed polynomials
   !> High-performance implementation with modern Fortran optimizations
   subroutine shsgc(nlat, nlon, mode, nt, g, idg, jdg, a, b, mdab, ndab, &
                    wshsgc, lshsgc, work, lwork, ierror)
      ! Input parameters
      integer(ip), intent(in) :: nlat, nlon, mode, nt, idg, jdg, mdab, ndab
      integer(ip), intent(in) :: lshsgc, lwork
      real(wp), intent(in) :: a(mdab, ndab, nt), b(mdab, ndab, nt)
      real(wp), intent(in) :: wshsgc(lshsgc)

      ! Output parameters
      real(wp), intent(out) :: g(idg, jdg, nt)
      real(wp), intent(inout) :: work(lwork)
      integer(ip), intent(out) :: ierror

      ! Local variables
      integer(ip) :: l, late, lat, l1, l2, ifft, ipmn

      ! Input validation with early returns
      if (nlat < 3) then
         ierror = 1; return
      end if
      if (nlon < 4) then
         ierror = 2; return
      end if
      if (mode < 0 .or. mode > 2) then
         ierror = 3; return
      end if
      if (nt < 1) then
         ierror = 4; return
      end if

      ! Set limit for m in a(m,n), b(m,n) computation
      l = min((nlon + 2) / 2, nlat)

      ! Set gaussian point nearest equator pointer
      late = (nlat + mod(nlat, 2)) / 2

      ! Set number of grid points for analysis/synthesis
      lat = merge(late, nlat, mode /= 0)

      ! Validate array dimensions
      if (idg < lat) then
         ierror = 5; return
      end if
      if (jdg < nlon) then
         ierror = 6; return
      end if
      if (mdab < l) then
         ierror = 7; return
      end if
      if (ndab < nlat) then
         ierror = 8; return
      end if

      ! Check workspace sizes
      l1 = l
      l2 = late

      if (lshsgc < nlat * (2 * l2 + 3 * l1 - 2) + 3 * l1 * (1 - l1) / 2 + nlon + 15) then
         ierror = 9; return
      end if

      if (mode == 0) then
         if (lwork < nlat * (nlon * nt + max(3 * l2, nlon))) then
            ierror = 10; return
         end if
      else
         if (lwork < l2 * (nlon * nt + max(3 * nlat, nlon))) then
            ierror = 10; return
         end if
      end if

      ierror = 0

      ! Calculate workspace pointers
      ifft = nlat + 2 * nlat * late + 3 * (l * (l - 1) / 2 + (nlat - l) * (l - 1)) + 1
      ipmn = lat * nlon * nt + 1

      ! Call optimized synthesis routine
      call shsgc1_optimized(nlat, nlon, l, lat, mode, g, idg, jdg, nt, a, b, mdab, ndab, &
                           wshsgc, wshsgc(ifft:), late, work(ipmn:), work)

   end subroutine shsgc

   !> @brief Optimized spherical harmonic synthesis implementation
   subroutine shsgc1_optimized(nlat, nlon, l, lat, mode, gs, idg, jdg, nt, a, b, mdab, &
                              ndab, w, wfft, late, pmn, g)
      ! Input parameters
      integer(ip), intent(in) :: nlat, nlon, l, lat, mode, idg, jdg, nt, mdab, ndab, late
      real(wp), intent(in) :: a(mdab, ndab, nt), b(mdab, ndab, nt)
      real(wp), intent(in) :: w(*), wfft(*)

      ! Working arrays
      real(wp), intent(inout) :: pmn(nlat, late, 3), g(lat, nlon, nt)

      ! Output
      real(wp), intent(out) :: gs(idg, jdg, nt)

      ! Local variables
      integer(ip) :: lm1, m, mp1, mp2, k, np1, i, is, nl2, km
      integer(ip) :: meo, ms, ns, lp1, j
      real(wp) :: t1, t2, t3, t4

      ! Initialize to zero with vectorization
      call initialize_grid_shsgc(g, lat, nlon, nt)

      lm1 = merge(l - 1, l, nlon == l + l - 2)

      if (mode == 0) then
         ! Full sphere case
         call process_full_sphere_shsgc(nlat, nlon, l, lat, late, lm1, nt, &
                                        a, b, mdab, ndab, w, pmn, g)
      else
         ! Half sphere case
         call process_half_sphere_shsgc(nlat, nlon, l, lat, late, lm1, nt, mode, &
                                        a, b, mdab, ndab, w, pmn, g)
      end if

      ! Apply inverse FFT with vectorization
      call apply_inverse_fft_shsgc(g, lat, nlon, nt, wfft, pmn)

      ! Scale and copy output with vectorization
      call scale_output_shsgc(gs, g, idg, jdg, lat, nlon, nt)

   end subroutine shsgc1_optimized

   !> @brief Initialize grid to zero with vectorization
   subroutine initialize_grid_shsgc(g, lat, nlon, nt)
      integer(ip), intent(in) :: lat, nlon, nt
      real(wp), intent(out) :: g(lat, nlon, nt)

      integer(ip) :: k, j, i

      do k = 1, nt
         do j = 1, nlon
            !DIR$ SIMD
            do i = 1, lat
               g(i, j, k) = ZERO
            end do
         end do
      end do

   end subroutine initialize_grid_shsgc

   !> @brief Process full sphere synthesis for shsgc
   subroutine process_full_sphere_shsgc(nlat, nlon, l, lat, late, lm1, nt, &
                                       a, b, mdab, ndab, w, pmn, g)
      integer(ip), intent(in) :: nlat, nlon, l, lat, late, lm1, nt, mdab, ndab
      real(wp), intent(in) :: a(mdab, ndab, nt), b(mdab, ndab, nt)
      real(wp), intent(in) :: w(*)
      real(wp), intent(inout) :: pmn(nlat, late, 3), g(lat, nlon, nt)

      integer(ip) :: m, mp1, mp2, k, np1, i, is, nl2, km, lp1
      real(wp) :: t1, t2, t3, t4

      ! Set first column in g (m = 0)
      m = 0
      call legin(0, l, nlat, m, w, pmn, km)

      do k = 1, nt
         ! n even
         !DIR$ SIMD
         do np1 = 1, nlat, 2
            !DIR$ SIMD
            do i = 1, late
               g(i, 1, k) = g(i, 1, k) + a(1, np1, k) * pmn(np1, i, km)
            end do
         end do

         ! n odd
         nl2 = nlat / 2
         !DIR$ SIMD
         do np1 = 2, nlat, 2
            !DIR$ SIMD
            do i = 1, nl2
               is = nlat - i + 1
               g(is, 1, k) = g(is, 1, k) + a(1, np1, k) * pmn(np1, i, km)
            end do
         end do

         ! Restore m=0 coefficients (reverse implicit even/odd reduction)
         !DIR$ SIMD
         do i = 1, nl2
            is = nlat - i + 1
            t1 = g(i, 1, k)
            t3 = g(is, 1, k)
            g(i, 1, k) = t1 + t3
            g(is, 1, k) = t1 - t3
         end do
      end do

      ! Sweep columns of g for which b is available
      do mp1 = 2, lm1
         m = mp1 - 1
         mp2 = m + 2
         call legin(0, l, nlat, m, w, pmn, km)

         do k = 1, nt
            ! For n-m even store synthesis
            !DIR$ SIMD
            do np1 = mp1, nlat, 2
               !DIR$ SIMD
               do i = 1, late
                  g(i, 2*m, k) = g(i, 2*m, k) + a(mp1, np1, k) * pmn(np1, i, km)
                  g(i, 2*m+1, k) = g(i, 2*m+1, k) + b(mp1, np1, k) * pmn(np1, i, km)
               end do
            end do

            ! For n-m odd store synthesis
            !DIR$ SIMD
            do np1 = mp2, nlat, 2
               !DIR$ SIMD
               do i = 1, nl2
                  is = nlat - i + 1
                  g(is, 2*m, k) = g(is, 2*m, k) + a(mp1, np1, k) * pmn(np1, i, km)
                  g(is, 2*m+1, k) = g(is, 2*m+1, k) + b(mp1, np1, k) * pmn(np1, i, km)
               end do
            end do

            ! Set fourier coefficients using even-odd reduction
            !DIR$ SIMD
            do i = 1, nl2
               is = nlat - i + 1
               t1 = g(i, 2*m, k)
               t2 = g(i, 2*m+1, k)
               t3 = g(is, 2*m, k)
               t4 = g(is, 2*m+1, k)
               g(i, 2*m, k) = t1 + t3
               g(i, 2*m+1, k) = t2 + t4
               g(is, 2*m, k) = t1 - t3
               g(is, 2*m+1, k) = t2 - t4
            end do
         end do
      end do

      ! Set last column (using a only) if necessary
      if (nlon == l + l - 2) then
         m = l - 1
         call legin(0, l, nlat, m, w, pmn, km)

         do k = 1, nt
            ! n-m even
            !DIR$ SIMD
            do np1 = l, nlat, 2
               !DIR$ SIMD
               do i = 1, late
                  g(i, nlon, k) = g(i, nlon, k) + TWO * a(l, np1, k) * pmn(np1, i, km)
               end do
            end do

            lp1 = l + 1
            ! n-m odd
            !DIR$ SIMD
            do np1 = lp1, nlat, 2
               !DIR$ SIMD
               do i = 1, nl2
                  is = nlat - i + 1
                  g(is, nlon, k) = g(is, nlon, k) + TWO * a(l, np1, k) * pmn(np1, i, km)
               end do
            end do

            !DIR$ SIMD
            do i = 1, nl2
               is = nlat - i + 1
               t1 = g(i, nlon, k)
               t3 = g(is, nlon, k)
               g(i, nlon, k) = t1 + t3
               g(is, nlon, k) = t1 - t3
            end do
         end do
      end if

   end subroutine process_full_sphere_shsgc

   !> @brief Process half sphere synthesis for shsgc
   subroutine process_half_sphere_shsgc(nlat, nlon, l, lat, late, lm1, nt, mode, &
                                       a, b, mdab, ndab, w, pmn, g)
      integer(ip), intent(in) :: nlat, nlon, l, lat, late, lm1, nt, mode, mdab, ndab
      real(wp), intent(in) :: a(mdab, ndab, nt), b(mdab, ndab, nt)
      real(wp), intent(in) :: w(*)
      real(wp), intent(inout) :: pmn(nlat, late, 3), g(lat, nlon, nt)

      integer(ip) :: m, mp1, k, np1, i, meo, ms, ns, km

      ! Set first column in g
      m = 0
      meo = merge(2, 1, mode == 1)
      ms = m + meo
      call legin(mode, l, nlat, m, w, pmn, km)

      do k = 1, nt
         !DIR$ SIMD
         do np1 = ms, nlat, 2
            !DIR$ SIMD
            do i = 1, late
               g(i, 1, k) = g(i, 1, k) + a(1, np1, k) * pmn(np1, i, km)
            end do
         end do
      end do

      ! Sweep interior columns of g
      do mp1 = 2, lm1
         m = mp1 - 1
         ms = m + meo
         call legin(mode, l, nlat, m, w, pmn, km)

         do k = 1, nt
            !DIR$ SIMD
            do np1 = ms, nlat, 2
               !DIR$ SIMD
               do i = 1, late
                  g(i, 2*m, k) = g(i, 2*m, k) + a(mp1, np1, k) * pmn(np1, i, km)
                  g(i, 2*m+1, k) = g(i, 2*m+1, k) + b(mp1, np1, k) * pmn(np1, i, km)
               end do
            end do
         end do
      end do

      if (nlon == l + l - 2) then
         ! Set last column
         m = l - 1
         call legin(mode, l, nlat, m, w, pmn, km)
         ns = merge(l + 1, l, mode == 1)

         do k = 1, nt
            !DIR$ SIMD
            do np1 = ns, nlat, 2
               !DIR$ SIMD
               do i = 1, late
                  g(i, nlon, k) = g(i, nlon, k) + TWO * a(l, np1, k) * pmn(np1, i, km)
               end do
            end do
         end do
      end if

   end subroutine process_half_sphere_shsgc

   !> @brief Apply inverse FFT with vectorization for shsgc
   subroutine apply_inverse_fft_shsgc(g, lat, nlon, nt, wfft, pmn)
      integer(ip), intent(in) :: lat, nlon, nt
      real(wp), intent(inout) :: g(lat, nlon, nt)
      real(wp), intent(in) :: wfft(*)
      real(wp), intent(inout) :: pmn(*)

      integer(ip) :: k

      do k = 1, nt
         call hrfftb(lat, nlon, g(1, 1, k), lat, wfft, pmn)
      end do

   end subroutine apply_inverse_fft_shsgc

   !> @brief Scale output with vectorization for shsgc
   subroutine scale_output_shsgc(gs, g, idg, jdg, lat, nlon, nt)
      integer(ip), intent(in) :: idg, jdg, lat, nlon, nt
      real(wp), intent(in) :: g(lat, nlon, nt)
      real(wp), intent(out) :: gs(idg, jdg, nt)

      integer(ip) :: k, j, i

      do k = 1, nt
         do j = 1, nlon
            !DIR$ SIMD
            do i = 1, lat
               gs(i, j, k) = HALF * g(i, j, k)
            end do
         end do
      end do

   end subroutine scale_output_shsgc

   !> @brief Initialize workspace for spherical harmonic synthesis
   subroutine shsgci(nlat, nlon, wshsgc, lshsgc, dwork, ldwork, ierror)
      ! Input parameters
      integer(ip), intent(in) :: nlat, nlon, lshsgc, ldwork
      real(wp), intent(inout) :: dwork(ldwork)

      ! Output parameters
      real(wp), intent(out) :: wshsgc(lshsgc)
      integer(ip), intent(out) :: ierror

      ! Local variables
      integer(ip) :: l, late, l1, l2
      integer(ip) :: i1, i2, i3, i4, i5, i6, i7, idth, idwts, iw

      ! Input validation
      if (nlat < 3) then
         ierror = 1; return
      end if
      if (nlon < 4) then
         ierror = 2; return
      end if

      ! Set triangular truncation limit
      l = min((nlon + 2) / 2, nlat)
      late = (nlat + mod(nlat, 2)) / 2
      l1 = l
      l2 = late

      ! Check workspace sizes
      if (lshsgc < nlat * (2 * l2 + 3 * l1 - 2) + 3 * l1 * (1 - l1) / 2 + nlon + 15) then
         ierror = 3; return
      end if

      if (ldwork < nlat * (nlat + 4)) then
         ierror = 4; return
      end if

      ierror = 0

      ! Set pointers
      i1 = 1
      i2 = i1 + nlat
      i3 = i2 + nlat * late
      i4 = i3 + nlat * late
      i5 = i4 + l * (l - 1) / 2 + (nlat - l) * (l - 1)
      i6 = i5 + l * (l - 1) / 2 + (nlat - l) * (l - 1)
      i7 = i6 + l * (l - 1) / 2 + (nlat - l) * (l - 1)

      ! Set indices for double precision Gaussian weights and points
      idth = 1
      idwts = idth + nlat
      iw = idwts + nlat

      call shsgci1_optimized(nlat, nlon, l, late, wshsgc(i1:), wshsgc(i2:), wshsgc(i3:), &
                            wshsgc(i4:), wshsgc(i5:), wshsgc(i6:), wshsgc(i7:), &
                            dwork(idth:), dwork(idwts:), dwork(iw:), ierror)

      if (ierror /= 0) ierror = 5

   end subroutine shsgci

   !> @brief Optimized initialization subroutine
   subroutine shsgci1_optimized(nlat, nlon, l, late, wts, p0n, p1n, abel, bbel, cbel, &
                                wfft, dtheta, dwts, work, ier)
      integer(ip), intent(in) :: nlat, nlon, l, late
      real(wp), intent(out) :: wts(nlat), p0n(nlat, late), p1n(nlat, late)
      real(wp), intent(out) :: abel(*), bbel(*), cbel(*), wfft(*)
      real(wp), intent(inout) :: dtheta(nlat), dwts(nlat)
      real(wp), intent(inout) :: work(*)
      integer(ip), intent(out) :: ier

      real(wp) :: pb
      integer(ip) :: np1, n, m, i, lw, imn, mlim

      ! Index functions
      integer(ip) :: indx, imndx
      indx(m, n) = (n - 1) * (n - 2) / 2 + m - 1
      imndx(m, n) = l * (l - 1) / 2 + (n - l - 1) * (l - 1) + m - 1

      call hrffti(nlon, wfft)

      ! Compute double precision Gaussian points and weights
      lw = nlat * (nlat + 2)
      call gaqd(nlat, dtheta, dwts, work, lw, ier)
      if (ier /= 0) return

      ! Store Gaussian weights in single precision
      !DIR$ SIMD
      do i = 1, nlat
         wts(i) = real(dwts(i), wp)
      end do

      ! Initialize p0n, p1n
      do np1 = 1, nlat
         !DIR$ SIMD
         do i = 1, late
            p0n(np1, i) = ZERO
            p1n(np1, i) = ZERO
         end do
      end do

      ! Compute m=n=0 Legendre polynomials
      np1 = 1
      n = 0
      m = 0
      call dnlfk(m, n, work)
      do i = 1, late
         call dnlft(m, n, dtheta(i), work, pb)
         p0n(1, i) = real(pb, wp)
      end do

      ! Compute p0n, p1n for all theta(i) when n > 0
      do np1 = 2, nlat
         n = np1 - 1
         m = 0
         call dnlfk(m, n, work)
         do i = 1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p0n(np1, i) = real(pb, wp)
         end do

         ! Compute m=1 Legendre polynomials
         m = 1
         call dnlfk(m, n, work)
         do i = 1, late
            call dnlft(m, n, dtheta(i), work, pb)
            p1n(np1, i) = real(pb, wp)
         end do
      end do

      ! Compute Swarztrauber recursion coefficients
      do n = 2, nlat
         mlim = min(n, l)
         do m = 2, mlim
            imn = indx(m, n)
            if (n >= l) imn = imndx(m, n)

            abel(imn) = sqrt(real((2*n+1)*(m+n-2)*(m+n-3), wp) / &
                            real((2*n-3)*(m+n-1)*(m+n), wp))
            bbel(imn) = sqrt(real((2*n+1)*(n-m-1)*(n-m), wp) / &
                            real((2*n-3)*(m+n-1)*(m+n), wp))
            cbel(imn) = sqrt(real((n-m+1)*(n-m+2), wp) / &
                            real((n+m-1)*(n+m), wp))
         end do
      end do

   end subroutine shsgci1_optimized

end module shsgc_mod
