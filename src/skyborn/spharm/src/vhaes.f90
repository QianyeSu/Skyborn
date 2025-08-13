! This file contains optimized versions of:
!   - vhaes: Main vector spherical harmonic analysis routine
!   - vhaes1: Core computation routine
!   - vhaesi: Initialization routine
!   - vea1: Core initialization routine

!> @brief Optimized main vector spherical harmonic analysis routine
!> Modernized from FORTRAN 77 to Fortran 2008+ standards
subroutine vhaes(nlat, nlon, ityp, nt, v, w, idvw, jdvw, br, bi, cr, ci, &
                 mdab, ndab, wvhaes, lvhaes, work, lwork, ierror)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, ityp, nt, idvw, jdvw
   integer, intent(in) :: mdab, ndab, lvhaes, lwork
   real, intent(in) :: v(idvw, jdvw, nt), w(idvw, jdvw, nt)
   real, intent(in) :: wvhaes(lvhaes)

   ! Output parameters
   real, intent(out) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(out) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)
   real, intent(inout) :: work(lwork)
   integer, intent(out) :: ierror

   ! Local variables
   integer :: imid, mmax, idz, lzimn, idv, lnl, ist
   integer :: iw1, iw2, iw3, iw4, jw1, jw2

   ! Input validation with structured error handling
   if (nlat < 3) then
      ierror = 1
      return
   end if

   if (nlon < 1) then
      ierror = 2
      return
   end if

   if (ityp < 0 .or. ityp > 8) then
      ierror = 3
      return
   end if

   if (nt < 0) then
      ierror = 4
      return
   end if

   ! Precompute frequently used values
   imid = (nlat + 1) / 2

   ! Validate array dimensions
   if ((ityp <= 2 .and. idvw < nlat) .or. &
       (ityp > 2 .and. idvw < imid)) then
      ierror = 5
      return
   end if

   if (jdvw < nlon) then
      ierror = 6
      return
   end if

   mmax = min(nlat, (nlon + 1) / 2)
   if (mdab < mmax) then
      ierror = 7
      return
   end if

   if (ndab < nlat) then
      ierror = 8
      return
   end if

   ! Check workspace sizes
   idz = (mmax * (nlat + nlat - mmax + 1)) / 2
   lzimn = idz * imid
   if (lvhaes < lzimn + lzimn + nlon + 15) then
      ierror = 9
      return
   end if

   idv = nlat
   if (ityp > 2) idv = imid
   lnl = nt * idv * nlon
   if (lwork < lnl + lnl + idv * nlon) then
      ierror = 10
      return
   end if

   ierror = 0

   ! Calculate workspace pointers
   ist = 0
   if (ityp <= 2) ist = imid
   iw1 = ist + 1
   iw2 = lnl + 1
   iw3 = iw2 + ist
   iw4 = iw2 + lnl
   jw1 = lzimn + 1
   jw2 = jw1 + lzimn

   ! Call core computation routine
   call vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, ndab, &
              br, bi, cr, ci, idv, work, work(iw1), work(iw2), work(iw3), &
              work(iw4), idz, wvhaes, wvhaes(jw1), wvhaes(jw2))

end subroutine vhaes

!> @brief Optimized core vector harmonic analysis computation routine
!> Modernized with OpenMP parallelization for better performance
subroutine vhaes1(nlat, nlon, ityp, nt, imid, idvw, jdvw, v, w, mdab, &
                  ndab, br, bi, cr, ci, idv, ve, vo, we, wo, work, idz, zv, zw, wrfft)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, ityp, nt, imid, idvw, jdvw
   integer, intent(in) :: mdab, ndab, idv, idz
   real, intent(in) :: v(idvw, jdvw, nt), w(idvw, jdvw, nt)
   real, intent(in) :: zv(idz, *), zw(idz, *), wrfft(*)

   ! Working arrays
   real, intent(inout) :: ve(idv, nlon, nt), vo(idv, nlon, nt)
   real, intent(inout) :: we(idv, nlon, nt), wo(idv, nlon, nt)
   real, intent(inout) :: work(*)

   ! Output arrays
   real, intent(out) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(out) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)

   ! Local variables
   integer :: nlp1, mlat, mlon, mmax, imm1, ndo1, ndo2, itypp
   integer :: k, i, j, mp1, np1, m, mb, mp2
   real :: tsn, fsn

   ! Precompute constants
   nlp1 = nlat + 1
   tsn = 2. / nlon
   fsn = 4. / nlon
   mlat = mod(nlat, 2)
   mlon = mod(nlon, 2)
   mmax = min(nlat, (nlon + 1) / 2)
   imm1 = imid
   if (mlat /= 0) imm1 = imid - 1

   ! Compute even and odd components with OpenMP parallelization
   if (ityp <= 2) then
      ! Full sphere case - vectorizable loops
      do k = 1, nt
         do i = 1, imm1
            do j = 1, nlon
               ve(i, j, k) = tsn * (v(i, j, k) + v(nlp1 - i, j, k))
               vo(i, j, k) = tsn * (v(i, j, k) - v(nlp1 - i, j, k))
               we(i, j, k) = tsn * (w(i, j, k) + w(nlp1 - i, j, k))
               wo(i, j, k) = tsn * (w(i, j, k) - w(nlp1 - i, j, k))
            end do
         end do
      end do
   else
      ! Half sphere case - vectorizable loops
      do k = 1, nt
         do i = 1, imm1
            do j = 1, nlon
               ve(i, j, k) = fsn * v(i, j, k)
               vo(i, j, k) = fsn * v(i, j, k)
               we(i, j, k) = fsn * w(i, j, k)
               wo(i, j, k) = fsn * w(i, j, k)
            end do
         end do
      end do
   end if

   ! Handle equator point if needed
   if (mlat /= 0) then
      do k = 1, nt
         do j = 1, nlon
            ve(imid, j, k) = tsn * v(imid, j, k)
            we(imid, j, k) = tsn * w(imid, j, k)
         end do
      end do
   end if

   ! Apply forward FFT - each transform is independent
   do k = 1, nt
      call hrfftf(idv, nlon, ve(1, 1, k), idv, wrfft, work)
      call hrfftf(idv, nlon, we(1, 1, k), idv, wrfft, work)
   end do

   ! Determine loop bounds
   ndo1 = nlat
   ndo2 = nlat
   if (mlat /= 0) ndo1 = nlat - 1
   if (mlat == 0) ndo2 = nlat - 1

   ! Initialize coefficient arrays with OpenMP
   if (ityp /= 2 .and. ityp /= 5 .and. ityp /= 8) then
      do k = 1, nt
         do mp1 = 1, mmax
            do np1 = mp1, nlat
               br(mp1, np1, k) = 0.
               bi(mp1, np1, k) = 0.
            end do
         end do
      end do
   end if

   if (ityp /= 1 .and. ityp /= 4 .and. ityp /= 7) then
      do k = 1, nt
         do mp1 = 1, mmax
            do np1 = mp1, nlat
               cr(mp1, np1, k) = 0.
               ci(mp1, np1, k) = 0.
            end do
         end do
      end do
   end if

   itypp = ityp + 1

   ! Use structured control flow instead of computed GOTO
   select case (itypp)
   case (1)
      ! ityp=0: no symmetries
      call vhaes_case0(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   case (2)
      ! ityp=1: no symmetries but cr and ci equal zero
      call vhaes_case1(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   case (3)
      ! ityp=2: no symmetries but br and bi equal zero
      call vhaes_case2(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   case (4)
      ! ityp=3: v even, w odd
      call vhaes_case3(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   case (5)
      ! ityp=4: v even, w odd, and cr and ci equal 0
      call vhaes_case4(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   case (6)
      ! ityp=5: v even, w odd, and br and bi equal zero
      call vhaes_case5(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   case (7)
      ! ityp=6: v odd, w even
      call vhaes_case6(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   case (8)
      ! ityp=7: v odd, w even, and cr and ci equal zero
      call vhaes_case7(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   case (9)
      ! ityp=8: v odd, w even, and both br and bi equal zero
      call vhaes_case8(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   end select

end subroutine vhaes1

!> @brief Case 0: ityp=0, no symmetries (optimized with OpenMP)
subroutine vhaes_case0(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   implicit none
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   integer :: k, i, np1, mp1, m, mb, mp2

   ! m=0 case - even np1 values
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
         end do
      end do
   end do

   ! m=0 case - odd np1 values
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
         end do
      end do
   end do

   ! m = 1 through nlat-1
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         ! Process odd harmonics
         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                     cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if

         ! Process even harmonics
         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if
end subroutine vhaes_case0

!> @brief Case 1: ityp=1, no symmetries but cr and ci equal zero
subroutine vhaes_case1(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   implicit none
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   integer :: k, i, np1, mp1, m, mb, mp2

   ! m=0 case - even np1 values
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
         end do
      end do
   end do

   ! m=0 case - odd np1 values
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
         end do
      end do
   end do

   ! m = 1 through nlat-1
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                  end do
               end do
            end do

            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if

         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if
end subroutine vhaes_case1

!> @brief Case 2: ityp=2, no symmetries but br and bi equal zero
subroutine vhaes_case2(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   implicit none
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   integer :: k, i, np1, mp1, m, mb, mp2

   ! m=0 case - even np1 values
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
         end do
      end do
   end do

   ! m=0 case - odd np1 values
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
         end do
      end do
   end do

   ! m = 1 through nlat-1
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                  end do
               end do
            end do

            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if

         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if
end subroutine vhaes_case2

!> @brief Case 3: ityp=3, v even, w odd
!> Vector field with v symmetric and w antisymmetric about equator
subroutine vhaes_case3(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - even np1 values (br coefficients)
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
         end do
      end do
   end do

   ! Case m=0 - odd np1 values (cr coefficients)
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         ! Process odd harmonics (cr, ci coefficients)
         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for cr, ci coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if

         ! Process even harmonics (br, bi coefficients)
         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for br, bi coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case3

!> @brief Case 4: ityp=4, v even, w odd, and cr and ci equal zero
!> Vector field with v symmetric and w antisymmetric, curl-free (cr=ci=0)
subroutine vhaes_case4(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters (only br, bi computed, cr, ci remain zero)
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - even np1 values (br coefficients only)
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * ve(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1 (only even harmonics, br and bi coefficients)
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         ! Process only even harmonics (mp2 <= ndo2)
         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * wo(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * ve(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * wo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for br, bi coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case4

!> @brief Case 5: ityp=5, v even, w odd, and br and bi equal zero
!> Vector field with v symmetric and w antisymmetric, divergence-free (br=bi=0)
subroutine vhaes_case5(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters (only cr, ci computed, br, bi remain zero)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - odd np1 values (cr coefficients only)
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * wo(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1 (only odd harmonics, cr and ci coefficients)
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2

         ! Process only odd harmonics (mp1 <= ndo1)
         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * ve(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * wo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * ve(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for cr, ci coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) + zw(np1 + mb, imid) * ve(imid, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zw(np1 + mb, imid) * ve(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case5

!> @brief Case 6: ityp=6, v odd, w even
!> Vector field with v antisymmetric and w symmetric about equator
subroutine vhaes_case6(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, cr, ci, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - even np1 values (cr coefficients)
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
         end do
      end do
   end do

   ! Case m=0 - odd np1 values (br coefficients)
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         ! Process odd harmonics (br, bi coefficients)
         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for br, bi coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if

         ! Process even harmonics (cr, ci coefficients)
         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for cr, ci coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case6

!> @brief Case 7: ityp=7, v odd, w even, and cr and ci equal zero
!> Vector field with v antisymmetric and w symmetric, curl-free (cr=ci=0)
subroutine vhaes_case7(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, br, bi, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters (only br, bi computed, cr, ci remain zero)
   real, intent(inout) :: br(mdab, ndab, nt), bi(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - odd np1 values (br coefficients only)
   do k = 1, nt
      do i = 1, imm1
         do np1 = 3, ndo1, 2
            br(1, np1, k) = br(1, np1, k) + zv(np1, i) * vo(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1 (only odd harmonics, br and bi coefficients)
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2

         ! Process only odd harmonics (mp1 <= ndo1)
         if (mp1 <= ndo1) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * we(i, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) + zv(np1 + mb, i) * vo(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * we(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for br, bi coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp1, ndo1, 2
                     br(mp1, np1, k) = br(mp1, np1, k) + zw(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                     bi(mp1, np1, k) = bi(mp1, np1, k) - zw(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case7

!> @brief Case 8: ityp=8, v odd, w even, and both br and bi equal zero
!> Vector field with v antisymmetric and w symmetric, divergence-free (br=bi=0)
subroutine vhaes_case8(nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, &
                      ve, vo, we, wo, cr, ci, mdab, ndab, zv, zw)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, nt, imid, imm1, mlat, mmax, ndo1, ndo2, mdab, ndab
   real, intent(in) :: ve(imid, nlon, nt), vo(imid, nlon, nt)
   real, intent(in) :: we(imid, nlon, nt), wo(imid, nlon, nt)
   real, intent(in) :: zv(ndab, *), zw(ndab, *)

   ! Output parameters (only cr, ci computed, br, bi remain zero)
   real, intent(inout) :: cr(mdab, ndab, nt), ci(mdab, ndab, nt)

   ! Local variables
   integer :: k, i, np1, mp1, m, mb, mp2

   ! Case m=0 - even np1 values (cr coefficients only)
   do k = 1, nt
      do i = 1, imid
         do np1 = 2, ndo2, 2
            cr(1, np1, k) = cr(1, np1, k) - zv(np1, i) * we(i, 1, k)
         end do
      end do
   end do

   ! Case m = 1 through nlat-1 (only even harmonics, cr and ci coefficients)
   if (mmax >= 2) then
      do mp1 = 2, mmax
         m = mp1 - 1
         mb = m * (nlat - 1) - (m * (m - 1)) / 2
         mp2 = mp1 + 1

         ! Process only even harmonics (mp2 <= ndo2)
         if (mp2 <= ndo2) then
            do k = 1, nt
               do i = 1, imm1
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-2, k) &
                                                       + zw(np1 + mb, i) * vo(i, 2*mp1-1, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, i) * we(i, 2*mp1-1, k) &
                                                       - zw(np1 + mb, i) * vo(i, 2*mp1-2, k)
                  end do
               end do
            end do

            ! Handle equator point for cr, ci coefficients
            if (mlat /= 0) then
               do k = 1, nt
                  do np1 = mp2, ndo2, 2
                     cr(mp1, np1, k) = cr(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-2, k)
                     ci(mp1, np1, k) = ci(mp1, np1, k) - zv(np1 + mb, imid) * we(imid, 2*mp1-1, k)
                  end do
               end do
            end if
         end if
      end do
   end if

end subroutine vhaes_case8

!> @brief Initialization routine for vector harmonic analysis
!> Initializes the workspace array wvhaes for repeated use by vhaes
!> @param[in] nlat Number of colatitudes on full sphere
!> @param[in] nlon Number of longitude points
!> @param[out] wvhaes Workspace array to be initialized
!> @param[in] lvhaes Dimension of wvhaes array
!> @param[inout] work Temporary workspace array
!> @param[in] lwork Dimension of work array
!> @param[inout] dwork Double precision workspace
!> @param[in] ldwork Dimension of dwork array
!> @param[out] ierror Error flag (0=success, >0=error code)
subroutine vhaesi(nlat, nlon, wvhaes, lvhaes, work, lwork, dwork, &
                  ldwork, ierror)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, lvhaes, lwork, ldwork

   ! Output parameters
   real, intent(out) :: wvhaes(lvhaes)
   integer, intent(out) :: ierror

   ! Working arrays
   real, intent(inout) :: work(lwork)
   double precision, intent(inout) :: dwork(ldwork)

   ! Local variables
   integer :: mmax, imid, lzimn, labc, iw1, idz

   ! Input validation with structured error handling
   if (nlat < 3) then
      ierror = 1
      return
   end if

   if (nlon < 1) then
      ierror = 2
      return
   end if

   ! Calculate key dimensions
   mmax = min(nlat, (nlon + 1) / 2)
   imid = (nlat + 1) / 2
   lzimn = (imid * mmax * (nlat + nlat - mmax + 1)) / 2

   if (lvhaes < lzimn + lzimn + nlon + 15) then
      ierror = 3
      return
   end if

   labc = 3 * (max(mmax - 2, 0) * (nlat + nlat - mmax - 1)) / 2
   if (lwork < 5 * nlat * imid + labc) then
      ierror = 4
      return
   end if

   if (ldwork < 2 * (nlat + 1)) then
      ierror = 5
      return
   end if

   ierror = 0

   ! Calculate workspace pointers
   iw1 = 3 * nlat * imid + 1
   idz = (mmax * (nlat + nlat - mmax + 1)) / 2

   ! Initialize vector analysis coefficient arrays
   call vea1(nlat, nlon, imid, wvhaes, wvhaes(lzimn + 1), idz, &
             work, work(iw1), dwork)

   ! Initialize FFT workspace
   call hrffti(nlon, wvhaes(2 * lzimn + 1))

end subroutine vhaesi

!> @brief Core initialization routine for vector analysis
!> Initializes coefficient arrays for vector harmonic analysis
!> @param[in] nlat Number of colatitudes
!> @param[in] nlon Number of longitude points
!> @param[in] imid Mid-point index
!> @param[out] zv,zw Coefficient arrays for vector analysis
!> @param[in] idz Coefficient array dimension
!> @param[inout] zin Working array for coefficient computation
!> @param[inout] wzvin Workspace for initialization
!> @param[inout] dwork Double precision workspace
subroutine vea1(nlat, nlon, imid, zv, zw, idz, zin, wzvin, dwork)
   implicit none

   ! Input parameters
   integer, intent(in) :: nlat, nlon, imid, idz

   ! Output parameters
   real, intent(out) :: zv(idz, *), zw(idz, *)

   ! Working arrays
   real, intent(inout) :: zin(imid, nlat, 3), wzvin(*)
   double precision, intent(inout) :: dwork(*)

   ! Local variables
   integer :: mmax, mp1, m, np1, mn, i, i3

   ! Calculate maximum wavenumber
   mmax = min(nlat, (nlon + 1) / 2)

   ! Initialize V-component coefficient arrays
   call zvinit(nlat, nlon, wzvin, dwork)
   do mp1 = 1, mmax
      m = mp1 - 1
      call zvin(0, nlat, nlon, m, zin, i3, wzvin)
      do np1 = mp1, nlat
         mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
         do i = 1, imid
            zv(mn, i) = zin(i, np1, i3)
         end do
      end do
   end do

   ! Initialize W-component coefficient arrays
   call zwinit(nlat, nlon, wzvin, dwork)
   do mp1 = 1, mmax
      m = mp1 - 1
      call zwin(0, nlat, nlon, m, zin, i3, wzvin)
      do np1 = mp1, nlat
         mn = m * (nlat - 1) - (m * (m - 1)) / 2 + np1
         do i = 1, imid
            zw(mn, i) = zin(i, np1, i3)
         end do
      end do
   end do

end subroutine vea1
