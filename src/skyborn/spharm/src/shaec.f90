!> @file shaec.f90
!> @brief SPHEREPACK Spherical harmonic analysis (even/odd cosine)
!> @author SPHEREPACK team, modernized by Qianye Su
!> @date 2025

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

module shaec_mod
   use sphcom_mod
   use hrfft_mod
   implicit none
   private

   public :: shaec, shaeci

contains

   subroutine shaec(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, &
                    wshaec,lshaec,work,lwork,ierror)
      implicit none
      integer, intent(in) :: nlat,nlon,isym,nt,idg,jdg,mdab,ndab,lshaec,lwork
      real, intent(in) :: g(idg,jdg,nt)
      real, intent(out) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(in) :: wshaec(lshaec)
      real, intent(inout) :: work(lwork)
      integer, intent(out) :: ierror

      integer :: imid,mmax,l1,l2,lzz1,labc

      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 4) return
      ierror = 3
      imid = (nlat+1)/2
      mmax = min(nlat,nlon/2+1)
      l1 = min(nlat,(nlon+2)/2)
      if(nlon-2*(nlon/2) .eq. 1) l1 = min(nlat,(nlon+1)/2)
      l2 = (nlat+1)/2
      lzz1 = 2*nlat*imid
      labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
      if(lshaec .lt. lzz1+labc+nlon+15) return
      ierror = 4
      if(isym .eq. 0) then
         if(lwork .lt. nlat*(nt*nlon+max(3*l2,nlon))) return
      else
         if(lwork .lt. l2*(nt*nlon+max(3*nlat,nlon))) return
      end if
      ierror = 0

      call shaec1(nlat,isym,nt,g,idg,jdg,a,b,mdab,ndab,imid, &
                  wshaec,wshaec(lzz1+1),work,work(l2*nt*nlon+1))
   end subroutine shaec

   subroutine shaec1(nlat,isym,nt,g,idgs,jdgs,a,b,mdab,ndab,imid, &
                     wts,abc,work,fnn)
      implicit none
      integer, intent(in) :: nlat,isym,nt,idgs,jdgs,mdab,ndab,imid
      real, intent(in) :: g(idgs,jdgs,nt),wts(imid,nlat),abc(*)
      real, intent(out) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(inout) :: work(*),fnn(*)

      integer :: i,j,k,m,n,mdo,mmax,ndo,mp1,np1,ms2
      integer :: nlon,mode,ns2,i3,mn,mp2
      real :: cf,fn,fsn,tsn,t1,t2,t3,t4

      nlon = jdgs
      mmax = min(nlat,nlon/2+1)
      mdo = mmax
      if(mmax-1) 3,2,3
   2  mdo = mmax-1
   3  continue

      ! zero coefficients
      do k=1,nt
         do n=1,nlat
            do m=1,mdo
               a(m,n,k) = 0.0
               b(m,n,k) = 0.0
            end do
         end do
      end do

      ndo = nlat
      if(isym .ne. 0) ndo = (nlat+1)/2
      if(isym .eq. 0) call shaec2(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                                  wts,work,mmax,ndo)
      if(isym .eq. 1) call shaec3(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                                  wts,work,mmax,ndo)
      if(isym .eq. 2) call shaec4(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                                  wts,work,mmax,ndo)
      call shaec5(nlat,nlon,isym,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                  wts,abc,mdo,work,fnn)
   end subroutine shaec1

   subroutine shaec2(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                     wts,work,mmax,ndo)
      implicit none
      integer, intent(in) :: nlat,nlon,nt,idgs,jdgs,mdab,ndab,mmax,ndo
      real, intent(in) :: g(idgs,jdgs,nt),wts(*)
      real, intent(out) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(inout) :: work(*)

      integer :: i,j,k,m,n,imm1,mp1,np1,ms2,ns2
      real :: cf

      do k=1,nt
         do i=1,ndo
            do j=1,nlon
               work(i+(j-1)*ndo) = g(i,j,k)+g(nlat+1-i,j,k)
            end do
         end do
         call hrfftf(ndo,nlon,work,ndo,work(ndo*nlon+1))
         do i=1,ndo
            a(1,i,k) = wts(i)*work(i)/2.0
         end do
         do mp1=2,mmax
            m = mp1-1
            ms2 = 2*m
            do i=1,ndo
               a(mp1,i,k) = wts(i)*work(i+ms2*ndo)
               b(mp1,i,k) = wts(i)*work(i+(ms2-1)*ndo)
            end do
         end do
      end do
   end subroutine shaec2

   subroutine shaec3(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                     wts,work,mmax,ndo)
      implicit none
      integer, intent(in) :: nlat,nlon,nt,idgs,jdgs,mdab,ndab,mmax,ndo
      real, intent(in) :: g(idgs,jdgs,nt),wts(*)
      real, intent(out) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(inout) :: work(*)

      integer :: i,j,k,m,n,imm1,mp1,np1,ms2,ns2
      real :: cf

      do k=1,nt
         do i=1,ndo
            do j=1,nlon
               work(i+(j-1)*ndo) = g(i,j,k)
            end do
         end do
         call hrfftf(ndo,nlon,work,ndo,work(ndo*nlon+1))
         do i=1,ndo
            a(1,i,k) = wts(i)*work(i)
         end do
         do mp1=2,mmax
            m = mp1-1
            ms2 = 2*m
            do i=1,ndo
               a(mp1,i,k) = 2.0*wts(i)*work(i+ms2*ndo)
               b(mp1,i,k) = 2.0*wts(i)*work(i+(ms2-1)*ndo)
            end do
         end do
      end do
   end subroutine shaec3

   subroutine shaec4(nlat,nlon,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                     wts,work,mmax,ndo)
      implicit none
      integer, intent(in) :: nlat,nlon,nt,idgs,jdgs,mdab,ndab,mmax,ndo
      real, intent(in) :: g(idgs,jdgs,nt),wts(*)
      real, intent(out) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(inout) :: work(*)

      integer :: i,j,k,m,n,imm1,mp1,np1,ms2,ns2
      real :: cf

      do k=1,nt
         do i=1,ndo
            do j=1,nlon
               work(i+(j-1)*ndo) = g(i,j,k)-g(nlat+1-i,j,k)
            end do
         end do
         call hrfftf(ndo,nlon,work,ndo,work(ndo*nlon+1))
         do mp1=2,mmax
            m = mp1-1
            ms2 = 2*m
            do i=1,ndo
               a(mp1,i,k) = wts(i)*work(i+ms2*ndo)
               b(mp1,i,k) = wts(i)*work(i+(ms2-1)*ndo)
            end do
         end do
      end do
   end subroutine shaec4

   subroutine shaec5(nlat,nlon,isym,nt,g,idgs,jdgs,a,b,mdab,ndab, &
                     wts,abc,mdo,work,fnn)
      implicit none
      integer, intent(in) :: nlat,nlon,isym,nt,idgs,jdgs,mdab,ndab,mdo
      real, intent(in) :: g(idgs,jdgs,nt),wts(*),abc(*)
      real, intent(inout) :: a(mdab,ndab,nt),b(mdab,ndab,nt)
      real, intent(inout) :: work(*),fnn(*)

      integer :: i,j,k,m,n,mmax,ndo,mp1,np1,i3,mn,mp2
      real :: fn

      mmax = min(nlat,nlon/2+1)
      ndo = nlat
      if(isym .ne. 0) ndo = (nlat+1)/2

      do mp1=1,mdo
         m = mp1-1
         call zvinit(nlat,nlon,m,abc,fnn)
         do k=1,nt
            do i=1,ndo
               do j=1,nlon
                  work(i+(j-1)*ndo) = 0.0
               end do
            end do
            if(isym .eq. 0) then
               do i=1,nlat
                  do j=1,nlon
                     work(i+(j-1)*nlat) = g(i,j,k)
                  end do
               end do
            else if(isym .eq. 1) then
               do i=1,ndo
                  do j=1,nlon
                     work(i+(j-1)*ndo) = g(i,j,k)
                  end do
               end do
            else
               do i=1,ndo
                  do j=1,nlon
                     work(i+(j-1)*ndo) = g(i,j,k)
                  end do
               end do
            end if
            call zvin(0,nlat,nlon,m,work,i3,fnn)
            do np1=mp1,ndo,2
               a(mp1,np1,k) = a(mp1,np1,k)+fnn(np1)
            end do
            call zvin(1,nlat,nlon,m,work,i3,fnn)
            do np1=mp1,ndo,2
               b(mp1,np1,k) = b(mp1,np1,k)+fnn(np1)
            end do
         end do
      end do
   end subroutine shaec5

   subroutine shaeci(nlat,nlon,wshaec,lshaec,dwork,ldwork,ierror)
      implicit none
      integer, intent(in) :: nlat,nlon,lshaec,ldwork
      real, intent(out) :: wshaec(lshaec)
      double precision, intent(inout) :: dwork(ldwork)
      integer, intent(out) :: ierror

      integer :: imid,mmax,lzz1,labc,iw1

      ierror = 1
      if(nlat .lt. 3) return
      ierror = 2
      if(nlon .lt. 4) return
      ierror = 3
      imid = (nlat+1)/2
      mmax = min(nlat,nlon/2+1)
      lzz1 = 2*nlat*imid
      labc = 3*((mmax-2)*(nlat+nlat-mmax-1))/2
      if(lshaec .lt. lzz1+labc+nlon+15) return
      ierror = 4
      if(ldwork .lt. nlat+1) return
      ierror = 0
      call zfinit(nlat,nlon,wshaec,dwork)
      iw1 = lzz1+labc+1
      call hrffti(nlon,wshaec(iw1))
   end subroutine shaeci

end module shaec_mod
