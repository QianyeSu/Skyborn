      subroutine onedtotwod(dataspec,a,b,nlat,nmdim,nt)
      integer nlat,nmdim,nt,nmstrt,i
      real a(nlat,nlat,nt),b(nlat,nlat,nt)
      complex dataspec(nmdim,nt)

      ntrunc = -1.5 + 0.5*sqrt(9.-8.*(1.-float(nmdim)))

      scale = 0.5

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         a(m,n,i) = real(dataspec(nm,i)/scale)
         b(m,n,i) = aimag(dataspec(nm,i)/scale)
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      enddo

      return
      end
