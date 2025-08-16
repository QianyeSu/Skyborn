      subroutine multsmoothfact(dataspec,dataspec_smooth,
     *                          smooth,nlat,nmdim,nt)
      integer nlat,ntrunc,nmdim
      complex dataspec(nmdim,nt)
      complex dataspec_smooth(nmdim,nt)
      real smooth(nlat)

      ntrunc = -1.5 + 0.5*sqrt(9.-8.*(1.-float(nmdim)))

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         dataspec_smooth(nm,i) = dataspec(nm,i)*smooth(n)
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      enddo

      return
      end
