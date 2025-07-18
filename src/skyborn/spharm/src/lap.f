      subroutine lap(dataspec,dataspec_lap,nmdim,nt,rsphere)
      integer nmdim
      complex dataspec(nmdim,nt),dataspec_lap(nmdim,nt)
      real rsphere

      ntrunc = -1.5 + 0.5*sqrt(9.-8.*(1.-float(nmdim)))

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         dataspec_lap(nm,i) =
     *   -(float(n)*float(n-1)/rsphere**2)*dataspec(nm,i)
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      enddo

      return
      end
