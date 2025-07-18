      subroutine invlap(dataspec,dataspec_ilap,nmdim,nt,rsphere)
      integer nmdim
      complex dataspec(nmdim,nt),dataspec_ilap(nmdim,nt)
      real rsphere

      ntrunc = -1.5 + 0.5*sqrt(9.-8.*(1.-float(nmdim)))

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      n1=m
      if (m .eq. 1) n1=2
      do n=n1,ntrunc+1
         nm = nmstrt + n - m + 1
         dataspec_ilap(nm,i) =
     *   -(rsphere**2/(float(n)*float(n-1)))*dataspec(nm,i)
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      dataspec_ilap(1,i) = (0.,0.)
      enddo

      return
      end
