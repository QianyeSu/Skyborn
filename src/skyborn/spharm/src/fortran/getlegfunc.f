      subroutine getlegfunc(legfunc,lat,ntrunc)
      integer ntrunc,m,n,nmstrt
      real lat,theta,pi,
     *     legfunc((ntrunc+1)*(ntrunc+2)/2),cp((ntrunc/2)+1)

      pi = 4.*atan(1.0)
      theta = 0.5*pi-(pi/180.)*lat
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         call alfk(n-1,m-1,cp)
         call lfpt(n-1,m-1,theta,cp,legfunc(nm))
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo

      return
      end
