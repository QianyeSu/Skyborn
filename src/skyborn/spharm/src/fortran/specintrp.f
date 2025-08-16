      subroutine specintrp(rlon,ntrunc,datnm,scrm,pnm,ob)

c given spectral coeffs of a function (datnm) truncated at T(ntrunc),
c longitude of pt in radians (rlon),
c and associated legendre polynomials for that pt (pnm), return
c interpolated function at that pt (ob).

      integer ntrunc,m,n,nmstrt
      real ob,rlon,pnm((ntrunc+1)*(ntrunc+2)/2)
      complex datnm((ntrunc+1)*(ntrunc+2)/2)
      complex scrm(ntrunc+1)

      mwaves = ntrunc + 1
      nmstrt = 0
      do m = 1, mwaves
         scrm(m) = cmplx(0.,0.)
         DO n = 1, mwaves-m+1
            nm = nmstrt + n
            scrm(m) = scrm(m)  +  datnm(nm) * pnm(NM)
         enddo
         nmstrt = nmstrt + mwaves-m+1
      enddo

      ob = scrm(1)
      do m=2,mwaves
         ob = ob +
     *        2.0*real(scrm(m))*cos(float(m-1)*rlon)-
     *        2.0*aimag(scrm(m))*sin(float(m-1)*rlon)
      enddo

      return
      end
