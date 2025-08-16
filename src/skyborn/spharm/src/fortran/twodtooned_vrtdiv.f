      subroutine twodtooned_vrtdiv(vrtspec,divspec,br,bi,cr,ci,
     *                             nlat,ntrunc,nt,rsphere)
      integer nlat,ntrunc,nt,nmstrt,i
      complex vrtspec((ntrunc+1)*(ntrunc+2)/2,nt),
     *        divspec((ntrunc+1)*(ntrunc+2)/2,nt)
      real br(nlat,nlat,nt),bi(nlat,nlat,nt),
     *     cr(nlat,nlat,nt),ci(nlat,nlat,nt),rsphere

      scale = 0.5

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         divspec(nm,i) = -(sqrt(float(n)*float(n-1))/rsphere)*
     *                   scale*cmplx(br(m,n,i),bi(m,n,i))
         vrtspec(nm,i) =  (sqrt(float(n)*float(n-1))/rsphere)*
     *                   scale*cmplx(cr(m,n,i),ci(m,n,i))
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      enddo

      return
      end
