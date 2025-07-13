      subroutine onedtotwod_vrtdiv(vrtspec,divspec,br,bi,cr,ci,
     *                             nlat,nmdim,nt,rsphere)
      integer nlat,nmdim,nt,nmstrt,i
      real br(nlat,nlat,nt),bi(nlat,nlat,nt),
     *     cr(nlat,nlat,nt),ci(nlat,nlat,nt),rsphere
      complex vrtspec(nmdim,nt),divspec(nmdim,nt)

      ntrunc = -1.5 + 0.5*sqrt(9.-8.*(1.-float(nmdim)))

      scale = 0.5

      do i=1,nt
      nmstrt = 0
      do m=1,ntrunc+1
      do n=m,ntrunc+1
         nm = nmstrt + n - m + 1
         br(m,n,i)=
     *   -(rsphere/sqrt(float(n)*float(n-1)))*real(divspec(nm,i)/scale)
         bi(m,n,i)=
     *   -(rsphere/sqrt(float(n)*float(n-1)))*aimag(divspec(nm,i)/scale)
         cr(m,n,i)=
     *   (rsphere/sqrt(float(n)*float(n-1)))*real(vrtspec(nm,i)/scale)
         ci(m,n,i)=
     *   (rsphere/sqrt(float(n)*float(n-1)))*aimag(vrtspec(nm,i)/scale)
      enddo
      nmstrt = nmstrt + ntrunc - m + 2
      enddo
      enddo

      return
      end
