










!
! $Header$
!
      subroutine scopy(n,sx,incx,sy,incy)
c
      IMPLICIT NONE
c
      integer n,incx,incy,ix,iy,i
      real sx((n-1)*incx+1),sy((n-1)*incy+1)
c
      iy=1
      ix=1
      do 10 i=1,n
         sy(iy)=sx(ix)
         ix=ix+incx
         iy=iy+incy
10    continue
c
      return
      end

      function ssum(n,sx,incx)
c
      IMPLICIT NONE
c
      integer n,incx,i,ix
      real ssum,sx((n-1)*incx+1)
c
      ssum=0.
      ix=1
      do 10 i=1,n
         ssum=ssum+sx(ix)
         ix=ix+incx
10    continue
c
      return
      end
