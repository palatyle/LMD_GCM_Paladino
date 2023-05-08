      subroutine scopyi(n,sx,incx,sy,incy)
c
      IMPLICIT NONE
c
      integer n,incx,incy,ix,iy,i
      integer sx((n-1)*incx+1),sy((n-1)*incy+1)
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
