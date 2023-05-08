      function isamax(n,sx,incx)
c
      IMPLICIT NONE
c
      INTEGER n,i,incx,isamax,ix
      real sx((n-1)*incx+1),sxmax,abssx
c
      ix=1
      isamax=1
      sxmax=ABS(sx(1))
      do 10 i=1,n-1
         ix=ix+incx
         abssx=ABS(sx(ix))
         if(abssx.gt.sxmax) then
             sxmax=abssx
             isamax=i+1
         endif
10    continue
c
      return
      end
