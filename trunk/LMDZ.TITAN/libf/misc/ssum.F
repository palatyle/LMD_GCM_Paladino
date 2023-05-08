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
