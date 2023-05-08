      FUNCTION sdot (n, x,incx, y,incy)
      IMPLICIT NONE
c
      integer n, incx, incy
      real sdot
      real x((n-1)*incx+1), y((n-1)*incy+1)
c
      integer i, ix, iy
c
      sdot = 0.0
      ix = 1
      iy = 1
c
      do 10 i = 1, n
         sdot = sdot + x(ix) * y(iy)
         ix = ix + incx
         iy = iy + incy
   10 continue
c
      return
      end
