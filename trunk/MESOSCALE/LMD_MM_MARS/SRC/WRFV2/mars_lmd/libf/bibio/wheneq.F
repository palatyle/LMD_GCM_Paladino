      SUBROUTINE wheneq(n,x,incx,target,index,nn)
      IMPLICIT NONE
c
      integer n, incx, nn
      integer index(n)
c
      integer target, x((n-1)*incx+1)
cccccc      real target, x((n-1)*incx+1)
c
      integer i, ix
c
      ix = 1
      nn = 0
      do 10 i = 1, n
         if (x(ix).EQ.target) then
            index(i) = ix
            nn = nn + 1
         endif
         ix = ix + incx
   10 continue
c
      return
      end
