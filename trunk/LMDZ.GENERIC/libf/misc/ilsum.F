      FUNCTION ilsum(n, logic, inc)
c
      IMPLICIT NONE
c
      integer ilsum
      integer n, inc
      logical logic((n-1)*inc+1)
c
      integer i, in
c
      ilsum = 0.
      in = 1
      do 10 i = 1, n
         if (logic(in)) ilsum = ilsum + 1
         in=in+inc
10    continue
c
      return
      end
