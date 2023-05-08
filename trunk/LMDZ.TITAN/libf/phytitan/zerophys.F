      SUBROUTINE zerophys(n,x)
      IMPLICIT NONE
      INTEGER n
      REAL x(n)
      INTEGER i

      DO i=1,n
         x(i)=0.
      ENDDO
      RETURN
      END
