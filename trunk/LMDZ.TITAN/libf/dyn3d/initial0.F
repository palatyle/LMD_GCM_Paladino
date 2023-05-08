      SUBROUTINE initial0(n,x)
      IMPLICIT NONE
      INTEGER n,i
      REAL x(n)
      DO 10 i=1,n
         x(i)=0.
10    CONTINUE
      RETURN
      END
