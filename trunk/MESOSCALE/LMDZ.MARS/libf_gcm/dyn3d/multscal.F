      SUBROUTINE multscal(n,x,lambda,y)
      IMPLICIT NONE
      INTEGER i,n
      REAL lambda,x(n),y(n)
      DO 10 i=1,n
         y(i)=x(i)*lambda
10    CONTINUE
      RETURN
      END
