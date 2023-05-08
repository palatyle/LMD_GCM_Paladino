      SUBROUTINE multipl(n,x1,x2,y)
      IMPLICIT NONE
c=======================================================================
c
c   multiplication de deux vecteurs
c
c=======================================================================
c
      INTEGER n,i
      REAL x1(n),x2(n),y(n)
c
      DO 10 i=1,n
         y(i)=x1(i)*x2(i)
10    CONTINUE
c
      RETURN
      END
