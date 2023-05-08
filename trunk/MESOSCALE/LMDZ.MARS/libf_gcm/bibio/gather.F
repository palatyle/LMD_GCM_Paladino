      SUBROUTINE GATHER(n,a,b,index)
c
      IMPLICIT NONE
C
      INTEGER n,index(n),i
      REAL a(n),b(n)
c
      DO 100 i=1,n
        a(i)=b(index(i))
100   CONTINUE
c
      RETURN
      END
c
c     voir aussi vec_gather(v,vindices,count,r)...p.11-14
c  
