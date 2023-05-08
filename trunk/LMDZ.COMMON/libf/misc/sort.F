!
! $Header$
!
C
C
      SUBROUTINE sort(n,d)
c
c     P.Le Van
c      
c...  cette routine met le tableau d  dans l'ordre croissant  ....
cc   ( pour avoir l'ordre decroissant,il suffit de remplacer l'instruc
c      tion  situee + bas  IF(d(j).LE.p)  THEN     par
c                           IF(d(j).GE.p)  THEN
c

      INTEGER n
      REAL d(n) , p
      INTEGER i,j,k

      DO i=1,n-1
        k=i
        p=d(i)
        DO j=i+1,n
         IF(d(j).LE.p) THEN
           k=j
           p=d(j)
         ENDIF
        ENDDO

       IF(k.ne.i) THEN
         d(k)=d(i)
         d(i)=p
       ENDIF
      ENDDO

       RETURN
       END
