!
! $Header$
!
      SUBROUTINE haut2bas(klon, klev, varB2H, varH2B)
      IMPLICIT NONE
c
      INTEGER klon, klev
      REAL varB2H(klon, klev), varH2B(klon, klev)
      INTEGER i, k, kinv
c
      DO k=1,klev  
       kinv=klev-k+1 
       DO i=1,klon
        varH2B(i,k)=varB2H(i,kinv)
       ENDDO
      ENDDO
c
      RETURN
      END
