!
! $Id: minmaxqfi.F90 1279 2009-12-10 09:02:56Z fairhead $
!
SUBROUTINE minmaxqfi(zq,qmin,qmax,comment)
  USE dimphy
  IMPLICIT NONE

! Entrees
  REAL,DIMENSION(klon,klev), INTENT(IN)   :: zq
  REAL,INTENT(IN)                         :: qmin,qmax
  CHARACTER(LEN=*),INTENT(IN)             :: comment

! Local  
  INTEGER,DIMENSION(klon)     :: jadrs 
  INTEGER                     :: i, jbad, k
  
  DO k = 1, klev
     jbad = 0
     DO i = 1, klon
        IF (zq(i,k).GT.qmax .OR. zq(i,k).LT.qmin) THEN
           jbad = jbad + 1
           jadrs(jbad) = i
        ENDIF
     ENDDO
     IF (jbad.GT.0) THEN
        WRITE(*,*)comment
        DO i = 1, jbad
           WRITE(*,*) "i,k,q=", jadrs(i),k,zq(jadrs(i),k)
        ENDDO
     ENDIF
  ENDDO
  
END SUBROUTINE minmaxqfi
