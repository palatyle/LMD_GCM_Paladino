      FUNCTION cvmgp(x1,x2,x3)
      IMPLICIT NONE

      REAL x1,x2,x3,cvmgp

      IF(x3.ge.0) then 
        cvmgp=x1
      ELSE
        cvmgp=x2
      ENDIF

      RETURN
      END
C
