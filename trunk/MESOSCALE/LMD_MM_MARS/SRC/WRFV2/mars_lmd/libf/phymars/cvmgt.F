      FUNCTION cvmgt(x1,x2,l)
      IMPLICIT NONE

      REAL x1,x2,cvmgt
      LOGICAL l

      IF(l) then 
        cvmgt=x1
      ELSE
        cvmgt=x2
      ENDIF

      RETURN
      END
C
