      SUBROUTINE addit(np,px,py,pz)
      IMPLICIT NONE

      INTEGER np
      REAL px(np),py(np),pz(np)

      INTEGER i

      DO i=1,np
	 pz(i)=px(i)+py(i)
      ENDDO

      RETURN
      END
