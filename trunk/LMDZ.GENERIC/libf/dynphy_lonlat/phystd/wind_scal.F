      SUBROUTINE wind_scal(pbaru,pbarv,us,vs)
c=======================================================================
c
c
c   Subject:
c   ------
c   On ramene les flux de masse /vents  aux points scalaires.
c
c=======================================================================
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c   Arguments:
c   ----------

      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
      REAL us(ip1jmp1,llm), vs (ip1jmp1,llm)

c   Local:
c   ------

      INTEGER ij,l

c-----------------------------------------------------------------------

c   transport zonal:
c   ----------------
      DO l=1,llm
	 DO ij=2,ip1jmp1
            us(ij,l)=.5*(pbaru(ij,l)+pbaru(ij-1,l))
	 ENDDO
      ENDDO
      CALL SCOPY(jjp1*llm,us(iip1,1),iip1,us(1,1),iip1)


c   Transport meridien:
c   -------------------
      DO l=1,llm
         DO ij=iip2,ip1jm
	    vs(ij,l)=.5*(pbarv(ij,l)+pbarv(ij-iip1,l))
	 ENDDO
	 DO ij=1,iip1
	    vs(ij,l)=0.
	    vs(ip1jm+ij,l)=0.
	 ENDDO
      ENDDO

      RETURN
      END
