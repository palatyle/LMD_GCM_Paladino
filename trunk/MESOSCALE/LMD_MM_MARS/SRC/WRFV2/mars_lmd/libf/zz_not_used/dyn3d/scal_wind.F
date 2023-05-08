      SUBROUTINE scal_wind(xus,xvs,xu,xv)
c=======================================================================
c
c
c   Subject:
c   ------
c On passe  les variable xus, xvs  aux points de vent u et v (xu et xv)
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

      REAL xu(iip1,jjp1,llm),xv(iip1,jjm,llm)
      REAL xus(iip1,jjp1,llm), xvs (iip1,jjp1,llm)

c   Local:
c   ------

      INTEGER i,j,l

c-----------------------------------------------------------------------

c   transport zonal:
c   ----------------
      DO l=1,llm
        Do j=1,jjp1
	      DO i=1,iim
            xu(i,j,l)=0.5*(xus(i,j,l)+xus(i+1,j,l))
	      ENDDO
          xu(iip1,j,l)=xu(1,j,l)
	    ENDDO
      ENDDO


c   Transport meridien:
c   -------------------
      DO l=1,llm
         DO j=1,jjm
           do i=1 ,iip1
	         xv(i,j,l)=.5*(xvs(i,j,l)+xvs(i,j+1,l))
           end do
	     ENDDO
	  ENDDO

      RETURN
      END
