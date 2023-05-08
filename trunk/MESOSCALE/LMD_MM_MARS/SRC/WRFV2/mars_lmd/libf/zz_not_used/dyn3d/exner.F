      SUBROUTINE exner(ngrid,pext,aire,pks,pksf)
      IMPLICIT NONE

c      .........    Version du 16/04/96     .........
c
c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "comconst.h"

      INTEGER ngrid
      REAL pext(ngrid),aire(ngrid),pks(ngrid),pksf(ngrid)

      INTEGER ig,iii

      EXTERNAL filtreg

c-----------------------------------------------------------------------

      DO 3 ig = 1,ngrid
         pks( ig )     = ( pext( ig )/ aire( ig ) ) ** kappa
         pksf( ig )    = pks( ig )
   3  CONTINUE

      CALL filtreg (pksf,jmp1,1,2,1,.true.,1)

      RETURN
      END
