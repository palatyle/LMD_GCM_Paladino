      SUBROUTINE laplacien ( klevel, teta, divgra )
c
c     P. Le Van
c
c   ************************************************************
c    ....     calcul de  (div( grad ))   de   teta  .....
c   ************************************************************
c     klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c
c    .........      variables  en arguments   ..............
c
      INTEGER klevel
      REAL teta( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
c
c    ............     variables  locales      ..............
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c    .......................................................

      EXTERNAL  SCOPY, grad, divergf, filtreg

c
      CALL SCOPY ( ip1jmp1 * klevel, teta, 1, divgra, 1 )

      CALL filtreg( divgra,  jjp1, klevel,  2, 1, .TRUE., 1 )
      CALL grad ( klevel,divgra,   ghx , ghy              )
      CALL divergf ( klevel, ghx , ghy  , divgra           )

      RETURN
      END
