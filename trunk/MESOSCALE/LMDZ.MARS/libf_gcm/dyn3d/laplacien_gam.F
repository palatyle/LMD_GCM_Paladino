      SUBROUTINE laplacien_gam ( klevel, cuvsga, cvusga, unsaigam ,
     *                        unsapolnga, unsapolsga, teta, divgra )

c  P. Le Van
c
c   ************************************************************
c
c      ....   calcul de  (div( grad ))   de   teta  .....
c   ************************************************************
c    klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

c
c    ............     variables  en arguments    ..........
c
      INTEGER klevel
      REAL teta( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
      REAL cuvsga(ip1jm) , cvusga( ip1jmp1 ),unsaigam(ip1jmp1),
     *     unsapolnga, unsapolsga
c
c    ...........    variables  locales    .................
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c    ......................................................

      EXTERNAL  filtreg
      EXTERNAL  SCOPY, grad, divergst
c
c
c   ...  cvuscugam  = ( cvu/ cu ) ** (- gamdissip )
c   ...  cuvscvgam  = ( cuv/ cv ) ** (- gamdissip )  calcules dans inigeom  ..
c   ...  unsairegam =  1. /  aire ** (- gamdissip )
c

      CALL SCOPY ( ip1jmp1 * klevel, teta, 1, divgra, 1 )
c
      CALL grad ( klevel, divgra, ghx, ghy )
c
      CALL diverg_gam ( klevel, cuvsga, cvusga,  unsaigam  ,
     *                 unsapolnga, unsapolsga, ghx , ghy , divgra )

c


      RETURN
      END
