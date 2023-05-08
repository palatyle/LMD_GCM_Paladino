      SUBROUTINE laplacien_gam_p ( klevel, cuvsga, cvusga, unsaigam ,
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
      USE parallel_lmdz
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

      INTEGER :: ijb,ije
      INTEGER :: l      
c
c
c   ...  cvuscugam  = ( cvu/ cu ) ** (- gamdissip )
c   ...  cuvscvgam  = ( cuv/ cv ) ** (- gamdissip )  calcules dans inigeom  ..
c   ...  unsairegam =  1. /  aire ** (- gamdissip )
c

c      CALL SCOPY ( ip1jmp1 * klevel, teta, 1, divgra, 1 )
      
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      if (pole_nord) ijb=ij_begin
      if (pole_sud ) ije=ij_end
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,klevel      
        divgra(ijb:ije,l)=teta(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT

c
      CALL   grad_p ( klevel, divgra, ghx, ghy )
c
      CALL  diverg_gam_p ( klevel, cuvsga, cvusga,  unsaigam  ,
     *                 unsapolnga, unsapolsga, ghx , ghy , divgra )

c


      RETURN
      END
