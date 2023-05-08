      SUBROUTINE laplacien_p ( klevel, teta, divgra )
c
c     P. Le Van
c
c   ************************************************************
c    ....     calcul de  (div( grad ))   de   teta  .....
c   ************************************************************
c     klevel et teta  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      USE parallel_lmdz
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
      INTEGER :: l
c
c    ............     variables  locales      ..............
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
c    .......................................................

      
      INTEGER :: ijb,ije,jjb,jje
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
      
      jjb=jj_begin-1
      jje=jj_end+1
      if (pole_nord) jjb=jj_begin
      if (pole_sud ) jje=jj_end
      
      CALL filtreg_p( divgra,jjb,jje,jjp1, klevel,  2, 1, .TRUE., 1 )
      CALL   grad_p ( klevel,divgra,   ghx , ghy              )
      CALL  divergf_p ( klevel, ghx , ghy  , divgra           )

      RETURN
      END
