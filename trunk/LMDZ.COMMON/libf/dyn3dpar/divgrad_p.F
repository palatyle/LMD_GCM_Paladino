      SUBROUTINE divgrad_p (klevel,h, lh, divgra_out )
      USE parallel_lmdz
      USE times
      IMPLICIT NONE
c
c=======================================================================
c
c  Auteur :   P. Le Van
c  ----------
c
c                              lh
c      calcul de  (div( grad ))   de h  .....
c      h  et lh  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
c=======================================================================
c
c   declarations:
c   -------------
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comdissipn.h"
c
      INTEGER klevel
      REAL h( ip1jmp1,klevel ), divgra_out( ip1jmp1,klevel )
      REAL,SAVE :: divgra( ip1jmp1,llm )

c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)

      INTEGER  l,ij,iter,lh
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY ( ip1jmp1*klevel,h,1,divgra,1 )
      
      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, klevel      
      divgra(ijb:ije,l)=h(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
c

c
      DO 10 iter = 1,lh

      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p ( divgra,jjb,jje,jjp1,klevel,2,1,.true.,1  )

c      call exchange_Hallo(divgra,ip1jmp1,llm,0,1)
c$OMP BARRIER
c$OMP MASTER      
      call suspend_timer(timer_dissip)
      call exchange_Hallo(divgra,ip1jmp1,llm,1,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER       
      CALL    grad_p (klevel,divgra, ghx  , ghy          )

c$OMP BARRIER
c$OMP MASTER   
      call suspend_timer(timer_dissip)
      call exchange_Hallo(ghy,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER            

      CALL  diverg_p (klevel,  ghx , ghy  , divgra       )

      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p( divgra,jjb,jje,jjp1,klevel,2,1,.true.,1)

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,klevel
      DO 4  ij = ijb, ije
      divgra_out( ij,l ) = - cdivh * divgra( ij,l )
   4  CONTINUE
   5  CONTINUE
c$OMP END DO NOWAIT
c
  10  CONTINUE
      RETURN
      END
