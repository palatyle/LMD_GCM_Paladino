      SUBROUTINE gradiv_p(klevel, xcov, ycov, ld, gdx_out, gdy_out )
c
c    Auteur :   P. Le Van
c
c   ***************************************************************
c
c                                ld
c       calcul  de  (grad (div) )   du vect. v ....
c
c     xcov et ycov etant les composant.covariantes de v
c   ****************************************************************
c    xcov , ycov et ld  sont des arguments  d'entree pour le s-prog
c     gdx   et  gdy     sont des arguments de sortie pour le s-prog
c
c     
      USE parallel_lmdz
      USE times
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comdissipn.h"

      INTEGER klevel
c
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL,SAVE :: gdx( ip1jmp1,llm ),   gdy( ip1jm,llm )

      REAL gdx_out( ip1jmp1,klevel ),   gdy_out( ip1jm,klevel )

      REAL,SAVE ::  div(ip1jmp1,llm)

      INTEGER l,ij,iter,ld
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY( ip1jmp1*klevel,xcov,1,gdx,1 )
c      CALL SCOPY( ip1jm*klevel,  ycov,1,gdy,1 )
      
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,klevel
        gdx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
      
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,klevel
        gdy(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT

c
      DO 10 iter = 1,ld

c$OMP BARRIER
c$OMP MASTER      
      call suspend_timer(timer_dissip)
      call exchange_Hallo(gdy,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER      
c$OMP BARRIER

      CALL  diverg_p( klevel,  gdx , gdy, div          )
      
      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p( div,jjb,jje, jjp1, klevel, 2,1, .true.,2 )
      
c      call exchange_Hallo(div,ip1jmp1,llm,0,1)

c$OMP BARRIER
c$OMP MASTER       
      call suspend_timer(timer_dissip)
      call exchange_Hallo(div,ip1jmp1,llm,1,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER
      
      CALL    grad_p( klevel,  div, gdx, gdy           )
c

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO 5  l = 1, klevel
      
      if(pole_sud) ije=ij_end
      DO 3 ij = ijb, ije
        gdx_out( ij,l ) = - gdx( ij,l ) * cdivu
   3  CONTINUE
   
      if(pole_sud) ije=ij_end-iip1
      DO 4 ij = ijb, ije
        gdy_out( ij,l ) = - gdy( ij,l ) * cdivu
   4  CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
c
  10  CONTINUE
      RETURN
      END
