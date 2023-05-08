      SUBROUTINE nxgrarot_p (klevel,xcov, ycov, lr, grx_out, gry_out )
c   ***********************************************************
c
c    Auteur :  P.Le Van  
c
c                                 lr
c      calcul de  ( nXgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      USE parallel_lmdz
      USE times
      USE write_field_p
      IMPLICIT NONE
c
c
#include "dimensions.h"
#include "paramet.h"
#include "comdissipn.h"
c
      INTEGER klevel
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL  grx_out( ip1jmp1,klevel ),  gry_out( ip1jm,klevel )
      REAL,SAVE ::  grx( ip1jmp1,llm ),  gry( ip1jm,llm )

c
      REAL,SAVE :: rot(ip1jm,llm)

      INTEGER l,ij,iter,lr
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY ( ip1jmp1*klevel, xcov, 1, grx, 1 )
c      CALL SCOPY (  ip1jm*klevel, ycov, 1, gry, 1 )
c
      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
      DO l = 1, klevel
        grx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO 
c$OMP END DO NOWAIT      

      if(pole_sud) ije=ij_end-iip1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1, klevel
        gry(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
      
      DO 10 iter = 1,lr
c$OMP BARRIER
c$OMP MASTER
      call suspend_timer(timer_dissip)
      call exchange_Hallo(grx,ip1jmp1,llm,0,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER

      CALL  rotat_p (klevel,grx, gry, rot )
c      call write_field3d_p('rot',reshape(rot,(/iip1,jjm,llm/)))
      
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      CALL filtreg_p( rot,jjb,jje, jjm, klevel, 2,1, .false.,2)

c$OMP BARRIER
c$OMP MASTER
      call suspend_timer(timer_dissip)
      call exchange_Hallo(rot,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER
      
      CALL nxgrad_p (klevel,rot, grx, gry )
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5  l = 1, klevel
      if(pole_sud) ije=ij_end-iip1
      DO 2 ij = ijb, ije
      gry_out( ij,l ) = - gry( ij,l ) * crot
   2  CONTINUE
      if(pole_sud) ije=ij_end
      DO 3 ij = ijb, ije
      grx_out( ij,l ) = - grx( ij,l ) * crot
   3  CONTINUE
   5  CONTINUE
c$OMP END DO NOWAIT
c      call write_field3d_p('grx',reshape(grx,(/iip1,jjp1,llm/)))
c      call write_field3d_p('gry',reshape(gry,(/iip1,jjm,llm/)))
c      stop
  10  CONTINUE
      RETURN
      END
