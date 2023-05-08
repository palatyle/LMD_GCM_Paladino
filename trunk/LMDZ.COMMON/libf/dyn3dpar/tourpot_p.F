      SUBROUTINE tourpot_p ( vcov, ucov, massebxy, vorpot )
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    .........      calcul du tourbillon potentiel             .........
c    *******************************************************************
c
c     vcov,ucov,fext et pbarxyfl sont des argum. d'entree pour le s-pg .
c             vorpot            est  un argum.de sortie pour le s-pg .
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      REAL  rot( ip1jm,llm )
      REAL vcov( ip1jm,llm ),ucov( ip1jmp1,llm )
      REAL massebxy( ip1jm,llm ),vorpot( ip1jm,llm )

      INTEGER l, ij ,ije,ijb,jje,jjb


      ijb=ij_begin-iip1
      ije=ij_end
      
      if (pole_nord) ijb=ij_begin
      
      
c  ... vorpot = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) /psbarxy ..



c    ........  Calcul du rotationnel du vent V  puis filtrage  ........
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5 l = 1,llm

      if (pole_sud)  ije=ij_end-iip1-1
      DO 2 ij = ijb, ije 
      rot( ij,l ) = vcov(ij+1,l)-vcov(ij,l)+ucov(ij+iip1,l)-ucov(ij,l)
   2  CONTINUE

c    ....  correction pour  rot( iip1,j,l )  .....
c    ....     rot(iip1,j,l) = rot(1,j,l)    .....

CDIR$ IVDEP

      if (pole_sud)  ije=ij_end-iip1
     
      DO 3 ij = ijb+iip1-1, ije, iip1
      rot( ij,l ) = rot( ij -iim, l )
   3  CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
      jjb=jj_begin-1
      jje=jj_end
      
      if (pole_nord) jjb=jjb+1
      if (pole_sud)  jje=jje-1
      CALL  filtreg_p( rot, jjb,jje,jjm, llm, 2, 1, .FALSE., 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)     
      DO 10 l = 1, llm
      
      if (pole_sud)  ije=ij_end-iip1-1  
      
      DO 6 ij = ijb, ije
      vorpot( ij,l ) = ( rot(ij,l) + fext(ij) ) / massebxy(ij,l)
   6  CONTINUE

c    ..... correction pour  vorpot( iip1,j,l)  .....
c    ....   vorpot(iip1,j,l)= vorpot(1,j,l) ....
CDIR$ IVDEP
      if (pole_sud)  ije=ij_end-iip1
      DO 8 ij = ijb+iip1-1, ije, iip1
      vorpot( ij,l ) = vorpot( ij -iim,l )
   8  CONTINUE

  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
