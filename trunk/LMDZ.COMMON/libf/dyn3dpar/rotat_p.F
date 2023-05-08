      SUBROUTINE rotat_p (klevel, x, y, rot )
c
c     Auteur : P.Le Van 
c**************************************************************
c.  calcule le rotationnel
c     a tous les niveaux d'1 vecteur de comp. x et y ..
c       x  et  y etant des composantes  covariantes  ...
c********************************************************************
c   klevel, x  et y   sont des arguments d'entree pour le s-prog
c        rot          est  un argument  de sortie pour le s-prog
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
c
c   .....  variables en arguments  ......
c
      INTEGER klevel
      REAL rot( ip1jm,klevel )
      REAL x( ip1jmp1,klevel ), y( ip1jm,klevel )
c
c  ...   variables  locales  ...
c
      INTEGER  l, ij
      INTEGER :: ijb,ije
c
c
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO  10 l = 1,klevel
c
        DO   ij = ijb, ije - 1
         rot( ij,l )  =    y( ij+1 , l )  -  y( ij,l )   +
     *                   x(ij +iip1, l )  -  x( ij,l )  
        ENDDO
c
c    .... correction pour rot( iip1,j,l)  ....
c    ....   rot(iip1,j,l)= rot(1,j,l) ...
CDIR$ IVDEP
        DO  ij = ijb+iip1-1, ije, iip1
         rot( ij,l ) = rot( ij -iim,l )
        ENDDO
c
  10  CONTINUE
c$OMP END DO NOWAIT
ccc        CALL filtreg( rot, jjm, klevel, 2, 2, .FALSE., 1 )
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO l = 1, klevel
          DO ij = ijb, ije
           rot(ij,l) = rot(ij,l) * unsairez(ij)
          ENDDO
        ENDDO
c$OMP END DO NOWAIT
c
c
      RETURN
      END
