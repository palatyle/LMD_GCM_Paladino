      SUBROUTINE  grad_p(klevel, pg,pgx,pgy )
c
c      P. Le Van
c
c    ******************************************************************
c     .. calcul des composantes covariantes en x et y du gradient de g
c
c    ******************************************************************
c             pg        est un   argument  d'entree pour le s-prog
c       pgx  et  pgy    sont des arguments de sortie pour le s-prog
c
      USE parallel_lmdz
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
      INTEGER klevel
      REAL  pg( ip1jmp1,klevel )
      REAL pgx( ip1jmp1,klevel ) , pgy( ip1jm,klevel )
      INTEGER  l,ij
      INTEGER :: ijb,ije,jjb,jje
c
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 6 l = 1,klevel
c
      ijb=ij_begin
      ije=ij_end
      DO 2  ij = ijb, ije - 1
        pgx( ij,l ) = pg( ij +1,l ) - pg( ij,l )
   2  CONTINUE
c
c    .... correction pour  pgx(ip1,j,l)  ....
c    ...    pgx(iip1,j,l)= pgx(1,j,l)  ....
CDIR$ IVDEP
      DO 3  ij = ijb+iip1-1, ije, iip1
        pgx( ij,l ) = pgx( ij -iim,l )
   3  CONTINUE
c
      ijb=ij_begin-iip1
      ije=ij_end
      if (pole_nord) ijb=ij_begin
      if (pole_sud)  ije=ij_end-iip1
      
      DO 4 ij = ijb,ije
        pgy( ij,l ) = pg( ij,l ) - pg( ij +iip1,l )
   4  CONTINUE
c
   6  CONTINUE
c$OMP END DO NOWAIT

      RETURN
      END
