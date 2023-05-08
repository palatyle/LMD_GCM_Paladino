      SUBROUTINE  grad(klevel, pg,pgx,pgy )
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
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
      INTEGER klevel
      REAL  pg( ip1jmp1,klevel )
      REAL pgx( ip1jmp1,klevel ) , pgy( ip1jm,klevel )
      INTEGER  l,ij
c
c
      DO 6 l = 1,klevel
c
      DO 2  ij = 1, ip1jmp1 - 1
      pgx( ij,l ) = pg( ij +1,l ) - pg( ij,l )
   2  CONTINUE
c
c    .... correction pour  pgx(ip1,j,l)  ....
c    ...    pgx(iip1,j,l)= pgx(1,j,l)  ....
CDIR$ IVDEP
      DO 3  ij = iip1, ip1jmp1, iip1
      pgx( ij,l ) = pgx( ij -iim,l )
   3  CONTINUE
c
      DO 4 ij = 1,ip1jm
      pgy( ij,l ) = pg( ij,l ) - pg( ij +iip1,l )
   4  CONTINUE
c
   6  CONTINUE
      RETURN
      END
