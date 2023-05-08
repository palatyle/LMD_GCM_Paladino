      SUBROUTINE rotatst (klevel,x, y, rot )
c
c  P. Le Van
c
c    *****************************************************************
c     .. calcule le rotationnel a tous les niveaux d'1 vecteur de comp. x et y ..
c         x  et  y etant des composantes  covariantes  .....
c    *****************************************************************
c        x  et y     sont des arguments d'entree pour le s-prog
c        rot          est  un argument  de sortie pour le s-prog
c
      IMPLICIT NONE
c
      INTEGER klevel
#include "dimensions.h"
#include "paramet.h"

      REAL rot( ip1jm,klevel )
      REAL x( ip1jmp1,klevel ), y( ip1jm,klevel )
      INTEGER  l, ij
c
c
      DO 5 l = 1,klevel
c
      DO 1 ij = 1, ip1jm - 1
      rot( ij,l )  =  (  y( ij+1 , l )  -  y( ij,l )   +
     *                 x(ij +iip1, l )  -  x( ij,l )  )
   1  CONTINUE
c
c    .... correction pour rot( iip1,j,l)  ....
c
c    ....   rot(iip1,j,l)= rot(1,j,l) ...
CDIR$ IVDEP
      DO 2 ij = iip1, ip1jm, iip1
      rot( ij,l ) = rot( ij -iim,l )
   2  CONTINUE
c
   5  CONTINUE
      RETURN
      END
