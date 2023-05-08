      SUBROUTINE divergst(klevel,x,y,div)
      IMPLICIT NONE
c
c     P. Le Van
c
c  ******************************************************************
c  ... calcule la divergence a tous les niveaux d'1 vecteur de compos. x et y...
c           x et y  etant des composantes contravariantes   ...
c  ****************************************************************
c      x  et  y  sont des arguments  d'entree pour le s-prog
c        div      est  un argument  de sortie pour le s-prog
c
c
c   -------------------------------------------------------------------
c
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      INTEGER klevel
      REAL x( ip1jmp1,klevel ),y( ip1jm,klevel ),div( ip1jmp1,klevel )
      INTEGER ij,l,i
      REAL aiy1( iip1 ) , aiy2( iip1 )
      REAL sumypn,sumyps

      REAL SSUM
      EXTERNAL SSUM
c
c
      DO 10 l = 1,klevel
c
      DO 1 ij = iip2, ip1jm - 1
      div( ij + 1, l ) = x(ij+1,l) - x(ij,l)+ y(ij-iim,l)-y(ij+1,l)
   1  CONTINUE
c
c     ....  correction pour  div( 1,j,l)  ......
c     ....   div(1,j,l)= div(iip1,j,l) ....
c
CDIR$ IVDEP
      DO 3 ij = iip2,ip1jm,iip1
      div( ij,l ) = div( ij + iim,l )
   3  CONTINUE
c
c     ....  calcul  aux poles  .....
c
c
      DO 5 i  = 1,iim
      aiy1(i)= y(i,l)
      aiy2(i)= y(i+ip1jmi1,l)
   5  CONTINUE
      sumypn = SSUM ( iim,aiy1,1 )
      sumyps = SSUM ( iim,aiy2,1 )
      DO 7 i = 1,iip1
      div(     i    , l ) = - sumypn/iim
      div( i + ip1jm, l ) =   sumyps/iim
   7  CONTINUE
c
  10  CONTINUE
      RETURN
      END
