      SUBROUTINE covcont (klevel,ucov, vcov, ucont, vcont )
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c  *********************************************************************
c    calcul des compos. contravariantes a partir des comp.covariantes
c  ********************************************************************
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      INTEGER klevel
      REAL ucov ( ip1jmp1,klevel ),  vcov( ip1jm,klevel )
      REAL ucont( ip1jmp1,klevel ), vcont( ip1jm,klevel )
      INTEGER   l,ij


      DO 10 l = 1,klevel

      DO 2  ij = iip2, ip1jm
      ucont( ij,l ) = ucov( ij,l ) * unscu2( ij )
   2  CONTINUE

      DO 4 ij = 1,ip1jm
      vcont( ij,l ) = vcov( ij,l ) * unscv2( ij )
   4  CONTINUE

  10  CONTINUE

      RETURN
      END
