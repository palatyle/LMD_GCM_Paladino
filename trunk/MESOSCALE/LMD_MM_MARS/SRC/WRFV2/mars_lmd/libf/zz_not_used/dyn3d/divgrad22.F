      SUBROUTINE divgrad2_mars (klevel,h, lh, divgra )
c
c  P. Le Van
c
c   ************************************************************
c
c                                lh
c        calcul de  (div( grad ))   de h  .....
c   ************************************************************
c      h  et lh  sont des arguments  d'entree pour le s-prog
c      divgra     est  un argument  de sortie pour le s-prog
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comdissipn.h"

      INTEGER klevel
      REAL h( ip1jmp1,klevel ), divgra( ip1jmp1,klevel )
c
      REAL ghy(ip1jm,llm), ghx(ip1jmp1,llm)
      REAL signe, nudivgrs

      INTEGER  l,ij,iter,lh

      EXTERNAL  filtreg
      EXTERNAL  SCOPY, grad, covcont, diverg, divergst


      signe    = (-1.)**lh
      nudivgrs = signe * cdivh

      CALL SCOPY ( ip1jmp1*klevel,h,1,divgra,1 )

      DO 2 iter = 1, lh -1
      CALL grad (klevel,divgra, ghx  , ghy          )
      CALL divergst (klevel,  ghx , ghy  , divgra       )
   2  CONTINUE

      CALL filtreg (divgra,jjp1,klevel,2,1,.true.,1)
      CALL grad  (klevel,divgra, ghx, ghy      )
      CALL covcont (klevel,ghx   , ghy, ghx ,ghy )
      CALL diverg  (klevel,ghx   , ghy,divgra    )

      CALL filtreg ( divgra,jjp1,klevel,2,1,.true.,1)

      DO 5 l = 1,klevel
      DO 4  ij = 1, ip1jmp1
      divgra( ij,l ) = nudivgrs * divgra( ij,l )
   4  CONTINUE
   5  CONTINUE

      RETURN
      END
