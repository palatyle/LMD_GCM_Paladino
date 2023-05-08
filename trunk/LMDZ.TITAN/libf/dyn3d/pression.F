      SUBROUTINE pression( ngrid, ap, bp, ps, p )
c

c      Auteurs : P. Le Van , Fr.Hourdin  .

c  ************************************************************************
c     Calcule la pression p(l) aux differents niveaux l = 1 ( niveau du
c     sol) a l = llm +1, ces niveaux correspondant aux interfaces des (llm) 
c     couches , avec  p(ij,llm +1) = 0.  et p(ij,1) = ps(ij)  .      
c  ************************************************************************
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
c
      INTEGER ngrid
      INTEGER l,ij
 
      REAL ap(llmp1), bp(llmp1), ps(ngrid), p(ngrid,llmp1) 
      
      DO    l    = 1, llmp1
        DO  ij   = 1, ngrid
         p(ij,l) = ap(l) + bp(l) * ps(ij)
        ENDDO
      ENDDO
   
      RETURN
      END
