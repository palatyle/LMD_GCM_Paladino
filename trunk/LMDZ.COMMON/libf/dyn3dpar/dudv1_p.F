      SUBROUTINE dudv1_p ( vorpot, pbaru, pbarv, du, dv )
      USE parallel_lmdz
      IMPLICIT NONE
c
c-----------------------------------------------------------------------
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c   calcul du terme de  rotation
c   ce terme est ajoute a  d(ucov)/dt et a d(vcov)/dt  ..
c   vorpot, pbaru et pbarv sont des arguments d'entree  pour le s-pg ..
c   du  et dv              sont des arguments de sortie pour le s-pg ..
c
c-----------------------------------------------------------------------

#include "dimensions.h"
#include "paramet.h"

      REAL vorpot( ip1jm,llm ) ,pbaru( ip1jmp1,llm ) ,
     *     pbarv( ip1jm,llm ) ,du( ip1jmp1,llm ) ,dv( ip1jm,llm )
      INTEGER  l,ij,ijb,ije
c
c
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO 10 l = 1,llm
c
      ijb=ij_begin
      ije=ij_end
      
      if (pole_nord) ijb=ij_begin+iip1
      if (pole_sud)  ije=ij_end-iip1
      
      DO 2  ij = ijb, ije-1 
      du( ij,l ) = 0.125 *(  vorpot(ij-iip1, l) + vorpot( ij, l)  ) *
     *                    (   pbarv(ij-iip1, l) + pbarv(ij-iim,  l) +
     *                        pbarv(   ij  , l) + pbarv(ij+ 1 ,  l)   )
   2  CONTINUE
   
 
c
      if (pole_nord) ijb=ij_begin
      
      DO 3 ij = ijb, ije-1 
      dv( ij+1,l ) = - 0.125 *(  vorpot(ij, l)  + vorpot(ij+1, l)  ) *
     *                        (   pbaru(ij, l)  +  pbaru(ij+1   , l) +
     *                       pbaru(ij+iip1, l)  +  pbaru(ij+iip2, l)  )
   3  CONTINUE
c
c    .... correction  pour  dv( 1,j,l )  .....
c    ....   dv(1,j,l)= dv(iip1,j,l) ....
c
CDIR$ IVDEP
      DO 4 ij = ijb, ije, iip1
      dv( ij,l ) = dv( ij + iim, l )
   4  CONTINUE
c
  10  CONTINUE
c$OMP END DO NOWAIT
      RETURN
      END
