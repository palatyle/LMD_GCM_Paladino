      SUBROUTINE convmas2_p ( convm )
c
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van , F. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   ********************************************************************
c   .... calcul de la convergence du flux de masse aux niveaux p ...
c   ********************************************************************
c
c
c     pbaru  et  pbarv  sont des arguments d'entree pour le s-pg  ....
c      .....  convm      est  un argument de sortie pour le s-pg  ....
c
c    le calcul se fait de haut en bas, 
c    la convergence de masse au niveau p(llm+1) est egale a 0. et
c    n'est pas stockee dans le tableau convm .
c
c
c=======================================================================
c
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"

      REAL pbaru( ip1jmp1,llm ),pbarv( ip1jm,llm )
      REAL :: convm(  ip1jmp1,llm )
      INTEGER   l,ij
      INTEGER ijb,ije,jjb,jje
 
c$OMP MASTER
c    integration de la convergence de masse de haut  en bas ......
       ijb=ij_begin
       ije=ij_end+iip1
       if (pole_sud) ije=ij_end
            
      DO      l      = llmm1, 1, -1
        DO    ij     = ijb, ije
         convm(ij,l) = convm(ij,l) + convm(ij,l+1)
        ENDDO
      ENDDO
c
c$OMP END MASTER
      RETURN
      END
