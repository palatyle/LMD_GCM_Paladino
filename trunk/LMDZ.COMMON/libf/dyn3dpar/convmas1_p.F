      SUBROUTINE convmas1_p (pbaru, pbarv, convm )
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
      REAL, target :: convm(  ip1jmp1,llm )
      INTEGER   l,ij

      INTEGER ijb,ije,jjb,jje
 
      
c-----------------------------------------------------------------------
c    ....  calcul de - (d(pbaru)/dx + d(pbarv)/dy ) ......

      CALL  convflu_p( pbaru, pbarv, llm, convm )

c-----------------------------------------------------------------------
c   filtrage:
c   ---------
       
       jjb=jj_begin
       jje=jj_end+1
       if (pole_sud) jje=jj_end
 
       CALL filtreg_p( convm, jjb, jje, jjp1, llm, 2, 2, .true., 1 )

c    integration de la convergence de masse de haut  en bas ......
c
      RETURN
      END
