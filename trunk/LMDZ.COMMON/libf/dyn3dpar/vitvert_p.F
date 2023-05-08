      SUBROUTINE vitvert_p ( convm , w )
c
      USE parallel_lmdz
      USE comvert_mod, ONLY: bp
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van , F. Hourdin .
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c  .... calcul de la vitesse verticale aux niveaux sigma  ....
c    *******************************************************************
c     convm   est un argument  d'entree pour le s-pg  ......
c       w     est un argument de sortie pour le s-pg  ......
c
c    la vitesse verticale est orientee de  haut en bas .
c    au sol, au niveau sigma(1),   w(i,j,1) = 0.
c    au sommet, au niveau sigma(llm+1) , la vit.verticale est aussi
c    egale a 0. et n'est pas stockee dans le tableau w  .
c
c
c=======================================================================

#include "dimensions.h"
#include "paramet.h"

      REAL w(ip1jmp1,llm),convm(ip1jmp1,llm)
      INTEGER   l, ij,ijb,ije


      ijb=ij_begin
      ije=ij_end+iip1
      
      if (pole_sud) ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO 2  l = 1,llmm1

      DO 1 ij = ijb,ije
      w( ij, l+1 ) = convm( ij, l+1 ) - bp(l+1) * convm( ij, 1 )
   1  CONTINUE

   2  CONTINUE
c$OMP END DO
c$OMP MASTER
      DO 5 ij  = ijb,ije
      w(ij,1)  = 0.
5     CONTINUE
c$OMP END MASTER
c$OMP BARRIER
      RETURN
      END
