      SUBROUTINE geopot_p ( ngrid, teta, pk, pks, phis, phi )
      USE parallel_lmdz
      IMPLICIT NONE
      
      
c=======================================================================
c
c   Auteur:  P. Le Van
c   -------
c
c   Objet:
c   ------
c
c    *******************************************************************
c    ....   calcul du geopotentiel aux milieux des couches    .....
c    *******************************************************************
c
c     ....   l'integration se fait de bas en haut  ....
c
c     .. ngrid,teta,pk,pks,phis sont des argum. d'entree pour le s-pg ..
c              phi               est un  argum. de sortie pour le s-pg .
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"

c   Arguments:
c   ----------
      INTEGER ngrid
      REAL teta(ngrid,llm),pks(ngrid),phis(ngrid),pk(ngrid,llm) ,
     *       phi(ngrid,llm)


c   Local:
c   ------
      
      INTEGER  l, ij,ijb,ije


c-----------------------------------------------------------------------
c     calcul de phi au niveau 1 pres du sol  .....
      ijb=ij_begin
      ije=ij_end+iip1
      
      IF (pole_sud)  ije=ij_end

      DO  ij  = ijb, ije
      phi( ij,1 ) = phis( ij ) + teta(ij,1) * ( pks(ij) - pk(ij,1) )
      ENDDO

c     calcul de phi aux niveaux superieurs  .......

      DO  l = 2,llm
        DO  ij    = ijb,ije
        phi(ij,l) = phi(ij,l-1) + 0.5 * ( teta(ij,l)  + teta(ij,l-1) ) 
     *                              *   (  pk(ij,l-1) -  pk(ij,l)    )
        ENDDO
      ENDDO

      RETURN
      END
