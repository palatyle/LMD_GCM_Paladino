!
! $Header$
!
      SUBROUTINE geopot (ngrid, teta, pk, pks, phis, phi )
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
c  This computation (with teta = cp T / pk !) is identical to 
c     delta phi = R/RMD T/p delta p         (r=R/RMD=cpp*kappa)
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

      INTEGER  l, ij


c-----------------------------------------------------------------------
c     calcul de phi au niveau 1 pres du sol  .....

      DO   1  ij  = 1, ngrid
      phi( ij,1 ) = phis( ij ) + teta(ij,1) * ( pks(ij) - pk(ij,1) )
   1  CONTINUE

c     calcul de phi aux niveaux superieurs  .......

      DO  l = 2,llm
        DO  ij    = 1,ngrid
        phi(ij,l) = phi(ij,l-1) + 0.5 * ( teta(ij,l)  + teta(ij,l-1) ) 
     *                              *   (  pk(ij,l-1) -  pk(ij,l)    )
        ENDDO
      ENDDO

      RETURN
      END
