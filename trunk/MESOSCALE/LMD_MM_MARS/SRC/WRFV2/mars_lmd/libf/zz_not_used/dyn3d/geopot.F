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
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comvert.h"

c   Arguments:
c   ----------

      INTEGER ngrid
      REAL teta(ngrid,llm),pks(ngrid),phis(ngrid),pk(ngrid,llm) ,
     *       phi(ngrid,llm)
      real zz(ngrid)


c   Local:
c   ------

      INTEGER  l, ij, count
      logical flag 


c-----------------------------------------------------------------------
c     calcul de phi au niveau 1 pres du sol  .....

      flag=.false.
      count=0
c     call dump2d(iip1,jjp1,teta(:,1),'TETA COUCHE 1')
c     zz(:)=pks(:)-pk(:,1)
c     call dump2d(iip1,jjp1,zz(:),'DPK SURF-1')

      DO   ij  = 1, ngrid
          phi( ij,1 ) = phis( ij ) + teta(ij,1) * ( pks(ij) - pk(ij,1) )
      ENDDO

c     calcul de phi aux niveaux superieurs  .......

      DO  l = 2,llm
        DO  ij    = 1,ngrid
        phi(ij,l) = phi(ij,l-1) + 0.5 * ( teta(ij,l)  + teta(ij,l-1) ) 
     *                              *   (  pk(ij,l-1) -  pk(ij,l)    )
        ENDDO
      ENDDO

      RETURN
      END
