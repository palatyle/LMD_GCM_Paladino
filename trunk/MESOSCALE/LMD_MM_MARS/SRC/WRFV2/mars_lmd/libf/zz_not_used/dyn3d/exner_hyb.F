      SUBROUTINE  exner_hyb ( ngrid, ps, p,beta, pks, pk, pkf )
c
c     Auteurs :  F. Forget , Y. Wanherdrick
c P.Le Van  , Fr. Hourdin  .
c    ..........
c
c    ....  ngrid, ps,p             sont des argum.d'entree  au sous-prog ...
c    .... beta, pks,pk,pkf   sont des argum.de sortie au sous-prog ...
c
c   ************************************************************************
c    Calcule la fonction d'Exner pk = Cp * (p/preff) ** kappa , aux milieux des 
c    couches .   Pk(l) sera calcule aux milieux  des couches l ,entre les
c    pressions p(l) et p(l+1) ,definis aux interfaces des llm couches .
c   ************************************************************************
c    .. N.B : Au sommet de l'atmosphere,  p(llm+1) = 0. , et ps et pks sont
c    la pression et la fonction d'Exner  au  sol  .
c
c     WARNING : CECI est une version speciale de exner_hyb originale
c               Utilis‰ dans la version martienne pour pouvoir 
c               tourner avec des coordonn‰es verticales complexe
c              => Il ne verifie PAS la condition la proportionalit‰ en 
c              ‰nergie totale/ interne / potentielle (F.Forget 2001)
c    ( voir note de Fr.Hourdin )  ,
c
      IMPLICIT NONE
c
#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comgeom.h"
#include "comvert.h"
#include "serre.h"

      INTEGER  ngrid
      REAL p(ngrid,llmp1),pk(ngrid,llm),pkf(ngrid,llm)
      REAL ps(ngrid),pks(ngrid), beta(ngrid,llm)

c    .... variables locales   ...

      INTEGER l, ij
      REAL dum1

      REAL ppn(iim),pps(iim)
      REAL xpn, xps
      REAL SSUM
      EXTERNAL filtreg, SSUM
      
c     -------------
c     Calcul de pks
c     -------------
   
      DO   ij  = 1, ngrid
        pks(ij) = cpp * ( ps(ij)/preff ) ** kappa
      ENDDO

      DO  ij   = 1, iim
        ppn(ij) = aire(   ij   ) * pks(  ij     )
        pps(ij) = aire(ij+ip1jm) * pks(ij+ip1jm )
      ENDDO
      xpn      = SSUM(iim,ppn,1) /apoln
      xps      = SSUM(iim,pps,1) /apols

      DO ij   = 1, iip1
        pks(   ij     )  =  xpn
        pks( ij+ip1jm )  =  xps
      ENDDO
c
c
c    .... Calcul de pk  pour la couche l 
c    --------------------------------------------
c
      dum1 = cpp * (2*preff)**(-kappa) 
      DO l = 1, llm-1
        DO   ij   = 1, ngrid
         pk(ij,l) = dum1 * (p(ij,l) + p(ij,l+1))**kappa
        ENDDO
      ENDDO

c    .... Calcul de pk  pour la couche l = llm ..
c    (on met la meme distance (en log pression)  entre Pk(llm)
c    et Pk(llm -1) qu'entre Pk(llm-1) et Pk(llm-2)

      DO   ij   = 1, ngrid
         pk(ij,llm) = pk(ij,llm-1)**2 / pk(ij,llm-2)
      ENDDO


c    calcul de pkf
c    -------------
      CALL SCOPY   ( ngrid * llm, pk, 1, pkf, 1 )
      CALL filtreg ( pkf, jmp1, llm, 2, 1, .TRUE., 1 )
      
c    EST-CE UTILE ?? : calcul de beta
c    --------------------------------
      DO l = 2, llm
        DO   ij   = 1, ngrid
          beta(ij,l) = pk(ij,l) / pk(ij,l-1)   
        ENDDO
      ENDDO

      RETURN
      END
