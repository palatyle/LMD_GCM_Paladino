      SUBROUTINE thermcell_env(ngrid,nlay,po,pt,pu,pv,pplay,  &
     &           pplev,zo,zh,zl,ztv,zthl,zu,zv,zpspsk,pqsat,lev_out)

!--------------------------------------------------------------
!thermcell_env: calcule les caracteristiques de l environnement
!necessaires au calcul des proprietes dans le thermique
!--------------------------------------------------------------

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"      
#include "iniprint.h"

      INTEGER ngrid,nlay
      REAL po(ngrid,nlay)
      REAL pt(ngrid,nlay)
      REAL pu(ngrid,nlay)
      REAL pv(ngrid,nlay)
      REAL pplay(ngrid,nlay)
      REAL pplev(ngrid,nlay+1)
      integer lev_out                           ! niveau pour les print

      REAL zo(ngrid,nlay)
      REAL zl(ngrid,nlay)
      REAL zh(ngrid,nlay)
      REAL ztv(ngrid,nlay)
      REAL zthl(ngrid,nlay)
      REAL zpspsk(ngrid,nlay)
      REAL zu(ngrid,nlay)
      REAL zv(ngrid,nlay)
      REAL pqsat(ngrid,nlay)

      INTEGER ig,ll

      real dqsat_dT
      real RLvCp

logical mask(ngrid,nlay)


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Initialisations :
!------------------

mask(:,:)=.true.
RLvCp = RLVTT/RCPD

!
! calcul des caracteristiques de l environnement
       DO  ll=1,nlay
         DO ig=1,ngrid
            zo(ig,ll)=po(ig,ll)
            zl(ig,ll)=0.
            zh(ig,ll)=pt(ig,ll)
         EndDO
       EndDO
!
!
! Condensation :
!---------------
! Calcul de l'humidite a saturation et de la condensation

call thermcell_qsat(ngrid*nlay,mask,pplev,pt,po,pqsat)
DO ll=1,nlay
   DO ig=1,ngrid
      zl(ig,ll) = max(0.,po(ig,ll)-pqsat(ig,ll))
      zh(ig,ll) = pt(ig,ll)+RLvCp*zl(ig,ll)         !   T = Tl + Lv/Cp ql
      zo(ig,ll) = po(ig,ll)-zl(ig,ll)
   ENDDO
ENDDO
!
!
!-----------------------------------------------------------------------

      if (prt_level.ge.1) print*,'0 OK convect8'

      DO ll=1,nlay
         DO ig=1,ngrid
             zpspsk(ig,ll)=(pplay(ig,ll)/100000.)**RKAPPA
             zu(ig,ll)=pu(ig,ll)
             zv(ig,ll)=pv(ig,ll)
!attention zh est maintenant le profil de T et plus le profil de theta !
! Quelle horreur ! A eviter.
!
!   T-> Theta
            ztv(ig,ll)=zh(ig,ll)/zpspsk(ig,ll)
!Theta_v
            ztv(ig,ll)=ztv(ig,ll)*(1.+RETV*(zo(ig,ll))-zl(ig,ll))
!Thetal
            zthl(ig,ll)=pt(ig,ll)/zpspsk(ig,ll)
!            
         ENDDO
      ENDDO
 
      RETURN
      END
