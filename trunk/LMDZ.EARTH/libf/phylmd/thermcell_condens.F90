subroutine thermcell_condens(klon,active,zpspsk,pplev,ztla,zqta,zqla)
implicit none

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"


!====================================================================
! DECLARATIONS
!====================================================================

! Arguments
INTEGER klon
REAL zpspsk(klon),pplev(klon)
REAL ztla(klon),zqta(klon),zqla(klon)
LOGICAL active(klon)

! Variables locales
INTEGER ig,iter
REAL Tbef(klon),DT(klon)
REAL tdelta,qsatbef,zcor,qlbef,zdelta,zcvm5,dqsat,num,denom,dqsat_dT
logical Zsat
REAL RLvCp
REAL, SAVE :: DDT0=.01
LOGICAL afaire(klon),tout_converge

!====================================================================
! INITIALISATIONS
!====================================================================

RLvCp = RLVTT/RCPD
tout_converge=.false.
afaire(:)=.false.
DT(:)=0.


!====================================================================
! Routine a vectoriser en copiant active dans converge et en mettant
! la boucle sur les iterations a l'exterieur est en mettant
! converge= false des que la convergence est atteinte.
!====================================================================

do ig=1,klon
   if (active(ig)) then
               Tbef(ig)=ztla(ig)*zpspsk(ig)
               zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
               qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
               qsatbef=MIN(0.5,qsatbef)
               zcor=1./(1.-retv*qsatbef)
               qsatbef=qsatbef*zcor
               qlbef=max(0.,zqta(ig)-qsatbef)
               DT(ig) = 0.5*RLvCp*qlbef
     endif
enddo

do iter=1,10
    afaire(:)=abs(DT(:)).gt.DDT0
    do ig=1,klon
               if (afaire(ig)) then
                 Tbef(ig)=Tbef(ig)+DT(ig)
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 qsatbef= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig)
                 qsatbef=MIN(0.5,qsatbef)
                 zcor=1./(1.-retv*qsatbef)
                 qsatbef=qsatbef*zcor
                 qlbef=zqta(ig)-qsatbef
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
                 zcor=1./(1.-retv*qsatbef)
                 dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef,zcor)
                 num=-Tbef(ig)+ztla(ig)*zpspsk(ig)+RLvCp*qlbef
                 denom=1.+RLvCp*dqsat_dT
                 zqla(ig) = max(0.,zqta(ig)-qsatbef) 
                 DT(ig)=num/denom
               endif
    enddo
enddo

return
end
