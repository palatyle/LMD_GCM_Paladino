       SUBROUTINE cloudth(ngrid,klev,ind2,  &
     &           ztv,po,zqta,fraca, & 
     &           qcloud,ctot,zpspsk,paprs,ztla,zthl, &
     &           ratqs,zqs,t)


      IMPLICIT NONE


!===========================================================================
! Auteur : Arnaud Octavio Jam (LMD/CNRS)
! Date : 25 Mai 2010
! Objet : calcule les valeurs de qc et rneb dans les thermiques
!===========================================================================


#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"
#include "thermcell.h"

      INTEGER itap,ind1,ind2
      INTEGER ngrid,klev,klon,l,ig 
      
      REAL ztv(ngrid,klev)
      REAL po(ngrid)
      REAL zqenv(ngrid)   
      REAL zqta(ngrid,klev)
          
      REAL fraca(ngrid,klev+1)
      REAL zpspsk(ngrid,klev)
      REAL paprs(ngrid,klev+1)
      REAL ztla(ngrid,klev)
      REAL zthl(ngrid,klev)

      REAL zqsatth(ngrid,klev)
      REAL zqsatenv(ngrid,klev)
      
      
      REAL sigma1(ngrid,klev)                                                         
      REAL sigma2(ngrid,klev)
      REAL qlth(ngrid,klev)
      REAL qlenv(ngrid,klev)
      REAL qltot(ngrid,klev) 
      REAL cth(ngrid,klev)  
      REAL cenv(ngrid,klev)   
      REAL ctot(ngrid,klev)
      REAL rneb(ngrid,klev)
      REAL t(ngrid,klev)                                                                  
      REAL qsatmmussig1,qsatmmussig2,sqrt2pi,pi
      REAL rdd,cppd,Lv
      REAL alth,alenv,ath,aenv
      REAL sth,senv,sigma1s,sigma2s,xth,xenv
      REAL Tbef,zdelta,qsatbef,zcor
      REAL alpha,qlbef  
      REAL ratqs(ngrid,klev) ! determine la largeur de distribution de vapeur
      
      REAL zpdf_sig(ngrid),zpdf_k(ngrid),zpdf_delta(ngrid)
      REAL zpdf_a(ngrid),zpdf_b(ngrid),zpdf_e1(ngrid),zpdf_e2(ngrid) 
      REAL zqs(ngrid), qcloud(ngrid)
      REAL erf





!      print*,ngrid,klev,ind1,ind2,ztv(ind1,ind2),po(ind1),zqta(ind1,ind2), &
!     &       fraca(ind1,ind2),zpspsk(ind1,ind2),paprs(ind1,ind2),ztla(ind1,ind2),zthl(ind1,ind2), &
!     &       'verif'


!      LOGICAL active(ngrid)   
      
!-----------------------------------------------------------------------------------------------------------------
! Initialisation des variables réelles
!-----------------------------------------------------------------------------------------------------------------
      sigma1(:,:)=0.
      sigma2(:,:)=0.
      qlth(:,:)=0.
      qlenv(:,:)=0.  
      qltot(:,:)=0.
      rneb(:,:)=0.
      qcloud(:)=0.
      cth(:,:)=0.
      cenv(:,:)=0.
      ctot(:,:)=0.
      qsatmmussig1=0.
      qsatmmussig2=0.
      rdd=287.04
      cppd=1005.7
      pi=3.14159 
      Lv=2.5e6
      sqrt2pi=sqrt(2.*pi)



!------------------------------------------------------------------------------------------------------------------
! Calcul de la fraction du thermique et des écart-types des distributions
!------------------------------------------------------------------------------------------------------------------                 
      do ind1=1,ngrid

      if ((ztv(ind1,1).gt.ztv(ind1,2)).and.(fraca(ind1,ind2).gt.1.e-10)) then 

      zqenv(ind1)=(po(ind1)-fraca(ind1,ind2)*zqta(ind1,ind2))/(1.-fraca(ind1,ind2))


!      zqenv(ind1)=po(ind1)
      Tbef=zthl(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef




      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 




      Tbef=ztla(ind1,ind2)*zpspsk(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatth(ind1,ind2)=qsatbef
            
      alth=(0.622*Lv*zqsatth(ind1,ind2))/(rdd*ztla(ind1,ind2)**2)   
      ath=1./(1.+(alth*Lv/cppd))
      sth=ath*(zqta(ind1,ind2)-zqsatth(ind1,ind2)) 
      
      

!-----------------------------------------------------------------------------------------------------------------
! Calcul des écart-types pour s
!-----------------------------------------------------------------------------------------------------------------

      sigma1s=(1.1**0.5)*(fraca(ind1,ind2)**0.6)/(1-fraca(ind1,ind2))*((sth-senv)**2)**0.5+ratqs(ind1,ind2)*po(ind1)
      sigma2s=0.11*((sth-senv)**2)**0.5/(fraca(ind1,ind2)+0.02)**0.4+0.002*zqta(ind1,ind2)  

 
!-----------------------------------------------------------------------------------------------------------------
! Calcul de l'eau condensée et de la couverture nuageuse
!-----------------------------------------------------------------------------------------------------------------
      sqrt2pi=sqrt(2.*pi)
      xth=sth/(sqrt(2.)*sigma2s)
      xenv=senv/(sqrt(2.)*sigma1s)
      cth(ind1,ind2)=0.5*(1.+1.*erf(xth))
      cenv(ind1,ind2)=0.5*(1.+1.*erf(xenv)) 
      ctot(ind1,ind2)=fraca(ind1,ind2)*cth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*cenv(ind1,ind2)   
!      ctot(ind1,ind2)=alpha*cth(ind1,ind2)+(1.-1.*alpha)*cenv(ind1,ind2) 



      qlth(ind1,ind2)=sigma2s*((exp(-1.*xth**2)/sqrt2pi)+xth*sqrt(2.)*cth(ind1,ind2))
      qlenv(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))   
      qltot(ind1,ind2)=fraca(ind1,ind2)*qlth(ind1,ind2)+(1.-1.*fraca(ind1,ind2))*qlenv(ind1,ind2)
!      qltot(ind1,ind2)=alpha*qlth(ind1,ind2)+(1.-1.*alpha)*qlenv(ind1,ind2)
     

!      print*,senv,sth,sigma1s,sigma2s,fraca(ind1,ind2),'senv et sth et sig1 et sig2 et alpha'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ctot(ind1,ind2).lt.1.e-10) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqs(ind1)

      endif                           
      
          
!     print*,sth,sigma2s,qlth(ind1,ind2),ctot(ind1,ind2),qltot(ind1,ind2),'verif'


      else  ! gaussienne environnement seule
      
      zqenv(ind1)=po(ind1)
      Tbef=t(ind1,ind2)
      zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
      qsatbef= R2ES * FOEEW(Tbef,zdelta)/paprs(ind1,ind2)
      qsatbef=MIN(0.5,qsatbef)
      zcor=1./(1.-retv*qsatbef)
      qsatbef=qsatbef*zcor
      zqsatenv(ind1,ind2)=qsatbef
      

!      qlbef=Max(po(ind1)-zqsatenv(ind1,ind2),0.)
      zthl(ind1,ind2)=t(ind1,ind2)*(101325/paprs(ind1,ind2))**(rdd/cppd)
      alenv=(0.622*Lv*zqsatenv(ind1,ind2))/(rdd*zthl(ind1,ind2)**2)  
      aenv=1./(1.+(alenv*Lv/cppd))
      senv=aenv*(po(ind1)-zqsatenv(ind1,ind2)) 
      

      sigma1s=ratqs(ind1,ind2)*zqenv(ind1)

      sqrt2pi=sqrt(2.*pi)
      xenv=senv/(sqrt(2.)*sigma1s)
      ctot(ind1,ind2)=0.5*(1.+1.*erf(xenv))
      qltot(ind1,ind2)=sigma1s*((exp(-1.*xenv**2)/sqrt2pi)+xenv*sqrt(2.)*cenv(ind1,ind2))
      
      if (ctot(ind1,ind2).lt.1.e-3) then
      ctot(ind1,ind2)=0.
      qcloud(ind1)=zqsatenv(ind1,ind2) 

      else   
                
      ctot(ind1,ind2)=ctot(ind1,ind2) 
      qcloud(ind1)=qltot(ind1,ind2)/ctot(ind1,ind2)+zqsatenv(ind1,ind2)

      endif    
 
 
 
 
 
 
      endif   
      enddo
     
      return
      end





                                                                            





