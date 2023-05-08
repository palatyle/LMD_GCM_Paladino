      SUBROUTINE cv3_vertmix(len,nd,iflag,plim1,plim2,p,ph,t,q,u,v
     :                     ,w,wi,nk,tmix,thmix,qmix,qsmix
     :                     ,umix,vmix,plcl)
***************************************************************
*                                                             *
* CV3_VERTMIX   Brassage adiabatique d'une couche d'epaisseur *
*               arbitraire.                                   *
*                                                             *
* written by   : Grandpeix Jean-Yves, 28/12/2001, 13.14.24    *
* modified by :  Filiberti M-A 06/2005 vectorisation          *
***************************************************************
*
       implicit none
C==============================================================
C
C vertmix : determine theta et r du melange obtenu en brassant
C adiabatiquement entre plim1 et plim2, avec une ponderation w.
C
C===============================================================

#include "cvthermo.h"
#include "YOETHF.h"
#include "YOMCST.h"
#include "FCTTRE.h"
c input :
      integer nd,len
      integer nk(len),iflag(len)
      real t(len,nd),q(len,nd),w(nd)
      real u(len,nd),v(len,nd)
      real p(len,nd),ph(len,nd+1)
      real plim1(len),plim2(len)
c output :
      real tmix(len),thmix(len),qmix(len),wi(len,nd)
      real umix(len),vmix(len)
      real qsmix(len)
      real plcl(len)
c internal variables :
      integer j1(len),j2(len),niflag7
      real A,B
      real ahm(len),dpw(len),coef(len)
      real p1(len,nd),p2(len,nd)
      real rdcp(len),a2(len),b2(len),pnk(len)
      real rh(len),chi(len)
      real cpn
      real x,y,p0,p0m1,zdelta,zcor

      integer i,j

      do j = 1,nd
        do i=1,len
          if (plim1(i).le.ph(i,j)) j1(i) = j
          if (plim2(i).ge.ph(i,j+1).and.plim2(i).lt.ph(i,j)) j2(i) = j
        enddo
      enddo
c
      do j=1,nd
        do i = 1,len
          wi(i,j) = 0.
        enddo
      enddo
      do i = 1,len
       ahm(i)=0.
       qmix(i)=0.
       umix(i)=0.
       vmix(i)=0.
       dpw(i) =0.
       a2(i)=0.0
       b2(i) = 0.
       pnk(i) = p(i,nk(i))
      enddo
c
      p0 = 1000.
      p0m1 = 1./p0
c
      do i=1,len
        coef(i) = 1./(plim1(i)-plim2(i))
      end do
c
      do  j=1,nd
        do i=1,len
          if (j.ge.j1(i).and.j.le.j2(i)) then
            p1(i,j) = min(ph(i,j),plim1(i))
            p2(i,j) = max(ph(i,j+1),plim2(i))
cCRtest:couplage thermiques: deja normalise
c             wi(i,j) = w(j)
c             print*,'wi',wi(i,j)
            wi(i,j) = w(j)*(p1(i,j)-p2(i,j))*coef(i)
            dpw(i) = dpw(i)+wi(i,j)
          endif
        end do
      end do
cCR:print
c      do i=1,len
c         print*,'plim',plim1(i),plim2(i)
c      enddo 
      do  j=1,nd
        do i=1,len
          if (j.ge.j1(i).and.j.le.j2(i)) then
            wi(i,j)=wi(i,j)/dpw(i)
            ahm(i)=ahm(i)+(cpd*(1.-q(i,j))+q(i,j)*cpv)*t(i,j)*wi(i,j)
            qmix(i)=qmix(i)+q(i,j)*wi(i,j)
            umix(i)=umix(i)+u(i,j)*wi(i,j)
            vmix(i)=vmix(i)+v(i,j)*wi(i,j)
          endif
        end do
      end do
c
      do i=1,len
         rdcp(i)=(rrd*(1.-qmix(i))+qmix(i)*rrv)/
     :            (cpd*(1.-qmix(i))+qmix(i)*cpv)
      end do
c

c
      do 20 j=1,nd
        do 18 i=1,len
          if (j.ge.j1(i).and.j.le.j2(i)) then
cc            x=(.5*(p1(i,j)+p2(i,j))*p0m1)**rdcp(i)
            y=(.5*(p1(i,j)+p2(i,j))/pnk(i))**rdcp(i)
cc            a2(i)=a2(i)+(cpd*(1.-qmix(i))+qmix(i)*cpv)*x*wi(i,j)
            b2(i)=b2(i)+(cpd*(1.-qmix(i))+qmix(i)*cpv)*y*wi(i,j)
          endif
   18  continue
   20  continue
c
       do i=1,len
         tmix(i) = ahm(i)/b2(i)
         thmix(i) =tmix(i)*(p0/pnk(i))**rdcp(i)
c         print*,'thmix ahm',ahm(i),b2(i)
c         print*,'thmix t',tmix(i),p0 
c         print*,'thmix p',pnk(i),rdcp(i)
c         print*,'thmix',thmix(i)
cc         thmix(i) = ahm(i)/a2(i)
cc         tmix(i)= thmix(i)*(pnk(i)*p0m1)**rdcp(i)
         zdelta=max(0.,sign(1.,rtt-tmix(i)))
         qsmix(i)= r2es*FOEEW(tmix(i),zdelta)/(pnk(i)*100.)
         qsmix(i)=min(0.5,qsmix(i))
         zcor=1./(1.-retv*qsmix(i))
         qsmix(i)=qsmix(i)*zcor
       end do
c
!-------------------------------------------------------------------
! --- Calculate lifted condensation level of air at parcel origin level
! --- (Within 0.2% of formula of Bolton, MON. WEA. REV.,1980)
!-------------------------------------------------------------------

       A = 1669.0 ! convect3
       B = 122.0  ! convect3


       niflag7=0
       do 260 i=1,len

        if (iflag(i).ne.7) then ! modif sb Jun7th 2002
c
        rh(i)=qmix(i)/qsmix(i)
        chi(i)=tmix(i)/(A-B*rh(i)-tmix(i)) ! convect3
c   ATTENTION, la LIGNE DESSOUS A ETE RAJOUTEE ARBITRAIREMENT ET
c  MASQUE UN PB POTENTIEL
        chi(i)=max(chi(i),0.)
        rh(i)=max(rh(i),0.)
        plcl(i)=pnk(i)*(rh(i)**chi(i))
        if(((plcl(i).lt.200.0).or.(plcl(i).ge.2000.0))
     &   .and.(iflag(i).eq.0))iflag(i)=8

        else

          niflag7=niflag7+1
          plcl(i)=plim2(i)
c
        endif ! iflag=7

c      print*,'NIFLAG7  =',niflag7

 260   continue

      return
      end

