      SUBROUTINE thermcell_2002(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      USE dimphy
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"

c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/
c$OMP THREADPRIVATE(idetr)

c   local:
c   ------

      INTEGER ig,k,l,lmax(klon,klev+1),lmaxa(klon),lmix(klon)
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev)
!      common/comtherm/thetath2,wth2

      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit
c$OMP THREADPRIVATE(isplit)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

      real wmax(klon,klev),wmaxa(klon)

      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

      character (len=2) :: str2
      character (len=10) :: str10

      character (len=20) :: modname='thermcell2002'
      character (len=80) :: abort_message

      LOGICAL vtest(klon),down

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
c$OMP THREADPRIVATE(ncorrec)

c
c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

      print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
            zo(ig,l)=po(ig,l)
            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
1015     CONTINUE
1010  CONTINUE

c     print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo

c     print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo

c     print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------


      do k=1,nlay-1
         do ig=1,ngrid
            wa(ig,k,k)=0.
            wa(ig,k,k+1)=2.*RG*(ztv(ig,k)-ztv(ig,k+1))/ztv(ig,k+1)
     s      *(zlev(ig,k+1)-zlev(ig,k))
         enddo
         do l=k+1,nlay-1
            do ig=1,ngrid
               wa(ig,k,l+1)=wa(ig,k,l)+
     s         2.*RG*(ztv(ig,k)-ztv(ig,l))/ztv(ig,l)
     s         *(zlev(ig,l+1)-zlev(ig,l))
            enddo
         enddo
         do ig=1,ngrid
            wa(ig,k,nlay+1)=0.
         enddo
      enddo

c     print*,'4 OK convect8'
c Calcul de la couche correspondant a la hauteur du thermique
      do k=1,nlay-1
         do ig=1,ngrid
            lmax(ig,k)=k
         enddo
         do l=nlay,k+1,-1
            do ig=1,ngrid
               if(wa(ig,k,l).le.1.e-10) lmax(ig,k)=l-1
            enddo
         enddo
      enddo

c     print*,'5 OK convect8'
c   Calcule du w max du thermique
      do k=1,nlay
      do ig=1,ngrid
         wmax(ig,k)=0.
      enddo
      enddo

      do k=1,nlay-1
         do l=k,nlay
            do ig=1,ngrid
               if (l.le.lmax(ig,k)) then
                  wa(ig,k,l)=sqrt(wa(ig,k,l))
                  wmax(ig,k)=max(wmax(ig,k),wa(ig,k,l))
               else
                  wa(ig,k,l)=0.
               endif
            enddo
         enddo
      enddo

      do k=1,nlay-1
         do ig=1,ngrid
             pu_therm(ig,k)=sqrt(wmax(ig,k))
             pv_therm(ig,k)=sqrt(wmax(ig,k))
         enddo
      enddo

c     print*,'6 OK convect8'
c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=500.
      enddo
c     print*,'LMAX LMAX LMAX '
      do k=1,nlay-1
         do  ig=1,ngrid
            zmax(ig)=max(zmax(ig),zlev(ig,lmax(ig,k))-zlev(ig,k))
         enddo
c     print*,k,lmax(1,k)
      enddo
c     print*,'ZMAX ZMAX ZMAX ',zmax
c      call dump2d(iim,jjm-1,zmax(2:ngrid-1),'ZMAX      ')

c   Calcul de l'entrainement.
c   Le rapport d'aspect relie la largeur de l'ascendance a l'epaisseur
c   de la couche d'alimentation en partant du principe que la vitesse
c   maximum dans l'ascendance est la vitesse d'entrainement horizontale.
      do k=1,nlay
         do ig=1,ngrid
            zzz=rho(ig,k)*wmax(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
     s    /(zmax(ig)*r_aspect)
            if(w2di.eq.2) then
               entr(ig,k)=entr(ig,k)+
     s         ptimestep*(zzz-entr(ig,k))/REAL(tho)
            else
               entr(ig,k)=zzz
            endif
            ztva(ig,k)=ztv(ig,k)
         enddo
      enddo

c     print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

c     print*,'8 OK convect8'
      do ig=1,ngrid
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo


      do l=1,nlay-2
         do ig=1,ngrid
c           if (zw2(ig,l).lt.1.e-10.and.ztv(ig,l).gt.ztv(ig,l+1)) then
c         print*,'COUCOU ',l,zw2(ig,l),ztv(ig,l),ztv(ig,l+1)
            if (zw2(ig,l).lt.1.e-10.and.ztv(ig,l).gt.ztv(ig,l+1)
     s       .and.entr(ig,l).gt.1.e-10) then
c        print*,'COUCOU cas 1'
c   Initialisation de l'ascendance
c              lmix(ig)=1
               ztva(ig,l)=ztv(ig,l)
               fmc(ig,l)=0.
               fmc(ig,l+1)=entr(ig,l)
               zw2(ig,l)=0.
c     if (.not.ztv(ig,l+1).gt.150.) then
c     print*,'ig,l+1,ztv(ig,l+1)'
c     print*, ig,l+1,ztv(ig,l+1)
c     endif
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s         *(zlev(ig,l+1)-zlev(ig,l))
               larg_detr(ig,l)=0.
            else if (zw2(ig,l).ge.1.e-10.and.
     .               fmc(ig,l)+entr(ig,l).gt.1.e-10) then
c   Incrementation...
               fmc(ig,l+1)=fmc(ig,l)+entr(ig,l)
c     if (.not.fmc(ig,l+1).gt.1.e-15) then
c     print*,'ig,l+1,fmc(ig,l+1)'
c     print*, ig,l+1,fmc(ig,l+1)
c     print*,'Fmc ',(fmc(ig,ll),ll=1,klev+1)
c     print*,'W2 ',(zw2(ig,ll),ll=1,klev+1)
c     print*,'Tv ',(ztv(ig,ll),ll=1,klev)
c     print*,'Entr ',(entr(ig,ll),ll=1,klev)
c     endif
               ztva(ig,l)=(fmc(ig,l)*ztva(ig,l-1)+entr(ig,l)*ztv(ig,l))
     s          /fmc(ig,l+1)
c  mise a jour de la vitesse ascendante (l'air entraine de la couche
c  consideree commence avec une vitesse nulle).
               zw2(ig,l+1)=zw2(ig,l)*(fmc(ig,l)/fmc(ig,l+1))**2+
     s         2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s         *(zlev(ig,l+1)-zlev(ig,l))
            endif
            if (zw2(ig,l+1).lt.0.) then
               zw2(ig,l+1)=0.
               lmaxa(ig)=l
            else
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
c        print*,'COUCOU cas 2 LMIX=',lmix(ig),wa_moy(ig,l+1),wmaxa(ig)
         enddo
      enddo

c     print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
            endif
         enddo
       enddo

c     print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   calcul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.

      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1.) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
               if(l.gt.lmix(ig)) then
                  xxx(ig,l)=(lmaxa(ig)+1.-l) / (lmaxa(ig)+1.-lmix(ig))
           if (idetr.eq.0) then
               fraca(ig,l)=fraca(ig,lmix(ig))
           else if (idetr.eq.1) then
               fraca(ig,l)=fraca(ig,lmix(ig))*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fraca(ig,lmix(ig))*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fraca(ig,lmix(ig))*xxx(ig,l)**2
           endif
               endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo

c     print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1) then
              abort_message = 'fracd trop petit'
              CALL abort_gcm (modname,abort_message,1) 
           else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

c     print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
c     print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
c    s         ,fm(ig,l+1)*ptimestep
c    s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(entr(ig,l)*ptimestep.gt.masse(ig,l)) then
c     print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
c    s         ,entr(ig,l)*ptimestep
c    s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
c     print*,'WARN!!! masse exagere ig=',ig,'   l=',l
c    s         ,'   M=',masse(ig,l)
c     print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
c    s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
c     print*,'zlev(ig,l+1),zlev(ig,l)'
c    s                ,zlev(ig,l+1),zlev(ig,l)
c     print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
c    s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
            endif
            if(.not.entr(ig,l).ge.0..or..not.entr(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/REAL(tho)
         entr0=entr0+ptimestep*(entr-entr0)/REAL(tho)
      else
         fm0=fm
         entr0=entr
      endif

      if (1.eq.1) then
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zh,zdhadj,zha)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zo,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

      do l=1,nlay
         do ig=1,ngrid
            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
         enddo
      enddo



c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

c     print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------

      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
                entr(ig,l)=entr(ig,l)-detr(ig,l)
                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
         enddo
      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1


c #define und
	goto 123
#ifdef und
      CALL writeg1d(1,nlay,wd,'wd      ','wd      ')
      CALL writeg1d(1,nlay,zwa,'wa      ','wa      ')
      CALL writeg1d(1,nlay,fracd,'fracd      ','fracd      ')
      CALL writeg1d(1,nlay,fraca,'fraca      ','fraca      ')
      CALL writeg1d(1,nlay,wa_moy,'wam         ','wam         ')
      CALL writeg1d(1,nlay,zla,'la      ','la      ')
      CALL writeg1d(1,nlay,zld,'ld      ','ld      ')
      CALL writeg1d(1,nlay,pt,'pt      ','pt      ')
      CALL writeg1d(1,nlay,zh,'zh      ','zh      ')
      CALL writeg1d(1,nlay,zha,'zha      ','zha      ')
      CALL writeg1d(1,nlay,zu,'zu      ','zu      ')
      CALL writeg1d(1,nlay,zv,'zv      ','zv      ')
      CALL writeg1d(1,nlay,zo,'zo      ','zo      ')
      CALL writeg1d(1,nlay,wh,'wh      ','wh      ')
      CALL writeg1d(1,nlay,wu,'wu      ','wu      ')
      CALL writeg1d(1,nlay,wv,'wv      ','wv      ')
      CALL writeg1d(1,nlay,wo,'w15uo     ','wXo     ')
      CALL writeg1d(1,nlay,zdhadj,'zdhadj      ','zdhadj      ')
      CALL writeg1d(1,nlay,pduadj,'pduadj      ','pduadj      ')
      CALL writeg1d(1,nlay,pdvadj,'pdvadj      ','pdvadj      ')
      CALL writeg1d(1,nlay,pdoadj,'pdoadj      ','pdoadj      ')
      CALL writeg1d(1,nlay,entr  ,'entr        ','entr        ')
      CALL writeg1d(1,nlay,detr  ,'detr        ','detr        ')
      CALL writeg1d(1,nlay,fm    ,'fm          ','fm          ')

      CALL writeg1d(1,nlay,pdtadj,'pdtadj    ','pdtadj    ')
      CALL writeg1d(1,nlay,pplay,'pplay     ','pplay     ')
      CALL writeg1d(1,nlay,pplev,'pplev     ','pplev     ')
c   recalcul des flux en diagnostique...
c     print*,'PAS DE TEMPS ',ptimestep
       call dt2F(pplev,pplay,pt,pdtadj,wh)
      CALL writeg1d(1,nlay,wh,'wh2     ','wh2     ')
#endif
123   continue
! #define troisD
#ifdef troisD
c       if (sorties) then
      print*,'Debut des wrgradsfi'

c      print*,'16 OK convect8'
         call wrgradsfi(1,nlay,wd,'wd        ','wd        ')
         call wrgradsfi(1,nlay,zwa,'wa        ','wa        ')
         call wrgradsfi(1,nlay,fracd,'fracd     ','fracd     ')
         call wrgradsfi(1,nlay,fraca,'fraca     ','fraca     ')
         call wrgradsfi(1,nlay,xxx,'xxx       ','xxx       ')
         call wrgradsfi(1,nlay,wa_moy,'wam       ','wam       ')
c      print*,'WA6 ',wa_moy
         call wrgradsfi(1,nlay,zla,'la        ','la        ')
         call wrgradsfi(1,nlay,zld,'ld        ','ld        ')
         call wrgradsfi(1,nlay,pt,'pt        ','pt        ')
         call wrgradsfi(1,nlay,zh,'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha,'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua,'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva,'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu,'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv,'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,wh,'wh        ','wh        ')
         call wrgradsfi(1,nlay,wu,'wu        ','wu        ')
         call wrgradsfi(1,nlay,wv,'wv        ','wv        ')
         call wrgradsfi(1,nlay,wo,'wo        ','wo        ')
         call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
         call wrgradsfi(1,nlay,zdhadj,'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj,'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj,'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj,'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr,'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr,'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm,'fm        ','fm        ')
         call wrgradsfi(1,nlay,fmc,'fmc       ','fmc       ')
         call wrgradsfi(1,nlay,zw2,'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,ztva,'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv,'ztv       ','ztv       ')

         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,larg_cons,'Lc        ','Lc        ')
         call wrgradsfi(1,nlay,larg_detr,'Ldetr     ','Ldetr     ')


c      print*,'17 OK convect8'

         do k=1,klev/10
            write(str2,'(i2.2)') k
            str10='wa'//str2
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=wa(ig,k,l)
               enddo
            enddo   
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=larg_part(ig,k,l)
               enddo
            enddo
            str10='la'//str2
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
         enddo


c     print*,'18 OK convect8'
c      endif
      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

c     print*,'19 OK convect8'
      return
      end

      SUBROUTINE thermcell_cld(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi,zlev,debut
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0,zqla,lmax
     s                  ,zmax_sec,wmax_sec,zw_sec,lmix_sec
     s                  ,ratqscth,ratqsdiff
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      USE dimphy
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"

c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/
c$OMP THREADPRIVATE(idetr)

c   local:
c   ------

      INTEGER ig,k,l,lmaxa(klon),lmix(klon)
      real zsortie1d(klon)
c CR: on remplace lmax(klon,klev+1)
      INTEGER lmax(klon),lmin(klon),lentr(klon)
      real linter(klon)
      real zmix(klon), fracazmix(klon)
      real alpha
      save alpha
      data alpha/1./
c$OMP THREADPRIVATE(alpha)

c RC 
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz
      real zmax_sec(klon)
      real zmax_sec2(klon)
      real zw_sec(klon,klev+1)
      INTEGER lmix_sec(klon)
      real w_est(klon,klev+1)
con garde le zmax du pas de temps precedent
c      real zmax0(klon)
c      save zmax0
c      real zmix0(klon)
c      save zmix0 
      REAL, SAVE, ALLOCATABLE :: zmax0(:), zmix0(:)
c$OMP THREADPRIVATE(zmax0, zmix0)

      real zlev(klon,klev+1),zlay(klon,klev)
      real deltaz(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      real zthl(klon,klev),zdthladj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      real zl(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zta(klon,klev)
      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev),wth3(klon,klev)
      real q2(klon,klev)
      real dtheta(klon,klev)
!      common/comtherm/thetath2,wth2
    
      real ratqscth(klon,klev)
      real sum
      real sumdiff
      real ratqsdiff(klon,klev)
      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit
c$OMP THREADPRIVATE(isplit)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

c     real wmax(klon,klev),wmaxa(klon)
      real wmax(klon),wmaxa(klon)
      real wmax_sec(klon)
      real wmax_sec2(klon)
      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real massetot(klon,klev)
      real detr0(klon,klev)
      real alim0(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

      real zcor,zdelta,zcvm5,qlbef
      real Tbef(klon),qsatbef(klon)
      real dqsat_dT,DT,num,denom
      REAL REPS,RLvCp,DDT0
      real ztla(klon,klev),zqla(klon,klev),zqta(klon,klev)
cCR niveau de condensation
      real nivcon(klon)
      real zcon(klon)
      real zqsat(klon,klev)
      real zqsatth(klon,klev) 
      PARAMETER (DDT0=.01)


cCR:nouvelles variables
      real f_star(klon,klev+1),entr_star(klon,klev)
      real detr_star(klon,klev)
      real alim_star_tot(klon),alim_star2(klon)
      real entr_star_tot(klon)
      real detr_star_tot(klon)
      real alim_star(klon,klev)
      real alim(klon,klev)
      real nu(klon,klev)
      real nu_e(klon,klev)
      real nu_min
      real nu_max
      real nu_r
      real f(klon)
c      real f(klon), f0(klon)
c     save f0
      REAL,SAVE, ALLOCATABLE :: f0(:)
c$OMP THREADPRIVATE(f0)

      real f_old
      real zlevinter(klon)
      logical, save :: first = .true.
c$OMP THREADPRIVATE(first)
c      data first /.false./
c      save first
      logical nuage
c      save nuage
      logical boucle
      logical therm
      logical debut
      logical rale
      integer test(klon)
      integer signe_zw2
cRC

      character*2 str2
      character*10 str10

      character (len=20) :: modname='thermcell_cld'
      character (len=80) :: abort_message

      LOGICAL vtest(klon),down
      LOGICAL Zsat(klon)

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
c$OMP THREADPRIVATE(ncorrec)

c

c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
      if (first) then
        allocate(zmix0(klon))
        allocate(zmax0(klon))
        allocate(f0(klon))
        first=.false.
      endif

       sorties=.false.
c     print*,'NOUVEAU DETR PLUIE '
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c Initialisation
      RLvCp = RLVTT/RCPD
      REPS  = RD/RV
cinitialisations de zqsat
      DO ll=1,nlay
         DO ig=1,ngrid
            zqsat(ig,ll)=0.
            zqsatth(ig,ll)=0.
         ENDDO
      ENDDO
c
con met le first a true pour le premier passage de la journée
      do ig=1,klon
            test(ig)=0
      enddo
      if (debut) then
         do ig=1,klon
            test(ig)=1
            f0(ig)=0.
            zmax0(ig)=0.
         enddo
      endif
      do ig=1,klon      
         if ((.not.debut).and.(f0(ig).lt.1.e-10)) then
            test(ig)=1
         endif
      enddo 
c     do ig=1,klon
c        print*,'test(ig)',test(ig),zmax0(ig)
c     enddo
      nuage=.false.
c-----------------------------------------------------------------------
cAM Calcul de T,q,ql a partir de Tl et qT
c   ---------------------------------------------------
c
c Pr Tprec=Tl calcul de qsat 
c Si qsat>qT T=Tl, q=qT
c Sinon DDT=(-Tprec+Tl+RLVCP (qT-qsat(T')) / (1+RLVCP dqsat/dt) 
c On cherche DDT < DDT0
c
c defaut
       DO  ll=1,nlay
         DO ig=1,ngrid
            zo(ig,ll)=po(ig,ll)
            zl(ig,ll)=0.
            zh(ig,ll)=pt(ig,ll)
         EndDO
       EndDO
       do ig=1,ngrid
          Zsat(ig)=.false.
       enddo
c
c
       DO ll=1,nlay
c les points insatures sont definitifs
         DO ig=1,ngrid
            Tbef(ig)=pt(ig,ll)
            zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
            qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,ll)
            qsatbef(ig)=MIN(0.5,qsatbef(ig))
            zcor=1./(1.-retv*qsatbef(ig))
            qsatbef(ig)=qsatbef(ig)*zcor
            Zsat(ig) = (max(0.,po(ig,ll)-qsatbef(ig)) .gt. 1.e-10)
         EndDO

         DO ig=1,ngrid
           if (Zsat(ig).and.(1.eq.1)) then
            qlbef=max(0.,po(ig,ll)-qsatbef(ig))
c si sature: ql est surestime, d'ou la sous-relax
            DT = 0.5*RLvCp*qlbef
c            write(18,*),'DT0=',DT
c on pourra enchainer 2 ou 3 calculs sans Do while
            do while (abs(DT).gt.DDT0)
c il faut verifier si c,a conserve quand on repasse en insature ...
              Tbef(ig)=Tbef(ig)+DT
              zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
              qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,ll)
              qsatbef(ig)=MIN(0.5,qsatbef(ig))
              zcor=1./(1.-retv*qsatbef(ig))
              qsatbef(ig)=qsatbef(ig)*zcor
c on veut le signe de qlbef
              qlbef=po(ig,ll)-qsatbef(ig)
              zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
              zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
              zcor=1./(1.-retv*qsatbef(ig))
              dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef(ig),zcor)
              num=-Tbef(ig)+pt(ig,ll)+RLvCp*qlbef
              denom=1.+RLvCp*dqsat_dT
              if (denom.lt.1.e-10) then
                  print*,'pb denom'
              endif
              DT=num/denom
            enddo
c on ecrit de maniere conservative (sat ou non)
            zl(ig,ll) = max(0.,qlbef)
c          T = Tl +Lv/Cp ql
            zh(ig,ll) = pt(ig,ll)+RLvCp*zl(ig,ll)
            zo(ig,ll) = po(ig,ll)-zl(ig,ll)
           endif
con ecrit zqsat 
            zqsat(ig,ll)=qsatbef(ig)     
         EndDO
       EndDO
cAM fin
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

c     print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
             zpspsk(ig,l)=(pplay(ig,l)/100000.)**RKAPPA
c            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
c            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
c            zo(ig,l)=po(ig,l)
c            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
cAM attention zh est maintenant le profil de T et plus le profil de theta !
c
c   T-> Theta
            ztv(ig,l)=zh(ig,l)/zpspsk(ig,l)
cAM Theta_v
            ztv(ig,l)=ztv(ig,l)*(1.+RETV*(zo(ig,l))
     s           -zl(ig,l))
cAM Thetal
            zthl(ig,l)=pt(ig,l)/zpspsk(ig,l)
c            
1015     CONTINUE
1010  CONTINUE

c     print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo
ccalcul de deltaz
      do l=1,nlay
         do ig=1,ngrid
            deltaz(ig,l)=zlev(ig,l+1)-zlev(ig,l)
         enddo
      enddo

c     print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
c            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
             rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*ztv(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo
cCr:ajout:calcul de la masse
      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo
c     print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------

cCR: ponderation entrainement des couches instables
cdef des alim_star tels que alim=f*alim_star      
      do l=1,klev
         do ig=1,ngrid 
            alim_star(ig,l)=0.
            alim(ig,l)=0.
         enddo
      enddo
c determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lentr(ig)=1
      enddo

con ne considere que les premieres couches instables
      therm=.false.
      do k=nlay-2,1,-1
         do ig=1,ngrid
            if (ztv(ig,k).gt.ztv(ig,k+1).and.
     s          ztv(ig,k+1).le.ztv(ig,k+2)) then
               lentr(ig)=k+1
               therm=.true.
            endif
          enddo
      enddo
c
c determination du lmin: couche d ou provient le thermique
      do ig=1,ngrid
         lmin(ig)=1
      enddo
      do ig=1,ngrid
         do l=nlay,2,-1
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.lt.lentr(ig)) then
cdef possibles pour alim_star: zdthetadz, dthetadz, zdtheta
             alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)
c     s                       *(zlev(ig,l+1)-zlev(ig,l))
     s                       *sqrt(zlev(ig,l+1))
c             alim_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
c     s                         /zlev(ig,lentr(ig)+2)))**(3./2.) 
            endif
         enddo
      enddo
      
c pas de thermique si couche 1 stable
      do ig=1,ngrid
c         if (lmin(ig).gt.1) then
cCRnouveau test
        if (alim_star(ig,1).lt.1.e-10) then 
            do l=1,klev
                alim_star(ig,l)=0.
            enddo
         endif
      enddo 
c calcul de l entrainement total
      do ig=1,ngrid
         alim_star_tot(ig)=0.
         entr_star_tot(ig)=0.
         detr_star_tot(ig)=0.
      enddo
      do ig=1,ngrid
         do k=1,klev
            alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,k)
         enddo
      enddo
c
c Calcul entrainement normalise
      do ig=1,ngrid 
         if (alim_star_tot(ig).gt.1.e-10) then
c         do l=1,lentr(ig)
          do l=1,klev
cdef possibles pour entr_star: zdthetadz, dthetadz, zdtheta 
            alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
         enddo
         endif
      enddo
       
c     print*,'fin calcul alim_star'

cAM:initialisations
      do k=1,nlay
         do ig=1,ngrid
            ztva(ig,k)=ztv(ig,k)
            ztla(ig,k)=zthl(ig,k)
            zqla(ig,k)=0.
            zqta(ig,k)=po(ig,k)
            Zsat(ig) =.false.
         enddo
      enddo 
      do k=1,klev
        do ig=1,ngrid
           detr_star(ig,k)=0.
           entr_star(ig,k)=0.
           detr(ig,k)=0.
           entr(ig,k)=0.
        enddo
      enddo
c     print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
cCR
            f_star(ig,k)=0.
cRC
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

cn     print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo

      nu_min=l_mix
      nu_max=1000.
c      do ig=1,ngrid
c      nu_max=wmax_sec(ig)
c      enddo
      do ig=1,ngrid
         do k=1,klev
            nu(ig,k)=0.
            nu_e(ig,k)=0.
         enddo
      enddo
cCalcul de l'excès de température du à la diffusion turbulente
      do ig=1,ngrid
         do l=1,klev
            dtheta(ig,l)=0.
         enddo
       enddo
      do ig=1,ngrid
         do l=1,lentr(ig)-1
      dtheta(ig,l)=sqrt(10.*0.4*zlev(ig,l+1)**2*1.
     s          *((ztv(ig,l+1)-ztv(ig,l))/(zlev(ig,l+1)-zlev(ig,l)))**2)
         enddo
       enddo
c      do l=1,nlay-2
      do l=1,klev-1
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.alim_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
cAM
ctest:on rajoute un excès de T dans couche alim
c               ztla(ig,l)=zthl(ig,l)+dtheta(ig,l)
               ztla(ig,l)=zthl(ig,l) 
ctest: on rajoute un excès de q dans la couche alim
c               zqta(ig,l)=po(ig,l)+0.001
               zqta(ig,l)=po(ig,l)
               zqla(ig,l)=zl(ig,l)
cAM
               f_star(ig,l+1)=alim_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               w_est(ig,l+1)=zw2(ig,l+1)
               larg_detr(ig,l)=0.
c     print*,'coucou boucle 1'
            else if ((zw2(ig,l).ge.1e-10).and.
     s         (f_star(ig,l)+alim_star(ig,l)).gt.1.e-10) then
c     print*,'coucou boucle 2'
cestimation du detrainement a partir de la geometrie du pas precedent
      if ((test(ig).eq.1).or.((.not.debut).and.(f0(ig).lt.1.e-10))) then
                  detr_star(ig,l)=0.
                  entr_star(ig,l)=0.
c     print*,'coucou test(ig)',test(ig),f0(ig),zmax0(ig)
             else
c     print*,'coucou debut detr'
ctests sur la definition du detr
        if (zqla(ig,l-1).gt.1.e-10) then
           nuage=.true.
        endif 

             w_est(ig,l+1)=zw2(ig,l)*
     s                   ((f_star(ig,l))**2)
     s                   /(f_star(ig,l)+alim_star(ig,l))**2+
     s                   2.*RG*(ztva(ig,l-1)-ztv(ig,l))/ztv(ig,l)
     s                   *(zlev(ig,l+1)-zlev(ig,l))
             if (w_est(ig,l+1).lt.0.) then
                w_est(ig,l+1)=zw2(ig,l)
             endif
      if (l.gt.2) then
             if ((w_est(ig,l+1).gt.w_est(ig,l)).and.
     s           (zlev(ig,l+1).lt.zmax_sec(ig)).and.
     s            (zqla(ig,l-1).lt.1.e-10)) then 
             detr_star(ig,l)=MAX(0.,(rhobarz(ig,l+1)
     s                *sqrt(w_est(ig,l+1))*sqrt(nu(ig,l)*zlev(ig,l+1))
     s       -rhobarz(ig,l)*sqrt(w_est(ig,l))*sqrt(nu(ig,l)*zlev(ig,l)))
     s       /(r_aspect*zmax_sec(ig)))
             else if ((zlev(ig,l+1).lt.zmax_sec(ig)).and.
     s                (zqla(ig,l-1).lt.1.e-10)) then
       detr_star(ig,l)=-f0(ig)*f_star(ig,lmix(ig))
     s /(rhobarz(ig,lmix(ig))*wmaxa(ig))*
     s (rhobarz(ig,l+1)*sqrt(w_est(ig,l+1))
     s *((zmax_sec(ig)-zlev(ig,l+1))/((zmax_sec(ig)-zlev(ig,lmix(ig)))))
     s **2.
     s -rhobarz(ig,l)*sqrt(w_est(ig,l))
     s *((zmax_sec(ig)-zlev(ig,l))/((zmax_sec(ig)-zlev(ig,lmix(ig)))))
     s **2.)
             else
       detr_star(ig,l)=0.002*f0(ig)*f_star(ig,l)
     s                *(zlev(ig,l+1)-zlev(ig,l))
             
             endif
        else
        detr_star(ig,l)=0.
        endif
       
         detr_star(ig,l)=detr_star(ig,l)/f0(ig)
         if (nuage) then
         entr_star(ig,l)=0.4*detr_star(ig,l)
         else
         entr_star(ig,l)=0.4*detr_star(ig,l)
         endif

             if ((detr_star(ig,l)).gt.f_star(ig,l)) then
              detr_star(ig,l)=f_star(ig,l)
c              entr_star(ig,l)=0.
              endif

             if ((l.lt.lentr(ig))) then
                 entr_star(ig,l)=0.
c                 detr_star(ig,l)=0.
             endif  

c           print*,'ok detr_star'
      endif
cprise en compte du detrainement dans le calcul du flux
             f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)
     s                      -detr_star(ig,l)
ctest
c             if (f_star(ig,l+1).lt.0.) then
c                f_star(ig,l+1)=0.
c                entr_star(ig,l)=0.
c                detr_star(ig,l)=f_star(ig,l)+alim_star(ig,l)
c             endif
ctest sur le signe de f_star
       if (f_star(ig,l+1).gt.1.e-10) then 
c                 then
ctest
c         if (((f_star(ig,l+1)+detr_star(ig,l)).gt.1.e-10)) then
cAM on melange Tl et qt du thermique
con rajoute un excès de T dans la couche alim
c               if (l.lt.lentr(ig)) then
c           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+
c     s     (alim_star(ig,l)+entr_star(ig,l))*(zthl(ig,l)+dtheta(ig,l)))
c     s     /(f_star(ig,l+1)+detr_star(ig,l))
c               else
               ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+
     s                    (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))
     s                 /(f_star(ig,l+1)+detr_star(ig,l))
c     s                    /(f_star(ig,l+1)) 
c               endif
con rajoute un excès de q dans la couche alim
c               if (l.lt.lentr(ig)) then
c               zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+
c     s           (alim_star(ig,l)+entr_star(ig,l))*(po(ig,l)+0.001))
c     s                 /(f_star(ig,l+1)+detr_star(ig,l))
c               else
               zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+
     s                    (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))
     s                 /(f_star(ig,l+1)+detr_star(ig,l))
c     s                   /(f_star(ig,l+1))
c               endif
cAM on en deduit thetav et ql du thermique
cCR test
c               Tbef(ig)=ztla(ig,l)*zpspsk(ig,l)
               Tbef(ig)=ztla(ig,l)*zpspsk(ig,l)
               zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
               qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,l)  
               qsatbef(ig)=MIN(0.5,qsatbef(ig))
               zcor=1./(1.-retv*qsatbef(ig))
               qsatbef(ig)=qsatbef(ig)*zcor
             Zsat(ig) = (max(0.,zqta(ig,l)-qsatbef(ig)) .gt. 1.e-10)

           if (Zsat(ig).and.(1.eq.1)) then
              qlbef=max(0.,zqta(ig,l)-qsatbef(ig))
              DT = 0.5*RLvCp*qlbef
c             write(17,*)'DT0=',DT
              do while (abs(DT).gt.DDT0)
c                 print*,'aie'
                 Tbef(ig)=Tbef(ig)+DT
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,l)
                 qsatbef(ig)=MIN(0.5,qsatbef(ig))
                 zcor=1./(1.-retv*qsatbef(ig))
                 qsatbef(ig)=qsatbef(ig)*zcor
                 qlbef=zqta(ig,l)-qsatbef(ig)

                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
                 zcor=1./(1.-retv*qsatbef(ig))
                 dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef(ig),zcor)
                 num=-Tbef(ig)+ztla(ig,l)*zpspsk(ig,l)+RLvCp*qlbef
                 denom=1.+RLvCp*dqsat_dT
                 if (denom.lt.1.e-10) then
                    print*,'pb denom'
                 endif
                 DT=num/denom
c                 write(17,*)'DT=',DT
              enddo
              zqla(ig,l) = max(0.,zqta(ig,l)-qsatbef(ig))
              zqla(ig,l) = max(0.,qlbef) 
c              zqla(ig,l)=0.
             endif
c             zqla(ig,l) = max(0.,zqta(ig,l)-qsatbef(ig))
c      
c on ecrit de maniere conservative (sat ou non)
c          T = Tl +Lv/Cp ql
cCR rq utilisation de humidite specifique ou rapport de melange?
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
con rajoute le calcul de zha pour diagnostiques (temp potentielle)
           zha(ig,l) = ztva(ig,l)
c           if (l.lt.lentr(ig)) then
c           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
c     s              -zqla(ig,l))-zqla(ig,l)) + 0.1
c           else
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
     s              -zqla(ig,l))-zqla(ig,l))
c           endif
c           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
c     s                 /(1.-retv*zqla(ig,l))
c           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
c           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
c     s                 /(1.-retv*zqta(ig,l))
c     s              -zqla(ig,l)/(1.-retv*zqla(ig,l)))
c     s              -zqla(ig,l)/(1.-retv*zqla(ig,l)))
c       write(13,*)zqla(ig,l),zqla(ig,l)/(1.-retv*zqla(ig,l))
con ecrit zqsat 
           zqsatth(ig,l)=qsatbef(ig)  
c        enddo
c        DO ig=1,ngrid
c           if (zw2(ig,l).ge.1.e-10.and.
c     s               f_star(ig,l)+entr_star(ig,l).gt.1.e-10) then
c  mise a jour de la vitesse ascendante (l'air entraine de la couche
c  consideree commence avec une vitesse nulle).
c
c            if (f_star(ig,l+1).gt.1.e-10) then
            zw2(ig,l+1)=zw2(ig,l)*
c     s                  ((f_star(ig,l)-detr_star(ig,l))**2)
c     s                  /f_star(ig,l+1)**2+
     s                   ((f_star(ig,l))**2)
     s                   /(f_star(ig,l+1)+detr_star(ig,l))**2+
c     s                    /(f_star(ig,l+1))**2+           
     s                   2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                   *(zlev(ig,l+1)-zlev(ig,l))
c     s                   *(f_star(ig,l)/f_star(ig,l+1))**2

            endif
        endif
c
            if (zw2(ig,l+1).lt.0.) then 
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
c              print*,'linter=',linter(ig)
c          else if ((zw2(ig,l+1).lt.1.e-10).and.(zw2(ig,l+1).ge.0.)) then
c               linter(ig)=l+1
c               print*,'linter=l',zw2(ig,l),zw2(ig,l+1)
            else
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
c            wa_moy(ig,l+1)=zw2(ig,l+1) 
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
      print*,'fin calcul zw2'
c
c Calcul de la couche correspondant a la hauteur du thermique
      do ig=1,ngrid
         lmax(ig)=lentr(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
            lmax(ig)=1
             lmin(ig)=1
             lentr(ig)=1
         endif
      enddo 
c    
c Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
                  print*,'pb2 zw2<0'
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
c calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*
     s    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)
     s    -zlev(ig,lmax(ig)))
cpour le cas ou on prend tjs lmin=1
c       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,1))
       zmax0(ig)=zmax(ig)
       write(11,*)'ig,lmax,linter',ig,lmax(ig),linter(ig)
       write(12,*)'ig,zlevinter,zmax',ig,zmax(ig),zlevinter(ig)
      enddo

cCalcul de zmax_sec et wmax_sec
      call fermeture_seche(ngrid,nlay
     s                  ,pplay,pplev,pphi,zlev,rhobarz,f0,zpspsk
     s                  ,alim,zh,zo,lentr,lmin,nu_min,nu_max,r_aspect
     s                  ,zmax_sec2,wmax_sec2)

      print*,'avant fermeture'
c Fermeture,determination de f
c en lmax f=d-e
      do ig=1,ngrid
c      entr_star(ig,lmax(ig))=0.
c      f_star(ig,lmax(ig)+1)=0.
c      detr_star(ig,lmax(ig))=f_star(ig,lmax(ig))+entr_star(ig,lmax(ig))
c     s                       +alim_star(ig,lmax(ig))
      enddo
c
      do ig=1,ngrid
         alim_star2(ig)=0.
      enddo
ccalcul de entr_star_tot
      do ig=1,ngrid
         do k=1,lmix(ig)
            entr_star_tot(ig)=entr_star_tot(ig)
c     s                        +entr_star(ig,k)
     s                        +alim_star(ig,k)
c     s                        -detr_star(ig,k)
            detr_star_tot(ig)=detr_star_tot(ig)
c     s                        +alim_star(ig,k)
     s                        -detr_star(ig,k)
     s                        +entr_star(ig,k)
         enddo
      enddo
      
      do ig=1,ngrid
         if (alim_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else   
c             do k=lmin(ig),lentr(ig)
             do k=1,lentr(ig)
                alim_star2(ig)=alim_star2(ig)+alim_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
             if ((zmax_sec(ig).gt.1.e-10).and.(1.eq.1)) then 
             f(ig)=wmax_sec(ig)/(max(500.,zmax_sec(ig))*r_aspect
     s             *alim_star2(ig))
             f(ig)=f(ig)+(f0(ig)-f(ig))*exp((-ptimestep/
     s                     zmax_sec(ig))*wmax_sec(ig))
             else
             f(ig)=wmax(ig)/(max(500.,zmax(ig))*r_aspect*alim_star2(ig))
            f(ig)=f(ig)+(f0(ig)-f(ig))*exp((-ptimestep/
     s                     zmax(ig))*wmax(ig))
             endif
         endif
         f0(ig)=f(ig)
      enddo
      print*,'apres fermeture'
c Calcul de l'entrainement
         do ig=1,ngrid 
            do k=1,klev
            alim(ig,k)=f(ig)*alim_star(ig,k)
         enddo
      enddo
cCR:test pour entrainer moins que la masse
c       do ig=1,ngrid
c          do l=1,lentr(ig)
c             if ((alim(ig,l)*ptimestep).gt.(0.9*masse(ig,l))) then
c                alim(ig,l+1)=alim(ig,l+1)+alim(ig,l)
c     s                       -0.9*masse(ig,l)/ptimestep
c                alim(ig,l)=0.9*masse(ig,l)/ptimestep
c             endif
c          enddo
c       enddo
c calcul du détrainement
         do ig=1,klon
             do k=1,klev
            detr(ig,k)=f(ig)*detr_star(ig,k)
            if (detr(ig,k).lt.0.) then
c               print*,'detr1<0!!!'
            endif
            enddo
            do k=1,klev
            entr(ig,k)=f(ig)*entr_star(ig,k)
            if (entr(ig,k).lt.0.) then
c               print*,'entr1<0!!!'
            endif
         enddo
      enddo
c
c       do ig=1,ngrid
c          do l=1,klev
c          if (((detr(ig,l)+entr(ig,l)+alim(ig,l))*ptimestep).gt.
c     s          (masse(ig,l))) then  
c      print*,'d2+e2+a2>m2','ig=',ig,'l=',l,'lmax(ig)=',lmax(ig),'d+e+a='
c     s,(detr(ig,l)+entr(ig,l)+alim(ig,l))*ptimestep,'m=',masse(ig,l)
c      endif
c      enddo
c      enddo
c Calcul des flux

      do ig=1,ngrid
         do l=1,lmax(ig)
c         do l=1,klev
c             fmc(ig,l+1)=f(ig)*f_star(ig,l+1)
            fmc(ig,l+1)=fmc(ig,l)+alim(ig,l)+entr(ig,l)-detr(ig,l)
c        print*,'??!!','ig=',ig,'l=',l,'lmax=',lmax(ig),'lmix=',lmix(ig),
c     s  'e=',entr(ig,l),'d=',detr(ig,l),'a=',alim(ig,l),'f=',fmc(ig,l),
c     s  'f+1=',fmc(ig,l+1)
          if (fmc(ig,l+1).lt.0.) then
               print*,'fmc1<0',l+1,lmax(ig),fmc(ig,l+1)
               fmc(ig,l+1)=fmc(ig,l)
               detr(ig,l)=alim(ig,l)+entr(ig,l)
c               fmc(ig,l+1)=0.
c               print*,'fmc1<0',l+1,lmax(ig),fmc(ig,l+1)
            endif
c       if ((fmc(ig,l+1).gt.fmc(ig,l)).and.(l.gt.lentr(ig))) then
c          f_old=fmc(ig,l+1)
c          fmc(ig,l+1)=fmc(ig,l)
c          detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l+1)
c       endif
       
c        if ((fmc(ig,l+1).gt.fmc(ig,l)).and.(l.gt.lentr(ig))) then
c          f_old=fmc(ig,l+1)
c          fmc(ig,l+1)=fmc(ig,l)
c          detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l)
c       endif
crajout du test sur alpha croissant
cif test
c       if (1.eq.0) then

       if (l.eq.klev) then
          print*,'THERMCELL PB ig=',ig,'   l=',l
          abort_message = 'THERMCELL PB'
          CALL abort_gcm (modname,abort_message,1)
       endif
!       if ((zw2(ig,l+1).gt.1.e-10).and.(zw2(ig,l).gt.1.e-10).and.
!     s     (l.ge.lentr(ig)).and.
       if ((zw2(ig,l+1).gt.1.e-10).and.(zw2(ig,l).gt.1.e-10).and.
     s         (l.ge.lentr(ig)) ) then
          if ( ((fmc(ig,l+1)/(rhobarz(ig,l+1)*zw2(ig,l+1))).gt.
     s     (fmc(ig,l)/(rhobarz(ig,l)*zw2(ig,l))))) then
           f_old=fmc(ig,l+1)
           fmc(ig,l+1)=fmc(ig,l)*rhobarz(ig,l+1)*zw2(ig,l+1)
     s                          /(rhobarz(ig,l)*zw2(ig,l))
           detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l+1) 
c           detr(ig,l)=(fmc(ig,l+1)-fmc(ig,l))/(0.4-1.)
c           entr(ig,l)=0.4*detr(ig,l)
c           entr(ig,l)=fmc(ig,l+1)-fmc(ig,l)+detr(ig,l)
        endif
        endif
        if ((fmc(ig,l+1).gt.fmc(ig,l)).and.(l.gt.lentr(ig))) then
          f_old=fmc(ig,l+1)
          fmc(ig,l+1)=fmc(ig,l)
          detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l+1)
       endif
       if (detr(ig,l).gt.fmc(ig,l)) then
               detr(ig,l)=fmc(ig,l)
               entr(ig,l)=fmc(ig,l+1)-alim(ig,l)
        endif
       if (fmc(ig,l+1).lt.0.) then
               detr(ig,l)=detr(ig,l)+fmc(ig,l+1)
               fmc(ig,l+1)=0.
               print*,'fmc2<0',l+1,lmax(ig)
            endif
            
ctest pour ne pas avoir f=0 et d=e/=0
c       if (fmc(ig,l+1).lt.1.e-10) then
c          detr(ig,l+1)=0.
c          entr(ig,l+1)=0.
c          zqla(ig,l+1)=0.
c          zw2(ig,l+1)=0.
c          lmax(ig)=l+1
c          zmax(ig)=zlev(ig,lmax(ig))
c       endif 
        if (zw2(ig,l+1).gt.1.e-10) then
       if ((((fmc(ig,l+1))/(rhobarz(ig,l+1)*zw2(ig,l+1))).gt.
     s      1.)) then
          f_old=fmc(ig,l+1)
          fmc(ig,l+1)=rhobarz(ig,l+1)*zw2(ig,l+1)
          zw2(ig,l+1)=0.
          zqla(ig,l+1)=0.
          detr(ig,l)=detr(ig,l)+f_old-fmc(ig,l+1)
         lmax(ig)=l+1
          zmax(ig)=zlev(ig,lmax(ig))
          print*,'alpha>1',l+1,lmax(ig)
       endif
        endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cendif test
c      endif
      enddo
      enddo
      do ig=1,ngrid
c         if (fmc(ig,lmax(ig)+1).ne.0.) then
         fmc(ig,lmax(ig)+1)=0.
         entr(ig,lmax(ig))=0.
         detr(ig,lmax(ig))=fmc(ig,lmax(ig))+entr(ig,lmax(ig))
     s                     +alim(ig,lmax(ig))
c         endif
      enddo
ctest sur le signe de fmc
       do ig=1,ngrid
         do l=1,klev+1
            if (fmc(ig,l).lt.0.) then
         print*,'fm1<0!!!','ig=',ig,'l=',l,'a=',alim(ig,l-1),'e='
     s ,entr(ig,l-1),'f=',fmc(ig,l-1),'d=',detr(ig,l-1),'f+1=',fmc(ig,l)
            endif
         enddo
       enddo
ctest de verification
      do ig=1,ngrid
       do l=1,lmax(ig)
       if ((abs(fmc(ig,l+1)-fmc(ig,l)-alim(ig,l)-entr(ig,l)+detr(ig,l)))
     s           .gt.1.e-4) then
c      print*,'pbcm!!','ig=',ig,'l=',l,'lmax=',lmax(ig),'lmix=',lmix(ig),
c     s  'e=',entr(ig,l),'d=',detr(ig,l),'a=',alim(ig,l),'f=',fmc(ig,l),
c     s  'f+1=',fmc(ig,l+1)
       endif
       if (detr(ig,l).lt.0.) then
          print*,'detrdemi<0!!!'
       endif
         enddo
      enddo
c
cRC
cCR def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1.) then
c test 
              if (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).gt.1e-10)
     s        then
c             
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))
     s        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
              else
              zmix(ig)=zlev(ig,lmix(ig))
              print*,'pb zmix'
              endif
          else 
              zmix(ig)=0.
          endif
ctest
         if ((zmax(ig)-zmix(ig)).le.0.) then
            zmix(ig)=0.9*zmax(ig)
c            print*,'pb zmix>zmax'
         endif
      enddo
      do ig=1,klon
         zmix0(ig)=zmix(ig)
      enddo
c
c calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,klev
            if (zmix(ig).ge.zlev(ig,l).and.
     s          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
c
cne devrait pas arriver!!!!!
      do ig=1,ngrid
       do l=1,klev
          if (detr(ig,l).gt.(fmc(ig,l)+alim(ig,l))+entr(ig,l)) then
             print*,'detr2>fmc2!!!','ig=',ig,'l=',l,'d=',detr(ig,l),
     s             'f=',fmc(ig,l),'lmax=',lmax(ig)
c             detr(ig,l)=fmc(ig,l)+alim(ig,l)+entr(ig,l)
c             entr(ig,l)=0.
c             fmc(ig,l+1)=0.
c             zw2(ig,l+1)=0.     
c             zqla(ig,l+1)=0.
           print*,'pb!fm=0 et f_star>0',l,lmax(ig)        
c             lmax(ig)=l
          endif
       enddo
      enddo
      do ig=1,ngrid
         do l=lmax(ig)+1,klev+1
c            fmc(ig,l)=0.
c            detr(ig,l)=0.
c            entr(ig,l)=0.
c            zw2(ig,l)=0.
c            zqla(ig,l)=0.
         enddo
      enddo

cCalcul du detrainement lors du premier passage
c     print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig).and.(test(ig).eq.1)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig).and.(test(ig).eq.1)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  if ((l_mix*zlev(ig,l)).lt.0.)then
                   print*,'pb l_mix*zlev<0'
                  endif
cCR: test: nouvelle def de lambda
c                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
                  if (zw2(ig,l).gt.1.e-10) then
                  larg_detr(ig,l)=sqrt((l_mix/zw2(ig,l))*zlev(ig,l))
                  else
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
                  endif
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
         endif
         enddo
       enddo

c     print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   cal1cul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.
c
c
      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1..and.(test(ig).eq.1)) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
c test
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo                  
cCR: calcul de fracazmix
       do ig=1,ngrid
          if (test(ig).eq.1) then
          fracazmix(ig)=(fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/
     s     (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig)
     s    +fraca(ig,lmix(ig))-zlev(ig,lmix(ig))*(fraca(ig,lmix(ig)+1)
     s    -fraca(ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
          endif
       enddo
c
       do l=2,nlay
          do ig=1,ngrid
             if(larg_cons(ig,l).gt.1..and.(test(ig).eq.1)) then
               if (l.gt.lmix(ig)) then
ctest
                 if (zmax(ig)-zmix(ig).lt.1.e-10) then
c                   print*,'pb xxx'
                    xxx(ig,l)=(lmax(ig)+1.-l)/(lmax(ig)+1.-lmix(ig))
                 else
                 xxx(ig,l)=(zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
                 endif
           if (idetr.eq.0) then
               fraca(ig,l)=fracazmix(ig)
           else if (idetr.eq.1) then
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)**2
           endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
             endif
            endif
         enddo
      enddo

      print*,'fin calcul fraca'
c     print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              if (test(ig).eq.1) then
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
cCR:test
              if (alim(ig,l-1).lt.1e-10.and.fm(ig,l).gt.fm(ig,l-1)
     s            .and.l.gt.lmix(ig)) then
                 fm(ig,l)=fm(ig,l-1)
c                 write(1,*)'ajustement fm, l',l
              endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cRC
              endif
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1.and.(test(ig).eq.1)) then
              abort_message = 'fracd trop petit'
              CALL abort_gcm (modname,abort_message,1)
            else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay+1
         do ig=1,ngrid
            if (test(ig).eq.0) then
              fm(ig,l)=fmc(ig,l)
            endif
         enddo
      enddo 
   
cfin du first
      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

      print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
      print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
     s         ,fm(ig,l+1)*ptimestep
     s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if((alim(ig,l)+entr(ig,l))*ptimestep.gt.masse(ig,l)) then
      print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
     s         ,(entr(ig,l)+alim(ig,l))*ptimestep
     s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
c     print*,'WARN!!! masse exagere ig=',ig,'   l=',l
c    s         ,'   M=',masse(ig,l)
c     print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
c    s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
c     print*,'zlev(ig,l+1),zlev(ig,l)'
c    s                ,zlev(ig,l+1),zlev(ig,l)
c     print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
c    s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
            endif
            if(.not.alim(ig,l).ge.0..or..not.alim(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

cCR:redefinition du entr
cCR:test:on ne change pas la def du entr mais la def du fm
       do l=1,nlay
         do ig=1,ngrid
            if (test(ig).eq.1) then
            detr(ig,l)=fm(ig,l)+alim(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
                fm(ig,l+1)=fm(ig,l)+alim(ig,l)
                detr(ig,l)=0.
c                write(11,*)'l,ig,entr',l,ig,entr(ig,l)
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
            endif
         enddo
      enddo
cRC

      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/REAL(tho)
         entr0=entr0+ptimestep*(alim+entr-entr0)/REAL(tho)
      else
         fm0=fm
         entr0=alim+entr
         detr0=detr
         alim0=alim
c         zoa=zqta
c         entr0=alim
      endif

      if (1.eq.1) then
c         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
c     .    ,zh,zdhadj,zha)
c         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
c     .    ,zo,pdoadj,zoa)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse,
     .                    zthl,zdthladj,zta)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse,
     .                   po,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

cCalcul des moments
c      do l=1,nlay
c         do ig=1,ngrid
c            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
c            zf2=zf/(1.-zf)
c            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
c            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
c         enddo
c      enddo






c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
c            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
           pdtadj(ig,l)=zdthladj(ig,l)*zpspsk(ig,l)  
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

      print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------
ccalcul de fraca pour les sorties
      do l=2,klev
         do ig=1,klon
            if (zw2(ig,l).gt.1.e-10) then
            fraca(ig,l)=fm(ig,l)/(rhobarz(ig,l)*zw2(ig,l))
            else
            fraca(ig,l)=0.
            endif
         enddo
      enddo
      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo
c CR calcul du niveau de condensation
c initialisation
      do ig=1,ngrid
         nivcon(ig)=0.
         zcon(ig)=0.
      enddo 
      do k=nlay,1,-1
         do ig=1,ngrid
            if (zqla(ig,k).gt.1e-10) then
               nivcon(ig)=k
               zcon(ig)=zlev(ig,k)
            endif
c            if (zcon(ig).gt.1.e-10) then
c               nuage=.true.
c            else 
c               nuage=.false.
c            endif
         enddo
      enddo
      
      do l=1,nlay
         do ig=1,ngrid
            zf=fraca(ig,l)
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l)/zpspsk(ig,l))**2
            wth2(ig,l)=zf2*(zw2(ig,l))**2
c           print*,'wth2=',wth2(ig,l)
            wth3(ig,l)=zf2*(1-2.*fraca(ig,l))/(1-fraca(ig,l))
     s                *zw2(ig,l)*zw2(ig,l)*zw2(ig,l)
            q2(ig,l)=zf2*(zqta(ig,l)*1000.-po(ig,l)*1000.)**2
ctest: on calcul q2/po=ratqsc
c            if (nuage) then
            ratqscth(ig,l)=sqrt(q2(ig,l))/(po(ig,l)*1000.)
c            else
c            ratqscth(ig,l)=0.
c            endif
         enddo
      enddo
ccalcul du ratqscdiff
      sum=0.
      sumdiff=0.
      ratqsdiff(:,:)=0.
      do ig=1,ngrid
         do l=1,lentr(ig)
            sum=sum+alim_star(ig,l)*zqta(ig,l)*1000.
         enddo
      enddo
      do ig=1,ngrid
          do l=1,lentr(ig)
          zf=fraca(ig,l)
          zf2=zf/(1.-zf)
       sumdiff=sumdiff+alim_star(ig,l)
     s           *(zqta(ig,l)*1000.-sum)**2
c      ratqsdiff=ratqsdiff+alim_star(ig,l)*
c     s          (zqta(ig,l)*1000.-po(ig,l)*1000.)**2
          enddo
      enddo
      do l=1,klev
      do ig=1,ngrid
      ratqsdiff(ig,l)=sqrt(sumdiff)/(po(ig,l)*1000.)   
c      write(11,*)'ratqsdiff=',ratqsdiff(ig,l)
      enddo
      enddo     
cdeja fait
c      do l=1,nlay
c         do ig=1,ngrid
c            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
c            if (detr(ig,l).lt.0.) then
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
c                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
c            endif
c         enddo
c      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1       


c #define und
	goto 123
#ifdef und
      CALL writeg1d(1,nlay,wd,'wd      ','wd      ')
      CALL writeg1d(1,nlay,zwa,'wa      ','wa      ')
      CALL writeg1d(1,nlay,fracd,'fracd      ','fracd      ')
      CALL writeg1d(1,nlay,fraca,'fraca      ','fraca      ')
      CALL writeg1d(1,nlay,wa_moy,'wam         ','wam         ')
      CALL writeg1d(1,nlay,zla,'la      ','la      ')
      CALL writeg1d(1,nlay,zld,'ld      ','ld      ')
      CALL writeg1d(1,nlay,pt,'pt      ','pt      ')
      CALL writeg1d(1,nlay,zh,'zh      ','zh      ')
      CALL writeg1d(1,nlay,zha,'zha      ','zha      ')
      CALL writeg1d(1,nlay,zu,'zu      ','zu      ')
      CALL writeg1d(1,nlay,zv,'zv      ','zv      ')
      CALL writeg1d(1,nlay,zo,'zo      ','zo      ')
      CALL writeg1d(1,nlay,wh,'wh      ','wh      ')
      CALL writeg1d(1,nlay,wu,'wu      ','wu      ')
      CALL writeg1d(1,nlay,wv,'wv      ','wv      ')
      CALL writeg1d(1,nlay,wo,'w15uo     ','wXo     ')
      CALL writeg1d(1,nlay,zdhadj,'zdhadj      ','zdhadj      ')
      CALL writeg1d(1,nlay,pduadj,'pduadj      ','pduadj      ')
      CALL writeg1d(1,nlay,pdvadj,'pdvadj      ','pdvadj      ')
      CALL writeg1d(1,nlay,pdoadj,'pdoadj      ','pdoadj      ')
      CALL writeg1d(1,nlay,entr  ,'entr        ','entr        ')
      CALL writeg1d(1,nlay,detr  ,'detr        ','detr        ')
      CALL writeg1d(1,nlay,fm    ,'fm          ','fm          ')

      CALL writeg1d(1,nlay,pdtadj,'pdtadj    ','pdtadj    ')
      CALL writeg1d(1,nlay,pplay,'pplay     ','pplay     ')
      CALL writeg1d(1,nlay,pplev,'pplev     ','pplev     ')

c   recalcul des flux en diagnostique...
c     print*,'PAS DE TEMPS ',ptimestep
       call dt2F(pplev,pplay,pt,pdtadj,wh)
      CALL writeg1d(1,nlay,wh,'wh2     ','wh2     ')
#endif
123   continue
! #define troisD
#ifdef troisD
c       if (sorties) then
      print*,'Debut des wrgradsfi'

c      print*,'16 OK convect8'
c         call wrgradsfi(1,nlay,wd,'wd        ','wd        ')
c         call wrgradsfi(1,nlay,zwa,'wa        ','wa        ')
         call wrgradsfi(1,nlay,fracd,'fracd     ','fracd     ')
         call wrgradsfi(1,nlay,fraca,'fraca     ','fraca     ')
c         call wrgradsfi(1,nlay,xxx,'xxx       ','xxx       ')
c         call wrgradsfi(1,nlay,wa_moy,'wam       ','wam       ')
c      print*,'WA6 ',wa_moy
c         call wrgradsfi(1,nlay,zla,'la        ','la        ')
c         call wrgradsfi(1,nlay,zld,'ld        ','ld        ')
         call wrgradsfi(1,nlay,pt,'pt        ','pt        ')
         call wrgradsfi(1,nlay,zh,'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha,'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua,'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva,'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu,'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv,'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,wh,'wh        ','wh        ')
         call wrgradsfi(1,nlay,wu,'wu        ','wu        ')
         call wrgradsfi(1,nlay,wv,'wv        ','wv        ')
         call wrgradsfi(1,nlay,wo,'wo        ','wo        ')
         call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
         call wrgradsfi(1,nlay,zdhadj,'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj,'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj,'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj,'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr,'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr,'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm,'fm        ','fm        ')
         call wrgradsfi(1,nlay,fmc,'fmc       ','fmc       ')
         call wrgradsfi(1,nlay,zw2,'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,w_est,'w_est      ','w_est      ')
con sort les moments
         call wrgradsfi(1,nlay,thetath2,'zh2       ','zh2       ')
         call wrgradsfi(1,nlay,wth2,'w2       ','w2       ')
         call wrgradsfi(1,nlay,wth3,'w3       ','w3       ')
         call wrgradsfi(1,nlay,q2,'q2       ','q2       ')
         call wrgradsfi(1,nlay,dtheta,'dT       ','dT       ')
c
         call wrgradsfi(1,nlay,zw_sec,'zw_sec       ','zw_sec       ')
         call wrgradsfi(1,nlay,ztva,'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv,'ztv       ','ztv       ')
         call wrgradsfi(1,nlay,nu,'nu       ','nu       ')

         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,zoa,'zoa        ','zoa        ')
c         call wrgradsfi(1,nlay,larg_cons,'Lc        ','Lc        ')
c         call wrgradsfi(1,nlay,larg_detr,'Ldetr     ','Ldetr     ')

cAM:nouveaux diagnostiques
         call wrgradsfi(1,nlay,zthl,'zthl        ','zthl        ')
         call wrgradsfi(1,nlay,zta,'zta        ','zta        ')
         call wrgradsfi(1,nlay,zl,'zl        ','zl        ')
         call wrgradsfi(1,nlay,zdthladj,'zdthladj    ',
     s        'zdthladj    ')
         call wrgradsfi(1,nlay,ztla,'ztla      ','ztla      ')
         call wrgradsfi(1,nlay,zqta,'zqta      ','zqta      ')
         call wrgradsfi(1,nlay,zqla,'zqla      ','zqla      ')
         call wrgradsfi(1,nlay,deltaz,'deltaz      ','deltaz      ')
cCR:nouveaux diagnostiques
      call wrgradsfi(1,nlay,entr_star  ,'entr_star   ','entr_star   ')
      call wrgradsfi(1,nlay,detr_star  ,'detr_star   ','detr_star   ')     
      call wrgradsfi(1,nlay,f_star    ,'f_star   ','f_star   ')
      call wrgradsfi(1,nlay,zqsat    ,'zqsat   ','zqsat   ')
      call wrgradsfi(1,nlay,zqsatth    ,'qsath   ','qsath   ')
      call wrgradsfi(1,nlay,alim_star    ,'alim_star   ','alim_star   ')
      call wrgradsfi(1,nlay,alim    ,'alim   ','alim   ')
      call wrgradsfi(1,1,f,'f      ','f      ')
      call wrgradsfi(1,1,alim_star_tot,'a_s_t      ','a_s_t      ')
      call wrgradsfi(1,1,alim_star2,'a_2      ','a_2      ')
      call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
      call wrgradsfi(1,1,zmax_sec,'z_sec      ','z_sec      ')
c      call wrgradsfi(1,1,zmax_sec2,'zz_se      ','zz_se      ')
      call wrgradsfi(1,1,zmix,'zmix      ','zmix      ') 
      call wrgradsfi(1,1,nivcon,'nivcon      ','nivcon      ')
      call wrgradsfi(1,1,zcon,'zcon      ','zcon      ')
      zsortie1d(:)=lmax(:)
      call wrgradsfi(1,1,zsortie1d,'lmax      ','lmax      ')
      call wrgradsfi(1,1,wmax,'wmax      ','wmax      ')
      call wrgradsfi(1,1,wmax_sec,'w_sec      ','w_sec      ')
      zsortie1d(:)=lmix(:)
      call wrgradsfi(1,1,zsortie1d,'lmix      ','lmix      ')
      zsortie1d(:)=lentr(:)
      call wrgradsfi(1,1,zsortie1d,'lentr      ','lentr     ')

c      print*,'17 OK convect8'

         do k=1,klev/10
            write(str2,'(i2.2)') k
            str10='wa'//str2
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=wa(ig,k,l)
               enddo
            enddo   
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=larg_part(ig,k,l)
               enddo
            enddo
           str10='la'//str2
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
         enddo


c     print*,'18 OK convect8'
c      endif
      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

      print*,'19 OK convect8'
      return
      end
      SUBROUTINE thermcell_eau(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      USE dimphy
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"

c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/
c$OMP THREADPRIVATE(idetr)

c   local:
c   ------

      INTEGER ig,k,l,lmaxa(klon),lmix(klon)
      real zsortie1d(klon)
c CR: on remplace lmax(klon,klev+1)
      INTEGER lmax(klon),lmin(klon),lentr(klon)
      real linter(klon)
      real zmix(klon), fracazmix(klon) 
c RC 
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      real zthl(klon,klev),zdthladj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      real zl(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zta(klon,klev)
      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev)
!      common/comtherm/thetath2,wth2

      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit
c$OMP THREADPRIVATE(isplit)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

c     real wmax(klon,klev),wmaxa(klon)
      real wmax(klon),wmaxa(klon)
      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

      real zcor,zdelta,zcvm5,qlbef
      real Tbef(klon),qsatbef(klon)
      real dqsat_dT,DT,num,denom
      REAL REPS,RLvCp,DDT0
      real ztla(klon,klev),zqla(klon,klev),zqta(klon,klev)
 
      PARAMETER (DDT0=.01)

cCR:nouvelles variables
      real f_star(klon,klev+1),entr_star(klon,klev)
      real entr_star_tot(klon),entr_star2(klon)
      real f(klon), f0(klon)
      real zlevinter(klon)
      logical first
      data first /.false./
      save first
c$OMP THREADPRIVATE(first)

cRC

      character*2 str2
      character*10 str10

      character (len=20) :: modname='thermcell_eau'
      character (len=80) :: abort_message

      LOGICAL vtest(klon),down
      LOGICAL Zsat(klon)

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
c$OMP THREADPRIVATE(ncorrec)

c

c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c Initialisation
      RLvCp = RLVTT/RCPD
      REPS  = RD/RV
c
c-----------------------------------------------------------------------
cAM Calcul de T,q,ql a partir de Tl et qT
c   ---------------------------------------------------
c
c Pr Tprec=Tl calcul de qsat 
c Si qsat>qT T=Tl, q=qT
c Sinon DDT=(-Tprec+Tl+RLVCP (qT-qsat(T')) / (1+RLVCP dqsat/dt) 
c On cherche DDT < DDT0
c
c defaut
       DO  ll=1,nlay
         DO ig=1,ngrid
            zo(ig,ll)=po(ig,ll)
            zl(ig,ll)=0.
            zh(ig,ll)=pt(ig,ll)
         EndDO
       EndDO
       do ig=1,ngrid
          Zsat(ig)=.false.
       enddo
c
c
       DO ll=1,nlay
c les points insatures sont definitifs
         DO ig=1,ngrid
            Tbef(ig)=pt(ig,ll)
            zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
            qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,ll)
            qsatbef(ig)=MIN(0.5,qsatbef(ig))
            zcor=1./(1.-retv*qsatbef(ig))
            qsatbef(ig)=qsatbef(ig)*zcor
            Zsat(ig) = (max(0.,po(ig,ll)-qsatbef(ig)) .gt. 0.00001)
         EndDO

         DO ig=1,ngrid
           if (Zsat(ig)) then
            qlbef=max(0.,po(ig,ll)-qsatbef(ig))
c si sature: ql est surestime, d'ou la sous-relax
            DT = 0.5*RLvCp*qlbef
c on pourra enchainer 2 ou 3 calculs sans Do while
            do while (DT.gt.DDT0)
c il faut verifier si c,a conserve quand on repasse en insature ...
              Tbef(ig)=Tbef(ig)+DT
              zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
              qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,ll)
              qsatbef(ig)=MIN(0.5,qsatbef(ig))
              zcor=1./(1.-retv*qsatbef(ig))
              qsatbef(ig)=qsatbef(ig)*zcor
c on veut le signe de qlbef
              qlbef=po(ig,ll)-qsatbef(ig)
c          dqsat_dT
              zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
              zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
              zcor=1./(1.-retv*qsatbef(ig))
              dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef(ig),zcor)
              num=-Tbef(ig)+pt(ig,ll)+RLvCp*qlbef
              denom=1.+RLvCp*dqsat_dT
              DT=num/denom
            enddo
c on ecrit de maniere conservative (sat ou non)
            zl(ig,ll) = max(0.,qlbef)
c          T = Tl +Lv/Cp ql
            zh(ig,ll) = pt(ig,ll)+RLvCp*zl(ig,ll)
            zo(ig,ll) = po(ig,ll)-zl(ig,ll)
           endif
         EndDO
       EndDO
cAM fin
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

      print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
c            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
c            zo(ig,l)=po(ig,l)
c            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
cAM attention zh est maintenant le profil de T et plus le profil de theta !
c
c   T-> Theta
            ztv(ig,l)=zh(ig,l)/zpspsk(ig,l)
cAM Theta_v
            ztv(ig,l)=ztv(ig,l)*(1.+RETV*(zo(ig,l))
     s           -zl(ig,l))
cAM Thetal
            zthl(ig,l)=pt(ig,l)/zpspsk(ig,l)
c            
1015     CONTINUE
1010  CONTINUE

c     print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo

c     print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
c            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
             rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*ztv(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo

c     print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------

cCR: ponderation entrainement des couches instables
cdef des entr_star tels que entr=f*entr_star      
      do l=1,klev
         do ig=1,ngrid 
            entr_star(ig,l)=0.
         enddo
      enddo
c determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lentr(ig)=1
      enddo

con ne considere que les premieres couches instables
      do k=nlay-1,1,-1
         do ig=1,ngrid
            if (ztv(ig,k).gt.ztv(ig,k+1).and.
     s          ztv(ig,k+1).lt.ztv(ig,k+2)) then
               lentr(ig)=k
            endif
          enddo
      enddo
    
c determination du lmin: couche d ou provient le thermique
      do ig=1,ngrid
         lmin(ig)=1
      enddo
      do ig=1,ngrid
         do l=nlay,2,-1
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.le.lentr(ig)) then 
                 entr_star(ig,l)=(ztv(ig,l)-ztv(ig,l+1))*
     s                          (zlev(ig,l+1)-zlev(ig,l))
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
            do l=1,klev
               entr_star(ig,l)=0.
            enddo
         endif
      enddo 
c calcul de l entrainement total
      do ig=1,ngrid
         entr_star_tot(ig)=0.
      enddo
      do ig=1,ngrid
         do k=1,klev
            entr_star_tot(ig)=entr_star_tot(ig)+entr_star(ig,k)
         enddo
      enddo
c
      do k=1,klev
         do ig=1,ngrid 
            ztva(ig,k)=ztv(ig,k)
         enddo
      enddo
cRC
cAM:initialisations
      do k=1,nlay
         do ig=1,ngrid
            ztva(ig,k)=ztv(ig,k)
            ztla(ig,k)=zthl(ig,k)
            zqla(ig,k)=0.
            zqta(ig,k)=po(ig,k)
            Zsat(ig) =.false.
         enddo
      enddo
c
c     print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
cCR
            f_star(ig,k)=0.
cRC
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

c     print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo

cCR: 
      do l=1,nlay-2
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.entr_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
cAM
               ztla(ig,l)=zthl(ig,l)
               zqta(ig,l)=po(ig,l)
               zqla(ig,l)=zl(ig,l)
cAM
               f_star(ig,l+1)=entr_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               larg_detr(ig,l)=0.
            else if ((zw2(ig,l).ge.1e-10).and.
     s               (f_star(ig,l)+entr_star(ig,l).gt.1.e-10)) then
               f_star(ig,l+1)=f_star(ig,l)+entr_star(ig,l)
c
cAM on melange Tl et qt du thermique
               ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+entr_star(ig,l)
     s                    *zthl(ig,l))/f_star(ig,l+1)
               zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+entr_star(ig,l)
     s                    *po(ig,l))/f_star(ig,l+1)
c
c               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)
c     s                    *ztv(ig,l))/f_star(ig,l+1)
c
cAM on en deduit thetav et ql du thermique
               Tbef(ig)=ztla(ig,l)*zpspsk(ig,l)
               zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
               qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,l)
               qsatbef(ig)=MIN(0.5,qsatbef(ig))
               zcor=1./(1.-retv*qsatbef(ig))
               qsatbef(ig)=qsatbef(ig)*zcor
               Zsat(ig) = (max(0.,zqta(ig,l)-qsatbef(ig)) .gt. 0.00001)
            endif
         enddo
         DO ig=1,ngrid
           if (Zsat(ig)) then
              qlbef=max(0.,zqta(ig,l)-qsatbef(ig))
              DT = 0.5*RLvCp*qlbef
              do while (DT.gt.DDT0)
                 Tbef(ig)=Tbef(ig)+DT
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 qsatbef(ig)= R2ES * FOEEW(Tbef(ig),zdelta)/pplev(ig,l)
                 qsatbef(ig)=MIN(0.5,qsatbef(ig))
                 zcor=1./(1.-retv*qsatbef(ig))
                 qsatbef(ig)=qsatbef(ig)*zcor
                 qlbef=zqta(ig,l)-qsatbef(ig)

                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef(ig)))
                 zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
                 zcor=1./(1.-retv*qsatbef(ig))
                 dqsat_dT=FOEDE(Tbef(ig),zdelta,zcvm5,qsatbef(ig),zcor)
                 num=-Tbef(ig)+ztla(ig,l)*zpspsk(ig,l)+RLvCp*qlbef
                 denom=1.+RLvCp*dqsat_dT
                 DT=num/denom
              enddo
              zqla(ig,l) = max(0.,zqta(ig,l)-qsatbef(ig))
           endif
c on ecrit de maniere conservative (sat ou non)
c          T = Tl +Lv/Cp ql
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)
     s              -zqla(ig,l))-zqla(ig,l))

        enddo
        DO ig=1,ngrid
           if (zw2(ig,l).ge.1.e-10.and.
     s               f_star(ig,l)+entr_star(ig,l).gt.1.e-10) then
c  mise a jour de la vitesse ascendante (l'air entraine de la couche
c  consideree commence avec une vitesse nulle).
c
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+
     s                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
c determination de zmax continu par interpolation lineaire
            if (zw2(ig,l+1).lt.0.) then
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmaxa(ig)=l
            else
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
c
c Calcul de la couche correspondant a la hauteur du thermique
      do ig=1,ngrid
         lmax(ig)=lentr(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
            lmax(ig)=1
            lmin(ig)=1
         endif
      enddo 
c    
c Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=500.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
c calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*
     s    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)
     s    -zlev(ig,lmax(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

c Fermeture,determination de f
      do ig=1,ngrid
         entr_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (entr_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else
             do k=lmin(ig),lentr(ig)
                entr_star2(ig)=entr_star2(ig)+entr_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
c Nouvelle fermeture
             f(ig)=wmax(ig)/(zmax(ig)*r_aspect*entr_star2(ig))
     s             *entr_star_tot(ig)
ctest
             if (first) then
             f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
     s             *wmax(ig))
             endif
         endif
         f0(ig)=f(ig)
         first=.true.
      enddo

c Calcul de l'entrainement
       do k=1,klev
         do ig=1,ngrid 
            entr(ig,k)=f(ig)*entr_star(ig,k)
         enddo
      enddo
c Calcul des flux
      do ig=1,ngrid
         do l=1,lmax(ig)-1
            fmc(ig,l+1)=fmc(ig,l)+entr(ig,l)
         enddo
      enddo

cRC


c     print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
            endif
         enddo
       enddo

c     print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   calcul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.
c
cCR def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1.) then
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))
     s        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
         else 
         zmix(ig)=0.
         endif
      enddo
c
c calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,klev
            if (zmix(ig).ge.zlev(ig,l).and.
     s          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
c
      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1.) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
c test
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo                  
cCR: calcul de fracazmix
       do ig=1,ngrid
          fracazmix(ig)=(fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/
     s     (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig)
     s    +fraca(ig,lmix(ig))-zlev(ig,lmix(ig))*(fraca(ig,lmix(ig)+1)
     s    -fraca(ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
       enddo
c
       do l=2,nlay
          do ig=1,ngrid
             if(larg_cons(ig,l).gt.1.) then
               if (l.gt.lmix(ig)) then
                 xxx(ig,l)=(zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
           if (idetr.eq.0) then
               fraca(ig,l)=fracazmix(ig)
           else if (idetr.eq.1) then
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)**2
           endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
             endif
            endif
         enddo
      enddo

c     print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
cCR:test
              if (entr(ig,l-1).lt.1e-10.and.fm(ig,l).gt.fm(ig,l-1)
     s            .and.l.gt.lmix(ig)) then
                 fm(ig,l)=fm(ig,l-1)
c                 write(1,*)'ajustement fm, l',l
              endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cRC
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1) then
              abort_message = 'fracd trop petit'
              CALL abort_gcm (modname,abort_message,1)
            else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

c     print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
c     print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
c    s         ,fm(ig,l+1)*ptimestep
c    s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(entr(ig,l)*ptimestep.gt.masse(ig,l)) then
c     print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
c    s         ,entr(ig,l)*ptimestep
c    s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
c     print*,'WARN!!! masse exagere ig=',ig,'   l=',l
c    s         ,'   M=',masse(ig,l)
c     print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
c    s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
c     print*,'zlev(ig,l+1),zlev(ig,l)'
c    s                ,zlev(ig,l+1),zlev(ig,l)
c     print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
c    s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
            endif
            if(.not.entr(ig,l).ge.0..or..not.entr(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/REAL(tho)
         entr0=entr0+ptimestep*(entr-entr0)/REAL(tho)
      else
         fm0=fm
         entr0=entr
      endif

      if (1.eq.1) then
c         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
c     .    ,zh,zdhadj,zha)
c         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
c     .    ,zo,pdoadj,zoa)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zthl,zdthladj,zta)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,po,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

      do l=1,nlay
         do ig=1,ngrid
            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
         enddo
      enddo



c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
c            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
           pdtadj(ig,l)=zdthladj(ig,l)*zpspsk(ig,l)  
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

c     print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------

      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
                entr(ig,l)=entr(ig,l)-detr(ig,l)
                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
         enddo
      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1


c #define und
	goto 123
#ifdef und
      CALL writeg1d(1,nlay,wd,'wd      ','wd      ')
      CALL writeg1d(1,nlay,zwa,'wa      ','wa      ')
      CALL writeg1d(1,nlay,fracd,'fracd      ','fracd      ')
      CALL writeg1d(1,nlay,fraca,'fraca      ','fraca      ')
      CALL writeg1d(1,nlay,wa_moy,'wam         ','wam         ')
      CALL writeg1d(1,nlay,zla,'la      ','la      ')
      CALL writeg1d(1,nlay,zld,'ld      ','ld      ')
      CALL writeg1d(1,nlay,pt,'pt      ','pt      ')
      CALL writeg1d(1,nlay,zh,'zh      ','zh      ')
      CALL writeg1d(1,nlay,zha,'zha      ','zha      ')
      CALL writeg1d(1,nlay,zu,'zu      ','zu      ')
      CALL writeg1d(1,nlay,zv,'zv      ','zv      ')
      CALL writeg1d(1,nlay,zo,'zo      ','zo      ')
      CALL writeg1d(1,nlay,wh,'wh      ','wh      ')
      CALL writeg1d(1,nlay,wu,'wu      ','wu      ')
      CALL writeg1d(1,nlay,wv,'wv      ','wv      ')
      CALL writeg1d(1,nlay,wo,'w15uo     ','wXo     ')
      CALL writeg1d(1,nlay,zdhadj,'zdhadj      ','zdhadj      ')
      CALL writeg1d(1,nlay,pduadj,'pduadj      ','pduadj      ')
      CALL writeg1d(1,nlay,pdvadj,'pdvadj      ','pdvadj      ')
      CALL writeg1d(1,nlay,pdoadj,'pdoadj      ','pdoadj      ')
      CALL writeg1d(1,nlay,entr  ,'entr        ','entr        ')
      CALL writeg1d(1,nlay,detr  ,'detr        ','detr        ')
      CALL writeg1d(1,nlay,fm    ,'fm          ','fm          ')

      CALL writeg1d(1,nlay,pdtadj,'pdtadj    ','pdtadj    ')
      CALL writeg1d(1,nlay,pplay,'pplay     ','pplay     ')
      CALL writeg1d(1,nlay,pplev,'pplev     ','pplev     ')

c   recalcul des flux en diagnostique...
c     print*,'PAS DE TEMPS ',ptimestep
       call dt2F(pplev,pplay,pt,pdtadj,wh)
      CALL writeg1d(1,nlay,wh,'wh2     ','wh2     ')
#endif
123   continue
! #define troisD
#ifdef troisD
c       if (sorties) then
      print*,'Debut des wrgradsfi'

c      print*,'16 OK convect8'
         call wrgradsfi(1,nlay,wd,'wd        ','wd        ')
         call wrgradsfi(1,nlay,zwa,'wa        ','wa        ')
         call wrgradsfi(1,nlay,fracd,'fracd     ','fracd     ')
         call wrgradsfi(1,nlay,fraca,'fraca     ','fraca     ')
         call wrgradsfi(1,nlay,xxx,'xxx       ','xxx       ')
         call wrgradsfi(1,nlay,wa_moy,'wam       ','wam       ')
c      print*,'WA6 ',wa_moy
         call wrgradsfi(1,nlay,zla,'la        ','la        ')
         call wrgradsfi(1,nlay,zld,'ld        ','ld        ')
         call wrgradsfi(1,nlay,pt,'pt        ','pt        ')
         call wrgradsfi(1,nlay,zh,'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha,'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua,'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva,'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu,'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv,'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,wh,'wh        ','wh        ')
         call wrgradsfi(1,nlay,wu,'wu        ','wu        ')
         call wrgradsfi(1,nlay,wv,'wv        ','wv        ')
         call wrgradsfi(1,nlay,wo,'wo        ','wo        ')
         call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
         call wrgradsfi(1,nlay,zdhadj,'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj,'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj,'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj,'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr,'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr,'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm,'fm        ','fm        ')
         call wrgradsfi(1,nlay,fmc,'fmc       ','fmc       ')
         call wrgradsfi(1,nlay,zw2,'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,ztva,'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv,'ztv       ','ztv       ')

         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,larg_cons,'Lc        ','Lc        ')
         call wrgradsfi(1,nlay,larg_detr,'Ldetr     ','Ldetr     ')

cAM:nouveaux diagnostiques
         call wrgradsfi(1,nlay,zthl,'zthl        ','zthl        ')
         call wrgradsfi(1,nlay,zta,'zta        ','zta        ')
         call wrgradsfi(1,nlay,zl,'zl        ','zl        ')
         call wrgradsfi(1,nlay,zdthladj,'zdthladj    ',
     s        'zdthladj    ')
         call wrgradsfi(1,nlay,ztla,'ztla      ','ztla      ')
         call wrgradsfi(1,nlay,zqta,'zqta      ','zqta      ')
         call wrgradsfi(1,nlay,zqla,'zqla      ','zqla      ')
cCR:nouveaux diagnostiques
      call wrgradsfi(1,nlay,entr_star  ,'entr_star   ','entr_star   ')     
      call wrgradsfi(1,nlay,f_star    ,'f_star   ','f_star   ')
      call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
      call wrgradsfi(1,1,zmix,'zmix      ','zmix      ') 
      zsortie1d(:)=lmax(:)
      call wrgradsfi(1,1,zsortie1d,'lmax      ','lmax      ')
      call wrgradsfi(1,1,wmax,'wmax      ','wmax      ')
      zsortie1d(:)=lmix(:)
      call wrgradsfi(1,1,zsortie1d,'lmix      ','lmix      ')
      zsortie1d(:)=lentr(:)
      call wrgradsfi(1,1,zsortie1d,'lentr      ','lentr     ')

c      print*,'17 OK convect8'

         do k=1,klev/10
            write(str2,'(i2.2)') k
            str10='wa'//str2
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=wa(ig,k,l)
               enddo
            enddo   
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=larg_part(ig,k,l)
               enddo
            enddo
            str10='la'//str2
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
         enddo


c     print*,'18 OK convect8'
c      endif
      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

c     print*,'19 OK convect8'
      return
      end

      SUBROUTINE thermcell(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      USE dimphy
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"

c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/
c$OMP THREADPRIVATE(idetr)

c   local:
c   ------

      INTEGER ig,k,l,lmaxa(klon),lmix(klon)
      real zsortie1d(klon)
c CR: on remplace lmax(klon,klev+1)
      INTEGER lmax(klon),lmin(klon),lentr(klon)
      real linter(klon)
      real zmix(klon), fracazmix(klon) 
c RC 
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev)
!      common/comtherm/thetath2,wth2

      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit
c$OMP THREADPRIVATE(isplit)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

c     real wmax(klon,klev),wmaxa(klon)
      real wmax(klon),wmaxa(klon)
      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

cCR:nouvelles variables
      real f_star(klon,klev+1),entr_star(klon,klev)
      real entr_star_tot(klon),entr_star2(klon)
      real f(klon), f0(klon)
      real zlevinter(klon)
      logical first
      data first /.false./
      save first
c$OMP THREADPRIVATE(first)
cRC

      character*2 str2
      character*10 str10

      character (len=20) :: modname='thermcell'
      character (len=80) :: abort_message

      LOGICAL vtest(klon),down

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
c$OMP THREADPRIVATE(ncorrec)
      
c
c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

       print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
            zo(ig,l)=po(ig,l)
            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
1015     CONTINUE
1010  CONTINUE

       print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo

c      print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo

c      print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------

cCR: ponderation entrainement des couches instables
cdef des entr_star tels que entr=f*entr_star      
      do l=1,klev
         do ig=1,ngrid 
            entr_star(ig,l)=0.
         enddo
      enddo
c determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lentr(ig)=1
      enddo

con ne considere que les premieres couches instables
      do k=nlay-2,1,-1
         do ig=1,ngrid
            if (ztv(ig,k).gt.ztv(ig,k+1).and.
     s          ztv(ig,k+1).le.ztv(ig,k+2)) then
               lentr(ig)=k
            endif
          enddo
      enddo
    
c determination du lmin: couche d ou provient le thermique
      do ig=1,ngrid
         lmin(ig)=1
      enddo
      do ig=1,ngrid
         do l=nlay,2,-1
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.le.lentr(ig)) then 
                 entr_star(ig,l)=(ztv(ig,l)-ztv(ig,l+1))*
     s                           (zlev(ig,l+1)-zlev(ig,l))
            endif
         enddo
      enddo
c pas de thermique si couches 1->5 stables
      do ig=1,ngrid
         if (lmin(ig).gt.5) then
            do l=1,klev
               entr_star(ig,l)=0.
            enddo
         endif
      enddo 
c calcul de l entrainement total
      do ig=1,ngrid
         entr_star_tot(ig)=0.
      enddo
      do ig=1,ngrid
         do k=1,klev
            entr_star_tot(ig)=entr_star_tot(ig)+entr_star(ig,k)
         enddo
      enddo
c
      print*,'fin calcul entr_star'
      do k=1,klev
         do ig=1,ngrid 
            ztva(ig,k)=ztv(ig,k)
         enddo
      enddo
cRC
c      print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
cCR
            f_star(ig,k)=0.
cRC
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

c      print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo

cCR: 
      do l=1,nlay-2
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.entr_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
               f_star(ig,l+1)=entr_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               larg_detr(ig,l)=0.
            else if ((zw2(ig,l).ge.1e-10).and.
     s               (f_star(ig,l)+entr_star(ig,l).gt.1.e-10)) then
               f_star(ig,l+1)=f_star(ig,l)+entr_star(ig,l)
               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)
     s                    *ztv(ig,l))/f_star(ig,l+1)
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+
     s                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
c determination de zmax continu par interpolation lineaire
            if (zw2(ig,l+1).lt.0.) then
ctest
               if (abs(zw2(ig,l+1)-zw2(ig,l)).lt.1e-10) then
                  print*,'pb linter'
               endif
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmaxa(ig)=l
            else
               if (zw2(ig,l+1).lt.0.) then
                  print*,'pb1 zw2<0'
               endif
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
      print*,'fin calcul zw2'
c
c Calcul de la couche correspondant a la hauteur du thermique
      do ig=1,ngrid
         lmax(ig)=lentr(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couches 1->5 stables
      do ig=1,ngrid
         if (lmin(ig).gt.5) then
            lmax(ig)=1
            lmin(ig)=1
         endif
      enddo 
c    
c Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
                  print*,'pb2 zw2<0'
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
c calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*
     s    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)
     s    -zlev(ig,lmax(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

      print*,'avant fermeture'
c Fermeture,determination de f
      do ig=1,ngrid
         entr_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (entr_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else
             do k=lmin(ig),lentr(ig)
                entr_star2(ig)=entr_star2(ig)+entr_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
c Nouvelle fermeture
             f(ig)=wmax(ig)/(max(500.,zmax(ig))*r_aspect
     s             *entr_star2(ig))*entr_star_tot(ig)
ctest
c             if (first) then
c             f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
c     s             *wmax(ig))
c             endif
         endif
c         f0(ig)=f(ig)
c         first=.true.
      enddo
      print*,'apres fermeture'

c Calcul de l'entrainement
       do k=1,klev
         do ig=1,ngrid 
            entr(ig,k)=f(ig)*entr_star(ig,k)
         enddo
      enddo
c Calcul des flux
      do ig=1,ngrid
         do l=1,lmax(ig)-1
            fmc(ig,l+1)=fmc(ig,l)+entr(ig,l)
         enddo
      enddo

cRC


c      print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  if ((l_mix*zlev(ig,l)).lt.0.)then
                   print*,'pb l_mix*zlev<0'
                  endif
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
            endif
         enddo
       enddo

c      print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   calcul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.
c
cCR def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1.) then
c test 
              if (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).gt.1e-10)
     s        then
c             
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))
     s        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
            else
            zmix(ig)=zlev(ig,lmix(ig))
            print*,'pb zmix'
            endif
         else 
         zmix(ig)=0.
         endif
ctest
         if ((zmax(ig)-zmix(ig)).lt.0.) then
            zmix(ig)=0.99*zmax(ig)
c            print*,'pb zmix>zmax'
         endif
      enddo
c
c calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,klev
            if (zmix(ig).ge.zlev(ig,l).and.
     s          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
c
      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1.) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
c test
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo                  
cCR: calcul de fracazmix
       do ig=1,ngrid
          fracazmix(ig)=(fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/
     s     (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig)
     s    +fraca(ig,lmix(ig))-zlev(ig,lmix(ig))*(fraca(ig,lmix(ig)+1)
     s    -fraca(ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
       enddo
c
       do l=2,nlay
          do ig=1,ngrid
             if(larg_cons(ig,l).gt.1.) then
               if (l.gt.lmix(ig)) then
ctest
                 if (zmax(ig)-zmix(ig).lt.1.e-10) then
c                   print*,'pb xxx'
                   xxx(ig,l)=(lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
                 else
                 xxx(ig,l)=(zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
                 endif
           if (idetr.eq.0) then
               fraca(ig,l)=fracazmix(ig)
           else if (idetr.eq.1) then
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)**2
           endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
             endif
            endif
         enddo
      enddo
      
      print*,'fin calcul fraca'
c      print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
cCR:test
              if (entr(ig,l-1).lt.1e-10.and.fm(ig,l).gt.fm(ig,l-1)
     s            .and.l.gt.lmix(ig)) then
                 fm(ig,l)=fm(ig,l-1)
c                 write(1,*)'ajustement fm, l',l
              endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cRC
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1) then
              abort_message = 'fracd trop petit'
              CALL abort_gcm (modname,abort_message,1)
            else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

      print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
c     print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
c    s         ,fm(ig,l+1)*ptimestep
c    s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(entr(ig,l)*ptimestep.gt.masse(ig,l)) then
c     print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
c    s         ,entr(ig,l)*ptimestep
c    s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
c     print*,'WARN!!! masse exagere ig=',ig,'   l=',l
c    s         ,'   M=',masse(ig,l)
c     print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
c    s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
c     print*,'zlev(ig,l+1),zlev(ig,l)'
c    s                ,zlev(ig,l+1),zlev(ig,l)
c     print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
c    s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
            endif
            if(.not.entr(ig,l).ge.0..or..not.entr(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

cCR:redefinition du entr
       do l=1,nlay
         do ig=1,ngrid
            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
                entr(ig,l)=entr(ig,l)-detr(ig,l)
                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
         enddo
      enddo
cRC
      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/REAL(tho)
         entr0=entr0+ptimestep*(entr-entr0)/REAL(tho)
      else
         fm0=fm
         entr0=entr
      endif

      if (1.eq.1) then
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zh,zdhadj,zha)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zo,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

      do l=1,nlay
         do ig=1,ngrid
            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
         enddo
      enddo



c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

      print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------

      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo

cdeja fait
c      do l=1,nlay
c         do ig=1,ngrid
c            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
c            if (detr(ig,l).lt.0.) then
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
c                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
c            endif
c         enddo
c      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1


c #define und
	goto 123
#ifdef und
      CALL writeg1d(1,nlay,wd,'wd      ','wd      ')
      CALL writeg1d(1,nlay,zwa,'wa      ','wa      ')
      CALL writeg1d(1,nlay,fracd,'fracd      ','fracd      ')
      CALL writeg1d(1,nlay,fraca,'fraca      ','fraca      ')
      CALL writeg1d(1,nlay,wa_moy,'wam         ','wam         ')
      CALL writeg1d(1,nlay,zla,'la      ','la      ')
      CALL writeg1d(1,nlay,zld,'ld      ','ld      ')
      CALL writeg1d(1,nlay,pt,'pt      ','pt      ')
      CALL writeg1d(1,nlay,zh,'zh      ','zh      ')
      CALL writeg1d(1,nlay,zha,'zha      ','zha      ')
      CALL writeg1d(1,nlay,zu,'zu      ','zu      ')
      CALL writeg1d(1,nlay,zv,'zv      ','zv      ')
      CALL writeg1d(1,nlay,zo,'zo      ','zo      ')
      CALL writeg1d(1,nlay,wh,'wh      ','wh      ')
      CALL writeg1d(1,nlay,wu,'wu      ','wu      ')
      CALL writeg1d(1,nlay,wv,'wv      ','wv      ')
      CALL writeg1d(1,nlay,wo,'w15uo     ','wXo     ')
      CALL writeg1d(1,nlay,zdhadj,'zdhadj      ','zdhadj      ')
      CALL writeg1d(1,nlay,pduadj,'pduadj      ','pduadj      ')
      CALL writeg1d(1,nlay,pdvadj,'pdvadj      ','pdvadj      ')
      CALL writeg1d(1,nlay,pdoadj,'pdoadj      ','pdoadj      ')
      CALL writeg1d(1,nlay,entr  ,'entr        ','entr        ')
      CALL writeg1d(1,nlay,detr  ,'detr        ','detr        ')
      CALL writeg1d(1,nlay,fm    ,'fm          ','fm          ')

      CALL writeg1d(1,nlay,pdtadj,'pdtadj    ','pdtadj    ')
      CALL writeg1d(1,nlay,pplay,'pplay     ','pplay     ')
      CALL writeg1d(1,nlay,pplev,'pplev     ','pplev     ')

c   recalcul des flux en diagnostique...
c     print*,'PAS DE TEMPS ',ptimestep
       call dt2F(pplev,pplay,pt,pdtadj,wh)
      CALL writeg1d(1,nlay,wh,'wh2     ','wh2     ')
#endif
123   continue
#define troisD
#ifdef troisD
c       if (sorties) then
      print*,'Debut des wrgradsfi'

c      print*,'16 OK convect8'
         call wrgradsfi(1,nlay,wd,'wd        ','wd        ')
         call wrgradsfi(1,nlay,zwa,'wa        ','wa        ')
         call wrgradsfi(1,nlay,fracd,'fracd     ','fracd     ')
         call wrgradsfi(1,nlay,fraca,'fraca     ','fraca     ')
         call wrgradsfi(1,nlay,xxx,'xxx       ','xxx       ')
         call wrgradsfi(1,nlay,wa_moy,'wam       ','wam       ')
c      print*,'WA6 ',wa_moy
         call wrgradsfi(1,nlay,zla,'la        ','la        ')
         call wrgradsfi(1,nlay,zld,'ld        ','ld        ')
         call wrgradsfi(1,nlay,pt,'pt        ','pt        ')
         call wrgradsfi(1,nlay,zh,'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha,'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua,'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva,'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu,'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv,'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,wh,'wh        ','wh        ')
         call wrgradsfi(1,nlay,wu,'wu        ','wu        ')
         call wrgradsfi(1,nlay,wv,'wv        ','wv        ')
         call wrgradsfi(1,nlay,wo,'wo        ','wo        ')
         call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
         call wrgradsfi(1,nlay,zdhadj,'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj,'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj,'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj,'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr,'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr,'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm,'fm        ','fm        ')
         call wrgradsfi(1,nlay,fmc,'fmc       ','fmc       ')
         call wrgradsfi(1,nlay,zw2,'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,ztva,'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv,'ztv       ','ztv       ')

         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,larg_cons,'Lc        ','Lc        ')
         call wrgradsfi(1,nlay,larg_detr,'Ldetr     ','Ldetr     ')

cCR:nouveaux diagnostiques
      call wrgradsfi(1,nlay,entr_star  ,'entr_star   ','entr_star   ')     
      call wrgradsfi(1,nlay,f_star    ,'f_star   ','f_star   ')
      call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
      call wrgradsfi(1,1,zmix,'zmix      ','zmix      ') 
      zsortie1d(:)=lmax(:)
      call wrgradsfi(1,1,zsortie1d,'lmax      ','lmax      ')
      call wrgradsfi(1,1,wmax,'wmax      ','wmax      ')
      zsortie1d(:)=lmix(:)
      call wrgradsfi(1,1,zsortie1d,'lmix      ','lmix      ')
      zsortie1d(:)=lentr(:)
      call wrgradsfi(1,1,zsortie1d,'lentr      ','lentr     ')

c      print*,'17 OK convect8'

         do k=1,klev/10
            write(str2,'(i2.2)') k
            str10='wa'//str2
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=wa(ig,k,l)
               enddo
            enddo   
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=larg_part(ig,k,l)
               enddo
            enddo
            str10='la'//str2
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
         enddo


c     print*,'18 OK convect8'
c      endif
      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

      print*,'19 OK convect8'
      return
      end

      subroutine dqthermcell(ngrid,nlay,ptimestep,fm,entr,
     .           masse,q,dq,qa)
      USE dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================

#include "dimensions.h"
cccc#include "dimphy.h"

      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),wqd(klon,klev+1)

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
ctest
            if (detr(ig,k).lt.0.) then
               entr(ig,k)=entr(ig,k)-detr(ig,k)
               detr(ig,k)=0.
c               print*,'detr2<0!!!','ig=',ig,'k=',k,'f=',fm(ig,k),
c     s         'f+1=',fm(ig,k+1),'e=',entr(ig,k),'d=',detr(ig,k)
            endif
            if (fm(ig,k+1).lt.0.) then
c               print*,'fm2<0!!!'
            endif
            if (entr(ig,k).lt.0.) then
c               print*,'entr2<0!!!'
            endif
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
         qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))
     s         /(fm(ig,k+1)+detr(ig,k))
            else
               qa(ig,k)=q(ig,k)
            endif
            if (qa(ig,k).lt.0.) then
c               print*,'qa<0!!!'
            endif
            if (q(ig,k).lt.0.) then
c               print*,'q<0!!!'
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
c             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
            wqd(ig,k)=fm(ig,k)*q(ig,k)
            if (wqd(ig,k).lt.0.) then
c               print*,'wqd<0!!!'
            endif
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo
     
      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(detr(ig,k)*qa(ig,k)-entr(ig,k)*q(ig,k)
     s               -wqd(ig,k)+wqd(ig,k+1))
     s               /masse(ig,k)
c            if (dq(ig,k).lt.0.) then
c               print*,'dq<0!!!'
c            endif
         enddo
      enddo

      return
      end
      subroutine dvthermcell(ngrid,nlay,ptimestep,fm,entr,masse
     .    ,fraca,larga
     .    ,u,v,du,dv,ua,va)
      USE dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================

#include "dimensions.h"
cccc#include "dimphy.h"

      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev)
      real wvd(klon,klev+1),wud(klon,klev+1)
      real gamma0,gamma(klon,klev+1)
      real dua,dva
      integer iter

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
c   On itère sur la valeur du coeff de freinage.
c              gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
               gamma0=masse(ig,k)
     s         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )
     s         *0.5/larga(ig)
c              gamma0=0.
c   la première fois on multiplie le coefficient de freinage
c   par le module du vent dans la couche en dessous.
               dua=ua(ig,k-1)-u(ig,k-1)
               dva=va(ig,k-1)-v(ig,k-1)
               do iter=1,5
                  gamma(ig,k)=gamma0*sqrt(dua**2+dva**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1)
     s               +(entr(ig,k)+gamma(ig,k))*u(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1)
     s               +(entr(ig,k)+gamma(ig,k))*v(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+gamma(ig,k))
c                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
                  dua=ua(ig,k)-u(ig,k)
                  dva=va(ig,k)-v(ig,k)
               enddo
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               gamma(ig,k)=0.
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*u(ig,k)
            wvd(ig,k)=fm(ig,k)*v(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*u(ig,k)
     s               -wud(ig,k)+wud(ig,k+1))
     s               /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*v(ig,k)
     s               -wvd(ig,k)+wvd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      subroutine dqthermcell2(ngrid,nlay,ptimestep,fm,entr,masse,frac
     .    ,q,dq,qa)
      USE dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================

#include "dimensions.h"
cccc#include "dimphy.h"

      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real entr(ngrid,nlay),frac(ngrid,nlay)
      real q(ngrid,nlay)
      real dq(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),wqd(klon,klev+1)
      real qe(klon,klev),zf,zf2

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         qa(ig,1)=q(ig,1)
         qe(ig,1)=q(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
               zf=0.5*(frac(ig,k)+frac(ig,k+1))
               zf2=1./(1.-zf)
               qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+zf2*entr(ig,k)*q(ig,k))
     s         /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2)
               qe(ig,k)=(q(ig,k)-zf*qa(ig,k))*zf2
            else
               qa(ig,k)=q(ig,k)
               qe(ig,k)=q(ig,k)
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
c             wqd(ig,k)=fm(ig,k)*0.5*(q(ig,k-1)+q(ig,k))
            wqd(ig,k)=fm(ig,k)*qe(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wqd(ig,1)=0.
         wqd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            dq(ig,k)=(detr(ig,k)*qa(ig,k)-entr(ig,k)*qe(ig,k)
     s               -wqd(ig,k)+wqd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      subroutine dvthermcell2(ngrid,nlay,ptimestep,fm,entr,masse
     .    ,fraca,larga
     .    ,u,v,du,dv,ua,va)
      use dimphy
      implicit none

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c   calcul du dq/dt une fois qu'on connait les ascendances
c
c=======================================================================

#include "dimensions.h"
cccc#include "dimphy.h"

      integer ngrid,nlay

      real ptimestep
      real masse(ngrid,nlay),fm(ngrid,nlay+1)
      real fraca(ngrid,nlay+1)
      real larga(ngrid)
      real entr(ngrid,nlay)
      real u(ngrid,nlay)
      real ua(ngrid,nlay)
      real du(ngrid,nlay)
      real v(ngrid,nlay)
      real va(ngrid,nlay)
      real dv(ngrid,nlay)

      real qa(klon,klev),detr(klon,klev),zf,zf2
      real wvd(klon,klev+1),wud(klon,klev+1)
      real gamma0,gamma(klon,klev+1)
      real ue(klon,klev),ve(klon,klev)
      real dua,dva
      integer iter

      integer ig,k

c   calcul du detrainement

      do k=1,nlay
         do ig=1,ngrid
            detr(ig,k)=fm(ig,k)-fm(ig,k+1)+entr(ig,k)
         enddo
      enddo

c   calcul de la valeur dans les ascendances
      do ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
         ue(ig,1)=u(ig,1)
         ve(ig,1)=v(ig,1)
      enddo

      do k=2,nlay
         do ig=1,ngrid
            if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.
     s         1.e-5*masse(ig,k)) then
c   On itère sur la valeur du coeff de freinage.
c              gamma0=rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k))
               gamma0=masse(ig,k)
     s         *sqrt( 0.5*(fraca(ig,k+1)+fraca(ig,k)) )
     s         *0.5/larga(ig)
     s         *1.
c    s         *0.5
c              gamma0=0.
               zf=0.5*(fraca(ig,k)+fraca(ig,k+1))
               zf=0.
               zf2=1./(1.-zf)
c   la première fois on multiplie le coefficient de freinage
c   par le module du vent dans la couche en dessous.
               dua=ua(ig,k-1)-u(ig,k-1)
               dva=va(ig,k-1)-v(ig,k-1)
               do iter=1,5
c   On choisit une relaxation lineaire.
                  gamma(ig,k)=gamma0
c   On choisit une relaxation quadratique.
                  gamma(ig,k)=gamma0*sqrt(dua**2+dva**2)
                  ua(ig,k)=(fm(ig,k)*ua(ig,k-1)
     s               +(zf2*entr(ig,k)+gamma(ig,k))*u(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2
     s                 +gamma(ig,k))
                  va(ig,k)=(fm(ig,k)*va(ig,k-1)
     s               +(zf2*entr(ig,k)+gamma(ig,k))*v(ig,k))
     s               /(fm(ig,k+1)+detr(ig,k)+entr(ig,k)*zf*zf2
     s                 +gamma(ig,k))
c                 print*,k,ua(ig,k),va(ig,k),u(ig,k),v(ig,k),dua,dva
                  dua=ua(ig,k)-u(ig,k)
                  dva=va(ig,k)-v(ig,k)
                  ue(ig,k)=(u(ig,k)-zf*ua(ig,k))*zf2
                  ve(ig,k)=(v(ig,k)-zf*va(ig,k))*zf2
               enddo
            else
               ua(ig,k)=u(ig,k)
               va(ig,k)=v(ig,k)
               ue(ig,k)=u(ig,k)
               ve(ig,k)=v(ig,k)
               gamma(ig,k)=0.
            endif
         enddo
      enddo

      do k=2,nlay
         do ig=1,ngrid
            wud(ig,k)=fm(ig,k)*ue(ig,k)
            wvd(ig,k)=fm(ig,k)*ve(ig,k)
         enddo
      enddo
      do ig=1,ngrid
         wud(ig,1)=0.
         wud(ig,nlay+1)=0.
         wvd(ig,1)=0.
         wvd(ig,nlay+1)=0.
      enddo

      do k=1,nlay
         do ig=1,ngrid
            du(ig,k)=((detr(ig,k)+gamma(ig,k))*ua(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*ue(ig,k)
     s               -wud(ig,k)+wud(ig,k+1))
     s               /masse(ig,k)
            dv(ig,k)=((detr(ig,k)+gamma(ig,k))*va(ig,k)
     s               -(entr(ig,k)+gamma(ig,k))*ve(ig,k)
     s               -wvd(ig,k)+wvd(ig,k+1))
     s               /masse(ig,k)
         enddo
      enddo

      return
      end
      SUBROUTINE thermcell_sec(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi,zlev
     s                  ,pu,pv,pt,po
     s                  ,pduadj,pdvadj,pdtadj,pdoadj
     s                  ,fm0,entr0
c    s                  ,pu_therm,pv_therm
     s                  ,r_aspect,l_mix,w2di,tho)

      use dimphy
      IMPLICIT NONE

c=======================================================================
c
c   Calcul du transport verticale dans la couche limite en presence
c   de "thermiques" explicitement representes
c
c   Réécriture à partir d'un listing papier à Habas, le 14/02/00
c
c   le thermique est supposé homogène et dissipé par mélange avec
c   son environnement. la longueur l_mix contrôle l'efficacité du
c   mélange
c
c   Le calcul du transport des différentes espèces se fait en prenant
c   en compte:
c     1. un flux de masse montant
c     2. un flux de masse descendant
c     3. un entrainement
c     4. un detrainement
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"

c   arguments:
c   ----------

      INTEGER ngrid,nlay,w2di,tho
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

      integer idetr
      save idetr
      data idetr/3/
c$OMP THREADPRIVATE(idetr)

c   local:
c   ------

      INTEGER ig,k,l,lmaxa(klon),lmix(klon)
      real zsortie1d(klon)
c CR: on remplace lmax(klon,klev+1)
      INTEGER lmax(klon),lmin(klon),lentr(klon)
      real linter(klon)
      real zmix(klon), fracazmix(klon) 
c RC 
      real zmax(klon),zw,zz,zw2(klon,klev+1),ztva(klon,klev),zzz

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
      real zwd(klon,klev+1)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zha(klon,klev)
      real wa_moy(klon,klev+1)
      real fraca(klon,klev+1)
      real fracc(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev)
!      common/comtherm/thetath2,wth2

      real count_time
      integer isplit,nsplit,ialt
      parameter (nsplit=10)
      data isplit/0/
      save isplit
c$OMP THREADPRIVATE(isplit)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev+1),masse(klon,klev)
      real zpspsk(klon,klev)

c     real wmax(klon,klev),wmaxa(klon)
      real wmax(klon),wmaxa(klon)
      real wa(klon,klev,klev+1)
      real wd(klon,klev+1)
      real larg_part(klon,klev,klev+1)
      real fracd(klon,klev+1)
      real xxx(klon,klev+1)
      real larg_cons(klon,klev+1)
      real larg_detr(klon,klev+1)
      real fm0(klon,klev+1),entr0(klon,klev),detr(klon,klev)
      real pu_therm(klon,klev),pv_therm(klon,klev)
      real fm(klon,klev+1),entr(klon,klev)
      real fmc(klon,klev+1)

cCR:nouvelles variables
      real f_star(klon,klev+1),entr_star(klon,klev)
      real entr_star_tot(klon),entr_star2(klon)
      real f(klon), f0(klon)
      real zlevinter(klon)
      logical first
      data first /.false./
      save first
c$OMP THREADPRIVATE(first)
cRC

      character*2 str2
      character*10 str10

      character (len=20) :: modname='thermcell_sec'
      character (len=80) :: abort_message

      LOGICAL vtest(klon),down

      EXTERNAL SCOPY

      integer ncorrec,ll
      save ncorrec
      data ncorrec/0/
c$OMP THREADPRIVATE(ncorrec)
      
c
c-----------------------------------------------------------------------
c   initialisation:
c   ---------------
c
       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
c
c-----------------------------------------------------------------------
c   incrementation eventuelle de tendances precedentes:
c   ---------------------------------------------------

c       print*,'0 OK convect8'

      DO 1010 l=1,nlay
         DO 1015 ig=1,ngrid
            zpspsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**RKAPPA
            zh(ig,l)=pt(ig,l)/zpspsk(ig,l)
            zu(ig,l)=pu(ig,l)
            zv(ig,l)=pv(ig,l)
            zo(ig,l)=po(ig,l)
            ztv(ig,l)=zh(ig,l)*(1.+0.61*zo(ig,l))
1015     CONTINUE
1010  CONTINUE

c       print*,'1 OK convect8'
c                       --------------------
c
c
c                       + + + + + + + + + + +
c
c
c  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
c  wh,wt,wo ...
c
c                       + + + + + + + + + + +  zh,zu,zv,zo,rho
c
c
c                       --------------------   zlev(1)
c                       \\\\\\\\\\\\\\\\\\\\
c
c

c-----------------------------------------------------------------------
c   Calcul des altitudes des couches
c-----------------------------------------------------------------------

      do l=2,nlay
         do ig=1,ngrid
            zlev(ig,l)=0.5*(pphi(ig,l)+pphi(ig,l-1))/RG
         enddo
      enddo
      do ig=1,ngrid
         zlev(ig,1)=0.
         zlev(ig,nlay+1)=(2.*pphi(ig,klev)-pphi(ig,klev-1))/RG
      enddo
      do l=1,nlay
         do ig=1,ngrid
            zlay(ig,l)=pphi(ig,l)/RG
         enddo
      enddo

c      print*,'2 OK convect8'
c-----------------------------------------------------------------------
c   Calcul des densites
c-----------------------------------------------------------------------

      do l=1,nlay
         do ig=1,ngrid
            rho(ig,l)=pplay(ig,l)/(zpspsk(ig,l)*RD*zh(ig,l))
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            rhobarz(ig,l)=0.5*(rho(ig,l)+rho(ig,l-1))
         enddo
      enddo

      do k=1,nlay
         do l=1,nlay+1
            do ig=1,ngrid
               wa(ig,k,l)=0.
            enddo
         enddo
      enddo

c      print*,'3 OK convect8'
c------------------------------------------------------------------
c   Calcul de w2, quarre de w a partir de la cape
c   a partir de w2, on calcule wa, vitesse de l'ascendance
c
c   ATTENTION: Dans cette version, pour cause d'economie de memoire,
c   w2 est stoke dans wa
c
c   ATTENTION: dans convect8, on n'utilise le calcule des wa
c   independants par couches que pour calculer l'entrainement
c   a la base et la hauteur max de l'ascendance.
c
c   Indicages:
c   l'ascendance provenant du niveau k traverse l'interface l avec
c   une vitesse wa(k,l).
c
c                       --------------------
c
c                       + + + + + + + + + + 
c
c  wa(k,l)   ----       --------------------    l
c             /\
c            /||\       + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||
c             ||        + + + + + + + + + + 
c             ||
c             ||        --------------------
c             ||__
c             |___      + + + + + + + + + +     k
c
c                       --------------------
c
c
c
c------------------------------------------------------------------

cCR: ponderation entrainement des couches instables
cdef des entr_star tels que entr=f*entr_star      
      do l=1,klev
         do ig=1,ngrid 
            entr_star(ig,l)=0.
         enddo
      enddo
c determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lentr(ig)=1
      enddo

con ne considere que les premieres couches instables
      do k=nlay-2,1,-1
         do ig=1,ngrid
            if (ztv(ig,k).gt.ztv(ig,k+1).and.
     s          ztv(ig,k+1).le.ztv(ig,k+2)) then
               lentr(ig)=k
            endif
          enddo
      enddo
    
c determination du lmin: couche d ou provient le thermique
      do ig=1,ngrid
         lmin(ig)=1
      enddo
      do ig=1,ngrid
         do l=nlay,2,-1
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.le.lentr(ig)) then 
                 entr_star(ig,l)=(ztv(ig,l)-ztv(ig,l+1))*
c     s                           (zlev(ig,l+1)-zlev(ig,l))
     s                           *sqrt(zlev(ig,l+1))
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
            do l=1,klev
               entr_star(ig,l)=0.
            enddo
         endif
      enddo 
c calcul de l entrainement total
      do ig=1,ngrid
         entr_star_tot(ig)=0.
      enddo
      do ig=1,ngrid
         do k=1,klev
            entr_star_tot(ig)=entr_star_tot(ig)+entr_star(ig,k)
         enddo
      enddo
c
c      print*,'fin calcul entr_star'
      do k=1,klev
         do ig=1,ngrid 
            ztva(ig,k)=ztv(ig,k)
         enddo
      enddo
cRC
c      print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            fmc(ig,k)=0.
cCR
            f_star(ig,k)=0.
cRC
            larg_cons(ig,k)=0.
            larg_detr(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

c      print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmaxa(ig)=1
         lmix(ig)=1
         wmaxa(ig)=0.
      enddo

cCR: 
      do l=1,nlay-2
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.entr_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
               f_star(ig,l+1)=entr_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               larg_detr(ig,l)=0.
            else if ((zw2(ig,l).ge.1e-10).and.
     s               (f_star(ig,l)+entr_star(ig,l).gt.1.e-10)) then
               f_star(ig,l+1)=f_star(ig,l)+entr_star(ig,l)
               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+entr_star(ig,l)
     s                    *ztv(ig,l))/f_star(ig,l+1)
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+
     s                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
c determination de zmax continu par interpolation lineaire
            if (zw2(ig,l+1).lt.0.) then
ctest
               if (abs(zw2(ig,l+1)-zw2(ig,l)).lt.1e-10) then
c                  print*,'pb linter'
               endif
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmaxa(ig)=l
            else
               if (zw2(ig,l+1).lt.0.) then
c                  print*,'pb1 zw2<0'
               endif
               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))
            endif
            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
c   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
c      print*,'fin calcul zw2'
c
c Calcul de la couche correspondant a la hauteur du thermique
      do ig=1,ngrid
         lmax(ig)=lentr(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.1) then
            lmax(ig)=1
            lmin(ig)=1
         endif
      enddo 
c    
c Determination de zw2 max
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
c                  print*,'pb2 zw2<0'
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

c   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
c calcul de zlevinter
          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*
     s    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)
     s    -zlev(ig,lmax(ig)))
       zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

c      print*,'avant fermeture'
c Fermeture,determination de f
      do ig=1,ngrid
         entr_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (entr_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else
             do k=lmin(ig),lentr(ig)
                entr_star2(ig)=entr_star2(ig)+entr_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
c Nouvelle fermeture
             f(ig)=wmax(ig)/(max(500.,zmax(ig))*r_aspect
     s             *entr_star2(ig))*entr_star_tot(ig)
ctest
c             if (first) then
c             f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
c     s             *wmax(ig))
c             endif
         endif
c         f0(ig)=f(ig)
c         first=.true.
      enddo
c      print*,'apres fermeture'

c Calcul de l'entrainement
       do k=1,klev
         do ig=1,ngrid 
            entr(ig,k)=f(ig)*entr_star(ig,k)
         enddo
      enddo
cCR:test pour entrainer moins que la masse
       do ig=1,ngrid
          do l=1,lentr(ig)
             if ((entr(ig,l)*ptimestep).gt.(0.9*masse(ig,l))) then
                entr(ig,l+1)=entr(ig,l+1)+entr(ig,l)
     s                       -0.9*masse(ig,l)/ptimestep
                entr(ig,l)=0.9*masse(ig,l)/ptimestep
             endif
          enddo
       enddo
cCR: fin test
c Calcul des flux
      do ig=1,ngrid
         do l=1,lmax(ig)-1
            fmc(ig,l+1)=fmc(ig,l)+entr(ig,l)
         enddo
      enddo

cRC


c      print*,'9 OK convect8'
c     print*,'WA1 ',wa_moy

c   determination de l'indice du debut de la mixed layer ou w decroit

c   calcul de la largeur de chaque ascendance dans le cas conservatif.
c   dans ce cas simple, on suppose que la largeur de l'ascendance provenant
c   d'une couche est égale à la hauteur de la couche alimentante.
c   La vitesse maximale dans l'ascendance est aussi prise comme estimation
c   de la vitesse d'entrainement horizontal dans la couche alimentante.

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
               zw=max(wa_moy(ig,l),1.e-10)
               larg_cons(ig,l)=zmax(ig)*r_aspect
     s         *fmc(ig,l)/(rhobarz(ig,l)*zw)
            endif
         enddo
      enddo

      do l=2,nlay
         do ig=1,ngrid
            if (l.le.lmaxa(ig)) then
c              if (idetr.eq.0) then
c  cette option est finalement en dur.
                  if ((l_mix*zlev(ig,l)).lt.0.)then
c                   print*,'pb l_mix*zlev<0'
                  endif
cCR: test: nouvelle def de lambda
c                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
                  if (zw2(ig,l).gt.1.e-10) then
                  larg_detr(ig,l)=sqrt((l_mix/zw2(ig,l))*zlev(ig,l))
                  else
                  larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
                  endif
cRC
c              else if (idetr.eq.1) then
c                 larg_detr(ig,l)=larg_cons(ig,l)
c    s            *sqrt(l_mix*zlev(ig,l))/larg_cons(ig,lmix(ig))
c              else if (idetr.eq.2) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *sqrt(wa_moy(ig,l))
c              else if (idetr.eq.4) then
c                 larg_detr(ig,l)=sqrt(l_mix*zlev(ig,l))
c    s            *wa_moy(ig,l)
c              endif
            endif
         enddo
       enddo

c      print*,'10 OK convect8'
c     print*,'WA2 ',wa_moy
c   calcul de la fraction de la maille concernée par l'ascendance en tenant
c   compte de l'epluchage du thermique.
c
cCR def de  zmix continu (profil parabolique des vitesses)
      do ig=1,ngrid
           if (lmix(ig).gt.1.) then
c test 
              if (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))).gt.1e-10)
     s        then
c             
            zmix(ig)=((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))
     s        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))
     s        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))
     s        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))
     s        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
            else
            zmix(ig)=zlev(ig,lmix(ig))
c            print*,'pb zmix'
            endif
         else 
         zmix(ig)=0.
         endif
ctest
         if ((zmax(ig)-zmix(ig)).lt.0.) then
            zmix(ig)=0.99*zmax(ig)
c            print*,'pb zmix>zmax'
         endif
      enddo
c
c calcul du nouveau lmix correspondant
      do ig=1,ngrid
         do l=1,klev
            if (zmix(ig).ge.zlev(ig,l).and.
     s          zmix(ig).lt.zlev(ig,l+1)) then
              lmix(ig)=l
             endif
          enddo
      enddo
c
      do l=2,nlay
         do ig=1,ngrid
            if(larg_cons(ig,l).gt.1.) then
c     print*,ig,l,lmix(ig),lmaxa(ig),larg_cons(ig,l),'  KKK'
               fraca(ig,l)=(larg_cons(ig,l)-larg_detr(ig,l))
     s            /(r_aspect*zmax(ig))
c test
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
            else
c              wa_moy(ig,l)=0.
               fraca(ig,l)=0.
               fracc(ig,l)=0.
               fracd(ig,l)=1.
            endif
         enddo
      enddo                  
cCR: calcul de fracazmix
       do ig=1,ngrid
          fracazmix(ig)=(fraca(ig,lmix(ig)+1)-fraca(ig,lmix(ig)))/
     s     (zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))*zmix(ig)
     s    +fraca(ig,lmix(ig))-zlev(ig,lmix(ig))*(fraca(ig,lmix(ig)+1)
     s    -fraca(ig,lmix(ig)))/(zlev(ig,lmix(ig)+1)-zlev(ig,lmix(ig)))
       enddo
c
       do l=2,nlay
          do ig=1,ngrid
             if(larg_cons(ig,l).gt.1.) then
               if (l.gt.lmix(ig)) then
ctest
                 if (zmax(ig)-zmix(ig).lt.1.e-10) then
c                   print*,'pb xxx'
                   xxx(ig,l)=(lmaxa(ig)+1.-l)/(lmaxa(ig)+1.-lmix(ig))
                 else
                 xxx(ig,l)=(zmax(ig)-zlev(ig,l))/(zmax(ig)-zmix(ig))
                 endif
           if (idetr.eq.0) then
               fraca(ig,l)=fracazmix(ig)
           else if (idetr.eq.1) then
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)
           else if (idetr.eq.2) then
               fraca(ig,l)=fracazmix(ig)*(1.-(1.-xxx(ig,l))**2)
           else
               fraca(ig,l)=fracazmix(ig)*xxx(ig,l)**2
           endif
c     print*,ig,l,lmix(ig),lmaxa(ig),xxx(ig,l),'LLLLLLL'
               fraca(ig,l)=max(fraca(ig,l),0.)
               fraca(ig,l)=min(fraca(ig,l),0.5)
               fracd(ig,l)=1.-fraca(ig,l)
               fracc(ig,l)=larg_cons(ig,l)/(r_aspect*zmax(ig))
             endif
            endif
         enddo
      enddo
      
c      print*,'fin calcul fraca'
c      print*,'11 OK convect8'
c     print*,'Ea3 ',wa_moy
c------------------------------------------------------------------
c   Calcul de fracd, wd
c   somme wa - wd = 0
c------------------------------------------------------------------


      do ig=1,ngrid
         fm(ig,1)=0.
         fm(ig,nlay+1)=0.
      enddo

      do l=2,nlay
           do ig=1,ngrid
              fm(ig,l)=fraca(ig,l)*wa_moy(ig,l)*rhobarz(ig,l)
cCR:test
              if (entr(ig,l-1).lt.1e-10.and.fm(ig,l).gt.fm(ig,l-1)
     s            .and.l.gt.lmix(ig)) then
                 fm(ig,l)=fm(ig,l-1)
c                 write(1,*)'ajustement fm, l',l
              endif
c              write(1,*)'ig,l,fm(ig,l)',ig,l,fm(ig,l)
cRC
           enddo
         do ig=1,ngrid
            if(fracd(ig,l).lt.0.1) then
              abort_message = 'fracd trop petit'
              CALL abort_gcm (modname,abort_message,1)
            else
c    vitesse descendante "diagnostique"
               wd(ig,l)=fm(ig,l)/(fracd(ig,l)*rhobarz(ig,l))
            endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
c           masse(ig,l)=rho(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
            masse(ig,l)=(pplev(ig,l)-pplev(ig,l+1))/RG
         enddo
      enddo

c      print*,'12 OK convect8'
c     print*,'WA4 ',wa_moy
cc------------------------------------------------------------------
c   calcul du transport vertical
c------------------------------------------------------------------

      go to 4444
c     print*,'XXXXXXXXXXXXXXX ptimestep= ',ptimestep
      do l=2,nlay-1
         do ig=1,ngrid
            if(fm(ig,l+1)*ptimestep.gt.masse(ig,l)
     s      .and.fm(ig,l+1)*ptimestep.gt.masse(ig,l+1)) then
c     print*,'WARN!!! FM>M ig=',ig,' l=',l,'  FM='
c    s         ,fm(ig,l+1)*ptimestep
c    s         ,'   M=',masse(ig,l),masse(ig,l+1)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(entr(ig,l)*ptimestep.gt.masse(ig,l)) then
c     print*,'WARN!!! E>M ig=',ig,' l=',l,'  E=='
c    s         ,entr(ig,l)*ptimestep
c    s         ,'   M=',masse(ig,l)
             endif
         enddo
      enddo

      do l=1,nlay
         do ig=1,ngrid
            if(.not.fm(ig,l).ge.0..or..not.fm(ig,l).le.10.) then
c     print*,'WARN!!! fm exagere ig=',ig,'   l=',l
c    s         ,'   FM=',fm(ig,l)
            endif
            if(.not.masse(ig,l).ge.1.e-10
     s         .or..not.masse(ig,l).le.1.e4) then
c     print*,'WARN!!! masse exagere ig=',ig,'   l=',l
c    s         ,'   M=',masse(ig,l)
c     print*,'rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)',
c    s                 rho(ig,l),pplay(ig,l),zpspsk(ig,l),RD,zh(ig,l)
c     print*,'zlev(ig,l+1),zlev(ig,l)'
c    s                ,zlev(ig,l+1),zlev(ig,l)
c     print*,'pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)'
c    s                ,pphi(ig,l-1),pphi(ig,l),pphi(ig,l+1)
            endif
            if(.not.entr(ig,l).ge.0..or..not.entr(ig,l).le.10.) then
c     print*,'WARN!!! entr exagere ig=',ig,'   l=',l
c    s         ,'   E=',entr(ig,l)
            endif
         enddo
      enddo

4444   continue

cCR:redefinition du entr
       do l=1,nlay
         do ig=1,ngrid
            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
            if (detr(ig,l).lt.0.) then
                entr(ig,l)=entr(ig,l)-detr(ig,l)
                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
            endif
         enddo
      enddo
cRC
      if (w2di.eq.1) then
         fm0=fm0+ptimestep*(fm-fm0)/REAL(tho)
         entr0=entr0+ptimestep*(entr-entr0)/REAL(tho)
      else
         fm0=fm
         entr0=entr
      endif

      if (1.eq.1) then
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zh,zdhadj,zha)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zo,pdoadj,zoa)
      else
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zh,zdhadj,zha)
         call dqthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse,fraca
     .    ,zo,pdoadj,zoa)
      endif

      if (1.eq.0) then
         call dvthermcell2(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,fraca,zmax
     .    ,zu,zv,pduadj,pdvadj,zua,zva)
      else
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zu,pduadj,zua)
         call dqthermcell(ngrid,nlay,ptimestep,fm0,entr0,masse
     .    ,zv,pdvadj,zva)
      endif

      do l=1,nlay
         do ig=1,ngrid
            zf=0.5*(fracc(ig,l)+fracc(ig,l+1))
            zf2=zf/(1.-zf)
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l))**2
            wth2(ig,l)=zf2*(0.5*(wa_moy(ig,l)+wa_moy(ig,l+1)))**2
         enddo
      enddo



c     print*,'13 OK convect8'
c     print*,'WA5 ',wa_moy
      do l=1,nlay
         do ig=1,ngrid
            pdtadj(ig,l)=zdhadj(ig,l)*zpspsk(ig,l)
         enddo
      enddo


c     do l=1,nlay
c        do ig=1,ngrid
c           if(abs(pdtadj(ig,l))*86400..gt.500.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdtadj=',pdtadj(ig,l)
c           endif
c           if(abs(pdoadj(ig,l))*86400..gt.1.) then
c     print*,'WARN!!! ig=',ig,'  l=',l
c    s         ,'   pdoadj=',pdoadj(ig,l)
c           endif
c        enddo
c      enddo

c      print*,'14 OK convect8'
c------------------------------------------------------------------
c   Calculs pour les sorties
c------------------------------------------------------------------

      if(sorties) then
      do l=1,nlay
         do ig=1,ngrid
            zla(ig,l)=(1.-fracd(ig,l))*zmax(ig)
            zld(ig,l)=fracd(ig,l)*zmax(ig)
            if(1.-fracd(ig,l).gt.1.e-10)
     s      zwa(ig,l)=wd(ig,l)*fracd(ig,l)/(1.-fracd(ig,l))
         enddo
      enddo

cdeja fait
c      do l=1,nlay
c         do ig=1,ngrid
c            detr(ig,l)=fm(ig,l)+entr(ig,l)-fm(ig,l+1)
c            if (detr(ig,l).lt.0.) then
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
c                detr(ig,l)=0.
c     print*,'WARNING !!! detrainement negatif ',ig,l
c            endif
c         enddo
c      enddo

c     print*,'15 OK convect8'

      isplit=isplit+1


c #define und
	goto 123
#ifdef und
      CALL writeg1d(1,nlay,wd,'wd      ','wd      ')
      CALL writeg1d(1,nlay,zwa,'wa      ','wa      ')
      CALL writeg1d(1,nlay,fracd,'fracd      ','fracd      ')
      CALL writeg1d(1,nlay,fraca,'fraca      ','fraca      ')
      CALL writeg1d(1,nlay,wa_moy,'wam         ','wam         ')
      CALL writeg1d(1,nlay,zla,'la      ','la      ')
      CALL writeg1d(1,nlay,zld,'ld      ','ld      ')
      CALL writeg1d(1,nlay,pt,'pt      ','pt      ')
      CALL writeg1d(1,nlay,zh,'zh      ','zh      ')
      CALL writeg1d(1,nlay,zha,'zha      ','zha      ')
      CALL writeg1d(1,nlay,zu,'zu      ','zu      ')
      CALL writeg1d(1,nlay,zv,'zv      ','zv      ')
      CALL writeg1d(1,nlay,zo,'zo      ','zo      ')
      CALL writeg1d(1,nlay,wh,'wh      ','wh      ')
      CALL writeg1d(1,nlay,wu,'wu      ','wu      ')
      CALL writeg1d(1,nlay,wv,'wv      ','wv      ')
      CALL writeg1d(1,nlay,wo,'w15uo     ','wXo     ')
      CALL writeg1d(1,nlay,zdhadj,'zdhadj      ','zdhadj      ')
      CALL writeg1d(1,nlay,pduadj,'pduadj      ','pduadj      ')
      CALL writeg1d(1,nlay,pdvadj,'pdvadj      ','pdvadj      ')
      CALL writeg1d(1,nlay,pdoadj,'pdoadj      ','pdoadj      ')
      CALL writeg1d(1,nlay,entr  ,'entr        ','entr        ')
      CALL writeg1d(1,nlay,detr  ,'detr        ','detr        ')
      CALL writeg1d(1,nlay,fm    ,'fm          ','fm          ')

      CALL writeg1d(1,nlay,pdtadj,'pdtadj    ','pdtadj    ')
      CALL writeg1d(1,nlay,pplay,'pplay     ','pplay     ')
      CALL writeg1d(1,nlay,pplev,'pplev     ','pplev     ')

c   recalcul des flux en diagnostique...
c     print*,'PAS DE TEMPS ',ptimestep
       call dt2F(pplev,pplay,pt,pdtadj,wh)
      CALL writeg1d(1,nlay,wh,'wh2     ','wh2     ')
#endif
123   continue
! #define troisD
#ifdef troisD
c       if (sorties) then
      print*,'Debut des wrgradsfi'

c      print*,'16 OK convect8'
         call wrgradsfi(1,nlay,wd,'wd        ','wd        ')
         call wrgradsfi(1,nlay,zwa,'wa        ','wa        ')
         call wrgradsfi(1,nlay,fracd,'fracd     ','fracd     ')
         call wrgradsfi(1,nlay,fraca,'fraca     ','fraca     ')
         call wrgradsfi(1,nlay,xxx,'xxx       ','xxx       ')
         call wrgradsfi(1,nlay,wa_moy,'wam       ','wam       ')
c      print*,'WA6 ',wa_moy
         call wrgradsfi(1,nlay,zla,'la        ','la        ')
         call wrgradsfi(1,nlay,zld,'ld        ','ld        ')
         call wrgradsfi(1,nlay,pt,'pt        ','pt        ')
         call wrgradsfi(1,nlay,zh,'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha,'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua,'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva,'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu,'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv,'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,wh,'wh        ','wh        ')
         call wrgradsfi(1,nlay,wu,'wu        ','wu        ')
         call wrgradsfi(1,nlay,wv,'wv        ','wv        ')
         call wrgradsfi(1,nlay,wo,'wo        ','wo        ')
         call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
         call wrgradsfi(1,nlay,zdhadj,'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj,'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj,'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj,'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr,'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr,'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm,'fm        ','fm        ')
         call wrgradsfi(1,nlay,fmc,'fmc       ','fmc       ')
         call wrgradsfi(1,nlay,zw2,'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,ztva,'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv,'ztv       ','ztv       ')

         call wrgradsfi(1,nlay,zo,'zo        ','zo        ')
         call wrgradsfi(1,nlay,larg_cons,'Lc        ','Lc        ')
         call wrgradsfi(1,nlay,larg_detr,'Ldetr     ','Ldetr     ')

cCR:nouveaux diagnostiques
      call wrgradsfi(1,nlay,entr_star  ,'entr_star   ','entr_star   ')     
      call wrgradsfi(1,nlay,f_star    ,'f_star   ','f_star   ')
      call wrgradsfi(1,1,zmax,'zmax      ','zmax      ')
      call wrgradsfi(1,1,zmix,'zmix      ','zmix      ') 
      zsortie1d(:)=lmax(:)
      call wrgradsfi(1,1,zsortie1d,'lmax      ','lmax      ')
      call wrgradsfi(1,1,wmax,'wmax      ','wmax      ')
      zsortie1d(:)=lmix(:)
      call wrgradsfi(1,1,zsortie1d,'lmix      ','lmix      ')
      zsortie1d(:)=lentr(:)
      call wrgradsfi(1,1,zsortie1d,'lentr      ','lentr     ')

c      print*,'17 OK convect8'

         do k=1,klev/10
            write(str2,'(i2.2)') k
            str10='wa'//str2
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=wa(ig,k,l)
               enddo
            enddo   
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
            do l=1,nlay
               do ig=1,ngrid
                  zsortie(ig,l)=larg_part(ig,k,l)
               enddo
            enddo
            str10='la'//str2
            CALL wrgradsfi(1,nlay,zsortie,str10,str10)
         enddo


c     print*,'18 OK convect8'
c      endif
      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

c      print*,'19 OK convect8'
      return
      end

