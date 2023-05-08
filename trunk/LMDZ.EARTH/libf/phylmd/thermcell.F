!
! $Id: thermcell.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE calcul_sec(ngrid,nlay,ptimestep
     s                  ,pplay,pplev,pphi,zlev
     s                  ,pu,pv,pt,po
     s                  ,zmax,wmax,zw2,lmix
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
      real zmax(klon),zw,zw2(klon,klev+1),ztva(klon,klev)

      real zlev(klon,klev+1),zlay(klon,klev)
      REAL zh(klon,klev),zdhadj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      REAL wh(klon,klev+1)
      real wu(klon,klev+1),wv(klon,klev+1),wo(klon,klev+1)
      real zla(klon,klev+1)
      real zwa(klon,klev+1)
      real zld(klon,klev+1)
!      real zwd(klon,klev+1)
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
!      integer isplit,nsplit
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
      real zalim(klon)
      integer lalim(klon)
      real norme(klon)
      real f(klon), f0(klon)
      real zlevinter(klon)
      logical therm
      logical first
      data first /.false./
      save first
c$OMP THREADPRIVATE(first)
cRC

      character*2 str2
      character*10 str10

      character (len=20) :: modname='calcul_sec'
      character (len=80) :: abort_message


!      LOGICAL vtest(klon),down

      EXTERNAL SCOPY

      integer ncorrec
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
climitation de la valeur du lentr
c      do ig=1,ngrid
c         lentr(ig)=min(5,lentr(ig))
c      enddo
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
cinitialisations
      do ig=1,ngrid
         zalim(ig)=0.
         norme(ig)=0.
         lalim(ig)=1
      enddo
      do k=1,klev-1
         do ig=1,ngrid
       zalim(ig)=zalim(ig)+zlev(ig,k)*MAX(0.,(ztv(ig,k)-ztv(ig,k+1))
     s          /(zlev(ig,k+1)-zlev(ig,k)))
c     s         *(zlev(ig,k+1)-zlev(ig,k))
       norme(ig)=norme(ig)+MAX(0.,(ztv(ig,k)-ztv(ig,k+1))
     s          /(zlev(ig,k+1)-zlev(ig,k)))
c    s          *(zlev(ig,k+1)-zlev(ig,k))
         enddo
       enddo
       do ig=1,ngrid
          if (norme(ig).gt.1.e-10) then
             zalim(ig)=max(10.*zalim(ig)/norme(ig),zlev(ig,2))
c             zalim(ig)=min(zalim(ig),zlev(ig,lentr(ig)))
          endif
       enddo
cdétermination du lalim correspondant
      do k=1,klev-1
         do ig=1,ngrid
      if ((zalim(ig).gt.zlev(ig,k)).and.(zalim(ig).le.zlev(ig,k+1))) 
     s   then
         lalim(ig)=k
      endif
         enddo
      enddo
c
c definition de l'entrainement des couches
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.lt.lentr(ig)) then 
                 entr_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)
c     s                           *(zlev(ig,l+1)-zlev(ig,l))
     s                           *sqrt(zlev(ig,l+1))
cautre def
c                entr_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
c     s                         /zlev(ig,lentr(ig)+2)))**(3./2.)
            endif
         enddo
      enddo
cnouveau test
c      if (therm) then
      do l=1,klev-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.
     s          l.ge.lmin(ig).and.l.le.lalim(ig)
     s          .and.zalim(ig).gt.1.e-10) then 
c            if (l.le.lentr(ig)) then 
c               entr_star(ig,l)=zlev(ig,l+1)*(1.-(zlev(ig,l+1)
c     s                         /zalim(ig)))**(3./2.)
c               write(10,*)zlev(ig,l),entr_star(ig,l)
            endif
         enddo
      enddo
c      endif
c pas de thermique si couche 1 stable
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
c Calcul entrainement normalise
      do ig=1,ngrid 
         if (entr_star_tot(ig).gt.1.e-10) then
c         do l=1,lentr(ig)
          do l=1,klev
cdef possibles pour entr_star: zdthetadz, dthetadz, zdtheta 
            entr_star(ig,l)=entr_star(ig,l)/entr_star_tot(ig)
         enddo
         endif
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
c         lmax(ig)=lalim(ig)
      enddo
      do ig=1,ngrid
         do l=nlay,lentr(ig)+1,-1
c         do l=nlay,lalim(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1
            endif
         enddo
      enddo
c pas de thermique si couche 1 stable
      do ig=1,ngrid
         if (lmin(ig).gt.5) then
            lmax(ig)=1
            lmin(ig)=1
            lentr(ig)=1
            lalim(ig)=1
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
      do ig=1,ngrid
c      write(8,*)zmax(ig),lmax(ig),lentr(ig),lmin(ig)
      enddo
con stope après les calculs de zmax et wmax
      RETURN

c      print*,'avant fermeture'
c Fermeture,determination de f
cAttention! entrainement normalisé ou pas?
      do ig=1,ngrid
         entr_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (entr_star_tot(ig).LT.1.e-10) then
            f(ig)=0.
         else
             do k=lmin(ig),lentr(ig)
c             do k=lmin(ig),lalim(ig) 
                entr_star2(ig)=entr_star2(ig)+entr_star(ig,k)**2
     s                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
c Nouvelle fermeture
             f(ig)=wmax(ig)/(max(500.,zmax(ig))*r_aspect
     s             *entr_star2(ig))
c     s            *entr_star_tot(ig)
ctest
c             if (first) then
             f(ig)=f(ig)+(f0(ig)-f(ig))*exp(-ptimestep/zmax(ig)
     s             *wmax(ig))
c             endif
         endif
         f0(ig)=f(ig)
c         first=.true.
      enddo
c      print*,'apres fermeture'
con stoppe après la fermeture
      RETURN
c Calcul de l'entrainement
       do k=1,klev
         do ig=1,ngrid 
            entr(ig,k)=f(ig)*entr_star(ig,k)
         enddo
      enddo
con stoppe après le calcul de entr
c      RETURN
cCR:test pour entrainer moins que la masse
c       do ig=1,ngrid
c          do l=1,lentr(ig)
c             if ((entr(ig,l)*ptimestep).gt.(0.9*masse(ig,l))) then
c                entr(ig,l+1)=entr(ig,l+1)+entr(ig,l)
c     s                       -0.9*masse(ig,l)/ptimestep
c                entr(ig,l)=0.9*masse(ig,l)/ptimestep
c             endif
c          enddo
c       enddo
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
c                entr(ig,l)=entr(ig,l)-detr(ig,l)
                fm(ig,l+1)=fm(ig,l)+entr(ig,l)
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
#define troisD
#ifdef troisD
c       if (sorties) then
c      print*,'Debut des wrgradsfi'

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
c      print*,'Fin des wrgradsfi'
#endif

      endif

c     if(wa_moy(1,4).gt.1.e-10) stop

c      print*,'19 OK convect8'
      return
      end

      SUBROUTINE fermeture_seche(ngrid,nlay
     s                ,pplay,pplev,pphi,zlev,rhobarz,f0,zpspsk
     s                ,alim_star,zh,zo,lentr,lmin,nu_min,nu_max,r_aspect
     s                ,zmax,wmax)

      USE dimphy
      IMPLICIT NONE

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"

      INTEGER ngrid,nlay
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)
      real zlev(klon,klev+1)
      real alim_star(klon,klev)
      real f0(klon)
      integer lentr(klon)
      integer lmin(klon)
      real zmax(klon)
      real wmax(klon)
      real nu_min
      real nu_max
      real r_aspect
      real rhobarz(klon,klev+1)
      REAL zh(klon,klev)
      real zo(klon,klev)
      real zpspsk(klon,klev)

      integer ig,l

      real f_star(klon,klev+1)
      real detr_star(klon,klev)
      real entr_star(klon,klev)
      real zw2(klon,klev+1)
      real linter(klon)
      integer lmix(klon)
      integer lmax(klon)
      real zlevinter(klon)
      real wa_moy(klon,klev+1)
      real wmaxa(klon)
      REAL ztv(klon,klev)
      REAL ztva(klon,klev)
      real nu(klon,klev)
!      real zmax0_sec(klon)
!      save zmax0_sec
       REAL, SAVE, ALLOCATABLE :: zmax0_sec(:)
c$OMP THREADPRIVATE(zmax0_sec)
      logical, save :: first = .true.
c$OMP THREADPRIVATE(first)

      if (first) then
        allocate(zmax0_sec(klon))
        first=.false.
      endif

      do l=1,nlay
         do ig=1,ngrid
      ztv(ig,l)=zh(ig,l)/zpspsk(ig,l)
      ztv(ig,l)=ztv(ig,l)*(1.+RETV*zo(ig,l))
         enddo
      enddo
      do l=1,nlay-2
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1)
     s         .and.alim_star(ig,l).gt.1.e-10
     s         .and.zw2(ig,l).lt.1e-10) then
               f_star(ig,l+1)=alim_star(ig,l)
ctest:calcul de dteta
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
     s                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
            else if ((zw2(ig,l).ge.1e-10).and.
     s         (f_star(ig,l)+alim_star(ig,l)).gt.1.e-10) then
cestimation du detrainement a partir de la geometrie du pas precedent
ctests sur la definition du detr
             nu(ig,l)=(nu_min+nu_max)/2.
     s         *(1.-(nu_max-nu_min)/(nu_max+nu_min)
     s  *tanh((((ztva(ig,l-1)-ztv(ig,l))/ztv(ig,l))/0.0005)))
         
             detr_star(ig,l)=rhobarz(ig,l)
     s                      *sqrt(zw2(ig,l)) 
     s                       /(r_aspect*zmax0_sec(ig))*
c     s                       /(r_aspect*zmax0(ig))*
     s                      (sqrt(nu(ig,l)*zlev(ig,l+1)
     s                /sqrt(zw2(ig,l)))
     s                     -sqrt(nu(ig,l)*zlev(ig,l)
     s                /sqrt(zw2(ig,l))))
         detr_star(ig,l)=detr_star(ig,l)/f0(ig)
         if ((detr_star(ig,l)).gt.f_star(ig,l)) then
              detr_star(ig,l)=f_star(ig,l)
         endif
         entr_star(ig,l)=0.9*detr_star(ig,l)
             if ((l.lt.lentr(ig))) then
                 entr_star(ig,l)=0.
c                 detr_star(ig,l)=0.
             endif 
c           print*,'ok detr_star'
cprise en compte du detrainement dans le calcul du flux
             f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)
     s                      -detr_star(ig,l)
ctest sur le signe de f_star
       if ((f_star(ig,l+1)+detr_star(ig,l)).gt.1.e-10) then 
cAM on melange Tl et qt du thermique
          ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+(entr_star(ig,l)
     s                    +alim_star(ig,l))
     s                    *ztv(ig,l))/(f_star(ig,l+1)+detr_star(ig,l))
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)
     s                     /(f_star(ig,l+1)+detr_star(ig,l)))**2+
     s                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
     s                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
        endif
c
            if (zw2(ig,l+1).lt.0.) then 
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))
     s           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
c              print*,'linter=',linter(ig)
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
c     print*,'fin calcul zw2'
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
c                 print*,'pb2 zw2<0'
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
       zmax0_sec(ig)=zmax(ig)
      enddo
       
      RETURN
      END
