!
! $Id: bilan_dyn.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE bilan_dyn (dt_app,dt_cum,
     s  ps,masse,pk,flux_u,flux_v,teta,phi,ucov,vcov,
     s  ducovdyn,ducovdis,ducovspg,ducovphy)
c si besoin des traceurs:
c      SUBROUTINE bilan_dyn (ntrac,dt_app,dt_cum,
c     s  ps,masse,pk,flux_u,flux_v,teta,phi,ucov,vcov,trac,
c     s  ducovdyn,ducovdis,ducovspg,ducovphy)

c   AFAIRE
c   Prevoir en champ nq+1 le diagnostique de l'energie
c   en faisant Qzon=Cv T + L * ...
c             vQ..A=Cp T + L * ...

#ifdef CPP_IOIPSL
      USE IOIPSL
#endif

      USE control_mod, ONLY: planet_type 
      USE cpdet_mod, only: tpot2t
      USE comvert_mod, ONLY: ap,bp,presnivs
      USE comconst_mod, ONLY: rad,omeg,pi
      USE temps_mod, ONLY: annee_ref,day_ref

      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "iniprint.h"

c====================================================================
c
c   Sous-programme consacre à des diagnostics dynamiques de base
c
c 
c   De facon generale, les moyennes des scalaires Q sont ponderees par
c   la masse.
c
c   Les flux de masse sont eux simplement moyennes.
c
c====================================================================

c   Arguments :
c   ===========

c      integer ntrac
      real dt_app,dt_cum
      real ps(iip1,jjp1)
      real masse(iip1,jjp1,llm),pk(iip1,jjp1,llm)
      real flux_u(iip1,jjp1,llm)
      real flux_v(iip1,jjm,llm)
      real teta(iip1,jjp1,llm)
      real phi(iip1,jjp1,llm)
      real ucov(iip1,jjp1,llm)
      real vcov(iip1,jjm,llm)
c      real trac(iip1,jjp1,llm,ntrac)
c Tendances en m/s2 :
      real ducovdyn(iip1,jjp1,llm)
      real ducovdis(iip1,jjp1,llm)
      real ducovspg(iip1,jjp1,llm)
      real ducovphy(iip1,jjp1,llm)

c   Local :
c   =======

      integer icum,ncum
      save icum,ncum

      integer i,j,l,iQ,num
      real zz,zqy,zfactv(jjm,llm),zfactw(jjm,llm)
      character*2 strd2
      real ww

      logical first
      save first
      data first/.true./

      integer i_sortie
      save i_sortie
      data i_sortie/1/

      real time
      integer itau
      save time,itau
      data time,itau/0.,0/

! facteur = -1. pour Venus
      real    fact_geovenus
      save    fact_geovenus

c   variables dynamiques intermédiaires
c -----------------------------------
      REAL vcont(iip1,jjm,llm),ucont(iip1,jjp1,llm)
      REAL ang(iip1,jjp1,llm),unat(iip1,jjp1,llm)
      REAL massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)
      REAL vorpot(iip1,jjm,llm)
      REAL w(iip1,jjp1,llm),ecin(iip1,jjp1,llm),convm(iip1,jjp1,llm)
      real temp(iip1,jjp1,llm)
      real dudyn(iip1,jjp1,llm)
      real dudis(iip1,jjp1,llm)
      real duspg(iip1,jjp1,llm)
      real duphy(iip1,jjp1,llm)

c CHAMPS SCALAIRES Q ADVECTES
c ----------------------------
      integer nQ
c avec tous les composes, ca fait trop.... Je les enleve
c     parameter (nQ=6+nqmx)
      parameter (nQ=6)

      character*6,save :: nom(nQ)
      character*6,save :: unites(nQ)

      integer itemp,igeop,iecin,iang,iu,iun
      save itemp,igeop,iecin,iang,iu,iun
      data itemp,igeop,iecin,iang,iu,iun/1,2,3,4,5,6/

c   champ contenant les scalaires advectés.
      real Q(iip1,jjp1,llm,nQ)
    
c   champs cumulés
      real ps_cum(iip1,jjp1)
      real masse_cum(iip1,jjp1,llm)
      real flux_u_cum(iip1,jjp1,llm)
      real flux_v_cum(iip1,jjm,llm)
      real flux_w_cum(iip1,jjp1,llm)
      real Q_cum(iip1,jjp1,llm,nQ)
      real flux_uQ_cum(iip1,jjp1,llm,nQ)
      real flux_vQ_cum(iip1,jjm,llm,nQ)
      real flux_wQ_cum(iip1,jjp1,llm,nQ)
      real dQ(iip1,jjp1,llm,nQ)

      save ps_cum,masse_cum,flux_u_cum,flux_v_cum
      save Q_cum,flux_uQ_cum,flux_vQ_cum
      save flux_w_cum,flux_wQ_cum

c   champs de transport en moyenne zonale
      integer ntr,itr
      parameter (ntr=5)

      character*10,save :: znom(ntr,nQ)
      character*20,save :: znoml(ntr,nQ)
      character*10,save :: zunites(ntr,nQ)
      character*10,save :: znom2(ntr,nQ)
      character*20,save :: znom2l(ntr,nQ)
      character*10,save :: zunites2(ntr,nQ)
      character*10,save :: znom3(nQ)
      character*20,save :: znom3l(nQ)
      character*10,save :: zunites3(nQ)

      integer iave,itot,immc,itrs,istn
      data iave,itot,immc,itrs,istn/1,2,3,4,5/
      character*3 ctrs(ntr)
      data ctrs/'  ','TOT','MMC','TRS','STN'/

      real zvQ(jjm,llm,ntr,nQ),zvQtmp(jjm,llm)
      real zwQ(jjm,llm,ntr,nQ),zwQtmp(jjm,llm)
      real zavQ(jjm,ntr,nQ),psiQ(jjm,llm+1,nQ)
      real zawQ(jjm,llm,ntr,nQ)
      real zdQ(jjm,llm,nQ)
      real zmasse(jjm,llm),zavmasse(jjm),zawmasse(llm)
      real psbar(jjm)

      real zv(jjm,llm),zw(jjp1,llm),psi(jjm,llm+1)

c TENDANCES POUR MOMENT CINETIQUE
c -------------------------------

      integer ntdc,itdc
      parameter (ntdc=4)

      integer itdcdyn,itdcdis,itdcspg,itdcphy
      data    itdcdyn,itdcdis,itdcspg,itdcphy/1,2,3,4/

      character*6,save :: nomtdc(ntdc)

c   champ contenant les tendances du moment cinetique.
      real    tdc(iip1,jjp1,llm,ntdc)
      real    ztdc(jjm,llm,ntdc)   ! moyenne zonale
    
c   champs cumulés
      real tdc_cum(iip1,jjp1,llm,ntdc)
      save tdc_cum

c   integrations completes
      real mtot,mctot,dmctot(ntdc)

c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      character*10 infile

      integer fileid
      integer thoriid, zvertiid
      save fileid

      integer ndex3d(jjm*llm)
      real    ztmp3d(jjm,llm)

C   Variables locales
C
      integer tau0
      real zjulian
      integer zan, dayref
C
      real rlong(jjm),rlatg(jjm)



c=====================================================================
c   Initialisation
c=====================================================================

      ndex3d=0

      if (first) then

        if (planet_type.eq."venus") then 
            fact_geovenus = -1.
        else
            fact_geovenus = 1.
        endif

        icum=0
c       initialisation des fichiers
        first=.false.
c   ncum est la frequence de stokage en pas de temps
        ncum=dt_cum/dt_app
        if (abs(ncum*dt_app-dt_cum).gt.1.e-2*dt_app) then
         if (abs((ncum+1)*dt_app-dt_cum).lt.1.e-2*dt_app) then
           ncum = ncum+1
         elseif (abs((ncum-1)*dt_app-dt_cum).lt.1.e-2*dt_app) then
           ncum = ncum-1
         else
           WRITE(lunout,*)
     .            'Pb : le pas de cumule doit etre multiple du pas'
           WRITE(lunout,*)'dt_app=',dt_app
           WRITE(lunout,*)'dt_cum=',dt_cum
           WRITE(lunout,*)'ncum*dt_app=',ncum*dt_app
           WRITE(lunout,*)'ncum=',ncum
           stop
         endif
        endif

c VARIABLES ADVECTEES:

        nom(itemp)='temp'
        nom(igeop)='gz'
        nom(iecin)='ecin'
        nom(iang)='ang'
        nom(iu)='u'
        nom(iun)='un'

        unites(itemp)='K'
        unites(igeop)='m2/s2'
        unites(iecin)='m2/s2'
        unites(iang)='ang'
        unites(iu)='m/s'
        unites(iun)='un'

c avec tous les composes, ca fait trop.... Je les enleve
c       do num=1,ntrac
c        write(strd2,'(i2.2)') num
c        nom(6+num)='trac'//strd2
c        unites(6+num)='kg/kg'
c       enddo

c TENDANCES MOMENT CIN:
        
        nomtdc(itdcdyn) ='dmcdyn' 
        nomtdc(itdcdis) ='dmcdis' 
        nomtdc(itdcspg) ='dmcspg' 
        nomtdc(itdcphy) ='dmcphy' 

c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      infile='dynzon'

      zan = annee_ref
      dayref = day_ref
      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
c     tau0 = itau_dyn
      tau0 = 0
      itau = tau0
      
      rlong=0.
      rlatg=rlatv*180./pi*fact_geovenus
       
      call histbeg(infile, 1, rlong, jjm, rlatg,
     .             1, 1, 1, jjm,
     .             tau0, zjulian, dt_cum, thoriid, fileid)

C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'presnivs', 'Niveaux sigma','mb',
     .              llm, presnivs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder

      do iQ=1,nQ
         do itr=1,ntr
            if(itr.eq.1) then
               znom(itr,iQ)    =nom(iQ)
               znoml(itr,iQ)   =nom(iQ)
               zunites(itr,iQ) =unites(iQ)
            else
           znom(itr,iQ)    =ctrs(itr)//'v'//nom(iQ)
           znoml(itr,iQ)   ='transport : v * '//nom(iQ)//' '//ctrs(itr)
           zunites(itr,iQ) ='m/s * '//unites(iQ)
           znom2(itr,iQ)   =ctrs(itr)//'w'//nom(iQ)
           znom2l(itr,iQ)  ='transport: w * '//nom(iQ)//' '//ctrs(itr)
           zunites2(itr,iQ)='Pa/s * '//unites(iQ)
            endif
         enddo
               znom3(iQ)='d'//nom(iQ)
               znom3l(iQ)='convergence: '//nom(iQ)
               zunites3(iQ)=unites(iQ)//' /s'
c          print*,'DEBUG:',znom3(iQ),znom3l(iQ),zunites3(iQ)
      enddo

c   Declarations des champs avec dimension verticale

      if (1.eq.0) then  ! on les sort, ou pas...

c     do iQ=1,nQ
c !!!! JE NE SORS ICI QUE temp et ang POUR CAUSE DE PLACE !
      do iQ=1,4,3
         do itr=1,ntr
      IF (prt_level > 5)
     . WRITE(lunout,*)'var ',itr,iQ
     .      ,znom(itr,iQ),znoml(itr,iQ),zunites(itr,iQ)
            call histdef(fileid,znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
     .        32,'ave(X)',dt_cum,dt_cum)
         enddo
c transport vertical:
         do itr=2,ntr
      IF (prt_level > 5)
     . WRITE(lunout,*)'var ',itr,iQ
     .      ,znom2(itr,iQ),znom2l(itr,iQ),zunites2(itr,iQ)
            call histdef(fileid,znom2(itr,iQ),znom2l(itr,iQ),
     .        zunites2(itr,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
     .        32,'ins(X)',dt_cum,dt_cum)
         enddo

c Declarations pour convergences
      IF (prt_level > 5)
     . WRITE(lunout,*)'var ',iQ
     .      ,znom3(iQ),znom3l(iQ),zunites3(iQ)
            call histdef(fileid,znom3(iQ),znom3l(iQ),
     .        zunites3(iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
     .        32,'ins(X)',dt_cum,dt_cum)

c   Declarations pour les fonctions de courant
c   Non sorties ici...
c          call histdef(fileid,'psi'//nom(iQ)
c     .      ,'stream fn. '//znoml(itot,iQ),
c     .      zunites(itot,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
c     .      32,'ave(X)',dt_cum,dt_cum)

      enddo ! iQ

      endif ! 1=1 sortie ou non...

c   Declarations pour les champs de transport d'air
      call histdef(fileid, 'masse', 'masse',
     .             'kg', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', dt_cum, dt_cum)
      call histdef(fileid, 'v', 'v',
     .             'm/s', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', dt_cum, dt_cum)
      call histdef(fileid, 'w', 'w',
     .             'Pa/s', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ins(X)', dt_cum, dt_cum)

c   Declarations pour la fonction de courant
          call histdef(fileid,'psi','stream fn. MMC ','mega t/s',
     .      1,jjm,thoriid,llm,1,llm,zvertiid,
     .      32,'ave(X)',dt_cum,dt_cum)


c   Declarations pour les tendances de moment cinetique
      do itdc=1,ntdc
      call histdef(fileid, nomtdc(itdc), nomtdc(itdc),
     .             'ang/s', 1, jjm, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ins(X)', dt_cum, dt_cum)
      enddo

c   Declaration des champs 1D de transport en latitude
c     do iQ=1,nQ
c !!!! JE NE SORS ICI QUE temp et ang POUR CAUSE DE PLACE !
      do iQ=1,4,3
         do itr=2,ntr
            call histdef(fileid,'a'//znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjm,thoriid,1,1,1,-99,
     .        32,'ave(X)',dt_cum,dt_cum)
c JE VIRE LE VERTICAL POUR L'INSTANT
c           call histdef(fileid,'a'//znom2(itr,iQ),znom2l(itr,iQ),
c    .        zunites2(itr,iQ),1,jjm,thoriid,llm,1,llm,zvertiid,
c    .        32,'ins(X)',dt_cum,dt_cum)
         enddo
      enddo

               CALL histend(fileid)


      endif  ! first


c=====================================================================
c   Calcul des champs dynamiques
c   ----------------------------

c   énergie cinétique
      ucont(:,:,:)=0
      CALL covcont(llm,ucov,vcov,ucont,vcont)
      CALL enercin(vcov,ucov,vcont,ucont,ecin)

c   moment cinétique et tendances
      dudyn = 0.
      dudis = 0.
      duspg = 0.
      duphy = 0.
      do l=1,llm
         unat(:,:,l)=ucont(:,:,l)*cu(:,:)
         dudyn(:,2:jjm,l) = ducovdyn(:,2:jjm,l)/cu(:,2:jjm)
         dudis(:,2:jjm,l) = ducovdis(:,2:jjm,l)/cu(:,2:jjm)
         duspg(:,2:jjm,l) = ducovspg(:,2:jjm,l)/cu(:,2:jjm)
         duphy(:,2:jjm,l) = ducovphy(:,2:jjm,l)/cu(:,2:jjm)
         do j=1,jjp1
          ang(:,j,l)= rad*cos(rlatu(j))*
     .     ( unat(:,j,l) + rad*cos(rlatu(j))*omeg )
          tdc(:,j,l,1) = rad*cos(rlatu(j))*dudyn(:,j,l)
          tdc(:,j,l,2) = rad*cos(rlatu(j))*dudis(:,j,l)
          tdc(:,j,l,3) = rad*cos(rlatu(j))*duspg(:,j,l)
          tdc(:,j,l,4) = rad*cos(rlatu(j))*duphy(:,j,l)
         enddo
      enddo
c Normalisation: 
      ang = ang / (2./3. *rad*rad*omeg)
      do itdc=1,ntdc
        tdc(:,:,:,itdc)=tdc(:,:,:,itdc) / (2./3. *rad*rad*omeg)
      enddo

! ADAPTATION GCM POUR CP(T)
      call tpot2t(ip1jmp1*llm,teta,temp,pk)
      Q(:,:,:,itemp) = temp(:,:,:)
      Q(:,:,:,igeop) =phi(:,:,:)
      Q(:,:,:,iecin) =ecin(:,:,:)
      Q(:,:,:,iang)  =ang(:,:,:)
      Q(:,:,:,iu)    =unat(:,:,:)
      Q(:,:,:,iun)   =1.
c avec tous les composes, ca fait trop.... Je les enleve
c     do num=1,ntrac
c      Q(:,:,:,6+num)=trac(:,:,:,num)
c     enddo

c   calcul du flux de masse vertical (+ vers le bas)
      call convmas(flux_u,flux_v,convm)
      CALL vitvert(convm,w)

c=====================================================================
c   Cumul
c=====================================================================
c
      if(icum.EQ.0) then
         ps_cum      = 0.
         masse_cum   = 0.
         flux_u_cum  = 0.
         flux_v_cum  = 0.
         flux_w_cum  = 0.
         Q_cum       = 0.
         flux_vQ_cum = 0.
         flux_uQ_cum = 0.
         flux_wQ_cum = 0.
         tdc_cum     = 0.
      endif

      IF (prt_level > 5)
     . WRITE(lunout,*)'dans bilan_dyn ',icum,'->',icum+1
      icum=icum+1

c   accumulation des flux de masse horizontaux
      ps_cum          = ps_cum     + ps
      masse_cum       = masse_cum  + masse
      flux_u_cum      = flux_u_cum + flux_u
      flux_v_cum      = flux_v_cum + flux_v
      flux_w_cum      = flux_w_cum + w
      do iQ=1,nQ
      Q_cum(:,:,:,iQ) = Q_cum(:,:,:,iQ) + Q(:,:,:,iQ)*masse(:,:,:)
      enddo
      do itdc=1,ntdc
      tdc_cum(:,:,:,itdc) =
     .       tdc_cum(:,:,:,itdc) + tdc(:,:,:,itdc)*masse(:,:,:)
      enddo

c=====================================================================
c  FLUX ET TENDANCES
c=====================================================================

c   Flux longitudinal
c   -----------------
      do iQ=1,nQ
         do l=1,llm
            do j=1,jjp1
               do i=1,iim
                  flux_uQ_cum(i,j,l,iQ)=flux_uQ_cum(i,j,l,iQ)
     s            +flux_u(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i+1,j,l,iQ))
               enddo
               flux_uQ_cum(iip1,j,l,iQ)=flux_uQ_cum(1,j,l,iQ)
            enddo
         enddo
      enddo

c    flux méridien
c    -------------
      do iQ=1,nQ
         do l=1,llm
            do j=1,jjm
               do i=1,iip1
                  flux_vQ_cum(i,j,l,iQ)=flux_vQ_cum(i,j,l,iQ)
     s            +flux_v(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i,j+1,l,iQ))
               enddo
            enddo
         enddo
      enddo

c   Flux vertical
c   -------------
      do iQ=1,nQ
         do l=2,llm
            do j=1,jjp1
               do i=1,iip1
                  flux_wQ_cum(i,j,l,iQ)=flux_wQ_cum(i,j,l,iQ)
     s            +w(i,j,l)*0.5*(Q(i,j,l-1,iQ)+Q(i,j,l,iQ))
               enddo
            enddo
         enddo
         flux_wQ_cum(:,:,1,iQ)=0.0e0
      enddo

c    tendances
c    ---------

c   convergence horizontale
      call  convflu(flux_uQ_cum,flux_vQ_cum,llm*nQ,dQ)

c   calcul de la vitesse verticale
      call convmas(flux_u_cum,flux_v_cum,convm)
      CALL vitvert(convm,w)

c  ajustement tendances (vitesse verticale)
      do iQ=1,nQ
         do l=1,llm-1
            do j=1,jjp1
               do i=1,iip1
                  ww=-0.5*w(i,j,l+1)*(Q(i,j,l,iQ)+Q(i,j,l+1,iQ))
                  dQ(i,j,l  ,iQ)=dQ(i,j,l  ,iQ)-ww
                  dQ(i,j,l+1,iQ)=dQ(i,j,l+1,iQ)+ww
               enddo
            enddo
         enddo
      enddo
      IF (prt_level > 5)
     . WRITE(lunout,*)'Apres les calculs fait a chaque pas'

c=====================================================================
c   PAS DE TEMPS D'ECRITURE
c=====================================================================
      if (icum.eq.ncum) then
c=====================================================================

      time=time+dt_cum
      itau=itau+1

      IF (prt_level > 5)
     . WRITE(lunout,*)'Pas d ecriture'

c   Normalisation
      do iQ=1,nQ
         Q_cum(:,:,:,iQ) = Q_cum(:,:,:,iQ)/masse_cum(:,:,:)
         dQ(:,:,:,iQ)    = dQ(:,:,:,iQ)   /masse_cum(:,:,:)
      enddo
      do itdc=1,ntdc
         tdc_cum(:,:,:,itdc) = tdc_cum(:,:,:,itdc)/masse_cum(:,:,:)
      enddo

      zz=1./REAL(ncum)
      ps_cum      = ps_cum      *zz
      masse_cum   = masse_cum   *zz
      flux_u_cum  = flux_u_cum  *zz
      flux_v_cum  = flux_v_cum  *zz
      flux_w_cum  = flux_w_cum  *zz
      flux_uQ_cum = flux_uQ_cum *zz
      flux_vQ_cum = flux_vQ_cum *zz
      flux_wQ_cum = flux_wQ_cum *zz

c Integration complete
      mtot  = 0.
      mctot  = 0.
      dmctot = 0.
      do l=1,llm
       do j=2,jjm
        do i=1,iim
          mtot  = mtot  + masse_cum(i,j,l)
          mctot = mctot + Q_cum(i,j,l,iang)*masse_cum(i,j,l)
        enddo
       enddo
      enddo
      mctot = mctot/mtot
      do itdc=1,ntdc
      do l=1,llm
       do j=2,jjm
        do i=1,iim
          dmctot(itdc) = dmctot(itdc) 
     .               + tdc_cum(i,j,l,itdc)*masse_cum(i,j,l)/mtot
        enddo
       enddo
      enddo
      enddo

c=====================================================================
c   Transport méridien
c=====================================================================

c   cumul zonal des masses des mailles
c   ----------------------------------
      zv=0.
      zw=0.
      zmasse=0.
      call massbar(masse_cum,massebx,masseby)

c moy zonale de la ps cumulee
         do j=1,jjm
            psbar(j)=0.
            do i=1,iim
               psbar(j)=psbar(j)+ps_cum(i,j)/iim
            enddo
         enddo

      do l=1,llm
         do j=1,jjm
            do i=1,iim
               zmasse(j,l)=zmasse(j,l)+masseby(i,j,l)
               zv(j,l)=zv(j,l)+flux_v_cum(i,j,l)
               zw(j,l)=zw(j,l)+flux_w_cum(i,j,l)
            enddo
            zfactv(j,l)=cv(1,j)/zmasse(j,l)
            zfactw(j,l)=((ap(l)-ap(l+1))+psbar(j)*(bp(l)-bp(l+1)))
     .                    /zmasse(j,l) 
         enddo
            do i=1,iim
               zw(jjp1,l)=zw(jjp1,l)+flux_w_cum(i,jjp1,l)
            enddo
      enddo

c     print*,'3OK'
c   --------------------------------------------------------------
c   calcul de la moyenne zonale du transport :
c   ------------------------------------------
c
c                                     --
c TOT : la circulation totale       [ vq ]
c
c                                      -     -
c MMC : mean meridional circulation [ v ] [ q ]
c
c                                     ----      --       - -
c TRS : transitoires                [ v'q'] = [ vq ] - [ v q ]
c
c                                     - * - *       - -       -     -
c STT : stationaires                [ v   q   ] = [ v q ] - [ v ] [ q ]
c
c                                              - -
c    on utilise aussi l'intermediaire TMP :  [ v q ]
c
c    la variable zfactv transforme un transport meridien cumule
c    en kg/s * unte-du-champ-transporte en m/s * unite-du-champ-transporte
c    la variable zfactw transforme un transport vertical cumule
c    en kg/s * unte-du-champ-transporte en Pa/s * unite-du-champ-transporte
c
c   --------------------------------------------------------------


c   ----------------------------------------
c   Transport dans le plan latitude-altitude
c   ----------------------------------------

      zvQ=0.
      zwQ=0.
      zdQ=0.
      psiQ=0.

      do iQ=1,nQ

c   transport meridien
         zvQtmp=0.
         do l=1,llm
            do j=1,jjm
c              print*,'j,l,iQ=',j,l,iQ
c   Calcul des moyennes zonales du transport total et de zvQtmp
               do i=1,iim
                  zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)
     s                            +flux_vQ_cum(i,j,l,iQ)
                  zqy=      0.5*(Q_cum(i,j,l,iQ)*masse_cum(i,j,l)+
     s                           Q_cum(i,j+1,l,iQ)*masse_cum(i,j+1,l))
                  zvQtmp(j,l)=zvQtmp(j,l)+flux_v_cum(i,j,l)*zqy
     s             /(0.5*(masse_cum(i,j,l)+masse_cum(i,j+1,l)))
                  zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)+zqy
               enddo
c              print*,'aOK'
c   Decomposition
               zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)/zmasse(j,l)
               zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)*zfactv(j,l)
               zvQtmp(j,l)=zvQtmp(j,l)*zfactv(j,l)
               zvQ(j,l,immc,iQ)=zv(j,l)*zvQ(j,l,iave,iQ)*zfactv(j,l)
               zvQ(j,l,itrs,iQ)=zvQ(j,l,itot,iQ)-zvQtmp(j,l)
               zvQ(j,l,istn,iQ)=zvQtmp(j,l)-zvQ(j,l,immc,iQ)
            enddo
         enddo
c   fonction de courant meridienne pour la quantite Q
         do l=llm,1,-1
            do j=1,jjm
             psiQ(j,l,iQ)=psiQ(j,l+1,iQ)+zvQ(j,l,itot,iQ)/zfactv(j,l)
            enddo
         enddo
!!      enddo

c   transport vertical
         zwQtmp=0.
         do l=1,llm
            do j=1,jjm
c              print*,'j,l,iQ=',j,l,iQ
c   Calcul des moyennes zonales du transport vertical total et de zwQtmp
               do i=1,iim
                  zwQ(j,l,itot,iQ)=zwQ(j,l,itot,iQ)
     s                            +flux_wQ_cum(i,j,l,iQ)
                  zqy=      0.5*(Q_cum(i,j,l,iQ)*masse_cum(i,j,l)+
     s                           Q_cum(i,j+1,l,iQ)*masse_cum(i,j+1,l))
                  zwQtmp(j,l)=zwQtmp(j,l)+flux_w_cum(i,j,l)*zqy
     s             /(0.5*(masse_cum(i,j,l)+masse_cum(i,j+1,l)))
                  zwQ(j,l,iave,iQ)=zwQ(j,l,iave,iQ)+zqy
               enddo
c   Decomposition
               zwQ(j,l,iave,iQ)=zwQ(j,l,iave,iQ)/zmasse(j,l)
               zwQ(j,l,itot,iQ)=zwQ(j,l,itot,iQ)*zfactw(j,l)
               zwQtmp(j,l)=zwQtmp(j,l)*zfactw(j,l)
               zwQ(j,l,immc,iQ)=zw(j,l)*zwQ(j,l,iave,iQ)*zfactw(j,l)
               zwQ(j,l,itrs,iQ)=zwQ(j,l,itot,iQ)-zwQtmp(j,l)
               zwQ(j,l,istn,iQ)=zwQtmp(j,l)-zwQ(j,l,immc,iQ)
            enddo
         enddo

c   convergence
c   Calcul moyenne zonale de la convergence totale
         do l=1,llm
            do j=1,jjm
c              print*,'j,l,iQ=',j,l,iQ
               do i=1,iim
                  zdQ(j,l,iQ)=zdQ(j,l,iQ) +
     .                   ( dQ(i,j,l,iQ)   * masse_cum(i,j,l)
     .                   + dQ(i,j+1,l,iQ) * masse_cum(i,j+1,l))
     .                 / ( masse_cum(i,j,l)+masse_cum(i,j+1,l))
               enddo
            enddo
         enddo
      enddo ! of do iQ=1,nQ

c   fonction de courant pour la circulation meridienne moyenne
      psi=0.
      do l=llm,1,-1
         do j=1,jjm
            psi(j,l)= psi(j,l+1)+zv(j,l)
            zv(j,l) = zv(j,l)*zfactv(j,l)
            zw(j,l) = 0.5*(zw(j,l)+zw(j+1,l))*zfactw(j,l)
         enddo
      enddo

c   Calcul moyenne zonale des tendances moment cin.
      ztdc=0.
      do itdc=1,ntdc
         do l=1,llm
            do j=1,jjm
               do i=1,iim
                  ztdc(j,l,itdc)=ztdc(j,l,itdc) +
     .            ( tdc_cum(i,j,l,itdc)   * masse_cum(i,j,l)
     .            + tdc_cum(i,j+1,l,itdc) * masse_cum(i,j+1,l))
     .          / ( masse_cum(i,j,l)+masse_cum(i,j+1,l))
               enddo
            enddo
         enddo
      enddo

c     print*,'4OK'

c--------------------------------------
c--------------------------------------
c   sorties proprement dites
c--------------------------------------
c--------------------------------------

      if (i_sortie.eq.1) then

c sortie des integrations completes dans le listing
      write(*,'(A12,5(1PE11.4,X))') "BILANMCDYN  ",mctot,dmctot

c sorties dans fichier dynzon

      if (1.eq.0) then  ! on les sort, ou pas...

c avec tous les composes, ca fait trop.... Je les enleve
c      do iQ=1,nQ
c      do iQ=1,6
c !!!! JE NE SORS ICI QUE temp et ang POUR CAUSE DE PLACE !
      do iQ=1,4,3

         ztmp3d(:,:)= zvQ(:,:,1,iQ) ! valeur moyenne
            call histwrite(fileid,znom(1,iQ),itau,ztmp3d
     s      ,jjm*llm,ndex3d)
       do itr=2,ntr
         ztmp3d(:,:)= zvQ(:,:,itr,iQ)*fact_geovenus ! transport horizontal
            call histwrite(fileid,znom(itr,iQ),itau,ztmp3d
     s      ,jjm*llm,ndex3d)
       enddo

       do itr=2,ntr
         ztmp3d(:,:)=zwQ(:,:,itr,iQ)
            call histwrite(fileid,znom2(itr,iQ),itau,ztmp3d
     s      ,jjm*llm,ndex3d)
       enddo
       
         ztmp3d(:,:)= zdQ(:,:,iQ)
            call histwrite(fileid,znom3(iQ),itau,ztmp3d
     s      ,jjm*llm,ndex3d)

c        ztmp3d(:,:)= psiQ(:,1:llm,iQ)*fact_geovenus
c        call histwrite(fileid,'psi'//nom(iQ),itau,ztmp3d
c    s      ,jjm*llm,ndex3d)
      enddo

      endif ! 1=1 sortie ou non...

      ztmp3d=zmasse
      call histwrite(fileid,'masse',itau,ztmp3d
     s   ,jjm*llm,ndex3d)
      
      ztmp3d= zv*fact_geovenus
      call histwrite(fileid,'v',itau,ztmp3d
     s   ,jjm*llm,ndex3d)
      ztmp3d(:,:)=zw(1:jjm,:)
      call histwrite(fileid,'w',itau,ztmp3d
     s   ,jjm*llm,ndex3d)
      ztmp3d= psi(:,1:llm)*1.e-9*fact_geovenus
      call histwrite(fileid,'psi',itau,ztmp3d,jjm*llm,ndex3d)

      do itdc=1,ntdc
         ztmp3d(:,:)= ztdc(:,:,itdc)
         call histwrite(fileid,nomtdc(itdc),itau,ztmp3d
     s    ,jjm*llm,ndex3d)
      enddo

      endif ! i_sortie


c   -----------------
c   Moyenne verticale
c   -----------------

      zavmasse=0.
      do l=1,llm
         zavmasse(:)=zavmasse(:)+zmasse(:,l)
      enddo
      zavQ=0.

c avec tous les composes, ca fait trop.... Je les enleve
c      do iQ=1,nQ
c      do iQ=1,6
c !!!! JE NE SORS ICI QUE temp et ang POUR CAUSE DE PLACE !
      do iQ=1,4,3
         do itr=2,ntr
            do l=1,llm
               zavQ(:,itr,iQ)=zavQ(:,itr,iQ)+zvQ(:,l,itr,iQ)*zmasse(:,l)
            enddo
            zavQ(:,itr,iQ)=zavQ(:,itr,iQ)/zavmasse(:)
      if (i_sortie.eq.1) then
         ztmp3d=0.0
         ztmp3d(:,1)= zavQ(:,itr,iQ)*fact_geovenus
         call histwrite(fileid,'a'//znom(itr,iQ),itau,ztmp3d
     .      ,jjm*llm,ndex3d)     
      endif
         enddo
      enddo

c   ------------------
c   Moyenne meridienne
c   ------------------

      zawmasse=0.
      do j=1,jjm
           do l=1,llm
         zawmasse(l)=zawmasse(l)+zmasse(j,l)
           enddo
      enddo
      zawQ=0.

c avec tous les composes, ca fait trop.... Je les enleve
c      do iQ=1,nQ
c      do iQ=1,6
c !!!! JE NE SORS ICI QUE temp et ang POUR CAUSE DE PLACE !
      do iQ=1,4,3
         do itr=2,ntr
           do l=1,llm
            do j=1,jjm
          zawQ(1,l,itr,iQ)=zawQ(1,l,itr,iQ)+zwQ(j,l,itr,iQ)*zmasse(j,l)
            enddo
            zawQ(1,l,itr,iQ)=zawQ(1,l,itr,iQ)/zawmasse(l)
           enddo
      if (i_sortie.eq.1) then
         ztmp3d=0.0
           do l=1,llm
         ztmp3d(1,l)=zawQ(1,l,itr,iQ)
           enddo
c JE VIRE LE VERTICAL POUR L'INSTANT
c        call histwrite(fileid,'a'//znom2(itr,iQ),itau,ztmp3d
c    .      ,jjm*llm,ndex3d)     
      endif
         enddo
      enddo

      call histsync(fileid)

c=====================================================================
c/////////////////////////////////////////////////////////////////////
      icum=0                  !///////////////////////////////////////
      endif ! icum.eq.ncum    !///////////////////////////////////////
c/////////////////////////////////////////////////////////////////////
c=====================================================================

      return
      end
