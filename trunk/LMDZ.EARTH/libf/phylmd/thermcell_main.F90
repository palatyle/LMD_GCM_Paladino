!
! $Id: thermcell_main.F90 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE thermcell_main(itap,ngrid,nlay,ptimestep  &
     &                  ,pplay,pplev,pphi,debut  &
     &                  ,pu,pv,pt,po  &
     &                  ,pduadj,pdvadj,pdtadj,pdoadj  &
     &                  ,fm0,entr0,detr0,zqta,zqla,lmax  &
     &                  ,ratqscth,ratqsdiff,zqsatth  &
     &                  ,r_aspect,l_mix,tau_thermals,iflag_thermals_ed &
     &                  ,Ale_bl,Alp_bl,lalim_conv,wght_th &
     &                  ,zmax0, f0,zw2,fraca,ztv &
     &                  ,zpspsk,ztla,zthl)

      USE dimphy
      USE comgeomphy , ONLY:rlond,rlatd
      IMPLICIT NONE

!=======================================================================
!   Auteurs: Frederic Hourdin, Catherine Rio, Anne Mathieu
!   Version du 09.02.07
!   Calcul du transport vertical dans la couche limite en presence
!   de "thermiques" explicitement representes avec processus nuageux
!
!   Reecriture a partir d'un listing papier a Habas, le 14/02/00
!
!   le thermique est suppose homogene et dissipe par melange avec
!   son environnement. la longueur l_mix controle l'efficacite du
!   melange
!
!   Le calcul du transport des differentes especes se fait en prenant
!   en compte:
!     1. un flux de masse montant
!     2. un flux de masse descendant
!     3. un entrainement
!     4. un detrainement
!
!=======================================================================

!-----------------------------------------------------------------------
!   declarations:
!   -------------

#include "dimensions.h"
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"

!   arguments:
!   ----------

!IM 140508
      INTEGER itap

      INTEGER ngrid,nlay,w2di
      real tau_thermals
      integer iflag_thermals_ed
      real ptimestep,l_mix,r_aspect
      REAL pt(ngrid,nlay),pdtadj(ngrid,nlay)
      REAL pu(ngrid,nlay),pduadj(ngrid,nlay)
      REAL pv(ngrid,nlay),pdvadj(ngrid,nlay)
      REAL po(ngrid,nlay),pdoadj(ngrid,nlay)
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      real pphi(ngrid,nlay)

!   local:
!   ------

      integer icount
      data icount/0/
      save icount
!$OMP THREADPRIVATE(icount)

      integer,save :: igout=1
!$OMP THREADPRIVATE(igout)
      integer,save :: lunout1=6
!$OMP THREADPRIVATE(lunout1)
      integer,save :: lev_out=10
!$OMP THREADPRIVATE(lev_out)

      INTEGER ig,k,l,ll
      real zsortie1d(klon)
      INTEGER lmax(klon),lmin(klon),lalim(klon)
      INTEGER lmix(klon)
      INTEGER lmix_bis(klon)
      real linter(klon)
      real zmix(klon)
      real zmax(klon),zw2(klon,klev+1),ztva(klon,klev),zw_est(klon,klev+1),ztva_est(klon,klev)
!      real fraca(klon,klev)

      real zmax_sec(klon)
!on garde le zmax du pas de temps precedent
      real zmax0(klon)
!FH/IM     save zmax0

      real lambda

      real zlev(klon,klev+1),zlay(klon,klev)
      real deltaz(klon,klev)
      REAL zh(klon,klev)
      real zthl(klon,klev),zdthladj(klon,klev)
      REAL ztv(klon,klev)
      real zu(klon,klev),zv(klon,klev),zo(klon,klev)
      real zl(klon,klev)
      real zsortie(klon,klev)
      real zva(klon,klev)
      real zua(klon,klev)
      real zoa(klon,klev)

      real zta(klon,klev)
      real zha(klon,klev)
      real fraca(klon,klev+1)
      real zf,zf2
      real thetath2(klon,klev),wth2(klon,klev),wth3(klon,klev)
      real q2(klon,klev)
! FH probleme de dimensionnement avec l'allocation dynamique
!     common/comtherm/thetath2,wth2
      real wq(klon,klev)
      real wthl(klon,klev)
      real wthv(klon,klev)
    
      real ratqscth(klon,klev)
      real var
      real vardiff
      real ratqsdiff(klon,klev)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev),masse(klon,klev)
      real zpspsk(klon,klev)

      real wmax(klon)
      real wmax_tmp(klon)
      real wmax_sec(klon)
      real fm0(klon,klev+1),entr0(klon,klev),detr0(klon,klev)
      real fm(klon,klev+1),entr(klon,klev),detr(klon,klev)

      real ztla(klon,klev),zqla(klon,klev),zqta(klon,klev)
!niveau de condensation
      integer nivcon(klon)
      real zcon(klon)
      REAL CHI
      real zcon2(klon)
      real pcon(klon)
      real zqsat(klon,klev)
      real zqsatth(klon,klev) 

      real f_star(klon,klev+1),entr_star(klon,klev)
      real detr_star(klon,klev)
      real alim_star_tot(klon)
      real alim_star(klon,klev)
      real alim_star_clos(klon,klev)
      real f(klon), f0(klon)
!FH/IM     save f0
      real zlevinter(klon)
      logical debut
       real seuil
      real csc(klon,klev)

!
      !nouvelles variables pour la convection
      real Ale_bl(klon)
      real Alp_bl(klon)
      real alp_int(klon)
      real ale_int(klon)
      integer n_int(klon)
      real fm_tot(klon)
      real wght_th(klon,klev)
      integer lalim_conv(klon)
!v1d     logical therm
!v1d     save therm

      character*2 str2
      character*10 str10

      character (len=20) :: modname='thermcell_main'
      character (len=80) :: abort_message

      EXTERNAL SCOPY
!

!-----------------------------------------------------------------------
!   initialisation:
!   ---------------
!

      seuil=0.25

      if (debut)  then
         fm0=0.
         entr0=0.
         detr0=0.


#undef wrgrads_thermcell
#ifdef wrgrads_thermcell
! Initialisation des sorties grads pour les thermiques.
! Pour l'instant en 1D sur le point igout.
! Utilise par thermcell_out3d.h
         str10='therm'
         call inigrads(1,1,rlond(igout),1.,-180.,180.,jjm, &
     &   rlatd(igout),-90.,90.,1.,llm,pplay(igout,:),1.,   &
     &   ptimestep,str10,'therm ')
#endif



      endif

      fm=0. ; entr=0. ; detr=0.


      icount=icount+1

!IM 090508 beg
!print*,'====================================================================='
!print*,'====================================================================='
!print*,' PAS ',icount,' PAS ',icount,' PAS ',icount,' PAS ',icount
!print*,'====================================================================='
!print*,'====================================================================='
!IM 090508 end

      if (prt_level.ge.1) print*,'thermcell_main V4'

       sorties=.true.
      IF(ngrid.NE.klon) THEN
         PRINT*
         PRINT*,'STOP dans convadj'
         PRINT*,'ngrid    =',ngrid
         PRINT*,'klon  =',klon
      ENDIF
!
!     write(lunout,*)'WARNING thermcell_main f0=max(f0,1.e-2)'
     do ig=1,klon
      if (prt_level.ge.20) then
       print*,'th_main ig f0',ig,f0(ig)
      endif
         f0(ig)=max(f0(ig),1.e-2)
         zmax0(ig)=max(zmax0(ig),40.)
!IMmarche pas ?!       if (f0(ig)<1.e-2) f0(ig)=1.e-2
     enddo

!-----------------------------------------------------------------------
! Calcul de T,q,ql a partir de Tl et qT dans l environnement
!   --------------------------------------------------------------------
!
      CALL thermcell_env(ngrid,nlay,po,pt,pu,pv,pplay,  &
     &           pplev,zo,zh,zl,ztv,zthl,zu,zv,zpspsk,zqsat,lev_out)
       
      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_env'

!------------------------------------------------------------------------
!                       --------------------
!
!
!                       + + + + + + + + + + +
!
!
!  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
!  wh,wt,wo ...
!
!                       + + + + + + + + + + +  zh,zu,zv,zo,rho
!
!
!                       --------------------   zlev(1)
!                       \\\\\\\\\\\\\\\\\\\\
!
!

!-----------------------------------------------------------------------
!   Calcul des altitudes des couches
!-----------------------------------------------------------------------

      do l=2,nlay
         zlev(:,l)=0.5*(pphi(:,l)+pphi(:,l-1))/RG
      enddo
         zlev(:,1)=0.
         zlev(:,nlay+1)=(2.*pphi(:,klev)-pphi(:,klev-1))/RG
      do l=1,nlay
         zlay(:,l)=pphi(:,l)/RG
      enddo
!calcul de l epaisseur des couches
      do l=1,nlay
         deltaz(:,l)=zlev(:,l+1)-zlev(:,l)
      enddo

!     print*,'2 OK convect8'
!-----------------------------------------------------------------------
!   Calcul des densites
!-----------------------------------------------------------------------

      do l=1,nlay
         rho(:,l)=pplay(:,l)/(zpspsk(:,l)*RD*ztv(:,l))
      enddo

!IM
     if (prt_level.ge.10)write(lunout,*)                                &
    &    'WARNING thermcell_main rhobarz(:,1)=rho(:,1)'
      rhobarz(:,1)=rho(:,1)

      do l=2,nlay
         rhobarz(:,l)=0.5*(rho(:,l)+rho(:,l-1))
      enddo

!calcul de la masse
      do l=1,nlay
         masse(:,l)=(pplev(:,l)-pplev(:,l+1))/RG
      enddo

      if (prt_level.ge.1) print*,'thermcell_main apres initialisation'

!------------------------------------------------------------------
!
!             /|\
!    --------  |  F_k+1 -------   
!                              ----> D_k
!             /|\              <---- E_k , A_k
!    --------  |  F_k --------- 
!                              ----> D_k-1
!                              <---- E_k-1 , A_k-1
!
!
!
!
!
!    ---------------------------
!
!    ----- F_lmax+1=0 ----------         \
!            lmax     (zmax)              |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |  E
!    ---------------------------          |  D
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------  \       |
!            lalim                 |      |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!                                  | A    |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!    lmin  (=1 pour le moment)     |      |
!    ----- F_lmin=0 ------------  /      /
!
!    ---------------------------
!    //////////////////////////
!
!
!=============================================================================
!  Calculs initiaux ne faisant pas intervenir les changements de phase
!=============================================================================

!------------------------------------------------------------------
!  1. alim_star est le profil vertical de l'alimentation a la base du
!     panache thermique, calcule a partir de la flotabilite de l'air sec
!  2. lmin et lalim sont les indices inferieurs et superieurs de alim_star
!------------------------------------------------------------------
!
      entr_star=0. ; detr_star=0. ; alim_star=0. ; alim_star_tot=0.
      lmin=1

!-----------------------------------------------------------------------------
!  3. wmax_sec et zmax_sec sont les vitesses et altitudes maximum d'un
!     panache sec conservatif (e=d=0) alimente selon alim_star 
!     Il s'agit d'un calcul de type CAPE
!     zmax_sec est utilise pour determiner la geometrie du thermique.
!------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!calcul du melange et des variables dans le thermique
!--------------------------------------------------------------------------------
!
      if (prt_level.ge.1) print*,'avant thermcell_plume ',lev_out
!IM 140508   CALL thermcell_plume(ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,  &

! Gestion temporaire de plusieurs appels à thermcell_plume au travers
! de la variable iflag_thermals

!      print*,'THERM thermcell_main iflag_thermals_ed=',iflag_thermals_ed
      if (iflag_thermals_ed<=9) then
!         print*,'THERM NOUVELLE/NOUVELLE Arnaud'
         CALL thermcell_plume(itap,ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,&
     &    zlev,pplev,pphi,zpspsk,alim_star,alim_star_tot,  &
     &    lalim,f0,detr_star,entr_star,f_star,csc,ztva,  &
     &    ztla,zqla,zqta,zha,zw2,zw_est,ztva_est,zqsatth,lmix,lmix_bis,linter &
     &    ,lev_out,lunout1,igout)

      elseif (iflag_thermals_ed>9) then
!        print*,'THERM RIO et al 2010, version d Arnaud'
         CALL thermcellV1_plume(itap,ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,&
     &    zlev,pplev,pphi,zpspsk,alim_star,alim_star_tot,  &
     &    lalim,f0,detr_star,entr_star,f_star,csc,ztva,  &
     &    ztla,zqla,zqta,zha,zw2,zw_est,ztva_est,zqsatth,lmix,lmix_bis,linter &
     &    ,lev_out,lunout1,igout)

      endif

      if (prt_level.ge.1) print*,'apres thermcell_plume ',lev_out

      call test_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_plum lalim ')
      call test_ltherm(ngrid,nlay,pplev,pplay,lmix ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_plum lmix  ')

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_plume'
      if (prt_level.ge.10) then
         write(lunout1,*) 'Dans thermcell_main 2'
         write(lunout1,*) 'lmin ',lmin(igout)
         write(lunout1,*) 'lalim ',lalim(igout)
         write(lunout1,*) ' ig l alim_star entr_star detr_star f_star '
         write(lunout1,'(i6,i4,4e15.5)') (igout,l,alim_star(igout,l),entr_star(igout,l),detr_star(igout,l) &
     &    ,f_star(igout,l+1),l=1,nint(linter(igout))+5)
      endif

!-------------------------------------------------------------------------------
! Calcul des caracteristiques du thermique:zmax,zmix,wmax
!-------------------------------------------------------------------------------
!
      CALL thermcell_height(ngrid,nlay,lalim,lmin,linter,lmix,zw2,  &
     &           zlev,lmax,zmax,zmax0,zmix,wmax,lev_out)
! Attention, w2 est transforme en sa racine carree dans cette routine
! Le probleme vient du fait que linter et lmix sont souvent égaux à 1.
      wmax_tmp=0.
      do  l=1,nlay
         wmax_tmp(:)=max(wmax_tmp(:),zw2(:,l))
      enddo
!     print*,"ZMAX ",lalim,lmin,linter,lmix,lmax,zmax,zmax0,zmix,wmax



      call test_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lalim ')
      call test_ltherm(ngrid,nlay,pplev,pplay,lmin ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmin  ')
      call test_ltherm(ngrid,nlay,pplev,pplay,lmix ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmix  ')
      call test_ltherm(ngrid,nlay,pplev,pplay,lmax ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmax  ')

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_height'

!-------------------------------------------------------------------------------
! Fermeture,determination de f
!-------------------------------------------------------------------------------
!
!
!!      write(lunout,*)'THERM NOUVEAU XXXXX'
      CALL thermcell_dry(ngrid,nlay,zlev,pphi,ztv,alim_star,  &
    &                      lalim,lmin,zmax_sec,wmax_sec,lev_out)

call test_ltherm(ngrid,nlay,pplev,pplay,lmin,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_dry  lmin  ')
call test_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_dry  lalim ')

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_dry'
      if (prt_level.ge.10) then
         write(lunout1,*) 'Dans thermcell_main 1b'
         write(lunout1,*) 'lmin ',lmin(igout)
         write(lunout1,*) 'lalim ',lalim(igout)
         write(lunout1,*) ' ig l alim_star entr_star detr_star f_star '
         write(lunout1,'(i6,i4,e15.5)') (igout,l,alim_star(igout,l) &
     &    ,l=1,lalim(igout)+4)
      endif




! Choix de la fonction d'alimentation utilisee pour la fermeture.
! Apparemment sans importance
      alim_star_clos(:,:)=alim_star(:,:)
      alim_star_clos(:,:)=entr_star(:,:)+alim_star(:,:)

! Appel avec la version seche
      CALL thermcell_closure(ngrid,nlay,r_aspect,ptimestep,rho,  &
     &   zlev,lalim,alim_star_clos,f_star,zmax_sec,wmax_sec,f,lev_out)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel avec les zmax et wmax tenant compte de la condensation
! Semble moins bien marcher
!     CALL thermcell_closure(ngrid,nlay,r_aspect,ptimestep,rho,  &
!    &   zlev,lalim,alim_star,f_star,zmax,wmax,f,lev_out)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(prt_level.ge.1)print*,'thermcell_closure apres thermcell_closure'

      if (tau_thermals>1.) then
         lambda=exp(-ptimestep/tau_thermals)
         f0=(1.-lambda)*f+lambda*f0
      else
         f0=f
      endif

! Test valable seulement en 1D mais pas genant
      if (.not. (f0(1).ge.0.) ) then
              abort_message = '.not. (f0(1).ge.0.)'
              CALL abort_gcm (modname,abort_message,1)
      endif

!-------------------------------------------------------------------------------
!deduction des flux
!-------------------------------------------------------------------------------

      CALL thermcell_flux2(ngrid,nlay,ptimestep,masse, &
     &       lalim,lmax,alim_star,  &
     &       entr_star,detr_star,f,rhobarz,zlev,zw2,fm,entr,  &
     &       detr,zqla,lev_out,lunout1,igout)
!IM 060508    &       detr,zqla,zmax,lev_out,lunout,igout)

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_flux'
      call test_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_flux lalim ')
      call test_ltherm(ngrid,nlay,pplev,pplay,lmax ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_flux lmax  ')

!------------------------------------------------------------------
!   On ne prend pas directement les profils issus des calculs precedents
!   mais on s'autorise genereusement une relaxation vers ceci avec
!   une constante de temps tau_thermals (typiquement 1800s).
!------------------------------------------------------------------

      if (tau_thermals>1.) then
         lambda=exp(-ptimestep/tau_thermals)
         fm0=(1.-lambda)*fm+lambda*fm0
         entr0=(1.-lambda)*entr+lambda*entr0
         detr0=(1.-lambda)*detr+lambda*detr0
      else
         fm0=fm
         entr0=entr
         detr0=detr
      endif

!c------------------------------------------------------------------
!   calcul du transport vertical
!------------------------------------------------------------------

      call thermcell_dq(ngrid,nlay,ptimestep,fm0,entr0,masse,  &
     &                    zthl,zdthladj,zta,lev_out)
      call thermcell_dq(ngrid,nlay,ptimestep,fm0,entr0,masse,  &
     &                   po,pdoadj,zoa,lev_out)

!------------------------------------------------------------------
! Calcul de la fraction de l'ascendance
!------------------------------------------------------------------
      do ig=1,klon
         fraca(ig,1)=0.
         fraca(ig,nlay+1)=0.
      enddo
      do l=2,nlay
         do ig=1,klon
            if (zw2(ig,l).gt.1.e-10) then
            fraca(ig,l)=fm(ig,l)/(rhobarz(ig,l)*zw2(ig,l))
            else
            fraca(ig,l)=0.
            endif
         enddo
      enddo
     
!------------------------------------------------------------------
!  calcul du transport vertical du moment horizontal
!------------------------------------------------------------------

!IM 090508  
      if (1.eq.1) then
!IM 070508 vers. _dq       
!     if (1.eq.0) then


! Calcul du transport de V tenant compte d'echange par gradient
! de pression horizontal avec l'environnement

         call thermcell_dv2(ngrid,nlay,ptimestep,fm0,entr0,masse  &
     &    ,fraca,zmax  &
     &    ,zu,zv,pduadj,pdvadj,zua,zva,lev_out)

      else

! calcul purement conservatif pour le transport de V
         call thermcell_dq(ngrid,nlay,ptimestep,fm0,entr0,masse  &
     &    ,zu,pduadj,zua,lev_out)
         call thermcell_dq(ngrid,nlay,ptimestep,fm0,entr0,masse  &
     &    ,zv,pdvadj,zva,lev_out)
      endif

!     print*,'13 OK convect8'
      do l=1,nlay
         do ig=1,ngrid
           pdtadj(ig,l)=zdthladj(ig,l)*zpspsk(ig,l)  
         enddo
      enddo

      if (prt_level.ge.1) print*,'14 OK convect8'
!------------------------------------------------------------------
!   Calculs de diagnostiques pour les sorties
!------------------------------------------------------------------
!calcul de fraca pour les sorties
      
      if (sorties) then
      if (prt_level.ge.1) print*,'14a OK convect8'
! calcul du niveau de condensation
! initialisation
      do ig=1,ngrid
         nivcon(ig)=0
         zcon(ig)=0.
      enddo 
!nouveau calcul
      do ig=1,ngrid
      CHI=zh(ig,1)/(1669.0-122.0*zo(ig,1)/zqsat(ig,1)-zh(ig,1))
      pcon(ig)=pplay(ig,1)*(zo(ig,1)/zqsat(ig,1))**CHI
      enddo
!IM   do k=1,nlay
      do k=1,nlay-1
         do ig=1,ngrid
         if ((pcon(ig).le.pplay(ig,k))  &
     &      .and.(pcon(ig).gt.pplay(ig,k+1))) then
            zcon2(ig)=zlay(ig,k)-(pcon(ig)-pplay(ig,k))/(RG*rho(ig,k))/100.
         endif
         enddo
      enddo
!IM
      do ig=1,ngrid
        if (pcon(ig).le.pplay(ig,nlay)) then 
           zcon2(ig)=zlay(ig,nlay)-(pcon(ig)-pplay(ig,nlay))/(RG*rho(ig,nlay))/100.
           abort_message = 'thermcellV0_main: les thermiques vont trop haut '
           CALL abort_gcm (modname,abort_message,1)
        endif
      enddo
      if (prt_level.ge.1) print*,'14b OK convect8'
      do k=nlay,1,-1
         do ig=1,ngrid
            if (zqla(ig,k).gt.1e-10) then
               nivcon(ig)=k
               zcon(ig)=zlev(ig,k)
            endif
         enddo
      enddo
      if (prt_level.ge.1) print*,'14c OK convect8'
!calcul des moments
!initialisation
      do l=1,nlay
         do ig=1,ngrid
            q2(ig,l)=0.
            wth2(ig,l)=0.
            wth3(ig,l)=0.
            ratqscth(ig,l)=0.
            ratqsdiff(ig,l)=0.
         enddo
      enddo      
      if (prt_level.ge.1) print*,'14d OK convect8'
      if (prt_level.ge.10)write(lunout,*)                                &
    &     'WARNING thermcell_main wth2=0. si zw2 > 1.e-10'
      do l=1,nlay
         do ig=1,ngrid
            zf=fraca(ig,l)
            zf2=zf/(1.-zf)
!
      if (prt_level.ge.10) print*,'14e OK convect8 ig,l,zf,zf2',ig,l,zf,zf2
!
      if (prt_level.ge.10) print*,'14f OK convect8 ig,l,zha zh zpspsk ',ig,l,zha(ig,l),zh(ig,l),zpspsk(ig,l)
            thetath2(ig,l)=zf2*(ztla(ig,l)-zthl(ig,l))**2
            if(zw2(ig,l).gt.1.e-10) then
             wth2(ig,l)=zf2*(zw2(ig,l))**2
            else
             wth2(ig,l)=0.
            endif
!           print*,'wth2=',wth2(ig,l)
            wth3(ig,l)=zf2*(1-2.*fraca(ig,l))/(1-fraca(ig,l))  &
     &                *zw2(ig,l)*zw2(ig,l)*zw2(ig,l)
      if (prt_level.ge.10) print*,'14g OK convect8 ig,l,po',ig,l,po(ig,l)
            q2(ig,l)=zf2*(zqta(ig,l)*1000.-po(ig,l)*1000.)**2
!test: on calcul q2/po=ratqsc
            ratqscth(ig,l)=sqrt(max(q2(ig,l),1.e-6)/(po(ig,l)*1000.))
         enddo
      enddo
!calcul des flux: q, thetal et thetav
      do l=1,nlay
         do ig=1,ngrid
      wq(ig,l)=fraca(ig,l)*zw2(ig,l)*(zqta(ig,l)*1000.-po(ig,l)*1000.)
      wthl(ig,l)=fraca(ig,l)*zw2(ig,l)*(ztla(ig,l)-zthl(ig,l))
      wthv(ig,l)=fraca(ig,l)*zw2(ig,l)*(ztva(ig,l)-ztv(ig,l))
         enddo
      enddo
!
!      print*,'avant calcul ale et alp' 
!calcul de ALE et ALP pour la convection
      Alp_bl(:)=0.
      Ale_bl(:)=0.
!          print*,'ALE,ALP ,l,zw2(ig,l),Ale_bl(ig),Alp_bl(ig)'
      do l=1,nlay
      do ig=1,ngrid
           Alp_bl(ig)=max(Alp_bl(ig),0.5*rhobarz(ig,l)*wth3(ig,l) )
           Ale_bl(ig)=max(Ale_bl(ig),0.5*zw2(ig,l)**2)
!          print*,'ALE,ALP',l,zw2(ig,l),Ale_bl(ig),Alp_bl(ig)
      enddo
      enddo

!     print*,'AAAAAAA ',Alp_bl,Ale_bl,lmix


! TEST. IL FAUT REECRIRE LES ALE et ALP
!     Ale_bl(:)=0.5*wmax(:)*wmax(:)
!     Alp_bl(:)=0.1*wmax(:)*wmax(:)*wmax(:)

!test:calcul de la ponderation des couches pour KE
!initialisations
!      print*,'ponderation'
      do ig=1,ngrid
           fm_tot(ig)=0.
      enddo
       do ig=1,ngrid
        do k=1,klev
           wght_th(ig,k)=1.
        enddo
       enddo
       do ig=1,ngrid
!         lalim_conv(ig)=lmix_bis(ig)
!la hauteur de la couche alim_conv = hauteur couche alim_therm
         lalim_conv(ig)=lalim(ig)
!         zentr(ig)=zlev(ig,lalim(ig))
      enddo
      do ig=1,ngrid
        do k=1,lalim_conv(ig)
           fm_tot(ig)=fm_tot(ig)+fm(ig,k)
        enddo
      enddo
      do ig=1,ngrid
        do k=1,lalim_conv(ig)
           if (fm_tot(ig).gt.1.e-10) then
!           wght_th(ig,k)=fm(ig,k)/fm_tot(ig)
           endif
!on pondere chaque couche par a*
             if (alim_star(ig,k).gt.1.e-10) then
             wght_th(ig,k)=alim_star(ig,k)
             else
             wght_th(ig,k)=1.
             endif
        enddo
      enddo
!      print*,'apres wght_th'
!test pour prolonger la convection
      do ig=1,ngrid
!v1d  if ((alim_star(ig,1).lt.1.e-10).and.(therm)) then
      if ((alim_star(ig,1).lt.1.e-10)) then
      lalim_conv(ig)=1
      wght_th(ig,1)=1.
!      print*,'lalim_conv ok',lalim_conv(ig),wght_th(ig,1)
      endif
      enddo

!calcul du ratqscdiff
      if (prt_level.ge.1) print*,'14e OK convect8'
      var=0.
      vardiff=0.
      ratqsdiff(:,:)=0.
      do ig=1,ngrid
         do l=1,lalim(ig)
            var=var+alim_star(ig,l)*zqta(ig,l)*1000.
         enddo
      enddo
      if (prt_level.ge.1) print*,'14f OK convect8'
      do ig=1,ngrid
          do l=1,lalim(ig)
          zf=fraca(ig,l)
          zf2=zf/(1.-zf)
          vardiff=vardiff+alim_star(ig,l)  &
     &           *(zqta(ig,l)*1000.-var)**2
!         ratqsdiff=ratqsdiff+alim_star(ig,l)*
!     s          (zqta(ig,l)*1000.-po(ig,l)*1000.)**2
          enddo
      enddo
      if (prt_level.ge.1) print*,'14g OK convect8'
      do l=1,nlay
         do ig=1,ngrid
            ratqsdiff(ig,l)=sqrt(vardiff)/(po(ig,l)*1000.)   
!           write(11,*)'ratqsdiff=',ratqsdiff(ig,l)
         enddo
      enddo 
!--------------------------------------------------------------------    
!
!ecriture des fichiers sortie
!     print*,'15 OK convect8'

#ifdef wrgrads_thermcell
      if (prt_level.ge.1) print*,'thermcell_main sorties 3D'
#include "thermcell_out3d.h"
#endif

      endif

      if (prt_level.ge.1) print*,'thermcell_main FIN  OK'

      return
      end

!-----------------------------------------------------------------------------

      subroutine test_ltherm(klon,klev,pplev,pplay,long,seuil,ztv,po,ztva,zqla,f_star,zw2,comment)
      IMPLICIT NONE
#include "iniprint.h"

      integer i, k, klon,klev
      real pplev(klon,klev+1),pplay(klon,klev)
      real ztv(klon,klev)
      real po(klon,klev)
      real ztva(klon,klev)
      real zqla(klon,klev)
      real f_star(klon,klev)
      real zw2(klon,klev)
      integer long(klon)
      real seuil
      character*21 comment

      if (prt_level.ge.1) THEN
       print*,'WARNING !!! TEST ',comment
      endif
      return

!  test sur la hauteur des thermiques ...
         do i=1,klon
!IMtemp           if (pplay(i,long(i)).lt.seuil*pplev(i,1)) then
           if (prt_level.ge.10) then
               print*,'WARNING ',comment,' au point ',i,' K= ',long(i)
               print*,'  K  P(MB)  THV(K)     Qenv(g/kg)THVA        QLA(g/kg)   F*        W2'
               do k=1,klev
                  write(6,'(i3,7f10.3)') k,pplay(i,k),ztv(i,k),1000*po(i,k),ztva(i,k),1000*zqla(i,k),f_star(i,k),zw2(i,k)
               enddo
           endif
         enddo


      return
      end

