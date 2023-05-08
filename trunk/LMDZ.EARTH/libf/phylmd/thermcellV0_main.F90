!
! $Id$
!
      SUBROUTINE thermcellV0_main(itap,ngrid,nlay,ptimestep  &
     &                  ,pplay,pplev,pphi,debut  &
     &                  ,pu,pv,pt,po  &
     &                  ,pduadj,pdvadj,pdtadj,pdoadj  &
     &                  ,fm0,entr0,detr0,zqta,zqla,lmax  &
     &                  ,ratqscth,ratqsdiff,zqsatth  &
     &                  ,r_aspect,l_mix,tau_thermals &
     &                  ,Ale_bl,Alp_bl,lalim_conv,wght_th &
     &                  ,zmax0, f0,zw2,fraca)

      USE dimphy
      USE comgeomphy , ONLY:rlond,rlatd
      IMPLICIT NONE

!=======================================================================
!   Auteurs: Frederic Hourdin, Catherine Rio, Anne Mathieu
!   Version du 09.02.07
!   Calcul du transport vertical dans la couche limite en presence
!   de "thermiques" explicitement representes avec processus nuageux
!
!   Réécriture à partir d'un listing papier à Habas, le 14/02/00
!
!   le thermique est supposé homogène et dissipé par mélange avec
!   son environnement. la longueur l_mix contrôle l'efficacité du
!   mélange
!
!   Le calcul du transport des différentes espèces se fait en prenant
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
      real zmax(klon),zw2(klon,klev+1),ztva(klon,klev),zw_est(klon,klev+1)
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
    
      real ratqscth(klon,klev)
      real var
      real vardiff
      real ratqsdiff(klon,klev)

      logical sorties
      real rho(klon,klev),rhobarz(klon,klev),masse(klon,klev)
      real zpspsk(klon,klev)

      real wmax(klon)
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
      real alim_star_tot(klon),alim_star2(klon)
      real alim_star(klon,klev)
      real f(klon), f0(klon)
!FH/IM     save f0
      real zlevinter(klon)
      logical debut
       real seuil

! Declaration uniquement pour les sorties dans thermcell_out3d.
! Inutilise en 3D
      real wthl(klon,klev)
      real wthv(klon,klev)
      real wq(klon,klev)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

      character (len=20) :: modname='thermcellV0_main'
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
!Initialisation
!
     if (prt_level.ge.10)write(lunout,*)                                &
    &     'WARNING thermcell_main f0=max(f0,1.e-2)'
     do ig=1,klon
         f0(ig)=max(f0(ig),1.e-2)
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
!  1. alim_star est le profil vertical de l'alimentation à la base du
!     panache thermique, calculé à partir de la flotabilité de l'air sec
!  2. lmin et lalim sont les indices inferieurs et superieurs de alim_star
!------------------------------------------------------------------
!
      entr_star=0. ; detr_star=0. ; alim_star=0. ; alim_star_tot=0.
      CALL thermcellV0_init(ngrid,nlay,ztv,zlay,zlev,  &
     &                    lalim,lmin,alim_star,alim_star_tot,lev_out)

call testV0_ltherm(ngrid,nlay,pplev,pplay,lmin,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_init lmin  ')
call testV0_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_init lalim ')


      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_init'
      if (prt_level.ge.10) then
         write(lunout1,*) 'Dans thermcell_main 1'
         write(lunout1,*) 'lmin ',lmin(igout)
         write(lunout1,*) 'lalim ',lalim(igout)
         write(lunout1,*) ' ig l alim_star thetav'
         write(lunout1,'(i6,i4,2e15.5)') (igout,l,alim_star(igout,l) &
     &   ,ztv(igout,l),l=1,lalim(igout)+4)
      endif

!v1d      do ig=1,klon
!v1d     if (alim_star(ig,1).gt.1.e-10) then
!v1d     therm=.true.
!v1d     endif
!v1d      enddo
!-----------------------------------------------------------------------------
!  3. wmax_sec et zmax_sec sont les vitesses et altitudes maximum d'un
!     panache sec conservatif (e=d=0) alimente selon alim_star 
!     Il s'agit d'un calcul de type CAPE
!     zmax_sec est utilisé pour déterminer la géométrie du thermique.
!------------------------------------------------------------------------------
!
      CALL thermcellV0_dry(ngrid,nlay,zlev,pphi,ztv,alim_star,  &
     &                      lalim,lmin,zmax_sec,wmax_sec,lev_out)

call testV0_ltherm(ngrid,nlay,pplev,pplay,lmin,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_dry  lmin  ')
call testV0_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_dry  lalim ')

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_dry'
      if (prt_level.ge.10) then
         write(lunout1,*) 'Dans thermcell_main 1b'
         write(lunout1,*) 'lmin ',lmin(igout)
         write(lunout1,*) 'lalim ',lalim(igout)
         write(lunout1,*) ' ig l alim_star entr_star detr_star f_star '
         write(lunout1,'(i6,i4,e15.5)') (igout,l,alim_star(igout,l) &
     &    ,l=1,lalim(igout)+4)
      endif



!---------------------------------------------------------------------------------
!calcul du melange et des variables dans le thermique
!--------------------------------------------------------------------------------
!
      if (prt_level.ge.1) print*,'avant thermcell_plume ',lev_out
!IM 140508   CALL thermcell_plume(ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,  &
      CALL thermcellV0_plume(itap,ngrid,nlay,ptimestep,ztv,zthl,po,zl,rhobarz,  &
     &           zlev,pplev,pphi,zpspsk,l_mix,r_aspect,alim_star,alim_star_tot,  &
     &           lalim,zmax_sec,f0,detr_star,entr_star,f_star,ztva,  &
     &           ztla,zqla,zqta,zha,zw2,zw_est,zqsatth,lmix,lmix_bis,linter &
     &            ,lev_out,lunout1,igout)
      if (prt_level.ge.1) print*,'apres thermcell_plume ',lev_out

      call testV0_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_plum lalim ')
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lmix ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_plum lmix  ')

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


      call testV0_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lalim ')
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lmin ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmin  ')
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lmix ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmix  ')
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lmax ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_heig lmax  ')

      if (prt_level.ge.1) print*,'thermcell_main apres thermcell_height'

!-------------------------------------------------------------------------------
! Fermeture,determination de f
!-------------------------------------------------------------------------------
!
!avant closure: on redéfinit lalim, alim_star_tot et alim_star
!       do ig=1,klon
!       do l=2,lalim(ig)
!       alim_star(ig,l)=entr_star(ig,l)
!       entr_star(ig,l)=0.
!       enddo
!       enddo

      CALL thermcellV0_closure(ngrid,nlay,r_aspect,ptimestep,rho,  &
     &   zlev,lalim,alim_star,alim_star_tot,zmax_sec,wmax_sec,zmax,wmax,f,lev_out)

      if(prt_level.ge.1)print*,'thermcell_closure apres thermcell_closure'

      if (tau_thermals>1.) then
         lambda=exp(-ptimestep/tau_thermals)
         f0=(1.-lambda)*f+lambda*f0
      else
         f0=f
      endif

! Test valable seulement en 1D mais pas genant
      if (.not. (f0(1).ge.0.) ) then
        abort_message = 'Dans thermcell_main f0(1).lt.0 '
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
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lalim,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_flux lalim ')
      call testV0_ltherm(ngrid,nlay,pplev,pplay,lmax ,seuil,ztv,po,ztva,zqla,f_star,zw2,'thermcell_flux lmax  ')

!------------------------------------------------------------------
!   On ne prend pas directement les profils issus des calculs precedents
!   mais on s'autorise genereusement une relaxation vers ceci avec
!   une constante de temps tau_thermals (typiquement 1800s).
!------------------------------------------------------------------

      if (tau_thermals>1.) then
         lambda=exp(-ptimestep/tau_thermals)
         fm0=(1.-lambda)*fm+lambda*fm0
         entr0=(1.-lambda)*entr+lambda*entr0
!        detr0=(1.-lambda)*detr+lambda*detr0
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
!IM 050508    &    ,zu,zv,pduadj,pdvadj,zua,zva,igout,lev_out)
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
            thetath2(ig,l)=zf2*(zha(ig,l)-zh(ig,l)/zpspsk(ig,l))**2
            if(zw2(ig,l).gt.1.e-10) then
             wth2(ig,l)=zf2*(zw2(ig,l))**2
            else
             wth2(ig,l)=0.
            endif
!           print*,'wth2=',wth2(ig,l)
            wth3(ig,l)=zf2*(1-2.*fraca(ig,l))/(1-fraca(ig,l))  &
     &                *zw2(ig,l)*zw2(ig,l)*zw2(ig,l)
            q2(ig,l)=zf2*(zqta(ig,l)*1000.-po(ig,l)*1000.)**2
!test: on calcul q2/po=ratqsc
            ratqscth(ig,l)=sqrt(max(q2(ig,l),1.e-6)/(po(ig,l)*1000.))
         enddo
      enddo

      if (prt_level.ge.10) then
          print*,'14e OK convect8 ig,l,zf,zf2',ig,l,zf,zf2
          ig=igout
          do l=1,nlay
             print*,'14f OK convect8 ig,l,zha zh zpspsk ',ig,l,zha(ig,l),zh(ig,l),zpspsk(ig,l)
          enddo
          do l=1,nlay
             print*,'14g OK convect8 ig,l,po',ig,l,po(ig,l)
          enddo
      endif

      do ig=1,ngrid
      alp_int(ig)=0.
      ale_int(ig)=0.
      n_int(ig)=0
      enddo
!
      do l=1,nlay
      do ig=1,ngrid
       if(l.LE.lmax(ig)) THEN
        alp_int(ig)=alp_int(ig)+0.5*rhobarz(ig,l)*wth3(ig,l)
        ale_int(ig)=ale_int(ig)+0.5*zw2(ig,l)**2
        n_int(ig)=n_int(ig)+1
       endif
      enddo
      enddo
!      print*,'avant calcul ale et alp' 
!calcul de ALE et ALP pour la convection
      do ig=1,ngrid
!      Alp_bl(ig)=0.5*rhobarz(ig,lmix_bis(ig))*wth3(ig,lmix(ig))
!          Alp_bl(ig)=0.5*rhobarz(ig,nivcon(ig))*wth3(ig,nivcon(ig))
!      Alp_bl(ig)=0.5*rhobarz(ig,lmix(ig))*wth3(ig,lmix(ig)) 
!     &           *0.1
!valeur integree de alp_bl * 0.5:
       if (n_int(ig).gt.0) then
       Alp_bl(ig)=0.5*alp_int(ig)/n_int(ig)
!       if (Alp_bl(ig).lt.0.) then
!       Alp_bl(ig)=0.
       endif
!       endif
!         write(18,*),'rhobarz,wth3,Alp',rhobarz(ig,nivcon(ig)),
!     s               wth3(ig,nivcon(ig)),Alp_bl(ig)
!            write(18,*),'ALP_BL',Alp_bl(ig),lmix(ig)
!      Ale_bl(ig)=0.5*zw2(ig,lmix_bis(ig))**2
!      if (nivcon(ig).eq.1) then
!       Ale_bl(ig)=0.
!       else
!valeur max de ale_bl:
       Ale_bl(ig)=0.5*zw2(ig,lmix(ig))**2 
!     & /2.
!     & *0.1
!        Ale_bl(ig)=0.5*zw2(ig,lmix_bis(ig))**2 
!       if (n_int(ig).gt.0) then
!       Ale_bl(ig)=ale_int(ig)/n_int(ig)
!        Ale_bl(ig)=4.
!       endif
!       endif
!            Ale_bl(ig)=0.5*wth2(ig,lmix_bis(ig))
!          Ale_bl(ig)=wth2(ig,nivcon(ig))
!          write(19,*),'wth2,ALE_BL',wth2(ig,nivcon(ig)),Ale_bl(ig)
      enddo
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

      if (prt_level.ge.1) print*,'thermcell_main sorties 3D'
#ifdef wrgrads_thermcell
#include "thermcell_out3d.h"
#endif

      endif

      if (prt_level.ge.1) print*,'thermcell_main FIN  OK'

!     if(icount.eq.501) stop'au pas 301 dans thermcell_main'
      return
      end

!-----------------------------------------------------------------------------

      subroutine testV0_ltherm(klon,klev,pplev,pplay,long,seuil,ztv,po,ztva,zqla,f_star,zw2,comment)
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

!==============================================================================
      SUBROUTINE thermcellV0_closure(ngrid,nlay,r_aspect,ptimestep,rho,  &
     &   zlev,lalim,alim_star,alim_star_tot,zmax_sec,wmax_sec,zmax,wmax,f,lev_out)

!-------------------------------------------------------------------------
!thermcell_closure: fermeture, determination de f
!-------------------------------------------------------------------------
      IMPLICIT NONE

#include "iniprint.h"
#include "thermcell.h"
      INTEGER ngrid,nlay
      INTEGER ig,k       
      REAL r_aspect,ptimestep
      integer lev_out                           ! niveau pour les print

      INTEGER lalim(ngrid)
      REAL alim_star(ngrid,nlay)
      REAL alim_star_tot(ngrid)
      REAL rho(ngrid,nlay)
      REAL zlev(ngrid,nlay)
      REAL zmax(ngrid),zmax_sec(ngrid)
      REAL wmax(ngrid),wmax_sec(ngrid)
      real zdenom

      REAL alim_star2(ngrid)

      REAL f(ngrid)

      character (len=20) :: modname='thermcellV0_main'
      character (len=80) :: abort_message

      do ig=1,ngrid
         alim_star2(ig)=0.
      enddo
      do ig=1,ngrid
         if (alim_star(ig,1).LT.1.e-10) then
            f(ig)=0.
         else   
             do k=1,lalim(ig)
                alim_star2(ig)=alim_star2(ig)+alim_star(ig,k)**2  &
     &                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
             enddo
             zdenom=max(500.,zmax(ig))*r_aspect*alim_star2(ig)
             if (zdenom<1.e-14) then
                print*,'ig=',ig
                print*,'alim_star2',alim_star2(ig)
                print*,'zmax',zmax(ig)
                print*,'r_aspect',r_aspect
                print*,'zdenom',zdenom
                print*,'alim_star',alim_star(ig,:)
                print*,'zmax_sec',zmax_sec(ig)
                print*,'wmax_sec',wmax_sec(ig)
                abort_message = 'zdenom<1.e-14'
                CALL abort_gcm (modname,abort_message,1)
             endif
             if ((zmax_sec(ig).gt.1.e-10).and.(iflag_thermals_ed.eq.0)) then 
             f(ig)=wmax_sec(ig)*alim_star_tot(ig)/(max(500.,zmax_sec(ig))*r_aspect  &
     &             *alim_star2(ig))
!            f(ig)=f(ig)+(f0(ig)-f(ig))*exp((-ptimestep/  &
!    &                     zmax_sec(ig))*wmax_sec(ig))
             if(prt_level.GE.10) write(lunout,*)'closure dry',f(ig),wmax_sec(ig),alim_star_tot(ig),zmax_sec(ig)
             else
             f(ig)=wmax(ig)*alim_star_tot(ig)/zdenom
!            f(ig)=f(ig)+(f0(ig)-f(ig))*exp((-ptimestep/  &
!     &                     zmax(ig))*wmax(ig))
             if(prt_level.GE.10) print*,'closure moist',f(ig),wmax(ig),alim_star_tot(ig),zmax(ig)
             endif
         endif
!         f0(ig)=f(ig)
      enddo
      if (prt_level.ge.1) print*,'apres fermeture'

!
      return
      end
!==============================================================================
      SUBROUTINE thermcellV0_plume(itap,ngrid,klev,ptimestep,ztv,zthl,po,zl,rhobarz,  &
     &           zlev,pplev,pphi,zpspsk,l_mix,r_aspect,alim_star,alim_star_tot,  &
     &           lalim,zmax_sec,f0,detr_star,entr_star,f_star,ztva,  &
     &           ztla,zqla,zqta,zha,zw2,w_est,zqsatth,lmix,lmix_bis,linter &
     &           ,lev_out,lunout1,igout)

!--------------------------------------------------------------------------
!thermcell_plume: calcule les valeurs de qt, thetal et w dans l ascendance
!--------------------------------------------------------------------------

      IMPLICIT NONE

#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
#include "iniprint.h"
#include "thermcell.h"

      INTEGER itap
      INTEGER lunout1,igout
      INTEGER ngrid,klev
      REAL ptimestep
      REAL ztv(ngrid,klev)
      REAL zthl(ngrid,klev)
      REAL po(ngrid,klev)
      REAL zl(ngrid,klev)
      REAL rhobarz(ngrid,klev)
      REAL zlev(ngrid,klev+1)
      REAL pplev(ngrid,klev+1)
      REAL pphi(ngrid,klev)
      REAL zpspsk(ngrid,klev)
      REAL alim_star(ngrid,klev)
      REAL zmax_sec(ngrid)
      REAL f0(ngrid)
      REAL l_mix
      REAL r_aspect
      INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print
      real zcon2(ngrid)
    
      real alim_star_tot(ngrid)

      REAL ztva(ngrid,klev)
      REAL ztla(ngrid,klev)
      REAL zqla(ngrid,klev)
      REAL zqla0(ngrid,klev)
      REAL zqta(ngrid,klev)
      REAL zha(ngrid,klev)

      REAL detr_star(ngrid,klev)
      REAL coefc
      REAL detr_stara(ngrid,klev)
      REAL detr_starb(ngrid,klev)
      REAL detr_starc(ngrid,klev)
      REAL detr_star0(ngrid,klev)
      REAL detr_star1(ngrid,klev)
      REAL detr_star2(ngrid,klev)

      REAL entr_star(ngrid,klev)
      REAL entr_star1(ngrid,klev)
      REAL entr_star2(ngrid,klev)
      REAL detr(ngrid,klev)
      REAL entr(ngrid,klev)

      REAL zw2(ngrid,klev+1)
      REAL w_est(ngrid,klev+1)
      REAL f_star(ngrid,klev+1)
      REAL wa_moy(ngrid,klev+1)

      REAL ztva_est(ngrid,klev)
      REAL zqla_est(ngrid,klev)
      REAL zqsatth(ngrid,klev)
      REAL zta_est(ngrid,klev)

      REAL linter(ngrid)
      INTEGER lmix(ngrid)
      INTEGER lmix_bis(ngrid)
      REAL    wmaxa(ngrid)

      INTEGER ig,l,k

      real zcor,zdelta,zcvm5,qlbef
      real Tbef,qsatbef
      real dqsat_dT,DT,num,denom
      REAL REPS,RLvCp,DDT0
      PARAMETER (DDT0=.01)
      logical Zsat
      REAL fact_gamma,fact_epsilon
      REAL c2(ngrid,klev)

      Zsat=.false.
! Initialisation
      RLvCp = RLVTT/RCPD
     
      if (iflag_thermals_ed==0) then
         fact_gamma=1.
         fact_epsilon=1.
      else if (iflag_thermals_ed==1)  then
         fact_gamma=1.
         fact_epsilon=1.
      else if (iflag_thermals_ed==2)  then
         fact_gamma=1.
         fact_epsilon=2.
      endif

      do l=1,klev
         do ig=1,ngrid
            zqla_est(ig,l)=0.
            ztva_est(ig,l)=ztva(ig,l)
            zqsatth(ig,l)=0.
         enddo
      enddo

!CR: attention test couche alim
!     do l=2,klev
!     do ig=1,ngrid
!        alim_star(ig,l)=0.
!     enddo
!     enddo
!AM:initialisations du thermique
      do k=1,klev
         do ig=1,ngrid
            ztva(ig,k)=ztv(ig,k)
            ztla(ig,k)=zthl(ig,k)
            zqla(ig,k)=0.
            zqta(ig,k)=po(ig,k)
!
            ztva(ig,k) = ztla(ig,k)*zpspsk(ig,k)+RLvCp*zqla(ig,k)
            ztva(ig,k) = ztva(ig,k)/zpspsk(ig,k)
            zha(ig,k) = ztva(ig,k)
!
         enddo
      enddo 
      do k=1,klev
        do ig=1,ngrid
           detr_star(ig,k)=0.
           entr_star(ig,k)=0.

           detr_stara(ig,k)=0.
           detr_starb(ig,k)=0.
           detr_starc(ig,k)=0.
           detr_star0(ig,k)=0.
           zqla0(ig,k)=0.
           detr_star1(ig,k)=0.
           detr_star2(ig,k)=0.
           entr_star1(ig,k)=0.
           entr_star2(ig,k)=0.

           detr(ig,k)=0.
           entr(ig,k)=0.
        enddo
      enddo
      if (prt_level.ge.1) print*,'7 OK convect8'
      do k=1,klev+1
         do ig=1,ngrid
            zw2(ig,k)=0.
            w_est(ig,k)=0.
            f_star(ig,k)=0.
            wa_moy(ig,k)=0.
         enddo
      enddo

      if (prt_level.ge.1) print*,'8 OK convect8'
      do ig=1,ngrid
         linter(ig)=1.
         lmix(ig)=1
         lmix_bis(ig)=2
         wmaxa(ig)=0.
      enddo

!-----------------------------------------------------------------------------------
!boucle de calcul de la vitesse verticale dans le thermique
!-----------------------------------------------------------------------------------
      do l=1,klev-1
         do ig=1,ngrid



! Calcul dans la premiere couche active du thermique (ce qu'on teste
! en disant que la couche est instable et que w2 en bas de la couche
! est nulle.

            if (ztv(ig,l).gt.ztv(ig,l+1)  &
     &         .and.alim_star(ig,l).gt.1.e-10  &
     &         .and.zw2(ig,l).lt.1e-10) then


! Le panache va prendre au debut les caracteristiques de l'air contenu
! dans cette couche.
               ztla(ig,l)=zthl(ig,l) 
               zqta(ig,l)=po(ig,l)
               zqla(ig,l)=zl(ig,l)
               f_star(ig,l+1)=alim_star(ig,l)

               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)  &
     &                     *(zlev(ig,l+1)-zlev(ig,l))  &
     &                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))
               w_est(ig,l+1)=zw2(ig,l+1)
!


            else if ((zw2(ig,l).ge.1e-10).and.  &
     &         (f_star(ig,l)+alim_star(ig,l)).gt.1.e-10) then
!estimation du detrainement a partir de la geometrie du pas precedent
!tests sur la definition du detr
!calcul de detr_star et entr_star



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH le test miraculeux de Catherine ? Le bout du tunel ?
!               w_est(ig,3)=zw2(ig,2)*  &
!    &                   ((f_star(ig,2))**2)  &
!    &                   /(f_star(ig,2)+alim_star(ig,2))**2+  &
!    &                   2.*RG*(ztva(ig,1)-ztv(ig,2))/ztv(ig,2)  &
!    &                   *(zlev(ig,3)-zlev(ig,2))
!     if (l.gt.2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! Premier calcul de la vitesse verticale a partir de la temperature
! potentielle virtuelle

! FH CESTQUOI CA ????
#define int1d2
!#undef int1d2
#ifdef int1d2
      if (l.ge.2) then
#else
      if (l.gt.2) then
#endif

      if (1.eq.1) then
          w_est(ig,3)=zw2(ig,2)* &
     &      ((f_star(ig,2))**2) &
     &      /(f_star(ig,2)+alim_star(ig,2))**2+ &
     &      2.*RG*(ztva(ig,2)-ztv(ig,2))/ztv(ig,2) &
!     &      *1./3. &
     &      *(zlev(ig,3)-zlev(ig,2))
       endif


!---------------------------------------------------------------------------
!calcul de l entrainement et du detrainement lateral
!---------------------------------------------------------------------------
!
!test:estimation de ztva_new_est sans entrainement

               Tbef=ztla(ig,l-1)*zpspsk(ig,l)
               zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
               qsatbef= R2ES * FOEEW(Tbef,zdelta)/pplev(ig,l)
               qsatbef=MIN(0.5,qsatbef)
               zcor=1./(1.-retv*qsatbef)
               qsatbef=qsatbef*zcor
               Zsat = (max(0.,zqta(ig,l-1)-qsatbef) .gt. 1.e-10)
               if (Zsat) then
               qlbef=max(0.,zqta(ig,l-1)-qsatbef)
               DT = 0.5*RLvCp*qlbef
               do while (abs(DT).gt.DDT0)
                 Tbef=Tbef+DT
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
                 qsatbef= R2ES * FOEEW(Tbef,zdelta)/pplev(ig,l)
                 qsatbef=MIN(0.5,qsatbef)
                 zcor=1./(1.-retv*qsatbef)
                 qsatbef=qsatbef*zcor
                 qlbef=zqta(ig,l-1)-qsatbef

                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
                 zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
                 zcor=1./(1.-retv*qsatbef)
                 dqsat_dT=FOEDE(Tbef,zdelta,zcvm5,qsatbef,zcor)
                 num=-Tbef+ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*qlbef
                 denom=1.+RLvCp*dqsat_dT
                 DT=num/denom
               enddo
                 zqla_est(ig,l) = max(0.,zqta(ig,l-1)-qsatbef) 
               endif
        ztva_est(ig,l) = ztla(ig,l-1)*zpspsk(ig,l)+RLvCp*zqla_est(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)/zpspsk(ig,l)
        zta_est(ig,l)=ztva_est(ig,l)
        ztva_est(ig,l) = ztva_est(ig,l)*(1.+RETV*(zqta(ig,l-1)  &
     &      -zqla_est(ig,l))-zqla_est(ig,l))

             w_est(ig,l+1)=zw2(ig,l)*  &
     &                   ((f_star(ig,l))**2)  &
     &                   /(f_star(ig,l)+alim_star(ig,l))**2+  &
     &                   2.*RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)  &
!     &                   *1./3. &
     &                   *(zlev(ig,l+1)-zlev(ig,l))
             if (w_est(ig,l+1).lt.0.) then
                w_est(ig,l+1)=zw2(ig,l)
             endif
!
!calcul du detrainement
!=======================

!CR:on vire les modifs
         if (iflag_thermals_ed==0) then

! Modifications du calcul du detrainement.
! Dans la version de la these de Catherine, on passe brusquement
! de la version seche a la version nuageuse pour le detrainement
! ce qui peut occasioner des oscillations.
! dans la nouvelle version, on commence par calculer un detrainement sec.
! Puis un autre en cas de nuages.
! Puis on combine les deux lineairement en fonction de la quantite d'eau.

#define int1d3
!#undef int1d3
#define RIO_TH
#ifdef RIO_TH
!1. Cas non nuageux
! 1.1 on est sous le zmax_sec et w croit
          if ((w_est(ig,l+1).gt.w_est(ig,l)).and.  &
     &       (zlev(ig,l+1).lt.zmax_sec(ig)).and.  &
#ifdef int1d3
     &       (zqla_est(ig,l).lt.1.e-10)) then 
#else
     &       (zqla(ig,l-1).lt.1.e-10)) then 
#endif
             detr_star(ig,l)=MAX(0.,(rhobarz(ig,l+1)  &
     &       *sqrt(w_est(ig,l+1))*sqrt(l_mix*zlev(ig,l+1))  &
     &       -rhobarz(ig,l)*sqrt(w_est(ig,l))*sqrt(l_mix*zlev(ig,l)))  &
     &       /(r_aspect*zmax_sec(ig)))
             detr_stara(ig,l)=detr_star(ig,l)

       if (prt_level.ge.20) print*,'coucou calcul detr 1: ig, l',ig,l

! 1.2 on est sous le zmax_sec et w decroit
          else if ((zlev(ig,l+1).lt.zmax_sec(ig)).and.  &
#ifdef int1d3
     &            (zqla_est(ig,l).lt.1.e-10)) then
#else
     &            (zqla(ig,l-1).lt.1.e-10)) then
#endif
             detr_star(ig,l)=-f0(ig)*f_star(ig,lmix(ig))  &
     &       /(rhobarz(ig,lmix(ig))*wmaxa(ig))*  &
     &       (rhobarz(ig,l+1)*sqrt(w_est(ig,l+1))  &
     &       *((zmax_sec(ig)-zlev(ig,l+1))/  &
     &       ((zmax_sec(ig)-zlev(ig,lmix(ig)))))**2.  &
     &       -rhobarz(ig,l)*sqrt(w_est(ig,l))  &
     &       *((zmax_sec(ig)-zlev(ig,l))/  &
     &       ((zmax_sec(ig)-zlev(ig,lmix(ig)))))**2.)
             detr_starb(ig,l)=detr_star(ig,l)

        if (prt_level.ge.20) print*,'coucou calcul detr 2: ig, l',ig,l

          else

! 1.3 dans les autres cas
             detr_star(ig,l)=0.002*f0(ig)*f_star(ig,l)  &
     &                      *(zlev(ig,l+1)-zlev(ig,l))
             detr_starc(ig,l)=detr_star(ig,l)

        if (prt_level.ge.20) print*,'coucou calcul detr 3 n: ig, l',ig, l
             
          endif

#else

! 1.1 on est sous le zmax_sec et w croit
          if ((w_est(ig,l+1).gt.w_est(ig,l)).and.  &
     &       (zlev(ig,l+1).lt.zmax_sec(ig)) ) then
             detr_star(ig,l)=MAX(0.,(rhobarz(ig,l+1)  &
     &       *sqrt(w_est(ig,l+1))*sqrt(l_mix*zlev(ig,l+1))  &
     &       -rhobarz(ig,l)*sqrt(w_est(ig,l))*sqrt(l_mix*zlev(ig,l)))  &
     &       /(r_aspect*zmax_sec(ig)))

       if (prt_level.ge.20) print*,'coucou calcul detr 1: ig, l', ig, l

! 1.2 on est sous le zmax_sec et w decroit
          else if ((zlev(ig,l+1).lt.zmax_sec(ig)) ) then
             detr_star(ig,l)=-f0(ig)*f_star(ig,lmix(ig))  &
     &       /(rhobarz(ig,lmix(ig))*wmaxa(ig))*  &
     &       (rhobarz(ig,l+1)*sqrt(w_est(ig,l+1))  &
     &       *((zmax_sec(ig)-zlev(ig,l+1))/  &
     &       ((zmax_sec(ig)-zlev(ig,lmix(ig)))))**2.  &
     &       -rhobarz(ig,l)*sqrt(w_est(ig,l))  &
     &       *((zmax_sec(ig)-zlev(ig,l))/  &
     &       ((zmax_sec(ig)-zlev(ig,lmix(ig)))))**2.)
       if (prt_level.ge.20) print*,'coucou calcul detr 1: ig, l', ig, l

          else
             detr_star=0.
          endif

! 1.3 dans les autres cas
          detr_starc(ig,l)=0.002*f0(ig)*f_star(ig,l)  &
     &                      *(zlev(ig,l+1)-zlev(ig,l))

          coefc=min(zqla(ig,l-1)/1.e-3,1.)
          if (zlev(ig,l+1).ge.zmax_sec(ig)) coefc=1.
          coefc=1.
! il semble qu'il soit important de baser le calcul sur
! zqla_est(ig,l-1) plutot que sur zqla_est(ig,l)
          detr_star(ig,l)=detr_starc(ig,l)*coefc+detr_star(ig,l)*(1.-coefc)

        if (prt_level.ge.20) print*,'coucou calcul detr 2: ig, l', ig, l

#endif


        if (prt_level.ge.20) print*,'coucou calcul detr 444: ig, l', ig, l
!IM 730508 beg
!        if(itap.GE.7200) THEN
!         print*,'th_plume ig,l,itap,zqla_est=',ig,l,itap,zqla_est(ig,l)
!        endif
!IM 730508 end
         
         zqla0(ig,l)=zqla_est(ig,l)
         detr_star0(ig,l)=detr_star(ig,l)
!IM 060508 beg
!         if(detr_star(ig,l).GT.1.) THEN
!          print*,'th_plumeBEF ig l detr_star detr_starc coefc',ig,l,detr_star(ig,l) &
!   &      ,detr_starc(ig,l),coefc
!         endif
!IM 060508 end
!IM 160508 beg
!IM 160508       IF (f0(ig).NE.0.) THEN
           detr_star(ig,l)=detr_star(ig,l)/f0(ig)
!IM 160508       ELSE IF(detr_star(ig,l).EQ.0.) THEN
!IM 160508        print*,'WARNING1  : th_plume f0=0, detr_star=0: ig, l, itap',ig,l,itap
!IM 160508       ELSE
!IM 160508        print*,'WARNING2  : th_plume f0=0, ig, l, itap, detr_star',ig,l,itap,detr_star(ig,l)
!IM 160508       ENDIF
!IM 160508 end
!IM 060508 beg
!        if(detr_star(ig,l).GT.1.) THEN
!         print*,'th_plumeAFT ig l detr_star f0 1/f0',ig,l,detr_star(ig,l),f0(ig), &
!   &     REAL(1)/f0(ig)
!        endif
!IM 060508 end
        if (prt_level.ge.20) print*,'coucou calcul detr 445: ig, l', ig, l
!
!calcul de entr_star

! #undef test2
! #ifdef test2
! La version test2 destabilise beaucoup le modele.
! Il semble donc que ca aide d'avoir un entrainement important sous
! le nuage.
!         if (zqla_est(ig,l-1).ge.1.e-10.and.l.gt.lalim(ig)) then
!          entr_star(ig,l)=0.4*detr_star(ig,l)
!         else
!          entr_star(ig,l)=0.
!         endif
! #else
!
! Deplacement du calcul de entr_star pour eviter d'avoir aussi
! entr_star > fstar.
! Redeplacer suite a la transformation du cas detr>f
! FH

        if (prt_level.ge.20) print*,'coucou calcul detr 446: ig, l', ig, l
#define int1d
!FH 070508 #define int1d4
!#undef int1d4
! L'option int1d4 correspond au choix dans le cas ou le detrainement
! devient trop grand.

#ifdef int1d

#ifdef int1d4
#else
       detr_star(ig,l)=min(detr_star(ig,l),f_star(ig,l))
!FH 070508 plus
       detr_star(ig,l)=min(detr_star(ig,l),1.)
#endif

       entr_star(ig,l)=max(0.4*detr_star(ig,l)-alim_star(ig,l),0.)

        if (prt_level.ge.20) print*,'coucou calcul detr 447: ig, l', ig, l
#ifdef int1d4
! Si le detrainement excede le flux en bas + l'entrainement, le thermique
! doit disparaitre.
       if (detr_star(ig,l)>f_star(ig,l)+entr_star(ig,l)) then
          detr_star(ig,l)=f_star(ig,l)+entr_star(ig,l)
          f_star(ig,l+1)=0.
          linter(ig)=l+1
          zw2(ig,l+1)=-1.e-10
       endif
#endif


#else

        if (prt_level.ge.20) print*,'coucou calcul detr 448: ig, l', ig, l
        if(l.gt.lalim(ig)) then
         entr_star(ig,l)=0.4*detr_star(ig,l)
        else

! FH :
! Cette ligne doit permettre de garantir qu'on a toujours un flux = 1
! en haut de la couche d'alimentation.
! A remettre en questoin a la premiere occasion mais ca peut aider a 
! ecrire un code robuste.
! Que ce soit avec ca ou avec l'ancienne facon de faire (e* = 0 mais
! d* non nul) on a une discontinuité de e* ou d* en haut de la couche
! d'alimentation, ce qui n'est pas forcement heureux.

        if (prt_level.ge.20) print*,'coucou calcul detr 449: ig, l', ig, l
#undef pre_int1c
#ifdef pre_int1c
         entr_star(ig,l)=max(detr_star(ig,l)-alim_star(ig,l),0.)
         detr_star(ig,l)=entr_star(ig,l)
#else
         entr_star(ig,l)=0.
#endif

        endif

#endif

        if (prt_level.ge.20) print*,'coucou calcul detr 440: ig, l', ig, l
        entr_star1(ig,l)=entr_star(ig,l)
        detr_star1(ig,l)=detr_star(ig,l)
!

#ifdef int1d
#else
        if (detr_star(ig,l).gt.f_star(ig,l)) then

!  Ce test est là entre autres parce qu'on passe par des valeurs
!  delirantes de detr_star.
!  ca vaut sans doute le coup de verifier pourquoi.

           detr_star(ig,l)=f_star(ig,l)
#ifdef pre_int1c
           if (l.gt.lalim(ig)+1) then
               entr_star(ig,l)=0.
               alim_star(ig,l)=0.
! FH ajout pour forcer a stoper le thermique juste sous le sommet
! de la couche (voir calcul de finter)
               zw2(ig,l+1)=-1.e-10
               linter(ig)=l+1
            else
               entr_star(ig,l)=0.4*detr_star(ig,l)
            endif
#else
           entr_star(ig,l)=0.4*detr_star(ig,l)
#endif
        endif
#endif

      else !l > 2
         detr_star(ig,l)=0.
         entr_star(ig,l)=0.
      endif

        entr_star2(ig,l)=entr_star(ig,l)
        detr_star2(ig,l)=detr_star(ig,l)
        if (prt_level.ge.20) print*,'coucou calcul detr 450: ig, l', ig, l

       endif  ! iflag_thermals_ed==0

!CR:nvlle def de entr_star et detr_star
      if (iflag_thermals_ed>=1) then
!      if (l.lt.lalim(ig)) then
!      if (l.lt.2) then 
!        entr_star(ig,l)=0.
!        detr_star(ig,l)=0.
!      else
!      if (0.001.gt.(RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l))/(2.*w_est(ig,l+1)))) then 
!         entr_star(ig,l)=0.001*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
!      else
!         entr_star(ig,l)=  &
!     &                f_star(ig,l)/(2.*w_est(ig,l+1))        &
!     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l))   &
!     &                *(zlev(ig,l+1)-zlev(ig,l))

 
         entr_star(ig,l)=MAX(0.*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l)),  &          
     &                f_star(ig,l)/(2.*w_est(ig,l+1))        &
     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)   &
     &                *(zlev(ig,l+1)-zlev(ig,l))) &
     &                +0.0001*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l))

        if (((ztva_est(ig,l)-ztv(ig,l)).gt.1.e-10).and.(l.le.lmix_bis(ig))) then
            alim_star_tot(ig)=alim_star_tot(ig)+entr_star(ig,l)
            lalim(ig)=lmix_bis(ig)
            if(prt_level.GE.10) print*,'alim_star_tot',alim_star_tot(ig),entr_star(ig,l)
        endif

        if (((ztva_est(ig,l)-ztv(ig,l)).gt.1.e-10).and.(l.le.lmix_bis(ig))) then
!        c2(ig,l)=2500000.*zqla_est(ig,l)/(1004.*zta_est(ig,l))
         c2(ig,l)=0.001
         detr_star(ig,l)=MAX(0.*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l)),  &
     &                c2(ig,l)*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l)) &
     &                -f_star(ig,l)/(2.*w_est(ig,l+1))       &
     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)       &
     &                *(zlev(ig,l+1)-zlev(ig,l)))                    &
     &                +0.0001*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l))

       else
!         c2(ig,l)=2500000.*zqla_est(ig,l)/(1004.*zta_est(ig,l))
          c2(ig,l)=0.003

         detr_star(ig,l)=MAX(0.*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l)),  &
     &                c2(ig,l)*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l)) &
     &                -f_star(ig,l)/(2.*w_est(ig,l+1))       &
     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)       &
     &                *(zlev(ig,l+1)-zlev(ig,l))) &
     &                +0.0002*f_star(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
       endif
         
           
!        detr_star(ig,l)=detr_star(ig,l)*3.
!        if (l.lt.lalim(ig)) then
!          entr_star(ig,l)=0.
!        endif
!        if (l.lt.2) then
!          entr_star(ig,l)=0.
!          detr_star(ig,l)=0.
!        endif


!      endif 
!      else if ((ztva_est(ig,l)-ztv(ig,l)).gt.1.e-10) then
!      entr_star(ig,l)=MAX(0.,0.8*f_star(ig,l)/(2.*w_est(ig,l+1))        &
!     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l))   &
!     &                *(zlev(ig,l+1)-zlev(ig,l))
!      detr_star(ig,l)=0.002*f_star(ig,l)                         &
!     &                *(zlev(ig,l+1)-zlev(ig,l))
!      else
!      entr_star(ig,l)=0.001*f_star(ig,l)                         &
!     &                *(zlev(ig,l+1)-zlev(ig,l))
!      detr_star(ig,l)=MAX(0.,-0.2*f_star(ig,l)/(2.*w_est(ig,l+1))       &
!     &                *RG*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l))       &
!     &                *(zlev(ig,l+1)-zlev(ig,l))                      &
!     &                +0.002*f_star(ig,l)                             &
!     &                *(zlev(ig,l+1)-zlev(ig,l))
!      endif

      endif   ! iflag_thermals_ed==1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH inutile si on conserve comme on l'a fait plus haut entr=detr
! dans la couche d'alimentation
!pas d entrainement dans la couche alim
!      if ((l.le.lalim(ig))) then
!           entr_star(ig,l)=0.
!      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!prise en compte du detrainement et de l entrainement dans le calcul du flux

      f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

!test sur le signe de f_star
        if (prt_level.ge.20) print*,'coucou calcul detr 451: ig, l', ig, l
       if (f_star(ig,l+1).gt.1.e-10) then 
!----------------------------------------------------------------------------
!calcul de la vitesse verticale en melangeant Tl et qt du thermique
!---------------------------------------------------------------------------
!
       Zsat=.false.
       ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*zthl(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
!
       zqta(ig,l)=(f_star(ig,l)*zqta(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*po(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
!  
               Tbef=ztla(ig,l)*zpspsk(ig,l)
               zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
               qsatbef= R2ES * FOEEW(Tbef,zdelta)/pplev(ig,l)               
               qsatbef=MIN(0.5,qsatbef)
               zcor=1./(1.-retv*qsatbef)
               qsatbef=qsatbef*zcor
               Zsat = (max(0.,zqta(ig,l)-qsatbef) .gt. 1.e-10)
               if (Zsat) then
               qlbef=max(0.,zqta(ig,l)-qsatbef)
               DT = 0.5*RLvCp*qlbef
               do while (abs(DT).gt.DDT0)
                 Tbef=Tbef+DT
                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
                 qsatbef= R2ES * FOEEW(Tbef,zdelta)/pplev(ig,l)
                 qsatbef=MIN(0.5,qsatbef)
                 zcor=1./(1.-retv*qsatbef)
                 qsatbef=qsatbef*zcor
                 qlbef=zqta(ig,l)-qsatbef

                 zdelta=MAX(0.,SIGN(1.,RTT-Tbef))
                 zcvm5=R5LES*(1.-zdelta) + R5IES*zdelta
                 zcor=1./(1.-retv*qsatbef)
                 dqsat_dT=FOEDE(Tbef,zdelta,zcvm5,qsatbef,zcor)
                 num=-Tbef+ztla(ig,l)*zpspsk(ig,l)+RLvCp*qlbef
                 denom=1.+RLvCp*dqsat_dT
                 DT=num/denom
              enddo
                 zqla(ig,l) = max(0.,qlbef) 
              endif
!    
        if (prt_level.ge.20) print*,'coucou calcul detr 4512: ig, l', ig, l
! on ecrit de maniere conservative (sat ou non)
!          T = Tl +Lv/Cp ql
           ztva(ig,l) = ztla(ig,l)*zpspsk(ig,l)+RLvCp*zqla(ig,l)
           ztva(ig,l) = ztva(ig,l)/zpspsk(ig,l)
!on rajoute le calcul de zha pour diagnostiques (temp potentielle)
           zha(ig,l) = ztva(ig,l)
           ztva(ig,l) = ztva(ig,l)*(1.+RETV*(zqta(ig,l)  &
     &              -zqla(ig,l))-zqla(ig,l))

!on ecrit zqsat 
           zqsatth(ig,l)=qsatbef  
!calcul de vitesse
           zw2(ig,l+1)=zw2(ig,l)*  &
     &                 ((f_star(ig,l))**2)  &
!  Tests de Catherine
!     &                 /(f_star(ig,l+1)+detr_star(ig,l))**2+             &
     &      /(f_star(ig,l+1)+detr_star(ig,l)-entr_star(ig,l)*(1.-fact_epsilon))**2+ &
     &                 2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
     &                 *fact_gamma &
     &                 *(zlev(ig,l+1)-zlev(ig,l))
!prise en compte des forces de pression que qd flottabilité<0
!              zw2(ig,l+1)=zw2(ig,l)*  &
!     &            1./(1.+2.*entr_star(ig,l)/f_star(ig,l)) + &        
!     &                 (f_star(ig,l))**2 &
!     &                 /(f_star(ig,l)+entr_star(ig,l))**2+ &
!     &                 (f_star(ig,l)-2.*entr_star(ig,l))**2/(f_star(ig,l)+2.*entr_star(ig,l))**2+  &        
!     &                 /(f_star(ig,l+1)+detr_star(ig,l)-entr_star(ig,l)*(1.-2.))**2+ &
!     &                 /(f_star(ig,l)**2+2.*2.*detr_star(ig,l)*f_star(ig,l)+2.*entr_star(ig,l)*f_star(ig,l))+ &
!     &                 2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
!     &                 *1./3. &
!     &                 *(zlev(ig,l+1)-zlev(ig,l))
          
!        write(30,*),l+1,zw2(ig,l+1)-zw2(ig,l), &
!     &              -2.*entr_star(ig,l)/f_star(ig,l)*zw2(ig,l), &
!     &               2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)*(zlev(ig,l+1)-zlev(ig,l))

 
!             zw2(ig,l+1)=zw2(ig,l)*  &
!     &                 (2.-2.*entr_star(ig,l)/f_star(ig,l)) &  
!     &                 -zw2(ig,l-1)+  &        
!     &                 2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
!     &                 *1./3. &
!     &                 *(zlev(ig,l+1)-zlev(ig,l))             

            endif
        endif
        if (prt_level.ge.20) print*,'coucou calcul detr 460: ig, l',ig, l
!
!initialisations pour le calcul de la hauteur du thermique, de l'inversion et de la vitesse verticale max 

            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
                print*,'On tombe sur le cas particulier de thermcell_plume'
                zw2(ig,l+1)=0.
                linter(ig)=l+1
            endif

!        if ((zw2(ig,l).gt.0.).and. (zw2(ig,l+1).le.0.)) then
        if (zw2(ig,l+1).lt.0.) then 
           linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &               -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
           zw2(ig,l+1)=0.
        endif

           wa_moy(ig,l+1)=sqrt(zw2(ig,l+1)) 

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
!on rajoute le calcul de lmix_bis
            if (zqla(ig,l).lt.1.e-10) then
               lmix_bis(ig)=l+1
            endif
            lmix(ig)=l+1
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
        enddo
      enddo

!on remplace a* par e* ds premiere couche
!      if (iflag_thermals_ed.ge.1) then
!       do ig=1,ngrid
!       do l=2,klev
!          if (l.lt.lalim(ig)) then
!             alim_star(ig,l)=entr_star(ig,l)
!          endif
!       enddo
!       enddo
!       do ig=1,ngrid
!          lalim(ig)=lmix_bis(ig)
!       enddo
!      endif
       if (iflag_thermals_ed.ge.1) then
          do ig=1,ngrid
             do l=2,lalim(ig)
                alim_star(ig,l)=entr_star(ig,l)
                entr_star(ig,l)=0.
             enddo
           enddo
       endif
        if (prt_level.ge.20) print*,'coucou calcul detr 470: ig, l', ig, l

!     print*,'thermcell_plume OK'

      return 
      end
!==============================================================================
       SUBROUTINE thermcellV0_dry(ngrid,nlay,zlev,pphi,ztv,alim_star,  &
     &                            lalim,lmin,zmax,wmax,lev_out)

!--------------------------------------------------------------------------
!thermcell_dry: calcul de zmax et wmax du thermique sec
!--------------------------------------------------------------------------
       IMPLICIT NONE
#include "YOMCST.h"       
#include "iniprint.h"
       INTEGER l,ig

       INTEGER ngrid,nlay
       REAL zlev(ngrid,nlay+1)
       REAL pphi(ngrid,nlay)
       REAl ztv(ngrid,nlay)
       REAL alim_star(ngrid,nlay)
       INTEGER lalim(ngrid)
      integer lev_out                           ! niveau pour les print

       REAL zmax(ngrid)
       REAL wmax(ngrid)

!variables locales
       REAL zw2(ngrid,nlay+1)
       REAL f_star(ngrid,nlay+1)
       REAL ztva(ngrid,nlay+1)
       REAL wmaxa(ngrid)
       REAL wa_moy(ngrid,nlay+1)
       REAL linter(ngrid),zlevinter(ngrid)
       INTEGER lmix(ngrid),lmax(ngrid),lmin(ngrid)

!initialisations
       do ig=1,ngrid
          do l=1,nlay+1
             zw2(ig,l)=0.
             wa_moy(ig,l)=0.
          enddo
       enddo
       do ig=1,ngrid
          do l=1,nlay
             ztva(ig,l)=ztv(ig,l)
          enddo
       enddo
       do ig=1,ngrid
          wmax(ig)=0.
          wmaxa(ig)=0.
       enddo
!calcul de la vitesse a partir de la CAPE en melangeant thetav

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A eliminer
! Ce if complique etait fait pour reperer la premiere couche instable
! Ici, c'est lmin.
!
!       do l=1,nlay-2
!         do ig=1,ngrid
!            if (ztv(ig,l).gt.ztv(ig,l+1)  &
!     &         .and.alim_star(ig,l).gt.1.e-10  &
!     &         .and.zw2(ig,l).lt.1e-10) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Calcul des F^*, integrale verticale de E^*
       f_star(:,1)=0.
       do l=1,nlay
          f_star(:,l+1)=f_star(:,l)+alim_star(:,l)
       enddo

! niveau (reel) auquel zw2 s'annule FH :n'etait pas initialise
       linter(:)=0.

! couche la plus haute concernee par le thermique. 
       lmax(:)=1

! Le niveau linter est une variable continue qui se trouve dans la couche
! lmax

       do l=1,nlay-2
         do ig=1,ngrid
            if (l.eq.lmin(ig).and.lalim(ig).gt.1) then

!------------------------------------------------------------------------
!  Calcul de la vitesse en haut de la premiere couche instable.
!  Premiere couche du panache thermique
!------------------------------------------------------------------------
               zw2(ig,l+1)=2.*RG*(ztv(ig,l)-ztv(ig,l+1))/ztv(ig,l+1)  &
     &                     *(zlev(ig,l+1)-zlev(ig,l))  &
     &                     *0.4*pphi(ig,l)/(pphi(ig,l+1)-pphi(ig,l))

!------------------------------------------------------------------------
! Tant que la vitesse en bas de la couche et la somme du flux de masse
! et de l'entrainement (c'est a dire le flux de masse en haut) sont
! positifs, on calcul
! 1. le flux de masse en haut  f_star(ig,l+1)
! 2. la temperature potentielle virtuelle dans la couche ztva(ig,l)
! 3. la vitesse au carré en haut zw2(ig,l+1)
!------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  A eliminer : dans cette version, si zw2 est > 0 on a un therique.
!  et donc, au dessus, f_star(ig,l+1) est forcement suffisamment 
!  grand puisque on n'a pas de detrainement.
!  f_star est une fonction croissante.
!  c'est donc vraiment sur zw2 uniquement qu'il faut faire le test.
!           else if ((zw2(ig,l).ge.1e-10).and.  &
!    &               (f_star(ig,l)+alim_star(ig,l).gt.1.e-10)) then
!              f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            else if (zw2(ig,l).ge.1e-10) then

               ztva(ig,l)=(f_star(ig,l)*ztva(ig,l-1)+alim_star(ig,l)  &
     &                    *ztv(ig,l))/f_star(ig,l+1)
               zw2(ig,l+1)=zw2(ig,l)*(f_star(ig,l)/f_star(ig,l+1))**2+  &
     &                     2.*RG*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)  &
     &                     *(zlev(ig,l+1)-zlev(ig,l))
            endif
! determination de zmax continu par interpolation lineaire
!------------------------------------------------------------------------

            if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
!               print*,'On tombe sur le cas particulier de thermcell_dry'
                zw2(ig,l+1)=0.
                linter(ig)=l+1
                lmax(ig)=l
            endif

            if (zw2(ig,l+1).lt.0.) then
               linter(ig)=(l*(zw2(ig,l+1)-zw2(ig,l))  &
     &           -zw2(ig,l))/(zw2(ig,l+1)-zw2(ig,l))
               zw2(ig,l+1)=0.
               lmax(ig)=l
            endif

               wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))

            if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
!   lmix est le niveau de la couche ou w (wa_moy) est maximum
               lmix(ig)=l+1
               wmaxa(ig)=wa_moy(ig,l+1)
            endif
         enddo
      enddo
       if (prt_level.ge.1) print*,'fin calcul zw2'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A eliminer :
! Ce calcul de lmax est fait en meme temps que celui de linter, plus haut
! Calcul de la couche correspondant a la hauteur du thermique
!      do ig=1,ngrid
!         lmax(ig)=lalim(ig)
!      enddo
!      do ig=1,ngrid
!         do l=nlay,lalim(ig)+1,-1
!            if (zw2(ig,l).le.1.e-10) then
!               lmax(ig)=l-1
!            endif
!         enddo
!      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    
! Determination de zw2 max
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

!   Longueur caracteristique correspondant a la hauteur des thermiques.
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo
      do  ig=1,ngrid
! calcul de zlevinter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH A eliminer
! Simplification
!          zlevinter(ig)=(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))*  &
!     &    linter(ig)+zlev(ig,lmax(ig))-lmax(ig)*(zlev(ig,lmax(ig)+1)  &
!     &    -zlev(ig,lmax(ig)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          zlevinter(ig)=zlev(ig,lmax(ig)) + &
     &    (linter(ig)-lmax(ig))*(zlev(ig,lmax(ig)+1)-zlev(ig,lmax(ig)))
           zmax(ig)=max(zmax(ig),zlevinter(ig)-zlev(ig,lmin(ig)))
      enddo

! Verification que lalim<=lmax
      do ig=1,ngrid
         if(lalim(ig)>lmax(ig)) then
           if ( prt_level > 1 ) THEN
            print*,'WARNING thermcell_dry ig=',ig,'  lalim=',lalim(ig),'  lmax(ig)=',lmax(ig)
           endif
           lmax(ig)=lalim(ig)
         endif
      enddo
      
      RETURN
      END
!==============================================================================
      SUBROUTINE thermcellV0_init(ngrid,nlay,ztv,zlay,zlev,  &
     &                  lalim,lmin,alim_star,alim_star_tot,lev_out)

!----------------------------------------------------------------------
!thermcell_init: calcul du profil d alimentation du thermique
!----------------------------------------------------------------------
      IMPLICIT NONE
#include "iniprint.h"
#include "thermcell.h"

      INTEGER l,ig
!arguments d entree
      INTEGER ngrid,nlay
      REAL ztv(ngrid,nlay)
      REAL zlay(ngrid,nlay)
      REAL zlev(ngrid,nlay+1)
!arguments de sortie
      INTEGER lalim(ngrid)
      INTEGER lmin(ngrid)
      REAL alim_star(ngrid,nlay)
      REAL alim_star_tot(ngrid)
      integer lev_out                           ! niveau pour les print
      
      REAL zzalim(ngrid)
!CR: ponderation entrainement des couches instables
!def des alim_star tels que alim=f*alim_star      

      do l=1,nlay
         do ig=1,ngrid 
            alim_star(ig,l)=0.
         enddo
      enddo
! determination de la longueur de la couche d entrainement
      do ig=1,ngrid
         lalim(ig)=1
      enddo

      if (iflag_thermals_ed.ge.1) then
!si la premiÃ¨re couche est instable, on declenche un thermique
         do ig=1,ngrid
            if (ztv(ig,1).gt.ztv(ig,2)) then
               lmin(ig)=1
               lalim(ig)=2
               alim_star(ig,1)=1.
               alim_star_tot(ig)=alim_star(ig,1)
               if(prt_level.GE.10) print*,'init',alim_star(ig,1),alim_star_tot(ig)
            else
                lmin(ig)=1
                lalim(ig)=1
                alim_star(ig,1)=0.
                alim_star_tot(ig)=0. 
            endif
         enddo
 
         else
!else iflag_thermals_ed=0 ancienne def de l alim 

!on ne considere que les premieres couches instables
      do l=nlay-2,1,-1
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1).and.  &
     &          ztv(ig,l+1).le.ztv(ig,l+2)) then
               lalim(ig)=l+1
            endif
          enddo
      enddo

! determination du lmin: couche d ou provient le thermique

      do ig=1,ngrid
! FH initialisation de lmin a nlay plutot que 1.
!        lmin(ig)=nlay
         lmin(ig)=1
      enddo
      do l=nlay,2,-1
         do ig=1,ngrid
            if (ztv(ig,l-1).gt.ztv(ig,l)) then
               lmin(ig)=l-1
            endif
         enddo
      enddo
!
      zzalim(:)=0.
      do l=1,nlay-1
         do ig=1,ngrid 
             if (l<lalim(ig)) then
                zzalim(ig)=zzalim(ig)+zlay(ig,l)*(ztv(ig,l)-ztv(ig,l+1))
             endif
          enddo
      enddo
      do ig=1,ngrid
          if (lalim(ig)>1) then
             zzalim(ig)=zlay(ig,1)+zzalim(ig)/(ztv(ig,1)-ztv(ig,lalim(ig)))
          else
             zzalim(ig)=zlay(ig,1)
          endif
      enddo

      if(prt_level.GE.10) print*,'ZZALIM LALIM ',zzalim,lalim,zlay(1,lalim(1))

! definition de l'entrainement des couches
      if (1.eq.1) then
      do l=1,nlay-1
         do ig=1,ngrid 
            if (ztv(ig,l).gt.ztv(ig,l+1).and.  &
     &          l.ge.lmin(ig).and.l.lt.lalim(ig)) then
!def possibles pour alim_star: zdthetadz, dthetadz, zdtheta
             alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1)) 
            endif
         enddo
      enddo
      else
      do l=1,nlay-1
         do ig=1,ngrid
            if (ztv(ig,l).gt.ztv(ig,l+1).and.  &
     &          l.ge.lmin(ig).and.l.lt.lalim(ig)) then
             alim_star(ig,l)=max(3.*zzalim(ig)-zlay(ig,l),0.) &
     &        *(zlev(ig,l+1)-zlev(ig,l))
            endif
         enddo
      enddo
      endif
      
! pas de thermique si couche 1 stable
      do ig=1,ngrid
!CRnouveau test
        if (alim_star(ig,1).lt.1.e-10) then 
            do l=1,nlay
                alim_star(ig,l)=0.
            enddo
            lmin(ig)=1
         endif
      enddo 
! calcul de l alimentation totale
      do ig=1,ngrid
         alim_star_tot(ig)=0.
      enddo
      do l=1,nlay
         do ig=1,ngrid
            alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
         enddo
      enddo
!
! Calcul entrainement normalise
      do l=1,nlay 
         do ig=1,ngrid 
            if (alim_star_tot(ig).gt.1.e-10) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo
       
!on remet alim_star_tot a 1
      do ig=1,ngrid 
         alim_star_tot(ig)=1.
      enddo

      endif
!endif iflag_thermals_ed
      return 
      end  
!==============================================================================
