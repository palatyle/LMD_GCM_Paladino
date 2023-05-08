










!
! $Id: calfis.F 1407 2010-07-07 10:31:52Z fairhead $
!
C
C
      SUBROUTINE calfis(lafin,
     $                  jD_cur, jH_cur,
     $                  pucov,
     $                  pvcov,
     $                  pteta,
     $                  pq,
     $                  pmasse,
     $                  pps,
     $                  pp,
     $                  ppk,
     $                  pphis,
     $                  pphi,
     $                  pducov,
     $                  pdvcov,
     $                  pdteta,
     $                  pdq,
     $                  flxw,
     $                  pdufi,
     $                  pdvfi,
     $                  pdhfi,
     $                  pdqfi,
     $                  pdpsfi)
c
c    Auteur :  P. Le Van, F. Hourdin 
c   .........
      USE infotrac, ONLY: nqtot, niadv, tname
      USE control_mod, ONLY: planet_type, nsplit_phys
      USE write_field
      USE cpdet_mod, only: t2tpot,tpot2t
      USE callphysiq_mod, ONLY: call_physiq

! used only for zonal averages
      USE moyzon_mod
      USE comvert_mod, ONLY: presnivs,preff
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,kappa,cpp,g,rad,pi
      USE logic_mod, ONLY: moyzon_ch,moyzon_mu

      IMPLICIT NONE
c=======================================================================
c
c   1. rearrangement des tableaux et transformation
c      variables dynamiques  >  variables physiques
c   2. calcul des termes physiques
c   3. retransformation des tendances physiques en tendances dynamiques
c
c   remarques:
c   ----------
c
c    - les vents sont donnes dans la physique par leurs composantes 
c      naturelles.
c    - la variable thermodynamique de la physique est une variable
c      intensive :   T 
c      pour la dynamique on prend    T * ( preff / p(l) ) **kappa
c    - les deux seules variables dependant de la geometrie necessaires
c      pour la physique sont la latitude pour le rayonnement et 
c      l'aire de la maille quand on veut integrer une grandeur 
c      horizontalement.
c    - les points de la physique sont les points scalaires de la 
c      la dynamique; numerotation:
c          1 pour le pole nord
c          (jjm-1)*iim pour l'interieur du domaine
c          ngridmx pour le pole sud
c      ---> ngridmx=2+(jjm-1)*iim
c
c     Input :
c     -------
c       pucov           covariant zonal velocity
c       pvcov           covariant meridional velocity 
c       pteta           potential temperature
c       pps             surface pressure
c       pmasse          masse d'air dans chaque maille
c       pts             surface temperature  (K)
c       callrad         clef d'appel au rayonnement
c
c    Output :
c    --------
c        pdufi          tendency for the natural zonal velocity (ms-1)
c        pdvfi          tendency for the natural meridional velocity 
c        pdhfi          tendency for the potential temperature (K/s)
c        pdtsfi         tendency for the surface temperature
c
c        pdtrad         radiative tendencies  \  both input
c        pfluxrad       radiative fluxes      /  and output
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------

      include "dimensions.h"
      include "paramet.h"

      INTEGER ngridmx
      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )

      include "comgeom2.h"
      include "iniprint.h"

c    Arguments :
c    -----------
      LOGICAL,INTENT(IN) ::  lafin ! .true. for the very last call to physics
      REAL,INTENT(IN) :: jD_cur, jH_cur
      REAL,INTENT(IN) :: pvcov(iip1,jjm,llm) ! covariant meridional velocity
      REAL,INTENT(IN) :: pucov(iip1,jjp1,llm) ! covariant zonal velocity
      REAL,INTENT(IN) :: pteta(iip1,jjp1,llm) ! potential temperature
      REAL,INTENT(IN) :: pmasse(iip1,jjp1,llm) ! mass in each cell ! not used
      REAL,INTENT(IN) :: pq(iip1,jjp1,llm,nqtot) ! tracers
      REAL,INTENT(IN) :: pphis(iip1,jjp1) ! surface geopotential
      REAL,INTENT(IN) :: pphi(iip1,jjp1,llm) ! geopotential

      REAL,INTENT(IN) :: pdvcov(iip1,jjm,llm) ! dynamical tendency on vcov
      REAL,INTENT(IN) :: pducov(iip1,jjp1,llm) ! dynamical tendency on ucov
      REAL,INTENT(IN) :: pdteta(iip1,jjp1,llm) ! dynamical tendency on teta
! commentaire SL: pdq ne sert que pour le calcul de pcvgq,
! qui lui meme ne sert a rien dans la routine telle qu'elle est
! ecrite, et que j'ai donc commente....
      REAL,INTENT(IN) :: pdq(iip1,jjp1,llm,nqtot) ! dynamical tendency on tracers
      ! NB: pdq is only used to compute pcvgq which is in fact not used...

      REAL,INTENT(IN) :: pps(iip1,jjp1) ! surface pressure (Pa)
      REAL,INTENT(IN) :: pp(iip1,jjp1,llmp1) ! pressure at mesh interfaces (Pa)
      REAL,INTENT(IN) :: ppk(iip1,jjp1,llm) ! Exner at mid-layer
      REAL,INTENT(IN) :: flxw(iip1,jjp1,llm)  ! Vertical mass flux on lower mesh interfaces (kg/s) (on llm because flxw(:,:,llm+1)=0)

      ! tendencies (in */s) from the physics
      REAL,INTENT(OUT) :: pdvfi(iip1,jjm,llm) ! tendency on covariant meridional wind
      REAL,INTENT(OUT) :: pdufi(iip1,jjp1,llm) ! tendency on covariant zonal wind
      REAL,INTENT(OUT) :: pdhfi(iip1,jjp1,llm) ! tendency on potential temperature (K/s)
      REAL,INTENT(OUT) :: pdqfi(iip1,jjp1,llm,nqtot) ! tendency on tracers
      REAL,INTENT(OUT) :: pdpsfi(iip1,jjp1) ! tendency on surface pressure (Pa/s)

c    Local variables :
c    -----------------

      INTEGER i,j,l,ig0,ig,iq,iiq
      REAL zpsrf(ngridmx)
      REAL zplev(ngridmx,llm+1),zplay(ngridmx,llm)
      REAL zphi(ngridmx,llm),zphis(ngridmx)

      REAL zrot(iip1,jjm,llm) ! AdlC May 2014
      REAL zufi(ngridmx,llm), zvfi(ngridmx,llm)
      REAL zrfi(ngridmx,llm) ! relative wind vorticity
      REAL ztfi(ngridmx,llm),zqfi(ngridmx,llm,nqtot)
! ADAPTATION GCM POUR CP(T)
      REAL zteta(ngridmx,llm)
      REAL zpk(ngridmx,llm)

! RQ SL 13/10/10:
! Ces calculs ne servent pas. 
! Si necessaire, decommenter ces variables et les calculs...
!      REAL pcvgu(ngridmx,llm), pcvgv(ngridmx,llm)
!      REAL pcvgt(ngridmx,llm), pcvgq(ngridmx,llm,2)

      REAL zdufi(ngridmx,llm),zdvfi(ngridmx,llm)
      REAL zdtfi(ngridmx,llm),zdqfi(ngridmx,llm,nqtot)
      REAL zdpsrf(ngridmx)

      REAL zdufic(ngridmx,llm),zdvfic(ngridmx,llm)
      REAL zdtfic(ngridmx,llm),zdqfic(ngridmx,llm,nqtot)
      REAL jH_cur_split,zdt_split
      LOGICAL debut_split,lafin_split
      INTEGER isplit

      REAL zsin(iim),zcos(iim),z1(iim)
      REAL zsinbis(iim),zcosbis(iim),z1bis(iim)
      REAL unskap, pksurcp
      save unskap

      REAL flxwfi(ngridmx,llm)  ! Flux de masse verticale sur la grille physiq
      
      REAL SSUM

      LOGICAL,SAVE :: firstcal=.true., debut=.true.
!      REAL rdayvrai

! For Titan only right now: 
! to allow for 2D computation of microphys and chemistry
      LOGICAL,save :: flag_moyzon
      REAL,allocatable,save :: tmpvar(:,:)
      REAL,allocatable,save :: tmpvarp1(:,:)
      REAL,allocatable,save :: tmpvarbar(:)
      REAL,allocatable,save :: tmpvarbarp1(:)
      real :: zz1,zz2

c-----------------------------------------------------------------------

c    1. Initialisations :
c    --------------------


      IF ( firstcal )  THEN
        debut = .TRUE.
        IF (ngridmx.NE.2+(jjm-1)*iim) THEN
         write(lunout,*) 'STOP dans calfis'
         write(lunout,*)
     &   'La dimension ngridmx doit etre egale a 2 + (jjm-1)*iim'
         write(lunout,*) '  ngridmx  jjm   iim   '
         write(lunout,*) ngridmx,jjm,iim
         STOP
        ENDIF

        unskap   = 1./ kappa

        flag_moyzon = .false.
        if(moyzon_ch.or.moyzon_mu) then
         flag_moyzon = .true.
         allocate(tmpvar(iip1,llm))
         allocate(tmpvarp1(iip1,llmp1))
         allocate(tmpvarbar(llm))
         allocate(tmpvarbarp1(llmp1))
        endif

        if (flag_moyzon) call moyzon_init(ngridmx,llm,nqtot)

c----------------------------------------------
c moyennes globales pour le profil de pression
       if(planet_type.eq."titan".or.planet_type.eq."venus") then
        ALLOCATE(plevmoy(llm+1))
        ALLOCATE(playmoy(llm))
        ALLOCATE(tmoy(llm))
        ALLOCATE(tetamoy(llm))
        ALLOCATE(pkmoy(llm))
        ALLOCATE(phimoy(0:llm))
        ALLOCATE(zlevmoy(llm+1))
        ALLOCATE(zlaymoy(llm))
        plevmoy=0.
        do l=1,llmp1
         do i=1,iip1
          do j=1,jjp1
            plevmoy(l)=plevmoy(l)+pp(i,j,l)/(iip1*jjp1)
          enddo
         enddo
        enddo
        tetamoy=0.
        pkmoy=0.
        phimoy=0.
        do i=1,iip1
         do j=1,jjp1
            phimoy(0)=phimoy(0)+pphis(i,j)/(iip1*jjp1)
         enddo
        enddo
        do l=1,llm
         do i=1,iip1
          do j=1,jjp1
            tetamoy(l)=tetamoy(l)+pteta(i,j,l)/(iip1*jjp1)
            pkmoy(l)=pkmoy(l)+ppk(i,j,l)/(iip1*jjp1)
            phimoy(l)=phimoy(l)+pphi(i,j,l)/(iip1*jjp1)
          enddo
         enddo
        enddo
        playmoy(:) = preff * (pkmoy(:)/cpp) ** unskap
        call tpot2t(llm,tetamoy,tmoy,pkmoy)
c SI ON TIENT COMPTE DE LA VARIATION DE G AVEC L'ALTITUDE:
        zlaymoy(1:llm) = g*rad*rad/(g*rad-phimoy(1:llm))-rad
        zlevmoy(1) = phimoy(0)/g
        DO l=2,llm
            zz1=(playmoy(l-1)+plevmoy(l))/(playmoy(l-1)-plevmoy(l))
            zz2=(plevmoy(l)  +playmoy(l))/(plevmoy(l)  -playmoy(l))
            zlevmoy(l)=(zz1*zlaymoy(l-1)+zz2*zlaymoy(l))/(zz1+zz2)
        ENDDO
        zlevmoy(llmp1)=zlaymoy(llm)+(zlaymoy(llm)-zlevmoy(llm))
c-------------------
c + lat index
        allocate(klat(ngridmx))
        klat=0
        klat(1)  = 1
        ig0  = 2
        DO j = 2,jjm
         do i=0,iim-1
          klat(ig0+i) = j
         enddo
         ig0 = ig0+iim
        ENDDO
        klat(ngridmx)  = jjp1
       endif   ! planet_type=titan
c----------------------------------------------
      ELSE
        debut = .FALSE.
      ENDIF ! of IF (firstcal)


c-----------------------------------------------------------------------
c   40. transformation des variables dynamiques en variables physiques:
c   ---------------------------------------------------------------

c   41. pressions au sol (en Pascals)
c   ----------------------------------
       
      zpsrf(1) = pps(1,1)

      ig0  = 2
      DO j = 2,jjm
         CALL SCOPY( iim,pps(1,j),1,zpsrf(ig0), 1 )
         ig0 = ig0+iim
      ENDDO

      zpsrf(ngridmx) = pps(1,jjp1)

c   42. pression intercouches et fonction d'Exner:

c   -----------------------------------------------------------------
c     .... zplev  definis aux (llm +1) interfaces des couches  ....
c     .... zplay  definis aux (  llm )    milieux des couches  .... 
c   -----------------------------------------------------------------

c    ...    Exner = cp * ( p(l) / preff ) ** kappa     ....

! ADAPTATION GCM POUR CP(T)
      DO l = 1, llm
        zpk(   1,l ) = ppk(1,1,l)
        zteta( 1,l ) = pteta(1,1,l)
        zplev( 1,l ) = pp(1,1,l)
        ig0 = 2
          DO j = 2, jjm
             DO i =1, iim
              zpk(   ig0,l ) = ppk(i,j,l)
              zteta( ig0,l ) = pteta(i,j,l)
              zplev( ig0,l ) = pp(i,j,l)
              ig0 = ig0 +1
             ENDDO
          ENDDO
        zpk(   ngridmx,l ) = ppk(1,jjp1,l)
        zteta( ngridmx,l ) = pteta(1,jjp1,l)
        zplev( ngridmx,l ) = pp(1,jjp1,l)
      ENDDO
        zplev( 1,llmp1 ) = pp(1,1,llmp1)
        ig0 = 2
          DO j = 2, jjm
             DO i =1, iim
              zplev( ig0,llmp1 ) = pp(i,j,llmp1)
              ig0 = ig0 +1
             ENDDO
          ENDDO
        zplev( ngridmx,llmp1 ) = pp(1,jjp1,llmp1)

      if (flag_moyzon) then
        tmpvarp1(:,:) = pp(:,1,:)
        call moyzon(llmp1,tmpvarp1,tmpvarbarp1)
        zplevbar(1,:) = tmpvarbarp1
        tmpvar(:,:) = ppk(:,1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zpkbar(1,:) = tmpvarbar
        tmpvar(:,:) = pteta(:,1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        ztetabar(1,:) = tmpvarbar
        call tpot2t(llm,ztetabar(1,:),ztfibar(1,:),zpkbar(1,:))
        ig0 = 2
         do j = 2, jjm
          tmpvarp1(:,:) = pp(:,j,:)
          call moyzon(llmp1,tmpvarp1,tmpvarbarp1)
          zplevbar(ig0,:) = tmpvarbarp1
          tmpvar(:,:) = ppk(:,j,:)
          call moyzon(llm,tmpvar,tmpvarbar)
          zpkbar(ig0,:) = tmpvarbar
          tmpvar(:,:) = pteta(:,j,:)
          call moyzon(llm,tmpvar,tmpvarbar)
          ztetabar(ig0,:) = tmpvarbar
          call tpot2t(llm,ztetabar(ig0,:),ztfibar(ig0,:),zpkbar(ig0,:))
          ig0 = ig0+1
          do i=2,iim
            zplevbar(ig0,:) = zplevbar(ig0-1,:)
            zpkbar(ig0,:)   = zpkbar(ig0-1,:)
            ztetabar(ig0,:) = ztetabar(ig0-1,:)
            ztfibar(ig0,:)  = ztfibar(ig0-1,:)
            ig0 = ig0+1
          enddo
         enddo
        tmpvarp1(:,:) = pp(:,jjp1,:)
        call moyzon(llmp1,tmpvarp1,tmpvarbarp1)
        zplevbar(ngridmx,:) = tmpvarbarp1
        tmpvar(:,:) = ppk(:,jjp1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zpkbar(ngridmx,:) = tmpvarbar
        tmpvar(:,:) = pteta(:,jjp1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        ztetabar(ngridmx,:) = tmpvarbar
        call tpot2t(llm,ztetabar(ngridmx,:),
     .                  ztfibar(ngridmx,:),zpkbar(ngridmx,:))
      endif

c   43. temperature naturelle (en K) et pressions milieux couches .
c   ---------------------------------------------------------------

! ADAPTATION GCM POUR CP(T)
         call tpot2t(ngridmx*llm,zteta,ztfi,zpk)

      DO l=1,llm

         pksurcp     =  ppk(1,1,l) / cpp
         zplay(1,l)  =  preff * pksurcp ** unskap
!         pcvgt(1,l)  =  pdteta(1,1,l) * pksurcp / pmasse(1,1,l)
         ig0         = 2

         DO j = 2, jjm
            DO i = 1, iim
              pksurcp        = ppk(i,j,l) / cpp
              zplay(ig0,l)   = preff * pksurcp ** unskap
!              pcvgt(ig0,l)   = pdteta(i,j,l) * pksurcp / pmasse(i,j,l)
              ig0            = ig0 + 1
            ENDDO
         ENDDO

         pksurcp       = ppk(1,jjp1,l) / cpp
         zplay(ig0,l)  = preff * pksurcp ** unskap
!         pcvgt(ig0,l)  = pdteta(1,jjp1,l) * pksurcp/ pmasse(1,jjp1,l)

      ENDDO

      if (flag_moyzon) then
        zplaybar(:,:) = preff * (zpkbar(:,:)/cpp)**unskap
      endif

c   43.bis traceurs (tous intensifs) 
c   ---------------

      DO iq=1,nqtot
         DO l=1,llm
            zqfi(1,l,iq) = pq(1,1,l,iq)
            ig0          = 2
            DO j=2,jjm
               DO i = 1, iim
                  zqfi(ig0,l,iq)  = pq(i,j,l,iq)
                  ig0             = ig0 + 1
               ENDDO
            ENDDO
            zqfi(ig0,l,iq) = pq(1,jjp1,l,iq)
         ENDDO
      ENDDO  ! boucle sur traceurs

      if (flag_moyzon) then
       DO iq=1,nqtot
! RQ: REVOIR A QUOI CA SERT... ET VERIFIER...
!       iiq=niadv(iq) 
! en fait, iiq=iq...
! FIN RQ
        tmpvar(:,:) = pq(:,1,:,iq)
        call moyzon(llm,tmpvar,tmpvarbar)
        zqfibar(1,:,iq) = tmpvarbar
        ig0 = 2
         do j = 2, jjm
          tmpvar(:,:) = pq(:,j,:,iq)
          call moyzon(llm,tmpvar,tmpvarbar)
          zqfibar(ig0,:,iq) = tmpvarbar
          ig0 = ig0+1
          do i=2,iim
            zqfibar(ig0,:,iq) = zqfibar(ig0-1,:,iq)
            ig0 = ig0+1
          enddo
         enddo
        tmpvar(:,:) = pq(:,jjp1,:,iq)
        call moyzon(llm,tmpvar,tmpvarbar)
        zqfibar(ngridmx,:,iq) = tmpvarbar
       ENDDO ! of DO iq=1,nqtot
      endif

! DEBUG
!     do ig0=1,ngridmx
!       write(*,'(6(e13.5,1x))') zqfibar(ig0,1,10),zqfi(ig0,1,10),
!    .                         zqfibar(ig0,llm/2,10),zqfi(ig0,llm/2,10), 
!    .                           zqfibar(ig0,llm,10),zqfi(ig0,llm,10) 
!     enddo
!     stop

!-----------------
! RQ SL 13/10/10:
! Ces calculs ne servent pas. 
! Si necessaire, decommenter ces variables et les calculs...
!
!   convergence dynamique pour les traceurs "EAU"
! Earth-specific treatment of first 2 tracers (water)
!      if (planet_type=="earth") then
!       DO iq=1,2
!        DO l=1,llm
!           pcvgq(1,l,iq)= pdq(1,1,l,iq) / pmasse(1,1,l)
!           ig0          = 2
!           DO j=2,jjm
!              DO i = 1, iim
!                 pcvgq(ig0,l,iq) = pdq(i,j,l,iq) / pmasse(i,j,l)
!                 ig0             = ig0 + 1
!              ENDDO
!           ENDDO
!           pcvgq(ig0,l,iq)= pdq(1,jjp1,l,iq) / pmasse(1,jjp1,l)
!        ENDDO
!       ENDDO
!      endif ! of if (planet_type=="earth")
!----------------

c   Geopotentiel calcule par rapport a la surface locale:
c   -----------------------------------------------------

      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,pphi,zphi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,pphis,zphis)
      DO l=1,llm
         DO ig=1,ngridmx
           zphi(ig,l)=zphi(ig,l)-zphis(ig)
         ENDDO
      ENDDO

      if (flag_moyzon) then
        tmpvar(:,1) = pphis(:,1)
        call moyzon(1,tmpvar(:,1),tmpvarbar(1))
        zphisbar(1) = tmpvarbar(1)
        tmpvar(:,:) = pphi(:,1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zphibar(1,:) = tmpvarbar
        ig0 = 2
         do j = 2, jjm
          tmpvar(:,1) = pphis(:,j)
          call moyzon(1,tmpvar(:,1),tmpvarbar(1))
          zphisbar(ig0) = tmpvarbar(1)
          tmpvar(:,:) = pphi(:,j,:)
          call moyzon(llm,tmpvar,tmpvarbar)
          zphibar(ig0,:) = tmpvarbar
          ig0 = ig0+1
          do i=2,iim
            zphisbar(ig0)  = zphisbar(ig0-1)
            zphibar(ig0,:) = zphibar(ig0-1,:)
            ig0 = ig0+1
          enddo
         enddo
        tmpvar(:,1) = pphis(:,jjp1)
        call moyzon(1,tmpvar(:,1),tmpvarbar(1))
        zphisbar(ngridmx) = tmpvarbar(1)
        tmpvar(:,:) = pphi(:,jjp1,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zphibar(ngridmx,:) = tmpvarbar
      endif

c   ....  Calcul de la vitesse  verticale  ( en Pa*m*s  ou Kg/s )  ....
c JG : ancien calcule de omega utilise dans physiq.F. Maintenant le flux 
c    de masse est calclue dans advtrac.F  
c      DO l=1,llm
c        pvervel(1,l)=pw(1,1,l) * g /apoln
c        ig0=2
c       DO j=2,jjm
c           DO i = 1, iim
c              pvervel(ig0,l) = pw(i,j,l) * g * unsaire(i,j)
c              ig0 = ig0 + 1
c           ENDDO
c       ENDDO
c        pvervel(ig0,l)=pw(1,jjp1,l) * g /apols
c      ENDDO

c
c   45. champ u:
c   ------------

      DO 50 l=1,llm

         DO 25 j=2,jjm
            ig0 = 1+(j-2)*iim
            zufi(ig0+1,l)= 0.5 * 
     $      ( pucov(iim,j,l)/cu(iim,j) + pucov(1,j,l)/cu(1,j) )
!            pcvgu(ig0+1,l)= 0.5 * 
!     $      ( pducov(iim,j,l)/cu(iim,j) + pducov(1,j,l)/cu(1,j) )
            DO 10 i=2,iim
               zufi(ig0+i,l)= 0.5 *
     $         ( pucov(i-1,j,l)/cu(i-1,j) + pucov(i,j,l)/cu(i,j) )
!               pcvgu(ig0+i,l)= 0.5 *
!     $         ( pducov(i-1,j,l)/cu(i-1,j) + pducov(i,j,l)/cu(i,j) )
10         CONTINUE
25      CONTINUE

50    CONTINUE


C  Alvaro de la Camara (May 2014)
C  46.1 Calcul de la vorticite et passage sur la grille physique
C  --------------------------------------------------------------
      DO l=1,llm
        do i=1,iim
          do j=1,jjm
            zrot(i,j,l) = (pvcov(i+1,j,l) - pvcov(i,j,l)
     $                   + pucov(i,j+1,l) - pucov(i,j,l)) 
     $                   / (cu(i,j)+cu(i,j+1)) 
     $                   / (cv(i+1,j)+cv(i,j)) *4
          enddo
        enddo
      ENDDO

c   46.champ v:
c   -----------

      DO l=1,llm
         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
               zvfi(ig0+i,l)= 0.5 *
     $         ( pvcov(i,j-1,l)/cv(i,j-1) + pvcov(i,j,l)/cv(i,j) )
c               pcvgv(ig0+i,l)= 0.5 *
c     $         ( pdvcov(i,j-1,l)/cv(i,j-1) + pdvcov(i,j,l)/cv(i,j) )
            ENDDO
               zrfi(ig0 + 1,l)= 0.25 *(zrot(iim,j-1,l)+zrot(iim,j,l)
     &                                +zrot(1,j-1,l)+zrot(1,j,l))
            DO i=2,iim
               zrfi(ig0 + i,l)= 0.25 *(zrot(i-1,j-1,l)+zrot(i-1,j,l)
     $                   +zrot(i,j-1,l)+zrot(i,j,l))   !  AdlC MAY 2014
            ENDDO
         ENDDO
      ENDDO


c   47. champs de vents aux pole nord   
c   ------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,1,l)/cv(1,1)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,1,l)/cv(1,1)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,1,l)/cv(i,1)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,1,l)/cv(i,1)
         ENDDO

         DO i=1,iim
            zcos(i)   = COS(rlonv(i))*z1(i)
            zcosbis(i)= COS(rlonv(i))*z1bis(i)
            zsin(i)   = SIN(rlonv(i))*z1(i)
            zsinbis(i)= SIN(rlonv(i))*z1bis(i)
         ENDDO

         zufi(1,l)  = SSUM(iim,zcos,1)/pi
!         pcvgu(1,l) = SSUM(iim,zcosbis,1)/pi
         zvfi(1,l)  = SSUM(iim,zsin,1)/pi
!         pcvgv(1,l) = SSUM(iim,zsinbis,1)/pi
         zrfi(1, l) = 0.
      ENDDO


c   48. champs de vents aux pole sud:
c   ---------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,jjm,l)/cv(1,jjm)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,jjm,l)/cv(1,jjm)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,jjm,l)/cv(i,jjm)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,jjm,l)/cv(i,jjm)
         ENDDO

         DO i=1,iim
            zcos(i)    = COS(rlonv(i))*z1(i)
            zcosbis(i) = COS(rlonv(i))*z1bis(i)
            zsin(i)    = SIN(rlonv(i))*z1(i)
            zsinbis(i) = SIN(rlonv(i))*z1bis(i)
         ENDDO

         zufi(ngridmx,l)  = SSUM(iim,zcos,1)/pi
!         pcvgu(ngridmx,l) = SSUM(iim,zcosbis,1)/pi
         zvfi(ngridmx,l)  = SSUM(iim,zsin,1)/pi
!         pcvgv(ngridmx,l) = SSUM(iim,zsinbis,1)/pi
         zrfi(ngridmx, l) = 0.
      ENDDO
c
c On change de grille, dynamique vers physiq, pour le flux de masse verticale
      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,flxw,flxwfi)

c-----------------------------------------------------------------------
c   Appel de la physique:
c   ---------------------

! Appel de la physique: pose probleme quand on tourne 
! SANS physique, car physiq.F est dans le repertoire phy[]... 
! Il faut une cle 1

! Le fait que les arguments de physiq soient differents selon les planetes 
! ne pose pas de probleme a priori.

!      write(lunout,*) 'PHYSIQUE AVEC NSPLIT_PHYS=',nsplit_phys
      zdt_split=dtphys/nsplit_phys
      zdufic(:,:)=0.
      zdvfic(:,:)=0.
      zdtfic(:,:)=0.
      zdqfic(:,:,:)=0.


      do isplit=1,nsplit_phys

         jH_cur_split=jH_cur+(isplit-1) * dtvr / (daysec *nsplit_phys)
         debut_split=debut.and.isplit==1
         lafin_split=lafin.and.isplit==nsplit_phys

        CALL call_physiq(ngridmx,llm,nqtot,tname,
     &                   debut_split,lafin_split,
     &                   jD_cur,jH_cur_split,zdt_split,
     &                   zplev,zplay,
     &                   zpk,zphi,zphis,
     &                   presnivs,
     &                   zufi,zvfi,zrfi,ztfi,zqfi,
     &                   flxwfi,pducov,
     &                   zdufi,zdvfi,zdtfi,zdqfi,zdpsrf)

!      if (planet_type.eq."earth") then
!         CALL physiq (ngridmx,
!     .             llm,
!     .             debut_split,
!     .             lafin_split,
!     .             jD_cur,
!     .             jH_cur_split,
!     .             zdt_split,
!     .             zplev,
!     .             zplay,
!     .             zphi,
!     .             zphis,
!     .             presnivs,
!     .             zufi,
!     .             zvfi,
!     .             ztfi,
!     .             zqfi,
!     .             flxwfi,
!     .             zdufi,
!     .             zdvfi,
!     .             zdtfi,
!     .             zdqfi,
!     .             zdpsrf,
!     .             pducov)
!
!      else if ( planet_type=="generic" ) then
!
!         CALL physiq (ngridmx,     !! ngrid
!     .             llm,            !! nlayer
!     .             nqtot,          !! nq
!     .             tname,          !! tracer names from dynamical core (given in infotrac)
!     .             debut_split,    !! firstcall 
!     .             lafin_split,    !! lastcall
!     .             jD_cur,         !! pday. see leapfrog
!     .             jH_cur_split,   !! ptime "fraction of day"
!     .             zdt_split,      !! ptimestep
!     .             zplev,          !! pplev
!     .             zplay,          !! pplay
!     .             zphi,           !! pphi
!     .             zufi,           !! pu
!     .             zvfi,           !! pv
!     .             ztfi,           !! pt
!     .             zqfi,           !! pq
!     .             flxwfi,         !! pw !! or 0. anyway this is for diagnostic. not used in physiq.
!     .             zdufi,          !! pdu
!     .             zdvfi,          !! pdv
!     .             zdtfi,          !! pdt
!     .             zdqfi,          !! pdq
!     .             zdpsrf,         !! pdpsrf
!     .             tracerdyn)      !! tracerdyn <-- utilite ???
!
!      else if ( planet_type=="mars" ) then
!
!        CALL physiq (ngridmx,    ! ngrid
!     .             llm,          ! nlayer
!     .             nqtot,        ! nq
!     .             debut_split,  ! firstcall
!     .             lafin_split,  ! lastcall
!     .             jD_cur,       ! pday
!     .             jH_cur_split, ! ptime
!     .             zdt_split,    ! ptimestep
!     .             zplev,        ! pplev
!     .             zplay,        ! pplay
!     .             zphi,         ! pphi
!     .             zufi,         ! pu
!     .             zvfi,         ! pv
!     .             ztfi,         ! pt
!     .             zqfi,         ! pq
!     .             flxwfi,       ! pw
!     .             zdufi,        ! pdu
!     .             zdvfi,        ! pdv
!     .             zdtfi,        ! pdt
!     .             zdqfi,        ! pdq
!     .             zdpsrf,       ! pdpsrf
!     .             tracerdyn)    ! tracerdyn (somewhat obsolete)
!
!      else if ((planet_type=="titan").or.(planet_type=="venus")) then
!
!         CALL physiq (ngridmx,
!     .             llm,
!     .             nqtot,
!     .             debut_split,
!     .             lafin_split,
!     .             jD_cur,
!     .             jH_cur_split,
!     .             zdt_split,
!     .             zplev,
!     .             zplay,
!     .             zpk,
!     .             zphi,
!     .             zphis,
!     .             presnivs,
!     .             zufi,
!     .             zvfi,
!     .             ztfi,
!     .             zqfi,
!     .             flxwfi,
!     .             zdufi,
!     .             zdvfi,
!     .             zdtfi,
!     .             zdqfi,
!     .             zdpsrf)
!
!      else ! unknown "planet_type"
!
!        write(lunout,*) "calfis_p: error, unknown planet_type: ",
!     &                  trim(planet_type)
!        stop
!
!      endif ! planet_type

         zufi(:,:)=zufi(:,:)+zdufi(:,:)*zdt_split
         zvfi(:,:)=zvfi(:,:)+zdvfi(:,:)*zdt_split
         ztfi(:,:)=ztfi(:,:)+zdtfi(:,:)*zdt_split
         zqfi(:,:,:)=zqfi(:,:,:)+zdqfi(:,:,:)*zdt_split

         zdufic(:,:)=zdufic(:,:)+zdufi(:,:)
         zdvfic(:,:)=zdvfic(:,:)+zdvfi(:,:)
         zdtfic(:,:)=zdtfic(:,:)+zdtfi(:,:)
         zdqfic(:,:,:)=zdqfic(:,:,:)+zdqfi(:,:,:)

      enddo ! of do isplit=1,nsplit_phys

! ATTENTION...
      if (flag_moyzon.and.(nsplit_phys.ne.1)) then
       print*,"WARNING ! flag_moyzon + nsplit_phys"
       print*,"zqfibar n'est pas implemente au cours des iterations"
       print*,"Donc a revoir..."
       stop
      endif

! #endif of #ifdef 1

      zdufi(:,:)=zdufic(:,:)/nsplit_phys
      zdvfi(:,:)=zdvfic(:,:)/nsplit_phys
      zdtfi(:,:)=zdtfic(:,:)/nsplit_phys
      zdqfi(:,:,:)=zdqfic(:,:,:)/nsplit_phys


500   CONTINUE

c-----------------------------------------------------------------------
c   transformation des tendances physiques en tendances dynamiques:
c   ---------------------------------------------------------------

c  tendance sur la pression :
c  -----------------------------------

      CALL gr_fi_dyn(1,ngridmx,iip1,jjp1,zdpsrf,pdpsfi)
c
c   62. enthalpie potentielle
c   ---------------------

! ADAPTATION GCM POUR CP(T)
      call t2tpot(ngridmx*llm,ztfi,zteta,zpk)

         DO i=1,iip1
      pdhfi(i,1,:) = (zteta(1,:) - pteta(i,1,:))/dtphys
      pdhfi(i,jjp1,:) = (zteta(ngridmx,:) - pteta(i,jjp1,:))/dtphys
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
      pdhfi(i,j,:) = (zteta(ig0+i,:) - pteta(i,j,:))/dtphys
            ENDDO
               pdhfi(iip1,j,:) =  pdhfi(1,j,:)
         ENDDO


c   62. humidite specifique
c   ---------------------
! Ehouarn: removed this useless bit: was overwritten at step 63 anyways
!      DO iq=1,nqtot
!         DO l=1,llm
!            DO i=1,iip1
!               pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)
!               pdqfi(i,jjp1,l,iq) = zdqfi(ngridmx,l,iq)
!            ENDDO
!            DO j=2,jjm
!               ig0=1+(j-2)*iim
!               DO i=1,iim
!                  pdqfi(i,j,l,iq) = zdqfi(ig0+i,l,iq)
!               ENDDO
!               pdqfi(iip1,j,l,iq) = pdqfi(1,j,l,iq)
!            ENDDO
!         ENDDO
!      ENDDO

c   63. traceurs (tous en intensifs)
c   ------------
C     initialisation des tendances
      pdqfi(:,:,:,:)=0.
C
       DO iq=1,nqtot
         iiq=niadv(iq)
         DO l=1,llm
            DO i=1,iip1
               pdqfi(i,1,l,iiq)    = zdqfi(1,l,iq)
               pdqfi(i,jjp1,l,iiq) = zdqfi(ngridmx,l,iq)
            ENDDO
            DO j=2,jjm
               ig0=1+(j-2)*iim
               DO i=1,iim
                  pdqfi(i,j,l,iiq) = zdqfi(ig0+i,l,iq)
               ENDDO
               pdqfi(iip1,j,l,iiq) = pdqfi(1,j,l,iq)
            ENDDO
         ENDDO
       ENDDO

c   65. champ u:
c   ------------

      DO l=1,llm

         DO i=1,iip1
            pdufi(i,1,l)    = 0.
            pdufi(i,jjp1,l) = 0.
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim-1
               pdufi(i,j,l)=
     $         0.5*(zdufi(ig0+i,l)+zdufi(ig0+i+1,l))*cu(i,j)
            ENDDO
            pdufi(iim,j,l)=
     $      0.5*(zdufi(ig0+1,l)+zdufi(ig0+iim,l))*cu(iim,j)
            pdufi(iip1,j,l)=pdufi(1,j,l)
         ENDDO

      ENDDO


c   67. champ v:
c   ------------

      DO l=1,llm

         DO j=2,jjm-1
            ig0=1+(j-2)*iim
            DO i=1,iim
               pdvfi(i,j,l)=
     $         0.5*(zdvfi(ig0+i,l)+zdvfi(ig0+i+iim,l))*cv(i,j)
            ENDDO
            pdvfi(iip1,j,l) = pdvfi(1,j,l)
         ENDDO
      ENDDO


c   68. champ v pres des poles:
c   ---------------------------
c      v = U * cos(long) + V * SIN(long)

      DO l=1,llm

         DO i=1,iim
            pdvfi(i,1,l)=
     $      zdufi(1,l)*COS(rlonv(i))+zdvfi(1,l)*SIN(rlonv(i))
            pdvfi(i,jjm,l)=zdufi(ngridmx,l)*COS(rlonv(i))
     $      +zdvfi(ngridmx,l)*SIN(rlonv(i))
            pdvfi(i,1,l)=
     $      0.5*(pdvfi(i,1,l)+zdvfi(i+1,l))*cv(i,1)
            pdvfi(i,jjm,l)=
     $      0.5*(pdvfi(i,jjm,l)+zdvfi(ngridmx-iip1+i,l))*cv(i,jjm)
          ENDDO

         pdvfi(iip1,1,l)  = pdvfi(1,1,l)
         pdvfi(iip1,jjm,l)= pdvfi(1,jjm,l)

      ENDDO

c-----------------------------------------------------------------------

700   CONTINUE
 
      firstcal = .FALSE.

! of #ifndef CPP_PARA
      END
