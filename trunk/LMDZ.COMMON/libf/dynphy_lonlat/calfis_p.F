!
! $Id: calfis_p.F 1407 2010-07-07 10:31:52Z fairhead $
!
C
C
      SUBROUTINE calfis_p(lafin,
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
#ifdef CPP_PHYS
! Ehouarn: if using (parallelized) physics
      USE dimphy
      USE mod_phys_lmdz_mpi_data, mpi_root_xx=>mpi_master
      USE mod_phys_lmdz_omp_data, ONLY: klon_omp, klon_omp_begin
      USE mod_const_mpi, ONLY: COMM_LMDZ
      USE mod_interface_dyn_phys
!      USE IOPHY
#endif
#ifdef CPP_PARA
      USE parallel_lmdz, ONLY: omp_chunk, using_mpi, AllGather_Field,
     &                         jjb_u,jje_u,jjb_v,jje_v,
     &                         jj_begin_dyn=>jj_begin,jj_end_dyn=>jj_end
      USE Write_Field
      Use Write_field_p
      USE Times
      USE cpdet_mod, only: tpot2t_p, t2tpot_p
! used only for zonal averages
      USE moyzon_mod
#endif
      USE infotrac, ONLY: nqtot, niadv, tname
      USE control_mod, ONLY: planet_type, nsplit_phys
      USE comvert_mod, ONLY: preff,presnivs
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,kappa,cpp,g,rad,pi
      USE logic_mod, ONLY: moyzon_ch,moyzon_mu
#ifdef CPP_PHYS
      USE callphysiq_mod, ONLY: call_physiq
#endif 

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
c       ecritphy        frequence d'ecriture (en jours)de histphy
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
#ifdef CPP_MPI
      include 'mpif.h'
#endif
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

#ifdef CPP_PARA
#ifdef CPP_PHYS
c    Local variables :
c    -----------------

      INTEGER i,j,l,ig0,ig,iq,iiq
      REAL,ALLOCATABLE,SAVE :: zpsrf(:)
      REAL,ALLOCATABLE,SAVE :: zplev(:,:),zplay(:,:)
      REAL,ALLOCATABLE,SAVE :: zphi(:,:),zphis(:)

!      REAL zrot(iip1,jjb_v:jje_v,llm) ! AdlC May 2014
      REAL :: zrot(iip1,jjm,llm)
      REAL,ALLOCATABLE,SAVE :: zufi(:,:), zvfi(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfi(:,:),zqfi(:,:,:)
! ADAPTATION GCM POUR CP(T)
      REAL,ALLOCATABLE,SAVE :: zteta(:,:)
      REAL,ALLOCATABLE,SAVE ::  zpk(:,:)

! Ces calculs ne servent pas. 
! Si necessaire, decommenter ces variables et les calculs...
!      REAL,ALLOCATABLE,SAVE :: pcvgu(:,:), pcvgv(:,:)
!      REAL,ALLOCATABLE,SAVE :: pcvgt(:,:), pcvgq(:,:,:)
c
      REAL,ALLOCATABLE,SAVE :: zdufi(:,:),zdvfi(:,:), zrfi(:,:)
      REAL,ALLOCATABLE,SAVE :: zdtfi(:,:),zdqfi(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zdpsrf(:)
      REAL,SAVE,ALLOCATABLE ::  flxwfi(:,:)     ! Flux de masse verticale sur la grille physiq

      REAL,ALLOCATABLE,SAVE :: zplev_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zplay_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zpk_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zphi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zphis_omp(:)
      REAL,ALLOCATABLE,SAVE :: presnivs_omp(:)
      REAL,ALLOCATABLE,SAVE :: zufi_omp(:,:) 
      REAL,ALLOCATABLE,SAVE :: zvfi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zrfi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zqfi_omp(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zdufi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdvfi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdtfi_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdqfi_omp(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zdpsrf_omp(:)
      REAL,SAVE,ALLOCATABLE ::  flxwfi_omp(:,:)     ! Flux de masse verticale sur la grille physiq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Introduction du splitting (FH)
! Question pour Yann :
! J'ai été surpris au début que les tableaux zufi_omp, zdufi_omp n'co soitent
! en SAVE. Je crois comprendre que c'est parce que tu voulais qu'il
! soit allocatable (plutot par exemple que de passer une dimension
! dépendant du process en argument des routines) et que, du coup,
! le SAVE évite d'avoir à refaire l'allocation à chaque appel.
! Tu confirmes ?
! J'ai suivi le même principe pour les zdufic_omp
! Mais c'est surement bien que tu controles.
! 

      REAL,ALLOCATABLE,SAVE :: zdufic_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdvfic_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdtfic_omp(:,:)
      REAL,ALLOCATABLE,SAVE :: zdqfic_omp(:,:,:)
      REAL jH_cur_split,zdt_split
      LOGICAL debut_split,lafin_split
      INTEGER isplit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c$OMP THREADPRIVATE(zplev_omp,zplay_omp,zpk_omp,zphi_omp,zphis_omp,
c$OMP+                 presnivs_omp,zufi_omp,zvfi_omp,ztfi_omp,
c$OMP+                 zrfi_omp,zqfi_omp,zdufi_omp,zdvfi_omp,
c$OMP+                 zdtfi_omp,zdqfi_omp,zdpsrf_omp,flxwfi_omp,
c$OMP+                 zdufic_omp,zdvfic_omp,zdtfic_omp,zdqfic_omp)       

      LOGICAL,SAVE :: first_omp=.true.
c$OMP THREADPRIVATE(first_omp)
      
      REAL zsin(iim),zcos(iim),z1(iim)
      REAL zsinbis(iim),zcosbis(iim),z1bis(iim)
      REAL unskap, pksurcp
      save unskap

      REAL SSUM

      LOGICAL,SAVE :: firstcal=.true., debut=.true.
c$OMP THREADPRIVATE(firstcal,debut)
      
      REAL,SAVE,dimension(1:iim,1:llm):: du_send,du_recv,dv_send,dv_recv
      INTEGER :: ierr
#ifdef CPP_MPI
      INTEGER,dimension(MPI_STATUS_SIZE,4) :: Status
#else
      INTEGER,dimension(1,4) :: Status
#endif
      INTEGER, dimension(4) :: Req
      REAL,ALLOCATABLE,SAVE:: zdufi2(:,:),zdvfi2(:,:)
      integer :: k,kstart,kend
      INTEGER :: offset  
      INTEGER :: jjb,jje

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


      klon=klon_mpi
      
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

c$OMP MASTER
        flag_moyzon = .false.
        if(moyzon_ch.or.moyzon_mu) then
         flag_moyzon = .true.
         allocate(tmpvar(iip1,llm))
         allocate(tmpvarp1(iip1,llmp1))
         allocate(tmpvarbar(llm))
         allocate(tmpvarbarp1(llmp1))
        endif

      ALLOCATE(zpsrf(klon))
      ALLOCATE(zplev(klon,llm+1),zplay(klon,llm))
      ALLOCATE(zphi(klon,llm),zphis(klon))
      ALLOCATE(zufi(klon,llm), zvfi(klon,llm))
      ALLOCATE(ztfi(klon,llm),zqfi(klon,llm,nqtot))
!      ALLOCATE(pcvgu(klon,llm), pcvgv(klon,llm))
!      ALLOCATE(pcvgt(klon,llm), pcvgq(klon,llm,2))
      ALLOCATE(zdufi(klon,llm),zdvfi(klon,llm),zrfi(klon,llm))
      ALLOCATE(zdtfi(klon,llm),zdqfi(klon,llm,nqtot))
      ALLOCATE(zdpsrf(klon))
      ALLOCATE(zdufi2(klon+iim,llm),zdvfi2(klon+iim,llm))
      ALLOCATE(flxwfi(klon,llm))
! ADAPTATION GCM POUR CP(T)
      ALLOCATE(zteta(klon,llm))
      ALLOCATE(zpk(klon,llm))

      ! zonal means. horizontal dimension should be iip1
      if (flag_moyzon) call moyzon_init(klon,llm,nqtot)

c------------------------------------------------------------------
c moyennes globales pour les profils de pression et de temperature
      if(planet_type.eq."titan".or.planet_type.eq."venus") then
        call AllGather_Field(pp,iip1*jjp1,llmp1)
        call AllGather_Field(pteta,iip1*jjp1,llm)
        call AllGather_Field(ppk,iip1*jjp1,llm)
        call AllGather_Field(pphi,iip1*jjp1,llm)
        call AllGather_Field(pphis,iip1*jjp1,1)
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
        call tpot2t_p(1,llm,tetamoy,tmoy,pkmoy)
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
        allocate(klat(klon))
        do ig0=1,klon
          j=index_j(ig0)
          klat(ig0)=j
        enddo
      endif   ! planet_type=titan
c------------------------------------------------------------------

c$OMP END MASTER
c$OMP BARRIER	  
      ELSE
          debut = .FALSE.
      ENDIF


c-----------------------------------------------------------------------
c   40. transformation des variables dynamiques en variables physiques:
c   ---------------------------------------------------------------

c   41. pressions au sol (en Pascals)
c   ----------------------------------

c$OMP MASTER
      call start_timer(timer_physic)
c$OMP END MASTER

c$OMP MASTER             
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
      do ig0=1,klon
        i=index_i(ig0)
        j=index_j(ig0)
        zpsrf(ig0)=pps(i,j)
      enddo
c$OMP END MASTER


c   42. pression intercouches et fonction d'Exner:

c   -----------------------------------------------------------------
c     .... zplev  definis aux (llm +1) interfaces des couches  ....
c     .... zplay  definis aux (  llm )    milieux des couches  .... 
c   -----------------------------------------------------------------

c    ...    Exner = cp * ( p(l) / preff ) ** kappa     ....

c      print *,omp_rank,'klon--->',klon
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llmp1
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
        do ig0=1,klon
          i=index_i(ig0)
          j=index_j(ig0)
          zplev( ig0,l ) = pp(i,j,l)
        enddo
      ENDDO
c$OMP END DO NOWAIT
      if (flag_moyzon) then
        call AllGather_Field(pp,iip1*jjp1,llmp1)
        j=index_j(1)
        tmpvarp1(:,:) = pp(:,j,:)
        call moyzon(llmp1,tmpvarp1,tmpvarbarp1)
        zplevbar_mpi(1,:) = tmpvarbarp1
        do ig0=2,klon
          j=index_j(ig0)
          if (j.ne.index_j(ig0-1)) then
            tmpvarp1(:,:) = pp(:,j,:)
            call moyzon(llmp1,tmpvarp1,tmpvarbarp1)
            zplevbar_mpi(ig0,:) = tmpvarbarp1
          else
            zplevbar_mpi(ig0,:) = zplevbar_mpi(ig0-1,:)
          endif
        enddo
      endif

! ADAPTATION GCM POUR CP(T)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
        do ig0=1,klon
          i=index_i(ig0)
          j=index_j(ig0)
          zpk(ig0,l)=ppk(i,j,l)
          zteta(ig0,l)=pteta(i,j,l)
        enddo
      ENDDO
c$OMP END DO NOWAIT
      if (flag_moyzon) then
        call AllGather_Field(pteta,iip1*jjp1,llm)
        call AllGather_Field(ppk,iip1*jjp1,llm)
        j=index_j(1)
        tmpvar(:,:) = pteta(:,j,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        ztetabar_mpi(1,:) = tmpvarbar
        tmpvar(:,:) = ppk(:,j,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zpkbar_mpi(1,:) = tmpvarbar
        call tpot2t_p(1,llm,ztetabar_mpi(1,:),ztfibar_mpi(1,:),
     &                      zpkbar_mpi(1,:))
        do ig0=2,klon
          j=index_j(ig0)
          if (j.ne.index_j(ig0-1)) then
            tmpvar(:,:) = pteta(:,j,:)
            call moyzon(llm,tmpvar,tmpvarbar)
            ztetabar_mpi(ig0,:) = tmpvarbar
            tmpvar(:,:) = ppk(:,j,:)
            call moyzon(llm,tmpvar,tmpvarbar)
            zpkbar_mpi(ig0,:) = tmpvarbar
            call tpot2t_p(1,llm,ztetabar_mpi(ig0,:),ztfibar_mpi(ig0,:),
     &                          zpkbar_mpi(ig0,:))
          else
            zpkbar_mpi(ig0,:)   = zpkbar_mpi(ig0-1,:)
            ztetabar_mpi(ig0,:) = ztetabar_mpi(ig0-1,:)
            ztfibar_mpi(ig0,:)  = ztfibar_mpi(ig0-1,:)
          endif
        enddo
      endif

c   43. temperature naturelle (en K) et pressions milieux couches .
c   ---------------------------------------------------------------

! ADAPTATION GCM POUR CP(T)
         call tpot2t_p(klon,llm,zteta,ztfi,zpk)

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
        do ig0=1,klon
          i=index_i(ig0)
          j=index_j(ig0)
          pksurcp        = ppk(i,j,l) / cpp
          zplay(ig0,l)   = preff * pksurcp ** unskap
        enddo
      ENDDO
c$OMP END DO NOWAIT
      if (flag_moyzon) then
        zplaybar_mpi(:,:) = preff * (zpkbar_mpi(:,:)/cpp)**unskap
      endif

c   43.bis traceurs (tous intensifs)
c   ---------------

      DO iq=1,nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
           do ig0=1,klon
             i=index_i(ig0)
             j=index_j(ig0)
             zqfi(ig0,l,iq)  = pq(i,j,l,iq)
           enddo
         ENDDO
c$OMP END DO NOWAIT	 
      ENDDO ! of DO iq=1,nqtot
      if (flag_moyzon) then
       DO iq=1,nqtot
         call AllGather_Field(pq(:,:,:,iq),iip1*jjp1,llm)
         j=index_j(1)
         tmpvar(:,:) = pq(:,j,:,iq)
         call moyzon(llm,tmpvar,tmpvarbar)
         zqfibar_mpi(1,:,iq) = tmpvarbar
         do ig0=2,klon
          j=index_j(ig0)
          if (j.ne.index_j(ig0-1)) then
            tmpvar(:,:) = pq(:,j,:,iq)
            call moyzon(llm,tmpvar,tmpvarbar)
            zqfibar_mpi(ig0,:,iq) = tmpvarbar
          else
            zqfibar_mpi(ig0,:,iq) = zqfibar_mpi(ig0-1,:,iq)
          endif
        enddo
       ENDDO ! of DO iq=1,nqtot
      endif


c   Geopotentiel calcule par rapport a la surface locale:
c   -----------------------------------------------------

      CALL gr_dyn_fi_p(llm,iip1,jjp1,klon,pphi,zphi)

      CALL gr_dyn_fi_p(1,iip1,jjp1,klon,pphis,zphis)

c$OMP BARRIER
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
	 DO ig=1,klon
	   zphi(ig,l)=zphi(ig,l)-zphis(ig)
	 ENDDO
      ENDDO
c$OMP END DO NOWAIT
      
      if (flag_moyzon) then
        call AllGather_Field(pphis,iip1*jjp1,1)
        call AllGather_Field(pphi,iip1*jjp1,llm)
        j=index_j(1)
        tmpvar(:,1) = pphis(:,j)
        call moyzon(1,tmpvar(:,1),tmpvarbar(1))
        zphisbar_mpi(1) = tmpvarbar(1)
        tmpvar(:,:) = pphi(:,j,:)
        call moyzon(llm,tmpvar,tmpvarbar)
        zphibar_mpi(1,:) = tmpvarbar-zphisbar_mpi(1)
        do ig0=2,klon
          j=index_j(ig0)
          if (j.ne.index_j(ig0-1)) then
            tmpvar(:,1) = pphis(:,j)
            call moyzon(1,tmpvar(:,1),tmpvarbar(1))
            zphisbar_mpi(ig0) = tmpvarbar(1)
            tmpvar(:,:) = pphi(:,j,:)
            call moyzon(llm,tmpvar,tmpvarbar)
            zphibar_mpi(ig0,:) = tmpvarbar-zphisbar_mpi(ig0)
          else
            zphisbar_mpi(ig0)  = zphisbar_mpi(ig0-1)
            zphibar_mpi(ig0,:) = zphibar_mpi(ig0-1,:)
          endif
        enddo
      endif

c
c   45. champ u:
c   ------------

      kstart=1
      kend=klon
      
      if (is_north_pole_dyn) kstart=2
      if (is_south_pole_dyn) kend=klon-1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
!CDIR SPARSE
        do ig0=kstart,kend
          i=index_i(ig0)
          j=index_j(ig0)
          if (i==1) then
            zufi(ig0,l)= 0.5 *(  pucov(iim,j,l)/cu(iim,j)
     $                         + pucov(1,j,l)/cu(1,j) )
          else
            zufi(ig0,l)= 0.5*(  pucov(i-1,j,l)/cu(i-1,j) 
     $                       + pucov(i,j,l)/cu(i,j) )
          endif
        enddo
      ENDDO
c$OMP END DO NOWAIT

c
C  Alvaro de la Camara (May 2014)
C  46.1 Calcul de la vorticite et passage sur la grille physique
C  --------------------------------------------------------------

      jjb=jj_begin_dyn-1
      jje=jj_end_dyn+1
      if (is_north_pole_dyn) jjb=1
      if (is_south_pole_dyn) jje=jjm

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)

      DO l=1,llm
        do i=1,iim
          do j=jjb,jje
            zrot(i,j,l) = (pvcov(i+1,j,l) - pvcov(i,j,l)
     $                   + pucov(i,j+1,l) - pucov(i,j,l))
     $                   / (cu(i,j)+cu(i,j+1))
     $                   / (cv(i+1,j)+cv(i,j)) *4
          enddo
        enddo
      ENDDO


c   46.2champ v:
c   -----------

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
        DO ig0=kstart,kend
          i=index_i(ig0)
          j=index_j(ig0)
          zvfi(ig0,l)= 0.5 *(  pvcov(i,j-1,l)/cv(i,j-1) 
     $                       + pvcov(i,j,l)/cv(i,j) )
    
          if (j==1 .OR. j==jjp1) then !  AdlC MAY 2014
            zrfi(ig0,l) = 0 !  AdlC MAY 2014
          else
            if(i==1)then
            zrfi(ig0,l)= 0.25 *(zrot(iim,j-1,l)+zrot(iim,j,l)
     $                   +zrot(1,j-1,l)+zrot(1,j,l))   !  AdlC MAY 2014
            else
            zrfi(ig0,l)= 0.25 *(zrot(i-1,j-1,l)+zrot(i-1,j,l)
     $                   +zrot(i,j-1,l)+zrot(i,j,l))   !  AdlC MAY 2014
            endif
          endif


         ENDDO
      ENDDO
c$OMP END DO NOWAIT

c   47. champs de vents aux pole nord   
c   ------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      if (is_north_pole_dyn) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO l=1,llm

           z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,1,l)/cv(1,1)
           DO i=2,iim
              z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,1,l)/cv(i,1)
           ENDDO
  
           DO i=1,iim
              zcos(i)   = COS(rlonv(i))*z1(i)
              zsin(i)   = SIN(rlonv(i))*z1(i)
           ENDDO
  
           zufi(1,l)  = SSUM(iim,zcos,1)/pi
           zvfi(1,l)  = SSUM(iim,zsin,1)/pi
           zrfi(1,l)  = 0.
  
        ENDDO
c$OMP END DO NOWAIT      
      endif


c   48. champs de vents aux pole sud:
c   ---------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      if (is_south_pole_dyn) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO l=1,llm
  
         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,jjm,l)/cv(1,jjm)
           DO i=2,iim
             z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,jjm,l)/cv(i,jjm)
	   ENDDO
  
           DO i=1,iim
              zcos(i)    = COS(rlonv(i))*z1(i)
              zsin(i)    = SIN(rlonv(i))*z1(i)
	   ENDDO
  
           zufi(klon,l)  = SSUM(iim,zcos,1)/pi
           zvfi(klon,l)  = SSUM(iim,zsin,1)/pi
           zrfi(klon,l)  = 0.
        ENDDO
c$OMP END DO NOWAIT       
      endif

c On change de grille, dynamique vers physiq, pour le flux de masse verticale
      CALL gr_dyn_fi_p(llm,iip1,jjp1,klon,flxw,flxwfi)

c-----------------------------------------------------------------------
c   Appel de la physique:
c   ---------------------

! Appel de la physique: pose probleme quand on tourne 
! SANS physique, car physiq.F est dans le repertoire phy[]... 
! Il faut une cle CPP_PHYS

! Le fait que les arguments de physiq soient differents selon les planetes 
! ne pose pas de probleme a priori.


c$OMP BARRIER
      if (first_omp) then
        klon=klon_omp

        allocate(zplev_omp(klon,llm+1))
        allocate(zplay_omp(klon,llm))
        allocate(zpk_omp(klon,llm))
        allocate(zphi_omp(klon,llm))
        allocate(zphis_omp(klon))
        allocate(presnivs_omp(llm))
        allocate(zufi_omp(klon,llm))
        allocate(zvfi_omp(klon,llm))
        allocate(zrfi_omp(klon,llm))  ! LG Ari 2014
        allocate(ztfi_omp(klon,llm))
        allocate(zqfi_omp(klon,llm,nqtot))
        allocate(zdufi_omp(klon,llm))
        allocate(zdvfi_omp(klon,llm))
        allocate(zdtfi_omp(klon,llm))
        allocate(zdqfi_omp(klon,llm,nqtot))
        allocate(zdufic_omp(klon,llm))
        allocate(zdvfic_omp(klon,llm))
        allocate(zdtfic_omp(klon,llm))
        allocate(zdqfic_omp(klon,llm,nqtot))
        allocate(zdpsrf_omp(klon))
        allocate(flxwfi_omp(klon,llm))

        if (flag_moyzon) call moyzon_init_omp(klon,llm,nqtot)

	first_omp=.false.
      endif
       
	   
      klon=klon_omp
      offset=klon_omp_begin-1
      
      do l=1,llm+1
        do i=1,klon
          zplev_omp(i,l)=zplev(offset+i,l)
	enddo 
      enddo
	  
       do l=1,llm
        do i=1,klon  
	  zplay_omp(i,l)=zplay(offset+i,l)
	enddo 
      enddo
	
       do l=1,llm
        do i=1,klon  
	  zpk_omp(i,l)=zpk(offset+i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zphi_omp(i,l)=zphi(offset+i,l)
	enddo 
      enddo
	
      do i=1,klon
	zphis_omp(i)=zphis(offset+i)
      enddo 
     
	
      do l=1,llm
        presnivs_omp(l)=presnivs(l)
      enddo 
	
      do l=1,llm
        do i=1,klon
	  zufi_omp(i,l)=zufi(offset+i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zvfi_omp(i,l)=zvfi(offset+i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zrfi_omp(i,l)=zrfi(offset+i,l)
	enddo 
      enddo
	
	
      do l=1,llm
        do i=1,klon
	  ztfi_omp(i,l)=ztfi(offset+i,l)
	enddo 
      enddo
	
      do iq=1,nqtot
        do l=1,llm
          do i=1,klon
            zqfi_omp(i,l,iq)=zqfi(offset+i,l,iq)
	  enddo
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zdufi_omp(i,l)=zdufi(offset+i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zdvfi_omp(i,l)=zdvfi(offset+i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
          zdtfi_omp(i,l)=zdtfi(offset+i,l)
	enddo 
      enddo
	
      do iq=1,nqtot
        do l=1,llm
          do i=1,klon
	    zdqfi_omp(i,l,iq)=zdqfi(offset+i,l,iq)
	  enddo 
        enddo
      enddo
      	
      do i=1,klon
	zdpsrf_omp(i)=zdpsrf(offset+i)
      enddo 

      do l=1,llm
        do i=1,klon
          flxwfi_omp(i,l)=flxwfi(offset+i,l)
	enddo 
      enddo

      if (flag_moyzon) then
       do l=1,llm+1
        do i=1,klon
          zplevbar(i,l)=zplevbar_mpi(offset+i,l)
	enddo 
       enddo
	  
       do l=1,llm
        do i=1,klon  
	  zplaybar(i,l)=zplaybar_mpi(offset+i,l)
	enddo 
       enddo
	
       do l=1,llm
        do i=1,klon
	  zphibar(i,l)=zphibar_mpi(offset+i,l)
	enddo 
       enddo
		
      do i=1,klon
	zphisbar(i)=zphisbar_mpi(offset+i)
      enddo 
     
      do l=1,llm
        do i=1,klon
	  ztfibar(i,l)=ztfibar_mpi(offset+i,l)
	enddo 
      enddo
	
      do iq=1,nqtot
        do l=1,llm
          do i=1,klon
            zqfibar(i,l,iq)=zqfibar_mpi(offset+i,l,iq)
	  enddo
	enddo 
       enddo
      endif

      
c$OMP BARRIER

!$OMP MASTER
!      write(lunout,*) 'PHYSIQUE AVEC NSPLIT_PHYS=',nsplit_phys
!$OMP END MASTER
      zdt_split=dtphys/nsplit_phys
      zdufic_omp(:,:)=0.
      zdvfic_omp(:,:)=0.
      zdtfic_omp(:,:)=0.
      zdqfic_omp(:,:,:)=0.

#ifdef CPP_PHYS
      do isplit=1,nsplit_phys

         jH_cur_split=jH_cur+(isplit-1) * dtvr / (daysec *nsplit_phys)
         debut_split=debut.and.isplit==1
         lafin_split=lafin.and.isplit==nsplit_phys

        CALL call_physiq(klon,llm,nqtot,tname,
     &                   debut_split,lafin_split,
     &                   jD_cur,jH_cur_split,zdt_split,
     &                   zplev_omp,zplay_omp,
     &                   zpk_omp,zphi_omp,zphis_omp,
     &                   presnivs_omp,
     &                   zufi_omp,zvfi_omp,zrfi_omp,ztfi_omp,zqfi_omp,
     &                   flxwfi_omp,pducov,
     &                   zdufi_omp,zdvfi_omp,zdtfi_omp,zdqfi_omp,
     &                   zdpsrf_omp)

!      if (planet_type=="earth") then
!        CALL physiq (klon,
!     .             llm,
!     .             debut_split,
!     .             lafin_split,
!     .             jD_cur,
!     .             jH_cur_split,
!     .             zdt_split,
!     .             zplev_omp,
!     .             zplay_omp,
!     .             zphi_omp,
!     .             zphis_omp,
!     .             presnivs_omp,
!     .             zufi_omp,
!     .             zvfi_omp,
!     .             ztfi_omp,
!     .             zqfi_omp,
!     .             flxwfi_omp,
!     .             zdufi_omp,
!     .             zdvfi_omp,
!     .             zdtfi_omp,
!     .             zdqfi_omp,
!     .             zdpsrf_omp,
!     .             pducov)
!
!      else if ( planet_type=="generic" ) then
!
!      CALL physiq (klon,     !! ngrid
!     .             llm,            !! nlayer
!     .             nqtot,          !! nq
!     .             tname,          !! tracer names from dynamical core (given in infotrac)
!     .             debut_split,    !! firstcall 
!     .             lafin_split,    !! lastcall
!     .             jD_cur,         !! pday. see leapfrog_p 
!     .             jH_cur_split,   !! ptime "fraction of day"
!     .             zdt_split,      !! ptimestep
!     .             zplev_omp,  !! pplev
!     .             zplay_omp,  !! pplay
!     .             zphi_omp,   !! pphi
!     .             zufi_omp,   !! pu
!     .             zvfi_omp,   !! pv
!     .             ztfi_omp,   !! pt
!     .             zqfi_omp,   !! pq
!     .             flxwfi_omp, !! pw !! or 0. anyway this is for diagnostic. not used in physiq.
!     .             zdufi_omp,  !! pdu
!     .             zdvfi_omp,  !! pdv
!     .             zdtfi_omp,  !! pdt
!     .             zdqfi_omp,  !! pdq
!     .             zdpsrf_omp, !! pdpsrf
!     .             tracerdyn)      !! tracerdyn <-- utilite ???
!
!      else if ( planet_type=="mars" ) then
!
!        CALL physiq (klon,       ! ngrid
!     .             llm,          ! nlayer
!     .             nqtot,        ! nq
!     .             debut_split,  ! firstcall
!     .             lafin_split,  ! lastcall
!     .             jD_cur,       ! pday
!     .             jH_cur_split, ! ptime
!     .             zdt_split,    ! ptimestep
!     .             zplev_omp,    ! pplev
!     .             zplay_omp,    ! pplay
!     .             zphi_omp,     ! pphi
!     .             zufi_omp,     ! pu
!     .             zvfi_omp,     ! pv
!     .             ztfi_omp,     ! pt
!     .             zqfi_omp,     ! pq
!     .             flxwfi_omp,   ! pw
!     .             zdufi_omp,    ! pdu
!     .             zdvfi_omp,    ! pdv
!     .             zdtfi_omp,    ! pdt
!     .             zdqfi_omp,    ! pdq
!     .             zdpsrf_omp,   ! pdpsrf
!     .             tracerdyn)    ! tracerdyn (somewhat obsolete)
!
!      else if ((planet_type=="titan").or.(planet_type=="venus")) then
!
!        CALL physiq (klon,
!     .             llm,
!     .             nqtot,
!     .             debut_split,
!     .             lafin_split,
!     .             jD_cur,
!     .             jH_cur_split,
!     .             zdt_split,
!     .             zplev_omp,
!     .             zplay_omp,
!     .             zpk_omp,
!     .             zphi_omp,
!     .             zphis_omp,
!     .             presnivs_omp,
!     .             zufi_omp,
!     .             zvfi_omp,
!     .             ztfi_omp,
!     .             zqfi_omp,
!     .             flxwfi_omp,
!     .             zdufi_omp,
!     .             zdvfi_omp,
!     .             zdtfi_omp,
!     .             zdqfi_omp,
!     .             zdpsrf_omp)
!
!      else ! unknown "planet_type"
!
!        write(lunout,*) "calfis_p: error, unknown planet_type: ",
!     &                  trim(planet_type)
!        stop
!
!      endif ! planet_type
         zufi_omp(:,:)=zufi_omp(:,:)+zdufi_omp(:,:)*zdt_split
         zvfi_omp(:,:)=zvfi_omp(:,:)+zdvfi_omp(:,:)*zdt_split
         ztfi_omp(:,:)=ztfi_omp(:,:)+zdtfi_omp(:,:)*zdt_split
         zqfi_omp(:,:,:)=zqfi_omp(:,:,:)+zdqfi_omp(:,:,:)*zdt_split

         zdufic_omp(:,:)=zdufic_omp(:,:)+zdufi_omp(:,:)
         zdvfic_omp(:,:)=zdvfic_omp(:,:)+zdvfi_omp(:,:)
         zdtfic_omp(:,:)=zdtfic_omp(:,:)+zdtfi_omp(:,:)
         zdqfic_omp(:,:,:)=zdqfic_omp(:,:,:)+zdqfi_omp(:,:,:)

      enddo

#endif
! of #ifdef CPP_PHYS

      zdufi_omp(:,:)=zdufic_omp(:,:)/nsplit_phys
      zdvfi_omp(:,:)=zdvfic_omp(:,:)/nsplit_phys
      zdtfi_omp(:,:)=zdtfic_omp(:,:)/nsplit_phys
      zdqfi_omp(:,:,:)=zdqfic_omp(:,:,:)/nsplit_phys

c$OMP BARRIER

      do l=1,llm+1
        do i=1,klon
          zplev(offset+i,l)=zplev_omp(i,l)
	enddo 
      enddo
	  
       do l=1,llm
        do i=1,klon  
	  zplay(offset+i,l)=zplay_omp(i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zphi(offset+i,l)=zphi_omp(i,l)
	enddo 
      enddo
	

      do i=1,klon
	zphis(offset+i)=zphis_omp(i)
      enddo 
     
	
      do l=1,llm
        presnivs(l)=presnivs_omp(l)
      enddo 
	
      do l=1,llm
        do i=1,klon
	  zufi(offset+i,l)=zufi_omp(i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zvfi(offset+i,l)=zvfi_omp(i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  ztfi(offset+i,l)=ztfi_omp(i,l)
	enddo 
      enddo
	
      do iq=1,nqtot
        do l=1,llm
          do i=1,klon
            zqfi(offset+i,l,iq)=zqfi_omp(i,l,iq)
	  enddo
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zdufi(offset+i,l)=zdufi_omp(i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
	  zdvfi(offset+i,l)=zdvfi_omp(i,l)
	enddo 
      enddo
	
      do l=1,llm
        do i=1,klon
          zdtfi(offset+i,l)=zdtfi_omp(i,l)
	enddo 
      enddo
	
      do iq=1,nqtot
        do l=1,llm
          do i=1,klon
	    zdqfi(offset+i,l,iq)=zdqfi_omp(i,l,iq)
	  enddo 
        enddo
      enddo
      	
      do i=1,klon
	zdpsrf(offset+i)=zdpsrf_omp(i)
      enddo 
      

      klon=klon_mpi
500   CONTINUE
c$OMP BARRIER

c$OMP MASTER
      call stop_timer(timer_physic)
c$OMP END MASTER

      IF (using_mpi) THEN
            
      if (MPI_rank>0) then

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
       DO l=1,llm      
        du_send(1:iim,l)=zdufi(1:iim,l)
        dv_send(1:iim,l)=zdvfi(1:iim,l)
       ENDDO
c$OMP END DO NOWAIT       

c$OMP BARRIER
#ifdef CPP_MPI 
c$OMP MASTER
!$OMP CRITICAL (MPI)
        call MPI_ISSEND(du_send,iim*llm,MPI_REAL8,MPI_Rank-1,401,
     &                   COMM_LMDZ,Req(1),ierr)
        call MPI_ISSEND(dv_send,iim*llm,MPI_REAL8,MPI_Rank-1,402,
     &                  COMM_LMDZ,Req(2),ierr)
!$OMP END CRITICAL (MPI)
c$OMP END MASTER
#endif
c$OMP BARRIER
     
      endif
   
      if (MPI_rank<MPI_Size-1) then
c$OMP BARRIER
#ifdef CPP_MPI 
c$OMP MASTER      
!$OMP CRITICAL (MPI)
        call MPI_IRECV(du_recv,iim*llm,MPI_REAL8,MPI_Rank+1,401,
     &                 COMM_LMDZ,Req(3),ierr)
        call MPI_IRECV(dv_recv,iim*llm,MPI_REAL8,MPI_Rank+1,402,
     &                 COMM_LMDZ,Req(4),ierr)
!$OMP END CRITICAL (MPI)
c$OMP END MASTER
#endif
      endif

c$OMP BARRIER


#ifdef CPP_MPI 
c$OMP MASTER    
!$OMP CRITICAL (MPI)
      if (MPI_rank>0 .and. MPI_rank< MPI_Size-1) then
        call MPI_WAITALL(4,Req(1),Status,ierr)
      else if (MPI_rank>0) then
        call MPI_WAITALL(2,Req(1),Status,ierr)
      else if (MPI_rank <MPI_Size-1) then
        call MPI_WAITALL(2,Req(3),Status,ierr)
      endif
!$OMP END CRITICAL (MPI)
c$OMP END MASTER
#endif

c$OMP BARRIER     

      ENDIF ! using_mpi
      
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
            
        zdufi2(1:klon,l)=zdufi(1:klon,l)
        zdufi2(klon+1:klon+iim,l)=du_recv(1:iim,l)
            
        zdvfi2(1:klon,l)=zdvfi(1:klon,l)
        zdvfi2(klon+1:klon+iim,l)=dv_recv(1:iim,l) 
  
         pdhfi(:,jj_begin,l)=0
         pdqfi(:,jj_begin,l,:)=0
         pdufi(:,jj_begin,l)=0
         pdvfi(:,jj_begin,l)=0
         
         if (.not. is_south_pole_dyn) then
           pdhfi(:,jj_end,l)=0
           pdqfi(:,jj_end,l,:)=0
           pdufi(:,jj_end,l)=0
           pdvfi(:,jj_end,l)=0
         endif
      
       ENDDO 
c$OMP END DO NOWAIT

c$OMP MASTER
       pdpsfi(:,jj_begin)=0    
       if (.not. is_south_pole_dyn) then
	 pdpsfi(:,jj_end)=0
       endif
c$OMP END MASTER
c-----------------------------------------------------------------------
c   transformation des tendances physiques en tendances dynamiques:
c   ---------------------------------------------------------------

c  tendance sur la pression :
c  -----------------------------------
      CALL gr_fi_dyn_p(1,klon,iip1,jjp1,zdpsrf,pdpsfi)
c
c   62. enthalpie potentielle
c   ---------------------

      kstart=1
      kend=klon

      if (is_north_pole_dyn) kstart=2
      if (is_south_pole_dyn)  kend=klon-1
      
! ADAPTATION GCM POUR CP(T)
      call t2tpot_p(klon,llm,ztfi,zteta,zpk)


c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
!cdir NODEP
        do ig0=kstart,kend
          i=index_i(ig0)
          j=index_j(ig0)
!          pdhfi(i,j,l) = cpp * zdtfi(ig0,l) / ppk(i,j,l)
          pdhfi(i,j,l) = (zteta(ig0,l) - pteta(i,j,l))/dtphys
!          if (i==1) pdhfi(iip1,j,l) =  cpp * zdtfi(ig0,l) / ppk(i,j,l)
          if (i==1) then
            pdhfi(iip1,j,l) = (zteta(ig0,l) - pteta(i,j,l))/dtphys
          endif
         enddo          

        if (is_north_pole_dyn) then
            DO i=1,iip1
!              pdhfi(i,1,l)    = cpp *  zdtfi(1,l)      / ppk(i, 1  ,l)
              pdhfi(i,1,l)    = (zteta(1,l) - pteta(i,1,l))/dtphys
            enddo
        endif
        
        if (is_south_pole_dyn) then
            DO i=1,iip1
!              pdhfi(i,jjp1,l) = cpp *  zdtfi(klon,l)/ ppk(i,jjp1,l)
              pdhfi(i,jjp1,l) = (zteta(klon,l) - pteta(i,jjp1,l))/dtphys
            ENDDO
        endif
      ENDDO
c$OMP END DO NOWAIT
      
c   62. humidite specifique
c   ---------------------
! Ehouarn: removed this useless bit: was overwritten at step 63 anyways
!      DO iq=1,nqtot
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
!         DO l=1,llm
!!!cdir NODEP 
!           do ig0=kstart,kend
!             i=index_i(ig0)
!             j=index_j(ig0)
!             pdqfi(i,j,l,iq) = zdqfi(ig0,l,iq) 
!             if (i==1) pdqfi(iip1,j,l,iq) = zdqfi(ig0,l,iq) 
!           enddo
!           
!           if (is_north_pole_dyn) then
!             do i=1,iip1
!               pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)             
!             enddo
!           endif
!           
!           if (is_south_pole_dyn) then
!             do i=1,iip1
!               pdqfi(i,jjp1,l,iq) = zdqfi(klon,l,iq) 
!             enddo
!           endif
!         ENDDO
!c$OMP END DO NOWAIT
!      ENDDO

c   63. traceurs (tous en intensifs)
c   ------------
C     initialisation des tendances

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        pdqfi(:,:,l,:)=0.
      ENDDO
c$OMP END DO NOWAIT	 

C
!cdir NODEP
      DO iq=1,nqtot
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
!cdir NODEP           
	     DO ig0=kstart,kend
              i=index_i(ig0)
              j=index_j(ig0)
              pdqfi(i,j,l,iq) = zdqfi(ig0,l,iq)
              if (i==1) pdqfi(iip1,j,l,iq) = zdqfi(ig0,l,iq)
            ENDDO
	    
	    IF (is_north_pole_dyn) then
	      DO i=1,iip1
                pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)
	      ENDDO
	    ENDIF
	    
	    IF (is_south_pole_dyn) then
	      DO i=1,iip1
                pdqfi(i,jjp1,l,iq) = zdqfi(klon,l,iq)
	      ENDDO
	    ENDIF
	    
         ENDDO
c$OMP END DO NOWAIT	 
      ENDDO
      
c   65. champ u:
c   ------------
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
!cdir NODEP
         do ig0=kstart,kend
           i=index_i(ig0)
           j=index_j(ig0)
           
           if (i/=iim) then
             pdufi(i,j,l)=0.5*(zdufi2(ig0,l)+zdufi2(ig0+1,l))*cu(i,j)
           endif
           
           if (i==1) then
              pdufi(iim,j,l)=0.5*(  zdufi2(ig0,l)
     $                            + zdufi2(ig0+iim-1,l))*cu(iim,j)
             pdufi(iip1,j,l)=0.5*(zdufi2(ig0,l)+zdufi2(ig0+1,l))*cu(i,j)
           endif
         
         enddo
         
         if (is_north_pole_dyn) then
           DO i=1,iip1
            pdufi(i,1,l)    = 0.
           ENDDO
         endif
         
         if (is_south_pole_dyn) then
           DO i=1,iip1
            pdufi(i,jjp1,l) = 0.
           ENDDO
         endif
         
      ENDDO
c$OMP END DO NOWAIT

c   67. champ v:
c   ------------

      kstart=1
      kend=klon

      if (is_north_pole_dyn) kstart=2
      if (is_south_pole_dyn)  kend=klon-1-iim
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l=1,llm
!CDIR ON_ADB(index_i)
!CDIR ON_ADB(index_j) 
!cdir NODEP
        do ig0=kstart,kend
           i=index_i(ig0)
           j=index_j(ig0)
           pdvfi(i,j,l)=0.5*(zdvfi2(ig0,l)+zdvfi2(ig0+iim,l))*cv(i,j)
           if (i==1) pdvfi(iip1,j,l) = 0.5*(zdvfi2(ig0,l)+
     $	                                    zdvfi2(ig0+iim,l))
     $				          *cv(i,j)
        enddo
         
      ENDDO
c$OMP END DO NOWAIT


c   68. champ v pres des poles:
c   ---------------------------
c      v = U * cos(long) + V * SIN(long)

      if (is_north_pole_dyn) then

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO l=1,llm

          DO i=1,iim
            pdvfi(i,1,l)=
     $      zdufi(1,l)*COS(rlonv(i))+zdvfi(1,l)*SIN(rlonv(i))
       
            pdvfi(i,1,l)=
     $      0.5*(pdvfi(i,1,l)+zdvfi(i+1,l))*cv(i,1)
          ENDDO

          pdvfi(iip1,1,l)  = pdvfi(1,1,l)

        ENDDO
c$OMP END DO NOWAIT

      endif    
      
      if (is_south_pole_dyn) then

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
         DO l=1,llm
  
           DO i=1,iim
              pdvfi(i,jjm,l)=zdufi(klon,l)*COS(rlonv(i))
     $        +zdvfi(klon,l)*SIN(rlonv(i))

              pdvfi(i,jjm,l)=
     $        0.5*(pdvfi(i,jjm,l)+zdvfi(klon-iip1+i,l))*cv(i,jjm)
           ENDDO

           pdvfi(iip1,jjm,l)= pdvfi(1,jjm,l)

        ENDDO
c$OMP END DO NOWAIT
     
      endif
c-----------------------------------------------------------------------

700   CONTINUE
 
      firstcal = .FALSE.

#else
      write(lunout,*)
     & "calfis_p: for now can only work with parallel physics"
      stop
#endif
! of #ifdef CPP_PHYS
#endif
! of #ifdef CPP_PARA
      END
