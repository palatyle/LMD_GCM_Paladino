C======================================================================
      PROGRAM newstart
c=======================================================================
c
c
c   Auteur:   Christophe Hourdin/Francois Forget/Yann Wanherdrick
c   ------
c             Derniere modif : 12/03
c
c
c   Objet:  Create or modify the initial state for the LMD Mars GCM
c   -----           (fichiers NetCDF start et startfi)
c
c
c=======================================================================

      use ioipsl_getincom, only: getin
      use mod_phys_lmdz_para, only: is_parallel, is_sequential,
     &                              is_mpi_root, is_omp_root,
     &                              is_master
      use infotrac, only: infotrac_init, nqtot, tname
      use tracer_mod, only: noms, mmol,
     &                      igcm_dust_number, igcm_dust_mass,
     &                      igcm_ccn_number, igcm_ccn_mass,
     &                      igcm_h2o_vap, igcm_h2o_ice, igcm_co2,
     &                      igcm_n2, igcm_ar, igcm_o2, igcm_co,
     &                      igcm_o, igcm_h2
      use surfdat_h, only: phisfi, z0, zmea, zstd, zsig, zgam, zthe,
     &                     albedodat, z0_default, qsurf, tsurf,
     &                     co2ice, emis, hmons, summit, base, watercap
      use comsoil_h, only: inertiedat, layer, mlayer, nsoilmx, tsoil
      use control_mod, only: day_step, iphysiq, anneeref, planet_type
      use geometry_mod, only: longitude,latitude,cell_area
      use phyetat0_mod, only: phyetat0
      use phyredem, only: physdem0, physdem1
      use iostart, only: open_startphy
      use dimradmars_mod, only: tauscaling, albedo
      use turb_mod, only: q2, wstar
      use co2cloud_mod, only: mem_Mccn_co2, mem_Mh2o_co2,
     &                        mem_Nccn_co2
      use filtreg_mod, only: inifilr
      USE mod_const_mpi, ONLY: COMM_LMDZ
      USE comvert_mod, ONLY: ap,bp,pa,preff
      USE comconst_mod, ONLY: lllm,daysec,dtphys,dtvr,
     .			cpp,kappa,rad,omeg,g,r,pi
      USE serre_mod, ONLY: alphax
      USE temps_mod, ONLY: day_ini,hour_ini
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0
      USE iniphysiq_mod, ONLY: iniphysiq
      USE exner_hyb_m, ONLY: exner_hyb

      implicit none

      include "dimensions.h"
      integer, parameter :: ngridmx = (2+(jjm-1)*iim - 1/jjm) 
      include "paramet.h"
      include "comgeom2.h"
      include "comdissnew.h"
      include "clesph0.h"
      include "netcdf.inc"
c=======================================================================
c   Declarations
c=======================================================================

c Variables dimension du fichier "start_archive"
c------------------------------------
      CHARACTER	relief*3

c et autres:
c----------

c Variables pour les lectures NetCDF des fichiers "start_archive" 
c--------------------------------------------------
      INTEGER nid_dyn, nid_fi,nid,nvarid
      INTEGER tab0

      REAL  date
      REAL p_rad,p_omeg,p_g,p_mugaz,p_daysec

c Variable histoire 
c------------------
      REAL vcov(iip1,jjm,llm),ucov(iip1,jjp1,llm) ! vents covariants
      REAL phis(iip1,jjp1)
      REAL,ALLOCATABLE :: q(:,:,:,:)               ! champs advectes

c autre variables dynamique nouvelle grille
c------------------------------------------
      REAL pks(iip1,jjp1)
      REAL w(iip1,jjp1,llm+1)
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm)
!      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
!      REAL dh(ip1jmp1,llm),dp(ip1jmp1)
      REAL phi(iip1,jjp1,llm)

      integer klatdat,klongdat
      PARAMETER (klatdat=180,klongdat=360)

c Physique sur grille scalaire 
c----------------------------
      real zmeaS(iip1,jjp1),zstdS(iip1,jjp1)
      real zsigS(iip1,jjp1),zgamS(iip1,jjp1),ztheS(iip1,jjp1)
      real hmonsS(iip1,jjp1)
      real summitS(iip1,jjp1)
      real baseS(iip1,jjp1)
      real zavgS(iip1,jjp1)
      real z0S(iip1,jjp1)

c variable physique
c------------------
      REAL tauscadyn(iip1,jjp1) ! dust conversion factor on the dynamics grid
      real alb(iip1,jjp1),albfi(ngridmx) ! albedos
      real ith(iip1,jjp1,nsoilmx),ithfi(ngridmx,nsoilmx) ! thermal inertia (3D)
      real surfith(iip1,jjp1),surfithfi(ngridmx) ! surface thermal inertia (2D)
!      REAL latfi(ngridmx),lonfi(ngridmx),airefi(ngridmx)

      INTEGER i,j,l,isoil,ig,idum
      real mugaz ! molar mass of the atmosphere

      integer ierr  !, nbetat

c Variables on the new grid along scalar points 
c------------------------------------------------------
!      REAL p(iip1,jjp1)
      REAL t(iip1,jjp1,llm)
      real phisold_newgrid(iip1,jjp1)
      REAL :: teta(iip1, jjp1, llm)
      REAL :: pk(iip1,jjp1,llm)
      REAL :: pkf(iip1,jjp1,llm)
      REAL :: ps(iip1, jjp1)
      REAL :: masse(iip1,jjp1,llm)
      REAL :: xpn,xps,xppn(iim),xpps(iim)
      REAL :: p3d(iip1, jjp1, llm+1)
!      REAL dteta(ip1jmp1,llm)

c Variable de l'ancienne grille 
c------------------------------
      real time
      real tab_cntrl(100)
      real tab_cntrl_bis(100)

c variables diverses
c-------------------
      real choix_1 ! ==0 : read start_archive file ; ==1: read start files
      character*80      fichnom
      integer Lmodif,iq
      integer flagthermo, flagh2o
      character modif*20
      real tsud,albsud,alb_bb,ith_bb,Tiso
      real ptoto,pcap,patm,airetot,ptotn,patmn
!      real ssum
      character*1 yes
      logical :: flagiso=.false. ,  flagps0=.false.
      real val, val2, val3 ! to store temporary variables
      real :: iceith=2000 ! thermal inertia of subterranean ice
      real :: iceithN,iceithS ! values of thermal inertias in N & S hemispheres
      integer iref,jref

      INTEGER :: itau
      
      INTEGER :: numvanle
      character(len=50) :: txt ! to store some text
      integer :: count
      real :: profile(llm+1) ! to store an atmospheric profile + surface value

! MONS data:
      real :: MONS_Hdn(iip1,jjp1) ! Hdn: %WEH=Mass fraction of H2O
      real :: MONS_d21(iip1,jjp1) ! ice table "depth" (in kg/m2)
      ! coefficient to apply to convert d21 to 'true' depth (m)
      real :: MONS_coeff
      real :: MONS_coeffS ! coeff for southern hemisphere
      real :: MONS_coeffN ! coeff for northern hemisphere
!      real,parameter :: icedepthmin=1.e-3 ! Ice begins at most at that depth
! Reference position for composition
      real :: latref,lonref,dlatmin,dlonmin
! Variable used to change composition
      real :: Svmr,Smmr,Smmr_old,Smmr_new,n,Sn
      real :: Mair_old,Mair_new,vmr_old,vmr_new
      real,allocatable :: coefvmr(:)  ! Correction coefficient when changing composition
      real :: maxq
      integer :: iloc(1), iqmax
! sub-grid cloud fraction
      real :: totcloudfrac(ngridmx)

c sortie visu pour les champs dynamiques
c---------------------------------------
!      INTEGER :: visuid
!      real :: time_step,t_ops,t_wrt
!      CHARACTER*80 :: visu_file

      cpp    = 744.499 ! for Mars, instead of 1004.70885 (Earth)
      preff  = 610.    ! for Mars, instead of 101325. (Earth)
      pa= 20           ! for Mars, instead of 500 (Earth)
      planet_type="mars"

! initialize "serial/parallel" related stuff:
! (required because we call tabfi() below, before calling iniphysiq)
      is_sequential=.true.
      is_parallel=.false.
      is_mpi_root=.true.
      is_omp_root=.true.
      is_master=.true.
      
! Load tracer number and names:
      call infotrac_init
! allocate arrays
      allocate(q(iip1,jjp1,llm,nqtot))
      allocate(coefvmr(nqtot))

c=======================================================================
c   Choice of the start file(s) to use
c=======================================================================

      write(*,*) 'From which kind of files do you want to create new',
     .  'start and startfi files'
      write(*,*) '    0 - from a file start_archive'
      write(*,*) '    1 - from files start and startfi'
 
c-----------------------------------------------------------------------
c   Open file(s) to modify (start or start_archive)
c-----------------------------------------------------------------------

      DO
         read(*,*,iostat=ierr) choix_1
         if ((choix_1 /= 0).OR.(choix_1 /=1)) EXIT
      ENDDO

c     Open start_archive
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (choix_1.eq.0) then

        write(*,*) 'Creating start files from:'
        write(*,*) './start_archive.nc'
        write(*,*)
        fichnom = 'start_archive.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Problem opening file:',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF
        tab0 = 50 
        Lmodif = 1

c     OR open start and startfi files
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else
        write(*,*) 'Creating start files from:'
        write(*,*) './start.nc and ./startfi.nc'
        write(*,*)
        fichnom = 'start.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid_dyn)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Problem opening file:',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF
 
        fichnom = 'startfi.nc'
        ierr = NF_OPEN (fichnom, NF_NOWRITE,nid_fi)
        IF (ierr.NE.NF_NOERR) THEN
          write(6,*)' Problem opening file:',fichnom
          write(6,*)' ierr = ', ierr
          CALL ABORT
        ENDIF

        tab0 = 0 
        Lmodif = 0

      endif

c-----------------------------------------------------------------------
c Lecture du tableau des parametres du run (pour la dynamique)
c-----------------------------------------------------------------------

      if (choix_1.eq.0) then

        write(*,*) 'reading tab_cntrl START_ARCHIVE'
c
        ierr = NF_INQ_VARID (nid, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
        ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
c
      else if (choix_1.eq.1) then

        write(*,*) 'reading tab_cntrl START'
c
        ierr = NF_INQ_VARID (nid_dyn, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid_dyn, nvarid, tab_cntrl)
#else
        ierr = NF_GET_VAR_REAL(nid_dyn, nvarid, tab_cntrl)
#endif
c
        write(*,*) 'reading tab_cntrl STARTFI'
c
        ierr = NF_INQ_VARID (nid_fi, "controle", nvarid)
#ifdef NC_DOUBLE
        ierr = NF_GET_VAR_DOUBLE(nid_fi, nvarid, tab_cntrl_bis)
#else
        ierr = NF_GET_VAR_REAL(nid_fi, nvarid, tab_cntrl_bis)
#endif
c
        do i=1,50
          tab_cntrl(i+50)=tab_cntrl_bis(i)
        enddo
      write(*,*) 'printing tab_cntrl', tab_cntrl
      do i=1,100
        write(*,*) i,tab_cntrl(i)
      enddo
      
      endif
c-----------------------------------------------------------------------
c		Initialisation des constantes dynamique
c-----------------------------------------------------------------------

      kappa = tab_cntrl(9) 
      etot0 = tab_cntrl(12)
      ptot0 = tab_cntrl(13)
      ztot0 = tab_cntrl(14)
      stot0 = tab_cntrl(15)
      ang0 = tab_cntrl(16)
      write(*,*) "Newstart: kappa,etot0,ptot0,ztot0,stot0,ang0"
      write(*,*) kappa,etot0,ptot0,ztot0,stot0,ang0

c-----------------------------------------------------------------------
c   Lecture du tab_cntrl et initialisation des constantes physiques
c  - pour start:  Lmodif = 0 => pas de modifications possibles
c                  (modif dans le tabfi de readfi + loin)
c  - pour start_archive:  Lmodif = 1 => modifications possibles
c-----------------------------------------------------------------------
      if (choix_1.eq.0) then
         ! tabfi requires that input file be first opened by open_startphy(fichnom)
         fichnom = 'start_archive.nc'
         call open_startphy(fichnom)
         call tabfi (nid,Lmodif,tab0,day_ini,lllm,p_rad,
     .            p_omeg,p_g,p_mugaz,p_daysec,time)
      else if (choix_1.eq.1) then
         fichnom = 'startfi.nc'
         call open_startphy(fichnom)
         call tabfi (nid_fi,Lmodif,tab0,day_ini,lllm,p_rad,
     .            p_omeg,p_g,p_mugaz,p_daysec,time)
      endif

      rad = p_rad
      omeg = p_omeg
      g = p_g
      mugaz = p_mugaz
      daysec = p_daysec


c=======================================================================
c  INITIALISATIONS DIVERSES
c=======================================================================

      day_step=180 !?! Note: day_step is a common in "control.h"
      CALL defrun_new( 99, .TRUE. )
      CALL iniconst 
      CALL inigeom
      idum=-1
      idum=0

! Initialize the physics
         CALL iniphysiq(iim,jjm,llm,
     &                  (jjm-1)*iim+2,comm_lmdz,
     &                  daysec,day_ini,dtphys,
     &                  rlatu,rlatv,rlonu,rlonv,
     &                  aire,cu,cv,rad,g,r,cpp,
     &                  1)

c=======================================================================
c   lecture topographie, albedo, inertie thermique, relief sous-maille
c=======================================================================

      if (choix_1.ne.1) then  ! pour ne pas avoir besoin du fichier 
                              ! surface.dat dans le cas des start

c do while((relief(1:3).ne.'mol').AND.(relief(1:3).ne.'pla'))
c       write(*,*)
c       write(*,*) 'choix du relief (mola,pla)'
c       write(*,*) '(Topographie MGS MOLA, plat)'
c       read(*,fmt='(a3)') relief
        relief="mola"
c     enddo

      CALL datareadnc(relief,phis,alb,surfith,z0S,
     &          zmeaS,zstdS,zsigS,zgamS,ztheS,
     &          hmonsS,summitS,baseS,zavgS)

      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,phis,phisfi)
!      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ith,ithfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,surfith,surfithfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,z0S,z0)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zmeaS,zmea)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zstdS,zstd)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zsigS,zsig)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zgamS,zgam)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ztheS,zthe)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,hmonsS,hmons)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,summitS,summit)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,baseS,base)
!      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,zavgS,zavg)

      endif ! of if (choix_1.ne.1)


c=======================================================================
c  Lecture des fichiers (start ou start_archive)
c=======================================================================

      if (choix_1.eq.0) then

        write(*,*) 'Reading file START_ARCHIVE'
        CALL lect_start_archive(ngridmx,llm,nqtot,
     &   date,tsurf,tsoil,emis,q2,
     &   t,ucov,vcov,ps,co2ice,teta,phisold_newgrid,q,qsurf,
     &   tauscaling,totcloudfrac,surfith,nid,watercap)
        write(*,*) "OK, read start_archive file"
	! copy soil thermal inertia
	ithfi(:,:)=inertiedat(:,:)
	
        ierr= NF_CLOSE(nid)

      else if (choix_1.eq.1) then !  c'est l'appel a tabfi de phyeta0 qui
                                  !  permet de changer les valeurs du 
                                  !  tab_cntrl Lmodif=1
        tab0=0
        Lmodif=1 ! Lmodif set to 1 to allow modifications in phyeta0                           
        write(*,*) 'Reading file START'
        fichnom = 'start.nc'
        CALL dynetat0(fichnom,vcov,ucov,teta,q,masse,
     &       ps,phis,time)

        write(*,*) 'Reading file STARTFI'
        fichnom = 'startfi.nc'
        CALL phyetat0 (fichnom,tab0,Lmodif,nsoilmx,ngridmx,llm,nqtot,
     &        day_ini,time,tsurf,tsoil,albedo,emis,
     &        q2,qsurf,co2ice,tauscaling,totcloudfrac,
     &        wstar,mem_Mccn_co2,mem_Nccn_co2,mem_Mh2o_co2,watercap)
        
        ! copy albedo and soil thermal inertia
        do i=1,ngridmx
          albfi(i) = albedodat(i)
	  do j=1,nsoilmx
           ithfi(i,j) = inertiedat(i,j)
	  enddo
        ! build a surfithfi(:) using 1st layer of ithfi(:), which might
        ! be neede later on if reinitializing soil thermal inertia
          surfithfi(i)=ithfi(i,1)
        enddo

      else 
        CALL exit(1)
      endif

      dtvr   = daysec/REAL(day_step)
      dtphys   = dtvr * REAL(iphysiq)

c=======================================================================
c 
c=======================================================================
! If tracer names follow 'old' convention (q01, q02, ...) then
! rename them
      count=0
      do iq=1,nqtot
        txt=" "
        write(txt,'(a1,i2.2)') 'q',iq
        if (txt.eq.tname(iq)) then
          count=count+1
        endif
      enddo ! of do iq=1,nqtot
      
      ! initialize tracer names noms(:) and indexes (igcm_co2, igcm_h2o_vap, ...)
      call initracer(ngridmx,nqtot,qsurf)
      
      if (count.eq.nqtot) then
        write(*,*) 'Newstart: updating tracer names'
        ! copy noms(:) to tname(:) to have matching tracer names in physics
        ! and dynamics
        tname(1:nqtot)=noms(1:nqtot)
      endif

c=======================================================================
c 
c=======================================================================

      do ! infinite loop on list of changes

      write(*,*)
      write(*,*)
      write(*,*) 'List of possible changes :'
      write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(*,*)
      write(*,*) 'flat         : no topography ("aquaplanet")'
      write(*,*) 'bilball      : uniform albedo and thermal inertia'
      write(*,*) 'z0           : set a uniform surface roughness length'
      write(*,*) 'coldspole    : cold subsurface and high albedo at
     $ S.Pole'
      write(*,*) 'qname        : change tracer name'
      write(*,*) 'q=0          : ALL tracer =zero'
      write(*,*) 'q=factor     : change tracer value by a multiplicative
     & factor'
      write(*,*) 'q=x          : give a specific uniform value to one
     $ tracer'
      write(*,*) 'q=profile    : specify a profile for a tracer'
      write(*,*) 'freedust     : rescale dust to a true value'
      write(*,*) 'ini_q        : tracers initialization for chemistry
     $ and water vapour'
      write(*,*) 'ini_q-h2o    : tracers initialization for chemistry
     $ only'
      write(*,*) 'composition  : change atm main composition: CO2,N2,Ar,
     $ O2,CO'
      write(*,*) 'ini_h2osurf  : reinitialize surface water ice '
      write(*,*) 'noglacier    : Remove tropical H2O ice if |lat|<45'
      write(*,*) 'watercapn    : H20 ice on permanent N polar cap '
      write(*,*) 'watercaps    : H20 ice on permanent S polar cap '
      write(*,*) 'wetstart     : start with a wet atmosphere'
      write(*,*) 'isotherm     : Isothermal Temperatures, wind set to
     $ zero'
      write(*,*) 'co2ice=0     : remove CO2 polar cap'
      write(*,*) 'ptot         : change total pressure'
      write(*,*) 'therm_ini_s  : set soil thermal inertia to reference
     $ surface values'
      write(*,*) 'subsoilice_n : put deep underground ice layer in
     $ northern hemisphere'
      write(*,*) 'subsoilice_s : put deep underground ice layer in
     $ southern hemisphere'
      write(*,*) 'mons_ice     : put underground ice layer according
     $ to MONS derived data'

        write(*,*)
        write(*,*) 'Change to perform ?'
        write(*,*) '   (enter keyword or return to end)'
        write(*,*)

        read(*,fmt='(a20)') modif
        if (modif(1:1) .eq. ' ') exit ! exit loop on changes

        write(*,*)
        write(*,*) trim(modif) , ' : '

c       'flat : no topography ("aquaplanet")'
c       -------------------------------------
        if (trim(modif) .eq. 'flat') then
c         set topo to zero 
          phis(:,:)=0
          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,phis,phisfi)
          write(*,*) 'topography set to zero.'
          write(*,*) 'WARNING : the subgrid topography parameters',
     &    ' were not set to zero ! => set calllott to F'                    

c        Choice for surface pressure
         yes=' '
         do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*) 'Do you wish to choose homogeneous surface',
     &                 'pressure (y) or let newstart interpolate ',
     &                 ' the previous field  (n)?'
             read(*,fmt='(a)') yes
         end do
         if (yes.eq.'y') then
           flagps0=.true.
           write(*,*) 'New value for ps (Pa) ?'
 201       read(*,*,iostat=ierr) patm
            if(ierr.ne.0) goto 201
             write(*,*)
             write(*,*) ' new ps everywhere (Pa) = ', patm
             write(*,*)
             do j=1,jjp1
               do i=1,iip1
                 ps(i,j)=patm
               enddo
             enddo
         end if

c       bilball : albedo, inertie thermique uniforme
c       --------------------------------------------
        else if (trim(modif) .eq. 'bilball') then
          write(*,*) 'constante albedo and iner.therm:'
          write(*,*) 'New value for albedo (ex: 0.25) ?'
 101      read(*,*,iostat=ierr) alb_bb
          if(ierr.ne.0) goto 101
          write(*,*)
          write(*,*) ' uniform albedo (new value):',alb_bb
          write(*,*)

          write(*,*) 'New value for thermal inertia (eg: 247) ?'
 102      read(*,*,iostat=ierr) ith_bb
          if(ierr.ne.0) goto 102
          write(*,*) 'uniform thermal inertia (new value):',ith_bb
          DO j=1,jjp1
             DO i=1,iip1
                alb(i,j) = alb_bb	! albedo
		do isoil=1,nsoilmx
                  ith(i,j,isoil) = ith_bb	! thermal inertia
		enddo
             END DO
          END DO
!          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,ith,ithfi)
          CALL gr_dyn_fi(nsoilmx,iip1,jjp1,ngridmx,ith,ithfi)
          CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)
        
         ! also reset surface roughness length to default value
         write(*,*) 'surface roughness length set to:',z0_default,' m'
         z0(:)=z0_default

!       z0 : set surface roughness length to a constant value
!       -----------------------------------------------------
        else if (trim(modif) .eq. 'z0') then
          write(*,*) 'set a uniform surface roughness length'
          write(*,*) ' value for z0_default (ex: ',z0_default,')?'
          ierr=1
          do while (ierr.ne.0)
            read(*,*,iostat=ierr) z0_default
          enddo
          z0(:)=z0_default

c       coldspole : sous-sol de la calotte sud toujours froid
c       -----------------------------------------------------
        else if (trim(modif) .eq. 'coldspole') then
          write(*,*)'new value for the subsurface temperature',
     &   ' beneath the permanent southern polar cap ? (eg: 141 K)'
 103      read(*,*,iostat=ierr) tsud
          if(ierr.ne.0) goto 103
          write(*,*)
          write(*,*) ' new value of the subsurface temperature:',tsud
c         nouvelle temperature sous la calotte permanente
          do l=2,nsoilmx
               tsoil(ngridmx,l) =  tsud
          end do


          write(*,*)'new value for the albedo',
     &   'of the permanent southern polar cap ? (eg: 0.75)'
 104      read(*,*,iostat=ierr) albsud
          if(ierr.ne.0) goto 104
          write(*,*)

c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         Option 1:  only the albedo of the pole is modified :    
          albfi(ngridmx)=albsud
          write(*,*) 'ig=',ngridmx,'   albedo perennial cap ',
     &    albfi(ngridmx)

c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
c          Option 2 A haute resolution : coordonnee de la vrai calotte ~    
c           DO j=1,jjp1
c             DO i=1,iip1
c                ig=1+(j-2)*iim +i
c                if(j.eq.1) ig=1
c                if(j.eq.jjp1) ig=ngridmx
c                if ((rlatu(j)*180./pi.lt.-84.).and.
c     &            (rlatu(j)*180./pi.gt.-91.).and.
c     &            (rlonv(i)*180./pi.gt.-91.).and.
c     &            (rlonv(i)*180./pi.lt.0.))         then
cc    albedo de la calotte permanente fixe a albsud
c                   alb(i,j)=albsud
c                   write(*,*) 'lat=',rlatu(j)*180./pi,
c     &                      ' lon=',rlonv(i)*180./pi
cc     fin de la condition sur les limites de la calotte permanente
c                end if
c             ENDDO
c          ENDDO
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c         CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,alb,albfi)


c       ptot : Modification of the total pressure: ice + current atmosphere 
c       -------------------------------------------------------------------
        else if (trim(modif) .eq. 'ptot') then

c         calcul de la pression totale glace + atm actuelle
          patm=0.
          airetot=0.
          pcap=0.
          DO j=1,jjp1
             DO i=1,iim
                ig=1+(j-2)*iim +i
                if(j.eq.1) ig=1
                if(j.eq.jjp1) ig=ngridmx
                patm = patm + ps(i,j)*aire(i,j)
                airetot= airetot + aire(i,j)
                pcap = pcap + aire(i,j)*co2ice(ig)*g
             ENDDO
          ENDDO
          ptoto = pcap + patm

          print*,'Current total pressure at surface (co2 ice + atm) ',
     &       ptoto/airetot

          print*,'new value?'
          read(*,*) ptotn
          ptotn=ptotn*airetot
          patmn=ptotn-pcap
          print*,'ptoto,patm,ptotn,patmn'
          print*,ptoto,patm,ptotn,patmn
          print*,'Mult. factor for pressure (atm only)', patmn/patm
          do j=1,jjp1
             do i=1,iip1
                ps(i,j)=ps(i,j)*patmn/patm
             enddo
          enddo

c        Correction pour la conservation des traceurs
         yes=' '
         do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*) 'Do you wish to conserve tracer total mass (y)',
     &              ' or tracer mixing ratio (n) ?'
             read(*,fmt='(a)') yes
         end do

         if (yes.eq.'y') then
           write(*,*) 'OK : conservation of tracer total mass'
           DO iq =1, nqtot
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=q(i,j,l,iq)*patm/patmn
                  ENDDO
               ENDDO
             ENDDO
           ENDDO
          else
            write(*,*) 'OK : conservation of tracer mixing ratio'
          end if

c       qname : change tracer name
c       --------------------------
        else if (trim(modif).eq.'qname') then
         yes='y'
         do while (yes.eq.'y')
          write(*,*) 'Which tracer name do you want to change ?'
          do iq=1,nqtot
            write(*,'(i3,a3,a20)')iq,' : ',trim(tname(iq))
          enddo
          write(*,'(a35,i3)')
     &            '(enter tracer number; between 1 and ',nqtot
          write(*,*)' or any other value to quit this option)'
          read(*,*) iq
          if ((iq.ge.1).and.(iq.le.nqtot)) then
            write(*,*)'Change tracer name ',trim(tname(iq)),' to ?'
            read(*,*) txt
            tname(iq)=txt
            write(*,*)'Do you want to change another tracer name (y/n)?'
            read(*,'(a)') yes 
          else
! inapropiate value of iq; quit this option
            yes='n'
          endif ! of if ((iq.ge.1).and.(iq.le.nqtot))
         enddo ! of do while (yes.ne.'y')

c       q=0 : set tracers to zero
c       -------------------------
        else if (trim(modif) .eq. 'q=0') then
c          mise a 0 des q (traceurs)
          write(*,*) 'Tracers set to 0 (1.E-30 in fact)'
           DO iq =1, nqtot
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=1.e-30
                  ENDDO
               ENDDO
             ENDDO
           ENDDO

c          set surface tracers to zero
           DO iq =1, nqtot
             DO ig=1,ngridmx
                 qsurf(ig,iq)=0.
             ENDDO
           ENDDO

c       q=factor : change value of tracer by a multiplicative factor
c       ------------------------------------------------------------
        else if (trim(modif) .eq. 'q=factor') then
             write(*,*) 'Which tracer do you want to modify ?'
             do iq=1,nqtot
               write(*,*)iq,' : ',trim(tname(iq))
             enddo
             write(*,*) '(choose between 1 and ',nqtot,')'
             read(*,*) iq 
             if ((iq.lt.1).or.(iq.gt.nqtot)) then
               ! wrong value for iq, go back to menu
               write(*,*) "wrong input value:",iq
               cycle
             endif
             write(*,*)"factor to multiply current mixing ratio by?"
             read(*,*) val
             
             q(1:iip1,1:jjp1,1:llm,iq)=q(1:iip1,1:jjp1,1:llm,iq)*val
             qsurf(1:ngridmx,iq)=qsurf(1:ngridmx,iq)*val

c       q=x : initialise tracer manually 
c       --------------------------------
        else if (trim(modif) .eq. 'q=x') then
             write(*,*) 'Which tracer do you want to modify ?'
             do iq=1,nqtot
               write(*,*)iq,' : ',trim(tname(iq))
             enddo
             write(*,*) '(choose between 1 and ',nqtot,')'
             read(*,*) iq 
             if ((iq.lt.1).or.(iq.gt.nqtot)) then
               ! wrong value for iq, go back to menu
               write(*,*) "wrong input value:",iq
               cycle
             endif
             write(*,*)'mixing ratio of tracer ',trim(tname(iq)),
     &                 ' ? (kg/kg)'
             read(*,*) val
             DO l=1,llm
               DO j=1,jjp1
                  DO i=1,iip1
                    q(i,j,l,iq)=val
                  ENDDO
               ENDDO
             ENDDO
             write(*,*) 'SURFACE value of tracer ',trim(tname(iq)),
     &                   ' ? (kg/m2)'
             read(*,*) val
             DO ig=1,ngridmx
                 qsurf(ig,iq)=val
             ENDDO

c       q=profile : initialize tracer with a given profile
c       --------------------------------------------------
        else if (trim(modif) .eq. 'q=profile') then
             write(*,*) 'Tracer profile will be sought in ASCII file'
             write(*,*) "'profile_tracer' where 'tracer' is tracer name"
             write(*,*) "(one value per line in file; starting with"
             write(*,*) "surface value, the 1st atmospheric layer"
             write(*,*) "followed by 2nd, etc. up to top of atmosphere)"
             write(*,*) 'Which tracer do you want to set?'
             do iq=1,nqtot
               write(*,*)iq,' : ',trim(tname(iq))
             enddo
             write(*,*) '(choose between 1 and ',nqtot,')'
             read(*,*) iq 
             if ((iq.lt.1).or.(iq.gt.nqtot)) then
               ! wrong value for iq, go back to menu
               write(*,*) "wrong input value:",iq
               cycle
             endif
             ! look for input file 'profile_tracer'
             txt="profile_"//trim(tname(iq))
             open(41,file=trim(txt),status='old',form='formatted',
     &            iostat=ierr)
             if (ierr.eq.0) then
               ! OK, found file 'profile_...', load the profile
               do l=1,llm+1
                 read(41,*,iostat=ierr) profile(l)
                 if (ierr.ne.0) then ! something went wrong
                   exit ! quit loop
                 endif
               enddo
               if (ierr.eq.0) then
                 ! initialize tracer values
                 qsurf(:,iq)=profile(1)
                 do l=1,llm
                   q(:,:,l,iq)=profile(l+1)
                 enddo
                 write(*,*)'OK, tracer ',trim(tname(iq)),
     &               ' initialized ','using values from file ',trim(txt)
               else
                 write(*,*)'problem reading file ',trim(txt),' !'
                 write(*,*)'No modifications to tracer ',trim(tname(iq))
               endif
             else
               write(*,*)'Could not find file ',trim(txt),' !'
               write(*,*)'No modifications to tracer ',trim(tname(iq))
             endif
             
c       convert dust from virtual to true values
c       --------------------------------------------------
        else if (trim(modif) .eq. 'freedust') then
         if (minval(tauscaling) .lt. 0) then
           write(*,*) 'WARNING conversion factor negative'
           write(*,*) 'This is probably because it was not present
     &in the file'
           write(*,*) 'A constant conversion is used instead.'
           tauscaling(:) = 1.e-3
         endif
         CALL gr_fi_dyn(1,ngridmx,iip1,jjp1,tauscaling,tauscadyn)
          do l=1,llm
            do j=1,jjp1
              do i=1,iip1
                if (igcm_dust_number .ne. 0) 
     &            q(i,j,l,igcm_dust_number) =
     &            q(i,j,l,igcm_dust_number) * tauscadyn(i,j)
                if (igcm_dust_mass .ne. 0) 
     &            q(i,j,l,igcm_dust_mass) =
     &            q(i,j,l,igcm_dust_mass) * tauscadyn(i,j)
                if (igcm_ccn_number .ne. 0) 
     &            q(i,j,l,igcm_ccn_number) =
     &            q(i,j,l,igcm_ccn_number) * tauscadyn(i,j)
                if (igcm_ccn_mass .ne. 0) 
     &            q(i,j,l,igcm_ccn_mass) =
     &            q(i,j,l,igcm_ccn_mass) * tauscadyn(i,j)
              end do
            end do
          end do

          tauscaling(:) = 1.

         ! We want to have the very same value at lon -180 and lon 180
          do l = 1,llm
             do j = 1,jjp1
                do iq = 1,nqtot
                   q(iip1,j,l,iq) = q(1,j,l,iq)
                end do
             end do
          end do

          write(*,*) 'done rescaling to true vale'

c       ini_q : Initialize tracers for chemistry
c       -----------------------------------------------
        else if (trim(modif) .eq. 'ini_q') then
          flagh2o    = 1
          flagthermo = 0
          yes=' '
c         For more than 32 layers, possible to initiate thermosphere only     
          if (llm.gt.32) then 
            do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*)'',
     &     'initialisation for thermosphere only? (y/n)'
            read(*,fmt='(a)') yes
            if (yes.eq.'y') then
            flagthermo=1 
            else
            flagthermo=0
            endif
            enddo  
          endif
          
          call inichim_newstart(ngridmx, nqtot, q, qsurf, ps, 
     &                          flagh2o, flagthermo)

         ! We want to have the very same value at lon -180 and lon 180
          do l = 1,llm
             do j = 1,jjp1
                do iq = 1,nqtot
                   q(iip1,j,l,iq) = q(1,j,l,iq)
                end do
             end do
          end do

          write(*,*) 'inichim_newstart: chemical species and
     $ water vapour initialised'

c       ini_q-h2o : as above except for the water vapour tracer 
c       ------------------------------------------------------
        else if (trim(modif) .eq. 'ini_q-h2o') then
          flagh2o    = 0
          flagthermo = 0
          yes=' '
          ! for more than 32 layers, possible to initiate thermosphere only     
          if(llm.gt.32) then
            do while ((yes.ne.'y').and.(yes.ne.'n'))
            write(*,*)'',
     &      'initialisation for thermosphere only? (y/n)'
            read(*,fmt='(a)') yes
            if (yes.eq.'y') then 
            flagthermo=1 
            else
            flagthermo=0
            endif
            enddo
          endif

          call inichim_newstart(ngridmx, nqtot, q, qsurf, ps, 
     &                          flagh2o, flagthermo)

         ! We want to have the very same value at lon -180 and lon 180
          do l = 1,llm
             do j = 1,jjp1
                do iq = 1,nqtot
                   q(iip1,j,l,iq) = q(1,j,l,iq)
                end do
             end do
          end do

          write(*,*) 'inichim_newstart: chemical species initialised
     $ (except water vapour)'

c      composition : change main composition: CO2,N2,Ar,O2,CO (FF 03/2014)
c      --------------------------------------------------------
       else if (trim(modif) .eq. 'composition') then
          write(*,*) "Lat (degN)  lon (degE) of the reference site ?"
          write(*,*) "e.g. MSL : -4.5  137.  "
 301      read(*,*,iostat=ierr) latref, lonref
          if(ierr.ne.0) goto 301


        !  Select GCM point close to reference site
          dlonmin =90.
          DO i=1,iip1-1
             if (abs(rlonv(i)*180./pi -lonref).lt.dlonmin)then
                iref=i
                dlonmin=abs(rlonv(i)*180./pi -lonref)
             end if   
          ENDDO
          dlatmin =45.
          DO j=1,jjp1
             if (abs(rlatu(j)*180./pi -latref).lt.dlatmin)then
                jref=j
                dlatmin=abs(rlatu(j)*180./pi -latref)
             end if   
          ENDDO
          write(*,*) "In GCM : lat= " ,  rlatu(jref)*180./pi
          write(*,*) "In GCM : lon= " ,  rlonv(iref)*180./pi
          write(*,*)

        ! Compute air molar mass at reference site
          Smmr=0.
          Sn = 0.
          write(*,*) 'igcm_co2 = ', igcm_co2
          write(*,*) 'igcm_n2 = ', igcm_n2
          write(*,*) 'igcm_ar = ', igcm_ar
          write(*,*) 'igcm_o2 = ', igcm_o2
          write(*,*) 'igcm_co = ', igcm_co
          write(*,*) noms
          do iq=1,nqtot 
             if ((iq.eq.igcm_co2).or.(iq.eq.igcm_n2)
     &      .or. (iq.eq.igcm_ar).or.(iq.eq.igcm_o2)
     &      .or. (iq.eq.igcm_co)) then
                 Smmr=Smmr+q(iref,jref,1,iq)
                 Sn=Sn+q(iref,jref,1,iq)/mmol(iq) 
             end if
          end do
        ! Special case : implicit non-co2 gases ! JN 11/2019
          if ((igcm_n2.eq.0) .or. (igcm_ar.eq.0)) then 
           write(*,*) "Warning : non-co2 gases are implicit :  "
           write(*,*) "At reference site :  "
       !    write(*,*) "q= ", q(iref, jref, 1,igcm_co2)
           write(*,*) "Sum of mass mix. ratio (ie MMR(co2))=",Smmr
           Mair_old = 44.0*Smmr  + 33.87226017157708*(1-Smmr) 
      
       !  33.87226017157708 is the 
       !   molar mass of non-co2 atmosphere measured by MSL at Ls ~184
       
          else
            ! Assume co2/n2/ar/o2/co are available
            Mair_old=(q(iref,jref,1,igcm_co2)*mmol(igcm_co2)
     &               +q(iref,jref,1,igcm_n2)*mmol(igcm_n2)
     &               +q(iref,jref,1,igcm_ar)*mmol(igcm_ar)
     &               +q(iref,jref,1,igcm_o2)*mmol(igcm_o2)
     &               +q(iref,jref,1,igcm_co)*mmol(igcm_co))/Smmr
          end if

          write(*,*)
     &      "Air molar mass (g/mol) at reference site= ",Mair_old

        ! Ask for new volume mixing ratio at reference site
          Svmr =0.
          Sn =0.
          coefvmr(igcm_co2)=1.

          do iq=1,nqtot 
           coefvmr(iq) = 1.
           if ((iq.eq.igcm_n2).or.(iq.eq.igcm_ar)
     &     .or. (iq.eq.igcm_o2).or.(iq.eq.igcm_co)) then

             vmr_old=q(iref,jref,1,iq)*Mair_old/mmol(iq)  
             write(*,*) "Previous vmr("//trim(tname(iq))//")= ", vmr_old

              if (iq.eq.igcm_n2) then
                write(*,*) "New vmr(n2)? (MSL: 2.8e-02 at Ls~180,",
     &           " Trainer et al. 2019)"
              endif
              if (iq.eq.igcm_ar) then
                write(*,*) "New vmr(ar)? (MSL: 2.1e-02 at Ls~180)"
              endif
              if (iq.eq.igcm_o2) then
                write(*,*) "New vmr(o2)? (MSL: 1.7e-03 at Ls~180)"
              endif
              if (iq.eq.igcm_co) then
                write(*,*) "New vmr(co)? (ACS: 1e-03 at Ls~180)"
              endif
 302          read(*,*,iostat=ierr) vmr_new
              if(ierr.ne.0) goto 302
              write(*,*) "New vmr("//trim(tname(iq))//")= ",vmr_new
              write(*,*) 
              coefvmr(iq) = vmr_new/vmr_old
              Svmr=Svmr+vmr_new
              Sn=Sn+vmr_new*mmol(iq)
           end if
          enddo ! of do iq=1,nqtot 

        ! Special case : implicit non-co2 gases JN 11/2019
          if ((igcm_n2.eq.0) .or. (igcm_ar.eq.0)) then 
            write(*,*) "Warning : non-co2 gases are implicit"
            vmr_old=q(iref,jref,1,igcm_co2)*Mair_old/mmol(igcm_co2)  
            write(*,*) "Previous vmr(co2)=", vmr_old
            write(*,*) "New vmr(co2) ? (MSL: 0.947 at Ls~180)",
     &                  " Trainer et al. 2019)"
 666          read(*,*,iostat=ierr) vmr_new
              if(ierr.ne.0) goto 666
              coefvmr(igcm_co2) = vmr_new/vmr_old
              Svmr=Svmr+vmr_new
              Sn=vmr_new*mmol(igcm_co2) + (1-vmr_new)
     &         *33.87226017157708 ! Molar mass of non-co2 atm from MSL
          end if
      !  Estimation of new Air molar mass at reference site (assuming vmr_co2 = 1-Svmr)
          Mair_new = Sn + (1-Svmr)*mmol(igcm_co2) 
      !  Estimation of new Air molar mass when non-co2 gases are implicit
          if ((igcm_n2.eq.0) .or. (igcm_ar.eq.0)) then 
              Mair_new=vmr_new*mmol(igcm_co2) + (1-vmr_new)
     &         *33.87226017157708 ! Molar mass of non-co2 atm from MSL
           write(*,*)
     &     "We consider non-co2 gases vmr measured from Curiosity"
          end if
          write(*,*)
     &     "NEW Air molar mass (g/mol) at reference site= ",Mair_new

        ! Compute mass mixing ratio changes  
          do iq=1,nqtot  
            if ((iq.eq.igcm_co2).or.(iq.eq.igcm_n2).or.(iq.eq.igcm_ar)
     &          .or. (iq.eq.igcm_o2).or.(iq.eq.igcm_co)) then
             write(*,*) "Everywhere mmr("//trim(tname(iq))//
     &        ") is multiplied by ",coefvmr(iq)*Mair_old/Mair_new
            end if
          end do

        ! Recompute mass mixing ratios everywhere, and adjust mmr of most abundant species
        ! to keep sum of mmr constant.
          do l=1,llm
           do j=1,jjp1
            do i=1,iip1
              Smmr_old = 0.
              Smmr_new = 0.
              do iq=1,nqtot  
                if ((iq.eq.igcm_co2).or.(iq.eq.igcm_n2)
     &          .or.(iq.eq.igcm_ar)
     &          .or. (iq.eq.igcm_o2).or.(iq.eq.igcm_co)
     &          .or. (iq.eq.igcm_o) .or. (iq.eq. igcm_h2) ) then
                   Smmr_old = Smmr_old + q(i,j,l,iq) ! sum of old mmr 
                   q(i,j,l,iq)=q(i,j,l,iq)*coefvmr(iq)*Mair_old/Mair_new
                   Smmr_new = Smmr_new + q(i,j,l,iq) ! sum of new mmr
                end if 
              enddo
              !iloc = maxloc(q(i,j,l,:))
              iqmax=0 ; maxq=0
              do iq=1,nqtot
                if ((iq.eq.igcm_co2).or.(iq.eq.igcm_n2)
     &          .or.(iq.eq.igcm_ar)
     &          .or. (iq.eq.igcm_o2).or.(iq.eq.igcm_co)
     &          .or. (iq.eq.igcm_o) .or. (iq.eq. igcm_h2) ) then
                  if (q(i,j,l,iq).gt.maxq) then
                    maxq=q(i,j,l,iq)
                    iqmax=iq
                  endif
                endif
              enddo
              !iqmax = iloc(1)
              q(i,j,l,iqmax) = q(i,j,l,iqmax) + Smmr_old - Smmr_new
            enddo
           enddo
          enddo

          write(*,*)
     &   "The most abundant species is modified everywhere to keep "//
     &   "sum of mmr constant"
          write(*,*) 'At reference site vmr(CO2)=', 
     &        q(iref,jref,1,igcm_co2)*Mair_new/mmol(igcm_co2)
          write(*,*) "Compared to MSL observation: vmr(CO2)= 0.947 "//
     &   "at Ls=180" 

          Sn = q(iref,jref,1,igcm_co2)*Mair_new/mmol(igcm_co2)
     &       + q(iref,jref,1,igcm_n2)*Mair_new/mmol(igcm_n2)
     &       + q(iref,jref,1,igcm_ar)*Mair_new/mmol(igcm_ar)
     &       + q(iref,jref,1,igcm_o2)*Mair_new/mmol(igcm_o2)
     &       + q(iref,jref,1,igcm_co)*Mair_new/mmol(igcm_co)

          write(*,*) 'Sum of volume mixing ratios = ', Sn

c      wetstart : wet atmosphere with a north to south gradient
c      --------------------------------------------------------
       else if (trim(modif) .eq. 'wetstart') then
        ! check that there is indeed a water vapor tracer
        if (igcm_h2o_vap.eq.0) then
          write(*,*) "No water vapour tracer! Can't use this option"
          stop
        endif
          DO l=1,llm
            DO j=1,jjp1
              DO i=1,iip1-1
                q(i,j,l,igcm_h2o_vap)=150.e-6 * (rlatu(j)+pi/2.) / pi
              ENDDO
              ! We want to have the very same value at lon -180 and lon 180
              q(iip1,j,l,igcm_h2o_vap) = q(1,j,l,igcm_h2o_vap)
            ENDDO
          ENDDO

         write(*,*) 'Water mass mixing ratio at north pole='
     *               ,q(1,1,1,igcm_h2o_vap)
         write(*,*) '---------------------------south pole='
     *               ,q(1,jjp1,1,igcm_h2o_vap)

c      ini_h2osurf : reinitialize surface water ice
c      --------------------------------------------------
        else if (trim(modif) .eq. 'ini_h2osurf') then
          write(*,*)'max surface ice left?(e.g. 0.2 kg/m2=200microns)'
 207      read(*,*,iostat=ierr) val
          if(ierr.ne.0) goto 207
          write(*,*)'also set negative values of surf ice to 0'
           do ig=1,ngridmx
              qsurf(ig,igcm_h2o_ice)=min(val,qsurf(ig,igcm_h2o_ice))
              qsurf(ig,igcm_h2o_ice)=max(0.,qsurf(ig,igcm_h2o_ice))
           end do

c      noglacier : remove tropical water ice (to initialize high res sim)
c      --------------------------------------------------
        else if (trim(modif) .eq. 'noglacier') then
           do ig=1,ngridmx
             j=(ig-2)/iim +2
              if(ig.eq.1) j=1
              write(*,*) 'OK: remove surface ice for |lat|<45'
              if (abs(rlatu(j)*180./pi).lt.45.) then
                   qsurf(ig,igcm_h2o_ice)=0.
              end if
           end do


c      watercapn : H20 ice on permanent northern cap
c      --------------------------------------------------
        else if (trim(modif) .eq. 'watercapn') then
           do ig=1,ngridmx
             j=(ig-2)/iim +2
              if(ig.eq.1) j=1
              if (rlatu(j)*180./pi.gt.80.) then
                   qsurf(ig,igcm_h2o_ice)=1.e5
                   write(*,*) 'ig=',ig,'    H2O ice mass (kg/m2)= ',
     &              qsurf(ig,igcm_h2o_ice)
                   write(*,*)'     ==> Ice mesh South boundary (deg)= ',
     &              rlatv(j)*180./pi
                end if
           enddo

c      watercaps : H20 ice on permanent southern cap
c      -------------------------------------------------
        else if (trim(modif) .eq. 'watercaps') then
           do ig=1,ngridmx
               j=(ig-2)/iim +2
               if(ig.eq.1) j=1
               if (rlatu(j)*180./pi.lt.-80.) then
                   qsurf(ig,igcm_h2o_ice)=1.e5
                   write(*,*) 'ig=',ig,'   H2O ice mass (kg/m2)= ',
     &              qsurf(ig,igcm_h2o_ice)
                   write(*,*)'     ==> Ice mesh North boundary (deg)= ',
     &              rlatv(j-1)*180./pi
               end if
           enddo

c       isotherm : Isothermal temperatures and no winds
c       ------------------------------------------------
        else if (trim(modif) .eq. 'isotherm') then

          write(*,*)'Isothermal temperature of the atmosphere, 
     &           surface and subsurface'
          write(*,*) 'Value of this temperature ? :'
 203      read(*,*,iostat=ierr) Tiso
          if(ierr.ne.0) goto 203

          do ig=1, ngridmx
            tsurf(ig) = Tiso
          end do
          do l=2,nsoilmx
            do ig=1, ngridmx
              tsoil(ig,l) = Tiso
            end do
          end do
          flagiso=.true.
          call initial0(llm*ip1jmp1,ucov)
          call initial0(llm*ip1jm,vcov)
          call initial0(ngridmx*(llm+1),q2)

c       co2ice=0 : remove CO2 polar ice caps'
c       ------------------------------------------------
        else if (trim(modif) .eq. 'co2ice=0') then
           do ig=1,ngridmx
              co2ice(ig)=0
              emis(ig)=emis(ngridmx/2)
           end do
        
!       therm_ini_s: (re)-set soil thermal inertia to reference surface values
!       ----------------------------------------------------------------------

	else if (trim(modif).eq.'therm_ini_s') then
!          write(*,*)"surfithfi(1):",surfithfi(1)
	  do isoil=1,nsoilmx
	    inertiedat(1:ngridmx,isoil)=surfithfi(1:ngridmx)
	  enddo
          write(*,*)'OK: Soil thermal inertia has been reset to referenc
     &e surface values'
!	  write(*,*)"inertiedat(1,1):",inertiedat(1,1)
	  ithfi(:,:)=inertiedat(:,:)
	 ! recast ithfi() onto ith()
	 call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,ithfi,ith)
! Check:
!         do i=1,iip1
!           do j=1,jjp1
!             do isoil=1,nsoilmx
!               write(77,*) i,j,isoil,"  ",ith(i,j,isoil)
!             enddo
!           enddo
!	 enddo

!       subsoilice_n: Put deep ice layer in northern hemisphere soil
!       ------------------------------------------------------------

	else if (trim(modif).eq.'subsoilice_n') then

         write(*,*)'From which latitude (in deg.), up to the north pole,
     &should we put subterranean ice?'
	 ierr=1
	 do while (ierr.ne.0)
	  read(*,*,iostat=ierr) val
	  if (ierr.eq.0) then ! got a value
	    ! do a sanity check
	    if((val.lt.0.).or.(val.gt.90)) then
	      write(*,*)'Latitude should be between 0 and 90 deg. !!!'
	      ierr=1
	    else ! find corresponding jref (nearest latitude)
	      ! note: rlatu(:) contains decreasing values of latitude
	      !       starting from PI/2 to -PI/2
	      do j=1,jjp1
	        if ((rlatu(j)*180./pi.ge.val).and.
     &              (rlatu(j+1)*180./pi.le.val)) then
		  ! find which grid point is nearest to val:
		  if (abs(rlatu(j)*180./pi-val).le.
     &                abs((rlatu(j+1)*180./pi-val))) then
		   jref=j
		  else
		   jref=j+1
		  endif
		 
		 write(*,*)'Will use nearest grid latitude which is:',
     &                     rlatu(jref)*180./pi
		endif
	      enddo ! of do j=1,jjp1
	    endif ! of if((val.lt.0.).or.(val.gt.90))
	  endif !of if (ierr.eq.0)
	 enddo ! of do while

         ! Build layers() (as in soil_settings.F)
	 val2=sqrt(mlayer(0)*mlayer(1))
	 val3=mlayer(1)/mlayer(0)
	 do isoil=1,nsoilmx
	   layer(isoil)=val2*(val3**(isoil-1))
	 enddo

         write(*,*)'At which depth (in m.) does the ice layer begin?'
         write(*,*)'(currently, the deepest soil layer extends down to:'
     &              ,layer(nsoilmx),')'
	 ierr=1
	 do while (ierr.ne.0)
	  read(*,*,iostat=ierr) val2
!	  write(*,*)'val2:',val2,'ierr=',ierr
	  if (ierr.eq.0) then ! got a value, but do a sanity check
	    if(val2.gt.layer(nsoilmx)) then
	      write(*,*)'Depth should be less than ',layer(nsoilmx)
	      ierr=1
	    endif
	    if(val2.lt.layer(1)) then
	      write(*,*)'Depth should be more than ',layer(1)
	      ierr=1
	    endif
	  endif
         enddo ! of do while
	 
	 ! find the reference index iref the depth corresponds to
!	 if (val2.lt.layer(1)) then
!	  iref=1
!	 else
	  do isoil=1,nsoilmx-1
	   if((val2.gt.layer(isoil)).and.(val2.lt.layer(isoil+1)))
     &       then
	     iref=isoil
	     exit
	   endif
	  enddo
!	 endif
	 
!	 write(*,*)'iref:',iref,'  jref:',jref
!	 write(*,*)'layer',layer
!	 write(*,*)'mlayer',mlayer
         
	 ! thermal inertia of the ice:
	 ierr=1
	 do while (ierr.ne.0)
          write(*,*)'What is the value of subterranean ice thermal inert
     &ia? (e.g.: 2000)'
          read(*,*,iostat=ierr)iceith
	 enddo ! of do while
	 
	 ! recast ithfi() onto ith()
	 call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,ithfi,ith)
	 
	 do j=1,jref
!	    write(*,*)'j:',j,'rlatu(j)*180./pi:',rlatu(j)*180./pi
	    do i=1,iip1 ! loop on longitudes
	     ! Build "equivalent" thermal inertia for the mixed layer
	     ith(i,j,iref+1)=sqrt((layer(iref+1)-layer(iref))/
     &                     (((val2-layer(iref))/(ith(i,j,iref)**2))+
     &                      ((layer(iref+1)-val2)/(iceith)**2)))
	     ! Set thermal inertia of lower layers
	     do isoil=iref+2,nsoilmx
	      ith(i,j,isoil)=iceith ! ice
	     enddo
	    enddo ! of do i=1,iip1 
	 enddo ! of do j=1,jjp1
	 

	 CALL gr_dyn_fi(nsoilmx,iip1,jjp1,ngridmx,ith,ithfi)

!         do i=1,nsoilmx
!	  write(*,*)'i:',i,'ithfi(1,i):',ithfi(1,i)
!	 enddo

	
!       subsoilice_s: Put deep ice layer in southern hemisphere soil
!       ------------------------------------------------------------

	else if (trim(modif).eq.'subsoilice_s') then

         write(*,*)'From which latitude (in deg.), down to the south pol
     &e, should we put subterranean ice?'
	 ierr=1
	 do while (ierr.ne.0)
	  read(*,*,iostat=ierr) val
	  if (ierr.eq.0) then ! got a value
	    ! do a sanity check
	    if((val.gt.0.).or.(val.lt.-90)) then
	      write(*,*)'Latitude should be between 0 and -90 deg. !!!'
	      ierr=1
	    else ! find corresponding jref (nearest latitude)
	      ! note: rlatu(:) contains decreasing values of latitude
	      !       starting from PI/2 to -PI/2
	      do j=1,jjp1
	        if ((rlatu(j)*180./pi.ge.val).and.
     &              (rlatu(j+1)*180./pi.le.val)) then
		  ! find which grid point is nearest to val:
		  if (abs(rlatu(j)*180./pi-val).le.
     &                abs((rlatu(j+1)*180./pi-val))) then
		   jref=j
		  else
		   jref=j+1
		  endif
		 
		 write(*,*)'Will use nearest grid latitude which is:',
     &                     rlatu(jref)*180./pi
		endif
	      enddo ! of do j=1,jjp1
	    endif ! of if((val.lt.0.).or.(val.gt.90))
	  endif !of if (ierr.eq.0)
	 enddo ! of do while

         ! Build layers() (as in soil_settings.F)
	 val2=sqrt(mlayer(0)*mlayer(1))
	 val3=mlayer(1)/mlayer(0)
	 do isoil=1,nsoilmx
	   layer(isoil)=val2*(val3**(isoil-1))
	 enddo

         write(*,*)'At which depth (in m.) does the ice layer begin?'
         write(*,*)'(currently, the deepest soil layer extends down to:'
     &              ,layer(nsoilmx),')'
	 ierr=1
	 do while (ierr.ne.0)
	  read(*,*,iostat=ierr) val2
!	  write(*,*)'val2:',val2,'ierr=',ierr
	  if (ierr.eq.0) then ! got a value, but do a sanity check
	    if(val2.gt.layer(nsoilmx)) then
	      write(*,*)'Depth should be less than ',layer(nsoilmx)
	      ierr=1
	    endif
	    if(val2.lt.layer(1)) then
	      write(*,*)'Depth should be more than ',layer(1)
	      ierr=1
	    endif
	  endif
         enddo ! of do while
	 
	 ! find the reference index iref the depth corresponds to
	  do isoil=1,nsoilmx-1
	   if((val2.gt.layer(isoil)).and.(val2.lt.layer(isoil+1)))
     &       then
	     iref=isoil
	     exit
	   endif
	  enddo
	 
!	 write(*,*)'iref:',iref,'  jref:',jref
         
	 ! thermal inertia of the ice:
	 ierr=1
	 do while (ierr.ne.0)
          write(*,*)'What is the value of subterranean ice thermal inert
     &ia? (e.g.: 2000)'
          read(*,*,iostat=ierr)iceith
	 enddo ! of do while
	 
	 ! recast ithfi() onto ith()
	 call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,ithfi,ith)
	 
	 do j=jref,jjp1
!	    write(*,*)'j:',j,'rlatu(j)*180./pi:',rlatu(j)*180./pi
	    do i=1,iip1 ! loop on longitudes
	     ! Build "equivalent" thermal inertia for the mixed layer
	     ith(i,j,iref+1)=sqrt((layer(iref+1)-layer(iref))/
     &                     (((val2-layer(iref))/(ith(i,j,iref)**2))+
     &                      ((layer(iref+1)-val2)/(iceith)**2)))
	     ! Set thermal inertia of lower layers
	     do isoil=iref+2,nsoilmx
	      ith(i,j,isoil)=iceith ! ice
	     enddo
	    enddo ! of do i=1,iip1 
	 enddo ! of do j=jref,jjp1
	 

	 CALL gr_dyn_fi(nsoilmx,iip1,jjp1,ngridmx,ith,ithfi)

c       'mons_ice' : use MONS data to build subsurface ice table
c       --------------------------------------------------------
        else if (trim(modif).eq.'mons_ice') then
        
       ! 1. Load MONS data
        call load_MONS_data(MONS_Hdn,MONS_d21)
        
        ! 2. Get parameters from user
        ierr=1
	do while (ierr.ne.0)
          write(*,*) "Coefficient to apply to MONS 'depth' in Northern",
     &               " Hemisphere?"
          write(*,*) " (should be somewhere between 3.2e-4 and 1.3e-3)"
          read(*,*,iostat=ierr) MONS_coeffN
        enddo
        ierr=1
	do while (ierr.ne.0)
          write(*,*) "Coefficient to apply to MONS 'depth' in Southern",
     &               " Hemisphere?"
          write(*,*) " (should be somewhere between 3.2e-4 and 1.3e-3)"
          read(*,*,iostat=ierr) MONS_coeffS
        enddo
        ierr=1
        do while (ierr.ne.0)
          write(*,*) "Value of subterranean ice thermal inertia ",
     &               " in Northern hemisphere?"
          write(*,*) " (e.g.: 2000, or perhaps 2290)"
!          read(*,*,iostat=ierr) iceith
          read(*,*,iostat=ierr) iceithN
        enddo
        ierr=1
        do while (ierr.ne.0)
          write(*,*) "Value of subterranean ice thermal inertia ",
     &               " in Southern hemisphere?"
          write(*,*) " (e.g.: 2000, or perhaps 2290)"
!          read(*,*,iostat=ierr) iceith
          read(*,*,iostat=ierr) iceithS
        enddo
        
        ! 3. Build subterranean thermal inertia
        
        ! initialise subsurface inertia with reference surface values
        do isoil=1,nsoilmx
          ithfi(1:ngridmx,isoil)=surfithfi(1:ngridmx)
        enddo
        ! recast ithfi() onto ith()
	call gr_fi_dyn(nsoilmx,ngridmx,iip1,jjp1,ithfi,ith)
        
        do i=1,iip1 ! loop on longitudes
          do j=1,jjp1 ! loop on latitudes
            ! set MONS_coeff
            if (rlatu(j).ge.0) then ! northern hemisphere
              ! N.B: rlatu(:) contains decreasing values of latitude
	      !       starting from PI/2 to -PI/2
              MONS_coeff=MONS_coeffN
              iceith=iceithN
            else ! southern hemisphere
              MONS_coeff=MONS_coeffS
              iceith=iceithS
            endif
            ! check if we should put subterranean ice
            if (MONS_Hdn(i,j).ge.14.0) then ! no ice if Hdn<14%
              ! compute depth at which ice lies:
              val=MONS_d21(i,j)*MONS_coeff
              ! compute val2= the diurnal skin depth of surface inertia
              ! assuming a volumetric heat cap. of C=1.e6 J.m-3.K-1
              val2=ith(i,j,1)*1.e-6*sqrt(88775./3.14159)
              if (val.lt.val2) then
                ! ice must be below the (surface inertia) diurnal skin depth
                val=val2
              endif
              if (val.lt.layer(nsoilmx)) then ! subterranean ice
                ! find the reference index iref that depth corresponds to
                iref=0
                do isoil=1,nsoilmx-1
                 if ((val.ge.layer(isoil)).and.(val.lt.layer(isoil+1)))
     &             then
	           iref=isoil
	           exit
	         endif
                enddo
                ! Build "equivalent" thermal inertia for the mixed layer
                ith(i,j,iref+1)=sqrt((layer(iref+1)-layer(iref))/
     &                     (((val-layer(iref))/(ith(i,j,iref+1)**2))+
     &                      ((layer(iref+1)-val)/(iceith)**2)))
	        ! Set thermal inertia of lower layers
                do isoil=iref+2,nsoilmx
                  ith(i,j,isoil)=iceith 
                enddo
              endif ! of if (val.lt.layer(nsoilmx))
            endif ! of if (MONS_Hdn(i,j).lt.14.0)
          enddo ! do j=1,jjp1
        enddo ! do i=1,iip1
        
! Check:
!         do i=1,iip1
!           do j=1,jjp1
!             do isoil=1,nsoilmx
!               write(77,*) i,j,isoil,"  ",ith(i,j,isoil)
!             enddo
!           enddo
!	 enddo

        ! recast ith() into ithfi()
        CALL gr_dyn_fi(nsoilmx,iip1,jjp1,ngridmx,ith,ithfi)
        
	else
          write(*,*) '       Unknown (misspelled?) option!!!'
        end if ! of if (trim(modif) .eq. '...') elseif ...
	
       enddo ! of do ! infinite loop on liste of changes

 999  continue

 
c=======================================================================
c   Correct pressure on the new grid (menu 0)
c=======================================================================

      if (choix_1.eq.0) then
        r = 1000.*8.31/mugaz

        do j=1,jjp1
          do i=1,iip1
             ps(i,j) = ps(i,j) * 
     .            exp((phisold_newgrid(i,j)-phis(i,j)) /
     .                                  (t(i,j,1) * r))
          end do
        end do
  
c periodicity of surface ps in longitude
        do j=1,jjp1
          ps(1,j) = ps(iip1,j)
        end do
      end if

c=======================================================================
c=======================================================================

c=======================================================================
c    Initialisation de la physique / ecriture de newstartfi :
c=======================================================================


      CALL inifilr 
      CALL pression(ip1jmp1, ap, bp, ps, p3d)

c-----------------------------------------------------------------------
c   Initialisation  pks:
c-----------------------------------------------------------------------

      CALL exner_hyb(ip1jmp1, ps, p3d, pks, pk, pkf)
! Calcul de la temperature potentielle teta

      if (flagiso) then
          DO l=1,llm
             DO j=1,jjp1
                DO i=1,iim
                   teta(i,j,l) = Tiso * cpp/pk(i,j,l)
                ENDDO
                teta (iip1,j,l)= teta (1,j,l)
             ENDDO
          ENDDO
      else if (choix_1.eq.0) then
         DO l=1,llm
            DO j=1,jjp1
               DO i=1,iim
                  teta(i,j,l) = t(i,j,l) * cpp/pk(i,j,l)
               ENDDO
               teta (iip1,j,l)= teta (1,j,l)
            ENDDO
         ENDDO
      endif

C Calcul intermediaire
c
      if (choix_1.eq.0) then
         CALL massdair( p3d, masse  )
c
         print *,' ALPHAX ',alphax

         DO  l = 1, llm
           DO  i    = 1, iim
             xppn(i) = aire( i, 1   ) * masse(  i     ,  1   , l )
             xpps(i) = aire( i,jjp1 ) * masse(  i     , jjp1 , l )
           ENDDO
             xpn      = SUM(xppn)/apoln
             xps      = SUM(xpps)/apols
           DO i   = 1, iip1
             masse(   i   ,   1     ,  l )   = xpn
             masse(   i   ,   jjp1  ,  l )   = xps
           ENDDO
         ENDDO
      endif
      phis(iip1,:) = phis(1,:)

c      CALL inidissip ( lstardis, nitergdiv, nitergrot, niterh,
c     *                tetagdiv, tetagrot , tetatemp  )
      itau=0
      if (choix_1.eq.0) then
         day_ini=int(date)
         hour_ini=date-int(date)
      endif
c
      CALL geopot  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )

      CALL caldyn0( itau,ucov,vcov,teta,ps,masse,pk,phis ,
     *                phi,w, pbaru,pbarv,day_ini+time )
c     CALL caldyn
c    $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
c    $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, day_ini )

      CALL dynredem0("restart.nc",day_ini,phis)
      CALL dynredem1("restart.nc",hour_ini,vcov,ucov,teta,q,
     .               masse,ps)
C
C Ecriture etat initial physique
C
      call physdem0("restartfi.nc",longitude,latitude,
     &              nsoilmx,ngridmx,llm,
     &              nqtot,dtphys,real(day_ini),0.0,cell_area,
     &              albfi,ithfi,zmea,zstd,zsig,zgam,zthe,
     &              hmons,summit,base)
      call physdem1("restartfi.nc",nsoilmx,ngridmx,llm,nqtot,
     &              dtphys,hour_ini,
     &              tsurf,tsoil,co2ice,albedo,emis,q2,qsurf,tauscaling,
     &              totcloudfrac,wstar,
     &              mem_Mccn_co2,mem_Nccn_co2,mem_Mh2o_co2,watercap)

c=======================================================================
c	Formats 
c=======================================================================

   1  FORMAT(//10x,'la valeur de im =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree iim =',i4//)
   2  FORMAT(//10x,'la valeur de jm =',i4,2x,'lue sur le fichier de dema
     *rrage est differente de la valeur parametree jjm =',i4//)
   3  FORMAT(//10x,'la valeur de lllm =',i4,2x,'lue sur le fichier demar
     *rage est differente de la valeur parametree llm =',i4//)

      write(*,*) "newstart: All is well that ends well."

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine load_MONS_data(MONS_Hdn,MONS_d21)
      
      use datafile_mod, only:datadir
      
      implicit none
      ! routine to load Benedicte Diez MONS dataset, fill in date in southern
      ! polar region, and interpolate the result onto the GCM grid
      include"dimensions.h"
      include"paramet.h"
      include"comgeom.h"
      ! arguments:
      real,intent(out) :: MONS_Hdn(iip1,jjp1) ! Hdn: %WEH=Mass fraction of H2O
      real,intent(out) :: MONS_d21(iip1,jjp1) ! ice table "depth" (in kg/m2)
      ! N.B MONS datasets should be of dimension (iip1,jjp1)
      ! local variables:
      character(len=88) :: filename="results_MONS_lat_lon_H_depth.txt"
      character(len=88) :: txt ! to store some text
      integer :: ierr,i,j
      integer,parameter :: nblon=180 ! number of longitudes of MONS datasets
      integer,parameter :: nblat=90 ! number of latitudes of MONS datasets
      real :: pi
      real :: longitudes(nblon) ! MONS dataset longitudes
      real :: latitudes(nblat)  ! MONS dataset latitudes
      ! MONS dataset: mass fraction of H2O where H is assumed to be in H2O
      real :: Hdn(nblon,nblat)
      real :: d21(nblon,nblat)! MONS dataset "depth" (g/cm2)

      ! Extended MONS dataset (for interp_horiz)
      real :: Hdnx(nblon+1,nblat)
      real :: d21x(nblon+1,nblat)
      real :: lon_bound(nblon+1) ! longitude boundaries
      real :: lat_bound(nblat-1) ! latitude boundaries

      ! 1. Initializations:

      write(*,*) "Loading MONS data"

      ! Open MONS datafile:
      open(42,file=trim(datadir)//"/"//trim(filename),
     &     status="old",iostat=ierr)
      if (ierr/=0) then
        write(*,*) "Error in load_MONS_data:"
        write(*,*) "Failed opening file ",
     &             trim(datadir)//"/"//trim(filename)
        write(*,*)'1) You can change this directory address in ',
     &            'callfis.def with'
        write(*,*)'   datadir=/path/to/datafiles'
        write(*,*)'2) If necessary ',trim(filename),
     &                 ' (and other datafiles)'
        write(*,*)'   can be obtained online at:'
        write(*,*)'http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
        CALL ABORT
      else ! skip first line of file (dummy read)
         read(42,*) txt
      endif

      pi=2.*asin(1.)
      
      !2. Load MONS data (on MONS grid)
      do j=1,nblat
        do i=1,nblon
        ! swap latitude index so latitudes go from north pole to south pole:
          read(42,*) latitudes(nblat-j+1),longitudes(i),
     &               Hdn(i,nblat-j+1),d21(i,nblat-j+1)
        ! multiply d21 by 10 to convert from g/cm2 to kg/m2
          d21(i,nblat-j+1)=d21(i,nblat-j+1)*10.0
        enddo
      enddo
      close(42)
      
      ! there is unfortunately no d21 data for latitudes -77 to -90
      ! so we build some by linear interpolation between values at -75
      ! and assuming d21=0 at the pole
      do j=84,90 ! latitudes(84)=-77 ; latitudes(83)=-75
        do i=1,nblon
          d21(i,j)=d21(i,83)*((latitudes(j)+90)/15.0)
        enddo
      enddo

      ! 3. Build extended MONS dataset & boundaries (for interp_horiz)
      ! longitude boundaries (in radians):
      do i=1,nblon
        ! NB: MONS data is every 2 degrees in longitude
        lon_bound(i)=(longitudes(i)+1.0)*pi/180.0
      enddo
      ! extra 'modulo' value
      lon_bound(nblon+1)=lon_bound(1)+2.0*pi
      
      ! latitude boundaries (in radians):
      do j=1,nblat-1
        ! NB: Mons data is every 2 degrees in latitude
        lat_bound(j)=(latitudes(j)-1.0)*pi/180.0
      enddo
      
      ! MONS datasets:
      do j=1,nblat
        Hdnx(1:nblon,j)=Hdn(1:nblon,j)
        Hdnx(nblon+1,j)=Hdnx(1,j)
        d21x(1:nblon,j)=d21(1:nblon,j)
        d21x(nblon+1,j)=d21x(1,j)
      enddo
      
      ! Interpolate onto GCM grid
      call interp_horiz(Hdnx,MONS_Hdn,nblon,nblat-1,iim,jjm,1,
     &                  lon_bound,lat_bound,rlonu,rlatv)
      call interp_horiz(d21x,MONS_d21,nblon,nblat-1,iim,jjm,1,
     &                  lon_bound,lat_bound,rlonu,rlatv)
      
      end subroutine
