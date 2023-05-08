      program rcm1d

! to use  'getin'
      use ioipsl_getincom, only: getin
      use dimphy, only : init_dimphy
      use mod_grid_phy_lmdz, only : regular_lonlat
      use infotrac, only: nqtot, tname
      use tracer_h, only: noms
      use surfdat_h, only: albedodat, phisfi, dryness, watercaptag,
     &                     zmea, zstd, zsig, zgam, zthe,
     &                     emissiv, emisice, iceradius,
     &                     dtemisice
      use comdiurn_h, only: sinlat, coslat, sinlon, coslon
      use comsoil_h, only: nsoilmx, layer, mlayer, inertiedat, volcapa
      use phyredem, only: physdem0,physdem1
      use geometry_mod, only: init_geometry
      use slab_ice_h, only: noceanmx
      use planete_mod, only: apoastr,periastr,year_day,peri_day,
     &         obliquit,nres,z0,lmixmin,emin_turb,coefvis,coefir,
     &         timeperi,e_elips,p_elips
      use comcstfi_mod, only: pi, cpp, rad, g, r,
     &                        mugaz, rcp, omeg
      use time_phylmdz_mod, only: daysec, dtphys, day_step, ecritphy,
     &                            nday, iphysiq
      use callkeys_mod, only: tracer,check_cpp_match,rings_shadow,
     &	                      specOLR,water,pceil,ok_slab_ocean,photochem
      USE comvert_mod, ONLY: ap,bp,aps,bps,pa,preff, sig,
     &                       presnivs,pseudoalt,scaleheight
      USE vertical_layers_mod, ONLY: init_vertical_layers
      USE logic_mod, ONLY: hybrid,autozlevs
      use regular_lonlat_mod, only: init_regular_lonlat
      use planete_mod, only: ini_planete_mod
      use physics_distribution_mod, only: init_physics_distribution
      use regular_lonlat_mod, only: init_regular_lonlat
      use mod_interface_dyn_phys, only: init_interface_dyn_phys
      use inifis_mod, only: inifis
      use phys_state_var_mod, only: phys_state_var_init
      use physiq_mod, only: physiq
      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Run the physics package of the universal model in a 1D column.
!     
!     It can be compiled with a command like (e.g. for 25 layers):
!                                  "makegcm -p std -d 25 rcm1d"
!     It requires the files "callphys.def", "z2sig.def",
!     "traceur.def", and "run.def" with a line "INCLUDEDEF=callphys.def"
!
!     Authors
!     -------
!     Frederic Hourdin
!     R. Fournier
!     F. Forget
!     F. Montmessin (water ice added)
!     R. Wordsworth
!     B. Charnay
!     A. Spiga
!     J. Leconte (2012)
!
!==================================================================

#include "dimensions.h"
#include "paramet.h"
!include "dimphys.h"
#include "netcdf.inc"
#include "comgeom.h"

c --------------------------------------------------------------
c  Declarations
c --------------------------------------------------------------
c
      INTEGER unitstart      ! unite d'ecriture de "startfi"
      INTEGER nlayer,nlevel,nsoil,ndt
      INTEGER ilayer,ilevel,isoil,idt,iq
      LOGICAl firstcall,lastcall
c
      INTEGER day0          ! date initial (sol ; =0 a Ls=0)
      REAL day              ! date durant le run
      REAL time             ! time (0<time<1 ; time=0.5 a midi)
      REAL play(llm)        ! Pressure at the middle of the layers (Pa)
      REAL plev(llm+1)      ! intermediate pressure levels (pa)
      REAL psurf,tsurf(1)      
      REAL u(llm),v(llm)    ! zonal, meridional wind
      REAL gru,grv          ! prescribed "geostrophic" background wind
      REAL temp(llm)        ! temperature at the middle of the layers
      REAL,ALLOCATABLE :: q(:,:)      ! tracer mixing ratio (e.g. kg/kg)
      REAL,ALLOCATABLE :: qsurf(:)    ! tracer surface budget (e.g. kg.m-2)
      REAL,ALLOCATABLE :: tsoil(:)    ! subsurface soil temperature (K)
!      REAL co2ice               ! co2ice layer (kg.m-2) !not used anymore
      integer :: i_co2_ice=0     ! tracer index of co2 ice
      integer :: i_h2o_ice=0     ! tracer index of h2o ice
      integer :: i_h2o_vap=0     ! tracer index of h2o vapor
      REAL emis(1)               ! surface layer
      REAL q2(llm+1)             ! Turbulent Kinetic Energy
      REAL zlay(llm)             ! altitude estimee dans les couches (km)

c    Physical and dynamical tandencies (e.g.  m.s-2, K/s, Pa/s)
      REAL du(llm),dv(llm),dtemp(llm)
      REAL dudyn(llm),dvdyn(llm),dtempdyn(llm)
      REAL dpsurf(1)    
      REAL,ALLOCATABLE :: dq(:,:)
      REAL,ALLOCATABLE :: dqdyn(:,:)

c   Various intermediate variables
!      INTEGER thermo
      REAL zls
      REAL phi(llm),h(llm),s(llm)
      REAL pks, ptif, w(llm)
      INTEGER ierr, aslun
      REAL tmp1(0:llm),tmp2(0:llm)
      integer :: nq !=1 ! number of tracers
 
      character*2 str2
      character (len=7) :: str7
      character(len=44) :: txt

      logical oldcompare, earthhack,saveprofile

!     added by RW for zlay computation
      real Hscale, Hmax, rho, dz

!     added by RW for autozlevs computation
      real nu, xx, pMIN, zlev, Htop
      real logplevs(llm)

!     added by BC
      REAL cloudfrac(1,llm)
      REAL hice(1),totcloudfrac(1)

!     added by BC for ocean
      real rnat(1)
      REAL tslab(1,noceanmx),tsea_ice(1),sea_ice(1)
      real pctsrf_sic(1)



!     added by AS to avoid the use of adv trac common
      character*30,allocatable :: nametmp(:)   !

      real :: latitude(1), longitude(1), cell_area(1)

c=======================================================================
c INITIALISATION
c=======================================================================

      ! read nq from traceur.def
      open(90,file='traceur.def',status='old',form='formatted',
     &       iostat=ierr)
      if (ierr.eq.0) then
        read(90,*,iostat=ierr) nq
      else
        nq=0
      endif
      close(90)
      
      ! Initialize dimphy module
      call init_dimphy(1,llm) 
      ! now initialize arrays using phys_state_var_init
      call phys_state_var_init(nq)
      
      saveprofile=.false.
      saveprofile=.true.

c ----------------------------------------
c  Default values  (corresponding to Mars)
c ----------------------------------------

      pi=2.E+0*asin(1.E+0)

c     Parametres Couche limite et Turbulence 
c     --------------------------------------
      z0 =  1.e-2                ! surface roughness (m) ~0.01 
      emin_turb = 1.e-6          ! energie minimale ~1.e-8
      lmixmin = 30               ! longueur de melange ~100
 
c     propriete optiques des calottes et emissivite du sol
c     ----------------------------------------------------
      emissiv= 0.95              ! Emissivite du sol martien ~.95
      emisice(1)=0.95            ! Emissivite calotte nord
      emisice(2)=0.95            ! Emissivite calotte sud  

      iceradius(1) = 100.e-6     ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6     ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 2.          ! time scale for snow metamorphism (north)
      dtemisice(2) = 2.          ! time scale for snow metamorphism (south
      hybrid=.false.

c ------------------------------------------------------
c  Load parameters from "run.def" and "gases.def" 
c ------------------------------------------------------

! check if 'rcm1d.def' file is around 
      open(90,file='rcm1d.def',status='old',form='formatted',
     &     iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) 'Cannot find required file "rcm1d.def"'
        write(*,*) 'which should contain some input parameters'
        write(*,*) ' ... might as well stop here ...'
        stop
      else
        close(90)
      endif

! now, run.def is needed anyway. so we create a dummy temporary one
! for ioipsl to work. if a run.def is already here, stop the
! program and ask the user to do a bit of cleaning
      open(90,file='run.def',status='old',form='formatted',
     &     iostat=ierr)
      if (ierr.eq.0) then
        close(90)
        write(*,*) 'There is already a run.def file.'
        write(*,*) 'This is not compatible with 1D runs.'
        write(*,*) 'Please remove the file and restart the run.'
        write(*,*) 'Runtime parameters are supposed to be in rcm1d.def'
        stop
      else
        call system('touch run.def')
        call system("echo 'INCLUDEDEF=callphys.def' >> run.def")
        call system("echo 'INCLUDEDEF=rcm1d.def' >> run.def")
      endif

! check if we are going to run with or without tracers
      write(*,*) "Run with or without tracer transport ?"
      tracer=.false. ! default value
      call getin("tracer",tracer)
      write(*,*) " tracer = ",tracer

! OK. now that run.def has been read once -- any variable is in memory.
! so we can dump the dummy run.def
!      call system("rm -rf run.def") ! Ehouarn: delay this to after inifis

! while we're at it, check if there is a 'traceur.def' file
! and preocess it, if necessary. Otherwise initialize tracer names
      if (tracer) then
      ! load tracer names from file 'traceur.def'
        open(90,file='traceur.def',status='old',form='formatted',
     &       iostat=ierr)
        if (ierr.eq.0) then
          write(*,*) "rcm1d: Reading file traceur.def"
          ! read number of tracers:
          read(90,*,iostat=ierr) nq
          nqtot=nq ! set value of nqtot (in infotrac module) as nq
          if (ierr.ne.0) then
            write(*,*) "rcm1d: error reading number of tracers"
            write(*,*) "   (first line of traceur.def) "
            stop
          endif
          if (nq>0) then
            allocate(tname(nq))
            allocate(noms(nq))
            allocate(q(llm,nq))
            allocate(qsurf(nq))
            allocate(dq(llm,nq))
            allocate(dqdyn(llm,nq))
          else
            write(*,*) "rcm1d: Error. You set tracer=.true."
            write(*,*) "       but # of tracers in traceur.def is ",nq
            stop
          endif
        
          do iq=1,nq
          ! minimal version, just read in the tracer names, 1 per line
            read(90,*,iostat=ierr) tname(iq)
            noms(iq)=tname(iq)
            if (ierr.ne.0) then
              write(*,*) 'rcm1d: error reading tracer names...'
              stop
            endif
          enddo !of do iq=1,nq
! check for co2_ice / h2o_ice tracers:
         i_co2_ice=0
         i_h2o_ice=0
         i_h2o_vap=0
         do iq=1,nq
           if (tname(iq)=="co2_ice") then
             i_co2_ice=iq
           elseif (tname(iq)=="h2o_ice") then
             i_h2o_ice=iq
           elseif (tname(iq)=="h2o_vap") then
             i_h2o_vap=iq
           endif
         enddo
        else
          write(*,*) 'Cannot find required file "traceur.def"'
          write(*,*) ' If you want to run with tracers, I need it'
          write(*,*) ' ... might as well stop here ...'
          stop
        endif
        close(90)


      else ! of if (tracer)
        nqtot=0
        nq=0
        ! still, make allocations for 1 dummy tracer
        allocate(tname(1)) 
        allocate(qsurf(1))
        allocate(q(llm,1))
        allocate(dq(llm,1))
      
       ! Check that tracer boolean is compliant with number of tracers
       ! -- otherwise there is an error (and more generally we have to be consistent)
       if (nq .ge. 1) then
          write(*,*) "------------------------------"
          write(*,*) "rcm1d: You set tracer=.false."
          write(*,*) " But set number of tracers to ",nq
          write(*,*) " > If you want tracers, set tracer=.true."
          write(*,*) "------------------------------"
          stop
       endif
      endif ! of if (tracer)

!!! We have to check that check_cpp_match is F for 1D computations
!!! We think this check is better than make a particular case for 1D in inifis or calc_cpp_mugaz
      check_cpp_match = .false.
      call getin("check_cpp_match",check_cpp_match)
      if (check_cpp_match) then
          print*,"In 1D modeling, check_cpp_match is supposed to be F"
          print*,"Please correct callphys.def"
          stop
      endif

!!! GEOGRAPHICAL INITIALIZATIONS
     !!! AREA. useless in 1D
      cell_area(1)=1.E+0 !JL+EM to have access to the area in the diagfi.nc files.
     !!! GEOPOTENTIAL. useless in 1D because control by surface pressure
      phisfi(1)=0.E+0
     !!! LATITUDE. can be set. 
      latitude=0 ! default value for latitude
      PRINT *,'latitude (in degrees) ?'
      call getin("latitude",latitude)
      write(*,*) " latitude = ",latitude
      latitude=latitude*pi/180.E+0
     !!! LONGITUDE. useless in 1D.
      longitude=0.E+0
      longitude=longitude*pi/180.E+0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PLANETARY CONSTANTS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      g = -99999.
      PRINT *,'GRAVITY in m s-2 ?'
      call getin("g",g)
      IF (g.eq.-99999.) THEN
          PRINT *,"STOP. I NEED g IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> g = ",g
      ENDIF

      rad = -99999.
      PRINT *,'PLANETARY RADIUS in m ?'
      call getin("rad",rad)
      ! Planetary  radius is needed to compute shadow of the rings
      IF (rad.eq.-99999. .and. rings_shadow .eqv. .true.) THEN
          PRINT *,"STOP. I NEED rad IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> rad = ",rad
      ENDIF

      daysec = -99999.
      PRINT *,'LENGTH OF A DAY in s ?'
      call getin("daysec",daysec)
      IF (daysec.eq.-99999.) THEN
          PRINT *,"STOP. I NEED daysec IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> daysec = ",daysec
      ENDIF
      omeg=4.*asin(1.)/(daysec)
      PRINT *,"OK. FROM THIS I WORKED OUT:"
      PRINT *,"--> omeg = ",omeg

      year_day = -99999.
      PRINT *,'LENGTH OF A YEAR in days ?'
      call getin("year_day",year_day)
      IF (year_day.eq.-99999.) THEN
          PRINT *,"STOP. I NEED year_day IN RCM1D.DEF."
          STOP
      ELSE 
          PRINT *,"--> year_day = ",year_day
      ENDIF

      periastr = -99999.
      PRINT *,'MIN DIST STAR-PLANET in AU ?'
      call getin("periastr",periastr)
      IF (periastr.eq.-99999.) THEN
          PRINT *,"STOP. I NEED periastr IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> periastr = ",periastr
      ENDIF

      apoastr = -99999.
      PRINT *,'MAX DIST STAR-PLANET in AU ?'
      call getin("apoastr",apoastr)
      IF (apoastr.eq.-99999.) THEN
          PRINT *,"STOP. I NEED apoastr IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> apoastr = ",apoastr
      ENDIF

      peri_day = -99999.
      PRINT *,'DATE OF PERIASTRON in days ?'
      call getin("peri_day",peri_day)
      IF (peri_day.eq.-99999.) THEN
          PRINT *,"STOP. I NEED peri_day IN RCM1D.DEF."
          STOP
      ELSE IF (peri_day.gt.year_day) THEN
          PRINT *,"STOP. peri_day.gt.year_day"
          STOP
      ELSE  
          PRINT *,"--> peri_day = ", peri_day
      ENDIF 

      obliquit = -99999. 
      PRINT *,'OBLIQUITY in deg ?'
      call getin("obliquit",obliquit)
      IF (obliquit.eq.-99999.) THEN
          PRINT *,"STOP. I NEED obliquit IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> obliquit = ",obliquit
      ENDIF 

      psurf = -99999.
      PRINT *,'SURFACE PRESSURE in Pa ?'
      call getin("psurf",psurf)
      IF (psurf.eq.-99999.) THEN
          PRINT *,"STOP. I NEED psurf IN RCM1D.DEF."
          STOP
      ELSE
          PRINT *,"--> psurf = ",psurf
      ENDIF
      !! we need reference pressures
      pa=psurf/30.
      preff=psurf ! these values are not needed in 1D  (are you sure JL12)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! END PLANETARY CONSTANTS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c  Date et heure locale du debut du run
c  ------------------------------------
c    Date (en sols depuis le solstice de printemps) du debut du run
      day0 = 0 ! default value for day0
      write(*,*) 'Initial date (in martian sols ; =0 at Ls=0)?'
      call getin("day0",day0)
      day=float(day0)
      write(*,*) " day0 = ",day0
c  Heure de demarrage
      time=0 ! default value for time
      write(*,*)'Initial local time (in hours, between 0 and 24)?'
      call getin("time",time)
      write(*,*)" time = ",time
      time=time/24.E+0 ! convert time (hours) to fraction of sol


c  Discretisation (Definition de la grille et des pas de temps)
c  --------------
c
      nlayer=llm
      nlevel=nlayer+1

      day_step=48 ! default value for day_step
      PRINT *,'Number of time steps per sol ?'
      call getin("day_step",day_step)
      write(*,*) " day_step = ",day_step

      iphysiq=1 ! in 1D model physics are called evry time step
      ecritphy=day_step ! default value for ecritphy
      PRINT *,'Nunber of steps between writediagfi ?'
      call getin("ecritphy",ecritphy)
      write(*,*) " ecritphy = ",ecritphy

      ndt=10 ! default value for ndt
      PRINT *,'Number of sols to run ?'
      call getin("ndt",ndt)
      write(*,*) " ndt = ",ndt
      nday=ndt

      ndt=ndt*day_step     
      dtphys=daysec/day_step  
      write(*,*)"-------------------------------------"
      write(*,*)"-------------------------------------"
      write(*,*)"--> Physical timestep is ",dtphys
      write(*,*)"-------------------------------------"
      write(*,*)"-------------------------------------"

      ! initializations, as with iniphysiq.F90 for the 3D GCM
      call init_physics_distribution(regular_lonlat,4,
     &                               1,1,1,nlayer,1)
      call init_interface_dyn_phys
      CALL init_regular_lonlat(1,1,longitude,latitude,
     &                   (/0.,0./),(/0.,0./))
      call init_geometry(1,longitude,latitude,
     &                   (/0.,0.,0.,0./),(/0.,0.,0.,0./),
     &                   cell_area)
! Ehouarn: init_vertial_layers called later (because disvert not called yet)
!      call init_vertical_layers(nlayer,preff,scaleheight,
!     &                      ap,bp,aps,bps,presnivs,pseudoalt)
!      call init_dimphy(1,nlayer) ! Initialize dimphy module
      call ini_planete_mod(nlayer,preff,ap,bp)

!!! CALL INIFIS
!!! - read callphys.def
!!! - calculate sine and cosine of longitude and latitude
!!! - calculate mugaz and cp 
!!! NB: some operations are being done dummily in inifis in 1D
!!! - physical constants: nevermind, things are done allright below
!!! - physical frequency: nevermind, in inifis this is a simple print
      cpp=-9999. ! dummy init for inifis, will be rewrite later on
      r=-9999.   ! dummy init for inifis, will be rewrite later on
      CALL inifis(1,llm,nq,day0,daysec,nday,dtphys,
     .            latitude,longitude,cell_area,rad,g,r,cpp)

      nsoil=nsoilmx
      allocate(tsoil(nsoilmx))
      !! those are defined in comsoil_h.F90
      IF (.not.ALLOCATED(layer)) ALLOCATE(layer(nsoilmx))
      IF (.not.ALLOCATED(mlayer)) ALLOCATE(mlayer(0:nsoilmx-1))
      IF (.not.ALLOCATED(inertiedat)) ALLOCATE(inertiedat(1,nsoilmx))

! At this point, both getin() and getin_p() functions have been used,
! and the run.def file can be removed.
      call system("rm -rf run.def")

!!! We check everything is OK.
      PRINT *,"CHECK"
      PRINT *,"--> mugaz = ",mugaz
      PRINT *,"--> cpp = ",cpp
      r = 8.314511E+0 * 1000.E+0 / mugaz
      rcp = r / cpp

c output spectrum?
      write(*,*)"Output spectral OLR?"
      specOLR=.false.
      call getin("specOLR",specOLR)
      write(*,*)" specOLR = ",specOLR

c
c  pour le schema d'ondes de gravite
c  ---------------------------------
      zmea(1)=0.E+0
      zstd(1)=0.E+0
      zsig(1)=0.E+0
      zgam(1)=0.E+0
      zthe(1)=0.E+0

c    Initialisation des traceurs
c    ---------------------------

      DO iq=1,nq
        DO ilayer=1,nlayer
           q(ilayer,iq) = 0.
        ENDDO
      ENDDO
     
      DO iq=1,nq
        qsurf(iq) = 0.
      ENDDO
      
      if (tracer) then ! Initialize tracers here. 
             
         write(*,*) "rcm1d : initializing tracers profiles"

         do iq=1,nq
         
            txt=""
            write(txt,"(a)") tname(iq)
            write(*,*)"  tracer:",trim(txt)
              
            ! CO2
            if (txt.eq."co2_ice") then
               q(:,iq)=0.   ! kg/kg of atmosphere
               qsurf(iq)=0. ! kg/m2 at the surface               
               ! Look for a "profile_co2_ice" input file
               open(91,file='profile_co2_ice',status='old',
     &         form='formatted',iostat=ierr)
               if (ierr.eq.0) then
                  read(91,*) qsurf(iq)
                  do ilayer=1,nlayer
                     read(91,*) q(ilayer,iq)
                  enddo
               else
                  write(*,*) "No profile_co2_ice file!"
               endif
               close(91)
            endif ! of if (txt.eq."co2")
          
            ! WATER VAPOUR
            if (txt.eq."h2o_vap") then
               q(:,iq)=0.   ! kg/kg of atmosphere
               qsurf(iq)=0. ! kg/m2 at the surface
               ! Look for a "profile_h2o_vap" input file   
               open(91,file='profile_h2o_vap',status='old',
     &         form='formatted',iostat=ierr)
               if (ierr.eq.0) then
                  read(91,*) qsurf(iq)
                  do ilayer=1,nlayer
                     read(91,*) q(ilayer,iq)
                  enddo
               else
                  write(*,*) "No profile_h2o_vap file!"
               endif
               close(91)
            endif ! of if (txt.eq."h2o_vap")
            
            ! WATER ICE
            if (txt.eq."h2o_ice") then
               q(:,iq)=0.   ! kg/kg of atmosphere
               qsurf(iq)=0. ! kg/m2 at the surface
               ! Look for a "profile_h2o_ice" input file
               open(91,file='profile_h2o_ice',status='old',
     &         form='formatted',iostat=ierr)
               if (ierr.eq.0) then
                  read(91,*) qsurf(iq)
                  do ilayer=1,nlayer
                     read(91,*) q(ilayer,iq)
                  enddo
               else
                  write(*,*) "No profile_h2o_ice file!"
               endif
               close(91)
            endif ! of if (txt.eq."h2o_ice")

         enddo ! of do iq=1,nq
         
      endif ! of tracer

c   Initialisation pour prendre en compte les vents en 1-D
c   ------------------------------------------------------
      ptif=2.E+0*omeg*sinlat(1)
 

c    vent geostrophique
      gru=10. ! default value for gru
      PRINT *,'zonal eastward component of the geostrophic wind (m/s) ?'
      call getin("u",gru)
      write(*,*) " u = ",gru
      grv=0. !default value for grv
      PRINT *,'meridional northward component of the geostrophic',
     &' wind (m/s) ?'
      call getin("v",grv)
      write(*,*) " v = ",grv

! To be clean, also set vertical winds to zero
      w(1:nlayer)=0

c     Initialisation des vents  au premier pas de temps
      DO ilayer=1,nlayer
         u(ilayer)=gru
         v(ilayer)=grv
      ENDDO

c     energie cinetique turbulente
      DO ilevel=1,nlevel
         q2(ilevel)=0.E+0
      ENDDO

c  emissivity / surface co2 ice ( + h2o ice??)
c  -------------------------------------------
      emis(1)=emissiv ! default value for emissivity
      PRINT *,'Emissivity of bare ground ?'
      call getin("emis",emis(1))
      write(*,*) " emis = ",emis(1)
      emissiv=emis(1) ! we do this so that condense_co2 sets things to the right
                   ! value if there is no snow

      if(i_co2_ice.gt.0)then
         qsurf(i_co2_ice)=0 ! default value for co2ice
         print*,'Initial CO2 ice on the surface (kg.m-2)'
         call getin("co2ice",qsurf(i_co2_ice))
         write(*,*) " co2ice = ",qsurf(i_co2_ice)
         IF (qsurf(i_co2_ice).ge.1.E+0) THEN
            ! if we have some CO2 ice on the surface, change emissivity
            if (latitude(1).ge.0) then ! northern hemisphere
              emis(1)=emisice(1)
            else ! southern hemisphere
              emis(1)=emisice(2)
            endif
         ENDIF
      endif

c  calcul des pressions et altitudes en utilisant les niveaux sigma
c  ----------------------------------------------------------------

c    Vertical Coordinates
c    """"""""""""""""""""
      hybrid=.true.
      PRINT *,'Hybrid coordinates ?'
      call getin("hybrid",hybrid)
      write(*,*) " hybrid = ", hybrid


      autozlevs=.false.
      PRINT *,'Auto-discretise vertical levels ?'
      call getin("autozlevs",autozlevs)
      write(*,*) " autozlevs = ", autozlevs

      pceil = psurf / 1000.0 ! Pascals
      PRINT *,'Ceiling pressure (Pa) ?'
      call getin("pceil",pceil)
      write(*,*) " pceil = ", pceil

! Test of incompatibility:
! if autozlevs used, cannot have hybrid too

      if (autozlevs.and.hybrid) then
         print*,'Cannot use autozlevs and hybrid together!'
         call abort
      endif

      if(autozlevs)then
            
         open(91,file="z2sig.def",form='formatted')
         read(91,*) Hscale
         DO ilayer=1,nlayer-2
            read(91,*) Hmax
         enddo
         close(91)
 
            
         print*,'Hmax = ',Hmax,' km'
         print*,'Auto-shifting Hscale to:'
!     Hscale = Hmax / log(psurf/100.0)
         Hscale = Hmax / log(psurf/pceil)
         print*,'Hscale = ',Hscale,' km'
         
! none of this matters if we dont care about zlay
         
      endif

      call disvert
      ! now that disvert has been called, initialize module vertical_layers_mod
      call init_vertical_layers(nlayer,preff,scaleheight,
     &                      ap,bp,aps,bps,presnivs,pseudoalt)

         if(.not.autozlevs)then
            ! we want only the scale height from z2sig, in order to compute zlay
            open(91,file="z2sig.def",form='formatted')
            read(91,*) Hscale
            close(91)
         endif

!         if(autozlevs)then
!            open(94,file="Hscale.temp",form='formatted')
!            read(94,*) Hscale
!            close(94)
!         endif

         DO ilevel=1,nlevel
            plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
         ENDDO
         
         DO ilayer=1,nlayer
            play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
         ENDDO
         


         DO ilayer=1,nlayer
!     zlay(ilayer)=-300.E+0 *r*log(play(ilayer)/plev(1))
!     &   /g
            zlay(ilayer)=-1000.0*Hscale*log(play(ilayer)/plev(1))
         ENDDO

!      endif

c  profil de temperature au premier appel
c  --------------------------------------
      pks=psurf**rcp

c altitude en km dans profile: on divise zlay par 1000
      tmp1(0)=0.E+0
      DO ilayer=1,nlayer
        tmp1(ilayer)=zlay(ilayer)/1000.E+0
      ENDDO
      call profile(nlayer+1,tmp1,tmp2)

      tsurf(1)=tmp2(0)
      DO ilayer=1,nlayer
        temp(ilayer)=tmp2(ilayer)
      ENDDO
      print*,"check"
      PRINT*,"INPUT SURFACE TEMPERATURE",tsurf(1)
      PRINT*,"INPUT TEMPERATURE PROFILE",temp

c  Initialisation albedo / inertie du sol
c  --------------------------------------
      albedodat(1)=0.2 ! default value for albedodat
      PRINT *,'Albedo of bare ground ?'
      call getin("albedo",albedodat(1))
      write(*,*) " albedo = ",albedodat(1)

      inertiedat(1,1)=400 ! default value for inertiedat
      PRINT *,'Soil thermal inertia (SI) ?'
      call getin("inertia",inertiedat(1,1))
      write(*,*) " inertia = ",inertiedat(1,1)

! Initialize soil properties and temperature
! ------------------------------------------
      volcapa=1.e6 ! volumetric heat capacity
      DO isoil=1,nsoil
         inertiedat(1,isoil)=inertiedat(1,1) ! soil thermal inertia
         tsoil(isoil)=tsurf(1)  ! soil temperature
      ENDDO

! Initialize depths
! -----------------
      do isoil=0,nsoil-1
        mlayer(isoil)=3.e-2*(2.**(isoil-0.5)) ! mid-layer depth
      enddo
      do isoil=1,nsoil
        layer(isoil)=3.e-2*(2.**(isoil-1)) ! layer depth
      enddo

! Initialize cloud fraction and oceanic ice
! -----------------------------------------
      hice=0.
      do ilayer=1,llm
        cloudfrac(1,ilayer)=0.
      enddo
      totcloudfrac=0.

! Initialize slab ocean
! -----------------
      rnat=1. ! default value for rnat
      if(inertiedat(1,1).GE.10000.)then
         rnat=0.
      endif
      if(ok_slab_ocean)then
      rnat=0.
      tslab(1,1)=tsurf(1)
      tslab(1,2)=tsurf(1)
      tsea_ice=tsurf
      pctsrf_sic=0.
      sea_ice=0.
      endif


! Initialize chemical species
! -----------------
#ifndef MESOSCALE
      if(tracer.and.photochem) then
           call initracer(1,nq,tname)
           allocate(nametmp(nq))
           nametmp(1:nq)=tname(1:nq)
           call inichim_1D(nq, q, qsurf, psurf, 0, 0)
           tname(1:nq)=nametmp(1:nq)
           noms(1:nq)=nametmp(1:nq)
      endif ! tracer and photochem
#endif


c  Write a "startfi" file
c  --------------------
c  This file will be read during the first call to "physiq".
c  It is needed to transfert physics variables to "physiq"...

      call physdem0("startfi.nc",longitude,latitude,nsoilmx,1,llm,nq,
     &              dtphys,real(day0),time,cell_area,
     &              albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe)
      call physdem1("startfi.nc",nsoilmx,1,llm,nq,
     &                dtphys,time,
     &                tsurf,tsoil,emis,q2,qsurf,
     &                cloudfrac,totcloudfrac,hice, 
     &                rnat,pctsrf_sic,tslab,tsea_ice,sea_ice)

c=======================================================================
c  BOUCLE TEMPORELLE DU MODELE 1D 
c=======================================================================

      firstcall=.true.
      lastcall=.false.

      DO idt=1,ndt
        IF (idt.eq.ndt) then       !test
         lastcall=.true.
         call stellarlong(day*1.0,zls)
!         write(103,*) 'Ls=',zls*180./pi
!         write(103,*) 'Lat=', latitude(1)*180./pi
!         write(103,*) 'RunEnd - Atmos. Temp. File'
!         write(103,*) 'RunEnd - Atmos. Temp. File'
!         write(104,*) 'Ls=',zls*180./pi
!         write(104,*) 'Lat=', latitude(1)
!         write(104,*) 'RunEnd - Atmos. Temp. File'
        ENDIF

c    calcul du geopotentiel 
c     ~~~~~~~~~~~~~~~~~~~~~


      DO ilayer=1,nlayer

!              if(autozlevs)then
!                 s(ilayer)=(play(ilayer)/psurf)**rcp
!              else
          s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
!              endif
              !s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
          h(ilayer)=cpp*temp(ilayer)/(pks*s(ilayer))
       ENDDO

!      DO ilayer=1,nlayer
!        s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
!        h(ilayer)=cpp*temp(ilayer)/(pks*s(ilayer))
!      ENDDO
      phi(1)=pks*h(1)*(1.E+0-s(1))
      DO ilayer=2,nlayer
         phi(ilayer)=phi(ilayer-1)+
     &               pks*(h(ilayer-1)+h(ilayer))*.5E+0
     &                  *(s(ilayer-1)-s(ilayer))

      ENDDO

c       appel de la physique
c       --------------------


      CALL physiq (1,llm,nq,
     .     tname,
     ,     firstcall,lastcall,
     ,     day,time,dtphys,
     ,     plev,play,phi,
     ,     u, v,temp, q,  
     ,     w,
C - sorties
     s     du, dv, dtemp, dq,dpsurf)


c       evolution du vent : modele 1D
c       -----------------------------
 
c       la physique calcule les derivees temporelles de u et v.
c       on y rajoute betement un effet Coriolis.
c
c       DO ilayer=1,nlayer
c          du(ilayer)=du(ilayer)+ptif*(v(ilayer)-grv)
c          dv(ilayer)=dv(ilayer)+ptif*(-u(ilayer)+gru)
c       ENDDO

c       Pour certain test : pas de coriolis a l'equateur
c       if(latitude(1).eq.0.) then
          DO ilayer=1,nlayer
             du(ilayer)=du(ilayer)+ (gru-u(ilayer))/1.e4
             dv(ilayer)=dv(ilayer)+ (grv-v(ilayer))/1.e4
          ENDDO
c       end if
c      
c
c       Calcul du temps au pas de temps suivant
c       ---------------------------------------
        firstcall=.false.
        time=time+dtphys/daysec
        IF (time.gt.1.E+0) then
            time=time-1.E+0
            day=day+1
        ENDIF

c       calcul des vitesses et temperature au pas de temps suivant
c       ----------------------------------------------------------

        DO ilayer=1,nlayer
           u(ilayer)=u(ilayer)+dtphys*du(ilayer)
           v(ilayer)=v(ilayer)+dtphys*dv(ilayer)
           temp(ilayer)=temp(ilayer)+dtphys*dtemp(ilayer)
        ENDDO

c       calcul des pressions au pas de temps suivant
c       ----------------------------------------------------------

           psurf=psurf+dtphys*dpsurf(1)   ! evolution de la pression de surface
           DO ilevel=1,nlevel
              plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
           ENDDO
           DO ilayer=1,nlayer
                 play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
           ENDDO

c       calcul traceur au pas de temps suivant
c       --------------------------------------

        DO iq = 1, nq
          DO ilayer=1,nlayer
             q(ilayer,iq)=q(ilayer,iq)+dtphys*dq(ilayer,iq)
          ENDDO
        END DO

c    ========================================================
c    GESTION DES SORTIE
c    ========================================================
      if(saveprofile)then
         OPEN(12,file='profile.out',form='formatted')
         write(12,*) tsurf
         DO ilayer=1,llm
            write(12,*) temp(ilayer) !, play(ilayer) !AS12 only temp so that iprofile=8 can be used
         ENDDO
         CLOSE(12)
      endif


      ENDDO   ! fin de la boucle temporelle

      write(*,*) "rcm1d: Everything is cool."

c    ========================================================
      end                       !rcm1d
 

