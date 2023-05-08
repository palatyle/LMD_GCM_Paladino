
      PROGRAM testphys1d
! to use  'getin'
      USE ioipsl_getincom, only: getin
      use dimphy, only : init_dimphy
      use mod_grid_phy_lmdz, only : regular_lonlat
      use infotrac, only: nqtot, tname
      use comsoil_h, only: volcapa, layer, mlayer, inertiedat, nsoilmx
      use comgeomfi_h, only: sinlat, ini_fillgeom
      use surfdat_h, only: albedodat, z0_default, emissiv, emisice,
     &                     albedice, iceradius, dtemisice, z0,
     &                     zmea, zstd, zsig, zgam, zthe, phisfi,
     &                     watercaptag, hmons, summit, base
      use slope_mod, only: theta_sl, psi_sl
      use phyredem, only: physdem0,physdem1
      use geometry_mod, only: init_geometry
      use planete_h, only: year_day, periheli, aphelie, peri_day,
     &                     obliquit, emin_turb, lmixmin
      use comcstfi_h, only: pi, rad, omeg, g, mugaz, rcp, r, cpp
      use time_phylmdz_mod, only: daysec, dtphys, day_step,
     &                            ecritphy, iphysiq
      use dimradmars_mod, only: tauscaling,tauvis,totcloudfrac
      use co2cloud_mod, only: mem_Mccn_co2, mem_Mh2o_co2,
     &                        mem_Nccn_co2
      USE comvert_mod, ONLY: ap,bp,aps,bps,pa,preff,sig,
     &                       presnivs,pseudoalt,scaleheight
      USE vertical_layers_mod, ONLY: init_vertical_layers
      USE logic_mod, ONLY: hybrid
      use physics_distribution_mod, only: init_physics_distribution
      use regular_lonlat_mod, only: init_regular_lonlat
      use mod_interface_dyn_phys, only: init_interface_dyn_phys
      USE phys_state_var_init_mod, ONLY: phys_state_var_init
      USE physiq_mod, ONLY: physiq
      IMPLICIT NONE

c=======================================================================
c   subject:
c   --------
c   PROGRAM useful to run physical part of the martian GCM in a 1D column
c       
c Can be compiled with a command like (e.g. for 25 layers)
c  "makegcm -p mars -d 25 testphys1d"
c It requires the files "testphys1d.def" "callphys.def"
c   and a 'run.def' file (containing a "INCLUDEDEF=callphys.def" line)
c      and a file describing the sigma layers (e.g. "z2sig.def")
c
c   author: Frederic Hourdin, R.Fournier,F.Forget
c   -------
c   
c   update: 12/06/2003 including chemistry (S. Lebonnois) 
c                            and water ice (F. Montmessin)
c 
c=======================================================================

      include "dimensions.h"
      integer, parameter :: ngrid = 1 !(2+(jjm-1)*iim - 1/jjm)
      integer, parameter :: nlayer = llm
!#include "dimradmars.h"
!#include "comgeomfi.h"
!#include "surfdat.h"
!#include "slope.h"
!#include "comsoil.h"
!#include "comdiurn.h"
      include "callkeys.h"
!#include "comsaison.h"
!#include "control.h"
      include "netcdf.inc"
      include "comg1d.h"
!#include "advtrac.h"

c --------------------------------------------------------------
c  Declarations
c --------------------------------------------------------------
c
      INTEGER unitstart      ! unite d'ecriture de "startfi"
      INTEGER nlevel,nsoil,ndt
      INTEGER ilayer,ilevel,isoil,idt,iq
      LOGICAl firstcall,lastcall
c
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
c
      INTEGER day0          ! date initial (sol ; =0 a Ls=0)
      REAL day           ! date durant le run
      REAL time             ! time (0<time<1 ; time=0.5 a midi)
      REAL play(nlayer)   ! Pressure at the middle of the layers (Pa)
      REAL plev(nlayer+1) ! intermediate pressure levels (pa)
      REAL psurf,tsurf(1)      
      REAL u(nlayer),v(nlayer)  ! zonal, meridional wind
      REAL gru,grv   ! prescribed "geostrophic" background wind
      REAL temp(nlayer)   ! temperature at the middle of the layers
      REAL,ALLOCATABLE :: q(:,:) ! tracer mixing ratio (e.g. kg/kg)
      REAL,ALLOCATABLE :: qsurf(:) ! tracer surface budget (e.g. kg.m-2)
      REAL tsoil(nsoilmx)   ! subsurface soik temperature (K)
      REAL co2ice(1)        ! co2ice layer (kg.m-2)
      REAL emis(1)          ! surface layer
      REAL albedo(1,1)      ! surface albedo
      REAL :: wstar(1)=0.    ! Thermals vertical velocity
      REAL q2(nlayer+1)   ! Turbulent Kinetic Energy
      REAL zlay(nlayer)   ! altitude estimee dans les couches (km)

c    Physical and dynamical tandencies (e.g.  m.s-2, K/s, Pa/s)
      REAL du(nlayer),dv(nlayer),dtemp(nlayer)
      REAL dudyn(nlayer),dvdyn(nlayer),dtempdyn(nlayer)
      REAL dpsurf(1)    
      REAL,ALLOCATABLE :: dq(:,:)
      REAL,ALLOCATABLE :: dqdyn(:,:)

c   Various intermediate variables
      INTEGER thermo
      REAL zls
      REAL phi(nlayer),h(nlayer),s(nlayer)
      REAL pks, ptif, w(nlayer)
      REAL qtotinit,qtot
      real,allocatable :: mqtot(:)
      INTEGER ierr, aslun
      REAL tmp1(0:nlayer),tmp2(0:nlayer)
      integer :: nq=1 ! number of tracers
      real :: latitude(1), longitude(1), cell_area(1)

      character*2 str2
      character (len=7) :: str7
      character(len=44) :: txt

c   New flag to compute paleo orbital configurations + few variables JN
      REAL halfaxe, excentric, Lsperi
      Logical paleomars

c=======================================================================

c=======================================================================
c INITIALISATION
c=======================================================================
! initialize "serial/parallel" related stuff
!      CALL init_phys_lmdz(iim,jjp1,llm,1,(/(jjm-1)*iim+2/))
!      CALL init_phys_lmdz(1,1,llm,1,(/1/))
!      call initcomgeomphy

c ------------------------------------------------------
c  Prescribed constants to be set here
c ------------------------------------------------------

      pi=2.E+0*asin(1.E+0)

c     Mars planetary constants
c     ----------------------------
      rad=3397200.               ! mars radius (m)  ~3397200 m
      daysec=88775.              ! length of a sol (s)  ~88775 s
      omeg=4.*asin(1.)/(daysec)  ! rotation rate (rad.s-1)
      g=3.72                     ! gravity (m.s-2) ~3.72  
      mugaz=43.49                ! atmosphere mola mass (g.mol-1) ~43.49
      rcp=.256793                ! = r/cp  ~0.256793
      r= 8.314511E+0 *1000.E+0/mugaz
      cpp= r/rcp
      year_day = 669             ! lenght of year (sols) ~668.6
      periheli = 206.66          ! minimum sun-mars distance (Mkm) ~206.66
      aphelie = 249.22           ! maximum sun-mars distance (Mkm) ~249.22
      halfaxe = 227.94           ! demi-grand axe de l'ellipse
      peri_day =  485.           ! perihelion date (sols since N. Spring)
      obliquit = 25.2            ! Obliquity (deg) ~25.2         
      excentric = 0.0934         ! Eccentricity (0.0934)         
 
c     Planetary Boundary Layer and Turbulence parameters 
c     --------------------------------------------------
      z0_default =  1.e-2        ! surface roughness (m) ~0.01 
      emin_turb = 1.e-6          ! minimal turbulent energy ~1.e-8
      lmixmin = 30               ! mixing length ~100
 
c     cap properties and surface emissivities
c     ----------------------------------------------------
      emissiv= 0.95              ! Bare ground emissivity ~.95
      emisice(1)=0.95            ! Northern cap emissivity
      emisice(2)=0.95            ! Southern cap emisssivity
      albedice(1)=0.5            ! Northern cap albedo
      albedice(2)=0.5            ! Southern cap albedo
      iceradius(1) = 100.e-6     ! mean scat radius of CO2 snow (north)
      iceradius(2) = 100.e-6     ! mean scat radius of CO2 snow (south)
      dtemisice(1) = 2.          ! time scale for snow metamorphism (north)
      dtemisice(2) = 2.          ! time scale for snow metamorphism (south

c     mesh surface (not a very usefull quantity in 1D)
c     ----------------------------------------------------
      cell_area(1)=1.E+0
      
c ------------------------------------------------------
c  Loading run parameters from "run.def" file
c ------------------------------------------------------


! check if 'run.def' file is around (otherwise reading parameters
! from callphys.def via getin() routine won't work.
      open(99,file='run.def',status='old',form='formatted',
     &     iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) 'Cannot find required file "run.def"'
        write(*,*) '  (which should contain some input parameters'
        write(*,*) '   along with the following line:'
        write(*,*) '   INCLUDEDEF=callphys.def'
        write(*,*) '   )'
        write(*,*) ' ... might as well stop here ...'
        stop
      else
        close(99)
      endif

! check if we are going to run with or without tracers
      write(*,*) "Run with or without tracer transport ?"
      tracer=.false. ! default value
      call getin("tracer",tracer)
      write(*,*) " tracer = ",tracer

! while we're at it, check if there is a 'traceur.def' file
! and process it.
      if (tracer) then
      ! load tracer names from file 'traceur.def'
        open(90,file='traceur.def',status='old',form='formatted',
     &       iostat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'Cannot find required file "traceur.def"'
          write(*,*) ' If you want to run with tracers, I need it'
          write(*,*) ' ... might as well stop here ...'
          stop
        else
          write(*,*) "testphys1d: Reading file traceur.def"
          ! read number of tracers:
          read(90,*,iostat=ierr) nq
          nqtot=nq ! set value of nqtot (in infotrac module) as nq
          if (ierr.ne.0) then
            write(*,*) "testphys1d: error reading number of tracers"
            write(*,*) "   (first line of traceur.def) "
            stop
          endif
        endif
        ! allocate arrays:
        allocate(tname(nq))
        allocate(q(nlayer,nq))
        allocate(qsurf(nq))
        allocate(dq(nlayer,nq))
        allocate(dqdyn(nlayer,nq))
        allocate(mqtot(nq))
        
        ! read tracer names from file traceur.def
        do iq=1,nq
          read(90,*,iostat=ierr) tname(iq)
          if (ierr.ne.0) then
            write(*,*) 'testphys1d: error reading tracer names...'
            stop
          endif
        enddo
        close(90)

        ! initialize tracers here:
        write(*,*) "testphys1d: initializing tracers"
        q(:,:)=0 ! default, set everything to zero
        qsurf(:)=0
        ! "smarter" initialization of some tracers
        ! (get values from "profile_*" files, if these are available)
        do iq=1,nq
          txt=""
          write(txt,"(a)") tname(iq)
          write(*,*)"  tracer:",trim(txt)
          ! CO2
          if (txt.eq."co2") then
            q(:,iq)=0.95   ! kg /kg of atmosphere
            qsurf(iq)=0. ! kg/m2 (not used for CO2)
            ! even better, look for a "profile_co2" input file
            open(91,file='profile_co2',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            endif
            close(91)
          endif ! of if (txt.eq."co2")
          ! Allow for an initial profile of argon
          ! Can also be used to introduce a decaying tracer
          ! in the 1D (TBD) to study thermals
          if (txt.eq."ar") then
            !look for a "profile_ar" input file
            open(91,file='profile_ar',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ar file!"
            endif
            close(91)
          endif ! of if (txt.eq."ar")

          ! WATER VAPOUR
          if (txt.eq."h2o_vap") then
            !look for a "profile_h2o_vap" input file
            open(91,file='profile_h2o_vap',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_h2o_vap file!"
            endif
            close(91)
          endif ! of if (txt.eq."h2o_ice")
          ! WATER ICE
          if (txt.eq."h2o_ice") then
            !look for a "profile_h2o_vap" input file
            open(91,file='profile_h2o_ice',status='old',
     &       form='formatted',iostat=ierr)
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
          ! DUST
          !if (txt(1:4).eq."dust") then
          !  q(:,iq)=0.4    ! kg/kg of atmosphere
          !  qsurf(iq)=100 ! kg/m2
          !endif
          ! DUST MMR
          if (txt.eq."dust_mass") then
            !look for a "profile_dust_mass" input file
            open(91,file='profile_dust_mass',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
!                write(*,*) "l=",ilayer," q(ilayer,iq)=",q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_dust_mass file!"
            endif
            close(91)
          endif ! of if (txt.eq."dust_mass")
          ! DUST NUMBER
          if (txt.eq."dust_number") then
            !look for a "profile_dust_number" input file
            open(91,file='profile_dust_number',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_dust_number file!"
            endif
            close(91)
          endif ! of if (txt.eq."dust_number")
          ! NB: some more initializations (chemistry) is done later
          ! CCN MASS
          if (txt.eq."ccn_mass") then
            !look for a "profile_ccn_mass" input file
            open(91,file='profile_ccn_mass',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ccn_mass file!"
            endif
            close(91)
          endif ! of if (txt.eq."ccn_mass")
          ! CCN NUMBER
          if (txt.eq."ccn_number") then
            !look for a "profile_ccn_number" input file
            open(91,file='profile_ccn_number',status='old',
     &       form='formatted',iostat=ierr)
            if (ierr.eq.0) then
              read(91,*) qsurf(iq)
              do ilayer=1,nlayer
                read(91,*) q(ilayer,iq)
              enddo
            else
              write(*,*) "No profile_ccn_number file!"
            endif
            close(91)
          endif ! of if (txt.eq."ccn_number")
        enddo ! of do iq=1,nq

      else
      ! we still need to set (dummy) tracer number and names for physdem1
        nq=1
        nqtot=nq ! set value of nqtot (in infotrac module) as nq
        ! allocate arrays:
        allocate(tname(nq))
        allocate(q(nlayer,nq))
        allocate(qsurf(nq))
        allocate(dq(nlayer,nq))
        allocate(dqdyn(nlayer,nq))
        allocate(mqtot(nq))
        do iq=1,nq
          write(str7,'(a1,i2.2)')'t',iq
          tname(iq)=str7
        enddo
      ! and just to be clean, also initialize tracers to zero for physdem1
        q(:,:)=0 
        qsurf(:)=0      
      endif ! of if (tracer)
      
      !write(*,*) "testphys1d q", q(1,:)
      !write(*,*) "testphys1d qsurf", qsurf

c  Date and local time at beginning of run
c  ---------------------------------------
c    Date (in sols since spring solstice) at beginning of run
      day0 = 0 ! default value for day0
      write(*,*) 'Initial date (in martian sols ; =0 at Ls=0)?'
      call getin("day0",day0)
      day=float(day0)
      write(*,*) " day0 = ",day0
c  Local time at beginning of run 
      time=0 ! default value for time
      write(*,*)'Initial local time (in hours, between 0 and 24)?'
      call getin("time",time)
      write(*,*)" time = ",time
      time=time/24.E+0 ! convert time (hours) to fraction of sol

c  Discretization (Definition of grid and time steps)
c  --------------
c
      nlevel=nlayer+1
      nsoil=nsoilmx

      day_step=48 ! default value for day_step
      PRINT *,'Number of time steps per sol ?'
      call getin("day_step",day_step)
      write(*,*) " day_step = ",day_step

      ecritphy=day_step ! default value for ecritphy, output every time step

      ndt=10 ! default value for ndt
      PRINT *,'Number of sols to run ?'
      call getin("ndt",ndt)
      write(*,*) " ndt = ",ndt

      ndt=ndt*day_step     
      dtphys=daysec/day_step  

c Imposed surface pressure
c ------------------------------------
c
      psurf=610. ! default value for psurf
      PRINT *,'Surface pressure (Pa) ?'
      call getin("psurf",psurf)
      write(*,*) " psurf = ",psurf
c Reference pressures
      pa=20.   ! transition pressure (for hybrid coord.)
      preff=610.      ! reference surface pressure
 
c Aerosol properties
c --------------------------------
      tauvis=0.2 ! default value for tauvis (dust opacity)
      write(*,'("Reference dust opacity at ",f4.0," Pa ?")')odpref
      call getin("tauvis",tauvis)
      write(*,*) " tauvis = ",tauvis

c Orbital parameters
c ------------------
      paleomars=.false. ! Default: no water ice reservoir
      call getin("paleomars",paleomars)
      if (paleomars==.true.) then
        write(*,*) "paleomars=", paleomars
        write(*,*) "Orbital parameters from callphys.def"
        write(*,*) "Enter eccentricity & Lsperi"
        print *, 'Martian eccentricity (0<e<1) ?'
        call getin('excentric ',excentric)
        write(*,*)"excentric =",excentric 
        print *, 'Solar longitude of perihelion (0<Ls<360) ?'
        call getin('Lsperi',Lsperi )
        write(*,*)"Lsperi=",Lsperi 
        Lsperi = Lsperi*pi/180.0 ! Put it in rad for peri_day
        periheli = halfaxe*(1-excentric)
        aphelie = halfaxe*(1+excentric)
        call call_dayperi(Lsperi,excentric,peri_day,year_day)
        write(*,*) "Corresponding orbital params for GCM"
        write(*,*) " periheli = ",periheli
        write(*,*) " aphelie = ",aphelie
        write(*,*) "date of perihelion (sol)",peri_day
      else
        write(*,*) "paleomars=", paleomars
        write(*,*) "Default present-day orbital parameters"
        write(*,*) "Unless specified otherwise"
        print *,'Min. distance Sun-Mars (Mkm)?'
        call getin("periheli",periheli)
        write(*,*) " periheli = ",periheli

        print *,'Max. distance Sun-Mars (Mkm)?'
        call getin("aphelie",aphelie)
        write(*,*) " aphelie = ",aphelie

        print *,'Day of perihelion?'
        call getin("periday",peri_day)
        write(*,*) " periday = ",peri_day

        print *,'Obliquity?'
        call getin("obliquit",obliquit)
        write(*,*) " obliquit = ",obliquit
      end if
 
c  latitude/longitude
c  ------------------
      latitude(1)=0 ! default value for latitude
      PRINT *,'latitude (in degrees) ?'
      call getin("latitude",latitude(1))
      write(*,*) " latitude = ",latitude
      latitude=latitude*pi/180.E+0
      longitude=0.E+0
      longitude=longitude*pi/180.E+0

!  some initializations (some of which have already been
!  done above!) and loads parameters set in callphys.def
!  and allocates some arrays 
!Mars possible matter with dtphys in input and include!!!
! Initializations below should mimick what is done in iniphysiq for 3D GCM
      call init_physics_distribution(regular_lonlat,4,
     &                               1,1,1,nlayer,1)
      call init_interface_dyn_phys
      call init_regular_lonlat(1,1,longitude,latitude,
     &                   (/0.,0./),(/0.,0./))
      call init_geometry(1,longitude,latitude,
     &                   (/0.,0.,0.,0./),(/0.,0.,0.,0./),
     &                   cell_area)
! Ehouarn: init_vertial_layers called later (because disvert not called yet)
!      call init_vertical_layers(nlayer,preff,scaleheight,
!     &                      ap,bp,aps,bps,presnivs,pseudoalt)
      call init_dimphy(1,nlayer) ! Initialize dimphy module
      call phys_state_var_init(1,llm,nq,tname,
     .          day0,time,daysec,dtphys,rad,g,r,cpp)
      call ini_fillgeom(1,latitude,longitude,(/1.0/))
      call conf_phys(1,llm,nq)

      ! in 1D model physics are called every time step
      ! ovverride iphysiq value that has been set by conf_phys
      if (iphysiq/=1) then
        write(*,*) "testphys1d: setting iphysiq=1"
        iphysiq=1
      endif 

c  Initialize albedo / soil thermal inertia
c  ----------------------------------------
c
      albedodat(1)=0.2 ! default value for albedodat
      PRINT *,'Albedo of bare ground ?'
      call getin("albedo",albedodat(1))
      write(*,*) " albedo = ",albedodat(1)
      albedo(1,1)=albedodat(1)

      inertiedat(1,1)=400 ! default value for inertiedat
      PRINT *,'Soil thermal inertia (SI) ?'
      call getin("inertia",inertiedat(1,1))
      write(*,*) " inertia = ",inertiedat(1,1)

      z0(1)=z0_default ! default value for roughness
      write(*,*) 'Surface roughness length z0 (m)?'
      call getin("z0",z0(1))
      write(*,*) " z0 = ",z0(1)

! Initialize local slope parameters (only matters if "callslope"
! is .true. in callphys.def)
      ! slope inclination angle (deg) 0: horizontal, 90: vertical
      theta_sl(1)=0.0 ! default: no inclination
      call getin("slope_inclination",theta_sl(1))
      ! slope orientation (deg)
      ! 0 == Northward, 90 == Eastward, 180 == Southward, 270 == Westward
      psi_sl(1)=0.0 ! default value
      call getin("slope_orientation",psi_sl(1))
      
c
c  for the gravity wave scheme
c  ---------------------------------
c
      zmea(1)=0.E+0
      zstd(1)=0.E+0
      zsig(1)=0.E+0
      zgam(1)=0.E+0
      zthe(1)=0.E+0
c
c  for the slope wind scheme
c  ---------------------------------
c  
      hmons(1)=0.E+0
      PRINT *,'hmons is initialized to ',hmons(1)
      summit(1)=0.E+0
      PRINT *,'summit is initialized to ',summit(1)
      base(1)=0.E+0
c
c  Default values initializing the coefficients calculated later
c  ---------------------------------
c  
      tauscaling(1)=1. ! calculated in aeropacity_mod.F
      totcloudfrac(1)=1. ! calculated in watercloud_mod.F      

c   Specific initializations for "physiq"
c   -------------------------------------
c   surface geopotential is not used (or useful) since in 1D
c   everything is controled by surface pressure
      phisfi(1)=0.E+0

c   Initialization to take into account prescribed winds
c   ------------------------------------------------------
      ptif=2.E+0*omeg*sinlat(1)
 
c    geostrophic wind
      gru=10. ! default value for gru
      PRINT *,'zonal eastward component of the geostrophic wind (m/s) ?'
      call getin("u",gru)
      write(*,*) " u = ",gru
      grv=0. !default value for grv
      PRINT *,'meridional northward component of the geostrophic',
     &' wind (m/s) ?'
      call getin("v",grv)
      write(*,*) " v = ",grv

c     Initialize winds  for first time step
      DO ilayer=1,nlayer
         u(ilayer)=gru
         v(ilayer)=grv
         w(ilayer)=0 ! default: no vertical wind
      ENDDO

c     Initialize turbulente kinetic energy
      DO ilevel=1,nlevel
         q2(ilevel)=0.E+0
      ENDDO

c  CO2 ice on the surface
c  -------------------
      co2ice(1)=0.E+0 ! default value for co2ice
      PRINT *,'Initial CO2 ice on the surface (kg.m-2)'
      call getin("co2ice",co2ice)
      write(*,*) " co2ice = ",co2ice
! Initialization for CO2 clouds (could be improved to read initial profiles)
      mem_Mccn_co2(:,:)=0
      mem_Mh2o_co2(:,:)=0
      mem_Nccn_co2(:,:)=0
c
c  emissivity
c  ----------
      emis=emissiv
      IF (co2ice(1).eq.1.E+0) THEN
         emis=emisice(1) ! northern hemisphere
         IF(latitude(1).LT.0) emis=emisice(2) ! southern hemisphere
      ENDIF

 

c  Compute pressures and altitudes of atmospheric levels 
c  ----------------------------------------------------------------

c    Vertical Coordinates
c    """"""""""""""""""""
      hybrid=.true.
      PRINT *,'Hybrid coordinates ?'
      call getin("hybrid",hybrid)
      write(*,*) " hybrid = ", hybrid

      CALL  disvert_noterre
      ! now that disvert has been called, initialize module vertical_layers_mod
      call init_vertical_layers(nlayer,preff,scaleheight,
     &                      ap,bp,aps,bps,presnivs,pseudoalt)

      DO ilevel=1,nlevel
        plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
      ENDDO

      DO ilayer=1,nlayer
        play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
      ENDDO

      DO ilayer=1,nlayer
        zlay(ilayer)=-200.E+0 *r*log(play(ilayer)/plev(1))
     &   /g
      ENDDO


c  Initialize temperature profile
c  --------------------------------------
      pks=psurf**rcp

c altitude in km in profile: divide zlay by 1000
      tmp1(0)=0.E+0
      DO ilayer=1,nlayer
        tmp1(ilayer)=zlay(ilayer)/1000.E+0
      ENDDO

      call profile(nlayer+1,tmp1,tmp2)

      tsurf=tmp2(0)
      DO ilayer=1,nlayer
        temp(ilayer)=tmp2(ilayer)
      ENDDO
      


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
        mlayer(isoil)=2.e-4*(2.**(isoil-0.5)) ! mid-layer depth
      enddo
      do isoil=1,nsoil
        layer(isoil)=2.e-4*(2.**(isoil-1)) ! layer depth
      enddo

c    Initialize traceurs
c    ---------------------------

      if (photochem.or.callthermos) then
         write(*,*) 'Initializing chemical species'
         ! thermo=0: initialize over all atmospheric layers
         thermo=0
         call inichim_newstart(q,psurf,sig,nq,latitude,longitude,
     $   cell_area,thermo,qsurf)
      endif

c Check if the surface is a water ice reservoir 
c --------------------------------------------------
      watercaptag(1)=.false. ! Default: no water ice reservoir
      print *,'Water ice cap on ground ?'
      call getin("watercaptag",watercaptag)
      write(*,*) " watercaptag = ",watercaptag
      

c    Initialization for GRADS outputs in "g1d.dat" and "g1d.ctl"
c    ----------------------------------------------------------------
c    (output done in "writeg1d", typically called by "physiq.F")

        g1d_nlayer=nlayer
        g1d_nomfich='g1d.dat'
        g1d_unitfich=40
        g1d_nomctl='g1d.ctl'
        g1d_unitctl=41
        g1d_premier=.true.
        g2d_premier=.true.

c  Write a "startfi" file
c  --------------------
c  This file will be read during the first call to "physiq".
c  It is needed to transfert physics variables to "physiq"...

      call physdem0("startfi.nc",longitude,latitude,nsoilmx,ngrid,llm,
     &              nq,dtphys,float(day0),0.,cell_area,
     &              albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe,
     &              hmons,summit,base)
      call physdem1("startfi.nc",nsoilmx,ngrid,llm,nq,
     &              dtphys,time,
     &              tsurf,tsoil,co2ice,albedo,emis,q2,qsurf,tauscaling,
     &              totcloudfrac,wstar,
     &              mem_Mccn_co2,mem_Nccn_co2,
     &              mem_Mh2o_co2)

c=======================================================================
c  1D MODEL TIME STEPPING LOOP
c=======================================================================
c
      firstcall=.true.
      lastcall=.false.

      DO idt=1,ndt
c        IF (idt.eq.ndt) lastcall=.true.
        IF (idt.eq.ndt-day_step-1) then       !test
         lastcall=.true.
         call solarlong(day*1.0,zls)
         write(103,*) 'Ls=',zls*180./pi
         write(103,*) 'Lat=', latitude(1)*180./pi
         write(103,*) 'Tau=', tauvis/odpref*psurf
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(103,*) 'RunEnd - Atmos. Temp. File'
         write(104,*) 'Ls=',zls*180./pi
         write(104,*) 'Lat=', latitude(1)
         write(104,*) 'Tau=', tauvis/odpref*psurf
         write(104,*) 'RunEnd - Atmos. Temp. File'
        ENDIF

c     compute geopotential
c     ~~~~~~~~~~~~~~~~~~~~~
      DO ilayer=1,nlayer
        s(ilayer)=(aps(ilayer)/psurf+bps(ilayer))**rcp
        h(ilayer)=cpp*temp(ilayer)/(pks*s(ilayer))
      ENDDO
      phi(1)=pks*h(1)*(1.E+0-s(1))
      DO ilayer=2,nlayer
         phi(ilayer)=phi(ilayer-1)+
     &               pks*(h(ilayer-1)+h(ilayer))*.5E+0
     &                  *(s(ilayer-1)-s(ilayer))

      ENDDO

c       call physics
c       --------------------
!      write(*,*) "testphys1d avant q", q(1,:)
      CALL physiq (1,llm,nq,
     ,     firstcall,lastcall,
     ,     day,time,dtphys,
     ,     plev,play,phi,
     ,     u, v,temp, q,  
     ,     w,
C - outputs
     s     du, dv, dtemp, dq,dpsurf)
!      write(*,*) "testphys1d apres q", q(1,:)


c       wind increment : specific for 1D
c       --------------------------------
 
c       The physics compute the tendencies on u and v,
c       here we just add Coriolos effect
c
c       DO ilayer=1,nlayer
c          du(ilayer)=du(ilayer)+ptif*(v(ilayer)-grv)
c          dv(ilayer)=dv(ilayer)+ptif*(-u(ilayer)+gru)
c       ENDDO

c       For some tests : No coriolis force at equator
c       if(latitude(1).eq.0.) then
          DO ilayer=1,nlayer
             du(ilayer)=du(ilayer)+ (gru-u(ilayer))/1.e4
             dv(ilayer)=dv(ilayer)+ (grv-v(ilayer))/1.e4
          ENDDO
c       end if
c      
c
c       Compute time for next time step
c       ---------------------------------------
        firstcall=.false.
        time=time+dtphys/daysec
        IF (time.gt.1.E+0) then
            time=time-1.E+0
            day=day+1
        ENDIF

c       compute winds and temperature for next time step
c       ----------------------------------------------------------

        DO ilayer=1,nlayer
           u(ilayer)=u(ilayer)+dtphys*du(ilayer)
           v(ilayer)=v(ilayer)+dtphys*dv(ilayer)
           temp(ilayer)=temp(ilayer)+dtphys*dtemp(ilayer)
        ENDDO

c       compute pressure for next time step 
c       ----------------------------------------------------------

           psurf=psurf+dtphys*dpsurf(1)   ! surface pressure change
           DO ilevel=1,nlevel
             plev(ilevel)=ap(ilevel)+psurf*bp(ilevel)
           ENDDO
           DO ilayer=1,nlayer
             play(ilayer)=aps(ilayer)+psurf*bps(ilayer)
           ENDDO

!       increment tracers
        DO iq = 1, nq
          DO ilayer=1,nlayer
             q(ilayer,iq)=q(ilayer,iq)+dtphys*dq(ilayer,iq)
          ENDDO
        ENDDO

      ENDDO   ! of idt=1,ndt ! end of time stepping loop

c    ========================================================
c    OUTPUTS
c    ========================================================

c    finalize and close grads files "g1d.dat" and "g1d.ctl"

c        CALL endg1d(1,nlayer,zphi/(g*1000.),ndt)
        CALL endg1d(1,nlayer,zlay/1000.,ndt)

      write(*,*) "testphys1d: Everything is cool."

      END
 
c***********************************************************************
c***********************************************************************
c     Dummy subroutines used only in 3D, but required to
c     compile testphys1d (to cleanly use writediagfi)

      subroutine gr_fi_dyn(nfield,ngrid,im,jm,pfi,pdyn)

      IMPLICIT NONE

      INTEGER im,jm,ngrid,nfield
      REAL pdyn(im,jm,nfield)
      REAL pfi(ngrid,nfield)
      
      if (ngrid.ne.1) then
        write(*,*) "gr_fi_dyn error: in 1D ngrid should be 1!!!"
        stop
      endif
      
      pdyn(1,1,1:nfield)=pfi(1,1:nfield)
      
      end
 
c***********************************************************************
c***********************************************************************

