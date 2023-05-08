MODULE inifis_mod
IMPLICIT NONE

CONTAINS

  SUBROUTINE inifis(ngrid,nlayer,nq, &
             day_ini,pdaysec,nday,ptimestep, &
             plat,plon,parea, &
             prad,pg,pr,pcpp)

  use init_print_control_mod, only: init_print_control
  use radinc_h, only: ini_radinc_h
  use datafile_mod
  use comdiurn_h, only: sinlat, coslat, sinlon, coslon
  use comgeomfi_h, only: totarea, totarea_planet
  use comsoil_h, only: ini_comsoil_h, nsoilmx, lay1_soil, alpha_soil
  use time_phylmdz_mod, only: ecritphy,day_step,iphysiq, &
                              init_time, daysec, dtphys
  use comcstfi_mod, only: rad, cpp, g, r, rcp, &
                          mugaz, pi, avocado, kbol
  use planete_mod, only: nres
  use planetwide_mod, only: planetwide_sumval
  use callkeys_mod
  use mod_phys_lmdz_para, only : is_parallel

!=======================================================================
!
!   purpose:
!   -------
!
!   Initialisation for the physical parametrisations of the LMD 
!   Generic Model.
!
!   author: Frederic Hourdin 15 / 10 /93
!   -------
!   modified: Sebastien Lebonnois 11/06/2003 (new callphys.def)
!             Ehouarn Millour (oct. 2008) tracers are now identified
!              by their names and may not be contiguously
!              stored in the q(:,:,:,:) array
!             E.M. (june 2009) use getin routine to load parameters
!
!
!   arguments:
!   ----------
!
!   input:
!   ------
!
!    ngrid                 Size of the horizontal grid.
!                          All internal loops are performed on that grid.
!    nlayer                Number of vertical layers.
!    pdayref               Day of reference for the simulation
!    pday                  Number of days counted from the North. Spring
!                          equinoxe.
!
!=======================================================================
!
!-----------------------------------------------------------------------
!   declarations:
!   -------------
  use ioipsl_getin_p_mod, only: getin_p
  IMPLICIT NONE



  REAL,INTENT(IN) :: prad,pg,pr,pcpp,pdaysec,ptimestep
  INTEGER,INTENT(IN) :: nday
  INTEGER,INTENT(IN) :: ngrid,nlayer,nq
  REAL,INTENT(IN) :: plat(ngrid),plon(ngrid),parea(ngrid)
  integer,intent(in) :: day_ini
  INTEGER ig,ierr
 
  EXTERNAL iniorbit,orbite
  EXTERNAL SSUM
  REAL SSUM
 
  ! Initialize flags lunout, prt_level, debug (in print_control_mod)
  CALL init_print_control

  ! initialize constants in comcstfi_mod
  rad=prad
  cpp=pcpp
  g=pg
  r=pr
  rcp=r/cpp
  mugaz=8.314*1000./pr ! dummy init
  pi=2.*asin(1.)
  avocado = 6.02214179e23   ! added by RW
  kbol = 1.38064852e-23  ! added by JVO for Titan chem

  ! Initialize some "temporal and calendar" related variables
  CALL init_time(day_ini,pdaysec,nday,ptimestep)

  ! read in some parameters from "run.def" for physics,
  ! or shared between dynamics and physics.
  call getin_p("ecritphy",ecritphy) ! frequency of outputs in physics,
                                    ! in dynamical steps
  call getin_p("day_step",day_step) ! number of dynamical steps per day
  call getin_p("iphysiq",iphysiq) ! call physics every iphysiq dyn step

  ! do we read a startphy.nc file? (default: .true.)
  call getin_p("startphy_file",startphy_file)
  
! --------------------------------------------------------------
!  Reading the "callphys.def" file controlling some key options
! --------------------------------------------------------------
     
  ! check that 'callphys.def' file is around
  OPEN(99,file='callphys.def',status='old',form='formatted',iostat=ierr)
  CLOSE(99)
  IF(ierr.EQ.0) iscallphys=.true. !iscallphys initialised as false in callkeys_mod module
      
!!!      IF(ierr.EQ.0) THEN
  IF(iscallphys) THEN
     PRINT*
     PRINT*
     PRINT*,'--------------------------------------------'
     PRINT*,' inifis: Parametres pour la physique (callphys.def)'
     PRINT*,'--------------------------------------------'

     write(*,*) "Directory where external input files are:"
     ! default 'datadir' is set in "datadir_mod"
     call getin_p("datadir",datadir) ! default path
     write(*,*) " datadir = ",trim(datadir)

     write(*,*) "Run with or without tracer transport ?"
     tracer=.false. ! default value
     call getin_p("tracer",tracer)
     write(*,*) " tracer = ",tracer

     write(*,*) "Run with or without atm mass update ", &
            " due to tracer evaporation/condensation?"
     mass_redistrib=.false. ! default value
     call getin_p("mass_redistrib",mass_redistrib)
     write(*,*) " mass_redistrib = ",mass_redistrib

     write(*,*) "Diurnal cycle ?"
     write(*,*) "(if diurnal=false, diurnal averaged solar heating)"
     diurnal=.true. ! default value
     call getin_p("diurnal",diurnal)
     write(*,*) " diurnal = ",diurnal

     write(*,*) "Seasonal cycle ?"
     write(*,*) "(if season=false, Ls stays constant, to value ", &
         "set in 'start'"
     season=.true. ! default value
     call getin_p("season",season)
     write(*,*) " season = ",season

     write(*,*) "Tidally resonant rotation ?"
     tlocked=.false. ! default value
     call getin_p("tlocked",tlocked)
     write(*,*) "tlocked = ",tlocked

     write(*,*) "Saturn ring shadowing ?"
     rings_shadow = .false.
     call getin_p("rings_shadow", rings_shadow)
     write(*,*) "rings_shadow = ", rings_shadow
         
     write(*,*) "Compute latitude-dependent gravity field?"
     oblate = .false.
     call getin_p("oblate", oblate)
     write(*,*) "oblate = ", oblate

     write(*,*) "Flattening of the planet (a-b)/a "
     flatten = 0.0
     call getin_p("flatten", flatten)
     write(*,*) "flatten = ", flatten         

     write(*,*) "Needed if oblate=.true.: J2"
     J2 = 0.0
     call getin_p("J2", J2)
     write(*,*) "J2 = ", J2
         
     write(*,*) "Needed if oblate=.true.: Planet mass (*1e24 kg)"
     MassPlanet = 0.0
     call getin_p("MassPlanet", MassPlanet)
     write(*,*) "MassPlanet = ", MassPlanet         

     write(*,*) "Needed if oblate=.true.: Planet mean radius (m)"
     Rmean = 0.0
     call getin_p("Rmean", Rmean)
     write(*,*) "Rmean = ", Rmean
     
     write(*,*) "Compute effective altitude-dependent gravity field?"
     eff_gz = .false.
     call getin_p("eff_gz", eff_gz)
     write(*,*) "eff_gz = ", eff_gz
     
     ! sanity check warning
     if (eff_gz) then
       print*,"WARNING : You run chemistry with effective altitude-dependent gravity field !!"
       print*,"You will have no coherence in your heating rates between physics and dynamics !!"
       print*,"I let you continue but you should rather set eff_gz =.false. ..."
     endif

         
! Test of incompatibility:
! if tlocked, then diurnal should be false
     if (tlocked.and.diurnal) then
       print*,'If diurnal=true, we should turn off tlocked.'
       stop
     endif

     write(*,*) "Tidal resonance ratio ?"
     nres=0          ! default value
     call getin_p("nres",nres)
     write(*,*) "nres = ",nres

     write(*,*) "Write some extra output to the screen ?"
     lwrite=.false. ! default value
     call getin_p("lwrite",lwrite)
     write(*,*) " lwrite = ",lwrite

     write(*,*) "Save statistics in file stats.nc ?"
     callstats=.true. ! default value
     call getin_p("callstats",callstats)
     write(*,*) " callstats = ",callstats

     write(*,*) "Test energy conservation of model physics ?"
     enertest=.false. ! default value
     call getin_p("enertest",enertest)
     write(*,*) " enertest = ",enertest

     write(*,*) "Check to see if cpp values used match gases.def ?"
     check_cpp_match=.true. ! default value
     call getin_p("check_cpp_match",check_cpp_match)
     write(*,*) " check_cpp_match = ",check_cpp_match

     write(*,*) "call radiative transfer ?"
     callrad=.true. ! default value
     call getin_p("callrad",callrad)
     write(*,*) " callrad = ",callrad

     write(*,*) "call correlated-k radiative transfer ?"
     corrk=.true. ! default value
     call getin_p("corrk",corrk)
     write(*,*) " corrk = ",corrk
     
     if (corrk) then
       ! default path is set in datadir
       write(*,*) "callcorrk: Correlated-k data base folder:",trim(datadir)
       call getin_p("corrkdir",corrkdir)
       write(*,*) " corrkdir = ",corrkdir
       
       write(*,*) "use correlated-k recombination instead of pre-mixed values ?"
       corrk_recombin=.false.! default value
       call getin_p("corrk_recombin",corrk_recombin)
       write(*,*) " corrk_recombin = ",corrk_recombin
     endif
     
     if (corrk .and. ngrid.eq.1) then
       write(*,*) "simulate global averaged conditions ?"
       global1d = .false. ! default value
       call getin_p("global1d",global1d)
       write(*,*) " global1d = ",global1d
       
       ! Test of incompatibility : if global1d is true, there should not be any diurnal cycle.
       if (global1d.and.diurnal) then
          write(*,*) "if global1d is true, diurnal must be set to false"
          stop
       endif

       if (global1d) then
          write(*,*) "Solar Zenith angle (deg.) ?"
          write(*,*) "(assumed for averaged solar flux S/4)"
          szangle=60.0  ! default value
          call getin_p("szangle",szangle)
          write(*,*) " szangle = ",szangle
       endif
     endif

     write(*,*) "prohibit calculations outside corrk T grid?"
     strictboundcorrk=.true. ! default value
     call getin_p("strictboundcorrk",strictboundcorrk)
     write(*,*) "strictboundcorrk = ",strictboundcorrk

     write(*,*) "call gaseous absorption in the visible bands?", &
                    "(matters only if callrad=T)"
     callgasvis=.false. ! default value
     call getin_p("callgasvis",callgasvis)
     write(*,*) " callgasvis = ",callgasvis
        
     write(*,*) "call continuum opacities in radiative transfer ?", &
                    "(matters only if callrad=T)"
     continuum=.true. ! default value
     call getin_p("continuum",continuum)
     write(*,*) " continuum = ",continuum
 
     write(*,*) "version for H2H2 CIA file ?"
     versH2H2cia=2011 ! default value (should be 2018 but retrocompatibility first)
     call getin_p("versH2H2cia",versH2H2cia)
     write(*,*) " versH2H2cia = ",versH2H2cia
     ! Sanity check
     if (versH2H2cia.ne.2011 .and. versH2H2cia.ne.2018) then
        print*,'Error: Choose a correct value (2011 or 2018) for versH2H2cia !'
        call abort
     endif 

     write(*,*) "call turbulent vertical diffusion ?"
     calldifv=.true. ! default value
     call getin_p("calldifv",calldifv)
     write(*,*) " calldifv = ",calldifv

     write(*,*) "use turbdiff instead of vdifc ?"
     UseTurbDiff=.true. ! default value
     call getin_p("UseTurbDiff",UseTurbDiff)
     write(*,*) " UseTurbDiff = ",UseTurbDiff

     write(*,*) "call convective adjustment ?"
     calladj=.true. ! default value
     call getin_p("calladj",calladj)
     write(*,*) " calladj = ",calladj

     write(*,*) "Radiative timescale for Newtonian cooling ?"
     tau_relax=30. ! default value
     call getin_p("tau_relax",tau_relax)
     write(*,*) " tau_relax = ",tau_relax
     tau_relax=tau_relax*24*3600 ! convert Earth days --> seconds

     write(*,*)"call thermal conduction in the soil ?"
     callsoil=.true. ! default value
     call getin_p("callsoil",callsoil)
     write(*,*) " callsoil = ",callsoil
         
     write(*,*)"Rad transfer is computed every iradia", &
                   " physical timestep"
     iradia=1 ! default value
     call getin_p("iradia",iradia)
     write(*,*)" iradia = ",iradia
       
     write(*,*)"Rayleigh scattering ?"
     rayleigh=.false.
     call getin_p("rayleigh",rayleigh)
     write(*,*)" rayleigh = ",rayleigh

     write(*,*) "Use blackbody for stellar spectrum ?"
     stelbbody=.false. ! default value
     call getin_p("stelbbody",stelbbody)
     write(*,*) " stelbbody = ",stelbbody

     write(*,*) "Stellar blackbody temperature ?"
     stelTbb=5800.0 ! default value
     call getin_p("stelTbb",stelTbb)
     write(*,*) " stelTbb = ",stelTbb

     write(*,*)"Output mean OLR in 1D?"
     meanOLR=.false.
     call getin_p("meanOLR",meanOLR)
     write(*,*)" meanOLR = ",meanOLR

     write(*,*)"Output spectral OLR in 3D?"
     specOLR=.false.
     call getin_p("specOLR",specOLR)
     write(*,*)" specOLR = ",specOLR

     write(*,*)"Output diagnostic optical thickness attenuation exp(-tau)?"
     diagdtau=.false.
     call getin_p("diagdtau",diagdtau)
     write(*,*)" diagdtau = ",diagdtau

     write(*,*)"Uniform absorption in radiative transfer?"
     graybody=.false.
     call getin_p("graybody",graybody)
     write(*,*)" graybody = ",graybody

! Chemistry

     write(*,*) "Run with or without chemistry?"
     callchim=.false. ! default value
     call getin_p("callchim",callchim)
     write(*,*) " callchim = ",callchim

     ! sanity check
     if (callchim.and.(.not.tracer)) then
       print*,"You are running chemistry without tracer"
       print*,"Please start again with tracer =.true."
       stop
     endif
     
     write(*,*)"Chemistry is computed every ichim", &
                   " physical timestep"
     ichim=1 ! default value
     call getin_p("ichim",ichim)
     write(*,*)" ichim = ",ichim

! Microphysics

     write(*,*) "Use haze seasonal model from Karkoschka 2016 ?"
     seashaze=.true. ! default value
     call getin_p("seashaze",seashaze)
     write(*,*)" seashaze = ",seashaze

     write(*,*) "Run with or without microphysics?"
     callmufi=.false. ! default value
     call getin_p("callmufi",callmufi)
     write(*,*)" callmufi = ",callmufi

     ! sanity check
     if (callmufi.and.(.not.tracer)) then
       print*,"You are running microphysics without tracer"
       print*,"Please start again with tracer =.true."
       stop
     endif

     write(*,*) "Compute clouds?"
     callclouds=.false. ! default value
     call getin_p("callclouds",callclouds)
     write(*,*)" callclouds = ",callclouds

     ! sanity check
     if (callclouds.and.(.not.callmufi)) then
       print*,"You are trying to make clouds without microphysics !"
       print*,"Please start again with callmufi =.true."
       stop
     endif

     write(*,*) "Disable the coupling of microphysics within rad. transf. ?"
     write(*,*) "If disabled we will assume a planetwide vert. profile of extinction ..."
     uncoupl_optic_haze=.true. ! default value - true as long as the microphysics is bugged
     call getin_p("uncoupl_optic_haze",uncoupl_optic_haze)
     write(*,*)" uncoupl_optic_haze = ",uncoupl_optic_haze

     write(*,*) "Pressure level of aer. production (Pa) ?"
     p_prod=1.0 ! default value
     call getin_p("p_prod",p_prod)
     write(*,*)" p_prod = ",p_prod
     
     write(*,*) "Aerosol production rate (kg.m-2.s-1) ?"
     tx_prod=3.5e-13 ! default value
     call getin_p("tx_prod",tx_prod)
     write(*,*)" tx_prod = ",tx_prod

     write(*,*) "Equivalent radius production (m) ?"
     rc_prod=2.0e-8 ! default value
     call getin_p("rc_prod",rc_prod)
     write(*,*)" rhc_prod = ",rc_prod

     write(*,*) "Radius of air (nitrogen) molecule (m) ?"
     air_rad=1.75e-10 ! default value
     call getin_p("air_rad",air_rad)
     write(*,*)" air_rad = ",air_rad

     write(*,*) "Path to microphys. config file ?"
     config_mufi='datagcm/microphysics/config.cfg' ! default value
     call getin_p("config_mufi",config_mufi)
     write(*,*)" config_mufi = ",config_mufi

! Soil model
     write(*,*)"Number of sub-surface layers for soil scheme?"
     ! default value of nsoilmx set in comsoil_h
     call getin_p("nsoilmx",nsoilmx)
     write(*,*)" nsoilmx=",nsoilmx
     
     write(*,*)"Thickness of topmost soil layer (m)?"
     ! default value of lay1_soil set in comsoil_h
     call getin_p("lay1_soil",lay1_soil)
     write(*,*)" lay1_soil = ",lay1_soil
     
     write(*,*)"Coefficient for soil layer thickness distribution?"
     ! default value of alpha_soil set in comsoil_h
     call getin_p("alpha_soil",alpha_soil)
     write(*,*)" alpha_soil = ",alpha_soil

     write(*,*)"Remove lower boundary?"
     nosurf=.false.
     call getin_p("nosurf",nosurf)
     write(*,*)" nosurf = ",nosurf

! Tests of incompatibility:
     if (nosurf.and.callsoil) then
       print*,'nosurf not compatible with soil scheme!'
       print*,'... got to make a choice!'
       call abort
     endif

     write(*,*)"Add an internal heat flux?", &
                   "... matters only if callsoil=F"
     intheat=0.
     call getin_p("intheat",intheat)
     write(*,*)" intheat = ",intheat

     write(*,*)"Use Newtonian cooling for radiative transfer?"
     newtonian=.false.
     call getin_p("newtonian",newtonian)
     write(*,*)" newtonian = ",newtonian

! Tests of incompatibility:
     if (newtonian.and.corrk) then
       print*,'newtonian not compatible with correlated-k!'
       call abort
     endif
     if (newtonian.and.calladj) then
       print*,'newtonian not compatible with adjustment!'
       call abort
     endif
     if (newtonian.and.calldifv) then
       print*,'newtonian not compatible with a boundary layer!'
       call abort
     endif

     write(*,*)"Test physics timescale in 1D?"
     testradtimes=.false.
     call getin_p("testradtimes",testradtimes)
     write(*,*)" testradtimes = ",testradtimes

! Test of incompatibility:
! if testradtimes used, we must be in 1D
     if (testradtimes.and.(ngrid.gt.1)) then
       print*,'testradtimes can only be used in 1D!'
       call abort
     endif

     write(*,*)"Which star?"
     startype=1 ! default value = Sol
     call getin_p("startype",startype)
     write(*,*)" startype = ",startype

     write(*,*)"Value of stellar flux at 1 AU?"
     Fat1AU=1356.0 ! default value = Sol today
     call getin_p("Fat1AU",Fat1AU)
     write(*,*)" Fat1AU = ",Fat1AU

     write(*,*) "Does user want to force cpp and mugaz?"
     force_cpp=.false. ! default value
     call getin_p("force_cpp",force_cpp)
     write(*,*) " force_cpp = ",force_cpp

     if (force_cpp) then
       mugaz = -99999.
       PRINT *,'MEAN MOLECULAR MASS in g mol-1 ?'
       call getin_p("mugaz",mugaz)
       IF (mugaz.eq.-99999.) THEN
           PRINT *, "mugaz must be set if force_cpp = T"
           STOP
       ELSE
           write(*,*) "mugaz=",mugaz
       ENDIF
       cpp = -99999.
       PRINT *,'SPECIFIC HEAT CAPACITY in J K-1 kg-1 ?'
       call getin_p("cpp",cpp)
       IF (cpp.eq.-99999.) THEN
           PRINT *, "cpp must be set if force_cpp = T"
           STOP
       ELSE
           write(*,*) "cpp=",cpp
       ENDIF
     endif ! of if (force_cpp)
     
     
     call su_gases(nlayer,tracer)     
     call calc_cpp_mugaz
     
     
     PRINT*,'--------------------------------------------'
     PRINT*
     PRINT*
  ELSE
     write(*,*)
     write(*,*) 'Cannot read file callphys.def. Is it here ?'
     stop
  ENDIF ! of IF(iscallphys)

  PRINT*
  PRINT*,'inifis: daysec',daysec
  PRINT*
  PRINT*,'inifis: The radiative transfer is computed:'
  PRINT*,'           each ',iradia,' physical time-step'
  PRINT*,'        or each ',iradia*dtphys,' seconds'
  PRINT*
  PRINT*,'inifis: Chemistry is computed:'
  PRINT*,'           each ',ichim,' physical time-step'
  PRINT*,'        or each ',ichim*dtphys,' seconds'
  PRINT*

!-----------------------------------------------------------------------
!     Some more initialization:
!     ------------------------

  ! Initializations for comgeomfi_h
  totarea=SSUM(ngrid,parea,1)
  call planetwide_sumval(parea,totarea_planet)

  !! those are defined in comdiurn_h.F90
  IF (.not.ALLOCATED(sinlat)) ALLOCATE(sinlat(ngrid))
  IF (.not.ALLOCATED(coslat)) ALLOCATE(coslat(ngrid))
  IF (.not.ALLOCATED(sinlon)) ALLOCATE(sinlon(ngrid))
  IF (.not.ALLOCATED(coslon)) ALLOCATE(coslon(ngrid))

  DO ig=1,ngrid
     sinlat(ig)=sin(plat(ig))
     coslat(ig)=cos(plat(ig))
     sinlon(ig)=sin(plon(ig))
     coslon(ig)=cos(plon(ig))
  ENDDO

  ! initialize variables in radinc_h
  call ini_radinc_h(nlayer)
 
  ! allocate "comsoil_h" arrays
  call ini_comsoil_h(ngrid)
      
  END SUBROUTINE inifis

END MODULE inifis_mod
