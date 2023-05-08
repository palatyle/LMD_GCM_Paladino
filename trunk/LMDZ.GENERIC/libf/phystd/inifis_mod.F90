MODULE inifis_mod
IMPLICIT NONE

CONTAINS

  SUBROUTINE inifis(ngrid,nlayer,nq, &
             day_ini,pdaysec,nday,ptimestep, &
             plat,plon,parea, &
             prad,pg,pr,pcpp)

  use init_print_control_mod, only: init_print_control
  use radinc_h, only: ini_radinc_h, naerkind
  use datafile_mod, only: datadir
  use comdiurn_h, only: sinlat, coslat, sinlon, coslon
  use comgeomfi_h, only: totarea, totarea_planet
  use comsoil_h, only: ini_comsoil_h, nsoilmx, lay1_soil, alpha_soil
  use time_phylmdz_mod, only: ecritphy,day_step,iphysiq, &
                              init_time, daysec, dtphys
  use comcstfi_mod, only: rad, cpp, g, r, rcp, &
                          mugaz, pi, avocado
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
  use datafile_mod, only: datadir
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

  ! Initialize some "temporal and calendar" related variables
#ifndef MESOSCALE
  CALL init_time(day_ini,pdaysec,nday,ptimestep)
#endif

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

     write(*,*) "prohibit calculations outside corrk T grid?"
     strictboundcorrk=.true. ! default value
     call getin_p("strictboundcorrk",strictboundcorrk)
     write(*,*) "strictboundcorrk = ",strictboundcorrk

     write(*,*) "Minimum atmospheric temperature for Planck function integration ?"
     tplanckmin=30.0 ! default value
     call getin_p("tplanckmin",tplanckmin)
     write(*,*) " tplanckmin = ",tplanckmin
 
     write(*,*) "Maximum atmospheric temperature for Planck function integration ?"
     tplanckmax=1500.0 ! default value
     call getin_p("tplanckmax",tplanckmax)
     write(*,*) " tplanckmax = ",tplanckmax
 
     write(*,*) "Temperature step for Planck function integration ?"
     dtplanck=0.1 ! default value
     call getin_p("dtplanck",dtplanck)
     write(*,*) " dtplanck = ",dtplanck
 
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

     write(*,*) "use analytic function for H2O continuum ?"
     H2Ocont_simple=.false. ! default value
     call getin_p("H2Ocont_simple",H2Ocont_simple)
     write(*,*) " H2Ocont_simple = ",H2Ocont_simple
 
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

     write(*,*) "call thermal plume model ?"
     calltherm=.false. ! default value
     call getin_p("calltherm",calltherm)
     write(*,*) " calltherm = ",calltherm

     write(*,*) "call CO2 condensation ?"
     co2cond=.false. ! default value
     call getin_p("co2cond",co2cond)
     write(*,*) " co2cond = ",co2cond
! Test of incompatibility
     if (co2cond.and.(.not.tracer)) then
        print*,'We need a CO2 ice tracer to condense CO2'
        call abort
     endif 
 
     write(*,*) "CO2 supersaturation level ?"
     co2supsat=1.0 ! default value
     call getin_p("co2supsat",co2supsat)
     write(*,*) " co2supsat = ",co2supsat

     write(*,*) "Radiative timescale for Newtonian cooling ?"
     tau_relax=30. ! default value
     call getin_p("tau_relax",tau_relax)
     write(*,*) " tau_relax = ",tau_relax
     tau_relax=tau_relax*24*3600 ! convert Earth days --> seconds

     write(*,*)"call thermal conduction in the soil ?"
     callsoil=.true. ! default value
     call getin_p("callsoil",callsoil)
     write(*,*) " callsoil = ",callsoil
     
     write(*,*)"Erupt volcano ?"
     callvolcano=.true. ! default value
     call getin_p("callvolcano",callvolcano)
     if (callvolcano) then ! Paladino
      call getin_p("volc_name",volc_name)
      ! Warning: There is no case for when input_key is empty
      call getin_p("input_key",input_key)
     endif 
     write(*,*) " callvolcano = ",callvolcano
    
     write(*,*)"output directory:" ! Paladino
     output_dir = "diagfi.nc" ! Default to where executable exists
     call getin_p("output_dir",output_dir)
     write(*,*) output_dir

     write(*,*)"Atmosphere Type?"
     atmos_type = "cd"
     call getin_p("atmos_type",atmos_type)
     write(*,*) atmos_type

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

     write(*,*)"Operate in kastprof mode?"
     kastprof=.false.
     call getin_p("kastprof",kastprof)
     write(*,*)" kastprof = ",kastprof

     write(*,*)"Uniform absorption in radiative transfer?"
     graybody=.false.
     call getin_p("graybody",graybody)
     write(*,*)" graybody = ",graybody

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

! Slab Ocean 
     write(*,*) "Use slab-ocean ?"
     ok_slab_ocean=.false.         ! default value
     call getin_p("ok_slab_ocean",ok_slab_ocean)
     write(*,*) "ok_slab_ocean = ",ok_slab_ocean
     ! Sanity check: for now slab oncean only works in serial mode
     if (ok_slab_ocean.and.is_parallel) then
       write(*,*) " Error: slab ocean should only be used in serial mode!"
       call abort
     endif

     write(*,*) "Use slab-sea-ice ?"
     ok_slab_sic=.true.         ! default value
     call getin_p("ok_slab_sic",ok_slab_sic)
     write(*,*) "ok_slab_sic = ",ok_slab_sic

     write(*,*) "Use heat transport for the ocean ?"
     ok_slab_heat_transp=.true.   ! default value
     call getin_p("ok_slab_heat_transp",ok_slab_heat_transp)
     write(*,*) "ok_slab_heat_transp = ",ok_slab_heat_transp

! Photochemistry and chemistry in the thermosphere

     write(*,*) "Use photochemistry ?"
     photochem=.false.         ! default value
     call getin_p("photochem",photochem)
     write(*,*) "photochem = ",photochem

     write(*,*)"Production of haze ?"
     haze=.false. ! default value
     call getin_p("haze",haze)
     write(*,*)" haze = ",haze


! Test of incompatibility:
! if kastprof used, we must be in 1D
     if (kastprof.and.(ngrid.gt.1)) then
       print*,'kastprof can only be used in 1D!'
       call abort
     endif

     write(*,*)"Stratospheric temperature for kastprof mode?"
     Tstrat=167.0
     call getin_p("Tstrat",Tstrat)
     write(*,*)" Tstrat = ",Tstrat

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

     write(*,*)"Default planetary temperature?"
     tplanet=215.0
     call getin_p("tplanet",tplanet)
     write(*,*)" tplanet = ",tplanet

     write(*,*)"Which star?"
     startype=1 ! default value = Sol
     call getin_p("startype",startype)
     write(*,*)" startype = ",startype

     write(*,*)"Value of stellar flux at 1 AU?"
     Fat1AU=1356.0 ! default value = Sol today
     call getin_p("Fat1AU",Fat1AU)
     write(*,*)" Fat1AU = ",Fat1AU


! TRACERS:

     write(*,*)"Varying H2O cloud fraction?"
     CLFvarying=.false.     ! default value
     call getin_p("CLFvarying",CLFvarying)
     write(*,*)" CLFvarying = ",CLFvarying

     write(*,*)"Value of fixed H2O cloud fraction?"
     CLFfixval=1.0                ! default value
     call getin_p("CLFfixval",CLFfixval)
     write(*,*)" CLFfixval = ",CLFfixval

     write(*,*)"fixed radii for Cloud particles?"
     radfixed=.false. ! default value
     call getin_p("radfixed",radfixed)
     write(*,*)" radfixed = ",radfixed

     if(kastprof)then
        radfixed=.true.
     endif  

     write(*,*)"Number mixing ratio of CO2 ice particles:"
     Nmix_co2=1.e6 ! default value
     call getin_p("Nmix_co2",Nmix_co2)
     write(*,*)" Nmix_co2 = ",Nmix_co2

!         write(*,*)"Number of radiatively active aerosols:"
!         naerkind=0. ! default value
!         call getin_p("naerkind",naerkind)
!         write(*,*)" naerkind = ",naerkind

     write(*,*)"Opacity of dust (if used):"
     dusttau=0. ! default value
     call getin_p("dusttau",dusttau)
     write(*,*)" dusttau = ",dusttau

     write(*,*)"Radiatively active CO2 aerosols?"
     aeroco2=.false.     ! default value
     call getin_p("aeroco2",aeroco2)
     write(*,*)" aeroco2 = ",aeroco2

     write(*,*)"Fixed CO2 aerosol distribution?"
     aerofixco2=.false.     ! default value
     call getin_p("aerofixco2",aerofixco2)
     write(*,*)" aerofixco2 = ",aerofixco2

     write(*,*)"Radiatively active water ice?"
     aeroh2o=.false.     ! default value
     call getin_p("aeroh2o",aeroh2o)
     write(*,*)" aeroh2o = ",aeroh2o

     write(*,*)"Fixed H2O aerosol distribution?"
     aerofixh2o=.false.     ! default value
     call getin_p("aerofixh2o",aerofixh2o)
     write(*,*)" aerofixh2o = ",aerofixh2o

     write(*,*)"Radiatively active sulfuric acid aerosols?"
     aeroh2so4=.false.     ! default value
     call getin_p("aeroh2so4",aeroh2so4)
     write(*,*)" aeroh2so4 = ",aeroh2so4
	 
!=================================

     write(*,*)"Radiatively active two-layer aerosols?"
     aeroback2lay=.false.     ! default value
     call getin_p("aeroback2lay",aeroback2lay)
     write(*,*)" aeroback2lay = ",aeroback2lay

     write(*,*)"Radiatively active ammonia cloud?"
     aeronh3=.false.     ! default value
     call getin_p("aeronh3",aeronh3)
     write(*,*)" aeronh3 = ",aeronh3

     write(*,*)"Radiatively active auroral aerosols?"
     aeroaurora=.false.     ! default value
     call getin_p("aeroaurora",aeroaurora)
     write(*,*)" aeroaurora = ",aeroaurora

     write(*,*)"TWOLAY AEROSOL: total optical depth ", &
                    "in the tropospheric layer (visible)"
     obs_tau_col_tropo=8.D0
     call getin_p("obs_tau_col_tropo",obs_tau_col_tropo)
     write(*,*)" obs_tau_col_tropo = ",obs_tau_col_tropo

     write(*,*)"TWOLAY AEROSOL: total optical depth ", &
                    "in the stratospheric layer (visible)"
     obs_tau_col_strato=0.08D0
     call getin_p("obs_tau_col_strato",obs_tau_col_strato)
     write(*,*)" obs_tau_col_strato = ",obs_tau_col_strato

     write(*,*)"TWOLAY AEROSOL: optprop_back2lay_vis?"
     optprop_back2lay_vis = 'optprop_saturn_vis_n20.dat'
     call getin_p("optprop_back2lay_vis",optprop_back2lay_vis)
     write(*,*)" optprop_back2lay_vis = ",optprop_back2lay_vis

     write(*,*)"TWOLAY AEROSOL: optprop_back2lay_ir?"
     optprop_back2lay_ir = 'optprop_saturn_ir_n20.dat'
     call getin_p("optprop_back2lay_ir",optprop_back2lay_ir)
     write(*,*)" optprop_back2lay_ir = ",optprop_back2lay_ir
     
     write(*,*)"TWOLAY AEROSOL: pres_bottom_tropo? in pa"
     pres_bottom_tropo=66000.0
     call getin_p("pres_bottom_tropo",pres_bottom_tropo)
     write(*,*)" pres_bottom_tropo = ",pres_bottom_tropo

     write(*,*)"TWOLAY AEROSOL: pres_top_tropo? in pa"
     pres_top_tropo=18000.0
     call getin_p("pres_top_tropo",pres_top_tropo)
     write(*,*)" pres_top_tropo = ",pres_top_tropo

     write(*,*)"TWOLAY AEROSOL: pres_bottom_strato? in pa"
     pres_bottom_strato=2000.0
     call getin_p("pres_bottom_strato",pres_bottom_strato)
     write(*,*)" pres_bottom_strato = ",pres_bottom_strato

     ! Sanity check
     if (pres_bottom_strato .gt. pres_top_tropo) then
       print*,'Error : TWOLAY AEROSOL, Please ensure that in callphys.def'
       print*,'you have pres_top_tropo > pres_bottom_strato !'
       stop
     endif

     write(*,*)"TWOLAY AEROSOL: pres_top_strato? in pa"
     pres_top_strato=100.0
     call getin_p("pres_top_strato",pres_top_strato)
     write(*,*)" pres_top_strato = ",pres_top_strato

     write(*,*)"TWOLAY AEROSOL: particle size in the ", &
                    "tropospheric layer, in meters"
     size_tropo=2.e-6
     call getin_p("size_tropo",size_tropo)
     write(*,*)" size_tropo = ",size_tropo

     write(*,*)"TWOLAY AEROSOL: particle size in the ", &
                    "stratospheric layer, in meters"
     size_strato=1.e-7
     call getin_p("size_strato",size_strato)
     write(*,*)" size_strato = ",size_strato

     write(*,*)"NH3 (thin) cloud: total optical depth"
     tau_nh3_cloud=7.D0
     call getin_p("tau_nh3_cloud",tau_nh3_cloud)
     write(*,*)" tau_nh3_cloud = ",tau_nh3_cloud

     write(*,*)"NH3 (thin) cloud pressure level"
     pres_nh3_cloud=7.D0
     call getin_p("pres_nh3_cloud",pres_nh3_cloud)
     write(*,*)" pres_nh3_cloud = ",pres_nh3_cloud

     write(*,*)"NH3 (thin) cloud: particle sizes"
     size_nh3_cloud=3.e-6
     call getin_p("size_nh3_cloud",size_nh3_cloud)
     write(*,*)" size_nh3_cloud = ",size_nh3_cloud

!=================================

     write(*,*)"Cloud pressure level (with kastprof only):"
     cloudlvl=0. ! default value
     call getin_p("cloudlvl",cloudlvl)
     write(*,*)" cloudlvl = ",cloudlvl

     write(*,*)"Is the variable gas species radiatively active?"
     Tstrat=167.0
     varactive=.false.
     call getin_p("varactive",varactive)
     write(*,*)" varactive = ",varactive

     write(*,*)"Is the variable gas species distribution set?"
     varfixed=.false.
     call getin_p("varfixed",varfixed)
     write(*,*)" varfixed = ",varfixed

     write(*,*)"What is the saturation % of the variable species?"
     satval=0.8
     call getin_p("satval",satval)
     write(*,*)" satval = ",satval


! Test of incompatibility:
! if varactive, then varfixed should be false
     if (varactive.and.varfixed) then
       print*,'if varactive, varfixed must be OFF!'
       stop
     endif

     write(*,*) "Gravitationnal sedimentation ?"
     sedimentation=.false. ! default value
     call getin_p("sedimentation",sedimentation)
     write(*,*) " sedimentation = ",sedimentation

     write(*,*) "Compute water cycle ?"
     water=.false. ! default value
     call getin_p("water",water)
     write(*,*) " water = ",water
         
! Test of incompatibility:
! if water is true, there should be at least a tracer
     if (water.and.(.not.tracer)) then
        print*,'if water is ON, tracer must be ON too!'
        stop
     endif

     write(*,*) "Include water condensation ?"
     watercond=.false. ! default value
     call getin_p("watercond",watercond)
     write(*,*) " watercond = ",watercond

! Test of incompatibility:
! if watercond is used, then water should be used too
     if (watercond.and.(.not.water)) then
        print*,'if watercond is used, water should be used too'
        stop
     endif

     write(*,*) "Include water precipitation ?"
     waterrain=.false. ! default value
     call getin_p("waterrain",waterrain)
     write(*,*) " waterrain = ",waterrain

     write(*,*) "Include surface hydrology ?"
     hydrology=.false. ! default value
     call getin_p("hydrology",hydrology)
     write(*,*) " hydrology = ",hydrology

     write(*,*) "Evolve surface water sources ?"
     sourceevol=.false. ! default value
     call getin_p("sourceevol",sourceevol)
     write(*,*) " sourceevol = ",sourceevol

     write(*,*) "Ice evolution timestep ?"
     icetstep=100.0 ! default value
     call getin_p("icetstep",icetstep)
     write(*,*) " icetstep = ",icetstep
         
     write(*,*) "Spectral Dependant albedo ?"
     albedo_spectral_mode=.false. ! default value
     call getin_p("albedo_spectral_mode",albedo_spectral_mode)
     write(*,*) " albedo_spectral_mode = ",albedo_spectral_mode

     write(*,*) "Snow albedo ?"
     write(*,*) "If albedo_spectral_mode=.true., then this "
     write(*,*) "corresponds to the 0.5 microns snow albedo."
     albedosnow=0.5         ! default value
     call getin_p("albedosnow",albedosnow)
     write(*,*) " albedosnow = ",albedosnow
         
     write(*,*) "CO2 ice albedo ?"
     albedoco2ice=0.5       ! default value
     call getin_p("albedoco2ice",albedoco2ice)
     write(*,*) " albedoco2ice = ",albedoco2ice

     write(*,*) "Maximum ice thickness ?"
     maxicethick=2.0         ! default value
     call getin_p("maxicethick",maxicethick)
     write(*,*) " maxicethick = ",maxicethick

     write(*,*) "Freezing point of seawater ?"
     Tsaldiff=-1.8          ! default value
     call getin_p("Tsaldiff",Tsaldiff)
     write(*,*) " Tsaldiff = ",Tsaldiff

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
     call su_gases
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


!-----------------------------------------------------------------------
!     Some more initialization:
!     ------------------------

  ! Initializations for comgeomfi_h
#ifndef MESOSCALE
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
#endif
  ! initialize variables in radinc_h
  call ini_radinc_h(nlayer,tplanckmin,tplanckmax,dtplanck)
 
  ! allocate "comsoil_h" arrays
  call ini_comsoil_h(ngrid)
    
  END SUBROUTINE inifis

END MODULE inifis_mod
