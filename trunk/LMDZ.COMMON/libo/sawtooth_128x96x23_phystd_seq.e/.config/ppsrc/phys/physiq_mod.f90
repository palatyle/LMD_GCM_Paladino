










      module physiq_mod
      
      implicit none
      
      contains
      
      subroutine physiq(ngrid,nlayer,nq,   &
                  nametrac,                & 
                  firstcall,lastcall,      &
                  pday,ptime,ptimestep,    &
                  pplev,pplay,pphi,        &
                  pu,pv,pt,pq,             &
                  flxw,                    &
                  pdu,pdv,pdt,pdq,pdpsrf)
 
      use radinc_h, only : L_NSPECTI,L_NSPECTV,naerkind
      use watercommon_h, only : RLVTT, Psat_water,epsi,su_watercycle, RV, T_h2o_ice_liq
      use thermcell_mod, only: init_thermcell_mod
      use gases_h, only: gnom, gfrac
      use radcommon_h, only: sigma, glat, grav, BWNV
      use radii_mod, only: h2o_reffrad, co2_reffrad
      use aerosol_mod, only: iaero_co2, iaero_h2o
      use surfdat_h, only: phisfi, zmea, zstd, zsig, zgam, zthe, &
                           dryness, watercaptag
      use datafile_mod, only: datadir
      use comdiurn_h, only: coslat, sinlat, coslon, sinlon
      use comsaison_h, only: mu0, fract, dist_star, declin, right_ascen
      use comsoil_h, only: nsoilmx, layer, mlayer, inertiedat
      use geometry_mod, only: latitude, longitude, latitude_deg, longitude_deg, cell_area
      USE comgeomfi_h, only: totarea, totarea_planet
      USE tracer_h, only: noms, mmol, radius, rho_q, qext, &
                          alpha_lift, alpha_devil, qextrhor, &
                          igcm_h2o_ice, igcm_h2o_vap, igcm_dustbin, &
                          igcm_co2_ice, igcm_volc_1, igcm_volc_2, &
                          igcm_volc_3, igcm_volc_4, igcm_volc_5, &
                          igcm_volc_6, igcm_h2so4
      use time_phylmdz_mod, only: ecritphy, iphysiq, nday
      use phyetat0_mod, only: phyetat0
      use phyredem, only: physdem0, physdem1
      use slab_ice_h, only: capcalocean, capcalseaice,capcalsno, &
                            noceanmx
      use ocean_slab_mod, only :ocean_slab_init, ocean_slab_ice, &
                                ini_surf_heat_transp_mod, &
                                ocean_slab_get_vars,ocean_slab_final
      use surf_heat_transp_mod,only :init_masquv
      use planetwide_mod, only: planetwide_minval,planetwide_maxval,planetwide_sumval
      use mod_phys_LMDZ_mpi_data, only: mpi_rank, mpi_size
      use print_control_mod, only: lunout
      use mod_phys_lmdz_para, only : is_master
      use planete_mod, only: apoastr, periastr, year_day, peri_day, &
                            obliquit, nres, z0
      use comcstfi_mod, only: pi, g, rcp, r, rad, mugaz, cpp
      use time_phylmdz_mod, only: daysec
      use callkeys_mod
      use conc_mod
      use phys_state_var_mod
      use callcorrk_mod, only: callcorrk
      use turb_mod, only : q2,sensibFlux,turb_resolved
      use vertical_layers_mod, only: presnivs, pseudoalt
      use mod_phys_lmdz_omp_data, ONLY: is_omp_master


      implicit none


!==================================================================
!     
!     Purpose
!     -------
!     Central subroutine for all the physics parameterisations in the
!     universal model. Originally adapted from the Mars LMDZ model.
!
!     The model can be run without or with tracer transport
!     depending on the value of "tracer" in file "callphys.def".
!
!
!   It includes:
!
!      I. Initialization :
!         I.1 Firstcall initializations.
!         I.2 Initialization for every call to physiq.
!
!      II. Compute radiative transfer tendencies (longwave and shortwave) :
!         II.a Option 1 : Call correlated-k radiative transfer scheme.
!         II.b Option 2 : Call Newtonian cooling scheme.
!         II.c Option 3 : Atmosphere has no radiative effect.
!
!      III. Vertical diffusion (turbulent mixing) :
!
!      IV. Convection :
!         IV.a Thermal plume model
!         IV.b Dry convective adjusment
!
!      V. Condensation and sublimation of gases (currently just CO2).
!
!      VI. Tracers
!         VI.1. Water and water ice.
!         VI.2  Photochemistry
!         VI.3. Aerosols and particles.
!         VI.4. Updates (pressure variations, surface budget).
!         VI.5. Slab Ocean.
!         VI.6. Surface Tracer Update.
!
!      VII. Surface and sub-surface soil temperature.
!
!      VIII. Perform diagnostics and write output files.
!
!
!   arguments
!   ---------
!
!   INPUT
!   -----
!
!    ngrid                 Size of the horizontal grid.
!    nlayer                Number of vertical layers.
!    nq                    Number of advected fields.
!    nametrac              Name of corresponding advected fields.
!
!    firstcall             True at the first call.
!    lastcall              True at the last call.
!
!    pday                  Number of days counted from the North. Spring equinoxe.
!    ptime                 Universal time (0<ptime<1): ptime=0.5 at 12:00 UT.
!    ptimestep             timestep (s).
!
!    pplay(ngrid,nlayer)   Pressure at the middle of the layers (Pa).
!    pplev(ngrid,nlayer+1) Intermediate pressure levels (Pa).
!    pphi(ngrid,nlayer)    Geopotential at the middle of the layers (m2.s-2).
!
!    pu(ngrid,nlayer)      u, zonal component of the wind (ms-1).
!    pv(ngrid,nlayer)      v, meridional component of the wind (ms-1).
!
!    pt(ngrid,nlayer)      Temperature (K).
!
!    pq(ngrid,nlayer,nq)   Advected fields.
!
!    pudyn(ngrid,nlayer)    !    pvdyn(ngrid,nlayer)     \ Dynamical temporal derivative for the
!    ptdyn(ngrid,nlayer)     / corresponding variables.
!    pqdyn(ngrid,nlayer,nq) /
!    flxw(ngrid,nlayer)      vertical mass flux (kg/s) at layer lower boundary
!
!   OUTPUT
!   ------
!
!    pdu(ngrid,nlayer)        !    pdv(ngrid,nlayer)         \  Temporal derivative of the corresponding
!    pdt(ngrid,nlayer)         /  variables due to physical processes.
!    pdq(ngrid,nlayer)        /
!    pdpsrf(ngrid)             /
!
!
!     Authors
!     -------
!           Frederic Hourdin        15/10/93
!           Francois Forget        1994
!           Christophe Hourdin        02/1997 
!           Subroutine completely rewritten by F. Forget (01/2000)
!           Water ice clouds: Franck Montmessin (update 06/2003)
!           Radiatively active tracers: J.-B. Madeleine (10/2008-06/2009)
!           New correlated-k radiative scheme: R. Wordsworth (2009)
!           Many specifically Martian subroutines removed: R. Wordsworth (2009)       
!           Improved water cycle: R. Wordsworth / B. Charnay (2010)
!           To F90: R. Wordsworth (2010)
!           New turbulent diffusion scheme: J. Leconte (2012)
!           Loops converted to F90 matrix format: J. Leconte (2012)
!           No more ngridmx/nqmx, F90 commons and adaptation to parallel: A. Spiga (2012)
!           Purge of the code : M. Turbet (2015)
!           Photochemical core developped by F. Lefevre: B. Charnay (2017)
!==================================================================


!    0.  Declarations :
!    ------------------

include "netcdf.inc"
include "dimensions.h"
! Arguments :
! -----------

!   INPUTS:
!   -------

      integer,intent(in) :: ngrid             ! Number of atmospheric columns.
      integer,intent(in) :: nlayer            ! Number of atmospheric layers.
      integer,intent(in) :: nq                ! Number of tracers.
      character*30,intent(in) :: nametrac(nq) ! Names of the tracers taken from dynamics.
      
      logical,intent(in) :: firstcall ! Signals first call to physics.
      logical,intent(in) :: lastcall  ! Signals last call to physics.
      
      real,intent(in) :: pday                  ! Number of elapsed sols since reference Ls=0.
      real,intent(in) :: ptime                 ! "Universal time", given as fraction of sol (e.g.: 0.5 for noon).
      real,intent(in) :: ptimestep             ! Physics timestep (s).
      real,intent(in) :: pplev(ngrid,nlayer+1) ! Inter-layer pressure (Pa).
      real,intent(in) :: pplay(ngrid,nlayer)   ! Mid-layer pressure (Pa).
      real,intent(in) :: pphi(ngrid,nlayer)    ! Geopotential at mid-layer (m2s-2).
      real,intent(in) :: pu(ngrid,nlayer)      ! Zonal wind component (m/s).
      real,intent(in) :: pv(ngrid,nlayer)      ! Meridional wind component (m/s).
      real,intent(in) :: pt(ngrid,nlayer)      ! Temperature (K).
      real,intent(in) :: pq(ngrid,nlayer,nq)   ! Tracers (kg/kg_of_air).
      real,intent(in) :: flxw(ngrid,nlayer)    ! Vertical mass flux (ks/s) at lower boundary of layer

!   OUTPUTS:
!   --------

!     Physical tendencies :

      real,intent(out) :: pdu(ngrid,nlayer)    ! Zonal wind tendencies (m/s/s).
      real,intent(out) :: pdv(ngrid,nlayer)    ! Meridional wind tendencies (m/s/s).
      real,intent(out) :: pdt(ngrid,nlayer)    ! Temperature tendencies (K/s).
      real,intent(out) :: pdq(ngrid,nlayer,nq) ! Tracer tendencies (kg/kg_of_air/s).
      real,intent(out) :: pdpsrf(ngrid)        ! Surface pressure tendency (Pa/s).

! Local saved variables:
! ----------------------
      integer,save :: day_ini                                      ! Initial date of the run (sol since Ls=0).
      integer,save :: icount                                       ! Counter of calls to physiq during the run.
!$OMP THREADPRIVATE(day_ini,icount)

! Local variables :
! -----------------

!     Aerosol (dust or ice) extinction optical depth  at reference wavelength 
!     for the "naerkind" optically active aerosols:

      real aerosol(ngrid,nlayer,naerkind) ! Aerosols.
      real zh(ngrid,nlayer)               ! Potential temperature (K).
      real pw(ngrid,nlayer)               ! Vertical velocity (m/s). (NOTE : >0 WHEN DOWNWARDS !!)
      real omega(ngrid,nlayer)            ! omega velocity (Pa/s, >0 when downward)

      integer l,ig,ierr,iq,nw,isoil
      
      real zls                       ! Solar longitude (radians).
      real zlss                      ! Sub solar point longitude (radians).
      real zday                      ! Date (time since Ls=0, calculated in sols).
      real zzlay(ngrid,nlayer)       ! Altitude at the middle of the atmospheric layers.
      real zzlev(ngrid,nlayer+1)     ! Altitude at the atmospheric layer boundaries.

! VARIABLES for the thermal plume model
      
      real f(ngrid)                       ! Mass flux norm
      real fm(ngrid,nlayer+1)             ! Mass flux
      real fm_bis(ngrid,nlayer)           ! Recasted fm
      real entr(ngrid,nlayer)             ! Entrainment
      real detr(ngrid,nlayer)             ! Detrainment
      real dqevap(ngrid,nlayer,nq)        ! water tracer mass mixing ratio variations due to evaporation
      real dtevap(ngrid,nlayer)           ! temperature variation due to evaporation
      real zqtherm(ngrid,nlayer,nq)       ! vapor mass mixing ratio after evaporation
      real zttherm(ngrid,nlayer)          ! temperature after evaporation
      real fraca(ngrid,nlayer+1)          ! Fraction of the surface that plumes occupies
      real zw2(ngrid,nlayer+1)            ! Vertical speed
      real zw2_bis(ngrid,nlayer)          ! Recasted zw2
      
! TENDENCIES due to various processes :

      ! For Surface Temperature : (K/s)
      real zdtsurf(ngrid)     ! Cumulated tendencies.
      real zdtsurfmr(ngrid)   ! Mass_redistribution routine.
      real zdtsurfc(ngrid)    ! Condense_co2 routine.
      real zdtsdif(ngrid)     ! Turbdiff/vdifc routines.
      real zdtsurf_hyd(ngrid) ! Hydrol routine.
            
      ! For Atmospheric Temperatures : (K/s)    
      real dtlscale(ngrid,nlayer)                             ! Largescale routine.
      real zdtc(ngrid,nlayer)                                 ! Condense_co2 routine.
      real zdtdif(ngrid,nlayer)                               ! Turbdiff/vdifc routines.
      real zdttherm(ngrid,nlayer)                             ! Calltherm routine.
      real zdtmr(ngrid,nlayer)                                ! Mass_redistribution routine.
      real zdtrain(ngrid,nlayer)                              ! Rain routine.
      real dtmoist(ngrid,nlayer)                              ! Moistadj routine.
      real dt_ekman(ngrid,noceanmx), dt_hdiff(ngrid,noceanmx) ! Slab_ocean routine.
      real zdtsw1(ngrid,nlayer), zdtlw1(ngrid,nlayer)         ! Callcorrk routine.
                              
      ! For Surface Tracers : (kg/m2/s)
      real dqsurf(ngrid,nq)                 ! Cumulated tendencies.
      real zdqsurfc(ngrid)                  ! Condense_co2 routine.
      real zdqsdif(ngrid,nq)                ! Turbdiff/vdifc routines.
      real zdqssed(ngrid,nq)                ! Callsedim routine.
      real zdqsurfmr(ngrid,nq)              ! Mass_redistribution routine.
      real zdqsrain(ngrid), zdqssnow(ngrid) ! Rain routine.
      real dqs_hyd(ngrid,nq)                ! Hydrol routine.
      real reevap_precip(ngrid)             ! re-evaporation flux of precipitation (integrated over the atmospheric column)
                  
      ! For Tracers : (kg/kg_of_air/s)
      real zdqc(ngrid,nlayer,nq)      ! Condense_co2 routine.
      real zdqadj(ngrid,nlayer,nq)    ! Convadj routine.
      real zdqdif(ngrid,nlayer,nq)    ! Turbdiff/vdifc routines.
      real zdqevap(ngrid,nlayer)      ! Turbdiff routine.
      real zdqtherm(ngrid,nlayer,nq)  ! Calltherm routine.
      real zdqsed(ngrid,nlayer,nq)    ! Callsedim routine.
      real zdqmr(ngrid,nlayer,nq)     ! Mass_redistribution routine.
      real zdqrain(ngrid,nlayer,nq)   ! Rain routine.
      real dqmoist(ngrid,nlayer,nq)   ! Moistadj routine.
      real dqvaplscale(ngrid,nlayer)  ! Largescale routine.
      real dqcldlscale(ngrid,nlayer)  ! Largescale routine.
      REAL zdqchim(ngrid,nlayer,nq)   ! Calchim_asis routine
      REAL zdqschim(ngrid,nq)         ! Calchim_asis routine
      REAL zdqvolc(ngrid,nlayer,nq)   ! Volcanism: Source     

      REAL array_zero1(ngrid)
      REAL array_zero2(ngrid,nlayer)
                        
      ! For Winds : (m/s/s)
      real zdvadj(ngrid,nlayer), zduadj(ngrid,nlayer)       ! Convadj routine.
      real zdutherm(ngrid,nlayer), zdvtherm(ngrid,nlayer)   ! Calltherm routine.
      real zdumr(ngrid,nlayer), zdvmr(ngrid,nlayer)         ! Mass_redistribution routine.
      real zdvdif(ngrid,nlayer), zdudif(ngrid,nlayer)       ! Turbdiff/vdifc routines.
      real zdhdif(ngrid,nlayer)                             ! Turbdiff/vdifc routines.
      real zdhadj(ngrid,nlayer)                             ! Convadj routine.
     
      ! For Pressure and Mass :
      real zdmassmr(ngrid,nlayer) ! Atmospheric Mass tendency for mass_redistribution (kg_of_air/m2/s).
      real zdmassmr_col(ngrid)    ! Atmospheric Column Mass tendency for mass_redistribution (kg_of_air/m2/s).
      real zdpsrfmr(ngrid)        ! Pressure tendency for mass_redistribution routine (Pa/s).

      
    
! Local variables for LOCAL CALCULATIONS:
! ---------------------------------------
      real zflubid(ngrid)
      real zplanck(ngrid),zpopsk(ngrid,nlayer)
      real ztim1,ztim2,ztim3, z1,z2
      real ztime_fin
      real zdh(ngrid,nlayer)
      real gmplanet
      real taux(ngrid),tauy(ngrid)


! local variables for DIAGNOSTICS : (diagfi & stat)
! -------------------------------------------------
      real ps(ngrid)                                     ! Surface Pressure.
      real zt(ngrid,nlayer)                              ! Atmospheric Temperature.
      real zu(ngrid,nlayer),zv(ngrid,nlayer)             ! Zonal and Meridional Winds.
      real zq(ngrid,nlayer,nq)                           ! Atmospheric Tracers.
      real zdtadj(ngrid,nlayer)                          ! Convadj Diagnostic.
      real zdtdyn(ngrid,nlayer)                          ! Dynamical Heating (K/s).
      real zdudyn(ngrid,nlayer)                          ! Dynamical Zonal Wind tendency (m.s-2).

      real reff(ngrid,nlayer)                       ! Effective dust radius (used if doubleq=T).
      real vmr(ngrid,nlayer)                        ! volume mixing ratio
      real time_phys
      
      real ISR,ASR,OLR,GND,DYN,GSR,Ts1,Ts2,Ts3,TsS ! for Diagnostic.
      
      real qcol(ngrid,nq) ! Tracer Column Mass (kg/m2).

!     included by RW for H2O Manabe scheme
      real rneb_man(ngrid,nlayer) ! H2O cloud fraction (moistadj).
      real rneb_lsc(ngrid,nlayer) ! H2O cloud fraction (large scale).


!     to test energy conservation (RW+JL)
      real mass(ngrid,nlayer),massarea(ngrid,nlayer)
      real dEtot, dEtots, AtmToSurf_TurbFlux
      real,save :: dEtotSW, dEtotsSW, dEtotLW, dEtotsLW
!$OMP THREADPRIVATE(dEtotSW, dEtotsSW, dEtotLW, dEtotsLW)
      real dEzRadsw(ngrid,nlayer),dEzRadlw(ngrid,nlayer),dEzdiff(ngrid,nlayer)
      real dEdiffs(ngrid),dEdiff(ngrid)
      real madjdE(ngrid), lscaledE(ngrid),madjdEz(ngrid,nlayer), lscaledEz(ngrid,nlayer)
      
!JL12 conservation test for mean flow kinetic energy has been disabled temporarily

      real dtmoist_max,dtmoist_min     
      real dItot, dItot_tmp, dVtot, dVtot_tmp


      real h2otot                      ! Total Amount of water. For diagnostic.
      real icesrf,liqsrf,icecol,vapcol ! Total Amounts of water (ice,liq,vap). For diagnostic.
      real dWtot, dWtot_tmp, dWtots, dWtots_tmp
      logical,save :: watertest
!$OMP THREADPRIVATE(watertest)

      real qsat(ngrid,nlayer) ! Water Vapor Volume Mixing Ratio at saturation (kg/kg_of_air).
      real RH(ngrid,nlayer)   ! Relative humidity.
      real H2Omaxcol(ngrid)   ! Maximum possible H2O column amount (at 100% saturation) (kg/m2).
      real psat_tmp
      
      logical clearsky ! For double radiative transfer call. By BC
      
      ! For Clear Sky Case.
      real fluxsurf_lw1(ngrid), fluxsurf_sw1(ngrid), fluxsurfabs_sw1(ngrid)  ! For SW/LW flux.
      real fluxtop_lw1(ngrid), fluxabs_sw1(ngrid)                            ! For SW/LW flux.
      real albedo_equivalent1(ngrid)                                         ! For Equivalent albedo calculation.
      real tau_col1(ngrid)                                                   ! For aerosol optical depth diagnostic.
      real OLR_nu1(ngrid,L_NSPECTI), OSR_nu1(ngrid,L_NSPECTV)                ! For Outgoing Radiation diagnostics.
      real int_dtaui1(ngrid,nlayer,L_NSPECTI),int_dtauv1(ngrid,nlayer,L_NSPECTV) ! For optical thickness diagnostics.
      real tf, ntf

      real nconsMAX, vdifcncons(ngrid), cadjncons(ngrid) ! Vdfic water conservation test. By RW

      real muvar(ngrid,nlayer+1) ! For Runaway Greenhouse 1D study. By RW

      real reffcol(ngrid,naerkind)

!     Sourceevol for 'accelerated ice evolution'. By RW
      real delta_ice,ice_tot
      integer num_run
      logical,save :: ice_update


      real :: tsurf2(ngrid)
      real :: flux_o(ngrid),flux_g(ngrid),fluxgrdocean(ngrid)
      real :: flux_sens_lat(ngrid)
      real :: qsurfint(ngrid,nq)

! VOLCANO 
      integer,save :: ivolc, i
      integer,save :: mpi_volc
      REAL :: lon_hs, lat_hs
      REAL :: dlon, dlat
      INTEGER, save :: inputstatus, nlines
      REAL, DIMENSION(:,:), SAVE, ALLOCATABLE  :: data_read_spr, data_read_sum, data_read_fall, data_read_win
      REAL, DIMENSION(:,:), SAVE, ALLOCATABLE  :: data_readt_spr, data_readt_sum, data_readt_fall, data_readt_win
      CHARACTER(LEN=300) :: filename_volc_spr, filename_volc_sum, filename_volc_fall, filename_volc_win
      REAL :: lon_volc 
      REAL :: lat_volc
 !     PLEASE DEFINE THE LOCATION OF THE ERUPTION BELOW :
 ! =============================================
 !     Coordinate of the volcano (degrees) :
 !       ex : Apollinaris Patera lon=174.4 and lat=-9.3
       ! ex : Elysium Mons lon=147 and lat=24.8
       ! ex : Cerberus lon=176.6 and lat=9.0
       ! ex : Olympus Mons lon=-133.9 and lat=18.7
       ! ex : Arsia Mons lon=-120.46 and lat=-9.14
       ! ex : Pavonis Mons lon=-112.85 and lat= 0.662
       ! ex : Ascraeus Mons lon=-104.37 and lat= 11.1
       ! ex : Syrtis Major lon=66.4 and lat= 9.85
       ! ex : Tyrrhenia Patera lon=106.55 and lat= -21.32
       ! ex : Hadriaca Patera lon=92.18 and lat= -30.44
       ! ex : Peneus Patera lon=60.76 and lat= -58.05
       ! ex : Alba Patera lon=-111.11 and lat= 39.19
       ! ex : Amphritites lon=52.66 and lat= -58.00
       ! ex : Hecates lon=150.08 and lat=31.68
       ! ex : Pityusa Patera lon=36.87 and lat=-66.77
       ! ex : Malea Patera lon=50.96 and lat=-63.09
       ! ex : Electris volcano lon=-173.21 and lat =-37.35


      ! Misc
      character*2 :: str2
!==================================================================================================

! -----------------
! I. INITIALISATION
! -----------------

! --------------------------------
! I.1   First Call Initialisation.
! --------------------------------
      if (firstcall) then
        ! Allocate saved arrays (except for 1D model, where this has already
        ! been done)
        if (ngrid>1) call phys_state_var_init(nq)

!        Variables set to 0
!        ~~~~~~~~~~~~~~~~~~
         dtrad(:,:) = 0.0
         fluxrad(:) = 0.0
         tau_col(:) = 0.0
         zdtsw(:,:) = 0.0
         zdtlw(:,:) = 0.0


!        Initialize aerosol indexes.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call iniaerosol()

      
!        Initialize tracer names, indexes and properties.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF (.NOT.ALLOCATED(noms)) ALLOCATE(noms(nq)) ! (because noms is an argument of physdem1 whether or not tracer is on)
         if (tracer) then
            call initracer(ngrid,nq,nametrac)
            if(photochem) then
              call ini_conc_mod(ngrid,nlayer)
            endif
         endif


!        Read 'startfi.nc' file.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call phyetat0(startphy_file,                                 &
                       ngrid,nlayer,"startfi.nc",0,0,nsoilmx,nq,      &
                       day_ini,time_phys,tsurf,tsoil,emis,q2,qsurf,   &
                       cloudfrac,totcloudfrac,hice,                   &
                       rnat,pctsrf_sic,tslab, tsea_ice,sea_ice)

         if (.not.startphy_file) then
           ! additionnal "academic" initialization of physics
           if (is_master) write(*,*) "Physiq: initializing tsurf(:) to pt(:,1) !!"
           tsurf(:)=pt(:,1)
           if (is_master) write(*,*) "Physiq: initializing tsoil(:) to pt(:,1) !!"
           do isoil=1,nsoilmx
             tsoil(1:ngrid,isoil)=tsurf(1:ngrid)
           enddo
           if (is_master) write(*,*) "Physiq: initializing day_ini to pdat !"
           day_ini=pday
         endif

         if (pday.ne.day_ini) then
           write(*,*) "ERROR in physiq.F90:"
           write(*,*) "bad synchronization between physics and dynamics"
           write(*,*) "dynamics day: ",pday
           write(*,*) "physics day:  ",day_ini
           stop
         endif

         write (*,*) 'In physiq day_ini =', day_ini

!        Initialize albedo calculation. 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         albedo(:,:)=0.0
          albedo_bareground(:)=0.0
          albedo_snow_SPECTV(:)=0.0
          albedo_co2_ice_SPECTV(:)=0.0
         call surfini(ngrid,nq,qsurf,albedo,albedo_bareground,albedo_snow_SPECTV,albedo_co2_ice_SPECTV) 
         
!        Initialize orbital calculation. 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call iniorbit(apoastr,periastr,year_day,peri_day,obliquit)


         if(tlocked)then
            print*,'Planet is tidally locked at resonance n=',nres
            print*,'Make sure you have the right rotation rate!!!'
         endif

!        Initialize oceanic variables.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (ok_slab_ocean)then

           call ocean_slab_init(ngrid,ptimestep, tslab, &
                                sea_ice, pctsrf_sic)

           call ini_surf_heat_transp_mod()
           
           knindex(:) = 0

           do ig=1,ngrid
              zmasq(ig)=1
              knindex(ig) = ig
              if (nint(rnat(ig)).eq.0) then
                 zmasq(ig)=0
              endif
           enddo

           CALL init_masquv(ngrid,zmasq)

         endif ! end of 'ok_slab_ocean'.


!        Initialize soil.
!        ~~~~~~~~~~~~~~~~
         if (callsoil) then
         
            call soil(ngrid,nsoilmx,firstcall,lastcall,inertiedat, &
                      ptimestep,tsurf,tsoil,capcal,fluxgrd)

            if (ok_slab_ocean) then
               do ig=1,ngrid
                  if (nint(rnat(ig)).eq.2) then
                     capcal(ig)=capcalocean
                     if (pctsrf_sic(ig).gt.0.5) then
                        capcal(ig)=capcalseaice
                        if (qsurf(ig,igcm_h2o_ice).gt.0.) then
                           capcal(ig)=capcalsno
                        endif
                     endif
                  endif
               enddo
            endif ! end of 'ok_slab_ocean'.

         else ! else of 'callsoil'.
         
            print*,'WARNING! Thermal conduction in the soil turned off'
            capcal(:)=1.e6
            fluxgrd(:)=intheat
            print*,'Flux from ground = ',intheat,' W m^-2'
            
         endif ! end of 'callsoil'.
         
         icount=1

!        Decide whether to update ice at end of run.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ice_update=.false.
         if(sourceevol)then
!$OMP MASTER
            open(128,file='num_run',form='formatted', &
                     status="old",iostat=ierr)
            if (ierr.ne.0) then
               write(*,*) "physiq: Error! No num_run file!"
               write(*,*) " (which is needed for sourceevol option)"
               stop
            endif
            read(128,*) num_run 
            close(128)
!$OMP END MASTER
!$OMP BARRIER
        
            if(num_run.ne.0.and.mod(num_run,2).eq.0)then
               print*,'Updating ice at end of this year!'
               ice_update=.true.
               ice_min(:)=1.e4
            endif
            
         endif ! end of 'sourceevol'.


         ! Here is defined the type of the surface : Continent or Ocean.
         ! BC2014 : This is now already done in newstart.
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (.not.ok_slab_ocean) then
         
           rnat(:)=1.
           do ig=1,ngrid
              if(inertiedat(ig,1).gt.1.E4)then
                 rnat(ig)=0
              endif
           enddo

           print*,'WARNING! Surface type currently decided by surface inertia'
           print*,'This should be improved e.g. in newstart.F'
           
         endif ! end of 'ok_slab_ocean'.


!        Initialize surface history variable.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         qsurf_hist(:,:)=qsurf(:,:)

!        Initialize variable for dynamical heating and zonal wind tendency diagnostic
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ztprevious(:,:)=pt(:,:)
         zuprevious(:,:)=pu(:,:)

!        Set temperature just above condensation temperature (for Early Mars)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if(nearco2cond) then
            write(*,*)' WARNING! Starting at Tcond+1K'
            do l=1, nlayer
               do ig=1,ngrid
                  pdt(ig,l)= ((-3167.8)/(log(.01*pplay(ig,l))-23.23)+4     &
                      -pt(ig,l)) / ptimestep
               enddo
            enddo
         endif

         if(meanOLR)then          
            call system('rm -f rad_bal.out') ! to record global radiative balance.           
            call system('rm -f tem_bal.out') ! to record global mean/max/min temperatures.            
            call system('rm -f h2o_bal.out') ! to record global hydrological balance.
         endif


         watertest=.false.        
         if(water)then ! Initialize variables for water cycle
            
            if(enertest)then
               watertest = .true.
            endif

            if(ice_update)then
               ice_initial(:)=qsurf(:,igcm_h2o_ice)
            endif

         endif
         
!        Set some parameters for the thermal plume model
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (calltherm) then
            CALL init_thermcell_mod(g, rcp, r, pi, T_h2o_ice_liq, RV)
         endif
         
         call su_watercycle ! even if we don't have a water cycle, we might
                            ! need epsi for the wvp definitions in callcorrk.F
                            ! or RETV, RLvCp for the thermal plume model
         if (ngrid.ne.1) then ! Note : no need to create a restart file in 1d.
            call physdem0("restartfi.nc",longitude,latitude,nsoilmx,ngrid,nlayer,nq, &
                          ptimestep,pday+nday,time_phys,cell_area,          &
                          albedo_bareground,inertiedat,zmea,zstd,zsig,zgam,zthe)
         endif
         
         ! XIOS outputs
         write(*,*) "physiq: end of firstcall"
      endif ! end of 'firstcall'

! ------------------------------------------------------
! I.2   Initializations done at every physical timestep:
! ------------------------------------------------------


      ! Initialize various variables
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      if ( .not.nearco2cond ) then
         pdt(1:ngrid,1:nlayer) = 0.0
      endif      
      zdtsurf(1:ngrid)      = 0.0
      pdq(1:ngrid,1:nlayer,1:nq) = 0.0
      dqsurf(1:ngrid,1:nq)= 0.0
      pdu(1:ngrid,1:nlayer) = 0.0
      pdv(1:ngrid,1:nlayer) = 0.0
      pdpsrf(1:ngrid)       = 0.0
      zflubid(1:ngrid)      = 0.0      
      flux_sens_lat(1:ngrid) = 0.0
      taux(1:ngrid) = 0.0
      tauy(1:ngrid) = 0.0

      zday=pday+ptime ! Compute time, in sols (and fraction thereof).

      ! Compute Stellar Longitude (Ls), and orbital parameters.
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (season) then
         call stellarlong(zday,zls)
      else
         call stellarlong(float(day_ini),zls)
      end if

      call orbite(zls,dist_star,declin,right_ascen)
      
      if (tlocked) then
              zlss=Mod(-(2.*pi*(zday/year_day)*nres - right_ascen),2.*pi)
      elseif (diurnal) then
              zlss=-2.*pi*(zday-.5)
      else if(diurnal .eqv. .false.) then
              zlss=9999.
      endif 


      ! Compute variations of g with latitude (oblate case).
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (oblate .eqv. .false.) then      
         glat(:) = g         
      else if (flatten .eq. 0.0 .or. J2 .eq. 0.0 .or. Rmean .eq. 0.0 .or. MassPlanet .eq. 0.0) then      
         print*,'I need values for flatten, J2, Rmean and MassPlanet to compute glat (else set oblate=.false.)'
         call abort
      else
         gmplanet = MassPlanet*grav*1e24
         do ig=1,ngrid
            glat(ig)= gmplanet/(Rmean**2) * (1.D0 + 0.75 *J2 - 2.0*flatten/3. + (2.*flatten - 15./4.* J2) * cos(2. * (pi/2. - latitude(ig)))) 
         end do
      endif

      ! Compute geopotential between layers.
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zzlay(1:ngrid,1:nlayer)=pphi(1:ngrid,1:nlayer)
      do l=1,nlayer         
         zzlay(1:ngrid,l)= zzlay(1:ngrid,l)/glat(1:ngrid)
      enddo

      zzlev(1:ngrid,1)=0.
      zzlev(1:ngrid,nlayer+1)=1.e7 ! Dummy top of last layer above 10000 km...

      do l=2,nlayer
         do ig=1,ngrid
            z1=(pplay(ig,l-1)+pplev(ig,l))/(pplay(ig,l-1)-pplev(ig,l))
            z2=(pplev(ig,l)+pplay(ig,l))/(pplev(ig,l)-pplay(ig,l))
            zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
         enddo
      enddo     

      ! Compute potential temperature
      ! Note : Potential temperature calculation may not be the same in physiq and dynamic...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do l=1,nlayer         
         do ig=1,ngrid
            zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
            mass(ig,l)  = (pplev(ig,l) - pplev(ig,l+1))/glat(ig)
            massarea(ig,l)=mass(ig,l)*cell_area(ig)
         enddo
      enddo

     ! Compute vertical velocity (m/s) from vertical mass flux
     ! w = F / (rho*area) and rho = P/(r*T)
     ! But first linearly interpolate mass flux to mid-layers
      do l=1,nlayer-1
         pw(1:ngrid,l)=0.5*(flxw(1:ngrid,l)+flxw(1:ngrid,l+1))
      enddo
      pw(1:ngrid,nlayer)=0.5*flxw(1:ngrid,nlayer) ! since flxw(nlayer+1)=0
      do l=1,nlayer
         pw(1:ngrid,l)=(pw(1:ngrid,l)*r*pt(1:ngrid,l)) /  &
                       (pplay(1:ngrid,l)*cell_area(1:ngrid))
      enddo
      ! omega in Pa/s
      do l=1,nlayer-1
         omega(1:ngrid,l)=0.5*(flxw(1:ngrid,l)+flxw(1:ngrid,l+1))
      enddo
      omega(1:ngrid,nlayer)=0.5*flxw(1:ngrid,nlayer) ! since flxw(nlayer+1)=0
      do l=1,nlayer
         omega(1:ngrid,l)=g*omega(1:ngrid,l)/cell_area(1:ngrid)
      enddo

      ! ----------------------------------------------------------------
      ! Compute mean mass, cp, and R
      ! --------------------------------
      if(photochem) then
         call concentrations(ngrid,nlayer,nq,pplay,pt,pdt,pq,pdq,ptimestep)
      endif 

!---------------------------------
! II. Compute radiative tendencies
!---------------------------------

      if (callrad) then
         if( mod(icount-1,iradia).eq.0.or.lastcall) then

          ! Compute local stellar zenith angles
          ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (tlocked) then
            ! JL14 corrects tidally resonant (and inclined) cases. nres=omega_rot/omega_orb
               ztim1=SIN(declin)
               ztim2=COS(declin)*COS(zlss)
               ztim3=COS(declin)*SIN(zlss)

               call stelang(ngrid,sinlon,coslon,sinlat,coslat,    &
                            ztim1,ztim2,ztim3,mu0,fract, flatten)

            elseif (diurnal) then
               ztim1=SIN(declin)
               ztim2=COS(declin)*COS(2.*pi*(zday-.5))
               ztim3=-COS(declin)*SIN(2.*pi*(zday-.5))

               call stelang(ngrid,sinlon,coslon,sinlat,coslat,    &
                            ztim1,ztim2,ztim3,mu0,fract, flatten)
            else if(diurnal .eqv. .false.) then

               call mucorr(ngrid,declin,latitude,mu0,fract,10000.,rad,flatten)
               ! WARNING: this function appears not to work in 1D

            endif 
           
            ! Eclipse incoming sunlight (e.g. Saturn ring shadowing).       
            if(rings_shadow) then
                call call_rings(ngrid, ptime, pday, diurnal)
            endif    


            if (corrk) then
            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! II.a Call correlated-k radiative transfer scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if(kastprof)then
                  print*,'kastprof should not = true here'
                  call abort
               endif
               if(water) then
                  muvar(1:ngrid,1:nlayer)=mugaz/(1.e0+(1.e0/epsi-1.e0)*pq(1:ngrid,1:nlayer,igcm_h2o_vap)) 
                  muvar(1:ngrid,nlayer+1)=mugaz/(1.e0+(1.e0/epsi-1.e0)*pq(1:ngrid,nlayer,igcm_h2o_vap)) 
                  ! take into account water effect on mean molecular weight
               else
                  muvar(1:ngrid,1:nlayer+1)=mugaz
               endif 
      

               if(ok_slab_ocean) then
                  tsurf2(:)=tsurf(:) 
                  do ig=1,ngrid
                     if (nint(rnat(ig))==0) then
                        tsurf(ig)=((1.-pctsrf_sic(ig))*tslab(ig,1)**4+pctsrf_sic(ig)*tsea_ice(ig)**4)**0.25
                     endif
                  enddo
               endif !(ok_slab_ocean)
                
               ! standard callcorrk
               clearsky=.false.
               call callcorrk(ngrid,nlayer,pq,nq,qsurf,                           &
                              albedo,albedo_equivalent,emis,mu0,pplev,pplay,pt,   &
                              tsurf,fract,dist_star,aerosol,muvar,                &
                              zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,                &
                              fluxsurfabs_sw,fluxtop_lw,                          &
                              fluxabs_sw,fluxtop_dn,OLR_nu,OSR_nu,                &
                              int_dtaui,int_dtauv,                                &
                              tau_col,cloudfrac,totcloudfrac,                     &
                              clearsky,firstcall,lastcall)     

                ! Option to call scheme once more for clear regions  
               if(CLFvarying)then

                  ! ---> PROBLEMS WITH ALLOCATED ARRAYS : temporary solution in callcorrk: do not deallocate if CLFvarying ...
                  clearsky=.true.
                  call callcorrk(ngrid,nlayer,pq,nq,qsurf,                           &
                                 albedo,albedo_equivalent1,emis,mu0,pplev,pplay,pt,   &
                                 tsurf,fract,dist_star,aerosol,muvar,                &
                                 zdtlw1,zdtsw1,fluxsurf_lw1,fluxsurf_sw1,            &
                                 fluxsurfabs_sw1,fluxtop_lw1,                        &
                                 fluxabs_sw1,fluxtop_dn,OLR_nu1,OSR_nu1,             &
                                 int_dtaui1,int_dtauv1,                              &
                                 tau_col1,cloudfrac,totcloudfrac,                    &
                                 clearsky,firstcall,lastcall)
                  clearsky = .false.  ! just in case.     

                  ! Sum the fluxes and heating rates from cloudy/clear cases
                  do ig=1,ngrid
                     tf=totcloudfrac(ig)
                     ntf=1.-tf                    
                     fluxsurf_lw(ig)       = ntf*fluxsurf_lw1(ig)       + tf*fluxsurf_lw(ig)
                     fluxsurf_sw(ig)       = ntf*fluxsurf_sw1(ig)       + tf*fluxsurf_sw(ig)
                     albedo_equivalent(ig) = ntf*albedo_equivalent1(ig) + tf*albedo_equivalent(ig)
                     fluxsurfabs_sw(ig)    = ntf*fluxsurfabs_sw1(ig)    + tf*fluxsurfabs_sw(ig)
                     fluxtop_lw(ig)        = ntf*fluxtop_lw1(ig)        + tf*fluxtop_lw(ig)
                     fluxabs_sw(ig)        = ntf*fluxabs_sw1(ig)        + tf*fluxabs_sw(ig)
                     tau_col(ig)           = ntf*tau_col1(ig)           + tf*tau_col(ig)
                    
                     zdtlw(ig,1:nlayer) = ntf*zdtlw1(ig,1:nlayer) + tf*zdtlw(ig,1:nlayer)
                     zdtsw(ig,1:nlayer) = ntf*zdtsw1(ig,1:nlayer) + tf*zdtsw(ig,1:nlayer)

                     OSR_nu(ig,1:L_NSPECTV) = ntf*OSR_nu1(ig,1:L_NSPECTV) + tf*OSR_nu(ig,1:L_NSPECTV)                       
                     OLR_nu(ig,1:L_NSPECTI) = ntf*OLR_nu1(ig,1:L_NSPECTI) + tf*OLR_nu(ig,1:L_NSPECTI)                       
                     int_dtauv(ig,:,1:L_NSPECTV) = ntf*int_dtauv1(ig,:,1:L_NSPECTV) + tf*int_dtauv(ig,:,1:L_NSPECTV)                       
                     int_dtaui(ig,:,1:L_NSPECTI) = ntf*int_dtaui1(ig,:,1:L_NSPECTI) + tf*int_dtaui(ig,:,1:L_NSPECTI)                       
                  enddo                               

               endif ! end of CLFvarying.

               if(ok_slab_ocean) then
                  tsurf(:)=tsurf2(:)
               endif


               ! Radiative flux from the sky absorbed by the surface (W.m-2).
               GSR=0.0
               fluxrad_sky(1:ngrid)=emis(1:ngrid)*fluxsurf_lw(1:ngrid)+fluxsurfabs_sw(1:ngrid)

                            !if(noradsurf)then ! no lower surface; SW flux just disappears
                            !   GSR = SUM(fluxsurf_sw(1:ngrid)*cell_area(1:ngrid))/totarea
                            !   fluxrad_sky(1:ngrid)=emis(1:ngrid)*fluxsurf_lw(1:ngrid)
                            !   print*,'SW lost in deep atmosphere = ',GSR,' W m^-2'
                            !endif

               ! Net atmospheric radiative heating rate (K.s-1)
               dtrad(1:ngrid,1:nlayer)=zdtsw(1:ngrid,1:nlayer)+zdtlw(1:ngrid,1:nlayer)
               
               ! Late initialization of the Ice Spectral Albedo. We needed the visible bands to do that ! 
               if (firstcall .and. albedo_spectral_mode) then
                  call spectral_albedo_calc(albedo_snow_SPECTV)
               endif

            elseif(newtonian)then
            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! II.b Call Newtonian cooling scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               call newtrelax(ngrid,nlayer,mu0,sinlat,zpopsk,pt,pplay,pplev,dtrad,firstcall)

               zdtsurf(1:ngrid) = +(pt(1:ngrid,1)-tsurf(1:ngrid))/ptimestep
               ! e.g. surface becomes proxy for 1st atmospheric layer ?

            else
            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! II.c Atmosphere has no radiative effect 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               fluxtop_dn(1:ngrid)  = fract(1:ngrid)*mu0(1:ngrid)*Fat1AU/dist_star**2
               if(ngrid.eq.1)then ! / by 4 globally in 1D case...
                  fluxtop_dn(1)  = fract(1)*Fat1AU/dist_star**2/2.0
               endif
               fluxsurf_sw(1:ngrid) = fluxtop_dn(1:ngrid)
               print*,'------------WARNING---WARNING------------' ! by MT2015.
               print*,'You are in corrk=false mode, '
               print*,'and the surface albedo is taken equal to the first visible spectral value'               
               
               fluxsurfabs_sw(1:ngrid) = fluxtop_dn(1:ngrid)*(1.-albedo(1:ngrid,1))
               fluxrad_sky(1:ngrid)    = fluxsurfabs_sw(1:ngrid)
               fluxtop_lw(1:ngrid)  = emis(1:ngrid)*sigma*tsurf(1:ngrid)**4

               dtrad(1:ngrid,1:nlayer)=0.0 ! no atmospheric radiative heating

            endif ! end of corrk

         endif ! of if(mod(icount-1,iradia).eq.0)
       

         ! Transformation of the radiative tendencies
         ! ------------------------------------------
         zplanck(1:ngrid)=tsurf(1:ngrid)*tsurf(1:ngrid)
         zplanck(1:ngrid)=emis(1:ngrid)*sigma*zplanck(1:ngrid)*zplanck(1:ngrid)
         fluxrad(1:ngrid)=fluxrad_sky(1:ngrid)-zplanck(1:ngrid)
         pdt(1:ngrid,1:nlayer)=pdt(1:ngrid,1:nlayer)+dtrad(1:ngrid,1:nlayer)
         
         ! Test of energy conservation
         !----------------------------
         if(enertest)then
            call planetwide_sumval(cpp*massarea(:,:)*zdtsw(:,:)/totarea_planet,dEtotSW)
            call planetwide_sumval(cpp*massarea(:,:)*zdtlw(:,:)/totarea_planet,dEtotLW)
            !call planetwide_sumval(fluxsurf_sw(:)*(1.-albedo_equivalent(:))*cell_area(:)/totarea_planet,dEtotsSW) !JL13 carefull, albedo can have changed since the last time we called corrk
            call planetwide_sumval(fluxsurfabs_sw(:)*cell_area(:)/totarea_planet,dEtotsSW) !JL13 carefull, albedo can have changed since the last time we called corrk
            call planetwide_sumval((fluxsurf_lw(:)*emis(:)-zplanck(:))*cell_area(:)/totarea_planet,dEtotsLW)
            dEzRadsw(:,:)=cpp*mass(:,:)*zdtsw(:,:)
            dEzRadlw(:,:)=cpp*mass(:,:)*zdtlw(:,:)
            if (is_master) then
               print*,'---------------------------------------------------------------'
               print*,'In corrk SW atmospheric heating       =',dEtotSW,' W m-2'
               print*,'In corrk LW atmospheric heating       =',dEtotLW,' W m-2'
               print*,'atmospheric net rad heating (SW+LW)   =',dEtotLW+dEtotSW,' W m-2'
               print*,'In corrk SW surface heating           =',dEtotsSW,' W m-2'
               print*,'In corrk LW surface heating           =',dEtotsLW,' W m-2'
               print*,'surface net rad heating (SW+LW)       =',dEtotsLW+dEtotsSW,' W m-2'
            endif
         endif ! end of 'enertest'

      endif ! of if (callrad)



!  --------------------------------------------
!  III. Vertical diffusion (turbulent mixing) :
!  --------------------------------------------

      if (calldifv) then
      
         zflubid(1:ngrid)=fluxrad(1:ngrid)+fluxgrd(1:ngrid)

         ! JL12 the following if test is temporarily there to allow us to compare the old vdifc with turbdiff.
         if (UseTurbDiff) then
         
            call turbdiff(ngrid,nlayer,nq,rnat,                  &
                          ptimestep,capcal,lwrite,               &
                          pplay,pplev,zzlay,zzlev,z0,            &
                          pu,pv,pt,zpopsk,pq,tsurf,emis,qsurf,   &
                          pdt,pdq,zflubid,                       &
                          zdudif,zdvdif,zdtdif,zdtsdif,          &
                          sensibFlux,q2,zdqdif,zdqevap,zdqsdif,  &
                          taux,tauy,lastcall)

         else
         
            zdh(1:ngrid,1:nlayer)=pdt(1:ngrid,1:nlayer)/zpopsk(1:ngrid,1:nlayer)
 
            call vdifc(ngrid,nlayer,nq,rnat,zpopsk,           &
                       ptimestep,capcal,lwrite,               &
                       pplay,pplev,zzlay,zzlev,z0,            &
                       pu,pv,zh,pq,tsurf,emis,qsurf,          &
                       zdh,pdq,zflubid,                       &
                       zdudif,zdvdif,zdhdif,zdtsdif,          &
                       sensibFlux,q2,zdqdif,zdqsdif,          &
                       taux,tauy,lastcall)

            zdtdif(1:ngrid,1:nlayer)=zdhdif(1:ngrid,1:nlayer)*zpopsk(1:ngrid,1:nlayer) ! for diagnostic only
            zdqevap(1:ngrid,1:nlayer)=0.

         end if !end of 'UseTurbDiff'

         zdtsurf(1:ngrid)=zdtsurf(1:ngrid)+zdtsdif(1:ngrid)

         !!! this is always done, except for turbulence-resolving simulations
         if (.not. turb_resolved) then
           pdv(1:ngrid,1:nlayer)=pdv(1:ngrid,1:nlayer)+zdvdif(1:ngrid,1:nlayer)
           pdu(1:ngrid,1:nlayer)=pdu(1:ngrid,1:nlayer)+zdudif(1:ngrid,1:nlayer)
           pdt(1:ngrid,1:nlayer)=pdt(1:ngrid,1:nlayer)+zdtdif(1:ngrid,1:nlayer)
         endif

         if(ok_slab_ocean)then
            flux_sens_lat(1:ngrid)=(zdtsdif(1:ngrid)*capcal(1:ngrid)-fluxrad(1:ngrid))
         endif


         if (tracer) then 
           pdq(1:ngrid,1:nlayer,1:nq)=pdq(1:ngrid,1:nlayer,1:nq)+ zdqdif(1:ngrid,1:nlayer,1:nq)
           dqsurf(1:ngrid,1:nq)=dqsurf(1:ngrid,1:nq) + zdqsdif(1:ngrid,1:nq)
         end if ! of if (tracer)


         ! test energy conservation
         !-------------------------
         if(enertest)then
         
            dEzdiff(:,:)=cpp*mass(:,:)*zdtdif(:,:)
            do ig = 1, ngrid
               dEdiff(ig)=SUM(dEzdiff (ig,:))+ sensibFlux(ig)! subtract flux to the ground
               dEzdiff(ig,1)= dEzdiff(ig,1)+ sensibFlux(ig)! subtract flux to the ground
            enddo
            
            call planetwide_sumval(dEdiff(:)*cell_area(:)/totarea_planet,dEtot)
            dEdiffs(:)=capcal(:)*zdtsdif(:)-zflubid(:)-sensibFlux(:)
            call planetwide_sumval(dEdiffs(:)*cell_area(:)/totarea_planet,dEtots)
            call planetwide_sumval(sensibFlux(:)*cell_area(:)/totarea_planet,AtmToSurf_TurbFlux)
            
            if (is_master) then
            
               if (UseTurbDiff) then
                         print*,'In TurbDiff sensible flux (atm=>surf) =',AtmToSurf_TurbFlux,' W m-2'
                         print*,'In TurbDiff non-cons atm nrj change   =',dEtot,' W m-2'
                  print*,'In TurbDiff (correc rad+latent heat) surf nrj change =',dEtots,' W m-2'
               else
                         print*,'In vdifc sensible flux (atm=>surf)    =',AtmToSurf_TurbFlux,' W m-2'
                         print*,'In vdifc non-cons atm nrj change      =',dEtot,' W m-2'
                         print*,'In vdifc (correc rad+latent heat) surf nrj change =',dEtots,' W m-2'
               end if
            endif ! end of 'is_master'
            
         ! JL12 : note that the black body radiative flux emitted by the surface has been updated by the implicit scheme but not given back elsewhere.
         endif ! end of 'enertest'


         ! Test water conservation.
         if(watertest.and.water)then
         
            call planetwide_sumval(massarea(:,:)*zdqdif(:,:,igcm_h2o_vap)*ptimestep/totarea_planet,dWtot_tmp)
            call planetwide_sumval(zdqsdif(:,igcm_h2o_vap)*cell_area(:)*ptimestep/totarea_planet,dWtots_tmp)
            do ig = 1, ngrid
               vdifcncons(ig)=SUM(mass(ig,:)*zdqdif(ig,:,igcm_h2o_vap))
            enddo
            call planetwide_sumval(massarea(:,:)*zdqdif(:,:,igcm_h2o_ice)*ptimestep/totarea_planet,dWtot)
            call planetwide_sumval(zdqsdif(:,igcm_h2o_ice)*cell_area(:)*ptimestep/totarea_planet,dWtots)
            dWtot = dWtot + dWtot_tmp
            dWtots = dWtots + dWtots_tmp
            do ig = 1, ngrid
               vdifcncons(ig)=vdifcncons(ig) + SUM(mass(ig,:)*zdqdif(ig,:,igcm_h2o_ice))
            enddo            
            call planetwide_maxval(vdifcncons(:),nconsMAX)

            if (is_master) then
               print*,'---------------------------------------------------------------'
               print*,'In difv atmospheric water change        =',dWtot,' kg m-2'
               print*,'In difv surface water change            =',dWtots,' kg m-2'
               print*,'In difv non-cons factor                 =',dWtot+dWtots,' kg m-2'
               print*,'In difv MAX non-cons factor             =',nconsMAX,' kg m-2 s-1'
            endif

         endif ! end of 'watertest'
         !-------------------------

      else ! calldifv

         if(.not.newtonian)then

            zdtsurf(1:ngrid) = zdtsurf(1:ngrid) + (fluxrad(1:ngrid) + fluxgrd(1:ngrid))/capcal(1:ngrid)

         endif

      endif ! end of 'calldifv'


!-------------------
!   IV. Convection :
!-------------------
      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
! IV.a Thermal plume model :
! ~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      IF (calltherm) THEN
         
! AB: We need to evaporate ice before calling thermcell_main.
         IF (water) THEN
            CALL evap(ngrid,nlayer,nq,ptimestep,pt,pq,pdq,pdt,dqevap,dtevap,zqtherm,zttherm)
         ELSE
            zttherm(:,:)   = pt(:,:)   + pdt(:,:)   * ptimestep
            zqtherm(:,:,:) = pq(:,:,:) + pdq(:,:,:) * ptimestep
         ENDIF
         
         CALL thermcell_main(ngrid, nlayer, nq, ptimestep, firstcall,            &
                             pplay, pplev, pphi, zpopsk,                         &
                             pu, pv, zttherm, zqtherm,                           &
                             zdutherm, zdvtherm, zdttherm, zdqtherm,             &
                             fm, entr, detr, zw2, fraca)
         
         pdu(1:ngrid,1:nlayer) = pdu(1:ngrid,1:nlayer) + zdutherm(1:ngrid,1:nlayer)
         pdv(1:ngrid,1:nlayer) = pdv(1:ngrid,1:nlayer) + zdvtherm(1:ngrid,1:nlayer)
         pdt(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer) + zdttherm(1:ngrid,1:nlayer)
         
         IF (tracer) THEN
            pdq(1:ngrid,1:nlayer,1:nq) = pdq(1:ngrid,1:nlayer,1:nq) + zdqtherm(1:ngrid,1:nlayer,1:nq)
         ENDIF
         
      ENDIF ! end of 'calltherm'
      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IV.b Dry convective adjustment :
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(calladj) then

         zdh(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer)/zpopsk(1:ngrid,1:nlayer)
         zduadj(1:ngrid,1:nlayer)=0.0
         zdvadj(1:ngrid,1:nlayer)=0.0
         zdhadj(1:ngrid,1:nlayer)=0.0


         call convadj(ngrid,nlayer,nq,ptimestep,            &
                      pplay,pplev,zpopsk,                   &
                      pu,pv,zh,pq,                          &
                      pdu,pdv,zdh,pdq,                      &
                      zduadj,zdvadj,zdhadj,                 &
                      zdqadj)

         pdu(1:ngrid,1:nlayer) = pdu(1:ngrid,1:nlayer) + zduadj(1:ngrid,1:nlayer)
         pdv(1:ngrid,1:nlayer) = pdv(1:ngrid,1:nlayer) + zdvadj(1:ngrid,1:nlayer)
         pdt(1:ngrid,1:nlayer)    = pdt(1:ngrid,1:nlayer) + zdhadj(1:ngrid,1:nlayer)*zpopsk(1:ngrid,1:nlayer)
         zdtadj(1:ngrid,1:nlayer) = zdhadj(1:ngrid,1:nlayer)*zpopsk(1:ngrid,1:nlayer) ! for diagnostic only

         if(tracer) then 
            pdq(1:ngrid,1:nlayer,1:nq) = pdq(1:ngrid,1:nlayer,1:nq) + zdqadj(1:ngrid,1:nlayer,1:nq) 
         end if

         ! Test energy conservation
         if(enertest)then
            call planetwide_sumval(cpp*massarea(:,:)*zdtadj(:,:)/totarea_planet,dEtot)
            if (is_master) print*,'In convadj atmospheric energy change  =',dEtot,' W m-2'
         endif

         ! Test water conservation
         if(watertest)then
            call planetwide_sumval(massarea(:,:)*zdqadj(:,:,igcm_h2o_vap)*ptimestep/totarea_planet,dWtot_tmp)
            do ig = 1, ngrid
               cadjncons(ig)=SUM(mass(ig,:)*zdqadj(ig,:,igcm_h2o_vap))
            enddo
            call planetwide_sumval(massarea(:,:)*zdqadj(:,:,igcm_h2o_ice)*ptimestep/totarea_planet,dWtot)
            dWtot = dWtot + dWtot_tmp
            do ig = 1, ngrid
               cadjncons(ig)=cadjncons(ig) + SUM(mass(ig,:)*zdqadj(ig,:,igcm_h2o_ice))
            enddo            
            call planetwide_maxval(cadjncons(:),nconsMAX)

            if (is_master) then
               print*,'In convadj atmospheric water change     =',dWtot,' kg m-2'
               print*,'In convadj MAX non-cons factor          =',nconsMAX,' kg m-2 s-1'
            endif
            
         endif ! end of 'watertest'
         
      endif ! end of 'calladj'
      
!-----------------------------------------------
!   V. Carbon dioxide condensation-sublimation :
!-----------------------------------------------

      if (co2cond) then
      
         if (.not.tracer) then
            print*,'We need a CO2 ice tracer to condense CO2'
            call abort
         endif
         call condense_co2(ngrid,nlayer,nq,ptimestep,                    &
                           capcal,pplay,pplev,tsurf,pt,                  &
                           pdt,zdtsurf,pq,pdq,                           &
                           qsurf,zdqsurfc,albedo,emis,                   &
                           albedo_bareground,albedo_co2_ice_SPECTV,      &
                           zdtc,zdtsurfc,pdpsrf,zdqc)

         pdt(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer)+zdtc(1:ngrid,1:nlayer)
         zdtsurf(1:ngrid)      = zdtsurf(1:ngrid) + zdtsurfc(1:ngrid)

         pdq(1:ngrid,1:nlayer,1:nq)   = pdq(1:ngrid,1:nlayer,1:nq)+ zdqc(1:ngrid,1:nlayer,1:nq)
         dqsurf(1:ngrid,igcm_co2_ice) = dqsurf(1:ngrid,igcm_co2_ice) + zdqsurfc(1:ngrid)

         ! test energy conservation
         if(enertest)then
            call planetwide_sumval(cpp*massarea(:,:)*zdtc(:,:)/totarea_planet,dEtot)
            call planetwide_sumval(capcal(:)*zdtsurfc(:)*cell_area(:)/totarea_planet,dEtots)
            if (is_master) then
               print*,'In co2cloud atmospheric energy change   =',dEtot,' W m-2'
               print*,'In co2cloud surface energy change       =',dEtots,' W m-2'
            endif
         endif

      endif  ! end of 'co2cond'


!---------------------------------------------
!   VI. Specific parameterizations for tracers 
!---------------------------------------------

      if (tracer) then
      
!        write(*,*) 'nlay or nlayer:', nlayer
!        write(*,*) 'zzlev', zzlev
!        -----------------------------------------
!        Volcanism : Ash source -MARS PALEOCLIMATE
!        -----------------------------------------
         if (callvolcano) then
         if (.not.tracer) then
            print*, 'We need volcanic tracers to run the volcano'
            call abort
         endif
         zdqvolc(:,:,:) = 0    
         IF (firstcall) THEN
            SELECT CASE (trim(volc_name))
               CASE ("Apollinaris_Patera")
                  lon_volc=174.4
                  lat_volc=-9.3
               CASE ("Elysium_Mons")
                  lon_volc=147.0
                  lat_volc=24.8
               CASE ("Cerberus")
                  lon_volc=176.6
                  lat_volc=9.0
               CASE ("Olympus_Mons")
                  lon_volc=-133.9
                  lat_volc=18.7
               CASE ("Arsia_Mons")
                  lon_volc=-120.46
                  lat_volc=-9.14
               CASE ("Pavonis_Mons")
                  lon_volc=-112.85
                  lat_volc=0.662
               CASE ("Ascraeus_Mons")
                  lon_volc=-104.37
                  lat_volc=11.1
               CASE ("Syrtis_Major")
                  lon_volc=66.4
                  lat_volc=9.85
               CASE ("Tyrrhenia_Patera")
                  lon_volc=106.55
                  lat_volc=-21.32
               CASE ("Hadriaca_Patera")
                  lon_volc=92.18
                  lat_volc=-30.44
               CASE ("Peneus_Patera")
                  lon_volc=60.76
                  lat_volc=-58.05
               CASE ("Alba_Patera")
                  lon_volc=-111.11
                  lat_volc=39.19
               CASE ("Amphritites")
                  lon_volc=52.66
                  lat_volc=-58.00
               CASE ("Hecates")
                  lon_volc=150.08
                  lat_volc=31.68
               CASE ("Pityusa_Patera")
                  lon_volc=36.87
                  lat_volc=-66.77
               CASE ("Malea_Patera")
                  lon_volc=50.96
                  lat_volc=-63.09
               CASE ("Electris")
                  lon_volc=-173.21
                  lat_volc=-37.35
               CASE ("Eden_Patera")
                  lon_volc=-11.1
                  lat_volc=33.6
               CASE ("Siloe_Patera")
                  lon_volc=6.6
                  lat_volc=35.2
               CASE ("Ismenia_Oxus")
                  lon_volc=0.0
                  lat_volc=38.5
               CASE Default
                  print*, 'Volcano location not found. Exiting...'
                  call abort
            END SELECT
            
            mpi_volc = mpi_size+1
            ! Lat and Lon halfspacing between gridpoints where jjm and iim are model dimensions
            lat_hs = (180./jjm)/2
            lon_hs = (360./iim)/2
            DO i=1,size(latitude_deg)
               ! Find lat,lon pair closest to grid point by comparing difference of volcano coords
               ! with half the spacing between grid points in both lat and lon (dependent on model dims)
               IF (abs(latitude_deg(i) - lat_volc) <= lat_hs .AND. &
                  abs (longitude_deg(i) - lon_volc) <= lon_hs) THEN
                  WRITE(lunout,*) 'Volcano active on proc: ', mpi_rank, "at lat ", &
                  latitude_deg(i), " lon ", longitude_deg(i), "ivolc = ", i
                  mpi_volc = mpi_rank
                  ivolc = i
                  write(lunout,*) "ivolc=",ivolc
               ELSE
                  IF (i == size(latitude_deg)) THEN
                     WRITE(lunout,*) 'Volcano inactive on proc: ', mpi_rank
                  END IF ! lat_deg
               END IF ! lat check
            END DO ! lat deg
            ! write(*,*) 'pdq(1,10,volc1)', pdq(1,10,igcm_volc_1)
            ! write(*,*) 'pdq(577,10,volc1)', pdq(577,10,igcm_volc_1)

            ! do l = 1,nlayer
            !    WRITE(*,*) mpi_rank, pdq(ivolc,l,:)
            ! end do


            ! do l = 1,nlayer
            !    WRITE(*,*) mpi_rank,'after', pdq(ivolc,l,:)
            ! end do
         
            IF (mpi_rank == mpi_volc) THEN

               WRITE(lunout,*) trim(volc_name), " selected at lat: ", lat_volc, " and lon: ", lon_volc
               filename_volc_spr = trim(datadir)//'/'//trim(volc_name)//'_spring_'//trim(atmos_type)//'_'//trim(input_key)//'.txt'
               filename_volc_sum = trim(datadir)//'/'//trim(volc_name)//'_summer_'//trim(atmos_type)//'_'//trim(input_key)//'.txt'
               filename_volc_fall = trim(datadir)//'/'//trim(volc_name)//'_fall_'//trim(atmos_type)//'_'//trim(input_key)//'.txt'
               filename_volc_win = trim(datadir)//'/'//trim(volc_name)//'_winter_'//trim(atmos_type)//'_'//trim(input_key)//'.txt'

               
               open(unit=11,file=filename_volc_spr,status='old',action='read')
               open(unit=12,file=filename_volc_sum,status='old',action='read')
               open(unit=13,file=filename_volc_fall,status='old',action='read')
               open(unit=14,file=filename_volc_win,status='old',action='read')

               nlines = 0
               ! Find length of file - Should be same for all files, so only doing once here
               DO
                  READ(11, *, IOSTAT = inputstatus)
                  IF (inputstatus > 0) STOP "*** Input error ***"
                  IF (inputstatus < 0) EXIT ! EOF
                  nlines = nlines + 1
               END DO

               nlines = nlines - 1
               REWIND 11
               ! Allocate with length of file. 5 corresponds to 1 column for height, 4 columns for tracers
               ALLOCATE(data_read_spr(5,nlines))
               ALLOCATE(data_readt_spr(nlines,5))

               ALLOCATE(data_read_sum(5,nlines))
               ALLOCATE(data_readt_sum(nlines,5))

               ALLOCATE(data_read_fall(5,nlines))
               ALLOCATE(data_readt_fall(nlines,5))

               ALLOCATE(data_read_win(5,nlines))
               ALLOCATE(data_readt_win(nlines,5))

               
               READ(11,*)
               READ(12,*)
               READ(13,*)
               READ(14,*)

               READ(11,*) data_read_spr
               data_readt_spr=transpose(data_read_spr)
               CLOSE(11)


               READ(12,*) data_read_sum
               data_readt_sum=transpose(data_read_sum)
               CLOSE(12)

               READ(13,*) data_read_fall
               data_readt_fall=transpose(data_read_fall)
               CLOSE(13)

               READ(14,*) data_read_win
               data_readt_win=transpose(data_read_win)
               CLOSE(14)
            END IF 
         zdqvolc(:,:,:) = 0
         END IF

         IF (mpi_rank == mpi_volc) THEN
            IF ((zls*180./pi) >= 0. .AND. (zls*180./pi) < 90.) THEN
               write(lunout,*) "Volcano erupting: Spring. Ls=",(zls*180/pi)
               call volcano(nq,ngrid,nlayer,ivolc,nlines,zzlev,data_readt_spr,zdqvolc)
            ELSEIF ((zls*180./pi) >= 90. .AND. (zls*180./pi) < 180.) THEN
               write(lunout,*) "Volcano erupting: Summer. Ls=",(zls*180/pi)
               call volcano(nq,ngrid,nlayer,ivolc,nlines,zzlev,data_readt_sum,zdqvolc)
            ELSEIF ((zls*180./pi) >= 180. .AND. (zls*180./pi) < 270.) THEN
               write(lunout,*) "Volcano erupting: Fall. Ls=",(zls*180/pi)
               call volcano(nq,ngrid,nlayer,ivolc,nlines,zzlev,data_readt_fall,zdqvolc)
            ELSEIF ((zls*180./pi) >= 270. .AND. (zls*180./pi) <= 360.) THEN
               write(lunout,*) "Volcano erupting: Winter. Ls=",(zls*180/pi)
               call volcano(nq,ngrid,nlayer,ivolc,nlines,zzlev,data_readt_win,zdqvolc)
            END IF
         END IF
         ! Add volcano tracers to main tracer array
         pdq(1:ngrid,1:nlayer,1:nq) = pdq(1:ngrid,1:nlayer,1:nq)       &
         + zdqvolc(1:ngrid,1:nlayer,1:nq)

           !do l = 1,nlayer
           !  WRITE(*,*) mpi_rank,'volc', zdqvolc(ivolc,l,:)
           !end do
!            pdq(1:ngrid,1:nlayer,igcm_volc_2) = pdq(1:ngrid,1:nlayer,igcm_volc_2)       &
!                                               +zdqvolc(1:ngrid,1:nlayer,igcm_volc_2)

!            pdq(1:ngrid,1:nlayer,igcm_h2o_vap) =pdq(1:ngrid,1:nlayer,igcm_h2o_vap)       &
!                                               +zdqvolc(1:ngrid,1:nlayer,igcm_h2o_vap)
 
!            write(*,*) 'pdq(1,10,volc1)', pdq(1,10,igcm_volc_1)
!            write(*,*) 'pdq(577,10,volc1)', pdq(577,10,igcm_volc_1)
!            DO iq=1,nq
!              DO ig=1, ngrid
!                DO l=1, nlayer
!                        pdq(ig,l,iq)=pdq(ig,l,iq) + zdqvolc(ig,l,iq)
!                ENDDO
!              ENDDO
!            ENDDO

!            WRITE(*,*)  'Volcano Erupting'
!            WRITE(*,*) 'zdqvolc(577,15,volc_1)', zdqvolc(577,10,2)
!            WRITE(*,*) 'zdqvolc(577,15,volc_1)',zdqvolc(577,10,igcm_volc_1)
!            WRITE(*,*) 'zdqvolc(577,15,h2o) ',zdqvolc(577,10,igcm_h2o_vap)
!            WRITE(*,*) 'pdqvolc(577,15,volc_1)', pdq(577,10,igcm_volc_1)
!          write(*,*) 'physiqvolclatlon',':LAT=',latitude(577)*180/pi
!          write(*,*)      'xLON=',longitude(577)*180/pi
         ENDIF


         !  do l = 1,nlayer
         !    WRITE(*,*) mpi_rank,'volc', zdqvolc(ivolc,l,:)
         !  end do


 ! ---------------------
  !   VI.1. Water and ice
  ! ---------------------
         if (water) then
            
            ! Water ice condensation in the atmosphere
            if(watercond.and.(RLVTT.gt.1.e-8))then
               
               if (.not.calltherm) then
                  dqmoist(1:ngrid,1:nlayer,1:nq)=0.
                  dtmoist(1:ngrid,1:nlayer)=0.
                  
                  ! Moist Convective Adjustment.
                  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call moistadj(ngrid,nlayer,nq,pt,pq,pdq,pplev,pplay,dtmoist,dqmoist,ptimestep,rneb_man)
                  
                 ! do l = 1,nlayer
                  !  WRITE(*,*) mpi_rank,'water', dqmoist(ivolc,l,:)
                  !end do

               
                  pdq(1:ngrid,1:nlayer,igcm_h2o_vap) = pdq(1:ngrid,1:nlayer,igcm_h2o_vap)     &
                                                     + dqmoist(1:ngrid,1:nlayer,igcm_h2o_vap)
                  pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = pdq(1:ngrid,1:nlayer,igcm_h2o_ice)     &
                                                     + dqmoist(1:ngrid,1:nlayer,igcm_h2o_ice)
                  pdt(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer)+dtmoist(1:ngrid,1:nlayer)
                   
                  ! Test energy conservation.
                  if(enertest)then
                     call planetwide_sumval(cpp*massarea(:,:)*dtmoist(:,:)/totarea_planet,dEtot)
                     call planetwide_maxval(dtmoist(:,:),dtmoist_max)
                     call planetwide_minval(dtmoist(:,:),dtmoist_min)
                     madjdEz(:,:)=cpp*mass(:,:)*dtmoist(:,:)
                     
                     do ig=1,ngrid
                        madjdE(ig) = cpp*SUM(mass(:,:)*dtmoist(:,:))
                     enddo
                     
                     if (is_master) then
                        print*,'In moistadj atmospheric energy change   =',dEtot,' W m-2'
                        print*,'In moistadj MAX atmospheric energy change   =',dtmoist_max*ptimestep,'K/step'
                        print*,'In moistadj MIN atmospheric energy change   =',dtmoist_min*ptimestep,'K/step'
                     endif
                     
                     call planetwide_sumval(massarea(:,:)*dqmoist(:,:,igcm_h2o_vap)*ptimestep/totarea_planet+        &
                                            massarea(:,:)*dqmoist(:,:,igcm_h2o_ice)*ptimestep/totarea_planet,dWtot)
                     if (is_master) print*,'In moistadj atmospheric water change    =',dWtot,' kg m-2'
                     
                  endif ! end of 'enertest'
               endif ! end of '.not.calltherm'
               
               ! Large scale condensation/evaporation.
               ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               call largescale(ngrid,nlayer,nq,ptimestep,pplev,pplay,pt,pq,pdt,pdq,dtlscale,dqvaplscale,dqcldlscale,rneb_lsc)

               pdt(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer)+dtlscale(1:ngrid,1:nlayer)
               pdq(1:ngrid,1:nlayer,igcm_h2o_vap) = pdq(1:ngrid,1:nlayer,igcm_h2o_vap)+dqvaplscale(1:ngrid,1:nlayer)
               pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = pdq(1:ngrid,1:nlayer,igcm_h2o_ice)+dqcldlscale(1:ngrid,1:nlayer)

               ! Test energy conservation.
               if(enertest)then
                  lscaledEz(:,:) = cpp*mass(:,:)*dtlscale(:,:)
                  do ig=1,ngrid
                     lscaledE(ig) = cpp*SUM(mass(:,:)*dtlscale(:,:))
                  enddo
                  call planetwide_sumval(cpp*massarea(:,:)*dtlscale(:,:)/totarea_planet,dEtot)

                  if (is_master) print*,'In largescale atmospheric energy change =',dEtot,' W m-2'

                  ! Test water conservation.
                  call planetwide_sumval(massarea(:,:)*dqvaplscale(:,:)*ptimestep/totarea_planet+        &
                                           SUM(massarea(:,:)*dqcldlscale(:,:))*ptimestep/totarea_planet,dWtot)
                        
                  if (is_master) print*,'In largescale atmospheric water change  =',dWtot,' kg m-2'
               endif ! end of 'enertest'

               ! Compute cloud fraction.
               do l = 1, nlayer
                  do ig=1,ngrid
                     cloudfrac(ig,l)=MAX(rneb_lsc(ig,l),rneb_man(ig,l)) 
                  enddo
               enddo

            endif ! end of 'watercond'
           
            ! Water ice / liquid precipitation.
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zdqrain(1:ngrid,1:nlayer,1:nq) = 0.0  !JL18 need to do that everytimestep if mass redis is on. 

            if(waterrain)then

               zdqsrain(1:ngrid)    = 0.0
               zdqssnow(1:ngrid)    = 0.0

               call rain(ngrid,nlayer,nq,ptimestep,pplev,pplay,pt,pdt,pq,pdq,            &
                         zdtrain,zdqrain,zdqsrain,zdqssnow,reevap_precip,cloudfrac)

               pdq(1:ngrid,1:nlayer,igcm_h2o_vap) = pdq(1:ngrid,1:nlayer,igcm_h2o_vap) &
                                                  + zdqrain(1:ngrid,1:nlayer,igcm_h2o_vap)
               pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = pdq(1:ngrid,1:nlayer,igcm_h2o_ice) &
                                                  + zdqrain(1:ngrid,1:nlayer,igcm_h2o_ice)
               pdt(1:ngrid,1:nlayer) = pdt(1:ngrid,1:nlayer)+zdtrain(1:ngrid,1:nlayer)
               
               dqsurf(1:ngrid,igcm_h2o_vap) = dqsurf(1:ngrid,igcm_h2o_vap)+zdqsrain(1:ngrid)
               dqsurf(1:ngrid,igcm_h2o_ice) = dqsurf(1:ngrid,igcm_h2o_ice)+zdqssnow(1:ngrid) 

               ! Test energy conservation.
               if(enertest)then
               
                  call planetwide_sumval(cpp*massarea(:,:)*zdtrain(:,:)/totarea_planet,dEtot)
                  if (is_master) print*,'In rain atmospheric T energy change       =',dEtot,' W m-2'
                  call planetwide_sumval(massarea(:,:)*zdqrain(:,:,igcm_h2o_ice)/totarea_planet*RLVTT/cpp,dItot_tmp)
                  call planetwide_sumval(cell_area(:)*zdqssnow(:)/totarea_planet*RLVTT/cpp,dItot)
                  dItot = dItot + dItot_tmp
                  call planetwide_sumval(massarea(:,:)*zdqrain(:,:,igcm_h2o_vap)*ptimestep/totarea_planet,dVtot_tmp)
                  call planetwide_sumval(cell_area(:)*zdqsrain(:)/totarea_planet*RLVTT/cpp,dVtot)
                  dVtot = dVtot + dVtot_tmp
                  dEtot = dItot + dVtot
                  
                  if (is_master) then
                     print*,'In rain dItot =',dItot,' W m-2'
                     print*,'In rain dVtot =',dVtot,' W m-2'
                     print*,'In rain atmospheric L energy change       =',dEtot,' W m-2'
                  endif

                  ! Test water conservation
                  call planetwide_sumval(massarea(:,:)*zdqrain(:,:,igcm_h2o_vap)*ptimestep/totarea_planet+        &
                          massarea(:,:)*zdqrain(:,:,igcm_h2o_ice)*ptimestep/totarea_planet,dWtot)
                  call planetwide_sumval((zdqsrain(:)+zdqssnow(:))*cell_area(:)*ptimestep/totarea_planet,dWtots)
                  
                  if (is_master) then
                          print*,'In rain atmospheric water change        =',dWtot,' kg m-2'
                          print*,'In rain surface water change            =',dWtots,' kg m-2'
                          print*,'In rain non-cons factor                 =',dWtot+dWtots,' kg m-2'
                  endif
                  
               endif ! end of 'enertest'

            end if ! enf of 'waterrain'
            
         end if ! end of 'water'

 ! -------------------------
  !   VI.2. Photochemistry
  ! -------------------------

         IF (photochem) then

             DO ig=1,ngrid
               array_zero1(ig)=0.0
               DO l=1,nlayer
                 array_zero2(ig,l)=0.
               ENDDO
             ENDDO

            call calchim_asis(ngrid,nlayer,nq,                          &
                        ptimestep,pplay,pplev,pt,pdt,dist_star,mu0,     &
                        fract,zzlev,zzlay,zday,pq,pdq,zdqchim,zdqschim, &
                        array_zero1,array_zero1,                        &
                        pu,pdu,pv,pdv,array_zero2,array_zero2)

           ! increment values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO l=1,nlayer
               DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+zdqchim(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO ! of DO iq=1,nq


           ! increment surface values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO ig=1,ngrid
!               dqsurf(ig,iq)=dqsurf(ig,iq)+zdqschim(ig,iq)
             ENDDO
           ENDDO ! of DO iq=1,nq

         END IF  ! of IF (photochem)


  ! -------------------------
  !   VI.3. Aerosol particles
  ! -------------------------

         ! Sedimentation.
         if (sedimentation) then
        
            zdqsed(1:ngrid,1:nlayer,1:nq) = 0.0
            zdqssed(1:ngrid,1:nq)  = 0.0

            if(watertest)then
           
               iq=igcm_h2o_ice
               call planetwide_sumval(massarea(:,:)*pq(:,:,iq)*ptimestep/totarea_planet,dWtot)
               call planetwide_sumval(massarea(:,:)*pdq(:,:,iq)*ptimestep/totarea_planet,dWtots)
               if (is_master) then
                        print*,'Before sedim pq  =',dWtot,' kg m-2'
                  print*,'Before sedim pdq =',dWtots,' kg m-2'
               endif
            endif
           
            call callsedim(ngrid,nlayer,ptimestep,           &
                          pplev,zzlev,pt,pdt,pq,pdq,        &
                          zdqsed,zdqssed,nq)

            if(watertest)then
               iq=igcm_h2o_ice
               call planetwide_sumval(massarea(:,:)*pq(:,:,iq)*ptimestep/totarea_planet,dWtot)
               call planetwide_sumval(massarea(:,:)*pdq(:,:,iq)*ptimestep/totarea_planet,dWtots)
               if (is_master) then
                        print*,'After sedim pq  =',dWtot,' kg m-2'
                        print*,'After sedim pdq =',dWtots,' kg m-2'
                 endif
            endif

            ! Whether it falls as rain or snow depends only on the surface temperature
            pdq(1:ngrid,1:nlayer,1:nq) = pdq(1:ngrid,1:nlayer,1:nq) + zdqsed(1:ngrid,1:nlayer,1:nq)
            dqsurf(1:ngrid,1:nq) = dqsurf(1:ngrid,1:nq) + zdqssed(1:ngrid,1:nq)

            ! Test water conservation
            if(watertest)then
               call planetwide_sumval(massarea(:,:)*(zdqsed(:,:,igcm_h2o_vap)+zdqsed(:,:,igcm_h2o_ice))*ptimestep/totarea_planet,dWtot)
               call planetwide_sumval((zdqssed(:,igcm_h2o_vap)+zdqssed(:,igcm_h2o_ice))*cell_area(:)*ptimestep/totarea_planet,dWtots)
               if (is_master) then
                        print*,'In sedim atmospheric ice change         =',dWtot,' kg m-2'
                        print*,'In sedim surface ice change             =',dWtots,' kg m-2'
                        print*,'In sedim non-cons factor                =',dWtot+dWtots,' kg m-2'
               endif
            endif

         end if ! end of 'sedimentation'


  ! ---------------
  !   VI.4. Updates
  ! ---------------

         ! Updating Atmospheric Mass and Tracers budgets.
         if(mass_redistrib) then

            zdmassmr(1:ngrid,1:nlayer) = mass(1:ngrid,1:nlayer) *     &
               (   zdqevap(1:ngrid,1:nlayer)                          &
                 + zdqrain(1:ngrid,1:nlayer,igcm_h2o_vap)             &
                 + dqmoist(1:ngrid,1:nlayer,igcm_h2o_vap)             &
                 + dqvaplscale(1:ngrid,1:nlayer) )

            do ig = 1, ngrid
               zdmassmr_col(ig)=SUM(zdmassmr(ig,1:nlayer))
            enddo
            
            call writediagfi(ngrid,"mass_evap","mass gain"," ",3,zdmassmr)
            call writediagfi(ngrid,"mass_evap_col","mass gain col"," ",2,zdmassmr_col)
            call writediagfi(ngrid,"mass","mass","kg/m2",3,mass)

            call mass_redistribution(ngrid,nlayer,nq,ptimestep,                     &
                                     rnat,capcal,pplay,pplev,pt,tsurf,pq,qsurf,     &
                                     pu,pv,pdt,zdtsurf,pdq,pdu,pdv,zdmassmr,        &
                                     zdtmr,zdtsurfmr,zdpsrfmr,zdumr,zdvmr,zdqmr,zdqsurfmr)
         
            pdq(1:ngrid,1:nlayer,1:nq) = pdq(1:ngrid,1:nlayer,1:nq) + zdqmr(1:ngrid,1:nlayer,1:nq)
            dqsurf(1:ngrid,1:nq)       = dqsurf(1:ngrid,1:nq) + zdqsurfmr(1:ngrid,1:nq)
            pdt(1:ngrid,1:nlayer)      = pdt(1:ngrid,1:nlayer) + zdtmr(1:ngrid,1:nlayer)
            pdu(1:ngrid,1:nlayer)      = pdu(1:ngrid,1:nlayer) + zdumr(1:ngrid,1:nlayer)
            pdv(1:ngrid,1:nlayer)      = pdv(1:ngrid,1:nlayer) + zdvmr(1:ngrid,1:nlayer)
            pdpsrf(1:ngrid)            = pdpsrf(1:ngrid) + zdpsrfmr(1:ngrid)
            zdtsurf(1:ngrid)           = zdtsurf(1:ngrid) + zdtsurfmr(1:ngrid)
            
         endif

  ! ------------------
  !   VI.5. Slab Ocean
  ! ------------------

         if (ok_slab_ocean)then

            do ig=1,ngrid
               qsurfint(:,igcm_h2o_ice)=qsurf(:,igcm_h2o_ice)
            enddo

            call ocean_slab_ice(ptimestep,                          &
                                ngrid, knindex, tsea_ice, fluxrad,  &
                                zdqssnow, qsurf(:,igcm_h2o_ice),    &
                              - zdqsdif(:,igcm_h2o_vap),            &
                                flux_sens_lat,tsea_ice, pctsrf_sic, &
                                taux,tauy,icount)


            call ocean_slab_get_vars(ngrid,tslab,      &
                                     sea_ice, flux_o,  &
                                     flux_g, dt_hdiff, &
                                     dt_ekman)
   
            do ig=1,ngrid
               if (nint(rnat(ig)).eq.1)then
                  tslab(ig,1) = 0.
                  tslab(ig,2) = 0.
                  tsea_ice(ig) = 0.
                  sea_ice(ig) = 0.
                  pctsrf_sic(ig) = 0.
                  qsurf(ig,igcm_h2o_ice) = qsurfint(ig,igcm_h2o_ice)
               endif
            enddo

         endif ! end of 'ok_slab_ocean'

  ! -----------------------------
  !   VI.6. Surface Tracer Update
  ! -----------------------------

         if(hydrology)then
            
            call hydrol(ngrid,nq,ptimestep,rnat,tsurf,qsurf,dqsurf,dqs_hyd, &
                        capcal,albedo,albedo_bareground,                    &
                        albedo_snow_SPECTV,albedo_co2_ice_SPECTV,           &
                        mu0,zdtsurf,zdtsurf_hyd,hice,pctsrf_sic,            &
                        sea_ice)

            zdtsurf(1:ngrid)     = zdtsurf(1:ngrid) + zdtsurf_hyd(1:ngrid)
            dqsurf(1:ngrid,1:nq) = dqsurf(1:ngrid,1:nq) + dqs_hyd(1:ngrid,1:nq)
            
            qsurf(1:ngrid,1:nq) = qsurf(1:ngrid,1:nq) + ptimestep*dqsurf(1:ngrid,1:nq)

            ! Test energy conservation
            if(enertest)then
               call planetwide_sumval(cell_area(:)*capcal(:)*zdtsurf_hyd(:)/totarea_planet,dEtots)
               if (is_master) print*,'In hydrol surface energy change     =',dEtots,' W m-2'
            endif

            ! test water conservation
            if(watertest)then
               call planetwide_sumval(dqs_hyd(:,igcm_h2o_ice)*cell_area(:)*ptimestep/totarea_planet,dWtots)
               if (is_master) print*,'In hydrol surface ice change            =',dWtots,' kg m-2'
                  call planetwide_sumval(dqs_hyd(:,igcm_h2o_vap)*cell_area(:)*ptimestep/totarea_planet,dWtots)
               if (is_master) then
                  print*,'In hydrol surface water change          =',dWtots,' kg m-2'
                  print*,'---------------------------------------------------------------'
               endif
            endif

         else ! of if (hydrology)

            qsurf(1:ngrid,1:nq) = qsurf(1:ngrid,1:nq) + ptimestep*dqsurf(1:ngrid,1:nq)

         end if ! of if (hydrology)

         ! Add qsurf to qsurf_hist, which is what we save in diagfi.nc. At the same time, we set the water 
         ! content of ocean gridpoints back to zero, in order to avoid rounding errors in vdifc, rain.
         qsurf_hist(:,:) = qsurf(:,:)

         if(ice_update)then
            ice_min(1:ngrid)=min(ice_min(1:ngrid),qsurf(1:ngrid,igcm_h2o_ice))
         endif

      endif! end of if 'tracer'


!------------------------------------------------
!   VII. Surface and sub-surface soil temperature 
!------------------------------------------------


      ! Increment surface temperature
      if(ok_slab_ocean)then
         do ig=1,ngrid
           if (nint(rnat(ig)).eq.0)then
             zdtsurf(ig)= (tslab(ig,1)              &
             + pctsrf_sic(ig)*(tsea_ice(ig)-tslab(ig,1))-tsurf(ig))/ptimestep
           endif
           tsurf(ig)=tsurf(ig)+ptimestep*zdtsurf(ig)
         enddo
   
      else
        tsurf(1:ngrid)=tsurf(1:ngrid)+ptimestep*zdtsurf(1:ngrid) 
      endif ! end of 'ok_slab_ocean'


      ! Compute soil temperatures and subsurface heat flux.
      if (callsoil) then
         call soil(ngrid,nsoilmx,.false.,lastcall,inertiedat,   &
                   ptimestep,tsurf,tsoil,capcal,fluxgrd)     
      endif


      if (ok_slab_ocean) then
      
         do ig=1,ngrid
         
            fluxgrdocean(ig)=fluxgrd(ig)
            if (nint(rnat(ig)).eq.0) then
               capcal(ig)=capcalocean
               fluxgrd(ig)=0.
               fluxgrdocean(ig)=pctsrf_sic(ig)*flux_g(ig)+(1-pctsrf_sic(ig))*(dt_hdiff(ig,1)+dt_ekman(ig,1))
               do iq=1,nsoilmx
                  tsoil(ig,iq)=tsurf(ig)
               enddo
               if (pctsrf_sic(ig).gt.0.5) then
                  capcal(ig)=capcalseaice
                  if (qsurf(ig,igcm_h2o_ice).gt.0.) then
                     capcal(ig)=capcalsno
                  endif
               endif               
            endif
            
         enddo
         
      endif !end of 'ok_slab_ocean'


      ! Test energy conservation
      if(enertest)then
         call planetwide_sumval(cell_area(:)*capcal(:)*zdtsurf(:)/totarea_planet,dEtots)         
         if (is_master) print*,'Surface energy change                 =',dEtots,' W m-2'
      endif


!---------------------------------------------------
!   VIII. Perform diagnostics and write output files
!---------------------------------------------------

      ! Note : For output only: the actual model integration is performed in the dynamics.


 
      ! Temperature, zonal and meridional winds.
      zt(1:ngrid,1:nlayer) = pt(1:ngrid,1:nlayer) + pdt(1:ngrid,1:nlayer)*ptimestep
      zu(1:ngrid,1:nlayer) = pu(1:ngrid,1:nlayer) + pdu(1:ngrid,1:nlayer)*ptimestep
      zv(1:ngrid,1:nlayer) = pv(1:ngrid,1:nlayer) + pdv(1:ngrid,1:nlayer)*ptimestep
      
      ! Recast thermal plume vertical velocity array for outputs
      IF (calltherm) THEN
         DO ig=1,ngrid
            DO l=1,nlayer
               zw2_bis(ig,l) = zw2(ig,l)
               fm_bis(ig,l) = fm(ig,l)
            ENDDO
         ENDDO
      ENDIF

      ! Diagnostic.
      zdtdyn(1:ngrid,1:nlayer)     = (pt(1:ngrid,1:nlayer)-ztprevious(1:ngrid,1:nlayer)) / ptimestep
      ztprevious(1:ngrid,1:nlayer) = zt(1:ngrid,1:nlayer)

      zdudyn(1:ngrid,1:nlayer)     = (pu(1:ngrid,1:nlayer)-zuprevious(1:ngrid,1:nlayer)) / ptimestep
      zuprevious(1:ngrid,1:nlayer) = zu(1:ngrid,1:nlayer)

      if(firstcall)then
         zdtdyn(1:ngrid,1:nlayer)=0.0
         zdudyn(1:ngrid,1:nlayer)=0.0
      endif

      ! Dynamical heating diagnostic.
      do ig=1,ngrid
         fluxdyn(ig)= SUM(zdtdyn(ig,:) *mass(ig,:))*cpp
      enddo

      ! Tracers.
      zq(1:ngrid,1:nlayer,1:nq) = pq(1:ngrid,1:nlayer,1:nq) + pdq(1:ngrid,1:nlayer,1:nq)*ptimestep

      ! Surface pressure.
      ps(1:ngrid) = pplev(1:ngrid,1) + pdpsrf(1:ngrid)*ptimestep



      ! Surface and soil temperature information
      call planetwide_sumval(cell_area(:)*tsurf(:)/totarea_planet,Ts1)
      call planetwide_minval(tsurf(:),Ts2)
      call planetwide_maxval(tsurf(:),Ts3)
      if(callsoil)then
         TsS = SUM(cell_area(:)*tsoil(:,nsoilmx))/totarea        ! mean temperature at bottom soil layer
         if (is_master) then
            print*,'          ave[Tsurf]             min[Tsurf]             max[Tsurf]             ave[Tdeep]'
            print*,Ts1,Ts2,Ts3,TsS
         end if
      else
         if (is_master) then
            print*,'          ave[Tsurf]             min[Tsurf]             max[Tsurf]'
            print*,Ts1,Ts2,Ts3
         endif
      end if


      ! Check the energy balance of the simulation during the run
      if(corrk)then

         call planetwide_sumval(cell_area(:)*fluxtop_dn(:)/totarea_planet,ISR)
         call planetwide_sumval(cell_area(:)*fluxabs_sw(:)/totarea_planet,ASR)
         call planetwide_sumval(cell_area(:)*fluxtop_lw(:)/totarea_planet,OLR)
         call planetwide_sumval(cell_area(:)*fluxgrd(:)/totarea_planet,GND)
         call planetwide_sumval(cell_area(:)*fluxdyn(:)/totarea_planet,DYN)
         do ig=1,ngrid
            if(fluxtop_dn(ig).lt.0.0)then
               print*,'fluxtop_dn has gone crazy'
               print*,'fluxtop_dn=',fluxtop_dn(ig)
               print*,'tau_col=',tau_col(ig)
               print*,'aerosol=',aerosol(ig,:,:)
               print*,'temp=   ',pt(ig,:)
               print*,'pplay=  ',pplay(ig,:)
               call abort
            endif
         end do
                     
         if(ngrid.eq.1)then
            DYN=0.0
         endif
         
         if (is_master) then
            print*,'  ISR            ASR            OLR            GND            DYN [W m^-2]'
            print*, ISR,ASR,OLR,GND,DYN
         endif

         if(enertest .and. is_master)then
            print*,'SW flux/heating difference SW++ - ASR = ',dEtotSW+dEtotsSW-ASR,' W m-2'
            print*,'LW flux/heating difference LW++ - OLR = ',dEtotLW+dEtotsLW+OLR,' W m-2'
            print*,'LW energy balance LW++ + ASR = ',dEtotLW+dEtotsLW+ASR,' W m-2'
         endif

         if(meanOLR .and. is_master)then
            if((ngrid.gt.1) .or. (mod(icount-1,ecritphy).eq.0))then
               ! to record global radiative balance
               open(92,file="rad_bal.out",form='formatted',position='append')
               write(92,*) zday,ISR,ASR,OLR
               close(92)
               open(93,file="tem_bal.out",form='formatted',position='append')
               if(callsoil)then
                         write(93,*) zday,Ts1,Ts2,Ts3,TsS
               else
                  write(93,*) zday,Ts1,Ts2,Ts3
               endif
               close(93)
            endif
         endif

      endif ! end of 'corrk'


      ! Diagnostic to test radiative-convective timescales in code.
      if(testradtimes)then
         open(38,file="tau_phys.out",form='formatted',position='append')
         ig=1
         do l=1,nlayer
            write(38,*) -1./pdt(ig,l),pt(ig,l),pplay(ig,l)
         enddo
         close(38)
         print*,'As testradtimes enabled,'
         print*,'exiting physics on first call'
         call abort
      endif


      ! Compute column amounts (kg m-2) if tracers are enabled.
      if(tracer)then
         qcol(1:ngrid,1:nq)=0.0
         do iq=1,nq
            do ig=1,ngrid
               qcol(ig,iq) = SUM( zq(ig,1:nlayer,iq) * mass(ig,1:nlayer))
            enddo
         enddo

         ! Generalised for arbitrary aerosols now. By LK
         reffcol(1:ngrid,1:naerkind)=0.0
         if(co2cond.and.(iaero_co2.ne.0))then
            call co2_reffrad(ngrid,nlayer,nq,zq,reffrad(1,1,iaero_co2))
            do ig=1,ngrid
               reffcol(ig,iaero_co2) = SUM(zq(ig,1:nlayer,igcm_co2_ice)*reffrad(ig,1:nlayer,iaero_co2)*mass(ig,1:nlayer))
            enddo
         endif
         if(water.and.(iaero_h2o.ne.0))then
            call h2o_reffrad(ngrid,nlayer,zq(1,1,igcm_h2o_ice),zt, &
                             reffrad(1,1,iaero_h2o),nueffrad(1,1,iaero_h2o))
            do ig=1,ngrid
               reffcol(ig,iaero_h2o) = SUM(zq(ig,1:nlayer,igcm_h2o_ice)*reffrad(ig,1:nlayer,iaero_h2o)*mass(ig,1:nlayer))
            enddo
         endif

      endif ! end of 'tracer'


      ! Test for water conservation.
      if(water)then

         call planetwide_sumval(cell_area(:)*qsurf_hist(:,igcm_h2o_ice)/totarea_planet,icesrf)
         call planetwide_sumval(cell_area(:)*qsurf_hist(:,igcm_h2o_vap)/totarea_planet,liqsrf)
         call planetwide_sumval(cell_area(:)*qcol(:,igcm_h2o_ice)/totarea_planet,icecol)
         call planetwide_sumval(cell_area(:)*qcol(:,igcm_h2o_vap)/totarea_planet,vapcol)

         h2otot = icesrf + liqsrf + icecol + vapcol
         
         if (is_master) then
            print*,' Total water amount [kg m^-2]: ',h2otot
            print*,' Surface ice    Surface liq.   Atmos. con.     Atmos. vap. [kg m^-2] '
            print*, icesrf,liqsrf,icecol,vapcol
         endif

         if(meanOLR .and. is_master)then
            if((ngrid.gt.1) .or. (mod(icount-1,ecritphy).eq.0))then
               ! to record global water balance
               open(98,file="h2o_bal.out",form='formatted',position='append')
               write(98,*) zday,icesrf,liqsrf,icecol,vapcol
               close(98)
            endif
         endif

      endif


      ! Calculate RH (Relative Humidity) for diagnostic.
      if(water)then
      
         do l = 1, nlayer
            do ig=1,ngrid
               call Psat_water(zt(ig,l),pplay(ig,l),psat_tmp,qsat(ig,l))
               RH(ig,l) = zq(ig,l,igcm_h2o_vap) / qsat(ig,l)
            enddo
         enddo

         ! Compute maximum possible H2O column amount (100% saturation).
         do ig=1,ngrid
            H2Omaxcol(ig) = SUM( qsat(ig,:) * mass(ig,:))
         enddo

      endif ! end of 'water'


      if (is_master) print*,'--> Ls =',zls*180./pi
      
      
!----------------------------------------------------------------------
!        Writing NetCDF file  "RESTARTFI" at the end of the run
!----------------------------------------------------------------------

!        Note: 'restartfi' is stored just before dynamics are stored
!              in 'restart'. Between now and the writting of 'restart',
!              there will have been the itau=itau+1 instruction and
!              a reset of 'time' (lastacll = .true. when itau+1= itaufin)
!              thus we store for time=time+dtvr



      if(lastcall) then
         ztime_fin = ptime + ptimestep/(float(iphysiq)*daysec)

         ! Update surface ice distribution to iterate to steady state if requested
         if(ice_update)then

            do ig=1,ngrid

               delta_ice = (qsurf(ig,igcm_h2o_ice)-ice_initial(ig))
               
               ! add multiple years of evolution
               qsurf_hist(ig,igcm_h2o_ice) = qsurf_hist(ig,igcm_h2o_ice) + delta_ice*icetstep 

               ! if ice has gone -ve, set to zero
               if(qsurf_hist(ig,igcm_h2o_ice).lt.0.0)then
                  qsurf_hist(ig,igcm_h2o_ice) = 0.0 
               endif

               ! if ice is seasonal, set to zero (NEW)
               if(ice_min(ig).lt.0.01)then
                  qsurf_hist(ig,igcm_h2o_ice) = 0.0 
               endif

            enddo

            ! enforce ice conservation
            ice_tot= SUM(qsurf_hist(:,igcm_h2o_ice)*cell_area(:) )/SUM(cell_area(:))
            qsurf_hist(:,igcm_h2o_ice) = qsurf_hist(:,igcm_h2o_ice)*(icesrf/ice_tot)

         endif
         
         if (ngrid.ne.1) then
            write(*,*)'PHYSIQ: for physdem ztime_fin =',ztime_fin
            
            call physdem1("restartfi.nc",nsoilmx,ngrid,nlayer,nq, &
                          ptimestep,ztime_fin,                    &
                          tsurf,tsoil,emis,q2,qsurf_hist,         &
                          cloudfrac,totcloudfrac,hice,            &
                          rnat,pctsrf_sic,tslab,tsea_ice,sea_ice)
         endif
         if(ok_slab_ocean) then
            call ocean_slab_final!(tslab, seaice)
         end if

    endif ! end of 'lastcall'


!-----------------------------------
!        Saving statistics :
!-----------------------------------

!    Note :("stats" stores and accumulates 8 key variables in file "stats.nc"
!           which can later be used to make the statistic files of the run:
!           "stats")          only possible in 3D runs !!!

         
      if (callstats) then

         call wstats(ngrid,"ps","Surface pressure","Pa",2,ps)
         call wstats(ngrid,"tsurf","Surface temperature","K",2,tsurf)
         call wstats(ngrid,"fluxsurf_lw",                               &
                     "Thermal IR radiative flux to surface","W.m-2",2,  &
                     fluxsurf_lw)
         call wstats(ngrid,"fluxtop_lw",                                &
                     "Thermal IR radiative flux to space","W.m-2",2,    &
                     fluxtop_lw)
                     
!            call wstats(ngrid,"fluxsurf_sw",                               &
!                        "Solar radiative flux to surface","W.m-2",2,       &
!                         fluxsurf_sw_tot)                     
!            call wstats(ngrid,"fluxtop_sw",                                &
!                        "Solar radiative flux to space","W.m-2",2,         &
!                        fluxtop_sw_tot)


         call wstats(ngrid,"ISR","incoming stellar rad.","W m-2",2,fluxtop_dn)
         call wstats(ngrid,"ASR","absorbed stellar rad.","W m-2",2,fluxabs_sw)
         call wstats(ngrid,"OLR","outgoing longwave rad.","W m-2",2,fluxtop_lw)
         !call wstats(ngrid,"ALB","Surface albedo"," ",2,albedo_equivalent)
         !call wstats(ngrid,"ALB_1st","First Band Surface albedo"," ",2,albedo(:,1))
         call wstats(ngrid,"p","Pressure","Pa",3,pplay)
         call wstats(ngrid,"temp","Atmospheric temperature","K",3,zt)
         call wstats(ngrid,"u","Zonal (East-West) wind","m.s-1",3,zu)
         call wstats(ngrid,"v","Meridional (North-South) wind","m.s-1",3,zv)
         call wstats(ngrid,"w","Vertical (down-up) wind","m.s-1",3,pw)
         call wstats(ngrid,"q2","Boundary layer eddy kinetic energy","m2.s-2",3,q2)

         if (tracer) then
            do iq=1,nq
               call wstats(ngrid,noms(iq),noms(iq),'kg/kg',3,zq(1,1,iq))
               call wstats(ngrid,trim(noms(iq))//'_surf',trim(noms(iq))//'_surf',  & 
                           'kg m^-2',2,qsurf(1,iq) )
               call wstats(ngrid,trim(noms(iq))//'_col',trim(noms(iq))//'_col',    & 
                          'kg m^-2',2,qcol(1,iq) )
                          
!              call wstats(ngrid,trim(noms(iq))//'_reff',                          & 
!                          trim(noms(iq))//'_reff',                                   & 
!                          'm',3,reffrad(1,1,iq))

            end do
            
            if (water) then
               vmr=zq(1:ngrid,1:nlayer,igcm_h2o_vap)*mugaz/mmol(igcm_h2o_vap)
               call wstats(ngrid,"vmr_h2ovapor",                             &
                           "H2O vapour volume mixing ratio","mol/mol",       & 
                           3,vmr)
            endif

         endif ! end of 'tracer' 

         if(watercond.and.CLFvarying)then
            call wstats(ngrid,"rneb_man","H2O cloud fraction (conv)"," ",3,rneb_man)
            call wstats(ngrid,"rneb_lsc","H2O cloud fraction (large scale)"," ",3,rneb_lsc)
            call wstats(ngrid,"CLF","H2O cloud fraction"," ",3,cloudfrac)
            call wstats(ngrid,"CLFt","H2O column cloud fraction"," ",2,totcloudfrac)
            call wstats(ngrid,"RH","relative humidity"," ",3,RH)
         endif

         if (ok_slab_ocean) then
            call wstats(ngrid,"dt_hdiff1","dt_hdiff1","K/s",2,dt_hdiff(:,1))
            call wstats(ngrid,"dt_hdiff2","dt_hdiff2","K/s",2,dt_hdiff(:,2))
            call wstats(ngrid,"dt_ekman1","dt_ekman1","K/s",2,dt_ekman(:,1))
            call wstats(ngrid,"dt_ekman2","dt_ekman2","K/s",2,dt_ekman(:,2))
            call wstats(ngrid,"tslab1","tslab1","K",2,tslab(:,1))
            call wstats(ngrid,"tslab2","tslab2","K",2,tslab(:,2))
            call wstats(ngrid,"pctsrf_sic","pct ice/sea","",2,pctsrf_sic)
            call wstats(ngrid,"tsea_ice","tsea_ice","K",2,tsea_ice)
            call wstats(ngrid,"sea_ice","sea ice","kg/m2",2,sea_ice)
            call wstats(ngrid,"rnat","nature of the surface","",2,rnat)
         endif

         if(lastcall) then
            write (*,*) "Writing stats..."
            call mkstats(ierr)
         endif
         
      endif ! end of 'callstats'

       
!-----------------------------------------------------------------------------------------------------
!           OUTPUT in netcdf file "DIAGFI.NC", containing any variable for diagnostic
!
!             Note 1 : output with  period "ecritphy", set in "run.def"
!
!             Note 2 : writediagfi can also be called from any other subroutine for any variable,
!                      but its preferable to keep all the calls in one place ...
!-----------------------------------------------------------------------------------------------------

      call writediagfi(ngrid,"Ls","solar longitude","deg",0,zls*180./pi)
      call writediagfi(ngrid,"Lss","sub solar longitude","deg",0,zlss*180./pi)
      call writediagfi(ngrid,"RA","right ascension","deg",0,right_ascen*180./pi)
      call writediagfi(ngrid,"Declin","solar declination","deg",0,declin*180./pi)
      call writediagfi(ngrid,"tsurf","Surface temperature","K",2,tsurf)
      call writediagfi(ngrid,"ps","Surface pressure","Pa",2,ps)
      call writediagfi(ngrid,"temp","temperature","K",3,zt)
      call writediagfi(ngrid,"teta","potential temperature","K",3,zh)
      call writediagfi(ngrid,"u","Zonal wind","m.s-1",3,zu)
      call writediagfi(ngrid,"v","Meridional wind","m.s-1",3,zv)
      call writediagfi(ngrid,"w","Vertical wind","m.s-1",3,pw)
      call writediagfi(ngrid,"p","Pressure","Pa",3,pplay)

!     Subsurface temperatures
!        call writediagsoil(ngrid,"tsurf","Surface temperature","K",2,tsurf)
!        call writediagsoil(ngrid,"temp","temperature","K",3,tsoil)

      ! Oceanic layers
      if(ok_slab_ocean) then
         call writediagfi(ngrid,"pctsrf_sic","pct ice/sea","",2,pctsrf_sic)
         call writediagfi(ngrid,"tsea_ice","tsea_ice","K",2,tsea_ice)
         call writediagfi(ngrid,"sea_ice","sea ice","kg/m2",2,sea_ice)
         call writediagfi(ngrid,"tslab1","tslab1","K",2,tslab(:,1))
         call writediagfi(ngrid,"tslab2","tslab2","K",2,tslab(:,2))
         call writediagfi(ngrid,"dt_hdiff1","dt_hdiff1","K/s",2,dt_hdiff(:,1))
         call writediagfi(ngrid,"dt_hdiff2","dt_hdiff2","K/s",2,dt_hdiff(:,2))
         call writediagfi(ngrid,"dt_ekman1","dt_ekman1","K/s",2,dt_ekman(:,1))
         call writediagfi(ngrid,"dt_ekman2","dt_ekman2","K/s",2,dt_ekman(:,2))
         call writediagfi(ngrid,"rnat","nature of the surface","",2,rnat)
         call writediagfi(ngrid,"sensibFlux","sensible heat flux","w.m^-2",2,sensibFlux)
         call writediagfi(ngrid,"latentFlux","latent heat flux","w.m^-2",2,zdqsdif(:,igcm_h2o_vap)*RLVTT)
      endif
      
      ! Thermal plume model
      if (calltherm) then
         call writediagfi(ngrid,'entr','Entrainment','kg m$^{-2}$ s$^{-1}$', 3, entr)
         call writediagfi(ngrid,'detr','Detrainment','kg m$^{-2}$ s$^{-1}$', 3, detr)
         call writediagfi(ngrid,'fm','Mass flux','kg m$^{-2}$ s$^{-1}$', 3, fm_bis)
         call writediagfi(ngrid,'w_plm','Squared vertical velocity','m s$^{-1}$', 3, zw2_bis)
         call writediagfi(ngrid,'fraca','Updraft fraction','', 3, fraca)
      endif
      
      ! Total energy balance diagnostics
      if(callrad.and.(.not.newtonian))then
      
         !call writediagfi(ngrid,"ALB","Surface albedo"," ",2,albedo_equivalent)
         !call writediagfi(ngrid,"ALB_1st","First Band Surface albedo"," ",2,albedo(:,1))
         call writediagfi(ngrid,"ISR","incoming stellar rad.","W m-2",2,fluxtop_dn)
         call writediagfi(ngrid,"ASR","absorbed stellar rad.","W m-2",2,fluxabs_sw)
         call writediagfi(ngrid,"OLR","outgoing longwave rad.","W m-2",2,fluxtop_lw)
         call writediagfi(ngrid,"shad","rings"," ", 2, fract)
         
!           call writediagfi(ngrid,"ASRcs","absorbed stellar rad (cs).","W m-2",2,fluxabs_sw1)
!           call writediagfi(ngrid,"OLRcs","outgoing longwave rad (cs).","W m-2",2,fluxtop_lw1)
!           call writediagfi(ngrid,"fluxsurfsw","sw surface flux.","W m-2",2,fluxsurf_sw)
!           call writediagfi(ngrid,"fluxsurflw","lw back radiation.","W m-2",2,fluxsurf_lw)
!           call writediagfi(ngrid,"fluxsurfswcs","sw surface flux (cs).","W m-2",2,fluxsurf_sw1)
!           call writediagfi(ngrid,"fluxsurflwcs","lw back radiation (cs).","W m-2",2,fluxsurf_lw1)

         if(ok_slab_ocean) then
            call writediagfi(ngrid,"GND","heat flux from ground","W m-2",2,fluxgrdocean)
         else
            call writediagfi(ngrid,"GND","heat flux from ground","W m-2",2,fluxgrd)
         endif
         
         call writediagfi(ngrid,"DYN","dynamical heat input","W m-2",2,fluxdyn)
         
      endif ! end of 'callrad'
        
      if(enertest) then
      
         if (calldifv) then
         
            call writediagfi(ngrid,"q2","turbulent kinetic energy","J.kg^-1",3,q2)
            call writediagfi(ngrid,"sensibFlux","sensible heat flux","w.m^-2",2,sensibFlux)
            
!             call writediagfi(ngrid,"dEzdiff","turbulent diffusion heating (-sensible flux)","w.m^-2",3,dEzdiff)
!             call writediagfi(ngrid,"dEdiff","integrated turbulent diffusion heating (-sensible flux)","w.m^-2",2,dEdiff)
!             call writediagfi(ngrid,"dEdiffs","In TurbDiff (correc rad+latent heat) surf nrj change","w.m^-2",2,dEdiffs)
             
         endif
          
         if (corrk) then
            call writediagfi(ngrid,"dEzradsw","radiative heating","w.m^-2",3,dEzradsw)
            call writediagfi(ngrid,"dEzradlw","radiative heating","w.m^-2",3,dEzradlw)
         endif
          
         if(watercond) then

            call writediagfi(ngrid,"lscaledE","heat from largescale","W m-2",2,lscaledE) 
            call writediagfi(ngrid,"madjdE","heat from moistadj","W m-2",2,madjdE)
            call writediagfi(ngrid,"qsatatm","atm qsat"," ",3,qsat)
             
!             call writediagfi(ngrid,"lscaledEz","heat from largescale","W m-2",3,lscaledEz) 
!             call writediagfi(ngrid,"madjdEz","heat from moistadj","W m-2",3,madjdEz)             
!             call writediagfi(ngrid,"h2o_max_col","maximum H2O column amount","kg.m^-2",2,H2Omaxcol)

         endif
          
      endif ! end of 'enertest'

        ! Diagnostics of optical thickness
        ! Warning this is exp(-tau), I let you postproc with -log to have tau itself - JVO 19
        if (diagdtau) then                
          do nw=1,L_NSPECTV
            write(str2,'(i2.2)') nw
            call writediagfi(ngrid,'dtauv'//str2,'Layer optical thickness attenuation in VI band '//str2,'',1,int_dtauv(:,nlayer:1:-1,nw))
          enddo
          do nw=1,L_NSPECTI
            write(str2,'(i2.2)') nw
            call writediagfi(ngrid,'dtaui'//str2,'Layer optical thickness attenuation in IR band '//str2,'',1,int_dtaui(:,nlayer:1:-1,nw))
          enddo
        endif


        ! Temporary inclusions for heating diagnostics.
        !call writediagfi(ngrid,"zdtsw","SW heating","T s-1",3,zdtsw)
        !call writediagfi(ngrid,"zdtlw","LW heating","T s-1",3,zdtlw)
        !call writediagfi(ngrid,"dtrad","radiative heating","K s-1",3,dtrad)
        ! call writediagfi(ngrid,"zdtdyn","Dyn. heating","T s-1",3,zdtdyn)
        
        ! For Debugging.
        !call writediagfi(ngrid,'rnat','Terrain type',' ',2,real(rnat))
        !call writediagfi(ngrid,'pphi','Geopotential',' ',3,pphi)
        

      ! Output aerosols.
      if (igcm_co2_ice.ne.0.and.iaero_co2.ne.0) &
         call writediagfi(ngrid,'CO2ice_reff','CO2ice_reff','m',3,reffrad(1,1,iaero_co2))
      if (igcm_h2o_ice.ne.0.and.iaero_h2o.ne.0) &
         call writediagfi(ngrid,'H2Oice_reff','H2Oice_reff','m',3,reffrad(:,:,iaero_h2o))
      if (igcm_co2_ice.ne.0.and.iaero_co2.ne.0) &
         call writediagfi(ngrid,'CO2ice_reffcol','CO2ice_reffcol','um kg m^-2',2,reffcol(1,iaero_co2))
      if (igcm_h2o_ice.ne.0.and.iaero_h2o.ne.0) &
         call writediagfi(ngrid,'H2Oice_reffcol','H2Oice_reffcol','um kg m^-2',2,reffcol(1,iaero_h2o))

      ! Output tracers.
      if (tracer) then
      
         do iq=1,nq
            call writediagfi(ngrid,noms(iq),noms(iq),'kg/kg',3,zq(1,1,iq))
            call writediagfi(ngrid,trim(noms(iq))//'_surf',trim(noms(iq))//'_surf',  & 
                             'kg m^-2',2,qsurf_hist(1,iq) )
            call writediagfi(ngrid,trim(noms(iq))//'_col',trim(noms(iq))//'_col',    & 
                            'kg m^-2',2,qcol(1,iq) )
                          
!          call writediagfi(ngrid,trim(noms(iq))//'_surf',trim(noms(iq))//'_surf',  & 
!                          'kg m^-2',2,qsurf(1,iq) )                          

            if(watercond.or.CLFvarying)then
               call writediagfi(ngrid,"rneb_man","H2O cloud fraction (conv)"," ",3,rneb_man)
               call writediagfi(ngrid,"rneb_lsc","H2O cloud fraction (large scale)"," ",3,rneb_lsc)
               call writediagfi(ngrid,"CLF","H2O cloud fraction"," ",3,cloudfrac)
               call writediagfi(ngrid,"CLFt","H2O column cloud fraction"," ",2,totcloudfrac)
               call writediagfi(ngrid,"RH","relative humidity"," ",3,RH)
            endif

            if(waterrain)then
               call writediagfi(ngrid,"rain","rainfall","kg m-2 s-1",2,zdqsrain)
               call writediagfi(ngrid,"snow","snowfall","kg m-2 s-1",2,zdqssnow)
               call writediagfi(ngrid,"reevap","reevaporation of precipitation","kg m-2 s-1",2,reevap_precip)
            endif

            if((hydrology).and.(.not.ok_slab_ocean))then
               call writediagfi(ngrid,"hice","oceanic ice height","m",2,hice)
            endif

            if(ice_update)then
               call writediagfi(ngrid,"ice_min","min annual ice","m",2,ice_min)
               call writediagfi(ngrid,"ice_ini","initial annual ice","m",2,ice_initial)
            endif

            do ig=1,ngrid
               if(tau_col(ig).gt.1.e3)then
                  print*,'WARNING: tau_col=',tau_col(ig)
                  print*,'at ig=',ig,'in PHYSIQ'
               endif         
            end do

            call writediagfi(ngrid,"tau_col","Total aerosol optical depth","[]",2,tau_col)

         enddo ! end of 'nq' loop
         
       endif ! end of 'tracer'


      ! Output spectrum.
      if(specOLR.and.corrk)then      
         call writediagspecIR(ngrid,"OLR3D","OLR(lon,lat,band)","W/m^2/cm^-1",3,OLR_nu)
         call writediagspecVI(ngrid,"OSR3D","OSR(lon,lat,band)","W/m^2/cm^-1",3,OSR_nu)
      endif


! XIOS outputs
      
      icount=icount+1
      
      end subroutine physiq
      
   end module physiq_mod
