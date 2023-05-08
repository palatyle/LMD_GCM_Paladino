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
 
      use radinc_h, only : L_NSPECTI,L_NSPECTV
      use radcommon_h, only: sigma, gzlat, gzlat_ig, Cmk, grav, BWNV
      use surfdat_h, only: phisfi, zmea, zstd, zsig, zgam, zthe
      use comchem_h, only: nkim, cnames, nlaykim_up, ykim_up, ykim_tot, botCH4
      use comdiurn_h, only: coslat, sinlat, coslon, sinlon
      use comsaison_h, only: mu0, fract, dist_star, declin, right_ascen
      use comsoil_h, only: nsoilmx, layer, mlayer, inertiedat
      use datafile_mod, only: datadir, corrkdir, banddir, haze_opt_file
      use geometry_mod, only: latitude, longitude, cell_area
      USE comgeomfi_h, only: totarea, totarea_planet
      USE tracer_h
      use time_phylmdz_mod, only: ecritphy, iphysiq, nday
      use phyetat0_mod, only: phyetat0
      use phyredem, only: physdem0, physdem1
      use planetwide_mod, only: planetwide_minval,planetwide_maxval,planetwide_sumval
      use mod_phys_lmdz_para, only : is_master
      use planete_mod, only: apoastr, periastr, year_day, peri_day, &
                            obliquit, nres, z0
      use comcstfi_mod, only: pi, g, rcp, r, rad, mugaz, cpp
      use time_phylmdz_mod, only: daysec
      use logic_mod, only: moyzon_ch
      use moyzon_mod, only:  zphibar, zphisbar, zplevbar, zplaybar, &
                             zzlevbar, zzlaybar, ztfibar, zqfibar
      use callkeys_mod
      use vertical_layers_mod, only: presnivs, pseudoalt
      use ioipsl_getin_p_mod, only: getin_p
      use mod_phys_lmdz_omp_data, ONLY: is_omp_master
#ifdef CPP_XIOS      
      use xios_output_mod, only: initialize_xios_output, &
                                 update_xios_timestep, &
                                 send_xios_field
      use wxios, only: wxios_context_init, xios_context_finalize
#endif
      use MMP_OPTICS
      use muphy_diag
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
!      IV. Dry Convective adjustment :
!
!      V. Tracers
!         V.1. Microphysics
!         V.2. Chemistry
!         V.3. Updates (pressure variations, surface budget).
!         V.4. Surface Tracer Update.
!
!      VI. Surface and sub-surface soil temperature.
!
!      VII. Perform diagnostics and write output files.
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
!    pudyn(ngrid,nlayer)    \
!    pvdyn(ngrid,nlayer)     \ Dynamical temporal derivative for the
!    ptdyn(ngrid,nlayer)     / corresponding variables.
!    pqdyn(ngrid,nlayer,nq) /
!    flxw(ngrid,nlayer)      vertical mass flux (kg/s) at layer lower boundary
!
!   OUTPUT
!   ------
!
!    pdu(ngrid,nlayer)        \
!    pdv(ngrid,nlayer)         \  Temporal derivative of the corresponding
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
!           Fork for Titan : J. Vatant d'Ollone (2017)
!                + clean of all too-generic (ocean, water, co2 ...) routines
!                + Titan's chemistry
!           Microphysical moment model - J.Burgalat / J.Vatant d'Ollone (2017-2018)
!============================================================================================


!    0.  Declarations :
!    ------------------

include "netcdf.inc"

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

      real, dimension(:),allocatable,save ::  tsurf                ! Surface temperature (K).
      real, dimension(:,:),allocatable,save ::  tsoil              ! Sub-surface temperatures (K).
      real, dimension(:,:),allocatable,save :: albedo              ! Surface Spectral albedo. By MT2015.
      real, dimension(:),allocatable,save :: albedo_equivalent     ! Spectral Mean albedo.
      
!$OMP THREADPRIVATE(tsurf,tsoil,albedo,albedo_equivalent)

      real,dimension(:),allocatable,save :: albedo_bareground ! Bare Ground Albedo. By MT 2015.
      
!$OMP THREADPRIVATE(albedo_bareground)

      real,dimension(:),allocatable,save :: emis        ! Thermal IR surface emissivity.
      real,dimension(:,:),allocatable,save :: dtrad     ! Net atmospheric radiative heating rate (K.s-1).
      real,dimension(:),allocatable,save :: fluxrad_sky ! Radiative flux from sky absorbed by surface (W.m-2).
      real,dimension(:),allocatable,save :: fluxrad     ! Net radiative surface flux (W.m-2).
      real,dimension(:),allocatable,save :: capcal      ! Surface heat capacity (J m-2 K-1).
      real,dimension(:),allocatable,save :: fluxgrd     ! Surface conduction flux (W.m-2).
      real,dimension(:,:),allocatable,save :: qsurf     ! Tracer on surface (e.g. kg.m-2).
      real,dimension(:,:),allocatable,save :: q2        ! Turbulent Kinetic Energy.
      
!$OMP THREADPRIVATE(emis,dtrad,fluxrad_sky,fluxrad,capcal,fluxgrd,qsurf,q2)


! Local variables :
! -----------------

      real zh(ngrid,nlayer)               ! Potential temperature (K).
      real pw(ngrid,nlayer)               ! Vertical velocity (m/s). (NOTE : >0 WHEN DOWNWARDS !!)

      integer l,ig,ierr,iq,nw,isoil
      
      ! FOR DIAGNOSTIC :
      
      real,dimension(:),allocatable,save :: fluxsurf_lw     ! Incident Long Wave (IR) surface flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxsurf_sw     ! Incident Short Wave (stellar) surface flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxsurfabs_sw  ! Absorbed Short Wave (stellar) flux by the surface (W.m-2).
      real,dimension(:),allocatable,save :: fluxtop_lw      ! Outgoing LW (IR) flux to space (W.m-2).
      real,dimension(:),allocatable,save :: fluxabs_sw      ! Absorbed SW (stellar) flux (W.m-2).
      real,dimension(:),allocatable,save :: fluxtop_dn      ! Incoming SW (stellar) radiation at the top of the atmosphere (W.m-2).
      real,dimension(:),allocatable,save :: fluxdyn         ! Horizontal heat transport by dynamics (W.m-2).
      real,dimension(:,:),allocatable,save :: OLR_nu        ! Outgoing LW radiation in each band (Normalized to the band width (W/m2/cm-1)).
      real,dimension(:,:),allocatable,save :: OSR_nu        ! Outgoing SW radiation in each band (Normalized to the band width (W/m2/cm-1)).
      real,dimension(:,:),allocatable,save :: zdtlw         ! LW heating tendencies (K/s).
      real,dimension(:,:),allocatable,save :: zdtsw         ! SW heating tendencies (K/s).
      real,dimension(:),allocatable,save :: sensibFlux      ! Turbulent flux given by the atmosphere to the surface (W.m-2).
      real,dimension(:,:,:),allocatable,save :: int_dtauv   ! VI optical thickness of layers within narrowbands for diags ().
      real,dimension(:,:,:),allocatable,save :: int_dtaui   ! IR optical thickness of layers within narrowbands for diags ().
     
!$OMP THREADPRIVATE(fluxsurf_lw,fluxsurf_sw,fluxsurfabs_sw,fluxtop_lw,fluxabs_sw,fluxtop_dn,fluxdyn,OLR_nu,OSR_nu,&        
        !$OMP zdtlw,zdtsw,sensibFlux,int_dtauv,int_dtaui))

      real zls                       ! Solar longitude (radians).
      real zlss                      ! Sub solar point longitude (radians).
      real zday                      ! Date (time since Ls=0, calculated in sols).
      real zzlay(ngrid,nlayer)       ! Altitude at the middle of the atmospheric layers (ref : local surf).
      real zzlev(ngrid,nlayer+1)     ! Altitude at the atmospheric layer boundaries (ref : local surf).
      real zzlay_eff(ngrid,nlayer)   ! Effective altitude at the middle of the atmospheric layers (ref : geoid ).
      real zzlev_eff(ngrid,nlayer+1) ! Effective altitude at the atmospheric layer boundaries ( ref : geoid ).
      
! TENDENCIES due to various processes :

      ! For Surface Temperature : (K/s)
      real zdtsurf(ngrid)     ! Cumulated tendencies.
      real zdtsurfmr(ngrid)   ! Mass_redistribution routine.
      real zdtsdif(ngrid)     ! Turbdiff/vdifc routines.
            
      ! For Atmospheric Temperatures : (K/s)    
      real zdtdif(ngrid,nlayer)                               ! Turbdiff/vdifc routines.
      real zdtmr(ngrid,nlayer)                                ! Mass_redistribution routine.
      real zdtsw1(ngrid,nlayer), zdtlw1(ngrid,nlayer)         ! Callcorrk routine.
                              
      ! For Surface Tracers : (kg/m2/s)
      real dqsurf(ngrid,nq)                 ! Cumulated tendencies.
      real zdqsdif(ngrid,nq)                ! Turbdiff/vdifc routines.
      real zdqsurfmr(ngrid,nq)              ! Mass_redistribution routine.
                  
      ! For Tracers : (kg/kg_of_air/s)
      real zdqadj(ngrid,nlayer,nq)    ! Convadj routine.
      real zdqdif(ngrid,nlayer,nq)    ! Turbdiff/vdifc routines.
      real zdqevap(ngrid,nlayer)      ! Turbdiff routine.
      real zdqmr(ngrid,nlayer,nq)     ! Mass_redistribution routine.
      
      real zdqchi(ngrid,nlayer,nq)    ! Chemical tendency ( chemistry routine ).
      
      real zdqmufi(ngrid,nlayer,nq)   ! Microphysical tendency.
      
      real zdqfibar(ngrid,nlayer,nq)   ! For 2D chemistry
      real zdqmufibar(ngrid,nlayer,nq) ! For 2D chemistry
                        
      ! For Winds : (m/s/s)
      real zdvadj(ngrid,nlayer),zduadj(ngrid,nlayer) ! Convadj routine.
      real zdumr(ngrid,nlayer),zdvmr(ngrid,nlayer)   ! Mass_redistribution routine.
      real zdvdif(ngrid,nlayer),zdudif(ngrid,nlayer) ! Turbdiff/vdifc routines.
      real zdhdif(ngrid,nlayer)                      ! Turbdiff/vdifc routines.
      real zdhadj(ngrid,nlayer)                      ! Convadj routine.
     
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
      real,allocatable,dimension(:,:),save :: ztprevious ! Previous loop Atmospheric Temperature (K)                                                         ! Useful for Dynamical Heating calculation.
      real,allocatable,dimension(:,:),save :: zuprevious ! Previous loop Zonal Wind (m.s-1)                                                                  ! Useful for Zonal Wind tendency calculation.
!$OMP THREADPRIVATE(ztprevious,zuprevious)

      real zhorizwind(ngrid,nlayer) ! Horizontal Wind ( sqrt(u**+v*v))

      real vmr(ngrid,nlayer)                        ! volume mixing ratio
      real time_phys
      
      real ISR,ASR,OLR,GND,DYN,GSR,Ts1,Ts2,Ts3,TsS ! for Diagnostic.
      
!     to test energy conservation (RW+JL)
      real mass(ngrid,nlayer),massarea(ngrid,nlayer)
      real dEtot, dEtots, AtmToSurf_TurbFlux
      real,save :: dEtotSW, dEtotsSW, dEtotLW, dEtotsLW
!$OMP THREADPRIVATE(dEtotSW, dEtotsSW, dEtotLW, dEtotsLW)
      real dEzRadsw(ngrid,nlayer),dEzRadlw(ngrid,nlayer),dEzdiff(ngrid,nlayer)
      real dEdiffs(ngrid),dEdiff(ngrid)
      
!JL12 conservation test for mean flow kinetic energy has been disabled temporarily

      real dItot, dItot_tmp, dVtot, dVtot_tmp
      
      real dWtot, dWtot_tmp, dWtots, dWtots_tmp
      
      
      ! For Clear Sky Case.
      real fluxsurf_lw1(ngrid), fluxsurf_sw1(ngrid), fluxsurfabs_sw1(ngrid)  ! For SW/LW flux.
      real fluxtop_lw1(ngrid), fluxabs_sw1(ngrid)                            ! For SW/LW flux.
      real albedo_equivalent1(ngrid)                                         ! For Equivalent albedo calculation.
      real tf, ntf    

      real,allocatable,dimension(:,:),save :: qsurf_hist
!$OMP THREADPRIVATE(qsurf_hist)
   
      ! Miscellaneous :
      character(len=10) :: tmp1
      character(len=10) :: tmp2
      
      character*2 :: str2

! Local variables for Titan chemistry and microphysics
! ----------------------------------------------------
 
      real,save :: ctimestep ! Chemistry timestep (s)
!$OMP THREADPRIVATE(ctimestep)
 
      ! Chemical tracers in molar fraction 
      real, dimension(ngrid,nlayer,nkim)          :: ychim ! (mol/mol)
      real, dimension(ngrid,nlayer,nkim)          :: ychimbar ! For 2D chemistry

      ! Molar fraction tendencies ( chemistry, condensation and evaporation ) for tracers (mol/mol/s)
      real, dimension(:,:,:), allocatable, save   :: dycchi         ! NB : Only for chem tracers. Saved since chemistry is not called every step.
!$OMP THREADPRIVATE(dycchi)
      real, dimension(ngrid,nlayer,nq)            :: dyccond        ! Condensation rate. NB : for all tracers, as we want to use indx on it.
      real, dimension(ngrid,nlayer,nq)            :: dyccondbar     ! For 2D chemistry
      real, dimension(ngrid)                      :: dycevapCH4     ! Surface "pseudo-evaporation" rate (forcing constant surface humidity).

      ! Saturation profiles
      real, dimension(ngrid,nlayer,nkim)          :: ysat ! (mol/mol)

      ! Surface methane
      real, dimension(:), allocatable, save       :: tankCH4    ! Depth of surface methane tank (m)
!$OMP THREADPRIVATE(tankCH4)

      real :: i2e(ngrid,nlayer)      ! int 2 ext factor ( X.kg-1 -> X.m-3 for diags )

      real,save,dimension(:,:,:), allocatable :: tpq ! Tracers for decoupled microphysical tests ( temporary in 01/18 )
!$OMP THREADPRIVATE(tpq)
      real,dimension(ngrid,nlayer,nq) :: dtpq ! (temporary in 01/18)

      logical file_ok

!-----------------------------------------------------------------------------
    ! Interface to calmufi 
    !   --> needed in order to pass assumed-shape arrays. Otherwise we must put calmufi in a module 
    !       (to have an explicit interface generated by the compiler).
    !   Or one can put calmufi in MMP_GCM module (in muphytitan).
    INTERFACE
      SUBROUTINE calmufi(dt, plev, zlev, play, zlay, g3d, temp, pq, zdqfi, zdq)
        REAL(kind=8), INTENT(IN)                 :: dt    !! Physics timestep (s).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: plev  !! Pressure levels (Pa).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: zlev  !! Altitude levels (m).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: play  !! Pressure layers (Pa).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: zlay  !! Altitude at the center of each layer (m).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: g3d   !! Latitude-Altitude depending gravitational acceleration (m.s-2).
        REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: temp  !! Temperature at the center of each layer (K).  
        REAL(kind=8), DIMENSION(:,:,:), INTENT(IN)  :: pq    !! Tracers (\(X.kg^{-1}}\)).
        REAL(kind=8), DIMENSION(:,:,:), INTENT(IN)  :: zdqfi !! Tendency from former processes for tracers (\(X.kg^{-1}}\)).
        REAL(kind=8), DIMENSION(:,:,:), INTENT(OUT) :: zdq   !! Microphysical tendency for tracers (\(X.kg^{-1}}\)).
      END SUBROUTINE calmufi
    END INTERFACE
      
!==================================================================================================

! -----------------
! I. INITIALISATION
! -----------------

! --------------------------------
! I.1   First Call Initialisation.
! --------------------------------
      if (firstcall) then
        allocate(tpq(ngrid,nlayer,nq))
        tpq(:,:,:) = pq(:,:,:)

        ! Initialisation of nmicro as well as tracers names, indexes ...
        if (ngrid.ne.1) then ! Already done in rcm1d
           call initracer2(nq,nametrac) ! WARNING JB (27/03/2018): should be wrapped in an OMP SINGLE statement (see module notes)
        endif

        ! Allocate saved arrays.
        ALLOCATE(tsurf(ngrid))
        ALLOCATE(tsoil(ngrid,nsoilmx))    
        ALLOCATE(albedo(ngrid,L_NSPECTV))
        ALLOCATE(albedo_equivalent(ngrid))        
        ALLOCATE(albedo_bareground(ngrid))                
        ALLOCATE(emis(ngrid))   
        ALLOCATE(dtrad(ngrid,nlayer))
        ALLOCATE(fluxrad_sky(ngrid))
        ALLOCATE(fluxrad(ngrid))    
        ALLOCATE(capcal(ngrid))    
        ALLOCATE(fluxgrd(ngrid))  
        ALLOCATE(qsurf(ngrid,nq))  
        ALLOCATE(q2(ngrid,nlayer+1)) 
        ALLOCATE(tankCH4(ngrid)) 
        ALLOCATE(ztprevious(ngrid,nlayer))
        ALLOCATE(zuprevious(ngrid,nlayer))
        ALLOCATE(qsurf_hist(ngrid,nq))
        ALLOCATE(fluxsurf_lw(ngrid))
        ALLOCATE(fluxsurf_sw(ngrid))
        ALLOCATE(fluxsurfabs_sw(ngrid))
        ALLOCATE(fluxtop_lw(ngrid))
        ALLOCATE(fluxabs_sw(ngrid))
        ALLOCATE(fluxtop_dn(ngrid))
        ALLOCATE(fluxdyn(ngrid))
        ALLOCATE(OLR_nu(ngrid,L_NSPECTI))
        ALLOCATE(OSR_nu(ngrid,L_NSPECTV))
        ALLOCATE(sensibFlux(ngrid))
        ALLOCATE(zdtlw(ngrid,nlayer))
        ALLOCATE(zdtsw(ngrid,nlayer))
        ALLOCATE(int_dtaui(ngrid,nlayer,L_NSPECTI))
        ALLOCATE(int_dtauv(ngrid,nlayer,L_NSPECTV))
 
        ! This is defined in comsaison_h
        ALLOCATE(mu0(ngrid))
        ALLOCATE(fract(ngrid))            
         ! This is defined in radcommon_h
        ALLOCATE(gzlat(ngrid,nlayer))
        ALLOCATE(gzlat_ig(nlayer))
        ALLOCATE(Cmk(nlayer))          

!        Variables set to 0
!        ~~~~~~~~~~~~~~~~~~
         dtrad(:,:) = 0.D0
         fluxrad(:) = 0.D0
         zdtsw(:,:) = 0.D0
         zdtlw(:,:) = 0.D0

!        Initialize setup for correlated-k radiative transfer
!        JVO 17 : Was in callcorrk firstcall, but we need spectral intervals for microphysics.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (corrk) then
         
           call system('rm -f surf_vals_long.out')

           write( tmp1, '(i3)' ) L_NSPECTI
           write( tmp2, '(i3)' ) L_NSPECTV
           banddir=trim(adjustl(tmp1))//'x'//trim(adjustl(tmp2))
           banddir=trim(adjustl(corrkdir))//'/'//trim(adjustl(banddir))

           call setspi            ! Basic infrared properties.
           call setspv            ! Basic visible properties.
           call sugas_corrk       ! Set up gaseous absorption properties.
         
           OLR_nu(:,:) = 0.D0
           OSR_nu(:,:) = 0.D0
           
           int_dtaui(:,:,:) = 0.D0
           int_dtauv(:,:,:) = 0.D0
           
           IF (callmufi .AND. (.NOT. uncoupl_optic_haze)) THEN
             haze_opt_file=trim(datadir)//'/HAZE_OPTIC_'//trim(adjustl(tmp1))//'x'//trim(adjustl(tmp2))//'.DAT'
             inquire(file=trim(haze_opt_file),exist=file_ok)
             if(.not.file_ok) then
               write(*,*) 'The file ',TRIM(haze_opt_file),' with the haze optical properties'
               write(*,*) 'was not found by optci.F90 ! Check in ', TRIM(datadir)
               write(*,*) 'that you have the one corresponding to the given spectral resolution !!'
               write(*,*) 'Meanwhile I abort ...'
               call abort
              endif
           ENDIF
           
         endif

!        Initialize names and timestep for chemistry
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if ( callchim ) then

            if ( moyzon_ch .and. ngrid.eq.1 ) then
               print *, "moyzon_ch=",moyzon_ch," and ngrid=1"
               print *, "Please desactivate zonal mean for 1D !"
               stop
            endif

            allocate(dycchi(ngrid,nlayer,nkim)) ! only for chemical tracers
            
            ! Chemistry timestep
            ctimestep = ptimestep*REAL(ichim)

         endif

!        Initialize microphysics.
!        ~~~~~~~~~~~~~~~~~~~~~~~~

         IF ( callmufi ) THEN
           ! WARNING JB (27/03/2018): inimufi should be wrapped in an OMP SINGLE statement.
           call inimufi(ptimestep)

           ! initialize microphysics diagnostics arrays.
           call ini_diag_arrays(ngrid,nlayer,nice)

         ENDIF

#ifdef CPP_XIOS
        ! Initialize XIOS context
        write(*,*) "physiq: call wxios_context_init"
        CALL wxios_context_init
#endif

!        Read 'startfi.nc' file.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call phyetat0(startphy_file,ngrid,nlayer,"startfi.nc",0,0,nsoilmx,nq,      &
                       day_ini,time_phys,tsurf,tsoil,emis,q2,qsurf,tankCH4)                       
         if (.not.startphy_file) then
           ! additionnal "academic" initialization of physics
           if (is_master) write(*,*) "Physiq: initializing tsurf(:) to pt(:,1) !!"
           tsurf(:)=pt(:,1)
           if (is_master) write(*,*) "Physiq: initializing tsoil(:) to pt(:,1) !!"
           do isoil=1,nsoilmx
             tsoil(:,isoil)=tsurf(:)
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
         albedo(:,:)=0.D0
          albedo_bareground(:)=0.D0
         call surfini(ngrid,nq,qsurf,albedo,albedo_bareground) 
         
!        Initialize orbital calculation. 
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call iniorbit(apoastr,periastr,year_day,peri_day,obliquit)

         if(tlocked)then
            print*,'Planet is tidally locked at resonance n=',nres
            print*,'Make sure you have the right rotation rate!!!'
         endif


!        Initialize soil.
!        ~~~~~~~~~~~~~~~~
         if (callsoil) then
         
            call soil(ngrid,nsoilmx,firstcall,lastcall,inertiedat, &
                      ptimestep,tsurf,tsoil,capcal,fluxgrd)

         else ! else of 'callsoil'.
         
            print*,'WARNING! Thermal conduction in the soil turned off'
            capcal(:)=1.e6
            fluxgrd(:)=intheat
            print*,'Flux from ground = ',intheat,' W m^-2'
            
         endif ! end of 'callsoil'.
         
         icount=1
           

!        Initialize surface history variable.
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         qsurf_hist(:,:)=qsurf(:,:)

!        Initialize variable for dynamical heating and zonal wind tendency diagnostic
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ztprevious(:,:)=pt(:,:)
         zuprevious(:,:)=pu(:,:)


         if(meanOLR)then          
            call system('rm -f rad_bal.out') ! to record global radiative balance.           
            call system('rm -f tem_bal.out') ! to record global mean/max/min temperatures.            
            call system('rm -f h2o_bal.out') ! to record global hydrological balance.
         endif

         if (ngrid.ne.1) then
            ! Note : no need to create a restart file in 1d.
            call physdem0("restartfi.nc",longitude,latitude,nsoilmx,ngrid,nlayer,nq, &
                          ptimestep,pday+nday,time_phys,cell_area,          &
                          albedo_bareground,inertiedat,zmea,zstd,zsig,zgam,zthe)
         endif
         
         ! XIOS outputs
#ifdef CPP_XIOS

         write(*,*) "physiq: call initialize_xios_output"
         call initialize_xios_output(pday,ptime,ptimestep,daysec, &
                                     presnivs,pseudoalt)
#endif
         write(*,*) "physiq: end of firstcall"
      endif ! end of 'firstcall'

! ------------------------------------------------------
! I.2   Initializations done at every physical timestep:
! ------------------------------------------------------

#ifdef CPP_XIOS      
      ! update XIOS time/calendar
      call update_xios_timestep
#endif      

      ! Initialize various variables
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pdt(:,:)    = 0.D0     
      zdtsurf(:)  = 0.D0
      pdq(:,:,:)  = 0.D0
      dqsurf(:,:) = 0.D0
      pdu(:,:)    = 0.D0
      pdv(:,:)    = 0.D0
      pdpsrf(:)   = 0.D0
      zflubid(:)  = 0.D0      
      taux(:)     = 0.D0
      tauy(:)     = 0.D0

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
         gzlat(:,:) = g         
      else if (flatten .eq. 0.0 .or. J2 .eq. 0.0 .or. Rmean .eq. 0.0 .or. MassPlanet .eq. 0.0) then      
         print*,'I need values for flatten, J2, Rmean and MassPlanet to compute gzlat (else set oblate=.false.)'
         call abort
      else
         gmplanet = MassPlanet*grav*1e24
         do ig=1,ngrid
            gzlat(ig,:)= gmplanet/(Rmean**2) * (1.D0 + 0.75 *J2 - 2.0*flatten/3. + (2.*flatten - 15./4.* J2) * cos(2. * (pi/2. - latitude(ig)))) 
         end do
      endif
      
      ! Compute altitudes with the geopotential coming from the dynamics. 
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (eff_gz .eqv. .false.) then
      
        do l=1,nlayer         
           zzlay(:,l) = pphi(:,l) / gzlat(:,l) ! Reference = local surface
        enddo
        
      else ! In this case we consider variations of g with altitude
      
        do l=1,nlayer
           zzlay(:,l) = g*rad*rad / ( g*rad - ( pphi(:,l) + phisfi(:) )) - rad
           gzlat(:,l) = g*rad*rad / ( rad + zzlay(:,l) )**2
        end do
      
      endif ! if eff_gz

      zzlev(:,1)=0.
      zzlev(:,nlayer+1)=1.e7 ! Dummy top of last layer above 10000 km...
      ! JVO 19 : This altitude is indeed dummy for the GCM and fits ptop=0
      ! but for upper chemistry that's a pb -> we anyway redefine it just after ..

      do l=2,nlayer
         do ig=1,ngrid
            z1=(pplay(ig,l-1)+pplev(ig,l))/(pplay(ig,l-1)-pplev(ig,l))
            z2=(pplev(ig,l)+pplay(ig,l))/(pplev(ig,l)-pplay(ig,l))
            zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
         enddo
      enddo     

      ! Effective altitudes ( eg needed for chemistry ) with correct g, and with reference to the geoid
      ! JVO 19 : We shall always have correct altitudes in chemistry no matter what's in physics
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (moyzon_ch) then ! Zonal averages
         
         zzlaybar(1,:)=g*rad*rad/(g*rad-(zphibar(1,:)+zphisbar(1)))-rad   ! reference = geoid 
         zzlevbar(1,1)=zphisbar(1)/g       
         DO l=2,nlayer
            z1=(zplaybar(1,l-1)+zplevbar(1,l))/(zplaybar(1,l-1)-zplevbar(1,l))
            z2=(zplevbar(1,l)  +zplaybar(1,l))/(zplevbar(1,l)  -zplaybar(1,l))
            zzlevbar(1,l)=(z1*zzlaybar(1,l-1)+z2*zzlaybar(1,l))/(z1+z2)
         ENDDO
         zzlevbar(1,nlayer+1)=zzlaybar(1,nlayer)+(zzlaybar(1,nlayer)-zzlevbar(1,nlayer))

         DO ig=2,ngrid
            if (latitude(ig).ne.latitude(ig-1)) then        
               DO l=1,nlayer   
                  zzlaybar(ig,l)=g*rad*rad/(g*rad-(zphibar(ig,l)+zphisbar(ig)))-rad
               ENDDO
               zzlevbar(ig,1)=zphisbar(ig)/g
               DO l=2,nlayer
                  z1=(zplaybar(ig,l-1)+zplevbar(ig,l))/ (zplaybar(ig,l-1)-zplevbar(ig,l))
                  z2=(zplevbar(ig,l)  +zplaybar(ig,l))/(zplevbar(ig,l)  -zplaybar(ig,l))
                  zzlevbar(ig,l)=(z1*zzlaybar(ig,l-1)+z2*zzlaybar(ig,l))/(z1+z2)
               ENDDO
               zzlevbar(ig,nlayer+1)=zzlaybar(ig,nlayer)+(zzlaybar(ig,nlayer)-zzlevbar(ig,nlayer))         
            else
               zzlaybar(ig,:)=zzlaybar(ig-1,:)
               zzlevbar(ig,:)=zzlevbar(ig-1,:)
            endif         
         ENDDO

      else !  if not moyzon
      
        DO ig=1,ngrid
          DO l=1,nlayer   
            zzlay_eff(ig,l)=g*rad*rad/(g*rad-(pphi(ig,l)+phisfi(ig)))-rad ! reference = geoid
          ENDDO
          zzlev_eff(ig,1)=phisfi(ig)/g
          DO l=2,nlayer
            z1=(pplay(ig,l-1)+pplev(ig,l))/ (pplay(ig,l-1)-pplev(ig,l))
            z2=(pplev(ig,l)  +pplay(ig,l))/(pplev(ig,l)  -pplay(ig,l))
            zzlev_eff(ig,l)=(z1*zzlay_eff(ig,l-1)+z2*zzlay_eff(ig,l))/(z1+z2)
          ENDDO
          zzlev_eff(ig,nlayer+1)=zzlay_eff(ig,nlayer)+(zzlay_eff(ig,nlayer)-zzlev_eff(ig,nlayer))
        ENDDO

      endif  ! moyzon

      ! -------------------------------------------------------------------------------------
      ! Compute potential temperature
      ! Note : Potential temperature calculation may not be the same in physiq and dynamic...
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do l=1,nlayer         
         do ig=1,ngrid
            zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
            mass(ig,l)  = (pplev(ig,l) - pplev(ig,l+1))/gzlat(ig,l)
            massarea(ig,l)=mass(ig,l)*cell_area(ig)
         enddo
      enddo

     ! Compute vertical velocity (m/s) from vertical mass flux
     ! w = F / (rho*area) and rho = P/(r*T)
     ! But first linearly interpolate mass flux to mid-layers
      do l=1,nlayer-1
         pw(:,l)=0.5*(flxw(:,l)+flxw(:,l+1))
      enddo
      pw(:,nlayer)=0.5*flxw(:,nlayer) ! since flxw(nlayer+1)=0
      do l=1,nlayer
         pw(:,l)=(pw(:,l)*r*pt(:,l)) / (pplay(:,l)*cell_area(:))
      enddo

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

               call call_profilgases(nlayer)

               ! standard callcorrk
               call callcorrk(ngrid,nlayer,pq,nq,qsurf,zday,                     &
                              albedo,albedo_equivalent,emis,mu0,pplev,pplay,pt,   &
                              tsurf,fract,dist_star,                              &
                              zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,                &
                              fluxsurfabs_sw,fluxtop_lw,                          &
                              fluxabs_sw,fluxtop_dn,OLR_nu,OSR_nu,                &
                              int_dtaui,int_dtauv,lastcall)

               ! Radiative flux from the sky absorbed by the surface (W.m-2).
               GSR=0.0
               fluxrad_sky(:)=emis(:)*fluxsurf_lw(:)+fluxsurfabs_sw(:)

                            !if(noradsurf)then ! no lower surface; SW flux just disappears
                            !   GSR = SUM(fluxsurf_sw(:)*cell_area(:))/totarea
                            !   fluxrad_sky(:)=emis(:)*fluxsurf_lw(:)
                            !   print*,'SW lost in deep atmosphere = ',GSR,' W m^-2'
                            !endif

               ! Net atmospheric radiative heating rate (K.s-1)
               dtrad(:,:)=zdtsw(:,:)+zdtlw(:,:)
               
            elseif(newtonian)then
            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! II.b Call Newtonian cooling scheme
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               call newtrelax(ngrid,nlayer,mu0,sinlat,zpopsk,pt,pplay,pplev,dtrad,firstcall)

               zdtsurf(:) = +(pt(:,1)-tsurf(:))/ptimestep
               ! e.g. surface becomes proxy for 1st atmospheric layer ?

            else
            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! II.c Atmosphere has no radiative effect 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               fluxtop_dn(:)  = fract(:)*mu0(:)*Fat1AU/dist_star**2
               if(ngrid.eq.1)then ! / by 4 globally in 1D case...
                  fluxtop_dn(1)  = fract(1)*Fat1AU/dist_star**2/2.0
               endif
               fluxsurf_sw(:) = fluxtop_dn(:)
               print*,'------------WARNING---WARNING------------' ! by MT2015.
               print*,'You are in corrk=false mode, '
               print*,'and the surface albedo is taken equal to the first visible spectral value'               
               
               fluxsurfabs_sw(:) = fluxtop_dn(:)*(1.-albedo(:,1))
               fluxrad_sky(:)    = fluxsurfabs_sw(:)
               fluxtop_lw(:)     = emis(:)*sigma*tsurf(:)**4

               dtrad(:,:)=0.D0 ! no atmospheric radiative heating

            endif ! end of corrk

         endif ! of if(mod(icount-1,iradia).eq.0)
       

         ! Transformation of the radiative tendencies
         ! ------------------------------------------
         zplanck(:)=tsurf(:)*tsurf(:)
         zplanck(:)=emis(:)*sigma*zplanck(:)*zplanck(:)
         fluxrad(:)=fluxrad_sky(:)-zplanck(:)
         pdt(:,:)=pdt(:,:)+dtrad(:,:)
         
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
      
         zflubid(:)=fluxrad(:)+fluxgrd(:)

         ! JL12 the following if test is temporarily there to allow us to compare the old vdifc with turbdiff.
         if (UseTurbDiff) then
         
            call turbdiff(ngrid,nlayer,nq,                       &
                          ptimestep,capcal,lwrite,               &
                          pplay,pplev,zzlay,zzlev,z0,            &
                          pu,pv,pt,zpopsk,pq,tsurf,emis,qsurf,   &
                          pdt,pdq,zflubid,                       &
                          zdudif,zdvdif,zdtdif,zdtsdif,          &
                          sensibFlux,q2,zdqdif,zdqsdif,          &
                          taux,tauy,lastcall)

         else
         
            zdh(:,:)=pdt(:,:)/zpopsk(:,:)
 
            call vdifc(ngrid,nlayer,nq,zpopsk,                &
                       ptimestep,capcal,lwrite,               &
                       pplay,pplev,zzlay,zzlev,z0,            &
                       pu,pv,zh,pq,tsurf,emis,qsurf,          &
                       zdh,pdq,zflubid,                       &
                       zdudif,zdvdif,zdhdif,zdtsdif,          &
                       sensibFlux,q2,zdqdif,zdqsdif,          &
                       taux,tauy,lastcall)

            zdtdif(:,:)=zdhdif(:,:)*zpopsk(:,:) ! for diagnostic only
            zdqevap(:,:)=0.

         end if !end of 'UseTurbDiff'


         pdv(:,:)=pdv(:,:)+zdvdif(:,:)
         pdu(:,:)=pdu(:,:)+zdudif(:,:)
         pdt(:,:)=pdt(:,:)+zdtdif(:,:)
         zdtsurf(:)=zdtsurf(:)+zdtsdif(:)

         if (tracer) then 
           pdq(:,:,:)=pdq(:,:,:)+ zdqdif(:,:,:)
           dqsurf(:,:)=dqsurf(:,:) + zdqsdif(:,:)
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

      else ! calldifv

         if(.not.newtonian)then

            zdtsurf(:) = zdtsurf(:) + (fluxrad(:) + fluxgrd(:))/capcal(:)

         endif

      endif ! end of 'calldifv'


!----------------------------------
!   IV. Dry convective adjustment :
!----------------------------------

      if(calladj) then

         zdh(:,:) = pdt(:,:)/zpopsk(:,:)
         zduadj(:,:)=0.D0
         zdvadj(:,:)=0.D0
         zdhadj(:,:)=0.D0


         call convadj(ngrid,nlayer,nq,ptimestep,            &
                      pplay,pplev,zpopsk,                   &
                      pu,pv,zh,pq,                          &
                      pdu,pdv,zdh,pdq,                      &
                      zduadj,zdvadj,zdhadj,                 &
                      zdqadj)

         pdu(:,:) = pdu(:,:) + zduadj(:,:)
         pdv(:,:) = pdv(:,:) + zdvadj(:,:)
         pdt(:,:)    = pdt(:,:) + zdhadj(:,:)*zpopsk(:,:)
         zdtadj(:,:) = zdhadj(:,:)*zpopsk(:,:) ! for diagnostic only

         if(tracer) then 
            pdq(:,:,:) = pdq(:,:,:) + zdqadj(:,:,:) 
         end if

         ! Test energy conservation
         if(enertest)then
            call planetwide_sumval(cpp*massarea(:,:)*zdtadj(:,:)/totarea_planet,dEtot)
            if (is_master) print*,'In convadj atmospheric energy change  =',dEtot,' W m-2'
         endif

         
      endif ! end of 'calladj'
      

!---------------------------------------------
!   V. Specific parameterizations for tracers 
!---------------------------------------------

      if (tracer) then

  ! -------------------
  !   V.1 Microphysics
  ! -------------------

         ! JVO 05/18 : We must call microphysics before chemistry, for condensation !
 
         if (callmufi) then

            zzlev(:,nlayer+1)=zzlay(:,nlayer)+(zzlay(:,nlayer)-zzlev(:,nlayer)) ! JVO 19 : We assume zzlev isn't reused later on (could be done cleaner)

#ifdef USE_QTEST
            dtpq(:,:,:) = 0.D0 ! we want tpq to go only through mufi
            call calmufi(ptimestep,pplev,zzlev,pplay,zzlay,gzlat,pt,tpq,dtpq,zdqmufi)
            tpq(:,:,:) = tpq(:,:,:) + zdqmufi(:,:,:)*ptimestep ! only manipulation of tpq->*ptimestep here
#else
            
            call calmufi(ptimestep,pplev,zzlev,pplay,zzlay,gzlat,pt,pq,pdq,zdqmufi)

            pdq(:,:,:) = pdq(:,:,:) + zdqmufi(:,:,:)

            ! Sanity check ( way safer to be done here rather than within YAMMS )
            ! Important : the sanity check intentionally include the former processes tendency !
            ! NB : Despite this sanity check there might be still some unphysical values going through :
            !        - Negatives, but harmless as it will be only for the output files
            !          just remove them in post-proc.
            !        - Weird unphysical ratio of m0 and m3, ok for now, but take care of them if
            !          you want to compute optics from radii.
            WHERE ( (pq(:,:,1)+pdq(:,:,1)*ptimestep < 0.D0) .OR. (pq(:,:,2)+pdq(:,:,2)*ptimestep < 0.D0) )
                    pdq(:,:,1) = (epsilon(1.0)-1.D0)*pq(:,:,1)/ptimestep 
                    pdq(:,:,2) = (epsilon(1.0)-1.D0)*pq(:,:,2)/ptimestep
            ENDWHERE 
            WHERE ( (pq(:,:,3)+pdq(:,:,3)*ptimestep < 0.D0) .OR. (pq(:,:,4)+pdq(:,:,4)*ptimestep < 0.D0) )
                    pdq(:,:,3) = (epsilon(1.0)-1.D0)*pq(:,:,3)/ptimestep 
                    pdq(:,:,4) = (epsilon(1.0)-1.D0)*pq(:,:,4)/ptimestep
            ENDWHERE
#endif

            ! Microphysics condensation for 2D fields to sent non-saturated fields to photochem
            if ( callclouds .and. moyzon_ch .and. mod(icount-1,ichim).eq.0 ) then
              zdqfibar(:,:,:) = 0.D0 ! We work in zonal average -> forget processes other than condensation
              call calmufi(ptimestep,zplevbar,zzlevbar,zplaybar,zzlaybar, &
                           gzlat,ztfibar,zqfibar,zdqfibar,zdqmufibar)
              ! TODO : Add a sanity check here !
            endif

         endif
      
  ! -----------------
  !   V.2. Chemistry
  ! -----------------
  !   NB : Must be call last ( brings fields back to an equilibrium )

  ! Known bug ? ( JVO 18 ) : If you try to use chimi_indx instead of iq+nmicro
  ! it leads to weird results / crash on dev mod ( ok in debug ) ... Why ? Idk ...

         if (callchim) then

            ! o. Convert updated tracers to molar fraction
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do iq = 1,nkim
               ychim(:,:,iq) = ( pq(:,:,iq+nmicro) + pdq(:,:,iq+nmicro) ) / rat_mmol(iq+nmicro)
            enddo

            ! JVO 05/18 : We update zonally averaged fields with condensation 
            ! as it is compulsory to have correct photochem production. But for other 
            ! processes ( convadj ... ) we miss them in any case as we work in zonally/diurnal 
            ! mean -> no fine diurnal/short time evolution, only seasonal evolution only.
            if ( moyzon_ch .and. mod(icount-1,ichim).eq. 0 ) then
              do iq = 1,nkim
                ychimbar(:,:,iq) =  zqfibar(:,:,iq+nmicro) / rat_mmol(iq+nmicro)
                  if ( callclouds ) then
                    ychimbar(:,:,iq) =  ychimbar(:,:,iq) + ( zdqmufibar(:,:,iq+nmicro) / rat_mmol(iq+nmicro) )
                  endif
              enddo
            endif

            ! i. Condensation of the 3D tracers after the transport
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            call calc_ysat(ngrid,nlayer,pplay/100.0,pt,ysat) ! Compute saturation profiles for every grid point

            dyccond(:,:,:) = 0.D0 ! Default value -> no condensation

            do iq=1,nkim
               where ( ychim(:,:,iq).gt.ysat(:,:,iq) )   &
                     dyccond(:,:,iq+nmicro) = ( -ychim(:,:,iq)+ysat(:,:,iq) ) / ptimestep
            enddo

            if ( callclouds ) dyccond(:,:,ices_indx) = 0.D0 ! Condensation have been calculated in the cloud microphysics 

            do iq=1,nkim
              ychim(:,:,iq) = ychim(:,:,iq) + dyccond(:,:,iq+nmicro) ! update molar ychim for following calchim

              pdq(:,:,iq+nmicro) = pdq(:,:,iq+nmicro) + dyccond(:,:,iq+nmicro)*rat_mmol(iq+nmicro) ! convert tendencies to mass mixing ratio
            enddo
            

            ! ii. 2D zonally averaged fields need to condense and evap before photochem
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if ( moyzon_ch .and. mod(icount-1,ichim).eq. 0 ) then

              call calc_ysat(ngrid,nlayer,zplaybar/100.0,ztfibar,ysat) ! Compute saturation profiles for every grid point for the zon-ave fields

              dyccondbar(:,:,:) = 0.D0 ! Default value -> no condensation
              
              do iq = 1,nkim
                 where ( ychimbar(:,:,iq).gt.ysat(:,:,iq) )  &
                       dyccondbar(:,:,iq+nmicro) = ( -ychimbar(:,:,iq)+ysat(:,:,iq) ) / ptimestep
              enddo

              if ( callclouds ) dyccondbar(:,:,ices_indx) = 0.D0 ! Condensation have been calculated in the cloud microphysics 

              do iq=1,nkim
                ychimbar(:,:,iq) = ychimbar(:,:,iq) + dyccondbar(:,:,iq+nmicro)
              enddo

              ! Pseudo-evap ( forcing constant surface humidity )
              do ig=1,ngrid
                 if ( ychimbar(ig,1,7) .lt. botCH4 ) ychimbar(ig,1,7) = botCH4
              enddo

            endif ! if ( moyzon_ch .and. mod(icount-1,ichim).eq. 0 )

            ! iii. Photochemistry ( must be call after condensation (and evap of 2D) )
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if( mod(icount-1,ichim).eq.0. ) then
               
               print *, "We enter in the chemistry ..."

               if (moyzon_ch) then ! 2D zonally averaged chemistry

                 ! Here we send zonal average fields ( corrected with cond ) from dynamics to chem. module
                 call calchim(ngrid,ychimbar,declin,ctimestep,ztfibar,zphibar,zphisbar,  &
                              zplaybar,zplevbar,zzlaybar,zzlevbar,dycchi)
               else ! 3D chemistry (or 1D run)
                 call calchim(ngrid,ychim,declin,ctimestep,pt,pphi,phisfi,  &
                              pplay,pplev,zzlay_eff,zzlev_eff,dycchi)
               endif ! if moyzon

            endif
            
            ! iv. Surface pseudo-evaporation
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            do ig=1,ngrid
               if ( (ychim(ig,1,7)+dycchi(ig,1,7)*ptimestep) .lt. botCH4 ) then ! +dycchi because ychim not yet updated
                  dycevapCH4(ig) = ( -ychim(ig,1,7)+botCH4 ) / ptimestep - dycchi(ig,1,7)
               else
                  dycevapCH4(ig) = 0.D0
               endif
            enddo

            pdq(:,1,7+nmicro) = pdq(:,1,7+nmicro) + dycevapCH4(:)*rat_mmol(7+nmicro)
            
            ! v. Updates and positivity check
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zdqchi(:,:,:)  = 0.D0 ! -> dycchi is saved but for the nmicro tracers we must update to 0 at each step

            do iq=1,nkim
              zdqchi(:,:,iq+nmicro) = dycchi(:,:,iq)*rat_mmol(iq+nmicro) ! convert tendencies to mass mixing ratio

              where( (pq(:,:,iq+nmicro) + ( pdq(:,:,iq+nmicro)+zdqchi(:,:,iq+nmicro) )*ptimestep ) .LT. 0.)  & ! When using zonal means we set the same tendency
                      zdqchi(:,:,iq+nmicro) = 1.D-30 - pdq(:,:,iq+nmicro) - pq(:,:,iq+nmicro)/ptimestep        ! everywhere in longitude -> could lead to negs !
            enddo

            pdq(:,:,:) = pdq(:,:,:) + zdqchi(:,:,:)
            
         endif ! end of 'callchim'

  ! ---------------
  !   V.3 Updates
  ! ---------------

         ! Updating Atmospheric Mass and Tracers budgets.
         if(mass_redistrib) then

            zdmassmr(:,:) = mass(:,:) * zdqevap(:,:)

            do ig = 1, ngrid
               zdmassmr_col(ig)=SUM(zdmassmr(ig,:))
            enddo
            
            call writediagfi(ngrid,"mass_evap","mass gain"," ",3,zdmassmr)
            call writediagfi(ngrid,"mass_evap_col","mass gain col"," ",2,zdmassmr_col)
            call writediagfi(ngrid,"mass","mass","kg/m2",3,mass)

            call mass_redistribution(ngrid,nlayer,nq,ptimestep,                     &
                                     capcal,pplay,pplev,pt,tsurf,pq,qsurf,          &
                                     pu,pv,pdt,zdtsurf,pdq,pdu,pdv,zdmassmr,        &
                                     zdtmr,zdtsurfmr,zdpsrfmr,zdumr,zdvmr,zdqmr,zdqsurfmr)
         
            pdq(:,:,:)  = pdq(:,:,:)  + zdqmr(:,:,:)
            dqsurf(:,:) = dqsurf(:,:) + zdqsurfmr(:,:)
            pdt(:,:)    = pdt(:,:)    + zdtmr(:,:)
            pdu(:,:)    = pdu(:,:)    + zdumr(:,:)
            pdv(:,:)    = pdv(:,:)    + zdvmr(:,:)
            pdpsrf(:)   = pdpsrf(:)   + zdpsrfmr(:)
            zdtsurf(:)  = zdtsurf(:)  + zdtsurfmr(:)
            
         endif

  ! -----------------------------
  !   V.4. Surface Tracer Update
  ! -----------------------------

        qsurf(:,:) = qsurf(:,:) + ptimestep*dqsurf(:,:)

         ! Add qsurf to qsurf_hist, which is what we save in diagfi.nc. At the same time, we set the water 
         ! content of ocean gridpoints back to zero, in order to avoid rounding errors in vdifc, rain.
         qsurf_hist(:,:) = qsurf(:,:)

      endif! end of if 'tracer'


!------------------------------------------------
!   VI. Surface and sub-surface soil temperature 
!------------------------------------------------


      ! Increment surface temperature

      tsurf(:)=tsurf(:)+ptimestep*zdtsurf(:) 

      ! Compute soil temperatures and subsurface heat flux.
      if (callsoil) then
         call soil(ngrid,nsoilmx,.false.,lastcall,inertiedat,   &
                   ptimestep,tsurf,tsoil,capcal,fluxgrd)     
      endif


      ! Test energy conservation
      if(enertest)then
         call planetwide_sumval(cell_area(:)*capcal(:)*zdtsurf(:)/totarea_planet,dEtots)         
         if (is_master) print*,'Surface energy change                 =',dEtots,' W m-2'
      endif


!---------------------------------------------------
!   VII. Perform diagnostics and write output files
!---------------------------------------------------

      ! Note : For output only: the actual model integration is performed in the dynamics.
 
      ! Temperature, zonal and meridional winds.
      zt(:,:) = pt(:,:) + pdt(:,:)*ptimestep
      zu(:,:) = pu(:,:) + pdu(:,:)*ptimestep
      zv(:,:) = pv(:,:) + pdv(:,:)*ptimestep

      ! Diagnostic.
      zdtdyn(:,:)     = (pt(:,:)-ztprevious(:,:)) / ptimestep
      ztprevious(:,:) = zt(:,:)

      zdudyn(:,:)     = (pu(:,:)-zuprevious(:,:)) / ptimestep
      zuprevious(:,:) = zu(:,:)

      if(firstcall)then
         zdtdyn(:,:)=0.D0
         zdudyn(:,:)=0.D0
      endif

      ! Horizotal wind
      zhorizwind(:,:) = sqrt( zu(:,:)*zu(:,:) + zv(:,:)*zv(:,:) )

      ! Dynamical heating diagnostic.
      do ig=1,ngrid
         fluxdyn(ig)= SUM(zdtdyn(ig,:) *mass(ig,:))*cpp
      enddo

      ! Tracers.
      zq(:,:,:) = pq(:,:,:) + pdq(:,:,:)*ptimestep

      ! Surface pressure.
      ps(:) = pplev(:,1) + pdpsrf(:)*ptimestep



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

         if (ngrid.ne.1) then
            write(*,*)'PHYSIQ: for physdem ztime_fin =',ztime_fin
            
            call physdem1("restartfi.nc",nsoilmx,ngrid,nlayer,nq, &
                          ptimestep,ztime_fin,                    &
                          tsurf,tsoil,emis,q2,qsurf_hist,tankCH4)
         endif
    endif ! end of 'lastcall'


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

      

      ! Total energy balance diagnostics
      if(callrad.and.(.not.newtonian))then
      
         !call writediagfi(ngrid,"ALB","Surface albedo"," ",2,albedo_equivalent)
         call writediagfi(ngrid,"ISR","incoming stellar rad.","W m-2",2,fluxtop_dn)
         call writediagfi(ngrid,"ASR","absorbed stellar rad.","W m-2",2,fluxabs_sw)
         call writediagfi(ngrid,"OLR","outgoing longwave rad.","W m-2",2,fluxtop_lw)
         
!           call writediagfi(ngrid,"ASRcs","absorbed stellar rad (cs).","W m-2",2,fluxabs_sw1)
!           call writediagfi(ngrid,"OLRcs","outgoing longwave rad (cs).","W m-2",2,fluxtop_lw1)
!           call writediagfi(ngrid,"fluxsurfsw","sw surface flux.","W m-2",2,fluxsurf_sw)
!           call writediagfi(ngrid,"fluxsurflw","lw back radiation.","W m-2",2,fluxsurf_lw)
!           call writediagfi(ngrid,"fluxsurfswcs","sw surface flux (cs).","W m-2",2,fluxsurf_sw1)
!           call writediagfi(ngrid,"fluxsurflwcs","lw back radiation (cs).","W m-2",2,fluxsurf_lw1)


         call writediagfi(ngrid,"GND","heat flux from ground","W m-2",2,fluxgrd)
         
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

        ! Temporary inclusions for winds diagnostics.
        call writediagfi(ngrid,"zdudif","Turbdiff tend. zon. wind","m s-2",3,zdudif)
        call writediagfi(ngrid,"zdudyn","Dyn. tend. zon. wind","m s-2",3,zdudyn)

        ! Temporary inclusions for heating diagnostics.
        call writediagfi(ngrid,"zdtsw","SW heating","T s-1",3,zdtsw)
        call writediagfi(ngrid,"zdtlw","LW heating","T s-1",3,zdtlw)
        call writediagfi(ngrid,"dtrad","radiative heating","K s-1",3,dtrad)
        call writediagfi(ngrid,"zdtdyn","Dyn. heating","T s-1",3,zdtdyn)
        
        ! For Debugging.
        call writediagfi(ngrid,'pphi','Geopotential',' ',3,pphi)
        

      ! Output tracers.
      if (tracer) then

         if (callmufi) then
           ! Microphysical tracers are expressed in unit/m3.
           ! convert X.kg-1 --> X.m-3 (whereas for optics was -> X.m-2)
           i2e(:,:) = ( pplev(:,1:nlayer)-pplev(:,2:nlayer+1) ) / gzlat(:,1:nlayer) /(zzlev(:,2:nlayer+1)-zzlev(:,1:nlayer))

#ifdef USE_QTEST
            ! Microphysical tracers passed through dyn+phys(except mufi)
            call writediagfi(ngrid,"mu_m0as_dp","Dynphys only spherical mode 0th order moment",'m-3',3,zq(:,:,micro_indx(1))*i2e)
            call writediagfi(ngrid,"mu_m3as_dp","Dynphys only spherical mode 3rd order moment",'m3/m3',3,zq(:,:,micro_indx(2))*i2e)
            call writediagfi(ngrid,"mu_m0af_dp","Dynphys only fractal mode 0th order moment",'m-3',3,zq(:,:,micro_indx(3))*i2e)
            call writediagfi(ngrid,"mu_m3af_dp","Dynphys only fractal mode 3rd order moment",'m3/m3',3,zq(:,:,micro_indx(4))*i2e)
            ! Microphysical tracers passed through mufi only
            call writediagfi(ngrid,"mu_m0as_mo","Mufi only spherical mode 0th order moment",'m-3',3,tpq(:,:,micro_indx(1))*i2e)
            call writediagfi(ngrid,"mu_m3as_mo","Mufi only spherical mode 3rd order moment",'m3/m3',3,tpq(:,:,micro_indx(2))*i2e)
            call writediagfi(ngrid,"mu_m0af_mo","Mufi only fractal mode 0th order moment",'m-3',3,tpq(:,:,micro_indx(3))*i2e)
            call writediagfi(ngrid,"mu_m3af_mo","Mufi only fractal mode 3rd order moment",'m3/m3',3,tpq(:,:,micro_indx(4))*i2e)
#else
            call writediagfi(ngrid,"mu_m0as","Spherical mode 0th order moment",'m-3',3,zq(:,:,micro_indx(1))*i2e)
            call writediagfi(ngrid,"mu_m3as","Spherical mode 3rd order moment",'m3/m3',3,zq(:,:,micro_indx(2))*i2e)
            call writediagfi(ngrid,"mu_m0af","Fractal mode 0th order moment",'m-3',3,zq(:,:,micro_indx(3))*i2e)
            call writediagfi(ngrid,"mu_m3af","Fractal mode 3rd order moment",'m3/m3',3,zq(:,:,micro_indx(4))*i2e)
#endif
            
            ! Microphysical diagnostics
            call writediagfi(ngrid,"mmd_aer_prec","Total aerosols precipitations",'m',2,mmd_aer_prec)
            call writediagfi(ngrid,"mmd_aer_s_flux","Spherical aerosols sedimentation flux",'kg.m-2.s-1',3,mmd_aer_s_flux)
            call writediagfi(ngrid,"mmd_aer_f_flux","Fractal aerosols sedimentation flux",'kg.m-2.s-1',3,mmd_aer_f_flux)
            call writediagfi(ngrid,"mmd_rc_sph","Spherical mode caracteristic radius",'m',3,mmd_rc_sph)
            call writediagfi(ngrid,"mmd_rc_fra","Fractal mode caracteristic radius",'m',3,mmd_rc_fra)

         endif ! end of 'callmufi'

         ! Chemical tracers
         if (callchim) then
           do iq=1,nkim
             call writediagfi(ngrid,cnames(iq),cnames(iq),'mol/mol',3,zq(:,:,iq+nmicro)/rat_mmol(iq+nmicro))
           enddo
           call writediagfi(ngrid,"evapCH4","Surface CH4 pseudo-evaporation rate",'mol/mol/s',2,dycevapCH4)
         endif

       endif ! end of 'tracer'

! XIOS outputs
#ifdef CPP_XIOS      
      ! Send fields to XIOS: (NB these fields must also be defined as
      ! <field id="..." /> in context_lmdz_physics.xml to be correctly used)
      CALL send_xios_field("ls",zls*180./pi)
      CALL send_xios_field("lss",zlss*180./pi)
      CALL send_xios_field("RA",right_ascen*180./pi)
      CALL send_xios_field("Declin",declin*180./pi)
      
      ! Total energy balance diagnostics
      if (callrad.and.(.not.newtonian)) then
         CALL send_xios_field("ISR_TOA",fluxtop_dn)
         CALL send_xios_field("OLR_TOA",fluxtop_lw)
      endif
      
      CALL send_xios_field("area",cell_area)
      CALL send_xios_field("pphi",pphi)
      CALL send_xios_field("pphis",phisfi)
      
      CALL send_xios_field("ps",ps)
      CALL send_xios_field("tsurf",tsurf)

      if(enertest) then
         if (calldifv) then
            CALL send_xios_field("sensibFlux",sensibFlux)
         endif
      endif

      CALL send_xios_field("temp",zt)
      CALL send_xios_field("teta",zh)
      CALL send_xios_field("u",zu)
      CALL send_xios_field("v",zv)
      CALL send_xios_field("w",pw)
      CALL send_xios_field("p",pplay)
      
      ! Winds diagnostics.
      CALL send_xios_field("dudif",zdudif)
      CALL send_xios_field("dudyn",zdudyn)

      CALL send_xios_field("horizwind",zhorizwind)

      ! Heating diagnostics.
      CALL send_xios_field("dtsw",zdtsw)
      CALL send_xios_field("dtlw",zdtlw)
      CALL send_xios_field("dtrad",dtrad)
      CALL send_xios_field("dtdyn",zdtdyn)
      CALL send_xios_field("dtdif",zdtdif)

      ! Chemical tracers
      if (callchim) then
      
        ! Advected fields
        do iq=1,nkim
          CALL send_xios_field(trim(cnames(iq)),zq(:,:,iq+nmicro)/rat_mmol(iq+nmicro)) ! kg/kg -> mol/mol
        enddo
        
        ! Upper chemistry fields
        do iq=1,nkim
          CALL send_xios_field(trim(cnames(iq))//"_up",ykim_up(iq,:,:)) ! mol/mol
        enddo
        
        ! Append fields in ykim_tot for output on the total vertical grid (0->1300km)
        do iq=1,nkim
          
          ! GCM levels
          do l=1,nlayer
            ykim_tot(iq,:,l) = zq(:,l,iq+nmicro)/rat_mmol(iq+nmicro)
          enddo
          ! Upper levels
          do l=1,nlaykim_up
            ykim_tot(iq,:,nlayer+l) = ykim_up(iq,:,l)
          enddo
          
          CALL send_xios_field(trim(cnames(iq))//"_tot",ykim_tot(iq,:,:)) ! mol/mol
          
        enddo
       
        ! Condensation tendencies ( mol/mol/s )
        CALL send_xios_field("dqcond_CH4",dyccond(:,:,7+nmicro))
        CALL send_xios_field("dqcond_C2H2",dyccond(:,:,10+nmicro))
        CALL send_xios_field("dqcond_C2H4",dyccond(:,:,12+nmicro))
        CALL send_xios_field("dqcond_C2H6",dyccond(:,:,14+nmicro))
        CALL send_xios_field("dqcond_C3H6",dyccond(:,:,17+nmicro))
        CALL send_xios_field("dqcond_C4H4",dyccond(:,:,21+nmicro))
        CALL send_xios_field("dqcond_CH3CCH",dyccond(:,:,23+nmicro))
        CALL send_xios_field("dqcond_C3H8",dyccond(:,:,24+nmicro))
        CALL send_xios_field("dqcond_C4H2",dyccond(:,:,25+nmicro))
        CALL send_xios_field("dqcond_C4H6",dyccond(:,:,26+nmicro))
        CALL send_xios_field("dqcond_C4H10",dyccond(:,:,27+nmicro))
        CALL send_xios_field("dqcond_AC6H6",dyccond(:,:,28+nmicro))
        CALL send_xios_field("dqcond_HCN",dyccond(:,:,35+nmicro))
        CALL send_xios_field("dqcond_CH3CN",dyccond(:,:,39+nmicro))
        CALL send_xios_field("dqcond_HC3N",dyccond(:,:,41+nmicro))
        CALL send_xios_field("dqcond_NCCN",dyccond(:,:,42+nmicro))
        CALL send_xios_field("dqcond_C4N2",dyccond(:,:,43+nmicro))

        ! Pseudo-evaporation flux (mol/mol/s)
        CALL send_xios_field("evapCH4",dycevapCH4(:))

      endif ! of 'if callchim'

      ! Microphysical tracers
      if (callmufi) then
        CALL send_xios_field("mu_m0as",zq(:,:,micro_indx(1))*i2e)
        CALL send_xios_field("mu_m3as",zq(:,:,micro_indx(2))*i2e)
        CALL send_xios_field("mu_m0af",zq(:,:,micro_indx(3))*i2e)
        CALL send_xios_field("mu_m3af",zq(:,:,micro_indx(4))*i2e)
      endif        

      if (lastcall.and.is_omp_master) then
        write(*,*) "physiq: call xios_context_finalize"
        call xios_context_finalize
      endif
#endif

      icount=icount+1

    end subroutine physiq
    
    end module physiq_mod
