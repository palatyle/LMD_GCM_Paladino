      MODULE physiq_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE physiq(
     $            ngrid,nlayer,nq
     $            ,firstcall,lastcall
     $            ,pday,ptime,ptimestep
     $            ,pplev,pplay,pphi
     $            ,pu,pv,pt,pq
     $            ,flxw
     $            ,pdu,pdv,pdt,pdq,pdpsrf)

      use watercloud_mod, only: watercloud, zdqcloud, zdqscloud
      use calchim_mod, only: calchim, ichemistry, zdqchim, zdqschim
      use watersat_mod, only: watersat
      use co2condens_mod, only: co2condens
      use co2cloud_mod, only: co2cloud, mem_Mccn_co2, mem_Mh2o_co2,
     &                        mem_Nccn_co2
      use callradite_mod, only: callradite
      use callsedim_mod, only: callsedim
      use rocketduststorm_mod, only: rocketduststorm, dustliftday
      use calcstormfract_mod, only: calcstormfract
      use topmons_mod, only: topmons,alpha_hmons
      use tracer_mod, only: noms, mmol, igcm_co2, igcm_n2, igcm_co2_ice,
     &                      igcm_co, igcm_o, igcm_h2o_vap, igcm_h2o_ice,
     &                      igcm_ccn_mass, igcm_ccn_number,
     &                      igcm_ccnco2_mass, igcm_ccnco2_number,
     &                      rho_ice_co2,nuiceco2_sed,nuiceco2_ref,
     &                      igcm_dust_mass, igcm_dust_number, igcm_h2o2,
     &                      nuice_ref, rho_ice, rho_dust, ref_r0,
     &                      igcm_he, igcm_stormdust_mass,
     &                      igcm_stormdust_number, igcm_topdust_mass,
     &                      igcm_topdust_number
      use comsoil_h, only: inertiedat, ! soil thermal inertia
     &                     tsoil, nsoilmx,!number of subsurface layers
     &                     mlayer,layer ! soil mid layer depths
      use geometry_mod, only: longitude, latitude, cell_area,
     &                        longitude_deg 
      use comgeomfi_h, only: sinlon, coslon, sinlat, coslat
      use surfdat_h, only: phisfi, albedodat, zmea, zstd, zsig, zgam,
     &                     zthe, z0, albedo_h2o_ice,
     &                     frost_albedo_threshold,
     &                     tsurf, co2ice, emis,
     &                     capcal, fluxgrd, qsurf,
     &                     hmons,summit,base,watercap
      use comsaison_h, only: dist_sol, declin, mu0, fract, local_time
      use slope_mod, only: theta_sl, psi_sl
      use conc_mod, only: rnew, cpnew, mmean
      use time_phylmdz_mod, only: iphysiq, day_step, ecritstart, daysec
      use dimradmars_mod, only: tauscaling, aerosol, totcloudfrac,
     &                          dtrad, fluxrad_sky, fluxrad, albedo,
     &                          naerkind, iaer_dust_doubleq, 
     &                          iaer_stormdust_doubleq
      use turb_mod, only: q2, wstar, ustar, sensibFlux, 
     &                    zmax_th, hfmax_th, turb_resolved
      use planete_h, only: aphelie, periheli, year_day, peri_day,
     &                     obliquit
      USE comcstfi_h, only: r, cpp, mugaz, g, rcp, pi, rad 
      USE calldrag_noro_mod, ONLY: calldrag_noro
      USE vdifc_mod, ONLY: vdifc
      use param_v4_h, only: nreact,n_avog, 
     &                      fill_data_thermos, allocate_param_thermos
      use iono_h, only: allocate_param_iono
      use compute_dtau_mod, only: compute_dtau
      use nonoro_gwd_ran_mod, only: nonoro_gwd_ran
#ifdef MESOSCALE
      use comsoil_h, only: mlayer,layer
      use surfdat_h, only: z0_default
      use comm_wrf
#else
      USE planetwide_mod, ONLY: planetwide_maxval, planetwide_minval,
     &                          planetwide_sumval
      use phyredem, only: physdem0, physdem1
      use phyetat0_mod, only: phyetat0
      use eofdump_mod, only: eofdump
      USE vertical_layers_mod, ONLY: ap,bp,aps,bps,presnivs,pseudoalt
      USE mod_phys_lmdz_omp_data, ONLY: is_omp_master
#endif

#ifdef CPP_XIOS      
      use xios_output_mod, only: initialize_xios_output,
     &                           update_xios_timestep,
     &                           send_xios_field
      use wxios, only: wxios_context_init, xios_context_finalize
#endif

      IMPLICIT NONE
c=======================================================================
c
c   subject:
c   --------
c
c   Organisation of the physical parametrisations of the LMD 
c   martian atmospheric general circulation model.
c
c   The GCM can be run without or with tracer transport
c   depending on the value of Logical "tracer" in file  "callphys.def"
c   Tracers may be water vapor, ice OR chemical species OR dust particles
c
c   SEE comments in initracer.F about numbering of tracer species...
c
c   It includes:
c
c      1. Initialization:
c      1.1 First call initializations
c      1.2 Initialization for every call to physiq
c      1.2.5 Compute mean mass and cp, R and thermal conduction coeff.
c      2. Compute radiative transfer tendencies
c         (longwave and shortwave) for CO2 and aerosols.
c      3. Gravity wave and subgrid scale topography drag :
c      4. Vertical diffusion (turbulent mixing):
c      5. Convective adjustment
c      6. Condensation and sublimation of carbon dioxide.
c      7.  TRACERS :
c       7a. water, water ice, co2 ice (clouds)
c       7b. call for photochemistry when tracers are chemical species
c       7c. other scheme for tracer (dust) transport (lifting, sedimentation)
c       7d. updates (CO2 pressure variations, surface budget)
c      8. Contribution to tendencies due to thermosphere
c      9. Surface and sub-surface temperature calculations
c     10. Write outputs :
c           - "startfi", "histfi" (if it's time)
c           - Saving statistics (if "callstats = .true.")
c           - Dumping eof (if "calleofdump = .true.")
c           - Output any needed variables in "diagfi" 
c     11. Diagnostic: mass conservation of tracers
c 
c   author: 
c   ------- 
c           Frederic Hourdin	15/10/93
c           Francois Forget		1994
c           Christophe Hourdin	02/1997 
c           Subroutine completly rewritten by F.Forget (01/2000)
c           Introduction of the photochemical module: S. Lebonnois (11/2002)
c           Introduction of the thermosphere module: M. Angelats i Coll (2002)
c           Water ice clouds: Franck Montmessin (update 06/2003)
c           Radiatively active tracers: J.-B. Madeleine (10/2008-06/2009)
c             Nb: See callradite.F for more information.
c           Mesoscale lines: Aymeric Spiga (2007 - 2011) -- check MESOSCALE flags
c           jul 2011 malv+fgg: Modified calls to NIR heating routine and 15 um cooling parameterization
c
c           10/16 J. Audouard: modifications for CO2 clouds scheme

c   arguments:
c   ----------
c
c   input:
c   ------
c    ecri                  period (in dynamical timestep) to write output
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    nq                    Number of advected fields
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c    ptime                 Universal time (0<ptime<1): ptime=0.5 at 12:00 UT
c    ptimestep             timestep (s)
c    pplay(ngrid,nlayer)   Pressure at the middle of the layers (Pa)
c    pplev(ngrid,nlayer+1) intermediate pressure levels (pa)
c    pphi(ngrid,nlayer)    Geopotential at the middle of the layers (m2s-2)
c    pu(ngrid,nlayer)      u component of the wind (ms-1)
c    pv(ngrid,nlayer)      v component of the wind (ms-1)
c    pt(ngrid,nlayer)      Temperature (K)
c    pq(ngrid,nlayer,nq)   Advected fields
c    pudyn(ngrid,nlayer)    | 
c    pvdyn(ngrid,nlayer)    | Dynamical temporal derivative for the
c    ptdyn(ngrid,nlayer)    | corresponding variables
c    pqdyn(ngrid,nlayer,nq) |
c    flxw(ngrid,nlayer)      vertical mass flux (kg/s) at layer lower boundary
c
c   output:
c   -------
c
c    pdu(ngrid,nlayer)       |
c    pdv(ngrid,nlayer)       |  Temporal derivative of the corresponding
c    pdt(ngrid,nlayer)       |  variables due to physical processes.
c    pdq(ngrid,nlayer,nq)    |
c    pdpsrf(ngrid)           |

c
c=======================================================================
c
c    0.  Declarations :
c    ------------------

      include "callkeys.h"
      include "comg1d.h"
      include "nlteparams.h"
      include "netcdf.inc"

c Arguments :
c -----------

c   inputs:
c   -------
      INTEGER,INTENT(in) :: ngrid ! number of atmospheric columns
      INTEGER,INTENT(in) :: nlayer ! number of atmospheric layers
      INTEGER,INTENT(in) :: nq ! number of tracers
      LOGICAL,INTENT(in) :: firstcall ! signals first call to physics
      LOGICAL,INTENT(in) :: lastcall ! signals last call to physics
      REAL,INTENT(in) :: pday ! number of elapsed sols since reference Ls=0
      REAL,INTENT(in) :: ptime ! "universal time", given as fraction of sol (e.g.: 0.5 for noon)
      REAL,INTENT(in) :: ptimestep ! physics timestep (s)
      REAL,INTENT(in) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pplay(ngrid,nlayer) ! mid-layer pressure (Pa)
      REAL,INTENT(IN) :: pphi(ngrid,nlayer) ! geopotential at mid-layer (m2s-2)
      REAL,INTENT(in) :: pu(ngrid,nlayer) ! zonal wind component (m/s)
      REAL,INTENT(in) :: pv(ngrid,nlayer) ! meridional wind component (m/s)
      REAL,INTENT(in) :: pt(ngrid,nlayer) ! temperature (K)
      REAL,INTENT(in) :: pq(ngrid,nlayer,nq) ! tracers (.../kg_of_air)
      REAL,INTENT(in) :: flxw(ngrid,nlayer) ! vertical mass flux (ks/s)
                                            ! at lower boundary of layer

c   outputs:
c   --------
c     physical tendencies
      REAL,INTENT(out) :: pdu(ngrid,nlayer) ! zonal wind tendency (m/s/s)
      REAL,INTENT(out) :: pdv(ngrid,nlayer) ! meridional wind tendency (m/s/s)
      REAL,INTENT(out) :: pdt(ngrid,nlayer) ! temperature tendency (K/s)
      REAL,INTENT(out) :: pdq(ngrid,nlayer,nq) ! tracer tendencies (../kg/s)
      REAL,INTENT(out) :: pdpsrf(ngrid) ! surface pressure tendency (Pa/s)

      

c Local saved variables:
c ----------------------
      INTEGER,SAVE :: day_ini ! Initial date of the run (sol since Ls=0)
      INTEGER,SAVE :: icount     ! counter of calls to physiq during the run.

#ifdef DUSTSTORM
      REAL pq_tmp(ngrid, nlayer, 2) ! To compute tendencies due the dust bomb
#endif
            
c     Variables used by the water ice microphysical scheme:
      REAL rice(ngrid,nlayer)    ! Water ice geometric mean radius (m)
      REAL nuice(ngrid,nlayer)   ! Estimated effective variance
                                     !   of the size distribution
      real rsedcloud(ngrid,nlayer) ! Cloud sedimentation radius
      real rhocloud(ngrid,nlayer)  ! Cloud density (kg.m-3)
      REAL inertiesoil(ngrid,nsoilmx) ! Time varying subsurface
                                      ! thermal inertia (J.s-1/2.m-2.K-1)
                                      ! (used only when tifeedback=.true.)
c     Variables used by the CO2 clouds microphysical scheme:
      DOUBLE PRECISION riceco2(ngrid,nlayer)   ! co2 ice geometric mean radius (m)
      real rsedcloudco2(ngrid,nlayer) !CO2 Cloud sedimentation radius
      real rhocloudco2(ngrid,nlayer) !co2 Cloud density (kg.m-3)
      real zdqssed_co2(ngrid)  ! CO2 flux at the surface  (kg.m-2.s-1)
      real zcondicea_co2microp(ngrid,nlayer)
c     Variables used by the photochemistry
      REAL surfdust(ngrid,nlayer) ! dust surface area (m2/m3, if photochemistry)
      REAL surfice(ngrid,nlayer)  !  ice surface area (m2/m3, if photochemistry)
c     Variables used by the slope model
      REAL sl_ls, sl_lct, sl_lat
      REAL sl_tau, sl_alb, sl_the, sl_psi
      REAL sl_fl0, sl_flu
      REAL sl_ra, sl_di0
      REAL sky

      REAL,PARAMETER :: stephan = 5.67e-08 ! Stephan Boltzman constant

c Local variables :
c -----------------

      REAL CBRT
      EXTERNAL CBRT

!      CHARACTER*80 fichier 
      INTEGER l,ig,ierr,igout,iq,tapphys,isoil

      REAL fluxsurf_lw(ngrid)      !incident LW (IR) surface flux (W.m-2)
      REAL fluxsurf_sw(ngrid,2)    !incident SW (solar) surface flux (W.m-2)
      REAL fluxtop_lw(ngrid)       !Outgoing LW (IR) flux to space (W.m-2)
      REAL fluxtop_sw(ngrid,2)     !Outgoing SW (solar) flux to space (W.m-2)
      REAL tauref(ngrid)           ! Reference column optical depth at odpref
c     rocket dust storm
      REAL totstormfract(ngrid)     ! fraction of the mesh where the dust storm is contained
      logical clearatm              ! clearatm used to calculate twice the radiative 
                                    ! transfer when rdstorm is active : 
                                    !            - in a mesh with stormdust and background dust (false)
                                    !            - in a mesh with background dust only (true)
c     entrainment by slope winds
      logical nohmons               ! nohmons used to calculate twice the radiative 
                                    ! transfer when slpwind is active : 
                                    !            - in a mesh with topdust and background dust (false)
                                    !            - in a mesh with background dust only (true)
      
      real,parameter :: odpref=610. ! DOD reference pressure (Pa)
      REAL tau(ngrid,naerkind)     ! Column dust optical depth at each point
                                   ! AS: TBD: this one should be in a module !
      REAL zls                       !  solar longitude (rad)
      REAL zday                      ! date (time since Ls=0, in martian days)
      REAL zzlay(ngrid,nlayer)     ! altitude at the middle of the layers
      REAL zzlev(ngrid,nlayer+1)   ! altitude at layer boundaries
!      REAL latvl1,lonvl1             ! Viking Lander 1 point (for diagnostic)

c     Tendancies due to various processes:
      REAL dqsurf(ngrid,nq)
      REAL zdtlw(ngrid,nlayer)     ! (K/s)
      REAL zdtsw(ngrid,nlayer)     ! (K/s)
      REAL pdqrds(ngrid,nlayer,nq) ! tendency for dust after rocketduststorm

      REAL zdtnirco2(ngrid,nlayer) ! (K/s)
      REAL zdtnlte(ngrid,nlayer)   ! (K/s)
      REAL zdtsurf(ngrid)            ! (K/s)
      REAL zdtcloud(ngrid,nlayer),zdtcloudco2(ngrid,nlayer)
      REAL zdvdif(ngrid,nlayer),zdudif(ngrid,nlayer)  ! (m.s-2)
      REAL zdhdif(ngrid,nlayer), zdtsdif(ngrid)         ! (K/s)
      REAL zdvadj(ngrid,nlayer),zduadj(ngrid,nlayer)  ! (m.s-2)
      REAL zdhadj(ngrid,nlayer)                           ! (K/s)
      REAL zdtgw(ngrid,nlayer)                            ! (K/s)
      REAL zdugw(ngrid,nlayer),zdvgw(ngrid,nlayer)    ! (m.s-2)
      REAL zdtc(ngrid,nlayer),zdtsurfc(ngrid)
      REAL zdvc(ngrid,nlayer),zduc(ngrid,nlayer)

      REAL zdqdif(ngrid,nlayer,nq), zdqsdif(ngrid,nq)
      REAL zdqsed(ngrid,nlayer,nq), zdqssed(ngrid,nq)
      REAL zdqdev(ngrid,nlayer,nq), zdqsdev(ngrid,nq)
      REAL zdqadj(ngrid,nlayer,nq)
      REAL zdqc(ngrid,nlayer,nq)
      REAL zdqcloudco2(ngrid,nlayer,nq)
      REAL zdqsc(ngrid,nq)

      REAL zdteuv(ngrid,nlayer)    ! (K/s)
      REAL zdtconduc(ngrid,nlayer) ! (K/s)
      REAL zdumolvis(ngrid,nlayer)
      REAL zdvmolvis(ngrid,nlayer)
      real zdqmoldiff(ngrid,nlayer,nq)

      REAL dwatercap(ngrid), dwatercap_dif(ngrid)     ! (kg/m-2)

c     Local variable for local intermediate calcul:
      REAL zflubid(ngrid)
      REAL zplanck(ngrid),zpopsk(ngrid,nlayer)
      REAL zdum1(ngrid,nlayer)
      REAL zdum2(ngrid,nlayer)
      REAL ztim1,ztim2,ztim3, z1,z2
      REAL ztime_fin
      REAL zdh(ngrid,nlayer)
      REAL zh(ngrid,nlayer)      ! potential temperature (K)
      REAL pw(ngrid,nlayer) ! vertical velocity (m/s) (>0 when downwards)
      INTEGER length
      PARAMETER (length=100)

c     Variables for the total dust for diagnostics
      REAL qdusttotal(ngrid,nlayer) !it equals to dust + stormdust

      INTEGER iaer

c local variables only used for diagnostic (output in file "diagfi" or "stats")
c -----------------------------------------------------------------------------
      REAL ps(ngrid), zt(ngrid,nlayer)
      REAL zu(ngrid,nlayer),zv(ngrid,nlayer)
      REAL zq(ngrid,nlayer,nq)

      REAL fluxtop_sw_tot(ngrid), fluxsurf_sw_tot(ngrid)
      character*2 str2
!      character*5 str5
      real zdtdif(ngrid,nlayer), zdtadj(ngrid,nlayer)
      real rdust(ngrid,nlayer) ! dust geometric mean radius (m)
      real rstormdust(ngrid,nlayer) ! stormdust geometric mean radius (m)
      real rtopdust(ngrid,nlayer)   ! topdust geometric mean radius (m)
      integer igmin, lmin
      logical tdiag

      real co2col(ngrid)        ! CO2 column
      ! pplev and pplay are dynamical inputs and must not be modified in the physics.
      ! instead, use zplay and zplev :
      REAL zplev(ngrid,nlayer+1),zplay(ngrid,nlayer) 
!      REAL zstress(ngrid),cd
      real tmean, zlocal(nlayer)
      real rho(ngrid,nlayer)  ! density
      real vmr(ngrid,nlayer)  ! volume mixing ratio
      real rhopart(ngrid,nlayer) ! number density of a given species
      real colden(ngrid,nq)     ! vertical column of tracers
      real mass(nq)             ! global mass of tracers (g)
      REAL mtot(ngrid)          ! Total mass of water vapor (kg/m2)
      REAL mstormdtot(ngrid)    ! Total mass of stormdust tracer (kg/m2)
      REAL mdusttot(ngrid)    ! Total mass of dust tracer (kg/m2)
      REAL icetot(ngrid)        ! Total mass of water ice (kg/m2)
      REAL mtotco2(ngrid)       ! Total mass of co2 vapor (kg/m2)
      REAL icetotco2(ngrid)     ! Total mass of co2 ice (kg/m2) 
      REAL Nccntot(ngrid)       ! Total number of ccn (nbr/m2)
      REAL Mccntot(ngrid)       ! Total mass of ccn (kg/m2)
      REAL rave(ngrid)          ! Mean water ice effective radius (m)
      REAL opTES(ngrid,nlayer)  ! abs optical depth at 825 cm-1
      REAL tauTES(ngrid)        ! column optical depth at 825 cm-1
      REAL Qabsice                ! Water ice absorption coefficient
      REAL taucloudtes(ngrid)  ! Cloud opacity at infrared
                               !   reference wavelength using
                               !   Qabs instead of Qext
                               !   (direct comparison with TES)
                               
      REAL dqdustsurf(ngrid) ! surface q dust flux (kg/m2/s)
      REAL dndustsurf(ngrid) ! surface n dust flux (number/m2/s)
      REAL ndust(ngrid,nlayer) ! true n dust (kg/kg)
      REAL qdust(ngrid,nlayer) ! true q dust (kg/kg)
      REAL nccn(ngrid,nlayer)  ! true n ccn (kg/kg)
      REAL qccn(ngrid,nlayer)  ! true q ccn (kg/kg)
c     definition tendancies of stormdust tracers
      REAL rdsdqdustsurf(ngrid) ! surface q stormdust flux (kg/m2/s)
      REAL rdsdndustsurf(ngrid) ! surface n stormdust flux (number/m2/s)
      REAL rdsndust(ngrid,nlayer) ! true n stormdust (kg/kg)
      REAL rdsqdust(ngrid,nlayer) ! true q stormdust (kg/kg)
      REAL wspeed(ngrid,nlayer+1) ! vertical velocity stormdust tracer 
      REAL dsodust(ngrid,nlayer) 
      REAL dsords(ngrid,nlayer)
      REAL dsotop(ngrid,nlayer)
      REAL wtop(ngrid,nlayer+1) ! vertical velocity topdust tracer

      REAL nccnco2(ngrid,nlayer)   ! true n ccnco2 (kg/kg)
      REAL qccnco2(ngrid,nlayer)  ! true q ccnco2 (kg/kg)

c Test 1d/3d scavenging
      real h2otot(ngrid)
      REAL satu(ngrid,nlayer)  ! satu ratio for output
      REAL zqsat(ngrid,nlayer) ! saturation
      REAL satuco2(ngrid,nlayer)  ! co2 satu ratio for output
      REAL zqsatco2(ngrid,nlayer) ! saturation co2
      REAL,SAVE :: time_phys

! Added for new NLTE scheme

      real co2vmr_gcm(ngrid,nlayer)
      real n2vmr_gcm(ngrid,nlayer)
      real ovmr_gcm(ngrid,nlayer)
      real covmr_gcm(ngrid,nlayer)
      integer ierr_nlte
      real*8  varerr

C  Non-oro GW drag & Calcul of Brunt-Vaisala freq. (BV2) 
      REAL ztetalev(ngrid,nlayer)
      real zdtetalev(ngrid,nlayer), zdzlev(ngrid,nlayer)
      REAL bv2(ngrid,nlayer)    ! BV2 at zlev   
c  Non-oro GW tendencies
      REAL d_u_hin(ngrid,nlayer), d_v_hin(ngrid,nlayer)
      REAL d_t_hin(ngrid,nlayer)
c  Diagnostics 2D of gw_nonoro
      REAL zustrhi(ngrid), zvstrhi(ngrid)
c Variables for PBL
      REAL zz1(ngrid)
      REAL lmax_th_out(ngrid)
      REAL pdu_th(ngrid,nlayer),pdv_th(ngrid,nlayer)
      REAL pdt_th(ngrid,nlayer),pdq_th(ngrid,nlayer,nq)
      INTEGER lmax_th(ngrid),dimout,n_out,n
      CHARACTER(50) zstring
      REAL dtke_th(ngrid,nlayer+1)
      REAL zcdv(ngrid), zcdh(ngrid)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: T_out
      REAL, ALLOCATABLE, DIMENSION(:,:) :: u_out ! Interpolated teta and u at z_out
      REAL u_out1(ngrid)
      REAL T_out1(ngrid)
      REAL, ALLOCATABLE, DIMENSION(:) :: z_out     ! height of interpolation between z0 and z1 [meters]
      REAL tstar(ngrid)  ! friction velocity and friction potential temp
      REAL L_mo(ngrid),vhf(ngrid),vvv(ngrid)
      real qdustrds0(ngrid,nlayer),qdustrds1(ngrid,nlayer)
      real qstormrds0(ngrid,nlayer),qstormrds1(ngrid,nlayer)   
      real qdusttotal0(ngrid),qdusttotal1(ngrid)

c sub-grid scale water ice clouds (A. Pottier 2013)
      logical clearsky
      ! flux for the part without clouds
      real zdtswclf(ngrid,nlayer)
      real zdtlwclf(ngrid,nlayer)
      real fluxsurf_lwclf(ngrid)     
      real fluxsurf_swclf(ngrid,2)  
      real fluxtop_lwclf(ngrid)
      real fluxtop_swclf(ngrid,2) 
      real taucloudtesclf(ngrid)
      real tf_clf, ntf_clf ! tf: fraction of clouds, ntf: fraction without clouds
      real rave2(ngrid), totrave2(ngrid) ! Mean water ice mean radius (m)

c entrainment by slope winds above sb-grid scale topography
      REAL pdqtop(ngrid,nlayer,nq) ! tendency for dust after topmons
      REAL hmax,hmin
      REAL hsummit(ngrid)

c when no startfi file is asked for init
      real alpha,lay1 ! coefficients for building layers
      integer iloop

!      LOGICAL startphy_file

c=======================================================================

c 1. Initialisation:
c -----------------

c  1.1   Initialisation only at first call
c  ---------------------------------------
      IF (firstcall) THEN

c        variables set to 0
c        ~~~~~~~~~~~~~~~~~~
         aerosol(:,:,:)=0
         dtrad(:,:)=0

#ifndef MESOSCALE
         fluxrad(:)=0
         wstar(:)=0.
#endif

#ifdef CPP_XIOS
         ! Initialize XIOS context
         write(*,*) "physiq: call wxios_context_init"
         CALL wxios_context_init
#endif

c        read startfi 
c        ~~~~~~~~~~~~
#ifndef MESOSCALE
! GCM. Read netcdf initial physical parameters.
         CALL phyetat0 ("startfi.nc",0,0,
     &         nsoilmx,ngrid,nlayer,nq,
     &         day_ini,time_phys,
     &         tsurf,tsoil,albedo,emis,
     &         q2,qsurf,co2ice,tauscaling,totcloudfrac,wstar,
     &         mem_Mccn_co2,mem_Nccn_co2,
     &         mem_Mh2o_co2,watercap)


#else
! MESOSCALE. Supposedly everything is already set in modules.
! So we just check. And we fill day_ini
      print*,"check: --- in physiq.F"
      print*,"check: rad,cpp,g,r,rcp,daysec"
      print*,rad,cpp,g,r,rcp,daysec
      PRINT*,'check: tsurf ',tsurf(1),tsurf(ngrid)
      PRINT*,'check: tsoil ',tsoil(1,1),tsoil(ngrid,nsoilmx)
      PRINT*,'check: inert ',inertiedat(1,1),inertiedat(ngrid,nsoilmx)
      PRINT*,'check: midlayer,layer ', mlayer(:),layer(:)
      PRINT*,'check: tracernames ', noms
      PRINT*,'check: emis ',emis(1),emis(ngrid)
      PRINT*,'check: q2 ',q2(1,1),q2(ngrid,nlayer+1)
      PRINT*,'check: qsurf ',qsurf(1,1),qsurf(ngrid,nq)
      PRINT*,'check: co2 ',co2ice(1),co2ice(ngrid)
      PRINT*,'check: watercap ',watercap(1),watercap(ngrid)
      !!!
      day_ini = pday
#endif
#ifndef MESOSCALE
         if (.not.startphy_file) then
           ! starting without startfi.nc and with callsoil
           ! is not yet possible as soildepth default is not defined
           if (callsoil) then                     
              ! default mlayer distribution, following a power law:
              !  mlayer(k)=lay1*alpha**(k-1/2)
              lay1=2.e-4
	      alpha=2
              do iloop=0,nsoilmx-1
	         mlayer(iloop)=lay1*(alpha**(iloop-0.5))
	      enddo
              lay1=sqrt(mlayer(0)*mlayer(1))
              alpha=mlayer(1)/mlayer(0)
              do iloop=1,nsoilmx
                 layer(iloop)=lay1*(alpha**(iloop-1))
              enddo
           endif
           ! additionnal "academic" initialization of physics
           write(*,*) "Physiq: initializing tsurf(:) to pt(:,1) !!"
           tsurf(:)=pt(:,1)
           write(*,*) "Physiq: initializing tsoil(:) to pt(:,1) !!"
           do isoil=1,nsoilmx
             tsoil(1:ngrid,isoil)=tsurf(1:ngrid)
           enddo
           write(*,*) "Physiq: initializing inertiedat !!"
           inertiedat(:,:)=400.
           write(*,*) "Physiq: initializing day_ini to pdat !"
           day_ini=pday
        endif
#endif
         if (pday.ne.day_ini) then
           write(*,*) "PHYSIQ: ERROR: bad synchronization between ",
     &                "physics and dynamics"
           write(*,*) "dynamics day [pday]: ",pday
           write(*,*) "physics day [day_ini]: ",day_ini
           call abort_physic("physiq","dynamics day /= physics day",1)
         endif

         write (*,*) 'In physiq day_ini =', day_ini

c        initialize tracers
c        ~~~~~~~~~~~~~~~~~~
         IF (tracer) THEN
            CALL initracer(ngrid,nq,qsurf)
         ENDIF  ! end tracer

c        Initialize albedo and orbital calculation
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL surfini(ngrid,co2ice,qsurf)
         CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)

c        initialize soil 
c        ~~~~~~~~~~~~~~~
         IF (callsoil) THEN
c           Thermal inertia feedback:
            IF (tifeedback) THEN
                CALL soil_tifeedback(ngrid,nsoilmx,qsurf,inertiesoil)
                CALL soil(ngrid,nsoilmx,firstcall,inertiesoil,
     s             ptimestep,tsurf,tsoil,capcal,fluxgrd)
            ELSE
                CALL soil(ngrid,nsoilmx,firstcall,inertiedat,
     s             ptimestep,tsurf,tsoil,capcal,fluxgrd)
            ENDIF ! of IF (tifeedback)
         ELSE
            PRINT*,
     &     'PHYSIQ WARNING! Thermal conduction in the soil turned off'
            DO ig=1,ngrid
               capcal(ig)=1.e5
               fluxgrd(ig)=0.
            ENDDO
         ENDIF
         icount=1

#ifndef MESOSCALE
c        Initialize thermospheric parameters
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (callthermos) then
            call fill_data_thermos
            call allocate_param_thermos(nlayer)
            call allocate_param_iono(nlayer,nreact)
            call param_read_e107
         endif
#endif
c        Initialize R and Cp as constant

         if (.not.callthermos .and. .not.photochem) then
                 do l=1,nlayer
                  do ig=1,ngrid
                   rnew(ig,l)=r
                   cpnew(ig,l)=cpp
                   mmean(ig,l)=mugaz
                   enddo
                  enddo  
         endif         

         if(callnlte.and.nltemodel.eq.2) call nlte_setup
         if(callnirco2.and.nircorr.eq.1) call NIR_leedat


        IF (tracer.AND.water.AND.(ngrid.NE.1)) THEN
          write(*,*)"physiq: water_param Surface water ice albedo:", 
     .                  albedo_h2o_ice
        ENDIF

#ifndef MESOSCALE

         if (ngrid.ne.1) then
	   ! no need to compute slopes when in 1D; it is an input
	   if (callslope) call getslopes(ngrid,phisfi)
	   ! no need to create a restart file in 1d
           call physdem0("restartfi.nc",longitude,latitude,
     &                   nsoilmx,ngrid,nlayer,nq,
     &                   ptimestep,pday,time_phys,cell_area,
     &                   albedodat,inertiedat,zmea,zstd,zsig,zgam,zthe,
     &                   hmons,summit,base)
         endif

c        Initialize mountain mesh fraction for the entrainment by slope wind param.
c        ~~~~~~~~~~~~~~~
        if (slpwind) then
          !! alpha_hmons calculation
          if (ngrid.gt.1) then
            call planetwide_maxval(hmons,hmax )
            call planetwide_minval(hmons,hmin )
            do ig=1,ngrid
              alpha_hmons(ig)= 0.5*(hmons(ig)-hmin)/(hmax-hmin)
            enddo
          else
            hmin=0.
            hmax=23162.1 !set here the height of the sub-grid scaled topography
            do ig=1,ngrid               
              alpha_hmons(ig)= (hmons(ig)-hmin)/(hmax-hmin) !0.1*(hmons(ig)-hmin)/(hmax-hmin)
              print*,"1D, hmons=",hmons(ig),"alpha=",alpha_hmons(ig)
            enddo
          endif ! (ngrid.gt.1)
        endif ! if (slpwind)

#endif
                 
#ifdef CPP_XIOS
        ! XIOS outputs
        write(*,*) "physiq firstcall: call initialize_xios_output"
        call initialize_xios_output(pday,ptime,ptimestep,daysec,
     &                              presnivs,pseudoalt)
#endif
      ENDIF        !  (end of "if firstcall")

c ---------------------------------------------------
c 1.2   Initializations done at every physical timestep:
c ---------------------------------------------------
c

#ifdef CPP_XIOS      
      ! update XIOS time/calendar
      call update_xios_timestep
#endif      
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      pdv(:,:)=0
      pdu(:,:)=0
      pdt(:,:)=0
      pdq(:,:,:)=0
      pdpsrf(:)=0
      zflubid(:)=0
      zdtsurf(:)=0
      dqsurf(:,:)=0
      dsodust(:,:)=0.
      dsords(:,:)=0.
      dsotop(:,:)=0.
      dwatercap(:)=0

#ifdef DUSTSTORM
      pq_tmp(:,:,:)=0
#endif
      igout=ngrid/2+1 


      zday=pday+ptime ! compute time, in sols (and fraction thereof)
      ! Compute local time at each grid point
      DO ig=1,ngrid
       local_time(ig)=modulo(1.+(zday-INT(zday))
     &                +(longitude_deg(ig)/15)/24,1.)
      ENDDO
c     Compute Solar Longitude (Ls) :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (season) then
         call solarlong(zday,zls)
      else
         call solarlong(float(day_ini),zls)
      end if

c     Initialize pressure levels
c     ~~~~~~~~~~~~~~~~~~
      zplev(:,:) = pplev(:,:)
      zplay(:,:) = pplay(:,:)
      ps(:) = pplev(:,1)

c     Compute geopotential at interlayers
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ponderation des altitudes au niveau des couches en dp/p

      DO l=1,nlayer
         DO ig=1,ngrid
            zzlay(ig,l)=pphi(ig,l)/g
         ENDDO
      ENDDO
      DO ig=1,ngrid
         zzlev(ig,1)=0.
         zzlev(ig,nlayer+1)=1.e7    ! dummy top of last layer above 10000 km...
      ENDDO
      DO l=2,nlayer
         DO ig=1,ngrid
            z1=(zplay(ig,l-1)+zplev(ig,l))/(zplay(ig,l-1)-zplev(ig,l))
            z2=(zplev(ig,l)+zplay(ig,l))/(zplev(ig,l)-zplay(ig,l))
            zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
         ENDDO
      ENDDO


!     Potential temperature calculation not the same in physiq and dynamic

c     Compute potential temperature
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO l=1,nlayer
         DO ig=1,ngrid 
            zpopsk(ig,l)=(zplay(ig,l)/zplev(ig,1))**rcp
            zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
         ENDDO
      ENDDO
#ifndef MESOSCALE
c-----------------------------------------------------------------------
c    1.2.5 Compute mean mass, cp, and R
c    --------------------------------

      if(photochem.or.callthermos) then
         call concentrations(ngrid,nlayer,nq,
     &                       zplay,pt,pdt,pq,pdq,ptimestep)
      endif
#endif

      ! Compute vertical velocity (m/s) from vertical mass flux
      ! w = F / (rho*area) and rho = P/(r*T)
      ! but first linearly interpolate mass flux to mid-layers
      do l=1,nlayer-1
       pw(1:ngrid,l)=0.5*(flxw(1:ngrid,l)+flxw(1:ngrid,l+1))
      enddo
      pw(1:ngrid,nlayer)=0.5*flxw(1:ngrid,nlayer) ! since flxw(nlayer+1)=0
      do l=1,nlayer
       pw(1:ngrid,l)=(pw(1:ngrid,l)*r*pt(1:ngrid,l)) / 
     &               (pplay(1:ngrid,l)*cell_area(1:ngrid))
       ! NB: here we use r and not rnew since this diagnostic comes
       ! from the dynamics
      enddo

c-----------------------------------------------------------------------
c    2. Compute radiative tendencies :
c------------------------------------

      IF (callrad) THEN

c       Local Solar zenith angle
c       ~~~~~~~~~~~~~~~~~~~~~~~~
        CALL orbite(zls,dist_sol,declin)

        IF (diurnal) THEN
            ztim1=SIN(declin)
            ztim2=COS(declin)*COS(2.*pi*(zday-.5))
            ztim3=-COS(declin)*SIN(2.*pi*(zday-.5))

            CALL solang(ngrid,sinlon,coslon,sinlat,coslat,
     &                  ztim1,ztim2,ztim3, mu0,fract)

        ELSE
            CALL mucorr(ngrid,declin,latitude,mu0,fract,10000.,rad)
        ENDIF ! of IF (diurnal)

         IF( MOD(icount-1,iradia).EQ.0) THEN

c          NLTE cooling from CO2 emission
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF(callnlte) then
              if(nltemodel.eq.0.or.nltemodel.eq.1) then
                 CALL nltecool(ngrid,nlayer,nq,zplay,pt,pq,zdtnlte)
              else if(nltemodel.eq.2) then
                co2vmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_co2)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_co2)
                n2vmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_n2)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_n2)
                covmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_co)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_co)
                ovmr_gcm(1:ngrid,1:nlayer)=
     &                      pq(1:ngrid,1:nlayer,igcm_o)*
     &                      mmean(1:ngrid,1:nlayer)/mmol(igcm_o)
                 
                 CALL nlte_tcool(ngrid,nlayer,zplay*9.869e-6,
     $                pt,zzlay,co2vmr_gcm, n2vmr_gcm, covmr_gcm, 
     $                ovmr_gcm,  zdtnlte,ierr_nlte,varerr )
                 if(ierr_nlte.gt.0) then
                    write(*,*)
     $                'WARNING: nlte_tcool output with error message',
     $                'ierr_nlte=',ierr_nlte,'varerr=',varerr
                    write(*,*)'I will continue anyway'
                 endif

                 zdtnlte(1:ngrid,1:nlayer)=
     &                             zdtnlte(1:ngrid,1:nlayer)/86400.
              endif
           ELSE
              zdtnlte(:,:)=0.
           ENDIF !end callnlte

c          Find number of layers for LTE radiation calculations
           IF(MOD(iphysiq*(icount-1),day_step).EQ.0)
     &          CALL nlthermeq(ngrid,nlayer,zplev,zplay)

c          rocketstorm : compute dust storm mesh fraction
           IF (rdstorm) THEN
              CALL calcstormfract(ngrid,nlayer,nq,pq,
     &                 totstormfract)
           ENDIF

c          Note: Dustopacity.F has been transferred to callradite.F

#ifdef DUSTSTORM
!! specific case: save the quantity of dust before adding perturbation
       if (firstcall) then
        pq_tmp(1:ngrid,1:nlayer,1)=pq(1:ngrid,1:nlayer,igcm_dust_mass)
        pq_tmp(1:ngrid,1:nlayer,2)=pq(1:ngrid,1:nlayer,igcm_dust_number)
       endif
#endif
         
c          Call main radiative transfer scheme
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          Transfer through CO2 (except NIR CO2 absorption)
c            and aerosols (dust and water ice)
           ! callradite for background dust
           clearatm=.true.
           !! callradite for background dust in the case of slope wind entrainment
           nohmons=.true.
c          Radiative transfer
c          ------------------
           ! callradite for the part with clouds
           clearsky=.false. ! part with clouds for both cases CLFvarying true/false 
           CALL callradite(icount,ngrid,nlayer,nq,zday,zls,pq,albedo,
     &     emis,mu0,zplev,zplay,pt,tsurf,fract,dist_sol,igout,
     &     zdtlw,zdtsw,fluxsurf_lw,fluxsurf_sw,fluxtop_lw,
     &     fluxtop_sw,tauref,tau,aerosol,dsodust,tauscaling,
     &     taucloudtes,rdust,rice,nuice,co2ice,rstormdust,rtopdust,
     &     totstormfract,clearatm,dsords,dsotop,alpha_hmons,nohmons,
     &     clearsky,totcloudfrac)

           ! case of sub-grid water ice clouds: callradite for the clear case
            IF (CLFvarying) THEN
               ! ---> PROBLEMS WITH ALLOCATED ARRAYS 
               ! (temporary solution in callcorrk: do not deallocate
               ! if
               ! CLFvarying ...) ?? AP ??
               clearsky=.true. ! 
               CALL callradite(icount,ngrid,nlayer,nq,zday,zls,pq,
     &              albedo,emis,mu0,zplev,zplay,pt,tsurf,fract,
     &              dist_sol,igout,zdtlwclf,zdtswclf,fluxsurf_lwclf,
     &              fluxsurf_swclf,fluxtop_lwclf,fluxtop_swclf,tauref,
     &              tau,aerosol,dsodust,tauscaling,taucloudtesclf,rdust,
     &              rice,nuice,co2ice,rstormdust,rtopdust,totstormfract,
     &              clearatm,dsords,dsotop,alpha_hmons,nohmons,
     &              clearsky,totcloudfrac)
               clearsky = .false.  ! just in case.
               ! Sum the fluxes and heating rates from cloudy/clear
               ! cases
               DO ig=1,ngrid
                  tf_clf=totcloudfrac(ig)
                  ntf_clf=1.-tf_clf
                  fluxsurf_lw(ig) = ntf_clf*fluxsurf_lwclf(ig) 
     &                                      + tf_clf*fluxsurf_lw(ig)
                  fluxsurf_sw(ig,1) = ntf_clf*fluxsurf_swclf(ig,1) 
     &                                      + tf_clf*fluxsurf_sw(ig,1)
                  fluxsurf_sw(ig,2) = ntf_clf*fluxsurf_swclf(ig,2) 
     &                                      + tf_clf*fluxsurf_sw(ig,2)
                  fluxtop_lw(ig)  = ntf_clf*fluxtop_lwclf(ig) 
     &                                      + tf_clf*fluxtop_lw(ig)
                  fluxtop_sw(ig,1)  = ntf_clf*fluxtop_swclf(ig,1)  
     &                                      + tf_clf*fluxtop_sw(ig,1)
                  fluxtop_sw(ig,2)  = ntf_clf*fluxtop_swclf(ig,2)  
     &                                      + tf_clf*fluxtop_sw(ig,2)
                  taucloudtes(ig) = ntf_clf*taucloudtesclf(ig) 
     &                                      + tf_clf*taucloudtes(ig)
                  zdtlw(ig,1:nlayer) = ntf_clf*zdtlwclf(ig,1:nlayer) 
     &                                      + tf_clf*zdtlw(ig,1:nlayer)
                  zdtsw(ig,1:nlayer) = ntf_clf*zdtswclf(ig,1:nlayer) 
     &                                      + tf_clf*zdtsw(ig,1:nlayer)
               ENDDO

            ENDIF ! (CLFvarying)
            
!============================================================================
           
#ifdef DUSTSTORM
!! specific case: compute the added quantity of dust for perturbation
       if (firstcall) then
         pdq(1:ngrid,1:nlayer,igcm_dust_mass)= 
     &      pdq(1:ngrid,1:nlayer,igcm_dust_mass)     
     &      - pq_tmp(1:ngrid,1:nlayer,1)
     &      + pq(1:ngrid,1:nlayer,igcm_dust_mass)
         pdq(1:ngrid,1:nlayer,igcm_dust_number)=
     &      pdq(1:ngrid,1:nlayer,igcm_dust_number)
     &      - pq_tmp(1:ngrid,1:nlayer,2)
     &      + pq(1:ngrid,1:nlayer,igcm_dust_number)
       endif
#endif

c          Outputs for basic check (middle of domain)
c          ------------------------------------------
           write(*,'("Ls =",f11.6," check lat =",f10.6,
     &               " lon =",f11.6)')
     &           zls*180./pi,latitude(igout)*180/pi,
     &                       longitude(igout)*180/pi

           write(*,'(" tauref(",f4.0," Pa) =",f9.6,
     &             " tau(",f4.0," Pa) =",f9.6)')
     &            odpref,tauref(igout),
     &            odpref,tau(igout,1)*odpref/zplev(igout,1)
c          ---------------------------------------------------------
c          Call slope parameterization for direct and scattered flux
c          ---------------------------------------------------------
           IF(callslope) THEN
            print *, 'Slope scheme is on and computing...'
            DO ig=1,ngrid  
              sl_the = theta_sl(ig)
              IF (sl_the .ne. 0.) THEN
                ztim1=fluxsurf_sw(ig,1)+fluxsurf_sw(ig,2)
                DO l=1,2
                 sl_lct = ptime*24. + 180.*longitude(ig)/pi/15.
                 sl_ra  = pi*(1.0-sl_lct/12.)
                 sl_lat = 180.*latitude(ig)/pi
                 sl_tau = tau(ig,1) !il faudrait iaerdust(iaer)
                 sl_alb = albedo(ig,l)
                 sl_psi = psi_sl(ig)
                 sl_fl0 = fluxsurf_sw(ig,l)
                 sl_di0 = 0.
                 if (mu0(ig) .gt. 0.) then
                  sl_di0 = mu0(ig)*(exp(-sl_tau/mu0(ig)))
                  sl_di0 = sl_di0*1370./dist_sol/dist_sol
                  sl_di0 = sl_di0/ztim1
                  sl_di0 = fluxsurf_sw(ig,l)*sl_di0
                 endif
                 ! you never know (roundup concern...)
                 if (sl_fl0 .lt. sl_di0) sl_di0=sl_fl0
                 !!!!!!!!!!!!!!!!!!!!!!!!!!
                 CALL param_slope( mu0(ig), declin, sl_ra, sl_lat, 
     &                             sl_tau, sl_alb, sl_the, sl_psi,
     &                             sl_di0, sl_fl0, sl_flu )
                 !!!!!!!!!!!!!!!!!!!!!!!!!!
                 fluxsurf_sw(ig,l) = sl_flu
                ENDDO
              !!! compute correction on IR flux as well
              sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
              fluxsurf_lw(ig)= fluxsurf_lw(ig)*sky
              ENDIF
            ENDDO
           ENDIF

c          CO2 near infrared absorption
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           zdtnirco2(:,:)=0
           if (callnirco2) then
              call nirco2abs (ngrid,nlayer,zplay,dist_sol,nq,pq,
     .                       mu0,fract,declin, zdtnirco2)
           endif

c          Radiative flux from the sky absorbed by the surface (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           DO ig=1,ngrid
               fluxrad_sky(ig)=emis(ig)*fluxsurf_lw(ig)
     $         +fluxsurf_sw(ig,1)*(1.-albedo(ig,1))
     $         +fluxsurf_sw(ig,2)*(1.-albedo(ig,2))
           ENDDO


c          Net atmospheric radiative heating rate (K.s-1)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           IF(callnlte) THEN
              CALL blendrad(ngrid, nlayer, zplay,
     &             zdtsw, zdtlw, zdtnirco2, zdtnlte, dtrad)
           ELSE
              DO l=1,nlayer
                 DO ig=1,ngrid
                    dtrad(ig,l)=zdtsw(ig,l)+zdtlw(ig,l)
     &                          +zdtnirco2(ig,l)
                  ENDDO
              ENDDO
           ENDIF

        ENDIF ! of if(mod(icount-1,iradia).eq.0)

c    Transformation of the radiative tendencies:
c    -------------------------------------------

c          Net radiative surface flux (W.m-2)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~
c
           DO ig=1,ngrid
               zplanck(ig)=tsurf(ig)*tsurf(ig)
               zplanck(ig)=emis(ig)*
     $         stephan*zplanck(ig)*zplanck(ig)
               fluxrad(ig)=fluxrad_sky(ig)-zplanck(ig)
               IF(callslope) THEN
                 sky= (1.+cos(pi*theta_sl(ig)/180.))/2.
                 fluxrad(ig)=fluxrad(ig)+(1.-sky)*zplanck(ig)
               ENDIF
           ENDDO

         DO l=1,nlayer
            DO ig=1,ngrid
               pdt(ig,l)=pdt(ig,l)+dtrad(ig,l)
            ENDDO
         ENDDO

      ENDIF ! of IF (callrad)

c     3.1 Rocket dust storm
c    -------------------------------------------
      IF (rdstorm) THEN
         clearatm=.false.
         pdqrds(:,:,:)=0.
         qdusttotal0(:)=0.
         qdusttotal1(:)=0.
         do ig=1,ngrid
           do l=1,nlayer
            zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l) ! updated potential
                                             ! temperature tendency
            ! for diagnostics
            qdustrds0(ig,l)=pq(ig,l,igcm_dust_mass)+
     &           pdq(ig,l,igcm_dust_mass)*ptimestep
            qstormrds0(ig,l)=pq(ig,l,igcm_stormdust_mass)+
     &           pdq(ig,l,igcm_stormdust_mass)*ptimestep
            qdusttotal0(ig)=qdusttotal0(ig)+(qdustrds0(ig,l)+
     &          qstormrds0(ig,l))*(zplev(ig,l)-
     &           zplev(ig,l+1))/g
           enddo
         enddo
         call writediagfi(ngrid,'qdustrds0','qdust before rds',
     &           'kg/kg ',3,qdustrds0)
         call writediagfi(ngrid,'qstormrds0','qstorm before rds',
     &           'kg/kg ',3,qstormrds0)

         CALL rocketduststorm(ngrid,nlayer,nq,ptime,ptimestep,
     &                      pq,pdq,pt,pdt,zplev,zplay,zzlev,
     &                      zzlay,zdtsw,zdtlw,
c               for radiative transfer
     &                      clearatm,icount,zday,zls,
     &                      tsurf,igout,totstormfract,
     &                      tauscaling,
c               input sub-grid scale cloud
     &                      clearsky,totcloudfrac,
c               input sub-grid scale topography
     &                      nohmons,alpha_hmons,
c               output
     &                      pdqrds,wspeed,dsodust,dsords,dsotop,
     &                      tauref)
                   
c      update the tendencies of both dust after vertical transport
         DO l=1,nlayer
          DO ig=1,ngrid
           pdq(ig,l,igcm_stormdust_mass)=
     &              pdq(ig,l,igcm_stormdust_mass)+ 
     &                      pdqrds(ig,l,igcm_stormdust_mass)
           pdq(ig,l,igcm_stormdust_number)=
     &              pdq(ig,l,igcm_stormdust_number)+
     &                      pdqrds(ig,l,igcm_stormdust_number)

           pdq(ig,l,igcm_dust_mass)=
     &       pdq(ig,l,igcm_dust_mass)+ pdqrds(ig,l,igcm_dust_mass)
           pdq(ig,l,igcm_dust_number)=
     &           pdq(ig,l,igcm_dust_number)+ 
     &                  pdqrds(ig,l,igcm_dust_number)

          ENDDO
         ENDDO
         do l=1,nlayer
           do ig=1,ngrid
             qdustrds1(ig,l)=pq(ig,l,igcm_dust_mass)+
     &           pdq(ig,l,igcm_dust_mass)*ptimestep
             qstormrds1(ig,l)=pq(ig,l,igcm_stormdust_mass)+
     &           pdq(ig,l,igcm_stormdust_mass)*ptimestep
             qdusttotal1(ig)=qdusttotal1(ig)+(qdustrds1(ig,l)+
     &          qstormrds1(ig,l))*(zplev(ig,l)-
     &           zplev(ig,l+1))/g
           enddo
         enddo

c        for diagnostics
         call writediagfi(ngrid,'qdustrds1','qdust after rds',
     &           'kg/kg ',3,qdustrds1)
         call writediagfi(ngrid,'qstormrds1','qstorm after rds',
     &           'kg/kg ',3,qstormrds1)
         
         call writediagfi(ngrid,'qdusttotal0','q sum before rds',
     &           'kg/kg ',2,qdusttotal0)
         call writediagfi(ngrid,'qdusttotal1','q sum after rds',
     &           'kg/kg ',2,qdusttotal1)

      ENDIF ! end of if(rdstorm)

c     3.2 Dust entrained from the PBL up to the top of sub-grid scale topography
c    -------------------------------------------
      IF (slpwind) THEN
         if (ngrid.gt.1) then
           hsummit(:)=summit(:)-phisfi(:)/g
         else
           hsummit(:)=14000.
         endif
         clearatm=.true. ! stormdust is not accounted in the extra heating on top of the mountains 
         nohmons=.false.
         pdqtop(:,:,:)=0.  
         CALL topmons(ngrid,nlayer,nq,ptime,ptimestep,
     &                pq,pdq,pt,pdt,zplev,zplay,zzlev,
     &                zzlay,zdtsw,zdtlw,
     &                icount,zday,zls,tsurf,igout,aerosol, 
     &                tauscaling,totstormfract,clearatm,
     &                clearsky,totcloudfrac,
     &                nohmons,hsummit,
     &                pdqtop,wtop,dsodust,dsords,dsotop,
     &                tauref)
     
                   
c      update the tendencies of both dust after vertical transport
         DO l=1,nlayer
          DO ig=1,ngrid
           pdq(ig,l,igcm_topdust_mass)=
     & pdq(ig,l,igcm_topdust_mass)+
     &                      pdqtop(ig,l,igcm_topdust_mass)
           pdq(ig,l,igcm_topdust_number)=
     & pdq(ig,l,igcm_topdust_number)+
     &                      pdqtop(ig,l,igcm_topdust_number)
           pdq(ig,l,igcm_dust_mass)=
     & pdq(ig,l,igcm_dust_mass)+ pdqtop(ig,l,igcm_dust_mass)
           pdq(ig,l,igcm_dust_number)=
     & pdq(ig,l,igcm_dust_number)+pdqtop(ig,l,igcm_dust_number)

          ENDDO
         ENDDO

      ENDIF ! end of if (slpwind)

c     3.3 Dust injection from the surface
c    -------------------------------------------
           if (dustinjection.gt.0) then

             CALL compute_dtau(ngrid,nlayer,
     &        			zday,pplev,tauref,
     &        			ptimestep,dustliftday,local_time)
           endif ! end of if (dustinjection.gt.0)

c-----------------------------------------------------------------------
c    4. Gravity wave and subgrid scale topography drag :
c    -------------------------------------------------


      IF(calllott)THEN

        CALL calldrag_noro(ngrid,nlayer,ptimestep,
     &                 zplay,zplev,pt,pu,pv,zdtgw,zdugw,zdvgw)
 
        DO l=1,nlayer
          DO ig=1,ngrid
            pdv(ig,l)=pdv(ig,l)+zdvgw(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdugw(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtgw(ig,l)
          ENDDO
        ENDDO
      ENDIF

c-----------------------------------------------------------------------
c    5. Vertical diffusion (turbulent mixing):
c    -----------------------------------------

      IF (calldifv) THEN

         DO ig=1,ngrid
            zflubid(ig)=fluxrad(ig)+fluxgrd(ig)
         ENDDO

         zdum1(:,:)=0
         zdum2(:,:)=0
         do l=1,nlayer
            do ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            enddo
         enddo

c ----------------------
c Treatment of a special case : using new surface layer (Richardson based)
c without using the thermals in gcm and mesoscale can yield problems in
c weakly unstable situations when winds are near to 0. For those cases, we add
c a unit subgrid gustiness. Remember that thermals should be used we using the 
c Richardson based surface layer model.
        IF ( .not.calltherm 
     .       .and. callrichsl 
     .       .and. .not.turb_resolved) THEN
          DO ig=1, ngrid
             IF (zh(ig,1) .lt. tsurf(ig)) THEN
               wstar(ig)=1.
               hfmax_th(ig)=0.2
             ELSE
               wstar(ig)=0.
               hfmax_th(ig)=0.
             ENDIF      
          ENDDO
        ENDIF
c ----------------------

         IF (tke_heat_flux .ne. 0.) THEN
             zz1(:)=(pt(:,1)+pdt(:,1)*ptimestep)*(r/g)*
     &                 (-alog(zplay(:,1)/zplev(:,1)))
             pdt(:,1)=pdt(:,1) + (tke_heat_flux/zz1(:))*zpopsk(:,1)
         ENDIF
c        Calling vdif (Martian version WITH CO2 condensation)
         CALL vdifc(ngrid,nlayer,nq,co2ice,zpopsk,
     $        ptimestep,capcal,lwrite,
     $        zplay,zplev,zzlay,zzlev,z0,
     $        pu,pv,zh,pq,tsurf,emis,qsurf,
     $        zdum1,zdum2,zdh,pdq,zflubid,
     $        zdudif,zdvdif,zdhdif,zdtsdif,q2,
     &        zdqdif,zdqsdif,wstar,zcdv,zcdh,hfmax_th,
     &        zcondicea_co2microp,sensibFlux,
     &        dustliftday,local_time,watercap,dwatercap_dif)

          DO ig=1,ngrid
             zdtsurf(ig)=zdtsurf(ig)+zdtsdif(ig)
             dwatercap(ig)=dwatercap(ig)+dwatercap_dif(ig)
          ENDDO

         IF (.not.turb_resolved) THEN
          DO l=1,nlayer
            DO ig=1,ngrid
               pdv(ig,l)=pdv(ig,l)+zdvdif(ig,l)
               pdu(ig,l)=pdu(ig,l)+zdudif(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhdif(ig,l)*zpopsk(ig,l)

               zdtdif(ig,l)=zdhdif(ig,l)*zpopsk(ig,l) ! for diagnostic only
            ENDDO
          ENDDO
          if (tracer) then
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdif(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
           DO iq=1, nq
              DO ig=1,ngrid
                 dqsurf(ig,iq)=dqsurf(ig,iq) + zdqsdif(ig,iq)
              ENDDO
           ENDDO
          end if ! of if (tracer)
         ELSE
           write (*,*) '******************************************'
           write (*,*) '** LES mode: the difv part is only used to'
           write (*,*) '**  - provide HFX and UST to the dynamics'
           write (*,*) '**  - update TSURF'
           write (*,*) '******************************************'
           !! Specific treatment for lifting in turbulent-resolving mode (AC)
           IF (lifting .and. doubleq) THEN
             !! lifted dust is injected in the first layer. 
             !! Sedimentation must be called after turbulent mixing, i.e. on next step, after WRF. 
             !! => lifted dust is not incremented before the sedimentation step.
             zdqdif(1:ngrid,1,1:nq)=0.
             zdqdif(1:ngrid,1,igcm_dust_number) = 
     .                  -zdqsdif(1:ngrid,igcm_dust_number)
             zdqdif(1:ngrid,1,igcm_dust_mass) = 
     .                  -zdqsdif(1:ngrid,igcm_dust_mass)
             zdqdif(1:ngrid,2:nlayer,1:nq) = 0.
             DO iq=1, nq
               IF ((iq .ne. igcm_dust_mass)
     &          .and. (iq .ne. igcm_dust_number)) THEN
                 zdqsdif(:,iq)=0.
               ENDIF
             ENDDO
           ELSE
             zdqdif(1:ngrid,1:nlayer,1:nq) = 0.
             zdqsdif(1:ngrid,1:nq) = 0.
           ENDIF
         ENDIF
      ELSE    
         DO ig=1,ngrid
            zdtsurf(ig)=zdtsurf(ig)+
     s        (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
         ENDDO
         IF (turb_resolved) THEN
            write(*,*) 'Turbulent-resolving mode !' 
            write(*,*) 'Please set calldifv to T in callphys.def'
            call abort_physic("physiq","turbulent-resolving mode",1)
         ENDIF
      ENDIF ! of IF (calldifv)

c-----------------------------------------------------------------------
c   6. Thermals :
c   -----------------------------

      if(calltherm .and. .not.turb_resolved) then
 
        call calltherm_interface(ngrid,nlayer,nq,
     $ tracer,igcm_co2,
     $ zzlev,zzlay,
     $ ptimestep,pu,pv,pt,pq,pdu,pdv,pdt,pdq,q2,
     $ zplay,zplev,pphi,zpopsk,
     $ pdu_th,pdv_th,pdt_th,pdq_th,lmax_th,zmax_th,
     $ dtke_th,zdhdif,hfmax_th,wstar,sensibFlux)
     
         DO l=1,nlayer
           DO ig=1,ngrid
              pdu(ig,l)=pdu(ig,l)+pdu_th(ig,l)
              pdv(ig,l)=pdv(ig,l)+pdv_th(ig,l)
              pdt(ig,l)=pdt(ig,l)+pdt_th(ig,l)
              q2(ig,l)=q2(ig,l)+dtke_th(ig,l)*ptimestep
           ENDDO
        ENDDO
 
        DO ig=1,ngrid
          q2(ig,nlayer+1)=q2(ig,nlayer+1)+dtke_th(ig,nlayer+1)*ptimestep
        ENDDO      
    
        if (tracer) then
        DO iq=1,nq
         DO l=1,nlayer
           DO ig=1,ngrid
             pdq(ig,l,iq)=pdq(ig,l,iq)+pdq_th(ig,l,iq)
           ENDDO
         ENDDO
        ENDDO
        endif

        lmax_th_out(:)=real(lmax_th(:))

      else   !of if calltherm
        lmax_th(:)=0
        wstar(:)=0.
        hfmax_th(:)=0.
        lmax_th_out(:)=0.
      end if
c-----------------------------------------------------------------------
c   7. Dry convective adjustment:
c   -----------------------------

      IF(calladj) THEN

         DO l=1,nlayer
            DO ig=1,ngrid
               zdh(ig,l)=pdt(ig,l)/zpopsk(ig,l)
            ENDDO
         ENDDO
         zduadj(:,:)=0
         zdvadj(:,:)=0
         zdhadj(:,:)=0

         CALL convadj(ngrid,nlayer,nq,ptimestep,
     $                zplay,zplev,zpopsk,lmax_th,
     $                pu,pv,zh,pq,
     $                pdu,pdv,zdh,pdq,
     $                zduadj,zdvadj,zdhadj,
     $                zdqadj)

         DO l=1,nlayer
            DO ig=1,ngrid
               pdu(ig,l)=pdu(ig,l)+zduadj(ig,l)
               pdv(ig,l)=pdv(ig,l)+zdvadj(ig,l)
               pdt(ig,l)=pdt(ig,l)+zdhadj(ig,l)*zpopsk(ig,l)

               zdtadj(ig,l)=zdhadj(ig,l)*zpopsk(ig,l) ! for diagnostic only
            ENDDO
         ENDDO

         if(tracer) then 
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqadj(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
         end if
      ENDIF ! of IF(calladj)

c-----------------------------------------------------
c    8. Non orographic Gravity waves :
c -------------------------------------------------

      IF (calllott_nonoro) THEN

         CALL nonoro_gwd_ran(ngrid,nlayer,ptimestep,zplay,
     &               zmax_th,                      ! max altitude reached by thermals (m)
     &               pt, pu, pv,
     &               pdt, pdu, pdv,
     &               zustrhi,zvstrhi,
     &               d_t_hin, d_u_hin, d_v_hin)

!  Update tendencies
         pdt(1:ngrid,1:nlayer)=pdt(1:ngrid,1:nlayer)
     &                         +d_t_hin(1:ngrid,1:nlayer)
!        d_t_hin(:,:)= d_t_hin(:,:)/ptimestep ! K/s (for outputs?)
         pdu(1:ngrid,1:nlayer)=pdu(1:ngrid,1:nlayer)
     &                         +d_u_hin(1:ngrid,1:nlayer) 
!        d_u_hin(:,:)= d_u_hin(:,:)/ptimestep ! (m/s)/s (for outputs?)
         pdv(1:ngrid,1:nlayer)=pdv(1:ngrid,1:nlayer)
     &                         +d_v_hin(1:ngrid,1:nlayer)
!        d_v_hin(:,:)= d_v_hin(:,:)/ptimestep ! (m/s)/s (for outputs?)

      ENDIF ! of IF (calllott_nonoro)

c-----------------------------------------------------------------------
c   9. Specific parameterizations for tracers 
c:   -----------------------------------------

      if (tracer) then 

c   9a. Water and ice
c     ---------------

c        ---------------------------------------
c        Water ice condensation in the atmosphere
c        ----------------------------------------
         IF (water) THEN

           call watercloud(ngrid,nlayer,ptimestep,
     &                zplev,zplay,pdpsrf,zzlay, pt,pdt,
     &                pq,pdq,zdqcloud,zdtcloud,
     &                nq,tau,tauscaling,rdust,rice,nuice,
     &                rsedcloud,rhocloud,totcloudfrac)
     
c Temperature variation due to latent heat release
           if (activice) then
               pdt(1:ngrid,1:nlayer) = 
     &          pdt(1:ngrid,1:nlayer) + 
     &          zdtcloud(1:ngrid,1:nlayer)
           endif
           
! increment water vapour and ice atmospheric tracers tendencies
           pdq(1:ngrid,1:nlayer,igcm_h2o_vap) = 
     &       pdq(1:ngrid,1:nlayer,igcm_h2o_vap) + 
     &       zdqcloud(1:ngrid,1:nlayer,igcm_h2o_vap)
           pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = 
     &       pdq(1:ngrid,1:nlayer,igcm_h2o_ice) + 
     &       zdqcloud(1:ngrid,1:nlayer,igcm_h2o_ice)

! increment dust and ccn masses and numbers
! We need to check that we have Nccn & Ndust > 0
! This is due to single precision rounding problems
           if (microphys) then
             pdq(1:ngrid,1:nlayer,igcm_ccn_mass) = 
     &         pdq(1:ngrid,1:nlayer,igcm_ccn_mass) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_ccn_mass)
             pdq(1:ngrid,1:nlayer,igcm_ccn_number) = 
     &         pdq(1:ngrid,1:nlayer,igcm_ccn_number) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_ccn_number)
             where (pq(:,:,igcm_ccn_mass) + 
     &       ptimestep*pdq(:,:,igcm_ccn_mass) < 0.)
               pdq(:,:,igcm_ccn_mass) = 
     &         - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_ccn_number) = 
     &         - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
             end where
             where (pq(:,:,igcm_ccn_number) + 
     &       ptimestep*pdq(:,:,igcm_ccn_number) < 0.)
               pdq(:,:,igcm_ccn_mass) = 
     &         - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_ccn_number) = 
     &         - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
             end where
           endif

           if (scavenging) then
             pdq(1:ngrid,1:nlayer,igcm_dust_mass) = 
     &         pdq(1:ngrid,1:nlayer,igcm_dust_mass) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_dust_mass)
             pdq(1:ngrid,1:nlayer,igcm_dust_number) = 
     &         pdq(1:ngrid,1:nlayer,igcm_dust_number) + 
     &         zdqcloud(1:ngrid,1:nlayer,igcm_dust_number)
             where (pq(:,:,igcm_dust_mass) + 
     &       ptimestep*pdq(:,:,igcm_dust_mass) < 0.)
               pdq(:,:,igcm_dust_mass) = 
     &         - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_dust_number) = 
     &         - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
             end where
             where (pq(:,:,igcm_dust_number) + 
     &       ptimestep*pdq(:,:,igcm_dust_number) < 0.)
               pdq(:,:,igcm_dust_mass) = 
     &         - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_dust_number) = 
     &         - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
             end where
           endif ! of if scavenging
                      
         END IF  ! of IF (water)

c   9a bis. CO2 clouds (CL & JA)
c        ---------------------------------------
c        CO2 ice cloud condensation in the atmosphere
c        ----------------------------------------
c flag needed in callphys.def:
c               co2clouds=.true. is mandatory (default is .false.)
c               co2useh2o=.true. if you want to allow co2 condensation 
c                                on water ice particles 
c               meteo_flux=.true. if you want to add a meteoritic
c                                 supply of CCN
c               CLFvaryingCO2=.true. if you want to have a sub-grid
c                                    temperature distribution 
c               spantCO2=integer (i.e. 3) amplitude of the sub-grid T disti
c               nuiceco2_sed=0.2 variance of the size distribution for the 
c                                sedimentation
c               nuiceco2_ref=0.2 variance of the size distribution for the 
c                                nucleation
c               imicroco2=50     micro-timestep is 1/50 of physical timestep
         
         IF (co2clouds ) THEN


            call co2cloud(ngrid,nlayer,ptimestep,
     &           zplev,zplay,pdpsrf,zzlay,pt,pdt,
     &           pq,pdq,zdqcloudco2,zdtcloudco2,
     &           nq,tau,tauscaling,rdust,rice,riceco2,nuice,
     &           rsedcloudco2,rhocloudco2,
     &           rsedcloud,rhocloud,zzlev,zdqssed_co2,
     &           pdu,pu,zcondicea_co2microp)


c Temperature variation due to latent heat release
c            if (activice) then  !Maybe create activice_co2 ?
               pdt(1:ngrid,1:nlayer) = 
     &              pdt(1:ngrid,1:nlayer) + 
     &              zdtcloudco2(1:ngrid,1:nlayer)! --> in co2condens
c            endif
            

! increment dust and ccn masses and numbers
! We need to check that we have Nccn & Ndust > 0
! This is due to single precision rounding problems
              
! increment dust tracers tendancies
               pdq(1:ngrid,1:nlayer,igcm_dust_mass) = 
     &              pdq(1:ngrid,1:nlayer,igcm_dust_mass) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_dust_mass)          
               pdq(1:ngrid,1:nlayer,igcm_dust_number) = 
     &              pdq(1:ngrid,1:nlayer,igcm_dust_number) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_dust_number)
               pdq(1:ngrid,1:nlayer,igcm_co2) = 
     &              pdq(1:ngrid,1:nlayer,igcm_co2) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_co2)
               pdq(1:ngrid,1:nlayer,igcm_co2_ice) = 
     &              pdq(1:ngrid,1:nlayer,igcm_co2_ice) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_co2_ice)
               pdq(1:ngrid,1:nlayer,igcm_ccnco2_mass) = 
     &              pdq(1:ngrid,1:nlayer,igcm_ccnco2_mass) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_ccnco2_mass)
               pdq(1:ngrid,1:nlayer,igcm_ccnco2_number) = 
     &              pdq(1:ngrid,1:nlayer,igcm_ccnco2_number) + 
     &              zdqcloudco2(1:ngrid,1:nlayer,igcm_ccnco2_number)
!Update water ice clouds values as well
             if (co2useh2o) then
                pdq(1:ngrid,1:nlayer,igcm_h2o_ice) = 
     &               pdq(1:ngrid,1:nlayer,igcm_h2o_ice) + 
     &               zdqcloudco2(1:ngrid,1:nlayer,igcm_h2o_ice)
                pdq(1:ngrid,1:nlayer,igcm_ccn_mass) = 
     &               pdq(1:ngrid,1:nlayer,igcm_ccn_mass) + 
     &               zdqcloudco2(1:ngrid,1:nlayer,igcm_ccn_mass)
                pdq(1:ngrid,1:nlayer,igcm_ccn_number) = 
     &               pdq(1:ngrid,1:nlayer,igcm_ccn_number) + 
     &               zdqcloudco2(1:ngrid,1:nlayer,igcm_ccn_number)
                where (pq(:,:,igcm_ccn_mass) + 
     &                 ptimestep*pdq(:,:,igcm_ccn_mass) < 0.)
                  pdq(:,:,igcm_ccn_mass) = 
     &               - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
                  pdq(:,:,igcm_ccn_number) = 
     &               - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
                end where
                where (pq(:,:,igcm_ccn_number) + 
     &                 ptimestep*pdq(:,:,igcm_ccn_number) < 0.)
                  pdq(:,:,igcm_ccn_mass) = 
     &              - pq(:,:,igcm_ccn_mass)/ptimestep + 1.e-30
                  pdq(:,:,igcm_ccn_number) = 
     &              - pq(:,:,igcm_ccn_number)/ptimestep + 1.e-30
                end where
             endif ! of if (co2useh2o)
c     Negative values?
             where (pq(:,:,igcm_ccnco2_mass) + 
     &              ptimestep*pdq(:,:,igcm_ccnco2_mass) < 0.)
               pdq(:,:,igcm_ccnco2_mass) = 
     &            - pq(:,:,igcm_ccnco2_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_ccnco2_number) = 
     &            - pq(:,:,igcm_ccnco2_number)/ptimestep + 1.e-30
             end where
             where (pq(:,:,igcm_ccnco2_number) + 
     &              ptimestep*pdq(:,:,igcm_ccnco2_number) < 0.)
               pdq(:,:,igcm_ccnco2_mass) = 
     &            - pq(:,:,igcm_ccnco2_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_ccnco2_number) = 
     &            - pq(:,:,igcm_ccnco2_number)/ptimestep + 1.e-30
             end where
       
c     Negative values?
             where (pq(:,:,igcm_dust_mass) + 
     &              ptimestep*pdq(:,:,igcm_dust_mass) < 0.)
               pdq(:,:,igcm_dust_mass) = 
     &           - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_dust_number) = 
     &           - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
             end where
             where (pq(:,:,igcm_dust_number) + 
     &              ptimestep*pdq(:,:,igcm_dust_number) < 0.)
               pdq(:,:,igcm_dust_mass) = 
     &           - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
               pdq(:,:,igcm_dust_number) = 
     &           - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
             end where
      
      END IF                    ! of IF (co2clouds)

c   9b. Aerosol particles
c     -------------------
c        ----------
c        Dust devil :
c        ----------
         IF(callddevil) then 
           call dustdevil(ngrid,nlayer,nq, zplev,pu,pv,pt, tsurf,q2,
     &                zdqdev,zdqsdev)
 
           if (dustbin.ge.1) then
              do iq=1,nq
                 DO l=1,nlayer
                    DO ig=1,ngrid
                       pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdev(ig,l,iq)
                    ENDDO
                 ENDDO
              enddo
              do iq=1,nq
                 DO ig=1,ngrid
                    dqsurf(ig,iq)= dqsurf(ig,iq) + zdqsdev(ig,iq)
                 ENDDO
              enddo
           endif  ! of if (dustbin.ge.1)

         END IF ! of IF (callddevil)

c        ------------- 
c        Sedimentation :   acts also on water ice
c        -------------
         IF (sedimentation) THEN 
           zdqsed(1:ngrid,1:nlayer,1:nq)=0
           zdqssed(1:ngrid,1:nq)=0

c Sedimentation for co2 clouds tracers are inside co2cloud microtimestep
c Zdqssed isn't
           call callsedim(ngrid,nlayer,ptimestep,
     &                zplev,zzlev,zzlay,pt,pdt,
     &                rdust,rstormdust,rtopdust,
     &                rice,rsedcloud,rhocloud,
     &                pq,pdq,zdqsed,zdqssed,nq, 
     &                tau,tauscaling)
c Flux at the surface of co2 ice computed in co2cloud microtimestep
           IF (co2clouds) THEN
              zdqssed(1:ngrid,igcm_co2_ice)=zdqssed_co2(1:ngrid)
           ENDIF

           IF (rdstorm) THEN
c             Storm dust cannot sediment to the surface
              DO ig=1,ngrid  
                 zdqsed(ig,1,igcm_stormdust_mass)=
     &             zdqsed(ig,1,igcm_stormdust_mass)+
     &             zdqssed(ig,igcm_stormdust_mass) /
     &             ((pplev(ig,1)-pplev(ig,2))/g)
                 zdqsed(ig,1,igcm_stormdust_number)=
     &             zdqsed(ig,1,igcm_stormdust_number)+
     &             zdqssed(ig,igcm_stormdust_number) /
     &               ((pplev(ig,1)-pplev(ig,2))/g) 
                 zdqssed(ig,igcm_stormdust_mass)=0.
                 zdqssed(ig,igcm_stormdust_number)=0.
              ENDDO
           ENDIF !rdstorm

           DO iq=1, nq
             DO l=1,nlayer
               DO ig=1,ngrid
                    pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqsed(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO
           DO iq=1, nq
             DO ig=1,ngrid
                dqsurf(ig,iq)= dqsurf(ig,iq) + zdqssed(ig,iq)
             ENDDO
           ENDDO

         END IF   ! of IF (sedimentation)

c Add lifted dust to tendancies after sedimentation in the LES (AC)
      IF (turb_resolved) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqdif(ig,l,iq)
              ENDDO
            ENDDO
           ENDDO
           DO iq=1, nq
              DO ig=1,ngrid
                 dqsurf(ig,iq)=dqsurf(ig,iq) + zdqsdif(ig,iq)
              ENDDO
           ENDDO
      ENDIF
c
c   9c. Chemical species
c     ------------------

#ifndef MESOSCALE
c        --------------
c        photochemistry :
c        --------------
         IF (photochem) then

           if (modulo(icount-1,ichemistry).eq.0) then
           ! compute chemistry every ichemistry physics step

!           dust and ice surface area
            call surfacearea(ngrid, nlayer, naerkind,
     $                       ptimestep, zplay, zzlay, 
     $                       pt, pq, pdq, nq, 
     $                       rdust, rice, tau, tauscaling, 
     $                       surfdust, surfice)
!           call photochemistry
            call calchim(ngrid,nlayer,nq,
     &                   ptimestep,zplay,zplev,pt,pdt,dist_sol,mu0,
     $                   zzlev,zzlay,zday,pq,pdq,zdqchim,zdqschim,
     $                   zdqcloud,zdqscloud,tau(:,1),co2ice,
     $                   pu,pdu,pv,pdv,surfdust,surfice)

            endif ! of if (modulo(icount-1,ichemistry).eq.0)

           ! increment values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO l=1,nlayer
               DO ig=1,ngrid
                 pdq(ig,l,iq)=pdq(ig,l,iq)+zdqchim(ig,l,iq)
               ENDDO
             ENDDO
           ENDDO ! of DO iq=1,nq
           
           ! add condensation tendency for H2O2
           if (igcm_h2o2.ne.0) then
             DO l=1,nlayer
               DO ig=1,ngrid
                 pdq(ig,l,igcm_h2o2)=pdq(ig,l,igcm_h2o2)
     &                                +zdqcloud(ig,l,igcm_h2o2)
               ENDDO
             ENDDO
           endif

           ! increment surface values of tracers:
           DO iq=1,nq ! loop on all tracers; tendencies for non-chemistry
                      ! tracers is zero anyways
             DO ig=1,ngrid
               dqsurf(ig,iq)=dqsurf(ig,iq)+zdqschim(ig,iq)
             ENDDO
           ENDDO ! of DO iq=1,nq

           ! add condensation tendency for H2O2
           if (igcm_h2o2.ne.0) then
             DO ig=1,ngrid
               dqsurf(ig,igcm_h2o2)=dqsurf(ig,igcm_h2o2)
     &                                +zdqscloud(ig,igcm_h2o2)
             ENDDO
           endif

         END IF  ! of IF (photochem)
#endif

      endif !  of if (tracer) 

#ifndef MESOSCALE
c-----------------------------------------------------------------------
c   10. THERMOSPHERE CALCULATION
c-----------------------------------------------------------------------

      if (callthermos) then
        call thermosphere(ngrid,nlayer,nq,zplev,zplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)

        DO l=1,nlayer
          DO ig=1,ngrid
            dtrad(ig,l)=dtrad(ig,l)+zdteuv(ig,l)
            pdt(ig,l)=pdt(ig,l)+zdtconduc(ig,l)
     &                         +zdteuv(ig,l)
            pdv(ig,l)=pdv(ig,l)+zdvmolvis(ig,l)
            pdu(ig,l)=pdu(ig,l)+zdumolvis(ig,l)
            DO iq=1, nq
              pdq(ig,l,iq)=pdq(ig,l,iq)+zdqmoldiff(ig,l,iq)
            ENDDO
          ENDDO
        ENDDO

      endif ! of if (callthermos)
#endif
c-----------------------------------------------------------------------
c   11. Carbon dioxide condensation-sublimation:
c     (should be the last atmospherical physical process to be computed)
c   -------------------------------------------

      IF (tituscap) THEN
         !!! get the actual co2 seasonal cap from Titus observations
         CALL geticecover( ngrid, 180.*zls/pi,
     .                  180.*longitude/pi, 180.*latitude/pi, co2ice )
         co2ice = co2ice * 10000.
      ENDIF
      
      
      pdpsrf(:) = 0

      IF (callcond) THEN
         CALL co2condens(ngrid,nlayer,nq,ptimestep,
     $              capcal,zplay,zplev,tsurf,pt,
     $              pphi,pdt,pdu,pdv,zdtsurf,pu,pv,pq,pdq,
     $              co2ice,albedo,emis,rdust,
     $              zdtc,zdtsurfc,pdpsrf,zduc,zdvc,zdqc,
     $              fluxsurf_sw,zls,
     $              zdqssed_co2,zcondicea_co2microp,
     &              zdtcloudco2,zdqsc)

         DO l=1,nlayer
           DO ig=1,ngrid
             pdt(ig,l)=pdt(ig,l)+zdtc(ig,l)
             pdv(ig,l)=pdv(ig,l)+zdvc(ig,l)
             pdu(ig,l)=pdu(ig,l)+zduc(ig,l)
           ENDDO
         ENDDO
         DO ig=1,ngrid
            zdtsurf(ig) = zdtsurf(ig) + zdtsurfc(ig)
         ENDDO

         IF (tracer) THEN
           DO iq=1, nq
            DO l=1,nlayer
              DO ig=1,ngrid
                pdq(ig,l,iq)=pdq(ig,l,iq)+ zdqc(ig,l,iq) 
              ENDDO
            ENDDO
           ENDDO
	   	   
	   DO iq=1, nq
            DO ig=1,ngrid
              dqsurf(ig,iq)=dqsurf(ig,iq)+zdqsc(ig,iq)
            ENDDO  ! (ig)
           ENDDO    ! (iq)
	   
         ENDIF ! of IF (tracer)

#ifndef MESOSCALE
        ! update surface pressure
        DO ig=1,ngrid
          ps(ig) = zplev(ig,1) + pdpsrf(ig)*ptimestep
        ENDDO
        
        ! update pressure levels
        DO l=1,nlayer
         DO ig=1,ngrid
          zplay(ig,l) = aps(l) + bps(l)*ps(ig)
          zplev(ig,l) = ap(l)  + bp(l)*ps(ig)
         ENDDO
        ENDDO
        zplev(:,nlayer+1) = 0.
        
        ! update layers altitude
        DO l=2,nlayer
          DO ig=1,ngrid
           z1=(zplay(ig,l-1)+zplev(ig,l))/(zplay(ig,l-1)-zplev(ig,l))
           z2=(zplev(ig,l)+zplay(ig,l))/(zplev(ig,l)-zplay(ig,l))
           zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
          ENDDO
        ENDDO
#endif
      
      ENDIF  ! of IF (callcond)

c-----------------------------------------------------------------------
c  Updating tracer budget on surface
c-----------------------------------------------------------------------	 
      IF (tracer) THEN
        DO iq=1, nq
          DO ig=1,ngrid

            qsurf(ig,iq)=qsurf(ig,iq)+ptimestep*dqsurf(ig,iq)

          ENDDO  ! (ig)
        ENDDO    ! (iq)
      ENDIF
      
c-----------------------------------------------------------------------
c   12. Surface  and sub-surface soil temperature
c-----------------------------------------------------------------------
c
c
c   12.1 Increment Surface temperature:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO ig=1,ngrid
         tsurf(ig)=tsurf(ig)+ptimestep*zdtsurf(ig) 
      ENDDO

c  Prescribe a cold trap at south pole (except at high obliquity !!)
c  Temperature at the surface is set there to be the temperature
c  corresponding to equilibrium temperature between phases of CO2


      IF (tracer.AND.water.AND.(ngrid.NE.1)) THEN
!#ifndef MESOSCALE 
!         if (caps.and.(obliquit.lt.27.)) then => now done in co2condens
           ! NB: Updated surface pressure, at grid point 'ngrid', is
           !     ps(ngrid)=zplev(ngrid,1)+pdpsrf(ngrid)*ptimestep
!           tsurf(ngrid)=1./(1./136.27-r/5.9e+5*alog(0.0095*
!     &                     (zplev(ngrid,1)+pdpsrf(ngrid)*ptimestep)))
!           tsurf(ngrid)=1./(1./136.27-r/5.9e+5*alog(0.0095*ps(ngrid)))
!         endif
!#endif
c       -------------------------------------------------------------
c       Change of surface albedo in case of ground frost
c       everywhere except on the north permanent cap and in regions
c       covered by dry ice. 
c              ALWAYS PLACE these lines after co2condens !!!
c       -------------------------------------------------------------
         do ig=1,ngrid
           if ((co2ice(ig).eq.0).and.
     &        (qsurf(ig,igcm_h2o_ice).gt.frost_albedo_threshold)) then
              albedo(ig,1) = albedo_h2o_ice
              albedo(ig,2) = albedo_h2o_ice
c              write(*,*) "frost thickness", qsurf(ig,igcm_h2o_ice)
c              write(*,*) "physiq.F frost :"
c     &        ,latitude(ig)*180./pi, longitude(ig)*180./pi
           endif
         enddo  ! of do ig=1,ngrid
      ENDIF  ! of IF (tracer.AND.water.AND.(ngrid.NE.1))

c
c   12.2 Compute soil temperatures and subsurface heat flux:
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (callsoil) THEN
c       Thermal inertia feedback
        IF (tifeedback) THEN
         CALL soil_tifeedback(ngrid,nsoilmx,qsurf,inertiesoil)
         CALL soil(ngrid,nsoilmx,.false.,inertiesoil,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
        ELSE
         CALL soil(ngrid,nsoilmx,.false.,inertiedat,
     s          ptimestep,tsurf,tsoil,capcal,fluxgrd)
        ENDIF
      ENDIF

c     To avoid negative values
      IF (rdstorm) THEN
           where (pq(:,:,igcm_stormdust_mass) + 
     &      ptimestep*pdq(:,:,igcm_stormdust_mass) < 0.)
             pdq(:,:,igcm_stormdust_mass) = 
     &       - pq(:,:,igcm_stormdust_mass)/ptimestep + 1.e-30
             pdq(:,:,igcm_stormdust_number) = 
     &       - pq(:,:,igcm_stormdust_number)/ptimestep + 1.e-30
           end where
           where (pq(:,:,igcm_stormdust_number) + 
     &      ptimestep*pdq(:,:,igcm_stormdust_number) < 0.)
             pdq(:,:,igcm_stormdust_mass) = 
     &       - pq(:,:,igcm_stormdust_mass)/ptimestep + 1.e-30
             pdq(:,:,igcm_stormdust_number) = 
     &       - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
           end where

            where (pq(:,:,igcm_dust_mass) + 
     &      ptimestep*pdq(:,:,igcm_dust_mass) < 0.)
             pdq(:,:,igcm_dust_mass) = 
     &       - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
             pdq(:,:,igcm_dust_number) = 
     &       - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
           end where
           where (pq(:,:,igcm_dust_number) + 
     &      ptimestep*pdq(:,:,igcm_dust_number) < 0.)
             pdq(:,:,igcm_dust_mass) = 
     &       - pq(:,:,igcm_dust_mass)/ptimestep + 1.e-30
             pdq(:,:,igcm_dust_number) = 
     &       - pq(:,:,igcm_dust_number)/ptimestep + 1.e-30
           end where
      ENDIF !(rdstorm)   

c-----------------------------------------------------------------------
c   J. Naar : Surface and sub-surface water ice
c-----------------------------------------------------------------------
c
c
c   Increment Watercap (surface h2o reservoirs):
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      DO ig=1,ngrid
         watercap(ig)=watercap(ig)+ptimestep*dwatercap(ig) 
      ENDDO

c-----------------------------------------------------------------------
c  13. Write output files
c  ----------------------

c    -------------------------------
c    Dynamical fields incrementation
c    -------------------------------
c (FOR OUTPUT ONLY : the actual model integration is performed in the dynamics)
      ! temperature, zonal and meridional wind
      DO l=1,nlayer
        DO ig=1,ngrid
          zt(ig,l)=pt(ig,l)  + pdt(ig,l)*ptimestep
          zu(ig,l)=pu(ig,l)  + pdu(ig,l)*ptimestep
          zv(ig,l)=pv(ig,l)  + pdv(ig,l)*ptimestep
        ENDDO
      ENDDO

      ! tracers
      DO iq=1, nq
        DO l=1,nlayer
          DO ig=1,ngrid
            zq(ig,l,iq)=pq(ig,l,iq) +pdq(ig,l,iq)*ptimestep
          ENDDO
        ENDDO
      ENDDO

      ! Density 
      DO l=1,nlayer
         DO ig=1,ngrid
            rho(ig,l) = zplay(ig,l)/(rnew(ig,l)*zt(ig,l))
         ENDDO
      ENDDO

      ! Potential Temperature

       DO ig=1,ngrid
          DO l=1,nlayer
              zh(ig,l) = zt(ig,l)*(zplev(ig,1)/zplay(ig,l))**rcp
          ENDDO
       ENDDO


c    Compute surface stress : (NB: z0 is a common in surfdat.h)
c     DO ig=1,ngrid
c        cd = (0.4/log(zzlay(ig,1)/z0(ig)))**2
c        zstress(ig) = rho(ig,1)*cd*(zu(ig,1)**2 + zv(ig,1)**2)
c     ENDDO

c     Sum of fluxes in solar spectral bands (for output only)
      DO ig=1,ngrid
             fluxtop_sw_tot(ig)=fluxtop_sw(ig,1) + fluxtop_sw(ig,2)
             fluxsurf_sw_tot(ig)=fluxsurf_sw(ig,1) + fluxsurf_sw(ig,2)
      ENDDO
c ******* TEST ******************************************************
      ztim1 = 999
      DO l=1,nlayer
        DO ig=1,ngrid
           if (pt(ig,l).lt.ztim1) then
               ztim1 = pt(ig,l)
               igmin = ig
               lmin = l 
           end if
        ENDDO
      ENDDO
      if(min(pt(igmin,lmin),zt(igmin,lmin)).lt.70.) then
        write(*,*) 'PHYSIQ: stability WARNING :'
        write(*,*) 'pt, zt Tmin = ', pt(igmin,lmin), zt(igmin,lmin),
     &              'ig l =', igmin, lmin
      end if
c *******************************************************************

c     ---------------------
c     Outputs to the screen 
c     ---------------------

      IF (lwrite) THEN
         PRINT*,'Global diagnostics for the physics'
         PRINT*,'Variables and their increments x and dx/dt * dt'
         WRITE(*,'(a6,a10,2a15)') 'Ts','dTs','ps','dps'
         WRITE(*,'(2f10.5,2f15.5)')
     s   tsurf(igout),zdtsurf(igout)*ptimestep,
     s   zplev(igout,1),pdpsrf(igout)*ptimestep
         WRITE(*,'(a4,a6,5a10)') 'l','u','du','v','dv','T','dT'
         WRITE(*,'(i4,6f10.5)') (l,
     s   pu(igout,l),pdu(igout,l)*ptimestep,
     s   pv(igout,l),pdv(igout,l)*ptimestep,
     s   pt(igout,l),pdt(igout,l)*ptimestep,
     s   l=1,nlayer)
      ENDIF ! of IF (lwrite)

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        INTERPOLATIONS IN THE SURFACE-LAYER
c        ----------------------------------------------------------
c        ----------------------------------------------------------

           n_out=0 ! number of elements in the z_out array.
                   ! for z_out=[3.,2.,1.,0.5,0.1], n_out must be set
                   ! to 5
           IF (n_out .ne. 0) THEN

           IF(.NOT. ALLOCATED(z_out)) ALLOCATE(z_out(n_out))
           IF(.NOT. ALLOCATED(T_out)) ALLOCATE(T_out(ngrid,n_out))
           IF(.NOT. ALLOCATED(u_out)) ALLOCATE(u_out(ngrid,n_out))

           z_out(:)=[3.,2.,1.,0.5,0.1]
           u_out(:,:)=0.
           T_out(:,:)=0.

           call pbl_parameters(ngrid,nlayer,ps,zplay,z0,
     & g,zzlay,zzlev,zu,zv,wstar,hfmax_th,zmax_th,tsurf,zh,z_out,n_out,
     & T_out,u_out,ustar,tstar,L_mo,vhf,vvv)
                   ! pourquoi ustar recalcule ici? fait dans vdifc.

#ifndef MESOSCALE
            IF (ngrid .eq. 1) THEN
               dimout=0
            ELSE
               dimout=2
            ENDIF
            DO n=1,n_out
               write(zstring, '(F8.6)') z_out(n)
               call WRITEDIAGFI(ngrid,'T_out_'//trim(zstring),
     &   'potential temperature at z_out','K',dimout,T_out(:,n))
               call WRITEDIAGFI(ngrid,'u_out_'//trim(zstring),
     &   'horizontal velocity norm at z_out','m/s',dimout,u_out(:,n))
            ENDDO
            call WRITEDIAGFI(ngrid,'u_star',
     &   'friction velocity','m/s',dimout,ustar)
            call WRITEDIAGFI(ngrid,'teta_star',
     &   'friction potential temperature','K',dimout,tstar)
!            call WRITEDIAGFI(ngrid,'L',
!     &   'Monin Obukhov length','m',dimout,L_mo)
            call WRITEDIAGFI(ngrid,'vvv',
     &   'Vertical velocity variance at zout','m',dimout,vvv)
            call WRITEDIAGFI(ngrid,'vhf',
     &   'Vertical heat flux at zout','m',dimout,vhf)
#else
         T_out1(:)=T_out(:,1)
         u_out1(:)=u_out(:,1)
#endif

         ENDIF

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        END OF SURFACE LAYER INTERPOLATIONS
c        ----------------------------------------------------------
c        ----------------------------------------------------------

      IF (ngrid.NE.1) THEN

#ifndef MESOSCALE
c        -------------------------------------------------------------------
c        Writing NetCDF file  "RESTARTFI" at the end of the run
c        -------------------------------------------------------------------
c        Note: 'restartfi' is stored just before dynamics are stored
c              in 'restart'. Between now and the writting of 'restart',
c              there will have been the itau=itau+1 instruction and
c              a reset of 'time' (lastacll = .true. when itau+1= itaufin)
c              thus we store for time=time+dtvr

         IF(   ((ecritstart.GT.0) .and. 
     .          (MOD(icount*iphysiq,ecritstart).EQ.0)) 
     .    .or. lastcall  ) THEN
         
          ztime_fin = pday + ptime  + ptimestep/(float(iphysiq)*daysec)
     .               - day_ini - time_phys
          print*, pday,ptime,day_ini, time_phys
          write(*,'(A,I7,A,F12.5)') 
     .         'PHYSIQ: Ecriture du fichier restartfi ; icount=',
     .          icount,' date=',ztime_fin
            
            
          call physdem1("restartfi.nc",nsoilmx,ngrid,nlayer,nq,
     .                ptimestep,ztime_fin,
     .                tsurf,tsoil,co2ice,albedo,emis,
     .                q2,qsurf,tauscaling,totcloudfrac,wstar,
     .                mem_Mccn_co2,mem_Nccn_co2,mem_Mh2o_co2,watercap)
          
         ENDIF
#endif

c        -------------------------------------------------------------------
c        Calculation of diagnostic variables written in both stats and
c          diagfi files
c        -------------------------------------------------------------------

         if (tracer) then
           ! Density-scaled opacities
              do ig=1,ngrid 
                dsodust(ig,:) = 
     &                dsodust(ig,:)*tauscaling(ig)
                dsords(ig,:) = 
     &                dsords(ig,:)*tauscaling(ig)
                dsotop(ig,:) =
     &                dsotop(ig,:)*tauscaling(ig)
              enddo
           if(doubleq) then
              do ig=1,ngrid 
                dqdustsurf(ig) = 
     &                zdqssed(ig,igcm_dust_mass)*tauscaling(ig)
                dndustsurf(ig) = 
     &                zdqssed(ig,igcm_dust_number)*tauscaling(ig)
                ndust(ig,:) =
     &                zq(ig,:,igcm_dust_number)*tauscaling(ig)
                qdust(ig,:) =
     &                zq(ig,:,igcm_dust_mass)*tauscaling(ig)
              enddo
              if (scavenging) then
                do ig=1,ngrid 
                  dqdustsurf(ig) = dqdustsurf(ig) + 
     &                     zdqssed(ig,igcm_ccn_mass)*tauscaling(ig)
                  dndustsurf(ig) = dndustsurf(ig) + 
     &                     zdqssed(ig,igcm_ccn_number)*tauscaling(ig)
                  nccn(ig,:) =
     &                     zq(ig,:,igcm_ccn_number)*tauscaling(ig)
                  qccn(ig,:) =
     &                     zq(ig,:,igcm_ccn_mass)*tauscaling(ig)
                enddo
              endif
           endif ! of (doubleq)

           if (rdstorm) then   ! diagnostics of stormdust tendancies for 1D and 3D
              mstormdtot(:)=0
              mdusttot(:)=0
              qdusttotal(:,:)=0
              do ig=1,ngrid 
                rdsdqdustsurf(ig) = 
     &                zdqssed(ig,igcm_stormdust_mass)*tauscaling(ig)
                rdsdndustsurf(ig) = 
     &                zdqssed(ig,igcm_stormdust_number)*tauscaling(ig)
                rdsndust(ig,:) =
     &                pq(ig,:,igcm_stormdust_number)*tauscaling(ig)
                rdsqdust(ig,:) =
     &                pq(ig,:,igcm_stormdust_mass)*tauscaling(ig)
                do l=1,nlayer
                    mstormdtot(ig) = mstormdtot(ig) + 
     &                      zq(ig,l,igcm_stormdust_mass) * 
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g
                    mdusttot(ig) = mdusttot(ig) + 
     &                      zq(ig,l,igcm_dust_mass) * 
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g
                    qdusttotal(ig,l) = qdust(ig,l)+rdsqdust(ig,l) !calculate total dust
                enddo
              enddo
           endif !(rdstorm)

           if (co2clouds) then
              do ig=1,ngrid 
                 nccnco2(ig,:) =
     &                zq(ig,:,igcm_ccnco2_number)*tauscaling(ig)
                 qccnco2(ig,:) =
     &                zq(ig,:,igcm_ccnco2_mass)*tauscaling(ig)
              enddo
c D. BARDET compute integrated CO2 vapor and ice content
              mtotco2(:)=0
              icetotco2(:)=0
              do ig=1,ngrid 
               do l=1,nlayer
                 mtotco2(ig) = mtotco2(ig) + 
     &                      zq(ig,l,igcm_co2) * 
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g
                 icetotco2(ig) = icetotco2(ig) + 
     &                        zq(ig,l,igcm_co2_ice) * 
     &                        (zplev(ig,l) - zplev(ig,l+1)) / g
               enddo
              enddo
           endif ! of if (co2clouds)
                              
           if (water) then
             mtot(:)=0
             icetot(:)=0
             rave(:)=0
             tauTES(:)=0
             do ig=1,ngrid 
               do l=1,nlayer
                 mtot(ig) = mtot(ig) + 
     &                      zq(ig,l,igcm_h2o_vap) * 
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g
                 icetot(ig) = icetot(ig) + 
     &                        zq(ig,l,igcm_h2o_ice) * 
     &                        (zplev(ig,l) - zplev(ig,l+1)) / g
c                Computing abs optical depth at 825 cm-1 in each
c                  layer to simulate NEW TES retrieval
                 Qabsice = min(
     &             max(0.4e6*rice(ig,l)*(1.+nuice_ref)-0.05 ,0.),1.2
     &                        )
                 opTES(ig,l)= 0.75 * Qabsice * 
     &             zq(ig,l,igcm_h2o_ice) *
     &             (zplev(ig,l) - zplev(ig,l+1)) / g
     &             / (rho_ice * rice(ig,l) * (1.+nuice_ref))
                 tauTES(ig)=tauTES(ig)+ opTES(ig,l) 
               enddo
c              rave(ig)=rave(ig)/max(icetot(ig),1.e-30)       ! mass weight
c               if (icetot(ig)*1e3.lt.0.01) rave(ig)=0.
             enddo

             call watersat(ngrid*nlayer,zt,zplay,zqsat)
             satu(:,:) = zq(:,:,igcm_h2o_vap)/zqsat(:,:)

             if (scavenging) then
               Nccntot(:)= 0
               Mccntot(:)= 0
               rave(:)=0
               do ig=1,ngrid 
                 do l=1,nlayer
                    Nccntot(ig) = Nccntot(ig) + 
     &              zq(ig,l,igcm_ccn_number)*tauscaling(ig)
     &              *(zplev(ig,l) - zplev(ig,l+1)) / g
                    Mccntot(ig) = Mccntot(ig) + 
     &              zq(ig,l,igcm_ccn_mass)*tauscaling(ig)
     &              *(zplev(ig,l) - zplev(ig,l+1)) / g
cccc Column integrated effective ice radius 
cccc is weighted by total ice surface area (BETTER than total ice mass) 
                    rave(ig) = rave(ig) + 
     &                      tauscaling(ig) *
     &                      zq(ig,l,igcm_ccn_number) *
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g * 
     &                      rice(ig,l) * rice(ig,l)*  (1.+nuice_ref)
                 enddo
               rave(ig)=(icetot(ig)/rho_ice+Mccntot(ig)/rho_dust)*0.75
     &                  /max(pi*rave(ig),1.e-30) ! surface weight
               if (icetot(ig)*1e3.lt.0.01) rave(ig)=0.
               enddo
             else ! of if (scavenging)
               rave(:)=0
               do ig=1,ngrid 
                 do l=1,nlayer
                 rave(ig) = rave(ig) + 
     &                      zq(ig,l,igcm_h2o_ice) *
     &                      (zplev(ig,l) - zplev(ig,l+1)) / g * 
     &                      rice(ig,l) * (1.+nuice_ref)
                 enddo
                 rave(ig) = max(rave(ig) / 
     &             max(icetot(ig),1.e-30),1.e-30) ! mass weight
               enddo
             endif ! of if (scavenging)

           !Alternative A. Pottier weighting
           rave2(:) = 0.
           totrave2(:) = 0.
           do ig=1,ngrid
              do l=1,nlayer
              rave2(ig) =rave2(ig)+ zq(ig,l,igcm_h2o_ice)*rice(ig,l)
              totrave2(ig) = totrave2(ig) + zq(ig,l,igcm_h2o_ice)
              end do
              rave2(ig)=max(rave2(ig)/max(totrave2(ig),1.e-30),1.e-30)
           end do

           endif ! of if (water)


        endif                   ! of if (tracer)

#ifndef MESOSCALE
c        -----------------------------------------------------------------
c        WSTATS: Saving statistics
c        -----------------------------------------------------------------
c        ("stats" stores and accumulates 8 key variables in file "stats.nc"
c        which can later be used to make the statistic files of the run:
c        "stats")          only possible in 3D runs !
         
       IF (callstats) THEN

        call wstats(ngrid,"ps","Surface pressure","Pa",2,ps)
        call wstats(ngrid,"tsurf","Surface temperature","K",2,tsurf)
        call wstats(ngrid,"co2ice","CO2 ice cover",
     &                "kg.m-2",2,co2ice)
        call wstats(ngrid,"watercap","H2O ice cover",
     &                "kg.m-2",2,watercap)
        call wstats(ngrid,"tauref","reference dod at 610 Pa","NU",
     &                2,tauref)
        call wstats(ngrid,"fluxsurf_lw",
     &                "Thermal IR radiative flux to surface","W.m-2",2,
     &                fluxsurf_lw)
        call wstats(ngrid,"fluxsurf_sw",
     &                "Solar radiative flux to surface","W.m-2",2,
     &                fluxsurf_sw_tot)
        call wstats(ngrid,"fluxtop_lw",
     &                "Thermal IR radiative flux to space","W.m-2",2,
     &                fluxtop_lw)
        call wstats(ngrid,"fluxtop_sw",
     &                "Solar radiative flux to space","W.m-2",2,
     &                fluxtop_sw_tot)
        call wstats(ngrid,"temp","Atmospheric temperature","K",3,zt)
        call wstats(ngrid,"u","Zonal (East-West) wind","m.s-1",3,zu)
        call wstats(ngrid,"v","Meridional (North-South) wind",
     &                "m.s-1",3,zv)
        call wstats(ngrid,"w","Vertical (down-up) wind",
     &                "m.s-1",3,pw)
        call wstats(ngrid,"rho","Atmospheric density","kg/m3",3,rho)
        call wstats(ngrid,"pressure","Pressure","Pa",3,zplay)
          call wstats(ngrid,"q2",
     &                "Boundary layer eddy kinetic energy",
     &                "m2.s-2",3,q2)
          call wstats(ngrid,"emis","Surface emissivity","w.m-1",2,
     &                emis)
c          call wstats(ngrid,"ssurf","Surface stress","N.m-2",
c    &                2,zstress)
c          call wstats(ngrid,"sw_htrt","sw heat.rate",
c    &                 "W.m-2",3,zdtsw)
c          call wstats(ngrid,"lw_htrt","lw heat.rate",
c    &                 "W.m-2",3,zdtlw)

          if (calltherm) then
            call wstats(ngrid,"zmax_th","Height of thermals",
     &                "m",2,zmax_th)
            call wstats(ngrid,"hfmax_th","Max thermals heat flux",
     &                "K.m/s",2,hfmax_th)
            call wstats(ngrid,"wstar",
     &                "Max vertical velocity in thermals",
     &                "m/s",2,wstar)
          endif

           if (tracer) then
             if (water) then
               vmr=zq(1:ngrid,1:nlayer,igcm_h2o_vap)
     &      *mmean(1:ngrid,1:nlayer)/mmol(igcm_h2o_vap)
               call wstats(ngrid,"vmr_h2ovap",
     &                    "H2O vapor volume mixing ratio","mol/mol",
     &                    3,vmr)
               vmr=zq(1:ngrid,1:nlayer,igcm_h2o_ice)
     &      *mmean(1:ngrid,1:nlayer)/mmol(igcm_h2o_ice)
               call wstats(ngrid,"vmr_h2oice",
     &                    "H2O ice volume mixing ratio","mol/mol",
     &                    3,vmr)
              ! also store vmr_ice*rice for better diagnostics of rice
               vmr(1:ngrid,1:nlayer)=vmr(1:ngrid,1:nlayer)*
     &                               rice(1:ngrid,1:nlayer)
               call wstats(ngrid,"vmr_h2oice_rice",
     &                "H2O ice mixing ratio times ice particule size",
     &                    "(mol/mol)*m",
     &                    3,vmr)
               vmr=zqsat(1:ngrid,1:nlayer)
     &      *mmean(1:ngrid,1:nlayer)/mmol(igcm_h2o_vap)
               call wstats(ngrid,"vmr_h2osat",
     &                    "saturation volume mixing ratio","mol/mol",
     &                    3,vmr)
               call wstats(ngrid,"h2o_ice_s",
     &                    "surface h2o_ice","kg/m2",
     &                    2,qsurf(1,igcm_h2o_ice))
               call wstats(ngrid,'albedo',
     &                         'albedo',
     &                         '',2,albedo(1,1))
               call wstats(ngrid,"mtot",
     &                    "total mass of water vapor","kg/m2",
     &                    2,mtot)
               call wstats(ngrid,"icetot",
     &                    "total mass of water ice","kg/m2",
     &                    2,icetot)
               call wstats(ngrid,"reffice",
     &                    "Mean reff","m",
     &                    2,rave)
              call wstats(ngrid,"Nccntot",
     &                    "condensation nuclei","Nbr/m2",
     &                    2,Nccntot)
              call wstats(ngrid,"Mccntot",
     &                    "condensation nuclei mass","kg/m2",
     &                    2,Mccntot)
              call wstats(ngrid,"rice",
     &                    "Ice particle size","m",
     &                    3,rice)
               if (.not.activice) then
                 call wstats(ngrid,"tauTESap",
     &                      "tau abs 825 cm-1","",
     &                      2,tauTES)
               else
                 call wstats(ngrid,'tauTES',
     &                   'tau abs 825 cm-1',
     &                  '',2,taucloudtes)
               endif

             endif ! of if (water)
             
             
           if (dustbin.ne.0) then
          
             call wstats(ngrid,'tau','taudust','SI',2,tau(1,1))
             
             if (doubleq) then
c            call wstats(ngrid,'qsurf','qsurf',
c     &                       'kg.m-2',2,qsurf(1,igcm_dust_mass))
c            call wstats(ngrid,'Nsurf','N particles',
c     &                       'N.m-2',2,qsurf(1,igcm_dust_number))
c            call wstats(ngrid,'dqsdev','ddevil lift',
c    &                       'kg.m-2.s-1',2,zdqsdev(1,1))
c            call wstats(ngrid,'dqssed','sedimentation',
c     &                       'kg.m-2.s-1',2,zdqssed(1,1))
c            call wstats(ngrid,'dqsdif','diffusion',
c     &                       'kg.m-2.s-1',2,zdqsdif(1,1))
               call wstats(ngrid,'dqsdust',
     &                        'deposited surface dust mass',
     &                        'kg.m-2.s-1',2,dqdustsurf)
               call wstats(ngrid,'dqndust',
     &                        'deposited surface dust number',
     &                        'number.m-2.s-1',2,dndustsurf)
               call wstats(ngrid,'reffdust','reffdust',
     &                        'm',3,rdust*ref_r0)
               call wstats(ngrid,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
               call wstats(ngrid,'dustN','Dust number',
     &                        'part/kg',3,ndust)
               if (rdstorm) then
                 call wstats(ngrid,'reffstormdust','reffdust',
     &                          'm',3,rstormdust*ref_r0)
                 call wstats(ngrid,'rdsdustq','Dust mass mr',
     &                          'kg/kg',3,rdsqdust)
                 call wstats(ngrid,'rdsdustN','Dust number',
     &                        'part/kg',3,rdsndust)
               end if
             else
               do iq=1,dustbin
                 write(str2(1:2),'(i2.2)') iq
                 call wstats(ngrid,'q'//str2,'mix. ratio',
     &                         'kg/kg',3,zq(1,1,iq))
                 call wstats(ngrid,'qsurf'//str2,'qsurf',
     &                         'kg.m-2',2,qsurf(1,iq))
               end do
             endif ! (doubleq)

             if (scavenging) then
               call wstats(ngrid,'ccnq','CCN mass mr',
     &                        'kg/kg',3,qccn)
               call wstats(ngrid,'ccnN','CCN number',
     &                        'part/kg',3,nccn)
             endif ! (scavenging)
          
           endif ! (dustbin.ne.0)

           if (photochem) then
              do iq=1,nq
                 if (noms(iq) .ne. "dust_mass" .and.
     $               noms(iq) .ne. "dust_number" .and.
     $               noms(iq) .ne. "ccn_mass" .and.
     $               noms(iq) .ne. "ccn_number" .and.
     $               noms(iq) .ne. "ccnco2_mass" .and.
     $               noms(iq) .ne. "ccnco2_number") then

!                   volume mixing ratio

                    vmr(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,iq)
     &                            *mmean(1:ngrid,1:nlayer)/mmol(iq)

                    call wstats(ngrid,"vmr_"//trim(noms(iq)),
     $                        "Volume mixing ratio","mol/mol",3,vmr)
                    if ((noms(iq).eq."o")
     $             .or. (noms(iq).eq."co2")
     $             .or. (noms(iq).eq."o3")
     $             .or. (noms(iq).eq."ar")
     $             .or. (noms(iq).eq."o2")
     $             .or. (noms(iq).eq."h2o_vap") ) then
                      call writediagfi(ngrid,"vmr_"//trim(noms(iq)),
     $                         "Volume mixing ratio","mol/mol",3,vmr)
                    end if

!                   number density (molecule.cm-3)

                    rhopart(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,iq)
     &                          *rho(1:ngrid,1:nlayer)*n_avog/
     &                           (1000*mmol(iq))

                   call wstats(ngrid,"rho_"//trim(noms(iq)),
     $                   "Number density","cm-3",3,rhopart)
                   call writediagfi(ngrid,"rho_"//trim(noms(iq)),
     $                  "Number density","cm-3",3,rhopart)

!                   vertical column (molecule.cm-2)

                    do ig = 1,ngrid
                       colden(ig,iq) = 0.                           
                    end do
                    do l=1,nlayer                                    
                       do ig=1,ngrid                                  
                          colden(ig,iq) = colden(ig,iq) + zq(ig,l,iq)
     $                                   *(zplev(ig,l)-zplev(ig,l+1)) 
     $                                   *6.022e22/(mmol(iq)*g)       
                       end do                                        
                    end do                                          

                    call wstats(ngrid,"c_"//trim(noms(iq)),           
     $                          "column","mol cm-2",2,colden(1,iq))  
                    call writediagfi(ngrid,"c_"//trim(noms(iq)),
     $                          "column","mol cm-2",2,colden(1,iq))

!                   global mass (g)
               
                    call planetwide_sumval(colden(:,iq)/6.022e23
     $                            *mmol(iq)*1.e4*cell_area(:),mass(iq))

                    call writediagfi(ngrid,"mass_"//trim(noms(iq)),
     $                              "global mass","g",0,mass(iq))

                 end if ! of if (noms(iq) .ne. "dust_mass" ...)
              end do ! of do iq=1,nq
           end if ! of if (photochem)

           end if ! of if (tracer)

           IF(lastcall) THEN
             write (*,*) "Writing stats..."
             call mkstats(ierr)
           ENDIF

         ENDIF !if callstats

c        (Store EOF for Mars Climate database software)
         IF (calleofdump) THEN
            CALL eofdump(ngrid, nlayer, zu, zv, zt, rho, ps)
         ENDIF
#endif
!endif of ifndef MESOSCALE

#ifdef MESOSCALE
      
      !! see comm_wrf. 
      !! not needed when an array is already in a shared module.
      !! --> example : hfmax_th, zmax_th

      CALL allocate_comm_wrf(ngrid,nlayer)

      !state  real  HR_SW     ikj   misc  1  -  h  "HR_SW"     "HEATING RATE SW"                 "K/s"
      comm_HR_SW(1:ngrid,1:nlayer) = zdtsw(1:ngrid,1:nlayer)
      !state  real  HR_LW     ikj   misc  1  -  h  "HR_LW"     "HEATING RATE LW"                 "K/s"
      comm_HR_LW(1:ngrid,1:nlayer) = zdtlw(1:ngrid,1:nlayer)
      !state  real  SWDOWNZ    ij   misc  1  -  h  "SWDOWNZ"   "DOWNWARD SW FLUX AT SURFACE"     "W m-2"
      comm_SWDOWNZ(1:ngrid) = fluxsurf_sw_tot(1:ngrid)
      !state  real  TAU_DUST   ij   misc  1  -  h  "TAU_DUST"  "REFERENCE VISIBLE DUST OPACITY"  ""
      comm_TAU_DUST(1:ngrid) = tauref(1:ngrid)
      !state  real  RDUST     ikj   misc  1  -  h  "RDUST"     "DUST RADIUS"                     "m"
      comm_RDUST(1:ngrid,1:nlayer) = rdust(1:ngrid,1:nlayer)
      !state  real  QSURFDUST  ij   misc  1  -  h  "QSURFDUST" "DUST MASS AT SURFACE"            "kg m-2"
      IF (igcm_dust_mass .ne. 0) THEN
        comm_QSURFDUST(1:ngrid) = qsurf(1:ngrid,igcm_dust_mass)
      ELSE
        comm_QSURFDUST(1:ngrid) = 0.
      ENDIF
      !state  real  MTOT       ij   misc  1  -  h  "MTOT"      "TOTAL MASS WATER VAPOR in pmic"  "pmic"
      comm_MTOT(1:ngrid) = mtot(1:ngrid) * 1.e6 / rho_ice
      !state  real  ICETOT     ij   misc  1  -  h  "ICETOT"    "TOTAL MASS WATER ICE"            "kg m-2"
      comm_ICETOT(1:ngrid) = icetot(1:ngrid) * 1.e6 / rho_ice
      !state  real  VMR_ICE   ikj   misc  1  -  h  "VMR_ICE"   "VOL. MIXING RATIO ICE"           "ppm"
      IF (igcm_h2o_ice .ne. 0) THEN
        comm_VMR_ICE(1:ngrid,1:nlayer) = 1.e6
     .    * zq(1:ngrid,1:nlayer,igcm_h2o_ice)
     .    * mmean(1:ngrid,1:nlayer) / mmol(igcm_h2o_ice)
      ELSE
        comm_VMR_ICE(1:ngrid,1:nlayer) = 0.
      ENDIF
      !state  real  TAU_ICE    ij   misc  1  -  h  "TAU_ICE"   "CLOUD OD at 825 cm-1 TES"        ""
      if (activice) then
        comm_TAU_ICE(1:ngrid) = taucloudtes(1:ngrid)
      else
        comm_TAU_ICE(1:ngrid) = tauTES(1:ngrid)
      endif
      !state  real  RICE      ikj   misc  1  -  h  "RICE"      "ICE RADIUS"                      "m"
      comm_RICE(1:ngrid,1:nlayer) = rice(1:ngrid,1:nlayer)
   
      !! calculate sensible heat flux in W/m2 for outputs
      !! -- the one computed in vdifc is not the real one
      !! -- vdifc must have been called
      if (.not.callrichsl) then
        sensibFlux(1:ngrid) = zflubid(1:ngrid)
     .         - capcal(1:ngrid)*zdtsdif(1:ngrid)
      else
        sensibFlux(1:ngrid) = 
     &   (pplay(1:ngrid,1)/(r*pt(1:ngrid,1)))*cpp
     &   *sqrt(pu(1:ngrid,1)*pu(1:ngrid,1)+pv(1:ngrid,1)*pv(1:ngrid,1)
     &         +(log(1.+0.7*wstar(1:ngrid) + 2.3*wstar(1:ngrid)**2))**2)
     &   *zcdh(1:ngrid)*(tsurf(1:ngrid)-zh(1:ngrid,1))
      endif
          
#else
#ifndef MESOINI

c        ==========================================================
c        WRITEDIAGFI: Outputs in netcdf file "DIAGFI", containing
c          any variable for diagnostic (output with period
c          "ecritphy", set in "run.def")
c        ==========================================================
c        WRITEDIAGFI can ALSO be called from any other subroutines
c        for any variables !!
c        call WRITEDIAGFI(ngrid,"emis","Surface emissivity","w.m-1",2,
c    &                  emis)
c        call WRITEDIAGFI(ngrid,"pplay","Pressure","Pa",3,zplay)
c        call WRITEDIAGFI(ngrid,"pplev","Pressure","Pa",3,zplev)
         call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",2,
     &                  tsurf)
         call WRITEDIAGFI(ngrid,"ps","surface pressure","Pa",2,ps)
         call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness"
     &                                         ,"kg.m-2",2,co2ice)
         call WRITEDIAGFI(ngrid,"watercap","Water ice thickness"
     &                                         ,"kg.m-2",2,watercap)

         call WRITEDIAGFI(ngrid,"temp7","temperature in layer 7",
     &                  "K",2,zt(1,7))
         call WRITEDIAGFI(ngrid,"fluxsurf_lw","fluxsurf_lw","W.m-2",2,
     &                  fluxsurf_lw)
         call WRITEDIAGFI(ngrid,"fluxsurf_sw","fluxsurf_sw","W.m-2",2,
     &                  fluxsurf_sw_tot)
         call WRITEDIAGFI(ngrid,"fluxtop_lw","fluxtop_lw","W.m-2",2,
     &                  fluxtop_lw)
         call WRITEDIAGFI(ngrid,"fluxtop_sw","fluxtop_sw","W.m-2",2,
     &                  fluxtop_sw_tot)
         call WRITEDIAGFI(ngrid,"temp","temperature","K",3,zt)
         call WRITEDIAGFI(ngrid,"Sols","Time","sols",0,zday)
         call WRITEDIAGFI(ngrid,"Ls","Solar longitude","deg",
     &                    0,zls*180./pi)
         call WRITEDIAGFI(ngrid,"u","Zonal wind","m.s-1",3,zu)
         call WRITEDIAGFI(ngrid,"v","Meridional wind","m.s-1",3,zv)
         call WRITEDIAGFI(ngrid,"w","Vertical wind","m.s-1",3,pw)
         call WRITEDIAGFI(ngrid,"rho","density","kg.m-3",3,rho)
c        call WRITEDIAGFI(ngrid,"q2","q2","kg.m-3",3,q2)
c        call WRITEDIAGFI(ngrid,'Teta','T potentielle','K',3,zh)
         call WRITEDIAGFI(ngrid,"pressure","Pressure","Pa",3,zplay)
c        call WRITEDIAGFI(ngrid,"ssurf","Surface stress","N.m-2",2,
c    &                  zstress)
c        call WRITEDIAGFI(ngrid,'sw_htrt','sw heat. rate',
c    &                   'w.m-2',3,zdtsw)
c        call WRITEDIAGFI(ngrid,'lw_htrt','lw heat. rate',
c    &                   'w.m-2',3,zdtlw)
 
            if (.not.activice) then
               CALL WRITEDIAGFI(ngrid,'tauTESap',
     &                         'tau abs 825 cm-1',
     &                         '',2,tauTES)
             else
               CALL WRITEDIAGFI(ngrid,'tauTES',
     &                         'tau abs 825 cm-1',
     &                         '',2,taucloudtes)
             endif

#else
     !!! this is to ensure correct initialisation of mesoscale model
        call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",2,
     &                  tsurf)
        call WRITEDIAGFI(ngrid,"ps","surface pressure","Pa",2,ps) 
        call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness","kg.m-2",2,
     &                  co2ice)
        call WRITEDIAGFI(ngrid,"temp","temperature","K",3,zt)
        call WRITEDIAGFI(ngrid,"u","Zonal wind","m.s-1",3,zu)
        call WRITEDIAGFI(ngrid,"v","Meridional wind","m.s-1",3,zv)
        call WRITEDIAGFI(ngrid,"emis","Surface emissivity","w.m-1",2,
     &                  emis)
        call WRITEDIAGFI(ngrid,"tsoil","Soil temperature",
     &                       "K",3,tsoil)
        call WRITEDIAGFI(ngrid,"inertiedat","Soil inertia",
     &                       "K",3,inertiedat)
#endif


c        ----------------------------------------------------------
c        Outputs of the CO2 cycle
c        ----------------------------------------------------------
            
         if (tracer.and.(igcm_co2.ne.0)) then
!          call WRITEDIAGFI(ngrid,"co2_l1","co2 mix. ratio in 1st layer",
!    &                     "kg/kg",2,zq(1,1,igcm_co2))
          call WRITEDIAGFI(ngrid,"co2","co2 mass mixing ratio",
     &                     "kg/kg",3,zq(1,1,igcm_co2))
          if (co2clouds) then
            CALL WRITEDIAGFI(ngrid,'mtotco2',
     &                       'total mass of CO2 vapor',
     &                       'kg/m2',2,mtotco2)
            CALL WRITEDIAGFI(ngrid,'zdtcloudco2',
     &                       'temperature variation of CO2 latent heat',
     &                       'K/s',3,zdtcloudco2)

            CALL WRITEDIAGFI(ngrid,'icetotco2',
     &                       'total mass of CO2 ice',
     &                       'kg/m2',2,icetotco2)
 
            call WRITEDIAGFI(ngrid,'ccnqco2','CCNco2 mass mr',
     &                       'kg/kg',3,qccnco2)
            call WRITEDIAGFI(ngrid,'ccnNco2','CCNco2 number',
     &                       'part/kg',3,nccnco2)
            call WRITEDIAGFI(ngrid,'co2_ice','co2_ice','kg/kg',
     &                       3,zq(:,:,igcm_co2_ice))
            call WRITEDIAGFI(ngrid,'precip_co2_ice',
     &                       'surface deposition of co2 ice',
     &                       'kg.m-2.s-1',2,
     &                       zdqssed(1:ngrid,igcm_co2_ice))
           endif ! of if (co2clouds)
         endif ! of if (tracer.and.(igcm_co2.ne.0))
        ! Output He tracer, if there is one
        if (tracer.and.(igcm_he.ne.0)) then
          call WRITEDIAGFI(ngrid,"he","helium mass mixing ratio",
     &                     "kg/kg",3,zq(1,1,igcm_he))
          vmr=zq(1:ngrid,1:nlayer,igcm_he)
     &       *mmean(1:ngrid,1:nlayer)/mmol(igcm_he)
          call WRITEDIAGFI(ngrid,'vmr_he','helium vol. mixing ratio',
     &                       'mol/mol',3,vmr)
        endif

c        ----------------------------------------------------------
c        Outputs of the water cycle
c        ----------------------------------------------------------
         if (tracer) then
           if (water) then

#ifdef MESOINI
            !!!! waterice = q01, voir readmeteo.F90
            call WRITEDIAGFI(ngrid,'q01',noms(igcm_h2o_ice),
     &                      'kg/kg',3,
     &                       zq(1:ngrid,1:nlayer,igcm_h2o_ice))
            !!!! watervapor = q02, voir readmeteo.F90
            call WRITEDIAGFI(ngrid,'q02',noms(igcm_h2o_vap),
     &                      'kg/kg',3,
     &                       zq(1:ngrid,1:nlayer,igcm_h2o_vap))
            !!!! surface waterice qsurf02 (voir readmeteo)
            call WRITEDIAGFI(ngrid,'qsurf02','surface tracer',
     &                      'kg.m-2',2,
     &                       qsurf(1:ngrid,igcm_h2o_ice))
#endif
            CALL WRITEDIAGFI(ngrid,'mtot',
     &                       'total mass of water vapor',
     &                       'kg/m2',2,mtot)
            CALL WRITEDIAGFI(ngrid,'icetot',
     &                       'total mass of water ice',
     &                       'kg/m2',2,icetot)
            vmr=zq(1:ngrid,1:nlayer,igcm_h2o_ice)
     &      *mmean(1:ngrid,1:nlayer)/mmol(igcm_h2o_ice)
            CALL WRITEDIAGFI(ngrid,'vmr_h2oice','h2o ice vmr',
     &                       'mol/mol',3,vmr)
            vmr=zq(1:ngrid,1:nlayer,igcm_h2o_vap)
     &      *mmean(1:ngrid,1:nlayer)/mmol(igcm_h2o_vap)
            CALL WRITEDIAGFI(ngrid,'vmr_h2ovap','h2o vap vmr',
     &                       'mol/mol',3,vmr)
            CALL WRITEDIAGFI(ngrid,'reffice',
     &                       'Mean reff',
     &                       'm',2,rave)

            call WRITEDIAGFI(ngrid,'h2o_ice','h2o_ice','kg/kg',
     &             3,zq(:,:,igcm_h2o_ice))
            call WRITEDIAGFI(ngrid,'h2o_vap','h2o_vap','kg/kg',
     &             3,zq(:,:,igcm_h2o_vap))


!A. Pottier
             CALL WRITEDIAGFI(ngrid,'rmoym',
     &                      'alternative reffice',
     &                      'm',2,rave2)
            call WRITEDIAGFI(ngrid,'saturation',
     &           'h2o vap saturation ratio','dimless',3,satu)
            if (scavenging) then
              CALL WRITEDIAGFI(ngrid,"Nccntot",
     &                    "condensation nuclei","Nbr/m2",
     &                    2,Nccntot)
              CALL WRITEDIAGFI(ngrid,"Mccntot",
     &                    "mass condensation nuclei","kg/m2",
     &                    2,Mccntot)
            endif
            call WRITEDIAGFI(ngrid,'rice','Ice particle size',
     &                       'm',3,rice)
            call WRITEDIAGFI(ngrid,'h2o_ice_s',
     &                       'surface h2o_ice',
     &                       'kg.m-2',2,qsurf(1,igcm_h2o_ice))
            CALL WRITEDIAGFI(ngrid,'albedo',
     &                         'albedo',
     &                         '',2,albedo(1,1))
              if (tifeedback) then
                 call WRITEDIAGSOIL(ngrid,"soiltemp",
     &                              "Soil temperature","K",
     &                              3,tsoil)
                 call WRITEDIAGSOIL(ngrid,'soilti',
     &                       'Soil Thermal Inertia',
     &                       'J.s-1/2.m-2.K-1',3,inertiesoil)
              endif
!A. Pottier
              if (CLFvarying) then !AP14 nebulosity
                 call WRITEDIAGFI(ngrid,'totcloudfrac',
     &                'Total cloud fraction',
     &                       ' ',2,totcloudfrac)
              endif !clf varying

           endif !(water)


           if (water.and..not.photochem) then
             iq=nq
c            write(str2(1:2),'(i2.2)') iq
c            call WRITEDIAGFI(ngrid,'dqs'//str2,'dqscloud',
c    &                       'kg.m-2',2,zdqscloud(1,iq))
c            call WRITEDIAGFI(ngrid,'dqch'//str2,'var chim',
c    &                       'kg/kg',3,zdqchim(1,1,iq))
c            call WRITEDIAGFI(ngrid,'dqd'//str2,'var dif',
c    &                       'kg/kg',3,zdqdif(1,1,iq))
c            call WRITEDIAGFI(ngrid,'dqa'//str2,'var adj',
c    &                       'kg/kg',3,zdqadj(1,1,iq))
c            call WRITEDIAGFI(ngrid,'dqc'//str2,'var c',
c    &                       'kg/kg',3,zdqc(1,1,iq))
           endif  !(water.and..not.photochem)
         endif

c        ----------------------------------------------------------
c        Outputs of the dust cycle
c        ----------------------------------------------------------

        call WRITEDIAGFI(ngrid,'tauref',
     &                    'Dust ref opt depth','NU',2,tauref)

         if (tracer.and.(dustbin.ne.0)) then

          call WRITEDIAGFI(ngrid,'tau','taudust','SI',2,tau(1,1))

#ifndef MESOINI
           if (doubleq) then
c            call WRITEDIAGFI(ngrid,'qsurf','qsurf',
c     &                       'kg.m-2',2,qsurf(1,igcm_dust_mass))
c            call WRITEDIAGFI(ngrid,'Nsurf','N particles',
c     &                       'N.m-2',2,qsurf(1,igcm_dust_number))
c            call WRITEDIAGFI(ngrid,'dqsdev','ddevil lift',
c    &                       'kg.m-2.s-1',2,zdqsdev(1,1))
c            call WRITEDIAGFI(ngrid,'dqssed','sedimentation',
c     &                       'kg.m-2.s-1',2,zdqssed(1,1))
c            call WRITEDIAGFI(ngrid,'dqsdif','diffusion',
c     &                       'kg.m-2.s-1',2,zdqsdif(1,1))
c             call WRITEDIAGFI(ngrid,'sedice','sedimented ice',
c     &                       'kg.m-2.s-1',2,zdqssed(:,igcm_h2o_ice))
c             call WRITEDIAGFI(ngrid,'subice','sublimated ice',
c     &                       'kg.m-2.s-1',2,zdqsdif(:,igcm_h2o_ice))
             call WRITEDIAGFI(ngrid,'dqsdust',
     &                        'deposited surface dust mass',
     &                        'kg.m-2.s-1',2,dqdustsurf)
             call WRITEDIAGFI(ngrid,'dqndust',
     &                        'deposited surface dust number',
     &                        'number.m-2.s-1',2,dndustsurf)
             call WRITEDIAGFI(ngrid,'reffdust','reffdust',
     &                        'm',3,rdust*ref_r0)
             call WRITEDIAGFI(ngrid,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
             call WRITEDIAGFI(ngrid,'dustN','Dust number',
     &                        'part/kg',3,ndust)
             
	     select case (trim(dustiropacity))
              case ("tes")
               call WRITEDIAGFI(ngrid,'dsodust',
     &  'density scaled extinction opacity of std dust at 9.3um(TES)',
     &         'm2.kg-1',3,dsodust)
               call WRITEDIAGFI(ngrid,'dso',
     &  'density scaled extinction opacity of all dust at 9.3um(TES)',
     &         'm2.kg-1',3,dsodust+dsords+dsotop)
              case ("mcs")
               call WRITEDIAGFI(ngrid,'dsodust',
     &  'density scaled extinction opacity of std dust at 21.6um(MCS)',
     &         'm2.kg-1',3,dsodust)
               call WRITEDIAGFI(ngrid,'dso',
     &  'density scaled extinction opacity of all dust at 21.6um(MCS)',
     &         'm2.kg-1',3,dsodust+dsords+dsotop)
             end select
           else ! (doubleq=.false.)
             do iq=1,dustbin
               write(str2(1:2),'(i2.2)') iq
               call WRITEDIAGFI(ngrid,'q'//str2,'mix. ratio',
     &                         'kg/kg',3,zq(1,1,iq))
               call WRITEDIAGFI(ngrid,'qsurf'//str2,'qsurf',
     &                         'kg.m-2',2,qsurf(1,iq))
             end do
           endif ! (doubleq)

           if (rdstorm) then  ! writediagfi tendencies stormdust tracers
             call WRITEDIAGFI(ngrid,'reffdust','reffdust',
     &                        'm',3,rdust*ref_r0)
             call WRITEDIAGFI(ngrid,'reffstormdust','reffstormdust',
     &                        'm',3,rstormdust*ref_r0)
             call WRITEDIAGFI(ngrid,'mstormdtot',
     &                        'total mass of stormdust only',
     &                        'kg.m-2',2,mstormdtot)
             call WRITEDIAGFI(ngrid,'mdusttot',
     &                        'total mass of dust only',
     &                        'kg.m-2',2,mdusttot)
             call WRITEDIAGFI(ngrid,'rdsdqsdust',
     &                        'deposited surface stormdust mass',
     &                        'kg.m-2.s-1',2,rdsdqdustsurf)
             call WRITEDIAGFI(ngrid,'rdsdustq','storm Dust mass mr',
     &                        'kg/kg',3,rdsqdust)
             call WRITEDIAGFI(ngrid,'rdsdustqmodel','storm Dust massmr',
     &                        'kg/kg',3,pq(:,:,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'rdsdustN','storm Dust number',
     &                        'part/kg',3,rdsndust)
             call WRITEDIAGFI(ngrid,"stormfract",
     &                "fraction of the mesh, with stormdust","none",
     &                                               2,totstormfract)
             call WRITEDIAGFI(ngrid,'qsurf',
     &                 'stormdust injection',
     &                 'kg.m-2',2,qsurf(:,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'pdqsurf',
     &                  'tendancy stormdust mass at surface',
     &                  'kg.m-2',2,dqsurf(:,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'wspeed','vertical speed stormdust',
     &                        'm/s',3,wspeed(:,1:nlayer))
             call WRITEDIAGFI(ngrid,'zdqsed_dustq'
     &          ,'sedimentation q','kg.m-2.s-1',3,
     &          zdqsed(:,:,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'zdqssed_dustq'
     &          ,'sedimentation q','kg.m-2.s-1',2,
     &          zdqssed(:,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'zdqsed_rdsq'
     &          ,'sedimentation q','kg.m-2.s-1',3,
     &          zdqsed(:,:,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'rdust','rdust',
     &                        'm',3,rdust)
             call WRITEDIAGFI(ngrid,'rstormdust','rstormdust',
     &                        'm',3,rstormdust)
             call WRITEDIAGFI(ngrid,'totaldustq','total dust mass',
     &                        'kg/kg',3,qdusttotal)
             
	     select case (trim(dustiropacity))
              case ("tes")
               call WRITEDIAGFI(ngrid,'dsords',
     &  'density scaled extinction opacity of stormdust at 9.3um(TES)',
     &         'm2.kg-1',3,dsords)
              case ("mcs")
               call WRITEDIAGFI(ngrid,'dsords',
     &  'density scaled extinction opacity of stormdust at 21.6um(MCS)',
     &         'm2.kg-1',3,dsords)
             end select
           endif ! (rdstorm)

           if (slpwind) then
	     select case (trim(dustiropacity))
              case ("tes")
               call WRITEDIAGFI(ngrid,'dsotop',
     &  'density scaled extinction opacity of topdust at 9.3um(TES)',
     &         'm2.kg-1',3,dsotop)
              case ("mcs")
               call WRITEDIAGFI(ngrid,'dsotop',
     &  'density scaled extinction opacity of topdust at 21.6um(MCS)',
     &         'm2.kg-1',3,dsotop)
             end select
           endif ! (slpwind)
	   
           if (scavenging) then
             call WRITEDIAGFI(ngrid,'ccnq','CCN mass mr',
     &                        'kg/kg',3,qccn)
             call WRITEDIAGFI(ngrid,'ccnN','CCN number',
     &                        'part/kg',3,nccn)
             call WRITEDIAGFI(ngrid,'surfccnq','Surf nuclei mass mr',
     &                        'kg.m-2',2,qsurf(1,igcm_ccn_mass))
             call WRITEDIAGFI(ngrid,'surfccnN','Surf nuclei number',
     &                        'kg.m-2',2,qsurf(1,igcm_ccn_number))
           endif ! (scavenging)

c          if (submicron) then
c            call WRITEDIAGFI(ngrid,'dustsubm','subm mass mr',
c    &                        'kg/kg',3,pq(1,1,igcm_dust_submicron))
c          endif ! (submicron)

#else
!     !!! to initialize mesoscale we need scaled variables
!     !!! because this must correspond to starting point for tracers
!             call WRITEDIAGFI(ngrid,'dustq','Dust mass mr',
!     &           'kg/kg',3,pq(1:ngrid,1:nlayer,igcm_dust_mass))
!             call WRITEDIAGFI(ngrid,'dustN','Dust number',
!     &           'part/kg',3,pq(1:ngrid,1:nlayer,igcm_dust_number))
!             call WRITEDIAGFI(ngrid,'ccn','Nuclei mass mr',
!     &           'kg/kg',3,pq(1:ngrid,1:nlayer,igcm_ccn_mass))
!             call WRITEDIAGFI(ngrid,'ccnN','Nuclei number',
!     &           'part/kg',3,pq(1:ngrid,1:nlayer,igcm_ccn_number))
      if (freedust) then
             call WRITEDIAGFI(ngrid,'dustq','Dust mass mr',
     &                        'kg/kg',3,qdust)
             call WRITEDIAGFI(ngrid,'dustN','Dust number',
     &                        'part/kg',3,ndust)
             call WRITEDIAGFI(ngrid,'ccn','CCN mass mr',
     &                        'kg/kg',3,qccn)
             call WRITEDIAGFI(ngrid,'ccnN','CCN number',
     &                        'part/kg',3,nccn)
      else
             call WRITEDIAGFI(ngrid,'dustq','Dust mass mr',
     &                        'kg/kg',3,pq(1,1,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'dustN','Dust number',
     &                        'part/kg',3,pq(1,1,igcm_dust_number))
             call WRITEDIAGFI(ngrid,'ccn','Nuclei mass mr',
     &                        'kg/kg',3,pq(1,1,igcm_ccn_mass))
             call WRITEDIAGFI(ngrid,'ccnN','Nuclei number',
     &                        'part/kg',3,pq(1,1,igcm_ccn_number))
      endif
#endif

         end if  ! (tracer.and.(dustbin.ne.0))

c        ----------------------------------------------------------
c        GW non-oro outputs
c        ----------------------------------------------------------

         if(calllott_nonoro) then 
           call WRITEDIAGFI(ngrid,"dugwno","GW non-oro dU","m/s2",
     $           3,d_u_hin/ptimestep)
           call WRITEDIAGFI(ngrid,"dvgwno","GW non-oro dV","m/s2",
     $           3,d_v_hin/ptimestep)
         endif                  !(calllott_nonoro)

c        ----------------------------------------------------------
c        Thermospheric outputs
c        ----------------------------------------------------------

         if(callthermos) then 

            call WRITEDIAGFI(ngrid,"q15um","15 um cooling","K/s",
     $           3,zdtnlte)
            call WRITEDIAGFI(ngrid,"quv","UV heating","K/s",
     $           3,zdteuv)
            call WRITEDIAGFI(ngrid,"cond","Thermal conduction","K/s",
     $           3,zdtconduc)
            call WRITEDIAGFI(ngrid,"qnir","NIR heating","K/s",
     $           3,zdtnirco2)

!            call wstats(ngrid,"q15um","15 um cooling","K/s",
!     $           3,zdtnlte)
!            call wstats(ngrid,"quv","UV heating","K/s",
!     $           3,zdteuv)
!            call wstats(ngrid,"cond","Thermal conduction","K/s",
!     $           3,zdtconduc)
!            call wstats(ngrid,"qnir","NIR heating","K/s",
!     $           3,zdtnirco2)

         endif  !(callthermos)

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        PBL OUTPUS
c        ----------------------------------------------------------
c        ----------------------------------------------------------

c        ----------------------------------------------------------
c        Outputs of thermals
c        ----------------------------------------------------------
         if (calltherm) then

!        call WRITEDIAGFI(ngrid,'dtke',
!     &              'tendance tke thermiques','m**2/s**2',
!     &                         3,dtke_th)
!        call WRITEDIAGFI(ngrid,'d_u_ajs',
!     &              'tendance u thermiques','m/s',
!     &                         3,pdu_th*ptimestep)
!        call WRITEDIAGFI(ngrid,'d_v_ajs',
!     &              'tendance v thermiques','m/s',
!     &                         3,pdv_th*ptimestep)
!        if (tracer) then
!        if (nq .eq. 2) then
!        call WRITEDIAGFI(ngrid,'deltaq_th',
!     &              'delta q thermiques','kg/kg',
!     &                         3,ptimestep*pdq_th(:,:,2))
!        endif
!        endif

        call WRITEDIAGFI(ngrid,'zmax_th',
     &              'hauteur du thermique','m',
     &                         2,zmax_th)
        call WRITEDIAGFI(ngrid,'hfmax_th',
     &              'maximum TH heat flux','K.m/s',
     &                         2,hfmax_th)
        call WRITEDIAGFI(ngrid,'wstar',
     &              'maximum TH vertical velocity','m/s',
     &                         2,wstar)

         endif

c        ----------------------------------------------------------
c        ----------------------------------------------------------
c        END OF PBL OUTPUS
c        ----------------------------------------------------------
c        ----------------------------------------------------------


c        ----------------------------------------------------------
c        Output in netcdf file "diagsoil.nc" for subterranean
c          variables (output every "ecritphy", as for writediagfi)
c        ----------------------------------------------------------

         ! Write soil temperature
!        call writediagsoil(ngrid,"soiltemp","Soil temperature","K",
!    &                     3,tsoil)
         ! Write surface temperature
!        call writediagsoil(ngrid,"tsurf","Surface temperature","K",
!    &                     2,tsurf)

c        ==========================================================
c        END OF WRITEDIAGFI
c        ==========================================================
#endif
! of ifdef MESOSCALE

      ELSE     ! if(ngrid.eq.1)

#ifndef MESOSCALE
         write(*,'("Ls =",f11.6," tauref(",f4.0," Pa) =",f9.6)')
     &    zls*180./pi,odpref,tauref
c      ----------------------------------------------------------------------
c      Output in grads file "g1d" (ONLY when using testphys1d)
c      (output at every X physical timestep)
c      ----------------------------------------------------------------------
c
c        CALL writeg1d(ngrid,1,fluxsurf_lw,'Fs_ir','W.m-2')
c         CALL writeg1d(ngrid,1,tsurf,'tsurf','K')
c         CALL writeg1d(ngrid,1,ps,'ps','Pa')
        
c         CALL writeg1d(ngrid,nlayer,zt,'T','K')
c        CALL writeg1d(ngrid,nlayer,pu,'u','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pv,'v','m.s-1')
c        CALL writeg1d(ngrid,nlayer,pw,'w','m.s-1')

! THERMALS STUFF 1D
         if(calltherm) then

        call WRITEDIAGFI(ngrid,'lmax_th',
     &              'hauteur du thermique','point',
     &                         0,lmax_th_out)
        call WRITEDIAGFI(ngrid,'zmax_th',
     &              'hauteur du thermique','m',
     &                         0,zmax_th)
        call WRITEDIAGFI(ngrid,'hfmax_th',
     &              'maximum TH heat flux','K.m/s',
     &                         0,hfmax_th)
        call WRITEDIAGFI(ngrid,'wstar',
     &              'maximum TH vertical velocity','m/s',
     &                         0,wstar)

         co2col(:)=0.
         if (tracer) then
         do l=1,nlayer
           do ig=1,ngrid
             co2col(ig)=co2col(ig)+
     &                  zq(ig,l,1)*(zplev(ig,l)-zplev(ig,l+1))/g
         enddo
         enddo

         end if
         call WRITEDIAGFI(ngrid,'co2col','integrated co2 mass'          &
     &                                      ,'kg/m-2',0,co2col)
         endif ! of if (calltherm)

         call WRITEDIAGFI(ngrid,'w','vertical velocity'                 &
     &                              ,'m/s',1,pw)
         call WRITEDIAGFI(ngrid,"q2","q2","kg.m-3",1,q2)
         call WRITEDIAGFI(ngrid,"tsurf","Surface temperature","K",0,
     &                  tsurf)
         call WRITEDIAGFI(ngrid,"u","u wind","m/s",1,zu)
         call WRITEDIAGFI(ngrid,"v","v wind","m/s",1,zv)

         call WRITEDIAGFI(ngrid,"pplay","Pressure","Pa",1,zplay)
         call WRITEDIAGFI(ngrid,"pplev","Pressure","Pa",1,zplev)
         call WRITEDIAGFI(ngrid,"rho","rho","kg.m-3",1,rho)
         call WRITEDIAGFI(ngrid,"dtrad","rad. heat. rate",              &
     &              "K.s-1",1,dtrad)
        call WRITEDIAGFI(ngrid,'sw_htrt','sw heat. rate',
     &                   'w.m-2',1,zdtsw)
        call WRITEDIAGFI(ngrid,'lw_htrt','lw heat. rate',
     &                   'w.m-2',1,zdtlw)
         call WRITEDIAGFI(ngrid,"co2ice","co2 ice thickness"
     &                                   ,"kg.m-2",0,co2ice)

         if (igcm_co2.ne.0) then
          call co2sat(ngrid*nlayer,zt,zplay,zqsatco2)
           do ig=1,ngrid 
            do l=1,nlayer
               satuco2(ig,l) = zq(ig,l,igcm_co2)* 
     &              (mmean(ig,l)/44.01)*zplay(ig,l)/zqsatco2(ig,l)
                  
c               write(*,*) "In PHYSIQMOD, pt,zt,time ",pt(ig,l)
c     &              ,zt(ig,l),ptime
            enddo
           enddo
          endif

c         CALL writeg1d(ngrid,nlayer,zt,'temp','K')
c         CALL writeg1d(ngrid,nlayer,riceco2,'riceco2','m')
c         CALL writeg1d(ngrid,nlayer,satuco2,'satuco2','satu')
         
         
c         call WRITEDIAGFI(ngrid,"satuco2","vap in satu","kg/kg",1,
c     &        satuco2)
c         call WRITEdiagfi(ngrid,"riceco2","ice radius","m"
c     &        ,1,riceco2)
! or output in diagfi.nc (for testphys1d)
         call WRITEDIAGFI(ngrid,'ps','Surface pressure','Pa',0,ps)
         call WRITEDIAGFI(ngrid,'temp','Temperature ',
     &                       'K JA',1,zt)
c         call WRITEDIAGFI(ngrid,'temp2','Temperature ',
c     &        'K JA2',1,pt)

         if(tracer) then
c           CALL writeg1d(ngrid,1,tau,'tau','SI')
            do iq=1,nq
c              CALL writeg1d(ngrid,nlayer,zq(1,1,iq),noms(iq),'kg/kg') 
               call WRITEDIAGFI(ngrid,trim(noms(iq)),
     &              trim(noms(iq)),'kg/kg',1,zq(1,1,iq))
            end do
           if (doubleq) then
             call WRITEDIAGFI(ngrid,'rdust','rdust',
     &                        'm',1,rdust)
           endif ! doubleq 1D
           if (rdstorm) then
             call writediagfi(1,'aerosol_dust','opacity of env. dust',''
     &                        ,1,aerosol(:,:,iaer_dust_doubleq))
             call writediagfi(1,'aerosol_stormdust',
     &                         'opacity of storm dust',''
     &                        ,1,aerosol(:,:,iaer_stormdust_doubleq))
             call WRITEDIAGFI(ngrid,'dqsdifdustq','diffusion',
     &                       'kg.m-2.s-1',0,zdqsdif(1,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'dqsdifrdsq','diffusion',
     &                  'kg.m-2.s-1',0,zdqsdif(1,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'mstormdtot',
     &                        'total mass of stormdust only',
     &                        'kg.m-2',0,mstormdtot)
             call WRITEDIAGFI(ngrid,'mdusttot',
     &                        'total mass of dust only',
     &                        'kg.m-2',0,mdusttot)
             call WRITEDIAGFI(ngrid,'tauref',
     &                    'Dust ref opt depth','NU',0,tauref)
             call WRITEDIAGFI(ngrid,'rdsdqsdust',
     &                        'deposited surface stormdust mass',
     &                        'kg.m-2.s-1',0,rdsdqdustsurf)
             call WRITEDIAGFI(ngrid,'rdsdustq','storm Dust mass mr',
     &                        'kg/kg',1,rdsqdust)
             call WRITEDIAGFI(ngrid,"stormfract", 
     &                         "fraction of the mesh,with stormdust",
     &                         "none",0,totstormfract)
             call WRITEDIAGFI(ngrid,'rdsqsurf',
     &                 'stormdust at surface',
     &                 'kg.m-2',0,qsurf(:,igcm_stormdust_mass))
             call WRITEDIAGFI(ngrid,'qsurf',
     &                  'dust mass at surface',
     &                  'kg.m-2',0,qsurf(:,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'wspeed','vertical speed stormdust',
     &                        'm/s',1,wspeed)
             call WRITEDIAGFI(ngrid,'totaldustq','total dust mass',
     &                        'kg/kg',1,qdusttotal)
             call WRITEDIAGFI(ngrid,'dsords',
     &                       'density scaled opacity of stormdust',
     &                       'm2.kg-1',1,dsords)
             call WRITEDIAGFI(ngrid,'zdqsed_dustq'
     &          ,'sedimentation q','kg.m-2.s-1',1,
     &          zdqsed(1,:,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'zdqsed_rdsq'
     &          ,'sedimentation q','kg.m-2.s-1',1,
     &          zdqsed(1,:,igcm_stormdust_mass))
           endif !(rdstorm 1D)

           if (water.AND.tifeedback) then
             call WRITEDIAGFI(ngrid,"soiltemp",
     &                              "Soil temperature","K",
     &                              1,tsoil)
             call WRITEDIAGFI(ngrid,'soilti',
     &                       'Soil Thermal Inertia',
     &                       'J.s-1/2.m-2.K-1',1,inertiesoil)
           endif
         end if
         
cccccccccccccccccc scavenging & water outputs 1D TN ccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF (water) THEN

           if (.not.activice) then

             tauTES=0
             do l=1,nlayer
               Qabsice = min(
     &             max(0.4e6*rice(1,l)*(1.+nuice_ref)-0.05 ,0.),1.2
     &                      )
               opTES(1,l)= 0.75 * Qabsice *
     &             zq(1,l,igcm_h2o_ice) *
     &             (zplev(1,l) - zplev(1,l+1)) / g
     &             / (rho_ice * rice(1,l) * (1.+nuice_ref))
               tauTES=tauTES+ opTES(1,l)
             enddo
             CALL WRITEDIAGFI(ngrid,'tauTESap',
     &                         'tau abs 825 cm-1',
     &                         '',0,tauTES)
           else

             CALL WRITEDIAGFI(ngrid,'tauTES',
     &                         'tau abs 825 cm-1',
     &                         '',0,taucloudtes)
           endif
     
           mtot = 0
           icetot = 0
           h2otot = qsurf(1,igcm_h2o_ice)

           do l=1,nlayer
             mtot = mtot +  zq(1,l,igcm_h2o_vap) 
     &                 * (zplev(1,l) - zplev(1,l+1)) / g
             icetot = icetot +  zq(1,l,igcm_h2o_ice) 
     &                 * (zplev(1,l) - zplev(1,l+1)) / g
           end do
           h2otot = h2otot+mtot+icetot

             CALL WRITEDIAGFI(ngrid,'h2otot',
     &                         'h2otot',
     &                         'kg/m2',0,h2otot)
             CALL WRITEDIAGFI(ngrid,'mtot',
     &                         'mtot',
     &                         'kg/m2',0,mtot)
             CALL WRITEDIAGFI(ngrid,'icetot',
     &                         'icetot',
     &                         'kg/m2',0,icetot)

           if (scavenging) then

             rave = 0
             do l=1,nlayer
cccc Column integrated effective ice radius
cccc is weighted by total ice surface area (BETTER) 
             rave = rave + tauscaling(1) *
     &              zq(1,l,igcm_ccn_number) *
     &              (zplev(1,l) - zplev(1,l+1)) / g * 
     &              rice(1,l) * rice(1,l)*  (1.+nuice_ref)
             enddo
             rave=icetot*0.75/max(rave*pi*rho_ice,1.e-30) ! surface weight

              Nccntot= 0
              call watersat(ngrid*nlayer,zt,zplay,zqsat)
               do l=1,nlayer
                   Nccntot = Nccntot + 
     &              zq(1,l,igcm_ccn_number)*tauscaling(1)
     &              *(zplev(1,l) - zplev(1,l+1)) / g
                   satu(1,l) = zq(1,l,igcm_h2o_vap)/zqsat(1,l)
                   satu(1,l) = (max(satu(1,l),float(1))-1)
!     &                      * zq(1,l,igcm_h2o_vap) *
!     &                      (zplev(1,l) - zplev(1,l+1)) / g
               enddo
               call WRITEDIAGFI(ngrid,"satu","vap in satu","kg/kg",1,
     &                    satu)
               CALL WRITEDIAGFI(ngrid,'Nccntot',
     &                         'Nccntot',
     &                         'nbr/m2',0,Nccntot)

             call WRITEDIAGFI(ngrid,'zdqsed_dustq'
     & ,'sedimentation q','kg.m-2.s-1',1,zdqsed(1,:,igcm_dust_mass))
             call WRITEDIAGFI(ngrid,'zdqsed_dustN'
     &,'sedimentation N','Nbr.m-2.s-1',1,
     &                                 zdqsed(1,:,igcm_dust_number))

           else ! of if (scavenging)

cccc Column integrated effective ice radius 
cccc is weighted by total ice mass         (LESS GOOD)
             rave = 0
             do l=1,nlayer
               rave = rave + zq(1,l,igcm_h2o_ice)
     &              * (zplev(1,l) - zplev(1,l+1)) / g
     &              *  rice(1,l) * (1.+nuice_ref)
             enddo
             rave=max(rave/max(icetot,1.e-30),1.e-30) ! mass weight
           endif ! of if (scavenging)


           CALL WRITEDIAGFI(ngrid,'reffice',
     &                      'reffice',
     &                      'm',0,rave)

          !Alternative A. Pottier weighting
           rave2 = 0.
           totrave2 = 0.
           do l=1,nlayer
              rave2 =rave2+ zq(1,l,igcm_h2o_ice)*rice(1,l)
              totrave2 = totrave2 + zq(1,l,igcm_h2o_ice)
           end do
           rave2=max(rave2/max(totrave2,1.e-30),1.e-30)
          CALL WRITEDIAGFI(ngrid,'rmoym',
     &                      'reffice',
     &                      'm',0,rave2)

           do iq=1,nq
             call WRITEDIAGFI(ngrid,trim(noms(iq))//'_s',
     &            trim(noms(iq))//'_s','kg/kg',0,qsurf(1,iq))
           end do
     
        call WRITEDIAGFI(ngrid,'zdqcloud_ice','cloud ice',
     &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_ice))
        call WRITEDIAGFI(ngrid,'zdqcloud_vap','cloud vap',
     &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_vap))
        call WRITEDIAGFI(ngrid,'zdqcloud','cloud ice',
     &            'kg.m-2.s-1',1,zdqcloud(1,:,igcm_h2o_ice)
     &                          +zdqcloud(1,:,igcm_h2o_vap))
     
        call WRITEDIAGFI(ngrid,"rice","ice radius","m",1,
     &                    rice)
              
              if (CLFvarying) then
                 call WRITEDIAGFI(ngrid,'totcloudfrac',
     &                'Total cloud fraction',
     &                       ' ',0,totcloudfrac)
              endif !clfvarying

         ENDIF ! of IF (water)
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         zlocal(1)=-log(zplay(1,1)/zplev(1,1))* Rnew(1,1)*zt(1,1)/g

         do l=2,nlayer-1
            tmean=zt(1,l)
            if(zt(1,l).ne.zt(1,l-1))
     &        tmean=(zt(1,l)-zt(1,l-1))/log(zt(1,l)/zt(1,l-1))
              zlocal(l)= zlocal(l-1)
     &        -log(zplay(1,l)/zplay(1,l-1))*rnew(1,l)*tmean/g
         enddo
         zlocal(nlayer)= zlocal(nlayer-1)-
     &                   log(zplay(1,nlayer)/zplay(1,nlayer-1))*
     &                   rnew(1,nlayer)*tmean/g
#endif

      END IF       ! if(ngrid.ne.1)

! XIOS outputs
#ifdef CPP_XIOS      
      ! Send fields to XIOS: (NB these fields must also be defined as
      ! <field id="..." /> in context_lmdz_physics.xml to be correctly used)
      CALL send_xios_field("ls",zls*180./pi)
      
      CALL send_xios_field("ps",ps)
      CALL send_xios_field("area",cell_area)
      
      CALL send_xios_field("temperature",zt)
      CALL send_xios_field("u",zu)
      CALL send_xios_field("v",zv)

      if (lastcall.and.is_omp_master) then
        write(*,*) "physiq lastcall: call xios_context_finalize"
        call xios_context_finalize
      endif
#endif


      icount=icount+1

      END SUBROUTINE physiq

      END MODULE physiq_mod
