!
! $Id: pbl_surface_mod.F90 1413 2010-07-13 07:55:22Z idelkadi $
!
MODULE pbl_surface_mod
!
! Planetary Boundary Layer and Surface module
!
! This module manage the calculation of turbulent diffusion in the boundary layer 
! and all interactions towards the differents sub-surfaces.
!
!
  USE dimphy
  USE mod_phys_lmdz_para,  ONLY : mpi_size
  USE ioipsl
  USE surface_data,        ONLY : type_ocean, ok_veget
  USE surf_land_mod,       ONLY : surf_land
  USE surf_landice_mod,    ONLY : surf_landice
  USE surf_ocean_mod,      ONLY : surf_ocean
  USE surf_seaice_mod,     ONLY : surf_seaice
  USE cpl_mod,             ONLY : gath2cpl
  USE climb_hq_mod,        ONLY : climb_hq_down, climb_hq_up
  USE climb_wind_mod,      ONLY : climb_wind_down, climb_wind_up
  USE coef_diff_turb_mod,  ONLY : coef_diff_turb
  USE control_mod


  IMPLICIT NONE

! Declaration of variables saved in restart file
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: qsol   ! water height in the soil (mm)
  !$OMP THREADPRIVATE(qsol)
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: fder   ! flux drift
  !$OMP THREADPRIVATE(fder)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: snow   ! snow at surface
  !$OMP THREADPRIVATE(snow)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: qsurf  ! humidity at surface
  !$OMP THREADPRIVATE(qsurf)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: evap   ! evaporation at surface
  !$OMP THREADPRIVATE(evap)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: rugos  ! rugosity at surface (m)
  !$OMP THREADPRIVATE(rugos)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE, SAVE   :: agesno ! age of snow at surface
  !$OMP THREADPRIVATE(agesno)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE :: ftsoil ! soil temperature
  !$OMP THREADPRIVATE(ftsoil)

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE pbl_surface_init(qsol_rst, fder_rst, snow_rst, qsurf_rst,&
       evap_rst, rugos_rst, agesno_rst, ftsoil_rst)

! This routine should be called after the restart file has been read.
! This routine initialize the restart variables and does some validation tests
! for the index of the different surfaces and tests the choice of type of ocean.

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "iniprint.h"
 
! Input variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(IN)                 :: qsol_rst
    REAL, DIMENSION(klon), INTENT(IN)                 :: fder_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: snow_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: qsurf_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: evap_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: rugos_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)          :: agesno_rst
    REAL, DIMENSION(klon, nsoilmx, nbsrf), INTENT(IN) :: ftsoil_rst

  
! Local variables
!****************************************************************************************
    INTEGER                       :: ierr
    CHARACTER(len=80)             :: abort_message
    CHARACTER(len = 20)           :: modname = 'pbl_surface_init'
    

!****************************************************************************************
! Allocate and initialize module variables with fields read from restart file.
!
!****************************************************************************************    
    ALLOCATE(qsol(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(fder(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(snow(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(qsurf(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(evap(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(rugos(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(agesno(klon,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)

    ALLOCATE(ftsoil(klon,nsoilmx,nbsrf), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm('pbl_surface_init', 'pb in allocation',1)


    qsol(:)       = qsol_rst(:)
    fder(:)       = fder_rst(:)
    snow(:,:)     = snow_rst(:,:)
    qsurf(:,:)    = qsurf_rst(:,:)
    evap(:,:)     = evap_rst(:,:)
    rugos(:,:)    = rugos_rst(:,:)
    agesno(:,:)   = agesno_rst(:,:)
    ftsoil(:,:,:) = ftsoil_rst(:,:,:)


!****************************************************************************************
! Test for sub-surface indices
!
!****************************************************************************************
    IF (is_ter /= 1) THEN 
      WRITE(lunout,*)" *** Warning ***"
      WRITE(lunout,*)" is_ter n'est pas le premier surface, is_ter = ",is_ter
      WRITE(lunout,*)"or on doit commencer par les surfaces continentales"
      abort_message="voir ci-dessus"
      CALL abort_gcm(modname,abort_message,1)
    ENDIF

    IF ( is_oce > is_sic ) THEN
      WRITE(lunout,*)' *** Warning ***'
      WRITE(lunout,*)' Pour des raisons de sequencement dans le code'
      WRITE(lunout,*)' l''ocean doit etre traite avant la banquise'
      WRITE(lunout,*)' or is_oce = ',is_oce, '> is_sic = ',is_sic
      abort_message='voir ci-dessus'
      CALL abort_gcm(modname,abort_message,1)
    ENDIF

    IF ( is_lic > is_sic ) THEN
      WRITE(lunout,*)' *** Warning ***'
      WRITE(lunout,*)' Pour des raisons de sequencement dans le code'
      WRITE(lunout,*)' la glace contineltalle doit etre traite avant la glace de mer'
      WRITE(lunout,*)' or is_lic = ',is_lic, '> is_sic = ',is_sic
      abort_message='voir ci-dessus'
      CALL abort_gcm(modname,abort_message,1)
    ENDIF

!****************************************************************************************
! Validation of ocean mode
!
!****************************************************************************************

    IF (type_ocean /= 'slab  ' .AND. type_ocean /= 'force ' .AND. type_ocean /= 'couple') THEN
       WRITE(lunout,*)' *** Warning ***'
       WRITE(lunout,*)'Option couplage pour l''ocean = ', type_ocean
       abort_message='option pour l''ocean non valable'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF

  END SUBROUTINE pbl_surface_init
!  
!****************************************************************************************
!  

  SUBROUTINE pbl_surface( &
       dtime,     date0,     itap,     jour,          &
       debut,     lafin,                              &
       rlon,      rlat,      rugoro,   rmu0,          &
       rain_f,    snow_f,    solsw_m,  sollw_m,       &
       t,         q,         u,        v,             &
       pplay,     paprs,     pctsrf,                  &
       ts,        alb1,      alb2,     u10m,   v10m,  &
       lwdown_m,  cdragh,    cdragm,   zu1,    zv1,   &
       alb1_m,    alb2_m,    zxsens,   zxevap,        &
       zxtsol,    zxfluxlat, zt2m,     qsat2m,        &
       d_t,       d_q,       d_u,      d_v,           & 
       zcoefh,    slab_wfbils,                        &
       qsol_d,    zq2m,      s_pblh,   s_plcl,        &
       s_capCL,   s_oliqCL,  s_cteiCL, s_pblT,        &
       s_therm,   s_trmb1,   s_trmb2,  s_trmb3,       &
       zxrugs,    zu10m,     zv10m,    fder_print,    &
       zxqsurf,   rh2m,      zxfluxu,  zxfluxv,       &
       rugos_d,   agesno_d,  sollw,    solsw,         &
       d_ts,      evap_d,    fluxlat,  t2m,           &
       wfbils,    wfbilo,    flux_t,   flux_u, flux_v,&
       dflux_t,   dflux_q,   zxsnow,                  &
       zxfluxt,   zxfluxq,   q2m,      flux_q, tke    )
!****************************************************************************************
! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
! Objet: interface de "couche limite" (diffusion verticale)
!
!AA REM:
!AA-----
!AA Tout ce qui a trait au traceurs est dans phytrac maintenant
!AA pour l'instant le calcul de la couche limite pour les traceurs
!AA se fait avec cltrac et ne tient pas compte de la differentiation
!AA des sous-fraction de sol.
!AA REM bis :
!AA----------
!AA Pour pouvoir extraire les coefficient d'echanges et le vent 
!AA dans la premiere couche, 3 champs supplementaires ont ete crees
!AA zcoefh, zu1 et zv1. Pour l'instant nous avons moyenne les valeurs
!AA de ces trois champs sur les 4 subsurfaces du modele. Dans l'avenir 
!AA si les informations des subsurfaces doivent etre prises en compte
!AA il faudra sortir ces memes champs en leur ajoutant une dimension, 
!AA c'est a dire nbsrf (nbre de subsurface).
!
! Arguments:
!
! dtime----input-R- interval du temps (secondes)
! itap-----input-I- numero du pas de temps
! date0----input-R- jour initial
! t--------input-R- temperature (K)
! q--------input-R- vapeur d'eau (kg/kg)
! u--------input-R- vitesse u
! v--------input-R- vitesse v
! ts-------input-R- temperature du sol (en Kelvin)
! paprs----input-R- pression a intercouche (Pa)
! pplay----input-R- pression au milieu de couche (Pa)
! rlat-----input-R- latitude en degree
! rugos----input-R- longeur de rugosite (en m)
!
! d_t------output-R- le changement pour "t"
! d_q------output-R- le changement pour "q"
! d_u------output-R- le changement pour "u"
! d_v------output-R- le changement pour "v"
! d_ts-----output-R- le changement pour "ts"
! flux_t---output-R- flux de chaleur sensible (CpT) J/m**2/s (W/m**2)
!                    (orientation positive vers le bas)
! tke---input/output-R- tke (kg/m**2/s)
! flux_q---output-R- flux de vapeur d'eau (kg/m**2/s)
! flux_u---output-R- tension du vent X: (kg m/s)/(m**2 s) ou Pascal
! flux_v---output-R- tension du vent Y: (kg m/s)/(m**2 s) ou Pascal
! dflux_t--output-R- derive du flux sensible
! dflux_q--output-R- derive du flux latent
! zu1------output-R- le vent dans la premiere couche
! zv1------output-R- le vent dans la premiere couche
! trmb1----output-R- deep_cape
! trmb2----output-R- inhibition 
! trmb3----output-R- Point Omega
! cteiCL---output-R- Critere d'instab d'entrainmt des nuages de CL
! plcl-----output-R- Niveau de condensation
! pblh-----output-R- HCL
! pblT-----output-R- T au nveau HCL
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, co2_send
    IMPLICIT NONE

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "iniprint.h"
    INCLUDE "FCTTRE.h"
    INCLUDE "clesphys.h"
    INCLUDE "compbl.h"
    INCLUDE "dimensions.h"
    INCLUDE "YOETHF.h"
    INCLUDE "temps.h"
! Input variables
!****************************************************************************************
    REAL,                         INTENT(IN)        :: dtime   ! time interval (s)
    REAL,                         INTENT(IN)        :: date0   ! initial day
    INTEGER,                      INTENT(IN)        :: itap    ! time step
    INTEGER,                      INTENT(IN)        :: jour    ! current day of the year
    LOGICAL,                      INTENT(IN)        :: debut   ! true if first run step
    LOGICAL,                      INTENT(IN)        :: lafin   ! true if last run step
    REAL, DIMENSION(klon),        INTENT(IN)        :: rlon    ! longitudes in degrees
    REAL, DIMENSION(klon),        INTENT(IN)        :: rlat    ! latitudes in degrees
    REAL, DIMENSION(klon),        INTENT(IN)        :: rugoro  ! rugosity length
    REAL, DIMENSION(klon),        INTENT(IN)        :: rmu0    ! cosine of solar zenith angle
    REAL, DIMENSION(klon),        INTENT(IN)        :: rain_f  ! rain fall
    REAL, DIMENSION(klon),        INTENT(IN)        :: snow_f  ! snow fall
    REAL, DIMENSION(klon),        INTENT(IN)        :: solsw_m ! net shortwave radiation at mean surface
    REAL, DIMENSION(klon),        INTENT(IN)        :: sollw_m ! net longwave radiation at mean surface
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: t       ! temperature (K)
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: q       ! water vapour (kg/kg)
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: u       ! u speed
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: v       ! v speed
    REAL, DIMENSION(klon,klev),   INTENT(IN)        :: pplay   ! mid-layer pression (Pa)
    REAL, DIMENSION(klon,klev+1), INTENT(IN)        :: paprs   ! pression between layers (Pa) 
    REAL, DIMENSION(klon, nbsrf), INTENT(IN)        :: pctsrf  ! sub-surface fraction

! Input/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: ts      ! temperature at surface (K)
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: alb1    ! albedo in visible SW interval
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: alb2    ! albedo in near infra-red SW interval
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: u10m    ! u speed at 10m
    REAL, DIMENSION(klon, nbsrf), INTENT(INOUT)     :: v10m    ! v speed at 10m
    REAL, DIMENSION(klon, klev+1, nbsrf), INTENT(INOUT) :: tke 

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon),        INTENT(OUT)       :: lwdown_m   ! Downcoming longwave radiation
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragh     ! drag coefficient for T and Q
    REAL, DIMENSION(klon),        INTENT(OUT)       :: cdragm     ! drag coefficient for wind
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zu1        ! u wind speed in first layer
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zv1        ! v wind speed in first layer
    REAL, DIMENSION(klon),        INTENT(OUT)       :: alb1_m     ! mean albedo in visible SW interval
    REAL, DIMENSION(klon),        INTENT(OUT)       :: alb2_m     ! mean albedo in near IR SW interval
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxsens     ! sensible heat flux at surface with inversed sign 
                                                                  ! (=> positive sign upwards)
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxevap     ! water vapour flux at surface, positiv upwards
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxtsol     ! temperature at surface, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxfluxlat  ! latent flux, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zt2m       ! temperature at 2m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: qsat2m
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_t        ! change in temperature 
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_q        ! change in water vapour
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_u        ! change in u speed
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: d_v        ! change in v speed
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: zcoefh     ! coef for turbulent diffusion of T and Q, mean for each grid point

! Output only for diagnostics
    REAL, DIMENSION(klon),        INTENT(OUT)       :: slab_wfbils! heat balance at surface only for slab at ocean points
    REAL, DIMENSION(klon),        INTENT(OUT)       :: qsol_d     ! water height in the soil (mm)
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zq2m       ! water vapour at 2m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblh     ! height of the planetary boundary layer(HPBL)
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_plcl     ! condensation level
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_capCL    ! CAPE of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_oliqCL   ! liquid water intergral of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_cteiCL   ! cloud top instab. crit. of PBL
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_pblT     ! temperature at PBLH
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_therm    ! thermal virtual temperature excess
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb1    ! deep cape, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb2    ! inhibition, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: s_trmb3    ! point Omega, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxrugs     ! rugosity at surface (m), mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zu10m      ! u speed at 10m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zv10m      ! v speed at 10m, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: fder_print ! fder for printing (=fder(i) + dflux_t(i) + dflux_q(i))
    REAL, DIMENSION(klon),        INTENT(OUT)       :: zxqsurf    ! humidity at surface, mean for each grid point
    REAL, DIMENSION(klon),        INTENT(OUT)       :: rh2m       ! relative humidity at 2m
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: zxfluxu    ! u wind tension, mean for each grid point
    REAL, DIMENSION(klon, klev),  INTENT(OUT)       :: zxfluxv    ! v wind tension, mean for each grid point
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: rugos_d    ! rugosity length (m)
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: agesno_d   ! age of snow at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: solsw      ! net shortwave radiation at surface 
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: sollw      ! net longwave radiation at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: d_ts       ! change in temperature at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: evap_d     ! evaporation at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: fluxlat    ! latent flux
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: t2m        ! temperature at 2 meter height
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfbils     ! heat balance at surface
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)       :: wfbilo     ! water balance at surface
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_t     ! sensible heat flux (CpT) J/m**2/s (W/m**2)
                                                                  ! positve orientation downwards
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_u     ! u wind tension (kg m/s)/(m**2 s) or Pascal
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_v     ! v wind tension (kg m/s)/(m**2 s) or Pascal

! Output not needed
    REAL, DIMENSION(klon),       INTENT(OUT)        :: dflux_t    ! change of sensible heat flux 
    REAL, DIMENSION(klon),       INTENT(OUT)        :: dflux_q    ! change of water vapour flux
    REAL, DIMENSION(klon),       INTENT(OUT)        :: zxsnow     ! snow at surface, mean for each grid point
    REAL, DIMENSION(klon, klev), INTENT(OUT)        :: zxfluxt    ! sensible heat flux, mean for each grid point
    REAL, DIMENSION(klon, klev), INTENT(OUT)        :: zxfluxq    ! water vapour flux, mean for each grid point
    REAL, DIMENSION(klon, nbsrf),INTENT(OUT)        :: q2m        ! water vapour at 2 meter height
    REAL, DIMENSION(klon, klev, nbsrf), INTENT(OUT) :: flux_q     ! water vapour flux(latent flux) (kg/m**2/s)


! Local variables with attribute SAVE
!****************************************************************************************
    INTEGER, SAVE                            :: nhoridbg, nidbg   ! variables for IOIPSL
!$OMP THREADPRIVATE(nhoridbg, nidbg)
    LOGICAL, SAVE                            :: debugindex=.FALSE.
!$OMP THREADPRIVATE(debugindex)
    LOGICAL, SAVE                            :: first_call=.TRUE.
!$OMP THREADPRIVATE(first_call)
    CHARACTER(len=8), DIMENSION(nbsrf), SAVE :: cl_surf
!$OMP THREADPRIVATE(cl_surf)

! Other local variables
!****************************************************************************************
    INTEGER                            :: i, k, nsrf 
    INTEGER                            :: knon, j
    INTEGER                            :: idayref
    INTEGER , DIMENSION(klon)          :: ni
    REAL                               :: zx_alf1, zx_alf2 !valeur ambiante par extrapola
    REAL                               :: amn, amx
    REAL                               :: f1 ! fraction de longeurs visibles parmi tout SW intervalle
    REAL, DIMENSION(klon)              :: r_co2_ppm     ! taux CO2 atmosphere
    REAL, DIMENSION(klon)              :: yts, yrugos, ypct, yz0_new
    REAL, DIMENSION(klon)              :: yalb, yalb1, yalb2
    REAL, DIMENSION(klon)              :: yu1, yv1
    REAL, DIMENSION(klon)              :: ysnow, yqsurf, yagesno, yqsol
    REAL, DIMENSION(klon)              :: yrain_f, ysnow_f
    REAL, DIMENSION(klon)              :: ysolsw, ysollw
    REAL, DIMENSION(klon)              :: yfder
    REAL, DIMENSION(klon)              :: yrugoro
    REAL, DIMENSION(klon)              :: yfluxlat
    REAL, DIMENSION(klon)              :: y_d_ts
    REAL, DIMENSION(klon)              :: y_flux_t1, y_flux_q1
    REAL, DIMENSION(klon)              :: y_dflux_t, y_dflux_q
    REAL, DIMENSION(klon)              :: y_flux_u1, y_flux_v1
    REAL, DIMENSION(klon)              :: yt2m, yq2m, yu10m
    REAL, DIMENSION(klon)              :: yustar
    REAL, DIMENSION(klon)              :: ywindsp
    REAL, DIMENSION(klon)              :: yt10m, yq10m
    REAL, DIMENSION(klon)              :: ypblh
    REAL, DIMENSION(klon)              :: ylcl
    REAL, DIMENSION(klon)              :: ycapCL
    REAL, DIMENSION(klon)              :: yoliqCL
    REAL, DIMENSION(klon)              :: ycteiCL
    REAL, DIMENSION(klon)              :: ypblT
    REAL, DIMENSION(klon)              :: ytherm
    REAL, DIMENSION(klon)              :: ytrmb1
    REAL, DIMENSION(klon)              :: ytrmb2
    REAL, DIMENSION(klon)              :: ytrmb3
    REAL, DIMENSION(klon)              :: uzon, vmer
    REAL, DIMENSION(klon)              :: tair1, qair1, tairsol
    REAL, DIMENSION(klon)              :: psfce, patm
    REAL, DIMENSION(klon)              :: qairsol, zgeo1
    REAL, DIMENSION(klon)              :: rugo1
    REAL, DIMENSION(klon)              :: yfluxsens
    REAL, DIMENSION(klon)              :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon)              :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon)              :: ypsref
    REAL, DIMENSION(klon)              :: yevap, ytsurf_new, yalb1_new, yalb2_new
    REAL, DIMENSION(klon)              :: ztsol
    REAL, DIMENSION(klon)              :: alb_m  ! mean albedo for whole SW interval
    REAL, DIMENSION(klon,klev)         :: y_d_t, y_d_q
    REAL, DIMENSION(klon,klev)         :: y_d_u, y_d_v
    REAL, DIMENSION(klon,klev)         :: y_flux_t, y_flux_q
    REAL, DIMENSION(klon,klev)         :: y_flux_u, y_flux_v
    REAL, DIMENSION(klon,klev)         :: ycoefh, ycoefm
    REAL, DIMENSION(klon)              :: ycdragh, ycdragm
    REAL, DIMENSION(klon,klev)         :: yu, yv
    REAL, DIMENSION(klon,klev)         :: yt, yq
    REAL, DIMENSION(klon,klev)         :: ypplay, ydelp
    REAL, DIMENSION(klon,klev)         :: delp
    REAL, DIMENSION(klon,klev+1)       :: ypaprs
    REAL, DIMENSION(klon,klev+1)       :: ytke
    REAL, DIMENSION(klon,nsoilmx)      :: ytsoil
    CHARACTER(len=80)                  :: abort_message
    CHARACTER(len=20)                  :: modname = 'pbl_surface'
    LOGICAL, PARAMETER                 :: zxli=.FALSE. ! utiliser un jeu de fonctions simples
    LOGICAL, PARAMETER                 :: check=.FALSE.

! For debugging with IOIPSL
    INTEGER, DIMENSION(iim*(jjm+1))    :: ndexbg
    REAL                               :: zjulian
    REAL, DIMENSION(klon)              :: tabindx
    REAL, DIMENSION(iim,jjm+1)         :: zx_lon, zx_lat
    REAL, DIMENSION(iim,jjm+1)         :: debugtab


    REAL, DIMENSION(klon,nbsrf)        :: pblh         ! height of the planetary boundary layer
    REAL, DIMENSION(klon,nbsrf)        :: plcl         ! condensation level
    REAL, DIMENSION(klon,nbsrf)        :: capCL
    REAL, DIMENSION(klon,nbsrf)        :: oliqCL
    REAL, DIMENSION(klon,nbsrf)        :: cteiCL
    REAL, DIMENSION(klon,nbsrf)        :: pblT
    REAL, DIMENSION(klon,nbsrf)        :: therm
    REAL, DIMENSION(klon,nbsrf)        :: trmb1        ! deep cape
    REAL, DIMENSION(klon,nbsrf)        :: trmb2        ! inhibition
    REAL, DIMENSION(klon,nbsrf)        :: trmb3        ! point Omega
    REAL, DIMENSION(klon,nbsrf)        :: zx_rh2m, zx_qsat2m
    REAL, DIMENSION(klon,nbsrf)        :: zx_t1
    REAL, DIMENSION(klon, nbsrf)       :: alb          ! mean albedo for whole SW interval
    REAL, DIMENSION(klon)              :: ylwdown      ! jg : temporary (ysollwdown)

    REAL                               :: zx_qs1, zcor1, zdelta1 

!****************************************************************************************
! Declarations specifiques pour le 1D. A reprendre 
  REAL  :: fsens,flat
  LOGICAL :: ok_flux_surf ! initialized during first_call below
  COMMON /flux_arp/fsens,flat,ok_flux_surf
!****************************************************************************************
! End of declarations
!****************************************************************************************


!****************************************************************************************
! 1) Initialisation and validation tests 
!    Only done first time entering this subroutine
!
!****************************************************************************************

    IF (first_call) THEN
       first_call=.FALSE.
      
       ! Initialize ok_flux_surf (for 1D model)
       if (klon>1) ok_flux_surf=.FALSE.
       
       ! Initilize debug IO
       IF (debugindex .AND. mpi_size==1) THEN 
          ! initialize IOIPSL output
          idayref = day_ini
          CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
          DO i = 1, iim
             zx_lon(i,1) = rlon(i+1)
             zx_lon(i,jjm+1) = rlon(i+1)
          ENDDO
          CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)
          CALL histbeg("sous_index", iim,zx_lon(:,1),jjm+1,zx_lat(1,:), &
               1,iim,1,jjm+1, &
               itau_phy,zjulian,dtime,nhoridbg,nidbg) 
          ! no vertical axis
          cl_surf(1)='ter'
          cl_surf(2)='lic'
          cl_surf(3)='oce'
          cl_surf(4)='sic'
          DO nsrf=1,nbsrf
             CALL histdef(nidbg, cl_surf(nsrf),cl_surf(nsrf), "-",iim, &
                  jjm+1,nhoridbg, 1, 1, 1, -99, 32, "inst", dtime,dtime) 
          END DO

          CALL histend(nidbg)
          CALL histsync(nidbg)

       END IF
       
    ENDIF
          
!****************************************************************************************
! Force soil water content to qsol0 if qsol0>0 and VEGET=F (use bucket
! instead of ORCHIDEE)
    IF (qsol0>0.) THEN
      PRINT*,'WARNING : On impose qsol=',qsol0
      qsol(:)=qsol0
    ENDIF
!****************************************************************************************

!****************************************************************************************
! 2) Initialization to zero 
!    Done for all local variables that will be compressed later
!    and argument with INTENT(OUT)
!****************************************************************************************
    cdragh = 0.0  ; cdragm = 0.0     ; dflux_t = 0.0   ; dflux_q = 0.0
    ypct = 0.0    ; yts = 0.0        ; ysnow = 0.0
    zv1 = 0.0     ; yqsurf = 0.0     ; yalb1 = 0.0     ; yalb2 = 0.0    
    yrain_f = 0.0 ; ysnow_f = 0.0    ; yfder = 0.0     ; ysolsw = 0.0    
    ysollw = 0.0  ; yrugos = 0.0     ; yu1 = 0.0    
    yv1 = 0.0     ; ypaprs = 0.0     ; ypplay = 0.0
    ydelp = 0.0   ; yu = 0.0         ; yv = 0.0        ; yt = 0.0         
    yq = 0.0      ; y_dflux_t = 0.0  ; y_dflux_q = 0.0 
    yrugoro = 0.0 ; ywindsp = 0.0   
    d_ts = 0.0    ; yfluxlat=0.0     ; flux_t = 0.0    ; flux_q = 0.0     
    flux_u = 0.0  ; flux_v = 0.0     ; d_t = 0.0       ; d_q = 0.0      
    d_u = 0.0     ; d_v = 0.0        ; yqsol = 0.0    
    ytherm = 0.0  ; ytke=0.
    
    zcoefh(:,:) = 0.0
    zcoefh(:,1) = 999999. ! zcoefh(:,k=1) should never be used
    ytsoil = 999999. 

    rh2m(:)        = 0.
    qsat2m(:)      = 0.
!****************************************************************************************
! 3) - Calculate pressure thickness of each layer
!    - Calculate the wind at first layer
!    - Mean calculations of albedo
!    - Calculate net radiance at sub-surface
!****************************************************************************************
    DO k = 1, klev
       DO i = 1, klon
          delp(i,k) = paprs(i,k)-paprs(i,k+1)
       ENDDO
    ENDDO

!****************************************************************************************
! Test for rugos........ from physiq.. A la fin plutot???
!
!****************************************************************************************

    zxrugs(:) = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          rugos(i,nsrf) = MAX(rugos(i,nsrf),0.000015)
          zxrugs(i) = zxrugs(i) + rugos(i,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO

! Mean calculations of albedo
!
! Albedo at sub-surface
! * alb1 : albedo in visible SW interval
! * alb2 : albedo in near infrared SW interval
! * alb  : mean albedo for whole SW interval
!
! Mean albedo for grid point
! * alb1_m : albedo in visible SW interval
! * alb2_m : albedo in near infrared SW interval
! * alb_m  : mean albedo at whole SW interval

    alb1_m(:) = 0.0
    alb2_m(:) = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          alb1_m(i) = alb1_m(i) + alb1(i,nsrf) * pctsrf(i,nsrf)
          alb2_m(i) = alb2_m(i) + alb2(i,nsrf) * pctsrf(i,nsrf)
       ENDDO
    ENDDO

! We here suppose the fraction f1 of incoming radiance of visible radiance 
! as a fraction of all shortwave radiance 
    f1 = 0.5 
!    f1 = 1    ! put f1=1 to recreate old calculations

    DO nsrf = 1, nbsrf
       DO i = 1, klon
          alb(i,nsrf) = f1*alb1(i,nsrf) + (1-f1)*alb2(i,nsrf)
       ENDDO
    ENDDO

    DO i = 1, klon
       alb_m(i) = f1*alb1_m(i) + (1-f1)*alb2_m(i)
    END DO

! Calculation of mean temperature at surface grid points
    ztsol(:) = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          ztsol(i) = ztsol(i) + ts(i,nsrf)*pctsrf(i,nsrf)
       ENDDO
    ENDDO

! Linear distrubution on sub-surface of long- and shortwave net radiance
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          sollw(i,nsrf) = sollw_m(i) + 4.0*RSIGMA*ztsol(i)**3 * (ztsol(i)-ts(i,nsrf))
          solsw(i,nsrf) = solsw_m(i) * (1.-alb(i,nsrf)) / (1.-alb_m(i))
       ENDDO
    ENDDO


! Downwelling longwave radiation at mean surface
    lwdown_m(:) = 0.0
    DO i = 1, klon
       lwdown_m(i) = sollw_m(i) + RSIGMA*ztsol(i)**4
    ENDDO

!****************************************************************************************
! 4) Loop over different surfaces
!
! Only points containing a fraction of the sub surface will be threated.
! 
!****************************************************************************************
   
    loop_nbsrf: DO nsrf = 1, nbsrf

! Search for index(ni) and size(knon) of domaine to treat
       ni(:) = 0
       knon  = 0
       DO i = 1, klon
          IF (pctsrf(i,nsrf) > 0.) THEN
             knon = knon + 1
             ni(knon) = i
          ENDIF
       ENDDO

       ! write index, with IOIPSL
       IF (debugindex .AND. mpi_size==1) THEN 
          tabindx(:)=0.
          DO i=1,knon
             tabindx(i)=REAL(i)
          END DO
          debugtab(:,:) = 0.
          ndexbg(:) = 0
          CALL gath2cpl(tabindx,debugtab,knon,ni)
          CALL histwrite(nidbg,cl_surf(nsrf),itap,debugtab,iim*(jjm+1), ndexbg)
       ENDIF
       
!****************************************************************************************
! 5) Compress variables 
!
!****************************************************************************************

       DO j = 1, knon
          i = ni(j)
          ypct(j)    = pctsrf(i,nsrf)
          yts(j)     = ts(i,nsrf)
          ysnow(j)   = snow(i,nsrf)
          yqsurf(j)  = qsurf(i,nsrf)
          yalb(j)    = alb(i,nsrf)
          yalb1(j)   = alb1(i,nsrf)
          yalb2(j)   = alb2(i,nsrf)
          yrain_f(j) = rain_f(i)
          ysnow_f(j) = snow_f(i)
          yagesno(j) = agesno(i,nsrf)
          yfder(j)   = fder(i)
          ysolsw(j)  = solsw(i,nsrf)
          ysollw(j)  = sollw(i,nsrf)
          yrugos(j)  = rugos(i,nsrf)
          yrugoro(j) = rugoro(i)
          yu1(j)     = u(i,1)
          yv1(j)     = v(i,1)
          ypaprs(j,klev+1) = paprs(i,klev+1)
          ywindsp(j) = SQRT(u10m(i,nsrf)**2 + v10m(i,nsrf)**2 )
       END DO

       DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             ypaprs(j,k) = paprs(i,k)
             ypplay(j,k) = pplay(i,k)
             ydelp(j,k)  = delp(i,k)
             ytke(j,k)   = tke(i,k,nsrf)
             yu(j,k) = u(i,k)
             yv(j,k) = v(i,k)
             yt(j,k) = t(i,k)
             yq(j,k) = q(i,k)
          ENDDO
       ENDDO
       
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ytsoil(j,k) = ftsoil(i,k,nsrf)
          END DO
       END DO
       
       ! qsol(water height in soil) only for bucket continental model
       IF ( nsrf .EQ. is_ter .AND. .NOT. ok_veget ) THEN 
          DO j = 1, knon
             i = ni(j)
             yqsol(j) = qsol(i)
          END DO
       ENDIF
       
!****************************************************************************************
! 6a) Calculate coefficients for turbulent diffusion at surface, cdragh et cdragm.
!
!****************************************************************************************

       CALL clcdrag( knon, nsrf, ypaprs, ypplay, &
            yu(:,1), yv(:,1), yt(:,1), yq(:,1), &
            yts, yqsurf, yrugos, &
            ycdragm, ycdragh )

!****************************************************************************************
! 6b) Calculate coefficients for turbulent diffusion in the atmosphere, ycoefm et ycoefm.
!
!****************************************************************************************

       CALL coef_diff_turb(dtime, nsrf, knon, ni,  &
            ypaprs, ypplay, yu, yv, yq, yt, yts, yrugos, yqsurf, ycdragm, &
            ycoefm, ycoefh, ytke)
       
!****************************************************************************************
! 
! 8) "La descente" - "The downhill"
!  
!  climb_hq_down and climb_wind_down calculate the coefficients
!  Ccoef_X et Dcoef_X for X=[H, Q, U, V].
!  Only the coefficients at surface for H and Q are returned.
!
!****************************************************************************************

! - Calculate the coefficients Ccoef_H, Ccoef_Q, Dcoef_H and Dcoef_Q 
       CALL climb_hq_down(knon, ycoefh, ypaprs, ypplay, &
            ydelp, yt, yq, dtime, &
            AcoefH, AcoefQ, BcoefH, BcoefQ)

! - Calculate the coefficients Ccoef_U, Ccoef_V, Dcoef_U and Dcoef_V
       CALL climb_wind_down(knon, dtime, ycoefm, ypplay, ypaprs, yt, ydelp, yu, yv, &
            AcoefU, AcoefV, BcoefU, BcoefV)
      

!****************************************************************************************
! 9) Small calculations
!
!****************************************************************************************

! - Reference pressure is given the values at surface level          
       ypsref(:) = ypaprs(:,1)  

! - CO2 field on 2D grid to be sent to ORCHIDEE
!   Transform to compressed field
       IF (carbon_cycle_cpl) THEN
          DO i=1,knon
             r_co2_ppm(i) = co2_send(ni(i))
          END DO
       ELSE
          r_co2_ppm(:) = co2_ppm     ! Constant field
       END IF

!****************************************************************************************
!
! Calulate t2m and q2m for the case of calculation at land grid points 
! t2m and q2m are needed as input to ORCHIDEE
!
!****************************************************************************************
       IF (nsrf == is_ter) THEN

          DO i = 1, knon
             zgeo1(i) = RD * yt(i,1) / (0.5*(ypaprs(i,1)+ypplay(i,1))) &
                  * (ypaprs(i,1)-ypplay(i,1))
          END DO

          ! Calculate the temperature et relative humidity at 2m and the wind at 10m 
          CALL stdlevvar(klon, knon, is_ter, zxli, &
               yu(:,1), yv(:,1), yt(:,1), yq(:,1), zgeo1, &
               yts, yqsurf, yrugos, ypaprs(:,1), ypplay(:,1), &
               yt2m, yq2m, yt10m, yq10m, yu10m, yustar)
          
       END IF

!****************************************************************************************
!
! 10) Switch selon current surface
!     It is necessary to start with the continental surfaces because the ocean
!     needs their run-off.
!
!****************************************************************************************
       SELECT CASE(nsrf)
     
       CASE(is_ter)
          ! ylwdown : to be removed, calculation is now done at land surface in surf_land
          ylwdown(:)=0.0
          DO i=1,knon
             ylwdown(i)=lwdown_m(ni(i))
          END DO
          CALL surf_land(itap, dtime, date0, jour, knon, ni,&
               rlon, rlat, &
               debut, lafin, ydelp(:,1), r_co2_ppm, ysolsw, ysollw, yalb, &
               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               AcoefH, AcoefQ, BcoefH, BcoefQ, & 
               AcoefU, AcoefV, BcoefU, BcoefV, & 
               ypsref, yu1, yv1, yrugoro, pctsrf, &
               ylwdown, yq2m, yt2m, &
               ysnow, yqsol, yagesno, ytsoil, &
               yz0_new, yalb1_new, yalb2_new, yevap, yfluxsens, yfluxlat, &
               yqsurf, ytsurf_new, y_dflux_t, y_dflux_q, &
               y_flux_u1, y_flux_v1 )
               
     
       CASE(is_lic)
          CALL surf_landice(itap, dtime, knon, ni, &
               ysolsw, ysollw, yts, ypplay(:,1), &
               ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, yrugoro, pctsrf, &
               ysnow, yqsurf, yqsol, yagesno, &
               ytsoil, yz0_new, yalb1_new, yalb2_new, yevap, yfluxsens, yfluxlat, &
               ytsurf_new, y_dflux_t, y_dflux_q, &
               y_flux_u1, y_flux_v1)
          
       CASE(is_oce)
          CALL surf_ocean(rlon, rlat, ysolsw, ysollw, yalb1, &
               yrugos, ywindsp, rmu0, yfder, yts, &
               itap, dtime, jour, knon, ni, &
               ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, yrugoro, pctsrf, &
               ysnow, yqsurf, yagesno, &
               yz0_new, yalb1_new, yalb2_new, yevap, yfluxsens, yfluxlat, &
               ytsurf_new, y_dflux_t, y_dflux_q, slab_wfbils, &
               y_flux_u1, y_flux_v1)
          
       CASE(is_sic)
          CALL surf_seaice( &
               rlon, rlat, ysolsw, ysollw, yalb1, yfder, &
               itap, dtime, jour, knon, ni, &
               lafin, &
               yts, ypplay(:,1), ycdragh, ycdragm, yrain_f, ysnow_f, yt(:,1), yq(:,1),&
               AcoefH, AcoefQ, BcoefH, BcoefQ, &
               AcoefU, AcoefV, BcoefU, BcoefV, &
               ypsref, yu1, yv1, yrugoro, pctsrf, &
               ysnow, yqsurf, yqsol, yagesno, ytsoil, &
               yz0_new, yalb1_new, yalb2_new, yevap, yfluxsens, yfluxlat, &
               ytsurf_new, y_dflux_t, y_dflux_q, &
               y_flux_u1, y_flux_v1)
          

       CASE DEFAULT
          WRITE(lunout,*) 'Surface index = ', nsrf
          abort_message = 'Surface index not valid'
          CALL abort_gcm(modname,abort_message,1)
       END SELECT


!****************************************************************************************
! 11) - Calcul the increment of surface temperature
!
!****************************************************************************************
       y_d_ts(1:knon)   = ytsurf_new(1:knon) - yts(1:knon)
 
!****************************************************************************************
!
! 12) "La remontee" - "The uphill"
!
!  The fluxes (y_flux_X) and tendancy (y_d_X) are calculated 
!  for X=H, Q, U and V, for all vertical levels.
!
!****************************************************************************************
! H and Q
       IF (ok_flux_surf) THEN
          PRINT *,'pbl_surface: fsens flat RLVTT=',fsens,flat,RLVTT
          y_flux_t1(:) =  fsens
          y_flux_q1(:) =  flat/RLVTT
          yfluxlat(:) =  flat
       ELSE
          y_flux_t1(:) =  yfluxsens(:)
          y_flux_q1(:) = -yevap(:)
       ENDIF

       CALL climb_hq_up(knon, dtime, yt, yq, &
            y_flux_q1, y_flux_t1, ypaprs, ypplay, &
            y_flux_q(:,:), y_flux_t(:,:), y_d_q(:,:), y_d_t(:,:))    
       

       CALL climb_wind_up(knon, dtime, yu, yv, y_flux_u1, y_flux_v1, &
            y_flux_u, y_flux_v, y_d_u, y_d_v)


       DO j = 1, knon
          y_dflux_t(j) = y_dflux_t(j) * ypct(j)
          y_dflux_q(j) = y_dflux_q(j) * ypct(j)
       ENDDO

!****************************************************************************************
! 13) Transform variables for output format : 
!     - Decompress
!     - Multiply with pourcentage of current surface
!     - Cumulate in global variable
!
!****************************************************************************************

       tke(:,:,nsrf) = 0.
       DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             y_d_t(j,k)  = y_d_t(j,k) * ypct(j)
             y_d_q(j,k)  = y_d_q(j,k) * ypct(j)
             y_d_u(j,k)  = y_d_u(j,k) * ypct(j)
             y_d_v(j,k)  = y_d_v(j,k) * ypct(j)

             flux_t(i,k,nsrf) = y_flux_t(j,k)
             flux_q(i,k,nsrf) = y_flux_q(j,k)
             flux_u(i,k,nsrf) = y_flux_u(j,k)
             flux_v(i,k,nsrf) = y_flux_v(j,k)

             tke(i,k,nsrf)    = ytke(j,k)

          ENDDO
       ENDDO

       evap(:,nsrf) = - flux_q(:,1,nsrf)
       
       alb1(:, nsrf) = 0.
       alb2(:, nsrf) = 0.
       snow(:, nsrf) = 0.
       qsurf(:, nsrf) = 0.
       rugos(:, nsrf) = 0.
       fluxlat(:,nsrf) = 0.
       DO j = 1, knon
          i = ni(j)
          d_ts(i,nsrf) = y_d_ts(j)
          alb1(i,nsrf) = yalb1_new(j)  
          alb2(i,nsrf) = yalb2_new(j)
          snow(i,nsrf) = ysnow(j)  
          qsurf(i,nsrf) = yqsurf(j)
          rugos(i,nsrf) = yz0_new(j)
          fluxlat(i,nsrf) = yfluxlat(j)
          agesno(i,nsrf) = yagesno(j)  
          cdragh(i) = cdragh(i) + ycdragh(j)*ypct(j)
          cdragm(i) = cdragm(i) + ycdragm(j)*ypct(j)
          dflux_t(i) = dflux_t(i) + y_dflux_t(j)
          dflux_q(i) = dflux_q(i) + y_dflux_q(j)
       END DO

       DO k = 2, klev
          DO j = 1, knon
             i = ni(j)
             zcoefh(i,k) = zcoefh(i,k) + ycoefh(j,k)*ypct(j)
          END DO
       END DO

       IF ( nsrf .EQ. is_ter ) THEN 
          DO j = 1, knon
             i = ni(j)
             qsol(i) = yqsol(j)
          END DO
       END IF
       
       ftsoil(:,:,nsrf) = 0.
       DO k = 1, nsoilmx
          DO j = 1, knon
             i = ni(j)
             ftsoil(i, k, nsrf) = ytsoil(j,k)
          END DO
       END DO
       
       
       DO k = 1, klev
          DO j = 1, knon
             i = ni(j)
             d_t(i,k) = d_t(i,k) + y_d_t(j,k)
             d_q(i,k) = d_q(i,k) + y_d_q(j,k)
             d_u(i,k) = d_u(i,k) + y_d_u(j,k)
             d_v(i,k) = d_v(i,k) + y_d_v(j,k)
          END DO
       END DO

!****************************************************************************************
! 14) Calculate the temperature et relative humidity at 2m and the wind at 10m 
!     Call HBTM
!
!****************************************************************************************
       t2m(:,nsrf)    = 0.
       q2m(:,nsrf)    = 0.
       u10m(:,nsrf)   = 0.
       v10m(:,nsrf)   = 0.

       pblh(:,nsrf)   = 0.        ! Hauteur de couche limite
       plcl(:,nsrf)   = 0.        ! Niveau de condensation de la CLA
       capCL(:,nsrf)  = 0.        ! CAPE de couche limite
       oliqCL(:,nsrf) = 0.        ! eau_liqu integree de couche limite
       cteiCL(:,nsrf) = 0.        ! cloud top instab. crit. couche limite
       pblt(:,nsrf)   = 0.        ! T a la Hauteur de couche limite
       therm(:,nsrf)  = 0.
       trmb1(:,nsrf)  = 0.        ! deep_cape
       trmb2(:,nsrf)  = 0.        ! inhibition 
       trmb3(:,nsrf)  = 0.        ! Point Omega

#undef T2m     
#define T2m     
#ifdef T2m
! Calculations of diagnostic t,q at 2m and u, v at 10m

       DO j=1, knon
          i = ni(j)
          uzon(j) = yu(j,1) + y_d_u(j,1)
          vmer(j) = yv(j,1) + y_d_v(j,1)
          tair1(j) = yt(j,1) + y_d_t(j,1)
          qair1(j) = yq(j,1) + y_d_q(j,1)
          zgeo1(j) = RD * tair1(j) / (0.5*(ypaprs(j,1)+ypplay(j,1))) &
               * (ypaprs(j,1)-ypplay(j,1))
          tairsol(j) = yts(j) + y_d_ts(j)
          rugo1(j) = yrugos(j)
          IF(nsrf.EQ.is_oce) THEN
             rugo1(j) = rugos(i,nsrf)
          ENDIF
          psfce(j)=ypaprs(j,1)
          patm(j)=ypplay(j,1)
          qairsol(j) = yqsurf(j)
       END DO
       

! Calculate the temperature et relative humidity at 2m and the wind at 10m 
       CALL stdlevvar(klon, knon, nsrf, zxli, &
            uzon, vmer, tair1, qair1, zgeo1, &
            tairsol, qairsol, rugo1, psfce, patm, &
            yt2m, yq2m, yt10m, yq10m, yu10m, yustar)

       DO j=1, knon
          i = ni(j)
          t2m(i,nsrf)=yt2m(j)
          q2m(i,nsrf)=yq2m(j)
          
          ! u10m, v10m : composantes du vent a 10m sans spirale de Ekman
          u10m(i,nsrf)=(yu10m(j) * uzon(j))/SQRT(uzon(j)**2+vmer(j)**2)
          v10m(i,nsrf)=(yu10m(j) * vmer(j))/SQRT(uzon(j)**2+vmer(j)**2)
       END DO

!IM Calcule de l'humidite relative a 2m (rh2m) pour diagnostique
!IM Ajoute dependance type surface
       IF (thermcep) THEN
          DO j = 1, knon
             i=ni(j)
             zdelta1 = MAX(0.,SIGN(1., rtt-yt2m(j) ))
             zx_qs1  = r2es * FOEEW(yt2m(j),zdelta1)/paprs(i,1)
             zx_qs1  = MIN(0.5,zx_qs1)
             zcor1   = 1./(1.-RETV*zx_qs1)
             zx_qs1  = zx_qs1*zcor1
             
             rh2m(i)   = rh2m(i)   + yq2m(j)/zx_qs1 * pctsrf(i,nsrf)
             qsat2m(i) = qsat2m(i) + zx_qs1  * pctsrf(i,nsrf)
          END DO
       END IF

       CALL HBTM(knon, ypaprs, ypplay, &
            yt2m,yt10m,yq2m,yq10m,yustar, &
            y_flux_t,y_flux_q,yu,yv,yt,yq, &
            ypblh,ycapCL,yoliqCL,ycteiCL,ypblT, &
            ytherm,ytrmb1,ytrmb2,ytrmb3,ylcl)
       
       DO j=1, knon
          i = ni(j)
          pblh(i,nsrf)   = ypblh(j)
          plcl(i,nsrf)   = ylcl(j)
          capCL(i,nsrf)  = ycapCL(j)
          oliqCL(i,nsrf) = yoliqCL(j)
          cteiCL(i,nsrf) = ycteiCL(j)
          pblT(i,nsrf)   = ypblT(j)
          therm(i,nsrf)  = ytherm(j)
          trmb1(i,nsrf)  = ytrmb1(j)
          trmb2(i,nsrf)  = ytrmb2(j)
          trmb3(i,nsrf)  = ytrmb3(j)
       END DO
       
#else 
! T2m not defined
! No calculation
       PRINT*,' Warning !!! No T2m calculation. Output is set to zero.'
#endif

!****************************************************************************************
! 15) End of loop over different surfaces
!
!****************************************************************************************
    END DO loop_nbsrf

!****************************************************************************************
! 16) Calculate the mean value over all sub-surfaces for som variables
!
!****************************************************************************************
    
    zxfluxt(:,:) = 0.0 ; zxfluxq(:,:) = 0.0
    zxfluxu(:,:) = 0.0 ; zxfluxv(:,:) = 0.0
    DO nsrf = 1, nbsrf
       DO k = 1, klev
          DO i = 1, klon
             zxfluxt(i,k) = zxfluxt(i,k) + flux_t(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxq(i,k) = zxfluxq(i,k) + flux_q(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxu(i,k) = zxfluxu(i,k) + flux_u(i,k,nsrf) * pctsrf(i,nsrf)
             zxfluxv(i,k) = zxfluxv(i,k) + flux_v(i,k,nsrf) * pctsrf(i,nsrf)
          END DO
       END DO
    END DO

    DO i = 1, klon
       zxsens(i)     = - zxfluxt(i,1) ! flux de chaleur sensible au sol
       zxevap(i)     = - zxfluxq(i,1) ! flux d'evaporation au sol
       fder_print(i) = fder(i) + dflux_t(i) + dflux_q(i)
    ENDDO
   
!
! Incrementer la temperature du sol
!
    zxtsol(:) = 0.0  ; zxfluxlat(:) = 0.0
    zt2m(:) = 0.0    ; zq2m(:) = 0.0 
    zu10m(:) = 0.0   ; zv10m(:) = 0.0
    s_pblh(:) = 0.0  ; s_plcl(:) = 0.0 
    s_capCL(:) = 0.0 ; s_oliqCL(:) = 0.0
    s_cteiCL(:) = 0.0; s_pblT(:) = 0.0
    s_therm(:) = 0.0 ; s_trmb1(:) = 0.0
    s_trmb2(:) = 0.0 ; s_trmb3(:) = 0.0
    
    
    DO nsrf = 1, nbsrf
       DO i = 1, klon          
          ts(i,nsrf) = ts(i,nsrf) + d_ts(i,nsrf)
          
          wfbils(i,nsrf) = ( solsw(i,nsrf) + sollw(i,nsrf) &
               + flux_t(i,1,nsrf) + fluxlat(i,nsrf) ) * pctsrf(i,nsrf)
          wfbilo(i,nsrf) = (evap(i,nsrf) - (rain_f(i) + snow_f(i))) * &
               pctsrf(i,nsrf)

          zxtsol(i)    = zxtsol(i)    + ts(i,nsrf)      * pctsrf(i,nsrf)
          zxfluxlat(i) = zxfluxlat(i) + fluxlat(i,nsrf) * pctsrf(i,nsrf)
          
          zt2m(i)  = zt2m(i)  + t2m(i,nsrf)  * pctsrf(i,nsrf)
          zq2m(i)  = zq2m(i)  + q2m(i,nsrf)  * pctsrf(i,nsrf)
          zu10m(i) = zu10m(i) + u10m(i,nsrf) * pctsrf(i,nsrf)
          zv10m(i) = zv10m(i) + v10m(i,nsrf) * pctsrf(i,nsrf)

          s_pblh(i)   = s_pblh(i)   + pblh(i,nsrf)  * pctsrf(i,nsrf)
          s_plcl(i)   = s_plcl(i)   + plcl(i,nsrf)  * pctsrf(i,nsrf)
          s_capCL(i)  = s_capCL(i)  + capCL(i,nsrf) * pctsrf(i,nsrf)
          s_oliqCL(i) = s_oliqCL(i) + oliqCL(i,nsrf)* pctsrf(i,nsrf)
          s_cteiCL(i) = s_cteiCL(i) + cteiCL(i,nsrf)* pctsrf(i,nsrf)
          s_pblT(i)   = s_pblT(i)   + pblT(i,nsrf)  * pctsrf(i,nsrf)
          s_therm(i)  = s_therm(i)  + therm(i,nsrf) * pctsrf(i,nsrf)
          s_trmb1(i)  = s_trmb1(i)  + trmb1(i,nsrf) * pctsrf(i,nsrf)
          s_trmb2(i)  = s_trmb2(i)  + trmb2(i,nsrf) * pctsrf(i,nsrf)
          s_trmb3(i)  = s_trmb3(i)  + trmb3(i,nsrf) * pctsrf(i,nsrf)
       END DO
    END DO

    IF (check) THEN
       amn=MIN(ts(1,is_ter),1000.)
       amx=MAX(ts(1,is_ter),-1000.)
       DO i=2, klon
          amn=MIN(ts(i,is_ter),amn)
          amx=MAX(ts(i,is_ter),amx)
       ENDDO
       PRINT*,' debut apres d_ts min max ftsol(ts)',itap,amn,amx
    ENDIF

!jg ?
!!$!
!!$! If a sub-surface does not exsist for a grid point, the mean value for all 
!!$! sub-surfaces is distributed.
!!$!
!!$    DO nsrf = 1, nbsrf
!!$       DO i = 1, klon
!!$          IF ((pctsrf_new(i,nsrf) .LT. epsfra) .OR. (t2m(i,nsrf).EQ.0.)) THEN
!!$             ts(i,nsrf)     = zxtsol(i)
!!$             t2m(i,nsrf)    = zt2m(i)
!!$             q2m(i,nsrf)    = zq2m(i)
!!$             u10m(i,nsrf)   = zu10m(i)
!!$             v10m(i,nsrf)   = zv10m(i)
!!$
!!$! Les variables qui suivent sont plus utilise, donc peut-etre pas la peine a les mettre ajour
!!$             pblh(i,nsrf)   = s_pblh(i)
!!$             plcl(i,nsrf)   = s_plcl(i)
!!$             capCL(i,nsrf)  = s_capCL(i)
!!$             oliqCL(i,nsrf) = s_oliqCL(i) 
!!$             cteiCL(i,nsrf) = s_cteiCL(i)
!!$             pblT(i,nsrf)   = s_pblT(i)
!!$             therm(i,nsrf)  = s_therm(i)
!!$             trmb1(i,nsrf)  = s_trmb1(i)
!!$             trmb2(i,nsrf)  = s_trmb2(i)
!!$             trmb3(i,nsrf)  = s_trmb3(i)
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO


    DO i = 1, klon
       fder(i) = - 4.0*RSIGMA*zxtsol(i)**3 
    ENDDO
    
    zxqsurf(:) = 0.0
    zxsnow(:)  = 0.0
    DO nsrf = 1, nbsrf
       DO i = 1, klon
          zxqsurf(i) = zxqsurf(i) + qsurf(i,nsrf) * pctsrf(i,nsrf)
          zxsnow(i)  = zxsnow(i)  + snow(i,nsrf)  * pctsrf(i,nsrf)
       END DO
    END DO

! Premier niveau de vent sortie dans physiq.F
    zu1(:) = u(:,1)
    zv1(:) = v(:,1)

! Some of the module declared variables are returned for printing in physiq.F
    qsol_d(:)     = qsol(:)
    evap_d(:,:)   = evap(:,:)
    rugos_d(:,:)  = rugos(:,:) 
    agesno_d(:,:) = agesno(:,:)


  END SUBROUTINE pbl_surface
!
!****************************************************************************************
!
  SUBROUTINE pbl_surface_final(qsol_rst, fder_rst, snow_rst, qsurf_rst, &
       evap_rst, rugos_rst, agesno_rst, ftsoil_rst)

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"

! Ouput variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)                 :: qsol_rst
    REAL, DIMENSION(klon), INTENT(OUT)                 :: fder_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: snow_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: qsurf_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: evap_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: rugos_rst
    REAL, DIMENSION(klon, nbsrf), INTENT(OUT)          :: agesno_rst
    REAL, DIMENSION(klon, nsoilmx, nbsrf), INTENT(OUT) :: ftsoil_rst

 
!****************************************************************************************
! Return module variables for writing to restart file
!
!****************************************************************************************    
    qsol_rst(:)       = qsol(:)
    fder_rst(:)       = fder(:)
    snow_rst(:,:)     = snow(:,:)
    qsurf_rst(:,:)    = qsurf(:,:)
    evap_rst(:,:)     = evap(:,:)
    rugos_rst(:,:)    = rugos(:,:)
    agesno_rst(:,:)   = agesno(:,:)
    ftsoil_rst(:,:,:) = ftsoil(:,:,:)

!****************************************************************************************
! Deallocate module variables
!
!****************************************************************************************
!   DEALLOCATE(qsol, fder, snow, qsurf, evap, rugos, agesno, ftsoil)
    IF (ALLOCATED(qsol)) DEALLOCATE(qsol)
    IF (ALLOCATED(fder)) DEALLOCATE(fder)
    IF (ALLOCATED(snow)) DEALLOCATE(snow)
    IF (ALLOCATED(qsurf)) DEALLOCATE(qsurf)
    IF (ALLOCATED(evap)) DEALLOCATE(evap)
    IF (ALLOCATED(rugos)) DEALLOCATE(rugos)
    IF (ALLOCATED(agesno)) DEALLOCATE(agesno)
    IF (ALLOCATED(ftsoil)) DEALLOCATE(ftsoil)

  END SUBROUTINE pbl_surface_final
!  
!****************************************************************************************
! 
  SUBROUTINE pbl_surface_newfrac(itime, pctsrf_new, pctsrf_old, tsurf, alb1, alb2, u10m, v10m, tke)

    ! Give default values where new fraction has appread

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "clesphys.h"
    INCLUDE "compbl.h"

! Input variables
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf_new, pctsrf_old

! InOutput variables
!****************************************************************************************
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: tsurf
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: alb1, alb2
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT)        :: u10m, v10m
    REAL, DIMENSION(klon,klev+1,nbsrf), INTENT(INOUT) :: tke

! Local variables
!****************************************************************************************
    INTEGER           :: nsrf, nsrf_comp1, nsrf_comp2, nsrf_comp3, i
    CHARACTER(len=80) :: abort_message
    CHARACTER(len=20) :: modname = 'pbl_surface_newfrac'
    INTEGER, DIMENSION(nbsrf) :: nfois=0, mfois=0, pfois=0
!
! All at once !! 
!****************************************************************************************
    
    DO nsrf = 1, nbsrf
       ! First decide complement sub-surfaces
       SELECT CASE (nsrf)
       CASE(is_oce)
          nsrf_comp1=is_sic
          nsrf_comp2=is_ter
          nsrf_comp3=is_lic
       CASE(is_sic)
          nsrf_comp1=is_oce
          nsrf_comp2=is_ter
          nsrf_comp3=is_lic
       CASE(is_ter)
          nsrf_comp1=is_lic
          nsrf_comp2=is_oce
          nsrf_comp3=is_sic
       CASE(is_lic)
          nsrf_comp1=is_ter
          nsrf_comp2=is_oce
          nsrf_comp3=is_sic
       END SELECT

       ! Initialize all new fractions
       DO i=1, klon
          IF (pctsrf_new(i,nsrf) > 0. .AND. pctsrf_old(i,nsrf) == 0.) THEN
             
             IF (pctsrf_old(i,nsrf_comp1) > 0.) THEN
                ! Use the complement sub-surface, keeping the continents unchanged
                qsurf(i,nsrf) = qsurf(i,nsrf_comp1)
                evap(i,nsrf)  = evap(i,nsrf_comp1)
                rugos(i,nsrf) = rugos(i,nsrf_comp1)
                tsurf(i,nsrf) = tsurf(i,nsrf_comp1)
                alb1(i,nsrf)  = alb1(i,nsrf_comp1)
                alb2(i,nsrf)  = alb2(i,nsrf_comp1)
                u10m(i,nsrf)  = u10m(i,nsrf_comp1)
                v10m(i,nsrf)  = v10m(i,nsrf_comp1)
                if (iflag_pbl > 1) then
                 tke(i,:,nsrf) = tke(i,:,nsrf_comp1)
                endif
                mfois(nsrf) = mfois(nsrf) + 1
             ELSE
                ! The continents have changed. The new fraction receives the mean sum of the existent fractions
                qsurf(i,nsrf) = qsurf(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + qsurf(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                evap(i,nsrf)  = evap(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + evap(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                rugos(i,nsrf) = rugos(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + rugos(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                tsurf(i,nsrf) = tsurf(i,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + tsurf(i,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                alb1(i,nsrf)  = alb1(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + alb1(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                alb2(i,nsrf)  = alb2(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + alb2(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                u10m(i,nsrf)  = u10m(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + u10m(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                v10m(i,nsrf)  = v10m(i,nsrf_comp2) *pctsrf_old(i,nsrf_comp2) + v10m(i,nsrf_comp3) *pctsrf_old(i,nsrf_comp3)
                if (iflag_pbl > 1) then
                 tke(i,:,nsrf) = tke(i,:,nsrf_comp2)*pctsrf_old(i,nsrf_comp2) + tke(i,:,nsrf_comp3)*pctsrf_old(i,nsrf_comp3)
                endif
            
                ! Security abort. This option has never been tested. To test, comment the following line.
!                abort_message='The fraction of the continents have changed!'
!                CALL abort_gcm(modname,abort_message,1)
                nfois(nsrf) = nfois(nsrf) + 1
             END IF
             snow(i,nsrf)     = 0.
             agesno(i,nsrf)   = 0.
             ftsoil(i,:,nsrf) = tsurf(i,nsrf)
          ELSE
             pfois(nsrf) = pfois(nsrf)+ 1
          END IF
       END DO
       
    END DO

  END SUBROUTINE pbl_surface_newfrac

!  
!****************************************************************************************
!  

END MODULE pbl_surface_mod
