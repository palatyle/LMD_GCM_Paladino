!*******************************************************************************
! PURPOSE: LMD_driver is the WRF mediation layer routine that provides the
! interface to LMD physics packages in the WRF model layer. For those familiar
! with the LMD GCM, the aim of this driver is to do part of the job of the
! calfis.F routine.
!*******************************************************************************
! AUTHOR: A. Spiga - January 2007
!*******************************************************************************
! UPDATES:  - included all final updates for the paper - March 2008
!           - general cleaning of code and comments - October 2008 
!           - additions for idealized cases - January 2009
!           - additions for new soil model in physics - January 2010 
!           - unified module_lmd_driver: old, new phys and LES - February 2011
!*******************************************************************************
MODULE module_lmd_driver
CONTAINS
    SUBROUTINE lmd_driver(id,max_dom,DT,ITIMESTEP,XLAT,XLONG, &
        IDS,IDE,JDS,JDE,KDS,KDE,IMS,IME,JMS,JME,KMS,KME, &
        i_start,i_end,j_start,j_end,kts,kte,num_tiles, &
        DX,DY, &
        MSFT,MSFU,MSFV, &
        GMT,JULYR,JULDAY, &
        P8W,DZ8W,T8W,Z,HT, &
        U,V,W,TH,T,P,EXNER,RHO, &
        PTOP, &
        RADT,CUDT, &
        TSK,PSFC, &
        RTHBLTEN,RUBLTEN,RVBLTEN, &
        num_3d_s,SCALAR, &
        MARS_MODE, &
        M_ALBEDO,M_TI,M_CO2ICE,M_EMISS, &
        M_H2OICE, &
        M_TSOIL, &
        M_Q2, &
        M_TSURF, &
#ifdef NEWPHYS
        M_FLUXRAD, &
        M_WSTAR, &
        M_ISOIL, &
        M_DSOIL, &
        M_Z0, &
        CST_Z0, &
#endif
        M_GW, &
        NUM_SOIL_LAYERS, &
        CST_AL, &
        CST_TI, &
        isfflx, &
        diff_opt, &
        km_opt, &
        HISTORY_INTERVAL, &
#ifndef NOPHYS
#include "module_lmd_driver_output1.inc"
#endif
        SLPX,SLPY,RESTART)
! NB: module_lmd_driver_output1.inc : output arguments generated from Registry





!==================================================================
! USES
!==================================================================
   USE module_model_constants
   USE module_wrf_error
! add new modules here, if needed ...

!==================================================================
  IMPLICIT NONE
!==================================================================

!==================================================================
! COMMON
!==================================================================

#ifndef NOPHYS
!! the only common needed is the one defining the physical grid
!! -- path is hardcoded, but the structure is not subject to change
!! -- please put # if needed by the pre-compilation process	
!
!
include "../mars_lmd/libf/phymars/dimphys.h"
!
!--to be commented because there are tests in the physics ?
!--TODO: get rid of the ...mx first in this routine and .inc

!
! INCLUDE AUTOMATIQUEMENT GENERE A PARTIR DU REGISTRY
!
include "../mars_lmd/libf/phymars/wrf_output_2d.h"
include "../mars_lmd/libf/phymars/wrf_output_3d.h"
#endif

!==================================================================
! VARIABLES
!==================================================================

! WRF Dimensions
INTEGER, INTENT(IN   )    ::   &
      ids,ide,jds,jde,kds,kde, &
      ims,ime,jms,jme,kms,kme, &
            kts,kte,num_tiles, &
              NUM_SOIL_LAYERS
INTEGER, DIMENSION(num_tiles), INTENT(IN) ::  &
      i_start,i_end,j_start,j_end   
! Scalars
INTEGER, INTENT(IN  ) :: JULDAY, itimestep, julyr,id,max_dom,MARS_MODE
INTEGER, INTENT(IN  ) :: isfflx,diff_opt,km_opt
REAL, INTENT(IN  ) :: GMT,dt,dx,dy,RADT,CUDT
REAL, INTENT(IN  ) :: CST_AL, CST_TI 
REAL, INTENT(IN  ) :: PTOP
! 2D arrays          
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   )  :: &
     MSFT,MSFU,MSFV, &
     XLAT,XLONG,HT,  &
     M_ALBEDO,M_TI,M_EMISS, &
     SLPX,SLPY
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT   )  :: &
     M_CO2ICE,M_H2OICE, &
     M_TSURF
! 3D arrays 
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: &
     dz8w,p8w,p,exner,t,t8w,rho,u,v,w,z,th
!REAL, DIMENSION( ims:ime, kms:kme+1, jms:jme ), INTENT(INOUT ) :: &
!     M_Q2
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT ) :: &
     M_Q2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, INTENT(IN   ) :: HISTORY_INTERVAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL, DIMENSION( ims:ime, NUM_SOIL_LAYERS, jms:jme ), INTENT(INOUT )  :: &
     M_TSOIL 
#ifdef NEWPHYS
REAL, INTENT(IN  ) :: CST_Z0
REAL, DIMENSION( ims:ime, NUM_SOIL_LAYERS, jms:jme ), INTENT(IN   )  :: &
     M_ISOIL, M_DSOIL         
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   )  :: &
     M_Z0
REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT   )  :: &
     M_FLUXRAD,M_WSTAR
#endif
REAL, DIMENSION( ims:ime, 5, jms:jme ), INTENT(IN   )  :: &
     M_GW
! 4D arrays
INTEGER, INTENT(IN ) :: num_3d_s
REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_3d_s ), INTENT(INOUT ) :: &
     scalar
! Logical
LOGICAL, INTENT(IN ) :: restart

!-------------------------------------------
! OUTPUT VARIABLES
!-------------------------------------------
!
! Generated from Registry
!
! default definitions :
! 2D : TSK, PSFC
! 3D : RTHBLTEN,RUBLTEN,RVBLTEN
#ifndef NOPHYS
#include "module_lmd_driver_output2.inc"
   REAL, DIMENSION(:,:), ALLOCATABLE :: output_tab2d
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: output_tab3d
#else
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: PSFC,TSK
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT)  :: RTHBLTEN,RUBLTEN,RVBLTEN
   REAL,DIMENSION(:),ALLOCATABLE,SAVE :: dtrad
#endif
!-------------------------------------------
! OUTPUT VARIABLES
!-------------------------------------------


!-------------------------------------------   
! LOCAL  VARIABLES
!-------------------------------------------
   INTEGER ::    i,j,k,its,ite,jts,jte,ij,kk
   INTEGER ::    subs,iii

   ! *** for Mars Mode 20 and 21 ***
   REAL ::    tau_decay, lct, zilctmin
   INTEGER, SAVE :: zipbl(500) !index of zi in file input_zipbl
   REAL, SAVE :: zilct(500) !corresponding local time in input_zipbl
   INTEGER :: lctindex,ziindex
   LOGICAL :: end_of_file
   ! *** ----------------------- ***

   ! *** for LMD physics
   ! ------> inputs:
   INTEGER :: ngrid,nlayer,nq,nsoil 
   REAL :: pday,ptime,MY  
   REAL :: aire_val,lat_val,lon_val
   REAL :: phisfi_val,albedodat_val,inertiedat_val
   REAL :: tsurf_val,co2ice_val,emis_val
   REAL :: zmea_val,zstd_val,zsig_val,zgam_val,zthe_val
   REAL :: theta_val, psi_val
   LOGICAL :: firstcall,lastcall,tracerdyn
   REAL,DIMENSION(:),ALLOCATABLE :: q2_val, qsurf_val, tsoil_val
#ifdef NEWPHYS
   REAL :: wstar_val,fluxrad_val
   REAL,DIMENSION(:),ALLOCATABLE :: isoil_val, dsoil_val
   REAL :: z0_val
#endif
   REAL,DIMENSION(:),ALLOCATABLE :: aire_vec,lat_vec,lon_vec
   REAL,DIMENSION(:),ALLOCATABLE :: walbedodat,winertiedat,wphisfi
   REAL,DIMENSION(:),ALLOCATABLE :: wzmea,wzstd,wzsig,wzgam,wzthe
   REAL,DIMENSION(:),ALLOCATABLE :: wtheta, wpsi
   ! v--- can they be modified ?
   REAL,DIMENSION(:),ALLOCATABLE :: wtsurf,wco2ice,wemis
   REAL,DIMENSION(:,:),ALLOCATABLE :: wq2,wqsurf,wtsoil 
#ifdef NEWPHYS
   REAL,DIMENSION(:),ALLOCATABLE :: wwstar,wfluxrad
   REAL,DIMENSION(:),ALLOCATABLE :: wz0tab
   REAL,DIMENSION(:,:),ALLOCATABLE :: wisoil,wdsoil
   CHARACTER*20,DIMENSION(:),ALLOCATABLE :: wtnom
#endif 
   ! ----------	
   REAL,DIMENSION(:,:),ALLOCATABLE :: pplev,pplay,pphi,pu,pv,pt,pw
   REAL,DIMENSION(:,:,:),ALLOCATABLE :: pq  

   ! <------ outputs:
   !     physical tendencies
   REAL,DIMENSION(:),ALLOCATABLE :: pdpsrf
   REAL,DIMENSION(:,:),ALLOCATABLE :: pdu,pdv,pdt
   REAL,DIMENSION(:,:,:),ALLOCATABLE :: pdq
   ! ... intermediate arrays
   REAL, DIMENSION(:), ALLOCATABLE  :: &
     dz8w_prof,p8w_prof,p_prof,t_prof,t8w_prof, &
     u_prof,v_prof,z_prof
!!   pi_prof, rho_prof, th_prof, &
#ifdef NEWPHYS
   REAL, DIMENSION(:,:), ALLOCATABLE :: q_prof
#else
   REAL, DIMENSION(:), ALLOCATABLE  :: &
     water_vapor_prof, water_ice_prof
#endif


   ! Additional control variables
   INTEGER :: sponge_top,relax,ips,ipe,jps,jpe,kps,kpe
   REAL :: elaps, ptimestep
   INTEGER :: wday_ini, test, test2
   REAL :: wappel_phys
   LOGICAL :: flag_LES
   LOGICAL, SAVE :: flag_first_restart
!**************************************************
! IMPORTANT: pour travailler avec grid nesting
! --- deux comportements distincts du save
! ... ex1: ferme planeto avec PGF+MPI: standard
! ... ex2: iDataPlex avec IFORT+MPI: SPECIAL_NEST_SAVE 
!**************************************************
#ifdef SPECIAL_NEST_SAVE
      !  une dimension supplementaire liee au nest
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             dp_save
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
             du_save, dv_save, dt_save
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
             dq_save      
#ifndef NORESTART
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
             save_tsoil_restart
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_tsurf_restart
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_co2ice_restart
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
             save_q2_restart
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
             save_qsurf_restart
#ifdef NEWPHYS
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_wstar_restart
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_fluxrad_restart
#endif
#endif
#else
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: &
             dp_save
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             du_save, dv_save, dt_save
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
             dq_save      
#ifndef NORESTART
!! FOR RESTART
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_tsoil_restart
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: &
             save_tsurf_restart
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: &
             save_co2ice_restart
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_q2_restart
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
             save_qsurf_restart
#ifdef NEWPHYS
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: &
             save_wstar_restart
      REAL, DIMENSION(:), ALLOCATABLE, SAVE :: &
             save_fluxrad_restart
#endif
#endif
#endif

!!!IDEALIZED IDEALIZED
      REAL :: lat_input, lon_input, ls_input, lct_input
!!!IDEALIZED IDEALIZED

!!!!!!!!!!!!!! TEST NaN or Inf : ne fonction qu'avec r4 i4
integer :: mask_nan = 2139095040
integer :: chuck_norris
!!!!!!!!!!!!!! TEST NaN or Inf : ne fonction qu'avec r4 i4

!REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
!     UMAX, UMIN, VMAX, VMIN, WMAX, WMIN, TMAX, TMIN



!==================================================================
! CODE
!==================================================================

IF (JULYR .ne. 9999) THEN
    flag_LES = .false.  ! "True" LES is not available in this version
    PRINT *, '*** REAL-CASE SIMULATION ***'
ELSE
     PRINT *, '*** IDEALIZED SIMULATION ***'
     IF ((diff_opt .eq. 2) .and. (km_opt .eq. 2)) THEN
          PRINT *, '*** type: LES ***'
          PRINT *, '*** diff_opt = 2 *** km_opt = 2'
          PRINT *, '*** forcing is isfflx = ',isfflx 
          flag_LES = .true.
          !! SPECIAL LES
     ELSE
          PRINT *, '*** type: not LES ***'
          PRINT *, '*** diff_opt = ',diff_opt
          PRINT *, '*** km_opt = ',km_opt
          flag_LES = .false.
          !! IDEALIZED, no LES
          !! cependant, ne veut-on pas pouvoir
          !! prescrire un flux en idealise ??
     ENDIF
ENDIF


print *,'** Mars ** DOMAIN',id

!-------------------------!
! TWEAK on WRF DIMENSIONS !
!-------------------------!
its = i_start(1)   ! define here one big tile
ite = i_end(num_tiles)
jts = j_start(1)
jte = j_end(num_tiles)
!!
IF (flag_LES .eqv. .false.) THEN
 relax=0
 sponge_top=0               ! another value than 0 triggers instabilities  
 IF (id .gt. 1) relax=2     ! fix to avoid noise in nesting simulations ; 1 >> too much noise ...
ENDIF
ips=its
ipe=ite
jps=jts
jpe=jte
IF (flag_LES .eqv. .false.) THEN
 IF (ips .eq. ids)   ips=its+relax !! IF tests necesary for parallel runs
 IF (ipe .eq. ide-1) ipe=ite-relax
 IF (jps .eq. jds)   jps=jts+relax
 IF (jpe .eq. jde-1) jpe=jte-relax
ENDIF
kps=kts         !! start at surface 
IF (flag_LES .eqv. .false.) THEN
 kpe=kte-sponge_top
ELSE
 PRINT *, '*** IDEALIZED SIMULATION: LES *** kpe=kte'
 kpe=kte !-sponge_top
ENDIF

!----------------------------!
! DIMENSIONS FOR THE PHYSICS !
!----------------------------!
wday_ini = JULDAY - 1      !! GCM convention
wappel_phys = RADT
ptimestep = dt*wappel_phys            ! physical timestep (s)
ngrid=(ipe-ips+1)*(jpe-jps+1)         ! size of the horizontal grid
nlayer = kpe-kps+1                    ! number of vertical layers: nlayermx
nsoil = NUM_SOIL_LAYERS               ! number of soil layers: nsoilmx
if (num_3d_s > 1) then                ! number of advected fields
        nq = num_3d_s-1               
else
        nq = 1
endif
! **** needed but hardcoded
lastcall = .false.
! **** needed but hardcoded

          !PRINT *, ips, ipe, jps, jpe
          !PRINT *, ngrid

elaps = (float(itimestep)-1.)*dt  ! elapsed seconds of simulation
!----------------------------------------------!
! what is done at the first step of simulation !
!----------------------------------------------!
IF (elaps .eq. 0.) THEN
   flag_first_restart = .false.
ELSE
   flag_first_restart=flag_first_restart.OR.(.NOT. ALLOCATED(dp_save))
ENDIF

IF ((elaps .eq. 0.).or.flag_first_restart) THEN
firstcall=.true.  !! for continuity with GCM, physics are always called at the first WRF timestep
  !firstcall=.false. !! just in case you'd want to get rid of the physics 
test=0
#ifdef SPECIAL_NEST_SAVE
PRINT *, 'ALLOCATED dp_save ???', ALLOCATED( dp_save ), id
IF( .NOT. ALLOCATED( dp_save  ) ) THEN
   PRINT *, '**** check **** OK I ALLOCATE one save ARRAY for all NESTS ', max_dom, id
   !! here are the arrays that would be useful to save physics tendencies
   ALLOCATE(dp_save(ngrid,max_dom))
   ALLOCATE(du_save(ngrid,nlayer,max_dom))
   ALLOCATE(dv_save(ngrid,nlayer,max_dom))
   ALLOCATE(dt_save(ngrid,nlayer,max_dom))
   ALLOCATE(dq_save(ngrid,nlayer,nq,max_dom))
   dp_save(:,:)=0.    !! initialize these arrays ...
   du_save(:,:,:)=0.
   dv_save(:,:,:)=0.
   dt_save(:,:,:)=0.
   dq_save(:,:,:,:)=0.
#ifndef NORESTART
   ! Restart save arrays
   ALLOCATE(save_tsoil_restart(ngrid,nsoil,max_dom))
   ALLOCATE(save_co2ice_restart(ngrid,max_dom))
   ALLOCATE(save_q2_restart(ngrid,nlayer+1,max_dom))
   ALLOCATE(save_qsurf_restart(ngrid,nq,max_dom))
   ALLOCATE(save_tsurf_restart(ngrid,max_dom))
   save_tsoil_restart(:,:,:)=0.
   save_co2ice_restart(:,:)=0.
   save_q2_restart(:,:,:)=0.
   save_qsurf_restart(:,:,:)=0.
   save_tsurf_restart(:,:)=0.
#ifdef NEWPHYS
   ALLOCATE(save_wstar_restart(ngrid,max_dom))
   ALLOCATE(save_fluxrad_restart(ngrid,max_dom))
   save_wstar_restart(:,:)=0.
   save_fluxrad_restart(:,:)=0.
#endif
#endif
ENDIF
IF (id .lt. max_dom) THEN
   flag_first_restart=.true.
ELSE
   flag_first_restart=.false.
ENDIF
#else
IF( .NOT. ALLOCATED( dp_save  ) ) THEN
ALLOCATE(dp_save(ngrid)) !! here are the arrays that would be useful to save physics tendencies
ALLOCATE(du_save(ngrid,nlayer))
ALLOCATE(dv_save(ngrid,nlayer))
ALLOCATE(dt_save(ngrid,nlayer))
ALLOCATE(dq_save(ngrid,nlayer,nq))
ENDIF
dp_save(:)=0.    !! initialize these arrays ...
du_save(:,:)=0.
dv_save(:,:)=0.
dt_save(:,:)=0.
dq_save(:,:,:)=0.
flag_first_restart=.false.
#ifndef NORESTART
! Restart save arrays
IF( .NOT. ALLOCATED( save_tsoil_restart  ) ) THEN
ALLOCATE(save_tsoil_restart(ngrid,nsoil))
ALLOCATE(save_co2ice_restart(ngrid))
ALLOCATE(save_q2_restart(ngrid,nlayer+1))
ALLOCATE(save_qsurf_restart(ngrid,nq))
ALLOCATE(save_tsurf_restart(ngrid))
ENDIF
save_tsoil_restart(:,:)=0.
save_co2ice_restart(:)=0.
save_q2_restart(:,:)=0.
save_qsurf_restart(:,:)=0.
save_tsurf_restart(:)=0.
#ifdef NEWPHYS
IF( .NOT. ALLOCATED( save_wstar_restart  ) ) THEN
ALLOCATE(save_wstar_restart(ngrid))
ALLOCATE(save_fluxrad_restart(ngrid))
ENDIF
save_wstar_restart(:)=0.
save_fluxrad_restart(:)=0.
#endif
#endif
#endif

!! put here some general information you'd like to print just once
    print *, 'TILES: ', i_start,i_end, j_start, j_end  ! numbers for simple runs, arrays for parallel runs
    print *, 'DOMAIN: ', ids, ide, jds, jde
    print *, 'MEMORY: ', ims, ime, jms, jme
    print *, 'ADVECTED TRACERS: ', num_3d_s-1
    print *, 'PHYSICS IS CALLED EACH...',wappel_phys
!! put here some general information you'd like to print just once
ELSE
!-------------------------------------------------!
! what is done for the other steps of simulation  !
!-------------------------------------------------!
IF (wappel_phys .gt. 0.) THEN
firstcall = .false.
test = MODULO(itimestep-1,int(wappel_phys)) ! WRF time is integrated time (itimestep=1 at the end of first step) 
                                            ! -- same strategy as diagfi in the LMD GCM 
                                            ! -- arguments of modulo must be of the same kind (here: integers)
ELSE
firstcall = .false.
test = 9999   
ENDIF
ENDIF 

!!!!! for 'subgrid' temporal diagnostics
!test2 = MODULO(elaps,history_interval*100.)

!!!******!!
!!! TIME !!
!!!******!!
IF (JULYR .ne. 9999) THEN
  !
  ! specified
  !
  ptime = (GMT + elaps/3700.) !! universal time (0<ptime<1): ptime=0.5 at 12:00 UT 
  ptime = MODULO(ptime,24.)   !! the two arguments of MODULO must be of the same type 
  ptime = ptime / 24.
  pday = (JULDAY - 1 + INT((3700*GMT+elaps)/88800))
  pday = MODULO(int(pday),669)
  MY = (JULYR-2000) + (88800.*(JULDAY - 1)+3700.*GMT+elaps)/59496000.
  MY = INT(MY)
ELSE
  !
  ! idealized
  !
  PRINT *,'** Mars ** IDEALIZED SIMULATION'
  PRINT *,'** Mars ** BEWARE: input_coord must be here'
  open(unit=14,file='input_coord',form='formatted',status='old')
  rewind(14)
  read(14,*) lon_input
  read(14,*) lat_input
  read(14,*) ls_input
  read(14,*) lct_input
  close(14)
  ptime = lct_input - lon_input / 15. + elaps/3700.
  ptime = MODULO(ptime,24.)
  ptime = ptime / 24.
  pday = floor(ls2sol(ls_input)) + INT((3700*(lct_input - lon_input / 15.) + elaps)/88800)
  pday = MODULO(int(pday),669)
  MY = 2024
  wday_ini = floor(ls2sol(ls_input))
ENDIF
print *,'** Mars ** TIME IS', pday, ptime*24.


!!****************!!
!! CHECK DYNAMICS !!
!!****************!!
!IF ((MAXVAL(t) > 500).OR.(MINVAL(t,MASK = t > 0) <= 50)) THEN
!        PRINT *,'******************   CRASH   *******************'
!        PRINT *,'Irrealistic temperature...', MAXLOC(t), MINLOC(t)
!PRINT *, t
!        PRINT *,'************************************************'
!        STOP
!ENDIF   

!print *, 'check dynamics'
!  !!! in some cases, weird values are displayed
!  !!! despite the fact that outputs are OK...
!  !print *, 'u', MAXVAL(u), MINVAL(u)
!  !print *, 'v', MAXVAL(v), MINVAL(v)
!  !print *, 'w', MAXVAL(w), MINVAL(w)
!  !print *, 't', MAXVAL(t), MINVAL(t, MASK = t > 0)
!print *, 'u', u(10,1,10), u(10,15,10)
!print *, 'v', v(10,1,10), v(10,15,10)
!print *, 'w', w(10,1,10), w(10,15,10) 
!print *, 't', t(10,1,10), t(10,15,10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!IF (test2.EQ.0) THEN
!  print *, 'compute stats'
!  print *, 'RESET'
!   uave  = uave*0.
!     vave  = vave*0.
!     tave  = tave*0.
!     wave  = wave*0.
!   ustd  = ustd*0.
!     vstd  = vstd*0.
!     tstd  = tstd*0.
!     wstd  = wstd*0.
!ENDIF
!  uave = uave + u  * dt / (float(history_interval)*100.)
!  vave = vave + v  * dt / (float(history_interval)*100.)
!  tave = tave + th * dt / (float(history_interval)*100.)
!  wave = wave + w  * dt / (float(history_interval)*100.)
!  ustd = ustd + u  * u  * dt / (float(history_interval)*100.)
!  vstd = vstd + v  * v  * dt / (float(history_interval)*100.)
!  tstd = tstd + th * th * dt / (float(history_interval)*100.)
!  wstd = wstd + w  * w  * dt / (float(history_interval)*100.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!
!!!!! PROBLEME avec 3 NESTS !!!
!!!!!
!chuck_norris = transfer( u(1,1,1),chuck_norris )        !! astuce J. Lefrere - test NaN or Inf 
!PRINT *,'check stability' 
!IF ( iand( chuck_norris,mask_nan ) == mask_nan ) THEN   !! astuce J. Lefrere - test NaN or Inf
!        PRINT *, u(1,1,1)
!        PRINT *,'******************   CRASH   *******************'
!        PRINT *,'************************************************'
!        PRINT *,'NaN or Inf appeared in the simulation ...'
!        PRINT *,'...this may be due to numerical or dynamical instability'
!        PRINT *,'************************************************'
!        PRINT *,'POSSIBLE SOLUTIONS:'
!        PRINT *,'>> IF nonhydrostatic mode,' 
!        PRINT *,'   --> check that smdiv, emdiv and epssm are not 0.'
!        PRINT *,'>> IF cfl is violated, '
!        PRINT *,'   --> try to lower the dynamical timestep'
!        PRINT *,'>> IF topographical gradients are high near specified bdy,'
!        PRINT *,'   --> try to redefine the domain'
!        PRINT *,'************************************************'
!        STOP
!ELSE
!        PRINT *,'OK OK OK OK'
!ENDIF
!ENDIF
        !IF (          ANY(isNaN(u)) &
        !         .OR. ANY(isNaN(v)) &
        !         .OR. ANY(isNaN(t)) ) THEN
        ! >>> ne marche qu'avec g95
!print *, 'check dynamics'
!print *, 'u', MAXVAL(u), MINVAL(u)
!print *, 'v', MAXVAL(v), MINVAL(v)
!print *, 't', MAXVAL(t), MINVAL(t, MASK = t > 0) 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------!
! ALLOCATE !
!----------!
ALLOCATE(pdpsrf(ngrid))
ALLOCATE(pdu(ngrid,nlayer))
ALLOCATE(pdv(ngrid,nlayer))
ALLOCATE(pdt(ngrid,nlayer))
ALLOCATE(pdq(ngrid,nlayer,nq))
!!!
!!! BIG LOOP : 1. no call for physics, used saved values
!!!
IF (test.NE.0) THEN
print *,'** Mars ** NO CALL FOR PHYSICS, go to next step...',test
#ifdef SPECIAL_NEST_SAVE
pdpsrf(:)=dp_save(:,id)
pdu(:,:)=du_save(:,:,id)
pdv(:,:)=dv_save(:,:,id)
pdt(:,:)=dt_save(:,:,id)
pdq(:,:,:)=dq_save(:,:,:,id)
#else
pdpsrf(:)=dp_save(:)
pdu(:,:)=du_save(:,:)
pdv(:,:)=dv_save(:,:)
pdt(:,:)=dt_save(:,:)
pdq(:,:,:)=dq_save(:,:,:)
#endif
!!!
!!! BIG LOOP : 2. calculate physical tendencies
!!!
ELSE
!----------!
! ALLOCATE !
!----------!
! inputs ...
IF (firstcall .EQV. .true.) THEN
  ALLOCATE(aire_vec(ngrid))!
  ALLOCATE(lon_vec(ngrid))!
  ALLOCATE(lat_vec(ngrid))!
  ALLOCATE(walbedodat(ngrid))!
  ALLOCATE(winertiedat(ngrid))!
  ALLOCATE(wphisfi(ngrid))!
  ALLOCATE(wzmea(ngrid))!
  ALLOCATE(wzstd(ngrid))!
  ALLOCATE(wzsig(ngrid))!
  ALLOCATE(wzgam(ngrid))!
  ALLOCATE(wzthe(ngrid))!
  ALLOCATE(wtheta(ngrid))!
  ALLOCATE(wpsi(ngrid))!
#ifdef NEWPHYS
  ALLOCATE(wz0tab(ngrid))!
#endif
ENDIF
ALLOCATE(wtsurf(ngrid))          !!!!!
#ifndef NOPHYS
ALLOCATE(output_tab2d(ngrid,n2d))
ALLOCATE(output_tab3d(ngrid,nlayer,n3d))
#endif
ALLOCATE(wco2ice(ngrid))         !!!!!
ALLOCATE(wemis(ngrid))           !!!!!
ALLOCATE(q2_val(nlayer+1))
ALLOCATE(qsurf_val(nq))
ALLOCATE(tsoil_val(nsoil))
#ifdef NEWPHYS
ALLOCATE(isoil_val(nsoil))
ALLOCATE(dsoil_val(nsoil))
#endif
ALLOCATE(wq2(ngrid,nlayer+1))    !!!!!
ALLOCATE(wqsurf(ngrid,nq))       !!!!!
ALLOCATE(wtsoil(ngrid,nsoil))    !!!!!
#ifdef NEWPHYS
ALLOCATE(wfluxrad(ngrid))
ALLOCATE(wwstar(ngrid))
ALLOCATE(wisoil(ngrid,nsoil))    !!!!!
ALLOCATE(wdsoil(ngrid,nsoil))    !!!!!
#endif
ALLOCATE(pplev(ngrid,nlayer+1))  !!!!!
ALLOCATE(pplay(ngrid,nlayer))    !!!!!
ALLOCATE(pphi(ngrid,nlayer))     !!!!!
ALLOCATE(pu(ngrid,nlayer))       !!!!!
ALLOCATE(pv(ngrid,nlayer))       !!!!!
ALLOCATE(pt(ngrid,nlayer))       !!!!!
ALLOCATE(pw(ngrid,nlayer))       !!!!!
ALLOCATE(pq(ngrid,nlayer,nq))    !!!!!
! interm
ALLOCATE(dz8w_prof(nlayer))
ALLOCATE(p8w_prof(nlayer))
ALLOCATE(p_prof(nlayer))
ALLOCATE(t_prof(nlayer))
ALLOCATE(t8w_prof(nlayer))
ALLOCATE(u_prof(nlayer))
ALLOCATE(v_prof(nlayer))
ALLOCATE(z_prof(nlayer))
!ALLOCATE(th_prof(nlayer))
!ALLOCATE(rho_prof(nlayer))
!ALLOCATE(pi_prof(nlayer))
#ifdef NEWPHYS
ALLOCATE(q_prof(nlayer,nq))
ALLOCATE(wtnom(nq))              !!!!!
#else
ALLOCATE(water_vapor_prof(nlayer))
ALLOCATE(water_ice_prof(nlayer))
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PREPARE PHYSICS INPUTS !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!SELECT CASE (MARS_MODE) !! ONLY ALLOW FOR MODES DEFINED IN Registry.EM
!   CASE(4-10,13-19,22-41,43:)      !! -- CHANGE THIS if YOU ADDED CASES in REGISTRY.EM
!   PRINT *, 'NOT SUPPORTED, to be done'
!   STOP
!END SELECT
!!!!!!!!!!!!!!!!!!! FOR REFERENCE ; FROM REGISTRY.EM
!package   default      mars==0                      -              -
!package   water        mars==1                      -              scalar:qh2o,qh2o_ice
!package   dust1        mars==2                      -              scalar:qdust
!package   dust2eq      mars==3                      -              scalar:qdust,qdustn
!package   newwater     mars==11                     -              scalar:qh2o,qh2o_ice,qdust,qdustn
!package   radioac      mars==20                     -              scalar:qtrac1
!package   radioac2     mars==21                     -              scalar:upward,downward
!package   photochem    mars==42                     -              scalar:qco2,chem_co,chem_o,chem_o1d,chem_o2,chem_o3,chem_h,chem_h2,chem_oh,chem_ho2,chem_h2o2,chem_ch4,chem_n2,chem_ar,qh2o_ice,qh2o,qdust,qdustn
!!!!!!!!!!!!!!!!!!! FOR REFERENCE

#ifdef NEWPHYS
!!! name of tracers to mimic entries in tracer.def 
!!! ----> IT IS IMPORTANT TO KEEP SAME ORDERING AS IN THE REGISTRY !!!!
!!!                               SAME NAMING   AS IN THE LMD PHYSICS !!!!
SELECT CASE (MARS_MODE)
    CASE(0,10)  
      wtnom(nq) = 'co2'
    CASE(1)
      wtnom(1)  = 'h2o_vap'
      wtnom(2)  = 'h2o_ice' 
    CASE(2)
      wtnom(1)  = 'dust01'      
    CASE(3)
      wtnom(1)  = 'dust_mass'
      wtnom(2)  = 'dust_number' 
    CASE(11) 
      wtnom(1)  = 'h2o_vap'
      wtnom(2)  = 'h2o_ice' 
      wtnom(3)  = 'dust_mass'
      wtnom(4)  = 'dust_number'
    CASE(12)
      wtnom(1)  = 'h2o_vap'
      wtnom(2)  = 'h2o_ice'
      wtnom(3)  = 'dust_mass'
      wtnom(4)  = 'dust_number'
      wtnom(5)  = 'ccn_mass'
      wtnom(6)  = 'ccn_number'
    CASE(20) 
      wtnom(1) = 'qtrac1'
    CASE(21)
      wtnom(1) = 'upward'
      wtnom(2) = 'downward' 
    CASE(42)
      wtnom(1)  = 'co2'
      wtnom(2)  = 'co'
      wtnom(3)  = 'o'
      wtnom(4)  = 'o1d'
      wtnom(5)  = 'o2'
      wtnom(6)  = 'o3'
      wtnom(7)  = 'h'
      wtnom(8)  = 'h2'
      wtnom(9)  = 'oh'
      wtnom(10)  = 'ho2'
      wtnom(11)  = 'h2o2'
      wtnom(12)  = 'ch4'
      wtnom(13)  = 'n2'
      wtnom(14)  = 'ar'
      wtnom(15)  = 'h2o_ice'
      wtnom(16)  = 'h2o_vap'
      wtnom(17)  = 'dust_mass'
      wtnom(18)  = 'dust_number'
END SELECT
#endif

!!*******************************************!!
!!*******************************************!!

DO j = jps,jpe
DO i = ips,ipe

!!*******************************************!!
!! FROM 3D WRF FIELDS TO 1D PHYSICS PROFILES !!
!!*******************************************!!

!-----------------------------------!
! 1D subscript for physics "cursor" !
!-----------------------------------!
subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

!--------------------------------------!
! 1D columns : inputs for the physics  !
! ... vert. coord., meteor. fields ... !
!--------------------------------------!
dz8w_prof(:) = dz8w(i,kps:kpe,j)   ! dz between full levels (m)   
p8w_prof(:) = p8w(i,kps:kpe,j)     ! pressure full level (Pa) >> pplev
p_prof(:) = p(i,kps:kpe,j)         ! pressure half level (Pa) >> pplay
t_prof(:) = t(i,kps:kpe,j)         ! temperature half level (K) >> pt
t8w_prof(:) = t8w(i,kps:kpe,j)     ! temperature full level (K)
u_prof(:) = u(i,kps:kpe,j)         ! zonal wind (A-grid: unstaggered) half level (m/s) >> pu  
v_prof(:) = v(i,kps:kpe,j)         ! meridional wind (A-grid: unstaggered) half level (m/s) >> pv
z_prof(:) = z(i,kps:kpe,j)         ! geopotential height half level (m) >> pphi/g
!pi_prof(:) = exner(i,kps:kpe,j)    ! exner function (dimensionless) half level 
!rho_prof(:) = rho(i,kps:kpe,j)     ! density half level  
!th_prof(:) = th(i,kps:kpe,j)       ! pot. temperature half level (K)

!--------------------------------!
! specific treatment for tracers !
!--------------------------------!
#ifndef NOPHYS
#ifdef NEWPHYS
IF (MARS_MODE .EQ. 0) THEN
    q_prof(:,1)=0.95
ELSE
    q_prof(:,1:nq) = SCALAR(i,kps:kpe,j,2:nq+1)  !! the names were set above !! one dummy tracer in WRF
ENDIF
  !!!! CAS DU CO2
  !DO iii=1,nq 
  ! IF ( wtnom(iii) .eq. 'co2' .and. (.not. restart)) q_prof(:,iii) = 0.95
  !ENDDO
  IF ((MARS_MODE .EQ. 20) .OR. (MARS_MODE .EQ. 21)) THEN
   IF (firstcall .EQV. .true. .and. (.not. restart)) THEN
      q_prof(:,:) = 0.95
   ENDIF
  ENDIF
#else
SELECT CASE (MARS_MODE)
    CASE(0)  !! NO TRACERS (mars=0)
    water_vapor_prof(:) = 0.
    water_ice_prof(:) = 0.
    CASE(1)  !! WATER CYCLE (mars=1)
    water_vapor_prof(:) = scalar(i,kps:kpe,j,2)  !! H2O vapor is tracer 1 in the Registry for mars=1
    water_ice_prof(:) = scalar(i,kps:kpe,j,3)    !! H2O ice is tracer 2 in the Registry for mars=1
    CASE(2)  !! DUST CYCLE (mars=2)
    water_vapor_prof(:) = 0.
    water_ice_prof(:) = 0.
    CASE(3:)
    print *, 'OOOPS... not ready yet.'
    STOP
END SELECT
#endif
#endif

!!**********************************************************!!
!! STATIC FIELDS FOR THE PHYSICS - NEEDED ONLY AT FIRSTCALL !!
!!**********************************************************!!
needed_ini_phys : IF (firstcall .EQV. .true.) THEN

!----------------------------------------!
! Surface of each part of the grid (m^2) !
!----------------------------------------!
!aire_val = dx*dy                           !! 1. idealized cases - computational grid
aire_val = (dx/msft(i,j))*(dy/msft(i,j))    !! 2. WRF map scale factors - assume that msfx=msfy (msf=covariance)
!aire_val=dx*dy/msfu(i,j)                   !! 3. special for Mercator GCM-like simulations

!---------------------------------------------!
! Mass-point latitude and longitude (radians) !
!---------------------------------------------!
IF (JULYR .ne. 9999) THEN
 lat_val = XLAT(i,j)*DEGRAD
 lon_val = XLONG(i,j)*DEGRAD
ELSE
 !!! IDEALIZED CASE
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION lat: ',lat_input
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION lon: ',lon_input
 lat_val = lat_input*DEGRAD
 lon_val = lon_input*DEGRAD
ENDIF

!-----------------------------------------!
! Gravity wave parametrization            !
! NB: usually 0 in mesoscale applications !
!-----------------------------------------!
IF (JULYR .ne. 9999) THEN
 zmea_val=M_GW(i,1,j)
 zstd_val=M_GW(i,2,j)
 zsig_val=M_GW(i,3,j)
 zgam_val=M_GW(i,4,j)
 zthe_val=M_GW(i,5,j)
ELSE
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION GWdrag OFF'
 zmea_val=0.
 zstd_val=0.
 zsig_val=0.
 zgam_val=0.
 zthe_val=0.
ENDIF

!---------------------------------!
! Ground albedo & Thermal Inertia !
!---------------------------------!
IF (JULYR .ne. 9999) THEN
 IF (CST_AL == 0) THEN 
 albedodat_val=M_ALBEDO(i,j) 
 ELSE 
 albedodat_val=CST_AL
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT ALBEDO ', albedodat_val
 ENDIF
 IF (CST_TI == 0) THEN 
 inertiedat_val=M_TI(i,j) 
 ELSE 
 inertiedat_val=CST_TI
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT THERMAL INERTIA ', inertiedat_val
 ENDIF
ELSE
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION albedo: ', CST_AL
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION inertia: ',CST_TI
 albedodat_val=CST_AL
 inertiedat_val=CST_TI
ENDIF

#ifdef NEWPHYS
!----------------------------!
! Variable surface roughness !
!----------------------------!
IF (JULYR .ne. 9999) THEN
 IF (CST_Z0 == 0) THEN 
   z0_val = M_Z0(i,j) 
 ELSE 
   z0_val = CST_Z0
   IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** SET CONSTANT SURF ROUGHNESS (m) ',CST_Z0
 ENDIF 
ELSE
 IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION z0 (m) ', CST_Z0
 z0_val=CST_Z0
ENDIF
!!!! ADDITIONAL SECURITY. THIS MIGHT HAPPEN WITH OLD INIT FILES.
IF (z0_val == 0.) THEN
   IF ( (i == ips) .AND. (j == jps) ) PRINT *, 'WELL, z0 is 0, this is no good. Setting to old defaults value 0.01 m'
   z0_val = 0.01
ENDIF
!!!! ADDITIONAL SECURITY. INTERP+SMOOTH IN GEOGRID MIGHT YIELD NEGATIVE Z0 !!!
IF (z0_val < 0.) THEN
   PRINT *, 'WELL, z0 is NEGATIVE, this is impossible. better stop here.'
   PRINT *, 'advice --> correct interpolation / smoothing of z0 in WPS'
   PRINT *, '           -- or check the constant value set in namelist.input'
   STOP
ENDIF
#endif

!---------------------------------------------------------!
! Ground geopotential                                     !
! (=g*HT(i,j), only used in the microphysics: newcondens) !
!---------------------------------------------------------! 
phisfi_val=g*(z_prof(1)-dz8w_prof(1)/2.)   !! NB: z_prof are half levels, so z_prof(1) is not the surface

!-----------------------------------------------!
! Ground temperature, emissivity, CO2 ice cover !
!-----------------------------------------------!
IF (M_TSURF(i,j) .gt. 0.) THEN
  tsurf_val=M_TSURF(i,j)
ELSE
  tsurf_val=TSK(i,j) ! retro-compatibility
ENDIF
emis_val=M_EMISS(i,j)
co2ice_val=M_CO2ICE(i,j)

!------------------------!
! Deep soil temperatures ! 
!------------------------!
IF (JULYR .ne. 9999) THEN
  IF (M_TSOIL(i,1,j) .gt. 0.) THEN
   tsoil_val(:)=M_TSOIL(i,:,j)
  ELSE
   tsoil_val = tsoil_val*0. + tsurf_val
  ENDIF
#ifdef NEWPHYS
  isoil_val(:)=M_ISOIL(i,:,j)
  dsoil_val(:)=M_DSOIL(i,:,j)
#endif
ELSE
  IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION tsoil is set to tsurf'
  do k=1,nsoil
   IF (.not.restart) THEN
     tsoil_val(k) = tsurf_val
   ELSE
     !this is a restart run. We must not set tsoil to tsurf in the init.
     !tsoil was saved in physiq.F under the name M_TSOIL in the restart file
     !(see Registry)
     tsoil_val(k)=M_TSOIL(i,k,j)
   ENDIF

#ifdef NEWPHYS
   IF ( nsoil .lt. 18 ) THEN
       PRINT *,'** Mars ** WRONG NUMBER OF SOIL LAYERS. SHOULD BE 18 AND IT IS ',nsoil
       STOP
   ENDIF
   IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION isoil and dsoil standard'
   isoil_val(k) = inertiedat_val
!  dsoil_val(k) = sqrt(887.75/3.14)*((2.**(k-0.5))-1.) * inertiedat_val / wvolcapa    ! old setting
   dsoil_val(k) = 2.E-4 * (2.**(k-0.5-1.))                                            ! new gcm settings
   IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** ISOIL DSOIL are ',isoil_val(k), dsoil_val(k)
#endif
  enddo
ENDIF

!-------------------!
! Tracer at surface !
!-------------------!
#ifndef NOPHYS
SELECT CASE (MARS_MODE)
    CASE(0)  
    qsurf_val(:)=0.
    CASE(1) 
    qsurf_val(1)=0. 
    qsurf_val(2)=M_H2OICE(i,j)    !! logique avec wtnom(2) = 'h2o_ice' defini ci-dessus
                                   !! ----- retrocompatible ancienne physique
                                   !! ----- [H2O ice is last tracer in qsurf in LMD physics]
    CASE(2)  
    qsurf_val(1)=0.                !! not coupled with lifting for the moment [non remobilise]
#ifdef NEWPHYS
    CASE(3)
    qsurf_val(1)=q_prof(1,1)                !!! temporaire, a definir       
    qsurf_val(2)=q_prof(1,2)     
    CASE(11) 
    qsurf_val(1)=0.
    qsurf_val(2)=M_H2OICE(i,j)    !! logique avec wtnom(2) = 'h2o_ice' defini ci-dessus
    qsurf_val(3)=0.                !! not coupled with lifting for the moment [non remobilise]
    qsurf_val(4)=0.
    CASE(12)    
    qsurf_val(1)=0.
    qsurf_val(2)=M_H2OICE(i,j)    !! logique avec wtnom(2) = 'h2o_ice' defini ci-dessus
    qsurf_val(3)=0.                !! not coupled with lifting for the moment [non remobilise]
    qsurf_val(4)=0.
    qsurf_val(5)=0.
    qsurf_val(6)=0.
    CASE(20)
    qsurf_val(:)=0.
    CASE(21)
    qsurf_val(:)=0.
#else
    CASE(3:)
    print *, 'OOOPS... not ready yet.'
    STOP
#endif

END SELECT
#endif

!-------------------!
! Slope inclination !
!-------------------!
IF (JULYR .ne. 9999) THEN
  theta_val=atan(sqrt( (1000.*SLPX(i,j))**2 + (1000.*SLPY(i,j))**2 ))
  theta_val=theta_val/DEGRAD
ELSE
  IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION slope inclination is 0'
  theta_val=0.
ENDIF

!-------------------------------------------!
! Slope orientation; 0 is north, 90 is east !
!-------------------------------------------!
IF (JULYR .ne. 9999) THEN
  psi_val=-90.*DEGRAD-atan(SLPY(i,j)/SLPX(i,j))
  if (SLPX(i,j) .ge. 0.) then
    psi_val=psi_val-180.*DEGRAD
  endif
  psi_val=360.*DEGRAD+psi_val
  psi_val=psi_val/DEGRAD
  psi_val = MODULO(psi_val+180.,360.)
ELSE
  IF ( (i == ips) .AND. (j == jps) ) PRINT *,'** Mars ** IDEALIZED SIMULATION slope orientation is 0 (well, whatever)'
  psi_val=0.
ENDIF

!-----------------!
! Optional output !
!-----------------!
IF ( (i == ips) .AND. (j == jps) ) THEN
PRINT *,'lat/lon ', lat_val/DEGRAD, lon_val/DEGRAD
PRINT *,'emiss ', emis_val
PRINT *,'albedo ', albedodat_val
PRINT *,'inertie ', inertiedat_val
PRINT *,'phi ',phisfi_val
PRINT *,'tsurf ',tsurf_val
PRINT *,'aire ',aire_val
PRINT *,'z_prof ',z_prof
PRINT *,'dz8w_prof',dz8w_prof
PRINT *,'p8w_prof ',p8w_prof
PRINT *,'p_prof ',p_prof
PRINT *,'t_prof ',t_prof
PRINT *,'t8w_prof ',t8w_prof
PRINT *,'u_prof ',u_prof
PRINT *,'v_prof ',v_prof
PRINT *,'tsoil ',tsoil_val 
PRINT *,'qsurf ',qsurf_val 
#ifdef NEWPHYS
PRINT *,'isoil ',isoil_val
PRINT *,'dsoil ',dsoil_val
PRINT *,'q_prof ',q_prof
PRINT *,'z0 ',z0_val
#endif
ENDIF

!-------------------------! 
!-------------------------!
!      PROVISOIRE         !
!-------------------------!
!-------------------------!
IF (.not. restart) THEN
   q2_val(:)=1.E-6      !PBL wind variance
#ifdef NEWPHYS
   fluxrad_val=0.
   wstar_val=0.
#endif
ELSE
   q2_val(:)=M_Q2(i,:,j)
#ifdef NEWPHYS
   fluxrad_val=M_FLUXRAD(i,j)
   wstar_val=M_WSTAR(i,j)
#endif
ENDIF

!-----------------!
! Fill the arrays !
!-----------------!
aire_vec(subs) = aire_val   !! NB: total area in square km is SUM(aire_vec)/1.0E6
lat_vec(subs) = lat_val
lon_vec(subs) = lon_val
wphisfi(subs) = phisfi_val
walbedodat(subs) = albedodat_val
winertiedat(subs) = inertiedat_val
wzmea(subs) = zmea_val
wzstd(subs) = zstd_val
wzsig(subs) = zsig_val
wzgam(subs) = zgam_val
wzthe(subs) = zthe_val
wtsurf(subs) = tsurf_val
wco2ice(subs) = co2ice_val
wemis(subs) = emis_val
wq2(subs,:) = q2_val(:)
wqsurf(subs,:) = qsurf_val(:)
wtsoil(subs,:) = tsoil_val(:)
#ifdef NEWPHYS
wfluxrad(subs) = fluxrad_val
wwstar(subs) = wstar_val
wisoil(subs,:) = isoil_val(:)
wdsoil(subs,:) = dsoil_val(:)
wz0tab(subs) = z0_val
#endif
wtheta(subs) = theta_val
wpsi(subs) = psi_val

ENDIF needed_ini_phys

!!*****************************!!
!! PREPARE "FOOD" FOR PHYSIQ.F !!
!!*****************************!!

!---------------------------------------------!
! in LMD physics, geopotential must be        !
! expressed with respect to the local surface !
!---------------------------------------------!
pphi(subs,:) = g*( z_prof(:)-(z_prof(1)-dz8w_prof(1)/2.) )

!--------------------------------!
! Dynamic fields for LMD physics !
!--------------------------------!
pplev(subs,1:nlayer) = p8w_prof(1:nlayer)  !! NB: last level: no data
pplay(subs,:) = p_prof(:)
pt(subs,:) = t_prof(:)
pu(subs,:) = u_prof(:)
pv(subs,:) = v_prof(:)
pw(subs,:) = 0   !! NB: not used in the physics, only diagnostic...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! for IDEALIZED CASES ONLY
IF (JULYR .eq. 9999) pplev(subs,nlayer+1)=0.  !! pplev(subs,nlayer+1)=ptop >> NO !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOTE:
! IF ( pplev(subs,nlayer+1) .le. 0 ) pplev(subs,nlayer+1)=ptop
! cree des diagnostics delirants et aleatoires dans le transfert radiatif 

!---------!
! Tracers !  
!---------!
#ifndef NOPHYS
#ifdef NEWPHYS
pq(subs,:,:) = q_prof(:,:)  !! traceurs generiques, seuls noms sont specifiques
#else
SELECT CASE (MARS_MODE)
    CASE(0)  !! NO TRACERS (mars=0)
    pq(subs,:,nq) = water_vapor_prof(:) !! NB: which is 0, actually (see above) 
    CASE(1)  !! WATER CYCLE (mars=1)
    pq(subs,:,nq) = water_vapor_prof(:)
    pq(subs,:,nq-1) = water_ice_prof(:)
    CASE(2)  !! DUST CYCLE (mars=2)
    pq(subs,:,nq) = water_vapor_prof(:) !! NB: which is 0, actually (see above)
    CASE(3:)
    print *, 'OOOPS... not ready yet.'
    STOP    
END SELECT
#endif
#endif
ENDDO
ENDDO

!!*****************!!
!! CLEAN THE PLACE !!
!!*****************!!
DEALLOCATE(q2_val)
DEALLOCATE(qsurf_val)
DEALLOCATE(tsoil_val)
#ifdef NEWPHYS
DEALLOCATE(isoil_val)
DEALLOCATE(dsoil_val)
#endif
DEALLOCATE(dz8w_prof)
DEALLOCATE(z_prof)
DEALLOCATE(p8w_prof)
DEALLOCATE(p_prof)
DEALLOCATE(t_prof)
DEALLOCATE(t8w_prof)
DEALLOCATE(u_prof)
DEALLOCATE(v_prof)
#ifdef NEWPHYS
DEALLOCATE(q_prof)
#else
DEALLOCATE(water_vapor_prof)
DEALLOCATE(water_ice_prof)
#endif
    !!! no use
    !DEALLOCATE(pi_prof)
    !DEALLOCATE(rho_prof)
    !DEALLOCATE(th_prof)


!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!! CALL LMD PHYSICS !!  
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!


!!***********************!!
!! INIFIS, AT FIRST CALL !!
!!***********************!!
IF (firstcall .EQV. .true.) THEN
#ifndef NOPHYS
print *, '** Mars ** LMD INITIALIZATION'
#include "../call_meso_inifis.inc"
!!! le # est important pour newphys
#endif
DEALLOCATE(aire_vec)
DEALLOCATE(lat_vec)
DEALLOCATE(lon_vec)
DEALLOCATE(walbedodat)
DEALLOCATE(winertiedat)
DEALLOCATE(wphisfi)
DEALLOCATE(wzmea)
DEALLOCATE(wzstd)
DEALLOCATE(wzsig)
DEALLOCATE(wzgam)
DEALLOCATE(wzthe)
DEALLOCATE(wtheta)
DEALLOCATE(wpsi)
#ifdef NEWPHYS
DEALLOCATE(wz0tab)
#endif
ENDIF

!!********!!
!! PHYSIQ !!
!!********!!
call_physics : IF (wappel_phys .ne. 0.) THEN

!-------------------------------------------------------------------------------!
! outputs:									!	
!    pdu(ngrid,nlayermx)        \						!
!    pdv(ngrid,nlayermx)         \  Temporal derivative of the corresponding	!
!    pdt(ngrid,nlayermx)         /  variables due to physical processes.	!
!    pdq(ngrid,nlayermx)        /						!
!    pdpsrf(ngrid)             /						!
!    tracerdyn                 call tracer in dynamical part of GCM ?		!
!-------------------------------------------------------------------------------!
pdpsrf(:)=0.
pdu(:,:)=0.
pdv(:,:)=0.
pdt(:,:)=0.
pdq(:,:,:)=0.
#ifndef NOPHYS
print *, '** Mars ** CALL TO LMD PHYSICS'
#include "../call_meso_physiq.inc" 
!!! le # est important pour newphys
#endif
DEALLOCATE(pplev)
DEALLOCATE(pplay)
DEALLOCATE(pphi)
DEALLOCATE(pu)
DEALLOCATE(pv)
DEALLOCATE(pt)
DEALLOCATE(pw)
DEALLOCATE(pq)
!DEALLOCATE(wtsurf)
!DEALLOCATE(wco2ice)
DEALLOCATE(wemis)
!DEALLOCATE(wq2)
!DEALLOCATE(wqsurf)
!DEALLOCATE(wtsoil)
!DEALLOCATE(wwstar)
!DEALLOCATE(wfluxrad)
#ifdef NEWPHYS
DEALLOCATE(wisoil)
DEALLOCATE(wdsoil)
DEALLOCATE(wtnom)
#endif


!-------------------------------!
! PHYSIQ OUTPUT IN THE WRF FILE !
!-------------------------------!
#ifndef NOPHYS
DO j = jps,jpe
DO i = ips,ipe
subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

!!!! SEMBLE T IL PROBLEME AVEC NEWPHYS .... qui marche tres bien sans !
!#ifndef NEWPHYS           
!            !!! sparadrap pour regler un probleme avec mpiifort en LES
!            !!! -- HFX apparaissait nul aux interfaces des tuiles
!            if ( output_tab2d(subs,ind_HFX) .le. 1.e-8 ) then
!               !print *, 'HFX is zero !!! ', i, j
!               !print *, 'I substituted the value right next to it ', output_tab2d(subs+1,ind_HFX) 
!               output_tab2d(subs,ind_HFX) = output_tab2d(subs+1,ind_HFX)
!            endif
!#endif

#include "module_lmd_driver_output3.inc"  
       !  ^-- generated from Registry
ENDDO
ENDDO
DEALLOCATE(output_tab2d)
DEALLOCATE(output_tab3d)
#endif

!!!!!! PATCH SPECIAL STORM
#ifdef NEWPHYS
#ifdef storm
#include "../mars_lmd/storm.inc"
#endif
#endif
!!!!!! PATCH SPECIAL STORM


!---------------------------------------------------------------------------------!
! PHYSIQ TENDENCIES ARE SAVED TO BE SPLIT WITHIN INTERMEDIATE DYNAMICAL TIMESTEPS !
!---------------------------------------------------------------------------------!
#ifdef SPECIAL_NEST_SAVE
dp_save(:,id)=pdpsrf(:)
du_save(:,:,id)=pdu(:,:)
dv_save(:,:,id)=pdv(:,:)
dt_save(:,:,id)=pdt(:,:)
dq_save(:,:,:,id)=pdq(:,:,:)
#ifndef NORESTART
save_tsoil_restart(:,:,id)=wtsoil(:,:)
save_co2ice_restart(:,id)=wco2ice(:)
save_q2_restart(:,:,id)=wq2(:,:)
save_qsurf_restart(:,:,id)=wqsurf(:,:)
save_tsurf_restart(:,id)=wtsurf(:)
#ifdef NEWPHYS
save_wstar_restart(:,id)=wwstar(:)
save_fluxrad_restart(:,id)=wfluxrad(:)
#endif
#endif
#else
dp_save(:)=pdpsrf(:)
du_save(:,:)=pdu(:,:)
dv_save(:,:)=pdv(:,:)
dt_save(:,:)=pdt(:,:)
dq_save(:,:,:)=pdq(:,:,:)
#ifndef NORESTART
save_tsoil_restart(:,:)=wtsoil(:,:)
save_co2ice_restart(:)=wco2ice(:)
save_q2_restart(:,:)=wq2(:,:)
save_qsurf_restart(:,:)=wqsurf(:,:)
save_tsurf_restart(:)=wtsurf(:)
#ifdef NEWPHYS
save_wstar_restart(:)=wwstar(:)
save_fluxrad_restart(:)=wfluxrad(:)
#endif
#endif
#endif
DEALLOCATE(wtsoil)
DEALLOCATE(wco2ice)
DEALLOCATE(wq2)
DEALLOCATE(wqsurf)
DEALLOCATE(wtsurf)
#ifdef NEWPHYS
DEALLOCATE(wfluxrad)
DEALLOCATE(wwstar)
#endif
ENDIF call_physics

ENDIF    ! end of BIG LOOP 2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!***************************!!
!! DEDUCE TENDENCIES FOR WRF !!
!!***************************!!
RTHBLTEN(ims:ime,kms:kme,jms:jme)=0.
RUBLTEN(ims:ime,kms:kme,jms:jme)=0.
RVBLTEN(ims:ime,kms:kme,jms:jme)=0.
PSFC(ims:ime,jms:jme)=p8w(ims:ime,kms,jms:jme) ! was done in surface driver in regular WRF 
!------------------------------------------------------------------!
! WRF adds the physical and dynamical tendencies                   !
! thus the tendencies must be passed as 'per second' tendencies    !
! --when physics is not called, the applied physical tendency      !
! --is the one calculated during the last call to physics          !
!------------------------------------------------------------------!

! PBL top index reading for MARS_MODE 21 :
IF (MARS_MODE .EQ. 21) THEN
      IF (firstcall .EQV. .true.) THEN
           open(unit=15,file='input_zipbl',form='formatted',status='old')
           rewind(15)
           end_of_file = .false.
           kk = 0
           do while (.not. end_of_file)
             read(15,*,end=101) zipbl(kk+1),zilct(kk+1)
             write(*,*) kk, zipbl(kk+1),zilct(kk+1)
             kk = kk+1
             go to 112
101          end_of_file = .true.
112          continue
           enddo
           close(unit=15,status = 'keep')
      ENDIF
      lct=lct_input + elaps/3700.
      zilctmin=MINVAL(ABS(zilct(:)-lct))
      lctindex=0
      DO kk=1,500
        IF(ABS(zilct(kk)-lct) .eq. zilctmin) lctindex=kk
      ENDDO
      IF(lctindex .eq. 0) print*, 'probleme index'
        IF ((zilct(lctindex) .gt. 12.) .and.( zilct(lctindex) .lt. 18)) THEN
            ziindex=MAX(0.,zipbl(lctindex)-7.)
        ELSE
            ziindex=0.
        ENDIF
ENDIF

#ifdef NOPHYS
!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
IF (firstcall .EQV. .true.) THEN
  ALLOCATE(dtrad(nlayer+1))
  DO k=kps,kpe+1
    !! the only solution to avoid 0 points in M_Q2
    !! -- in each domain decomposition case
    dtrad(k) = MAXVAL(M_Q2(:,k,:)) + MINVAL(M_Q2(:,k,:))
  ENDDO
  print*,"HEATING RATE",dtrad(kps:kpe+1)
ENDIF
DO i= 1,ngrid
  pdt(i,kps:kpe) = dtrad(kps:kpe)
ENDDO
!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!
#endif

DO j = jps,jpe
DO i = ips,ipe
subs = (j-jps)*(ipe-ips+1)+(i-ips+1)

!------------!
! zonal wind !
!------------!
RUBLTEN(i,kps:kpe,j) = pdu(subs,kps:kpe)                        

!-----------------!
! meridional wind !
!-----------------!
RVBLTEN(i,kps:kpe,j) = pdv(subs,kps:kpe)                        

!-----------------------!
! potential temperature !
!-----------------------!
! (dT = dtheta * exner for isobaric coordinates or if pressure variations are negligible)
RTHBLTEN(i,kps:kpe,j) = pdt(subs,kps:kpe) / exner(i,kps:kpe,j)   

!---------------------------!
! update surface pressure   !
! (cf CO2 cycle in physics) !
!---------------------------!
PSFC(i,j)=PSFC(i,j)+pdpsrf(subs)*dt    !!! here dt is needed

!------------------------------------!
! Save key variables for restart !  
!------------------------------------!
#ifndef NORESTART
#ifdef SPECIAL_NEST_SAVE
M_TSOIL(i,:,j)=save_tsoil_restart(subs,:,id)
M_CO2ICE(i,j)=save_co2ice_restart(subs,id)
M_Q2(i,kps:kpe+1,j)=save_q2_restart(subs,:,id)
SELECT CASE (MARS_MODE)
   CASE (1,11,12)
     M_H2OICE(i,j)=save_qsurf_restart(subs,2,id)  !! see above Tracer at surface
END SELECT
M_TSURF(i,j)=save_tsurf_restart(subs,id)
#ifdef NEWPHYS
M_WSTAR(i,j)=save_wstar_restart(subs,id)
M_FLUXRAD(i,j)=save_fluxrad_restart(subs,id)
#endif
#else
M_TSOIL(i,:,j)=save_tsoil_restart(subs,:)
M_CO2ICE(i,j)=save_co2ice_restart(subs)
M_Q2(i,kps:kpe+1,j)=save_q2_restart(subs,:)
SELECT CASE (MARS_MODE)
   CASE (1,11,12)
     M_H2OICE(i,j)=save_qsurf_restart(subs,2)  !! see above Tracer at surface
END SELECT
M_TSURF(i,j)=save_tsurf_restart(subs)
#ifdef NEWPHYS
M_WSTAR(i,j)=save_wstar_restart(subs)
M_FLUXRAD(i,j)=save_fluxrad_restart(subs)
#endif
#endif
#endif

!---------!
! Tracers !  
!---------!
#ifndef NOPHYS
#ifdef NEWPHYS
SCALAR(i,kps:kpe,j,1)=0.

SELECT CASE (MARS_MODE)
CASE(0)
        SCALAR(i,kps:kpe,j,:)=0.
CASE(20)
   !! Mars mode 20 : add a passive tracer with radioactive-like decay
   IF ( (i == ips) .AND. (j == jps) )   print *, 'RADIOACTIVE-LIKE TRACER WITH SOURCE AT SURFACE LAYER.'
   tau_decay=60.*10. !! why not make it a namelist argument?
   SCALAR(i,kps:kpe,j,2) = SCALAR(i,kps:kpe,j,2)*exp(-dt/tau_decay)
   SCALAR(i,1,j,2) = SCALAR(i,1,j,2) + 1. !! this tracer is emitted in the surface layer
CASE(21)
   !! Mars mode 21 : add a passive tracer with radioactive-like decay at surface and pbl top
   IF ( (i == ips) .AND. (j == jps) )   print *, 'RADIOACTIVE-LIKE TRACER WITH SOURCE AT SURFACE LAYER AND PBL TOP.'
   tau_decay=60.*10. !! why not make it a namelist argument?
   SCALAR(i,kps:kpe,j,2) = SCALAR(i,kps:kpe,j,2)*exp(-dt/tau_decay)
   SCALAR(i,kps:kpe,j,3) = SCALAR(i,kps:kpe,j,3)*exp(-dt/tau_decay)
   SCALAR(i,1,j,2) = SCALAR(i,1,j,2) + 1. !! this tracer is emitted in the surface layer
   IF ( (i == ips) .AND. (j == jps) )   print *, 'index of pbl emission: ',ziindex
   IF (ziindex .NE. 0) THEN
   SCALAR(i,ziindex,j,3) = SCALAR(i,ziindex,j,3) + 1. !! this tracer is emitted at pbl top
   ENDIF
CASE DEFAULT
   SCALAR(i,kps:kpe,j,2:nq+1)=SCALAR(i,kps:kpe,j,2:nq+1)+pdq(subs,kps:kpe,1:nq)*dt  !!! here dt is needed
END SELECT
#else
SELECT CASE (MARS_MODE)
CASE(0)
        SCALAR(i,kps:kpe,j,:)=0.
CASE(1)
        !!! Water vapor
        SCALAR(i,kps:kpe,j,2)=SCALAR(i,kps:kpe,j,2)+pdq(subs,kps:kpe,nq)*dt  !!! here dt is needed
        !!! Water ice
        SCALAR(i,kps:kpe,j,3)=SCALAR(i,kps:kpe,j,3)+pdq(subs,kps:kpe,nq-1)*dt  !!! here dt is needed
CASE(2)
        !!! Dust 
        SCALAR(i,kps:kpe,j,2)=SCALAR(i,kps:kpe,j,2)+pdq(subs,kps:kpe,nq)*dt  !!! here dt is needed
CASE(3:)
        print *, 'OOOPS... not ready yet.'
        STOP    
END SELECT
#endif
#endif
        !!TODO: check if adding the whole tendency once, and set the
        !!TODO: following tendencies to 0 until physics is called again
        !!TODO: is a good strategy ?
        !RUBLTEN(i,kps:kpe,j) = pdu(subs,kps:kpe)*ptimestep/dt
        !RVBLTEN(i,kps:kpe,j) = pdv(subs,kps:kpe)*ptimestep/dt
        !RTHBLTEN(i,kps:kpe,j) = pdt(subs,kps:kpe)*ptimestep/dt
        !RTHBLTEN(i,kps:kpe,j) = RTHBLTEN(i,kps:kpe,j)/exner(i,kps:kpe,j) 
        !PSFC(i,j)=PSFC(i,j)+pdpsrf(subs)*ptimestep/dt      
        !SELECT CASE (MARS_MODE)
        !CASE(0)
        !SCALAR(i,kps:kpe,j,:)=0.
        !CASE(1)
        !!!! Water vapor
        !SCALAR(i,kps:kpe,j,2)=SCALAR(i,kps:kpe,j,2)+pdq(subs,kps:kpe,nq)*ptimestep/dt
        !!!! Water ice
        !SCALAR(i,kps:kpe,j,3)=SCALAR(i,kps:kpe,j,3)+pdq(subs,kps:kpe,nq-1)*ptimestep/dt
        !CASE(2)
        !!!! Dust 
        !SCALAR(i,kps:kpe,j,2)=SCALAR(i,kps:kpe,j,2)+pdq(subs,kps:kpe,nq)*ptimestep/dt
        !END SELECT

ENDDO
ENDDO
DEALLOCATE(pdpsrf)
DEALLOCATE(pdu)
DEALLOCATE(pdv)
DEALLOCATE(pdt)
DEALLOCATE(pdq)

!!---------!
!! display !
!!---------!
!PRINT *, '** Mars ** Results from LMD physics'
!PRINT *, 'u non-zero tendencies'
!!PRINT *, 'max',MAXVAL(RUBLTEN, MASK=RUBLTEN/=0.),&
!!        ' at',MAXLOC(RUBLTEN, MASK=RUBLTEN/=0.)
!!PRINT *, 'min',MINVAL(RUBLTEN, MASK=RUBLTEN/=0.),&
!!        ' at',MINLOC(RUBLTEN, MASK=RUBLTEN/=0.)
!PRINT *, RUBLTEN(10,1,10), RUBLTEN(10,15,10)
!PRINT *, 'v non-zero tendencies'
!!PRINT *, 'max',MAXVAL(RVBLTEN, MASK=RVBLTEN/=0.),&
!!        ' at',MAXLOC(RVBLTEN, MASK=RVBLTEN/=0.)
!!PRINT *, 'min',MINVAL(RVBLTEN, MASK=RVBLTEN/=0.),&
!!        ' at',MINLOC(RVBLTEN, MASK=RVBLTEN/=0.)
!PRINT *, RVBLTEN(10,1,10), RVBLTEN(10,15,10)
!!! STOP IF CRASH
!IF (MAXVAL(RUBLTEN, MASK=RUBLTEN/=0.) == 0.) STOP
!IF (MAXVAL(RVBLTEN, MASK=RVBLTEN/=0.) == 0.) STOP

!!*****!!
!! END !!
!!*****!!
print *,'** Mars ** END LMD PHYSICS'
!----------------------------------------------------------------!
! use debug (see solve_em) if you wanna display some message ... !
!----------------------------------------------------------------!
END SUBROUTINE lmd_driver

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function ls2sol(ls)

!c  Returns solar longitude, Ls (in deg.), from day number (in sol),
!c  where sol=0=Ls=0 at the northern hemisphere spring equinox

      implicit none

!c  Arguments:
      real ls

!c  Local:
      double precision xref,zx0,zteta,zz
!c	xref: mean anomaly, zteta: true anomaly, zx0: eccentric anomaly
      double precision year_day 
      double precision peri_day,timeperi,e_elips
      double precision pi,degrad 
      parameter (year_day=668.6d0) ! number of sols in a amartian year
!c      data peri_day /485.0/
      parameter (peri_day=485.35d0) ! date (in sols) of perihelion
!c  timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
      parameter (timeperi=1.90258341759902d0)
      parameter (e_elips=0.0934d0)  ! eccentricity of orbit
      parameter (pi=3.14159265358979d0)
      parameter (degrad=57.2957795130823d0)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            ls2sol = 0.0
         else
            ls2sol = year_day
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      ls2sol = zz*year_day + peri_day
      if (ls2sol.lt.0.0) ls2sol = ls2sol + year_day
      if (ls2sol.ge.year_day) ls2sol = ls2sol - year_day

      return
      end function ls2sol
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

END MODULE module_lmd_driver
