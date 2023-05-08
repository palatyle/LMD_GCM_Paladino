!
MODULE cpl_mod
!
! This module excahanges and transforms all fields that should be recieved or sent to 
! coupler. The transformation of the fields are done from the grid 1D-array in phylmd 
! to the regular 2D grid accepted by the coupler. Cumulation of the fields for each 
! timestep is done in here. 
!
! Each type of surface that recevie fields from the coupler have a subroutine named 
! cpl_receive_XXX_fields and each surface that have fields to be sent to the coupler 
! have a subroutine named cpl_send_XXX_fields.
!
!*************************************************************************************

! Use statements
!*************************************************************************************
  USE dimphy, ONLY : klon
  USE mod_phys_lmdz_para
  USE ioipsl
  USE iophy

! The module oasis is always used. Without the cpp key CPP_COUPLE only the parameters 
! in the module are compiled and not the subroutines.
  USE oasis
  USE write_field_phy
  USE control_mod

  
! Global attributes
!*************************************************************************************
  IMPLICIT NONE
  PRIVATE

  ! All subroutine are public except cpl_send_all
  PUBLIC :: cpl_init, cpl_receive_frac, cpl_receive_ocean_fields, cpl_receive_seaice_fields, &
       cpl_send_ocean_fields, cpl_send_seaice_fields, cpl_send_land_fields, &
       cpl_send_landice_fields, gath2cpl
  

! Declaration of module variables
!*************************************************************************************
! variable for coupling period
  INTEGER, SAVE :: nexca
  !$OMP THREADPRIVATE(nexca)

! variables for cumulating fields during a coupling periode :
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_sols, cpl_nsol, cpl_rain
  !$OMP THREADPRIVATE(cpl_sols,cpl_nsol,cpl_rain)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_snow, cpl_evap, cpl_tsol
  !$OMP THREADPRIVATE(cpl_snow,cpl_evap,cpl_tsol)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_fder, cpl_albe, cpl_taux, cpl_tauy
  !$OMP THREADPRIVATE(cpl_fder,cpl_albe,cpl_taux,cpl_tauy)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_windsp
  !$OMP THREADPRIVATE(cpl_windsp)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_taumod
  !$OMP THREADPRIVATE(cpl_taumod)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_atm_co2
  !$OMP THREADPRIVATE(cpl_atm_co2)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_rriv2D, cpl_rcoa2D, cpl_rlic2D
  !$OMP THREADPRIVATE(cpl_rriv2D,cpl_rcoa2D,cpl_rlic2D)

! variables read from coupler :
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sst     ! sea surface temperature
  !$OMP THREADPRIVATE(read_sst)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sit     ! sea ice temperature
  !$OMP THREADPRIVATE(read_sit)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sic     ! sea ice fraction
  !$OMP THREADPRIVATE(read_sic)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_alb_sic ! albedo at sea ice
  !$OMP THREADPRIVATE(read_alb_sic)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_u0, read_v0 ! ocean surface current
  !$OMP THREADPRIVATE(read_u0,read_v0)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_co2     ! ocean co2 flux 
  !$OMP THREADPRIVATE(read_co2)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: unity
  !$OMP THREADPRIVATE(unity)
  INTEGER, SAVE                             :: nidct, nidcs
  !$OMP THREADPRIVATE(nidct,nidcs)

! variables to be sent to the coupler
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_sols2D, cpl_nsol2D, cpl_rain2D
  !$OMP THREADPRIVATE(cpl_sols2D, cpl_nsol2D, cpl_rain2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_snow2D, cpl_evap2D, cpl_tsol2D
  !$OMP THREADPRIVATE(cpl_snow2D, cpl_evap2D, cpl_tsol2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_fder2D, cpl_albe2D
  !$OMP THREADPRIVATE(cpl_fder2D, cpl_albe2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_taux2D, cpl_tauy2D
  !$OMP THREADPRIVATE(cpl_taux2D, cpl_tauy2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_taumod2D
  !$OMP THREADPRIVATE(cpl_taumod2D)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_windsp2D
  !$OMP THREADPRIVATE(cpl_windsp2D)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_atm_co22D
  !$OMP THREADPRIVATE(cpl_atm_co22D)

CONTAINS
!
!************************************************************************************
!
  SUBROUTINE cpl_init(dtime, rlon, rlat)
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, fco2_ocn_day
    USE surface_data

    INCLUDE "dimensions.h"
    INCLUDE "indicesol.h"
    INCLUDE "temps.h"
    INCLUDE "iniprint.h"

! Input arguments
!*************************************************************************************
    REAL, INTENT(IN)                  :: dtime
    REAL, DIMENSION(klon), INTENT(IN) :: rlon, rlat

! Local variables
!*************************************************************************************
    INTEGER                           :: error, sum_error, ig, i
    INTEGER                           :: jf, nhoridct
    INTEGER                           :: nhoridcs
    INTEGER                           :: idtime
    INTEGER                           :: idayref
    INTEGER                           :: npas ! only for OASIS2
    REAL                              :: zjulian
    REAL, DIMENSION(iim,jjm+1)        :: zx_lon, zx_lat
    CHARACTER(len = 20)               :: modname = 'cpl_init'
    CHARACTER(len = 80)               :: abort_message
    CHARACTER(len=80)                 :: clintocplnam, clfromcplnam

!*************************************************************************************
! Calculate coupling period
!
!*************************************************************************************
     
    npas = itaufin/ iphysiq
    nexca = 86400 / dtime
    WRITE(lunout,*)' ##### Ocean couple #####'
    WRITE(lunout,*)' Valeurs des pas de temps'
    WRITE(lunout,*)' npas = ', npas
    WRITE(lunout,*)' nexca = ', nexca
    
!*************************************************************************************
! Allocate variables
!
!*************************************************************************************
    error = 0
    sum_error = 0

    ALLOCATE(unity(klon), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_sols(klon,2), stat = error) 
    sum_error = sum_error + error
    ALLOCATE(cpl_nsol(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rain(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_snow(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_evap(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_tsol(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_fder(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_albe(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_taux(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_tauy(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_windsp(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_taumod(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rriv2D(iim,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rcoa2D(iim,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rlic2D(iim,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(read_sst(iim, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_sic(iim, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_sit(iim, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_alb_sic(iim, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_u0(iim, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_v0(iim, jj_nb), stat = error)
    sum_error = sum_error + error

    IF (carbon_cycle_cpl) THEN
       ALLOCATE(read_co2(iim, jj_nb), stat = error)
       sum_error = sum_error + error
       ALLOCATE(cpl_atm_co2(klon,2), stat = error)
       sum_error = sum_error + error

! Allocate variable in carbon_cycle_mod
       ALLOCATE(fco2_ocn_day(klon), stat = error)
       sum_error = sum_error + error
    END IF

    IF (sum_error /= 0) THEN
       abort_message='Pb allocation variables couplees'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
!*************************************************************************************
! Initialize the allocated varaibles
!
!*************************************************************************************
    DO ig = 1, klon
       unity(ig) = ig
    ENDDO

!*************************************************************************************
! Initialize coupling
!
!*************************************************************************************
    idtime = INT(dtime)
#ifdef CPP_COUPLE
    CALL inicma
#endif

!*************************************************************************************
! initialize NetCDF output
!
!*************************************************************************************
    IF (is_sequential) THEN
       idayref = day_ini
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
       CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlon,zx_lon)
       DO i = 1, iim
          zx_lon(i,1) = rlon(i+1)
          zx_lon(i,jjm+1) = rlon(i+1)
       ENDDO
       CALL gr_fi_ecrit(1,klon,iim,jjm+1,rlat,zx_lat)
       clintocplnam="cpl_atm_tauflx"
       CALL histbeg(clintocplnam, iim,zx_lon(:,1),jjm+1,zx_lat(1,:),&
            1,iim,1,jjm+1, itau_phy,zjulian,dtime,nhoridct,nidct) 
! no vertical axis
       CALL histdef(nidct, 'tauxe','tauxe', &
            "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tauyn','tauyn', &
            "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lon','tmp_lon', &
            "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lat','tmp_lat', &
            "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       DO jf=1,maxsend
         IF (infosend(i)%action) THEN
             CALL histdef(nidct, infosend(i)%name ,infosend(i)%name , &
                "-",iim, jjm+1, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
         ENDIF
       END DO
       CALL histend(nidct)
       CALL histsync(nidct)
       
       clfromcplnam="cpl_atm_sst"
       CALL histbeg(clfromcplnam, iim,zx_lon(:,1),jjm+1,zx_lat(1,:),1,iim,1,jjm+1, &
            0,zjulian,dtime,nhoridcs,nidcs) 
! no vertical axis
       DO jf=1,maxrecv
         IF (inforecv(i)%action) THEN
             CALL histdef(nidcs,inforecv(i)%name ,inforecv(i)%name , &
                "-",iim, jjm+1, nhoridcs, 1, 1, 1, -99, 32, "inst", dtime,dtime)
         ENDIF
       END DO
       CALL histend(nidcs)
       CALL histsync(nidcs)

    ENDIF    ! is_sequential
    

!*************************************************************************************
! compatibility test
!
!*************************************************************************************
    IF (carbon_cycle_cpl .AND. version_ocean=='opa8') THEN
       abort_message='carbon_cycle_cpl does not work with opa8'
       CALL abort_gcm(modname,abort_message,1)
    END IF

  END SUBROUTINE cpl_init
  
!
!*************************************************************************************
!
 
  SUBROUTINE cpl_receive_frac(itime, dtime, pctsrf, is_modified)
! This subroutine receives from coupler for both ocean and seaice
! 4 fields : read_sst, read_sic, read_sit and read_alb_sic. 
! The new sea-ice-land-landice fraction is returned. The others fields 
! are stored in this module.
    USE surface_data
    USE phys_state_var_mod, ONLY : rlon, rlat
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
    
    INCLUDE "indicesol.h"
    INCLUDE "temps.h"
    INCLUDE "iniprint.h"
    INCLUDE "YOMCST.h"
    INCLUDE "dimensions.h"

! Arguments
!************************************************************************************
    INTEGER, INTENT(IN)                        :: itime
    REAL, INTENT(IN)                           :: dtime
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: pctsrf
    LOGICAL, INTENT(OUT)                       :: is_modified

! Local variables
!************************************************************************************
    INTEGER                                 :: j, i, time_sec
    INTEGER                                 :: itau_w
    INTEGER, DIMENSION(iim*(jjm+1))         :: ndexcs
    CHARACTER(len = 20)                     :: modname = 'cpl_receive_frac'
    CHARACTER(len = 80)                     :: abort_message
    REAL, DIMENSION(klon)                   :: read_sic1D
    REAL, DIMENSION(iim,jj_nb,maxrecv)      :: tab_read_flds
    REAL, DIMENSION(klon,nbsrf)             :: pctsrf_old
    REAL, DIMENSION(klon_mpi)               :: rlon_mpi, rlat_mpi
    REAL, DIMENSION(iim, jj_nb)             :: tmp_lon, tmp_lat
    REAL, DIMENSION(iim, jj_nb)             :: tmp_r0

!*************************************************************************************
! Start calculation
! Get fields from coupler
!
!*************************************************************************************

    is_modified=.FALSE.

! Check if right moment to receive from coupler
    IF (MOD(itime, nexca) == 1) THEN
       is_modified=.TRUE.
 
       time_sec=(itime-1)*dtime
#ifdef CPP_COUPLE
!$OMP MASTER
    CALL fromcpl(time_sec, tab_read_flds)
!$OMP END MASTER
#endif
    
! NetCDF output of received fields
       IF (is_sequential) THEN
          ndexcs(:) = 0
          itau_w = itau_phy + itime
          DO i = 1, maxrecv
            IF (inforecv(i)%action) THEN
                CALL histwrite(nidcs,inforecv(i)%name,itau_w,tab_read_flds(:,:,i),iim*(jjm+1),ndexcs)
            ENDIF
          END DO
       ENDIF


! Save each field in a 2D array. 
!$OMP MASTER
       read_sst(:,:)     = tab_read_flds(:,:,idr_sisutw)  ! Sea surface temperature
       read_sic(:,:)     = tab_read_flds(:,:,idr_icecov)  ! Sea ice concentration
       read_alb_sic(:,:) = tab_read_flds(:,:,idr_icealw)  ! Albedo at sea ice
       read_sit(:,:)     = tab_read_flds(:,:,idr_icetem)  ! Sea ice temperature
!$OMP END MASTER

       IF (cpl_current) THEN

! Transform the longitudes and latitudes on 2D arrays
          CALL gather_omp(rlon,rlon_mpi)
          CALL gather_omp(rlat,rlat_mpi)
!$OMP MASTER
          CALL Grid1DTo2D_mpi(rlon_mpi,tmp_lon)
          CALL Grid1DTo2D_mpi(rlat_mpi,tmp_lat)

! Transform the currents from cartesian to spheric coordinates
! tmp_r0 should be zero
          CALL geo2atm(iim, jj_nb, tab_read_flds(:,:,idr_curenx), &
             tab_read_flds(:,:,idr_cureny), tab_read_flds(:,:,idr_curenz), &
               tmp_lon, tmp_lat, &
               read_u0(:,:), read_v0(:,:), tmp_r0(:,:))
!$OMP END MASTER

      ELSE
          read_u0(:,:) = 0.
          read_v0(:,:) = 0.
      ENDIF

       IF (carbon_cycle_cpl) THEN
!$OMP MASTER
           read_co2(:,:) = tab_read_flds(:,:,idr_oceco2) ! CO2 flux
!$OMP END MASTER
       ENDIF

!*************************************************************************************
!  Transform seaice fraction (read_sic : ocean-seaice mask) into global 
!  fraction (pctsrf : ocean-seaice-land-landice mask)
!
!*************************************************************************************
       CALL cpl2gath(read_sic, read_sic1D, klon, unity)

       pctsrf_old(:,:) = pctsrf(:,:)
       DO i = 1, klon
          ! treatment only of points with ocean and/or seaice
          ! old land-ocean mask can not be changed
          IF (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic) > 0.) THEN
             pctsrf(i,is_sic) = (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic)) &
                  * read_sic1D(i)
             pctsrf(i,is_oce) = (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic)) &
                  - pctsrf(i,is_sic)
          ENDIF
       ENDDO

    END IF ! if time to receive

  END SUBROUTINE cpl_receive_frac

!
!*************************************************************************************
!

  SUBROUTINE cpl_receive_ocean_fields(knon, knindex, tsurf_new, u0_new, v0_new)
!
! This routine returns the field for the ocean that has been read from the coupler
! (done earlier with cpl_receive_frac). The field is the temperature.
! The temperature is transformed into 1D array with valid points from index 1 to knon.
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, fco2_ocn_day
    INCLUDE "indicesol.h"

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)      :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: u0_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: v0_new

! Local variables
!*************************************************************************************
    INTEGER                  :: i
    INTEGER, DIMENSION(klon) :: index
    REAL, DIMENSION(klon)    :: sic_new

!*************************************************************************************
! Transform read_sst into compressed 1D variable tsurf_new
!
!*************************************************************************************
    CALL cpl2gath(read_sst, tsurf_new, knon, knindex)
    CALL cpl2gath(read_sic, sic_new, knon, knindex)
    CALL cpl2gath(read_u0, u0_new, knon, knindex)
    CALL cpl2gath(read_v0, v0_new, knon, knindex)

!*************************************************************************************
! Transform read_co2 into uncompressed 1D variable fco2_ocn_day added directly in 
! the module carbon_cycle_mod
!
!*************************************************************************************
    IF (carbon_cycle_cpl) THEN
       DO i=1,klon
          index(i)=i
       END DO
       CALL cpl2gath(read_co2, fco2_ocn_day, klon, index)
    END IF

!*************************************************************************************
! The fields received from the coupler have to be weighted with the fraction of ocean 
! in relation to the total sea-ice+ocean
!
!*************************************************************************************
    DO i=1, knon
       tsurf_new(i) = tsurf_new(i)/(1. - sic_new(i))
    END DO

  END SUBROUTINE cpl_receive_ocean_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_receive_seaice_fields(knon, knindex, &
       tsurf_new, alb_new, u0_new, v0_new)
!
! This routine returns the fields for the seaice that have been read from the coupler
! (done earlier with cpl_receive_frac). These fields are the temperature and 
! albedo at sea ice surface and fraction of sea ice.
! The fields are transformed into 1D arrays with valid points from index 1 to knon. 
!

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)      :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: alb_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: u0_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: v0_new

! Local variables
!*************************************************************************************
    INTEGER               :: i
    REAL, DIMENSION(klon) :: sic_new

!*************************************************************************************
! Transform fields read from coupler from 2D into compressed 1D variables
!
!*************************************************************************************
    CALL cpl2gath(read_sit, tsurf_new, knon, knindex)
    CALL cpl2gath(read_alb_sic, alb_new, knon, knindex)
    CALL cpl2gath(read_sic, sic_new, knon, knindex)
    CALL cpl2gath(read_u0, u0_new, knon, knindex)
    CALL cpl2gath(read_v0, v0_new, knon, knindex)

!*************************************************************************************
! The fields received from the coupler have to be weighted with the sea-ice 
! concentration (in relation to the total sea-ice + ocean).
!
!*************************************************************************************
    DO i= 1, knon
       tsurf_new(i) = tsurf_new(i) / sic_new(i)
       alb_new(i)   = alb_new(i)   / sic_new(i)
    END DO

  END SUBROUTINE cpl_receive_seaice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_ocean_fields(itime, knon, knindex, &
       swdown, lwdown, fluxlat, fluxsens, &
       precip_rain, precip_snow, evap, tsurf, fder, albsol, taux, tauy, windsp)
!
! This subroutine cumulates some fields for each time-step during a coupling 
! period. At last time-step in a coupling period the fields are transformed to the 
! grid accepted by the coupler. No sending to the coupler will be done from here 
! (it is done in cpl_send_seaice_fields).
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, co2_send
    INCLUDE "indicesol.h"
    INCLUDE "dimensions.h"

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, DIMENSION(klon), INTENT(IN)       :: swdown, lwdown 
    REAL, DIMENSION(klon), INTENT(IN)       :: fluxlat, fluxsens
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: evap, tsurf, fder, albsol
    REAL, DIMENSION(klon), INTENT(IN)       :: taux, tauy, windsp

! Local variables
!*************************************************************************************
    INTEGER                                 :: cpl_index, ig 
    INTEGER                                 :: error, sum_error
    CHARACTER(len = 25)                     :: modname = 'cpl_send_ocean_fields'
    CHARACTER(len = 80)                     :: abort_message

!*************************************************************************************
! Start calculation
! The ocean points are saved with second array index=1
!
!*************************************************************************************
    cpl_index = 1

!*************************************************************************************
! Reset fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
       cpl_sols(1:knon,cpl_index) = 0.0
       cpl_nsol(1:knon,cpl_index) = 0.0
       cpl_rain(1:knon,cpl_index) = 0.0
       cpl_snow(1:knon,cpl_index) = 0.0
       cpl_evap(1:knon,cpl_index) = 0.0
       cpl_tsol(1:knon,cpl_index) = 0.0
       cpl_fder(1:knon,cpl_index) = 0.0
       cpl_albe(1:knon,cpl_index) = 0.0
       cpl_taux(1:knon,cpl_index) = 0.0
       cpl_tauy(1:knon,cpl_index) = 0.0
       cpl_windsp(1:knon,cpl_index) = 0.0
       cpl_taumod(1:knon,cpl_index) = 0.0
       IF (carbon_cycle_cpl) cpl_atm_co2(1:knon,cpl_index) = 0.0
    ENDIF
       
!*************************************************************************************
! Cumulate at each time-step
!
!*************************************************************************************    
    DO ig = 1, knon
       cpl_sols(ig,cpl_index) = cpl_sols(ig,cpl_index) + &
            swdown(ig)      / REAL(nexca)
       cpl_nsol(ig,cpl_index) = cpl_nsol(ig,cpl_index) + &
            (lwdown(ig) + fluxlat(ig) +fluxsens(ig)) / REAL(nexca)
       cpl_rain(ig,cpl_index) = cpl_rain(ig,cpl_index) + &
            precip_rain(ig) / REAL(nexca)
       cpl_snow(ig,cpl_index) = cpl_snow(ig,cpl_index) + &
            precip_snow(ig) / REAL(nexca)
       cpl_evap(ig,cpl_index) = cpl_evap(ig,cpl_index) + &
            evap(ig)        / REAL(nexca)
       cpl_tsol(ig,cpl_index) = cpl_tsol(ig,cpl_index) + &
            tsurf(ig)       / REAL(nexca)
       cpl_fder(ig,cpl_index) = cpl_fder(ig,cpl_index) + &
            fder(ig)        / REAL(nexca)
       cpl_albe(ig,cpl_index) = cpl_albe(ig,cpl_index) + &
            albsol(ig)      / REAL(nexca)
       cpl_taux(ig,cpl_index) = cpl_taux(ig,cpl_index) + &
            taux(ig)        / REAL(nexca)
       cpl_tauy(ig,cpl_index) = cpl_tauy(ig,cpl_index) + &
            tauy(ig)        / REAL(nexca)      
       cpl_windsp(ig,cpl_index) = cpl_windsp(ig,cpl_index) + &
            windsp(ig)      / REAL(nexca)
       cpl_taumod(ig,cpl_index) =   cpl_taumod(ig,cpl_index) + &
          SQRT ( taux(ig)*taux(ig)+tauy(ig)*tauy(ig) ) / REAL (nexca)

       IF (carbon_cycle_cpl) THEN
          cpl_atm_co2(ig,cpl_index) = cpl_atm_co2(ig,cpl_index) + &
               co2_send(knindex(ig))/ REAL(nexca) 
       END IF
     ENDDO

!*************************************************************************************
! If the time-step corresponds to the end of coupling period the 
! fields are transformed to the 2D grid. 
! No sending to the coupler (it is done from cpl_send_seaice_fields).
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 0) THEN

       IF (.NOT. ALLOCATED(cpl_sols2D)) THEN
          sum_error = 0
          ALLOCATE(cpl_sols2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_nsol2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_rain2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_snow2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_evap2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tsol2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_fder2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_albe2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taux2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tauy2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_windsp2D(iim,jj_nb), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taumod2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          
          IF (carbon_cycle_cpl) THEN
             ALLOCATE(cpl_atm_co22D(iim,jj_nb), stat=error)
             sum_error = sum_error + error
          END IF

          IF (sum_error /= 0) THEN
             abort_message='Pb allocation variables couplees pour l''ecriture'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       

       CALL gath2cpl(cpl_sols(:,cpl_index), cpl_sols2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_nsol(:,cpl_index), cpl_nsol2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_rain(:,cpl_index), cpl_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_snow(:,cpl_index), cpl_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_evap(:,cpl_index), cpl_evap2D(:,:,cpl_index), &
            knon, knindex)

! cpl_tsol2D(:,:,:) not used!
       CALL gath2cpl(cpl_tsol(:,cpl_index), cpl_tsol2D(:,:, cpl_index), &
            knon, knindex)

! cpl_fder2D(:,:,1) not used, only cpl_fder(:,:,2)!
       CALL gath2cpl(cpl_fder(:,cpl_index), cpl_fder2D(:,:,cpl_index), &
            knon, knindex)

! cpl_albe2D(:,:,:) not used!
       CALL gath2cpl(cpl_albe(:,cpl_index), cpl_albe2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taux(:,cpl_index), cpl_taux2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_tauy(:,cpl_index), cpl_tauy2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_windsp(:,cpl_index), cpl_windsp2D(:,:), &
            knon, knindex)

       CALL gath2cpl(cpl_taumod(:,cpl_index), cpl_taumod2D(:,:,cpl_index), &
            knon, knindex)

       IF (carbon_cycle_cpl) &
            CALL gath2cpl(cpl_atm_co2(:,cpl_index), cpl_atm_co22D(:,:), knon, knindex)
   ENDIF

  END SUBROUTINE cpl_send_ocean_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_seaice_fields(itime, dtime, knon, knindex, &
       pctsrf, lafin, rlon, rlat, &
       swdown, lwdown, fluxlat, fluxsens, &
       precip_rain, precip_snow, evap, tsurf, fder, albsol, taux, tauy)
!
! This subroutine cumulates some fields for each time-step during a coupling 
! period. At last time-step in a coupling period the fields are transformed to the 
! grid accepted by the coupler. All fields for all types of surfaces are sent to
! the coupler.
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
    INCLUDE "indicesol.h"
    INCLUDE "dimensions.h"

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon), INTENT(IN)       :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)       :: swdown, lwdown 
    REAL, DIMENSION(klon), INTENT(IN)       :: fluxlat, fluxsens
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: evap, tsurf, fder
    REAL, DIMENSION(klon), INTENT(IN)       :: albsol, taux, tauy
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf
    LOGICAL, INTENT(IN)                     :: lafin

! Local variables
!*************************************************************************************
    INTEGER                                 :: cpl_index, ig 
    INTEGER                                 :: error, sum_error
    CHARACTER(len = 25)                     :: modname = 'cpl_send_seaice_fields'
    CHARACTER(len = 80)                     :: abort_message
    REAL, DIMENSION(klon)                   :: cpl_fder_tmp

!*************************************************************************************
! Start calulation
! The sea-ice points are saved with second array index=2
!
!*************************************************************************************
    cpl_index = 2

!*************************************************************************************
! Reset fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
       cpl_sols(1:knon,cpl_index) = 0.0
       cpl_nsol(1:knon,cpl_index) = 0.0
       cpl_rain(1:knon,cpl_index) = 0.0
       cpl_snow(1:knon,cpl_index) = 0.0
       cpl_evap(1:knon,cpl_index) = 0.0
       cpl_tsol(1:knon,cpl_index) = 0.0
       cpl_fder(1:knon,cpl_index) = 0.0
       cpl_albe(1:knon,cpl_index) = 0.0
       cpl_taux(1:knon,cpl_index) = 0.0
       cpl_tauy(1:knon,cpl_index) = 0.0
       cpl_taumod(1:knon,cpl_index) = 0.0
    ENDIF
       
!*************************************************************************************
! Cumulate at each time-step
!
!*************************************************************************************    
    DO ig = 1, knon
       cpl_sols(ig,cpl_index) = cpl_sols(ig,cpl_index) + &
            swdown(ig)      / REAL(nexca)
       cpl_nsol(ig,cpl_index) = cpl_nsol(ig,cpl_index) + &
            (lwdown(ig) + fluxlat(ig) +fluxsens(ig)) / REAL(nexca)
       cpl_rain(ig,cpl_index) = cpl_rain(ig,cpl_index) + &
            precip_rain(ig) / REAL(nexca)
       cpl_snow(ig,cpl_index) = cpl_snow(ig,cpl_index) + &
            precip_snow(ig) / REAL(nexca)
       cpl_evap(ig,cpl_index) = cpl_evap(ig,cpl_index) + &
            evap(ig)        / REAL(nexca)
       cpl_tsol(ig,cpl_index) = cpl_tsol(ig,cpl_index) + &
            tsurf(ig)       / REAL(nexca)
       cpl_fder(ig,cpl_index) = cpl_fder(ig,cpl_index) + &
            fder(ig)        / REAL(nexca)
       cpl_albe(ig,cpl_index) = cpl_albe(ig,cpl_index) + &
            albsol(ig)      / REAL(nexca)
       cpl_taux(ig,cpl_index) = cpl_taux(ig,cpl_index) + &
            taux(ig)        / REAL(nexca)
       cpl_tauy(ig,cpl_index) = cpl_tauy(ig,cpl_index) + &
            tauy(ig)        / REAL(nexca)     
       cpl_taumod(ig,cpl_index) = cpl_taumod(ig,cpl_index) + &
            SQRT ( taux(ig)*taux(ig)+tauy(ig)*tauy(ig) ) / REAL(nexca) 
    ENDDO

!*************************************************************************************
! If the time-step corresponds to the end of coupling period the 
! fields are transformed to the 2D grid and all fields are sent to coupler.
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 0) THEN
       IF (.NOT. ALLOCATED(cpl_sols2D)) THEN
          sum_error = 0
          ALLOCATE(cpl_sols2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_nsol2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_rain2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_snow2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_evap2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tsol2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_fder2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_albe2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taux2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tauy2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_windsp2D(iim,jj_nb), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taumod2D(iim,jj_nb,2), stat=error)
          sum_error = sum_error + error

          IF (carbon_cycle_cpl) THEN
             ALLOCATE(cpl_atm_co22D(iim,jj_nb), stat=error)
             sum_error = sum_error + error
          END IF

          IF (sum_error /= 0) THEN
             abort_message='Pb allocation variables couplees pour l''ecriture'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF

       CALL gath2cpl(cpl_sols(:,cpl_index), cpl_sols2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_nsol(:,cpl_index), cpl_nsol2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_rain(:,cpl_index), cpl_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_snow(:,cpl_index), cpl_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_evap(:,cpl_index), cpl_evap2D(:,:,cpl_index), &
            knon, knindex)

! cpl_tsol2D(:,:,:) not used!
       CALL gath2cpl(cpl_tsol(:,cpl_index), cpl_tsol2D(:,:, cpl_index), &
            knon, knindex)

       ! Set default value and decompress before gath2cpl
       cpl_fder_tmp(:) = -20.
       DO ig = 1, knon
          cpl_fder_tmp(knindex(ig))=cpl_fder(ig,cpl_index)
       END DO
       CALL gath2cpl(cpl_fder_tmp(:), cpl_fder2D(:,:,cpl_index), &
            klon, unity)

! cpl_albe2D(:,:,:) not used!
       CALL gath2cpl(cpl_albe(:,cpl_index), cpl_albe2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taux(:,cpl_index), cpl_taux2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_tauy(:,cpl_index), cpl_tauy2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taumod(:,cpl_index), cpl_taumod2D(:,:,cpl_index), &
            knon, knindex)

       ! Send all fields
       CALL cpl_send_all(itime, dtime, pctsrf, lafin, rlon, rlat)
    ENDIF

  END SUBROUTINE cpl_send_seaice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_land_fields(itime, knon, knindex, rriv_in, rcoa_in)
!
! This subroutine cumulates some fields for each time-step during a coupling 
! period. At last time-step in a coupling period the fields are transformed to the 
! grid accepted by the coupler. No sending to the coupler will be done from here 
! (it is done in cpl_send_seaice_fields).
!
    INCLUDE "dimensions.h"

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                       :: itime
    INTEGER, INTENT(IN)                       :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex
    REAL, DIMENSION(klon), INTENT(IN)         :: rriv_in
    REAL, DIMENSION(klon), INTENT(IN)         :: rcoa_in

! Local variables
!*************************************************************************************
    REAL, DIMENSION(iim,jj_nb)             :: rriv2D
    REAL, DIMENSION(iim,jj_nb)             :: rcoa2D

!*************************************************************************************
! Rearrange fields in 2D variables 
! First initialize to zero to avoid unvalid points causing problems
!
!*************************************************************************************
!$OMP MASTER
    rriv2D(:,:) = 0.0
    rcoa2D(:,:) = 0.0
!$OMP END MASTER
    CALL gath2cpl(rriv_in, rriv2D, knon, knindex)
    CALL gath2cpl(rcoa_in, rcoa2D, knon, knindex)

!*************************************************************************************
! Reset cumulated fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
!$OMP MASTER
       cpl_rriv2D(:,:) = 0.0
       cpl_rcoa2D(:,:) = 0.0
!$OMP END MASTER
    ENDIF

!*************************************************************************************
! Cumulate : Following fields should be cumulated at each time-step
!
!*************************************************************************************    
!$OMP MASTER
    cpl_rriv2D(:,:) = cpl_rriv2D(:,:) + rriv2D(:,:) / REAL(nexca)
    cpl_rcoa2D(:,:) = cpl_rcoa2D(:,:) + rcoa2D(:,:) / REAL(nexca)
!$OMP END MASTER

  END SUBROUTINE cpl_send_land_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_landice_fields(itime, knon, knindex, rlic_in)
! This subroutine cumulates the field for melting ice for each time-step 
! during a coupling period. This routine will not send to coupler. Sending 
! will be done in cpl_send_seaice_fields.
!

    INCLUDE "dimensions.h"

! Input varibales
!*************************************************************************************
    INTEGER, INTENT(IN)                       :: itime
    INTEGER, INTENT(IN)                       :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex
    REAL, DIMENSION(klon), INTENT(IN)         :: rlic_in

! Local varibales
!*************************************************************************************
    REAL, DIMENSION(iim,jj_nb)             :: rlic2D

!*************************************************************************************
! Rearrange field in a 2D variable 
! First initialize to zero to avoid unvalid points causing problems
!
!*************************************************************************************
!$OMP MASTER
    rlic2D(:,:) = 0.0
!$OMP END MASTER
    CALL gath2cpl(rlic_in, rlic2D, knon, knindex)

!*************************************************************************************
! Reset field to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
!$OMP MASTER
       cpl_rlic2D(:,:) = 0.0
!$OMP END MASTER
    ENDIF

!*************************************************************************************
! Cumulate : Melting ice should be cumulated at each time-step
!
!*************************************************************************************    
!$OMP MASTER
    cpl_rlic2D(:,:) = cpl_rlic2D(:,:) + rlic2D(:,:) / REAL(nexca)
!$OMP END MASTER

  END SUBROUTINE cpl_send_landice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_all(itime, dtime, pctsrf, lafin, rlon, rlat)
! This routine will send fields for all different surfaces to the coupler.
! This subroutine should be executed after calculations by the last surface(sea-ice),
! all calculations at the different surfaces have to be done before. 
!    
    USE surface_data
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
! Some includes
!*************************************************************************************
    INCLUDE "indicesol.h"
    INCLUDE "temps.h"
    INCLUDE "dimensions.h"
    
! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                                  :: itime
    REAL, INTENT(IN)                                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)                    :: rlon, rlat
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)              :: pctsrf
    LOGICAL, INTENT(IN)                                  :: lafin
    
! Local variables
!*************************************************************************************
    INTEGER                                              :: error, sum_error, j
    INTEGER                                              :: itau_w
    INTEGER                                              :: time_sec
    INTEGER, DIMENSION(iim*(jjm+1))                      :: ndexct
    REAL                                                 :: Up, Down
    REAL, DIMENSION(iim, jj_nb)                          :: tmp_lon, tmp_lat
    REAL, DIMENSION(iim, jj_nb, 4)                       :: pctsrf2D
    REAL, DIMENSION(iim, jj_nb)                          :: deno
    CHARACTER(len = 20)                                  :: modname = 'cpl_send_all'
    CHARACTER(len = 80)                                  :: abort_message
   
! Variables with fields to coupler
    REAL, DIMENSION(iim, jj_nb)                          :: tmp_taux
    REAL, DIMENSION(iim, jj_nb)                          :: tmp_tauy
    REAL, DIMENSION(iim, jj_nb)                          :: tmp_calv
! Table with all fields to send to coupler
    REAL, DIMENSION(iim, jj_nb, maxsend)                 :: tab_flds
    REAL, DIMENSION(klon_mpi)                            :: rlon_mpi, rlat_mpi

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                  :: status
#endif

! End definitions
!*************************************************************************************
    


!*************************************************************************************
! All fields are stored in a table tab_flds(:,:,:)
! First store the fields which are already on the right format
!
!*************************************************************************************
!$OMP MASTER
    tab_flds(:,:,ids_windsp) = cpl_windsp2D(:,:)
    tab_flds(:,:,ids_shfice) = cpl_sols2D(:,:,2)
    tab_flds(:,:,ids_nsfice) = cpl_nsol2D(:,:,2)
    tab_flds(:,:,ids_dflxdt) = cpl_fder2D(:,:,2)
    
    IF (version_ocean=='nemo') THEN
       tab_flds(:,:,ids_liqrun) = cpl_rriv2D(:,:) + cpl_rcoa2D(:,:)
       IF (carbon_cycle_cpl) tab_flds(:,:,ids_atmco2)=cpl_atm_co22D(:,:)
    ELSE IF (version_ocean=='opa8') THEN
       tab_flds(:,:,ids_shfoce) = cpl_sols2D(:,:,1)
       tab_flds(:,:,ids_nsfoce) = cpl_nsol2D(:,:,1)
       tab_flds(:,:,ids_icevap) = cpl_evap2D(:,:,2)
       tab_flds(:,:,ids_ocevap) = cpl_evap2D(:,:,1)
       tab_flds(:,:,ids_runcoa) = cpl_rcoa2D(:,:)
       tab_flds(:,:,ids_rivflu) = cpl_rriv2D(:,:)
    END IF

!*************************************************************************************
! Transform the fraction of sub-surfaces from 1D to 2D array
!
!*************************************************************************************
    pctsrf2D(:,:,:) = 0.
!$OMP END MASTER
    CALL gath2cpl(pctsrf(:,is_oce), pctsrf2D(:,:,is_oce), klon, unity)
    CALL gath2cpl(pctsrf(:,is_sic), pctsrf2D(:,:,is_sic), klon, unity)
    CALL gath2cpl(pctsrf(:,is_lic), pctsrf2D(:,:,is_lic), klon, unity)

!*************************************************************************************
! Calculate the average calving per latitude
! Store calving in tab_flds(:,:,19)
! 
!*************************************************************************************      
    IF (is_omp_root) THEN

      DO j = 1, jj_nb
         tmp_calv(:,j) = DOT_PRODUCT (cpl_rlic2D(1:iim,j), &
              pctsrf2D(1:iim,j,is_lic)) / REAL(iim)
      ENDDO
    
    
      IF (is_parallel) THEN
         IF (.NOT. is_north_pole) THEN
#ifdef CPP_MPI
            CALL MPI_RECV(Up,1,MPI_REAL_LMDZ,mpi_rank-1,1234,COMM_LMDZ_PHY,status,error)
            CALL MPI_SEND(tmp_calv(1,1),1,MPI_REAL_LMDZ,mpi_rank-1,1234,COMM_LMDZ_PHY,error)
#endif
         ENDIF
       
         IF (.NOT. is_south_pole) THEN
#ifdef CPP_MPI
            CALL MPI_SEND(tmp_calv(1,jj_nb),1,MPI_REAL_LMDZ,mpi_rank+1,1234,COMM_LMDZ_PHY,error)
            CALL MPI_RECV(down,1,MPI_REAL_LMDZ,mpi_rank+1,1234,COMM_LMDZ_PHY,status,error)
#endif
         ENDIF
         
         IF (.NOT. is_north_pole .AND. ii_begin /=1) THEN
            Up=Up+tmp_calv(iim,1)
            tmp_calv(:,1)=Up
         ENDIF
         
         IF (.NOT. is_south_pole .AND. ii_end /= iim) THEN
            Down=Down+tmp_calv(1,jj_nb)
            tmp_calv(:,jj_nb)=Down	 
         ENDIF
      ENDIF
      
      tab_flds(:,:,ids_calvin) = tmp_calv(:,:)

!*************************************************************************************
! Calculate total flux for snow, rain and wind with weighted addition using the 
! fractions of ocean and seaice.
!
!*************************************************************************************    
       ! fraction oce+seaice
       deno =  pctsrf2D(:,:,is_oce) + pctsrf2D(:,:,is_sic) 

       IF (version_ocean=='nemo') THEN
          tab_flds(:,:,ids_shftot)  = 0.0
          tab_flds(:,:,ids_nsftot) = 0.0
          tab_flds(:,:,ids_totrai) = 0.0
          tab_flds(:,:,ids_totsno) = 0.0
          tab_flds(:,:,ids_toteva) = 0.0
          tab_flds(:,:,ids_taumod) = 0.0
  
          tmp_taux(:,:)    = 0.0
          tmp_tauy(:,:)    = 0.0
          ! For all valid grid cells containing some fraction of ocean or sea-ice
          WHERE ( deno(:,:) /= 0 )
             tmp_taux = cpl_taux2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taux2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tmp_tauy = cpl_tauy2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_tauy2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)

             tab_flds(:,:,ids_shftot) = cpl_sols2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_sols2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_nsftot) = cpl_nsol2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_nsol2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totrai) = cpl_rain2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_rain2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totsno) = cpl_snow2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_snow2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_toteva) = cpl_evap2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_evap2D(:,:,2)  * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_taumod) = cpl_taumod2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taumod2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             
         ENDWHERE

          tab_flds(:,:,ids_icevap) = cpl_evap2D(:,:,2) 
          
       ELSE IF (version_ocean=='opa8') THEN
          ! Store fields for rain and snow in tab_flds(:,:,15) and tab_flds(:,:,16)
          tab_flds(:,:,ids_totrai) = 0.0
          tab_flds(:,:,ids_totsno) = 0.0
          tmp_taux(:,:)    = 0.0
          tmp_tauy(:,:)    = 0.0
          ! For all valid grid cells containing some fraction of ocean or sea-ice
          WHERE ( deno(:,:) /= 0 )
             tab_flds(:,:,ids_totrai) = cpl_rain2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_rain2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totsno) = cpl_snow2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_snow2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             
             tmp_taux = cpl_taux2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taux2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tmp_tauy = cpl_tauy2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_tauy2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
          ENDWHERE
       END IF

    ENDIF ! is_omp_root
  
!*************************************************************************************
! Transform the wind components from local atmospheric 2D coordinates to geocentric 
! 3D coordinates. 
! Store the resulting wind components in tab_flds(:,:,1:6)
!*************************************************************************************

! Transform the longitudes and latitudes on 2D arrays
    
    CALL gather_omp(rlon,rlon_mpi)
    CALL gather_omp(rlat,rlat_mpi)
!$OMP MASTER
    CALL Grid1DTo2D_mpi(rlon_mpi,tmp_lon)
    CALL Grid1DTo2D_mpi(rlat_mpi,tmp_lat)
!$OMP END MASTER    

    IF (is_sequential) THEN
       IF (is_north_pole) tmp_lon(:,1)     = tmp_lon(:,2)
       IF (is_south_pole) tmp_lon(:,jjm+1) = tmp_lon(:,jjm)
    ENDIF
      
! NetCDF output of the wind before transformation of coordinate system
    IF (is_sequential) THEN
       ndexct(:) = 0
       itau_w = itau_phy + itime
       CALL histwrite(nidct,'tauxe',itau_w,tmp_taux,iim*(jjm+1),ndexct)
       CALL histwrite(nidct,'tauyn',itau_w,tmp_tauy,iim*(jjm+1),ndexct)
       CALL histwrite(nidct,'tmp_lon',itau_w,tmp_lon,iim*(jjm+1),ndexct)
       CALL histwrite(nidct,'tmp_lat',itau_w,tmp_lat,iim*(jjm+1),ndexct)
    ENDIF

! Transform the wind from spherical atmospheric 2D coordinates to geocentric
! cartesian 3D coordinates 
!$OMP MASTER
    CALL atm2geo (iim, jj_nb, tmp_taux, tmp_tauy, tmp_lon, tmp_lat, &
         tab_flds(:,:,ids_tauxxu), tab_flds(:,:,ids_tauyyu), tab_flds(:,:,ids_tauzzu) )
    
    tab_flds(:,:,ids_tauxxv)  = tab_flds(:,:,ids_tauxxu)
    tab_flds(:,:,ids_tauyyv)  = tab_flds(:,:,ids_tauyyu)
    tab_flds(:,:,ids_tauzzv)  = tab_flds(:,:,ids_tauzzu)
!$OMP END MASTER

!*************************************************************************************
! NetCDF output of all fields just before sending to coupler.
!
!*************************************************************************************
    IF (is_sequential) THEN
        DO j=1,maxsend
          IF (infosend(j)%action) CALL histwrite(nidct,infosend(j)%name, itau_w, &
             tab_flds(:,:,j),iim*(jjm+1),ndexct)
        ENDDO
    ENDIF
!*************************************************************************************
! Send the table of all fields
!
!*************************************************************************************
    time_sec=(itime-1)*dtime
#ifdef CPP_COUPLE
!$OMP MASTER
    CALL intocpl(time_sec, lafin, tab_flds(:,:,:))
!$OMP END MASTER
#endif

!*************************************************************************************
! Finish with some dellocate
!
!*************************************************************************************  
    sum_error=0
    DEALLOCATE(cpl_sols2D, cpl_nsol2D, cpl_rain2D, cpl_snow2D, stat=error )
    sum_error = sum_error + error
    DEALLOCATE(cpl_evap2D, cpl_tsol2D, cpl_fder2D, cpl_albe2D, stat=error )
    sum_error = sum_error + error
    DEALLOCATE(cpl_taux2D, cpl_tauy2D, cpl_windsp2D, cpl_taumod2D, stat=error )
    sum_error = sum_error + error
    
    IF (carbon_cycle_cpl) THEN
       DEALLOCATE(cpl_atm_co22D, stat=error )
       sum_error = sum_error + error
    END IF

    IF (sum_error /= 0) THEN
       abort_message='Pb in deallocation of cpl_xxxx2D coupling variables'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    
  END SUBROUTINE cpl_send_all
!
!*************************************************************************************
!
  SUBROUTINE cpl2gath(champ_in, champ_out, knon, knindex)
  USE mod_phys_lmdz_para
! Cette routine transforme un champs de la grille 2D recu du coupleur sur la grille 
! 'gathered' (la grille physiq comprime).
!
! 
! input:         
!   champ_in     champ sur la grille 2D
!   knon         nombre de points dans le domaine a traiter
!   knindex      index des points de la surface a traiter
!
! output:
!   champ_out    champ sur la grille 'gatherd'
!
    INCLUDE "dimensions.h"

! Input
    INTEGER, INTENT(IN)                       :: knon
    REAL, DIMENSION(iim,jj_nb), INTENT(IN)    :: champ_in
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex

! Output
    REAL, DIMENSION(klon_mpi), INTENT(OUT)        :: champ_out

! Local
    INTEGER                                   :: i, ig
    REAL, DIMENSION(klon_mpi)                 :: temp_mpi
    REAL, DIMENSION(klon)                     :: temp_omp

!*************************************************************************************
!
    

! Transform from 2 dimensions (iim,jj_nb) to 1 dimension (klon)
!$OMP MASTER 
    CALL Grid2Dto1D_mpi(champ_in,temp_mpi)
!$OMP END MASTER

    CALL scatter_omp(temp_mpi,temp_omp)
    
! Compress from klon to knon
    DO i = 1, knon
       ig = knindex(i)
       champ_out(i) = temp_omp(ig)
    ENDDO

  END SUBROUTINE cpl2gath
!
!*************************************************************************************
!
  SUBROUTINE gath2cpl(champ_in, champ_out, knon, knindex)
  USE mod_phys_lmdz_para
! Cette routine ecrit un champ 'gathered' sur la grille 2D pour le passer
! au coupleur.
!
! input:         
!   champ_in     champ sur la grille gathere        
!   knon         nombre de points dans le domaine a traiter
!   knindex      index des points de la surface a traiter
!
! output:
!   champ_out    champ sur la grille 2D
!
    INCLUDE "dimensions.h"
    
! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                    :: knon
    REAL, DIMENSION(klon), INTENT(IN)      :: champ_in
    INTEGER, DIMENSION(klon), INTENT(IN)   :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(iim,jj_nb), INTENT(OUT) :: champ_out

! Local variables
!*************************************************************************************
    INTEGER                                :: i, ig
    REAL, DIMENSION(klon)                  :: temp_omp
    REAL, DIMENSION(klon_mpi)              :: temp_mpi
!*************************************************************************************

! Decompress from knon to klon
    temp_omp = 0.
    DO i = 1, knon
       ig = knindex(i)
       temp_omp(ig) = champ_in(i)
    ENDDO

! Transform from 1 dimension (klon) to 2 dimensions (iim,jj_nb)
    CALL gather_omp(temp_omp,temp_mpi)

!$OMP MASTER    
    CALL Grid1Dto2D_mpi(temp_mpi,champ_out)
    
    IF (is_north_pole) champ_out(:,1)=temp_mpi(1)
    IF (is_south_pole) champ_out(:,jj_nb)=temp_mpi(klon)
!$OMP END MASTER
    
  END SUBROUTINE gath2cpl
!
!*************************************************************************************
!
END MODULE cpl_mod

