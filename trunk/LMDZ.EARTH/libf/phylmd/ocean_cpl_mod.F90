!
MODULE ocean_cpl_mod
!
! This module is used both for the sub-surface ocean and sea-ice for the case of a 
! coupled model configuration, ocean=couple. 
!

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ocean_cpl_init, ocean_cpl_noice, ocean_cpl_ice

!****************************************************************************************
!
CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE ocean_cpl_init(dtime, rlon, rlat)
!
! Allocate fields for this module and initailize the module mod_cpl
!
    USE dimphy,           ONLY : klon
    USE cpl_mod

! Input arguments
!*************************************************************************************
    REAL, INTENT(IN)                  :: dtime
    REAL, DIMENSION(klon), INTENT(IN) :: rlon, rlat

! Local variables
!*************************************************************************************
    INTEGER              :: error
    CHARACTER (len = 80) :: abort_message
    CHARACTER (len = 20) :: modname = 'ocean_cpl_init'

! Initialize module cpl_init
    CALL cpl_init(dtime, rlon, rlat)
    
  END SUBROUTINE ocean_cpl_init
!
!****************************************************************************************
!
  SUBROUTINE ocean_cpl_noice( &
       swnet, lwnet, alb1, &
       windsp, fder_old, &
       itime, dtime, knon, knindex, &
       p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, &
       radsol, snow, agesno, &
       qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l)

!
! This subroutine treats the "open ocean", all grid points that are not entierly covered
! by ice. The subroutine first receives fields from coupler, then some calculations at 
! surface is done and finally it sends some fields to the coupler.
!
    USE dimphy,           ONLY : klon
    USE cpl_mod
    USE calcul_fluxs_mod

    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"
!    
! Input arguments  
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: windsp
    REAL, DIMENSION(klon), INTENT(IN)        :: fder_old
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)     :: radsol
    REAL, DIMENSION(klon), INTENT(INOUT)     :: snow
    REAL, DIMENSION(klon), INTENT(INOUT)     :: agesno
  
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      

! Local variables
!****************************************************************************************
    INTEGER               :: i
    INTEGER, DIMENSION(1) :: iloc
    REAL, DIMENSION(klon) :: cal, beta, dif_grnd
    REAL, DIMENSION(klon) :: fder_new
    REAL, DIMENSION(klon) :: tsurf_cpl
    REAL, DIMENSION(klon) :: u0_cpl, v0_cpl
    REAL, DIMENSION(klon) :: u1_lay, v1_lay
    LOGICAL               :: check=.FALSE.

! End definitions
!****************************************************************************************

    IF (check) WRITE(*,*)' Entering ocean_cpl_noice'

!****************************************************************************************
! Receive sea-surface temperature(tsurf_cpl) from coupler
!
!****************************************************************************************
    CALL cpl_receive_ocean_fields(knon, knindex, tsurf_cpl, u0_cpl, v0_cpl)

!****************************************************************************************
! Calculate fluxes at surface
!
!****************************************************************************************
    cal = 0.
    beta = 1.
    dif_grnd = 0.
    agesno(:) = 0.

    DO i = 1, knon
       u1_lay(i) = u1(i) - u0_cpl(i)
       v1_lay(i) = v1(i) - v0_cpl(i)
    END DO

    CALL calcul_fluxs(knon, is_oce, dtime, &
         tsurf_cpl, p1lay, cal, beta, cdragh, ps, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
         AcoefH, AcoefQ, BcoefH, BcoefQ, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)
    
! - Flux calculation at first modele level for U and V
    CALL calcul_flux_wind(knon, dtime, &
         u0_cpl, v0_cpl, u1, v1, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)  

!****************************************************************************************
! Calculate fder : flux derivative (sensible and latente)
!
!****************************************************************************************
    fder_new(:) = fder_old(:) + dflux_s(:) + dflux_l(:)
    
    iloc = MAXLOC(fder_new(1:klon))
    IF (check .AND. fder_new(iloc(1))> 0.) THEN
       WRITE(*,*)'**** Debug fder****'
       WRITE(*,*)'max fder(',iloc(1),') = ',fder_new(iloc(1))
       WRITE(*,*)'fder_old, dflux_s, dflux_l',fder_old(iloc(1)), &
            dflux_s(iloc(1)), dflux_l(iloc(1))
    ENDIF

!****************************************************************************************
! Send and cumulate fields to the coupler
!
!****************************************************************************************

    CALL cpl_send_ocean_fields(itime, knon, knindex, &
         swnet, lwnet, fluxlat, fluxsens, &
         precip_rain, precip_snow, evap, tsurf_new, fder_new, alb1, flux_u1, flux_v1, windsp)
    

  END SUBROUTINE ocean_cpl_noice
!
!****************************************************************************************
!
  SUBROUTINE ocean_cpl_ice( &
       rlon, rlat, swnet, lwnet, alb1, &
       fder_old, &
       itime, dtime, knon, knindex, &
       lafin, &
       p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, pctsrf, &
       radsol, snow, qsurf, &
       alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l)
!
! This subroutine treats the ocean where there is ice. The subroutine first receives 
! fields from coupler, then some calculations at surface is done and finally sends 
! some fields to the coupler.
!    
    USE dimphy,           ONLY : klon
    USE cpl_mod
    USE calcul_fluxs_mod

    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    LOGICAL, INTENT(IN)                      :: lafin
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: fder_old
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)  :: pctsrf

! In/output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)     :: radsol
    REAL, DIMENSION(klon), INTENT(INOUT)     :: snow

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new, alb2_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      

! Local variables
!****************************************************************************************
    INTEGER                 :: i
    INTEGER, DIMENSION(1)   :: iloc
    LOGICAL                 :: check=.FALSE.
    REAL, PARAMETER         :: t_grnd=271.35
    REAL, DIMENSION(klon)   :: cal, beta, dif_grnd
    REAL, DIMENSION(klon)   :: tsurf_cpl, fder_new
    REAL, DIMENSION(klon)   :: alb_cpl
    REAL, DIMENSION(klon)   :: u0, v0
    REAL, DIMENSION(klon)   :: u1_lay, v1_lay

! End definitions
!****************************************************************************************
    
    IF (check) WRITE(*,*)'Entering surface_seaice, knon=',knon 

!****************************************************************************************
! Receive ocean temperature(tsurf_cpl) and albedo(alb_new) from coupler
!
!****************************************************************************************

    CALL cpl_receive_seaice_fields(knon, knindex, &
         tsurf_cpl, alb_cpl, u0, v0)

    alb1_new(1:knon) = alb_cpl(1:knon)
    alb2_new(1:knon) = alb_cpl(1:knon)    

    
!****************************************************************************************
! Calculate fluxes at surface
!
!****************************************************************************************
    cal = 0.
    dif_grnd = 0.
    beta = 1.0
    
    DO i = 1, knon
       u1_lay(i) = u1(i) - u0(i)
       v1_lay(i) = v1(i) - v0(i)
    END DO

    CALL calcul_fluxs(knon, is_sic, dtime, &
         tsurf_cpl, p1lay, cal, beta, cdragh, ps, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
         AcoefH, AcoefQ, BcoefH, BcoefQ, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)


! - Flux calculation at first modele level for U and V
    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)  

!****************************************************************************************
! Calculate fder : flux derivative (sensible and latente)
!
!****************************************************************************************
    fder_new(:) = fder_old(:) + dflux_s(:) + dflux_l(:)
    
    iloc = MAXLOC(fder_new(1:klon))
    IF (check .AND. fder_new(iloc(1))> 0.) THEN
       WRITE(*,*)'**** Debug fder ****'
       WRITE(*,*)'max fder(',iloc(1),') = ',fder_new(iloc(1))
       WRITE(*,*)'fder_old, dflux_s, dflux_l',fder_old(iloc(1)), &
            dflux_s(iloc(1)), dflux_l(iloc(1))
    ENDIF

!****************************************************************************************
! Send and cumulate fields to the coupler
!
!****************************************************************************************

    CALL cpl_send_seaice_fields(itime, dtime, knon, knindex, &
       pctsrf, lafin, rlon, rlat, &
       swnet, lwnet, fluxlat, fluxsens, &
       precip_rain, precip_snow, evap, tsurf_new, fder_new, alb1, flux_u1, flux_v1)
 

  END SUBROUTINE ocean_cpl_ice
!  
!****************************************************************************************
!
END MODULE ocean_cpl_mod
