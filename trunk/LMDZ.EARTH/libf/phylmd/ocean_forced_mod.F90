!
MODULE ocean_forced_mod
!
! This module is used for both the sub-surfaces ocean and sea-ice for the case of a 
! forced ocean,  "ocean=force".
!
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE ocean_forced_noice( &
       itime, dtime, jour, knon, knindex, &
       p1lay, cdragh, cdragm, precip_rain, precip_snow, &
       temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, &
       radsol, snow, agesno, & 
       qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l)
!
! This subroutine treats the "open ocean", all grid points that are not entierly covered
! by ice.
! The routine receives data from climatologie file limit.nc and does some calculations at the 
! surface. 
!
    USE dimphy
    USE calcul_fluxs_mod
    USE limit_read_mod
    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    REAL, INTENT(IN)                         :: dtime
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
    REAL, DIMENSION(klon), INTENT(INOUT)     :: agesno !? put to 0 in ocean
  
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l

! Local variables
!****************************************************************************************
    INTEGER                     :: i
    REAL, DIMENSION(klon)       :: cal, beta, dif_grnd
    REAL, DIMENSION(klon)       :: alb_neig, tsurf_lim, zx_sl
    REAL, DIMENSION(klon)       :: u0, v0
    REAL, DIMENSION(klon)       :: u1_lay, v1_lay
    LOGICAL                     :: check=.FALSE.

!****************************************************************************************
! Start calculation
!****************************************************************************************
    IF (check) WRITE(*,*)' Entering ocean_forced_noice'
    
!****************************************************************************************
! 1)    
! Read sea-surface temperature from file limit.nc
!
!****************************************************************************************
    CALL limit_read_sst(knon,knindex,tsurf_lim)

!****************************************************************************************
! 2)
! Flux calculation
!
!****************************************************************************************
! Set some variables for calcul_fluxs
    cal = 0.
    beta = 1.
    dif_grnd = 0.
    alb_neig(:) = 0.
    agesno(:) = 0.
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

! Calcul de tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l and qsurf
    CALL calcul_fluxs(knon, is_oce, dtime, &
         tsurf_lim, p1lay, cal, beta, cdragh, ps, &
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

  END SUBROUTINE ocean_forced_noice
!
!***************************************************************************************
!
  SUBROUTINE ocean_forced_ice( &
       itime, dtime, jour, knon, knindex, &
       tsurf_in, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, &
       radsol, snow, qsol, agesno, tsoil, &
       qsurf, alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l)
!
! This subroutine treats the ocean where there is ice. 
! The routine reads data from climatologie file and does flux calculations at the 
! surface.
!
    USE dimphy
    USE calcul_fluxs_mod
    USE surface_data,     ONLY : calice, calsno, tau_gl
    USE limit_read_mod
    USE fonte_neige_mod,  ONLY : fonte_neige

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "clesphys.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex
    REAL, INTENT(IN)                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf_in
    REAL, DIMENSION(klon), INTENT(IN)    :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)    :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)    :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)    :: ps
    REAL, DIMENSION(klon), INTENT(IN)    :: u1, v1

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: radsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)            :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb1_new  ! new albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)            :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)            :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)            :: dflux_s, dflux_l      

! Local variables
!****************************************************************************************
    LOGICAL                     :: check=.FALSE.
    INTEGER                     :: i
    REAL                        :: zfra
    REAL, PARAMETER             :: t_grnd=271.35
    REAL, DIMENSION(klon)       :: cal, beta, dif_grnd, capsol
    REAL, DIMENSION(klon)       :: alb_neig, tsurf_tmp
    REAL, DIMENSION(klon)       :: soilcap, soilflux
    REAL, DIMENSION(klon)       :: u0, v0
    REAL, DIMENSION(klon)       :: u1_lay, v1_lay

!****************************************************************************************
! Start calculation
!****************************************************************************************
    IF (check) WRITE(*,*)'Entering surface_seaice, knon=',knon 

!****************************************************************************************
! 1) 
! Flux calculation : tsurf_new, evap, fluxlat, fluxsens, flux_u1, flux_v1
!                    dflux_s, dflux_l and qsurf
!****************************************************************************************
    tsurf_tmp(:) = tsurf_in(:)

! calculate the parameters cal, beta, capsol and dif_grnd
    CALL calbeta(dtime, is_sic, knon, snow, qsol, beta, capsol, dif_grnd)

    
    IF (soil_model) THEN 
! update tsoil and calculate soilcap and soilflux
       CALL soil(dtime, is_sic, knon, snow, tsurf_tmp, tsoil,soilcap, soilflux)
       cal(1:knon) = RCPD / soilcap(1:knon)
       radsol(1:knon) = radsol(1:knon)  + soilflux(1:knon)
       dif_grnd = 1.0 / tau_gl
    ELSE 
       dif_grnd = 1.0 / tau_gl
       cal = RCPD * calice
       WHERE (snow > 0.0) cal = RCPD * calsno 
    ENDIF

    beta = 1.0
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)
    CALL calcul_fluxs(knon, is_sic, dtime, &
         tsurf_tmp, p1lay, cal, beta, cdragh, ps, &
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
! 2)
! Calculations due to snow and runoff
!
!****************************************************************************************
    CALL fonte_neige( knon, is_sic, knindex, dtime, &
         tsurf_tmp, precip_rain, precip_snow, &
         snow, qsol, tsurf_new, evap)
    
! Calculation of albedo at snow (alb_neig) and update the age of snow (agesno)
! 
    CALL albsno(klon, knon, dtime, agesno(:), alb_neig(:), precip_snow(:))  

    WHERE (snow(1:knon) .LT. 0.0001) agesno(1:knon) = 0.

    alb1_new(:) = 0.0
    DO i=1, knon
       zfra = MAX(0.0,MIN(1.0,snow(i)/(snow(i)+10.0)))
       alb1_new(i) = alb_neig(i) * zfra +  0.6 * (1.0-zfra)
    ENDDO

    alb2_new(:) = alb1_new(:)

  END SUBROUTINE ocean_forced_ice
!
!****************************************************************************************
!
END MODULE ocean_forced_mod






