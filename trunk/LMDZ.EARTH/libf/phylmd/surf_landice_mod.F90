!
MODULE surf_landice_mod
  
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE surf_landice(itime, dtime, knon, knindex, &
       swnet, lwnet, tsurf, p1lay, &
       cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, rugoro, pctsrf, &
       snow, qsurf, qsol, agesno, &
       tsoil, z0_new, alb1, alb2, evap, fluxsens, fluxlat, &
       tsurf_new, dflux_s, dflux_l, &
       flux_u1, flux_v1)

    USE dimphy
    USE surface_data,     ONLY : type_ocean, calice, calsno
    USE fonte_neige_mod,  ONLY : fonte_neige, run_off_lic
    USE cpl_mod,          ONLY : cpl_send_landice_fields
    USE calcul_fluxs_mod
    USE phys_output_var_mod

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"
    INCLUDE "clesphys.h"

! Input variables 
!****************************************************************************************
    INTEGER, INTENT(IN)                           :: itime, knon
    INTEGER, DIMENSION(klon), INTENT(in)          :: knindex
    REAL, INTENT(in)                              :: dtime
    REAL, DIMENSION(klon), INTENT(IN)             :: swnet ! net shortwave radiance
    REAL, DIMENSION(klon), INTENT(IN)             :: lwnet ! net longwave radiance
    REAL, DIMENSION(klon), INTENT(IN)             :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)             :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)             :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)             :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)             :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)             :: AcoefH, AcoefQ
    REAL, DIMENSION(klon), INTENT(IN)             :: BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)             :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)             :: ps
    REAL, DIMENSION(klon), INTENT(IN)             :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)             :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)       :: pctsrf

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)            :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)            :: z0_new
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb1  ! new albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: alb2  ! new albedo in near IR interval
    REAL, DIMENSION(klon), INTENT(OUT)            :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)            :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)            :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)            :: flux_u1, flux_v1

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon)    :: soilcap, soilflux
    REAL, DIMENSION(klon)    :: cal, beta, dif_grnd
    REAL, DIMENSION(klon)    :: zfra, alb_neig
    REAL, DIMENSION(klon)    :: radsol
    REAL, DIMENSION(klon)    :: u0, v0, u1_lay, v1_lay
    INTEGER                  :: i,j

! End definition
!****************************************************************************************
!
! Initialize output variables
    alb2(:) = 999999.
    alb1(:) = 999999.

!****************************************************************************************
! Calculate total absorbed radiance at surface
!
!****************************************************************************************
    radsol(:) = 0.0
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

!****************************************************************************************
! Soil calculations
! 
!****************************************************************************************
    IF (soil_model) THEN 
       CALL soil(dtime, is_lic, knon, snow, tsurf, tsoil, soilcap, soilflux)
       cal(1:knon) = RCPD / soilcap(1:knon)
       radsol(1:knon)  = radsol(1:knon) + soilflux(1:knon)
    ELSE 
       cal = RCPD * calice
       WHERE (snow > 0.0) cal = RCPD * calsno
    ENDIF


!****************************************************************************************
! Calulate fluxes
!
!****************************************************************************************
    beta(:) = 1.0
    dif_grnd(:) = 0.0

! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

    CALL calcul_fluxs(knon, is_lic, dtime, &
         tsurf, p1lay, cal, beta, cdragh, ps, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
         AcoefH, AcoefQ, BcoefH, BcoefQ, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)

    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)

!****************************************************************************************
! Calculate snow height, age, run-off,..
!    
!****************************************************************************************
    CALL fonte_neige( knon, is_lic, knindex, dtime, &
         tsurf, precip_rain, precip_snow, &
         snow, qsol, tsurf_new, evap)


!****************************************************************************************
! Calculate albedo
!
!****************************************************************************************
    CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
    WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
    zfra(1:knon) = MAX(0.0,MIN(1.0,snow(1:knon)/(snow(1:knon)+10.0)))
    alb1(1:knon) = alb_neig(1:knon)*zfra(1:knon) + &
         0.6 * (1.0-zfra(1:knon))
!
!IM: plusieurs choix/tests sur l'albedo des "glaciers continentaux"
!       alb1(1 : knon)  = 0.6 !IM cf FH/GK 
!       alb1(1 : knon)  = 0.82
!       alb1(1 : knon)  = 0.77 !211003 Ksta0.77
!       alb1(1 : knon)  = 0.8 !KstaTER0.8 & LMD_ARMIP5
!IM: KstaTER0.77 & LMD_ARMIP6    

! Attantion: alb1 and alb2 are the same!
    alb1(1:knon)  = 0.77
    alb2(1:knon)  = alb1(1:knon)


!****************************************************************************************
! Rugosity
!
!****************************************************************************************
    z0_new(:) = MAX(1.E-3,rugoro(:))

!****************************************************************************************
! Send run-off on land-ice to coupler if coupled ocean.
! run_off_lic has been calculated in fonte_neige
!
!****************************************************************************************
    IF (type_ocean=='couple') THEN
       CALL cpl_send_landice_fields(itime, knon, knindex, run_off_lic)
    ENDIF
  
!****************************************************************************************
       snow_o=0.
       zfra_o = 0.
       DO j = 1, knon
           i = knindex(j)
           snow_o(i) = snow(j) 
           zfra_o(i) = zfra(j)
       ENDDO

!****************************************************************************************
       snow_o=0.
       zfra_o = 0.
       DO j = 1, knon
           i = knindex(j)
           snow_o(i) = snow(j)
           zfra_o(i) = zfra(j)
       ENDDO


  END SUBROUTINE surf_landice
!
!****************************************************************************************
!
END MODULE surf_landice_mod



