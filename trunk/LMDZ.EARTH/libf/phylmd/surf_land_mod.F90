!
MODULE surf_land_mod
  
  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!  
  SUBROUTINE surf_land(itime, dtime, date0, jour, knon, knindex, &
       rlon, rlat, &
       debut, lafin, zlev, ccanopy, swnet, lwnet, albedo, &
       tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, & 
       pref, u1, v1, rugoro, pctsrf, &
       lwdown_m, q2m, t2m, &
       snow, qsol, agesno, tsoil, &
       z0_new, alb1_new, alb2_new, evap, fluxsens, fluxlat, &
       qsurf, tsurf_new, dflux_s, dflux_l, &
       flux_u1, flux_v1 ) 

    USE dimphy
    USE surface_data, ONLY    : ok_veget

#ifdef ORCHIDEE_NOOPENMP
    USE surf_land_orchidee_noopenmp_mod
#else
    USE surf_land_orchidee_mod
#endif
    USE surf_land_bucket_mod
    USE calcul_fluxs_mod

    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"

! Input variables  
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, INTENT(IN)                        :: date0
    REAL, DIMENSION(klon), INTENT(IN)       :: rlon, rlat
    LOGICAL, INTENT(IN)                     :: debut, lafin
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon), INTENT(IN)       :: zlev, ccanopy
    REAL, DIMENSION(klon), INTENT(IN)       :: swnet, lwnet
    REAL, DIMENSION(klon), INTENT(IN)       :: albedo  ! albedo for whole short-wave interval
    REAL, DIMENSION(klon), INTENT(IN)       :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)       :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)       :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)       :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)       :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)       :: pref   ! pressure reference
    REAL, DIMENSION(klon), INTENT(IN)       :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)       :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf
    REAL, DIMENSION(klon), INTENT(IN)       :: lwdown_m  ! downwelling longwave radiation at mean surface
                                                         ! corresponds to previous sollwdown
    REAL, DIMENSION(klon), INTENT(IN)       :: q2m, t2m

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new ! albdeo for shortwave interval 1(visible)
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new ! albedo for shortwave interval 2(near infrared)
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap
    REAL, DIMENSION(klon), INTENT(OUT)       :: fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1  ! flux for U and V at first model level

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon) :: p1lay_tmp
    REAL, DIMENSION(klon) :: pref_tmp
    REAL, DIMENSION(klon) :: swdown     ! downwelling shortwave radiation at land surface
    REAL, DIMENSION(klon) :: lwdown     ! downwelling longwave radiation at land surface
    REAL, DIMENSION(klon) :: epot_air           ! potential air temperature
    REAL, DIMENSION(klon) :: tsol_rad, emis_new ! output from interfsol not used
    REAL, DIMENSION(klon) :: u0, v0     ! surface speed
    INTEGER               :: i


!**************************************************************************************** 
! Choice between call to vegetation model (ok_veget=true) or simple calculation below
!
!****************************************************************************************
   IF (ok_veget) THEN
!****************************************************************************************
!  Call model sechiba in model ORCHIDEE
!
!****************************************************************************************
       p1lay_tmp(:)      = 0.0
       pref_tmp(:)       = 0.0
       p1lay_tmp(1:knon) = p1lay(1:knon)/100.
       pref_tmp(1:knon)  = pref(1:knon)/100.
! 
!* Calculate incoming flux for SW and LW interval: swdown, lwdown
!
       swdown(:) = 0.0
       lwdown(:) = 0.0
       DO i = 1, knon
          swdown(i) = swnet(i)/(1-albedo(i))
          lwdown(i) = lwnet(i) + RSIGMA*tsurf(i)**4
       END DO
!
!* Calculate potential air temperature
!
       epot_air(:) = 0.0
       DO i = 1, knon
          epot_air(i) = RCPD*temp_air(i)*(pref(i)/p1lay(i))**RKAPPA
       END DO

       ! temporary for keeping same results using lwdown_m instead of lwdown
       CALL surf_land_orchidee(itime, dtime, date0, knon, &
            knindex, rlon, rlat, pctsrf, &
            debut, lafin, &
            zlev,  u1, v1, temp_air, spechum, epot_air, ccanopy, & 
            cdragh, AcoefH, AcoefQ, BcoefH, BcoefQ, &
            precip_rain, precip_snow, lwdown_m, swnet, swdown, &
            pref_tmp, q2m, t2m, &
            evap, fluxsens, fluxlat, &              
            tsol_rad, tsurf_new, alb1_new, alb2_new, &
            emis_new, z0_new, qsurf)       

!  
!* Add contribution of relief to surface roughness
!  
       DO i=1,knon
          z0_new(i) = MAX(1.5e-05,SQRT(z0_new(i)**2 + rugoro(i)**2))
       ENDDO

    ELSE  ! not ok_veget
!****************************************************************************************
! No extern vegetation model choosen, call simple bucket calculations instead.
!
!****************************************************************************************
       CALL surf_land_bucket(itime, jour, knon, knindex, debut, dtime,&
            tsurf, p1lay, cdragh, precip_rain, precip_snow, temp_air, &
            spechum, AcoefH, AcoefQ, BcoefH, BcoefQ, pref, &
            u1, v1, rugoro, swnet, lwnet, &
            snow, qsol, agesno, tsoil, &
            qsurf, z0_new, alb1_new, alb2_new, evap, &
            fluxsens, fluxlat, tsurf_new, dflux_s, dflux_l)

    ENDIF ! ok_veget

!****************************************************************************************
! Calculation for all land models
! - Flux calculation at first modele level for U and V
!****************************************************************************************
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    CALL calcul_flux_wind(knon, dtime, &
         u0, v0, u1, v1, cdragm, &
         AcoefU, AcoefV, BcoefU, BcoefV, &
         p1lay, temp_air, &
         flux_u1, flux_v1)
    
  END SUBROUTINE surf_land
!
!****************************************************************************************
!  
END MODULE surf_land_mod
!
!****************************************************************************************
!  
