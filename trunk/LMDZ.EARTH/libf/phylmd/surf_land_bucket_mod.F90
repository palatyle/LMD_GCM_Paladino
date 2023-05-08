!
MODULE surf_land_bucket_mod
!
! Surface land bucket module
!
! This module is used when no external land model is choosen.
!
  IMPLICIT NONE

CONTAINS

  SUBROUTINE surf_land_bucket(itime, jour, knon, knindex, debut, dtime,&
       tsurf, p1lay, tq_cdrag, precip_rain, precip_snow, temp_air, &
       spechum, petAcoef, peqAcoef, petBcoef, peqBcoef, pref, &
       u1, v1, rugoro, swnet, lwnet, &
       snow, qsol, agesno, tsoil, &
       qsurf, z0_new, alb1_new, alb2_new, evap, &
       fluxsens, fluxlat, tsurf_new, dflux_s, dflux_l)

    USE limit_read_mod
    USE surface_data
    USE fonte_neige_mod
    USE calcul_fluxs_mod
    USE cpl_mod
    USE dimphy
    USE mod_grid_phy_lmdz
    USE mod_phys_lmdz_para
!****************************************************************************************
! Bucket calculations for surface. 
!
    INCLUDE "clesphys.h"
    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"
    INCLUDE "YOMCST.h"

! Input variables  
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    LOGICAL, INTENT(IN)                     :: debut
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon), INTENT(IN)       :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)       :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)       :: tq_cdrag
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)       :: petAcoef, peqAcoef
    REAL, DIMENSION(klon), INTENT(IN)       :: petBcoef, peqBcoef
    REAL, DIMENSION(klon), INTENT(IN)       :: pref
    REAL, DIMENSION(klon), INTENT(IN)       :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)       :: rugoro
    REAL, DIMENSION(klon), INTENT(IN)       :: swnet, lwnet

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new, alb2_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon) :: soilcap, soilflux
    REAL, DIMENSION(klon) :: cal, beta, dif_grnd, capsol
    REAL, DIMENSION(klon) :: alb_neig, alb_lim
    REAL, DIMENSION(klon) :: zfra
    REAL, DIMENSION(klon) :: radsol       ! total net radiance at surface
    REAL, DIMENSION(klon) :: u0, v0, u1_lay, v1_lay
    REAL, DIMENSION(klon) :: dummy_riverflow,dummy_coastalflow 
    INTEGER               :: i
!
!****************************************************************************************


!
!* Read from limit.nc : albedo(alb_lim), length of rugosity(z0_new)
!
    CALL limit_read_rug_alb(itime, dtime, jour,&
         knon, knindex, &
         z0_new, alb_lim)
!
!* Calcultaion of fluxes 
!

! calculate total absorbed radiance at surface
       radsol(:) = 0.0
       radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

! calculate constants
    CALL calbeta(dtime, is_ter, knon, snow, qsol, beta, capsol, dif_grnd)
       
! calculate temperature, heat capacity and conduction flux in soil
    IF (soil_model) THEN 
       CALL soil(dtime, is_ter, knon, snow, tsurf, tsoil, soilcap, soilflux)
       DO i=1, knon
          cal(i) = RCPD / soilcap(i)
          radsol(i) = radsol(i)  + soilflux(i)
       END DO
    ELSE 
       cal(:) = RCPD * capsol(:)
    ENDIF
    
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

    CALL calcul_fluxs(knon, is_ter, dtime, &
         tsurf, p1lay, cal, beta, tq_cdrag, pref, &
         precip_rain, precip_snow, snow, qsurf,  &
         radsol, dif_grnd, temp_air, spechum, u1_lay, v1_lay, &
         petAcoef, peqAcoef, petBcoef, peqBcoef, &
         tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)
    
!
!* Calculate snow height, run_off, age of snow
!      
    CALL fonte_neige( knon, is_ter, knindex, dtime, &
         tsurf, precip_rain, precip_snow, &
         snow, qsol, tsurf_new, evap)
!
!* Calculate the age of snow
!
    CALL albsno(klon,knon,dtime,agesno(:),alb_neig(:), precip_snow(:))  
    
    WHERE (snow(1 : knon) .LT. 0.0001) agesno(1 : knon) = 0.
    
    DO i=1, knon
       zfra(i) = MAX(0.0,MIN(1.0, snow(i)/(snow(i)+10.0)))
       alb_lim(i)  = alb_neig(i) *zfra(i) + alb_lim(i)*(1.0-zfra(i))
    END DO

!
!* Return albedo : 
!    alb1_new and alb2_new are here given the same values
!
    alb1_new(:) = 0.0
    alb2_new(:) = 0.0
    alb1_new(1:knon) = alb_lim(1:knon)
    alb2_new(1:knon) = alb_lim(1:knon)
       
!
!* Calculate the rugosity
!
    DO i = 1, knon
       z0_new(i) = MAX(1.5e-05,SQRT(z0_new(i)**2 + rugoro(i)**2))
    END DO

!* Send to coupler
!  The run-off from river and coast are not calculated in the bucket modele.
!  For testing purpose of the coupled modele we put the run-off to zero.
    IF (type_ocean=='couple') THEN
       dummy_riverflow(:)   = 0.0
       dummy_coastalflow(:) = 0.0
       CALL cpl_send_land_fields(itime, knon, knindex, &
            dummy_riverflow, dummy_coastalflow)
    ENDIF

!
!* End
!
  END SUBROUTINE surf_land_bucket
!
!****************************************************************************************
!
END MODULE surf_land_bucket_mod
