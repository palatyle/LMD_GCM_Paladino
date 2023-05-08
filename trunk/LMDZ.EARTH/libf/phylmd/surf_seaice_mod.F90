!
MODULE surf_seaice_mod

  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE surf_seaice( & 
       rlon, rlat, swnet, lwnet, alb1, fder, &
       itime, dtime, jour, knon, knindex, &
       lafin, &
       tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, rugoro, pctsrf, &
       snow, qsurf, qsol, agesno, tsoil, &
       z0_new, alb1_new, alb2_new, evap, fluxsens, fluxlat, &
       tsurf_new, dflux_s, dflux_l, &
       flux_u1, flux_v1)

  USE dimphy
  USE surface_data
  USE ocean_forced_mod, ONLY : ocean_forced_ice
  USE ocean_cpl_mod, ONLY    : ocean_cpl_ice

!
! This subroutine will make a call to ocean_XXX_ice according to the ocean mode (force, 
! slab or couple). The calculation of rugosity for the sea-ice surface is also done
! in here because it is the same calculation for the different modes of ocean.
!
    INCLUDE "indicesol.h"
    INCLUDE "dimsoil.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    LOGICAL, INTENT(IN)                      :: lafin
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet  ! net shortwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet  ! net longwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: fder
    REAL, DIMENSION(klon), INTENT(IN)        :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)        :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)  :: pctsrf

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)          :: snow, qsurf, qsol
    REAL, DIMENSION(klon), INTENT(INOUT)          :: agesno
    REAL, DIMENSION(klon, nsoilmx), INTENT(INOUT) :: tsoil

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new  ! new albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1

! Local arguments
!****************************************************************************************
    REAL, DIMENSION(klon)  :: radsol

!
! End definitions
!****************************************************************************************


!****************************************************************************************
! Calculate total net radiance at surface
!
!****************************************************************************************
    radsol(:) = 0.0
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

!****************************************************************************************
! Switch according to type of ocean (couple, slab or forced)
!
!****************************************************************************************
    IF (type_ocean == 'couple') THEN
       
       CALL ocean_cpl_ice( &
            rlon, rlat, swnet, lwnet, alb1, & 
            fder, & 
            itime, dtime, knon, knindex, &
            lafin,&
            p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, pctsrf, &
            radsol, snow, qsurf, &
            alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)
       
    ELSE IF (type_ocean == 'force' .OR. (type_ocean == 'slab' .AND. version_ocean=='sicOBS')) THEN
       CALL ocean_forced_ice( &
            itime, dtime, jour, knon, knindex, &
            tsurf, p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, &
            radsol, snow, qsol, agesno, tsoil, &
            qsurf, alb1_new, alb2_new, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)

    ELSE IF (type_ocean == 'slab') THEN
!!$       CALL ocean_slab_ice( & 
!!$          itime, dtime, jour, knon, knindex, &
!!$          debut, &
!!$          tsurf, p1lay, cdragh, precip_rain, precip_snow, temp_air, spechum,&
!!$          AcoefH, AcoefQ, BcoefH, BcoefQ, &
!!$          ps, u1, v1, pctsrf, &
!!$          radsol, snow, qsurf, qsol, agesno, tsoil, &
!!$          alb1_new, alb2_new, evap, fluxsens, fluxlat, &
!!$          tsurf_new, dflux_s, dflux_l)

    END IF

!****************************************************************************************
! Calculate rugosity
!
!****************************************************************************************
    z0_new = 0.002
    z0_new = SQRT(z0_new**2+rugoro**2)

  END SUBROUTINE surf_seaice
!
!****************************************************************************************
!
END MODULE surf_seaice_mod

