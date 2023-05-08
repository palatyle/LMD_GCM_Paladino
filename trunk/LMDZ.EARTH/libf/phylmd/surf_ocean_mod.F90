!
MODULE surf_ocean_mod

  IMPLICIT NONE

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE surf_ocean(rlon, rlat, swnet, lwnet, alb1, &
       rugos, windsp, rmu0, fder, tsurf_in, &
       itime, dtime, jour, knon, knindex, &
       p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, rugoro, pctsrf, &
       snow, qsurf, agesno, &
       z0_new, alb1_new, alb2_new, evap, fluxsens, fluxlat, &
       tsurf_new, dflux_s, dflux_l, lmt_bils, &
       flux_u1, flux_v1)

  USE dimphy
  USE surface_data, ONLY     : type_ocean
  USE ocean_forced_mod, ONLY : ocean_forced_noice
  USE ocean_slab_mod, ONLY   : ocean_slab_noice
  USE ocean_cpl_mod, ONLY    : ocean_cpl_noice
!
! This subroutine will make a call to ocean_XXX_noice according to the ocean mode (force, 
! slab or couple). The calculations of albedo and rugosity for the ocean surface are 
! done in here because they are identical for the different modes of ocean. 
!
    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

! Input variables
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: itime, jour, knon
    INTEGER, DIMENSION(klon), INTENT(IN)     :: knindex
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon), INTENT(IN)        :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)        :: swnet  ! net shortwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: lwnet  ! net longwave radiation at surface  
    REAL, DIMENSION(klon), INTENT(IN)        :: alb1   ! albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(IN)        :: rugos
    REAL, DIMENSION(klon), INTENT(IN)        :: windsp
    REAL, DIMENSION(klon), INTENT(IN)        :: rmu0  
    REAL, DIMENSION(klon), INTENT(IN)        :: fder
    REAL, DIMENSION(klon), INTENT(IN)        :: tsurf_in
    REAL, DIMENSION(klon), INTENT(IN)        :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragh
    REAL, DIMENSION(klon), INTENT(IN)        :: cdragm
    REAL, DIMENSION(klon), INTENT(IN)        :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)        :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)        :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)        :: ps
    REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)        :: rugoro
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)  :: pctsrf

! In/Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT)     :: snow
    REAL, DIMENSION(klon), INTENT(INOUT)     :: qsurf
    REAL, DIMENSION(klon), INTENT(INOUT)     :: agesno

! Output variables
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: z0_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb1_new  ! new albedo in visible SW interval
    REAL, DIMENSION(klon), INTENT(OUT)       :: alb2_new  ! new albedo in near IR interval
    REAL, DIMENSION(klon), INTENT(OUT)       :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)       :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)       :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)       :: lmt_bils
    REAL, DIMENSION(klon), INTENT(OUT)       :: flux_u1, flux_v1

! Local variables
!****************************************************************************************
    INTEGER               :: i
    REAL                  :: tmp
    REAL, PARAMETER       :: cepdu2=(0.1)**2
    REAL, DIMENSION(klon) :: alb_eau
    REAL, DIMENSION(klon) :: radsol

! End definition
!****************************************************************************************


!****************************************************************************************
! Calculate total net radiance at surface
!
!****************************************************************************************
    radsol(:) = 0.0
    radsol(1:knon) = swnet(1:knon) + lwnet(1:knon)

!****************************************************************************************
! Switch according to type of ocean (couple, slab or forced)
!****************************************************************************************
    SELECT CASE(type_ocean)
    CASE('couple')
       CALL ocean_cpl_noice( &
            swnet, lwnet, alb1, &
            windsp, fder, & 
            itime, dtime, knon, knindex, &
            p1lay, cdragh, cdragm, precip_rain, precip_snow,temp_air,spechum,& 
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, &
            radsol, snow, agesno, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)

    CASE('slab')
       CALL ocean_slab_noice( &
            itime, dtime, jour, knon, knindex, &
            p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum,&
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, tsurf_in, &
            radsol, snow, agesno, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l, lmt_bils)
       
    CASE('force')
       CALL ocean_forced_noice( &
            itime, dtime, jour, knon, knindex, &
            p1lay, cdragh, cdragm, precip_rain, precip_snow, &
            temp_air, spechum, &
            AcoefH, AcoefQ, BcoefH, BcoefQ, &
            AcoefU, AcoefV, BcoefU, BcoefV, &
            ps, u1, v1, &
            radsol, snow, agesno, &
            qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
            tsurf_new, dflux_s, dflux_l)
    END SELECT

!****************************************************************************************
! Calculate albedo
!
!****************************************************************************************
    IF ( MINVAL(rmu0) == MAXVAL(rmu0) .AND. MINVAL(rmu0) == -999.999 ) THEN
       CALL alboc(REAL(jour),rlat,alb_eau)
    ELSE  ! diurnal cycle
       CALL alboc_cd(rmu0,alb_eau)
    ENDIF

    DO i =1, knon
       alb1_new(i) = alb_eau(knindex(i))
    ENDDO
    alb2_new(1:knon) = alb1_new(1:knon)

!****************************************************************************************
! Calculate the rugosity
!
!****************************************************************************************
    DO i = 1, knon
       tmp = MAX(cepdu2,u1(i)**2+v1(i)**2)
       z0_new(i) = 0.018*cdragm(i) * (u1(i)**2+v1(i)**2)/RG  &
            +  0.11*14e-6 / SQRT(cdragm(i) * tmp)
       z0_new(i) = MAX(1.5e-05,z0_new(i))
    ENDDO   
!
!****************************************************************************************
!    
  END SUBROUTINE surf_ocean
!
!****************************************************************************************
!
END MODULE surf_ocean_mod
