!
MODULE ocean_slab_mod
!
! This module is used for both surface ocean and sea-ice when using the slab ocean,
! "ocean=slab".
!
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ocean_slab_frac, ocean_slab_noice

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_frac(itime, dtime, jour, pctsrf, is_modified)

    USE dimphy
    USE limit_read_mod
    USE surface_data
    INCLUDE "indicesol.h"
!    INCLUDE "clesphys.h"

! Arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                        :: itime   ! numero du pas de temps courant
    INTEGER, INTENT(IN)                        :: jour    ! jour a lire dans l'annee
    REAL   , INTENT(IN)                        :: dtime   ! pas de temps de la physique (en s)
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: pctsrf  ! sub-surface fraction
    LOGICAL, INTENT(OUT)                       :: is_modified ! true if pctsrf is modified at this time step

! Local variables
!****************************************************************************************
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'ocean_slab_frac'


    IF (version_ocean == 'sicOBS') THEN   
       CALL limit_read_frac(itime, dtime, jour, pctsrf, is_modified)
    ELSE
       abort_message='Ocean slab model without forced sea-ice fractions has to be rewritten!!!'
       CALL abort_gcm(modname,abort_message,1)
! Here should sea-ice/ocean fraction either be calculated or returned if saved as a module varaiable 
! (in the case the new fractions are calculated in ocean_slab_ice or ocean_slab_noice subroutines).  
    END IF

  END SUBROUTINE ocean_slab_frac
!
!****************************************************************************************
!
  SUBROUTINE ocean_slab_noice( & 
       itime, dtime, jour, knon, knindex, &
       p1lay, cdragh, cdragm, precip_rain, precip_snow, temp_air, spechum, &
       AcoefH, AcoefQ, BcoefH, BcoefQ, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       ps, u1, v1, tsurf_in, &
       radsol, snow, agesno, &
       qsurf, evap, fluxsens, fluxlat, flux_u1, flux_v1, &
       tsurf_new, dflux_s, dflux_l, lmt_bils)
    
    USE dimphy
    USE calcul_fluxs_mod
  
    INCLUDE "indicesol.h"
    INCLUDE "iniprint.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: itime
    INTEGER, INTENT(IN)                  :: jour
    INTEGER, INTENT(IN)                  :: knon
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex
    REAL, INTENT(IN)                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: p1lay
    REAL, DIMENSION(klon), INTENT(IN)    :: cdragh, cdragm
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)    :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefH, AcoefQ, BcoefH, BcoefQ
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)    :: ps
    REAL, DIMENSION(klon), INTENT(IN)    :: u1, v1
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf_in

! In/Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT) :: radsol
    REAL, DIMENSION(klon), INTENT(INOUT) :: snow
    REAL, DIMENSION(klon), INTENT(INOUT) :: agesno
    
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)   :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)   :: evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)   :: flux_u1, flux_v1
    REAL, DIMENSION(klon), INTENT(OUT)   :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)   :: dflux_s, dflux_l      
    REAL, DIMENSION(klon), INTENT(OUT)   :: lmt_bils

! Local variables
!****************************************************************************************
    INTEGER               :: i
    REAL, DIMENSION(klon) :: cal, beta, dif_grnd
    REAL, DIMENSION(klon) :: lmt_bils_oce, lmt_foce, diff_sst
    REAL, DIMENSION(klon) :: u0, v0
    REAL, DIMENSION(klon) :: u1_lay, v1_lay
    REAL                  :: calc_bils_oce, deltat
    REAL, PARAMETER       :: cyang=50.0 * 4.228e+06 ! capacite calorifique volumetrique de l'eau J/(m2 K)

!****************************************************************************************
! 1) Flux calculation
!
!****************************************************************************************
    cal(:)      = 0.
    beta(:)     = 1.
    dif_grnd(:) = 0.
    agesno(:)   = 0.
    
! Suppose zero surface speed
    u0(:)=0.0
    v0(:)=0.0
    u1_lay(:) = u1(:) - u0(:)
    v1_lay(:) = v1(:) - v0(:)

    CALL calcul_fluxs(knon, is_oce, dtime, &
         tsurf_in, p1lay, cal, beta, cdragh, ps, &
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
! 2) Get global variables lmt_bils and lmt_foce from file limit_slab.nc
!
!****************************************************************************************
    CALL limit_slab(itime, dtime, jour, lmt_bils, lmt_foce, diff_sst)  ! global pour un processus

    lmt_bils_oce(:) = 0.
    WHERE (lmt_foce > 0.) 
       lmt_bils_oce = lmt_bils / lmt_foce ! global 
    END WHERE

!****************************************************************************************
! 3) Recalculate new temperature
!
!****************************************************************************************
    DO i = 1, knon
       calc_bils_oce = radsol(i) + fluxsens(i) + fluxlat(i)
       deltat        = (calc_bils_oce - lmt_bils_oce(knindex(i)))*dtime/cyang +diff_sst(knindex(i))
       tsurf_new(i)  = tsurf_in(i) + deltat
    END DO

  END SUBROUTINE ocean_slab_noice
!
!****************************************************************************************
!
END MODULE ocean_slab_mod
