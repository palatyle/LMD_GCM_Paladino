!
! $Header$
!
MODULE change_srf_frac_mod

  IMPLICIT NONE

CONTAINS
! 
! Change Surface Fractions
! Author J Ghattas 2008

  SUBROUTINE change_srf_frac(itime, dtime, jour, &
       pctsrf, alb1, alb2, tsurf, u10m, v10m, pbl_tke)
!
! This subroutine is called from physiq.F at each timestep. 
! 1- For each type of ocean (force, slab, couple) receive new fractions only if
!    it's time to modify (is_modified=true). Otherwise nothing is done (is_modified=false).   
! If received new fraction :
! 2- Tests and ajustements are done on the fractions 
! 3- Initialize variables where a new fraction(new or melted ice) has appered, 
!

    USE dimphy 
    USE surface_data, ONLY : type_ocean
    USE limit_read_mod
    USE pbl_surface_mod, ONLY : pbl_surface_newfrac
    USE cpl_mod, ONLY : cpl_receive_frac
    USE ocean_slab_mod, ONLY : ocean_slab_frac

    INCLUDE "iniprint.h"
    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: itime   ! current time step
    INTEGER, INTENT(IN)                     :: jour    ! day of the year
    REAL,    INTENT(IN)                     :: dtime   ! length of time step (s)
  
! In-Output arguments
!****************************************************************************************
   
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: pctsrf ! sub-surface fraction
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: alb1   ! albedo first interval in SW spektrum
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: alb2   ! albedo second interval in SW spektrum
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: tsurf
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: u10m
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: v10m
    REAL, DIMENSION(klon,klev+1,nbsrf), INTENT(INOUT) :: pbl_tke

! Loccal variables
!****************************************************************************************
    INTEGER                        :: i, nsrf
    LOGICAL                        :: is_modified   ! true if pctsrf is modified at this time step
    LOGICAL                        :: test_sum=.FALSE.
    LOGICAL, DIMENSION(klon,nbsrf) :: new_surf
    REAL, DIMENSION(klon,nbsrf)    :: pctsrf_old    ! fraction from previous time-step
    REAL                           :: tmpsum

    pctsrf_old(:,:) = pctsrf(:,:)
!****************************************************************************************
! 1) 
! For each type of ocean (force, slab, couple) receive new fractions only if it's time  
! to modify (is_modified=true). Otherwise nothing is done (is_modified=false).   
!****************************************************************************************
    SELECT CASE (type_ocean)
    CASE ('force')
       ! Read fraction from limit.nc
       CALL limit_read_frac(itime, dtime, jour, pctsrf, is_modified)
    CASE ('slab')
       ! Get fraction from slab module
       CALL ocean_slab_frac(itime, dtime, jour, pctsrf, is_modified)
    CASE ('couple')
       ! Get fraction from the coupler
       CALL cpl_receive_frac(itime, dtime, pctsrf, is_modified)
    END SELECT


!****************************************************************************************
! 2) 
! Tests and ajustements on the new fractions :
! - Put to zero fractions that are too small
! - Test total fraction sum is one for each grid point
!
!****************************************************************************************
    IF (is_modified) THEN
  
! Test and exit if a fraction is negative
       IF (MINVAL(pctsrf(:,:)) < 0.) THEN
          WRITE(lunout,*)'Warning! One or several fractions are negative, itime=',itime
          WRITE(lunout,*)'at point = ',MINLOC(pctsrf(:,:))
          WRITE(lunout,*)'value = ',MINVAL(pctsrf(:,:)) 
          CALL abort_gcm('change_srf_frac','Negative fraction',1)
       END IF

! Optional test on the incoming fraction 
       IF (test_sum) THEN
          DO i= 1, klon
             tmpsum = SUM(pctsrf(i,:))
             IF (ABS(1. - tmpsum) > 0.05) CALL abort_gcm('change_srf_frac','Total fraction not equal 1.',1)
          END DO
       END IF

! Test for too small fractions of the sum land+landice and ocean+sea-ice
       WHERE ((pctsrf(:,is_ter) + pctsrf(:,is_lic)) < 2*EPSFRA)
          pctsrf(:,is_ter) = 0.
          pctsrf(:,is_lic) = 0.
       END WHERE

       WHERE ((pctsrf(:,is_oce) + pctsrf(:,is_sic)) < 2*EPSFRA)
          pctsrf(:,is_oce) = 0.
          pctsrf(:,is_sic) = 0.
       END WHERE

! Normalize to force total fraction to be equal one
       DO i= 1, klon
          tmpsum = SUM(pctsrf(i,:))
          DO nsrf = 1, nbsrf
             pctsrf(i,nsrf) = pctsrf(i,nsrf) / tmpsum
          END DO
       END DO

! Test for too small fractions at each sub-surface
       WHERE (pctsrf(:,is_ter) < EPSFRA)
          pctsrf(:,is_lic) = pctsrf(:,is_ter) + pctsrf(:,is_lic)
          pctsrf(:,is_ter) = 0.
       END WHERE

       WHERE (pctsrf(:,is_lic) < EPSFRA)
          pctsrf(:,is_ter) = pctsrf(:,is_ter) + pctsrf(:,is_lic)
          pctsrf(:,is_lic) = 0.
       END WHERE

       WHERE (pctsrf(:,is_oce) < EPSFRA)
          pctsrf(:,is_sic) = pctsrf(:,is_oce) + pctsrf(:,is_sic)
          pctsrf(:,is_oce) = 0.
       END WHERE

       WHERE (pctsrf(:,is_sic) < EPSFRA)
          pctsrf(:,is_oce) = pctsrf(:,is_oce) + pctsrf(:,is_sic)
          pctsrf(:,is_sic) = 0.
       END WHERE

!****************************************************************************************
! 3)
! Initialize variables where a new fraction has appered, 
! i.e. where new sea ice has been formed
! or where ice free ocean has appread in a grid cell
! 
!****************************************************************************************
       CALL pbl_surface_newfrac(itime, pctsrf, pctsrf_old, tsurf, alb1, alb2, u10m, v10m, pbl_tke)

    ELSE
       ! No modifcation should be done
       pctsrf(:,:) = pctsrf_old(:,:)

    END IF ! is_modified

  END SUBROUTINE change_srf_frac


END MODULE change_srf_frac_mod
