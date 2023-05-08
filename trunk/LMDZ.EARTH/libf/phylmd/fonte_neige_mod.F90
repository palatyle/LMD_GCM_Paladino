!
! $Header$
!
MODULE fonte_neige_mod
!
! This module will treat the process of snow, melting, accumulating, calving, in 
! case of simplified soil model.
!
!****************************************************************************************
  USE dimphy, ONLY : klon

  IMPLICIT NONE
  SAVE

! run_off_ter and run_off_lic are the runoff at the compressed grid knon for 
! land and land-ice respectively
! Note: run_off_lic is used in mod_landice and therfore not private
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE    :: run_off_ter
  !$OMP THREADPRIVATE(run_off_ter)
  REAL, ALLOCATABLE, DIMENSION(:)             :: run_off_lic
  !$OMP THREADPRIVATE(run_off_lic)

! run_off_lic_0 is the runoff at land-ice a time-step earlier, on the global 1D array grid
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE    :: run_off_lic_0
  !$OMP THREADPRIVATE(run_off_lic_0)
  
  REAL, PRIVATE                               :: tau_calv  
  !$OMP THREADPRIVATE(tau_calv)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE  :: ffonte_global
  !$OMP THREADPRIVATE(ffonte_global)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE  :: fqfonte_global
  !$OMP THREADPRIVATE(fqfonte_global)
  REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE  :: fqcalving_global
  !$OMP THREADPRIVATE(fqcalving_global)

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE fonte_neige_init(restart_runoff)

! This subroutine allocates and initialize variables in the module. 
! The variable run_off_lic_0 is initialized to the field read from
! restart file. The other variables are initialized to zero.
!
    INCLUDE "indicesol.h"
!****************************************************************************************
! Input argument
    REAL, DIMENSION(klon), INTENT(IN) :: restart_runoff 

! Local variables
    INTEGER                           :: error
    CHARACTER (len = 80)              :: abort_message 
    CHARACTER (len = 20)              :: modname = 'fonte_neige_init'


!****************************************************************************************
! Allocate run-off at landice and initilize with field read from restart
!
!****************************************************************************************

    ALLOCATE(run_off_lic_0(klon), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation run_off_lic'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    run_off_lic_0(:) = restart_runoff(:) 

!****************************************************************************************
! Allocate other variables and initilize to zero
!
!****************************************************************************************
    ALLOCATE(run_off_ter(klon), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation run_off_ter'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    run_off_ter(:) = 0.
    
    ALLOCATE(run_off_lic(klon), stat = error)
    IF (error /= 0) THEN
       abort_message='Pb allocation run_off_lic'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    run_off_lic(:) = 0.
    
    ALLOCATE(ffonte_global(klon,nbsrf))
    IF (error /= 0) THEN
       abort_message='Pb allocation ffonte_global'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    ffonte_global(:,:) = 0.0

    ALLOCATE(fqfonte_global(klon,nbsrf))
    IF (error /= 0) THEN
       abort_message='Pb allocation fqfonte_global'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    fqfonte_global(:,:) = 0.0

    ALLOCATE(fqcalving_global(klon,nbsrf))
    IF (error /= 0) THEN
       abort_message='Pb allocation fqcalving_global'
       CALL abort_gcm(modname,abort_message,1)
    ENDIF
    fqcalving_global(:,:) = 0.0

!****************************************************************************************
! Read tau_calv
!
!****************************************************************************************
    CALL conf_interface(tau_calv)


  END SUBROUTINE fonte_neige_init
!
!****************************************************************************************
!
  SUBROUTINE fonte_neige( knon, nisurf, knindex, dtime, &
       tsurf, precip_rain, precip_snow, &
       snow, qsol, tsurf_new, evap)
        
! Routine de traitement de la fonte de la neige dans le cas du traitement
! de sol simplifie!
! LF 03/2001
! input:
!   knon         nombre de points a traiter
!   nisurf       surface a traiter
!   knindex      index des mailles valables pour surface a traiter
!   dtime        
!   tsurf        temperature de surface
!   precip_rain  precipitations liquides
!   precip_snow  precipitations solides
!
! input/output:
!   snow         champs hauteur de neige
!   qsol         hauteur d'eau contenu dans le sol
!   tsurf_new    temperature au sol
!   evap
!
  INCLUDE "indicesol.h"
  INCLUDE "dimensions.h"
  INCLUDE "YOETHF.h"
  INCLUDE "YOMCST.h"
  INCLUDE "FCTTRE.h"
  INCLUDE "clesphys.h"

! Input variables
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: knon
    INTEGER, INTENT(IN)                  :: nisurf
    INTEGER, DIMENSION(klon), INTENT(IN) :: knindex
    REAL   , INTENT(IN)                  :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_snow
    
! Input/Output variables
!****************************************************************************************

    REAL, DIMENSION(klon), INTENT(INOUT) :: snow
    REAL, DIMENSION(klon), INTENT(INOUT) :: qsol
    REAL, DIMENSION(klon), INTENT(INOUT) :: tsurf_new
    REAL, DIMENSION(klon), INTENT(INOUT) :: evap

! Local variables
!****************************************************************************************

    INTEGER               :: i, j
    REAL                  :: fq_fonte
    REAL                  :: coeff_rel
    REAL, PARAMETER       :: snow_max=3000.
    REAL, PARAMETER       :: max_eau_sol = 150.0
!! PB temporaire en attendant mieux pour le modele de neige
! REAL, parameter :: chasno = RLMLT/(2.3867E+06*0.15)
    REAL, PARAMETER       :: chasno = 3.334E+05/(2.3867E+06*0.15)
!IM cf JLD/ GKtest
    REAL, PARAMETER       :: chaice = 3.334E+05/(2.3867E+06*0.15)
! fin GKtest
    REAL, DIMENSION(klon) :: ffonte
    REAL, DIMENSION(klon) :: fqcalving, fqfonte
    REAL, DIMENSION(klon) :: d_ts
    REAL, DIMENSION(klon) :: bil_eau_s, snow_evap

    LOGICAL               :: neige_fond

!****************************************************************************************
! Start calculation
! - Initialization
!
!****************************************************************************************
    coeff_rel = dtime/(tau_calv * rday)
    
    bil_eau_s(:) = 0.

!****************************************************************************************
! - Increment snow due to precipitation and evaporation
! - Calculate the water balance due to precipitation and evaporation (bil_eau_s)
!
!****************************************************************************************
    WHERE (precip_snow > 0.) 
       snow = snow + (precip_snow * dtime)
    END WHERE

    snow_evap = 0.
    WHERE (evap > 0. ) 
       snow_evap = MIN (snow / dtime, evap) 
       snow = snow - snow_evap * dtime
       snow = MAX(0.0, snow)
    END WHERE
    
    bil_eau_s(:) = (precip_rain(:) * dtime) - (evap(:) - snow_evap(:)) * dtime


!****************************************************************************************
! - Calculate melting snow
! - Calculate calving and decrement snow, if there are to much snow
! - Update temperature at surface
!
!****************************************************************************************

    ffonte(:) = 0.0
    fqcalving(:) = 0.0
    fqfonte(:) = 0.0
    DO i = 1, knon
       ! Y'a-t-il fonte de neige?
       neige_fond = ((snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
            .AND. tsurf_new(i) >= RTT)
       IF (neige_fond) THEN
          fq_fonte     = MIN( MAX((tsurf_new(i)-RTT )/chasno,0.0),snow(i))
          ffonte(i)    = fq_fonte * RLMLT/dtime
          fqfonte(i)   = fq_fonte/dtime
          snow(i)      = MAX(0., snow(i) - fq_fonte)
          bil_eau_s(i) = bil_eau_s(i) + fq_fonte 
          tsurf_new(i) = tsurf_new(i) - fq_fonte * chasno  

!IM cf JLD OK     
!IM cf JLD/ GKtest fonte aussi pour la glace
          IF (nisurf == is_sic .OR. nisurf == is_lic ) THEN
             fq_fonte = MAX((tsurf_new(i)-RTT )/chaice,0.0)
             ffonte(i) = ffonte(i) + fq_fonte * RLMLT/dtime
             IF ( ok_lic_melt ) THEN
                fqfonte(i) = fqfonte(i) + fq_fonte/dtime
                bil_eau_s(i) = bil_eau_s(i) + fq_fonte
             ENDIF
             tsurf_new(i) = RTT
          ENDIF
          d_ts(i) = tsurf_new(i) - tsurf(i)
       ENDIF

       ! s'il y a une hauteur trop importante de neige, elle s'coule
       fqcalving(i) = MAX(0., snow(i) - snow_max)/dtime
       snow(i)=MIN(snow(i),snow_max)
    END DO


    IF (nisurf == is_ter) THEN
       DO i = 1, knon
          qsol(i) = qsol(i) + bil_eau_s(i)
          run_off_ter(i) = run_off_ter(i) + MAX(qsol(i) - max_eau_sol, 0.0)
          qsol(i) = MIN(qsol(i), max_eau_sol) 
       END DO
    ELSE IF (nisurf == is_lic) THEN
       DO i = 1, knon
          j = knindex(i)
          run_off_lic(i)   = (coeff_rel *  fqcalving(i)) + &
               (1. - coeff_rel) * run_off_lic_0(j)
          run_off_lic_0(j) = run_off_lic(i)
          run_off_lic(i)   = run_off_lic(i) + fqfonte(i)/dtime
       END DO
    ENDIF
    
!****************************************************************************************
! Save ffonte, fqfonte and fqcalving in global arrays for each 
! sub-surface separately
!
!****************************************************************************************
    DO i = 1, knon
       ffonte_global(knindex(i),nisurf)    = ffonte(i)
       fqfonte_global(knindex(i),nisurf)   = fqfonte(i)
       fqcalving_global(knindex(i),nisurf) = fqcalving(i)
    ENDDO

  END SUBROUTINE fonte_neige
!
!****************************************************************************************
!
  SUBROUTINE fonte_neige_final(restart_runoff)
!
! This subroutine returns run_off_lic_0 for later writing to restart file.
!
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT) :: restart_runoff

!****************************************************************************************
! Set the output variables
    restart_runoff(:) = run_off_lic_0(:)

! Deallocation of all varaibles in the module
!   DEALLOCATE(run_off_lic_0, run_off_ter, run_off_lic, ffonte_global, &
!        fqfonte_global, fqcalving_global)

    IF (ALLOCATED(run_off_lic_0)) DEALLOCATE(run_off_lic_0)
    IF (ALLOCATED(run_off_ter)) DEALLOCATE(run_off_ter)
    IF (ALLOCATED(run_off_lic)) DEALLOCATE(run_off_lic)
    IF (ALLOCATED(ffonte_global)) DEALLOCATE(ffonte_global)
    IF (ALLOCATED(fqfonte_global)) DEALLOCATE(fqfonte_global)
    IF (ALLOCATED(fqcalving_global)) DEALLOCATE(fqcalving_global)

  END SUBROUTINE fonte_neige_final
!
!****************************************************************************************
!
  SUBROUTINE fonte_neige_get_vars(pctsrf, fqcalving_out, &
       fqfonte_out, ffonte_out)

! Cumulate ffonte, fqfonte and fqcalving respectively for
! all type of surfaces according to their fraction.
!
! This routine is called from physiq.F before histwrite.

    INCLUDE "indicesol.h"
!****************************************************************************************
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf

    REAL, DIMENSION(klon), INTENT(OUT)      :: fqcalving_out
    REAL, DIMENSION(klon), INTENT(OUT)      :: fqfonte_out
    REAL, DIMENSION(klon), INTENT(OUT)      :: ffonte_out

    INTEGER   :: nisurf
!****************************************************************************************

    ffonte_out(:)    = 0.0
    fqfonte_out(:)   = 0.0
    fqcalving_out(:) = 0.0

    DO nisurf = 1, nbsrf
       ffonte_out(:) = ffonte_out(:) + ffonte_global(:,nisurf)*pctsrf(:,nisurf)
       fqfonte_out(:) = fqfonte_out(:) + fqfonte_global(:,nisurf)*pctsrf(:,nisurf)
       fqcalving_out(:) = fqcalving_out(:) + fqcalving_global(:,nisurf)*pctsrf(:,nisurf)
    ENDDO

  END SUBROUTINE fonte_neige_get_vars
!
!****************************************************************************************
!
END MODULE fonte_neige_mod



