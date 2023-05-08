!
! $Id: $
!
MODULE callphysiq_mod

IMPLICIT NONE

 !! ----------------------------------------------------------
 !! ---- CALL TO VENUS PHYSICS
 !! ----------------------------------------------------------

CONTAINS

SUBROUTINE call_physiq(planet_type, klon,llm,nqtot,                       &
                       debut_split,lafin_split)

  USE variables_mod
  USE physiq_mod, ONLY: physiq
  USE module_model_constants, only : p0,rcp,cp
  IMPLICIT NONE

  character(len=10),INTENT(IN) :: planet_type ! planet type ('earth','mars',...)

  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
  LOGICAL,INTENT(IN) :: debut_split ! .true. if very first call to physics
  LOGICAL,INTENT(IN) :: lafin_split ! .true. if last call to physics
  REAL*8 :: zplevmoy(llm+1) ! planet-averaged mean pressure (Pa) at interfaces
  REAL*8 :: ztmoy(llm)

!  ! Local variables
!  CHARACTER(len=11) :: modname="call_physiq"
!  LOGICAL,SAVE :: firstcall=.true.
!!$OMP THREADPRIVATE(firstcall)
!
!! Sanity check on physics package type
!  IF (firstcall) THEN
!    IF (planet_type.ne."mars") THEN
!      CALL abort_gcm(modname,"wrong planet_type for this physics package",1)
!    ENDIF
!    firstcall=.false.
!  ENDIF

! Sanity check on physics package type
 IF (debut_split) THEN
   IF (planet_type.ne."venus") THEN
     PRINT*,"wrong planet_type for this physics package"
     STOP
   ENDIF
 ENDIF

! Set dummy variables for Mars to zero (additional and prob useless security)
! NB: tname already filled with tracers' names (though not used here)
  zpk_omp(1:klon,1:llm)=cp*((zplay_omp(1:klon,1:llm)/p0)**rcp)
  !print*,'zpk_omp',zpk_omp(1,:)
  zphis_omp(1:klon)=zphi_omp(1:klon,1)
  presnivs_omp(:)=0.
  zrfi_omp(:,:)=0.
  ztmoy(:)=ztfi_omp(1,:)
  zplevmoy(:)=zplev_omp(1,:)
! Call physics package with required inputs/outputs
  CALL physiq(klon,           &
              llm,            &
              nqtot,          &
              debut_split,    &
              lafin_split,    &
              jD_cur,         &
              jH_cur_split,   &
              zdt_split,      &
              zplev_omp,      &
              zplay_omp,      &
              zpk_omp,        &
              zphi_omp,       &
              zphis_omp,      &
              presnivs_omp,   &
              zufi_omp,       &
              zvfi_omp,       &
              ztfi_omp,       &
              zqfi_omp,       &
              flxwfi_omp,     &
              zdufi_omp,      &
              zdvfi_omp,      &
              zdtfi_omp,      &
              zdqfi_omp,      &
              zdpsrf_omp)
END SUBROUTINE call_physiq

END MODULE callphysiq_mod
