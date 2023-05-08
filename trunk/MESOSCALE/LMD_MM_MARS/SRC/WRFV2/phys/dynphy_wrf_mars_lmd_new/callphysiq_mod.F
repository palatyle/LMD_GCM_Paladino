!
! $Id: $
!
MODULE callphysiq_mod

IMPLICIT NONE

 !! ----------------------------------------------------------
 !! ---- CALL TO MARS PHYSICS
 !! ----------------------------------------------------------

CONTAINS

SUBROUTINE call_physiq(planet_type, klon,llm,nqtot,                       &
                       debut_split,lafin_split)

  USE variables_mod
  USE physiq_mod, ONLY: physiq
  IMPLICIT NONE

  character(len=10),INTENT(IN) :: planet_type ! planet type ('earth','mars',...)

  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
  LOGICAL,INTENT(IN) :: debut_split ! .true. if very first call to physics
  LOGICAL,INTENT(IN) :: lafin_split ! .true. if last call to physics

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
   IF (planet_type.ne."mars") THEN
     PRINT*,"wrong planet_type for this physics package"
     STOP
   ENDIF
 ENDIF

! Set dummy variables for Mars to zero (additional and prob useless security)
  zpk_omp(:,:)=0.
  zphis_omp(:)=0.
  presnivs_omp(:)=0.
  zrfi_omp(:,:)=0.
 
! Call physics package with required inputs/outputs
  CALL physiq(klon,           & ! ngrid
              llm,            & ! nlayer
              nqtot,          & ! nq
              debut_split,    & ! firstcall
              lafin_split,    & ! lastcall
              jD_cur,         & ! pday
              jH_cur_split,   & ! ptime
              zdt_split,      & ! ptimestep
              zplev_omp,      & ! pplev
              zplay_omp,      & ! pplay
              zphi_omp,       & ! pphi
              zufi_omp,       & ! pu
              zvfi_omp,       & ! pv
              ztfi_omp,       & ! pt
              zqfi_omp,       & ! pq
              flxwfi_omp,     & ! flxw
              zdufi_omp,      & ! pdu
              zdvfi_omp,      & ! pdv
              zdtfi_omp,      & ! pdt
              zdqfi_omp,      & ! pdq
              zdpsrf_omp)       ! pdpsrf


END SUBROUTINE call_physiq

END MODULE callphysiq_mod
