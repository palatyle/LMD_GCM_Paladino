!
! $Id: $
!
MODULE callphysiq_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE call_physiq(planet_type,klon,llm,nqtot,                  &
                       debut_split,lafin_split)
 
  USE variables_mod
  USE physiq_mod, ONLY: physiq
  use phys_state_var_mod , only : qsurf
  use comm_wrf, only : allocate_comm_wrf
  use tracer_h, only: noms
 
  IMPLICIT NONE
  character(len=10),INTENT(IN) :: planet_type ! planet type ('earth','mars',...)
  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
  LOGICAL,INTENT(IN) :: debut_split ! .true. if very first call to physics
  LOGICAL,INTENT(IN) :: lafin_split ! .true. if last call to physics

  ! Local variables
!  CHARACTER(len=11) :: modname="call_physiq"

! Sanity check on physics package type
  IF (debut_split) THEN
    IF (planet_type.ne."generic") THEN
     PRINT*,"wrong planet_type for this physics package"
     STOP
    ENDIF
  ENDIF

  qsurf(:,:)=zqfi_omp(:,1,:)
  print*,'zdt_split',zdt_split

  call allocate_comm_wrf(klon,llm)

! Call physics package with required inputs/outputs
  CALL physiq(klon,           & ! ngrid
              llm,            & ! nlayer
              nqtot,          & ! nq
              noms,          & ! nametrac
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
