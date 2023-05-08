!
! $Id: $
!
MODULE callphysiq_mod

IMPLICIT NONE

 !! ----------------------------------------------------------
 !! ---- CALL TO NO PHYSICS
 !! ----------------------------------------------------------

CONTAINS

SUBROUTINE call_physiq(planet_type, klon,llm,nqtot,                       &
                       debut_split,lafin_split)

  USE variables_mod

  IMPLICIT NONE

  character(len=10),INTENT(IN) :: planet_type ! planet type ('earth','mars',...)
  INTEGER,INTENT(IN) :: klon ! (local) number of atmospheric columns
  INTEGER,INTENT(IN) :: llm  ! number of atmospheric layers
  INTEGER,INTENT(IN) :: nqtot ! number of tracers
  LOGICAL,INTENT(IN) :: debut_split ! .true. if very first call to physics
  LOGICAL,INTENT(IN) :: lafin_split ! .true. if last call to physics

END SUBROUTINE call_physiq

END MODULE callphysiq_mod
