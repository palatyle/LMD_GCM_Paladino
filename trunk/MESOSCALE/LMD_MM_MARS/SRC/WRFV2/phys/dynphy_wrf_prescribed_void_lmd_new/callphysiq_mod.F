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
  ! data file
  INTEGER :: i,k
  REAL :: hr_sw(llm)
  REAL :: hr_lw(llm)
  REAL :: hr_dyn(llm)

  open(17,file='prescribed_sw.txt',form='formatted',status='old')
  open(18,file='prescribed_lw.txt',form='formatted',status='old')
  open(19,file='prescribed_dyn.txt',form='formatted',status='old')
  rewind(17)
  rewind(18)
  rewind(19)
  DO k=1,llm
    read(17,*) hr_sw(k)
    read(18,*) hr_lw(k)
    read(19,*) hr_dyn(k)
  ENDDO
  close(17)
  close(18)
  close(19)

  DO i=1,klon
    zdtfi_omp(i,:)=hr_sw(:)+hr_lw(:)+hr_dyn(:)
  ENDDO


END SUBROUTINE call_physiq

END MODULE callphysiq_mod
