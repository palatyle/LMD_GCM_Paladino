
! $Id: $

MODULE infotrac_phy

! Infotrac for physics; contains the same information as infotrac for
! the dynamics
  IMPLICIT NONE

! iflag_trac: ==1 if running with tracers
  INTEGER,SAVE :: iflag_trac
!$OMP THREADPRIVATE(iflag_trac)

! nqtot : total number of tracers
  INTEGER,SAVE :: nqtot
!$OMP THREADPRIVATE(nqtot)

! tracer names
  CHARACTER(len=20),ALLOCATABLE,DIMENSION(:),SAVE :: tname
  CHARACTER(len=23),ALLOCATABLE,DIMENSION(:),SAVE :: ttext ! tracer long name for diagnostics
!$OMP THREADPRIVATE(tname,ttext)

CONTAINS

  SUBROUTINE init_infotrac_phy(iflag_trac_,nqtot_,tname_,ttext_)
  ! Initialize module variables
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: iflag_trac_ ! ==1 if running with tracers
  INTEGER,INTENT(IN) :: nqtot_ ! total number of tracers
  CHARACTER(LEN=*),INTENT(IN) :: tname_(nqtot_)
  CHARACTER(LEN=*),INTENT(IN) :: ttext_(nqtot_)
  
  INTEGER :: iq
  
  iflag_trac=iflag_trac_
  nqtot=nqtot_
  
  ALLOCATE(tname(nqtot))
  ALLOCATE(ttext(nqtot))
  DO iq=1,nqtot
    tname(iq)=tname_(iq)
    ttext(iq)=ttext_(iq)
  ENDDO
  
  END SUBROUTINE init_infotrac_phy

END MODULE infotrac_phy
