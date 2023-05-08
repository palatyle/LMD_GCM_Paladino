
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
  CHARACTER(len=30),ALLOCATABLE,DIMENSION(:),SAVE :: tname
  CHARACTER(len=33),ALLOCATABLE,DIMENSION(:),SAVE :: ttext ! tracer long name for diagnostics
!$OMP THREADPRIVATE(tname,ttext)

CONTAINS

  SUBROUTINE init_infotrac_phy(iflag_trac_,nqtot_,tname_,ttext_)
  ! Initialize module variables
  
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: iflag_trac_ ! ==1 if running with tracers
  INTEGER,INTENT(IN) :: nqtot_ ! total number of tracers
  CHARACTER(LEN=*),INTENT(IN) :: tname_(nqtot_)
  CHARACTER(LEN=*),INTENT(IN) :: ttext_(nqtot_)
  
  INTEGER :: iq
  CHARACTER(LEN=30) :: modname="init_infotrac_phy"
  CHARACTER(LEN=50) :: abort_message
  
  iflag_trac=iflag_trac_
  nqtot=nqtot_
  
  ALLOCATE(tname(nqtot))
  ALLOCATE(ttext(nqtot))
  DO iq=1,nqtot
    if (len_trim(tname_(iq)).gt.len(tname(iq))) then
      abort_message="init_infotrac_phy error; tname_() too long!"
      call abort_physic(modname,abort_message,1)
    endif
    tname(iq)=tname_(iq)
    ttext(iq)=ttext_(iq)
  ENDDO
  
  END SUBROUTINE init_infotrac_phy

END MODULE infotrac_phy
