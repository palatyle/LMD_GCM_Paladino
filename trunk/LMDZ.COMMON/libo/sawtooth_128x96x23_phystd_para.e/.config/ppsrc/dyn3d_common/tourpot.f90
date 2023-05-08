










SUBROUTINE tourpot ( vcov, ucov, massebxy, vorpot )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van.
!-------------------------------------------------------------------------------
! Purpose: Compute potential vorticity.
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: vcov    (ip1jm,  llm)
  REAL, INTENT(IN)  :: ucov    (ip1jmp1,llm)
  REAL, INTENT(IN)  :: massebxy(ip1jm,  llm)
  REAL, INTENT(OUT) :: vorpot  (ip1jm,  llm)
!===============================================================================
! Method used:
!   vorpot = ( Filtre( d(vcov)/dx - d(ucov)/dy ) + fext ) /psbarxy
!===============================================================================
! Local variables:
  INTEGER :: l, ij
  REAL    :: rot(ip1jm,llm)
!===============================================================================

!--- Wind vorticity ; correction: rot(iip1,j,l) = rot(1,j,l)
  DO l=1,llm
    DO ij=1,ip1jm-1
      rot(ij,l)=vcov(ij+1,l)-vcov(ij,l)+ucov(ij+iip1,l)-ucov(ij,l)
    END DO
    DO ij=iip1,ip1jm,iip1; rot(ij,l)=rot(ij-iim,l); END DO
  END DO

!--- Filter
  CALL  filtreg(rot,jjm,llm,2,1,.FALSE.,1)

!--- Potential vorticity ; correction: rot(iip1,j,l) = rot(1,j,l)
  DO l=1,llm
    DO ij=1,ip1jm-1
      vorpot(ij,l)=(rot(ij,l)+fext(ij))/massebxy(ij,l)
    END DO
    DO ij=iip1,ip1jm,iip1; vorpot(ij,l)=vorpot(ij-iim,l); END DO
  END DO

END SUBROUTINE tourpot
