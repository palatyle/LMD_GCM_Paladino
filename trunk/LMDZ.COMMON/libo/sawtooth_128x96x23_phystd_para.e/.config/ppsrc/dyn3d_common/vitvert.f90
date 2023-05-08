










!
! $Header$
!
SUBROUTINE vitvert (convm, w)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute vertical speed at sigma levels.
  USE comvert_mod, ONLY: bp
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: convm(ip1jmp1,llm)
  REAL, INTENT(OUT) :: w    (ip1jmp1,llm)
!===============================================================================
! Notes: Vertical speed is oriented from bottom to top.
!   * At ground - level sigma(1):     w(i,j,1) = 0.
!   * At top    - level sigma(llm+1): w(i,j,l) = 0. (not stored in w)
!===============================================================================
! Local variables:
  INTEGER :: l
!===============================================================================
  DO l=1,llmm1; w(:,l+1)=convm(:,l+1)-bp(l+1)*convm(:,1); END DO
  w(:,1)=0.

END SUBROUTINE vitvert

