SUBROUTINE convmas (pbaru, pbarv, convm)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux convergence at p levels.
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: pbaru(ip1jmp1,llm)
  REAL, INTENT(IN)  :: pbarv(ip1jm  ,llm)
  REAL, INTENT(OUT) :: convm(ip1jmp1,llm)
!===============================================================================
! Method used:   Computation from top to bottom.
!   Mass convergence at level llm is equal to zero and is not stored in convm.
!===============================================================================
! Local variables:
  INTEGER :: l
!===============================================================================

!--- Computation of - (d(pbaru)/dx + d(pbarv)/dy )
  CALL convflu( pbaru, pbarv, llm, convm )

!--- Filter
  CALL filtreg( convm, jjp1, llm, 2, 2, .TRUE., 1 )

!--- Mass convergence is integrated from top to bottom
  DO l=llmm1,1,-1
    convm(:,l) = convm(:,l) + convm(:,l+1)
  END DO

END SUBROUTINE convmas
