SUBROUTINE pression( ngrid, ap, bp, ps, p )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr.Hourdin 
!-------------------------------------------------------------------------------
! Purpose: Compute pressure p(l) at different levels from l = 1 (ground level)
!          to l = llm +1. Those levels correspond to the llm layers interfaces,
!          with p(ij,llm+1) = 0. and  p(ij,1) = ps(ij)  .   
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
!===============================================================================
! Arguments:
  INTEGER, INTENT(IN)  :: ngrid                 !--- NUMBER OF GRID POINTS
  REAL,    INTENT(IN)  :: ap(llmp1), bp(llmp1)  !--- HYBRID COEFFICIENTS
  REAL,    INTENT(IN)  :: ps(ngrid)             !--- SURFACE PRESSURE
  REAL,    INTENT(OUT) :: p(ngrid,llmp1)        !--- 3D PRESSURE FIELD 
!===============================================================================
! Local variables:
  INTEGER :: l
!===============================================================================
  DO l=1,llmp1;  p(:,l) = ap(l) + bp(l) * ps(:);  END DO

END SUBROUTINE pression


