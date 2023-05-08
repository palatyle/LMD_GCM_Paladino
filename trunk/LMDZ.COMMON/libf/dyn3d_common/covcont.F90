SUBROUTINE covcont (klevel,ucov, vcov, ucont, vcont )
!
!-------------------------------------------------------------------------------
! Author: P. Le Van
!-------------------------------------------------------------------------------
! Purpose: Compute contravariant components from covariant components.
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  INTEGER, INTENT(IN)  :: klevel                    !--- VERTICAL LEVELS NUMBER
  REAL,    INTENT(IN)  :: ucov ( ip1jmp1,klevel )   !--- U COVARIANT WIND
  REAL,    INTENT(IN)  :: vcov ( ip1jm  ,klevel )   !--- V COVARIANT WIND
  REAL,    INTENT(OUT) :: ucont( ip1jmp1,klevel )   !--- U CONTRAVAR WIND
  REAL,    INTENT(OUT) :: vcont( ip1jm  ,klevel )   !--- V CONTRAVAR WIND
!===============================================================================
!   Local variables:
  INTEGER :: l
!===============================================================================
  DO l=1,klevel
    ucont(iip2:ip1jm,l)=ucov(iip2:ip1jm,l) * unscu2(iip2:ip1jm)
    vcont(   1:ip1jm,l)=vcov(   1:ip1jm,l) * unscv2(   1:ip1jm)
  END DO

END SUBROUTINE covcont

