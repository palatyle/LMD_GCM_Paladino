SUBROUTINE massbarxy(masse,massebxy)
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute air mass mean along X and Y in each cell.
! See iniconst for more details.
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: masse   (ip1jmp1,llm)
  REAL, INTENT(OUT) :: massebxy(ip1jm  ,llm)
!===============================================================================
! Local variables:
  INTEGER :: ij, l
!===============================================================================
  DO l=1,llm
    DO ij=1,ip1jm-1
      massebxy(ij,l)=masse(ij     ,l)*alpha2(ij     ) + &
                     masse(ij+1   ,l)*alpha3(ij+1   ) + &
                     masse(ij+iip1,l)*alpha1(ij+iip1) + &
                     masse(ij+iip2,l)*alpha4(ij+iip2)
    END DO
    DO ij=iip1,ip1jm,iip1; massebxy(ij,l)=massebxy(ij-iim,l); END DO
  END DO

END SUBROUTINE massbarxy
