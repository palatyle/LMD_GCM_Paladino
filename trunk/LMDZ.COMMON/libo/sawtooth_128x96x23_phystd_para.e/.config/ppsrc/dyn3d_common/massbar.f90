










SUBROUTINE massbar(masse,massebx,masseby)
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
  REAL, INTENT(IN)  :: masse  (ip1jmp1,llm)
  REAL, INTENT(OUT) :: massebx(ip1jmp1,llm)
  REAL, INTENT(OUT) :: masseby(ip1jm  ,llm)
!-------------------------------------------------------------------------------
! Method used. Each scalar point is associated to 4 area coefficients:
!    * alpha1(i,j) at point ( i+1/4,j-1/4 )
!    * alpha2(i,j) at point ( i+1/4,j+1/4 )
!    * alpha3(i,j) at point ( i-1/4,j+1/4 )
!    * alpha4(i,j) at point ( i-1/4,j-1/4 )
! where alpha1(i,j) = aire(i+1/4,j-1/4)/ aire(i,j)
!
!   alpha4 .         . alpha1    . alpha4
!    (i,j)             (i,j)       (i+1,j)
!
!             P .        U .          . P
!           (i,j)       (i,j)         (i+1,j)
!
!   alpha3 .         . alpha2    .alpha3 
!    (i,j)              (i,j)     (i+1,j)
!
!             V .        Z .          . V
!           (i,j)
!
!   alpha4 .         . alpha1    .alpha4
!   (i,j+1)            (i,j+1)   (i+1,j+1) 
!
!             P .        U .          . P
!          (i,j+1)                    (i+1,j+1)
!
!
!    massebx(i,j) = masse(i  ,j) * ( alpha1(i  ,j) + alpha2(i,j))   +
!                   masse(i+1,j) * ( alpha3(i+1,j) + alpha4(i+1,j) )
!     localized at point  ... U (i,j) ...
!
!    masseby(i,j) = masse(i,j  ) * ( alpha2(i,j  ) + alpha3(i,j  )  +
!                   masse(i,j+1) * ( alpha1(i,j+1) + alpha4(i,j+1)  
!     localized at point  ... V (i,j) ...
!===============================================================================
! Local variables:
  INTEGER :: ij, l
!===============================================================================
  DO l=1,llm
    DO ij=1,ip1jmp1-1
      massebx(ij,l)=masse(ij,l)*alpha1p2(ij)+masse(ij+1   ,l)*alpha3p4(ij+1)
    END DO
    DO ij=iip1,ip1jmp1,iip1; massebx(ij,l)=massebx(ij-iim,l); END DO
    DO ij=1,ip1jm
      masseby(ij,l)=masse(ij,l)*alpha2p3(ij)+masse(ij+iip1,l)*alpha1p4(ij+iip1)
    END DO
  END DO

END SUBROUTINE massbar

