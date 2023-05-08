SUBROUTINE enercin ( vcov, ucov, vcont, ucont, ecin )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van.
!-------------------------------------------------------------------------------
! Purpose: Compute kinetic energy at sigma levels.
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: vcov    (ip1jm,  llm)
  REAL, INTENT(IN)  :: ucov    (ip1jmp1,llm)
  REAL, INTENT(IN)  :: vcont   (ip1jm,  llm)
  REAL, INTENT(IN)  :: ucont   (ip1jmp1,llm)
  REAL, INTENT(OUT) :: ecin    (ip1jmp1,llm)
!===============================================================================
! Notes:
!                 . V
!                i,j-1
!
!      alpha4 .       . alpha1
!
!
!        U .      . P     . U
!       i-1,j    i,j      i,j
!
!      alpha3 .       . alpha2
!
!
!                 . V
!                i,j
!
! Kinetic energy at scalar point P(i,j) (excluding poles) is:
!       Ecin = 0.5 * U(i-1,j)**2 *( alpha3 + alpha4 )  +
!              0.5 * U(i  ,j)**2 *( alpha1 + alpha2 )  +
!              0.5 * V(i,j-1)**2 *( alpha1 + alpha4 )  +
!              0.5 * V(i,  j)**2 *( alpha2 + alpha3 )
!===============================================================================
! Local variables:
  INTEGER :: l, ij, i
  REAL    :: ecinni(iip1), ecinsi(iip1), ecinpn, ecinps
!===============================================================================
  DO l=1,llm
    DO ij = iip2, ip1jm -1
      ecin(ij+1,l)=0.5*(ucov(ij    ,l)*ucont(ij    ,l)*alpha3p4(ij +1)          &
                      + ucov(ij+1  ,l)*ucont(ij+1  ,l)*alpha1p2(ij +1)          &
                      + vcov(ij-iim,l)*vcont(ij-iim,l)*alpha1p4(ij +1)          &
                      + vcov(ij+1  ,l)*vcont(ij+1  ,l)*alpha2p3(ij +1) )
    END DO
    !--- Correction: ecin(1,j,l)= ecin(iip1,j,l)
    DO ij=iip2,ip1jm,iip1; ecin(ij,l) = ecin(ij+iim,l); END DO

    !--- North pole
    DO i=1,iim
      ecinni(i) = vcov(i,l)*vcont(i,l)*aire(i)
    END DO
    ecinpn = 0.5*SUM(ecinni(1:iim))/apoln
    DO ij=1,iip1; ecin(ij,l)=ecinpn; END DO

    !--- South pole
    DO i=1,iim
      ecinsi(i) = vcov(i+ip1jmi1,l)*vcont(i+ip1jmi1,l)*aire(i+ip1jm)
    END DO
    ecinps = 0.5*SUM(ecinsi(1:iim))/apols
    DO ij=1,iip1; ecin(ij+ip1jm,l)=ecinps; END DO
  END DO

END SUBROUTINE enercin

