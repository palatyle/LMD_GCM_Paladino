SUBROUTINE flumass (massebx,masseby, vcont, ucont, pbaru, pbarv )
!
!-------------------------------------------------------------------------------
! Authors: P. Le Van , Fr. Hourdin.
!-------------------------------------------------------------------------------
! Purpose: Compute mass flux at s levels.
  IMPLICIT NONE
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
!===============================================================================
! Arguments:
  REAL, INTENT(IN)  :: massebx(ip1jmp1,llm)
  REAL, INTENT(IN)  :: masseby(ip1jm  ,llm)
  REAL, INTENT(IN)  :: vcont  (ip1jm  ,llm)
  REAL, INTENT(IN)  :: ucont  (ip1jmp1,llm)
  REAL, INTENT(OUT) :: pbaru  (ip1jmp1,llm)
  REAL, INTENT(OUT) :: pbarv  (ip1jm  ,llm)
!===============================================================================
! Method used:   A 2 equations system is solved.
!   * 1st one describes divergence computation at pole point nr. i (i=1 to im):
!     (0.5*(pbaru(i)-pbaru(i-1))-pbarv(i))/aire(i) = - SUM(pbarv(n))/aire pole
!   * 2nd one specifies that mean mass flux at pole is equal to 0:
!     SUM(pbaru(n)*local_area(n))=0
! This way, we determine additive constant common to pbary elements representing
!   pbaru(0,j,l) in divergence computation equation for point i=1. (i=1 to im)
!===============================================================================
! Local variables:
  REAL    :: sairen, saireun, ctn, ctn0, apbarun(iip1)
  REAL    :: saires, saireus, cts, cts0, apbarus(iip1)
  INTEGER :: l, i
!===============================================================================
  DO l=1,llm
    pbaru(iip2:ip1jm,l)=massebx(iip2:ip1jm,l)*ucont(iip2:ip1jm,l)
    pbarv(   1:ip1jm,l)=masseby(   1:ip1jm,l)*vcont(   1:ip1jm,l)
  END DO

  !--- NORTH POLE
  sairen =SUM(aire (1:iim))
  saireun=SUM(aireu(1:iim))
  DO l = 1,llm
    ctn=SUM(pbarv(1:iim,l))/sairen
    pbaru(1,l)= pbarv(1,l)-ctn*aire(1)
    DO i=2,iim
      pbaru(i,l)=pbaru(i-1,l)+pbarv(i,l)-ctn*aire(i)
    END DO
    DO i=1,iim
      apbarun(i)=aireu(i)*pbaru(i,l)
    END DO
    ctn0 = -SUM(apbarun(1:iim))/saireun
    DO i = 1,iim
      pbaru(i,l)=2.*(pbaru(i,l)+ctn0)
    END DO
    pbaru(iip1,l)=pbaru(1,l)
  END DO

  !--- SOUTH POLE
  saires =SUM(aire (ip1jm+1:ip1jmp1-1))
  saireus=SUM(aireu(ip1jm+1:ip1jmp1-1))
  DO l = 1,llm
    cts=SUM(pbarv(ip1jmi1+1:ip1jm-1,l))/saires
    pbaru(1+ip1jm,l)=-pbarv(1+ip1jmi1,l)+cts*aire(1+ip1jm)
    DO i=2,iim
      pbaru(i+ip1jm,l)=pbaru(i-1+ip1jm,l)-pbarv(i+ip1jmi1,l)+cts*aire(i+ip1jm)
    END DO
    DO i=1,iim
      apbarus(i)=aireu(i+ip1jm)*pbaru(i+ip1jm,l)
    END DO
    cts0 = -SUM(apbarus(1:iim))/saireus
    DO i = 1,iim
      pbaru(i+ip1jm,l)=2.*(pbaru(i+ip1jm,l)+cts0)
    END DO
    pbaru(ip1jmp1,l)=pbaru(1+ip1jm,l)
  END DO

END SUBROUTINE flumass
