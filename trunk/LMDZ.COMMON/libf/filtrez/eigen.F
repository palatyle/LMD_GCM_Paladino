!
! $Header$
!
      SUBROUTINE eigen( e,d)
      IMPLICIT NONE
#include "dimensions.h"
      real :: e( iim,iim ), d( iim )
      real :: asm( iim )
      integer :: im,i,j
      im=iim
c
      DO 48 i = 1,im
	 asm( i ) = d( im-i+1 )
 48   CONTINUE
      DO 49 i = 1,iim
	 d( i ) = asm( i )
 49   CONTINUE
c
c     PRINT 70,d
 70   FORMAT(5x,'Valeurs propres',/,8(1x,8f10.4,/),/)
		print *
c
      DO 51 i = 1,im
	 DO 52 j = 1,im
            asm( j ) = e( i , im-j+1 )
 52      CONTINUE
	 DO 50 j = 1,im
	    e( i,j ) = asm( j )
 50      CONTINUE
 51   CONTINUE

      RETURN
      END
