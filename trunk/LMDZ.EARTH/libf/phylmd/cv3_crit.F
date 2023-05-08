        SUBROUTINE CV3_CRIT (nloc,ncum,nd,icb,inb,p,ph,pzero
     $             ,v,threshold,kcrit,pcrit)
***************************************************************
*                                                             *
* CV3_CRIT   Find pressure level where vertical profile of    *
*            variable 'v' intersects 'threshold'              *
*                                                             *
* written by   : FROHWIRTH Julie, 13/08/2003, 21.55.12        *
* modified by :                                               *
***************************************************************
*

#include "cv3param.h"

c input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real pzero(nloc)
      real v(nloc,nd),threshold

c output:
      integer kcrit(nloc)
      real pcrit(nloc)

c local variables
      integer i,j,k,il
      logical ok(nloc)

      do il = 1,ncum
        ok(il) = .true.
        pcrit(il) = -1.
        kcrit(il) = 0
      enddo
c
      DO i = 1,nl
        DO il = 1,ncum
        IF (i .GT. icb(il) .AND. i .LE. inb(il)) THEN
        IF (P(il,i) .LE. Pzero(il) .AND. ok(il)) THEN
          IF ( (v(il,i)-threshold)*(v(il,i-1)-threshold) .LT. 0.) THEN
            pcrit(il) =
     $       ((threshold-v(il,i))*P(il,i-1)-
     $        (threshold-v(il,i-1))*P(il,i))
     $           /(v(il,i-1)-v(il,i))
           IF (pcrit(il) .gt. Pzero(il)) THEN
            pcrit(il) = -1.
           ELSE
            ok(il) = .false.
            kcrit(il) = i
            IF (pcrit(il) .LT. PH(il,i)) kcrit(il) = kcrit(il)+1
           ENDIF
          ENDIF  ! end IF (v(i) ...
        ENDIF    ! end IF (P(i) ...
        ENDIF    ! end IF (icb+1 le i le inb)
        ENDDO
      ENDDO
125   CONTINUE


      RETURN
      END
