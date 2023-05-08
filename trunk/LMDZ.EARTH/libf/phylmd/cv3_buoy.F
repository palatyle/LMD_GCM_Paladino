        SUBROUTINE CV3_BUOY (nloc,ncum,nd,icb,inb
     :                      ,pbase,plcl,p,ph,Ale,Cin
     :                      ,tv,tvp
     :                      ,buoy )
***************************************************************
*                                                             *
* CV3_BUOY                                                    *
*         Buoyancy corrections to account for ALE             *
*                                                             *
* written by   : MOREAU Cecile, 07/08/2003, 15.55.48          *
* modified by :                                               *
***************************************************************
*
      implicit none

#include "cvthermo.h"
#include "cv3param.h"

c input:
      integer ncum, nd, nloc
      integer icb(nloc), inb(nloc)
      real pbase(nloc),plcl(nloc)
      real p(nloc,nd), ph(nloc,nd+1)
      real Ale(nloc), Cin(nloc)
      real tv(nloc,nd), tvp(nloc,nd)

c output:
      real buoy(nloc,nd)

c local variables:
      integer il, k
      integer kmx(nloc)
      real bll(nloc), bmx(nloc)
      real gamma(nloc)
      real dgamma
      real buoymin
      logical ok(nloc)

      data dgamma /2.e-03/ !dgamma gamma
      data buoymin /2./

      logical fixed_bll
      SAVE fixed_bll
      data fixed_bll /.TRUE./
c$OMP THREADPRIVATE(fixed_bll)


c      print *,' Ale+cin ',ale(1)+cin(1)
c--------------------------------------------------------------
c      Recompute buoyancies
c--------------------------------------------------------------
      DO k = 1,nl
        DO il = 1,ncum
           buoy(il,k) = tvp(il,k) - tv(il,k)
        ENDDO
      ENDDO

c -------------------------------------------------------------
c -- Compute low level buoyancy ( function of Ale+Cin )
c -------------------------------------------------------------
      IF (fixed_bll) THEN
c
      do il = 1,ncum
        bll(il) = 0.5
      end DO
      else

      do il = 1,ncum
       IF (Ale(il)+Cin(il) .GT. 0.) THEN
        gamma(il) = 4.*buoy(il,icb(il))**2
     :           + 8.*dgamma*(Ale(il)+Cin(il))*tv(il,icb(il))/grav
        gamma(il) = max(gamma(il),1.e-10)
       ENDIF
      end do

      do il = 1,ncum
       IF (Ale(il)+Cin(il) .GT. 0.) THEN
        bll(il) = 4.*dgamma*(Ale(il)+Cin(il))*tv(il,icb(il))
     :         /(grav*(abs(buoy(il,icb(il))+0.5*sqrt(gamma(il)))))
       ENDIF
      end do

      do il = 1,ncum
       IF (Ale(il)+Cin(il) .GT. 0.) THEN
        bll(il) = min(bll(il),buoymin)
       ENDIF
      end DO
c
      ENDIF     !(fixed_bll)


c -------------------------------------------------------------
c --Get highest buoyancy among levels below LCL-200hPa
c -------------------------------------------------------------

      do il = 1,ncum
       bmx(il) =-1000.
       kmx(il) = icb(il)
       ok(il) = .true.
      end do

      do k = 1,nl
       do il = 1,ncum
        IF (Ale(il)+Cin(il) .GT. 0. .AND. ok(il)) THEN
        IF (k .GT. icb(il) .AND. k .LE. inb(il)) THEN
cc         print *,'k,p(il,k),plcl(il)-200. ', k,p(il,k),plcl(il)-200.
         IF (P(il,k) .GT. plcl(il)-200.) THEN
          IF (buoy(il,k) .GT. bmx(il)) THEN
           bmx(il) = buoy(il,k)
           kmx(il) = k
           IF (bmx(il) .GE. bll(il)) ok(il)=.false.
          ENDIF
         ENDIF
        ENDIF
        ENDIF
       end do
      end do

c      print *,' ==cv3_buoy== bll(1),bmx(1),icb(1),kmx(1) '
c     $       ,bll(1),bmx(1),icb(1),kmx(1)

c -------------------------------------------------------------
c --Calculate modified buoyancies
c -------------------------------------------------------------

      do il = 1,ncum
       IF (Ale(il)+Cin(il) .GT. 0.) THEN
        bll(il) = min(bll(il),bmx(il))
       ENDIF
      end do

      do k = 1,nl
       do il = 1,ncum
        IF (Ale(il)+Cin(il) .GT. 0.) THEN
         IF (k .GE. icb(il) .AND. k .LE. kmx(il)-1) THEN
           buoy(il,k) = bll(il)
         ENDIF
        ENDIF
       end do
      end do



      RETURN
      END
