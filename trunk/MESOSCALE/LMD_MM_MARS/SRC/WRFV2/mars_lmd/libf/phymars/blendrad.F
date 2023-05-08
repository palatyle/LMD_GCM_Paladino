      subroutine blendrad(ngrid, nlayer, pplay,
     &                    zdtsw, zdtlw, zdtnirco2, zdtnlte, dtrad)
c
c  Combine radiative tendencies.  LTE contributions (zdtsw and zdtlw)
c  have been calculated for the first NLAYLTE layers, zdtnirco2 and
c  zdtnlte have been calculated for all nlayer layers (but zdtnlte may
c  be zero low down).  zdtlw is phased out in favour of zdtnlte with
c  height; zdtsw is also phased out to remove possible spurious heating
c  at low pressures.  The pressure at which the transition occurs and
c  the scale over which this happens are set in the nlteparams.h file.
c  Above layer NLAYLTE the tendency is purely the sum of NLTE contributions.
c  (Note : nlaylte is calculated by "nlthermeq" and stored in common "yomlw.h")
c  Stephen Lewis 6/2000 FF 
c
      implicit none
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "nlteparams.h"
#include "yomlw.h"

c     Input:
      integer ngrid, nlayer
      real pplay(ngrid, nlayer)
      real zdtlw(ngrid, nlayer)
      real zdtsw(ngrid, nlayer)
      real zdtnirco2(ngrid, nlayer)
      real zdtnlte(ngrid, nlayer)
c
c     Output:
      real dtrad(ngrid, nlayer)
c
c     Local:
      integer l, ig
      real alpha
c
c     This is split into two loops to minimize number of calculations,
c     but for vector machines it may be faster to perform one big
c     loop from 1 to nlayer and remove the second loop.
c
c     Loop over layers for which zdtsw/lw have been calculated.
      do l = 1,nlaylte 
         do ig = 1, ngrid
c           alpha is actually 0.5*(1+tanh((z-ztrans)/zw))
c           written here in a simpler form, with z=-ln(p) and zwi=2/zw
            alpha = 1./(1.+(pplay(ig,l)/ptrans)**zwi)
            dtrad(ig,l) = (1.-alpha)*(zdtsw(ig,l)+zdtlw(ig,l))
     &                  + zdtnirco2(ig,l) + alpha*zdtnlte(ig,l)
         enddo
      enddo
c
c     Faster loop over any remaining layers.
      do l = nlaylte+1, nlayer
         do ig = 1, ngrid
            dtrad(ig,l) = zdtnirco2(ig,l) + zdtnlte(ig,l)
         enddo
      enddo
c
      return
      end
