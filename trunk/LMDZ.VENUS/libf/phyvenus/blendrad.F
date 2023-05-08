      subroutine blendrad(nlon, nlev, pplay, heat,
     &                     cool, pdtnirco2,zdtnlte, dtsw,dtlw )
c
c  Combine radiative tendencies.  LTE contributions (heat and cool)
c  have been calculated for the first NLAYLTE layers, zdtnirco2 and
c  zdtnlte have been calculated for all nlev layers (but zdtnlte may
c  be zero low down).  cool is phased out in favour of zdtnlte with
c  height; heat is also phased out to remove possible spurious heating
c  at low pressures.  The pressure at which the transition occurs and
c  the scale over which this happens are set in the nlteparams.h file.
c  Above layer NLAYLTE the tendency is purely the sum of NLTE contributions.
c  (Note : nlaylte is calculated by "nlthermeq" and stored in common "yomlw.h")
c  Stephen Lewis 6/2000 FF 

      use dimphy
      implicit none
c#include "dimradmars.h"
#include "nlteparams.h"
c#include "yomlw.h"
#include "YOMCST.h"

c     Input:
      integer nlon, nlev
      real pplay(nlon, nlev)
      real cool(nlon, nlev)
      real heat(nlon, nlev)
      real pdtnirco2(nlon, nlev)
      real zdtnlte(nlon, nlev)
c
c     Output:
c      real dtrad(nlon, nlev)
      real dtlw(nlon, nlev)
      real dtsw(nlon, nlev)
c
c     Local:
      integer l, ig
      real alpha, alpha2
      real, parameter :: p_lowup = 10.e2

c
c     This is split into two loops to minimize number of calculations,
c     but for vector machines it may be faster to perform one big
c     loop from 1 to nlev and remove the second loop.
c

c          print*, '--- NLAYTE value is: ---'
c          print*, nlaylte

c     Loop over layers for which heat/lw have been calculated.
       do l = 1,nlaylte !defini dans nlthermeq 
         do ig = 1, nlon
c           alpha is actually 0.5*(1+tanh((z-ztrans)/zw))
c           written here in a simpler form, with z=-ln(p) and zwi=2/zw
            alpha  = 1./(1.+(pplay(ig,l)/ptrans)**zwi)
            alpha2 = 1./(1.+(pplay(ig,l)/p_lowup)**zwi)
            
c     This formula is used in the Martian routines 
c            dtrad(ig,l) = (1.-alpha)*(heat(ig,l)+cool(ig,l))
c     &                  + pdtnirco2(ig,l) + alpha*zdtnlte(ig,l)

            dtlw(ig,l) = (1.-alpha)*(-cool(ig,l)) 
     &			+ alpha*zdtnlte(ig,l)
            dtsw(ig,l) = (1-alpha2)*(heat(ig,l)) 
     &			+  alpha2*pdtnirco2(ig,l)

         enddo
      enddo

c 
c     Faster loop over any remaining layers.
      do l = nlaylte+1, nlev
         do ig = 1, nlon

c           dtrad(ig,l) = pdtnirco2(ig,l) + zdtnlte(ig,l)

            dtlw(ig,l) = zdtnlte(ig,l)
            dtsw(ig,l) = pdtnirco2(ig,l)
         enddo
      enddo

      return
      end
