      subroutine nlthermeq(ngrid, nlayer, pplev, pplay)
c
c  Compute the number of layers nlaylte (stored in common yomlw.h)
c  over which local thermodynamic equilibrium
c  radiation scheme should be run to be sure of covering at least to a
c  height greater than (pressure lower than) p=pminte, set in nlteparams.h.
c  The maximum layer needed is found for the worst possible case.
c  Stephen Lewis 6/2000
c  Modified Y. Wanherdrick/ F. Forget 09/2000

      implicit none
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "nlteparams.h"
#include "yomlw.h"
#include "callkeys.h"

c
c     Input:
      integer ngrid, nlayer
      real pplev(ngrid, nlayer+1)
      real pplay(ngrid, nlayer)
c
c     Local:
      integer igpmax, ismax
      logical firstcall
      data firstcall /.true./
      save firstcall, igpmax
c
      if(firstcall) then
c     Find the location of maximum surface pressure.
c     Location won't vary much so only do it at the start;
c     with no topography location would vary, but this is only
c     needed for an estimate so any point would do in that case.
         igpmax = ismax(ngrid, pplev, 1)
         write(*, 10) ptrans
         write(*, 20) zw
         write(*, 30) pminte
         firstcall = .false.
      endif
c
      IF(callnlte) then
c       Find first layer above pminte at this location
        do nlaylte = nlayer, 1, -1
           if (pplay(igpmax, nlaylte).gt.pminte)  go to 100
        enddo
      ELSE
        nlaylte=nlayer
      END IF
  100 write(*,*) 'LTE rad. calculations up to layer ',  nlaylte
c
      return
c
   10 format(' nlthermeq: transition to NLTE centred at ',f6.2,'Pa')
   20 format('               half-width (scale heights) ',f6.2)
   30 format('          suggested LTE coverage at least ',f6.2,'Pa')
      end
