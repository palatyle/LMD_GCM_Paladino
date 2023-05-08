      subroutine nlthermeq(nlon, nlayer, pplev, pplay)
c
c  Compute the number of layers nlaylte (stored in common yomlw.h)
c  over which local thermodynamic equilibrium
c  radiation scheme should be run to be sure of covering at least to a
c  height greater than (pressure lower than) p=pminte, set in nlteparams.h.
c  The maximum layer needed is found for the worst possible case.
c  Stephen Lewis 6/2000
c  Modified Y. Wanherdrick/ F. Forget 09/2000
      use dimphy
      implicit none
c#include "dimradmars.h"
#include "nlteparams.h"
c#include "yomlw.h"
#include "clesphys.h"

c
c     Input:
      integer nlon, nlayer
      real pplev(nlon, nlayer+1)
      real pplay(nlon, nlayer)
c
c     Local:
      integer igpmax, ismax
      logical firstcall
      data firstcall /.true./
      save firstcall, igpmax


      INTEGER i,ix
      real sxmax

ccc      
      if(firstcall) then
c     Find the location of maximum surface pressure.
c     Location won't vary much so only do it at the start;
c     with no topography location would vary, but this is only
c     needed for an estimate so any point would do in that case.
      ismax=1  
      sxmax=pplev(1,1)
      ix=1
        do i=1,nlon-1
         if(pplev(i,ix).gt.sxmax) then
           sxmax=pplev(i,ix)
           ismax=i+1
         endif
       enddo


         igpmax = ismax            ! longitude/ latitude where pression is maximum
         write(*, 10) ptrans
         write(*, 20) zw
         write(*, 30) pminte
         firstcall = .false.
      endif

      IF(callnlte .or. callnirco2) THEN
c       Find first layer above pminte at this location
        do nlaylte = nlayer, 1, -1
      if (pplay(igpmax, nlaylte).gt.pminte)  go to 100
        enddo
      ELSE
        nlaylte=nlayer        
      END IF

c
 100    return
c

   10 format(' nlthermeq: transition to NLTE centred at ',f6.2,'Pa')
   20 format('               half-width (scale heights) ',f6.2)
   30 format('          suggested LTE coverage at least ',f6.2,'Pa')
   40 format(' nlthermeq: purely NLTE contribution over (nlayer) ',f6.4)

      end
