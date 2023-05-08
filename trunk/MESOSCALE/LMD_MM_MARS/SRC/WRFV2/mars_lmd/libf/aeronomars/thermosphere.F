      subroutine thermosphere(pplev,pplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)

      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "comdiurn.h"
#include "param.h"
#include "param_v3.h"
#include "chimiedata.h"
#include "conc.h"


      INTEGER l,ig

      REAL pplay(ngridmx,nlayermx)
      real pplev(ngridmx,nlayermx+1)
      REAL zzlay(ngridmx,nlayermx)
      real zzlev(ngridmx,nlayermx+1)
      REAL pt(ngridmx,nlayermx)
      real zday
      REAL dist_sol
      real mu0(ngridmx)
      real pq(ngridmx,nlayermx,nqmx)
      real ptimestep
      real ptime
      real tsurf(ngridmx)
      REAL pu(ngridmx,nlayermx),pv(ngridmx,nlayermx)
      REAL pdt(ngridmx,nlayermx),pdq(ngridmx,nlayermx,nqmx)

      REAL zdteuv(ngridmx,nlayermx)
      REAL zdtconduc(ngridmx,nlayermx)
      REAL zdumolvis(ngridmx,nlayermx)
      REAL zdvmolvis(ngridmx,nlayermx)
      real zdqmoldiff(ngridmx,nlayermx,nqmx)

      logical firstcall
      save firstcall
      data firstcall /.true./

      if (firstcall) then
        if (.not. tracer) then
          do l=1,nlayermx
            do ig=1,ngridmx
              rnew(ig,l)=r
              cpnew(ig,l)=cpp
            enddo
          enddo
        endif
        firstcall= .false.
      endif

      if (calleuv) then
        call zerophys(ngridmx*nlayermx,zdteuv)
        call euvheat(pt,pdt,pplev,pplay,zzlay,dist_sol,
     $               mu0,ptimestep,ptime,zday,pq,pdq,zdteuv)
      endif

      if (callconduct) THEN
        call zerophys(ngridmx*nlayermx,zdtconduc)
        call conduction(ptimestep,pplay,pplev,pt,zdteuv,
     $                   tsurf,zzlev,zzlay,zdtconduc)
      endif

      if (callmolvis) THEN
        call zerophys(ngridmx*nlayermx,zdumolvis)
        call molvis(ptimestep,pplay,pplev,pt,zdteuv,zdtconduc,pu,
     $                   tsurf,zzlev,zzlay,zdumolvis)
        call zerophys(ngridmx*nlayermx,zdvmolvis)
        call molvis(ptimestep,pplay,pplev,pt,zdteuv,zdtconduc,pv,
     $                   tsurf,zzlev,zzlay,zdvmolvis)
      endif

      if (callmoldiff) THEN
        call zerophys(ngridmx*nlayermx*nqmx,zdqmoldiff)
        call moldiff(pplay,pplev,pt,pdt,pq,pdq,ptimestep,
     &                   zzlay,zdteuv,zdtconduc,zdqmoldiff)
      endif

      return
      end


