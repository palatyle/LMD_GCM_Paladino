      subroutine thermosphere(ngrid,nlayer,nq,
     &     pplev,pplay,dist_sol,
     $     mu0,ptimestep,ptime,zday,tsurf,zzlev,zzlay,
     &     pt,pq,pu,pv,pdt,pdq,
     $     zdteuv,zdtconduc,zdumolvis,zdvmolvis,zdqmoldiff)

      use conc_mod, only: rnew, cpnew
      USE comcstfi_h, only: r, cpp
      implicit none

#include "callkeys.h"

      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of advected tracers
      REAL,INTENT(in) :: pplay(ngrid,nlayer)
      REAL,INTENT(in) :: pplev(ngrid,nlayer+1)
      REAL,INTENT(in) :: zzlay(ngrid,nlayer)
      REAL,INTENT(in) :: zzlev(ngrid,nlayer+1)
      REAL,INTENT(in) :: pt(ngrid,nlayer)
      REAL,INTENT(in) :: zday
      REAL,INTENT(in) :: dist_sol
      REAL,INTENT(in) :: mu0(ngrid)
      REAL,INTENT(in) :: pq(ngrid,nlayer,nq)
      REAL,INTENT(in) :: ptimestep
      REAL,INTENT(in) :: ptime
      REAL,INTENT(in) :: tsurf(ngrid)
      REAL,INTENT(in) :: pu(ngrid,nlayer),pv(ngrid,nlayer)
      REAL,INTENT(in) :: pdt(ngrid,nlayer),pdq(ngrid,nlayer,nq)

      REAL,INTENT(out) :: zdteuv(ngrid,nlayer)
      REAL,INTENT(out) :: zdtconduc(ngrid,nlayer)
      REAL,INTENT(out) :: zdumolvis(ngrid,nlayer)
      REAL,INTENT(out) :: zdvmolvis(ngrid,nlayer)
      REAL,INTENT(out) :: zdqmoldiff(ngrid,nlayer,nq)

      INTEGER :: l,ig
      logical,save :: firstcall=.true.

      if (firstcall) then
        if (.not. tracer) then
          do l=1,nlayer
            do ig=1,ngrid
              rnew(ig,l)=r
              cpnew(ig,l)=cpp
            enddo
          enddo
        endif
        firstcall= .false.
      endif

      ! initialize tendencies to zero in all cases
      ! (tendencies are added later on, even if parametrization is not called)
      zdteuv(1:ngrid,1:nlayer)=0
      zdtconduc(1:ngrid,1:nlayer)=0
      zdumolvis(1:ngrid,1:nlayer)=0
      zdvmolvis(1:ngrid,1:nlayer)=0
      zdqmoldiff(1:ngrid,1:nlayer,1:nq)=0
      
      if (calleuv) then
        call euvheat(ngrid,nlayer,nq,pt,pdt,pplev,pplay,zzlay,
     $               mu0,ptimestep,ptime,zday,pq,pdq,zdteuv)
      endif

      if (callconduct) THEN
        call conduction(ngrid,nlayer,ptimestep,pplay,pplev,pt,zdteuv,
     $                   tsurf,zzlev,zzlay,zdtconduc)
      endif

      if (callmolvis) THEN
        call molvis(ngrid,nlayer,ptimestep,pplay,pplev,pt,
     &                zdteuv,zdtconduc,pu,
     $                   tsurf,zzlev,zzlay,zdumolvis)
        call molvis(ngrid,nlayer,ptimestep,pplay,pplev,pt,
     &                zdteuv,zdtconduc,pv,
     $                   tsurf,zzlev,zzlay,zdvmolvis)
      endif

      if (callmoldiff) THEN
        call moldiff_red(ngrid,nlayer,nq,
     &                   pplay,pplev,pt,pdt,pq,pdq,ptimestep,
     &                   zzlay,zdteuv,zdtconduc,zdqmoldiff)
      endif

      end


