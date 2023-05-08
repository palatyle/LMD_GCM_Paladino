      subroutine inistats(ierr)

      implicit none

#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comvert.h"
#include "comconst.h"
#include "statto.h"
#include "netcdf.inc"

      integer,intent(out) :: ierr
      integer :: nid
      integer :: l,nsteppd
      real, dimension(llm) ::  sig_s
      integer :: idim_lat,idim_lon,idim_llm,idim_llmp1,idim_time
      real, dimension(istime) :: lt
      integer :: nvarid

      write (*,*) 
      write (*,*) '                        || STATS ||'
      write (*,*) 
      write (*,*) 'daysec',daysec
      write (*,*) 'dtphys',dtphys
      nsteppd=nint(daysec/dtphys)
      write (*,*) 'nsteppd=',nsteppd
      if (abs(float(nsteppd)-daysec/dtphys).gt.1.e-8*daysec)
     ,   stop'Dans Instat:  1jour .ne. n pas physiques'

      if(mod(nsteppd,istime).ne.0)
     ,   stop'Dans Instat:  1jour .ne. n*istime pas physiques'

      istats=nsteppd/istime
      write (*,*) 'istats=',istats
      write (*,*) 'Storing ',istime,'times per day'
      write (*,*) 'thus every ',istats,'physical timestep '
      write (*,*) 

      do l= 1, llm
         sig_s(l)=((ap(l)+ap(l+1))/preff+bp(l)+bp(l+1))/2.
         pseudoalt(l)=-10.*log(presnivs(l)/preff)   
      enddo

      ierr = NF_CREATE("stats.nc",NF_CLOBBER,nid)
      if (ierr.ne.NF_NOERR) then
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif

      ierr = NF_DEF_DIM (nid, "latitude", jjp1, idim_lat)
      ierr = NF_DEF_DIM (nid, "longitude", iip1, idim_lon)
      ierr = NF_DEF_DIM (nid, "altitude", llm, idim_llm)
      ierr = NF_DEF_DIM (nid, "llmp1", llm+1, idim_llmp1)
      ierr = NF_DEF_DIM (nid, "Time", NF_UNLIMITED, idim_time)

      ierr = NF_ENDDEF(nid)
      call def_var(nid,"Time","Time",
     .            "days since 0000-00-0 00:00:00",1,
     .            idim_time,nvarid,ierr)
! Time is initialised later by mkstats subroutine

      call def_var(nid,"latitude","latitude","degrees_north",1,
     .            idim_lat,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu/pi*180)
#endif
      call def_var(nid,"longitude","East longitude","degrees_east",1,
     .            idim_lon,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv/pi*180)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv/pi*180)
#endif

! Niveaux verticaux, aps et bps
      ierr = NF_REDEF (nid)
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"altitude", NF_DOUBLE, 1,idim_llm,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"altitude", NF_FLOAT, 1,idim_llm,nvarid)
#endif
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,"long_name",8,"altitude")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'units',2,"km")
      ierr = NF_PUT_ATT_TEXT (nid,nvarid,'positive',2,"up")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,pseudoalt)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,pseudoalt)
#endif 
      call def_var(nid,"aps","hybrid pressure at midlayers"," ",
     .            1,idim_llm,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aps)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aps)
#endif

      call def_var(nid,"bps","hybrid sigma at midlayers"," ",
     .            1,idim_llm,nvarid,ierr)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bps)
#else 
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bps)
#endif

      ierr=NF_CLOSE(nid)

      end
