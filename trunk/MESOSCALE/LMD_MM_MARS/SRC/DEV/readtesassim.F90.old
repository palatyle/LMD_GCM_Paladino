subroutine readtesassim(ngrid,nlayer,zday,pplev,tauref)

! Reading of the dust assimilation file

implicit none

include "dimensions.h"
include "dimphys.h"
include "comgeomfi.h"
include "netcdf.inc"
include "datafile.h"

integer, intent(in) :: ngrid,nlayer
real, intent(in) :: zday ! date in martian days
real, dimension(ngrid,nlayer+1), intent(in) :: pplev
real, dimension(ngrid), intent(out) :: tauref

! Local variables

real :: realday
integer nid,nvarid,ierr
integer ltloop,lsloop,iloop,jloop,varloop,ig
real, dimension(2) :: taubuf
real tau1(4),tau
real alt(4)
integer latp(4),lonp(4)
integer yinf,ysup,xinf,xsup,tinf,tsup
real latinf,latsup,loninf,lonsup
real latintmp,lonintmp,latdeg,londeg
real colat,dlat,dlon,colattmp
logical, save :: firstcall=.true.
logical :: timeflag
real radeg,pi
integer :: timedim,londim,latdim
real, dimension(:), allocatable, save :: lat,lon,time
real, dimension(:,:,:), allocatable, save :: tautes
integer, save :: timelen,lonlen,latlen
INTEGER, external:: lnblnk

pi=acos(-1.)
radeg=180/pi
realday=mod(zday,669.)

if (firstcall) then
   firstcall=.false.

   ! Note: datafile() is defined in "datafile.h"
   ierr=NF_OPEN (datafile(1:lnblnk(datafile))//"/dust_tes.nc",NF_NOWRITE,nid)
   IF (ierr.NE.NF_NOERR) THEN
     write(*,*)'Problem opening dust.nc (in phymars/readtesassim.F90)'
     write(*,*)'It should be in :',datafile(1:lnblnk(datafile)),'/'
     write(*,*)'1) You can change this directory address in '
     write(*,*)'   file phymars/datafile.h'
     write(*,*)'2) If necessary, dust.nc (and other datafiles)'
     write(*,*)'   can be obtained online on:'
     write(*,*)'   http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
     CALL ABORT
   ENDIF

   ierr=NF_INQ_DIMID(nid,"Time",timedim)
   ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
   ierr=NF_INQ_DIMID(nid,"latitude",latdim)
   ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
   ierr=NF_INQ_DIMID(nid,"longitude",londim)
   ierr=NF_INQ_DIMLEN(nid,londim,lonlen)

   allocate(tautes(lonlen,latlen,timelen))
   allocate(lat(latlen), lon(lonlen), time(timelen))

   ierr = NF_INQ_VARID (nid, "dustop",nvarid)
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tautes)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, tautes)
#endif
   IF (ierr .NE. NF_NOERR) THEN
      PRINT*, "Error: Readtesassim <dustop> not found"
      stop
   ENDIF
   ierr = NF_INQ_VARID (nid, "Time",nvarid)
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif
   IF (ierr .NE. NF_NOERR) THEN
      PRINT*, "Error: Readtesassim <Time> not found"
      stop
   ENDIF

   ierr = NF_INQ_VARID (nid, "latitude",nvarid)
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, lat)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, lat)
#endif
   IF (ierr .NE. NF_NOERR) THEN
      PRINT*, "Error: Readtesassim <latitude> not found"
      stop
   ENDIF

   ierr = NF_INQ_VARID (nid, "longitude",nvarid)
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, lon)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, lon)
#endif
   IF (ierr .NE. NF_NOERR) THEN
      PRINT*, "Error: Readtesassim <longitude> not found"
      stop
   ENDIF

endif ! of if (firstcall)

do ig=1,ngrid

! Find the four nearest points, arranged as follows:
!                               1 2
!                               3 4

   latdeg=lati(ig)*radeg  ! latitude, in degrees
   londeg=long(ig)*radeg  ! longitude, in degrees east
   colat=90-latdeg        ! colatitude, in degrees

 ! Find encompassing latitudes
   if (colat<(90-lat(1))) then ! between north pole and lat(1)
      ysup=1
      yinf=1
   else if (colat>=90-(lat(latlen))) then ! between south pole and lat(laten)
      ysup=latlen
      yinf=latlen
   else ! general case
      do iloop=2,latlen
         if(colat<(90-lat(iloop))) then
           ysup=iloop-1
           yinf=iloop
           exit
         endif
      enddo
   endif

   latinf=lat(yinf)
   latsup=lat(ysup)


 ! Find encompassing longitudes
   ! Note: in input file, lon(1)=-180.
   if (londeg>lon(lonlen)) then
      xsup=1
      xinf=lonlen
      loninf=lon(xsup)
      lonsup=180.0 ! since lon(1)=-180.0
   else
      do iloop=2,lonlen
         if(londeg<lon(iloop)) then
              xsup=iloop
              xinf=iloop-1
              exit
         endif
      enddo
      loninf=lon(xinf)
      lonsup=lon(xsup)
   endif

   if ((xsup.gt.lonlen).OR.(yinf.gt.latlen).OR.(xinf.lt.1)&
     .OR.(ysup.lt.1)) then
      write (*,*) "Readtesassim: SYSTEM ERROR on x or y in ts_gcm"
      write (*,*) "xinf: ",xinf
      write (*,*) "xsup: ",xsup
      write (*,*) "yinf: ",yinf
      write (*,*) "ysup: ",ysup
      stop
   endif

!   loninf=lon(xinf)
!   lonsup=lon(xsup)
!   latinf=lat(yinf)
!   latsup=lat(ysup)

! The four neighbouring points are arranged as follows:
!                               1 2
!                               3 4

   latp(1)=ysup
   latp(2)=ysup
   latp(3)=yinf
   latp(4)=yinf

   lonp(1)=xinf
   lonp(2)=xsup
   lonp(3)=xinf
   lonp(4)=xsup

! Linear interpolation on time, for all four neighbouring points

   if ((realday<time(1)).or.(realday>time(timelen))) then
      tinf=timelen
      tsup=1
      timeflag=.true.
   else 
      do iloop=2,timelen
         if (realday<time(iloop)) then
            tinf=iloop-1
            tsup=iloop
            exit
         endif
      enddo
   endif

! Bilinear interpolation on the four nearest points

   if ((colat<(90-lat(1))).OR.(colat>(90-lat(latlen))).OR.(latsup==latinf)) then
      dlat=0
   else
      dlat=((90-latinf)-colat)/((90-latinf)-(90-latsup))
   endif

   if (lonsup==loninf) then
      dlon=0
   else
      dlon=(londeg-loninf)/(lonsup-loninf)
   endif

   do iloop=1,4
      taubuf(1)=tautes(lonp(iloop),latp(iloop),tinf)
      taubuf(2)=tautes(lonp(iloop),latp(iloop),tsup)
      if (timeflag) then
         if (realday>time(timelen)) then
            tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*(realday-time(tinf))/(time(tsup)+(669-time(tinf))) 
         else
            tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*realday/(time(tsup)+(669-time(tinf)))
         endif
      else
         tau1(iloop)=taubuf(1)+(taubuf(2)-taubuf(1))*(realday-time(tinf))/(time(tsup)-time(tinf))
      endif
      if (tau1(iloop)<0) then
          write (*,*) "Readtesassim: SYSTEM ERROR on tau"
          write (*,*) "utime ",realday
          write (*,*) "time(tinf) ",time(tinf)
          write (*,*) "time(tsup) ",time(tsup)
          write (*,*) "tau1 ",taubuf(1)
          write (*,*) "tau2 ",taubuf(2)
          write (*,*) "tau ",tau1(iloop)
          stop
      endif
   enddo
   timeflag=.false.

   if ((dlat>1).OR.(dlon>1) .OR. (dlat<0) .OR. (dlon<0)) then
      write (*,*) "Readtesassim: SYSTEM ERROR on dlat or dlon in ts_gcm"
      write (*,*) "dlat: ",dlat
      write (*,*) "lat: ",latdeg
      write (*,*) "dlon: ",dlon
      write (*,*) "lon: ",londeg
      stop
   endif

   tau= dlat*(dlon*(tau1(2)+tau1(3)-tau1(1)-tau1(4))+tau1(1)-tau1(3)) +dlon*(tau1(4)-tau1(3))+tau1(3)

!   tauref(ig)=tau*700/pplev(ig,1)*0.8
    tauref(ig)=tau*0.825
enddo

end
