










subroutine mkstats(ierr)

 
!
!  This program writes a stats.nc file from sums and sums of squares
!  to means and standard deviations and also writes netcdf style
!  file so that the data can be viewed easily.  The data file is
!  overwritten in place.  
!  SRL  21 May 1996
!  Yann W. july 2003

use statto_mod, only: istime,count
use mod_phys_lmdz_para, only : is_master
use mod_grid_phy_lmdz, only : nbp_lon, nbp_lat, nbp_lev, klon_glo

implicit none

include "netcdf.inc"

integer :: ierr,nid,nbvar,i,ndims,lt,nvarid
integer, dimension(4) :: id,varid,start,size
integer, dimension(5) :: dimids
character (len=50) :: name,nameout,units,title
real,allocatable :: sum3d(:,:,:),square3d(:,:,:),mean3d(:,:,:),sd3d(:,:,:)
real,allocatable :: sum2d(:,:),square2d(:,:),mean2d(:,:),sd2d(:,:)
real, dimension(istime) :: time
real, dimension(nbp_lat) :: lat
real,allocatable :: lon(:)
real, dimension(nbp_lev) :: alt
logical :: lcopy=.true.
!integer :: latid,lonid,altid,timeid
integer :: meanid,sdid
!integer, dimension(4) :: dimout

! Incrementation of count for the last step, which is not done in wstats
count(istime)=count(istime)+1

if (is_master) then
! only the master needs do this
if (klon_glo>1) then
  allocate(lon(nbp_lon+1))
  allocate(sum3d(nbp_lon+1,nbp_lat,nbp_lev))
  allocate(square3d(nbp_lon+1,nbp_lat,nbp_lev))
  allocate(mean3d(nbp_lon+1,nbp_lat,nbp_lev))
  allocate(sd3d(nbp_lon+1,nbp_lat,nbp_lev))
  allocate(sum2d(nbp_lon+1,nbp_lat))
  allocate(square2d(nbp_lon+1,nbp_lat))
  allocate(mean2d(nbp_lon+1,nbp_lat))
  allocate(sd2d(nbp_lon+1,nbp_lat))
else ! 1D model case
  allocate(lon(1))
  allocate(sum3d(1,1,nbp_lev))
  allocate(square3d(1,1,nbp_lev))
  allocate(mean3d(1,1,nbp_lev))
  allocate(sd3d(1,1,nbp_lev))
  allocate(sum2d(1,1))
  allocate(square2d(1,1))
  allocate(mean2d(1,1))
  allocate(sd2d(1,1))
endif

ierr = NF_OPEN("stats.nc",NF_WRITE,nid)

! We catch the id of dimensions of the stats file

ierr= NF_INQ_DIMID(nid,"latitude",id(1))
ierr= NF_INQ_DIMID(nid,"longitude",id(2))
ierr= NF_INQ_DIMID(nid,"altitude",id(3))
ierr= NF_INQ_DIMID(nid,"Time",id(4))

ierr= NF_INQ_VARID(nid,"latitude",varid(1))
ierr= NF_INQ_VARID(nid,"longitude",varid(2))
ierr= NF_INQ_VARID(nid,"altitude",varid(3))
ierr= NF_INQ_VARID(nid,"Time",varid(4))

! Time initialisation

do i=1,istime
   time(i)=i*24./istime
   ierr= NF_PUT_VARA_DOUBLE(nid,varid(4),i,1,time(i))
enddo

! We catche the values of the variables

         ierr = NF_GET_VAR_DOUBLE(nid,varid(1),lat)
         ierr = NF_GET_VAR_DOUBLE(nid,varid(2),lon)
         ierr = NF_GET_VAR_DOUBLE(nid,varid(3),alt)

! We catch the number of variables in the stats file
ierr = NF_INQ_NVARS(nid,nbvar)

! to catche the "real" number of variables (without the "additionnal variables")
nbvar=(nbvar-4)/2 

do i=1,nbvar
   varid=(i-1)*2+5

   ! What's the variable's name?
   ierr=NF_INQ_VARNAME(nid,varid,name)
   write(*,*) "OK variable ",name
   ! Its units?
   units=" "
   ierr=NF_GET_ATT_TEXT(nid,varid,"units",units)
   ! Its title?
   title=" "
   ierr=NF_GET_ATT_TEXT(nid,varid,"title",title)
   ! Its number of dimensions?   
   ierr=NF_INQ_VARNDIMS(nid,varid,ndims)
   ! Its values?

   if(ndims==4) then ! lat, lon, alt & time

!      dimout(1)=lonid
!      dimout(2)=latid
!      dimout(3)=altid
!      dimout(4)=timeid

      if (klon_glo>1) then ! general case
        size=(/nbp_lon+1,nbp_lat,nbp_lev,1/)
      else ! 1D model
        size=(/1,1,nbp_lev,1/)
      endif
      do lt=1,istime
         start=(/1,1,1,lt/)
         ! Extraction of the "source" variables
         ierr = NF_GET_VARA_DOUBLE(nid,varid,start,size,sum3d)
         ierr = NF_GET_VARA_DOUBLE(nid,varid+1,start,size,square3d)
         ! Calculation of these variables
         mean3d=sum3d/count(lt)
         sd3d=sqrt(max(0.,square3d/count(lt)-mean3d**2))
         ! Writing of the variables
         ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,size,mean3d)
         ierr = NF_PUT_VARA_DOUBLE(nid,varid+1,start,size,sd3d)
      enddo

    else if (ndims.eq.3) then

!      dimout(1)=lonid
!      dimout(2)=latid
!      dimout(3)=timeid

      if (klon_glo>1) then ! general case
        size=(/nbp_lon+1,nbp_lat,1,0/)
      else
        size=(/1,1,1,0/)
      endif
      do lt=1,istime
         start=(/1,1,lt,0/)
         ! Extraction of the "source" variables
         ierr = NF_GET_VARA_DOUBLE(nid,varid,start,size,sum2d)
         ierr = NF_GET_VARA_DOUBLE(nid,varid+1,start,size,square2d)
         ! Calculation of these variables
         mean2d=sum2d/count(lt)
         sd2d=sqrt(max(0.,square2d/count(lt)-mean2d**2))
         ! Writing of the variables
         ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,size,mean2d)
         ierr = NF_PUT_VARA_DOUBLE(nid,varid+1,start,size,sd2d)
      enddo

    endif 
enddo

ierr= NF_CLOSE(nid)

endif ! of if (is_master)

end
