










subroutine iniwritesoil(nid,ngrid,inertia,area,nbplon,nbplat)

! initialization routine for 'writediagoil'. Here we create/define
! dimensions (longitude, latitude, depth and time) and other fixed
! (time-independent) parameters.

use comsoil_h, only: mlayer, nsoilmx
USE comcstfi_mod, only: pi
USE regular_lonlat_mod, ONLY: lon_reg, lat_reg
use mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat

implicit none

include"netcdf.inc"

! Arguments:
integer,intent(in) :: ngrid
integer,intent(in) :: nid ! NetCDF output file ID
real,intent(in) :: inertia(nbplon,nbplat,nsoilmx)
real,intent(in) :: area(nbplon,nbp_lat) ! mesh area (m2)
integer,intent(in) :: nbplon,nbplat ! sizes of area

! Local variables:

! NetCDF stuff:
integer :: ierr ! NetCDF routines return code
integer :: idim_rlatu ! ID of the 'latitude' dimension
integer :: idim_rlonv ! ID of the 'longitude' dimension
integer :: idim_depth ! ID of the 'depth' dimension
integer :: idim_time  ! ID of the 'time' dimension
integer :: varid ! to store NetCDF ID of a variable
integer,dimension(3) :: dimids ! to store IDs of dimensions of a variable
character(len=60) :: text ! to store some text
real,dimension(nbp_lon+1,nbp_lat,nsoilmx) :: data3 ! to store 3D data
integer :: i,j,l,ig0
real,allocatable :: lon_reg_ext(:) ! extended longitudes


if (nbp_lon*nbp_lat==1) then
  ! 1D model
  allocate(lon_reg_ext(1))
else
  ! 3D model
  allocate(lon_reg_ext(nbp_lon+1))
endif

! 1. Define the dimensions
! Switch to NetCDF define mode
ierr=NF_REDEF(nid)

! Define the dimensions
if (nbp_lon*nbp_lat==1) then
  ierr=NF_DEF_DIM(nid,"longitude",1,idim_rlonv)
else
  ierr=NF_DEF_DIM(nid,"longitude",nbp_lon+1,idim_rlonv)
endif
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define longitude dimension"
endif
ierr=NF_DEF_DIM(nid,"latitude",nbp_lat,idim_rlatu)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define latitude dimension"
endif
ierr=NF_DEF_DIM(nid,"depth",nsoilmx,idim_depth)
! nsoilmx known from comsoil_h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth dimension"
endif
ierr=NF_DEF_DIM(nid,"time",NF_UNLIMITED,idim_time)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define time dimension"
endif

! Switch out of NetCDF define mode
ierr=NF_ENDDEF(nid)

! 2. Define (as variables) and write dimensions, as well as their attributes
! 2.1. Longitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
ierr=NF_DEF_VAR(nid,"longitude",NF_DOUBLE,1,idim_rlonv,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define longitude variable"
endif

! Longitude attributes
text="East longitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_east"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

lon_reg_ext(1:nbp_lon)=lon_reg(1:nbp_lon)
if (nbp_lon*nbp_lat/=1) then ! in 3D only:
  ! add extra redundant point (180 degrees, since lon_reg starts at -180
  lon_reg_ext(nbp_lon+1)=-lon_reg_ext(1)
endif

! Write longitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
ierr=NF_PUT_VAR_DOUBLE(nid,varid,lon_reg_ext*(180./pi))
! Note: rlonv is known from comgeom.h and pi from comcstfi.h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write longitude variable"
endif

! 2.2. Latitude
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
ierr=NF_DEF_VAR(nid,"latitude",NF_DOUBLE,1,idim_rlatu,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define latitude variable"
endif

! Latitude attributes
text="North latitude"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="degrees_north"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write latitude to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
ierr=NF_PUT_VAR_DOUBLE(nid,varid,lat_reg*(180./pi))
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write longitude variable"
endif

! 2.3. Depth
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
ierr=NF_DEF_VAR(nid,"depth",NF_DOUBLE,1,idim_depth,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth variable"
endif

! Depth attributes
text="Soil mid-layer depth"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="m"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)
text="down"
ierr=NF_PUT_ATT_TEXT(nid,varid,"positive",len_trim(text),text)

! Write depth to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
ierr=NF_PUT_VAR_DOUBLE(nid,varid,mlayer)
! Note mlayer(0:nsoilmx-1) known from comsoil_h
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write depth variable"
endif

! 2.4. Time
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
ierr=NF_DEF_VAR(nid,"time",NF_DOUBLE,1,idim_time,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define depth variable"
endif

! time attributes
text="Time"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="days since 0000-01-01 00:00:00"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Note no need to write time variable here; it is done in writediagsoil.

! 3. Other variables to be included

! 3.1 mesh area surrounding each horizontal point
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
dimids(1)=idim_rlonv ! ID of the 'longitude' dimension
dimids(2)=idim_rlatu ! ID of the 'latitude' dimension
ierr=NF_DEF_VAR(nid,"area",NF_DOUBLE,2,dimids,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define area variable"
endif

! Area attributes
text="Mesh area"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="m2"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write area to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
ierr=NF_PUT_VAR_DOUBLE(nid,varid,area)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write area variable"
endif

! 3.2 Thermal inertia
ierr=NF_REDEF(nid) ! switch to NetCDF define mode

! Define the variable
dimids(1)=idim_rlonv ! ID of the 'longitude' dimension
dimids(2)=idim_rlatu ! ID of the 'latitude' dimension
dimids(3)=idim_depth ! ID of the 'depth' dimension
ierr=NF_DEF_VAR(nid,"th_inertia",NF_DOUBLE,3,dimids,varid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not define th_inertia variable"
endif

! Attributes
text="Thermal inertia"
ierr=NF_PUT_ATT_TEXT(nid,varid,"long_name",len_trim(text),text)
text="J.s-1/2.m-2.K-1"
ierr=NF_PUT_ATT_TEXT(nid,varid,"units",len_trim(text),text)

! Write data2 to file
ierr=NF_ENDDEF(nid) ! switch out of NetCDF define mode
! Write
ierr=NF_PUT_VAR_DOUBLE(nid,varid,inertia)
if (ierr.ne.NF_NOERR) then
  write(*,*)"iniwritesoil: Error, could not write th_inertia variable"
endif

end subroutine iniwritesoil
