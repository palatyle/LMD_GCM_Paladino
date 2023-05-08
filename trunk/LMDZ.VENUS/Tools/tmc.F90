program tmc

! SL 01/2010:
! This program reads 4D (lon-lat-alt-time) fields recast in log P coordinates
!
! it computes angular momentum transport from high-frequency outputs:
!
! totvang -- 2D -- Meridional transport of angular momentum, total (m3 s-2)
! totwang -- 2D --   Vertical transport of angular momentum, total (m3 s-2)
! mmcvang -- 2D -- Meridional transport of angular momentum, by MMC (m3 s-2)
! mmcwang -- 2D --   Vertical transport of angular momentum, by MMC (m3 s-2)
! trsvang -- 2D -- Meridional transport of angular momentum, transients (m3 s-2)
! trswang -- 2D --   Vertical transport of angular momentum, transients (m3 s-2)
! stnvang -- 2D -- Meridional transport of angular momentum, stationaries (m3 s-2)
! stnwang -- 2D --   Vertical transport of angular momentum, stationaries (m3 s-2)
! dmass   -- 2D -- Mass in each cell (dmassmeanbar)
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - pressure vertical coordinate
! - surface pressure
! - zonal, meridional and vertical (Pa/s) winds
! - altitude above areoid
!
! Convention: qbar  <=> zonal average    / qstar = q - qbar
!             qmean <=> temporal average / qprim = q - qmean
!
!  Therefore: ((qv)mean)bar                         (total)
!                          =  qmeanbar *  vmeanbar  (mmc)
!                      + (qstarmean * vstarmean)bar (stn)
!                         + (qprim * vprim)meanbar  (trs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (name_P.nc)
character (len=128) :: outfile ! output file name

character (len=64) :: text ! to store some text
integer infid ! NetCDF input file ID
integer outfid ! NetCDF output file ID
integer lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid
integer,dimension(2) :: datashape2d ! shape of 3D datasets
integer,dimension(3) :: datashape3d ! shape of 3D datasets
integer,dimension(4) :: datashape4d ! shape of 4D datasets

real :: miss_val ! special "missing value" to specify missing data
real,parameter :: miss_val_def=-9.99e+33 ! default value for "missing value"
real :: pi
real,dimension(:),allocatable :: lon ! longitude
integer lonlength ! # of grid points along longitude
real,dimension(:),allocatable :: lat ! latitude
real,dimension(:),allocatable :: coslat ! cos of latitude
integer latlength ! # of grid points along latitude
real,dimension(:),allocatable :: plev ! Pressure levels (Pa)
integer altlength ! # of grid point along altitude (of input datasets)
real,dimension(:),allocatable :: time ! time
integer timelength ! # of points along time
real,dimension(:,:,:),allocatable :: ps ! surface pressure
real,dimension(:,:,:,:),allocatable :: vitu ! zonal wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitw ! vertical wind (in Pa/s, then converted in m/s)
real,dimension(:,:,:,:),allocatable :: za ! above areoid levels (m)

!!! output variables
real,dimension(:,:),allocatable :: totvang ! merid transport of ang momentum
real,dimension(:,:),allocatable :: totwang ! verti transport of ang momentum
real,dimension(:,:),allocatable :: mmcvang ! merid transport of ang momentum
real,dimension(:,:),allocatable :: mmcwang ! verti transport of ang momentum
real,dimension(:,:),allocatable :: trsvang ! merid transport of ang momentum
real,dimension(:,:),allocatable :: trswang ! verti transport of ang momentum
real,dimension(:,:),allocatable :: stnvang ! merid transport of ang momentum
real,dimension(:,:),allocatable :: stnwang ! verti transport of ang momentum
real,dimension(:,:),allocatable :: dmassmeanbar 

! local variables
real :: deltalat,deltalon ! lat and lon intervals in radians
real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)
real,dimension(:,:,:),allocatable   :: dmassmean 

real,dimension(:,:,:,:),allocatable :: osam ! planetary rotation specific ang (m2/s)
real,dimension(:,:,:,:),allocatable :: rsam ! zonal wind specific ang (m2/s)
real,dimension(:,:,:,:),allocatable :: tsam ! total specific ang (m2/s)

real,dimension(:,:,:,:),allocatable :: vang ! v * specific ang (m3/s)
real,dimension(:,:,:,:),allocatable :: wang ! w * specific ang (m3/s)
real,dimension(:,:,:,:),allocatable :: vstar
real,dimension(:,:,:,:),allocatable :: wstar
real,dimension(:,:,:,:),allocatable :: angstar 
real,dimension(:,:,:,:),allocatable :: vprim
real,dimension(:,:,:,:),allocatable :: wprim
real,dimension(:,:,:,:),allocatable :: angprim 
real,dimension(:,:,:,:),allocatable :: angprimvprim
real,dimension(:,:,:,:),allocatable :: angprimwprim

! lon,lat,alt
real,dimension(:,:,:),allocatable :: vangmean
real,dimension(:,:,:),allocatable :: wangmean
real,dimension(:,:,:),allocatable :: vmean 
real,dimension(:,:,:),allocatable :: wmean 
real,dimension(:,:,:),allocatable :: angmean 
real,dimension(:,:,:),allocatable :: angprimvprimmean 
real,dimension(:,:,:),allocatable :: angprimwprimmean 
real,dimension(:,:,:),allocatable :: vstarmean
real,dimension(:,:,:),allocatable :: wstarmean
real,dimension(:,:,:),allocatable :: angstarmean 
real,dimension(:,:,:),allocatable :: angstarmeanvstarmean
real,dimension(:,:,:),allocatable :: angstarmeanwstarmean
real,dimension(:,:,:),allocatable :: v3d    ! intermediate for vbar
real,dimension(:,:,:),allocatable :: w3d    ! intermediate for wbar
real,dimension(:,:,:),allocatable :: ang3d  ! intermediate for angbar

! lat,alt,time
real,dimension(:,:,:),allocatable :: vbar 
real,dimension(:,:,:),allocatable :: wbar 
real,dimension(:,:,:),allocatable :: angbar

!lat,alt
real,dimension(:,:),allocatable :: vmeanbar
real,dimension(:,:),allocatable :: wmeanbar
real,dimension(:,:),allocatable :: angmeanbar
real,dimension(:,:),allocatable :: vbar2d     ! intermediate for vbar
real,dimension(:,:),allocatable :: wbar2d     ! intermediate for wbar
real,dimension(:,:),allocatable :: angbar2d   ! intermediate for qbar

integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
logical :: lmdflag

include "planet.h"

!===============================================================================
! 1. Input parameters
!===============================================================================

pi = 2.*asin(1.)
miss_val = miss_val_def

write(*,*) ""
write(*,*) "You are working on the atmosphere of ",planet

!===============================================================================
! 1.1 Input file
!===============================================================================

write(*,*) ""
write(*,*) "Program valid for files with pressure axis (*_P.nc)"
write(*,*) "Enter input file name:"

read(*,'(a128)') infile
write(*,*) ""

! open input file

ierr = NF_OPEN(infile,NF_NOWRITE,infid)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening file ',trim(infile)
   stop ""
endif

!===============================================================================
! 1.2 Get grids in lon,lat,alt(pressure),time
!===============================================================================

call get_iddim(infid,lat_varid,latlength,lon_varid,lonlength,&
                     alt_varid,altlength,time_varid,timelength,lmdflag )

allocate(lon(lonlength))
ierr=NF_GET_VAR_REAL(infid,lon_varid,lon)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading longitude"

allocate(lat(latlength))
ierr=NF_GET_VAR_REAL(infid,lat_varid,lat)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading lat"

allocate(coslat(latlength))
! Beware of rounding problems at poles...
coslat(:) = max(0.,cos(lat(:)*pi/180.))

! Lat, lon pressure intervals
deltalat = abs(lat(2)-lat(1))*pi/180.
deltalon = abs(lon(2)-lon(1))*pi/180.

allocate(plev(altlength))
ierr=NF_GET_VAR_REAL(infid,alt_varid,plev)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading altitude (ie pressure levels)"

allocate(time(timelength))
ierr=NF_GET_VAR_REAL(infid,time_varid,time)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading time"

!===============================================================================
! 1.3 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
outfile=infile(1:len_trim(infile)-3)//"_TMC.nc"
write(*,*) "Output file name is: "//trim(outfile)



!===============================================================================
! 2.1 Store needed fields 
!===============================================================================

!===============================================================================
! 2.1.1 Surface pressure
!===============================================================================
allocate(ps(lonlength,latlength,timelength))

text="ps"
call get_var3d(infid,lonlength,latlength,timelength,text,ps,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for psol instead... "
  text="psol"
  call get_var3d(infid,lonlength,latlength,timelength,text,ps,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) stop "Error: Failed to get psol ID"
endif
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading surface pressure"

!===============================================================================
! 2.1.3 Winds
!===============================================================================
allocate(vitu(lonlength,latlength,altlength,timelength))
allocate(vitv(lonlength,latlength,altlength,timelength))
allocate(vitw(lonlength,latlength,altlength,timelength))

! zonal wind vitu (in m/s)
text="vitu"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitu,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitu ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zonal wind"

! meridional wind vitv (in m/s)
text="vitv"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitv,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitv ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading meridional wind"

! vertical wind vitw (in Pa/s)
text="vitw"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitw,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitw ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading vertical wind"

!===============================================================================
! 2.1.4 Altitude above areoide
!===============================================================================
! Only needed if g(z) on Titan...

! allocate(za(lonlength,latlength,altlength,timelength))

! text="zareoid"
! call get_var4d(infid,lonlength,latlength,altlength,timelength,text,za,miss_val,ierr1,ierr2)
! if (ierr1.ne.NF_NOERR) stop "Error: Failed to get za ID"
! if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zareoid"

!===============================================================================
! 2.2 Computations 
!===============================================================================

!===============================================================================
! 2.2.2 Mass in cells
!===============================================================================
allocate(rayon(lonlength,latlength,altlength,timelength))
allocate( grav(lonlength,latlength,altlength,timelength))
allocate(dmass(lonlength,latlength,altlength,timelength))
allocate(dmassmean(lonlength,latlength,altlength))
allocate(dmassmeanbar(latlength,altlength))

do itim=1,timelength
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
! Need to be consistent with GCM computations
    if (vitu(ilon,ilat,ilev,itim).ne.miss_val) then
     rayon(ilon,ilat,ilev,itim) = a0
!     rayon(ilon,ilat,ilev,itim) = za(ilon,ilat,ilev,itim) + a0
      grav(ilon,ilat,ilev,itim) = g0*a0*a0 &
                 /(rayon(ilon,ilat,ilev,itim)*rayon(ilon,ilat,ilev,itim))
    else
     rayon(ilon,ilat,ilev,itim) = miss_val
      grav(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo
enddo ! timelength

call cellmass(infid,latlength,lonlength,altlength,timelength,lmdflag, &
              miss_val,deltalon,deltalat,coslat,plev,ps,grav,rayon, &
              dmass )

call moytim(lonlength,latlength,altlength,timelength,miss_val,dmass,dmassmean)
call moyzon2(lonlength,latlength,altlength,miss_val,dmassmean,dmassmeanbar)

print*,"OK dmass"

!===============================================================================
! 2.2.3 Specific angular momentum
!===============================================================================
allocate(osam(lonlength,latlength,altlength,timelength))
allocate(rsam(lonlength,latlength,altlength,timelength))
allocate(tsam(lonlength,latlength,altlength,timelength))

do itim=1,timelength
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      osam(ilon,ilat,ilev,itim) = &
         rayon(ilon,ilat,ilev,itim)*rayon(ilon,ilat,ilev,itim) &
       * coslat(ilat)*coslat(ilat) &
       * omega
      rsam(ilon,ilat,ilev,itim) = vitu(ilon,ilat,ilev,itim) &
       * rayon(ilon,ilat,ilev,itim)*coslat(ilat)
      tsam(ilon,ilat,ilev,itim) = osam(ilon,ilat,ilev,itim)&
                                + rsam(ilon,ilat,ilev,itim)
    else
      osam(ilon,ilat,ilev,itim) = miss_val
      rsam(ilon,ilat,ilev,itim) = miss_val
      tsam(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo
enddo ! timelength

print*,"debut tprt ang"

!===============================================================================
! 2.2.4 Angular momentum transport
!===============================================================================
! allocations
!-------------
allocate(vang(lonlength,latlength,altlength,timelength))
allocate(wang(lonlength,latlength,altlength,timelength))
allocate(  vstar(lonlength,latlength,altlength,timelength))
allocate(  wstar(lonlength,latlength,altlength,timelength))
allocate(angstar(lonlength,latlength,altlength,timelength))
allocate(  vprim(lonlength,latlength,altlength,timelength))
allocate(  wprim(lonlength,latlength,altlength,timelength))
allocate(angprim(lonlength,latlength,altlength,timelength))
allocate(angprimvprim(lonlength,latlength,altlength,timelength))
allocate(angprimwprim(lonlength,latlength,altlength,timelength))

! lon,lat,alt
allocate(vangmean(lonlength,latlength,altlength))
allocate(wangmean(lonlength,latlength,altlength))
allocate(  vmean(lonlength,latlength,altlength))
allocate(  wmean(lonlength,latlength,altlength))
allocate(angmean(lonlength,latlength,altlength))
allocate(angprimvprimmean(lonlength,latlength,altlength))
allocate(angprimwprimmean(lonlength,latlength,altlength))
allocate(  vstarmean(lonlength,latlength,altlength))
allocate(  wstarmean(lonlength,latlength,altlength))
allocate(angstarmean(lonlength,latlength,altlength))
allocate(angstarmeanvstarmean(lonlength,latlength,altlength))
allocate(angstarmeanwstarmean(lonlength,latlength,altlength))
allocate(  v3d(lonlength,latlength,altlength))
allocate(  w3d(lonlength,latlength,altlength))
allocate(ang3d(lonlength,latlength,altlength))

! lat,alt,time
allocate(  vbar(latlength,altlength,timelength))
allocate(  wbar(latlength,altlength,timelength))
allocate(angbar(latlength,altlength,timelength))

!lat,alt
allocate(  vmeanbar(latlength,altlength))
allocate(  wmeanbar(latlength,altlength))
allocate(angmeanbar(latlength,altlength))
allocate(  vbar2d(latlength,altlength))
allocate(  wbar2d(latlength,altlength))
allocate(angbar2d(latlength,altlength))

allocate(totvang(latlength,altlength))
allocate(totwang(latlength,altlength))
allocate(mmcvang(latlength,altlength))
allocate(mmcwang(latlength,altlength))
allocate(trsvang(latlength,altlength))
allocate(trswang(latlength,altlength))
allocate(stnvang(latlength,altlength))
allocate(stnwang(latlength,altlength))

! intermediates
!-----------------

do itim=1,timelength
   v3d(:,:,:) = vitv(:,:,:,itim)
   w3d(:,:,:) = vitw(:,:,:,itim)
 ang3d(:,:,:) = tsam(:,:,:,itim)
 call moyzon2(lonlength,latlength,altlength,miss_val,  v3d,  vbar2d)
 call moyzon2(lonlength,latlength,altlength,miss_val,  w3d,  wbar2d)
 call moyzon2(lonlength,latlength,altlength,miss_val,ang3d,angbar2d)
   vbar(:,:,itim) =   vbar2d(:,:)
   wbar(:,:,itim) =   wbar2d(:,:)
 angbar(:,:,itim) = angbar2d(:,:)
enddo ! timelength

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
   do itim=1,timelength
    if ((vitv(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (vbar(ilat,ilev,itim)     .ne.miss_val)) then
   vstar(ilon,ilat,ilev,itim) = vitv(ilon,ilat,ilev,itim)-  vbar(ilat,ilev,itim)
    else
   vstar(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((vitw(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (wbar(ilat,ilev,itim)     .ne.miss_val)) then
   wstar(ilon,ilat,ilev,itim) = vitw(ilon,ilat,ilev,itim)-  wbar(ilat,ilev,itim)
    else
   wstar(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((  tsam(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (angbar(ilat,ilev,itim)     .ne.miss_val)) then
 angstar(ilon,ilat,ilev,itim) = tsam(ilon,ilat,ilev,itim)-angbar(ilat,ilev,itim)
    else
 angstar(ilon,ilat,ilev,itim) = miss_val
    endif
   enddo
  enddo
 enddo
enddo ! lonlength
call moytim(lonlength,latlength,altlength,timelength,miss_val,  vstar,  vstarmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,  wstar,  wstarmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,angstar,angstarmean)
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if ((angstarmean(ilon,ilat,ilev).ne.miss_val).and. &
        (  vstarmean(ilon,ilat,ilev).ne.miss_val)) then
angstarmeanvstarmean(ilon,ilat,ilev) = angstarmean(ilon,ilat,ilev)*vstarmean(ilon,ilat,ilev)
    else
angstarmeanvstarmean(ilon,ilat,ilev) = miss_val
    endif
    if ((angstarmean(ilon,ilat,ilev).ne.miss_val).and. &
        (  wstarmean(ilon,ilat,ilev).ne.miss_val)) then
angstarmeanwstarmean(ilon,ilat,ilev) = angstarmean(ilon,ilat,ilev)*wstarmean(ilon,ilat,ilev)
    else
angstarmeanwstarmean(ilon,ilat,ilev) = miss_val
    endif
  enddo
 enddo
enddo ! lonlength

call moytim(lonlength,latlength,altlength,timelength,miss_val,vitv,  vmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,vitw,  wmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,tsam,angmean)
call moyzon2(lonlength,latlength,altlength,miss_val,  vmean,  vmeanbar)
call moyzon2(lonlength,latlength,altlength,miss_val,  wmean,  wmeanbar)
call moyzon2(lonlength,latlength,altlength,miss_val,angmean,angmeanbar)

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
   do itim=1,timelength
    if ((vitv(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (tsam(ilon,ilat,ilev,itim).ne.miss_val)) then
vang(ilon,ilat,ilev,itim) = vitv(ilon,ilat,ilev,itim)*tsam(ilon,ilat,ilev,itim)
    else
vang(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((vitw(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (tsam(ilon,ilat,ilev,itim).ne.miss_val)) then
wang(ilon,ilat,ilev,itim) = vitw(ilon,ilat,ilev,itim)*tsam(ilon,ilat,ilev,itim)
    else
wang(ilon,ilat,ilev,itim) = miss_val
    endif
   enddo
  enddo
 enddo
enddo ! lonlength
call moytim(lonlength,latlength,altlength,timelength,miss_val,vang,vangmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,wang,wangmean)

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
   do itim=1,timelength
    if ((vitv(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (vmean(ilon,ilat,ilev)    .ne.miss_val)) then
  vprim(ilon,ilat,ilev,itim) = vitv(ilon,ilat,ilev,itim)-  vmean(ilon,ilat,ilev)
    else
  vprim(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((vitw(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (wmean(ilon,ilat,ilev)    .ne.miss_val)) then
  wprim(ilon,ilat,ilev,itim) = vitw(ilon,ilat,ilev,itim)-  wmean(ilon,ilat,ilev)
    else
  wprim(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((tsam(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (angmean(ilon,ilat,ilev)  .ne.miss_val)) then
angprim(ilon,ilat,ilev,itim) = tsam(ilon,ilat,ilev,itim)-angmean(ilon,ilat,ilev)
    else
angprim(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((angprim(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (  vprim(ilon,ilat,ilev,itim).ne.miss_val)) then
angprimvprim(ilon,ilat,ilev,itim) = angprim(ilon,ilat,ilev,itim)*vprim(ilon,ilat,ilev,itim)
    else
angprimvprim(ilon,ilat,ilev,itim) = miss_val
    endif
    if ((angprim(ilon,ilat,ilev,itim).ne.miss_val).and. &
        (  wprim(ilon,ilat,ilev,itim).ne.miss_val)) then
angprimwprim(ilon,ilat,ilev,itim) = angprim(ilon,ilat,ilev,itim)*wprim(ilon,ilat,ilev,itim)
    else
angprimwprim(ilon,ilat,ilev,itim) = miss_val
    endif
   enddo
  enddo
 enddo
enddo ! lonlength
call moytim(lonlength,latlength,altlength,timelength,miss_val,&
                  angprimvprim,angprimvprimmean)
call moytim(lonlength,latlength,altlength,timelength,miss_val,&
                  angprimwprim,angprimwprimmean)

! ang transport terms
!----------------------

call moyzon2(lonlength,latlength,altlength,miss_val,vangmean,totvang)
call moyzon2(lonlength,latlength,altlength,miss_val,wangmean,totwang)

do ilat=1,latlength
 do ilev=1,altlength
    if ((angmeanbar(ilat,ilev).ne.miss_val).and. &
        (  vmeanbar(ilat,ilev).ne.miss_val)) then
mmcvang(ilat,ilev) = angmeanbar(ilat,ilev)*vmeanbar(ilat,ilev)
    else
mmcvang(ilat,ilev) = miss_val
    endif
    if ((angmeanbar(ilat,ilev).ne.miss_val).and. &
        (  wmeanbar(ilat,ilev).ne.miss_val)) then
mmcwang(ilat,ilev) = angmeanbar(ilat,ilev)*wmeanbar(ilat,ilev)
    else
mmcwang(ilat,ilev) = miss_val
    endif
 enddo
enddo

call moyzon2(lonlength,latlength,altlength,miss_val,angprimvprimmean,trsvang)
call moyzon2(lonlength,latlength,altlength,miss_val,angprimwprimmean,trswang)

call moyzon2(lonlength,latlength,altlength,miss_val,angstarmeanvstarmean,stnvang)
call moyzon2(lonlength,latlength,altlength,miss_val,angstarmeanwstarmean,stnwang)


print*,"End of computations"

!===============================================================================
! 3. Create output file 
!===============================================================================

! Create output file
ierr=NF_CREATE(outfile,NF_CLOBBER,outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile
  stop
endif

!===============================================================================
! 3.1. Define and write dimensions
!===============================================================================

call write_dim(outfid,lonlength,latlength,altlength,timelength, &
    lon,lat,plev,time,lon_dimid,lat_dimid,alt_dimid,time_dimid)

!===============================================================================
! 3.2. Define and write variables
!===============================================================================

datashape2d(1)=lat_dimid
datashape2d(2)=alt_dimid

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "totvang   ", "tot hor trpt of ang ","m3 s-2    ",miss_val,&
                 totvang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "totwang   ", "tot ver trpt of ang ","m3 s-2    ",miss_val,&
                 totwang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "mmcvang   ", "MMC hor trpt of ang ","m3 s-2    ",miss_val,&
                 mmcvang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "mmcwang   ", "MMC ver trpt of ang ","m3 s-2    ",miss_val,&
                 mmcwang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "trsvang   ", "trs hor trpt of ang ","m3 s-2    ",miss_val,&
                 trsvang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "trswang   ", "trs ver trpt of ang ","m3 s-2    ",miss_val,&
                 trswang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "stnvang   ", "stn hor trpt of ang ","m3 s-2    ",miss_val,&
                 stnvang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "stnwang   ", "stn ver trpt of ang ","m3 s-2    ",miss_val,&
                 stnwang )

call write_var2d(outfid,datashape2d,latlength,altlength,&
                "dmass     ", "mass in each cell   ","kg        ",miss_val,&
                 dmassmeanbar )


!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile


end program
