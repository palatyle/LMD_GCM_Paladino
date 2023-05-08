program streamfunction

! SL 12/2009:
! This program reads 4D (lon-lat-alt-time) fields directly from LMD outputs
! without regrid : histmth OR from files recast in log P coordinates (_P)
!
! it computes:
! dmass -- 4D -- mass of each cell
! psi   -- 3D -- Stream function
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - surface pressure and surface geopotential
! - meridional winds

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (name_P.nc)
character (len=128) :: outfile ! output file name

character (len=64) :: text ! to store some text
integer infid ! NetCDF input file ID
integer outfid ! NetCDF output file ID
integer lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid
integer              :: datashape1d ! shape of 1D datasets
integer,dimension(2) :: datashape2d ! shape of 2D datasets
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
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)
real,dimension(:,:,:,:),allocatable :: za   ! areoid levels (m)

real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)
real,dimension(:,:,:,:),allocatable :: vm ! meridional mass flux
real,dimension(:,:,:),allocatable :: psi ! stream function

integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
logical :: lmdflag

real :: deltalat,deltalon ! lat and lon intervals in radians

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
write(*,*) " Program valid for Venus or Titan LMD output files"
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

! Time axis IN PLANET DAYS

time=time/localday

!===============================================================================
! 1.3 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
outfile=infile(1:len_trim(infile)-3)//"_PSI.nc"
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
allocate(vitv(lonlength,latlength,altlength,timelength))

! meridional wind vitv (in m/s)
text="vitv"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitv,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitv ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading meridional wind"

!===============================================================================
! 2.1.4 Altitude above areoide
!===============================================================================

allocate(za(lonlength,latlength,altlength,timelength))
! present only in _P regrided files
! For others, using g0*a0*a0/(g0*a0-geop)-a0

text="zareoid"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,za,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for geop instead... "
  text="geop"
  call get_var4d(infid,lonlength,latlength,altlength,timelength,text,za,miss_val,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) stop "Error: Failed to get geop ID"
do itim=1,timelength
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (za(ilon,ilat,ilev,itim).ne.miss_val) then
        za(ilon,ilat,ilev,itim) = (g0*a0*a0/(g0*a0-za(ilon,ilat,ilev,itim))-a0)/1000. ! in km
    else
        za(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo
enddo
endif
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zareoid/geop"

!===============================================================================
! 2.2 Computations 
!===============================================================================

!===============================================================================
! 2.2.2 Mass in cells
!===============================================================================
allocate(rayon(lonlength,latlength,altlength,timelength))
allocate(grav(lonlength,latlength,altlength,timelength))
allocate(dmass(lonlength,latlength,altlength,timelength))

do itim=1,timelength
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
! Need to be consistent with GCM computations
    if (za(ilon,ilat,ilev,itim).ne.miss_val) then
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

!===============================================================================
! 2.2.8 Stream function
!===============================================================================
allocate(vm(lonlength,latlength,altlength,timelength))
allocate(psi(latlength,altlength,timelength))

do itim=1,timelength

do ilat=1,latlength
 do ilev=1,altlength
  do ilon=1,lonlength
    if (dmass(ilon,ilat,ilev,itim).ne.miss_val) then
      vm(ilon,ilat,ilev,itim) = vitv(ilon,ilat,ilev,itim)  &
                              * dmass(ilon,ilat,ilev,itim) &
                              / (rayon(ilon,ilat,ilev,itim)*deltalat)
    else
      vm(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo


do ilat=1,latlength
  psi(ilat,altlength,itim) = 0.
  do ilon=1,lonlength
    if (vm(ilon,ilat,altlength,itim).ne.miss_val) then
      psi(ilat,altlength,itim) = psi(ilat,altlength,itim) &
           + vm(ilon,ilat,altlength,itim)
    endif
  enddo
 do ilev=altlength-1,1,-1
  psi(ilat,ilev,itim) = psi(ilat,ilev+1,itim)
  do ilon=1,lonlength
    if (vm(ilon,ilat,ilev,itim).ne.miss_val) then
      psi(ilat,ilev,itim) = psi(ilat,ilev,itim) &
           + vm(ilon,ilat,ilev,itim)
    endif
  enddo
 enddo
enddo

enddo ! timelength

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

! 1D Variables

datashape1d=time_dimid
 
! 3D variables

datashape3d(1)=lon_dimid
datashape3d(2)=lat_dimid
datashape3d(3)=time_dimid

call write_var3d(outfid,datashape3d,lonlength,latlength,timelength,&
                "ps        ", "Surface pressure    ","Pa        ",miss_val,&
                 ps )

datashape3d(1)=lat_dimid
datashape3d(2)=alt_dimid
datashape3d(3)=time_dimid

call write_var3d(outfid,datashape3d,lonlength,latlength,timelength,&
                "psi       ", "Stream function     ","kg/s      ",miss_val,&
                 psi )

! 4D variables

datashape4d(1)=lon_dimid
datashape4d(2)=lat_dimid
datashape4d(3)=alt_dimid
datashape4d(4)=time_dimid

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "dmass     ", "Mass                ","kg        ",miss_val,&
                 dmass )

!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
endif


end program
