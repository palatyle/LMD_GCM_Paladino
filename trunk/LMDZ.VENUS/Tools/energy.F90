program energy

! SL 12/2009:
! This program reads 4D (lon-lat-alt-time) fields directly from LMD outputs without regrid : histmth 
!
! it computes:
! dmass -- 4D -- mass of each cell
! sek   -- 4D -- specific kinetic energy
! ek    -- 1D -- integrated kinetic energy
! sep   -- 4D -- specific potential energy
! ep    -- 1D -- integrated potential energy
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - surface pressure
! - atmospheric temperature
! - zonal and meridional winds

! VERTICAL WIND SPEED IS NEGLECTED IN KINETIC ENERGY

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
integer altlength ! # of grid point along presnivs (of input datasets)
real,dimension(:),allocatable :: time ! time
integer timelength ! # of points along time
real,dimension(:,:,:),allocatable :: ps ! surface pressure
real,dimension(:,:,:,:),allocatable :: temp ! atmospheric temperature
real,dimension(:,:,:,:),allocatable :: vitu ! zonal wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)

real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)
real,dimension(:,:,:,:),allocatable :: sek ! specific kinetic energy
real,dimension(:),allocatable :: ek ! total kinetic energy
real,dimension(:,:,:,:),allocatable :: sep ! specific potential energy
real,dimension(:),allocatable :: ep ! total potential energy

integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops

real :: deltalat,deltalon ! lat and lon intervals in radians
real,dimension(:,:,:,:),allocatable :: deltap ! pressure thickness of each layer (Pa)
real :: tmpp ! temporary pressure
real :: dz ! altitude diff
real :: signe ! orientation of lon axis for mountain torque computation
logical :: lmdflag

real :: cpdet

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
write(*,*) " Program valid for Venus or Titan LMD, or Venus CAM output files"
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
if(lon(1).gt.lon(2)) then
  signe=-1.
else
  signe=1.
endif

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
if (ierr.ne.NF_NOERR) stop "Error: Failed reading presnivs (ie pressure levels)"
if(.not.lmdflag) then
! in CAM files, pressure in mbar and reversed...
  call reverselev(altlength,plev)
  plev=plev*100.  ! convertion to Pa
endif

allocate(time(timelength))
ierr=NF_GET_VAR_REAL(infid,time_varid,time)
if (ierr.ne.NF_NOERR) stop "Error: Failed reading time"

! Time axis IN PLANET DAYS

if(.not.lmdflag) then
! in CAM files, time in Earth days...
!   => seconds
  time=time*86400.
endif
time=time/localday

!===============================================================================
! 1.3 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
outfile=infile(1:len_trim(infile)-3)//"_NRG.nc"
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
if((.not.lmdflag).and.(planet.eq."Venus")) call reverse3d(lonlength,latlength,timelength,ps)

!===============================================================================
! 2.1.2 Atmospheric temperature
!===============================================================================
allocate(temp(lonlength,latlength,altlength,timelength))

if(lmdflag) then
  text="temp"
else
  text="T"
endif
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for t instead... "
  text="t"
  call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) stop "Error: Failed to get temperature ID"
endif
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading temperature"

if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,temp)

!===============================================================================
! 2.1.3 Winds
!===============================================================================
allocate(vitu(lonlength,latlength,altlength,timelength))
allocate(vitv(lonlength,latlength,altlength,timelength))

! zonal wind vitu (in m/s)
if(lmdflag) then
  text="vitu"
else
  text="U"
endif
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitu,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitu ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zonal wind"

if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,vitu)

! meridional wind vitv (in m/s)
if(lmdflag) then
  text="vitv"
else
  text="V"
endif
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitv,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get vitv ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading meridional wind"

if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,vitv)

!===============================================================================
! 2.1.4 Altitude above areoide
!===============================================================================
! Only needed if g(z) on Titan...

!allocate(za(lonlength,latlength,altlength,timelength))

!text="zareoid"
!call get_var4d(infid,lonlength,latlength,altlength,timelength,text,za,miss_val,ierr1,ierr2)
!if (ierr1.ne.NF_NOERR) stop "Error: Failed to get za ID"
!if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zareoid"

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
!    if (za(ilon,ilat,ilev,itim).ne.miss_val) then
     rayon(ilon,ilat,ilev,itim) = a0
!     rayon(ilon,ilat,ilev,itim) = za(ilon,ilat,ilev,itim) + a0
      grav(ilon,ilat,ilev,itim) = g0*a0*a0 &
                 /(rayon(ilon,ilat,ilev,itim)*rayon(ilon,ilat,ilev,itim))
!    else
!     rayon(ilon,ilat,ilev,itim) = miss_val
!      grav(ilon,ilat,ilev,itim) = miss_val
!    endif
  enddo
 enddo
enddo
enddo ! timelength

call cellmass(infid,latlength,lonlength,altlength,timelength,lmdflag, &
              miss_val,deltalon,deltalat,coslat,plev,ps,grav,rayon, &
              dmass )

!===============================================================================
! 2.2.6 Specific energies
!===============================================================================
allocate(sek(lonlength,latlength,altlength,timelength))
allocate(sep(lonlength,latlength,altlength,timelength))

do itim=1,timelength

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      if ((vitu(ilon,ilat,ilev,itim).lt.miss_val) &
     .and.(vitv(ilon,ilat,ilev,itim).lt.miss_val)) then
      sek(ilon,ilat,ilev,itim) = 0.5 * &
       ( vitu(ilon,ilat,ilev,itim)*vitu(ilon,ilat,ilev,itim) &
       + vitv(ilon,ilat,ilev,itim)*vitv(ilon,ilat,ilev,itim) )
       else
      sek(ilon,ilat,ilev,itim) = miss_val
       endif
      if (temp(ilon,ilat,ilev,itim).ne.miss_val) then
      sep(ilon,ilat,ilev,itim) = temp(ilon,ilat,ilev,itim) &
                             * cpdet(temp(ilon,ilat,ilev,itim))
       else
      sep(ilon,ilat,ilev,itim) = miss_val
       endif
    else
      sek(ilon,ilat,ilev,itim) = miss_val
      sep(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo

enddo ! timelength

!===============================================================================
! 2.2.7 Total energies
!===============================================================================
allocate(ek(timelength))
allocate(ep(timelength))

do itim=1,timelength

ek(itim) = 0.
ep(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (sek(ilon,ilat,ilev,itim).ne.miss_val) then
      ek(itim) = ek(itim) &
       + sek(ilon,ilat,ilev,itim) * dmass(ilon,ilat,ilev,itim)
    endif
    if (sep(ilon,ilat,ilev,itim).ne.miss_val) then
      ep(itim) = ep(itim) &
       + sep(ilon,ilat,ilev,itim) * dmass(ilon,ilat,ilev,itim)
    endif
  enddo
 enddo
enddo
if (ek(itim).eq.0.) then
  ek(itim) = miss_val
  ep(itim) = miss_val
endif

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
 
call write_var1d(outfid,datashape1d,timelength,&
                "ek        ", "Total kinetic energy","J         ",miss_val,&
                 ek )

call write_var1d(outfid,datashape1d,timelength,&
                "ep        ", "Total pot energy    ","J         ",miss_val,&
                 ep )

! 3D variables

datashape3d(1)=lon_dimid
datashape3d(2)=lat_dimid
datashape3d(3)=time_dimid

call write_var3d(outfid,datashape3d,lonlength,latlength,timelength,&
                "ps        ", "Surface pressure    ","Pa        ",miss_val,&
                 ps )

! 4D variables

datashape4d(1)=lon_dimid
datashape4d(2)=lat_dimid
datashape4d(3)=alt_dimid
datashape4d(4)=time_dimid

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "dmass     ", "Mass                ","kg        ",miss_val,&
                 dmass )

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "sek       ", "Specific kin energy ","J/kg      ",miss_val,&
                 sek )

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "sep       ", "Specific pot energy ","J/kg      ",miss_val,&
                 sep )


!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
endif


end program
