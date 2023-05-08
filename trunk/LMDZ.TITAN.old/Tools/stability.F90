program stability

! SL 12/2009:
! This program reads 4D (lon-lat-alt-time) fields recast in log P coordinates
!
! it computes:
! stab  -- 4D -- stability
! Ri    -- 4D -- Richardson number
! deqc  -- 3D -- distance to cyclostrophic equilibrium
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - surface pressure and surface geopotential
! - atmospheric temperature
! - zonal and meridional winds
! - altitude above areoid

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (name_P.nc)
character (len=128) :: outfile ! output file name

character (len=64) :: text ! to store some text
integer infid ! NetCDF input file ID
integer outfid ! NetCDF output file ID
integer lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid
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
real,dimension(:,:,:,:),allocatable :: temp ! atmospheric temperature
real,dimension(:,:,:,:),allocatable :: vitu ! zonal wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)
real,dimension(:,:,:,:),allocatable :: za ! above areoid levels (m)

!!! output variables
real,dimension(:,:,:,:),allocatable :: stab ! stability (K/km)
real,dimension(:,:,:,:),allocatable :: Ri ! Richardson number
real,dimension(:,:,:),allocatable :: deqc ! distance to cyclostrophic equilibrium

! variables prepared for computation (4D)
real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)

! variables prepared for computation inside timeloop
real,dimension(:,:,:),allocatable :: t3d ! temp
real,dimension(:,:,:),allocatable :: u3d ! zonal wind
real,dimension(:,:,:),allocatable :: v3d ! merid wind
! variables obtained from computation inside timeloop
real,dimension(:,:,:),allocatable :: stab3d ! stability (K/km)
real,dimension(:,:,:),allocatable :: Ri3d ! Richardson number
real,dimension(:,:),allocatable :: deqc2d ! distance to cyclostrophic equilibrium

real,dimension(:,:),allocatable :: t2d   ! temp
real,dimension(:,:),allocatable :: u2d   ! zonal wind
real,dimension(:,:),allocatable :: v2d   ! merid wind
real,dimension(:,:),allocatable :: dtsdp ! d(temp)/d(plev)
real,dimension(:,:),allocatable :: dusdp ! du/d(plev)
real,dimension(:,:),allocatable :: dvsdp ! dv/d(plev)

integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
logical :: lmdflag

real :: deltalat,deltalon ! lat and lon intervals in radians
real :: fac1,ecden,ecnum ! for cyclo eq.

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
outfile=infile(1:len_trim(infile)-3)//"_STA.nc"
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
! 2.1.2 Atmospheric temperature
!===============================================================================
allocate(temp(lonlength,latlength,altlength,timelength))

text="temp"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for t instead... "
  text="t"
  call get_var4d(infid,lonlength,latlength,altlength,timelength,text,temp,miss_val,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) stop "Error: Failed to get temperature ID"
endif
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading temperature"

!===============================================================================
! 2.1.3 Winds
!===============================================================================
allocate(vitu(lonlength,latlength,altlength,timelength))
allocate(vitv(lonlength,latlength,altlength,timelength))

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
!!! Allocations before timeloop 
!===============================================================================

! latlength correspond a jjm+1
! mais lonlength correspond a iim
! pour boucler en longitude, on a besoin du point iim+1 (= 1)

allocate(rayon(lonlength+1,latlength,altlength,timelength))
allocate(grav(lonlength+1,latlength,altlength,timelength))
allocate(dmass(lonlength+1,latlength,altlength,timelength))

allocate(t3d(lonlength+1,latlength,altlength))
allocate(u3d(lonlength+1,latlength,altlength))
allocate(v3d(lonlength+1,latlength,altlength))

allocate(t2d(latlength,altlength))
allocate(u2d(latlength,altlength))
allocate(v2d(latlength,altlength))
allocate(dtsdp(latlength,altlength))
allocate(dusdp(latlength,altlength))
allocate(dvsdp(latlength,altlength))

allocate(stab(lonlength,latlength,altlength,timelength))
allocate(Ri(lonlength,latlength,altlength,timelength))
allocate(deqc(latlength,altlength,timelength))

allocate(stab3d(lonlength+1,latlength,altlength))
allocate(Ri3d(lonlength+1,latlength,altlength))
allocate(deqc2d(latlength,altlength))

!===============================================================================
! 2.2.2 Mass in cells
!===============================================================================

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

rayon(lonlength+1,:,:,:) = rayon(1,:,:,:)
 grav(lonlength+1,:,:,:) =  grav(1,:,:,:)

call cellmass(infid,latlength,lonlength+1,altlength,timelength,lmdflag, &
              miss_val,deltalon,deltalat,coslat,plev,ps,grav,rayon, &
              dmass )

!===============================================================================
!!! GLOBAL TIME LOOP !!!
!===============================================================================
do itim=1,timelength

!===============================================================================
! 2.2 Computations 
!===============================================================================

!===============================================================================
! 2.2.3 Init of 3D variables
!===============================================================================

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
      t3d(ilon,ilat,ilev) = temp(ilon,ilat,ilev,itim)
      u3d(ilon,ilat,ilev) = vitu(ilon,ilat,ilev,itim)
      v3d(ilon,ilat,ilev) = vitv(ilon,ilat,ilev,itim)
  enddo
 enddo
enddo

 t3d(lonlength+1,:,:) =  t3d(1,:,:)
 u3d(lonlength+1,:,:) =  u3d(1,:,:)
 v3d(lonlength+1,:,:) =  v3d(1,:,:)

!===============================================================================
! 2.2.4 Stability
!===============================================================================

do ilon=1,lonlength+1
 t2d(:,:) =  t3d(ilon,:,:)
 call dx_dp(latlength,altlength,miss_val,plev,t2d,dtsdp)
 do ilat=1,latlength
  do ilev=1,altlength
    if ((grav(ilon,ilat,ilev,itim).ne.miss_val).and. &
        ( t3d(ilon,ilat,ilev).ne.miss_val) ) then
      stab3d(ilon,ilat,ilev) = grav(ilon,ilat,ilev,itim)* &
                (  1./cpdet(t2d(ilat,ilev))          &
                 - plev(ilev)*dtsdp(ilat,ilev)/(R0*t2d(ilat,ilev)) )
      stab3d(ilon,ilat,ilev) = stab3d(ilon,ilat,ilev)*1000. ! passage en K/km
    else
      stab3d(ilon,ilat,ilev) = miss_val
    endif
  enddo
 enddo
enddo

!===============================================================================
! 2.2.5 Richardson number
!===============================================================================

do ilon=1,lonlength+1
 u2d(:,:) =  u3d(ilon,:,:)
 v2d(:,:) =  v3d(ilon,:,:)
 call dx_dp(latlength,altlength,miss_val,plev,u2d,dusdp)
 call dx_dp(latlength,altlength,miss_val,plev,v2d,dvsdp)
 do ilat=1,latlength
  do ilev=1,altlength
    if ((grav(ilon,ilat,ilev,itim).ne.miss_val).and. &
        ( u3d(ilon,ilat,ilev).ne.miss_val).and. &
        ( v3d(ilon,ilat,ilev).ne.miss_val).and. &
        ( t3d(ilon,ilat,ilev).ne.miss_val) ) then
      Ri3d(ilon,ilat,ilev) =  & ! attention, transfo a cause de du/dp au lieu de du/dz
          stab3d(ilon,ilat,ilev)*t3d(ilon,ilat,ilev)*R0*R0  &
      / (grav(ilon,ilat,ilev,itim)*plev(ilev)*plev(ilev))        &
      / (dusdp(ilat,ilev)*dusdp(ilat,ilev)+dvsdp(ilat,ilev)*dvsdp(ilat,ilev))
    else
      Ri3d(ilon,ilat,ilev) = miss_val
    endif
  enddo
 enddo
enddo

!===============================================================================
! 2.2.6 Distance to cyclostrophic equilibrium
!===============================================================================

call moyzon(lonlength,latlength,altlength,miss_val,u3d,u2d)
call moyzon(lonlength,latlength,altlength,miss_val,t3d,t2d)
call dx_dp(latlength,altlength,miss_val,plev,u2d,dusdp)

do ilat=2,latlength-1
   if (tan(lat(ilat)*pi/180.).ne.0.) then
     fac1 = R0/tan(lat(ilat)*pi/180.)
   else
     fac1 = miss_val
   endif
  do ilev=1,altlength
   if ((dusdp(ilat,ilev).ne.miss_val).and. &
       (  u2d(ilat,ilev).ne.miss_val).and. &
       (            fac1.ne.miss_val).and. &
       (  t2d(ilat,ilev).ne.miss_val) ) then
    ecden = dusdp(ilat,ilev)*(2.*u2d(ilat,ilev)*plev(ilev))
    ecnum = ((t2d(ilat+1,ilev)-t2d(ilat-1,ilev))/(2.*deltalat)*fac1-ecden)*100.
    deqc2d(ilat,ilev) = ecnum/ecden
   else
    deqc2d(ilat,ilev) = miss_val
   endif
  enddo
enddo
do ilev=1,altlength
    deqc2d(1,ilev)         = miss_val
    deqc2d(latlength,ilev) = miss_val
enddo

!===============================================================================
! 2.2.7 Building 3D+time variables
!===============================================================================

    deqc(:,:,itim)   = deqc2d(:,:)
    stab(:,:,:,itim) = stab3d(1:lonlength,:,:)
      Ri(:,:,:,itim) =   Ri3d(1:lonlength,:,:)


enddo ! timelength
!===============================================================================
!!! END GLOBAL TIME LOOP !!!
!===============================================================================

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

! 3D Variables

datashape3d(1)=lat_dimid
datashape3d(2)=alt_dimid
datashape3d(3)=time_dimid

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "deqc      ", "Distance to cyclo eq","per cent  ",miss_val,&
                 deqc )

! 4D Variables

datashape4d(1)=lon_dimid
datashape4d(2)=lat_dimid
datashape4d(3)=alt_dimid
datashape4d(4)=time_dimid

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "stab      ", "Stability           ","K/km      ",miss_val,&
                 stab )

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "Ri        ", "Richardson number   ","          ",miss_val,&
                 Ri )


!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
endif


end program
