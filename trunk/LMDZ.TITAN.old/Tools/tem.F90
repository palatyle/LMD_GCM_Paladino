program tem

! SL 01/2010:
! This program reads 4D (lon-lat-alt-time) fields recast in log P coordinates
! Developed from the tool built by Audrey Crespin during her PhD.
!
! it computes TransEulerianMean variables:
!
! vtem   -- 3D -- Residual meridional speed (m s-1)
! wtem   -- 3D -- Residual   vertical speed (Pa s-1)
! psitem -- 3D -- Residual stream function (kg s-1)
! epfy   -- 3D -- meridional component of Eliassen-Palm flux
! epfz   -- 3D -- vertical component of Eliassen-Palm flux
! divepf -- 3D -- Divergence of Eliassen-Palm flux
! ammctem - 3D -- Acc due to residual MMC
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - pressure vertical coordinate
! - surface pressure
! - atmospheric temperature
! - zonal, meridional and vertical (Pa/s) winds
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
real,dimension(:,:,:,:),allocatable :: vitw ! vertical wind (in Pa/s, then converted in m/s)
real,dimension(:,:,:,:),allocatable :: za ! above areoid levels (m)

!!! output variables
real,dimension(:,:,:),allocatable :: epy ! merid component of EP flux
real,dimension(:,:,:),allocatable :: epz ! verti component of EP flux
real,dimension(:,:,:),allocatable :: divep ! divergence of EP flux
real,dimension(:,:,:),allocatable :: ammctem ! acc by residual mmc
real,dimension(:,:,:),allocatable :: uzon ! mean zonal wind
real,dimension(:,:,:),allocatable :: vtem ! residual merid wind
real,dimension(:,:,:),allocatable :: wtem ! residual verti wind
real,dimension(:,:,:),allocatable :: psitem ! residual stream function

! variables prepared for computation (4D)
real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)

! variables prepared for computation inside timeloop
real,dimension(:,:,:),allocatable :: r3d ! distance to center (m)
real,dimension(:,:,:),allocatable :: rsurg ! rayon/grav
real,dimension(:,:,:),allocatable :: t3d ! temp
real,dimension(:,:,:),allocatable :: u3d ! zonal wind
real,dimension(:,:,:),allocatable :: v3d ! merid wind
real,dimension(:,:,:),allocatable :: w3d ! verti wind
real,dimension(:,:,:),allocatable :: pk3d ! Exner function
real,dimension(:,:,:),allocatable :: teta ! potential temp
! variables obtained from computation inside timeloop
real,dimension(:,:),allocatable :: epy2d ! merid component of EP flux
real,dimension(:,:),allocatable :: epz2d ! verti component of EP flux
real,dimension(:,:),allocatable :: div2d ! divergence of EP flux
real,dimension(:,:),allocatable :: ammc2d ! acc by residual mmc
real,dimension(:,:),allocatable :: ubar   ! mean zonal wind
real,dimension(:,:),allocatable :: vtem2d ! residual merid wind
real,dimension(:,:),allocatable :: wtem2d ! residual verti wind

real,dimension(:,:),allocatable :: rbar   ! distance to center (zonal ave)
real,dimension(:,:),allocatable :: rsurgbar ! rayon/grav
real,dimension(:,:),allocatable :: vm   ! merid mass flux (zonal ave)
real,dimension(:,:),allocatable :: psi  ! residual stream function
real :: deltalat,deltalon ! lat and lon intervals in radians
real,dimension(:,:,:),allocatable :: deltap ! pressure thickness of each layer (Pa)

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
outfile=infile(1:len_trim(infile)-3)//"_TEM.nc"
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
!!! Allocations before timeloop 
!===============================================================================

! latlength correspond a jjm+1
! mais lonlength correspond a iim
! pour boucler en longitude, on a besoin du point iim+1 (= 1)

allocate(rayon(lonlength+1,latlength,altlength,timelength))
allocate(grav(lonlength+1,latlength,altlength,timelength))
allocate(dmass(lonlength+1,latlength,altlength,timelength))

allocate(r3d(lonlength+1,latlength,altlength))
allocate(rsurg(lonlength+1,latlength,altlength))
allocate(rbar(latlength,altlength))
allocate(rsurgbar(latlength,altlength))

allocate(t3d(lonlength+1,latlength,altlength))
allocate(u3d(lonlength+1,latlength,altlength))
allocate(v3d(lonlength+1,latlength,altlength))
allocate(w3d(lonlength+1,latlength,altlength))
allocate(pk3d(lonlength+1,latlength,altlength))
allocate(teta(lonlength+1,latlength,altlength))

allocate(epy(latlength,altlength,timelength))
allocate(epz(latlength,altlength,timelength))
allocate(divep(latlength,altlength,timelength))
allocate(ammctem(latlength,altlength,timelength))
allocate(uzon(latlength,altlength,timelength))
allocate(vtem(latlength,altlength,timelength))
allocate(wtem(latlength,altlength,timelength))

allocate(epy2d(latlength,altlength))
allocate(epz2d(latlength,altlength))
allocate(div2d(latlength,altlength))
allocate(ammc2d(latlength,altlength))
allocate(ubar(latlength,altlength))
allocate(vtem2d(latlength,altlength))
allocate(wtem2d(latlength,altlength))

allocate(vm(latlength,altlength))
allocate(psi(latlength,altlength))
allocate(psitem(latlength,altlength,timelength))

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
      r3d(ilon,ilat,ilev) = rayon(ilon,ilat,ilev,itim)
      t3d(ilon,ilat,ilev) = temp(ilon,ilat,ilev,itim)
      u3d(ilon,ilat,ilev) = vitu(ilon,ilat,ilev,itim)
      v3d(ilon,ilat,ilev) = vitv(ilon,ilat,ilev,itim)
      w3d(ilon,ilat,ilev) = vitw(ilon,ilat,ilev,itim)
     pk3d(ilon,ilat,ilev) = cp0*(plev(ilev)/psref)**(R0/cp0)
  enddo
 enddo
enddo

 t3d(lonlength+1,:,:) =  t3d(1,:,:)
 u3d(lonlength+1,:,:) =  u3d(1,:,:)
 v3d(lonlength+1,:,:) =  v3d(1,:,:)
 w3d(lonlength+1,:,:) =  w3d(1,:,:)
pk3d(lonlength+1,:,:) = pk3d(1,:,:)

call t2tpot((lonlength+1)*latlength*altlength,t3d,teta,pk3d)

!===============================================================================
! 2.2.4 TEM and Eliassen-Palm
!===============================================================================

print*,"eliasflu_meridien",itim

call moyzon(lonlength,latlength,altlength,miss_val,r3d,rbar)
call moyzon(lonlength,latlength,altlength,miss_val,u3d,ubar)

call epflux(lonlength+1,latlength,altlength,miss_val,lat,rbar &
            ,teta,u3d,v3d,w3d,plev &
            ,epy2d,epz2d,div2d,vtem2d,wtem2d,ammc2d &
!           ,vpupbar2d,wpupbar2d,vpvpbar2d,wpvpbar2d &
!           ,vptetapbar2d,wptetapbar2d &
           )

!===============================================================================
! 2.2.5 Stream function
!===============================================================================

do ilon=1,lonlength+1
 do ilat=1,latlength
  do ilev=1,altlength
    if (dmass(ilon,ilat,ilev,itim).ne.miss_val) then
! rsurg: r*dp/g = dm/(r cos(lat) dlon dlat) !!!
     rsurg(ilon,ilat,ilev) = dmass(ilon,ilat,ilev,itim) &
          / (r3d(ilon,ilat,ilev)*coslat(ilat)*deltalon*deltalat)
    else
     rsurg(ilon,ilat,ilev) = miss_val
    endif
  enddo
 enddo
enddo

call moyzon(lonlength,latlength,altlength,miss_val,rsurg,rsurgbar)

do ilat=1,latlength
 do ilev=1,altlength
    if (  (vtem2d(ilat,ilev).ne.miss_val).and. &
        (rsurgbar(ilat,ilev).ne.miss_val) ) then
      vm(ilat,ilev) = vtem2d(ilat,ilev) &
            * 2.*pi*rsurgbar(ilat,ilev)*coslat(ilat)
    else
      vm(ilat,ilev) = miss_val
    endif
 enddo
enddo


do ilat=1,latlength
  psi(ilat,altlength) = 0.
    if (vm(ilat,altlength).ne.miss_val) then
      psi(ilat,altlength) = psi(ilat,altlength) &
           + vm(ilat,altlength)
    endif
 do ilev=altlength-1,1,-1
  psi(ilat,ilev) = psi(ilat,ilev+1)
    if (vm(ilat,ilev).ne.miss_val) then
      psi(ilat,ilev) = psi(ilat,ilev) &
           + vm(ilat,ilev)
    endif
 enddo
enddo

!===============================================================================
! 2.2.6 Building 2D+time variables
!===============================================================================

    epy(:,:,itim) = epy2d(:,:)
    epz(:,:,itim) = epz2d(:,:)
  divep(:,:,itim) = div2d(:,:)
ammctem(:,:,itim) = ammc2d(:,:)
   uzon(:,:,itim) =   ubar(:,:)
   vtem(:,:,itim) = vtem2d(:,:)
   wtem(:,:,itim) = wtem2d(:,:)
 psitem(:,:,itim) = psi(:,:)


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

datashape3d(1)=lat_dimid
datashape3d(2)=alt_dimid
datashape3d(3)=time_dimid

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "epy       ", "EP flux on lat      ","m3 s-2    ",miss_val,&
                 epy )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "epz       ", "EP flux on press    ","m3 s-2    ",miss_val,&
                 epz )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "divep     ", "Div of EP flux      ","m s-2     ",miss_val,&
                 divep )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "ammctem   ", "acc by residual mmc ","m s-2     ",miss_val,&
                 ammctem )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "uzon      ", "Mean zonal wind     ","m s-1     ",miss_val,&
                 uzon )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "vtem      ", "Resid TEM merid wind","m s-1     ",miss_val,&
                 vtem )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "wtem      ", "Resid TEM verti wind","Pa s-1    ",miss_val,&
                 wtem )

call write_var3d(outfid,datashape3d,latlength,altlength,timelength,&
                "psitem    ", "Resid stream funct  ","kg s-1    ",miss_val,&
                 psitem )


!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) write(*,*) 'Error, failed to close output file ',outfile


end program
