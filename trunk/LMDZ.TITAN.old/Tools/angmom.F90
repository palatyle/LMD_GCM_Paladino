program angmom

! SL 12/2009:
! This program reads 4D (lon-lat-alt-time) fields directly from LMD (or CAM) outputs
! without regrid : histmth OR from files recast in log P coordinates (_P)
! (+ dynzon for LMD if present, but beware of recast consistency)
!
! it computes:
! dmass -- 4D -- mass of each cell
! osam  -- 4D -- specific angular momentum (omega term)
! rsam  -- 4D -- specific angular momentum (zonal wind term)
! oaam  -- 1D -- integrated angular momentum (omega term)
! raam  -- 1D -- integrated angular momentum (zonal wind term)
! tmou  -- 1D -- Mountain torque
! tbls  -- 1D -- Surface friction torque IF duvdf is present 
!                                       or if simple friction
! tdyn  -- 1D -- Dynamics torque IF dudyn is present 
! tajs  -- 1D -- Torque from convective adjustment IF dudyn is present 
! tgwo  -- 1D -- Orographic GW torque IF dugwo is present 
! tgwno -- 1D -- Non-Orographic GW torque IF dugwno is present 
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - surface pressure and surface geopotential
! - zonal wind
! Optional: dudyn, duvdf, duajs, dugwo, dugwno (acceleration terms from physiq param)

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name
character (len=128) :: dzfile ! input dynzon file name
character (len=128) :: outfile ! output file name

character (len=64) :: text ! to store some text
integer infid ! NetCDF input file ID
integer dzfid ! NetCDF input dynzon file ID
integer outfid ! NetCDF output file ID
integer lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs
integer latdz_dimid,altdz_dimid,timedz_dimid ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid, tmpvarid
integer latdz_varid,altdz_varid,timedz_varid
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
real,dimension(:,:,:),allocatable :: phis3d ! surface geopotential
real,dimension(:,:),allocatable :: phis ! surface geopotential
real,dimension(:,:,:,:),allocatable :: temp ! atmospheric temperature
real,dimension(:,:,:,:),allocatable :: vitu ! zonal wind (in m/s)
real,dimension(:,:,:,:),allocatable :: vitv ! meridional wind (in m/s)

real,dimension(:,:,:,:),allocatable :: duvdf  ! Friction in BL
real,dimension(:,:,:,:),allocatable :: dudyn  ! Dynamics 
real,dimension(:,:,:,:),allocatable :: duajs  ! Convective adjustment
real,dimension(:,:,:,:),allocatable :: dugwo  ! Orographic Gravity Waves
real,dimension(:,:,:,:),allocatable :: dugwno ! Non-Orographic Gravity Waves

real,dimension(:,:,:),allocatable :: dmcdyn  ! Dynamics dAM from dynzon
real,dimension(:,:,:),allocatable :: dmcdis  ! Dissip dAM from dynzon
real,dimension(:,:,:),allocatable :: dmcspg  ! Sponge dAM from dynzon
real,dimension(:,:,:),allocatable :: dmcphy  ! Physics dAM from dynzon

real,dimension(:,:,:,:),allocatable :: rayon ! distance to center (m)
real,dimension(:,:,:,:),allocatable :: grav ! gravity field (m s-2)
real,dimension(:,:,:,:),allocatable :: dmass ! mass in cell (kg)
real,dimension(:,:,:,:),allocatable :: osam ! planetary rotation specific ang (m2/s)
real,dimension(:,:,:,:),allocatable :: rsam ! zonal wind specific ang (m2/s)
real,dimension(:),allocatable :: oaam ! planetary rotation total ang (kg m2/s)
real,dimension(:),allocatable :: raam ! zonal wind total ang (kg m2/s)
real,dimension(:),allocatable :: tmou ! mountain torque (kg m2/s2)
real,dimension(:),allocatable :: tdyn ! dynamics torque (kg m2/s2)
real,dimension(:),allocatable :: tajs ! convective adjustment torque (kg m2/s2)
real,dimension(:),allocatable :: tbls ! friction torque (kg m2/s2)
real,dimension(:),allocatable :: tgwo ! oro GW torque (kg m2/s2)
real,dimension(:),allocatable :: tgwno! non-oro GW torque (kg m2/s2)
real,dimension(:),allocatable :: tdyndz ! dynamics torque (kg m2/s2) from dynzon
real,dimension(:),allocatable :: tdisdz ! dissip torque (kg m2/s2) from dynzon
real,dimension(:),allocatable :: tspgdz ! sponge torque (kg m2/s2) from dynzon
real,dimension(:),allocatable :: tphydz ! physics torque (kg m2/s2) from dynzon

! for angular momentum budget normalisation
real,parameter :: hadley=1.e18
real,parameter :: hadday=1.e25

integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
integer idlsurf ! for option ideal surface
logical :: flag_duvdf,flag_dudyn,flag_duajs,flag_dugwo,flag_dugwno,lmdflag,dzflag

real :: deltalat,deltalon ! lat and lon intervals in radians
real :: tmpp ! temporary pressure
real :: dz ! altitude diff
real :: signe ! orientation of lon axis for mountain torque computation
real :: norm ! for dynzon

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
if (ierr.ne.NF_NOERR) stop "Error: Failed reading pressure levels"
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
! 1.3 dynzon file if present
!===============================================================================

! RAJOUTER UN TEST COHERENCE _P...

dzflag=.false.

if(lmdflag) then

write(*,*) "Enter dynzon file name (<return> if not present):"

read(*,'(a128)') dzfile
write(*,*) ""

if(dzfile.ne."") then 

! open dynzon file
ierr = NF_OPEN(dzfile,NF_NOWRITE,dzfid)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening file ',trim(dzfile)
else
   dzflag=.true.
endif

endif ! dzfile.ne.""
endif ! lmdflag

!===============================================================================
! 1.4 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
outfile=infile(1:len_trim(infile)-3)//"_GAM.nc"
write(*,*) "Output file name is: "//trim(outfile)



!===============================================================================
! 2.1 Store needed fields 
!===============================================================================

!===============================================================================
! 2.1.1 Surface pressure and geopotential
!===============================================================================
allocate(ps(lonlength,latlength,timelength))
allocate(phis3d(lonlength,latlength,timelength))
allocate(phis(lonlength,latlength))

text="PS"
call get_var3d(infid,lonlength,latlength,timelength,text,ps,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "  looking for ps instead... "
  text="ps"
  call get_var3d(infid,lonlength,latlength,timelength,text,ps,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) then
    write(*,*) "  looking for psol instead... "
    text="psol"
    call get_var3d(infid,lonlength,latlength,timelength,text,ps,ierr1,ierr2)
    if (ierr1.ne.NF_NOERR) stop "Error: Failed to get psol ID"
  endif
endif
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading surface pressure"
if((.not.lmdflag).and.(planet.eq."Venus")) call reverse3d(lonlength,latlength,timelength,ps)

text="PHIS"
call get_var3d(infid,lonlength,latlength,timelength,text,phis3d,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) " Failed to get PHIS ID (3d), looking for phis (2d) instead... "
  text="phis"
  call get_var2d(infid,lonlength,latlength,text,phis,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) stop "Error: Failed to get phis ID"
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading surface geopotential"
else
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading surface geopotential"
  phis(:,:)=phis3d(:,:,1)
  if((.not.lmdflag).and.(planet.eq."Venus")) call reverse2d(lonlength,latlength,phis)
endif

!===============================================================================
! 2.1.3 Winds
!===============================================================================
allocate(vitu(lonlength,latlength,altlength,timelength))

! zonal wind vitu / U (in m/s)
if(lmdflag) then
  text="vitu"
else
  text="U"
endif

call get_var4d(infid,lonlength,latlength,altlength,timelength,text,vitu,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get U ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading zonal wind"

if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,vitu)

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
! 2.1.5 Friction in Boundary layer
!===============================================================================
allocate(duvdf(lonlength,latlength,altlength,timelength))

if(lmdflag) then
  text="duvdf"
else
  text="DUVDF"
endif
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,duvdf,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Failed to get duvdf ID"
  flag_duvdf = .false.

if(.not.lmdflag) then
! IDEAL FRICTION ?
write(*,*) ""
write(*,*) " Is the friction at surface ideal ? 0 for no, 1 for yes"
write(*,*) " duvdf = -u/(86400*30) in surface layer"
write(*,*) " timescale hard-coded: 30 Edays"
read(*,'(i1)') idlsurf
!write(*,*) ""
!write(*,*) " ASSUMED YES ! "
!idlsurf=1
write(*,*) ""
else
idlsurf=0
endif

if (idlsurf.eq.1) then
  flag_duvdf = .true.
  duvdf = 0.
  ilev=1
do ilon=1,lonlength
 do ilat=1,latlength
  do itim=1,timelength
   duvdf(ilon,ilat,ilev,itim) = vitu(ilon,ilat,ilev,itim) * (-1.)/(86400.*30.)
  enddo
 enddo
enddo
endif

else !err1
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading duvdf"
  if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,duvdf)
  flag_duvdf = .true.
endif !err1

!===============================================================================
! 2.1.6 Orographic and Non-orographic gravity waves
!===============================================================================
allocate(dugwo(lonlength,latlength,altlength,timelength))
allocate(dugwno(lonlength,latlength,altlength,timelength))

if(lmdflag) then

text="dugwo"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,dugwo,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Failed to get dugwo ID"
  flag_dugwo = .false.
else
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dugwo"
  flag_dugwo = .true.
endif

text="dugwno"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,dugwno,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Failed to get dugwno ID"
  flag_dugwno = .false.
else
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dugwno"
  flag_dugwno = .true.
endif

else   ! lmdflag
  print*,"dugwo and dugwno not in CAM simulations"
  flag_dugwo = .false.
  flag_dugwno = .false.
endif  ! lmdflag

!===============================================================================
! 2.1.65 Accelerations from convective adjustment
!===============================================================================
allocate(duajs(lonlength,latlength,altlength,timelength))

if(lmdflag) then

text="duajs"
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,duajs,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Failed to get duajs ID"
  flag_duajs = .false.
else
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading duajs"
  flag_duajs = .true.
endif

else   ! lmdflag
  print*,"duajs not in CAM simulations"
  flag_duajs = .false.
endif  ! lmdflag

!===============================================================================
! 2.1.7 Dynamics (includes the mountain torque...)
!===============================================================================
allocate(dudyn(lonlength,latlength,altlength,timelength))

if(lmdflag) then
  text="dudyn"
else
  text="DUDYN"
endif
call get_var4d(infid,lonlength,latlength,altlength,timelength,text,dudyn,miss_val,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Failed to get dudyn ID"
  flag_dudyn = .false.
else
  if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dudyn"
  if(.not.lmdflag) call reverse4d(lonlength,latlength,altlength,timelength,dudyn)
  flag_dudyn = .true.
endif

!===============================================================================
! 2.1.8 Dynzon dAM/dt
!===============================================================================
if(dzflag) then

allocate(dmcdyn(latlength-1,altlength,timelength))
allocate(dmcdis(latlength-1,altlength,timelength))
allocate(dmcspg(latlength-1,altlength,timelength))
allocate(dmcphy(latlength-1,altlength,timelength))

  text="dmcdyn"
call get_var3d(dzfid,latlength-1,altlength,timelength,text,dmcdyn,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get dmcdyn ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dmcdyn"

  text="dmcdis"
call get_var3d(dzfid,latlength-1,altlength,timelength,text,dmcdis,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get dmcdis ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dmcdis"

  text="dmcspg"
call get_var3d(dzfid,latlength-1,altlength,timelength,text,dmcspg,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get dmcspg ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dmcspg"

  text="dmcphy"
call get_var3d(dzfid,latlength-1,altlength,timelength,text,dmcphy,ierr1,ierr2)
if (ierr1.ne.NF_NOERR) stop "Error: Failed to get dmcphy ID"
if (ierr2.ne.NF_NOERR) stop "Error: Failed reading dmcphy"

endif ! dzflag

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

!===============================================================================
! 2.2.3 Specific angular momentum
!===============================================================================
allocate(osam(lonlength,latlength,altlength,timelength))
allocate(rsam(lonlength,latlength,altlength,timelength))

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
    else
      osam(ilon,ilat,ilev,itim) = miss_val
      rsam(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo

enddo ! timelength

!===============================================================================
! 2.2.4 Total angular momentum
!===============================================================================
allocate(oaam(timelength))
allocate(raam(timelength))

do itim=1,timelength

oaam(itim) = 0.
raam(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      oaam(itim) = oaam(itim) &
       + osam(ilon,ilat,ilev,itim)/ hadday * dmass(ilon,ilat,ilev,itim)
      raam(itim) = raam(itim) &
       + rsam(ilon,ilat,ilev,itim)/ hadday * dmass(ilon,ilat,ilev,itim)
    endif
  enddo
 enddo
enddo
if (oaam(itim).eq.0.) then
  oaam(itim) = miss_val
  raam(itim) = miss_val
endif

enddo ! timelength

!===============================================================================
! 2.2.5 Mountain, friction, convective adjustment, dynamics and GW torques
!===============================================================================
allocate(tmou(timelength))
if (flag_dudyn)  allocate(tdyn(timelength))
if (flag_duajs)  allocate(tajs(timelength))
if (flag_duvdf)  allocate(tbls(timelength))
if (flag_dugwo)  allocate(tgwo(timelength))
if (flag_dugwno) allocate(tgwno(timelength))

do itim=1,timelength

tmou(itim) = 0.
do ilon=1,lonlength
 do ilat=2,latlength-1
      if (ilon.eq.lonlength) then
         dz = (phis(   1,ilat)     -phis(ilon,ilat))/g0
       tmpp = (ps(     1,ilat,itim)  +ps(ilon,ilat,itim))/2.
      else
         dz = (phis(ilon+1,ilat)   -phis(ilon,ilat))/g0
       tmpp = (ps(ilon+1,ilat,itim)  +ps(ilon,ilat,itim))/2.
      endif
      tmou(itim) = tmou(itim) + a0*a0* deltalat*coslat(ilat) &
       * signe*dz*tmpp &
       / hadley
 enddo
enddo

if (flag_dudyn) then
tdyn(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      tdyn(itim) = tdyn(itim) + dudyn(ilon,ilat,ilev,itim) &
       * rayon(ilon,ilat,ilev,itim)*coslat(ilat)      &
       * dmass(ilon,ilat,ilev,itim) &
       / hadley
    endif
  enddo
 enddo
enddo
  if (tdyn(itim).eq.0.) tdyn(itim) = miss_val
endif

if (flag_duajs) then
tajs(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      tajs(itim) = tajs(itim) + duajs(ilon,ilat,ilev,itim) &
       * rayon(ilon,ilat,ilev,itim)*coslat(ilat)      &
       * dmass(ilon,ilat,ilev,itim) &
       / hadley
    endif
  enddo
 enddo
enddo
  if (tajs(itim).eq.0.) tajs(itim) = miss_val
endif

if (flag_duvdf) then
tbls(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      tbls(itim) = tbls(itim) + duvdf(ilon,ilat,ilev,itim) &
       * rayon(ilon,ilat,ilev,itim)*coslat(ilat)      &
       * dmass(ilon,ilat,ilev,itim) &
       / hadley
    endif
  enddo
 enddo
enddo
  if (tbls(itim).eq.0.) tbls(itim) = miss_val
endif

if (flag_dugwo) then
tgwo(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      tgwo(itim) = tgwo(itim) + dugwo(ilon,ilat,ilev,itim) &
       * rayon(ilon,ilat,ilev,itim)*coslat(ilat)      &
       * dmass(ilon,ilat,ilev,itim) &
       / hadley
    endif
  enddo
 enddo
enddo
  if (tgwo(itim).eq.0.) tgwo(itim) = miss_val
endif

if (flag_dugwno) then
tgwno(itim) = 0.
do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if (rayon(ilon,ilat,ilev,itim).ne.miss_val) then
      tgwno(itim) = tgwno(itim) + dugwno(ilon,ilat,ilev,itim) &
         * rayon(ilon,ilat,ilev,itim)*coslat(ilat)       &
         * dmass(ilon,ilat,ilev,itim)  &
       / hadley
    endif
  enddo
 enddo
enddo
  if (tgwno(itim).eq.0.) tgwno(itim) = miss_val
endif

enddo ! timelength

!===============================================================================
! 2.2.6 Torques from dynzon
!===============================================================================
if(dzflag) then

allocate(tdyndz(timelength))
allocate(tdisdz(timelength))
allocate(tspgdz(timelength))
allocate(tphydz(timelength))
norm=2./3.*a0*a0*omega

do itim=1,timelength

tdyndz(itim) = 0.
tdisdz(itim) = 0.
tspgdz(itim) = 0.
tphydz(itim) = 0.
do ilon=1,lonlength
 do ilat=2,latlength
  do ilev=1,altlength
      tdyndz(itim) = tdyndz(itim) + & 
     (dmcdyn(ilat-1,ilev,itim)+dmcdyn(ilat,ilev,itim))/(2*lonlength) &
             * norm * dmass(ilon,ilat,ilev,itim) / hadley 
      tdisdz(itim) = tdisdz(itim) + &
     (dmcdis(ilat-1,ilev,itim)+dmcdis(ilat,ilev,itim))/(2*lonlength) &
             * norm * dmass(ilon,ilat,ilev,itim) / hadley 
      tspgdz(itim) = tspgdz(itim) + &
     (dmcspg(ilat-1,ilev,itim)+dmcspg(ilat,ilev,itim))/(2*lonlength) &
             * norm * dmass(ilon,ilat,ilev,itim) / hadley 
      tphydz(itim) = tphydz(itim) + &
     (dmcphy(ilat-1,ilev,itim)+dmcphy(ilat,ilev,itim))/(2*lonlength) &
             * norm * dmass(ilon,ilat,ilev,itim) / hadley 
  enddo
 enddo
enddo
  if (tdyndz(itim).eq.0.) tdyndz(itim) = miss_val
  if (tdisdz(itim).eq.0.) tdisdz(itim) = miss_val
  if (tspgdz(itim).eq.0.) tspgdz(itim) = miss_val
  if (tphydz(itim).eq.0.) tphydz(itim) = miss_val

enddo ! timelength

endif ! dzflag

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

! Check variables to output

do itim=1,timelength
  if (flag_dudyn .and.( tdyn(itim).eq.miss_val)) flag_dudyn =.false.
  if (flag_duajs .and.( tajs(itim).eq.miss_val)) flag_duajs =.false.
  if (flag_duvdf .and.( tbls(itim).eq.miss_val)) flag_duvdf =.false.
  if (flag_dugwo .and.( tgwo(itim).eq.miss_val)) flag_dugwo =.false.
  if (flag_dugwno.and.(tgwno(itim).eq.miss_val)) flag_dugwno=.false.
enddo ! timelength

if(dzflag) then
do itim=1,timelength
  if (tdyndz(itim).eq.miss_val) dzflag=.false.
  if (tdisdz(itim).eq.miss_val) dzflag=.false.
  if (tspgdz(itim).eq.miss_val) dzflag=.false.
  if (tphydz(itim).eq.miss_val) dzflag=.false.
enddo ! timelength
endif ! dzflag


! 1D Variables

datashape1d=time_dimid
 
call write_var1d(outfid,datashape1d,timelength,&
                "oaam      ", "Total rotation ang  ","E25kgm2s-1",miss_val,&
                 oaam )

call write_var1d(outfid,datashape1d,timelength,&
                "raam      ", "Total wind ang      ","E25kgm2s-1",miss_val,&
                 raam )

call write_var1d(outfid,datashape1d,timelength,&
                "tmou      ", "Mountain torque     ","E18kgm2s-2",miss_val,&
                 tmou )

if (flag_dudyn) then
call write_var1d(outfid,datashape1d,timelength,&
                "tdyn      ", "Dynamics torque     ","E18kgm2s-2",miss_val,&
                 tdyn )
endif

if (flag_duajs) then
call write_var1d(outfid,datashape1d,timelength,&
                "tajs      ", "Dynamics torque     ","E18kgm2s-2",miss_val,&
                 tajs )
endif

if (flag_duvdf) then
call write_var1d(outfid,datashape1d,timelength,&
                "tbls      ", "Friction torque     ","E18kgm2s-2",miss_val,&
                 tbls )
endif

if (flag_dugwo) then
call write_var1d(outfid,datashape1d,timelength,&
                "tgwo      ", "Orographic GW torque","E18kgm2s-2",miss_val,&
                 tgwo )
endif

if (flag_dugwno) then
call write_var1d(outfid,datashape1d,timelength,&
                "tgwno     ", "Non-orogr. GW torque","E18kgm2s-2",miss_val,&
                 tgwno )
endif

if(dzflag) then

call write_var1d(outfid,datashape1d,timelength,&
                "tdyndz    ", "Dynamics torque DZ  ","E18kgm2s-2",miss_val,&
                 tdyndz )

call write_var1d(outfid,datashape1d,timelength,&
                "tdisdz    ", "Dissip torque DZ    ","E18kgm2s-2",miss_val,&
                 tdisdz )

call write_var1d(outfid,datashape1d,timelength,&
                "tspgdz    ", "Sponge torque DZ    ","E18kgm2s-2",miss_val,&
                 tspgdz )

call write_var1d(outfid,datashape1d,timelength, &
                "tphydz    ", "Physics torque DZ   ","E18kgm2s-2",miss_val,&
                 tphydz )

endif ! dzflag

! 2D variables

datashape2d(1)=lon_dimid
datashape2d(2)=lat_dimid

call write_var2d(outfid,datashape2d,lonlength,latlength,&
                "phis      ", "Surface geopot      ","m2 s-2    ",miss_val,&
                 phis )

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
                "osam      ", "Specific rotat ang  ","E25m2s-1  ",miss_val,&
                 osam )

call write_var4d(outfid,datashape4d,lonlength,latlength,altlength,timelength,&
                "rsam      ", "Specific wind ang   ","E25m2s-1  ",miss_val,&
                 rsam )

!!!! Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
endif


end program
