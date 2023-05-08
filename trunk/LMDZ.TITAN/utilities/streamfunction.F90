module planet_const
! planetary constants (set via init_planet_const routine)
implicit none
real :: a0 ! Mean planetary radius
real :: g0 ! gravity
real :: Rmean ! reduced gaz constant
real :: omega ! rotation rate
end module planet_const

program streamfunction

! ----------------------------------------------------------------------------
! Program to calculate the stream function and the total angular momentum
! from stats and diagfi files
! input : diagfi.nc  / concat.nc / stats.nc kind of files
! author: F. Gonzalez-Galindo
!
! Updates :
! 
! - 02/2016 (M. Turbet) : The code was adapted to work for any planet.
!
! ----------------------------------------------------------------------------

use planet_const
implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (diagfi.nc or stats.nc format)
character (len=128) :: infile2 ! second input file (may be needed for 'aire' and 'cv')
character (len=64) :: text ! to store some text

character (len=50), dimension(:), allocatable :: var
! var(): name(s) of variable(s) that will be concatenated
character (len=100) :: filename
! filename(): output file name
integer :: nid,ierr,miss
! nid: [netcdf] file ID #
! ierr: [netcdf] subroutine returned error code
! miss: [netcdf] subroutine returned error code
integer :: tmpvarid
! varid: temporary store a variable ID #
integer tmpdimid ! temporarily store a dimension ID
integer infid ! NetCDF input file ID (of diagfi.nc or stats.nc format)
integer infid2 ! NetCDF input file which contains 'phisini' dataset (diagfi.nc)
integer nbvarinfile ! # of variables in input file
real, dimension(:), allocatable:: lat,lon,alt,time
! lat(): array, stores latitude coordinates
! lon(): array, stores longitude coordinates
! alt(): array, stores altitude coordinates
! time(): array, stores time coordinates
integer :: ilat,ilon,ialt,it,i
! ilat: index for loops on latitude
! ilon: index for loops on longitude 
! ialt: index for loops on altitude
! it: # index for loops on time
! i: index for loops
integer :: nout,latdimout,londimout,altdimout,timedimout,lonvarout,timevarout
integer :: psivarout,momvarout,uvarout,vvarout,rhovarout,tempvarout
integer :: psvarout, phisinitvarout
! nout: [netcdf] output file ID
! latdimout: [netcdf] output latitude (dimension) ID
! londimout: [netcdf] output longitude (dimension) ID
! altdimout: [netcdf] output altitude (dimension) ID
! timedimout: [netcdf] output time (dimension) ID
! lonvarout: [netcdf] ID of output "Time" variable
integer lonlength ! # of grid points along longitude
integer latlength ! # of grid points along latitude
integer altlength ! # of grid point along altitude (of input datasets)
integer timelength ! # of points along time
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:),allocatable :: sigma ! sigma levels
real,dimension(:),allocatable :: lon_fake ! Fake longitude to be written
real,dimension(:,:),allocatable :: aire ! Area of the box
real,dimension(:,:),allocatable :: cv ! Conversion between natural and covariant v
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
real,dimension(:,:,:),allocatable :: ps ! GCM surface pressure
real,dimension(:,:,:,:),allocatable :: press ! GCM atmospheric pressure 
real,dimension(:,:,:,:),allocatable :: temp ! GCM atmospheric temperature
real,dimension(:,:,:,:),allocatable :: u ! GCM zonal wind
real,dimension(:,:,:,:),allocatable :: v ! GCM meridional wind
real,dimension(:,:,:,:),allocatable :: rho ! GCM atmospheric density
real,dimension(:,:,:),allocatable :: vcont ! Covariant meridional wind
real,dimension(:,:),allocatable :: vcontcum ! Zonal mean covariant meridional wind
real,dimension(:,:,:),allocatable :: plev ! Pressure from hybrid coordinates
real,dimension(:,:,:),allocatable :: mass ! Mass in a grid box
real,dimension(:,:,:),allocatable :: dmass ! Mass variation
real,dimension(:,:),allocatable :: psi ! Stream function
real,dimension(:,:,:),allocatable :: mom ! Angular momentum
real,dimension(:,:),allocatable :: momave ! Zonally averaged angular momentum
real,dimension(:,:,:),allocatable :: ucum ! Temporally averaged zonal wind
real,dimension(:,:,:),allocatable :: vcum ! Temporally averaged meridional wind
real,dimension(:,:,:),allocatable :: rhocum ! Temporally averaged density
real,dimension(:,:,:),allocatable :: tempcum ! Temporally averaged zonal wind
real,dimension(:,:),allocatable :: pscum ! Temporally averaged zonal wind
real,dimension(:,:),allocatable :: uzm ! Zonally averaged zonal wind
real,dimension(:,:),allocatable :: vzm ! Zonally averaged meridional wind
real,dimension(:,:),allocatable :: rhozm ! Zonally averaged density
real,dimension(:,:),allocatable :: tempzm ! Zonally averaged temperature
real,dimension(:),allocatable :: pszm ! Zonally averaged surface pressure
real,dimension(:),allocatable :: phisinitzm ! Zonally averaged ground geopotential

!===============================================================================
! 1.1 Input file
!===============================================================================

write(*,*) ""
write(*,*) " Program valid for diagfi.nc, concatnc.nc and stats.nc files"
write(*,*) "Enter input file name:"

read(*,'(a128)') infile
write(*,*) ""

! open input file

ierr = NF_OPEN(infile,NF_NOWRITE,infid)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening input file'
   stop ""
endif

! load planet constants (radius, gravity, ...)

call init_planet_const(infid)

!===============================================================================
! 1.2 Get grids in lon,lat,alt,time,
!     as well as hybrid coordinates aps() and bps() (or sigma levels sigma())
!     aire() and cv() from input file
!===============================================================================

! 1.2.1 longitude, latitude, altitude and time

! latitude
ierr=NF_INQ_DIMID(infid,"latitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get latitude dimension ID"
else
  ierr=NF_INQ_VARID(infid,"latitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get latitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,latlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get latitude length"
    else
      allocate(lat(latlength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,lat)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading latitude"
      endif
    endif
  endif
endif

! longitude
ierr=NF_INQ_DIMID(infid,"longitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get longitude dimension ID"
else
  ierr=NF_INQ_VARID(infid,"longitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get longitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,lonlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get longitude length"
    else
      allocate(lon(lonlength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,lon)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading longitude"
      endif
    endif
  endif
endif

! time
ierr=NF_INQ_DIMID(infid,"Time",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get Time dimension ID"
else
  ierr=NF_INQ_VARID(infid,"Time",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get Time ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get Time length"
    else
      allocate(time(timelength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,time)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading Time"
      endif
    endif
  endif
endif

! altlength
ierr=NF_INQ_DIMID(infid,"altitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get altitude dimension ID"
else
  ierr=NF_INQ_VARID(infid,"altitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
     stop "Error: Failed to get altitude length"
  else
     ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
     if (ierr.ne.NF_NOERR) then
        stop "Error: Failed to get altitude length"
     else   
        allocate(alt(altlength))
        ierr=NF_GET_VAR_REAL(infid,tmpvarid,alt)
        if (ierr.ne.NF_NOERR) then
           stop "Error: Failed reading Altitude"
        endif
     endif
  endif
endif

! 1.2.2 Get hybrid coordinates 

! look for hybrid coordinates

! hybrid coordinate aps
ierr=NF_INQ_VARID(infid,"aps",tmpvarid)
if (ierr.ne.NF_NOERR) then
   stop "Error: Failed to get aps ID"
else
   allocate(aps(altlength))
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading aps"
   endif
endif

! hybrid coordinate bps
ierr=NF_INQ_VARID(infid,"bps",tmpvarid)
if (ierr.ne.NF_NOERR) then
   stop "Error: Failed to get bps ID"
else
   allocate(bps(altlength))
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading bps"
   endif
endif

!surface pressure
ierr=NF_INQ_VARID(infid,"ps",tmpvarid)
if (ierr.ne.NF_NOERR) then
   stop "Error: Failed to get ps ID"
else
   allocate(ps(lonlength,latlength,timelength))
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,ps)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading ps"
   endif
endif


!===============================================================================
! 1.3 Reads variables in input file
!===============================================================================

ierr=NF_INQ_NVARS(infid,nbvarinfile)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'ERROR: Failed geting number of variables from file'
  stop
endif

!1.3.1 Zonal wind, meridional wind

!Zonal wind
ierr=NF_INQ_VARID(infid,"u",tmpvarid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get u ID"
else
  allocate(u(lonlength,latlength,altlength,timelength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,u)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading zonal wind"
  endif
endif

!Meridional wind
ierr=NF_INQ_VARID(infid,"v",tmpvarid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get v ID"
else
  allocate(v(lonlength,latlength,altlength,timelength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,v)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading zonal wind"
  endif
endif


! 1.3.2 aire

allocate(aire(lonlength,latlength))
! look for 'aire' in current file
ierr=NF_INQ_VARID(infid,"aire",tmpvarid) 
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get aire ID from file ",trim(infile)
  infile2="diagfi.nc"
  write(*,*) "         Trying file ",trim(infile2)
  ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
  if (ierr.ne.NF_NOERR) then
     infile2="diagfi1.nc"
     write(*,*) "         Trying file ",trim(infile2)
     ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
     if (ierr.ne.NF_NOERR) then
        write(*,*) "Problem: Could not find/open these files"
        stop "Might as well stop here"
     endif
  endif

   ! Get ID for aire
   ierr=NF_INQ_VARID(infid2,"aire",tmpvarid)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get aire ID"
   endif
   ! Get aire
   ierr=NF_GET_VAR_REAL(infid2,tmpvarid,aire)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading aire"
   endif
   ! Close file
   write(*,*) 'OK, got aire'
   ierr=NF_CLOSE(infid2)
else
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,aire)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading aire"
   endif
endif


! 1.3.3 phisinit

allocate(phisinit(lonlength,latlength))
! look for '' in current file
ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid) 
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get phisinit ID from file ",trim(infile)
  infile2="diagfi.nc"
  write(*,*) "         Trying file ",trim(infile2)
  ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
  if (ierr.ne.NF_NOERR) then
     infile2="diagfi1.nc"
     write(*,*) "         Trying file ",trim(infile2)
     ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
     if (ierr.ne.NF_NOERR) then
        write(*,*) "Problem: Could not find/open these files"
        stop "Might as well stop here"
     endif
  endif


   ! Get ID for phisinit
   ierr=NF_INQ_VARID(infid2,"phisinit",tmpvarid)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get phisinit ID"
   endif
   ! Get phisinit
   ierr=NF_GET_VAR_REAL(infid2,tmpvarid,phisinit)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading phisinit"
   endif
   ! Close file
   write(*,*) 'OK, got phisinit'
   ierr=NF_CLOSE(infid2)
else
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,phisinit)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading phisinit"
   endif
endif


! 1.3.4 cv

allocate(cv(lonlength,latlength))
! look for 'cv' in current file
ierr=NF_INQ_VARID(infid,"cv",tmpvarid) 
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get cv ID from file ",trim(infile)
  infile2="diagfi.nc"
  write(*,*) "         Trying file ",trim(infile2)
  ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
  if (ierr.ne.NF_NOERR) then
     infile2="diagfi1.nc"
     write(*,*) "         Trying file ",trim(infile2)
     ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
     if (ierr.ne.NF_NOERR) then
        write(*,*) "Problem: Could not find/open these files"
        stop "Might as well stop here"
     endif
  endif


   ! Get ID for cv
   ierr=NF_INQ_VARID(infid2,"cv",tmpvarid)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get cv ID"
   endif
   ! Get cv
   ierr=NF_GET_VAR_REAL(infid2,tmpvarid,cv)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading cv"
   endif
   ! Close file
   write(*,*) 'OK, got cv'
   ierr=NF_CLOSE(infid2)
else
   ierr=NF_GET_VAR_REAL(infid,tmpvarid,cv)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading cv"
   endif
endif


!Other variables: rho, temp

ierr=NF_INQ_VARID(infid,"temp",tmpvarid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get temp ID"
else
  allocate(temp(lonlength,latlength,altlength,timelength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,temp)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading temperature"
  endif
endif

ierr=NF_INQ_VARID(infid,"rho",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Failed to get rho ID"
  write(*,*) "Let's compute it from temperature"
  allocate(rho(lonlength,latlength,altlength,timelength))
  
  do it=1,timelength
    do ilat=1,latlength
      do ilon=1,lonlength
        do ialt=1,altlength
          rho(ilon,ilat,ialt,it)= (aps(ialt)+bps(ialt)*ps(ilon,ilat,it))/(Rmean*temp(ilon,ilat,ialt,it))
        enddo
      enddo
    enddo
  enddo
else
  allocate(rho(lonlength,latlength,altlength,timelength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,rho)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading mass density"
  endif
endif

!==============================================================================
! 1.4. Output file name
!==============================================================================
filename=infile(1:len_trim(infile)-3)//"_stream.nc"
write(*,*)'The output file has name:'
write(*,*) filename


!==============================================================================
! 2.1. Temporal averages
!==============================================================================

allocate(ucum(lonlength,latlength,altlength)) !Time mean zonal wind
allocate(vcum(lonlength,latlength,altlength)) !Time mean meridional wind
allocate(tempcum(lonlength,latlength,altlength)) !Time mean temperature
allocate(rhocum(lonlength,latlength,altlength)) !Time mean density
allocate(pscum(lonlength,latlength)) !Time mean surface pressure
do ilat=1,latlength
   do ilon=1,lonlength
      pscum(ilon,ilat)=0.
      do ialt=1,altlength
            ucum(ilon,ilat,ialt)=0.
            vcum(ilon,ilat,ialt)=0.
            tempcum(ilon,ilat,ialt)=0.
            rhocum(ilon,ilat,ialt)=0.
            do it=1,timelength
               ucum(ilon,ilat,ialt)=ucum(ilon,ilat,ialt)+u(ilon,ilat,ialt,it)
               vcum(ilon,ilat,ialt)=vcum(ilon,ilat,ialt)+v(ilon,ilat,ialt,it)
               tempcum(ilon,ilat,ialt)=tempcum(ilon,ilat,ialt)+temp(ilon,ilat,ialt,it)
               rhocum(ilon,ilat,ialt)=rhocum(ilon,ilat,ialt)+rho(ilon,ilat,ialt,it)
            enddo
            ucum(ilon,ilat,ialt)=ucum(ilon,ilat,ialt)/timelength
            vcum(ilon,ilat,ialt)=vcum(ilon,ilat,ialt)/timelength
            tempcum(ilon,ilat,ialt)=tempcum(ilon,ilat,ialt)/timelength
            rhocum(ilon,ilat,ialt)=rhocum(ilon,ilat,ialt)/timelength
      enddo
      do it=1,timelength
         pscum(ilon,ilat)=pscum(ilon,ilat)+ps(ilon,ilat,it)
      enddo
      pscum(ilon,ilat)=pscum(ilon,ilat)/timelength
   enddo
enddo

!==============================================================================
! 2.2. Calculation of the stream function
!==============================================================================

!Contravariant meridional wind
allocate(vcont(lonlength,latlength,altlength))
do ilon=1,lonlength
   do ialt=1,altlength
      do ilat=1,latlength-1
         vcont(ilon,ilat,ialt)=0.5*&
              (vcum(ilon,ilat,ialt)+vcum(ilon,ilat+1,ialt))/cv(ilon,ilat)
      enddo
      vcont(ilon,latlength,ialt)=0.!vcont(ilon,latlength-1,ialt)
   enddo
enddo

!Pressure from hybrid levels
allocate(Plev(lonlength,latlength,altlength))
do ilon=1,lonlength
   do ilat=1,latlength
      do ialt=1,altlength
         Plev(ilon,ilat,ialt)=aps(ialt)+bps(ialt)*pscum(ilon,ilat)
      enddo
   enddo
enddo

!Mass in the box
allocate(mass(lonlength,latlength,altlength))
do ilon=1,lonlength
   do ilat=1,latlength
      do ialt=1,altlength
         mass(ilon,ilat,ialt)=aire(ilon,ilat)*&
              Plev(ilon,ilat,ialt)/g0
      enddo
   enddo
enddo

!Mass variation in the box
allocate(dmass(lonlength,latlength,altlength))
do ilon=1,lonlength
   do ilat=1,latlength
      do ialt=1,altlength-1
         dmass(ilon,ilat,ialt)=mass(ilon,ilat,ialt)-&
              mass(ilon,ilat,ialt+1)
      enddo
      dmass(ilon,ilat,altlength)=0.!dmass(ilon,ilat,altlength-1)
   enddo
enddo

!Stream function
allocate(psi(latlength,altlength))
allocate(vcontcum(latlength,altlength))
do ilat=1,latlength-1
   psi(ilat,altlength)=0.
   do ialt=altlength-1,1,-1
      psi(ilat,ialt)=psi(ilat,ialt+1)
      vcontcum(ilat,ialt)=0.
      do ilon=1,lonlength
         vcontcum(ilat,ialt)=vcontcum(ilat,ialt)+vcont(ilon,ilat,ialt)/lonlength
         psi(ilat,ialt)=psi(ilat,ialt)+(vcont(ilon,ilat,ialt)&
              *(dmass(ilon,ilat,ialt)+dmass(ilon,ilat+1,ialt)))
      enddo
      psi(latlength,ialt)=0.!psi(latlength-1,ialt)
   enddo
enddo

!==============================================================================
! 2.3. Calculation of the total angular momentum
!==============================================================================

!Angular momentum
allocate(mom(lonlength,latlength,altlength))
allocate(momave(latlength,altlength))
do ilat=1,latlength
   do ialt=1,altlength
      momave(ilat,ialt)=0.
      do ilon=1,lonlength
!            mom(ilon,ilat,ialt,it)=dmass(ilon,ilat,ialt,it)*a0*&
         mom(ilon,ilat,ialt)=a0*cos(lat(ilat)*3.141592/180.)*&
              (omega*a0*cos(lat(ilat)*3.141592/180.)+ucum(ilon,ilat,ialt))
         
         momave(ilat,ialt)=momave(ilat,ialt)+mom(ilon,ilat,ialt)
      enddo
      momave(ilat,ialt)=momave(ilat,ialt)/lonlength
   enddo
enddo

!==============================================================================
! 2.4. Zonal mean winds, temperature and density
!==============================================================================

allocate(uzm(latlength,altlength)) !Zonal mean zonal wind
allocate(vzm(latlength,altlength)) !Zonal mean meridional wind
allocate(tempzm(latlength,altlength)) !Zonal mean temperature
allocate(rhozm(latlength,altlength)) !Zonal mean density
allocate(pszm(latlength)) !Zonal mean surface pressure
allocate(phisinitzm(latlength)) !Zonal mean ground geopotential
do ilat=1,latlength
   phisinitzm(ilat)=0.
   pszm(ilat)=0.
   do ialt=1,altlength
      uzm(ilat,ialt)=0.
      vzm(ilat,ialt)=0.
      tempzm(ilat,ialt)=0.
      rhozm(ilat,ialt)=0.
      do ilon=1,lonlength
         uzm(ilat,ialt)=uzm(ilat,ialt)+ucum(ilon,ilat,ialt)
         vzm(ilat,ialt)=vzm(ilat,ialt)+vcum(ilon,ilat,ialt)
         tempzm(ilat,ialt)=tempzm(ilat,ialt)+tempcum(ilon,ilat,ialt)
         rhozm(ilat,ialt)=rhozm(ilat,ialt)+rhocum(ilon,ilat,ialt)
      enddo
      uzm(ilat,ialt)=uzm(ilat,ialt)/lonlength
      vzm(ilat,ialt)=vzm(ilat,ialt)/lonlength
      tempzm(ilat,ialt)=tempzm(ilat,ialt)/lonlength
      rhozm(ilat,ialt)=rhozm(ilat,ialt)/lonlength
   enddo
   do ilon=1,lonlength
      pszm(ilat)=pszm(ilat)+pscum(ilon,ilat)
      phisinitzm(ilat)=phisinitzm(ilat)+phisinit(ilon,ilat)
   enddo
   pszm(ilat)=pszm(ilat)/lonlength
   phisinitzm(ilat)=phisinitzm(ilat)/lonlength  
enddo


!==============================================================================
! 2.3. Recalculation of angular momentum using zonally averaged wind
!==============================================================================

do ilat=1,latlength
   do ialt=1,altlength
      momave(ilat,ialt)=a0*cos(lat(ilat)*3.141592/180.)*&
           (omega*a0*cos(lat(ilat)*3.141592/180.)+uzm(ilat,ialt))
   enddo
enddo


!==============================================================================
! 3.1 Dimensions in output file
!==============================================================================

! 3.1.1 Initialize output file's lat,lon,alt and time dimensions
call initiate (filename,lat,alt,time,nout,&
     latdimout,londimout,altdimout,timedimout,lonvarout,timevarout)
! Initialize output file's aps,bps variables
call init2(infid,lonlength,latlength,altlength,&
     nout,londimout,latdimout,altdimout)

! 3.1.2 New longitude dimension/value in output file
allocate(lon_fake(1))
lon_fake(1)=0.
#ifdef NC_DOUBLE
        ierr= NF_PUT_VARA_DOUBLE(nout,lonvarout,1,1,lon_fake)
#else
        ierr= NF_PUT_VARA_REAL(nout,lonvarout,1,1,lon_fake)
#endif

! 3.1.3 New time dimension/value in output file
#ifdef NC_DOUBLE
        ierr= NF_PUT_VARA_DOUBLE(nout,timevarout,1,1,lon_fake)
#else
        ierr= NF_PUT_VARA_REAL(nout,timevarout,1,1,lon_fake)
#endif


!==============================================================================
! 3.2 Write variables
!==============================================================================

!Define the dimensions of the variables to be written

!!3.2.1 Stream function
call def_var(nout,"psi","Stream function","",4,&
             (/londimout,latdimout,altdimout,timedimout/),psivarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable psi'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,psivarout,psi)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable psi'
   stop
endif

!3.2.2 Momentum
call def_var(nout,"momave","Angular momentum","",4,&
             (/londimout,latdimout,altdimout,timedimout/),momvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable momave'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,momvarout,momave)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable momave'
   stop
endif


!3.2.3 Zonal mean zonal wind
call def_var(nout,"u","Zonal mean zonal wind","m/s",4,&
             (/londimout,latdimout,altdimout,timedimout/),uvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable u'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,uvarout,uzm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable u'
   stop
endif


!3.2.4 Zonal mean meridional wind
call def_var(nout,"v","Zonal mean meridional wind","m/s",4,&
             (/londimout,latdimout,altdimout,timedimout/),vvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable v'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,vvarout,vzm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable v'
   stop
endif

!3.2.5 Zonal mean density
call def_var(nout,"rho","Zonal mean atmospheric density","",4,&
             (/londimout,latdimout,altdimout,timedimout/),rhovarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable rho'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,rhovarout,rhozm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable rho'
   stop
endif

!3.2.6 Zonal mean temperature
call def_var(nout,"temp","Zonal mean temperature","K",4,&
             (/londimout,latdimout,altdimout,timedimout/),tempvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable temp'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,tempvarout,tempzm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable temp'
   stop
endif

!3.2.7 Zonal mean surface pressure
call def_var(nout,"ps","Zonal mean surface pressure","Pa",3,&
             (/londimout,latdimout,timedimout/),psvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable ps'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,psvarout,pszm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable ps'
   stop
endif


!3.2.8 Zonal mean geopotential
call def_var(nout,"phisinit","Zonal mean initial geopotential","",2,&
             (/londimout,latdimout/),phisinitvarout,ierr)
if (ierr.ne.NF_NOERR) then
   write(*,*) 'Error, could not define variable phisinit'
   stop ""
endif

!Write in the output file
ierr=NF_PUT_VAR_REAL(nout,phisinitvarout,phisinitzm)
if (ierr.ne.NF_NOERR) then
   write(*,*)'Error, Failed to write variable phisinit'
   stop
endif

! Close input file
ierr=nf_close(nid)

! Close output file
ierr=NF_CLOSE(nout)

contains

!******************************************************************************
Subroutine initiate (filename,lat,alt,time,&
                     nout,latdimout,londimout,altdimout,timedimout,&
                     lonvarout,timevarout)
!==============================================================================
! Purpose:
! Create and initialize a data file (NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file (created in this subroutine) remains open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
character (len=*), intent(in):: filename
! filename(): the file's name
real, dimension(:), intent(in):: lat
! lat(): latitude
real, dimension(:), intent(in):: alt
! alt(): altitude
real, dimension(:), intent(in):: time
! time(): Time
integer, intent(out):: nout
! nout: [netcdf] file ID
integer, intent(out):: latdimout
! latdimout: [netcdf] lat() (i.e.: latitude)  ID
integer, intent(out):: londimout
! londimout: [netcdf] lon()  ID
integer, intent(out):: altdimout
! altdimout: [netcdf] alt()  ID
integer, intent(out):: timedimout
! timedimout: [netcdf] "Time"  ID
integer, intent(out):: lonvarout
! timevarout: [netcdf] Longiture (considered as a variable) ID
integer, intent(out):: timevarout
! timevarout: [netcdf] Time (considered as a variable) ID

!==============================================================================
! Local variables:
!==============================================================================
!integer :: latdim,londim,altdim,timedim
integer :: nvarid,ierr
! nvarid: [netcdf] ID of a variable
! ierr: [netcdf]  return error code (from called subroutines)

!==============================================================================
! 1. Create (and open) output file
!==============================================================================
write(*,*) "creating "//trim(adjustl(filename))//'...'
ierr = NF_CREATE(filename,NF_CLOBBER,nout)
! NB: setting NF_CLOBBER mode means that it's OK to overwrite an existing file
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'ERROR: Impossible to create the file.'
   stop ""
endif

!==============================================================================
! 2. Define/write "dimensions" and get their IDs
!==============================================================================

ierr = NF_DEF_DIM(nout, "latitude", size(lat), latdimout)
!ierr = NF_DEF_DIM(nout, "longitude", NF_UNLIMITED, londimout)
ierr = NF_DEF_DIM(nout, "longitude", 1, londimout)
ierr = NF_DEF_DIM(nout, "altitude", size(alt), altdimout)
ierr = NF_DEF_DIM(nout, "Time", 1, timedimout)

! End netcdf define mode
ierr = NF_ENDDEF(nout)

!==============================================================================
! 3. Write "Time" (attributes)
!==============================================================================

call def_var(nout,"Time","Time","years since 0000-00-0 00:00:00",1,&
             (/timedimout/),timevarout,ierr)

!==============================================================================
! 4. Write "latitude" (data and attributes)
!==============================================================================

call def_var(nout,"latitude","latitude","degrees_north",1,&
             (/latdimout/),nvarid,ierr)

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,lat)
#else
ierr = NF_PUT_VAR_REAL (nout,nvarid,lat)
#endif

!==============================================================================
! 4. Write "longitude" (attributes)
!==============================================================================

call def_var(nout,"longitude","East longitude","degrees_east",1,&
             (/londimout/),lonvarout,ierr)


!==============================================================================
! 4. Write "altitude" (data and attributes)
!==============================================================================

! Switch to netcdf define mode

call def_var(nout,"altitude","Altitude","km",1,&
             (/altdimout/),nvarid,ierr)

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,alt)
#else
ierr = NF_PUT_VAR_REAL (nout,nvarid,alt)
#endif 

end Subroutine initiate
!******************************************************************************
subroutine init2(infid,lonlen,latlen,altlen, &
                 outfid,londimout,latdimout,altdimout)
!==============================================================================
! Purpose:
! Copy aps() and bps() from input file to outpout file
!==============================================================================
! Remarks:
! The NetCDF files must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: infid  ! NetCDF output file ID
integer, intent(in) :: lonlen ! # of grid points along longitude
integer, intent(in) :: latlen ! # of grid points along latitude
integer, intent(in) :: altlen ! # of grid points along latitude
integer, intent(in) :: outfid ! NetCDF output file ID
integer, intent(in) :: londimout ! longitude dimension ID
integer, intent(in) :: latdimout ! latitude dimension ID
integer, intent(in) :: altdimout ! altitude dimension ID
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
integer :: apsid,bpsid
integer :: ierr
integer :: tmpvarid ! temporary variable ID
logical :: aps_ok, bps_ok ! are "phisinit" "aps" "bps" available ?


!==============================================================================
! 1. Read data from input file
!==============================================================================

! hybrid coordinate aps
  allocate(aps(altlen))
ierr=NF_INQ_VARID(infid,"aps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "oops Failed to get aps ID. OK"
  aps_ok=.false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
  if (ierr.ne.NF_NOERR) then
   stop "error: Failed reading aps"
  endif
  aps_ok=.true.
endif

! hybrid coordinate bps
  allocate(bps(altlen))
ierr=NF_INQ_VARID(infid,"bps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "oops: Failed to get bps ID. OK"
  bps_ok=.false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading bps"
  endif
  bps_ok=.true.
endif

!==============================================================================
! 2. Write
!==============================================================================

!==============================================================================
! 2.2. Hybrid coordinates aps() and bps()
!==============================================================================

IF (aps_ok) then 
   ! define aps
   call def_var(nout,"aps","hybrid pressure at midlayers"," ",1,&
        (/altdimout/),apsid,ierr)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to def_var aps"
   endif

   ! write aps
#ifdef NC_DOUBLE
   ierr=NF_PUT_VAR_DOUBLE(outfid,apsid,aps)
#else
   ierr=NF_PUT_VAR_REAL(outfid,apsid,aps)
#endif
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to write aps"
   endif
ENDIF 

IF (bps_ok) then 
   ! define bps
   call def_var(nout,"bps","hybrid sigma at midlayers"," ",1,&
        (/altdimout/),bpsid,ierr)
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to def_var bps"
   endif

   ! write bps
#ifdef NC_DOUBLE
   ierr=NF_PUT_VAR_DOUBLE(outfid,bpsid,bps)
#else
   ierr=NF_PUT_VAR_REAL(outfid,bpsid,bps)
#endif
   if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to write bps"
   endif
ENDIF



! Cleanup
deallocate(aps)
deallocate(bps)

end subroutine init2
!******************************************************************************
subroutine def_var(nid,name,title,units,nbdim,dim,nvarid,ierr)
!==============================================================================
! Purpose: Write a variable (i.e: add a variable to a dataset)
! called "name"; along with its attributes "title", "units"...
! to a file (following the NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nid
! nid: [netcdf] file ID #
character (len=*), intent(in) :: name
! name(): [netcdf] variable's name
character (len=*), intent(in) :: title
! title(): [netcdf] variable's "title" attribute
character (len=*), intent(in) :: units
! unit(): [netcdf] variable's "units" attribute
integer, intent(in) :: nbdim
! nbdim: number of dimensions of the variable
integer, dimension(nbdim), intent(in) :: dim
! dim(nbdim): [netcdf] dimension(s) ID(s)
integer, intent(out) :: nvarid
! nvarid: [netcdf] ID # of the variable
integer, intent(out) :: ierr
! ierr: [netcdf] subroutines returned error code

! Switch to netcdf define mode
ierr=NF_REDEF(nid)

! Insert the definition of the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,adjustl(name),NF_DOUBLE,nbdim,dim,nvarid)
#else
ierr=NF_DEF_VAR(nid,adjustl(name),NF_FLOAT,nbdim,dim,nvarid)
#endif

! Write the attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",len_trim(adjustl(title)),adjustl(title))
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",len_trim(adjustl(units)),adjustl(units))

! End netcdf define mode
ierr=NF_ENDDEF(nid)

end subroutine def_var

end program streamfunction

subroutine init_planet_const(infid)
! initialize planetary constants using the "controle" array in the file
! if "controle" array not found in file, look for it in "diagfi.nc"
use planet_const
use netcdf
implicit none

integer,intent(in) :: infid ! input file ID

! local variables
character(len=100) :: varname
integer :: varid ! variable ID
integer :: status
integer :: infid2 ! for an alternate input file
real :: controle(100) ! to store "controle" array

! look for "controle" in input file
infid2=infid ! initialization
varname="controle" 
status=nf90_inq_varid(infid,trim(varname),varid)
if (status.ne.nf90_noerr) then
  write(*,*) "init_planet_const: Failed to find ",trim(varname)
  write(*,*) " looking for it in file diagfi.nc"
  status=nf90_open("diagfi.nc",NF90_NOWRITE,infid2)
  if (status.ne.nf90_noerr) then
    write(*,*) " no diafi.nc file ... looking for diagfi1.nc"
    status=nf90_open("diagfi1.nc",NF90_NOWRITE,infid2)
    if (status.ne.nf90_noerr) then
      write(*,*) "might as well stop here..."
      stop
    endif
  endif
  status=nf90_inq_varid(infid2,trim(varname),varid)
  if (status.ne.nf90_noerr) then
    write(*,*) " Failed to find ",trim(varname)," in file!!"
    stop
  endif
  write(*,*) "OK, found ",trim(varname)
endif
status=nf90_get_var(infid2,varid,controle)
if (status.ne.nf90_noerr) then
  write(*,*) "init_planet_const: Failed to load ",trim(varname)
  stop
endif

! now that we have "controle"; extract relevent informations
a0=controle(5) ! = radius of planet (m)
g0=controle(7) ! = gravity (m.s-2) at a0
Rmean=1000.*8.3144598/controle(8) ! controle(8)=mugaz = molar mass (g.mol-1) of atm.
omega=controle(6) ! rotation rate (rad.s-1)

end subroutine init_planet_const
