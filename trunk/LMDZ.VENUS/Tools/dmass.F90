subroutine cellmass(infid,latlength,lonlength,altlength,timelength,lmdflag, &
                     miss_val,deltalon,deltalat,coslat,plev,ps,grav,rayon, &
                     dmass )

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer infid ! NetCDF input file ID
integer lonlength ! # of grid points along longitude
integer latlength ! # of grid points along latitude
integer altlength ! # of grid point along altitude (of input datasets)
integer timelength ! # of points along time
logical lmdflag ! true=LMD, false=CAM
real :: miss_val ! missing value
real :: deltalat,deltalon ! lat and lon intervals in radians
real,dimension(latlength) :: coslat ! cos of latitudes
real,dimension(altlength) :: plev ! Pressure levels (Pa)
real,dimension(lonlength,latlength,timelength),intent(in) :: ps ! surface pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: grav,rayon
real,dimension(lonlength,latlength,altlength,timelength),intent(out) :: dmass

!local
character (len=64) :: text ! to store some text
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:,:,:,:),allocatable :: press ! GCM atmospheric pressure 
real,dimension(:),allocatable :: pp ! Pressure levels (Pa)
real,dimension(:,:,:,:),allocatable :: deltap ! pressure thickness of each layer (Pa)
integer ierr,ierr1,ierr2 ! NetCDF routines return codes
integer i,j,ilon,ilat,ilev,itim ! for loops
real :: tmpp ! temporary pressure
real :: p0 ! GCM reference pressure
integer tmpvarid

include "planet.h"

!===============================================================================
! 2.2.1 Pressure intervals
!===============================================================================
allocate(deltap(lonlength,latlength,altlength,timelength))
allocate(press(lonlength,latlength,altlength,timelength))
allocate(pp(altlength))

if(lmdflag) then  ! LMD vs CAM
  text="pres"
  call get_var4d(infid,lonlength,latlength,altlength,timelength,text,press,miss_val,ierr1,ierr2)
  if (ierr1.ne.NF_NOERR) then
! assume we are using _P files 
do itim=1,timelength
 do ilon=1,lonlength
  do ilat=1,latlength
     press(ilon,ilat,:,itim) = plev
  enddo
 enddo
enddo   
  else
    if (ierr2.ne.NF_NOERR) stop "Error: Failed reading pres"
  endif

else  ! LMD vs CAM
  ierr=NF_INQ_VARID(infid,"P0",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Failed to get P0 ID, used psref=",psref
    p0=psref
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,p0)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Failed reading reference pressure, used psref=",psref
      p0=psref
    endif
  endif
  ! hybrid coordinate aps/hyam
  text="hyam"
  ierr=NF_INQ_VARID(infid,text,tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not find aps/hyam coordinates"
  else
    allocate(aps(altlength))
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading aps/hyam"
    endif
    call reverselev(altlength,aps)
  endif

  ! hybrid coordinate hybm
  text="hybm"
  ierr=NF_INQ_VARID(infid,text,tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get bps/hybm ID"
  else
    allocate(bps(altlength))
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading bps/hybm"
    endif
    call reverselev(altlength,bps)
  endif

  aps=aps*p0

  do itim=1,timelength
    do ilev=1,altlength
      do ilat=1,latlength
        do ilon=1,lonlength
          press(ilon,ilat,ilev,itim)=aps(ilev)+bps(ilev)*ps(ilon,ilat,itim)
        enddo
      enddo
    enddo
  enddo

endif  ! LMD vs CAM

do itim=1,timelength

do ilon=1,lonlength
 do ilat=1,latlength
  tmpp=ps(ilon,ilat,itim)
! Beware when miss_val is negative...
  do ilev=1,altlength
    if (press(ilon,ilat,ilev,itim).eq.miss_val) then
      pp(ilev)=1.e33
    else
      pp(ilev)=press(ilon,ilat,ilev,itim)
    endif
  enddo
if (tmpp.ge.pp(1)) then
  deltap(ilon,ilat,1,itim)= tmpp - exp((log(pp(1))+log(pp(2)))/2.) ! initialization, rectified with ps below
else
  deltap(ilon,ilat,1,itim)= miss_val
endif
do ilev=2,altlength-1
  deltap(ilon,ilat,ilev,itim)= exp((log(pp(ilev-1))+log(pp(ilev)))/2.)-&
                     exp((log(pp(ilev))+log(pp(ilev+1)))/2.)
enddo
deltap(ilon,ilat,altlength,itim)=exp((log(pp(altlength-1))+log(pp(altlength)))/2.)
 enddo
enddo

do ilon=1,lonlength
 do ilat=1,latlength
  tmpp=ps(ilon,ilat,itim)
! Beware when miss_val is negative...
  do ilev=1,altlength
    if (press(ilon,ilat,ilev,itim).eq.miss_val) then
      pp(ilev)=1.e33
    else
      pp(ilev)=press(ilon,ilat,ilev,itim)
    endif
  enddo
  do ilev=altlength,2,-1
    if ((pp(ilev).le.tmpp).and.(pp(ilev-1).gt.tmpp)) then
      deltap(ilon,ilat,ilev,itim)= tmpp - &
               exp((log(pp(ilev))+log(pp(ilev+1)))/2.)
    elseif (pp(ilev).gt.tmpp) then
      deltap(ilon,ilat,ilev,itim)=miss_val
    endif
  enddo
 enddo
enddo

enddo ! timelength

!===============================================================================
! 2.2.2 Mass in cells
!===============================================================================

do itim=1,timelength

do ilon=1,lonlength
 do ilat=1,latlength
  do ilev=1,altlength
    if ((deltap(ilon,ilat,ilev,itim).ne.miss_val) &
   .and. (rayon(ilon,ilat,ilev,itim).ne.miss_val) &
   .and.  (grav(ilon,ilat,ilev,itim).ne.miss_val)) then
     dmass(ilon,ilat,ilev,itim) = &
               rayon(ilon,ilat,ilev,itim)*coslat(ilat)*deltalon &
             * rayon(ilon,ilat,ilev,itim)*deltalat &
             * deltap(ilon,ilat,ilev,itim)/grav(ilon,ilat,ilev,itim)
    else
     dmass(ilon,ilat,ilev,itim) = miss_val
    endif
  enddo
 enddo
enddo

enddo ! timelength

end subroutine cellmass
