










subroutine writediagsoil(ngrid,name,title,units,dimpx,px)

! Write variable 'name' to NetCDF file 'diagsoil.nc'.
! The variable may be 3D (lon,lat,depth) subterranean field,
! a 2D (lon,lat) surface field, or a simple scalar (0D variable).
!
! Calls to 'writediagsoil' can originate from anywhere in the program;
! An initialisation of variable 'name' is done if it is the first time
! that this routine is called with given 'name'; otherwise data is appended
! (yielding the sought time series of the variable)

! Modifs: Aug.2010 Ehouarn: enforce outputs to be real*4

use comsoil_h, only: nsoilmx, inertiedat
use geometry_mod, only: cell_area
use time_phylmdz_mod, only: ecritphy, day_step, iphysiq
use mod_phys_lmdz_para, only : is_mpi_root, is_master, gather
use mod_grid_phy_lmdz, only : klon_glo, Grid1Dto2D_glo, &
                              nbp_lon, nbp_lat

implicit none

include"netcdf.inc"

! Arguments:
integer,intent(in) :: ngrid ! number of (horizontal) points of physics grid
! i.e. ngrid = 2+(jjm-1)*iim - 1/jjm
character(len=*),intent(in) :: name ! 'name' of the variable
character(len=*),intent(in) :: title ! 'long_name' attribute of the variable
character(len=*),intent(in) :: units ! 'units' attribute of the variable
integer,intent(in) :: dimpx ! dimension of the variable (3,2 or 0)
real,dimension(ngrid,nsoilmx),intent(in) :: px ! variable

! Local variables:
real*4,dimension(nbp_lon+1,nbp_lat,nsoilmx) :: data3 ! to store 3D data
real*4,dimension(nbp_lon+1,nbp_lat) :: data2 ! to store 2D data
real*4 :: data0 ! to store 0D data
real*4 :: data3_1d(1,nsoilmx) ! to store a profile in 1D model
real*4 :: data2_1d ! to store surface value with 1D model
integer :: i,j,l ! for loops
integer :: ig0

real*4,save :: date ! time counter (in elapsed days)

real :: inertia((nbp_lon+1),nbp_lat,nsoilmx)
real :: area((nbp_lon+1),nbp_lat)

real :: inertiafi_glo(klon_glo,nsoilmx)
real :: areafi_glo(klon_glo)

integer,save :: isample ! sample rate at which data is to be written to output
integer,save :: ntime=0 ! counter to internally store time steps
character(len=20),save :: firstname="1234567890"
integer,save :: zitau=0
!$OMP THREADPRIVATE(date,isample,ntime,firstname,zitau)

character(len=30) :: filename="diagsoil.nc"

! NetCDF stuff:
integer :: nid ! NetCDF output file ID
integer :: varid ! NetCDF ID of a variable
integer :: ierr ! NetCDF routines return code
integer,dimension(4) :: id ! NetCDF IDs of the dimensions of the variable
integer,dimension(4) :: edges,corners


! 1. Initialization step
if (firstname.eq."1234567890") then
  ! Store 'name' as 'firstname'
  firstname=name
  ! From now on, if 'name'.eq.'firstname', then it is a new time cycle

  ! just to be sure, check that firstnom is large enough to hold nom
  if (len_trim(firstname).lt.len_trim(name)) then
    write(*,*) "writediagsoil: Error !!!"
    write(*,*) "   firstname string not long enough!!"
    write(*,*) "   increase its size to at least ",len_trim(name)
    stop
  endif
  
  ! Set output sample rate
  isample=int(ecritphy) ! same as for diagfi outputs
  ! Note ecritphy is known from control.h
  
  ! Create output NetCDF file
  if (is_master) then
   ierr=NF_CREATE(filename,NF_CLOBBER,nid)
   if (ierr.ne.NF_NOERR) then
    write(*,*)'writediagsoil: Error, failed creating file '//trim(filename)
    stop
   endif

         inertiafi_glo(:,:)=inertiedat(:,:)
         areafi_glo(:)=cell_area(:)

   ! build inertia() and area()
   if (klon_glo>1) then
    do i=1,nbp_lon+1 ! poles
     inertia(i,1,1:nsoilmx)=inertiafi_glo(1,1:nsoilmx)
     inertia(i,nbp_lat,1:nsoilmx)=inertiafi_glo(klon_glo,1:nsoilmx)
     ! for area, divide at the poles by nbp_lon
     area(i,1)=areafi_glo(1)/nbp_lon
     area(i,nbp_lat)=areafi_glo(klon_glo)/nbp_lon
    enddo
    do j=2,nbp_lat-1
     ig0= 1+(j-2)*nbp_lon
     do i=1,nbp_lon
        inertia(i,j,1:nsoilmx)=inertiafi_glo(ig0+i,1:nsoilmx)
        area(i,j)=areafi_glo(ig0+i)
     enddo
     ! handle redundant point in longitude
     inertia(nbp_lon+1,j,1:nsoilmx)=inertia(1,j,1:nsoilmx)
     area(nbp_lon+1,j)=area(1,j)
    enddo
   endif
   
   ! write "header" of file (longitudes, latitudes, geopotential, ...)
   if (klon_glo>1) then ! general 3D case
     call iniwritesoil(nid,ngrid,inertia,area,nbp_lon+1,nbp_lat)
   else ! 1D model
     call iniwritesoil(nid,ngrid,inertiafi_glo(1,:),areafi_glo(1),1,1)
   endif

  endif ! of if (is_master)
  
  ! set zitau to -1 to be compatible with zitau incrementation step below
  zitau=-1
  
else
  ! If not an initialization call, simply open the NetCDF file
  if (is_master) then
   ierr=NF_OPEN(filename,NF_WRITE,nid)
  endif
endif ! of if (firstname.eq."1234567890")

! 2. Increment local time counter, if necessary
if (name.eq.firstname) then
  ! if we run across 'firstname', then it is a new time step
  zitau=zitau+iphysiq
  ! Note iphysiq is known from control.h
endif

! 3. Write data, if the time index matches the sample rate
if (mod(zitau+1,isample).eq.0) then

! 3.1 If first call at this date, update 'time' variable
  if (name.eq.firstname) then
    ntime=ntime+1
    date=float(zitau+1)/float(day_step)
    ! Note: day_step is known from control.h
    
    if (is_master) then
     ! Get NetCDF ID for "time"
     ierr=NF_INQ_VARID(nid,"time",varid)
     ! Add the current value of date to the "time" array
!#ifdef 1
!     ierr=NF_PUT_VARA_DOUBLE(nid,varid,ntime,1,date)
!#else
     ierr=NF_PUT_VARA_REAL(nid,varid,ntime,1,date)
!#endif
     if (ierr.ne.NF_NOERR) then
      write(*,*)"writediagsoil: Failed writing date to time variable"
      stop 
     endif
    endif ! of if (is_master)
  endif ! of if (name.eq.firstname)

! 3.2 Write the variable to the NetCDF file
if (dimpx.eq.3) then ! Case of a 3D variable
  ! A. Recast data along 'dynamics' grid
  if (klon_glo>1) then ! General case
   do l=1,nsoilmx
    ! handle the poles
    do i=1,nbp_lon+1
      data3(i,1,l)=px(1,l)
      data3(i,nbp_lat,l)=px(ngrid,l)
    enddo
    ! rest of the grid
    do j=2,nbp_lat-1
      ig0=1+(j-2)*nbp_lon
      do i=1,nbp_lon
        data3(i,j,l)=px(ig0+i,l)
      enddo
      data3(nbp_lon+1,j,l)=data3(1,j,l) ! extra (modulo) longitude
    enddo
   enddo
  else ! 1D model case
   data3_1d(1,1:nsoilmx)=px(1,1:nsoilmx)
  endif
  
  ! B. Write (append) the variable to the NetCDF file
  if (is_master) then
  ! B.1. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"longitude",id(1))
    ierr=NF_INQ_DIMID(nid,"latitude",id(2))
    ierr=NF_INQ_DIMID(nid,"depth",id(3))
    ierr=NF_INQ_DIMID(nid,"time",id(4))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagsoil: creating variable "//trim(name)
    call def_var(nid,name,title,units,4,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)
  
  ! B.2. Prepare things to be able to write/append the variable
  corners(1)=1
  corners(2)=1
  corners(3)=1
  corners(4)=ntime
  
  if (klon_glo==1) then
    edges(1)=1
  else
    edges(1)=nbp_lon+1
  endif
  edges(2)=nbp_lat
  edges(3)=nsoilmx
  edges(4)=1
  
  ! B.3. Write the slab of data
!#ifdef 1
!  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,data3)
!#else
  if (klon_glo>1) then
    ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,data3)
  else
    ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,data3_1d)
  endif
!#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagsoil: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif
  endif ! of if (is_master)

elseif (dimpx.eq.2) then ! Case of a 2D variable

  ! A. Recast data along 'dynamics' grid
  if (klon_glo>1) then ! general case
    ! handle the poles
    do i=1,nbp_lon+1
      data2(i,1)=px(1,1)
      data2(i,nbp_lat)=px(ngrid,1)
    enddo
    ! rest of the grid
    do j=2,nbp_lat-1
      ig0=1+(j-2)*nbp_lon
      do i=1,nbp_lon
        data2(i,j)=px(ig0+i,1)
      enddo
      data2(nbp_lon+1,j)=data2(1,j) ! extra (modulo) longitude
    enddo
  else ! 1D model case
    data2_1d=px(1,1)
  endif

  ! B. Write (append) the variable to the NetCDF file
  if (is_master) then
  ! B.1. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"longitude",id(1))
    ierr=NF_INQ_DIMID(nid,"latitude",id(2))
    ierr=NF_INQ_DIMID(nid,"time",id(3))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagsoil: creating variable "//trim(name)
    call def_var(nid,name,title,units,3,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)

  ! B.2. Prepare things to be able to write/append the variable
  corners(1)=1
  corners(2)=1
  corners(3)=ntime
  
  if (klon_glo==1) then
    edges(1)=1
  else
    edges(1)=nbp_lon+1
  endif
  edges(2)=nbp_lat
  edges(3)=1
  
  ! B.3. Write the slab of data
!#ifdef 1
!  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,data2)
!#else
  if (klon_glo>1) then ! General case
    ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,data2)
  else
    ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,data2_1d)
  endif
!#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagsoil: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif
  endif ! of if (is_master)

elseif (dimpx.eq.0) then ! Case of a 0D variable
  ! A. Copy data value
  data0=px(1,1)

  ! B. Write (append) the variable to the NetCDF file
  ! B.1. Get the ID of the variable
  ierr=NF_INQ_VARID(nid,name,varid)
  if (ierr.ne.NF_NOERR) then
    ! If we failed geting the variable's ID, we assume it is because
    ! the variable doesn't exist yet and must be created.
    ! Start by obtaining corresponding dimensions IDs
    ierr=NF_INQ_DIMID(nid,"time",id(1))
    ! Tell the world about it
    write(*,*) "====================="
    write(*,*) "writediagsoil: creating variable "//trim(name)
    call def_var(nid,name,title,units,1,id,varid,ierr)
  endif ! of if (ierr.ne.NF_NOERR)

  ! B.2. Prepare things to be able to write/append the variable
  corners(1)=ntime
  
  edges(1)=1

  ! B.3. Write the data
!#ifdef 1
!  ierr=NF_PUT_VARA_DOUBLE(nid,varid,corners,edges,data0)
!#else
  ierr=NF_PUT_VARA_REAL(nid,varid,corners,edges,data0)
!#endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) "writediagsoil: Error: Failed writing "//trim(name)//&
               " to file "//trim(filename)//" at time",date
  endif

endif ! of if (dimpx.eq.3) elseif (dimpx.eq.2) ...
endif ! of if (mod(zitau+1,isample).eq.0)

! 4. Close the NetCDF file
if (is_master) then
  ierr=NF_CLOSE(nid)
endif

end subroutine writediagsoil
