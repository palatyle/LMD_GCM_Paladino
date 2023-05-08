PROGRAM	extrapol_icefield

!--------------------------------------------------------------------------------------------------
! This program is a tool to accelerate the calculation of ice fieds evolution.
!
! It uses data files (diagfi.nc) to extrapolate surface 
! physical fields (ice fields typically) in time.
!
! 1. We load data file(s) 'diagfi.nc' and get dimensions (longitude,latitude,altitude,time).
! 2. We get a surface field from the loaded 'diagfi.nc' file.
! 3. We make the extrapolation.
! 4. We load a start file 'startfi.nc' and copy it into a new start file 'startfi_extrapolated.nc'.
! 5. We modify the 'startfi_extrapolated.nc' according to the extrapolation calculations.
!
! --> 'startfi_extrapolated.nc' is the interpolated new start file that can be used to run
!
!
! Author : M. Turbet (2016) [Adapted from E. Millour previous work]
!
!
!--------------------------------------------------------------------------------------------------

use netcdf
 
implicit none

integer :: j,ig,t,flag

! NetCDF variables
integer :: status_d,inid,varid_d
integer :: lat_dimid,lon_dimid,alt_dimid,time_dimid
integer :: status_s,fileid,varid_s
integer :: physical_points_dimid

character(len=128) :: diag_file=""
character(len=128) :: start_file_in=""
character(len=128) :: start_file_out="startfi_extrapolated.nc"
character(len=128) :: varname_d=""
character(len=128) :: varname_s=""
character(len=128) :: varname_area="area"
 
real,dimension(:),allocatable :: longitude,latitude,altitude,time ! # dimensions
integer lonlen,latlen,altlen,timelen ! # of grid points along longitude/latitude/altitude/time

integer :: physical_points ! number of atmospheric columns (created)
integer :: physical_points_s ! number of atmospheric columns (extracted from startfi file)

real,dimension(:,:,:),allocatable :: field_3D ! LON x LAT x TIME field
real,dimension(:,:),allocatable :: field_phys_3D ! Physical_points x Time field
real,dimension(:),allocatable :: field_phys_3D_reduced ! Physical_points field
real,allocatable :: surf_field_3D(:) ! Physical_points field
real,allocatable :: area(:) ! to store the 1D field of the area

real :: field_ini
real :: field_fin

integer :: nb_field
integer :: datashape(3)

integer :: n_years

integer :: extrapolation_mode

write(*,*) ""
write(*,*) "Welcome in extrapol_icefield routine !"
write(*,*) "This tool was designed to extrapolate ice field evolution."
write(*,*) ""


!===============================================================================
!
! I. INPUT FILES
!
!===============================================================================

write(*,*) ""
write(*,*) "Enter nectdf (diagfi.nc file type) data file name:"
read(*,'(a128)') diag_file
write(*,*) ""
 
write(*,*) ""
write(*,*) "Enter nectdf (startfi.nc file type) input file name:"
read(*,'(a128)') start_file_in
write(*,*) ""

write(*,*) "Output file name is now ", start_file_out
write(*,*) ""

! We copy the startfi (input) file and work now with the startfi (output) file.
CALL system ("cp "//trim(adjustl(start_file_in))//" "//trim(adjustl(start_file_out)))

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! I.1. Open diagfi.nc file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


status_d=nf90_open(diag_file,NF90_NOWRITE,inid)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to open datafile ",trim(diag_file)
  write(*,*)trim(nf90_strerror(status_d))
 stop
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! I.2. Get grids for dimensions longitude,latitude,altitude and time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Get latitudes.
status_d=nf90_inq_dimid(inid,"latitude",lat_dimid)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find latitude dimension"
  write(*,*)trim(nf90_strerror(status_d))
 stop
endif
status_d=nf90_inquire_dimension(inid,lat_dimid,len=latlen)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find latitude length"
  write(*,*)trim(nf90_strerror(status_d))
endif
allocate(latitude(latlen))
status_d=nf90_inq_varid(inid,"latitude",varid_d)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to find latitude ID"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif
status_d=nf90_get_var(inid,varid_d,latitude)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to load latitude"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif

! Get longitudes.
status_d=nf90_inq_dimid(inid,"longitude",lon_dimid)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find longitude dimension"
  write(*,*)trim(nf90_strerror(status_d))
 stop
endif
status_d=nf90_inquire_dimension(inid,lon_dimid,len=lonlen)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find longitude length"
  write(*,*)trim(nf90_strerror(status_d))
endif
allocate(longitude(lonlen))
status_d=nf90_inq_varid(inid,"longitude",varid_d)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to find longitude ID"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif
status_d=nf90_get_var(inid,varid_d,longitude)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to load longitude"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif

! Get times.
status_d=nf90_inq_dimid(inid,"Time",time_dimid)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find Time dimension"
  write(*,*)trim(nf90_strerror(status_d))
 stop
endif
status_d=nf90_inquire_dimension(inid,time_dimid,len=timelen)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find Time length"
  write(*,*)trim(nf90_strerror(status_d))
endif
allocate(time(timelen))
status_d=nf90_inq_varid(inid,"Time",varid_d)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to find Time ID"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif
status_d=nf90_get_var(inid,varid_d,time)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to load Time"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif

! Get altitudes.
status_d=nf90_inq_dimid(inid,"altitude",alt_dimid)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find altitude dimension"
  write(*,*)trim(nf90_strerror(status_d))
 stop
endif
status_d=nf90_inquire_dimension(inid,alt_dimid,len=altlen)
if (status_d.ne.nf90_noerr) then
  write(*,*)"Failed to find altitude length"
  write(*,*)trim(nf90_strerror(status_d))
endif
allocate(altitude(altlen))
status_d=nf90_inq_varid(inid,"altitude",varid_d)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to find altitude ID"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif
status_d=nf90_get_var(inid,varid_d,altitude)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to load altitude"
  write(*,*)trim(nf90_strerror(status_d))
  stop
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! I.3. Get the data field.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 write(*,*) ""
 write(*,*) "Enter surface tracer field (in diagfi.nc) to extrapolate"
 read(*,'(a128)') varname_d
 write(*,*) ""
 
 write(*,*) ""
 write(*,*) "Enter corresponding surface tracer field (in startfi.nc) to extrapolate"
 read(*,'(a128)') varname_s
 write(*,*) ""


status_d=nf90_inq_varid(inid,trim(varname_d),varid_d)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to find variable ",trim(varname_d)," in ",trim(diag_file)
  write(*,*)trim(nf90_strerror(status_d))
  write(*,*) " Might as well stop here."
  stop
endif

! Sanity checks on the variable.
status_d=nf90_inquire_variable(inid,varid_d,ndims=nb_field,dimids=datashape)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to obtain information on variable ",trim(varname_d)
  write(*,*)trim(nf90_strerror(status_d))
  write(*,*) " Might as well stop here."
  stop
else
  ! check that it is a 3D variable.
  if (nb_field.ne.3) then
    write(*,*) "Error, expected a 3D (lon-lat-time) variable!"
    stop
  endif
  ! check that its dimensions are indeed lon,lat,time.
  if (datashape(1).ne.lon_dimid) then
    write(*,*) "Error, expected first dimension to be longitude!"
    write(*,*) "Please check your diagfi.nc file"
    stop
  endif
  if (datashape(2).ne.lat_dimid) then
    write(*,*) "Error, expected second dimension to be latitude!"
    write(*,*) "Please check your diagfi.nc file"
    stop
  endif
  if (datashape(3).ne.time_dimid) then
    write(*,*) "Error, expected third dimension to be time!"
    write(*,*) "Please check your diagfi.nc file"
    stop
  endif
endif

allocate(field_3D(lonlen,latlen,timelen))

! Load the data field.
status_d=nf90_get_var(inid,varid_d,field_3D)
if (status_d.ne.nf90_noerr) then
  write(*,*) "Failed to load ",trim(varname_d)
  write(*,*)trim(nf90_strerror(status_d))
  stop
else
  write(*,*) "Loaded ",trim(varname_d)
endif


!===============================================================================
!
! II. MANIPULATION OF DATA FIELDS
!
!===============================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! II.1. Transfer the 2D data onto the "physics" 1-dimensional grid.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the physical 1D grid.
physical_points=2+(latlen-2)*(lonlen-1)
write(*,*) " lonlen=",lonlen," latlen=",latlen
write(*,*) " physical_points=",physical_points

allocate(field_phys_3D(physical_points,timelen))
allocate(field_phys_3D_reduced(physical_points))

! North pole :
field_phys_3D(1,1:timelen)=field_3D(1,1,1:timelen)
! South pole :
field_phys_3D(physical_points,1:timelen)=field_3D(lonlen,latlen,1:timelen)
! Rest of the world :
do j=2,latlen-1
  ig=2+(j-2)*(lonlen-1)
  field_phys_3D(ig:ig+(lonlen-2),1:timelen)=field_3D(1:lonlen-1,j,1:timelen)
enddo

field_phys_3D_reduced(1:physical_points)=0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! II.2. Extrapolation of data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 write(*,*) ""
 write(*,*) "Extrapolation/Speed Factor ?"
 read(*,*) n_years
 write(*,*) ""
 
 write(*,*) ""
 write(*,*) "Extrapolation Modes :"
 write(*,*) "Mode 1 : We extrapolate using only the first value and the last value."
 write(*,*) "Mode 2 : We extrapolate using the first value and the mean value."
 write(*,*) ""
 write(*,*) "Write the number of the mode you want :"
 read(*,*) extrapolation_mode
 write(*,*) ""


if (extrapolation_mode .eq. 1) then

   do ig = 1, physical_points 
      flag = 0
      do t = 1,timelen
         if (field_phys_3D(ig,t) .eq. 0) flag = 1 ! Here, we remove the "seasonal" ice deposit.
      enddo
      if (flag .eq. 0) then
         field_phys_3D_reduced(ig)=(field_phys_3D(ig,timelen)-field_phys_3D(ig,1))*n_years
         if(field_phys_3D_reduced(ig) .lt. 0) field_phys_3D_reduced(ig) = 0   
      else
         field_phys_3D_reduced(ig)=0. ! Here, we remove the "seasonal" ice deposit.
      endif
   enddo
   
else if (extrapolation_mode .eq. 2) then

   do ig = 1, physical_points 
      flag = 0
      do t = 1,timelen
         if (field_phys_3D(ig,t) .eq. 0) flag = 1 ! Here, we remove the "seasonal" ice deposit.
      enddo
      if (flag .eq. 0) then
         do t = 1,timelen
            field_phys_3D_reduced(ig)=field_phys_3D_reduced(ig)+field_phys_3D(ig,t)
         enddo
         field_phys_3D_reduced(ig)=field_phys_3D_reduced(ig)/timelen
         field_phys_3D_reduced(ig)=(field_phys_3D_reduced(ig)-field_phys_3D(ig,1))*2*n_years
         if(field_phys_3D_reduced(ig) .lt. 0) field_phys_3D_reduced(ig) = 0   
      else
         field_phys_3D_reduced(ig)=0.
      endif
   enddo
else
    write(*,*) "Please check your extrapolation mode because "
    write(*,*) "it does not match the possibilities !"
    stop
endif


!===============================================================================
!
! III. CREATION OF THE EXTRAPOLATED STARTFI FILE
!
!===============================================================================


! Open output file in read/write mode.
status_s=nf90_open(start_file_out,NF90_WRITE,fileid)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to open output file ",trim(start_file_out)
  write(*,*)trim(nf90_strerror(status_s))
 stop
endif
write(*,*) "Reading output file: ",trim(start_file_out)

! Get value of physical_points.
status_s=nf90_inq_dimid(fileid,"physical_points",physical_points_dimid)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to find physical_points dimension"
  write(*,*)trim(nf90_strerror(status_s))
 stop
endif

status_s=nf90_inquire_dimension(fileid,physical_points_dimid,len=physical_points_s)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to find physical_points value"
  write(*,*)trim(nf90_strerror(status_s))
 stop
else
  write(*,*) " physical_points = ",physical_points_s
endif

if (physical_points_s .ne. physical_points) then
   write(*,*) "the number of physical points in diagfi and startfi are not equal!"
   stop
endif

! Load surface field
allocate(surf_field_3D(physical_points))
allocate(area(physical_points))

! Get id of the area
status_s=nf90_inq_varid(fileid,varname_area,varid_s)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to find ",trim(varname_area)," varid_s"
  write(*,*)"Check that the area is included in your files"
  write(*,*)trim(nf90_strerror(status_s))
  stop
endif

! Read variable from file
status_s=nf90_get_var(fileid,varid_s,area)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to load ",trim(varname_area)
  write(*,*)trim(nf90_strerror(status_s))
  stop
else
  write(*,*)"Loaded ",trim(varname_s)
endif

! Get id of variable
status_s=nf90_inq_varid(fileid,varname_s,varid_s)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to find ",trim(varname_s)," varid_s"
  write(*,*)trim(nf90_strerror(status_s))
  stop
endif

! Read variable from file
status_s=nf90_get_var(fileid,varid_s,surf_field_3D)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to load ",trim(varname_s)
  write(*,*)trim(nf90_strerror(status_s))
  stop
else
  write(*,*)"Loaded ",trim(varname_s)
endif


! Here, we calculate the final ice field.
! In particular, we do some calculations to conserve the total amount of ice.

field_ini=0.
do ig=1,physical_points
   field_ini = field_ini + area(ig)*surf_field_3D(ig)
enddo

surf_field_3D(1:physical_points)= surf_field_3D(1:physical_points) &
   + field_phys_3D_reduced(1:physical_points)

field_fin=0.
do ig=1,physical_points
   field_fin = field_fin + area(ig)*surf_field_3D(ig)
enddo

do ig=1, physical_points
   surf_field_3D(ig) = surf_field_3D(ig)*field_ini/field_fin
enddo


! 3. Write the data to the startfi file
status_s=nf90_put_var(fileid,varid_s,surf_field_3D)
if (status_s.ne.nf90_noerr) then
  write(*,*)"Failed to write ",trim(varname_s)
  write(*,*)trim(nf90_strerror(status_s))
  stop
else
  write(*,*)"Wrote ",trim(varname_s)
endif

! 4. Close file
status_s=nf90_close(fileid)

END PROGRAM extrapol_icefield
