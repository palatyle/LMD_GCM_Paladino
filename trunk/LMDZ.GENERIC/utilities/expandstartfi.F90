program expandstartfi

! program which takes a GCM startfi.nc file and puts it on the
! corresponding lonxlat grid (so it can be displayed using Grads, Ferret, etc.)
! usage:
! expandstartfi  [infile.nc] [outfile.nc]
!     if outfile is not specified it is built as "infile_ex.nc"
!     if infile is not specified, "startfi.nc" is used as default

use netcdf
implicit none

! Input and output files:
character(len=256) :: infile="startfi.nc"      ! default input file
character(len=256) :: outfile="startfi_ex.nc"  ! default output file

! NetCDF stuff
integer :: status ! NetCDF return code
integer :: inid,outid ! NetCDF file IDs
integer :: varid ! to store the ID of a variable
integer :: datashape(4) ! to store dimension IDs of a given dataset
integer :: corner(3),edges(3) ! to read data with a time axis
character(len=90) :: varname ! name of a variable
character(len=90) :: varatt ! name of attribute of a variable
character(len=90) :: varattcontent ! content of the attribute
! Input file dimension IDs:
integer :: physical_points_dimid
integer :: subsurface_layers_dimid
integer :: nlayer_plus_1_dimid
integer :: number_of_advected_fields_dimid
integer :: time_dimid
integer :: nbindims ! number of dimensions in input file
integer :: nbinvars ! number of variables in input file
integer :: invarid ! to store ID of an input variable
!Input file variables
integer :: physical_points
integer :: subsurface_layers
integer :: nlayer_plus_1
integer :: number_of_advected_fields
integer :: timelen
real,allocatable :: surf_field(:) ! to store a 1D field of physical_points elements
real,allocatable :: subsurf_field(:,:) ! to store subsurface (2D field)

! Output file dimensions:
integer :: latlen
integer :: lonlen
! Output file variables:
real,allocatable :: longitude(:)
real,allocatable :: latitude(:)
real,allocatable :: depth(:)
real,allocatable :: out_surf_field(:,:)
real,allocatable :: out_subsurf_field(:,:,:)
! Output file dimension IDs
integer :: lon_dimid
integer :: lat_dimid
integer :: depth_dimid
! IDs of Output file variables
integer :: lon_varid
integer :: lat_varid
integer :: depth_varid

integer :: i,j,k,ig0,ivar
integer :: nbdim,nbatt,shape(4)
integer :: nbarg ! # of program arguments
character(len=256) :: arg ! to store a program argument
real :: pi ! 3.14159...

pi=2.*asin(1.)

! 0. Preliminary stuff, check arguments (input and output files)
! get number of arguments this program was called with
nbarg=command_argument_count()

if (nbarg.ge.1) then
  call get_command_argument(1,arg) ! get argument # 1
  infile=trim(arg)
  if (nbarg.eq.2) then
    call get_command_argument(2,arg) ! get argument # 2
    outfile=trim(arg)
  else
   ! build outfile from infile (replace ".nc" with "_ex.nc"
   outfile=trim(infile(1:index(infile,".nc",back=.true.)-1))//"_ex.nc"
  endif
  if (nbarg.ge.3) then
    write(*,*) ' Warning: Too many arguments...'
    write(*,*) '         will only use the first 2 '
  endif
endif

! 1. open input file
status=nf90_open(infile,NF90_NOWRITE,inid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to open datafile ",trim(infile)
  write(*,*)trim(nf90_strerror(status))
 stop
endif
write(*,*) "Reading input file: ",trim(infile)

! 1.2 Identify/load dimensions in input file
status=nf90_inq_dimid(inid,"physical_points",physical_points_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find physical_points dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,physical_points_dimid,len=physical_points)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find physical_points value"
  write(*,*)trim(nf90_strerror(status))
 stop
else
  write(*,*) " physical_points = ",physical_points
endif

status=nf90_inq_dimid(inid,"subsurface_layers",subsurface_layers_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find subsurface_layers dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,subsurface_layers_dimid,len=subsurface_layers)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find subsurface_layers value"
  write(*,*)trim(nf90_strerror(status))
 stop
else
  write(*,*) " subsurface_layers = ",subsurface_layers
endif

status=nf90_inq_dimid(inid,"nlayer_plus_1",nlayer_plus_1_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find nlayer_plus_1 dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,nlayer_plus_1_dimid,len=nlayer_plus_1)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find nlayer_plus_1 value"
  write(*,*)trim(nf90_strerror(status))
 stop
else
  write(*,*) " nlayer_plus_1 = ",nlayer_plus_1
endif

status=nf90_inq_dimid(inid,"number_of_advected_fields",number_of_advected_fields_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find number_of_advected_fields dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,number_of_advected_fields_dimid,len=number_of_advected_fields)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find number_of_advected_fields value"
  write(*,*)trim(nf90_strerror(status))
 stop
else
  write(*,*) " number_of_advected_fields = ",number_of_advected_fields
endif

status=nf90_inq_dimid(inid,"Time",time_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find Time dimension"
  write(*,*)trim(nf90_strerror(status))
  timelen = 0
else
  status=nf90_inquire_dimension(inid,time_dimid,len=timelen)
  if (status.ne.nf90_noerr) then
    write(*,*)"Failed to read Time dimension"
    write(*,*)trim(nf90_strerror(status))
   stop
  else
    write(*,*) " time length = ",timelen
  endif
endif

! 1.3 Allocate memory for input fields
allocate(surf_field(physical_points))
allocate(subsurf_field(physical_points,subsurface_layers))

! 2.1. Create output file
status=nf90_create(outfile,NF90_CLOBBER,outid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to create output file: ",trim(outfile)
  write(*,*)trim(nf90_strerror(status))
  stop
endif
write(*,*) " "
write(*,*) "Created output file: ",trim(outfile)

! 2.2. Build dimensions for output file

! 2.2.1 Compute lonlen and latlen from  physical_points
! load "longitude(physical_points)" from input file
status=nf90_inq_varid(inid,"longitude",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find longitude ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read longitude
status=nf90_get_var(inid,varid,surf_field)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load longitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! count the number of points before longitude(i) wraps around
i=3
lonlen=1
!write(*,*) "longitude(2)=",surf_field(2)
do while (surf_field(i).ne.surf_field(2))
!write(*,*) "i:",i,"surf_field(i)=",surf_field(i)
 i=i+1
 lonlen=lonlen+1
enddo
! and add 1 because we want a redundant lon=180 point
lonlen=lonlen+1
write(*,*) " => lonlen=",lonlen

allocate(longitude(lonlen))
! fill longitude(1:lonlen)
longitude(1:lonlen-1)=surf_field(2:lonlen)
longitude(lonlen)=-longitude(1) ! redundant +Pi/2
! convert to degrees
longitude(:)=longitude(:)*(180./pi)

! knowing lonlen, compute latlen
latlen=2+(physical_points-2)/(lonlen-1)
write(*,*) " => latlen=",latlen

allocate(latitude(latlen))
! get field "latitude(physical_points)" from infile
status=nf90_inq_varid(inid,"latitude",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find latitude ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read latitude
status=nf90_get_var(inid,varid,surf_field)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load latitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

latitude(1)=surf_field(1)
!write(*,*) "latitude(1)=",latitude(1)
do i=2,latlen-1
  latitude(i)=surf_field((i-1)*(lonlen-1))
!  write(*,*) "i:",i,"latitude(i)=",latitude(i)
enddo
latitude(latlen)=surf_field(physical_points)
!write(*,*) "latitude(latlen)=",latitude(latlen)
! convert to degrees
latitude(:)=latitude(:)*(180./pi)

! depth:
allocate(depth(subsurface_layers))
! load "soildepth(subsurface_layers)" from input file
status=nf90_inq_varid(inid,"soildepth",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find soildepth ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read longitude
status=nf90_get_var(inid,varid,depth)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load soildepth"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! 2.2.2 define dimensions to output file
! longitude:
status=nf90_def_dim(outid,"longitude",lonlen,lon_dimid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating longitude dimension"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! latitude:
status=nf90_def_dim(outid,"latitude",latlen,lat_dimid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating longitude dimension"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! depth:
status=nf90_def_dim(outid,"depth",subsurface_layers,depth_dimid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating depth dimension"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

!2.3. write variables to output file
! 2.3.1. define 1D variables
! longitude
datashape(1)=lon_dimid
status=nf90_def_var(outid,"longitude",NF90_REAL,lon_dimid,lon_varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating longitude variable"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! longitude attributes
status=nf90_put_att(outid,lon_varid,"long_name","East longitude")
status=nf90_put_att(outid,lon_varid,"units","degrees_east")

!latitude
datashape(2)=lat_dimid
status=nf90_def_var(outid,"latitude",NF90_REAL,lat_dimid,lat_varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating latitude variable"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! latitude attributes
status=nf90_put_att(outid,lat_varid,"long_name","North latitude")
status=nf90_put_att(outid,lat_varid,"units","degrees_north")

! depth
status=nf90_def_var(outid,"depth",NF90_REAL,depth_dimid,depth_varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed creating depth variable"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
!depth atributes
status=nf90_put_att(outid,depth_varid,"long_name","Soil mid-layer depth")
status=nf90_put_att(outid,depth_varid,"units","m")
status=nf90_put_att(outid,depth_varid,"positive","down")

! 2.3.2 write 1D variable

! swich out of NetCDF define mode
status=nf90_enddef(outid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to swich out of define mode"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! write longitude
status=nf90_put_var(outid,lon_varid,longitude)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed writing longitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! write latitude
status=nf90_put_var(outid,lat_varid,latitude)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed writing latitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! write depth
status=nf90_put_var(outid,depth_varid,depth)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed writing depth"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

! 2.3. 2D (surface) variables
! First find out how many variables there are in input file
status=nf90_inquire(inid,nbindims,nbinvars)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed obtaining nbindims and nbinvars from input file"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

allocate(out_surf_field(lonlen,latlen))

shape(:) = 0
do ivar=1,nbinvars ! loop on all input variables
  ! find out what dimensions are linked to this variable
  status=nf90_inquire_variable(inid,ivar,name=varname,ndims=nbdim,&
                               dimids=shape,natts=nbatt)
  if (((nbdim==1).and.(shape(1)==physical_points_dimid))&
  .or.((nbdim==2).and.(shape(1)==physical_points_dimid)&
                 .and.(shape(2)==time_dimid))) then
  
    corner(1) = 1
    corner(2) = timelen
    edges(1)  = physical_points
    edges(2)  = 1
    
    ! skip "longitude" and "latitude"
    if (trim(varname)=="longitude") cycle
    if (trim(varname)=="latitude") cycle
    
    write(*,*) " processing: ",trim(varname)

    ! load input data:
    status=nf90_inq_varid(inid,varname,invarid)
    status=nf90_get_var(inid,invarid,surf_field,corner,edges)
    
    ! switch output file to to define mode
    status=nf90_redef(outid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to swich to define mode"
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
    !define the variable
    status=nf90_def_var(outid,trim(varname),NF90_REAL,&
                         (/lon_dimid,lat_dimid/),varid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed creating variable ",trim(varname)
      write(*,*)trim(nf90_strerror(status))
      stop
    endif

    ! variable attributes
    do k=1,nbatt
      status=nf90_inq_attname(inid,invarid,k,varatt)
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed getting attribute number",k," for ",trim(varname)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif
      status=nf90_get_att(inid,invarid,trim(varatt),varattcontent)
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed loading attribute ",trim(varatt)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif

      ! write the attribute and its value to output
      status=nf90_put_att(outid,varid,trim(varatt),trim(varattcontent))
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed writing attribute ",trim(varatt)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif
    enddo ! do k=1,nbatt

    ! swich out of NetCDF define mode
    status=nf90_enddef(outid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to swich out of define mode"
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
    
    ! "convert" from surf_field(:) to out_surf_field(:,:)
    do i=1,lonlen
      out_surf_field(i,1)=surf_field(1) ! North Pole
      out_surf_field(i,latlen)=surf_field(physical_points) ! South Pole
    enddo
    do j=2,latlen-1
      ig0=1+(j-2)*(lonlen-1)
      do i=1,lonlen-1
        out_surf_field(i,j)=surf_field(ig0+i)
      enddo
      out_surf_field(lonlen,j)=out_surf_field(1,j) ! redundant lon=180
    enddo
    
    ! write the variable
    status=nf90_put_var(outid,varid,out_surf_field)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to write ",trim(varname)
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
  endif ! of if ((nbdim==1).and.(shape(1)==physical_points_dimid))
enddo ! of do i=1,nbinvars 

! 2.4. 3D (subsurface) variables
allocate(out_subsurf_field(lonlen,latlen,subsurface_layers))

do ivar=1,nbinvars ! loop on all input variables
  ! find out what dimensions are linked to this variable
  status=nf90_inquire_variable(inid,ivar,name=varname,ndims=nbdim,&
                               dimids=shape,natts=nbatt)
  if (((nbdim==2).and.(shape(1)==physical_points_dimid)&
                 .and.(shape(2)==subsurface_layers_dimid))&
  .or.((nbdim==3).and.(shape(1)==physical_points_dimid)&
                 .and.(shape(2)==subsurface_layers_dimid)&
                 .and.(shape(3)==time_dimid))) then
    
    corner(1) = 1
    corner(2) = 1
    corner(3) = timelen
    edges(1)  = physical_points
    edges(2)  = subsurface_layers
    edges(3)  = 1
    
    write(*,*) " processing: ",trim(varname)

    ! load input data:
    status=nf90_inq_varid(inid,varname,invarid)
    status=nf90_get_var(inid,invarid,subsurf_field,corner,edges)
    
    ! switch output file to to define mode
    status=nf90_redef(outid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to swich to define mode"
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
    !define the variable
    status=nf90_def_var(outid,trim(varname),NF90_REAL,&
                         (/lon_dimid,lat_dimid,depth_dimid/),varid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed creating variable ",trim(varname)
      write(*,*)trim(nf90_strerror(status))
      stop
    endif

    ! variable attributes
    do k=1,nbatt
      status=nf90_inq_attname(inid,invarid,k,varatt)
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed getting attribute number",k," for ",trim(varname)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif
      status=nf90_get_att(inid,invarid,trim(varatt),varattcontent)
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed loading attribute ",trim(varatt)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif

      ! write the attribute and its value to output
      status=nf90_put_att(outid,varid,trim(varatt),trim(varattcontent))
      if (status.ne.nf90_noerr) then
        write(*,*) "Failed writing attribute ",trim(varatt)
        write(*,*)trim(nf90_strerror(status))
      stop
      endif
    enddo ! do k=1,nbatt

    ! swich out of NetCDF define mode
    status=nf90_enddef(outid)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to swich out of define mode"
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
    
    ! "convert" from subsurf_field(:,:) to out_subsurf_field(:,:,:)
    do i=1,lonlen
      out_subsurf_field(i,1,:)=subsurf_field(1,:) ! North Pole
      out_subsurf_field(i,latlen,:)=subsurf_field(physical_points,:) ! South Pole
    enddo
    do j=2,latlen-1
      ig0=1+(j-2)*(lonlen-1)
      do i=1,lonlen-1
        out_subsurf_field(i,j,:)=subsurf_field(ig0+i,:)
      enddo
      out_subsurf_field(lonlen,j,:)=out_subsurf_field(1,j,:) ! redundant lon=180
    enddo
    
    ! write the variable
    status=nf90_put_var(outid,varid,out_subsurf_field)
    if (status.ne.nf90_noerr) then
      write(*,*) "Failed to write ",trim(varname)
      write(*,*)trim(nf90_strerror(status))
      stop
    endif
  endif ! of if ((nbdim==1).and.(shape(1)==physical_points_dimid)...)
enddo ! of do i=1,nbinvars 


! 3. Finish things and cleanup
! Close output file
status=nf90_close(outid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to close file: ",trim(outfile)
  write(*,*)trim(nf90_strerror(status))
  stop
endif

end program expandstartfi
