

program hrecast

! This program reads fields from GCM output files
! (diagfi.nc, stats.nc, concat.nc, ...)
! and recasts it on a new horizontal grid specified by the user.
! The output file name is automatically generated from input file name:
! if input file name is 'input.nc', then output file will be 'input_h.nc' 
!
! NB: It is OK if the output grid is not like typical GCM outputs grids,
!     but longitudes range must be something like from -180 to 180 (although
!     you may omit either of these endpoints), and latitudes must range
!     from 90 to -90 (again, you may omit the endpoints).
! EM 09/2009
! TN 01/2013 : Adapted for large output files with at least 2 variables > 2 GiB

implicit none

include "netcdf.inc" ! NetCDF definitions

character(len=128) :: infile ! input file name (diagfi.nc or stats.nc format)
character(len=128) :: outfile ! output file name
character(len=64) :: text ! to store some text
character(len=64) :: tmpvarname ! temporarily store a variable name
integer :: infid ! NetCDF input file ID
integer :: outfid ! NetCDF output file ID
integer :: tmpvarid ! temporarily store a variable ID
integer :: tmpdimid ! temporarily store a dimension ID
integer :: tmpndims ! temporarily store # of dimensions of a variable
integer :: nbvarinfile ! # of variables in input file
integer :: nbattr ! # of attributes of a given variable in input file
integer :: nbvar3dinfile ! # of 3D (lon,lat,time) variables in input file
integer :: nbvar4dinfile ! # of 4D (lon,lat,alt,time) variables in input file

integer :: i,j
integer :: lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs
integer :: lon_varid,lat_varid,alt_varid,time_varid ! NetCDF variable IDs
!integer gcm_layers_dimid ! NetCDF dimension ID for # of layers in GCM
integer :: sigma_varid,aps_varid,bps_varid
integer :: phisinit_varid
integer,dimension(4) :: datashape ! shape of 4D datasets
integer :: ierr ! NetCDF routines return code
character (len=64), dimension(:), allocatable :: var
! var(): names of variables that will be processed
integer nbvar ! # of variables to process
integer,dimension(:),allocatable :: var_id ! IDs of variables var() (in outfile)
real :: miss_val=-9.99e+33 ! special "missing value" to specify missing data
real,parameter :: miss_val_def=-9.99e+33 ! default value for "missing value"
real,dimension(:),allocatable :: inlat ! input latitude
integer :: inlatlength  ! # of elements in input latitude
real,dimension(:),allocatable :: inlon ! input longitude
integer :: inlonlength ! # of elements in input longitude
real,dimension(:),allocatable :: alt ! altitude
integer :: altlength ! # of elements in altitude
real,dimension(:),allocatable :: time ! input time
integer :: timelength ! # of elements in time(:)
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:),allocatable :: sigma ! sigma levels
real,dimension(:,:),allocatable :: inphisinit ! input ground geopotential
real,dimension(:),allocatable :: lon ! output longitude
integer :: lonlength ! # of elements in lon()
real,dimension(:),allocatable :: wklon ! work longitude (with modulo element)
integer :: wklonlength ! # of elements in wklon()
real,dimension(:),allocatable :: lat ! output latitude
integer :: latlength ! # of elements in lat()
real,dimension(:),allocatable :: wklat ! work latitude (includes poles)
integer :: wklatlength ! # of elements in wklat()
!real,dimension(:),allocatable :: lon_bound ! output longitude boundaries
!real,dimension(:),allocatable :: lat_bound ! output latitude boundaries
real,dimension(:),allocatable :: in_lon_bound ! input longitude boundaries
real,dimension(:),allocatable :: in_lat_bound ! input latitude boundaries
real,dimension(:),allocatable :: wk_lon_bound ! work longitude boundaries
real,dimension(:),allocatable :: wk_lat_bound ! work latitude boundaries
real,dimension(:,:),allocatable :: in_2d_data ! input 2D (lon-lat) dataset
real,dimension(:,:),allocatable :: wk_2d_data ! work 2D dataset
real,dimension(:,:),allocatable :: out_2d_data ! output 2D dataset
real,dimension(:,:,:),allocatable :: in_3d_data ! intput 3D dataset
real,dimension(:,:,:),allocatable :: wk_3d_data ! work 3D dataset
real,dimension(:,:,:),allocatable :: out_3d_data ! output 3D dataset
real,dimension(:,:,:,:),allocatable :: in_4d_data ! intput 4D dataset
real,dimension(:,:,:,:),allocatable :: wk_4D_data ! work 4D dataset
real,dimension(:,:,:,:),allocatable :: out_4d_data ! output 4D dataset

real :: pi ! =3.14...
logical :: have_sigma ! Flag: true if sigma levels are known (false if hybrid
                      !       coordinates are used)
logical :: have_geopot ! Flag: true if input file contains ground geopotential
                       ! phisinit()
logical :: out_mod_lon ! Flag: true if output grid has modulo longitude (ie: 
                       ! first and last point are in fact at same longitude)
logical :: out_has_poles ! Flag: true if output grid includes North and South
                         ! poles 
                         
integer, dimension(4) :: edges,corner ! needed to write variables for big files

!===============================================================================
! 1. Input parameters
!===============================================================================
pi=2.*asin(1.)

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
   write(*,*) 'ERROR: Pb opening file ',trim(infile)
   stop ""
endif

!===============================================================================
! 1.2 Get # and names of variables in input file
!===============================================================================

ierr=NF_INQ_NVARS(infid,nbvarinfile)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'ERROR: Failed geting number of variables from file'
  stop
endif

write(*,*)" The following variables have been found:"
nbvar3dinfile=0
nbvar4dinfile=0
do i=1,nbvarinfile
  ! get name of variable # i
  ierr=NF_INQ_VARNAME(infid,i,tmpvarname)
  ! check if it is a 3D variable
  ierr=NF_INQ_VARNDIMS(infid,i,tmpndims)
  if (tmpndims.eq.3) then
    nbvar3dinfile=nbvar3dinfile+1
    write(*,*) trim(tmpvarname)
  endif
  ! check if it is a 4D variable
  ierr=NF_INQ_VARNDIMS(infid,i,tmpndims)
  if (tmpndims.eq.4) then
    nbvar4dinfile=nbvar4dinfile+1
    write(*,*) trim(tmpvarname)
  endif
enddo

allocate(var(nbvar3dinfile+nbvar4dinfile))

write(*,*) ""
write(*,*) "Which variable do you want to keep?"
write(*,*) "all or list of <variables> (separated by <Return>s)"
write(*,*) "(an empty line , i.e: just <Return>, implies end of list)"
nbvar=0
read(*,'(a64)') tmpvarname
do while ((tmpvarname.ne.' ').and.(trim(tmpvarname).ne.'all'))
  ! check if tmpvarname is valid
  ierr=NF_INQ_VARID(infid,tmpvarname,tmpvarid)
  if (ierr.eq.NF_NOERR) then ! valid name
    nbvar=nbvar+1
    var(nbvar)=tmpvarname
  else ! invalid name
    write(*,*) 'Error: ',trim(tmpvarname),' is not a valid name'
    write(*,*) '       (we''ll skip that one)'
  endif
  read(*,'(a64)') tmpvarname
enddo

! handle "all" case
if (tmpvarname.eq.'all') then
  nbvar=0
  do i=1,nbvarinfile
    ! look for 4D variables
    ierr=NF_INQ_VARNDIMS(infid,i,tmpndims)
    if (tmpndims.eq.4) then
      nbvar=nbvar+1
      ! get the corresponding name
      ierr=NF_INQ_VARNAME(infid,i,tmpvarname)
      var(nbvar)=tmpvarname
    endif
  enddo
endif

! Check that there is at least 1 variable to process
if (nbvar.eq.0) then
  write(*,*) 'No variables to process !?'
  write(*,*) 'Might as well stop here'
  stop ""
else
  write(*,*) ""
  write(*,*) 'OK, the following variables will be processed:'
  do i=1,nbvar
    write(*,*) var(i)
  enddo
endif

!===============================================================================
! 1.3 Get input grids in lon,lat,alt,time,
!     as well as hybrid coordinates aps() and bps() (or sigma levels sigma())
!     and eventually phisinit() from input file
!===============================================================================

! 1.3.1 input longitude, latitude, altitude and time

! latitude
ierr=NF_INQ_DIMID(infid,"latitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get latitude dimension ID"
else
  ierr=NF_INQ_VARID(infid,"latitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get latitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,inlatlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get latitude length"
    else
      allocate(inlat(inlatlength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,inlat)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading latitude"
      endif
    endif
  endif
endif

! check that these latitudes are 'GCM' latitudes (i.e. poles are included)
if ((abs(inlat(1)-90.0)).gt.0.001) then
  write(*,*) "Error: Input latitudes should include north pole, but"
  write(*,*) " lat(1)=",inlat(1)
  stop
endif
if ((abs(inlat(inlatlength)+90.0)).gt.0.001) then
  write(*,*) "Error: Input latitudes should include south pole, but"
  write(*,*) " lat(inlatlength)=",inlat(inlatlength)
  stop
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
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,inlonlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get longitude length"
    else
      allocate(inlon(inlonlength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,inlon)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading longitude"
      endif
    endif
  endif
endif

! check that these longitudes are 'GCM' longitudes (i.e range from -180 to 180)
if (abs(180.+inlon(1)).gt.0.001) then
  write(*,*) "Error: Input latitudes should start at -180, but"
  write(*,*) " lon(1)=",inlon(1)
  stop
endif
if (abs(inlon(inlonlength)-180).gt.0.001) then
  write(*,*) "Error: Input latitudes should end at 180, but"
  write(*,*) " lon(inlonlength)=",inlon(inlonlength)
  stop
endif


! altitude
ierr=NF_INQ_DIMID(infid,"altitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get altitude dimension ID"
else
  ierr=NF_INQ_VARID(infid,"altitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get altitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get altitude length"
    else
      allocate(alt(altlength))
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,alt)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading altitude"
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

! 1.3.2 Get hybrid coordinates (or sigma levels)

! start by looking for sigma levels
ierr=NF_INQ_VARID(infid,"sigma",tmpvarid)
if (ierr.ne.NF_NOERR) then
  have_sigma=.false.
  write(*,*) "Could not find sigma levels... will look for hybrid coordinates"
else
  have_sigma=.true.
  allocate(sigma(altlength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,sigma)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading sigma"
  endif
endif

! if no sigma levels, look for hybrid coordinates
if (.not.have_sigma) then
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
endif !of if (.not.have_sigma)

! 1.3.3 Get ground geopotential phisinit, if available

! look for 'phisinit' in current file
ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid) 
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get phisinit ID from file ",trim(infile)
  write(*,*) "  ...will not store geopotential in output... "
  have_geopot=.false.
else
  have_geopot=.true.
  ! Get input physinit
  allocate(inphisinit(inlonlength,inlatlength))
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,inphisinit)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading phisinit"
  endif
endif

! 1.3.4 Create input longitude and latitude boundaries

! build input longitude boundaries (in radians)
allocate(in_lon_bound(inlonlength))
do i=1,inlonlength-1
  in_lon_bound(i)=0.5*(inlon(i+1)+inlon(i))*pi/180.0
enddo
! we have inlon(1)=inlon(inlonlength) modulo 360
! and thus in_lon_bound(inlonlength)=in_lon_bound(1)+360
in_lon_bound(inlonlength)=2.*pi+in_lon_bound(1)
!do i=1,inlonlength
!  write(*,*) "i=",i,"180/pi*in_lon_bound(i)=",(180./pi)*in_lon_bound(i)
!enddo

! build input latitude boundaries (in radians)
allocate(in_lat_bound(inlatlength-1))
do i=1,inlatlength-1
  in_lat_bound(i)=0.5*(inlat(i)+inlat(i+1))*pi/180.0
enddo
!do i=1,inlatlength-1
!  write(*,*) "i=",i,"180/pi*in_lat_bound(i)",(180./pi)*in_lat_bound(i)
!enddo

!===============================================================================
! 1.4 Get output longitude and latitude coordinates
!===============================================================================

write(*,*) ""
write(*,*) "Output horizontal grid:"
write(*,*) "Number of grid points in longitude?"
read(*,*) lonlength
write(*,*) "Enter longitudes (degrees, in [-180:180]), "
write(*,*) " in increasing order (one per line):"
allocate(lon(lonlength))
do i=1,lonlength
  read(*,*) lon(i)
enddo

!! build 'work' longitude (which must be a modulo axis; i.e. first and
!! last points e.g. -180 and 180 are included)
if (abs((lon(1)+360.-lon(lonlength))).le.0.01) then
  ! the axis already has modulo endpoints
  out_mod_lon=.true.
  wklonlength=lonlength
  allocate(wklon(wklonlength))
  wklon(1:lonlength)=lon(1:lonlength)
else
  ! add an extra point
  out_mod_lon=.false.
  wklonlength=lonlength+1
  allocate(wklon(wklonlength))
  wklon(1:lonlength)=lon(1:lonlength)
  wklon(wklonlength)=wklon(1)+360.0
endif

! build work longitude boundaries (in radians)
allocate(wk_lon_bound(wklonlength))
do i=1,lonlength-1
  wk_lon_bound(i)=0.5*(wklon(i+1)+wklon(i))*pi/180.0
enddo
if (out_mod_lon) then
  ! we have lon(1)=lon(lonlength) modulo 360
  ! and thus lon_bound(lonlength)=lon_bound(1)+360
  wk_lon_bound(wklonlength)=2.*pi+wk_lon_bound(1)
else
  wk_lon_bound(wklonlength-1)=0.5*(wklon(wklonlength)+wklon(wklonlength-1))*pi/180.0
  wk_lon_bound(wklonlength)=2.*pi+wk_lon_bound(1)
endif


write(*,*) "Number of grid points in latitude?"
read(*,*) latlength
write(*,*) "Enter latitudes (degrees), in decreasing order (from northernmost"
write(*,*) " to southernmost), one per line:"
allocate(lat(latlength))
do i=1,latlength
  read(*,*) lat(i)
enddo

! build 'work' latitude (which must include poles, ie lat=90 and -90)
if (abs(lat(1)-90.0).le.0.001) then
  out_has_poles=.true.
  wklatlength=latlength
  allocate(wklat(wklatlength))
  wklat(1:latlength)=lat(1:latlength)
else
  out_has_poles=.false.
  ! add poles
  wklatlength=latlength+2
  allocate(wklat(wklatlength))
  wklat(1)=90
  wklat(2:latlength+1)=lat(1:latlength)
  wklat(wklatlength)=-90
endif

! build work latitude boundaries (in radians)
allocate(wk_lat_bound(wklatlength-1))
if (out_has_poles) then
  do i=1,wklatlength-1
    wk_lat_bound(i)=0.5*(wklat(i)+wklat(i+1))*pi/180.0
  enddo
else
  ! put northermost boundary near pole
  wk_lat_bound(1)=(90-0.01*(90.-lat(1)))*pi/180.0
  do i=2,wklatlength-2
    wk_lat_bound(i)=0.5*(wklat(i)+wklat(i+1))*pi/180.0
  enddo
  ! put southernmost boundary near pole
  wk_lat_bound(wklatlength-1)=(-90.0-0.01*(-90.-lat(latlength)))*pi/180.0
endif

!do i=1,wklatlength-1
!  write(*,*) "i=",i,"180/pi*wk_lat_bound(i)",(180./pi)*wk_lat_bound(i)
!enddo

!===============================================================================
! 1.5 Output file
!===============================================================================
write(*,*) ""
outfile=infile(1:len_trim(infile)-3)//"_h.nc"
write(*,*) "Output file name is: ",trim(outfile)


!===============================================================================
! 2. Create output file and initialize definitions of variables and dimensions
!===============================================================================

!===============================================================================
! 2.1. Output file
!===============================================================================

! Create output file
ierr=NF_CREATE(outfile,IOR(NF_CLOBBER,NF_64BIT_OFFSET),outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile
  stop
endif

!===============================================================================
! 2.2. Define dimensions
!===============================================================================
! longitude
ierr=NF_DEF_DIM(outfid,"longitude",lonlength,lon_dimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define longitude dimension"
endif

! latitude
ierr=NF_DEF_DIM(outfid,"latitude",latlength,lat_dimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define latitude dimension"
endif

! altitude
ierr=NF_DEF_DIM(outfid,"altitude",altlength,alt_dimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define altitude dimension"
endif

! time
ierr=NF_DEF_DIM(outfid,"Time",NF_UNLIMITED,time_dimid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define latitude dimension"
endif

!===============================================================================
! 2.3. Define variables and their attributes
!===============================================================================

! 2.3.1 Define 1D variables

! longitude 
datashape(1)=lon_dimid
ierr=NF_DEF_VAR(outfid,"longitude",NF_REAL,1,datashape(1),lon_varid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define longitude variable"
endif

! longitude attributes
text='east longitude'
ierr=NF_PUT_ATT_TEXT(outfid,lon_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing long_name for longitude"
endif
text='degrees_east'
ierr=NF_PUT_ATT_TEXT(outfid,lon_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing units for longitude"
endif

! latitude
datashape(2)=lat_dimid
ierr=NF_DEF_VAR(outfid,"latitude",NF_REAL,1,datashape(2),lat_varid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define latitude variable"
endif

! latitude attributes
text='north latitude'
ierr=NF_PUT_ATT_TEXT(outfid,lat_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing long_name for latitude"
endif
text='degrees_north'
ierr=NF_PUT_ATT_TEXT(outfid,lat_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing units for latitude"
endif

! altitude
datashape(3)=alt_dimid
ierr=NF_DEF_VAR(outfid,"altitude",NF_REAL,1,datashape(3),alt_varid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define altitude variable"
endif

! altitude attributes
! preliminary stuff, get input file altitude variable ID 
ierr=NF_INQ_VARID(infid,"altitude",tmpvarid)

! look for a 'long_name' attribute
text=" "
ierr=NF_GET_ATT_TEXT(infid,tmpvarid,"long_name",text)
if (ierr.eq.NF_NOERR) then
  ! found the attribute; write it to output file 
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
endif

! look for a 'unit' attribute
text=" "
ierr=NF_GET_ATT_TEXT(infid,tmpvarid,"units",text)
if (ierr.eq.NF_NOERR) then
  ! found the attribute; write it to output file 
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
endif

! look for a 'positive' attribute
text=" "
ierr=NF_GET_ATT_TEXT(infid,tmpvarid,"positive",text)
if (ierr.eq.NF_NOERR) then
  ! found the attribute; write it to output file 
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
endif

! sigma levels or hybrid coordinates
if (have_sigma) then
  ierr=NF_DEF_VAR(outfid,"sigma",NF_REAL,1,alt_dimid,sigma_varid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not define sigma variable"
  endif
else ! hybrid coordinates
  ierr=NF_DEF_VAR(outfid,"aps",NF_REAL,1,alt_dimid,aps_varid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not define aps variable"
  endif
  ierr=NF_DEF_VAR(outfid,"bps",NF_REAL,1,alt_dimid,bps_varid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not define bps variable"
  endif
endif

! sigma levels (or hybrid coordinates) attributes
if (have_sigma) then
  text="sigma levels"
  ierr=NF_PUT_ATT_TEXT(outfid,sigma_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for sigma"
  endif
else ! hybrid coordinates
  text="hybrid pressure at midlayers"
  ierr=NF_PUT_ATT_TEXT(outfid,aps_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for aps"
  endif
  text="hybrid sigma at midlayers"
  ierr=NF_PUT_ATT_TEXT(outfid,bps_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for bps"
  endif
endif ! of if (have_sigma)

! time
datashape(4)=time_dimid
ierr=NF_DEF_VAR(outfid,"Time",NF_REAL,1,datashape(4),time_varid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not define Time variable"
endif

! time attributes
! preliminary stuff, get input file time variable ID 
ierr=NF_INQ_VARID(infid,"Time",tmpvarid)
! look for a 'unit' attribute
text=" "
ierr=NF_GET_ATT_TEXT(infid,tmpvarid,"units",text)
if (ierr.eq.NF_NOERR) then
  ! found the attribute; write it to output file 
  ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'units',len_trim(text),text)
else
  ! write something not too stupid
  text='days since 0000-01-1 00:00:00'
  ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for Time"
  endif
endif

! 2.3.2 Define 2D (lon-lat, time independent) variables

! ground geopotential
if (have_geopot) then
  if (have_geopot) then
    ierr=NF_DEF_VAR(outfid,"phisinit",NF_REAL,2,datashape,phisinit_varid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Could not define phisinit variable"
    endif
  endif
  ! ground geopotential attributes
  text='Geopotential at the surface'
  ierr=NF_PUT_ATT_TEXT(outfid,phisinit_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for phisinit"
  endif
endif ! of if (have_geopot)

! 2.3.3 Define 3D&4D variables (ie: surface/volume+time variables)

! variables requested by user
allocate(var_id(nbvar))
do i=1,nbvar
  ! Get the (input file) ID of the variable (stored in tmpvarid)
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to get ID for input variable ',trim(var(i))
    stop ""
  endif
  
  ! Get # of dimensions of the variable
  ! and define the variable in output file
  write(*,*) ""
  write(*,*) "Creating ",trim(var(i))
  ierr=NF_INQ_VARNDIMS(infid,tmpvarid,tmpndims)
  if (tmpndims.eq.3) then
    datashape(3)=time_dimid
    ierr=NF_DEF_VAR(outfid,var(i),NF_REAL,3,datashape,var_id(i))
  else
    if (tmpndims.eq.4) then
      datashape(3)=alt_dimid
      datashape(4)=time_dimid
      ierr=NF_DEF_VAR(outfid,var(i),NF_REAL,4,datashape,var_id(i))
    else
      write(*,*) "Error: number of dimensions of input variable ",trim(var(i))
      write(*,*) "       is ",tmpndims
      stop
    endif
  endif
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, could not define variable ',trim(var(i))
    stop ""
  endif
  
  ! Get the (input file) ID for the variable and
  ! the # of attributes associated to that variable
  ierr=NF_INQ_VARNATTS(infid,tmpvarid,nbattr)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, could not get number of attributes for variable ',&
               trim(var(i))
    stop ""
  endif
  ! inititialize j == number of attributes written to output
  j=0
  
  ! look for a "long_name" attribute (or eventually a 'title' attribute)
  text='   '
  ierr=NF_GET_ATT_TEXT(infid,tmpvarid,'long_name',text)
  if (ierr.ne.NF_NOERR) then ! no long_name attribute
    ! try to find an equivalent 'title' attribute
    text='   '
    ierr=NF_GET_ATT_TEXT(infid,tmpvarid,'title',text)
    if (ierr.eq.NF_NOERR) then ! found 'title' attribute
      write(*,*) "Found title ",trim(text)
      j=j+1
      ! write it as a 'long_name' attribute
      ierr=NF_PUT_ATT_TEXT(outfid,var_id(i),'long_name',len_trim(text),text)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "Error failed to copy title attribute:",trim(text)
      stop ""
      endif
    endif
  else ! found long_name; write it to outfile
    write(*,*) "Found long_name ",trim(text)
    j=j+1
    ierr=NF_PUT_ATT_TEXT(outfid,var_id(i),'long_name',len_trim(text),text)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error failed to copy long_name attribute:",trim(text)
      stop""
    endif
  endif
  
  ! look for a "units" attribute
  text='   '
  ierr=NF_GET_ATT_TEXT(infid,tmpvarid,'units',text)
  if (ierr.eq.NF_NOERR) then ! found 'units' attribute
    write(*,*) "Found units ",trim(text)
    j=j+1
    ! write it to output
    ierr=NF_PUT_ATT_TEXT(outfid,var_id(i),'units',len_trim(text),text)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error failed to copy units attribute:",trim(text)
      stop""
    endif
  endif
  
  ! look for a "missing_value" attribute
  ierr=NF_GET_ATT_REAL(infid,tmpvarid,"missing_value",miss_val)
  if (ierr.eq.NF_NOERR) then ! found 'missing_value' attribute
    write(*,*) "Found missing_value ",miss_val
    j=j+1
  else ! no 'missing_value' attribute, set miss_val to default
    miss_val=miss_val_def
  endif
  
  ! write the missing_value attribute to output
  ierr=NF_PUT_ATT_REAL(outfid,var_id(i),'missing_value',NF_REAL,1,miss_val)
  if (ierr.ne.NF_NOERR) then
    stop "Error, failed to write missing_value attribute"
  endif
  
  ! warn if some attributes were missed
  if (j.ne.nbattr) then
    write(*,*)'Warning, it seems some attributes of variable ',trim(var(i))
    write(*,*)"were not transfered to the new file"
    write(*,*)'nbattr:',nbattr,' j:',j
  endif

enddo ! of do i=1,nbvar

!===============================================================================
! 2.4. Write dimensions (and time-independent variables)
!===============================================================================
! Switch out of NetCDF define mode
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not switch out of define mode"
endif

! Write longitude
ierr=NF_PUT_VAR_REAL(outfid,lon_varid,lon)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not write longitude data to output file"
endif

! Write latitude
ierr=NF_PUT_VAR_REAL(outfid,lat_varid,lat)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not write latitude data to output file"
endif

! write altitude
ierr=NF_PUT_VAR_REAL(outfid,alt_varid,alt)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not write altitude data to output file"
endif

! Write sigma levels (or hybrid coordinates)
if (have_sigma) then
  ierr=NF_PUT_VAR_REAL(outfid,sigma_varid,sigma)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not write sigma data to output file"
  endif
else ! hybrid coordinates
  ierr=NF_PUT_VAR_REAL(outfid,aps_varid,aps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not write aps data to output file"
  endif
  ierr=NF_PUT_VAR_REAL(outfid,bps_varid,bps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not write bps data to output file"
  endif
endif

! write time
ierr=NF_PUT_VARA_REAL(outfid,time_varid,1,timelength,time)
if (ierr.ne.NF_NOERR) then
  stop "Error: Could not write Time data to output file"
endif

!===============================================================================
! 3. Interpolate and write 2D variables
!===============================================================================
write(*,*) "interpolate 2D variables"
allocate(in_2d_data(inlonlength,inlatlength))
allocate(wk_2d_data(wklonlength,wklatlength))
allocate(out_2d_data(lonlength,latlength))

! ground geopotential
if (have_geopot) then
  ! load input geopotential: get ID for input phisinit
  ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get phisinit ID"
  endif
  ! Get physinit
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,in_2d_data)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading input phisinit"
  endif
  
  ! interpolate onto new grid
  call interp_horiz(in_2d_data,wk_2d_data,inlonlength-1,inlatlength-1,&
                    wklonlength-1,wklatlength-1,1,&
                    in_lon_bound,in_lat_bound,wk_lon_bound,wk_lat_bound)
  
  ! copy (and possibly reshape) data to output grid
  if (out_has_poles) then
    out_2d_data(1:lonlength,1:latlength)=wk_2d_data(1:lonlength,1:latlength)
  else
    out_2d_data(1:lonlength,1:latlength)=wk_2d_data(1:lonlength,2:latlength+1)
  endif
  
  ! write interpolated phisinit to output
  ierr=NF_PUT_VAR_REAL(outfid,phisinit_varid,out_2d_data)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Could not write phisinit data to output file"
  endif
endif ! of if (have_geopot)

!===============================================================================
! 4. Interpolate and write 3D (surface+time) variables
!===============================================================================
write(*,*) "interpolate 3D variables"
allocate(in_3d_data(inlonlength,inlatlength,timelength))
allocate(wk_3d_data(wklonlength,wklatlength,timelength))
allocate(out_3d_data(lonlength,latlength,timelength))

do i=1,nbvar ! loop on all selected 3D&4D variables
  ! get input file ID for this variable
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Failed to get input ID for ",trim(var(i))
    stop
  endif
  ! get # of dimensions for this variable
  ierr=NF_INQ_VARNDIMS(infid,tmpvarid,tmpndims)

  if (tmpndims.eq.3) then ! if it is indeed at 3D variable
    ! get the data
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,in_3d_data)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Failed reading input ",trim(var(i))
      stop
    endif
    ! interpolate data
    do j=1,timelength
      call interp_horiz(in_3d_data(1,1,j),wk_3d_data(1,1,j),&
                        inlonlength-1,inlatlength-1,&
                        wklonlength-1,wklatlength-1,1,&
                        in_lon_bound,in_lat_bound,wk_lon_bound,wk_lat_bound)
      ! copy (and possibly reshape) data to output grid
      if (out_has_poles) then
        out_3d_data(1:lonlength,1:latlength,j)=&
                        wk_3d_data(1:lonlength,1:latlength,j)
      else
        out_3d_data(1:lonlength,1:latlength,j)=&
                        wk_3d_data(1:lonlength,2:latlength+1,j)
      endif
    enddo
    ! write interpolated data to output
    corner(:)=1
    edges(1)=lonlength
    edges(2)=latlength
    edges(3)=timelength
    ierr=NF_PUT_VARA_REAL(outfid,var_id(i),corner(1:3),edges(1:3),out_3d_data)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write ",trim(var(i))," data to output file"
      stop
    else
      write(*,*) "wrote ",trim(var(i))
    endif
  endif ! of if (tmpndims.eq.3)

enddo ! of do i=1,nbvar

!===============================================================================
! 5. Interpolate and write 4D variables
!===============================================================================
allocate(in_4d_data(inlonlength,inlatlength,altlength,timelength))
allocate(wk_4d_data(wklonlength,wklatlength,altlength,timelength))
allocate(out_4d_data(lonlength,latlength,altlength,timelength))

do i=1,nbvar ! loop on all selected 3D&4D variables
  ! get input file ID for this variable
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Failed to get input ID for ",trim(var(i))
    stop
  endif
  ! get # of dimensions for this variable
  ierr=NF_INQ_VARNDIMS(infid,tmpvarid,tmpndims)

  if (tmpndims.eq.4) then ! if it is indeed at 4D variable
    ! get the data
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,in_4d_data)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Failed reading input ",trim(var(i))
      stop
    endif
    ! interpolate data
    do j=1,timelength
      call interp_horiz(in_4d_data(1,1,1,j),wk_4d_data(1,1,1,j),&
                        inlonlength-1,inlatlength-1,&
                        wklonlength-1,wklatlength-1,altlength,&
                        in_lon_bound,in_lat_bound,wk_lon_bound,wk_lat_bound)
      ! copy (and possibly reshape) data to output grid
      if (out_has_poles) then
        out_4d_data(1:lonlength,1:latlength,1:altlength,j)=&
                        wk_4d_data(1:lonlength,1:latlength,1:altlength,j)
      else
        out_4d_data(1:lonlength,1:latlength,1:altlength,j)=&
                        wk_4d_data(1:lonlength,2:latlength+1,1:altlength,j)
      endif
    enddo
    ! write interpolated data to output
    corner(:)=1
    edges(1)=lonlength
    edges(2)=latlength
    edges(3)=altlength
    edges(4)=timelength
    ierr=NF_PUT_VARA_REAL(outfid,var_id(i),corner,edges,out_4d_data)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write ",trim(var(i))," data to output file"
      stop
    else
      write(*,*) "wrote ",trim(var(i))
    endif
  endif ! of if (tmpndims.eq.3)

enddo ! of do i=1,nbvar

!===============================================================================
! 6. Close output file
!===============================================================================
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
endif

end program hrecast



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine interp_horiz (varo,varn,imo,jmo,imn,jmn,lm,            &
     &  rlonuo,rlatvo,rlonun,rlatvn)  

!===========================================================
!  Interpolation Horizontales des variables d'une grille LMDZ
! (des points SCALAIRES au point SCALAIRES)
!  dans une autre grille LMDZ en conservant la quantite
!  totale pour les variables intensives (/m2) : ex : Pression au sol
!
! Francois Forget (01/1995)
!===========================================================

      IMPLICIT NONE 

!   Declarations:
! ==============
!
!  ARGUMENTS
!  """""""""
        
       INTEGER imo, jmo ! dimensions ancienne grille (input)
       INTEGER imn,jmn  ! dimensions nouvelle grille (input)

       REAL rlonuo(imo+1)     !  Latitude et
       REAL rlatvo(jmo)       !  longitude des
       REAL rlonun(imn+1)     !  bord des 
       REAL rlatvn(jmn)     !  cases "scalaires" (input)

       INTEGER lm ! dimension verticale (input)
       REAL varo (imo+1, jmo+1,lm) ! var dans l'ancienne grille (input)
       REAL varn (imn+1,jmn+1,lm) ! var dans la nouvelle grille (output)

! Autres variables
! """"""""""""""""
       INTEGER imnmx2,jmnmx2
!       parameter (imnmx2=190,jmnmx2=100)
       parameter (imnmx2=360,jmnmx2=190)
       REAL airetest(imnmx2+1,jmnmx2+1)
       INTEGER ii,jj,l

       REAL airen ((imnmx2+1)*(jmnmx2+1)) ! aire dans la nouvelle grille
       REAL airentotn	! aire totale pole nord dans la nouvelle grille
       REAL airentots	! aire totale pole sud dans la nouvelle grille
!    Info sur les ktotal intersection entre les cases new/old grille

! kmax: le nombre  max  d'intersections entre les 2 grilles horizontales
! On fixe kmax a la taille de la grille des donnees martiennes (360x179) 
! + des pouiemes (cas ou une maille est a cheval sur 2 ou 4 mailles)
!  Il y a un test dans iniinterp_h pour s'assurer que ktotal < kmax
       INTEGER kmax, k, ktotal
       parameter (kmax = 360*179 + 200000)
!      parameter (kmax = 360*179 + 40000)

       INTEGER iik(kmax), jjk(kmax),jk(kmax),ik(kmax)
       REAL intersec(kmax)
       REAL R
       REAL totn, tots
       integer prev_sumdim
       save prev_sumdim
       data prev_sumdim /0/
       

       logical firsttest, aire_ok
       save firsttest
       data firsttest /.true./
       data aire_ok /.true./

       integer imoS,jmoS,imnS,jmnS
       save imoS,jmoS,imnS,jmnS
       save ktotal,iik,jjk,jk,ik,intersec,airen
!       REAL pi

! Test dimensions imnmx2 jmnmx2
!""""""""""""""""""""""""""""""
! test dimensionnement tableau airetest
      if (imn.GT.imnmx2.OR.jmn.GT.jmnmx2) then
         write(*,*) 'STOP pb dimensionnement tableau airetest'
         write(*,*) 'il faut imn < imnmx2 et jmn < jmnmx2'
         write(*,*) 'imn imnmx2', imn,imnmx2
         write(*,*) 'jmn jmnmx2', jmn,jmnmx2
         call exit(1)
      endif

! initialisation
! --------------
! Si c'est le premier appel,  on prepare l'interpolation
! en calculant pour chaque case autour d'un point scalaire de la
! nouvelle grille, la surface  de intersection avec chaque
!    case de l'ancienne grille.

!  This must also be done if we change the dimension
      if (imo+jmo+imn+jmn.ne.prev_sumdim) then
          firsttest=.true.
          prev_sumdim=imo+jmo+imn+jmn
      end if       

      if (firsttest) then 
        call iniinterp_h(imo,jmo,imn,jmn ,kmax,                         &
     &       rlonuo,rlatvo,rlonun,rlatvn,                               &
     &          ktotal,iik,jjk,jk,ik,intersec,airen)
       imoS=imo
       jmoS=jmo
       imnS=imn
       jmnS=jmn
      else
       if(imo.NE.imoS.OR.jmo.NE.jmoS.OR.imn.NE.imnS.OR.jmn.NE.jmnS) then
        call iniinterp_h(imo,jmo,imn,jmn ,kmax,                         &
     &       rlonuo,rlatvo,rlonun,rlatvn,                               &
     &          ktotal,iik,jjk,jk,ik,intersec,airen)
       imoS=imo
       jmoS=jmo
       imnS=imn
       jmnS=jmn
       end if
      end if

      do l=1,lm
       do jj =1 , jmn+1
        do ii=1, imn+1
          varn(ii,jj,l) =0.
        end do
       end do
      end do 
       
! Interpolation horizontale
! -------------------------
! boucle sur toute les ktotal intersections entre les cases
! de l'ancienne et la  nouvelle grille
!
     
      do k=1,ktotal
        do l=1,lm
         varn(iik(k),jjk(k),l) = varn(iik(k),jjk(k),l)                  &
     &   + varo(ik(k), jk(k),l)*intersec(k)/airen(iik(k)                &
     &   +(jjk(k)-1)*(imn+1))
        end do
      end do

! Une seule valeur au pole pour les variables ! :
! -----------------------------------------------
      DO l=1, lm
         totn =0.
         tots =0.


! moyenne du champ au poles (ponderee par les aires)
!"""""""""""""""""""""""""""""""
         airentotn=0.
         airentots=0.

         do ii =1, imn+1
            totn = totn + varn(ii,1,l)*airen(ii)
            tots = tots + varn (ii,jmn+1,l)*airen(jmn*(imn+1)+ii)
            airentotn=airentotn + airen(ii)
            airentots=airentots + airen(jmn*(imn+1)+ii)
         end do 

         do ii =1, imn+1
            varn(ii,1,l) = totn/airentotn
            varn(ii,jmn+1,l) = tots/airentots
         end do 

      ENDDO
           

!---------------------------------------------------------------
!  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST  TEST 
      if (firsttest) then
!      pi=2.*asin(1.)
      firsttest = .false.
!      write (*,*) 'INTERP. HORIZ. : TEST SUR LES AIRES:'

      do jj =1 , jmn+1
        do ii=1, imn+1
          airetest(ii,jj) =0.
        end do
      end do 
      do k=1,ktotal
         airetest(iik(k),jjk(k))= airetest(iik(k),jjk(k)) +intersec(k) 
      end do
      do jj =1 , jmn+1
       do ii=1, imn+1
         r = airen(ii+(jj-1)*(imn+1))/airetest(ii,jj)
         if ((r.gt.1.001).or.(r.lt.0.999)) then
             write (*,*) '********** PROBLEME D'' AIRES !!!',           &
     &                   ' DANS L''INTERPOLATION HORIZONTALE'
             write(*,*)'ii,jj,airen,airetest',                          &
     &          ii,jj,airen(ii+(jj-1)*(imn+1)),airetest(ii,jj)
             aire_ok = .false.
         end if
       end do
      end do
!      if (aire_ok) write(*,*) 'INTERP. HORIZ. : AIRES OK'
      endif

! FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST  FIN TEST
! --------------------------------------------------------------------


        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine iniinterp_h (imo,jmo,imn,jmn ,kmax,                    &
     &       rlonuo,rlatvo,rlonun,rlatvn,                               &
     &       ktotal,iik,jjk,jk,ik,intersec,airen)
   
      implicit none



! ---------------------------------------------------------
! Prepare l' interpolation des variables d'une grille LMDZ
!  dans une autre grille LMDZ en conservant la quantite
!  totale pour les variables intensives (/m2) : ex : Pression au sol
!
!   (Pour chaque case autour d'un point scalaire de la nouvelle
!    grille, on calcule la surface (en m2)en intersection avec chaque
!    case de l'ancienne grille , pour la future interpolation)
!
! on calcule aussi l' aire dans la nouvelle grille 
!
!
!   Auteur:  F.Forget 01/1995
!   -------
!
! ---------------------------------------------------------
!   Declarations:
! ==============
!
!  ARGUMENTS
!  """""""""
! INPUT
       integer imo, jmo ! dimensions ancienne grille
       integer imn,jmn  ! dimensions nouvelle grille
       integer kmax ! taille du tableau des intersections
       real rlonuo(imo+1)     !  Latitude et
       real rlatvo(jmo)       !  longitude des
       real rlonun(imn+1)     !  bord des
       real rlatvn(jmn)     !  cases "scalaires" (input)

! OUTPUT
       integer ktotal ! nombre totale d''intersections reperees
       integer iik(kmax), jjk(kmax),jk(kmax),ik(kmax)
       real intersec(kmax)  ! surface des intersections (m2)
       real airen (imn+1,jmn+1) ! aire dans la nouvelle grille


       
 
! Autres variables
! """"""""""""""""
       integer i,j,ii,jj
       integer imomx,jmomx,imnmx1,jmnmx1
!       parameter (imomx=361,jmomx=180,imnmx1=190,jmnmx1=100)
       parameter (imomx=361,jmomx=180,imnmx1=360,jmnmx1=190)
       real a(imomx+1),b(imomx+1),c(jmomx+1),d(jmomx+1)
       real an(imnmx1+1),bn(imnmx1+1)
       real cn(jmnmx1+1),dn(jmnmx1+1)
       real aa, bb,cc,dd
       real pi

       pi      = 2.*ASIN( 1. )

! Test dimensions imnmx1 jmnmx1
!""""""""""""""""""""""""""""""
! test dimensionnement tableau airetest
      if (imn.GT.imnmx1.OR.jmn.GT.jmnmx1) then
        write(*,*) 'STOP pb dimensionnement'
        write(*,*) 'il faut imn < imnmx1 et jmn < jmnmx1'
        write(*,*) 'imn imnmx1', imn,imnmx1
        write(*,*) 'jmn jmnmx1', jmn,jmnmx1
        call exit(1)
      endif

      if (imo.GT.imomx.OR.jmo.GT.jmomx) then
        write(*,*) 'STOP pb dimensionnement'
        write(*,*) 'il faut imo < imomx et jmo < jmomx'
        write(*,*) 'imo imomx', imo,imomx
        write(*,*) 'jmo jmomx', jmo,jmomx
        call exit(1)
      endif

! On repere les frontieres des cases :
! =================================== 
!
! Attention, on ruse avec des latitudes = 90 deg au pole.


!  Ancienne grile
!  """"""""""""""
      a(1) =   - rlonuo(imo+1)
      b(1) = rlonuo(1)
      do i=2,imo+1
         a(i) = rlonuo(i-1)
         b(i) =  rlonuo(i)
      end do

      d(1) = pi/2 
      do j=1,jmo
         c(j) = rlatvo(j) 
         d(j+1) = rlatvo(j)
      end do
      c(jmo+1) = -pi/2 
      

!  Nouvelle grille
!  """""""""""""""
      an(1) =  - rlonun(imn+1)
      bn(1) = rlonun(1)
      do i=2,imn+1
         an(i) = rlonun(i-1)
         bn(i) =  rlonun(i)
      end do

      dn(1) = pi/2 
      do j=1,jmn
         cn(j) = rlatvn(j)
         dn(j+1) = rlatvn(j)
      end do
      cn(jmn+1) = -pi/2 

! Calcul de la surface des cases scalaires de la nouvelle grille
! ==============================================================
      do ii=1,imn + 1
        do jj = 1,jmn+1
           airen(ii,jj) = (bn(ii)-an(ii))*(sin(dn(jj))-sin(cn(jj)))
        end do
      end do

! Calcul de la surface des intersections
! ======================================

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! test dimenssion kmax (calcul de ktotal)
!"""""""""""""""""""""""""""""""""""""""
! Calcul de ktotal, mais ralonge beaucoup en temps => pour debug
!      write(*,*) 
!      write(*,*) 'DEBUT DU TEST KMAX'
      ktotal = 0
      do jj = 1,jmn+1
       do j=1, jmo+1
          if((cn(jj).lt.d(j)).and.(dn(jj).gt.c(j)))then
              do ii=1,imn + 1
                do i=1, imo +1
                    if (  ((an(ii).lt.b(i)).and.(bn(ii).gt.a(i)))       &
     &        .or. ((an(ii).lt.b(i)-2*pi).and.(bn(ii).gt.a(i)-2*pi)     &
     &             .and.(b(i)-2*pi.lt.-pi) )                            &
     &        .or. ((an(ii).lt.b(i)+2*pi).and.(bn(ii).gt.a(i)+2*pi)     &
     &             .and.(a(i)+2*pi.gt.pi) )                             &
     &                     )then
                      ktotal = ktotal +1
                     end if
                enddo
              enddo
           end if
        enddo
      enddo

      if (kmax.LT.ktotal) then
         write(*,*)
         write(*,*) '******** ATTENTION ********' 
         write(*,*) 'kmax =',kmax 
         write(*,*) 'ktotal =',ktotal 
         write(*,*) 'Changer la valeur de kmax dans interp_horiz.F ' 
         write(*,*) 'avec kmax >= ktotal' 
         write(*,*) 'EXIT dans iniinterp_h'
         call exit(1)
      else
!         write(*,*) 'kmax =',kmax 
!         write(*,*) 'ktotal =',ktotal
      end if
!      write(*,*) 'FIN DU TEST KMAX'
!      write(*,*) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     boucle sur la nouvelle grille
!     """"""""""""""""""""""""""""
      ktotal = 0
      do jj = 1,jmn+1
       do j=1, jmo+1
          if((cn(jj).lt.d(j)).and.(dn(jj).gt.c(j)))then
              do ii=1,imn + 1
                do i=1, imo +1
                    if (  ((an(ii).lt.b(i)).and.(bn(ii).gt.a(i)))       &
     &        .or. ((an(ii).lt.b(i)-2*pi).and.(bn(ii).gt.a(i)-2*pi)     &
     &             .and.(b(i)-2*pi.lt.-pi) )                            &
     &        .or. ((an(ii).lt.b(i)+2*pi).and.(bn(ii).gt.a(i)+2*pi)     &
     &             .and.(a(i)+2*pi.gt.pi) )                             &
     &                     )then
                      ktotal = ktotal +1
                      iik(ktotal) =ii
                      jjk(ktotal) =jj
                      ik(ktotal) =i
                      jk(ktotal) =j

                      dd = min(d(j), dn(jj))
                      cc = cn(jj)
                      if (cc.lt. c(j))cc=c(j)
                      if((an(ii).lt.b(i)-2*pi).and.                     &
     &                  (bn(ii).gt.a(i)-2*pi)) then 
                          bb = min(b(i)-2*pi,bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)-2*pi) aa=a(i)-2*pi
                      else if((an(ii).lt.b(i)+2*pi).and.                &
     &                       (bn(ii).gt.a(i)+2*pi)) then
                          bb = min(b(i)+2*pi,bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)+2*pi) aa=a(i)+2*pi
                      else 
                          bb = min(b(i),bn(ii))
                          aa = an(ii)
                          if (aa.lt.a(i)) aa=a(i)
                      end if
                      intersec(ktotal)=(bb-aa)*(sin(dd)-sin(cc))
                     end if
                end do
               end do
             end if
         end do
       end do       


!     TEST  INFO
!     DO k=1,ktotal 
!      ii = iik(k) 
!      jj = jjk(k)
!      i = ik(k)
!      j = jk(k)
!      if ((ii.eq.10).and.(jj.eq.10).and.(i.eq.10).and.(j.eq.10))then
!      if (jj.eq.2.and.(ii.eq.1))then
!          write(*,*) '**************** jj=',jj,'ii=',ii
!          write(*,*) 'i,j =',i,j
!          write(*,*) 'an,bn,cn,dn', an(ii), bn(ii), cn(jj),dn(jj)
!          write(*,*) 'a,b,c,d', a(i), b(i), c(j),d(j)
!          write(*,*) 'intersec(k)',intersec(k)
!          write(*,*) 'airen(ii,jj)=',airen(ii,jj)
!      end if
!     END DO 

      return
      end
