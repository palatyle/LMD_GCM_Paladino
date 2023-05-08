program extract

! program to extract (ie: interpolates) pointwise values of an atmospheric 
! variable from a 'zrecast'ed diagfi file (works if altitude is geometrical 
! height or a pressure vertical coordinates)
! user has to specify:
! - name of input file
! - date (in sols) offset wrt the input file (e.g. if the input file "begins"
!   at Ls=0, then the offset is 0; if the input file begins at Ls=30, the
!   offset date corresponding to the first 3 months is 61+66+66=193 sols, etc.)
! - the "extraction mode": 
!      1: extract individual values; user will specify values of
!         lon lat alt Ls LT (all on a same line)
!         on as many lines as there are sought values
!      2: extract a profile: user will specify on a first line the values of
!         lon lat Ls LT (all on a same line)
!         and then only specify values of altitudes (m or Pa depending on the
!         coordinate in the input file), one per line, at which values are
!         sought
! - output values are sent to (ASCII) output file 'infile_var_.dat', where
!   'infile' is the input file name (without trailing '.nc') and
!   'var' is the sought variable, for extraction mode 1 as
!   lines of "lon lat alt Ls LT value" and for a profile (extraction mode 2)
!   as lines of "alt value"
!
!  NB: If there is no data to do an appropriate interpolation to extract
!      the sought value, then a "missing_value" (taken from the variable's
!      attribute in the input file, most likely -9.99E33) is returned.
!
! EM. Sept. 2011

use netcdf
implicit none

! Input file:
character(len=256) :: infile
character(len=256) :: outfile

character (len=256) :: text ! to store some text
character (len=64) :: varname ! to store the name of the variable to retreive

! NetCDF stuff
integer :: status ! NetCDF return code
integer :: inid ! NetCDF file IDs
integer :: varid ! to store the ID of a variable
integer :: lat_dimid,lon_dimid,alt_dimid,time_dimid
integer :: datashape(4)

real,dimension(:),allocatable :: longitude ! longitude
integer lonlen ! # of grid points along longitude
real,dimension(:),allocatable :: latitude ! latitude
integer latlen ! # of grid points along latitude
real,dimension(:),allocatable :: altitude ! can be geometric heights or pressure
integer altlen ! # of grid point along altitude
real,dimension(:),allocatable :: time ! time
integer timelen ! # of points along time
character :: alttype ! altitude coord. type:'z' (altitude, m) 'p' (pressure, Pa)
real,dimension(:,:,:,:),allocatable :: field
real :: missing_value ! value to denote non-existant data
real :: starttimeoffset ! offset (in sols) wrt Ls=0 of sol 0 in file
integer :: extract_mode ! 1: point-by-point extraction 2: extract a profile  

! point at which data is sought:
real :: lon,lat,alt,Ls,LT,value
real :: sol ! sol GCM date corresponding to sought Ls and LT

integer :: nb

!===============================================================================
! 1.1 Input file
!===============================================================================

write(*,*) ""
write(*,*) " Program valid for diagfi.nc or concatnc.nc files"
write(*,*) "  processed by zrecast "
write(*,*) " Enter input file name:"

read(*,'(a)') infile
write(*,*) ""

! open input file
status=nf90_open(infile,NF90_NOWRITE,inid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to open datafile ",trim(infile)
  write(*,*)trim(nf90_strerror(status))
 stop
endif

write(*,*) " Beginning date of the file?"
write(*,*) " (i.e. number of sols since Ls=0 that the Time=0.0 in the input"
write(*,*) "  file corresponds to)"
read(*,*) starttimeoffset 

write(*,*) " Extraction mode?"
write(*,*) " ( 1: pointwise extraction , 2: profile extraction)"
read(*,*,iostat=status) extract_mode
if ((status.ne.0).or.(extract_mode.lt.1) &
                 .or.(extract_mode.gt.2)) then
  write(*,*) "Error: invalid extraction mode:",extract_mode
  stop
endif

!===============================================================================
! 1.2 Input variable to extract
!===============================================================================

write(*,*) "Enter variable to extract:"
read(*,*) varname
! check that input file contains that variable
status=nf90_inq_varid(inid,trim(varname),varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find variable ",trim(varname)," in ",trim(infile)
  write(*,*)trim(nf90_strerror(status))
  write(*,*) " Might as well stop here."
  stop
endif

!===============================================================================
! 1.3 Get grids for dimensions lon,lat,alt,time
!===============================================================================

! latitude
status=nf90_inq_dimid(inid,"latitude",lat_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find latitude dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,lat_dimid,len=latlen)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find latitude length"
  write(*,*)trim(nf90_strerror(status))
endif
allocate(latitude(latlen))
status=nf90_inq_varid(inid,"latitude",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find latitude ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read latitude
status=nf90_get_var(inid,varid,latitude)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load latitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

!longitude
status=nf90_inq_dimid(inid,"longitude",lon_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find longitude dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,lon_dimid,len=lonlen)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find longitude length"
  write(*,*)trim(nf90_strerror(status))
endif
allocate(longitude(lonlen))
status=nf90_inq_varid(inid,"longitude",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find longitude ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read longitude
status=nf90_get_var(inid,varid,longitude)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load longitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif

!time
status=nf90_inq_dimid(inid,"Time",time_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find Time dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,time_dimid,len=timelen)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find Time length"
  write(*,*)trim(nf90_strerror(status))
endif
allocate(time(timelen))
status=nf90_inq_varid(inid,"Time",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find Time ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read Time
status=nf90_get_var(inid,varid,time)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load Time"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! add the offset to time(:)
time(:)=time(:)+starttimeoffset

!altitude
status=nf90_inq_dimid(inid,"altitude",alt_dimid)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find altitude dimension"
  write(*,*)trim(nf90_strerror(status))
 stop
endif
status=nf90_inquire_dimension(inid,alt_dimid,len=altlen)
if (status.ne.nf90_noerr) then
  write(*,*)"Failed to find altitude length"
  write(*,*)trim(nf90_strerror(status))
endif
allocate(altitude(altlen))
status=nf90_inq_varid(inid,"altitude",varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find altitude ID"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! read altitude
status=nf90_get_var(inid,varid,altitude)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load altitude"
  write(*,*)trim(nf90_strerror(status))
  stop
endif
! check altitude attribute "units" to find out altitude type
status=nf90_get_att(inid,varid,"units",text)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load altitude units attribute"
  write(*,*)trim(nf90_strerror(status))
  stop
else
  if (trim(text).eq."Pa") then
    alttype="p" ! pressure coordinate
  else if (trim(text).eq."m") then
    alttype="z" ! altitude coordinate
  else
    write(*,*)" I do not understand this unit ",trim(text)," for altitude!"
    stop
  endif
endif

!===============================================================================
! 1.3 Get input dataset
!===============================================================================
status=nf90_inq_varid(inid,trim(varname),varid)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to find variable ",trim(varname)," in ",trim(infile)
  write(*,*)trim(nf90_strerror(status))
  write(*,*) " Might as well stop here."
  stop
endif

! sanity checks on the variable
status=nf90_inquire_variable(inid,varid,ndims=nb,dimids=datashape)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to obtain information on variable ",trim(varname)
  write(*,*)trim(nf90_strerror(status))
  write(*,*) " Might as well stop here."
  stop
else
  ! check that it is a 4D variable
  if (nb.ne.4) then
    write(*,*) "Error, expected a 4D (lon-lat-alt-time) variable!"
    stop
  endif
  ! check that its dimensions are indeed lon,lat,alt,time
  if (datashape(1).ne.lon_dimid) then
    write(*,*) "Error, expected first dimension to be longitude!"
    stop
  endif
  if (datashape(2).ne.lat_dimid) then
    write(*,*) "Error, expected second dimension to be latitude!"
    stop
  endif
  if (datashape(3).ne.alt_dimid) then
    write(*,*) "Error, expected third dimension to be altitude!"
    stop
  endif
  if (datashape(4).ne.time_dimid) then
    write(*,*) "Error, expected fourth dimension to be time!"
    stop
  endif
endif

allocate(field(lonlen,latlen,altlen,timelen))

! load dataset
status=nf90_get_var(inid,varid,field)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load ",trim(varname)
  write(*,*)trim(nf90_strerror(status))
  stop
else
  write(*,*) "Loaded ",trim(varname)
endif
! get dataset's missing_value attribute
status=nf90_get_att(inid,varid,"missing_value",missing_value)
if (status.ne.nf90_noerr) then
  write(*,*) "Failed to load missing_value attribute"
  write(*,*)trim(nf90_strerror(status))
  stop
else
  write(*,'("  with missing_value attribute : ",1pe12.5)') missing_value
endif

!===============================================================================
! 2. Process and interpolate
!===============================================================================

!===============================================================================
! 2.1 Create output file
!===============================================================================

outfile=trim(infile(1:index(infile,".nc",back=.true.)-1))//"_"//&
        trim(varname)//".dat"
open(42,file=outfile,form="formatted")
write(*,*) "Output file is: ",trim(outfile)

!===============================================================================
! 2.2 Extract values
!===============================================================================

if (extract_mode==1) then ! pointwise extraction
 write(*,*) " Enter values of: lon lat alt Ls LT (all on the same line,"
 write(*,*) "                    using as many lines as sought data values)"
 write(*,*) "    (enter anything else, e.g. blank line, nonesense,... to quit)"
else ! extract_mode==2 , profile extraction
 write(*,*) " Enter values of: lon lat Ls LT (all on the same line)"
 read(*,*) lon,lat,Ls,LT
 write(*,*) " Enter values of altitude (one per line)"
 write(*,*) "    (enter anything else, e.g. blank line, nonesense,... to quit)"
endif
do 
  ! 2.1 read coordinates and do some sanity checks
  text="" !initialize text to empty string
  read(*,'(a)') text ! store input as text (to spot empty input line)
  if (len_trim(adjustl(text)).eq.0) exit
  
  if (extract_mode==1) then
    read(text,*,iostat=status) lon,lat,alt,Ls,LT
    ! Ls in degrees, LT (local true solar time) in hours, i.e.  in [0:24]
  else ! extract_mode==2 , read alt only
    read(text,*,iostat=status) alt
  endif
  if (status.ne.0) exit

  if ((lon.lt.-360.).or.(lon.gt.360.)) then
    write(*,*) "Unexpected value for lon: ",lon
    stop
  endif
  ! we want lon in [-180:180]
  if (lon.lt.-180.) lon=lon+360.
  if (lon.gt.180.) lon=lon-360.
  
  if ((lat.lt.-90.).or.(lat.gt.90.)) then
    write(*,*) "Unexpected value for lat: ",lat
    stop
  endif
  
  if ((Ls.lt.0.).or.(Ls.gt.360.)) then
    write(*,*) "Unexpected value for Ls: ",Ls
    stop
  endif
  
  if ((LT.lt.0.).or.(LT.gt.24.)) then
    write(*,*) "Unexpected value for LT: ",LT
    stop
  endif
  
  ! 2.2 compute GCM sol date corresponding to sought  Ls and LT
  call ls2sol(Ls,sol)
  !shift 'sol' decimal part to ensure a compatible LT with the desired one
  sol=floor(sol)+(LT-lon/15.)/24.
  ! handle bordeline cases:
  if (sol.gt.669) sol=sol-669
  if (sol.lt.0) sol=sol+669
  
!  write(*,*) " Ls=",Ls," LT=",LT," => sol=",sol
  
  ! 2.3 do the interpolation
  call extraction(lon,lat,alt,sol,&
                  lonlen,latlen,altlen,timelen,&
                  longitude,latitude,altitude,time,&
                  field,missing_value,alttype,varname,value)
  ! 2.4 Write value to output
  if (extract_mode==1) then ! pointwise extraction
    write(42,'(6(1x,1pe12.5))')lon,lat,alt,Ls,LT,value
  else ! profile 
    write(42,'(6(1x,1pe12.5))')alt,value
  endif
enddo ! of do while

! close output file
close(42)

end program extract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extraction(lon,lat,alt,sol,&
                  lonlen,latlen,altlen,timelen,&
                  longitude,latitude,altitude,time,&
                  field,missing_value,alttype,varname,value)

implicit none
! Arguments:
real,intent(in) :: lon  ! sought longitude (deg, in [-180:180])
real,intent(in) :: lat  ! sought latitude (deg, in [-90:90])
real,intent(in) :: alt  ! sought altitude (m or Pa)
real,intent(in) :: sol  ! sought date (sols)
integer,intent(in) :: lonlen
integer,intent(in) :: latlen
integer,intent(in) :: altlen
integer,intent(in) :: timelen
real,intent(in) :: longitude(lonlen)
real,intent(in) :: latitude(latlen)
real,intent(in) :: altitude(altlen)
real,intent(in) :: time(timelen)
real,intent(in) :: field(lonlen,latlen,altlen,timelen)
real,intent(in) :: missing_value ! default value in GCM file for "no data"
character,intent(in) :: alttype ! altitude coord. type:'z' (altitude, m) 
                                !                      'p' (pressure, Pa)
character(len=*),intent(in) :: varname ! variable name (in GCM file)
real,intent(out) :: value

! Local variables:
real,save :: prev_lon=-99 ! previous value of 'lon' routine was called with
real,save :: prev_lat=-99 ! previous value of 'lat' routine was called with
real,save :: prev_alt=-9.e20 ! ! previous value of 'alt'
real,save :: prev_sol=-99 ! previous value of 'sol' routine was called with

! encompasing indexes:
integer,save :: ilon_inf=-1,ilon_sup=-1 ! along longitude
integer,save :: ilat_inf=-1,ilat_sup=-1 ! along latitude
integer,save :: ialt_inf=-1,ialt_sup=-1 ! along altitude
integer,save :: itim_inf=-1,itim_sup=-1 ! along time

! intermediate interpolated values
real :: t_interp(2,2,2) ! after time interpolation
real :: zt_interp(2,2) ! after altitude interpolation
real :: yzt_interp(2) ! after latitude interpolation
real :: coeff ! interpolation coefficient

integer :: i

! By default, set value to missing_value
value=missing_value

! 1. Find encompassing indexes
if (lon.ne.prev_lon) then
  do i=1,lonlen-1
    if (longitude(i).le.lon) then
      ilon_inf=i
    else
      exit
    endif
  enddo
  ilon_sup=ilon_inf+1
endif ! of if (lon.ne.prev_lon)
!write(*,*) 'ilon_inf=',ilon_inf," longitude(ilon_inf)=",longitude(ilon_inf)

if (lat.ne.prev_lat) then
  ! recall that latitudes start from north pole
  do i=1,latlen-1
    if (latitude(i).ge.lat) then
      ilat_inf=i
    else
      exit
    endif
  enddo
  ilat_sup=ilat_inf+1
endif ! of if (lat.ne.prev_lat)
!write(*,*) 'ilat_inf=',ilat_inf," latitude(ilat_inf)=",latitude(ilat_inf)

if (alt.ne.prev_alt) then
  if (alttype.eq.'p') then ! pressures are ordered from max to min
    !handle special case for alt not in the altitude(1:altlen) interval
    if ((alt.gt.altitude(1)).or.(alt.lt.altitude(altlen))) then
      ialt_inf=-1
      ialt_sup=-1
      ! return to main program (with value=missing_value)
      return
    else ! general case
      do i=1,altlen-1
        if (altitude(i).ge.alt) then
          ialt_inf=i
        else
          exit
        endif
      enddo
      ialt_sup=ialt_inf+1
    endif ! of if ((alt.gt.altitude(1)).or.(alt.lt.altitude(altlen)))
  else ! altitudes (m) are ordered from min to max
    !handle special case for alt not in the altitude(1:altlen) interval
    if ((alt.lt.altitude(1)).or.(alt.gt.altitude(altlen))) then
      ialt_inf=-1
      ialt_sup=-1
      ! return to main program (with value=missing_value)
      return
    else ! general case
      do i=1,altlen-1
        if (altitude(i).le.alt) then
          ialt_inf=i
        else
          exit
        endif
      enddo
      ialt_sup=ialt_inf+1
    endif ! of if ((alt.lt.altitude(1)).or.(alt.gt.altitude(altlen)))
  endif ! of if (alttype.eq.'p') 
endif ! of if (alt.ne.prev_alt)
!write(*,*) 'ialt_inf=',ialt_inf," altitude(ialt_inf)=",altitude(ialt_inf)

if (sol.ne.prev_sol) then
  !handle special case for sol not in the time(1):time(timenlen) interval
  if ((sol.lt.time(1)).or.(sol.gt.time(timelen))) then
    itim_inf=-1
    itim_sup=-1
    ! return to main program (with value=missing_value)
    return
  else ! general case
    do i=1,timelen-1
      if (time(i).le.sol) then
        itim_inf=i
      else
        exit
      endif
    enddo
    itim_sup=itim_inf+1
  endif ! of if ((sol.lt.time(1)).or.(sol.gt.time(timenlen)))
endif ! of if (sol.ne.prev_sol)
!write(*,*) 'itim_inf=',itim_inf," time(itim_inf)=",time(itim_inf)
!write(*,*) 'itim_sup=',itim_sup," time(itim_sup)=",time(itim_sup)

! 2. Interpolate
! check that there are no "missing_value" in the field() elements we need
! otherwise return to main program (with value=missing_value)
if (field(ilon_inf,ilat_inf,ialt_inf,itim_inf).eq.missing_value) return
if (field(ilon_inf,ilat_inf,ialt_inf,itim_sup).eq.missing_value) return
if (field(ilon_inf,ilat_inf,ialt_sup,itim_inf).eq.missing_value) return
if (field(ilon_inf,ilat_inf,ialt_sup,itim_sup).eq.missing_value) return
if (field(ilon_inf,ilat_sup,ialt_inf,itim_inf).eq.missing_value) return
if (field(ilon_inf,ilat_sup,ialt_inf,itim_sup).eq.missing_value) return
if (field(ilon_inf,ilat_sup,ialt_sup,itim_inf).eq.missing_value) return
if (field(ilon_inf,ilat_sup,ialt_sup,itim_sup).eq.missing_value) return
if (field(ilon_sup,ilat_inf,ialt_inf,itim_inf).eq.missing_value) return
if (field(ilon_sup,ilat_inf,ialt_inf,itim_sup).eq.missing_value) return
if (field(ilon_sup,ilat_inf,ialt_sup,itim_inf).eq.missing_value) return
if (field(ilon_sup,ilat_inf,ialt_sup,itim_sup).eq.missing_value) return
if (field(ilon_sup,ilat_sup,ialt_inf,itim_inf).eq.missing_value) return
if (field(ilon_sup,ilat_sup,ialt_inf,itim_sup).eq.missing_value) return
if (field(ilon_sup,ilat_sup,ialt_sup,itim_inf).eq.missing_value) return
if (field(ilon_sup,ilat_sup,ialt_sup,itim_sup).eq.missing_value) return

! 2.1 Linear interpolation in time
coeff=(sol-time(itim_inf))/(time(itim_sup)-time(itim_inf))
t_interp(1,1,1)=field(ilon_inf,ilat_inf,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_inf,ialt_inf,itim_sup)- &
                       field(ilon_inf,ilat_inf,ialt_inf,itim_inf))
t_interp(1,1,2)=field(ilon_inf,ilat_inf,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_inf,ialt_sup,itim_sup)- &
                       field(ilon_inf,ilat_inf,ialt_sup,itim_inf))
t_interp(1,2,1)=field(ilon_inf,ilat_sup,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_sup,ialt_inf,itim_sup)- &
                       field(ilon_inf,ilat_sup,ialt_inf,itim_inf))
t_interp(1,2,2)=field(ilon_inf,ilat_sup,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_inf,ilat_sup,ialt_sup,itim_sup)- &
                       field(ilon_inf,ilat_sup,ialt_sup,itim_inf))
t_interp(2,1,1)=field(ilon_sup,ilat_inf,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_inf,ialt_inf,itim_sup)- &
                       field(ilon_sup,ilat_inf,ialt_inf,itim_inf))
t_interp(2,1,2)=field(ilon_sup,ilat_inf,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_inf,ialt_sup,itim_sup)- &
                       field(ilon_sup,ilat_inf,ialt_sup,itim_inf))
t_interp(2,2,1)=field(ilon_sup,ilat_sup,ialt_inf,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_sup,ialt_inf,itim_sup)- &
                       field(ilon_sup,ilat_sup,ialt_inf,itim_inf))
t_interp(2,2,2)=field(ilon_sup,ilat_sup,ialt_sup,itim_inf)+ &
                coeff*(field(ilon_sup,ilat_sup,ialt_sup,itim_sup)- &
                       field(ilon_sup,ilat_sup,ialt_sup,itim_inf))

! 2.2 Vertical interpolation
if (((varname=='rho').or.(varname=='pressure')).and.(alttype=='z')) then
  ! do the interpolation on the log of the quantity
  coeff=(alt-altitude(ialt_inf))/(altitude(ialt_sup)-altitude(ialt_inf))
  zt_interp(1,1)=log(t_interp(1,1,1))+coeff* &
                             (log(t_interp(1,1,2))-log(t_interp(1,1,1)))
  zt_interp(1,2)=log(t_interp(1,2,1))+coeff* &
                             (log(t_interp(1,2,2))-log(t_interp(1,2,1)))
  zt_interp(2,1)=log(t_interp(2,1,1))+coeff* &
                             (log(t_interp(2,1,2))-log(t_interp(2,1,1)))
  zt_interp(2,2)=log(t_interp(2,2,1))+coeff* &
                             (log(t_interp(2,2,2))-log(t_interp(2,2,1)))
  zt_interp(1:2,1:2)=exp(zt_interp(1:2,1:2))
else ! general case
  coeff=(alt-altitude(ialt_inf))/(altitude(ialt_sup)-altitude(ialt_inf))
  zt_interp(1,1)=t_interp(1,1,1)+coeff*(t_interp(1,1,2)-t_interp(1,1,1))
  zt_interp(1,2)=t_interp(1,2,1)+coeff*(t_interp(1,2,2)-t_interp(1,2,1))
  zt_interp(2,1)=t_interp(2,1,1)+coeff*(t_interp(2,1,2)-t_interp(2,1,1))
  zt_interp(2,2)=t_interp(2,2,1)+coeff*(t_interp(2,2,2)-t_interp(2,2,1))
endif

! 2.3 Latitudinal interpolation
coeff=(lat-latitude(ilat_inf))/(latitude(ilat_sup)-latitude(ilat_inf))
yzt_interp(1)=zt_interp(1,1)+coeff*(zt_interp(1,2)-zt_interp(1,1))
yzt_interp(2)=zt_interp(2,1)+coeff*(zt_interp(2,2)-zt_interp(2,1))

! 2.4 longitudinal interpolation
coeff=(lon-longitude(ilon_inf))/(longitude(ilon_sup)-longitude(ilon_inf))
value=yzt_interp(1)+coeff*(yzt_interp(2)-yzt_interp(1))

end subroutine extraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ls2sol(ls,sol)

implicit none
!  Arguments:
real,intent(in) :: ls
real,intent(out) :: sol

!  Local:
double precision xref,zx0,zteta,zz
!xref: mean anomaly, zteta: true anomaly, zx0: eccentric anomaly
double precision year_day 
double precision peri_day,timeperi,e_elips
double precision pi,degrad 
parameter (year_day=668.6d0) ! number of sols in a martian year
parameter (peri_day=485.35d0) ! date (in sols) of perihelion
!  timeperi: 2*pi*( 1 - Ls(perihelion)/ 360 ); Ls(perihelion)=250.99
parameter (timeperi=1.90258341759902d0)
parameter (e_elips=0.0934d0)  ! eccentricity of orbit
parameter (pi=3.14159265358979d0)
parameter (degrad=57.2957795130823d0)

      if (abs(ls).lt.1.0e-5) then
         if (ls.ge.0.0) then
            sol = 0.0
         else
            sol = real(year_day)
         end if
         return
      end if

      zteta = ls/degrad + timeperi
      zx0 = 2.0*datan(dtan(0.5*zteta)/dsqrt((1.+e_elips)/(1.-e_elips)))
      xref = zx0-e_elips*dsin(zx0)
      zz = xref/(2.*pi)
      sol = real(zz*year_day + peri_day)
      if (sol.lt.0.0) sol = sol + real(year_day)
      if (sol.ge.year_day) sol = sol - real(year_day)


end subroutine ls2sol

