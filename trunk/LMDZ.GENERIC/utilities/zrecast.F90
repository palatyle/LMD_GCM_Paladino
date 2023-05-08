module planet_const
! planetary constants (set via init_planet_const routine)
implicit none
real :: a0 ! Mean planetary radius (a0=3396.E3 for Mars)
real :: g0 ! gravity at a0  (Mars: g0=3.7257964 ; Lemoine et al. 2001)
real :: Rmean ! reduced gaz constant (Mars: R=191)
end module planet_const

program zrecast

! This program reads 4D (lon-lat-alt-time) fields from GCM output files
! (ie: diagfi.nc time series or concat.nc or stats.nc files) and, by
! integrating the hydrostatic equation, recasts data along the vertical
! direction.
! The vertical coordinate can be either 1) pressure, 2) above areoid
! altitudes, 3) above local surface altitudes  or 4) distance to center of
! the planet. Some interpolation along the vertical direction is also
! done, following instructions given by user (levels may be specified
! or given as minimu,maximum and number of levels).
! For "above areoid altitudes" output, Atmospheric pressure is added to
! output dataset; for "pressure coordinate" outputs, the above areoid
! altitude of pressure is added to output dataset.
!
! Minimal requirements and dependencies:
! The dataset must include the following data:
! - surface pressure
! - atmospheric temperature
! - hybrid coordinates aps() and bps(), or sigma levels() (see section 1.3.2)
! - ground geopotential (in input file; if not found, it is sought
!   in a 'diagfi.nc' file. If not found there, it is then sought in
!   a 'phisinit.nc' file  (see section 1.3.3 of program)
!
! - When integration the hydrostatic equation, we assume that R, the molecular
!   Gas Constant, may not be constant, so it is computed as
!      R=P/(rho*T) (P=Pressure, rho=density, T=temperature)
!   If 'rho' is not available, then we use a constant R (see section 2.2)
!
! WARNING: Asking for many points along the vertical direction quickly
!          leads to HUGE output files.
!
! EM 01/2006 : Corrected a bug in vertical (log) interpolation for pressure
!              and density
! EM 10/2006 : Modified program so that it can now process 'stats.nc'
!              files obtained from British GCM (ie: vertical coordinate
!              given as sigma levels and geopotential read from file
!              'phisinit.nc')
! EM 02/2007 : Changed behavior for "altitude above surface" case
!              (for MCD RMS computation). Number of levels is then set as 
!              number of levels in initial file,
!              and the new set of above surface levels follow a more elaborate
!              distribution (see build_zs.F90 routine).
! EM 08/2009 : User may now specify values of each vertical level,
!              or just min,max and number of levels (as before)
! EM 01/2010 : Corrected bug in 'zs_coord_interp' to correctly handle the case
!              when interpolating log-wise and the density field is not
!              available.
! EM 03/2011 : Add possibility to have outputs as distance to center
!              of planet
! EM 11/2012 : Adapted so it can be used on "generic" model outputs; planet
!              constants (radius, R, etc.) are now read from file
! TN 01/2013 : Adapted for large output files with at least 2 variables > 2 GiB
! TN 10/2013 : Read and write controle field, if available
!
implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=128) :: infile ! input file name (diagfi.nc or stats.nc format)
character (len=128) :: infile2 ! second input file (may be needed for 'phisini')
character (len=128) :: outfile ! output file name

character (len=64) :: text ! to store some text
character (len=64) :: tmpvarname ! temporarily store a variable name
integer tmpvarid ! temporarily store a variable ID
integer tmpdimid ! temporarily store a dimension ID
integer tmpndims ! temporarily store # of dimensions of a variable
integer infid ! NetCDF input file ID (of diagfi.nc or stats.nc format)
integer infid2 ! NetCDF input file which contains 'phisini' dataset (diagfi.nc)
integer nbvarinfile ! # of variables in input file
integer nbattr ! # of attributes of a given variable in input file
integer nbvar4dinfile ! # of 4D (lon,lat,alt,time) variables in input file
integer outfid ! NetCDF output file ID
integer lon_dimid,lat_dimid,alt_dimid,time_dimid,ctl_dimid ! NetCDF dimension IDs
integer lon_varid,lat_varid,alt_varid,time_varid,ctl_varid
integer gcm_layers_dimid ! NetCDF dimension ID for # of layers in GCM
integer sigma_varid,aps_varid,bps_varid
integer za_varid,p_varid ! above areoid and pressure data IDs



integer ps_varid ! surface pressure data ID
integer,dimension(4) :: datashape ! shape of 4D datasets
integer,dimension(3) :: surfdatashape ! shape of 3D (surface+time) datasets

real :: miss_val=-9.99e+33 ! special "missing value" to specify missing data
real,parameter :: miss_val_def=-9.99e+33 ! default value for "missing value"
character (len=64), dimension(:), allocatable :: var
! var(): names of variables that will be processed
integer nbvar ! # of variables to process
integer,dimension(:),allocatable :: var_id ! IDs of variables var() (in outfile)
real,dimension(:),allocatable :: lon ! longitude
integer lonlength ! # of grid points along longitude
real,dimension(:),allocatable :: lat ! latitude
integer latlength ! # of grid points along latitude
integer altlength ! # of grid point along altitude (of input datasets)
real,dimension(:),allocatable :: ctl ! controle
integer ctllength ! # of grid points along controle
real,dimension(:),allocatable :: time ! time
integer timelength ! # of points along time
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:),allocatable :: sigma ! sigma levels
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
real,dimension(:,:,:),allocatable :: ps ! GCM surface pressure
real,dimension(:,:,:,:),allocatable :: press ! GCM atmospheric pressure 
real,dimension(:,:,:,:),allocatable :: temp ! GCM atmospheric temperature
real,dimension(:,:,:,:),allocatable :: teta ! GCM atmospheric potential temperature
real,dimension(:,:,:,:),allocatable :: rho ! GCM atmospheric density
real,dimension(:,:,:,:),allocatable :: za_gcm ! GCM above areoid levels (m)
real,dimension(:,:,:,:),allocatable :: zs_gcm ! GCM above surface heights (m)
real,dimension(:,:,:,:),allocatable :: ra_gcm ! radial distance (m) to
                                          ! center of Mars, (of GCM grid points)
real,dimension(:,:,:,:),allocatable :: indata ! to store a GCM 4D dataset
real,dimension(:,:,:,:),allocatable :: outdata ! to store a 4D dataset
integer ierr ! NetCDF routines return code
integer i,j,ilon,ilat,ilev,itim ! for loops

integer ztype ! Flag for vertical coordinate of output
! ztype=1: pressure  ztype=2: above areoid  ztype=3: above local surface
! ztype=4: distance to center of planet
integer nblev ! # of levels (along vertical coordinate) for output data
real pmin,pmax ! min and max values for output pressure coordinate
real,dimension(:),allocatable :: plevel ! Pressure levels for output data
real zamin,zamax ! min and max values for output above areoid coordinate
real,dimension(:),allocatable :: zareoid ! Above areoid heights for output data
real,dimension(:),allocatable :: zsurface ! Above surface heights for output
real,dimension(:),allocatable :: zradius ! distance to center of planet
logical :: have_rho ! Flag: true if density 'rho' is available
logical :: have_sigma ! Flag: true if sigma levels are known (false if hybrid
                      !       coordinates are used)
logical :: auto_vert_levels ! Flag: true if the positions of vertical levels
                            ! has to be computed; false if these are given
                            ! by the user
logical :: auto_mcd_levels=.false. ! Flag: specific case for MCD automatic 
                           ! above local surface levels
integer,dimension(4) :: edges,corner ! needed to write variables for big files

!===============================================================================
! 1. Input parameters
!===============================================================================

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
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! load planet constants (radius, gravity, ...)

call init_planet_const(infid)
 
!===============================================================================
! 1.2 Get # and names of variables in input file
!===============================================================================

ierr=NF_INQ_NVARS(infid,nbvarinfile)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'ERROR: Failed geting number of variables from file'
  write(*,*) NF_STRERROR(ierr)
  stop
endif

write(*,*)" The following variables have been found:"
nbvar4dinfile=0
do i=1,nbvarinfile
  ! get name of variable # i
  ierr=NF_INQ_VARNAME(infid,i,tmpvarname)
  ! check if it is a 4D variable
  ierr=NF_INQ_VARNDIMS(infid,i,tmpndims)
  if (tmpndims.eq.4) then
    nbvar4dinfile=nbvar4dinfile+1
    write(*,*) trim(tmpvarname)
  endif
enddo

allocate(var(nbvar4dinfile),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error! failed memory allocation of var(nbvar4dinfile)"
  stop
endif

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
    ! also check that it is indeed a 4D variable
    ierr=NF_INQ_VARNDIMS(infid,tmpvarid,tmpndims)
    if (tmpndims.eq.4) then
      nbvar=nbvar+1
      var(nbvar)=tmpvarname    
    else
      write(*,*) 'Error: ',trim(tmpvarname),' is not a 4D variable'
      write(*,*) '       (we''ll skip that one)'
    endif
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
! 1.3 Get grids in lon,lat,alt,time,
!     as well as hybrid coordinates aps() and bps() (or sigma levels sigma())
!     and phisinit() from input file
!===============================================================================

! 1.3.1 longitude, latitude, altitude and time

! latitude
ierr=NF_INQ_DIMID(infid,"latitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get latitude dimension ID"
  write(*,*) "  looking for lat dimension instead... "
  ierr=NF_INQ_DIMID(infid,"lat",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get lat dimension ID"
  else
    ierr=NF_INQ_VARID(infid,"lat",tmpvarid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get lat ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,latlength)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed to get lat length"
      else
        allocate(lat(latlength),stat=ierr)
        if (ierr.ne.0) then
          write(*,*) "Error: Failed to allocate lat(latlength)"
          write(*,*) "     latlength=",latlength
          stop
        endif
        ierr=NF_GET_VAR_REAL(infid,tmpvarid,lat)
        if (ierr.ne.NF_NOERR) then
          stop "Error: Failed reading lat"
        endif
      endif
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"latitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get latitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,latlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get latitude length"
    else
      allocate(lat(latlength),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate lat(latlength)"
        write(*,*) "     latlength=",latlength
        stop
      endif
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
  write(*,*) "Could not get longitude dimension ID"
  write(*,*) "  looking for lon dimension instead... "
  ierr=NF_INQ_DIMID(infid,"lon",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get lon dimension ID"
  else
    ierr=NF_INQ_VARID(infid,"lon",tmpvarid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get lon ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,lonlength)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed to get lon length"
      else
        allocate(lon(lonlength),stat=ierr)
        if (ierr.ne.0) then
          write(*,*) "Error: Failed to allocate lon(lonlength)"
          write(*,*) "     lonlength=",lonlength
          stop
        endif
        ierr=NF_GET_VAR_REAL(infid,tmpvarid,lon)
        if (ierr.ne.NF_NOERR) then
          stop "Error: Failed reading longitude"
        endif
      endif
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"longitude",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get longitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,lonlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get longitude length"
    else
      allocate(lon(lonlength),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate lon(lonlength)"
        write(*,*) "     lonlength=",lonlength
        stop
      endif
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
  write(*,*) "Could not get Time dimension ID"
  write(*,*) "  looking for time dimension instead... "
  ierr=NF_INQ_DIMID(infid,"time",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get lon dimension ID"
  else
    ierr=NF_INQ_VARID(infid,"time",tmpvarid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get time ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed to get time length"
      else
        allocate(time(timelength),stat=ierr)
        if (ierr.ne.0) then
          write(*,*) "Error: Failed to allocate time(timelength)"
          write(*,*) "     timelength=",timelength
          stop
        endif
        ierr=NF_GET_VAR_REAL(infid,tmpvarid,time)
        if (ierr.ne.NF_NOERR) then
          stop "Error: Failed reading time"
        endif
      endif
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"Time",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get Time ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get Time length"
    else
      allocate(time(timelength),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate time(timelength)"
        write(*,*) "     timelength=",timelength
        stop
      endif
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
  write(*,*) "Could not get altitude dimension ID"
  write(*,*) "  looking for sigma dimension instead... "
  ierr=NF_INQ_DIMID(infid,"sigma",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get sigma dimension ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get altitude length"
    endif
  endif
else
  ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
  if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get altitude length"
  endif
endif

! controle
ierr=NF_INQ_DIMID(infid,"index",tmpdimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get controle dimension ID"
  ctllength = 0
else
  ierr=NF_INQ_VARID(infid,"controle",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get controle ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,ctllength)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get controle length"
    else
      allocate(ctl(ctllength),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate ctl(ctllength)"
        write(*,*) "     ctllength=",ctllength
        stop
      endif
      ierr=NF_GET_VAR_REAL(infid,tmpvarid,ctl)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed reading controle"
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
  allocate(sigma(altlength),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "Error: Failed to allocate sigma(altlength)"
    write(*,*) "     altlength=",altlength
    stop
  endif
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
    allocate(aps(altlength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate aps(altlength)"
      write(*,*) "     altlength=",altlength
      stop
    endif
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
    allocate(bps(altlength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate bps(altlength)"
      write(*,*) "     altlength=",altlength
      stop
    endif
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading bps"
    endif
  endif
endif !of if (.not.have_sigma)

! 1.3.3 ground geopotential phisinit

allocate(phisinit(lonlength,latlength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Failed allocation of phisinit(lonlength,latlength) !!!"
  write(*,*) "lonlength=",lonlength," latlength=",latlength
  stop
endif
! look for 'phisinit' in current file
ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid) 
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get phisinit ID from file ",trim(infile)
  infile2="diagfi.nc"
  write(*,*) "         Trying file ",trim(infile2)
  ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Problem: Could not find/open that file"
    infile2="diagfi1.nc"
    write(*,*) "         Trying file ",trim(infile2)
    ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "Problem: Could not find/open that file"
      infile2="phisinit.nc"
      write(*,*) "         Trying file ",trim(infile2)
      ierr=NF_OPEN(infile2,NF_NOWRITE,infid2)
      if (ierr.ne.NF_NOERR) then
        write(*,*) "Error: Could not open that file either"
        stop "Might as well stop here"
      endif
    endif
  endif

  ! Get ID for phisinit
  ierr=NF_INQ_VARID(infid2,"phisinit",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get phisinit ID"
  endif
  ! Get physinit
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

!===============================================================================
! 1.4 Choose and build the new vertical coordinate
!===============================================================================

write(*,*) ""
write(*,*) "Which vertical coordinate should the output be in?"
ztype=0
do while ((ztype.lt.1).or.(ztype.gt.4))
  write(*,*) "(1: pressure, 2: altitude above areoid, "
  write(*,*) " 3: above local surface, 4: distance to center of Mars)"
  read(*,*)ztype
enddo

text="dummy" ! dummy initialization
auto_vert_levels=.true. ! dummy initialization to get rid of compiler warning
do while ((trim(text).ne."yes").and.(trim(text).ne."no").and.(trim(text).ne."mcd"))
  write(*,*) ""
  write(*,*) "Automatic generation of vertical levels distribution? (yes/no)"
  write(*,*) "(yes: you only provide min, max and number of levels)"
  write(*,*) "(no: you provide values for each level)"
  write(*,*) "(mcd: special case for above local surface levels for MCD)"
  read(*,'(a64)') text
  if (trim(text).eq."yes") then
    auto_vert_levels=.true.
    auto_mcd_levels=.false.
  elseif(trim(text).eq."no") then
    auto_vert_levels=.false.
    auto_mcd_levels=.false.
  elseif(trim(text).eq."mcd") then
    auto_vert_levels=.true.
    auto_mcd_levels=.true.
  endif
enddo

if (auto_vert_levels) then
  ! ask for # of points and end values for pressure or above areoid cases
  ! except in MCD case
  if (.not.auto_mcd_levels) then
    write(*,*) ""
    write(*,*) "Enter min and max of vertical coordinate (Pa or m)"
    write(*,*) " (in that order and on the same line)"
    if (ztype.eq.1) then ! pressure coordinate
      read(*,*) pmin,pmax
    else ! altitude coordinate, except in MCD case
      read(*,*) zamin,zamax
    endif
  endif

  ! Build corresponding vertical coordinates
  if (ztype.eq.1) then ! pressure coordinate
    write(*,*) "Number of levels along vertical coordinate?"
    read(*,*) nblev
    allocate(plevel(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate plevel(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    if (nblev.eq.1) then ! in case only one level is asked for
      plevel(nblev)=pmin
    else
      do i=1,nblev
        !    build exponentially spread layers
        plevel(i)=exp(log(pmax)+(log(pmin)-log(pmax))* &
                      ((real(i)-1.0)/(real(nblev)-1.0)))
      enddo
    endif
  else if (ztype.eq.2) then ! above areoid heights
    write(*,*) "Number of levels along vertical coordinate?"
    read(*,*) nblev
    allocate(zareoid(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zareoid(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    if (nblev.eq.1) then ! in case only one level is asked for
      zareoid(nblev)=zamin
    else 
      do i=1,nblev
        zareoid(i)=zamin+(real(i)-1.0)*((zamax-zamin)/(real(nblev)-1.0))
      enddo
    endif
  else if (ztype.eq.3) then ! above local surface
    if (auto_mcd_levels) then
      ! set nblev to # of levels from input data files
      nblev=altlength
      allocate(zsurface(nblev),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate zsurface(nblev)"
        write(*,*) "     nblev=",nblev
        stop
      endif
      ! build specific above local surface altitudes
      call build_zs(nblev,have_sigma,sigma,aps,bps,zsurface)
    else
      write(*,*) "Number of levels along vertical coordinate?"
      read(*,*) nblev
      allocate(zsurface(nblev),stat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error: Failed to allocate zsurface(nblev)"
        write(*,*) "     nblev=",nblev
        stop
      endif
      if (nblev.eq.1) then ! in case only one level is asked for
        zsurface(nblev)=zamin
      else 
        do i=1,nblev
          zsurface(i)=zamin+(real(i)-1.0)*((zamax-zamin)/(real(nblev)-1.0))
        enddo
      endif
    endif ! of if (auto_mcd_levels)
  else ! distance to center of planet
   write(*,*) "Number of levels along vertical coordinate?"
    read(*,*) nblev
    allocate(zradius(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zradius(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    if (nblev.eq.1) then ! in case only one level is asked for
      zradius(nblev)=zamin
    else 
      do i=1,nblev
        zradius(i)=zamin+(real(i)-1.0)*((zamax-zamin)/(real(nblev)-1.0))
      enddo
    endif
  endif
else ! auto_vert_levels=.false. ; user provides values
  ! ask for # of points along the vertical
  write(*,*) ""
  write(*,*) "Number of levels along vertical coordinate?"
  read(*,*) nblev
  if (ztype.eq.1) then ! pressure coordinate
    allocate(plevel(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate plevel(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    write(*,*) "Enter Pressure (Pa) of levels, ordered"
    write(*,*) " from max (near-surface) to min (top of atmosphere),"
    write(*,*) " (one value per line)"
    do i=1,nblev
      read(*,*) plevel(i)
    enddo
  else if (ztype.eq.2) then ! above areoid heights
    allocate(zareoid(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zareoid(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    write(*,*) "Enter altitude (m) above areoid of levels, ordered"
    write(*,*) " from min to max (one value per line)"
    do i=1,nblev
      read(*,*) zareoid(i)
    enddo
  else if (ztype.eq.3) then ! above local surface
    allocate(zsurface(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zsurface(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    write(*,*) "Enter altitude (m) above surface of levels, ordered"
    write(*,*) " from min to max (one value per line)"
    do i=1,nblev
      read(*,*) zsurface(i)
    enddo
  else ! distance to center of planet
    allocate(zradius(nblev),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zradius(nblev)"
      write(*,*) "     nblev=",nblev
      stop
    endif
    write(*,*) "Enter distance (m) to center of Mars of levels, ordered"
    write(*,*) " from min to max (one value per line)"
    do i=1,nblev
      read(*,*) zradius(i)
    enddo
  endif
endif ! of if (auto_vert_levels)

!===============================================================================
! 1.5 Get output file name
!===============================================================================
write(*,*) ""
!write(*,*) "Enter output file name"
!read(*,*) outfile
if (ztype.eq.1) then ! pressure coordinate
  outfile=infile(1:len_trim(infile)-3)//"_P.nc"
else if (ztype.eq.2) then ! above areoid coordinate
  outfile=infile(1:len_trim(infile)-3)//"_A.nc"
else if (ztype.eq.3) then ! above local surface
  outfile=infile(1:len_trim(infile)-3)//"_S.nc"
else ! radial distance to center of Mars
  outfile=infile(1:len_trim(infile)-3)//"_R.nc"
endif
write(*,*) "Output file name is: "//trim(outfile)

!===============================================================================
! 2.1 Build/store GCM fields which will be used later
!===============================================================================

!===============================================================================
! 2.1.1 Surface pressure
!===============================================================================
ierr=NF_INQ_VARID(infid,"ps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  stop "Error: Failed to get ps ID"
else
  allocate(ps(lonlength,latlength,timelength),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "Error: Failed to allocate ps(lonlength,latlength,timelength)"
    write(*,*) "       lonlength=",lonlength," latlength=",latlength
    write(*,*) "       timelength=",timelength
    stop
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,ps)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading surface pressure"
  endif
endif

!===============================================================================
! 2.1.2 Atmospheric pressure
!===============================================================================
allocate(press(lonlength,latlength,altlength,timelength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate press(lonlength,latlength,altlength,timelength)"
  write(*,*) "       lonlength=",lonlength," latlength=",latlength
  write(*,*) "       altlength=",altlength," timelength=",timelength
  stop
endif

if (have_sigma) then ! sigma coordinate
  do itim=1,timelength
    do ilev=1,altlength
      do ilat=1,latlength
        do ilon=1,lonlength
          press(ilon,ilat,ilev,itim)=sigma(ilev)*ps(ilon,ilat,itim)
        enddo
      enddo
    enddo
  enddo
else ! hybrid coordinates
  do itim=1,timelength
    do ilev=1,altlength
      do ilat=1,latlength
        do ilon=1,lonlength
          press(ilon,ilat,ilev,itim)=aps(ilev)+bps(ilev)*ps(ilon,ilat,itim)
        enddo
      enddo
    enddo
  enddo
endif

!===============================================================================
! 2.1.3 Atmospheric temperature
!===============================================================================
allocate(temp(lonlength,latlength,altlength,timelength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate temp(lonlength,latlength,altlength,timelength)"
  write(*,*) "       lonlength=",lonlength," latlength=",latlength
  write(*,*) "       altlength=",altlength," timelength=",timelength
  stop
endif

allocate(teta(lonlength,latlength,altlength,timelength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate teta(lonlength,latlength,altlength,timelength)"
  write(*,*) "       lonlength=",lonlength," latlength=",latlength
  write(*,*) "       altlength=",altlength," timelength=",timelength
  stop
endif


ierr=NF_INQ_VARID(infid,"temp",tmpvarid)
if (ierr.ne.NF_NOERR) then
  ! stop "Error: Failed to get temp ID"
  ! try "t" for temperature
  ierr=NF_INQ_VARID(infid,"t",tmpvarid)
  if (ierr.ne.NF_NOERR) then
     ierr=NF_INQ_VARID(infid,"teta",tmpvarid)
     if (ierr.ne.NF_NOERR) then
       stop "Error: Failed to get t or teta ID"
     endif
       ierr=NF_GET_VAR_REAL(infid,tmpvarid,teta)
       if (ierr.ne.NF_NOERR) then
         stop "Error: Failed reading atmospheric temperature"
       endif
     
  do itim=1,timelength
    do ilev=1,altlength
      do ilat=1,latlength
        do ilon=1,lonlength
          temp(ilon,ilat,ilev,itim)=teta(ilon,ilat,ilev,itim)*(press(ilon,ilat,ilev,itim)/ps(ilon,ilat,itim))**(.256793)
        enddo
      enddo
    enddo
  enddo

  endif
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,temp)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading atmospheric temperature"
  endif
endif

!===============================================================================
! 2.1.4 Atmospheric density
!===============================================================================

ierr=NF_INQ_VARID(infid,"rho",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Warning: Failed to get rho ID"
  have_rho=.false.
else
  have_rho=.true.
  allocate(rho(lonlength,latlength,altlength,timelength),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "Error: Failed to allocate rho(lonlength,latlength,altlength,timelength)"
    write(*,*) "       lonlength=",lonlength," latlength=",latlength
    write(*,*) "       altlength=",altlength," timelength=",timelength
    stop
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,rho)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed reading atmospheric density"
  endif
endif

!===============================================================================
! 2.2 Build GCM Above areoid (or above surface) altitudes of GCM nodes
!===============================================================================

if (have_rho) then
  if (ztype.le.2) then ! above areoid altitudes (also needed for ztype=1)
    allocate(za_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate za_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call build_gcm_za(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,rho,za_gcm)
  else if (ztype.eq.3) then! above local surface altitudes
    allocate(zs_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zs_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call build_gcm_zs(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,rho,zs_gcm)
  else ! radial distance to center of planet
    allocate(za_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate za_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    allocate(ra_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate ra_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call build_gcm_za(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,rho,za_gcm)
    ! and now convert the altitudes above areoid to radial distances
    call build_gcm_ra(lonlength,latlength,altlength,timelength, &
                      lon,lat,za_gcm,ra_gcm)
   ! free za_gcm
   deallocate(za_gcm)
  endif
else
  write(*,*)"Warning: Using constant R to integrate hydrostatic equation"
  if (ztype.le.2) then ! above areoid altitudes (also needed for ztype=1)
    allocate(za_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate za_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call crude_gcm_za(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,za_gcm)
  else if (ztype.eq.3) then ! above local surface altitudes
    allocate(zs_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate zs_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call crude_gcm_zs(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,zs_gcm)
  else ! radial distance to center of planet
    allocate(za_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate za_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    allocate(ra_gcm(lonlength,latlength,altlength,timelength),stat=ierr)
    if (ierr.ne.0) then
      write(*,*) "Error: Failed to allocate ra_gcm(lonlength,latlength,altlength,timelength)"
      write(*,*) "       lonlength=",lonlength," latlength=",latlength
      write(*,*) "       altlength=",altlength," timelength=",timelength
      stop
    endif
    call crude_gcm_za(lonlength,latlength,altlength,timelength, &
                      phisinit,ps,press,temp,za_gcm)
    ! and now convert the altitudes above areoid to radial distances
    call build_gcm_ra(lonlength,latlength,altlength,timelength, &
                      lon,lat,za_gcm,ra_gcm)
   ! free za_gcm
   deallocate(za_gcm)
  endif
endif ! of if (have_rho)

!===============================================================================
! 3. Create output file and initialize definitions of variables and dimensions
!===============================================================================

!===============================================================================
! 3.1. Output file
!===============================================================================

! Create output file
ierr=NF_CREATE(outfile,IOR(NF_CLOBBER,NF_64BIT_OFFSET),outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"Error: could not create file ",outfile
  write(*,*) NF_STRERROR(ierr)
  stop
endif

!===============================================================================
! 3.2. Define dimensions
!===============================================================================
! longitude
ierr=NF_DEF_DIM(outfid,"longitude",lonlength,lon_dimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define longitude dimension"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! latitude
ierr=NF_DEF_DIM(outfid,"latitude",latlength,lat_dimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define latitude dimension"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! altitude
ierr=NF_DEF_DIM(outfid,"altitude",nblev,alt_dimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define altitude dimension"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! time
ierr=NF_DEF_DIM(outfid,"Time",NF_UNLIMITED,time_dimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define Time dimension"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! controle
if (ctllength .ne. 0) then
  ierr=NF_DEF_DIM(outfid,"index",ctllength,ctl_dimid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define controle dimension"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif

! GCM layers (for sigma or aps and bps)
ierr=NF_DEF_DIM(outfid,"GCM_layers",altlength,gcm_layers_dimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define GCM_layers dimension"
  write(*,*) NF_STRERROR(ierr)
  stop
endif


!===============================================================================
! 3.3. Define variables and their attributes
!===============================================================================

! 3.3.1 Define 1D variables

! longitude 
datashape(1)=lon_dimid
ierr=NF_DEF_VAR(outfid,"longitude",NF_REAL,1,datashape(1),lon_varid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define longitude variable"
  write(*,*) NF_STRERROR(ierr)
  stop
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
  write(*,*) "Error: Could not define latitude variable"
  write(*,*) NF_STRERROR(ierr)
  stop
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
  write(*,*) "Error: Could not define altitude variable"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

!altitude attributes
if (ztype.eq.1) then ! pressure vertical coordinate
  text='Pressure levels'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for altitude"
  endif
  text='Pa'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for altitude"
  endif
  text='down'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing positive for altitude"
  endif
else if (ztype.eq.2) then ! above areoid vertical coordinate
  text='Altitude above areoid'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for altitude"
  endif
  text='m'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for altitude"
  endif
  text='up'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing positive for altitude"
  endif
else if (ztype.eq.3) then ! above surface vertical coordinate
  text='Altitude above local surface'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for altitude"
  endif
  text='m'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for altitude"
  endif
  text='up'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing positive for altitude"
  endif
else ! radial distance to center of planet
  text='Distance to center of Mars'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for altitude"
  endif
  text='m'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for altitude"
  endif
  text='up'
  ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing positive for altitude"
  endif
endif ! of if (have_sigma)

! controle
if (ctllength .ne. 0) then
  ierr=NF_DEF_VAR(outfid,"controle",NF_REAL,1,ctl_dimid,ctl_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define controle variable"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif

  ! controle attributes
  text='Control parameters'
  ierr=NF_PUT_ATT_TEXT(outfid,ctl_varid,'title',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing title for controle"
  endif
endif

! GCM_layers
!ierr=NF_DEF_VAR(outfid,"gcm_layers",NF_REAL,1,,sigma_varid)

! sigma levels or hybrid coordinates
if (have_sigma) then
  ierr=NF_DEF_VAR(outfid,"sigma",NF_REAL,1,gcm_layers_dimid,sigma_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define sigma variable"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
else ! hybrid coordinates
  ierr=NF_DEF_VAR(outfid,"aps",NF_REAL,1,gcm_layers_dimid,aps_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define aps variable"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ierr=NF_DEF_VAR(outfid,"bps",NF_REAL,1,gcm_layers_dimid,bps_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define bps variable"
    write(*,*) NF_STRERROR(ierr)
    stop
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
  write(*,*) "Error: Could not define Time variable"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! time attributes
text='Time'
ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing long_name for Time"
endif
text='days since 0000-01-1 00:00:00'
ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing units for Time"
endif

! 3.3.2 Define 3D variables (ie: surface+time variables)

! Surface pressure
surfdatashape(1)=lon_dimid
surfdatashape(2)=lat_dimid
surfdatashape(3)=time_dimid
ierr=NF_DEF_VAR(outfid,"ps",NF_REAL,3,surfdatashape,ps_varid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not define ps variable"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! Surface pressure attributes
text='Surface pressure'
ierr=NF_PUT_ATT_TEXT(outfid,ps_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing long_name for surface pressure"
endif
text='Pa'
ierr=NF_PUT_ATT_TEXT(outfid,ps_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) then
  stop "Error: Problem writing units for surface pressure"
endif

! 3.3.3 Define 4D variables

! add pressure or zareoid
if (ztype.eq.1) then ! pressure vertical coordinate
  ! zareoid dataset
  ierr=NF_DEF_VAR(outfid,"zareoid",NF_REAL,4,datashape,za_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define zareoid variable"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ! zareoid attributes
  text='altitude above areoid'
  ierr=NF_PUT_ATT_TEXT(outfid,za_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for zareoid"
  endif
  text='m'
  ierr=NF_PUT_ATT_TEXT(outfid,za_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for zareoid"
  endif
  ! zareoid missing value
  ierr=NF_PUT_ATT_REAL(outfid,za_varid,'missing_value',NF_REAL,1,miss_val)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing missing_value for zareoid"
  endif
else ! above areoid or above local surface vertical coordinate
  ! pressure dataset
  ierr=NF_DEF_VAR(outfid,"pressure",NF_REAL,4,datashape,p_varid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not define pressure variable"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ! pressure attributes
  text='Atmospheric pressure'
  ierr=NF_PUT_ATT_TEXT(outfid,p_varid,'long_name',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing long_name for pressure"
  endif
  text='Pa'
  ierr=NF_PUT_ATT_TEXT(outfid,p_varid,'units',len_trim(text),text)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing units for pressure"
  endif
  ! pressure missing value
  ierr=NF_PUT_ATT_REAL(outfid,p_varid,'missing_value',NF_REAL,1,miss_val)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Problem writing missing_value for pressure"
  endif
endif

! add zs_gcm
if (ztype.eq.3) then
endif

! variables requested by user
allocate(var_id(nbvar),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate var_id(nbvar)"
  write(*,*) "       nbvar=",nbvar
  stop
endif
do i=1,nbvar
  write(*,*) ""
  write(*,*) "Creating ",trim(var(i))
  ! define the variable
  ierr=NF_DEF_VAR(outfid,var(i),NF_REAL,4,datashape,var_id(i))
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, could not define variable ',trim(var(i))
    write(*,*) NF_STRERROR(ierr)
    stop 
  endif

  ! Get the (input file) ID for the variable and
  ! the # of attributes associated to that variable
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to get ID for variable ',trim(var(i))
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ierr=NF_INQ_VARNATTS(infid,tmpvarid,nbattr)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, could not get number of attributes for variable ',&
               trim(var(i))
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ! inititialize j == number of attributes written to output
  j=0

  ! look for a "long_name" attribute
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
    else ! no 'title' attribute
      ! try to find a "Physics_diagnostic" attribute (UK GCM outputs)
      text='   '
      ierr=NF_GET_ATT_TEXT(infid,tmpvarid,'Physics_diagnostic',text)
      if (ierr.eq.NF_NOERR) then ! found 'Physics_diagnostic' attribute
        write(*,*) "Found Physics_diagnostic ",trim(text)
        j=j+1
        ! write it as a 'long_name' attribute
        ierr=NF_PUT_ATT_TEXT(outfid,var_id(i),'long_name',len_trim(text),text)
        if (ierr.ne.NF_NOERR) then
          write(*,*) "Error failed to copy Physics_diagnostic attribute:",trim(text)
          stop
        endif
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
    
enddo ! of do=1,nbvar


!===============================================================================
! 3.4. Write dimensions (and 1D varaiables)
!===============================================================================
! Switch out of NetCDF define mode
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not switch out of define mode"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! Write longitude
ierr=NF_PUT_VAR_REAL(outfid,lon_varid,lon)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not write longitude data to output file"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! Write latitude
ierr=NF_PUT_VAR_REAL(outfid,lat_varid,lat)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not write latitude data to output file"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

! Write altitude
if (ztype.eq.1) then ! pressure vertical coordinate
  ierr=NF_PUT_VAR_REAL(outfid,alt_varid,plevel)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write altitude data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
else if (ztype.eq.2) then ! above areoid altitude
  ierr=NF_PUT_VAR_REAL(outfid,alt_varid,zareoid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write altitude data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
else if (ztype.eq.3) then ! above local surface
  ierr=NF_PUT_VAR_REAL(outfid,alt_varid,zsurface)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write altitude data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
else ! radial distance to center of Mars
  ierr=NF_PUT_VAR_REAL(outfid,alt_varid,zradius)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write altitude data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif

! Write sigma levels (or hybrid coordinates)
if (have_sigma) then
  ierr=NF_PUT_VAR_REAL(outfid,sigma_varid,sigma)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write sigma data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
else ! hybrid coordinates
  ierr=NF_PUT_VAR_REAL(outfid,aps_varid,aps)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write aps data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  ierr=NF_PUT_VAR_REAL(outfid,bps_varid,bps)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write bps data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif

! Write controle
if (ctllength .ne. 0) then
  ierr=NF_PUT_VAR_REAL(outfid,ctl_varid,ctl)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error: Could not write controle data to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif

! Write time
ierr=NF_PUT_VARA_REAL(outfid,time_varid,1,timelength,time)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not write Time data to output file"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

!===============================================================================
! 3.5 Write 3D variables
!===============================================================================

! Write surface pressure
corner(:)=1
edges(1)=lonlength
edges(2)=latlength
edges(3)=timelength
ierr=NF_PUT_VARA_REAL(outfid,ps_varid,corner(1:3),edges(1:3),ps)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Error: Could not write ps data to output file"
  write(*,*) NF_STRERROR(ierr)
  stop
endif

!===============================================================================
! 4. Interpolate and write 4D variables
!===============================================================================

! 4.0 Allocations
!indata() to store input (on GCM grid) data
allocate(indata(lonlength,latlength,altlength,timelength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate indata(lonlength,latlength,altlength,timelength)"
  write(*,*) "       lonlength=",lonlength," latlength=",latlength
  write(*,*) "       altlength=",altlength," timelength=",timelength
  stop
endif
! outdata() to store output (on new vertical grid) data
allocate(outdata(lonlength,latlength,nblev,timelength),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "Error: Failed to allocate outdata(lonlength,latlength,nblev,timelength)"
  write(*,*) "       lonlength=",lonlength," latlength=",latlength
  write(*,*) "       nblev=",nblev," timelength=",timelength
  stop
endif

! 4.1 If output is in pressure coordinate
if (ztype.eq.1) then
  do i=1,nbvar ! loop on 4D variable to process
  ! identify and read a dataset
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to get ID for variable ',var(i)
    stop
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,indata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to load variable ',var(i)
    stop
  endif
  
  ! interpolate dataset onto new grid
  call p_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,ps,press,indata,plevel,outdata)
  
  ! write new dataset to output file
  ierr=NF_PUT_VAR_REAL(outfid,var_id(i),outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*)'Error, Failed to write variable ',var(i)
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
  
  enddo ! of do i=1,nbvar
  
  ! interpolate zareoid onto new grid
  call p_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,ps,press,za_gcm,plevel,outdata)
  ! write result to output file
  corner(:)=1
  edges(1)=lonlength
  edges(2)=latlength
  edges(3)=nblev
  edges(4)=timelength
  ierr=NF_PUT_VARA_REAL(outfid,za_varid,corner,edges,outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error, Failed to write zareoid to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif ! of if (ztype.eq.1)

! 4.2 If output is in above areoid altitude
if (ztype.eq.2) then
  do i=1,nbvar ! loop on 4D variable to process
  ! identify and read a dataset
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to get ID for variable ',var(i)
    stop
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,indata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to load variable ',var(i)
    stop
  endif
  
  ! interpolate dataset onto new grid
  ! check if variable is "rho" (to set flag for interpolation below)
  if (var(i).eq.'rho') then
    j=1
  else
    j=0
  endif
  
  call z_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,za_gcm,indata,j,zareoid,outdata)

  ! write new dataset to output file
  ierr=NF_PUT_VAR_REAL(outfid,var_id(i),outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*)'Error, Failed to write variable ',var(i)
    write(*,*) NF_STRERROR(ierr)
    stop
  endif

  enddo ! of do i=1,nbvar
  
  ! interpolate pressure onto new grid
  write(*,*) ""
  write(*,*) "Processing pressure"
  call z_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,za_gcm,press,1,zareoid,outdata)

  ! write result to output file
  ierr=NF_PUT_VAR_REAL(outfid,p_varid,outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error, Failed to write pressure to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif ! of if (ztype.eq.2)

! 4.3 If output is in above local surface altitude
if (ztype.eq.3) then
  do i=1,nbvar ! loop on 4D variable to process
    write(*,*) " "
    write(*,*) "Processing "//trim(var(i))
    ! identify and read a dataset
    ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) 'Error, failed to get ID for variable ',var(i)
      stop
    endif
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,indata)
    if (ierr.ne.NF_NOERR) then
      write(*,*) 'Error, failed to load variable ',var(i)
      stop
    endif
  
    ! interpolate dataset onto new grid
    ! check if variable is "rho" (to set flag for interpolation below)
    if (var(i).eq.'rho') then
      j=1
    else
      j=0
    endif
  
    call zs_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                        miss_val,zs_gcm,indata,phisinit,press,temp,rho,&
                        have_rho,j,zsurface,outdata)

    ! write new dataset to output file
    ierr=NF_PUT_VAR_REAL(outfid,var_id(i),outdata)
    if (ierr.ne.NF_NOERR) then
      write(*,*)'Error, Failed to write variable ',var(i)
      write(*,*) NF_STRERROR(ierr)
      stop
    endif
  enddo ! of do i=1,nbvar
  
  ! interpolate pressure onto new grid
  write(*,*) ""
  write(*,*) "Processing pressure"
  call zs_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,zs_gcm,press,phisinit,press,temp,rho,&
                      have_rho,1,zsurface,outdata)

  ! write result to output file
  ierr=NF_PUT_VAR_REAL(outfid,p_varid,outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error, Failed to write pressure to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif

endif ! of if (ztype.eq.3)

! 4.4 If output is in radial distance to center of planet
if (ztype.eq.4) then
  do i=1,nbvar ! loop on 4D variable to process
  ! identify and read a dataset
  ierr=NF_INQ_VARID(infid,var(i),tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to get ID for variable ',var(i)
    stop
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,indata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) 'Error, failed to load variable ',var(i)
    stop
  endif
  
  ! interpolate dataset onto new grid
  ! check if variable is "rho" (to set flag for interpolation below)
  if (var(i).eq.'rho') then
    j=1
  else
    j=0
  endif
  
  call z_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,ra_gcm,indata,j,zradius,outdata)

  ! write new dataset to output file
  ierr=NF_PUT_VAR_REAL(outfid,var_id(i),outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*)'Error, Failed to write variable ',var(i)
    write(*,*) NF_STRERROR(ierr)
    stop
  endif

  enddo ! of do i=1,nbvar
  
  ! interpolate pressure onto new grid
  write(*,*) ""
  write(*,*) "Processing pressure"
  call z_coord_interp(lonlength,latlength,altlength,timelength,nblev, &
                      miss_val,ra_gcm,press,1,zradius,outdata)

  ! write result to output file
  ierr=NF_PUT_VAR_REAL(outfid,p_varid,outdata)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Error, Failed to write pressure to output file"
    write(*,*) NF_STRERROR(ierr)
    stop
  endif
endif ! of if (ztype.eq.4)

! 4.5 Close output file
ierr=NF_CLOSE(outfid)
if (ierr.ne.NF_NOERR) then
  write(*,*) 'Error, failed to close output file ',outfile
  write(*,*) NF_STRERROR(ierr)
endif

end program

!===============================================================================

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
Rmean=1000.*8.314511/controle(8) ! controle(8)=mugaz = molar mass (g.mol-1) of atm.

end subroutine init_planet_const

!===============================================================================

subroutine build_gcm_zs(lonlength,latlength,altlength,timelength, &
                         phis,ps,press,temp,rho,zs_gcm)
!==============================================================================
! Purpose: Integrate hydrostatic equation in order to build the "above local
!          surface altitudes" corresponding to GCM atmospheric levels
!==============================================================================
use planet_const, only : a0,g0
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlength ! # of points along longitude
integer,intent(in) :: latlength ! # of points along latitude
integer,intent(in) :: altlength ! # of points along altitude
integer,intent(in) :: timelength ! # of stored time steps
real,dimension(lonlength,latlength),intent(in) :: phis
! phis(:,:) is the ground geopotential
real,dimension(lonlength,latlength,timelength),intent(in) :: ps
! ps(:,:) is the surface pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: press
! press(:,:,:,:) is atmospheric pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: temp
! temp(:,:,:,:) is atmospheric temperature
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: rho
! rho(:,:,:,:) is atmospheric density

real,dimension(lonlength,latlength,altlength,timelength),intent(out) :: zs_gcm
! zs_gcm(:,:,:,:) are the above local surface altitudes of GCM levels

!===============================================================================
! Local variables:
!===============================================================================
real,dimension(:),allocatable :: sigma ! GCM sigma levels
real,dimension(:,:,:,:),allocatable :: R ! molecular gas constant
!real,dimension(:,:,:,:),allocatable :: zlocal ! local above surface altitude
real :: Rmean ! mean value of R for a given level
real :: Tmean ! "mean" value of temperature for a given level
integer iloop,jloop,kloop,tloop

! Parameters needed to integrate hydrostatic equation:
!real,parameter :: g0=3.7257964
!g0: exact mean gravity at radius=3396.km (Lemoine et al. 2001)
!real,parameter :: a0=3396.E3
!a0: 'mean' Mars radius=3396.km
real :: gz 
! gz: gravitational acceleration at a given (zareoid) altitude

!===============================================================================
! 1. Various initialisations
!===============================================================================
allocate(sigma(altlength))
allocate(R(lonlength,latlength,altlength,timelength))
!allocate(zlocal(lonlength,latlength,altlength,timelength))

!==============================================================================
! 2. Compute Molecular Gas Constant (J kg-1 K-1) at every grid point 
!   --later used to integrate the hydrostatic equation--
!==============================================================================
do tloop=1,timelength
  do kloop=1,altlength
    do jloop=1,latlength
      do iloop=1,lonlength
        R(iloop,jloop,kloop,tloop)=press(iloop,jloop,kloop,tloop)/ &
                                            (rho(iloop,jloop,kloop,tloop)* &
                                             temp(iloop,jloop,kloop,tloop))
      enddo
    enddo
  enddo
enddo

!===============================================================================
! 3. Integrate hydrostatic equation to compute zlocal and za_gcm
!===============================================================================
do tloop=1,timelength
  do jloop=1,latlength
    do iloop=1,lonlength
      ! handle case of first altitude level
      sigma(1)=press(iloop,jloop,1,tloop)/ps(iloop,jloop,tloop)
      zs_gcm(iloop,jloop,1,tloop)=-log(sigma(1))*R(iloop,jloop,1,tloop)* &
                                                 temp(iloop,jloop,1,tloop)/g0
!      za_gcm(iloop,jloop,1,tloop)=zlocal(iloop,jloop,1,tloop)+ &
!                                  phis(iloop,jloop)/g0
      do kloop=2,altlength
        ! compute sigma level of layer
        sigma(kloop)=press(iloop,jloop,kloop,tloop)/ps(iloop,jloop,tloop)
        
        ! compute "mean" temperature of the layer
        if(temp(iloop,jloop,kloop,tloop).eq.  &
           temp(iloop,jloop,kloop-1,tloop)) then
          Tmean=temp(iloop,jloop,kloop,tloop)
        else
          Tmean=(temp(iloop,jloop,kloop,tloop)-      &
                 temp(iloop,jloop,kloop-1,tloop))/   &
                 log(temp(iloop,jloop,kloop,tloop)/  &
                     temp(iloop,jloop,kloop-1,tloop))
        endif
        
        ! compute mean value of R of the layer
        Rmean=0.5*(R(iloop,jloop,kloop,tloop)+R(iloop,jloop,kloop-1,tloop))
        
        ! compute gravitational acceleration (at altitude zaeroid(kloop-1))
        ! NB: zareoid=zsurface+phis/g0
        gz=g0*(a0**2)/ &
           (a0+zs_gcm(iloop,jloop,kloop-1,tloop)+phis(iloop,jloop)/g0)**2
        
        ! compute zs_gcm(iloop,jloop,kloop,tloop)
        zs_gcm(iloop,jloop,kloop,tloop)=zs_gcm(iloop,jloop,kloop-1,tloop)- &
                log(sigma(kloop)/sigma(kloop-1))*Rmean*Tmean/gz
        
        ! compute za_gcm(kloop)
!        za_gcm(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop,tloop)+ &
!                                        phis(iloop,jloop)/g0
      enddo ! kloop
    enddo ! iloop
  enddo ! jloop
enddo ! tloop

! Cleanup
deallocate(sigma)
deallocate(R)
!deallocate(zlocal)

end subroutine build_gcm_zs

!===============================================================================

subroutine build_gcm_za(lonlength,latlength,altlength,timelength, &
                         phis,ps,press,temp,rho,za_gcm)
!==============================================================================
! Purpose: Integrate hydrostatic equation in order to build the "above areoid
!          altitudes" corresponding to GCM atmospheric levels
!==============================================================================
use planet_const, only : a0,g0
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlength ! # of points along longitude
integer,intent(in) :: latlength ! # of points along latitude
integer,intent(in) :: altlength ! # of points along altitude
integer,intent(in) :: timelength ! # of stored time steps
real,dimension(lonlength,latlength),intent(in) :: phis
! phis(:,:) is the ground geopotential
real,dimension(lonlength,latlength,timelength),intent(in) :: ps
! ps(:,:) is the surface pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: press
! press(:,:,:,:) is atmospheric pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: temp
! temp(:,:,:,:) is atmospheric temperature
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: rho
! rho(:,:,:,:) is atmospheric density

real,dimension(lonlength,latlength,altlength,timelength),intent(out) :: za_gcm
! za_gcm(:,:,:,:) are the above aroid heights of GCM levels

!===============================================================================
! Local variables:
!===============================================================================
real,dimension(:),allocatable :: sigma ! GCM sigma levels
real,dimension(:,:,:,:),allocatable :: R ! molecular gas constant
real,dimension(:,:,:,:),allocatable :: zlocal ! local above surface altitude
real :: Rmean ! mean value of R for a given level
real :: Tmean ! "mean" value of temperature for a given level
integer iloop,jloop,kloop,tloop

! Parameters needed to integrate hydrostatic equation:
!real,parameter :: g0=3.7257964
!g0: exact mean gravity at radius=3396.km (Lemoine et al. 2001)
!real,parameter :: a0=3396.E3
!a0: 'mean' Mars radius=3396.km
real :: gz 
! gz: gravitational acceleration at a given (zareoid) altitude

!===============================================================================
! 1. Various initialisations
!===============================================================================
allocate(sigma(altlength))
allocate(R(lonlength,latlength,altlength,timelength))
allocate(zlocal(lonlength,latlength,altlength,timelength))

!==============================================================================
! 2. Compute Molecular Gas Constant (J kg-1 K-1) at every grid point 
!   --later used to integrate the hydrostatic equation--
!==============================================================================
do tloop=1,timelength
  do kloop=1,altlength
    do jloop=1,latlength
      do iloop=1,lonlength
        R(iloop,jloop,kloop,tloop)=press(iloop,jloop,kloop,tloop)/ &
                                            (rho(iloop,jloop,kloop,tloop)* &
                                             temp(iloop,jloop,kloop,tloop))
      enddo
    enddo
  enddo
enddo

!===============================================================================
! 3. Integrate hydrostatic equation to compute zlocal and za_gcm
!===============================================================================
do tloop=1,timelength
  do jloop=1,latlength
    do iloop=1,lonlength
      ! handle case of first altitude level
      sigma(1)=press(iloop,jloop,1,tloop)/ps(iloop,jloop,tloop)
      zlocal(iloop,jloop,1,tloop)=-log(sigma(1))*R(iloop,jloop,1,tloop)* &
                                                 temp(iloop,jloop,1,tloop)/g0
      za_gcm(iloop,jloop,1,tloop)=zlocal(iloop,jloop,1,tloop)+ &
                                  phis(iloop,jloop)/g0
      do kloop=2,altlength
        ! compute sigma level of layer
        sigma(kloop)=press(iloop,jloop,kloop,tloop)/ps(iloop,jloop,tloop)
        
        ! compute "mean" temperature of the layer
        if(temp(iloop,jloop,kloop,tloop).eq.  &
           temp(iloop,jloop,kloop-1,tloop)) then
          Tmean=temp(iloop,jloop,kloop,tloop)
        else
          Tmean=(temp(iloop,jloop,kloop,tloop)-      &
                 temp(iloop,jloop,kloop-1,tloop))/   &
                 log(temp(iloop,jloop,kloop,tloop)/  &
                     temp(iloop,jloop,kloop-1,tloop))
        endif
        
        ! compute mean value of R of the layer
        Rmean=0.5*(R(iloop,jloop,kloop,tloop)+R(iloop,jloop,kloop-1,tloop))
        
        ! compute gravitational acceleration (at altitude zaeroid(kloop-1))
        gz=g0*(a0**2)/(a0+za_gcm(iloop,jloop,kloop-1,tloop))**2
        
        ! compute zlocal(kloop)
        zlocal(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop-1,tloop)- &
                log(sigma(kloop)/sigma(kloop-1))*Rmean*Tmean/gz
        
        ! compute za_gcm(kloop)
        za_gcm(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop,tloop)+ &
                                        phis(iloop,jloop)/g0
      enddo ! kloop
    enddo ! iloop
  enddo ! jloop
enddo ! tloop

! Cleanup
deallocate(sigma)
deallocate(R)
deallocate(zlocal)

end subroutine build_gcm_za

!===============================================================================

subroutine build_zs(altlength,have_sigma,sigma,aps,bps,zsurface)
!==============================================================================
! Build generic above surface altitudes, using either sigma levels
! or hybrid vertical coordinates.
! In order to do so, we need to use scale heights that vary with altitude.
! The scale heights distribution is then set to vary from a min. to a max.
! following a tanh() law over an imposed transition range (see below).
! Similarily, a reference mean surface pressure (610 Pa) is used to
! compute generic above surface altitudes.
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: altlength ! # of points along altitude
logical,intent(in) :: have_sigma ! true if sigma(:) are known
! have_sigma is false if aps() and bps(:) are given instead
real,dimension(altlength),intent(in) :: sigma ! sigma levels
real,dimension(altlength),intent(in) :: aps ! hybrid pressure levels
real,dimension(altlength),intent(in) :: bps ! hybrid sigma levels
real,dimension(altlength),intent(out) :: zsurface ! altitudes (m)

!===============================================================================
! Local variables:
!===============================================================================
real,dimension(:),allocatable :: H ! scale heights
real,parameter :: H_low=9650   ! scale height at low altitudes
real,parameter :: H_high=15000 ! scale height at high altitudes
real,parameter :: trans_window=10 ! # of levels over which H(:) goes
                                  ! from H_low to H_high
real,parameter :: lev_trans=32+trans_window/2
! lev_trans: level at which H(lev_trans)=(H_low+H_high)/2
! N.B.: keep lev_trans and lev_trans as reals to avoid truncation issues

real,parameter :: P_ref=610 ! reference surface pressure used to build zsurface
integer i 

! 1. Build scale heights
allocate(H(altlength))
do i=1,altlength
  H(i)=H_low+(H_high-H_low)*0.5*(1.0+tanh(6.*(i-lev_trans)/trans_window))
enddo

! 2. Compute zsurface(:)
if (have_sigma) then ! use sigma levels
  do i=1,altlength
    zsurface(i)=-H(i)*log(sigma(i)*P_ref)
  enddo
else ! use hybrid coordinates
  do i=1,altlength
    zsurface(i)=-H(i)*log((aps(i)/P_ref)+bps(i))
  enddo
endif

! Cleanup
deallocate(H)

end subroutine build_zs

!===============================================================================

subroutine crude_gcm_zs(lonlength,latlength,altlength,timelength, &
                         phis,ps,press,temp,zs_gcm)
!==============================================================================
! Purpose: Integrate hydrostatic equation in order to build the "above local
!          surface altitudes" corresponding to GCM atmospheric levels
!==============================================================================
use planet_const, only : a0,g0,Rmean
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlength ! # of points along longitude
integer,intent(in) :: latlength ! # of points along latitude
integer,intent(in) :: altlength ! # of points along altitude
integer,intent(in) :: timelength ! # of stored time steps
real,dimension(lonlength,latlength),intent(in) :: phis
! phis(:,:) is the ground geopotential
real,dimension(lonlength,latlength,timelength),intent(in) :: ps
! ps(:,:) is the surface pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: press
! press(:,:,:,:) is atmospheric pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: temp
! temp(:,:,:,:) is atmospheric temperature

real,dimension(lonlength,latlength,altlength,timelength),intent(out) :: zs_gcm
! zs_gcm(:,:,:,:) are the above local surface altitudes of GCM levels

!===============================================================================
! Local variables:
!===============================================================================
real,dimension(:),allocatable :: sigma ! GCM sigma levels
!real,parameter :: R=191 ! molecular gas constant
!real,dimension(:,:,:,:),allocatable :: zlocal ! local above surface altitude
real :: Tmean ! "mean" value of temperature for a given level
integer iloop,jloop,kloop,tloop

! Parameters needed to integrate hydrostatic equation:
!real,parameter :: g0=3.7257964
!g0: exact mean gravity at radius=3396.km (Lemoine et al. 2001)
!real,parameter :: a0=3396.E3
!a0: 'mean' Mars radius=3396.km
real :: gz 
! gz: gravitational acceleration at a given (zareoid) altitude

!===============================================================================
! 1. Various initialisations
!===============================================================================
allocate(sigma(altlength))
!allocate(zlocal(lonlength,latlength,altlength,timelength))

!===============================================================================
! 2. Integrate hydrostatic equation to compute zlocal and za_gcm
!===============================================================================
do tloop=1,timelength
  do jloop=1,latlength
    do iloop=1,lonlength
      ! handle case of first altitude level
      sigma(1)=press(iloop,jloop,1,tloop)/ps(iloop,jloop,tloop)
      zs_gcm(iloop,jloop,1,tloop)=-log(sigma(1))* &
                                           Rmean*temp(iloop,jloop,1,tloop)/g0
!      za_gcm(iloop,jloop,1,tloop)=zlocal(iloop,jloop,1,tloop)+ &
!                                  phis(iloop,jloop)/g0
      do kloop=2,altlength
        ! compute sigma level of layer
        sigma(kloop)=press(iloop,jloop,kloop,tloop)/ps(iloop,jloop,tloop)
        
        ! compute "mean" temperature of the layer
        if(temp(iloop,jloop,kloop,tloop).eq.  &
           temp(iloop,jloop,kloop-1,tloop)) then
          Tmean=temp(iloop,jloop,kloop,tloop)
        else
          Tmean=(temp(iloop,jloop,kloop,tloop)-      &
                 temp(iloop,jloop,kloop-1,tloop))/   &
                 log(temp(iloop,jloop,kloop,tloop)/  &
                     temp(iloop,jloop,kloop-1,tloop))
        endif
                
        ! compute gravitational acceleration (at altitude zaeroid(kloop-1))
        ! NB: zareoid=zsurface+phis/g0
        gz=g0*(a0**2)/ &
           (a0+zs_gcm(iloop,jloop,kloop-1,tloop)+phis(iloop,jloop)/g0)**2
        
        ! compute zs_gcm(iloop,jloop,kloop,tloop)
        zs_gcm(iloop,jloop,kloop,tloop)=zs_gcm(iloop,jloop,kloop-1,tloop)- &
                log(sigma(kloop)/sigma(kloop-1))*Rmean*Tmean/gz
        
        ! compute za_gcm(kloop)
!        za_gcm(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop,tloop)+ &
!                                        phis(iloop,jloop)/g0
      enddo ! kloop
    enddo ! iloop
  enddo ! jloop
enddo ! tloop

! Cleanup
deallocate(sigma)
!deallocate(zlocal)

end subroutine crude_gcm_zs

!===============================================================================

subroutine crude_gcm_za(lonlength,latlength,altlength,timelength, &
                         phis,ps,press,temp,za_gcm)
!==============================================================================
! Purpose: Integrate hydrostatic equation in order to build the "above areoid
!          altitudes" corresponding to GCM atmospheric levels
!==============================================================================
use planet_const, only : a0,g0,Rmean
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlength ! # of points along longitude
integer,intent(in) :: latlength ! # of points along latitude
integer,intent(in) :: altlength ! # of points along altitude
integer,intent(in) :: timelength ! # of stored time steps
real,dimension(lonlength,latlength),intent(in) :: phis
! phis(:,:) is the ground geopotential
real,dimension(lonlength,latlength,timelength),intent(in) :: ps
! ps(:,:) is the surface pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: press
! press(:,:,:,:) is atmospheric pressure
real,dimension(lonlength,latlength,altlength,timelength),intent(in) :: temp
! temp(:,:,:,:) is atmospheric temperature

real,dimension(lonlength,latlength,altlength,timelength),intent(out) :: za_gcm
! za_gcm(:,:,:,:) are the above aroid heights of GCM levels

!===============================================================================
! Local variables:
!===============================================================================
real,dimension(:),allocatable :: sigma ! GCM sigma levels
!real,parameter :: R=191 ! molecular gas constant
real,dimension(:,:,:,:),allocatable :: zlocal ! local above surface altitude
real :: Tmean ! "mean" value of temperature for a given level
integer iloop,jloop,kloop,tloop

! Parameters needed to integrate hydrostatic equation:
!real,parameter :: g0=3.7257964
!g0: exact mean gravity at radius=3396.km (Lemoine et al. 2001)
!real,parameter :: a0=3396.E3
!a0: 'mean' Mars radius=3396.km
real :: gz 
! gz: gravitational acceleration at a given (zareoid) altitude

!===============================================================================
! 1. Various initialisations
!===============================================================================
allocate(sigma(altlength))
allocate(zlocal(lonlength,latlength,altlength,timelength))

!===============================================================================
! 2. Integrate hydrostatic equation to compute zlocal and za_gcm
!===============================================================================
do tloop=1,timelength
  do jloop=1,latlength
    do iloop=1,lonlength
      ! handle case of first altitude level
      sigma(1)=press(iloop,jloop,1,tloop)/ps(iloop,jloop,tloop)
      zlocal(iloop,jloop,1,tloop)=-log(sigma(1))* &
                                            Rmean*temp(iloop,jloop,1,tloop)/g0
      za_gcm(iloop,jloop,1,tloop)=zlocal(iloop,jloop,1,tloop)+ &
                                  phis(iloop,jloop)/g0
      do kloop=2,altlength
        ! compute sigma level of layer
        sigma(kloop)=press(iloop,jloop,kloop,tloop)/ps(iloop,jloop,tloop)
        
        ! compute "mean" temperature of the layer
        if(temp(iloop,jloop,kloop,tloop).eq.  &
           temp(iloop,jloop,kloop-1,tloop)) then
          Tmean=temp(iloop,jloop,kloop,tloop)
        else
          Tmean=(temp(iloop,jloop,kloop,tloop)-      &
                 temp(iloop,jloop,kloop-1,tloop))/   &
                 log(temp(iloop,jloop,kloop,tloop)/  &
                     temp(iloop,jloop,kloop-1,tloop))
        endif
                
        ! compute gravitational acceleration (at altitude zaeroid(kloop-1))
        gz=g0*(a0**2)/(a0+za_gcm(iloop,jloop,kloop-1,tloop))**2
        
        ! compute zlocal(kloop)
        zlocal(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop-1,tloop)- &
                log(sigma(kloop)/sigma(kloop-1))*Rmean*Tmean/gz
        
        ! compute za_gcm(kloop)
        za_gcm(iloop,jloop,kloop,tloop)=zlocal(iloop,jloop,kloop,tloop)+ &
                                        phis(iloop,jloop)/g0
      enddo ! kloop
    enddo ! iloop
  enddo ! jloop
enddo ! tloop

! Cleanup
deallocate(sigma)
deallocate(zlocal)

end subroutine crude_gcm_za

!===============================================================================

subroutine p_coord_interp(lonlen,latlen,altlen,tlen,newaltlen, &
                          missing,ps,press,gcmdata,plevels,newdata)
!==============================================================================
! Purpose:
! Recast a 4D (spatio-temporal) GCM dataset, which has a vertical coordinate
! in pseudo-altitude, into a dataset which has a vertical coordinate at given
! pressure levels
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlen ! # of points along longitude (GCM dataset)
integer,intent(in) :: latlen ! # of points along latitude (GCM dataset)
integer,intent(in) :: altlen ! # of points along altitude (GCM dataset)
integer,intent(in) :: tlen   ! # of stored time steps (GCM dataset)
integer,intent(in) :: newaltlen ! # of points along altitude
real,intent(in) :: missing ! missing value
real,dimension(lonlen,latlen,tlen),intent(in) :: ps ! GCM surface pressure
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: press ! GCM pressure
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: gcmdata ! GCM dataset
real,dimension(newaltlen),intent(in) :: plevels
! plevels(:) pressure levels at which gcmdata is to be interpolated

real,dimension(lonlen,latlen,newaltlen,tlen),intent(out) :: newdata
! newdata(:,:,:,:) gcmdata recasted along vertical coordinate
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: lnp
! lnp(:): used to store log(pressure) values
real,dimension(:),allocatable :: q
! q(:): used to store values along vertical direction
real :: x,y
! x,y: temporary variables
integer :: iloop,jloop,kloop,tloop

allocate(lnp(altlen))
allocate(q(altlen))

do tloop=1,tlen
  do jloop=1,latlen
    do iloop=1,lonlen
      ! build lnp() and q() along vertical coordinate
      do kloop=1,altlen
        lnp(kloop)=-log(press(iloop,jloop,kloop,tloop))
        q(kloop)=gcmdata(iloop,jloop,kloop,tloop)
      enddo
      
      ! Interpolation  
      do kloop=1,newaltlen
        ! compute, by interpolation, value at pressure level plevels(kloop)
        if ((plevels(kloop).le.ps(iloop,jloop,tloop)).and. &
            (plevels(kloop).ge.press(iloop,jloop,altlen,tloop))) then
          x=-log(plevels(kloop))
          call interpolf(x,y,missing,lnp,q,altlen)
          newdata(iloop,jloop,kloop,tloop) = y
        else ! if plevels(kloop) is out of range,
          ! assign a "missing_value" at this node
          newdata(iloop,jloop,kloop,tloop) = missing
        endif
      enddo
      
    enddo !iloop
  enddo !jloop
enddo !tloop

! Cleanup
deallocate(lnp)
deallocate(q)

end subroutine p_coord_interp

!===============================================================================

subroutine z_coord_interp(lonlen,latlen,altlen,tlen,newaltlen, &
                          missing,z_gcm,gcmdata,flag,z_new,newdata)
!==============================================================================
! Purpose:
! Recast a 4D (spatio-temporal) GCM dataset 'gcmdata', for which corresponding
! grid values of altitude are known 'z_gcm', onto a new altitude grid 'z_new'.
! "Altitudes" can be above areoid or above local surface altitudes, as long as
! both 'z_gcm' and 'znew' refer to altitudes of a same type.
! Note: If altitudes in 'znew' fall outside of the range of altitudes
!       in 'z_gcm' then corresponding 'newdata' value is set to 'missing'.
!       'z_gcm' and 'znew' altitudes must be monotone increasing sequences.
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlen ! # of points along longitude (GCM dataset)
integer,intent(in) :: latlen ! # of points along latitude (GCM dataset)
integer,intent(in) :: altlen ! # of points along altitude (GCM dataset)
integer,intent(in) :: tlen   ! # of stored time steps (GCM dataset)
integer,intent(in) :: newaltlen ! # of points along altitude
real,intent(in) :: missing ! missing value
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: z_gcm
!z_gcm(:,:,:,) GCM grid points altitude (above areoid or above surface)
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: gcmdata ! GCM dataset
integer,intent(in) :: flag ! flag (==1 for 'log' interpolation)
! flag==0 (standard linear interpolation)
real,dimension(newaltlen),intent(in) :: z_new
! z_new(:) altitudes (above areoid or surface) at which data must be recast

real,dimension(lonlen,latlen,newaltlen,tlen),intent(out) :: newdata
! newdata(:,:,:,:) gcmdata recasted along vertical coordinate

!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: za,q
real,dimension(:),allocatable :: logq
! za(:): used to store z_gcm at a given location
! q(:): used to store field values along the vertical direction
real :: x,y ! temporary variables
integer :: iloop,jloop,kloop,tloop

allocate(za(altlen))
allocate(q(altlen))
allocate(logq(altlen))

if (flag.eq.0) then
 do tloop=1,tlen
  do jloop=1,latlen
    do iloop=1,lonlen
      do kloop=1,altlen
        ! extract the vertical coordinates
        za(kloop)=z_gcm(iloop,jloop,kloop,tloop)
        ! store values along altitude
        q(kloop)=gcmdata(iloop,jloop,kloop,tloop)
      enddo !kloop
      
      ! Interpolation
      do kloop=1,newaltlen
      ! Check if z_new(kloop) is within range
      if ((z_new(kloop).ge.z_gcm(iloop,jloop,1,tloop)).and. &
          (z_new(kloop).le.z_gcm(iloop,jloop,altlen,tloop))) then
        ! z_new(kloop) is within range
        x=z_new(kloop)
        call interpolf(x,y,missing,za,q,altlen)
        newdata(iloop,jloop,kloop,tloop)=y
      else ! z_new(kloop) is out of range
        newdata(iloop,jloop,kloop,tloop)=missing
      endif 
      enddo !kloop
    enddo !iloop
  enddo !jloop
 enddo !tloop 
else ! case when flag=1 (i.e.: rho)
 do tloop=1,tlen
  do jloop=1,latlen
    do iloop=1,lonlen
      do kloop=1,altlen
        ! extract the vertical coordinates 
        za(kloop)=z_gcm(iloop,jloop,kloop,tloop)
        ! store log values along altitude
        logq(kloop)=log(gcmdata(iloop,jloop,kloop,tloop))
      enddo !kloop
      
      ! Interpolation
      do kloop=1,newaltlen
      ! Check if z_new(kloop) is within range
      if ((z_new(kloop).ge.z_gcm(iloop,jloop,1,tloop)).and. &
          (z_new(kloop).le.z_gcm(iloop,jloop,altlen,tloop))) then
        ! z_new(kloop) is within range
        x=z_new(kloop)
        call interpolf(x,y,missing,za,logq,altlen)
        newdata(iloop,jloop,kloop,tloop)=exp(y)
      else  ! z_new(kloop) is out of range
        newdata(iloop,jloop,kloop,tloop)=missing
      endif
      enddo !kloop
    enddo !iloop
  enddo !jloop
 enddo !tloop 
endif

! Cleanup
deallocate(za)
deallocate(q)
deallocate(logq)

end subroutine z_coord_interp

!===============================================================================

subroutine zs_coord_interp(lonlen,latlen,altlen,tlen,newaltlen, &
                          missing,z_gcm,gcmdata,phis,press,temp,rho, &
                          have_rho,flag,z_new,newdata)
!==============================================================================
! Purpose:
! Recast a 4D (spatio-temporal) GCM dataset 'gcmdata', for which corresponding
! grid values of altitude are known 'z_gcm', onto a new altitude grid 'z_new'.
! "Altitudes" must be above local surface altitudes.
! Notes: 
!       'z_gcm' and 'znew' altitudes must be monotone increasing sequences.
!       If altitudes in 'znew(i)' fall below those in 'z_gcm(:,:,1,:)', then
!       'newdata(:,:,i,:)' is set to constant value 'gcmdata(:,:,1,:)' if
!       flag=0, and extrapolated (exponentially) if flag=1.
!       If altitudes in 'znew(i)' are above those in 'z_gcm(:,:,altlen,:)', then
!       'newdata(:,:,i,:)' is set to constant value 'gcmdata(:,:,altlen,:)' 
!       if flag=0, and extrapolated (exponentially) if flag=1.
!==============================================================================
use planet_const, only : a0,g0,Rmean
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlen ! # of points along longitude (GCM dataset)
integer,intent(in) :: latlen ! # of points along latitude (GCM dataset)
integer,intent(in) :: altlen ! # of points along altitude (GCM dataset)
integer,intent(in) :: tlen   ! # of stored time steps (GCM dataset)
integer,intent(in) :: newaltlen ! # of points along altitude
real,intent(in) :: missing ! missing value
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: z_gcm
!z_gcm(:,:,:,) GCM grid points altitude (above areoid or above surface)
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: gcmdata ! GCM dataset
real,dimension(lonlen,latlen),intent(in) :: phis
! phis(:,:) is the ground geopotential
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: press
! press(:,:,:,:) is atmospheric pressure on GCM levels
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: temp
! temp(:,:,:,:) is atmospheric temperature on GCM levels
real,dimension(lonlen,latlen,altlen,tlen),intent(in) :: rho
! rho(:,:,:,:) is atmospheric density on GCM levels
logical,intent(in) :: have_rho ! trueif we have density at hand
integer,intent(in) :: flag ! flag (==1 for 'log' interpolation)
! flag==0 (standard linear interpolation)
real,dimension(newaltlen),intent(in) :: z_new
! z_new(:) altitudes (above areoid or surface) at which data must be recast

real,dimension(lonlen,latlen,newaltlen,tlen),intent(out) :: newdata
! newdata(:,:,:,:) gcmdata recasted along vertical coordinate

!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: z,q
real,dimension(:),allocatable :: logq
! z(:): used to store z_gcm at a given location
! q(:): used to store field values along the vertical direction
real,dimension(:,:,:),allocatable :: Rbottom,Rtop
! values of R (gas constant) below and above GCM layers
real :: x,y ! temporary variables
integer :: iloop,jloop,kloop,tloop
!real,parameter :: g0=3.7257964
!g0: exact mean gravity at radius=3396.km (Lemoine et al. 2001)
!real,parameter :: a0=3396.E3
!a0: 'mean' Mars radius=3396.km
!real,parameter :: Rmean=191 ! molecular gas constant
real :: gz 
! gz: gravitational acceleration at a given (above areoid) altitude

allocate(z(altlen))
allocate(q(altlen))
allocate(logq(altlen))

! 1. Build Rbottom and Rtop (only necessary if flag=1)
if (flag.eq.1) then
  allocate(Rbottom(lonlen,latlen,tlen))
  allocate(Rtop(lonlen,latlen,tlen))
  if (have_rho) then
   do iloop=1,lonlen
    do jloop=1,latlen
      do tloop=1,tlen
        Rbottom(iloop,jloop,tloop)=press(iloop,jloop,1,tloop)/ &
                                            (rho(iloop,jloop,1,tloop)* &
                                             temp(iloop,jloop,1,tloop))
        Rtop(iloop,jloop,tloop)=press(iloop,jloop,altlen,tloop)/ &
                                            (rho(iloop,jloop,altlen,tloop)* &
                                             temp(iloop,jloop,altlen,tloop))
      enddo
    enddo
   enddo
  else ! we don't have density at hand; use mean molecular gas constant value
    Rbottom(:,:,:)=Rmean
    Rtop(:,:,:)=Rmean
  endif
endif ! if (flag.eq.1)

! 2. Interpolation
if (flag.eq.0) then
 do tloop=1,tlen
  do jloop=1,latlen
    do iloop=1,lonlen
      ! preliminary stuff
      do kloop=1,altlen
        ! extract the vertical coordinates
        z(kloop)=z_gcm(iloop,jloop,kloop,tloop)
        ! store values along altitude
        q(kloop)=gcmdata(iloop,jloop,kloop,tloop)
      enddo !kloop
      
      ! Interpolation
      do kloop=1,newaltlen
        if (z_new(kloop).lt.z_gcm(iloop,jloop,1,tloop)) then
          ! z_new(kloop) is below lowest GCM level
          newdata(iloop,jloop,kloop,tloop)=gcmdata(iloop,jloop,1,tloop)
        else if (z_new(kloop).gt.z_gcm(iloop,jloop,altlen,tloop)) then
          ! z_new(kloop) is above highest GCM level
          newdata(iloop,jloop,kloop,tloop)=gcmdata(iloop,jloop,altlen,tloop)
        else ! z_new(kloop) is within range
          x=z_new(kloop)
          call interpolf(x,y,missing,z,q,altlen)
          newdata(iloop,jloop,kloop,tloop)=y
        endif
      enddo !kloop
    enddo !iloop
  enddo !jloop
 enddo !tloop 
else ! case when flag=1 (i.e.: rho or pressure)
 do tloop=1,tlen
  do jloop=1,latlen
    do iloop=1,lonlen
      ! preliminary stuff
      do kloop=1,altlen
        ! extract the vertical coordinates 
        z(kloop)=z_gcm(iloop,jloop,kloop,tloop)
        ! store log values along altitude
        logq(kloop)=log(gcmdata(iloop,jloop,kloop,tloop))
      enddo !kloop
      
      ! Interpolation
      do kloop=1,newaltlen
        if (z_new(kloop).lt.z_gcm(iloop,jloop,1,tloop)) then
          ! z_new(kloop) is below lowest GCM level
          newdata(iloop,jloop,kloop,tloop)=gcmdata(iloop,jloop,1,tloop)* &
                   exp(((z_gcm(iloop,jloop,1,tloop)-z_new(kloop))*g0)/ &
                       (Rbottom(iloop,jloop,tloop)* &
                        temp(iloop,jloop,1,tloop)))
        else if (z_new(kloop).gt.z_gcm(iloop,jloop,altlen,tloop)) then
          ! z_new(kloop) is above highest GCM level
          ! NB: zareoid=zsurface+phis/g0
          gz=g0/(1.+((z_new(kloop)+phis(iloop,jloop)/g0)/a0))/(1+(z_gcm(iloop,jloop,altlen,tloop)+phis(iloop,jloop)/g0)/a0) !JL
          newdata(iloop,jloop,kloop,tloop)=gcmdata(iloop,jloop,altlen,tloop)* &
                   exp(((z_gcm(iloop,jloop,altlen,tloop)-z_new(kloop))*gz)/ &
                       (Rtop(iloop,jloop,tloop)* &
                        temp(iloop,jloop,altlen,tloop)))
        else ! z_new(kloop) is within range
          x=z_new(kloop)
          call interpolf(x,y,missing,z,logq,altlen)
          newdata(iloop,jloop,kloop,tloop)=exp(y)
        endif
      enddo !kloop
    enddo !iloop
  enddo !jloop
 enddo !tloop 
endif

! Cleanup
deallocate(z)
deallocate(q)
deallocate(logq)
if (flag.eq.1) then
  deallocate(Rbottom)
  deallocate(Rtop)
endif

end subroutine zs_coord_interp

!===============================================================================

subroutine interpolf(x,y,missing,xd,yd,nd)
!==============================================================================
! Purpose:
! Yield y=f(x), where xd() end yd() are arrays of known values,
! using linear interpolation
! If x is not included in the interval spaned by xd(), then y is set
! to a default value 'missing'
! Note:
! Array xd() should contain ordered (either increasing or decreasing) abscissas
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
real,intent(in) :: x ! abscissa at which interpolation is to be done
real,intent(in) :: missing ! missing value (if no interpolation is performed)
integer :: nd ! size of arrays
real,dimension(nd),intent(in) :: xd ! array of known absissas
real,dimension(nd),intent(in) :: yd ! array of correponding values

real,intent(out) :: y ! interpolated value
!==============================================================================
! Local variables:
!==============================================================================
integer :: i

! default: set y to 'missing'
y=missing

   do i=1,nd-1
     if (((x.ge.xd(i)).and.(x.le.xd(i+1))).or.&
          ((x.le.xd(i)).and.(x.ge.xd(i+1)))) then
        y=yd(i)+(x-xd(i))*(yd(i+1)-yd(i))/(xd(i+1)-xd(i))
        exit
     endif
   enddo


end subroutine interpolf

!===============================================================================

subroutine build_gcm_ra(lonlength,latlength,altlength,timelength, &
                         lon,lat,za_gcm,ra_gcm)
!==============================================================================
! Purpose:  Build, from the "above areoid altitudes" of GCM levels,
!           the corresponding "radial distance to center of Mars" values
!==============================================================================
implicit none
!==============================================================================
! Arguments:
!==============================================================================
integer,intent(in) :: lonlength ! # of points along longitude
integer,intent(in) :: latlength ! # of points along latitude
integer,intent(in) :: altlength ! # of points along altitude
integer,intent(in) :: timelength ! # of stored time steps
real,intent(in) :: lon(lonlength) ! GCM longitudes
real,intent(in) :: lat(latlength) ! GCM latitudes
!za_gcm: above areoid altitudes of GCM levels
real,intent(in) :: za_gcm(lonlength,latlength,altlength,timelength)
!ra_gcm: radial distance to center of planet of GCM levels
real,intent(out) :: ra_gcm(lonlength,latlength,altlength,timelength)

!===============================================================================
! Local variables:
!===============================================================================
real :: rareoid(lonlength,latlength) ! radial position (ie: distance to the
! center of Mars) of the reference areoid for GCM longitudes and latitudes
integer :: iloop,jloop,kloop,tloop
double precision :: dlon,dlat,dareoid ! double precision temp variables

! common data shared with 'areoid_ini' and 'geoid' routines
integer,parameter :: ndeg=90,ndeg2=2*ndeg,nd2p3=ndeg2+3
! data structure of gravity field
        double precision v0,omega,ae,gm
	double precision clm(0:ndeg,0:ndeg),slm(0:ndeg,0:ndeg)
	common /gmm1/v0,omega,ae,gm,clm,slm
        integer lmin,lmax
        double precision root(nd2p3)
        double precision requator
	common /sqr/ lmin,lmax,root,requator

!===============================================================================
! 1. Various initialisations
!===============================================================================
! initalize coefficients of areoid model
call areoid_ini()

! compute rareoid for GCM longitude/latitude grid
do iloop=1,lonlength
 dlon=lon(iloop) ! double precision version of lon
 do jloop=1,latlength
   dlat=lat(jloop) ! double precision version of lat
   call geoid(dlon,dlat,dareoid)
   rareoid(iloop,jloop)=dareoid
 enddo
enddo

!===============================================================================
! 1. Compute ra_gcm
!===============================================================================
!!! test:
!rareoid(:,:)=0

do tloop=1,timelength
  do kloop=1,altlength
    ra_gcm(:,:,kloop,tloop)=za_gcm(:,:,kloop,tloop)+rareoid(:,:)
  enddo
enddo

end subroutine build_gcm_ra



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines corresponding to Lemoine et al. 2001 implementation.
! Based on "areoid.f" program by G. Neumann, which is available at
! ftp://ltpftp.gsfc.nasa.gov/projects/tharsis/MOLA/SOFTWARE/
! uses the gravity field coefficients file 'mgm1025' available in
! the same ftp directory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine areoid_ini()
implicit none
!character(len=250) :: mgm="/u/emlmd/MCD_v4.3/data/mgm1025"

!character(len=64) :: title ! first line of file
!character(len=6) :: gcoef
! data structure of gravity field
         integer ndeg,ndeg2,nd2p3
	  parameter (ndeg=90,ndeg2=2*ndeg,nd2p3=ndeg2+3)
	  double precision c(0:ndeg,0:ndeg),s(0:ndeg,0:ndeg)
	  double precision  v0,omega,ae,gm !,clm,slm
	  double precision  plm(ndeg)
	  common /gmm1/ v0,omega,ae,gm,c,s
	  integer lmin,lmax
          double precision root(nd2p3)
          double precision r
          common /sqr/ lmin,lmax, root,r
!
!        integer lcmax,mmax
        integer k,l,m
!        double precision coef
        double precision x,xi,sum
        
	  ae= 3396000.d0
	  gm =42828.37d9
	  omega=0.70882181d-4
! this value makes the equatorial mean radius equal to 3396 km.
	do k=1,nd2p3
	 root(k)=sqrt(dble(k))
	enddo
!  initialize coefficients c and s
c(:,:)=0
s(:,:)=0

! fill in the values of c() and s() (extracted from mgm1025 file)
c( 2, 0)=   -0.87450441550858D-03
c( 3, 0)=   -0.11889205541698D-04
c( 4, 0)=    0.51231389950069D-05
c( 5, 0)=   -0.17246183479555D-05
c( 6, 0)=    0.13442413561338D-05
c( 7, 0)=    0.10568085417745D-05
c( 8, 0)=    0.14344394468535D-06
c( 9, 0)=   -0.28623030841535D-06
c(10, 0)=    0.72517371380309D-06
c(11, 0)=   -0.26515503657960D-06
c(12, 0)=    0.25791004374137D-06
c(13, 0)=   -0.48987732989967D-06
c(14, 0)=    0.30149567154029D-06
c(15, 0)=    0.49309755579596D-06
c(16, 0)=    0.62303904221200D-06
c(17, 0)=   -0.53940568680542D-07
c(18, 0)=   -0.31286409396935D-06
c(19, 0)=    0.66506982323697D-08
c(20, 0)=    0.29800506646171D-06
c(21, 0)=    0.12382554869535D-06
c(22, 0)=   -0.11079454891980D-06
c(23, 0)=    0.46651763615791D-07
c(24, 0)=    0.26561879988827D-06
c(25, 0)=   -0.14012911779311D-06
c(26, 0)=   -0.17890980799686D-06
c(27, 0)=   -0.22473934735542D-06
c(28, 0)=    0.35017911372558D-07
c(29, 0)=    0.96416933029589D-07
c(30, 0)=    0.13454167912854D-07
c(31, 0)=    0.11415389076391D-06
c(32, 0)=    0.14309722973488D-06
c(33, 0)=   -0.93613803420433D-07
c(34, 0)=   -0.11634597651051D-06
c(35, 0)=    0.17390073509388D-07
c(36, 0)=    0.17980990446829D-06
c(37, 0)=   -0.58813753768476D-07
c(38, 0)=   -0.15639701438663D-06
c(39, 0)=    0.67338858450205D-07
c(40, 0)=    0.10015015531020D-06
c(41, 0)=    0.65659089401047D-07
c(42, 0)=   -0.32828754909310D-07
c(43, 0)=    0.55320735747522D-07
c(44, 0)=    0.10719796171993D-06
c(45, 0)=   -0.24466821475646D-07
c(46, 0)=    0.53066022003621D-07
c(47, 0)=    0.18183459032733D-08
c(48, 0)=   -0.90336249640171D-07
c(49, 0)=    0.10172719809900D-07
c(50, 0)=    0.14157647949966D-06
c(51, 0)=   -0.13459597844273D-07
c(52, 0)=   -0.10770557778749D-06
c(53, 0)=    0.41656444281046D-07
c(54, 0)=    0.23491811928996D-07
c(55, 0)=    0.39514870435300D-07
c(56, 0)=    0.17750185036483D-07
c(57, 0)=   -0.16389790284882D-07
c(58, 0)=    0.15924966924563D-07
c(59, 0)=   -0.52269191539423D-07
c(60, 0)=    0.22651535424788D-07
c(61, 0)=   -0.22088875934791D-07
c(62, 0)=   -0.89365975311230D-07
c(63, 0)=   -0.22314891242664D-06
c(64, 0)=   -0.10126090956466D-06
c(65, 0)=    0.53472579669405D-07
c(66, 0)=   -0.79239790651784D-07
c(67, 0)=   -0.93403592836398D-07
c(68, 0)=   -0.92958885757065D-07
c(69, 0)=   -0.23268566991842D-06
c(70, 0)=   -0.61104272230081D-07
c(71, 0)=   -0.15227130308936D-06
c(72, 0)=    0.44022260088073D-07
c(73, 0)=   -0.18296622662986D-06
c(74, 0)=   -0.21563284613746D-07
c(75, 0)=   -0.14408830353913D-06
c(76, 0)=   -0.31971691137664D-07
c(77, 0)=   -0.96143919013694D-07
c(78, 0)=   -0.14039664820363D-07
c(79, 0)=   -0.56173045626134D-07
c(80, 0)=   -0.42587385931103D-07
c( 2, 1)=    0.94588820913966D-10
c( 3, 1)=    0.38007291643504D-05
c( 4, 1)=    0.42113282449449D-05
c( 5, 1)=    0.48295248273193D-06
c( 6, 1)=    0.17989300452431D-05
c( 7, 1)=    0.13721377780493D-05
c( 8, 1)=   -0.13172502269848D-06
c( 9, 1)=    0.41968715869939D-06
c(10, 1)=    0.92089843305934D-06
c(11, 1)=   -0.81438583356007D-06
c(12, 1)=   -0.11226169861948D-05
c(13, 1)=   -0.56756560374700D-06
c(14, 1)=    0.75140494090486D-06
c(15, 1)=    0.23398249720387D-06
c(16, 1)=   -0.23315134603590D-06
c(17, 1)=   -0.37162884465144D-06
c(18, 1)=   -0.18094168688736D-06
c(19, 1)=    0.46181639714625D-07
c(20, 1)=    0.18344606355547D-06
c(21, 1)=   -0.14791536512739D-06
c(22, 1)=   -0.40126140681542D-06
c(23, 1)=    0.24408804946295D-06
c(24, 1)=    0.18984659529516D-06
c(25, 1)=   -0.77644847446303D-07
c(26, 1)=   -0.18954371508813D-06
c(27, 1)=    0.58215212507928D-07
c(28, 1)=    0.15662107285173D-06
c(29, 1)=    0.69848912908713D-07
c(30, 1)=   -0.14594550925334D-06
c(31, 1)=   -0.60978169516080D-07
c(32, 1)=   -0.69857883539046D-07
c(33, 1)=   -0.40910599890025D-07
c(34, 1)=    0.69727112840257D-07
c(35, 1)=    0.69051821175597D-07
c(36, 1)=   -0.89123393751658D-07
c(37, 1)=   -0.10308616478520D-06
c(38, 1)=   -0.17953527802480D-07
c(39, 1)=    0.13357163900313D-06
c(40, 1)=    0.32997141656491D-07
c(41, 1)=   -0.42770439863940D-07
c(42, 1)=   -0.33815214882692D-07
c(43, 1)=    0.10228125797069D-06
c(44, 1)=    0.83199427365936D-07
c(45, 1)=   -0.26716533044488D-08
c(46, 1)=   -0.21968332747552D-07
c(47, 1)=    0.26051927992496D-07
c(48, 1)=   -0.20011329156204D-07
c(49, 1)=   -0.46000606601472D-07
c(50, 1)=   -0.50465401285307D-08
c(51, 1)=   -0.21435559836597D-07
c(52, 1)=   -0.79586940432080D-07
c(53, 1)=   -0.21721909216360D-08
c(54, 1)=    0.11211857004569D-07
c(55, 1)=   -0.14572008323947D-07
c(56, 1)=   -0.76165251339758D-07
c(57, 1)=   -0.33685497326864D-07
c(58, 1)=    0.29095996986081D-07
c(59, 1)=    0.11183726327720D-07
c(60, 1)=   -0.29984787445527D-07
c(61, 1)=    0.32521236160127D-07
c(62, 1)=    0.81714334232117D-08
c(63, 1)=   -0.44106457563710D-08
c(64, 1)=   -0.11522745649087D-07
c(65, 1)=    0.24452822525910D-07
c(66, 1)=    0.39582110377419D-08
c(67, 1)=   -0.16048574963929D-07
c(68, 1)=   -0.10635293597363D-08
c(69, 1)=   -0.25348757056472D-08
c(70, 1)=   -0.13990070916062D-07
c(71, 1)=   -0.21113532040349D-07
c(72, 1)=    0.21388767403742D-07
c(73, 1)=   -0.45295107975469D-07
c(74, 1)=    0.11061952017768D-08
c(75, 1)=   -0.20070701321975D-07
c(76, 1)=   -0.81623263977223D-08
c(77, 1)=   -0.48608889440242D-08
c(78, 1)=    0.48605157923197D-08
c(79, 1)=   -0.72265799247382D-08
c(80, 1)=    0.12959728064204D-07
c( 2, 2)=   -0.84585751018071D-04
c( 3, 2)=   -0.15933741663502D-04
c( 4, 2)=   -0.95107786023564D-06
c( 5, 2)=   -0.42918839549229D-05
c( 6, 2)=    0.85990806408567D-06
c( 7, 2)=    0.28084893201890D-05
c( 8, 2)=    0.18057628945403D-05
c( 9, 2)=    0.11356212264308D-05
c(10, 2)=    0.73405458338506D-08
c(11, 2)=   -0.31264422685096D-06
c(12, 2)=   -0.97454248928457D-07
c(13, 2)=    0.19140550616679D-06
c(14, 2)=    0.44883687385942D-06
c(15, 2)=   -0.35213615909880D-06
c(16, 2)=   -0.48079087360350D-06
c(17, 2)=   -0.20190086133281D-06
c(18, 2)=   -0.10432508924443D-06
c(19, 2)=   -0.52238899978893D-07
c(20, 2)=    0.69247343387462D-07
c(21, 2)=   -0.93994254289733D-07
c(22, 2)=    0.10295574560935D-06
c(23, 2)=    0.19547803118600D-06
c(24, 2)=    0.19407504646311D-07
c(25, 2)=   -0.17219187149635D-06
c(26, 2)=   -0.20135410709838D-06
c(27, 2)=   -0.28273103642271D-08
c(28, 2)=    0.30690987897186D-06
c(29, 2)=    0.79628272547828D-07
c(30, 2)=   -0.18900571342100D-06
c(31, 2)=    0.83337491764168D-07
c(32, 2)=    0.13315305225107D-06
c(33, 2)=    0.58675965618782D-07
c(34, 2)=   -0.10511750363253D-06
c(35, 2)=   -0.41561363195631D-07
c(36, 2)=    0.10803644254118D-06
c(37, 2)=    0.57900788447265D-07
c(38, 2)=   -0.10674339261125D-07
c(39, 2)=    0.10267822284623D-08
c(40, 2)=   -0.39409739046749D-07
c(41, 2)=    0.48560442200681D-07
c(42, 2)=    0.44295705013232D-08
c(43, 2)=   -0.23863651808965D-07
c(44, 2)=   -0.40529365741889D-07
c(45, 2)=   -0.72544669649420D-07
c(46, 2)=   -0.30840437495884D-07
c(47, 2)=   -0.10134355839673D-07
c(48, 2)=   -0.74066104086251D-07
c(49, 2)=   -0.42575092353871D-07
c(50, 2)=    0.38002915272366D-07
c(51, 2)=    0.30850031943561D-08
c(52, 2)=   -0.38896999892434D-07
c(53, 2)=   -0.20811133913029D-08
c(54, 2)=    0.60282320372300D-07
c(55, 2)=    0.66842630552045D-08
c(56, 2)=   -0.15054903643693D-07
c(57, 2)=   -0.16272429639032D-07
c(58, 2)=   -0.77062433134402D-08
c(59, 2)=    0.22152915027638D-07
c(60, 2)=    0.33331254807163D-07
c(61, 2)=    0.47393585702229D-07
c(62, 2)=    0.24590842039888D-07
c(63, 2)=   -0.25081905036699D-07
c(64, 2)=    0.85705073676948D-08
c(65, 2)=    0.60565718726133D-09
c(66, 2)=    0.13996463097736D-07
c(67, 2)=    0.36801218056836D-07
c(68, 2)=   -0.17804762229462D-07
c(69, 2)=    0.15279165045386D-07
c(70, 2)=    0.11881995727714D-07
c(71, 2)=   -0.23621731481701D-07
c(72, 2)=    0.16890385593438D-07
c(73, 2)=   -0.20283358613253D-07
c(74, 2)=    0.64980567456639D-08
c(75, 2)=   -0.23260084626746D-07
c(76, 2)=    0.99374517344100D-08
c(77, 2)=   -0.32363754554167D-07
c(78, 2)=    0.10856848087610D-07
c(79, 2)=   -0.27919914518472D-07
c(80, 2)=    0.16116790335365D-07
c( 3, 3)=    0.35023491711606D-04
c( 4, 3)=    0.64490745811941D-05
c( 5, 3)=    0.33076885084810D-05
c( 6, 3)=    0.95401012517447D-06
c( 7, 3)=    0.87863220403753D-06
c( 8, 3)=   -0.12040062428688D-05
c( 9, 3)=   -0.99094051068927D-06
c(10, 3)=   -0.29834384339389D-06
c(11, 3)=   -0.12986012000905D-05
c(12, 3)=   -0.14371536216267D-05
c(13, 3)=    0.21515088011757D-06
c(14, 3)=    0.42474467046937D-06
c(15, 3)=   -0.55907478206270D-06
c(16, 3)=   -0.53391027929080D-06
c(17, 3)=   -0.45940051707924D-07
c(18, 3)=   -0.29341890127158D-06
c(19, 3)=   -0.52823002208831D-07
c(20, 3)=    0.42381732878213D-06
c(21, 3)=   -0.15864628605041D-06
c(22, 3)=   -0.37694486993635D-06
c(23, 3)=    0.10153362064588D-06
c(24, 3)=    0.33326669401862D-06
c(25, 3)=    0.10923543811362D-06
c(26, 3)=   -0.15114438271846D-06
c(27, 3)=    0.12204330344988D-06
c(28, 3)=    0.19906682949500D-06
c(29, 3)=    0.10324442774105D-06
c(30, 3)=   -0.35165969895604D-07
c(31, 3)=   -0.37695255894017D-07
c(32, 3)=   -0.15441100504594D-06
c(33, 3)=   -0.12566971141268D-06
c(34, 3)=    0.86388483626716D-07
c(35, 3)=    0.73898333015289D-07
c(36, 3)=   -0.11654986577645D-06
c(37, 3)=   -0.12337255642019D-06
c(38, 3)=   -0.93452853607312D-07
c(39, 3)=    0.15866384390936D-07
c(40, 3)=    0.22621856963862D-07
c(41, 3)=   -0.82261088157716D-07
c(42, 3)=   -0.11509240581167D-06
c(43, 3)=   -0.51209031157573D-07
c(44, 3)=    0.87755823722094D-07
c(45, 3)=   -0.32736783199659D-07
c(46, 3)=   -0.44469896535062D-07
c(47, 3)=    0.53202827198090D-07
c(48, 3)=    0.31840012781999D-07
c(49, 3)=   -0.39439367510789D-07
c(50, 3)=   -0.24919100537213D-07
c(51, 3)=    0.21683500372067D-07
c(52, 3)=   -0.33140152676634D-08
c(53, 3)=   -0.18205948135764D-07
c(54, 3)=    0.43727885184963D-07
c(55, 3)=    0.14256829781278D-07
c(56, 3)=   -0.33010270260415D-07
c(57, 3)=    0.13439525530206D-07
c(58, 3)=    0.19680334210526D-07
c(59, 3)=    0.49353515402636D-07
c(60, 3)=   -0.32647840865208D-08
c(61, 3)=    0.81100619477578D-08
c(62, 3)=    0.43104153603246D-07
c(63, 3)=    0.41737189429633D-08
c(64, 3)=   -0.76395988054220D-08
c(65, 3)=    0.18729578607082D-07
c(66, 3)=    0.17754380139907D-07
c(67, 3)=   -0.84496500704376D-08
c(68, 3)=    0.16529568807175D-07
c(69, 3)=    0.16310377740880D-07
c(70, 3)=    0.35280078353266D-08
c(71, 3)=    0.18307764248237D-07
c(72, 3)=    0.31724949328092D-08
c(73, 3)=    0.11816052480633D-07
c(74, 3)=    0.17494675215840D-07
c(75, 3)=    0.26118248858656D-08
c(76, 3)=    0.30023266822151D-07
c(77, 3)=   -0.65712049623256D-08
c(78, 3)=    0.24493993741607D-07
c(79, 3)=   -0.81982981193597D-09
c(80, 3)=    0.14472986941000D-07
c( 4, 4)=    0.30959206658010D-06
c( 5, 4)=   -0.46336405889077D-05
c( 6, 4)=    0.10065686648503D-05
c( 7, 4)=    0.24640560610782D-05
c( 8, 4)=    0.15844031607287D-05
c( 9, 4)=    0.29405976756683D-06
c(10, 4)=   -0.12101666736715D-05
c(11, 4)=   -0.15688434008850D-05
c(12, 4)=   -0.10108862643404D-06
c(13, 4)=    0.17054914228551D-06
c(14, 4)=    0.22728427855596D-06
c(15, 4)=   -0.71711294666294D-06
c(16, 4)=   -0.77468414495595D-06
c(17, 4)=    0.45639058452256D-06
c(18, 4)=    0.87264101737834D-06
c(19, 4)=    0.13800527907950D-06
c(20, 4)=   -0.35839832070397D-06
c(21, 4)=   -0.39829410795150D-06
c(22, 4)=    0.19549557765958D-06
c(23, 4)=    0.44405614702060D-06
c(24, 4)=    0.57394522076625D-07
c(25, 4)=   -0.61925818643782D-06
c(26, 4)=   -0.37209534817516D-06
c(27, 4)=    0.34334229006640D-06
c(28, 4)=    0.40345026906393D-06
c(29, 4)=   -0.18728923925709D-06
c(30, 4)=   -0.26245016435872D-06
c(31, 4)=   -0.31150685284646D-07
c(32, 4)=    0.20037736081893D-06
c(33, 4)=    0.91242269592696D-07
c(34, 4)=   -0.29153815376169D-07
c(35, 4)=   -0.17807633359855D-06
c(36, 4)=   -0.15303133769982D-06
c(37, 4)=    0.10039017742705D-06
c(38, 4)=    0.17157873926809D-06
c(39, 4)=   -0.14539420941324D-07
c(40, 4)=   -0.89280013487568D-07
c(41, 4)=    0.37666735858964D-07
c(42, 4)=    0.13264285074014D-06
c(43, 4)=    0.51205577350747D-07
c(44, 4)=    0.42197628209560D-07
c(45, 4)=    0.37851610019160D-07
c(46, 4)=   -0.75655211066113D-08
c(47, 4)=    0.90026651213639D-07
c(48, 4)=    0.93009448368906D-07
c(49, 4)=   -0.23184611520205D-07
c(50, 4)=   -0.20719026291426D-07
c(51, 4)=    0.34422439963500D-07
c(52, 4)=    0.21952192711957D-07
c(53, 4)=   -0.21041597028035D-07
c(54, 4)=   -0.24228687045471D-07
c(55, 4)=    0.37895570193813D-08
c(56, 4)=   -0.19352611895629D-07
c(57, 4)=    0.16368584240792D-07
c(58, 4)=    0.69419737541584D-07
c(59, 4)=   -0.69539185815161D-08
c(60, 4)=   -0.18799793363085D-07
c(61, 4)=    0.14374232953337D-07
c(62, 4)=   -0.24130754111132D-07
c(63, 4)=   -0.36015801612887D-08
c(64, 4)=    0.40250447914756D-08
c(65, 4)=    0.10342517649973D-07
c(66, 4)=   -0.37158460287582D-07
c(67, 4)=   -0.35013020730218D-08
c(68, 4)=   -0.26839843046241D-09
c(69, 4)=    0.76971770252262D-09
c(70, 4)=    0.17192596869177D-07
c(71, 4)=    0.26611715724886D-08
c(72, 4)=    0.51616950532301D-08
c(73, 4)=    0.17438582702884D-07
c(74, 4)=    0.11010503509019D-07
c(75, 4)=    0.61518542523363D-08
c(76, 4)=    0.39330429681511D-08
c(77, 4)=    0.70796199616500D-08
c(78, 4)=   -0.78748564738901D-08
c(79, 4)=    0.77449061439055D-08
c(80, 4)=   -0.25351771157557D-08
c( 5, 5)=   -0.44434141649300D-05
c( 6, 5)=    0.16548025737626D-05
c( 7, 5)=   -0.19096934491812D-06
c( 8, 5)=   -0.27850787272597D-05
c( 9, 5)=   -0.22656267960896D-05
c(10, 5)=    0.42032262399202D-06
c(11, 5)=    0.13520771564056D-05
c(12, 5)=    0.71777774077014D-06
c(13, 5)=    0.61450262542325D-08
c(14, 5)=   -0.45276921372495D-06
c(15, 5)=   -0.93462951523554D-06
c(16, 5)=   -0.34377637201001D-06
c(17, 5)=    0.69789569926980D-06
c(18, 5)=    0.20629370256866D-06
c(19, 5)=   -0.28358145151180D-06
c(20, 5)=    0.16848241831807D-06
c(21, 5)=    0.32529795139056D-06
c(22, 5)=    0.72679019916679D-07
c(23, 5)=   -0.31654147385163D-06
c(24, 5)=   -0.34631841874185D-06
c(25, 5)=    0.28611650642088D-06
c(26, 5)=    0.20360981554297D-06
c(27, 5)=    0.96293497515718D-07
c(28, 5)=   -0.47930164566393D-07
c(29, 5)=   -0.12691726734393D-06
c(30, 5)=   -0.93858451381326D-07
c(31, 5)=    0.92431999857506D-07
c(32, 5)=    0.93351803773724D-07
c(33, 5)=   -0.15400246903553D-06
c(34, 5)=   -0.14119500703571D-06
c(35, 5)=    0.78736519330474D-07
c(36, 5)=    0.83329148294896D-07
c(37, 5)=    0.73174483987068D-07
c(38, 5)=   -0.55893430456265D-07
c(39, 5)=   -0.11157629707788D-06
c(40, 5)=    0.54238258568734D-07
c(41, 5)=    0.11263379495700D-06
c(42, 5)=    0.41202421137030D-07
c(43, 5)=   -0.59032315787153D-07
c(44, 5)=    0.26883253020711D-07
c(45, 5)=    0.51890763975186D-07
c(46, 5)=    0.13952729851129D-07
c(47, 5)=    0.83317364677177D-08
c(48, 5)=   -0.13102674596668D-07
c(49, 5)=   -0.21274731423140D-07
c(50, 5)=   -0.54059556337886D-07
c(51, 5)=   -0.18854105667823D-07
c(52, 5)=    0.22621073263953D-07
c(53, 5)=   -0.19805104368574D-07
c(54, 5)=   -0.19550955358517D-07
c(55, 5)=   -0.70267276107722D-08
c(56, 5)=   -0.30278638914856D-07
c(57, 5)=   -0.80729568890359D-08
c(58, 5)=   -0.33069309873432D-08
c(59, 5)=    0.20020433625436D-07
c(60, 5)=   -0.25864119973521D-07
c(61, 5)=   -0.39757987838197D-07
c(62, 5)=    0.64985926165149D-08
c(63, 5)=    0.19392998715642D-07
c(64, 5)=   -0.12325977875120D-07
c(65, 5)=   -0.23164237582529D-07
c(66, 5)=    0.12950766375034D-07
c(67, 5)=    0.28954128332616D-09
c(68, 5)=    0.14496023097164D-07
c(69, 5)=    0.35344064539468D-07
c(70, 5)=    0.11948634893181D-07
c(71, 5)=    0.50615198739125D-08
c(72, 5)=    0.29805375871243D-08
c(73, 5)=   -0.67955847102354D-08
c(74, 5)=   -0.68758792996420D-08
c(75, 5)=    0.10638093088951D-07
c(76, 5)=   -0.14851134350421D-07
c(77, 5)=    0.19291099080097D-07
c(78, 5)=   -0.16916285273421D-07
c(79, 5)=    0.17689701598385D-07
c(80, 5)=   -0.38699679010122D-08
c( 6, 6)=    0.27571599805567D-05
c( 7, 6)=   -0.55805315254086D-06
c( 8, 6)=   -0.91195769387441D-06
c( 9, 6)=    0.80627407047627D-06
c(10, 6)=    0.66759500742877D-06
c(11, 6)=   -0.24393303094978D-06
c(12, 6)=   -0.40190828647835D-06
c(13, 6)=    0.43978424085633D-07
c(14, 6)=   -0.66013890282255D-07
c(15, 6)=   -0.62099703560131D-07
c(16, 6)=    0.14317921213443D-06
c(17, 6)=    0.63476326900737D-06
c(18, 6)=    0.22350772799262D-06
c(19, 6)=   -0.27396734006137D-06
c(20, 6)=    0.20104741582389D-06
c(21, 6)=    0.62410418046131D-07
c(22, 6)=   -0.26852512844608D-06
c(23, 6)=    0.38219653909366D-06
c(24, 6)=    0.31595908792017D-06
c(25, 6)=   -0.15880959481453D-06
c(26, 6)=   -0.32691684424342D-06
c(27, 6)=    0.27157517394775D-06
c(28, 6)=    0.45910522899934D-06
c(29, 6)=   -0.90360646208109D-07
c(30, 6)=   -0.14350225793397D-06
c(31, 6)=    0.58517088199898D-07
c(32, 6)=    0.82448774379966D-08
c(33, 6)=   -0.14039660982872D-06
c(34, 6)=   -0.10319846464728D-07
c(35, 6)=    0.16365404957830D-06
c(36, 6)=   -0.11439331698275D-06
c(37, 6)=   -0.17957148759724D-06
c(38, 6)=    0.99909642061381D-07
c(39, 6)=    0.73601014779209D-07
c(40, 6)=   -0.11129691323560D-06
c(41, 6)=   -0.50363372367066D-07
c(42, 6)=    0.48970749788252D-07
c(43, 6)=   -0.21003915377398D-07
c(44, 6)=   -0.99990391085428D-07
c(45, 6)=    0.33708262897175D-07
c(46, 6)=    0.12638028137008D-07
c(47, 6)=   -0.91631664850957D-08
c(48, 6)=   -0.88122254176770D-08
c(49, 6)=    0.40675867295780D-07
c(50, 6)=    0.52310566613043D-07
c(51, 6)=   -0.13272662627539D-07
c(52, 6)=   -0.73904952552898D-08
c(53, 6)=   -0.14993717253759D-07
c(54, 6)=   -0.44698503381972D-07
c(55, 6)=    0.17039893002879D-07
c(56, 6)=    0.15830410793048D-07
c(57, 6)=   -0.14936960045191D-07
c(58, 6)=    0.60206055639023D-08
c(59, 6)=    0.46368464974162D-07
c(60, 6)=    0.40230875207353D-07
c(61, 6)=   -0.15769084210498D-07
c(62, 6)=   -0.21166941664643D-07
c(63, 6)=    0.18703755701939D-07
c(64, 6)=   -0.27691990214752D-08
c(65, 6)=   -0.33352911189790D-07
c(66, 6)=   -0.35608979285730D-07
c(67, 6)=   -0.70033417115897D-08
c(68, 6)=   -0.19185441063505D-07
c(69, 6)=   -0.32700471757675D-07
c(70, 6)=    0.12564201099496D-07
c(71, 6)=   -0.15126602104398D-07
c(72, 6)=    0.13317694048512D-07
c(73, 6)=    0.13006332176297D-07
c(74, 6)=   -0.71969446649702D-08
c(75, 6)=    0.28709756345817D-07
c(76, 6)=   -0.19623113496586D-07
c(77, 6)=    0.18195733917348D-07
c(78, 6)=   -0.21640275443293D-07
c(79, 6)=    0.11458350447296D-07
c(80, 6)=   -0.12170377130640D-07
c( 7, 7)=    0.43995886876871D-06
c( 8, 7)=   -0.47261950038642D-06
c( 9, 7)=   -0.61651917619159D-06
c(10, 7)=    0.31718853806594D-06
c(11, 7)=    0.65859011010243D-06
c(12, 7)=    0.40292072266072D-06
c(13, 7)=   -0.76808712383584D-06
c(14, 7)=   -0.44557079395778D-06
c(15, 7)=    0.92204367511949D-06
c(16, 7)=    0.38659794764365D-06
c(17, 7)=    0.33821902436418D-06
c(18, 7)=   -0.19710160923205D-06
c(19, 7)=   -0.48845474893250D-06
c(20, 7)=    0.68730006673845D-07
c(21, 7)=    0.24155063400434D-06
c(22, 7)=    0.26175060260014D-06
c(23, 7)=   -0.29231102856546D-06
c(24, 7)=   -0.40727652556414D-06
c(25, 7)=    0.11980582758862D-06
c(26, 7)=    0.12312238978368D-06
c(27, 7)=    0.24812236262542D-08
c(28, 7)=   -0.17273809020356D-06
c(29, 7)=   -0.25850308737583D-06
c(30, 7)=   -0.84399266886281D-07
c(31, 7)=    0.24130381448124D-06
c(32, 7)=    0.16317990329089D-06
c(33, 7)=   -0.89111302464619D-07
c(34, 7)=   -0.17537996827503D-06
c(35, 7)=    0.36799721551061D-07
c(36, 7)=    0.10372056536401D-06
c(37, 7)=    0.19570404665206D-07
c(38, 7)=   -0.55267439226573D-07
c(39, 7)=   -0.79573733294900D-07
c(40, 7)=    0.23027537491822D-07
c(41, 7)=    0.69887735189113D-07
c(42, 7)=    0.72321482138545D-07
c(43, 7)=   -0.10565869255953D-07
c(44, 7)=   -0.10151552325077D-06
c(45, 7)=    0.35931009870257D-07
c(46, 7)=    0.75225923562369D-08
c(47, 7)=   -0.22752655202284D-07
c(48, 7)=    0.89441389869944D-08
c(49, 7)=    0.37220465089298D-07
c(50, 7)=   -0.85192448435310D-08
c(51, 7)=   -0.40081109799497D-07
c(52, 7)=    0.23073769709323D-07
c(53, 7)=    0.45667290438653D-07
c(54, 7)=   -0.17445936216360D-07
c(55, 7)=   -0.15369846125184D-07
c(56, 7)=    0.18826249672075D-07
c(57, 7)=   -0.31417790708496D-07
c(58, 7)=   -0.24984816683585D-07
c(59, 7)=    0.37629174651362D-07
c(60, 7)=    0.56899536588022D-08
c(61, 7)=   -0.10496967073751D-07
c(62, 7)=    0.28100413921520D-07
c(63, 7)=    0.45442449578000D-08
c(64, 7)=   -0.11540662509264D-07
c(65, 7)=   -0.25710679867463D-07
c(66, 7)=   -0.19698525852962D-07
c(67, 7)=   -0.13303329727190D-07
c(68, 7)=    0.10982345162823D-07
c(69, 7)=   -0.21402856717440D-08
c(70, 7)=   -0.64088891211600D-08
c(71, 7)=   -0.13008094605731D-07
c(72, 7)=   -0.86969840430305D-08
c(73, 7)=   -0.22626636302946D-07
c(74, 7)=    0.54862757485682D-08
c(75, 7)=   -0.83509452181555D-08
c(76, 7)=   -0.12666008338105D-07
c(77, 7)=   -0.76562284824962D-08
c(78, 7)=   -0.17724487244370D-07
c(79, 7)=   -0.36520982438646D-08
c(80, 7)=   -0.16859220330869D-07
c( 8, 8)=   -0.30973350203001D-06
c( 9, 8)=    0.12047847138462D-05
c(10, 8)=    0.53865780719934D-06
c(11, 8)=   -0.11785529186851D-05
c(12, 8)=   -0.15982435716482D-05
c(13, 8)=   -0.23194739260953D-06
c(14, 8)=    0.81535531994532D-06
c(15, 8)=    0.11565702748603D-05
c(16, 8)=    0.23447441593593D-06
c(17, 8)=   -0.36982763604740D-06
c(18, 8)=    0.61118874148442D-07
c(19, 8)=    0.42698157940794D-06
c(20, 8)=    0.31067740467773D-06
c(21, 8)=    0.40705531212846D-07
c(22, 8)=   -0.37201577463817D-06
c(23, 8)=   -0.21605680473365D-06
c(24, 8)=    0.18543841568708D-06
c(25, 8)=    0.24648179571309D-06
c(26, 8)=    0.86142134273613D-07
c(27, 8)=   -0.34810645313143D-06
c(28, 8)=   -0.24819931440283D-06
c(29, 8)=   -0.41746703183486D-07
c(30, 8)=    0.53557629619836D-07
c(31, 8)=   -0.31189662495483D-08
c(32, 8)=   -0.23684136584399D-06
c(33, 8)=   -0.11135804639784D-06
c(34, 8)=    0.21017986319691D-07
c(35, 8)=    0.18144509931141D-06
c(36, 8)=    0.68093513349118D-07
c(37, 8)=   -0.24994894975000D-06
c(38, 8)=   -0.78471090833385D-07
c(39, 8)=    0.13239300047411D-06
c(40, 8)=    0.54880323402336D-07
c(41, 8)=    0.13430167979081D-07
c(42, 8)=    0.56339230993143D-07
c(43, 8)=    0.28565133925337D-08
c(44, 8)=   -0.42488219823149D-07
c(45, 8)=   -0.73817974736043D-08
c(46, 8)=    0.69193895041098D-07
c(47, 8)=    0.20359281649319D-07
c(48, 8)=   -0.47559903158974D-07
c(49, 8)=    0.11527250499104D-07
c(50, 8)=    0.31827816790805D-07
c(51, 8)=   -0.62623816165509D-07
c(52, 8)=   -0.35933127317171D-07
c(53, 8)=    0.19711901123882D-07
c(54, 8)=   -0.77427771586048D-08
c(55, 8)=   -0.41303575042150D-07
c(56, 8)=   -0.15072982855835D-07
c(57, 8)=   -0.20451938702364D-08
c(58, 8)=   -0.37398071041519D-07
c(59, 8)=    0.23820342762430D-07
c(60, 8)=    0.25502369443082D-07
c(61, 8)=   -0.15910304881318D-07
c(62, 8)=    0.13795120944057D-07
c(63, 8)=    0.32525287884179D-07
c(64, 8)=    0.13809423551226D-08
c(65, 8)=   -0.10812715236635D-09
c(66, 8)=    0.67634725831510D-09
c(67, 8)=   -0.16676212050011D-07
c(68, 8)=    0.12322542466311D-07
c(69, 8)=   -0.87952271938768D-08
c(70, 8)=    0.70662206807584D-08
c(71, 8)=    0.97306837784429D-08
c(72, 8)=   -0.10409665918863D-07
c(73, 8)=    0.22526300287359D-07
c(74, 8)=   -0.11929360719090D-07
c(75, 8)=    0.56724970381373D-08
c(76, 8)=    0.39296322296042D-08
c(77, 8)=    0.40729954548785D-09
c(78, 8)=   -0.14579810396259D-07
c(79, 8)=    0.33817539135385D-08
c(80, 8)=   -0.12114249017571D-07
c( 9, 9)=   -0.11685977561307D-05
c(10, 9)=   -0.14529152650329D-05
c(11, 9)=   -0.41376264440044D-06
c(12, 9)=    0.70595661731973D-06
c(13, 9)=    0.11036467020707D-05
c(14, 9)=    0.20920331569049D-06
c(15, 9)=   -0.30020191696171D-06
c(16, 9)=   -0.35633738393987D-06
c(17, 9)=   -0.19460945904606D-06
c(18, 9)=    0.69551342431549D-07
c(19, 9)=    0.21017949476738D-06
c(20, 9)=    0.58506560488758D-07
c(21, 9)=   -0.17704293895129D-06
c(22, 9)=   -0.79873680455104D-09
c(23, 9)=    0.23863912649051D-07
c(24, 9)=    0.11269490152068D-06
c(25, 9)=    0.80663857082553D-07
c(26, 9)=   -0.52332864057279D-07
c(27, 9)=   -0.40077218712098D-07
c(28, 9)=    0.92385377802570D-07
c(29, 9)=    0.89560884551153D-07
c(30, 9)=   -0.74814774016663D-07
c(31, 9)=   -0.13348919881375D-06
c(32, 9)=   -0.46162437050820D-08
c(33, 9)=    0.86098107102499D-07
c(34, 9)=    0.35831903511112D-07
c(35, 9)=   -0.38733132408446D-07
c(36, 9)=   -0.43099792134629D-07
c(37, 9)=   -0.10107397896816D-06
c(38, 9)=   -0.34748056946680D-07
c(39, 9)=    0.10009476196909D-06
c(40, 9)=    0.60783352459747D-07
c(41, 9)=   -0.70182088303399D-07
c(42, 9)=   -0.18055961756054D-07
c(43, 9)=    0.57484699400282D-07
c(44, 9)=    0.10553683342488D-07
c(45, 9)=   -0.42685258985362D-07
c(46, 9)=    0.49309943027769D-07
c(47, 9)=    0.22793754906480D-07
c(48, 9)=   -0.23998217097360D-07
c(49, 9)=    0.10914696187724D-07
c(50, 9)=    0.32781688335377D-07
c(51, 9)=   -0.62506643838801D-08
c(52, 9)=   -0.64979961127411D-07
c(53, 9)=    0.51096919194121D-07
c(54, 9)=    0.29356063667930D-07
c(55, 9)=   -0.17802875815950D-07
c(56, 9)=   -0.33534152808031D-07
c(57, 9)=   -0.11201661074046D-07
c(58, 9)=   -0.17348653702844D-07
c(59, 9)=    0.30891923032634D-08
c(60, 9)=    0.31455792355522D-07
c(61, 9)=   -0.15441400708974D-07
c(62, 9)=   -0.48101211929862D-08
c(63, 9)=   -0.67324611707600D-08
c(64, 9)=   -0.23034248154764D-07
c(65, 9)=   -0.87353984493366D-08
c(66, 9)=    0.29023336917157D-07
c(67, 9)=   -0.12812531346711D-07
c(68, 9)=   -0.74928735796058D-08
c(69, 9)=    0.18566108048567D-07
c(70, 9)=   -0.88858650319843D-08
c(71, 9)=   -0.44843642059847D-08
c(72, 9)=    0.14316683400634D-07
c(73, 9)=   -0.50703872752401D-09
c(74, 9)=    0.80343938717078D-08
c(75, 9)=    0.24856787548675D-07
c(76, 9)=    0.10960775781968D-08
c(77, 9)=    0.25742924795831D-07
c(78, 9)=    0.84686476194344D-08
c(79, 9)=    0.76200924005946D-08
c(80, 9)=    0.12750704377667D-07
c(10,10)=   -0.27418592728775D-06
c(11,10)=    0.33152053925569D-06
c(12,10)=    0.48642163810406D-06
c(13,10)=   -0.15737288606159D-06
c(14,10)=    0.25566932292312D-07
c(15,10)=   -0.15979532465480D-06
c(16,10)=   -0.51356597033173D-06
c(17,10)=   -0.56809357230694D-06
c(18,10)=    0.66632156103219D-06
c(19,10)=    0.60140162556967D-06
c(20,10)=   -0.42013699476283D-06
c(21,10)=   -0.28400079019097D-06
c(22,10)=   -0.25361012940757D-06
c(23,10)=   -0.89752221004529D-07
c(24,10)=   -0.12256105013166D-06
c(25,10)=    0.21528872721283D-06
c(26,10)=    0.14690117401347D-06
c(27,10)=   -0.31119309740966D-06
c(28,10)=   -0.82984621401518D-07
c(29,10)=    0.13876448127514D-06
c(30,10)=    0.67547141689097D-07
c(31,10)=   -0.80968956240002D-07
c(32,10)=   -0.25061516036457D-07
c(33,10)=    0.48032933510956D-08
c(34,10)=   -0.20507365083276D-07
c(35,10)=   -0.55221010597376D-07
c(36,10)=    0.75914839445849D-07
c(37,10)=   -0.34872065650435D-08
c(38,10)=   -0.52058651133105D-07
c(39,10)=    0.84942358720890D-07
c(40,10)=    0.11866701531355D-07
c(41,10)=   -0.76495579940316D-07
c(42,10)=   -0.17727920207917D-07
c(43,10)=    0.61463584510654D-07
c(44,10)=    0.33270487334672D-07
c(45,10)=   -0.87602100996217D-07
c(46,10)=    0.26372866601973D-07
c(47,10)=    0.30349490906782D-07
c(48,10)=   -0.27923311429486D-07
c(49,10)=   -0.13302497453275D-07
c(50,10)=   -0.25820559835020D-07
c(51,10)=   -0.23374649723249D-07
c(52,10)=    0.82265363238744D-08
c(53,10)=    0.46393708341138D-07
c(54,10)=   -0.15072675889249D-07
c(55,10)=   -0.34347539063150D-07
c(56,10)=    0.20485320872554D-07
c(57,10)=    0.42233478341289D-07
c(58,10)=    0.50401390086354D-07
c(59,10)=   -0.10171279820395D-07
c(60,10)=   -0.35575784347796D-07
c(61,10)=    0.42127433292415D-07
c(62,10)=   -0.10118118317467D-07
c(63,10)=    0.91817545738375D-08
c(64,10)=    0.35421298579918D-07
c(65,10)=    0.12182613264688D-07
c(66,10)=   -0.64310568888630D-08
c(67,10)=    0.46161869451084D-08
c(68,10)=   -0.29013715978110D-07
c(69,10)=    0.19689845799256D-07
c(70,10)=   -0.10549308957055D-08
c(71,10)=    0.17281557475445D-07
c(72,10)=   -0.11057461689040D-07
c(73,10)=   -0.97301930678866D-08
c(74,10)=   -0.81094277056530D-08
c(75,10)=   -0.37042327999666D-09
c(76,10)=   -0.15641894574580D-07
c(77,10)=    0.11854403909782D-07
c(78,10)=   -0.54294732577521D-08
c(79,10)=    0.78724534170822D-08
c(80,10)=    0.42317866507287D-08
c(11,11)=   -0.44684635932806D-07
c(12,11)=    0.86209182518877D-06
c(13,11)=    0.75993437561620D-06
c(14,11)=   -0.92572155775018D-06
c(15,11)=   -0.77135611592181D-06
c(16,11)=   -0.66555587400325D-07
c(17,11)=    0.35050392061217D-06
c(18,11)=    0.27473842288581D-06
c(19,11)=    0.44343867592765D-07
c(20,11)=   -0.16158145338117D-06
c(21,11)=   -0.27987346427304D-06
c(22,11)=   -0.86898719785704D-07
c(23,11)=    0.33488893134310D-06
c(24,11)=   -0.35108546304128D-08
c(25,11)=   -0.42067068385092D-06
c(26,11)=   -0.46190939130907D-06
c(27,11)=   -0.19330288476283D-06
c(28,11)=    0.21188939514542D-06
c(29,11)=    0.29489368172660D-06
c(30,11)=    0.13455548251108D-06
c(31,11)=   -0.19411432800647D-06
c(32,11)=   -0.59424865170610D-07
c(33,11)=    0.17417689865144D-06
c(34,11)=    0.12434837415707D-06
c(35,11)=   -0.56095234006330D-07
c(36,11)=   -0.78614167217231D-07
c(37,11)=    0.76271352504045D-07
c(38,11)=    0.44466538643800D-07
c(39,11)=    0.49697346366543D-07
c(40,11)=    0.95342048484995D-07
c(41,11)=   -0.22060919870347D-07
c(42,11)=   -0.10627362439600D-06
c(43,11)=    0.64554255720816D-08
c(44,11)=    0.58642333698728D-07
c(45,11)=   -0.26890924444448D-07
c(46,11)=   -0.35211105586852D-07
c(47,11)=   -0.10464139564664D-08
c(48,11)=   -0.18663789655521D-07
c(49,11)=    0.36421059241949D-08
c(50,11)=    0.39530294301080D-08
c(51,11)=   -0.46538292713927D-08
c(52,11)=   -0.73516296354175D-07
c(53,11)=    0.71800656224743D-08
c(54,11)=    0.19838048143693D-07
c(55,11)=   -0.27800604308575D-07
c(56,11)=   -0.23549993100805D-07
c(57,11)=    0.23355590833905D-07
c(58,11)=    0.12048348807616D-07
c(59,11)=   -0.35020125986908D-07
c(60,11)=    0.13173194177869D-07
c(61,11)=    0.50660334043555D-08
c(62,11)=   -0.93139166850523D-08
c(63,11)=    0.10367939361458D-07
c(64,11)=   -0.20853446969930D-07
c(65,11)=    0.15681558651406D-07
c(66,11)=   -0.14422762490918D-07
c(67,11)=    0.19077120844381D-08
c(68,11)=   -0.14617125298418D-07
c(69,11)=   -0.27984412638652D-07
c(70,11)=    0.12427476941900D-08
c(71,11)=    0.14146068560767D-09
c(72,11)=    0.10794013713742D-07
c(73,11)=   -0.33633551806085D-08
c(74,11)=    0.73440033006451D-08
c(75,11)=   -0.19989334493278D-07
c(76,11)=    0.18759121278194D-07
c(77,11)=   -0.19463700567277D-07
c(78,11)=    0.10058388000952D-07
c(79,11)=   -0.82265261612786D-08
c(80,11)=    0.23951018907818D-08
c(12,12)=   -0.64231718972073D-08
c(13,12)=   -0.14518332762979D-05
c(14,12)=   -0.44667389507501D-06
c(15,12)=    0.97595914392990D-06
c(16,12)=    0.40830279631007D-06
c(17,12)=    0.20089686057968D-06
c(18,12)=    0.73988670101664D-07
c(19,12)=   -0.11269182077243D-06
c(20,12)=   -0.47171328165674D-06
c(21,12)=    0.25397802768283D-06
c(22,12)=    0.37922079263538D-06
c(23,12)=    0.74687070585526D-07
c(24,12)=   -0.29603322643162D-06
c(25,12)=   -0.32919144778195D-06
c(26,12)=    0.22776493253620D-06
c(27,12)=    0.59673534864034D-07
c(28,12)=    0.89182305751498D-07
c(29,12)=   -0.14013340779859D-06
c(30,12)=   -0.30325844797736D-06
c(31,12)=   -0.83515417019970D-07
c(32,12)=    0.69433465327683D-07
c(33,12)=    0.17567234786505D-06
c(34,12)=   -0.70733077990148D-07
c(35,12)=   -0.15850573423284D-06
c(36,12)=   -0.37550527011772D-07
c(37,12)=    0.78044879320411D-07
c(38,12)=   -0.10785624108895D-07
c(39,12)=   -0.77044658435730D-07
c(40,12)=   -0.43096549983760D-10
c(41,12)=   -0.79751941321668D-07
c(42,12)=   -0.27619361101061D-07
c(43,12)=    0.86397321275642D-07
c(44,12)=    0.39211827033959D-07
c(45,12)=   -0.10868277448053D-06
c(46,12)=   -0.13572773598541D-06
c(47,12)=    0.67101554496870D-07
c(48,12)=    0.10473935783814D-06
c(49,12)=   -0.44412396872399D-07
c(50,12)=   -0.58344190655476D-07
c(51,12)=    0.34682409657115D-07
c(52,12)=    0.44490481568391D-07
c(53,12)=   -0.60142679556923D-07
c(54,12)=   -0.41393872946390D-07
c(55,12)=    0.24090662325969D-07
c(56,12)=    0.19446206959575D-07
c(57,12)=    0.19876152940267D-07
c(58,12)=   -0.80326007048200D-08
c(59,12)=   -0.11205579119766D-07
c(60,12)=    0.16730592865468D-07
c(61,12)=    0.30935461684387D-07
c(62,12)=    0.63355437631645D-08
c(63,12)=   -0.80774330138626D-08
c(64,12)=   -0.54235063151270D-08
c(65,12)=   -0.90466913004512D-08
c(66,12)=   -0.18000612328190D-07
c(67,12)=    0.36176316060519D-07
c(68,12)=   -0.14836236343949D-07
c(69,12)=   -0.10115866668775D-07
c(70,12)=    0.61568005209749D-09
c(71,12)=   -0.63328475344205D-08
c(72,12)=   -0.13228133600474D-07
c(73,12)=    0.11675415257145D-07
c(74,12)=   -0.88084298812367D-08
c(75,12)=   -0.16421185729936D-08
c(76,12)=   -0.28705724964639D-08
c(77,12)=    0.76454526940509D-08
c(78,12)=   -0.97983270912013D-08
c(79,12)=    0.93094685912931D-08
c(80,12)=   -0.14952267692181D-07
c(13,13)=    0.42601241621664D-06
c(14,13)=    0.75129408112081D-06
c(15,13)=    0.17845595134875D-06
c(16,13)=   -0.64644274826991D-07
c(17,13)=   -0.44363357576703D-06
c(18,13)=   -0.29893151037383D-06
c(19,13)=   -0.20081378465993D-06
c(20,13)=   -0.38974995460658D-09
c(21,13)=    0.26038531389766D-08
c(22,13)=   -0.16993532987053D-06
c(23,13)=   -0.22669233943721D-06
c(24,13)=    0.77346455693605D-07
c(25,13)=   -0.18200920903674D-07
c(26,13)=    0.51923621297858D-07
c(27,13)=   -0.69534755150605D-07
c(28,13)=   -0.13011181784465D-06
c(29,13)=    0.17478512424802D-06
c(30,13)=    0.97172404864323D-07
c(31,13)=    0.29060801404821D-07
c(32,13)=    0.49415299254759D-07
c(33,13)=   -0.63664146670117D-09
c(34,13)=   -0.42098806139770D-07
c(35,13)=   -0.63199371645429D-07
c(36,13)=    0.14925139462609D-06
c(37,13)=    0.77616737374773D-07
c(38,13)=   -0.49196687691755D-07
c(39,13)=   -0.21612344963652D-07
c(40,13)=    0.80541885220445D-08
c(41,13)=    0.43986031382998D-07
c(42,13)=    0.53831179685585D-09
c(43,13)=   -0.48153297585660D-07
c(44,13)=   -0.37336817198338D-07
c(45,13)=   -0.72522229447116D-07
c(46,13)=    0.52595884906078D-07
c(47,13)=    0.11224855683888D-07
c(48,13)=   -0.84039896860964D-08
c(49,13)=   -0.10002921042261D-07
c(50,13)=    0.20111357105294D-07
c(51,13)=   -0.37326778289764D-07
c(52,13)=   -0.49297223414475D-07
c(53,13)=    0.24567721109946D-07
c(54,13)=    0.36389985182291D-07
c(55,13)=   -0.33801457283125D-07
c(56,13)=   -0.45368225288631D-07
c(57,13)=    0.97773492056029D-08
c(58,13)=    0.32131377937284D-07
c(59,13)=    0.26674049173607D-07
c(60,13)=    0.88098981589739D-08
c(61,13)=    0.22112055523099D-07
c(62,13)=    0.26945186382454D-08
c(63,13)=   -0.37446774164464D-07
c(64,13)=    0.24823107737634D-07
c(65,13)=    0.16925015623816D-07
c(66,13)=   -0.16280148223969D-07
c(67,13)=    0.13036214805850D-08
c(68,13)=   -0.57340326176229D-08
c(69,13)=    0.29455152151528D-07
c(70,13)=   -0.14619484493505D-07
c(71,13)=    0.90805605878318D-08
c(72,13)=    0.70049802313740D-08
c(73,13)=   -0.17737979012555D-07
c(74,13)=    0.10487531174343D-07
c(75,13)=   -0.24750045534933D-07
c(76,13)=    0.15308268186489D-07
c(77,13)=   -0.13983560080127D-07
c(78,13)=    0.15885401503704D-07
c(79,13)=   -0.27785701635593D-08
c(80,13)=    0.11282900601153D-07
c(14,14)=    0.47021295278236D-07
c(15,14)=    0.22219457857795D-06
c(16,14)=   -0.16944993310769D-06
c(17,14)=   -0.35316177030083D-06
c(18,14)=   -0.36574949577333D-06
c(19,14)=    0.15370674221947D-06
c(20,14)=    0.28054717681758D-06
c(21,14)=    0.13450955702606D-06
c(22,14)=   -0.14565214032387D-06
c(23,14)=   -0.83807732025708D-07
c(24,14)=    0.33039361477979D-08
c(25,14)=   -0.76090567719229D-07
c(26,14)=    0.12905418648105D-06
c(27,14)=    0.17337430654230D-06
c(28,14)=    0.80857002198244D-07
c(29,14)=    0.70110690240570D-07
c(30,14)=   -0.13733277347911D-06
c(31,14)=   -0.21562170698916D-07
c(32,14)=    0.21759747165202D-07
c(33,14)=    0.61333748496406D-07
c(34,14)=    0.38829485595747D-07
c(35,14)=   -0.64086563137960D-07
c(36,14)=   -0.10153809047224D-06
c(37,14)=    0.32039980307240D-08
c(38,14)=    0.36304249021626D-07
c(39,14)=   -0.27257848252935D-07
c(40,14)=    0.53269190848855D-08
c(41,14)=   -0.51651425829049D-07
c(42,14)=   -0.57946053082766D-07
c(43,14)=    0.45498065486750D-08
c(44,14)=    0.47767471675605D-07
c(45,14)=   -0.55594018064780D-08
c(46,14)=   -0.10999591526210D-06
c(47,14)=    0.18923122289816D-07
c(48,14)=    0.61960720902588D-07
c(49,14)=   -0.60598476698241D-07
c(50,14)=   -0.67740190638225D-07
c(51,14)=    0.25780522864084D-07
c(52,14)=    0.63555181447774D-07
c(53,14)=   -0.12055036245110D-07
c(54,14)=   -0.53230826865350D-07
c(55,14)=    0.85375653220950D-08
c(56,14)=   -0.18483807864022D-07
c(57,14)=   -0.25455211194517D-09
c(58,14)=    0.12082251417562D-07
c(59,14)=    0.93215676397531D-08
c(60,14)=   -0.22104077649682D-07
c(61,14)=    0.23226932687473D-07
c(62,14)=    0.36234693962067D-07
c(63,14)=    0.11914217891800D-07
c(64,14)=   -0.17905878561972D-07
c(65,14)=    0.11619420494160D-07
c(66,14)=    0.57319890444538D-08
c(67,14)=    0.19784210676846D-07
c(68,14)=   -0.17912811632289D-07
c(69,14)=   -0.80003248139404D-08
c(70,14)=    0.23602690947874D-07
c(71,14)=   -0.94580549890193D-08
c(72,14)=    0.14153322210681D-07
c(73,14)=   -0.94356315047272D-08
c(74,14)=    0.28953059753020D-08
c(75,14)=   -0.25026524449757D-07
c(76,14)=    0.91427701028891D-08
c(77,14)=   -0.23200692086839D-07
c(78,14)=    0.72605453361048D-08
c(79,14)=   -0.12468131635802D-07
c(80,14)=    0.48904006092653D-08
c(15,15)=   -0.33826575599910D-06
c(16,15)=   -0.85495972949958D-06
c(17,15)=   -0.36408938792661D-06
c(18,15)=    0.35819799150966D-06
c(19,15)=    0.82529172013831D-06
c(20,15)=    0.78475364933386D-07
c(21,15)=   -0.32932122376644D-06
c(22,15)=   -0.41965832774604D-06
c(23,15)=   -0.27719475994770D-06
c(24,15)=    0.21532172296135D-06
c(25,15)=    0.44144073998089D-06
c(26,15)=    0.94034486508637D-07
c(27,15)=   -0.50527553014782D-07
c(28,15)=   -0.87669846139116D-07
c(29,15)=    0.41888714845582D-07
c(30,15)=    0.12166056752935D-06
c(31,15)=    0.88299403590952D-07
c(32,15)=    0.46689356080411D-07
c(33,15)=   -0.12797747227092D-06
c(34,15)=   -0.11305342583159D-06
c(35,15)=    0.79647332221832D-07
c(36,15)=    0.25153748542570D-06
c(37,15)=    0.43048298307130D-07
c(38,15)=   -0.69563590983603D-07
c(39,15)=   -0.46081191209894D-07
c(40,15)=   -0.98718472252429D-07
c(41,15)=    0.94967876174956D-07
c(42,15)=    0.97997180838959D-07
c(43,15)=   -0.53778195033675D-08
c(44,15)=   -0.42357934669160D-07
c(45,15)=   -0.12123269508550D-07
c(46,15)=    0.11349722928787D-06
c(47,15)=    0.11067875281430D-07
c(48,15)=   -0.13053510212507D-07
c(49,15)=   -0.14638954657636D-07
c(50,15)=    0.49823097696749D-07
c(51,15)=    0.56629472909396D-07
c(52,15)=   -0.27650738137218D-07
c(53,15)=   -0.10838656591511D-07
c(54,15)=    0.15888344512037D-07
c(55,15)=    0.86772124063542D-08
c(56,15)=   -0.80460750479536D-09
c(57,15)=    0.20227305510350D-08
c(58,15)=    0.31307506718066D-07
c(59,15)=   -0.21600268000035D-07
c(60,15)=    0.20834762794785D-07
c(61,15)=    0.19649645104048D-07
c(62,15)=   -0.10319384192510D-07
c(63,15)=   -0.28311597849162D-07
c(64,15)=   -0.72646727376410D-07
c(65,15)=    0.23316763995366D-07
c(66,15)=   -0.26509433663467D-08
c(67,15)=   -0.26642380677339D-07
c(68,15)=    0.30153600914930D-07
c(69,15)=   -0.34775082982131D-07
c(70,15)=   -0.15669476712055D-07
c(71,15)=    0.58140866847644D-08
c(72,15)=   -0.15813244179083D-07
c(73,15)=   -0.21739180642430D-07
c(74,15)=    0.51335559578317D-08
c(75,15)=   -0.14775847459000D-07
c(76,15)=    0.45298066015886D-08
c(77,15)=    0.12785646702054D-08
c(78,15)=   -0.16274874220888D-07
c(79,15)=    0.72392959718484D-08
c(80,15)=   -0.23845399951142D-07
c(16,16)=    0.20940180005886D-06
c(17,16)=    0.10696458697612D-05
c(18,16)=    0.41157939324289D-06
c(19,16)=   -0.61844626633270D-08
c(20,16)=    0.46528363403318D-07
c(21,16)=   -0.28509034304189D-06
c(22,16)=   -0.20032041830870D-06
c(23,16)=    0.43647221573199D-06
c(24,16)=    0.60074270247024D-06
c(25,16)=   -0.88444347559313D-07
c(26,16)=   -0.50086802979082D-06
c(27,16)=   -0.10370177302162D-06
c(28,16)=    0.17209430787186D-06
c(29,16)=    0.14527833834916D-06
c(30,16)=    0.67423619085705D-07
c(31,16)=   -0.13012023367532D-06
c(32,16)=   -0.18209698545153D-06
c(33,16)=    0.29330141545847D-07
c(34,16)=    0.23040995993151D-06
c(35,16)=    0.11137405742133D-06
c(36,16)=   -0.14322365665332D-06
c(37,16)=   -0.13393323975087D-06
c(38,16)=   -0.58158427107805D-07
c(39,16)=    0.66403824552258D-07
c(40,16)=    0.13326888331571D-06
c(41,16)=   -0.36830244580688D-07
c(42,16)=   -0.14702764073939D-06
c(43,16)=   -0.76430545168527D-07
c(44,16)=    0.11263069331412D-06
c(45,16)=    0.63796262698772D-07
c(46,16)=   -0.87686464926144D-07
c(47,16)=   -0.97964463362824D-08
c(48,16)=    0.79365091027521D-07
c(49,16)=    0.47651494645736D-07
c(50,16)=   -0.24066271282983D-07
c(51,16)=    0.56223103817103D-07
c(52,16)=    0.47985593815186D-07
c(53,16)=   -0.27226622049035D-07
c(54,16)=   -0.31213765163699D-07
c(55,16)=    0.27525580536844D-07
c(56,16)=    0.40704550615353D-07
c(57,16)=    0.14539391121680D-07
c(58,16)=   -0.13130877228086D-07
c(59,16)=    0.14699653481615D-07
c(60,16)=   -0.91637998362305D-08
c(61,16)=    0.81179669745568D-08
c(62,16)=    0.17151092609237D-07
c(63,16)=   -0.74429374406359D-08
c(64,16)=    0.67003182535719D-09
c(65,16)=   -0.35789700257476D-07
c(66,16)=   -0.22138650131817D-07
c(67,16)=   -0.30092279766274D-07
c(68,16)=   -0.18685560165448D-07
c(69,16)=   -0.26432788454162D-07
c(70,16)=    0.50424290960415D-08
c(71,16)=   -0.28808648811523D-07
c(72,16)=   -0.11382697449659D-07
c(73,16)=   -0.14815829721152D-08
c(74,16)=   -0.35217346892953D-07
c(75,16)=   -0.11363933421473D-07
c(76,16)=   -0.29903695425997D-07
c(77,16)=   -0.10568165816869D-07
c(78,16)=   -0.19737582568469D-07
c(79,16)=   -0.19551425947931D-07
c(80,16)=   -0.16322989519910D-08
c(17,17)=    0.54701387168645D-07
c(18,17)=    0.71135256682422D-07
c(19,17)=   -0.11665317418912D-06
c(20,17)=   -0.48568162906700D-06
c(21,17)=   -0.26328670190302D-06
c(22,17)=    0.33024327781270D-06
c(23,17)=    0.31284561172707D-06
c(24,17)=   -0.10503059994605D-06
c(25,17)=   -0.35143134866898D-06
c(26,17)=   -0.24326830441430D-06
c(27,17)=    0.22100890791401D-06
c(28,17)=    0.16019813922073D-06
c(29,17)=   -0.28819396574269D-07
c(30,17)=   -0.13171119227208D-06
c(31,17)=   -0.47817763167130D-07
c(32,17)=    0.12225078912847D-07
c(33,17)=    0.42623845131959D-07
c(34,17)=   -0.83790307497045D-07
c(35,17)=   -0.14399968540526D-06
c(36,17)=   -0.63350654236970D-08
c(37,17)=    0.61044104459295D-07
c(38,17)=    0.61949450540370D-07
c(39,17)=    0.44177767550151D-07
c(40,17)=   -0.71427176328117D-07
c(41,17)=    0.16646511539951D-07
c(42,17)=    0.30889394313482D-07
c(43,17)=   -0.59842571662947D-09
c(44,17)=    0.13614720720525D-07
c(45,17)=    0.11726912716675D-07
c(46,17)=    0.15012106749067D-07
c(47,17)=    0.25102121647612D-09
c(48,17)=   -0.13071448281777D-07
c(49,17)=   -0.34637835249062D-07
c(50,17)=   -0.32933308128326D-09
c(51,17)=    0.38939007315884D-07
c(52,17)=    0.13581147762298D-08
c(53,17)=   -0.93264032646212D-08
c(54,17)=    0.96767939090607D-08
c(55,17)=    0.36659358735424D-07
c(56,17)=    0.37491483020846D-08
c(57,17)=   -0.33698632494570D-08
c(58,17)=   -0.23932301146387D-08
c(59,17)=    0.42250810294220D-08
c(60,17)=    0.50024857363612D-08
c(61,17)=   -0.15884841625309D-07
c(62,17)=    0.16657425735481D-07
c(63,17)=    0.24035489426911D-07
c(64,17)=   -0.14758055967310D-07
c(65,17)=    0.22469107130498D-07
c(66,17)=   -0.46668229709955D-07
c(67,17)=   -0.30884264246269D-07
c(68,17)=   -0.28919828125255D-08
c(69,17)=    0.28658969794463D-08
c(70,17)=    0.29255472177867D-07
c(71,17)=    0.37611588847544D-08
c(72,17)=    0.84076306333866D-08
c(73,17)=    0.86650572592168D-08
c(74,17)=    0.15425043157405D-07
c(75,17)=   -0.13780132543492D-07
c(76,17)=    0.21776929781402D-07
c(77,17)=   -0.18096281174059D-07
c(78,17)=    0.80052085590603D-08
c(79,17)=   -0.16578597468070D-07
c(80,17)=   -0.13357084307772D-08
c(18,18)=    0.46656982881352D-06
c(19,18)=   -0.72203974438337D-06
c(20,18)=    0.35213761763877D-07
c(21,18)=    0.24334025394976D-06
c(22,18)=   -0.11671446266874D-06
c(23,18)=   -0.13241514416301D-06
c(24,18)=    0.36035680514514D-06
c(25,18)=    0.12798344905849D-06
c(26,18)=   -0.27199670154689D-06
c(27,18)=   -0.32523985662353D-07
c(28,18)=    0.27683282308885D-06
c(29,18)=    0.28581207181270D-07
c(30,18)=   -0.71344708930735D-07
c(31,18)=    0.30547665378528D-07
c(32,18)=    0.76591740848835D-07
c(33,18)=   -0.12586254793145D-06
c(34,18)=   -0.27422892767733D-07
c(35,18)=    0.11475323877322D-06
c(36,18)=   -0.33482636659461D-07
c(37,18)=   -0.11596384891875D-06
c(38,18)=   -0.25844182306635D-07
c(39,18)=    0.65617178872133D-07
c(40,18)=    0.58005502905905D-07
c(41,18)=   -0.11349233087356D-06
c(42,18)=   -0.72573977433579D-07
c(43,18)=    0.16836709714452D-07
c(44,18)=    0.58205565043271D-07
c(45,18)=    0.67135385454244D-07
c(46,18)=   -0.42278595866491D-07
c(47,18)=   -0.69071860661362D-07
c(48,18)=    0.27022559026842D-07
c(49,18)=    0.54608578248484D-07
c(50,18)=   -0.51026025277597D-07
c(51,18)=   -0.36922133821202D-07
c(52,18)=   -0.21243464334118D-08
c(53,18)=   -0.63537119606441D-08
c(54,18)=    0.35311423634546D-07
c(55,18)=   -0.64441967139307D-09
c(56,18)=   -0.63934099899106D-08
c(57,18)=   -0.17517136269818D-07
c(58,18)=   -0.12479573707963D-07
c(59,18)=    0.79659850821093D-08
c(60,18)=   -0.18721276925731D-07
c(61,18)=    0.17670204724543D-07
c(62,18)=    0.43540870673701D-08
c(63,18)=    0.16298063385003D-07
c(64,18)=    0.65394780694313D-08
c(65,18)=   -0.14960967688679D-07
c(66,18)=   -0.11950207010582D-07
c(67,18)=    0.16994741177372D-08
c(68,18)=    0.61450979011810D-08
c(69,18)=   -0.24717849546371D-07
c(70,18)=    0.21767138894177D-07
c(71,18)=   -0.29289376420573D-07
c(72,18)=   -0.14163742552877D-07
c(73,18)=    0.45636403434965D-08
c(74,18)=   -0.48803178398436D-07
c(75,18)=    0.48548247736960D-08
c(76,18)=   -0.19930049143878D-07
c(77,18)=   -0.55322134805482D-08
c(78,18)=   -0.21643701140752D-07
c(79,18)=    0.41623489747571D-08
c(80,18)=   -0.20210577532281D-07
c(19,19)=   -0.36488623031146D-06
c(20,19)=    0.17498191632079D-06
c(21,19)=    0.18234844161543D-06
c(22,19)=    0.41049548405926D-06
c(23,19)=    0.41382769125263D-06
c(24,19)=   -0.28283847306805D-06
c(25,19)=   -0.67171858587578D-06
c(26,19)=   -0.13281287624361D-07
c(27,19)=    0.73831317520578D-06
c(28,19)=    0.31697019030015D-06
c(29,19)=   -0.19570027691728D-06
c(30,19)=   -0.25865913112589D-06
c(31,19)=    0.66189798286727D-07
c(32,19)=    0.21940030137740D-06
c(33,19)=    0.21730441828881D-06
c(34,19)=    0.11070102602743D-07
c(35,19)=   -0.10875568103679D-06
c(36,19)=   -0.98549961802302D-07
c(37,19)=    0.11773175276404D-06
c(38,19)=    0.12373075793390D-06
c(39,19)=   -0.80672746880080D-07
c(40,19)=   -0.10499627369546D-06
c(41,19)=    0.97463003969700D-08
c(42,19)=    0.77667013884904D-07
c(43,19)=   -0.25448050543306D-08
c(44,19)=   -0.84328415252217D-08
c(45,19)=    0.19663700560579D-07
c(46,19)=   -0.57579508769271D-07
c(47,19)=    0.62220402666774D-07
c(48,19)=    0.80456306652161D-07
c(49,19)=   -0.10942297172432D-07
c(50,19)=   -0.59666734494807D-07
c(51,19)=   -0.29821353916073D-07
c(52,19)=    0.52731420460540D-07
c(53,19)=    0.16921809036510D-07
c(54,19)=   -0.66440525931939D-07
c(55,19)=   -0.77889674318441D-08
c(56,19)=    0.25265777071564D-07
c(57,19)=   -0.19961563366707D-08
c(58,19)=    0.34843747062013D-08
c(59,19)=    0.26913236517537D-08
c(60,19)=   -0.21092730784263D-07
c(61,19)=   -0.27368683972922D-07
c(62,19)=    0.10436698163346D-07
c(63,19)=    0.34728429204960D-08
c(64,19)=   -0.21415899983517D-07
c(65,19)=   -0.29172479955474D-07
c(66,19)=    0.20484328937738D-07
c(67,19)=   -0.77003006269710D-08
c(68,19)=    0.68077110643878D-08
c(69,19)=    0.27576319284877D-07
c(70,19)=   -0.17063671948984D-07
c(71,19)=    0.35652582279772D-08
c(72,19)=    0.51331374560066D-08
c(73,19)=   -0.16726262326640D-07
c(74,19)=    0.17879784259480D-07
c(75,19)=    0.68971937496097D-08
c(76,19)=   -0.61963011248938D-08
c(77,19)=    0.90614875079607D-08
c(78,19)=   -0.43906604006747D-08
c(79,19)=   -0.38438442046203D-08
c(80,19)=   -0.65132030997370D-08
c(20,20)=   -0.47388247284749D-06
c(21,20)=   -0.12339000794407D-06
c(22,20)=   -0.43652709220027D-06
c(23,20)=   -0.49050613929173D-06
c(24,20)=   -0.81235739566793D-07
c(25,20)=    0.56499626215347D-06
c(26,20)=    0.12278930428748D-06
c(27,20)=   -0.24818904834720D-06
c(28,20)=   -0.26539606385896D-06
c(29,20)=   -0.12786594632194D-06
c(30,20)=    0.12291704860035D-06
c(31,20)=    0.40161300329303D-06
c(32,20)=    0.15962202179924D-06
c(33,20)=   -0.16312204203802D-06
c(34,20)=   -0.17246098146981D-06
c(35,20)=    0.14159663364171D-06
c(36,20)=    0.18338783821004D-06
c(37,20)=    0.65842534620178D-07
c(38,20)=    0.17607805933065D-07
c(39,20)=    0.16630905489238D-08
c(40,20)=    0.47066879402688D-07
c(41,20)=    0.31993670191674D-07
c(42,20)=    0.24112045504797D-07
c(43,20)=    0.29386948733193D-08
c(44,20)=   -0.63550823151513D-07
c(45,20)=    0.50430699402188D-08
c(46,20)=    0.41103919023346D-07
c(47,20)=   -0.24645662834438D-07
c(48,20)=    0.53048318352021D-08
c(49,20)=   -0.16932289464472D-07
c(50,20)=   -0.41841399284243D-07
c(51,20)=   -0.11495118139598D-07
c(52,20)=    0.36232278648845D-07
c(53,20)=    0.15835361411647D-07
c(54,20)=   -0.17140178972305D-07
c(55,20)=   -0.60642737543764D-08
c(56,20)=    0.13380551924738D-08
c(57,20)=   -0.30162880909350D-08
c(58,20)=    0.19626397354767D-08
c(59,20)=    0.41233270172503D-07
c(60,20)=    0.21408192120022D-07
c(61,20)=   -0.78802683142367D-08
c(62,20)=   -0.17958417108302D-07
c(63,20)=   -0.56071100273366D-07
c(64,20)=    0.21343376042478D-08
c(65,20)=    0.54170022859117D-09
c(66,20)=   -0.31859741637542D-07
c(67,20)=    0.53342732542924D-07
c(68,20)=   -0.21166186689024D-07
c(69,20)=   -0.28514542148274D-08
c(70,20)=    0.27951518506446D-07
c(71,20)=   -0.10307834534860D-07
c(72,20)=    0.15107322121888D-07
c(73,20)=    0.33817694479123D-09
c(74,20)=    0.26777201832760D-07
c(75,20)=    0.12375624289778D-07
c(76,20)=    0.30069089333642D-07
c(77,20)=    0.34309026025203D-08
c(78,20)=    0.15330082776951D-07
c(79,20)=    0.14757179279055D-07
c(80,20)=   -0.38003438288289D-08
c(21,21)=    0.66801656361563D-06
c(22,21)=    0.22695041961075D-06
c(23,21)=    0.44090545303246D-06
c(24,21)=   -0.37895350129349D-07
c(25,21)=   -0.68730646069251D-06
c(26,21)=   -0.26279832988690D-06
c(27,21)=    0.28469163064165D-06
c(28,21)=    0.21734973853847D-06
c(29,21)=   -0.17343242283494D-06
c(30,21)=   -0.32809610520172D-07
c(31,21)=   -0.12733505527466D-07
c(32,21)=   -0.67932447668075D-07
c(33,21)=   -0.47019512161784D-07
c(34,21)=    0.31471635002656D-07
c(35,21)=    0.49048388496592D-07
c(36,21)=    0.83329503886604D-07
c(37,21)=   -0.49282623952816D-08
c(38,21)=   -0.12925156961751D-07
c(39,21)=   -0.31963257671986D-07
c(40,21)=    0.66590167520988D-08
c(41,21)=    0.98817058942399D-07
c(42,21)=    0.37605383750701D-07
c(43,21)=   -0.27830621504435D-07
c(44,21)=   -0.40418029820324D-07
c(45,21)=    0.91859461635090D-09
c(46,21)=   -0.19730249807615D-07
c(47,21)=    0.11851628452068D-07
c(48,21)=    0.14639167924279D-07
c(49,21)=   -0.21754507960563D-07
c(50,21)=   -0.24067375212203D-07
c(51,21)=    0.28599563504452D-08
c(52,21)=    0.35681365916782D-07
c(53,21)=    0.97529859432951D-08
c(54,21)=   -0.75983263026241D-07
c(55,21)=    0.29918164273647D-07
c(56,21)=    0.31950356789542D-07
c(57,21)=   -0.18132609646459D-07
c(58,21)=   -0.29257205779927D-07
c(59,21)=   -0.37200365621068D-08
c(60,21)=    0.38205597192254D-07
c(61,21)=    0.11272521808560D-07
c(62,21)=   -0.11664113135659D-07
c(63,21)=   -0.73057985316199D-08
c(64,21)=   -0.27677612090994D-08
c(65,21)=    0.36040573095370D-08
c(66,21)=    0.25787279693229D-07
c(67,21)=   -0.86957805251398D-08
c(68,21)=   -0.55649103460140D-08
c(69,21)=    0.15265548801377D-07
c(70,21)=    0.35448812479815D-07
c(71,21)=   -0.68589306764546D-08
c(72,21)=    0.31834576620702D-07
c(73,21)=   -0.49374229123815D-08
c(74,21)=    0.16439779910660D-09
c(75,21)=    0.31044128392424D-07
c(76,21)=   -0.30113345576158D-08
c(77,21)=    0.39170181976378D-07
c(78,21)=   -0.59917985487534D-08
c(79,21)=    0.23178362595262D-07
c(80,21)=   -0.12786737597268D-07
c(22,22)=   -0.11411133475898D-06
c(23,22)=    0.17110677302269D-06
c(24,22)=   -0.10587573763555D-06
c(25,22)=    0.33521479699086D-06
c(26,22)=    0.40227056541479D-06
c(27,22)=   -0.23272077782554D-06
c(28,22)=   -0.37341707266588D-06
c(29,22)=    0.78100415541802D-07
c(30,22)=    0.20280238119689D-06
c(31,22)=    0.13784175106419D-06
c(32,22)=   -0.12420910000420D-06
c(33,22)=   -0.63631476345334D-07
c(34,22)=   -0.65546742937671D-07
c(35,22)=    0.72321494313189D-07
c(36,22)=    0.51436658429770D-07
c(37,22)=    0.39640401394404D-07
c(38,22)=   -0.10409129234201D-06
c(39,22)=   -0.65354707774874D-07
c(40,22)=    0.42634055748588D-07
c(41,22)=    0.65460387769598D-07
c(42,22)=   -0.18109241916822D-07
c(43,22)=   -0.65632458727538D-07
c(44,22)=    0.68983050838698D-08
c(45,22)=    0.55442014524560D-07
c(46,22)=    0.39443958851056D-07
c(47,22)=   -0.57185005444180D-07
c(48,22)=   -0.34467977713549D-07
c(49,22)=   -0.34596063089222D-07
c(50,22)=    0.13892672880796D-07
c(51,22)=    0.44209988486808D-07
c(52,22)=    0.10896735238304D-07
c(53,22)=   -0.17661974329943D-07
c(54,22)=   -0.53634574743683D-07
c(55,22)=   -0.17248553870907D-07
c(56,22)=    0.30432449984706D-07
c(57,22)=   -0.15796666457769D-07
c(58,22)=   -0.87713858685771D-08
c(59,22)=   -0.20939165130728D-07
c(60,22)=    0.19574575636057D-07
c(61,22)=   -0.30975426084544D-07
c(62,22)=   -0.10973232656097D-07
c(63,22)=    0.40264187905897D-07
c(64,22)=   -0.17368952968311D-07
c(65,22)=    0.43794519119160D-08
c(66,22)=    0.13406416109298D-07
c(67,22)=    0.37085940789819D-07
c(68,22)=    0.10617163289153D-07
c(69,22)=    0.46762767510145D-08
c(70,22)=    0.30719518942792D-07
c(71,22)=   -0.11985554100366D-07
c(72,22)=    0.32837674655834D-07
c(73,22)=    0.18015628046874D-07
c(74,22)=    0.28297092215668D-07
c(75,22)=    0.25034784697646D-07
c(76,22)=    0.34219864341146D-07
c(77,22)=    0.22944444993649D-07
c(78,22)=    0.35818440943536D-07
c(79,22)=    0.21635353557979D-07
c(80,22)=    0.33894259709796D-07
c(23,23)=   -0.55132911635500D-06
c(24,23)=   -0.58926412579959D-06
c(25,23)=   -0.51315458270389D-06
c(26,23)=   -0.26261888686332D-06
c(27,23)=    0.14576740229800D-06
c(28,23)=    0.38585456798044D-06
c(29,23)=    0.66995964161193D-07
c(30,23)=   -0.16251214080634D-06
c(31,23)=   -0.12238386862490D-06
c(32,23)=   -0.42458438711385D-07
c(33,23)=   -0.18067172599805D-07
c(34,23)=    0.13009099109376D-06
c(35,23)=    0.31651928110492D-07
c(36,23)=    0.11978687288709D-09
c(37,23)=   -0.12629933344705D-06
c(38,23)=   -0.37130646665923D-07
c(39,23)=    0.86734371358331D-07
c(40,23)=    0.80611044030297D-07
c(41,23)=   -0.44577951142381D-07
c(42,23)=   -0.54936752461265D-07
c(43,23)=   -0.11297736606126D-06
c(44,23)=   -0.83269881447660D-08
c(45,23)=    0.49991835706267D-07
c(46,23)=   -0.28159240767130D-07
c(47,23)=   -0.68223851181465D-07
c(48,23)=   -0.39901537310458D-07
c(49,23)=    0.48584499142797D-07
c(50,23)=    0.14446762547444D-08
c(51,23)=   -0.50994668094026D-07
c(52,23)=   -0.17651577544801D-07
c(53,23)=    0.15463684721894D-07
c(54,23)=    0.17832090824119D-07
c(55,23)=   -0.42173396292640D-08
c(56,23)=    0.95349773671057D-08
c(57,23)=   -0.52830499167672D-07
c(58,23)=   -0.38825709417983D-07
c(59,23)=    0.17027700586923D-07
c(60,23)=    0.33326808633737D-07
c(61,23)=    0.86978868777464D-08
c(62,23)=   -0.23749175485776D-07
c(63,23)=   -0.30714494159292D-07
c(64,23)=    0.13103887643843D-07
c(65,23)=    0.16155192522717D-07
c(66,23)=    0.32984699484048D-07
c(67,23)=    0.56546884026074D-08
c(68,23)=    0.82254206219378D-09
c(69,23)=   -0.13077802587556D-07
c(70,23)=   -0.41592406679985D-07
c(71,23)=    0.13722687204235D-07
c(72,23)=   -0.39592855730196D-07
c(73,23)=    0.37758246583069D-08
c(74,23)=   -0.24256660738867D-07
c(75,23)=   -0.26000421236925D-07
c(76,23)=    0.83798077844146D-08
c(77,23)=   -0.48513087131603D-07
c(78,23)=    0.15559617982592D-07
c(79,23)=   -0.55293533765553D-07
c(80,23)=    0.20134258862280D-07
c(24,24)=    0.27360023822351D-06
c(25,24)=    0.32738384010408D-06
c(26,24)=    0.72863473401760D-06
c(27,24)=    0.16093057770313D-06
c(28,24)=   -0.29589926512110D-06
c(29,24)=   -0.10991726774138D-06
c(30,24)=    0.11254291701804D-06
c(31,24)=    0.17859073686593D-07
c(32,24)=   -0.30011726079240D-07
c(33,24)=    0.13915220352998D-06
c(34,24)=    0.11409795602060D-06
c(35,24)=   -0.13798663590998D-06
c(36,24)=   -0.13157126113546D-06
c(37,24)=    0.43616254907620D-07
c(38,24)=    0.53242550915175D-07
c(39,24)=   -0.15872770405754D-07
c(40,24)=    0.49378347054287D-07
c(41,24)=   -0.23370839834228D-07
c(42,24)=   -0.61613101099447D-07
c(43,24)=   -0.19440049821244D-07
c(44,24)=    0.10079463392137D-06
c(45,24)=   -0.14693885308818D-07
c(46,24)=   -0.57817890838879D-07
c(47,24)=   -0.57254998897843D-07
c(48,24)=    0.15034522215740D-07
c(49,24)=    0.19077209840090D-07
c(50,24)=   -0.95208092046963D-08
c(51,24)=    0.10282440206743D-07
c(52,24)=   -0.46351326670846D-07
c(53,24)=   -0.10863182835033D-07
c(54,24)=   -0.74793148216455D-09
c(55,24)=    0.93620343855328D-08
c(56,24)=    0.29068723779422D-08
c(57,24)=   -0.15419884256602D-07
c(58,24)=   -0.60395360286350D-08
c(59,24)=    0.39581769608064D-07
c(60,24)=   -0.26508552555905D-08
c(61,24)=   -0.18491414748868D-07
c(62,24)=   -0.17172922202434D-07
c(63,24)=   -0.29424835187963D-08
c(64,24)=    0.19662057722625D-07
c(65,24)=    0.24066937436565D-08
c(66,24)=    0.17953685293943D-07
c(67,24)=    0.95222343380676D-08
c(68,24)=    0.21946265479023D-08
c(69,24)=    0.23364113366357D-08
c(70,24)=   -0.20275669559082D-07
c(71,24)=   -0.29967567672062D-08
c(72,24)=   -0.53160274878373D-08
c(73,24)=    0.31489670100519D-07
c(74,24)=   -0.14624284525145D-08
c(75,24)=    0.28207731919307D-07
c(76,24)=   -0.13192466585347D-07
c(77,24)=    0.35850722943136D-07
c(78,24)=   -0.13251015428247D-07
c(79,24)=    0.29557122352244D-07
c(80,24)=   -0.19331113253332D-07
c(25,25)=    0.13827677523287D-06
c(26,25)=    0.27625280762581D-06
c(27,25)=   -0.96240679400536D-07
c(28,25)=   -0.50080049447013D-07
c(29,25)=    0.10837288324825D-06
c(30,25)=   -0.18429883698663D-06
c(31,25)=    0.32165873057388D-08
c(32,25)=    0.10573172690568D-06
c(33,25)=    0.50146203682854D-07
c(34,25)=   -0.11281578875297D-06
c(35,25)=   -0.54004629034111D-07
c(36,25)=    0.28592437724991D-07
c(37,25)=    0.73009647804798D-07
c(38,25)=   -0.24905250419934D-07
c(39,25)=   -0.42476689752654D-07
c(40,25)=   -0.72060583598888D-07
c(41,25)=    0.20897240094878D-07
c(42,25)=    0.16619890321859D-07
c(43,25)=    0.34629351316775D-07
c(44,25)=   -0.29466152701850D-07
c(45,25)=    0.55075897055936D-07
c(46,25)=   -0.15017910423017D-08
c(47,25)=   -0.34302255994406D-07
c(48,25)=    0.61340681274207D-08
c(49,25)=    0.29036834928664D-07
c(50,25)=   -0.81096146062157D-09
c(51,25)=   -0.42957902037937D-07
c(52,25)=   -0.94697605834800D-08
c(53,25)=   -0.17030074324832D-07
c(54,25)=    0.77264677944891D-08
c(55,25)=   -0.11146747060136D-07
c(56,25)=   -0.94023525814262D-08
c(57,25)=    0.49020310280130D-07
c(58,25)=    0.20124054899184D-07
c(59,25)=   -0.17995425121176D-07
c(60,25)=   -0.56268220366364D-08
c(61,25)=    0.33779675058183D-07
c(62,25)=    0.18273518403616D-07
c(63,25)=   -0.27665906223905D-07
c(64,25)=    0.13771301241879D-07
c(65,25)=   -0.14506820228763D-07
c(66,25)=    0.15878668828642D-08
c(67,25)=    0.17622768093748D-07
c(68,25)=    0.29476695214199D-07
c(69,25)=    0.26577367028634D-07
c(70,25)=    0.11677120879497D-07
c(71,25)=   -0.37909563759775D-07
c(72,25)=    0.25008418617247D-07
c(73,25)=   -0.28380435542435D-08
c(74,25)=    0.10649666747561D-07
c(75,25)=    0.24119782709760D-07
c(76,25)=   -0.54017796573998D-08
c(77,25)=    0.19111003800308D-07
c(78,25)=   -0.31631113712616D-07
c(79,25)=    0.22478056142447D-07
c(80,25)=   -0.39976017400591D-07
c(26,26)=   -0.15768781906645D-06
c(27,26)=   -0.28716714447585D-06
c(28,26)=   -0.37837486086759D-06
c(29,26)=   -0.12189316103849D-07
c(30,26)=    0.11694968385426D-06
c(31,26)=    0.17530693411769D-06
c(32,26)=    0.32497432436030D-07
c(33,26)=   -0.23119735238602D-07
c(34,26)=    0.74636012461474D-07
c(35,26)=    0.15395549829320D-08
c(36,26)=    0.20013793882874D-08
c(37,26)=   -0.42502753577933D-07
c(38,26)=    0.43313630950955D-07
c(39,26)=   -0.17243272467394D-07
c(40,26)=    0.16182982494337D-07
c(41,26)=   -0.24319674509892D-07
c(42,26)=    0.11703675729424D-07
c(43,26)=    0.40833570563849D-07
c(44,26)=   -0.33754590797787D-09
c(45,26)=   -0.15309772857120D-07
c(46,26)=   -0.11040468272219D-07
c(47,26)=    0.39941386102419D-08
c(48,26)=    0.37957000702940D-07
c(49,26)=    0.37206048721622D-07
c(50,26)=    0.33144313356696D-07
c(51,26)=   -0.19472123432831D-08
c(52,26)=   -0.36716371681056D-07
c(53,26)=    0.12402559308799D-07
c(54,26)=    0.40800283942106D-07
c(55,26)=    0.18524878935786D-07
c(56,26)=   -0.27461362896829D-07
c(57,26)=   -0.53435381799615D-07
c(58,26)=    0.74740768845563D-08
c(59,26)=    0.21871523310841D-07
c(60,26)=   -0.19021180048423D-07
c(61,26)=   -0.26836181875340D-08
c(62,26)=   -0.18952070737695D-07
c(63,26)=    0.19662374227625D-07
c(64,26)=    0.53486771190449D-08
c(65,26)=   -0.27400639652929D-08
c(66,26)=   -0.16533488030242D-07
c(67,26)=    0.13897698716283D-07
c(68,26)=   -0.52270927859002D-08
c(69,26)=   -0.15466373293019D-07
c(70,26)=    0.22100777349175D-07
c(71,26)=   -0.18112639846698D-07
c(72,26)=    0.34923579670438D-07
c(73,26)=    0.95812667845168D-08
c(74,26)=    0.14320091388089D-07
c(75,26)=    0.35129192097712D-08
c(76,26)=    0.21728683861324D-07
c(77,26)=   -0.12469401072512D-07
c(78,26)=    0.25525630309300D-07
c(79,26)=   -0.27617997939105D-08
c(80,26)=    0.16711466880933D-07
c(27,27)=    0.49922745314867D-06
c(28,27)=    0.51307032689722D-07
c(29,27)=    0.47342079945331D-06
c(30,27)=    0.31273964601630D-06
c(31,27)=   -0.93208680613626D-07
c(32,27)=   -0.49604553155496D-07
c(33,27)=   -0.52592873912420D-07
c(34,27)=   -0.10971440832705D-06
c(35,27)=   -0.94295297611662D-07
c(36,27)=    0.66563298481375D-07
c(37,27)=    0.14317673509148D-06
c(38,27)=   -0.17144489535359D-07
c(39,27)=   -0.94722990492321D-07
c(40,27)=    0.43786436520367D-08
c(41,27)=    0.12572169813857D-06
c(42,27)=    0.46379869838391D-07
c(43,27)=   -0.46811398187818D-07
c(44,27)=   -0.67977000624105D-07
c(45,27)=   -0.39846423793425D-07
c(46,27)=    0.18625166975218D-07
c(47,27)=    0.20793584765986D-07
c(48,27)=   -0.16229061557653D-07
c(49,27)=   -0.41254092172753D-07
c(50,27)=   -0.11133624841479D-08
c(51,27)=   -0.95346373542353D-08
c(52,27)=    0.38083286492211D-09
c(53,27)=   -0.32630560970930D-09
c(54,27)=   -0.17581841044420D-08
c(55,27)=   -0.40120425280532D-07
c(56,27)=   -0.53321988305913D-08
c(57,27)=    0.62440969854788D-07
c(58,27)=    0.27465497133955D-07
c(59,27)=   -0.29902266610927D-07
c(60,27)=   -0.22972128426252D-07
c(61,27)=    0.33610641602112D-07
c(62,27)=    0.24584322369889D-07
c(63,27)=   -0.24389544807378D-07
c(64,27)=   -0.62052229565562D-07
c(65,27)=    0.20781774132365D-07
c(66,27)=    0.24924741317856D-07
c(67,27)=   -0.40707642790285D-07
c(68,27)=    0.21050431435801D-09
c(69,27)=   -0.16447647711241D-07
c(70,27)=   -0.36497401984733D-08
c(71,27)=    0.11138804780997D-07
c(72,27)=    0.20916302338556D-07
c(73,27)=    0.36258873465065D-07
c(74,27)=    0.14782292380843D-07
c(75,27)=    0.36989493200005D-07
c(76,27)=    0.43421887974999D-08
c(77,27)=    0.41521520832111D-07
c(78,27)=   -0.24515004678780D-08
c(79,27)=    0.52012055408943D-07
c(80,27)=   -0.22821743159857D-07
c(28,28)=   -0.13689973617542D-06
c(29,28)=   -0.36388931964975D-07
c(30,28)=   -0.22231722481905D-06
c(31,28)=   -0.38150827821589D-06
c(32,28)=   -0.58287455217108D-07
c(33,28)=   -0.74350715266180D-07
c(34,28)=    0.81814502729351D-07
c(35,28)=    0.15080464160424D-06
c(36,28)=   -0.75436706016174D-08
c(37,28)=   -0.17476378746928D-06
c(38,28)=   -0.11431195964044D-06
c(39,28)=    0.19545174111152D-06
c(40,28)=    0.16305356191226D-06
c(41,28)=   -0.93078528856659D-07
c(42,28)=   -0.12897344777906D-06
c(43,28)=   -0.77923049190565D-07
c(44,28)=   -0.97268960102326D-08
c(45,28)=    0.74039577562575D-07
c(46,28)=   -0.45808062259145D-08
c(47,28)=   -0.69851857908157D-07
c(48,28)=   -0.38470248635592D-07
c(49,28)=    0.49105275517028D-07
c(50,28)=    0.30297015546051D-07
c(51,28)=   -0.34991597752385D-07
c(52,28)=   -0.51839954778737D-07
c(53,28)=   -0.13214778203740D-07
c(54,28)=    0.22477820102354D-07
c(55,28)=   -0.13908576065961D-07
c(56,28)=   -0.12162546127106D-07
c(57,28)=   -0.15722553589095D-07
c(58,28)=   -0.33487420235316D-07
c(59,28)=    0.47044156737570D-08
c(60,28)=    0.28198960167965D-07
c(61,28)=    0.73303752073192D-08
c(62,28)=   -0.18370041251986D-07
c(63,28)=    0.12349279004626D-08
c(64,28)=    0.35007731749392D-08
c(65,28)=   -0.23224498552509D-07
c(66,28)=    0.13970133302314D-07
c(67,28)=   -0.89120626545945D-09
c(68,28)=   -0.79376188101558D-08
c(69,28)=    0.35387467400759D-08
c(70,28)=   -0.13220441921215D-07
c(71,28)=   -0.32123127538395D-08
c(72,28)=    0.44176267253961D-08
c(73,28)=    0.67659583418974D-08
c(74,28)=    0.14056659707927D-07
c(75,28)=    0.20844639281022D-07
c(76,28)=    0.17272049921388D-07
c(77,28)=    0.22371913041315D-07
c(78,28)=    0.87387875836527D-08
c(79,28)=    0.19380306095852D-07
c(80,28)=    0.24314259472560D-08
c(29,29)=   -0.15397394928722D-06
c(30,29)=    0.17163119431251D-07
c(31,29)=   -0.75620640053598D-07
c(32,29)=    0.42868420782768D-07
c(33,29)=    0.15799589034123D-06
c(34,29)=    0.10232700385077D-06
c(35,29)=   -0.12829305964441D-06
c(36,29)=   -0.14293188494881D-06
c(37,29)=    0.82381952665583D-07
c(38,29)=    0.12901950578175D-06
c(39,29)=   -0.72815490635225D-07
c(40,29)=   -0.72543507067565D-07
c(41,29)=   -0.30558594101122D-07
c(42,29)=    0.48650112950610D-08
c(43,29)=    0.28832565940202D-09
c(44,29)=   -0.20017089623151D-07
c(45,29)=   -0.12708975149658D-07
c(46,29)=   -0.32574403218797D-07
c(47,29)=    0.34663095824221D-08
c(48,29)=    0.31860141170627D-08
c(49,29)=   -0.31449739233846D-08
c(50,29)=   -0.28462300498113D-07
c(51,29)=   -0.35729047384193D-07
c(52,29)=    0.23769432318297D-08
c(53,29)=    0.46161104813680D-07
c(54,29)=   -0.94831705493119D-08
c(55,29)=   -0.56340321229327D-07
c(56,29)=   -0.45657689848635D-07
c(57,29)=    0.29094738968036D-07
c(58,29)=    0.36711216012873D-07
c(59,29)=    0.15306448205221D-08
c(60,29)=   -0.14249087220871D-07
c(61,29)=    0.28649671796700D-07
c(62,29)=   -0.11416285183053D-07
c(63,29)=   -0.35195295471856D-07
c(64,29)=    0.19563969323385D-07
c(65,29)=    0.21661845663115D-07
c(66,29)=    0.10551643257578D-08
c(67,29)=    0.14730370872549D-07
c(68,29)=   -0.13561004099641D-07
c(69,29)=   -0.20476363454595D-07
c(70,29)=   -0.67750299933167D-09
c(71,29)=   -0.35005990687414D-07
c(72,29)=   -0.65643894247056D-08
c(73,29)=    0.21692423049340D-08
c(74,29)=   -0.22267588946526D-07
c(75,29)=   -0.17568510128030D-07
c(76,29)=   -0.25388068203378D-08
c(77,29)=   -0.22367038245313D-07
c(78,29)=    0.27866557535431D-08
c(79,29)=   -0.24300705562891D-07
c(80,29)=    0.11905917970883D-07
c(30,30)=    0.83186343568163D-07
c(31,30)=    0.23812610706667D-07
c(32,30)=    0.16269246219649D-06
c(33,30)=    0.12954915838060D-06
c(34,30)=    0.36293181471303D-07
c(35,30)=    0.10376351794984D-06
c(36,30)=   -0.99078257625299D-07
c(37,30)=   -0.93028550280416D-07
c(38,30)=   -0.12621296536070D-06
c(39,30)=    0.12295135568984D-06
c(40,30)=    0.16096277489877D-06
c(41,30)=   -0.27673871698105D-07
c(42,30)=   -0.10023188068847D-06
c(43,30)=    0.34525964033799D-07
c(44,30)=    0.10836442675597D-06
c(45,30)=    0.42389965745143D-07
c(46,30)=   -0.23482682177555D-07
c(47,30)=   -0.58934570584881D-07
c(48,30)=   -0.36157018585747D-07
c(49,30)=    0.47357700819730D-07
c(50,30)=    0.63047674075323D-07
c(51,30)=   -0.28469935204260D-08
c(52,30)=   -0.27982829983927D-07
c(53,30)=   -0.32706828381456D-07
c(54,30)=    0.21802672309297D-07
c(55,30)=   -0.25004109572926D-07
c(56,30)=    0.10570594075491D-07
c(57,30)=    0.20347369544623D-07
c(58,30)=   -0.13128767510778D-07
c(59,30)=   -0.20356843497471D-07
c(60,30)=   -0.96912610461841D-09
c(61,30)=    0.11291487103604D-07
c(62,30)=    0.67482106273810D-08
c(63,30)=   -0.89895652575622D-08
c(64,30)=    0.39542250836942D-07
c(65,30)=    0.90599837235250D-08
c(66,30)=   -0.22073492422253D-07
c(67,30)=   -0.41919069815473D-07
c(68,30)=    0.28348949353530D-07
c(69,30)=    0.36447692780697D-07
c(70,30)=   -0.70187849639355D-08
c(71,30)=    0.16788401195784D-07
c(72,30)=   -0.23089232450709D-07
c(73,30)=    0.14948522539138D-07
c(74,30)=   -0.93865279422310D-09
c(75,30)=   -0.18807804688485D-07
c(76,30)=    0.36938784201911D-08
c(77,30)=   -0.15504731737758D-07
c(78,30)=   -0.86966446791838D-08
c(79,30)=   -0.37600224133010D-08
c(80,30)=   -0.10561182421300D-07
c(31,31)=   -0.19971160085949D-06
c(32,31)=    0.10908816421837D-06
c(33,31)=   -0.16423185813899D-06
c(34,31)=   -0.27831701396602D-06
c(35,31)=   -0.60764412378810D-07
c(36,31)=   -0.11192741536182D-06
c(37,31)=    0.15141181019967D-06
c(38,31)=    0.30323807895755D-06
c(39,31)=    0.43022525943111D-07
c(40,31)=   -0.21100144591716D-06
c(41,31)=   -0.59251904915924D-07
c(42,31)=    0.94353917088774D-07
c(43,31)=    0.13777426702155D-06
c(44,31)=    0.23134424874716D-08
c(45,31)=   -0.90228954967828D-07
c(46,31)=   -0.12517418053937D-07
c(47,31)=    0.35862897194960D-07
c(48,31)=    0.94098714309386D-07
c(49,31)=    0.42723201776514D-07
c(50,31)=   -0.12902144742638D-07
c(51,31)=   -0.17588938173710D-07
c(52,31)=    0.13165185712010D-07
c(53,31)=    0.61512672049867D-07
c(54,31)=    0.21734143962000D-07
c(55,31)=   -0.32452109194103D-07
c(56,31)=    0.48910497644945D-08
c(57,31)=    0.33392689067592D-07
c(58,31)=   -0.38347685177227D-08
c(59,31)=   -0.44165074507785D-09
c(60,31)=   -0.17841906031136D-07
c(61,31)=   -0.27569789564735D-07
c(62,31)=    0.13380484859283D-07
c(63,31)=   -0.75272826765232D-08
c(64,31)=   -0.45360521108785D-08
c(65,31)=    0.18960524919855D-09
c(66,31)=    0.56356045681231D-08
c(67,31)=    0.30023740484481D-07
c(68,31)=    0.84992532727179D-08
c(69,31)=   -0.59599812310000D-08
c(70,31)=    0.30639656550984D-09
c(71,31)=    0.13947996581212D-07
c(72,31)=    0.16382973370411D-07
c(73,31)=   -0.11097442635374D-08
c(74,31)=    0.28179014724023D-07
c(75,31)=    0.42297861194098D-08
c(76,31)=    0.94844315150209D-08
c(77,31)=    0.24016512447752D-07
c(78,31)=   -0.59973313546418D-08
c(79,31)=    0.26764883203518D-07
c(80,31)=   -0.48425109191681D-08
c(32,32)=    0.87589187448183D-07
c(33,32)=    0.42337514614445D-07
c(34,32)=    0.11983455625456D-06
c(35,32)=    0.22254793177251D-06
c(36,32)=    0.20223389264696D-06
c(37,32)=   -0.17352387443074D-07
c(38,32)=   -0.12309263008740D-06
c(39,32)=   -0.77038893054049D-07
c(40,32)=    0.10648215450483D-06
c(41,32)=    0.89715620444346D-07
c(42,32)=   -0.11317612632298D-07
c(43,32)=   -0.29920121867975D-07
c(44,32)=   -0.26650675104411D-07
c(45,32)=   -0.33637956681839D-07
c(46,32)=    0.22088137286031D-07
c(47,32)=    0.60744666629907D-07
c(48,32)=   -0.16552897447257D-07
c(49,32)=   -0.55388812557025D-07
c(50,32)=    0.30675953304965D-07
c(51,32)=    0.27637435278865D-07
c(52,32)=   -0.78186216052303D-08
c(53,32)=    0.15699836962564D-07
c(54,32)=    0.36391188947904D-07
c(55,32)=   -0.16401826604118D-07
c(56,32)=   -0.70290099949870D-08
c(57,32)=    0.50425652873629D-07
c(58,32)=    0.46358498149771D-08
c(59,32)=   -0.23250058492860D-07
c(60,32)=   -0.33095184207838D-07
c(61,32)=    0.31577510302841D-07
c(62,32)=   -0.17001480715715D-07
c(63,32)=   -0.22013692833670D-07
c(64,32)=   -0.24269775264514D-07
c(65,32)=   -0.39725406675138D-08
c(66,32)=    0.17237680450052D-07
c(67,32)=    0.96163156520467D-09
c(68,32)=   -0.29832574133303D-07
c(69,32)=    0.44200002774657D-07
c(70,32)=   -0.79566261925044D-08
c(71,32)=   -0.13613272106422D-07
c(72,32)=    0.27937227740734D-08
c(73,32)=   -0.93001102663037D-08
c(74,32)=   -0.22116499900137D-07
c(75,32)=   -0.68369804873006D-08
c(76,32)=   -0.22934832888321D-07
c(77,32)=   -0.24642366270291D-07
c(78,32)=   -0.14311989857057D-07
c(79,32)=   -0.29557459108403D-07
c(80,32)=   -0.61874821463050D-08
c(33,33)=    0.58475398205840D-07
c(34,33)=   -0.38256469157865D-06
c(35,33)=   -0.82770977697643D-07
c(36,33)=    0.93631906354654D-07
c(37,33)=   -0.23063691191785D-06
c(38,33)=    0.79175188405789D-07
c(39,33)=    0.15217645575778D-06
c(40,33)=   -0.90179486378835D-07
c(41,33)=   -0.14315529062634D-06
c(42,33)=    0.12455438073450D-06
c(43,33)=    0.13091780059606D-06
c(44,33)=   -0.10659899187999D-06
c(45,33)=   -0.10237156989509D-06
c(46,33)=    0.22138057630585D-07
c(47,33)=    0.46071502629398D-07
c(48,33)=    0.57083405416726D-08
c(49,33)=    0.21112663812878D-07
c(50,33)=   -0.40987304628857D-07
c(51,33)=    0.88179277427948D-08
c(52,33)=    0.20761917487268D-07
c(53,33)=   -0.21153986027314D-07
c(54,33)=   -0.16223554683602D-07
c(55,33)=   -0.13943481594437D-07
c(56,33)=    0.14425864738890D-07
c(57,33)=    0.12394851169921D-07
c(58,33)=   -0.35004298789192D-07
c(59,33)=    0.20642131984743D-07
c(60,33)=    0.14431777555028D-07
c(61,33)=    0.58679802259372D-09
c(62,33)=   -0.63721540402709D-08
c(63,33)=   -0.23478906602072D-07
c(64,33)=   -0.33635494239517D-07
c(65,33)=   -0.50026523423955D-08
c(66,33)=    0.34255035230904D-07
c(67,33)=    0.14400063700947D-07
c(68,33)=    0.26176792756362D-08
c(69,33)=   -0.40449361853318D-08
c(70,33)=   -0.52362825164384D-08
c(71,33)=    0.10143959437298D-07
c(72,33)=   -0.22030297152042D-07
c(73,33)=    0.67502982264843D-08
c(74,33)=    0.98139417034179D-09
c(75,33)=    0.22460790031195D-08
c(76,33)=    0.38068505390521D-08
c(77,33)=   -0.15652074355205D-08
c(78,33)=   -0.21939585571664D-08
c(79,33)=    0.22810602352599D-08
c(80,33)=   -0.45829111183152D-08
c(34,34)=    0.51278731778469D-08
c(35,34)=    0.14088901232734D-06
c(36,34)=   -0.90304181898181D-07
c(37,34)=   -0.13093622779341D-06
c(38,34)=   -0.66317165476431D-07
c(39,34)=   -0.37316803541154D-07
c(40,34)=    0.36525759389135D-07
c(41,34)=    0.14011876832441D-06
c(42,34)=   -0.72015594055707D-08
c(43,34)=   -0.18591326784172D-06
c(44,34)=   -0.48810626144485D-07
c(45,34)=    0.44297320015628D-07
c(46,34)=    0.72841636235943D-07
c(47,34)=    0.96218257796547D-08
c(48,34)=   -0.29381533458598D-07
c(49,34)=   -0.34294171926132D-07
c(50,34)=   -0.22011655861774D-07
c(51,34)=   -0.71549393377881D-08
c(52,34)=    0.33432039816191D-07
c(53,34)=   -0.30479096304569D-07
c(54,34)=   -0.13037337002687D-07
c(55,34)=    0.11644996746542D-07
c(56,34)=   -0.24663480766096D-08
c(57,34)=   -0.15438293532936D-07
c(58,34)=   -0.97851886351454D-08
c(59,34)=    0.14521584137480D-07
c(60,34)=    0.53921040095297D-08
c(61,34)=   -0.11477146333923D-08
c(62,34)=    0.52975867054377D-08
c(63,34)=    0.56175520709770D-08
c(64,34)=    0.10829415430945D-07
c(65,34)=   -0.30524661947070D-07
c(66,34)=    0.77871891153686D-08
c(67,34)=   -0.50818167995811D-09
c(68,34)=    0.30032772156280D-08
c(69,34)=   -0.99777924570013D-09
c(70,34)=    0.98682117429245D-09
c(71,34)=   -0.76516873233813D-08
c(72,34)=   -0.96833533366251D-08
c(73,34)=   -0.10819444526080D-07
c(74,34)=    0.10514207644296D-08
c(75,34)=    0.30582673848055D-08
c(76,34)=    0.39220734316602D-08
c(77,34)=   -0.16301991141856D-07
c(78,34)=    0.23167286395576D-07
c(79,34)=   -0.39324638318927D-07
c(80,34)=    0.36060892185870D-07
c(35,35)=    0.10875503401032D-06
c(36,35)=    0.19141905282068D-06
c(37,35)=    0.21205540234216D-06
c(38,35)=    0.12835760211068D-06
c(39,35)=    0.12759888990090D-06
c(40,35)=   -0.86444568960824D-07
c(41,35)=   -0.19167486373494D-06
c(42,35)=   -0.11961658002801D-06
c(43,35)=    0.11232157498663D-06
c(44,35)=    0.81707928150757D-07
c(45,35)=   -0.34076762837328D-07
c(46,35)=   -0.27400516712106D-07
c(47,35)=    0.34995215911193D-07
c(48,35)=   -0.12038435759851D-07
c(49,35)=   -0.40104793213950D-08
c(50,35)=   -0.16886724287599D-07
c(51,35)=    0.65342450503192D-07
c(52,35)=    0.21523009086677D-07
c(53,35)=   -0.62022321859295D-07
c(54,35)=   -0.28373306059328D-07
c(55,35)=    0.63117998229223D-08
c(56,35)=    0.16876051399852D-07
c(57,35)=    0.25868272424565D-07
c(58,35)=    0.70500786720980D-08
c(59,35)=   -0.12251863005692D-07
c(60,35)=   -0.93566015604113D-08
c(61,35)=   -0.84019414365513D-08
c(62,35)=    0.26742010109669D-07
c(63,35)=    0.76057211425958D-08
c(64,35)=   -0.22491965600775D-07
c(65,35)=   -0.25038314711882D-07
c(66,35)=   -0.34041829506929D-07
c(67,35)=    0.16062129005521D-07
c(68,35)=   -0.10548112606566D-07
c(69,35)=   -0.88880184970768D-08
c(70,35)=    0.31750901623981D-08
c(71,35)=    0.42044182945617D-08
c(72,35)=   -0.11988809888886D-07
c(73,35)=   -0.69820172758214D-10
c(74,35)=    0.91115164241102D-08
c(75,35)=   -0.11609547895818D-07
c(76,35)=    0.11516958940068D-07
c(77,35)=   -0.21074821878740D-07
c(78,35)=    0.34835485227049D-07
c(79,35)=   -0.21629803665303D-07
c(80,35)=    0.34293942445702D-07
c(36,36)=   -0.70916055549077D-07
c(37,36)=   -0.29552977164803D-06
c(38,36)=   -0.10590143858948D-06
c(39,36)=   -0.35195335816752D-07
c(40,36)=   -0.76751593371798D-08
c(41,36)=    0.30079574120527D-07
c(42,36)=    0.70733237972441D-07
c(43,36)=    0.51232856265250D-08
c(44,36)=    0.12026768299596D-07
c(45,36)=    0.42192228750503D-07
c(46,36)=   -0.15831748222013D-07
c(47,36)=   -0.33434991659727D-07
c(48,36)=   -0.42199573302386D-07
c(49,36)=   -0.56358286053357D-08
c(50,36)=    0.23883933279818D-07
c(51,36)=    0.26370392520937D-08
c(52,36)=   -0.15649192586338D-07
c(53,36)=    0.95648409726559D-08
c(54,36)=   -0.15123091187138D-07
c(55,36)=    0.20918159294584D-07
c(56,36)=   -0.70802809683391D-08
c(57,36)=   -0.36492181914188D-07
c(58,36)=   -0.27325960701348D-07
c(59,36)=    0.37477801504421D-07
c(60,36)=   -0.14643991148462D-07
c(61,36)=   -0.95999339767453D-08
c(62,36)=   -0.14701469821659D-07
c(63,36)=    0.18978276808263D-07
c(64,36)=   -0.15637269426543D-07
c(65,36)=    0.54689016064712D-08
c(66,36)=    0.64229928471165D-08
c(67,36)=    0.23849761413833D-07
c(68,36)=   -0.22910566563468D-07
c(69,36)=    0.32871975932697D-07
c(70,36)=    0.12911775666932D-07
c(71,36)=   -0.34242915583168D-08
c(72,36)=    0.11392274884802D-07
c(73,36)=   -0.75910173414701D-08
c(74,36)=   -0.91071088434135D-08
c(75,36)=   -0.20378647678455D-07
c(76,36)=    0.27905963048534D-08
c(77,36)=   -0.16130379168336D-07
c(78,36)=    0.13527489004341D-07
c(79,36)=   -0.25215869054319D-07
c(80,36)=    0.35975273158689D-07
c(37,37)=   -0.85045436427460D-07
c(38,37)=    0.24140798605316D-07
c(39,37)=   -0.62842977823986D-07
c(40,37)=   -0.36861072623085D-07
c(41,37)=    0.72215177795741D-07
c(42,37)=    0.12829880604610D-07
c(43,37)=    0.47597533262804D-07
c(44,37)=    0.73830430107757D-07
c(45,37)=   -0.11811188081236D-06
c(46,37)=   -0.12248817840432D-06
c(47,37)=    0.60659863233175D-07
c(48,37)=    0.84979491220722D-07
c(49,37)=    0.70853251954797D-08
c(50,37)=   -0.16687012678264D-07
c(51,37)=    0.14630482125012D-07
c(52,37)=    0.19236863718465D-07
c(53,37)=   -0.22491318978291D-07
c(54,37)=   -0.95007356401955D-08
c(55,37)=   -0.40253426495364D-08
c(56,37)=   -0.16646635097155D-07
c(57,37)=    0.82754017205748D-08
c(58,37)=    0.13998509257893D-07
c(59,37)=   -0.32988321707420D-07
c(60,37)=   -0.41060253450443D-08
c(61,37)=   -0.29000352325951D-07
c(62,37)=   -0.17074531725920D-07
c(63,37)=    0.11390580203762D-08
c(64,37)=    0.70220783291552D-08
c(65,37)=    0.46659847901810D-08
c(66,37)=   -0.17328699515012D-07
c(67,37)=    0.55327166917024D-09
c(68,37)=    0.53167674136827D-07
c(69,37)=   -0.66684357806224D-08
c(70,37)=    0.40490753477360D-07
c(71,37)=    0.76012426198717D-09
c(72,37)=    0.15572657080339D-07
c(73,37)=   -0.82052508194241D-08
c(74,37)=   -0.71427980308794D-08
c(75,37)=    0.15040299820799D-07
c(76,37)=   -0.69944438257929D-08
c(77,37)=    0.16974340870149D-07
c(78,37)=   -0.24435107865475D-07
c(79,37)=    0.29149038795301D-07
c(80,37)=   -0.35587716469847D-07
c(38,38)=    0.13347014331245D-06
c(39,38)=    0.24047024641936D-06
c(40,38)=    0.10798059966788D-06
c(41,38)=    0.39114089924291D-07
c(42,38)=    0.28125967035620D-07
c(43,38)=   -0.66988229940973D-07
c(44,38)=   -0.11179251119048D-06
c(45,38)=   -0.23819711611369D-07
c(46,38)=    0.42118175811780D-07
c(47,38)=    0.44074049448197D-08
c(48,38)=    0.10270838687362D-07
c(49,38)=   -0.12385049538435D-07
c(50,38)=   -0.89625479508237D-08
c(51,38)=   -0.42945144119609D-07
c(52,38)=   -0.44067169866758D-07
c(53,38)=    0.26040844656555D-07
c(54,38)=   -0.36358336030660D-09
c(55,38)=    0.40094057796099D-07
c(56,38)=    0.19688765599997D-07
c(57,38)=   -0.47801829982398D-07
c(58,38)=   -0.51146761471813D-07
c(59,38)=    0.47886017097438D-07
c(60,38)=    0.27525666149222D-07
c(61,38)=    0.42213139440541D-08
c(62,38)=   -0.17814299187379D-07
c(63,38)=    0.29031052627425D-07
c(64,38)=    0.35779061759296D-07
c(65,38)=   -0.27282827507959D-07
c(66,38)=    0.14753903190742D-07
c(67,38)=   -0.16579971933734D-08
c(68,38)=   -0.57436727657306D-08
c(69,38)=    0.36564824924530D-07
c(70,38)=   -0.32991506827445D-07
c(71,38)=    0.16195339671985D-07
c(72,38)=   -0.35675973184680D-07
c(73,38)=   -0.26513674991594D-07
c(74,38)=   -0.31202448371751D-07
c(75,38)=   -0.48994367964207D-07
c(76,38)=   -0.34750817055747D-07
c(77,38)=   -0.67185675417512D-07
c(78,38)=   -0.11320947600947D-07
c(79,38)=   -0.90812483442622D-07
c(80,38)=    0.36964635550072D-07
c(39,39)=   -0.12498518201365D-07
c(40,39)=   -0.13249419606610D-06
c(41,39)=   -0.51104602604330D-07
c(42,39)=   -0.14888223222037D-07
c(43,39)=   -0.81923363205969D-07
c(44,39)=    0.47012056247530D-07
c(45,39)=    0.94195726890302D-07
c(46,39)=   -0.28033365586129D-07
c(47,39)=    0.72678982196798D-07
c(48,39)=    0.32175051925217D-07
c(49,39)=   -0.54587681320482D-07
c(50,39)=   -0.35906860360574D-07
c(51,39)=    0.35743456986384D-07
c(52,39)=    0.45527074760768D-07
c(53,39)=    0.33476912203896D-07
c(54,39)=    0.21777990528633D-07
c(55,39)=    0.23226980951241D-07
c(56,39)=   -0.43380426291089D-07
c(57,39)=   -0.45781622243592D-08
c(58,39)=    0.54269052639050D-07
c(59,39)=    0.13474235123675D-07
c(60,39)=   -0.31128483256094D-07
c(61,39)=   -0.15497264978552D-07
c(62,39)=   -0.10038525100273D-08
c(63,39)=    0.48634086698810D-08
c(64,39)=    0.95461834588078D-08
c(65,39)=    0.15479454147057D-07
c(66,39)=   -0.14549195991840D-08
c(67,39)=    0.25879413188807D-08
c(68,39)=   -0.41548366603569D-08
c(69,39)=   -0.38736734081212D-07
c(70,39)=    0.82514376192889D-08
c(71,39)=    0.14486363767362D-07
c(72,39)=    0.15291206955104D-07
c(73,39)=    0.13652057707307D-07
c(74,39)=    0.88101634633753D-08
c(75,39)=    0.74584143811033D-10
c(76,39)=    0.54102488893796D-08
c(77,39)=    0.11294842282177D-08
c(78,39)=   -0.16518184430928D-07
c(79,39)=    0.86319238945332D-08
c(80,39)=   -0.21355768760937D-07
c(40,40)=   -0.28282625688147D-06
c(41,40)=    0.71688952427537D-07
c(42,40)=   -0.41047158588585D-07
c(43,40)=   -0.13949634065030D-08
c(44,40)=    0.97186473366003D-07
c(45,40)=    0.12331926439725D-07
c(46,40)=   -0.67038586580410D-07
c(47,40)=   -0.23701097252138D-07
c(48,40)=   -0.27951691934561D-07
c(49,40)=   -0.10094153922806D-08
c(50,40)=    0.63522211342458D-07
c(51,40)=    0.81590490393668D-07
c(52,40)=   -0.52115468198460D-07
c(53,40)=   -0.86155361382134D-07
c(54,40)=    0.84227815420622D-08
c(55,40)=    0.18715156588866D-07
c(56,40)=   -0.67972245257240D-08
c(57,40)=   -0.24573951818749D-07
c(58,40)=   -0.22061955816573D-08
c(59,40)=    0.12130463609119D-07
c(60,40)=    0.21499168174430D-07
c(61,40)=   -0.12609775667312D-07
c(62,40)=    0.11832773538245D-07
c(63,40)=   -0.77427790428475D-08
c(64,40)=    0.19520948878565D-07
c(65,40)=    0.30749231518606D-07
c(66,40)=   -0.10000105914878D-07
c(67,40)=   -0.88965453695869D-08
c(68,40)=   -0.20541193801204D-07
c(69,40)=    0.75723105050778D-08
c(70,40)=   -0.12781738343996D-07
c(71,40)=   -0.24448621383419D-07
c(72,40)=    0.11948259811746D-07
c(73,40)=   -0.21398505614222D-07
c(74,40)=    0.82732733199092D-08
c(75,40)=    0.13044211664951D-07
c(76,40)=   -0.15609725733120D-07
c(77,40)=    0.16822678035570D-07
c(78,40)=   -0.23368023000363D-07
c(79,40)=    0.19885107540687D-07
c(80,40)=   -0.26675411298633D-07
c(41,41)=    0.66269264403178D-07
c(42,41)=    0.99485360418763D-07
c(43,41)=    0.93705529606588D-07
c(44,41)=   -0.19515115240073D-07
c(45,41)=   -0.68356932414673D-07
c(46,41)=   -0.13908385672117D-08
c(47,41)=   -0.28006145007972D-07
c(48,41)=   -0.47709422194968D-08
c(49,41)=   -0.14037285758001D-07
c(50,41)=   -0.15162775455965D-07
c(51,41)=   -0.49284491497797D-07
c(52,41)=   -0.27703139523507D-07
c(53,41)=    0.62703453874591D-07
c(54,41)=    0.16410550192500D-07
c(55,41)=   -0.66846435652823D-07
c(56,41)=   -0.82250945164372D-08
c(57,41)=    0.19896704641883D-08
c(58,41)=    0.22316449818608D-07
c(59,41)=   -0.27565696280757D-07
c(60,41)=   -0.17661700212443D-07
c(61,41)=    0.12602239023716D-07
c(62,41)=   -0.36581893103834D-08
c(63,41)=    0.38730066928903D-09
c(64,41)=   -0.27027514350908D-07
c(65,41)=    0.28227954080531D-07
c(66,41)=    0.25878266572920D-07
c(67,41)=    0.29373762962172D-07
c(68,41)=   -0.82278483119246D-08
c(69,41)=    0.90009074627914D-08
c(70,41)=   -0.10582409007455D-08
c(71,41)=   -0.93034890642953D-08
c(72,41)=    0.10047844986013D-07
c(73,41)=    0.23913638286220D-07
c(74,41)=   -0.90193920166047D-08
c(75,41)=    0.57658148193283D-08
c(76,41)=   -0.67469936719169D-08
c(77,41)=   -0.32940382954878D-08
c(78,41)=   -0.27510094448350D-07
c(79,41)=    0.14175735259075D-07
c(80,41)=   -0.44034849277347D-07
c(42,42)=    0.92084705347499D-07
c(43,42)=   -0.13577783088227D-06
c(44,42)=    0.66054997633427D-07
c(45,42)=   -0.55789097880559D-07
c(46,42)=   -0.90848449752462D-07
c(47,42)=   -0.25582225100060D-08
c(48,42)=    0.11413175256446D-06
c(49,42)=    0.22357274453811D-07
c(50,42)=   -0.62716570436004D-08
c(51,42)=   -0.74847824791414D-09
c(52,42)=   -0.69677744138901D-07
c(53,42)=   -0.49174179537698D-07
c(54,42)=    0.22202981900828D-07
c(55,42)=    0.66096468373953D-07
c(56,42)=   -0.27445455755823D-07
c(57,42)=    0.19616284961736D-08
c(58,42)=    0.28918259082415D-07
c(59,42)=   -0.22802650617594D-07
c(60,42)=   -0.88137084110375D-08
c(61,42)=    0.31929626108723D-07
c(62,42)=   -0.14807107399537D-08
c(63,42)=   -0.26965637758467D-07
c(64,42)=   -0.15718904895240D-07
c(65,42)=    0.30913064719717D-08
c(66,42)=    0.14991231218581D-07
c(67,42)=   -0.42132291844845D-08
c(68,42)=   -0.46684718650143D-08
c(69,42)=   -0.22935192150385D-07
c(70,42)=    0.16358324993681D-08
c(71,42)=   -0.13419448164149D-07
c(72,42)=   -0.25778244251830D-08
c(73,42)=   -0.25467697514425D-08
c(74,42)=    0.18493268392668D-09
c(75,42)=   -0.13303291742363D-07
c(76,42)=   -0.17476449335409D-07
c(77,42)=   -0.29458611799455D-07
c(78,42)=   -0.11589069751728D-07
c(79,42)=   -0.35071126898184D-07
c(80,42)=    0.22849369637846D-08
c(43,43)=   -0.11442187755961D-06
c(44,43)=   -0.18708421375589D-07
c(45,43)=    0.61895517531317D-08
c(46,43)=    0.57517781137047D-07
c(47,43)=    0.13652878643041D-06
c(48,43)=    0.19555034319472D-07
c(49,43)=   -0.17818932260147D-07
c(50,43)=   -0.75669723906277D-07
c(51,43)=   -0.40222254874503D-07
c(52,43)=    0.60811877672993D-07
c(53,43)=    0.11049247406639D-06
c(54,43)=    0.48427088590006D-07
c(55,43)=   -0.78582760780145D-07
c(56,43)=   -0.68780516844109D-07
c(57,43)=   -0.21773151329278D-07
c(58,43)=    0.17147740562778D-07
c(59,43)=   -0.29014997969292D-08
c(60,43)=    0.76737763879791D-08
c(61,43)=    0.16512977570143D-07
c(62,43)=    0.25950586898876D-07
c(63,43)=    0.21419877306605D-07
c(64,43)=   -0.11142692166081D-07
c(65,43)=   -0.15667713090307D-07
c(66,43)=   -0.22924332340415D-07
c(67,43)=    0.26159412835283D-08
c(68,43)=    0.64597941450654D-08
c(69,43)=   -0.74105375701351D-08
c(70,43)=   -0.57090236353885D-08
c(71,43)=   -0.89446579989596D-08
c(72,43)=    0.14337002102460D-07
c(73,43)=   -0.82877772512053D-08
c(74,43)=   -0.81212139359209D-08
c(75,43)=    0.14519859525652D-07
c(76,43)=    0.24246669613919D-08
c(77,43)=    0.18897270806945D-07
c(78,43)=    0.24274786938760D-07
c(79,43)=    0.89619348006897D-08
c(80,43)=    0.19895050638066D-07
c(44,44)=    0.15034955386023D-07
c(45,44)=    0.29116782959401D-07
c(46,44)=   -0.17559317091249D-07
c(47,44)=   -0.18860928689845D-07
c(48,44)=   -0.10326269581702D-07
c(49,44)=   -0.37636204272158D-07
c(50,44)=    0.12467321887837D-07
c(51,44)=    0.93317756497594D-07
c(52,44)=    0.32832347845314D-07
c(53,44)=   -0.42941186913527D-07
c(54,44)=   -0.31413798106672D-07
c(55,44)=   -0.14791229842655D-07
c(56,44)=    0.88408079487054D-08
c(57,44)=   -0.19625091623408D-07
c(58,44)=    0.51543244270834D-07
c(59,44)=   -0.12919504987766D-07
c(60,44)=   -0.13264958891254D-07
c(61,44)=    0.20677079745552D-07
c(62,44)=   -0.65029971790238D-07
c(63,44)=    0.98131242873875D-08
c(64,44)=    0.25244962626010D-07
c(65,44)=   -0.13905280467247D-07
c(66,44)=    0.61066655051075D-07
c(67,44)=   -0.43741884460756D-07
c(68,44)=   -0.19657163308199D-07
c(69,44)=    0.44381263132876D-07
c(70,44)=   -0.10264324063237D-06
c(71,44)=    0.14578568501114D-07
c(72,44)=    0.45166455610180D-07
c(73,44)=   -0.58415182390034D-07
c(74,44)=    0.33722342847193D-07
c(75,44)=    0.23189641608263D-07
c(76,44)=   -0.24921693853061D-07
c(77,44)=   -0.38291464174011D-07
c(78,44)=    0.34358181778215D-07
c(79,44)=   -0.77493084150848D-07
c(80,44)=    0.34141243790682D-07
c(45,45)=    0.14265682580741D-06
c(46,45)=    0.15292837930747D-07
c(47,45)=   -0.21445895294116D-07
c(48,45)=    0.20341197365456D-07
c(49,45)=   -0.89532491837208D-07
c(50,45)=   -0.57507075673080D-07
c(51,45)=   -0.33405316514918D-07
c(52,45)=   -0.41281866726819D-07
c(53,45)=    0.49242933525228D-07
c(54,45)=    0.53303213636083D-07
c(55,45)=   -0.34108191093865D-07
c(56,45)=   -0.62113473636493D-07
c(57,45)=    0.76198176751553D-07
c(58,45)=    0.57042231324263D-07
c(59,45)=   -0.20561091210655D-07
c(60,45)=   -0.32671481821249D-07
c(61,45)=    0.20167512122980D-07
c(62,45)=    0.30195507431502D-07
c(63,45)=    0.95719760972733D-08
c(64,45)=    0.36321603100966D-07
c(65,45)=    0.56831220959255D-08
c(66,45)=   -0.83012494997678D-08
c(67,45)=    0.18359905888740D-07
c(68,45)=    0.28494515715017D-08
c(69,45)=   -0.10647112912967D-07
c(70,45)=   -0.88631691215427D-08
c(71,45)=   -0.21502271205195D-09
c(72,45)=    0.47858427168376D-08
c(73,45)=    0.38933163129512D-07
c(74,45)=    0.25025007992658D-07
c(75,45)=    0.14841842949944D-07
c(76,45)=    0.12748094097083D-07
c(77,45)=    0.17726062828373D-07
c(78,45)=    0.93737960399375D-08
c(79,45)=    0.37839699704908D-07
c(80,45)=    0.47862829358012D-08
c(46,46)=   -0.14124569438264D-06
c(47,46)=    0.12506343805339D-08
c(48,46)=   -0.39402913988663D-07
c(49,46)=    0.64742134407597D-07
c(50,46)=    0.77238528885076D-07
c(51,46)=    0.48117475078108D-07
c(52,46)=    0.65934142496894D-08
c(53,46)=   -0.27167517524648D-07
c(54,46)=   -0.52970215357704D-07
c(55,46)=   -0.89538854076998D-08
c(56,46)=    0.72924442906164D-07
c(57,46)=    0.35820881327144D-07
c(58,46)=   -0.79457083115507D-07
c(59,46)=   -0.88463103236001D-08
c(60,46)=    0.48770163670525D-07
c(61,46)=    0.68898336547070D-08
c(62,46)=    0.56198496682383D-08
c(63,46)=    0.76879764196924D-08
c(64,46)=    0.28406084271576D-07
c(65,46)=    0.51038041852556D-08
c(66,46)=   -0.18592227300313D-07
c(67,46)=    0.31702157827593D-07
c(68,46)=   -0.30109462601346D-08
c(69,46)=   -0.21910454839202D-07
c(70,46)=    0.28039934921293D-07
c(71,46)=   -0.37140052867618D-07
c(72,46)=    0.11943303499421D-08
c(73,46)=   -0.18540140679958D-07
c(74,46)=   -0.30158913772799D-08
c(75,46)=    0.85438958070824D-08
c(76,46)=    0.40976066867798D-08
c(77,46)=    0.32634408061461D-07
c(78,46)=   -0.90090093293770D-08
c(79,46)=    0.68672076020273D-07
c(80,46)=   -0.19362365837174D-07
c(47,47)=    0.14104959906971D-07
c(48,47)=   -0.56545643851864D-07
c(49,47)=   -0.38175920466332D-07
c(50,47)=   -0.76892238928920D-07
c(51,47)=   -0.78104320931309D-08
c(52,47)=    0.19157556304477D-07
c(53,47)=    0.20290009300674D-07
c(54,47)=    0.29724165786259D-07
c(55,47)=   -0.31083742120195D-08
c(56,47)=   -0.43551560493261D-07
c(57,47)=    0.77487657315157D-08
c(58,47)=    0.34053742630077D-07
c(59,47)=   -0.34855336573877D-09
c(60,47)=   -0.30695640747836D-07
c(61,47)=    0.31835869509673D-08
c(62,47)=    0.36686279312462D-08
c(63,47)=   -0.22215743048625D-07
c(64,47)=   -0.50574708035914D-07
c(65,47)=   -0.20035066392557D-07
c(66,47)=    0.81522797029569D-08
c(67,47)=    0.20022103306298D-07
c(68,47)=    0.15718793088804D-08
c(69,47)=   -0.28678333511384D-08
c(70,47)=   -0.26607709402548D-07
c(71,47)=    0.41520598257674D-08
c(72,47)=   -0.11803899095500D-07
c(73,47)=    0.36826455634889D-08
c(74,47)=    0.58832581399847D-08
c(75,47)=    0.20245747297436D-07
c(76,47)=    0.96155820778275D-08
c(77,47)=    0.17867757804120D-07
c(78,47)=    0.34855665927062D-07
c(79,47)=    0.28893007548978D-07
c(80,47)=    0.38683398105126D-07
c(48,48)=    0.13183786875222D-06
c(49,48)=    0.10087023421292D-07
c(50,48)=    0.81723772372626D-07
c(51,48)=    0.26523736965109D-07
c(52,48)=   -0.75389297257786D-07
c(53,48)=   -0.63839178227448D-07
c(54,48)=   -0.27153787256252D-07
c(55,48)=   -0.13374928239397D-07
c(56,48)=    0.25304555690120D-07
c(57,48)=    0.26820681196692D-08
c(58,48)=   -0.49324278938910D-07
c(59,48)=   -0.31186651188467D-07
c(60,48)=    0.33338897659264D-07
c(61,48)=    0.20864299491934D-07
c(62,48)=    0.95898773671841D-08
c(63,48)=   -0.12398556095698D-07
c(64,48)=   -0.10627481591800D-07
c(65,48)=   -0.92468866785723D-08
c(66,48)=   -0.94963352921816D-08
c(67,48)=   -0.20438806970465D-08
c(68,48)=   -0.17983398277993D-07
c(69,48)=   -0.30340242280293D-07
c(70,48)=   -0.21676501819313D-07
c(71,48)=   -0.35822028310007D-07
c(72,48)=    0.14350716357890D-07
c(73,48)=   -0.43515657192583D-08
c(74,48)=   -0.12839602474258D-07
c(75,48)=    0.52903331325952D-08
c(76,48)=   -0.12157397860099D-07
c(77,48)=    0.87611039507312D-08
c(78,48)=   -0.27391677973795D-07
c(79,48)=    0.95069694460083D-08
c(80,48)=   -0.25340201833080D-07
c(49,49)=   -0.83893778506803D-07
c(50,49)=    0.27064741304335D-07
c(51,49)=    0.59735502015113D-09
c(52,49)=    0.30218896960283D-07
c(53,49)=    0.78701606316708D-07
c(54,49)=   -0.46582548862016D-08
c(55,49)=   -0.34068124375426D-09
c(56,49)=    0.48514433884466D-07
c(57,49)=   -0.16947411657877D-07
c(58,49)=   -0.54934651340397D-08
c(59,49)=    0.24107767614798D-07
c(60,49)=   -0.25104549570907D-08
c(61,49)=   -0.57972944552729D-07
c(62,49)=    0.42758974997715D-07
c(63,49)=   -0.64445010884034D-08
c(64,49)=    0.86496732192467D-08
c(65,49)=    0.11232181522790D-07
c(66,49)=   -0.37864499506033D-07
c(67,49)=    0.19317408091003D-07
c(68,49)=    0.17397051050412D-07
c(69,49)=    0.80501430900921D-08
c(70,49)=    0.16955407897328D-07
c(71,49)=    0.38286819337312D-07
c(72,49)=   -0.10128367419507D-07
c(73,49)=    0.20312064206339D-07
c(74,49)=    0.24389096991361D-07
c(75,49)=    0.19831553378102D-07
c(76,49)=    0.16471155022668D-07
c(77,49)=    0.27719829269071D-07
c(78,49)=    0.27858782388775D-07
c(79,49)=   -0.10207394251253D-07
c(80,49)=    0.37372005385793D-07
c(50,50)=   -0.78655304961773D-07
c(51,50)=   -0.76718587761635D-07
c(52,50)=   -0.89359338973498D-07
c(53,50)=   -0.63253751970261D-07
c(54,50)=    0.96908722492314D-08
c(55,50)=    0.25257480614650D-08
c(56,50)=    0.65906394738951D-08
c(57,50)=   -0.64046436732057D-08
c(58,50)=   -0.24866048819383D-07
c(59,50)=   -0.56934326330242D-08
c(60,50)=    0.23431408908595D-08
c(61,50)=    0.42512188872396D-07
c(62,50)=    0.84098611284314D-08
c(63,50)=   -0.21477523224098D-07
c(64,50)=   -0.21031569730425D-07
c(65,50)=    0.41137055754030D-07
c(66,50)=   -0.21919227094243D-07
c(67,50)=   -0.57077805469915D-08
c(68,50)=   -0.11200861055523D-07
c(69,50)=   -0.31850325078367D-07
c(70,50)=   -0.39860295573898D-08
c(71,50)=   -0.39638944627352D-07
c(72,50)=    0.27386483149971D-07
c(73,50)=   -0.33781222070150D-08
c(74,50)=    0.99179260059346D-08
c(75,50)=    0.14836368757641D-08
c(76,50)=    0.73951694721136D-08
c(77,50)=    0.63976287061156D-08
c(78,50)=   -0.37047176831862D-08
c(79,50)=    0.12257327135528D-07
c(80,50)=   -0.65840102450554D-09
c(51,51)=    0.71649831664665D-07
c(52,51)=    0.16063833192547D-07
c(53,51)=    0.31663668748921D-07
c(54,51)=   -0.96481014932878D-09
c(55,51)=   -0.15373391967241D-07
c(56,51)=    0.31492453053932D-07
c(57,51)=   -0.11570563810668D-09
c(58,51)=   -0.28931597834427D-08
c(59,51)=    0.41435796906410D-07
c(60,51)=    0.13070549173929D-07
c(61,51)=   -0.56997555582700D-07
c(62,51)=   -0.31280643852365D-07
c(63,51)=    0.39321121551819D-07
c(64,51)=   -0.88565493270598D-08
c(65,51)=   -0.22640084711286D-07
c(66,51)=    0.18325345120306D-07
c(67,51)=    0.18139330298322D-07
c(68,51)=    0.20664491841579D-07
c(69,51)=    0.37650603546457D-07
c(70,51)=    0.32543754149647D-07
c(71,51)=    0.26063198733519D-07
c(72,51)=    0.32848204100458D-07
c(73,51)=    0.16593720400871D-07
c(74,51)=    0.22770517219112D-07
c(75,51)=    0.21571734230323D-07
c(76,51)=   -0.44162812278580D-08
c(77,51)=    0.28007202452446D-07
c(78,51)=   -0.22122967357936D-07
c(79,51)=    0.19328668658253D-07
c(80,51)=   -0.29308862022510D-07
c(52,52)=   -0.11917850450220D-07
c(53,52)=    0.10178324076334D-06
c(54,52)=    0.58929881779205D-07
c(55,52)=    0.34829088506542D-07
c(56,52)=    0.16503415265138D-08
c(57,52)=   -0.32575071175232D-08
c(58,52)=    0.47579921152148D-08
c(59,52)=    0.11888762468345D-07
c(60,52)=   -0.17557956528633D-07
c(61,52)=   -0.10272219459391D-07
c(62,52)=   -0.56820460819980D-08
c(63,52)=   -0.35737699480519D-07
c(64,52)=   -0.16699373339017D-07
c(65,52)=   -0.20041144180385D-07
c(66,52)=    0.46210866588471D-07
c(67,52)=   -0.31160658888733D-07
c(68,52)=    0.12946214273948D-07
c(69,52)=   -0.23434559023731D-07
c(70,52)=   -0.83481726529725D-08
c(71,52)=    0.92231887371209D-08
c(72,52)=    0.31895493452823D-08
c(73,52)=    0.22255625341153D-07
c(74,52)=    0.21680626238655D-07
c(75,52)=    0.26803795127612D-08
c(76,52)=    0.11705238582242D-07
c(77,52)=    0.28329758737764D-07
c(78,52)=    0.21223322726206D-07
c(79,52)=    0.16061455486414D-07
c(80,52)=    0.22122360251331D-07
c(53,53)=   -0.68906317081098D-07
c(54,53)=   -0.60139427876214D-07
c(55,53)=    0.15391096176413D-07
c(56,53)=   -0.33805643074270D-07
c(57,53)=   -0.50938541055426D-07
c(58,53)=    0.16348008538981D-07
c(59,53)=    0.49039386929627D-08
c(60,53)=    0.44413621448349D-08
c(61,53)=   -0.32225383168232D-07
c(62,53)=   -0.89356017798621D-08
c(63,53)=    0.94239635687397D-08
c(64,53)=   -0.87572493735644D-08
c(65,53)=    0.13209339403168D-07
c(66,53)=    0.27687647322889D-09
c(67,53)=    0.11123610748244D-07
c(68,53)=    0.38325877465571D-07
c(69,53)=    0.27840788777681D-07
c(70,53)=   -0.88360645883774D-08
c(71,53)=    0.13562389883347D-07
c(72,53)=   -0.66732701659572D-08
c(73,53)=    0.41933879004779D-08
c(74,53)=    0.27542979980989D-08
c(75,53)=   -0.26217514680132D-08
c(76,53)=    0.10740428081553D-08
c(77,53)=   -0.25936431338946D-07
c(78,53)=    0.12336307137725D-08
c(79,53)=   -0.58450093285968D-07
c(80,53)=    0.23015764407575D-07
c(54,54)=    0.25146386421733D-07
c(55,54)=   -0.41312492468982D-07
c(56,54)=    0.13899542785488D-07
c(57,54)=   -0.27243615989484D-07
c(58,54)=    0.16963677848444D-07
c(59,54)=    0.30057873584805D-07
c(60,54)=    0.90897022651945D-08
c(61,54)=    0.14362518907253D-07
c(62,54)=    0.21614788154163D-07
c(63,54)=    0.67266755587643D-08
c(64,54)=    0.80274054593026D-08
c(65,54)=   -0.14779647252088D-07
c(66,54)=    0.26377341781642D-08
c(67,54)=    0.13119241916865D-07
c(68,54)=   -0.36825921062611D-07
c(69,54)=    0.15505522683942D-07
c(70,54)=    0.15039615228945D-07
c(71,54)=   -0.16702220423832D-08
c(72,54)=   -0.15348632944715D-07
c(73,54)=   -0.24829612157751D-07
c(74,54)=    0.78502933547880D-08
c(75,54)=   -0.67078194155254D-08
c(76,54)=    0.14748767487734D-07
c(77,54)=   -0.13101636965934D-07
c(78,54)=    0.13844005003076D-07
c(79,54)=    0.72503899731971D-08
c(80,54)=    0.25821120091196D-07
c(55,55)=   -0.98051494924926D-08
c(56,55)=    0.14126817716778D-06
c(57,55)=   -0.13467101173673D-07
c(58,55)=    0.82297465056759D-08
c(59,55)=    0.29225038907521D-07
c(60,55)=   -0.64093659413411D-07
c(61,55)=   -0.47913699208843D-07
c(62,55)=   -0.10389687527535D-07
c(63,55)=    0.49810414202377D-08
c(64,55)=   -0.22994709008642D-07
c(65,55)=    0.17242577871404D-07
c(66,55)=    0.21061388375067D-07
c(67,55)=   -0.68768515944441D-08
c(68,55)=    0.45347017660204D-08
c(69,55)=    0.47904452762682D-07
c(70,55)=    0.18448559507244D-07
c(71,55)=   -0.21244351613833D-07
c(72,55)=   -0.13105884046395D-08
c(73,55)=    0.33496923412302D-08
c(74,55)=   -0.13093431332956D-07
c(75,55)=   -0.55420200956180D-08
c(76,55)=   -0.16499505957784D-07
c(77,55)=    0.20582488512699D-08
c(78,55)=   -0.44222802265208D-07
c(79,55)=   -0.77237356795639D-08
c(80,55)=   -0.45631690754147D-07
c(56,56)=   -0.40123407345616D-07
c(57,56)=   -0.48078870485104D-07
c(58,56)=    0.19069834981467D-07
c(59,56)=    0.26403983071577D-07
c(60,56)=   -0.20851446900971D-07
c(61,56)=   -0.38687964516932D-07
c(62,56)=    0.40027130893790D-07
c(63,56)=    0.56000185804383D-07
c(64,56)=   -0.47512071486663D-08
c(65,56)=    0.29803130463985D-08
c(66,56)=   -0.11097075764010D-07
c(67,56)=   -0.98398333106898D-09
c(68,56)=   -0.89812932768307D-08
c(69,56)=   -0.27444091868882D-07
c(70,56)=    0.71394078207807D-08
c(71,56)=    0.80232100103553D-08
c(72,56)=    0.88993436969500D-08
c(73,56)=    0.15966538088772D-07
c(74,56)=    0.38662741665250D-08
c(75,56)=    0.89148213228963D-08
c(76,56)=   -0.13968865312220D-07
c(77,56)=   -0.44765328926921D-08
c(78,56)=   -0.19401767715379D-07
c(79,56)=   -0.39704639568240D-08
c(80,56)=   -0.22012530731357D-07
c(57,57)=    0.63662142430206D-07
c(58,57)=   -0.61005760103547D-07
c(59,57)=   -0.34773117369159D-07
c(60,57)=   -0.13499097103984D-07
c(61,57)=    0.17126155435825D-07
c(62,57)=    0.64963959932636D-07
c(63,57)=    0.54259025484029D-07
c(64,57)=   -0.52186134100488D-07
c(65,57)=   -0.22130732676743D-07
c(66,57)=   -0.54542285960446D-08
c(67,57)=   -0.12538641918910D-07
c(68,57)=    0.14420623833334D-07
c(69,57)=    0.21102228863191D-07
c(70,57)=    0.26963740902659D-08
c(71,57)=   -0.16541292896987D-07
c(72,57)=    0.13496408458053D-07
c(73,57)=    0.47681230402104D-08
c(74,57)=    0.96128372730776D-08
c(75,57)=    0.40469998860768D-08
c(76,57)=    0.25058187955984D-07
c(77,57)=    0.48739736631978D-08
c(78,57)=    0.85817486067055D-08
c(79,57)=    0.91923573652914D-08
c(80,57)=    0.21580207304601D-07
c(58,58)=    0.30404168717792D-07
c(59,58)=    0.36637732444766D-07
c(60,58)=   -0.48187228617950D-07
c(61,58)=   -0.18450902782091D-07
c(62,58)=    0.23639019012449D-07
c(63,58)=   -0.31789808240154D-07
c(64,58)=   -0.51495679858299D-07
c(65,58)=   -0.22612778319004D-07
c(66,58)=    0.33350776016663D-08
c(67,58)=   -0.30328810291125D-08
c(68,58)=   -0.12511875027624D-07
c(69,58)=   -0.18952692994231D-07
c(70,58)=    0.16125248447427D-07
c(71,58)=   -0.12079424125372D-07
c(72,58)=    0.11802305307720D-07
c(73,58)=   -0.98850212842383D-08
c(74,58)=    0.64901372670500D-08
c(75,58)=   -0.65035971316783D-10
c(76,58)=    0.12462278158387D-08
c(77,58)=   -0.11518357795713D-07
c(78,58)=    0.17927360960238D-07
c(79,58)=   -0.20106153109405D-07
c(80,58)=    0.23847471373295D-07
c(59,59)=   -0.54197484848084D-07
c(60,59)=    0.13195419173402D-07
c(61,59)=    0.45503356796075D-07
c(62,59)=   -0.18488950216058D-07
c(63,59)=   -0.13960575345108D-07
c(64,59)=   -0.36197787334473D-07
c(65,59)=   -0.33970411942424D-07
c(66,59)=    0.31570141904188D-07
c(67,59)=   -0.22622266641646D-07
c(68,59)=   -0.23545771401746D-07
c(69,59)=    0.26108287400461D-07
c(70,59)=   -0.10133573729510D-07
c(71,59)=   -0.25833239151630D-08
c(72,59)=   -0.13609862257131D-07
c(73,59)=    0.12741460641084D-08
c(74,59)=   -0.48779453360635D-07
c(75,59)=    0.34181500763575D-08
c(76,59)=   -0.11496318806930D-07
c(77,59)=    0.18047284149487D-07
c(78,59)=    0.17360033625305D-07
c(79,59)=   -0.93953597834784D-10
c(80,59)=    0.38800928176602D-07
c(60,60)=    0.28593707509687D-07
c(61,60)=   -0.75616956756271D-07
c(62,60)=   -0.55700540200929D-08
c(63,60)=   -0.22091231827130D-07
c(64,60)=    0.28948480397681D-07
c(65,60)=    0.62989187973610D-07
c(66,60)=   -0.71384896383999D-08
c(67,60)=   -0.45844593434291D-07
c(68,60)=   -0.19217972358172D-07
c(69,60)=   -0.19635478036358D-07
c(70,60)=   -0.11957403033925D-07
c(71,60)=   -0.12261960577181D-07
c(72,60)=    0.65877363610953D-08
c(73,60)=   -0.17113095123597D-07
c(74,60)=   -0.27127106512783D-07
c(75,60)=    0.83263511549959D-08
c(76,60)=    0.29824364947494D-08
c(77,60)=    0.19921028326643D-07
c(78,60)=    0.46332146449539D-07
c(79,60)=    0.35514037943076D-07
c(80,60)=    0.48865245673134D-07
c(61,61)=    0.52319776819198D-07
c(62,61)=    0.15023349539380D-07
c(63,61)=   -0.70027422021410D-07
c(64,61)=   -0.34638400597222D-07
c(65,61)=    0.14985893846734D-08
c(66,61)=   -0.29791974365691D-07
c(67,61)=   -0.90679449958987D-08
c(68,61)=    0.12857956782801D-07
c(69,61)=    0.53845517212863D-07
c(70,61)=   -0.17458835856659D-07
c(71,61)=    0.36191908982051D-07
c(72,61)=    0.25837294831375D-07
c(73,61)=   -0.12606088590518D-07
c(74,61)=   -0.25171077476141D-08
c(75,61)=   -0.35504806941465D-07
c(76,61)=    0.21289043925316D-10
c(77,61)=   -0.66199539992973D-07
c(78,61)=   -0.93599768844641D-08
c(79,61)=   -0.96226191599110D-07
c(80,61)=   -0.32834066503532D-07
c(62,62)=   -0.60237139592494D-08
c(63,62)=   -0.20625629252480D-08
c(64,62)=    0.30807889241947D-07
c(65,62)=    0.19055891691252D-07
c(66,62)=    0.23197555972264D-07
c(67,62)=   -0.41975845620426D-07
c(68,62)=   -0.49644119375293D-09
c(69,62)=   -0.51833992488221D-08
c(70,62)=    0.67550663105708D-08
c(71,62)=    0.10323538786660D-07
c(72,62)=    0.15241950767748D-07
c(73,62)=    0.23488887949895D-08
c(74,62)=   -0.74696335029352D-08
c(75,62)=   -0.10698619291426D-07
c(76,62)=   -0.28652339175347D-07
c(77,62)=    0.24374250151510D-08
c(78,62)=   -0.33061642688537D-07
c(79,62)=    0.14350862995235D-07
c(80,62)=   -0.51276446311965D-07
c(63,63)=   -0.22756202330126D-07
c(64,63)=   -0.66785580121160D-08
c(65,63)=    0.14062329305521D-07
c(66,63)=    0.57591865682234D-07
c(67,63)=    0.38626292276263D-07
c(68,63)=    0.15931686122385D-07
c(69,63)=    0.19912556468545D-07
c(70,63)=    0.39780213948434D-07
c(71,63)=   -0.64204764298785D-07
c(72,63)=    0.17881688771298D-07
c(73,63)=   -0.26418452436418D-07
c(74,63)=    0.18081288199856D-07
c(75,63)=    0.26357469243573D-07
c(76,63)=    0.35714299208398D-07
c(77,63)=    0.31287884421194D-07
c(78,63)=    0.19831070181322D-07
c(79,63)=    0.71195365286147D-07
c(80,63)=   -0.17057535049470D-07
c(64,64)=    0.16001141443898D-07
c(65,64)=    0.26376404806720D-07
c(66,64)=   -0.44667259109842D-07
c(67,64)=   -0.33013127101211D-07
c(68,64)=   -0.24688575424350D-07
c(69,64)=   -0.16536702551016D-07
c(70,64)=    0.33616196659677D-07
c(71,64)=    0.82053664939067D-08
c(72,64)=    0.53768906241930D-08
c(73,64)=    0.22133367848765D-07
c(74,64)=    0.34709204101415D-07
c(75,64)=    0.11497386113963D-07
c(76,64)=    0.28299933040159D-07
c(77,64)=    0.42666219442628D-07
c(78,64)=    0.89976141970644D-08
c(79,64)=    0.61504928080787D-07
c(80,64)=    0.46337885997585D-08
c(65,65)=   -0.37923340594247D-07
c(66,65)=    0.57908038036215D-07
c(67,65)=    0.39678249190651D-07
c(68,65)=   -0.32980329401940D-08
c(69,65)=   -0.12573178107348D-07
c(70,65)=    0.21840245123388D-07
c(71,65)=   -0.15649518181756D-07
c(72,65)=    0.37446874360709D-07
c(73,65)=   -0.36013087655695D-07
c(74,65)=    0.26242447015256D-07
c(75,65)=    0.65993616223791D-08
c(76,65)=   -0.89022251743086D-08
c(77,65)=    0.25086514838982D-07
c(78,65)=   -0.13654570164780D-07
c(79,65)=    0.41828803742268D-07
c(80,65)=   -0.21721222732775D-07
c(66,66)=   -0.42645750642653D-07
c(67,66)=   -0.35890545057486D-07
c(68,66)=    0.48737014786927D-07
c(69,66)=    0.89928827974056D-08
c(70,66)=    0.45315771776624D-08
c(71,66)=   -0.67900004779097D-08
c(72,66)=   -0.22716964517560D-07
c(73,66)=   -0.17907817841494D-07
c(74,66)=   -0.22089058870163D-07
c(75,66)=   -0.22542089210875D-07
c(76,66)=   -0.24111771883861D-07
c(77,66)=   -0.32975407909178D-07
c(78,66)=    0.18032489207429D-07
c(79,66)=   -0.67365391398799D-07
c(80,66)=    0.27165847045611D-07
c(67,67)=    0.13746247493082D-07
c(68,67)=   -0.46840571421664D-07
c(69,67)=   -0.46082937833652D-07
c(70,67)=   -0.17875821555122D-07
c(71,67)=   -0.13324915725140D-07
c(72,67)=   -0.34047797587000D-07
c(73,67)=   -0.46711868852491D-07
c(74,67)=   -0.35449274763391D-08
c(75,67)=   -0.20397141981446D-07
c(76,67)=   -0.11108360493265D-07
c(77,67)=    0.11381321576696D-07
c(78,67)=   -0.63488647439425D-08
c(79,67)=   -0.48538968483360D-07
c(80,67)=   -0.71983453260028D-07
c(68,68)=    0.16680960733166D-07
c(69,68)=   -0.52083592951271D-08
c(70,68)=    0.32325341968542D-07
c(71,68)=   -0.16496832468463D-07
c(72,68)=   -0.20643422798399D-09
c(73,68)=    0.21679139252239D-07
c(74,68)=    0.39106999354482D-07
c(75,68)=    0.18143968667046D-07
c(76,68)=   -0.68931580196856D-08
c(77,68)=    0.33085534802760D-07
c(78,68)=   -0.11037676026853D-07
c(79,68)=    0.31808918194418D-07
c(80,68)=   -0.10952551120856D-07
c(69,69)=   -0.12245535537515D-07
c(70,69)=    0.10287607959453D-07
c(71,69)=    0.48907175215348D-07
c(72,69)=   -0.66055000909424D-08
c(73,69)=   -0.18998722014041D-07
c(74,69)=    0.50721545495667D-07
c(75,69)=    0.13625539452155D-07
c(76,69)=    0.20921948458561D-07
c(77,69)=   -0.80026826707920D-08
c(78,69)=    0.10498380995640D-08
c(79,69)=   -0.26003036261122D-08
c(80,69)=   -0.18329141763226D-07
c(70,70)=    0.26129165325730D-08
c(71,70)=    0.39445680827134D-07
c(72,70)=   -0.36440945745060D-07
c(73,70)=   -0.34304469545119D-07
c(74,70)=    0.90666251366812D-08
c(75,70)=    0.98845829432444D-08
c(76,70)=    0.20848195762778D-07
c(77,70)=    0.65496386810649D-08
c(78,70)=    0.13420363140932D-07
c(79,70)=    0.24097504159691D-07
c(80,70)=   -0.61345329316370D-07
c(71,71)=   -0.23979906811712D-07
c(72,71)=    0.81589711992019D-09
c(73,71)=   -0.12038754975651D-07
c(74,71)=   -0.33769050389210D-08
c(75,71)=    0.22626329292769D-07
c(76,71)=    0.34705206656060D-08
c(77,71)=   -0.47557144208902D-08
c(78,71)=   -0.12837520837664D-07
c(79,71)=   -0.47143600956230D-08
c(80,71)=   -0.13673596061185D-07
c(72,72)=    0.20935785099319D-07
c(73,72)=    0.15457694715164D-07
c(74,72)=    0.43118438437429D-07
c(75,72)=    0.44347445470690D-08
c(76,72)=    0.86935013976379D-08
c(77,72)=    0.14238230603767D-07
c(78,72)=   -0.75862730206882D-08
c(79,72)=    0.18694963412678D-07
c(80,72)=   -0.22644291755096D-08
c(73,73)=    0.10665385729218D-07
c(74,73)=    0.10893726176099D-07
c(75,73)=   -0.39069862760945D-07
c(76,73)=   -0.68738761190501D-07
c(77,73)=   -0.25753617667176D-07
c(78,73)=   -0.42049642957662D-07
c(79,73)=    0.96057051523230D-08
c(80,73)=   -0.54754957010108D-09
c(74,74)=    0.35533577546674D-07
c(75,74)=    0.24381875264283D-08
c(76,74)=   -0.16956003331522D-07
c(77,74)=    0.14996059285111D-07
c(78,74)=    0.15017114043667D-07
c(79,74)=    0.28157653772661D-07
c(80,74)=   -0.10857895553919D-07
c(75,75)=   -0.20557935405623D-07
c(76,75)=   -0.11445288933741D-07
c(77,75)=   -0.13998374279198D-07
c(78,75)=   -0.86094853130988D-07
c(79,75)=   -0.44073733643326D-07
c(80,75)=   -0.97888515027293D-07
c(76,76)=    0.16947695586139D-07
c(77,76)=    0.48271839875412D-07
c(78,76)=    0.87890863155167D-07
c(79,76)=   -0.47066652207509D-07
c(80,76)=   -0.18420259011703D-07
c(77,77)=    0.22222570253966D-07
c(78,77)=   -0.35567341235073D-07
c(79,77)=    0.15074521892097D-07
c(80,77)=   -0.65806552524787D-07
c(78,78)=    0.42469008564636D-07
c(79,78)=   -0.75948714508049D-08
c(80,78)=   -0.16243257517502D-07
c(79,79)=   -0.12617301383201D-07
c(80,79)=    0.16728015999556D-07
c(80,80)=    0.54813675815846D-08
s( 2, 1)=   -0.11401847225989D-09
s( 3, 1)=    0.25155150744332D-04
s( 4, 1)=    0.37588917182210D-05
s( 5, 1)=    0.21199078107223D-05
s( 6, 1)=   -0.15154712831579D-05
s( 7, 1)=   -0.22693254012337D-06
s( 8, 1)=    0.74869335197712D-06
s( 9, 1)=   -0.48887826939649D-06
s(10, 1)=    0.22250477035440D-06
s(11, 1)=   -0.22078265839368D-06
s(12, 1)=   -0.51290190208753D-07
s(13, 1)=    0.22519671490044D-06
s(14, 1)=    0.30949508730742D-06
s(15, 1)=    0.25961987086161D-06
s(16, 1)=   -0.55346995650064D-06
s(17, 1)=   -0.63687648573341D-06
s(18, 1)=   -0.57862275066301D-06
s(19, 1)=    0.50580889644975D-06
s(20, 1)=    0.55159048065934D-06
s(21, 1)=   -0.42466628430160D-06
s(22, 1)=   -0.11861574241659D-06
s(23, 1)=    0.74492590496240D-07
s(24, 1)=    0.16898337729283D-06
s(25, 1)=   -0.16782838780194D-07
s(26, 1)=   -0.38583254338416D-06
s(27, 1)=   -0.42826946644146D-07
s(28, 1)=    0.10245418978571D-06
s(29, 1)=    0.20097501534675D-06
s(30, 1)=   -0.18646620419890D-08
s(31, 1)=   -0.82686085187458D-07
s(32, 1)=   -0.18316285527029D-06
s(33, 1)=   -0.93637017303617D-07
s(34, 1)=   -0.51407090820137D-07
s(35, 1)=    0.15318194346565D-06
s(36, 1)=   -0.28160803509340D-07
s(37, 1)=   -0.18555020643732D-06
s(38, 1)=   -0.12915239202711D-07
s(39, 1)=    0.11449391594282D-06
s(40, 1)=   -0.12753203712807D-07
s(41, 1)=   -0.67949821258339D-07
s(42, 1)=   -0.16396289398253D-07
s(43, 1)=   -0.43771996408463D-07
s(44, 1)=   -0.12376736507648D-08
s(45, 1)=    0.53832073271757D-07
s(46, 1)=    0.15683302078249D-07
s(47, 1)=   -0.95298601138481D-07
s(48, 1)=   -0.59657889079683D-07
s(49, 1)=    0.77816630739978D-07
s(50, 1)=    0.39362579443894D-07
s(51, 1)=   -0.50147019425721D-07
s(52, 1)=   -0.76175189595045D-08
s(53, 1)=    0.44470821222391D-07
s(54, 1)=    0.58422934124297D-07
s(55, 1)=   -0.27531333298184D-07
s(56, 1)=    0.47425144454705D-07
s(57, 1)=    0.12692149563230D-07
s(58, 1)=   -0.62104080669803D-07
s(59, 1)=   -0.78129455156601D-09
s(60, 1)=    0.19250839858214D-07
s(61, 1)=    0.29756964767289D-07
s(62, 1)=   -0.36635669397708D-07
s(63, 1)=    0.39026727963702D-07
s(64, 1)=   -0.11666343243320D-08
s(65, 1)=   -0.39507755076019D-08
s(66, 1)=    0.13315627623450D-07
s(67, 1)=   -0.90964859573059D-08
s(68, 1)=    0.31922474126431D-07
s(69, 1)=   -0.89309440163424D-08
s(70, 1)=   -0.24655505132494D-08
s(71, 1)=   -0.42602408922089D-08
s(72, 1)=   -0.12030889796115D-07
s(73, 1)=   -0.10741257587395D-07
s(74, 1)=   -0.55499974865562D-08
s(75, 1)=   -0.35174435080564D-08
s(76, 1)=    0.15678335042096D-07
s(77, 1)=   -0.28372251356007D-08
s(78, 1)=    0.20084489528831D-07
s(79, 1)=    0.22045740795526D-08
s(80, 1)=    0.15356372022802D-07
s( 2, 2)=    0.48905551075504D-04
s( 3, 2)=    0.83539192403728D-05
s( 4, 2)=   -0.89701777491190D-05
s( 5, 2)=   -0.11641453008893D-05
s( 6, 2)=    0.14667268956886D-05
s( 7, 2)=   -0.62842100807606D-06
s( 8, 2)=    0.50596700143418D-06
s( 9, 2)=    0.38142891566534D-06
s(10, 2)=   -0.11124225282454D-05
s(11, 2)=   -0.99100594960733D-06
s(12, 2)=    0.53128493253623D-06
s(13, 2)=    0.85028075079998D-06
s(14, 2)=   -0.52161616518213D-06
s(15, 2)=   -0.78050858035765D-06
s(16, 2)=   -0.87214358332511D-07
s(17, 2)=    0.14579826763514D-06
s(18, 2)=    0.43798807886695D-07
s(19, 2)=    0.87664275315519D-07
s(20, 2)=   -0.23739329968468D-06
s(21, 2)=   -0.55617022532929D-06
s(22, 2)=    0.74634801430613D-07
s(23, 2)=    0.47058080251370D-06
s(24, 2)=   -0.13833943024497D-06
s(25, 2)=   -0.36206673055862D-06
s(26, 2)=    0.13771086763274D-07
s(27, 2)=    0.26908130569571D-06
s(28, 2)=    0.25184065503879D-06
s(29, 2)=    0.29109346693879D-08
s(30, 2)=   -0.44165883077743D-07
s(31, 2)=    0.13908167846995D-07
s(32, 2)=    0.18305835915308D-06
s(33, 2)=    0.20487537947298D-06
s(34, 2)=    0.47320171392541D-07
s(35, 2)=   -0.58882619512619D-07
s(36, 2)=   -0.12311353910300D-06
s(37, 2)=   -0.21919769425805D-07
s(38, 2)=    0.11709291424090D-06
s(39, 2)=    0.92699311123633D-07
s(40, 2)=   -0.90368996680850D-07
s(41, 2)=   -0.11102320608435D-06
s(42, 2)=    0.85441307939682D-07
s(43, 2)=    0.15019286440103D-06
s(44, 2)=   -0.15758773652609D-07
s(45, 2)=   -0.56546577893551D-07
s(46, 2)=    0.46414838712203D-07
s(47, 2)=    0.69440803211037D-07
s(48, 2)=   -0.41543319320956D-07
s(49, 2)=   -0.25453153619611D-07
s(50, 2)=    0.61069806618382D-07
s(51, 2)=   -0.21431871272039D-07
s(52, 2)=   -0.27493738239856D-07
s(53, 2)=    0.40499074786701D-07
s(54, 2)=    0.16409989015642D-07
s(55, 2)=   -0.46733116241997D-07
s(56, 2)=   -0.49935750762261D-07
s(57, 2)=    0.38916218999179D-07
s(58, 2)=   -0.49768112095755D-08
s(59, 2)=   -0.28668365214044D-07
s(60, 2)=    0.10106627158974D-07
s(61, 2)=   -0.15637139662540D-07
s(62, 2)=   -0.20915792684197D-07
s(63, 2)=   -0.54227075337162D-07
s(64, 2)=    0.33432699268793D-07
s(65, 2)=    0.10034791702217D-07
s(66, 2)=    0.10133225197098D-07
s(67, 2)=   -0.22905770890077D-07
s(68, 2)=   -0.22301580854256D-07
s(69, 2)=   -0.92177350430440D-08
s(70, 2)=   -0.18780903598666D-08
s(71, 2)=   -0.75542782412997D-09
s(72, 2)=    0.49577307486526D-08
s(73, 2)=   -0.10826784880676D-07
s(74, 2)=   -0.23224199983071D-07
s(75, 2)=    0.19561774315674D-07
s(76, 2)=   -0.30216464519069D-07
s(77, 2)=    0.17749568424107D-07
s(78, 2)=   -0.33399757719077D-07
s(79, 2)=    0.18010662015649D-07
s(80, 2)=   -0.32871462541030D-07
s( 3, 3)=    0.25551471726737D-04
s( 4, 3)=   -0.19274827646727D-06
s( 5, 3)=    0.27140581800371D-06
s( 6, 3)=    0.33222302672742D-06
s( 7, 3)=   -0.39609564442357D-06
s( 8, 3)=   -0.13382050309136D-05
s( 9, 3)=   -0.10048343370863D-05
s(10, 3)=    0.44208842620839D-06
s(11, 3)=    0.76239768915178D-06
s(12, 3)=    0.31335143763597D-06
s(13, 3)=    0.18130002960504D-06
s(14, 3)=   -0.31104921120362D-06
s(15, 3)=   -0.29856958137367D-06
s(16, 3)=    0.42481073247226D-06
s(17, 3)=    0.19122120603783D-06
s(18, 3)=   -0.46117693571270D-06
s(19, 3)=   -0.20023496082877D-06
s(20, 3)=    0.20997193488831D-06
s(21, 3)=    0.24698772474013D-06
s(22, 3)=    0.14526401968646D-07
s(23, 3)=   -0.13094420881434D-06
s(24, 3)=   -0.28724262482356D-06
s(25, 3)=   -0.14705252644312D-06
s(26, 3)=    0.12595891565708D-06
s(27, 3)=    0.31207723482287D-06
s(28, 3)=   -0.27650809086232D-06
s(29, 3)=   -0.30671076129865D-06
s(30, 3)=    0.84155662713398D-08
s(31, 3)=    0.54069483045578D-07
s(32, 3)=    0.37773579574019D-07
s(33, 3)=   -0.11854623883255D-06
s(34, 3)=   -0.16136229886968D-06
s(35, 3)=   -0.59238139036980D-07
s(36, 3)=    0.52705208366376D-07
s(37, 3)=    0.61587862484921D-07
s(38, 3)=   -0.11594308999457D-07
s(39, 3)=   -0.21618306258924D-07
s(40, 3)=   -0.14792502688220D-07
s(41, 3)=    0.62343085836212D-07
s(42, 3)=    0.25134346110112D-07
s(43, 3)=    0.24424514891413D-07
s(44, 3)=   -0.26996983351598D-07
s(45, 3)=   -0.72534663934446D-08
s(46, 3)=    0.29618360439917D-07
s(47, 3)=    0.52654231743027D-08
s(48, 3)=   -0.10583447514475D-07
s(49, 3)=    0.37626191938029D-07
s(50, 3)=    0.25596290670976D-07
s(51, 3)=   -0.23932391975686D-07
s(52, 3)=   -0.49857747681036D-07
s(53, 3)=    0.95895058624665D-08
s(54, 3)=    0.71293575106710D-07
s(55, 3)=    0.14572395176144D-07
s(56, 3)=   -0.46910828885870D-07
s(57, 3)=    0.35071316888314D-07
s(58, 3)=    0.31409118702607D-07
s(59, 3)=   -0.14417701951018D-07
s(60, 3)=    0.13876650780031D-07
s(61, 3)=   -0.35334741851824D-08
s(62, 3)=   -0.30324625813371D-07
s(63, 3)=   -0.32787634088217D-07
s(64, 3)=    0.12694055876807D-07
s(65, 3)=    0.22406210975383D-07
s(66, 3)=   -0.18141161946838D-07
s(67, 3)=    0.54521038956287D-09
s(68, 3)=   -0.38138689619787D-08
s(69, 3)=   -0.16093220344950D-07
s(70, 3)=   -0.37036403872396D-07
s(71, 3)=   -0.47186949728098D-08
s(72, 3)=   -0.27956746206030D-07
s(73, 3)=    0.85634852321390D-09
s(74, 3)=   -0.38325087199963D-07
s(75, 3)=    0.15398137020464D-07
s(76, 3)=   -0.39023342002934D-07
s(77, 3)=    0.19635939217429D-07
s(78, 3)=   -0.35504514917566D-07
s(79, 3)=    0.21761683080775D-07
s(80, 3)=   -0.31625978448866D-07
s( 4, 4)=   -0.12857823170070D-04
s( 5, 4)=   -0.33773263335749D-05
s( 6, 4)=    0.26339507884398D-05
s( 7, 4)=   -0.42144144853044D-06
s( 8, 4)=    0.14777413181852D-06
s( 9, 4)=    0.16194710093435D-05
s(10, 4)=   -0.69313507253819D-07
s(11, 4)=   -0.64233300320774D-06
s(12, 4)=    0.12301690995295D-06
s(13, 4)=    0.68736009234055D-06
s(14, 4)=   -0.95770079832859D-06
s(15, 4)=   -0.83528061339268D-06
s(16, 4)=   -0.54260378146251D-07
s(17, 4)=   -0.93929405238052D-08
s(18, 4)=    0.13191461217734D-06
s(19, 4)=    0.25226344496632D-06
s(20, 4)=    0.20621490254119D-06
s(21, 4)=   -0.14043749966466D-06
s(22, 4)=   -0.37858016062249D-06
s(23, 4)=    0.73495913921712D-07
s(24, 4)=    0.11779230311764D-06
s(25, 4)=    0.15527691778753D-06
s(26, 4)=    0.21058477382704D-06
s(27, 4)=    0.89681836558224D-07
s(28, 4)=    0.19622213253189D-07
s(29, 4)=    0.47678131095151D-07
s(30, 4)=    0.21380065806948D-06
s(31, 4)=    0.75620818268705D-08
s(32, 4)=    0.69970608136749D-07
s(33, 4)=    0.31985403196957D-07
s(34, 4)=    0.98787731189489D-07
s(35, 4)=    0.79432952777437D-07
s(36, 4)=   -0.22350132896334D-07
s(37, 4)=    0.66789039516572D-08
s(38, 4)=    0.64998554245920D-07
s(39, 4)=    0.36097324429760D-07
s(40, 4)=    0.53798620841362D-07
s(41, 4)=   -0.60633572489982D-08
s(42, 4)=    0.43908209925831D-07
s(43, 4)=    0.69304474175146D-07
s(44, 4)=    0.17136856044562D-07
s(45, 4)=   -0.21172262663520D-07
s(46, 4)=    0.15964270630663D-07
s(47, 4)=    0.68319629051684D-07
s(48, 4)=    0.22307463806329D-07
s(49, 4)=   -0.22098747481583D-07
s(50, 4)=   -0.13416105611403D-07
s(51, 4)=   -0.28727639863575D-08
s(52, 4)=    0.61362062327015D-08
s(53, 4)=    0.28029179020574D-07
s(54, 4)=    0.20493852005795D-07
s(55, 4)=   -0.34076081538132D-07
s(56, 4)=   -0.50027878148709D-07
s(57, 4)=    0.25275768023372D-07
s(58, 4)=    0.32984416637604D-07
s(59, 4)=    0.46748302094913D-08
s(60, 4)=   -0.28731083679469D-07
s(61, 4)=   -0.12185426332213D-08
s(62, 4)=    0.24196923372757D-07
s(63, 4)=   -0.38300851302597D-07
s(64, 4)=    0.23385333301277D-07
s(65, 4)=    0.10699297135166D-07
s(66, 4)=    0.24635886833335D-07
s(67, 4)=    0.16525365641931D-08
s(68, 4)=    0.11620271286385D-08
s(69, 4)=    0.21844011370929D-07
s(70, 4)=   -0.21954322860600D-07
s(71, 4)=   -0.29656501717783D-08
s(72, 4)=    0.11819742766063D-07
s(73, 4)=   -0.33898434165158D-07
s(74, 4)=    0.10552495877755D-07
s(75, 4)=   -0.60657709996744D-08
s(76, 4)=   -0.15629828556566D-07
s(77, 4)=   -0.20599003025979D-08
s(78, 4)=   -0.38749943891530D-08
s(79, 4)=   -0.12713374207514D-07
s(80, 4)=    0.22321713775651D-09
s( 5, 5)=    0.37742550130119D-05
s( 6, 5)=    0.16199010884387D-05
s( 7, 5)=   -0.13557075249867D-05
s( 8, 5)=   -0.16260186706174D-05
s( 9, 5)=   -0.15457469241661D-05
s(10, 5)=   -0.10600209030637D-05
s(11, 5)=    0.89636346823903D-06
s(12, 5)=    0.99814527374599D-06
s(13, 5)=   -0.40334461020521D-06
s(14, 5)=   -0.96150455430912D-06
s(15, 5)=   -0.47286535275867D-06
s(16, 5)=    0.58730941024010D-06
s(17, 5)=    0.31152083265288D-06
s(18, 5)=   -0.61825374203407D-06
s(19, 5)=   -0.44472448343507D-06
s(20, 5)=    0.14586438492309D-07
s(21, 5)=    0.16628168439163D-06
s(22, 5)=    0.31241402527141D-07
s(23, 5)=    0.18854458337510D-06
s(24, 5)=   -0.16489726638894D-06
s(25, 5)=   -0.31845421512104D-06
s(26, 5)=    0.22721950714351D-06
s(27, 5)=    0.28299545652663D-06
s(28, 5)=   -0.37389057002835D-07
s(29, 5)=   -0.17226603663226D-06
s(30, 5)=    0.14710390042110D-06
s(31, 5)=    0.17561089328596D-06
s(32, 5)=   -0.65422322261088D-07
s(33, 5)=    0.54998698313070D-07
s(34, 5)=   -0.58162952622315D-07
s(35, 5)=   -0.21460608956987D-06
s(36, 5)=    0.10197746868170D-07
s(37, 5)=    0.12520063761894D-06
s(38, 5)=   -0.35506011727971D-07
s(39, 5)=   -0.19201531312844D-06
s(40, 5)=   -0.12269200040578D-06
s(41, 5)=    0.89778390439825D-07
s(42, 5)=   -0.41913786741451D-07
s(43, 5)=   -0.91520006488999D-07
s(44, 5)=   -0.44920044599259D-07
s(45, 5)=    0.11643541732227D-07
s(46, 5)=   -0.26360725053581D-08
s(47, 5)=   -0.44760884506842D-08
s(48, 5)=   -0.21592819756822D-07
s(49, 5)=   -0.60937408379074D-08
s(50, 5)=    0.27401279526375D-08
s(51, 5)=    0.15851055218607D-07
s(52, 5)=   -0.39887968690482D-07
s(53, 5)=   -0.18343774352404D-07
s(54, 5)=    0.62372930417058D-08
s(55, 5)=    0.20924085053186D-07
s(56, 5)=   -0.30622186799619D-07
s(57, 5)=   -0.11748953359707D-07
s(58, 5)=    0.23880590247310D-07
s(59, 5)=    0.20034896605191D-07
s(60, 5)=   -0.16471279024280D-07
s(61, 5)=   -0.11357411549884D-07
s(62, 5)=   -0.18049010749233D-07
s(63, 5)=   -0.65094139439041D-07
s(64, 5)=   -0.20737977692137D-07
s(65, 5)=   -0.29791542438761D-07
s(66, 5)=    0.78550958031711D-09
s(67, 5)=    0.98305589918783D-08
s(68, 5)=    0.15360752910204D-07
s(69, 5)=    0.17205978100653D-07
s(70, 5)=    0.82258548517498D-08
s(71, 5)=    0.85106428921503D-08
s(72, 5)=    0.11577667889131D-07
s(73, 5)=    0.16308111802428D-08
s(74, 5)=    0.14229372953148D-07
s(75, 5)=   -0.82577367551857D-08
s(76, 5)=    0.16419910665001D-07
s(77, 5)=   -0.76141933727246D-08
s(78, 5)=    0.15663078015234D-07
s(79, 5)=   -0.38076118797541D-08
s(80, 5)=    0.19338680156276D-07
s( 6, 6)=    0.82048914571916D-06
s( 7, 6)=   -0.18976252698758D-05
s( 8, 6)=   -0.17859907392665D-05
s( 9, 6)=    0.57635343927424D-06
s(10, 6)=    0.11141482531981D-05
s(11, 6)=    0.20520417204917D-07
s(12, 6)=   -0.16241518863413D-05
s(13, 6)=   -0.92898699384914D-06
s(14, 6)=   -0.83030911917673D-07
s(15, 6)=    0.63029479290680D-06
s(16, 6)=    0.22377218327274D-06
s(17, 6)=   -0.30803335024641D-06
s(18, 6)=   -0.67899952282254D-06
s(19, 6)=   -0.15250901091529D-06
s(20, 6)=    0.37831534677241D-06
s(21, 6)=    0.27709441869198D-06
s(22, 6)=   -0.52996707240534D-06
s(23, 6)=   -0.49676668701636D-06
s(24, 6)=    0.36332436606584D-07
s(25, 6)=    0.35015834406465D-06
s(26, 6)=    0.58170241349480D-07
s(27, 6)=   -0.32085152929534D-06
s(28, 6)=   -0.26174056751663D-06
s(29, 6)=    0.74912443471519D-08
s(30, 6)=    0.36062403088968D-06
s(31, 6)=    0.11445836997106D-06
s(32, 6)=   -0.21026389477025D-06
s(33, 6)=   -0.26176185843027D-06
s(34, 6)=    0.27649159609722D-07
s(35, 6)=    0.18659652364132D-06
s(36, 6)=    0.64867610890407D-08
s(37, 6)=   -0.12178155734875D-06
s(38, 6)=   -0.11711628132855D-06
s(39, 6)=   -0.54345146341479D-08
s(40, 6)=    0.71059805395230D-07
s(41, 6)=    0.46535261677561D-07
s(42, 6)=   -0.19077964613060D-07
s(43, 6)=   -0.61725941246726D-07
s(44, 6)=    0.56827990615115D-07
s(45, 6)=    0.49360626933475D-07
s(46, 6)=   -0.44437585732807D-07
s(47, 6)=   -0.16590532201585D-07
s(48, 6)=    0.52118679119891D-07
s(49, 6)=   -0.39270934787392D-07
s(50, 6)=   -0.29200578663369D-07
s(51, 6)=    0.71930752953130D-08
s(52, 6)=   -0.55013112384309D-08
s(53, 6)=   -0.41179828878416D-07
s(54, 6)=    0.27775828195230D-07
s(55, 6)=    0.40970764709541D-07
s(56, 6)=   -0.93882334391395D-08
s(57, 6)=   -0.19406697176474D-07
s(58, 6)=    0.33469412443438D-07
s(59, 6)=    0.12935720779420D-07
s(60, 6)=   -0.11381909184848D-07
s(61, 6)=    0.94081584158453D-08
s(62, 6)=    0.47789557270427D-07
s(63, 6)=    0.17622793297102D-07
s(64, 6)=   -0.18345308909856D-07
s(65, 6)=    0.26940425956765D-08
s(66, 6)=   -0.12905439356089D-07
s(67, 6)=   -0.19649164082368D-07
s(68, 6)=    0.15817487307388D-07
s(69, 6)=    0.98315960599906D-08
s(70, 6)=   -0.36651513112467D-08
s(71, 6)=   -0.31455261139310D-08
s(72, 6)=    0.23442550631311D-07
s(73, 6)=   -0.91447673426722D-08
s(74, 6)=    0.12777811795930D-07
s(75, 6)=   -0.38604892131951D-08
s(76, 6)=    0.15707324950716D-07
s(77, 6)=   -0.54637645120108D-08
s(78, 6)=    0.12652007604582D-07
s(79, 6)=    0.76513109485040D-08
s(80, 6)=    0.80948836322871D-08
s( 7, 7)=   -0.17718148905347D-05
s( 8, 7)=    0.16408397974219D-05
s( 9, 7)=    0.86647021106569D-06
s(10, 7)=   -0.61835547459555D-06
s(11, 7)=   -0.86232464665112D-06
s(12, 7)=   -0.11248859759827D-06
s(13, 7)=    0.39634821776235D-06
s(14, 7)=    0.24045824377199D-06
s(15, 7)=    0.36039156989872D-06
s(16, 7)=    0.31060939091896D-06
s(17, 7)=   -0.24009029696849D-06
s(18, 7)=   -0.16932094993760D-06
s(19, 7)=    0.18821371044436D-06
s(20, 7)=    0.40126716006775D-07
s(21, 7)=   -0.30263479586734D-06
s(22, 7)=    0.12697718689314D-07
s(23, 7)=    0.28287048180980D-06
s(24, 7)=    0.50713252130312D-06
s(25, 7)=    0.26061153429560D-06
s(26, 7)=   -0.23148709970149D-06
s(27, 7)=   -0.11121282776403D-06
s(28, 7)=    0.19604577315001D-06
s(29, 7)=    0.16686812676659D-06
s(30, 7)=    0.30214879537025D-07
s(31, 7)=   -0.13377877566804D-06
s(32, 7)=   -0.16610438475314D-06
s(33, 7)=    0.28746228761159D-06
s(34, 7)=    0.22764174789091D-06
s(35, 7)=   -0.10676698900774D-06
s(36, 7)=   -0.11823947886321D-06
s(37, 7)=   -0.18761169597968D-07
s(38, 7)=    0.15498239194539D-06
s(39, 7)=   -0.39950543926411D-08
s(40, 7)=   -0.43229921049524D-07
s(41, 7)=    0.11412007712194D-06
s(42, 7)=   -0.14466968510796D-07
s(43, 7)=   -0.50395367209663D-07
s(44, 7)=    0.10038828070189D-07
s(45, 7)=    0.89846934942603D-07
s(46, 7)=    0.26646513249887D-07
s(47, 7)=   -0.10731318307303D-06
s(48, 7)=    0.42121798444175D-08
s(49, 7)=    0.75577709045872D-07
s(50, 7)=   -0.42621912089417D-07
s(51, 7)=   -0.16964770982555D-07
s(52, 7)=    0.22999036083795D-07
s(53, 7)=   -0.20887288050776D-07
s(54, 7)=   -0.44108867497604D-07
s(55, 7)=    0.13930811021967D-08
s(56, 7)=   -0.83764704949305D-08
s(57, 7)=   -0.69153676523682D-08
s(58, 7)=    0.18007382119490D-08
s(59, 7)=    0.33336892845039D-08
s(60, 7)=   -0.41896555386511D-08
s(61, 7)=   -0.17039520923528D-07
s(62, 7)=   -0.16712720904900D-07
s(63, 7)=    0.37298748429177D-07
s(64, 7)=   -0.25654109191059D-07
s(65, 7)=   -0.13749720601388D-07
s(66, 7)=    0.26676806808541D-07
s(67, 7)=    0.16476181922272D-07
s(68, 7)=    0.83007982073902D-08
s(69, 7)=    0.26847074426587D-07
s(70, 7)=    0.18339215465452D-07
s(71, 7)=   -0.12096493390590D-07
s(72, 7)=    0.30873049918127D-07
s(73, 7)=   -0.22948720255903D-08
s(74, 7)=    0.42303559984739D-09
s(75, 7)=   -0.31468500927299D-09
s(76, 7)=   -0.78036646136261D-08
s(77, 7)=   -0.15461825889251D-07
s(78, 7)=    0.35576129993784D-08
s(79, 7)=   -0.87002735760482D-08
s(80, 7)=   -0.34389487185570D-09
s( 8, 8)=   -0.24959849099076D-06
s( 9, 8)=   -0.14580693512087D-06
s(10, 8)=    0.81805036262635D-06
s(11, 8)=    0.76433996425429D-06
s(12, 8)=   -0.38702672712394D-06
s(13, 8)=   -0.64465263597708D-08
s(14, 8)=    0.37300526843875D-06
s(15, 8)=    0.47222387124156D-06
s(16, 8)=    0.38955608515904D-06
s(17, 8)=   -0.20383440752989D-06
s(18, 8)=   -0.26151194759556D-06
s(19, 8)=   -0.13178247598978D-06
s(20, 8)=    0.31649702590755D-06
s(21, 8)=    0.31101517689284D-06
s(22, 8)=   -0.21168255118607D-06
s(23, 8)=   -0.28466023855760D-06
s(24, 8)=    0.11719052550134D-06
s(25, 8)=    0.15843866010930D-06
s(26, 8)=    0.23810692528317D-07
s(27, 8)=    0.53647006877697D-07
s(28, 8)=   -0.31906307271979D-08
s(29, 8)=   -0.10391985530998D-07
s(30, 8)=    0.95109571880006D-07
s(31, 8)=    0.57365947063787D-07
s(32, 8)=   -0.40518346014117D-07
s(33, 8)=   -0.15499078054985D-06
s(34, 8)=    0.37313770945502D-07
s(35, 8)=    0.93476114898862D-07
s(36, 8)=   -0.45262218231854D-07
s(37, 8)=   -0.98638451420277D-07
s(38, 8)=    0.69215062222509D-07
s(39, 8)=    0.30032636960184D-07
s(40, 8)=   -0.75098783816941D-07
s(41, 8)=   -0.63812778667225D-07
s(42, 8)=   -0.57734182666462D-08
s(43, 8)=   -0.43605753090059D-07
s(44, 8)=    0.53239834915912D-07
s(45, 8)=    0.29151322995370D-07
s(46, 8)=   -0.19230450020761D-07
s(47, 8)=   -0.39920862610279D-07
s(48, 8)=    0.27925566188812D-07
s(49, 8)=    0.58064949635228D-07
s(50, 8)=   -0.71988778892891D-08
s(51, 8)=    0.78900516452218D-08
s(52, 8)=    0.28429706686802D-07
s(53, 8)=    0.76661527952012D-08
s(54, 8)=    0.72627505164607D-08
s(55, 8)=   -0.67541296281969D-08
s(56, 8)=    0.17380975529252D-07
s(57, 8)=   -0.64317567606288D-08
s(58, 8)=    0.11915145445347D-07
s(59, 8)=    0.12196090940292D-07
s(60, 8)=   -0.27987670013253D-07
s(61, 8)=   -0.71249626874001D-08
s(62, 8)=   -0.51867829783309D-08
s(63, 8)=   -0.19468386912618D-07
s(64, 8)=    0.23098695913554D-07
s(65, 8)=    0.28411589627994D-07
s(66, 8)=    0.27127785981138D-07
s(67, 8)=    0.20265388073624D-07
s(68, 8)=   -0.69432862125621D-08
s(69, 8)=    0.51173112821441D-08
s(70, 8)=   -0.43588065245559D-08
s(71, 8)=    0.61109017794536D-08
s(72, 8)=    0.21118182460781D-07
s(73, 8)=   -0.94860333532873D-09
s(74, 8)=    0.43499469247456D-08
s(75, 8)=    0.15440798733987D-07
s(76, 8)=    0.42930395860814D-08
s(77, 8)=   -0.10648979146836D-07
s(78, 8)=    0.23716172947373D-07
s(79, 8)=   -0.23048343056173D-07
s(80, 8)=    0.16544765807105D-07
s( 9, 9)=   -0.65627326888320D-06
s(10, 9)=   -0.14595950764968D-05
s(11, 9)=   -0.41836765661321D-06
s(12, 9)=    0.47512601540373D-06
s(13, 9)=    0.80247505155953D-06
s(14, 9)=    0.11498533744602D-05
s(15, 9)=   -0.37235910487197D-06
s(16, 9)=   -0.93489733864836D-06
s(17, 9)=   -0.32805356533853D-06
s(18, 9)=    0.22858126868659D-06
s(19, 9)=    0.39400348404035D-06
s(20, 9)=    0.18072704192137D-06
s(21, 9)=   -0.10183574008031D-06
s(22, 9)=   -0.30654676245482D-07
s(23, 9)=    0.10420719046342D-06
s(24, 9)=    0.33274868234587D-06
s(25, 9)=    0.21290437080171D-06
s(26, 9)=   -0.51156668516244D-06
s(27, 9)=   -0.38743075186968D-06
s(28, 9)=    0.42432381937433D-06
s(29, 9)=    0.33265775129990D-06
s(30, 9)=   -0.14915562556088D-07
s(31, 9)=   -0.41861144470828D-07
s(32, 9)=   -0.99150362320510D-07
s(33, 9)=    0.10711465229193D-06
s(34, 9)=    0.12194543304261D-06
s(35, 9)=   -0.46902631272439D-07
s(36, 9)=   -0.66049672182928D-07
s(37, 9)=   -0.52095734962230D-07
s(38, 9)=    0.13573875181403D-06
s(39, 9)=    0.18417451602229D-07
s(40, 9)=   -0.11934084112825D-06
s(41, 9)=   -0.18995682305052D-07
s(42, 9)=    0.54219013565157D-07
s(43, 9)=    0.35057697093395D-07
s(44, 9)=    0.53305704921528D-09
s(45, 9)=    0.64794407176129D-07
s(46, 9)=    0.27331745831916D-07
s(47, 9)=   -0.58988474395239D-07
s(48, 9)=    0.17461117306505D-07
s(49, 9)=    0.87270646852327D-07
s(50, 9)=   -0.71342940132626D-08
s(51, 9)=   -0.28435611154970D-07
s(52, 9)=    0.11805683008097D-07
s(53, 9)=    0.89143404074597D-08
s(54, 9)=   -0.26555392074377D-07
s(55, 9)=    0.33828855511918D-07
s(56, 9)=    0.42546390083977D-07
s(57, 9)=   -0.36734873423041D-07
s(58, 9)=   -0.22395650708275D-07
s(59, 9)=    0.24315885012474D-07
s(60, 9)=    0.79944980951610D-08
s(61, 9)=   -0.38296575561497D-07
s(62, 9)=    0.70571132839722D-08
s(63, 9)=    0.25038246911884D-07
s(64, 9)=    0.34757865716423D-07
s(65, 9)=    0.71465665945082D-08
s(66, 9)=    0.96220861771692D-08
s(67, 9)=    0.14907466291611D-07
s(68, 9)=   -0.12459707429573D-07
s(69, 9)=    0.25196834012392D-09
s(70, 9)=   -0.24635603469987D-07
s(71, 9)=   -0.34947072691407D-08
s(72, 9)=   -0.40986826279927D-07
s(73, 9)=   -0.53285000800155D-08
s(74, 9)=   -0.21864929212364D-07
s(75, 9)=    0.31655289123553D-08
s(76, 9)=   -0.71972613348107D-08
s(77, 9)=   -0.11715970978303D-07
s(78, 9)=   -0.91581498788497D-08
s(79, 9)=   -0.33289555113209D-08
s(80, 9)=   -0.69618551966212D-08
s(10,10)=    0.75110049124197D-06
s(11,10)=    0.19602479146368D-05
s(12,10)=    0.13735593875557D-05
s(13,10)=   -0.64528317351519D-06
s(14,10)=   -0.15806738741427D-05
s(15,10)=   -0.15651979227771D-06
s(16,10)=    0.34092895016214D-06
s(17,10)=    0.26625026138228D-06
s(18,10)=   -0.66632869044241D-07
s(19,10)=   -0.78558502011907D-07
s(20,10)=   -0.27558540599626D-07
s(21,10)=    0.10612069522761D-06
s(22,10)=    0.27759528560520D-06
s(23,10)=   -0.14762936279683D-06
s(24,10)=   -0.30649417359411D-06
s(25,10)=   -0.24934120573227D-06
s(26,10)=    0.45214442081102D-07
s(27,10)=    0.24920251008875D-07
s(28,10)=    0.11146872495122D-07
s(29,10)=   -0.26060426643668D-07
s(30,10)=   -0.23781098460740D-06
s(31,10)=   -0.22249872682332D-06
s(32,10)=    0.10939022689204D-06
s(33,10)=   -0.99274702786069D-08
s(34,10)=   -0.11057200235541D-06
s(35,10)=   -0.71666688130456D-07
s(36,10)=    0.34024516192873D-07
s(37,10)=    0.55363294504361D-07
s(38,10)=    0.35642895088882D-08
s(39,10)=    0.10576871280035D-07
s(40,10)=   -0.11697411911326D-06
s(41,10)=   -0.73996793185219D-07
s(42,10)=    0.38678791991800D-07
s(43,10)=    0.24828811765688D-07
s(44,10)=   -0.42505224541213D-07
s(45,10)=   -0.87400153715382D-07
s(46,10)=   -0.52257551226251D-08
s(47,10)=   -0.15139873000105D-07
s(48,10)=   -0.30765583056167D-07
s(49,10)=    0.37035367656631D-07
s(50,10)=   -0.25016564928436D-07
s(51,10)=   -0.35771827517194D-07
s(52,10)=   -0.33478355753907D-09
s(53,10)=    0.20219187490712D-07
s(54,10)=    0.38102648478186D-08
s(55,10)=   -0.34616352014803D-07
s(56,10)=   -0.51970170086516D-08
s(57,10)=    0.46109370443295D-07
s(58,10)=   -0.62288697447250D-08
s(59,10)=    0.16408443629483D-07
s(60,10)=    0.17080610897797D-07
s(61,10)=   -0.20496638449623D-07
s(62,10)=    0.88434962009344D-08
s(63,10)=   -0.45568632645146D-07
s(64,10)=   -0.32009339658894D-07
s(65,10)=    0.22650298617807D-08
s(66,10)=    0.11883874117333D-07
s(67,10)=   -0.29624723139076D-07
s(68,10)=    0.91697466414316D-08
s(69,10)=   -0.40427796463915D-08
s(70,10)=   -0.78382603351011D-08
s(71,10)=    0.62882121043255D-08
s(72,10)=    0.10388060801758D-07
s(73,10)=   -0.10887005067750D-07
s(74,10)=    0.93326955954817D-08
s(75,10)=   -0.35375814440110D-08
s(76,10)=    0.11440477826086D-07
s(77,10)=    0.11280465372383D-07
s(78,10)=    0.80737354952143D-09
s(79,10)=    0.13169751291096D-07
s(80,10)=   -0.89384411088562D-08
s(11,11)=   -0.32243874465540D-06
s(12,11)=   -0.16488455183527D-05
s(13,11)=   -0.82051095635799D-06
s(14,11)=    0.15963288344652D-06
s(15,11)=    0.57043780534665D-06
s(16,11)=   -0.31162197414919D-06
s(17,11)=    0.51284823319059D-07
s(18,11)=   -0.30113265411810D-06
s(19,11)=   -0.39066464008848D-06
s(20,11)=    0.19716210205580D-07
s(21,11)=   -0.17102747377246D-06
s(22,11)=   -0.16030018161286D-06
s(23,11)=   -0.31284845289221D-06
s(24,11)=    0.30924886178203D-06
s(25,11)=    0.40328090713224D-06
s(26,11)=   -0.22008090349919D-06
s(27,11)=   -0.22263962027705D-06
s(28,11)=    0.28201322279299D-06
s(29,11)=    0.58131493898922D-07
s(30,11)=   -0.11375026639270D-06
s(31,11)=    0.18978094742895D-06
s(32,11)=    0.81541810295914D-07
s(33,11)=   -0.15861779452979D-06
s(34,11)=   -0.21333796420055D-06
s(35,11)=    0.33398335918711D-07
s(36,11)=    0.37795393184265D-07
s(37,11)=   -0.47045914061136D-07
s(38,11)=   -0.96238569513032D-07
s(39,11)=   -0.26332968760336D-07
s(40,11)=   -0.57618669052149D-07
s(41,11)=   -0.12609690722743D-07
s(42,11)=    0.92647921190134D-07
s(43,11)=    0.15554619277272D-07
s(44,11)=   -0.14729797978674D-06
s(45,11)=   -0.35970928874431D-07
s(46,11)=    0.63535970160862D-07
s(47,11)=    0.36411397717139D-07
s(48,11)=   -0.96957422611593D-07
s(49,11)=   -0.30295526388977D-08
s(50,11)=    0.71644391986672D-07
s(51,11)=    0.24307941537757D-08
s(52,11)=   -0.40647107856712D-07
s(53,11)=   -0.21995970977921D-07
s(54,11)=    0.10271836091161D-07
s(55,11)=    0.19536869823667D-07
s(56,11)=    0.28264297233668D-07
s(57,11)=    0.59982895719428D-08
s(58,11)=    0.52151716282265D-08
s(59,11)=   -0.21249449590778D-07
s(60,11)=    0.16590631640645D-07
s(61,11)=    0.30957968332036D-07
s(62,11)=    0.12664703332618D-07
s(63,11)=    0.58286283310618D-08
s(64,11)=   -0.15894064666470D-07
s(65,11)=   -0.10004824485461D-07
s(66,11)=    0.18897384981152D-08
s(67,11)=    0.92202405384102D-08
s(68,11)=   -0.41312478365989D-07
s(69,11)=   -0.47410459394041D-08
s(70,11)=   -0.13962678469304D-07
s(71,11)=   -0.29120421461712D-07
s(72,11)=   -0.33846955715306D-08
s(73,11)=   -0.21463688511713D-07
s(74,11)=   -0.15861238274070D-07
s(75,11)=    0.40683442664938D-08
s(76,11)=   -0.18831197923717D-07
s(77,11)=    0.12897496781262D-07
s(78,11)=   -0.89125352021773D-08
s(79,11)=   -0.13556648439713D-08
s(80,11)=   -0.44640224309969D-08
s(12,12)=   -0.88920779082167D-07
s(13,12)=   -0.42925246963553D-06
s(14,12)=   -0.30381120209259D-06
s(15,12)=    0.75705454254370D-06
s(16,12)=    0.66822060944807D-06
s(17,12)=   -0.31080678295820D-08
s(18,12)=   -0.10087201702018D-06
s(19,12)=   -0.21529973577778D-07
s(20,12)=   -0.17140589471271D-06
s(21,12)=    0.77850994715886D-07
s(22,12)=   -0.69979530788204D-07
s(23,12)=   -0.93144143314537D-08
s(24,12)=   -0.31991542759975D-06
s(25,12)=   -0.23780819412334D-06
s(26,12)=   -0.57158660489078D-07
s(27,12)=   -0.62059213801999D-08
s(28,12)=   -0.15907547414198D-07
s(29,12)=    0.10171806953978D-06
s(30,12)=    0.10491439130093D-06
s(31,12)=   -0.45676262519329D-07
s(32,12)=    0.25768458653127D-07
s(33,12)=   -0.35446527289965D-07
s(34,12)=   -0.47116393079744D-07
s(35,12)=    0.45383703782855D-07
s(36,12)=    0.58386526496614D-08
s(37,12)=   -0.85162673042796D-07
s(38,12)=   -0.91882765548383D-07
s(39,12)=   -0.36025149381032D-07
s(40,12)=    0.14595185410301D-07
s(41,12)=    0.17593887978747D-07
s(42,12)=    0.92166455685116D-08
s(43,12)=    0.14563806235588D-07
s(44,12)=    0.63545789496187D-08
s(45,12)=   -0.17181064007162D-07
s(46,12)=    0.67254550365132D-07
s(47,12)=   -0.88853639772640D-08
s(48,12)=   -0.39625401175776D-07
s(49,12)=    0.33547808860120D-07
s(50,12)=   -0.27539465725800D-07
s(51,12)=   -0.72009659384084D-07
s(52,12)=   -0.65160563283952D-08
s(53,12)=    0.58297585674734D-07
s(54,12)=    0.16743434242418D-08
s(55,12)=   -0.31503739387471D-07
s(56,12)=    0.22844849793649D-07
s(57,12)=    0.12620994119193D-07
s(58,12)=   -0.28295430174200D-07
s(59,12)=   -0.23511171040329D-07
s(60,12)=    0.22816487987095D-07
s(61,12)=    0.17900096358017D-07
s(62,12)=   -0.19359220793239D-07
s(63,12)=   -0.20436015960463D-07
s(64,12)=    0.29421441916531D-07
s(65,12)=    0.25063647453935D-09
s(66,12)=   -0.23434053627949D-07
s(67,12)=    0.47996977707156D-08
s(68,12)=   -0.26536851961145D-07
s(69,12)=   -0.24253678258790D-07
s(70,12)=   -0.29390882380831D-08
s(71,12)=   -0.16333925954900D-08
s(72,12)=   -0.28490175642328D-07
s(73,12)=    0.60768142527698D-08
s(74,12)=   -0.72357878393461D-08
s(75,12)=    0.47286146894435D-08
s(76,12)=   -0.29640799063399D-07
s(77,12)=    0.77994725347343D-08
s(78,12)=   -0.25082417039953D-07
s(79,12)=    0.44238135506005D-08
s(80,12)=   -0.29369750400863D-07
s(13,13)=    0.87090151746402D-06
s(14,13)=    0.20603291432190D-05
s(15,13)=    0.70781194137062D-06
s(16,13)=   -0.59505211895115D-06
s(17,13)=   -0.79988519244723D-06
s(18,13)=   -0.26781919637356D-06
s(19,13)=   -0.39638127126335D-07
s(20,13)=    0.34308208665972D-06
s(21,13)=    0.10176613624133D-06
s(22,13)=   -0.17303916297882D-06
s(23,13)=   -0.97853601308965D-07
s(24,13)=   -0.23487959545685D-07
s(25,13)=    0.27930133682935D-06
s(26,13)=   -0.81776277087869D-08
s(27,13)=   -0.11040107048514D-06
s(28,13)=    0.35518665380933D-07
s(29,13)=    0.28024876055130D-07
s(30,13)=    0.62729815192089D-07
s(31,13)=    0.14212575392451D-06
s(32,13)=    0.18847039668974D-07
s(33,13)=   -0.18720676208875D-07
s(34,13)=   -0.12328174691030D-06
s(35,13)=    0.13471841102089D-07
s(36,13)=    0.18504594338397D-06
s(37,13)=    0.41332784858941D-07
s(38,13)=   -0.90322029947740D-07
s(39,13)=   -0.89696583112129D-07
s(40,13)=    0.31621240917421D-07
s(41,13)=   -0.10732345514734D-08
s(42,13)=    0.51671318931160D-08
s(43,13)=    0.18086307630789D-08
s(44,13)=   -0.49976630665191D-07
s(45,13)=   -0.28555583722465D-07
s(46,13)=    0.72434265587148D-07
s(47,13)=    0.11841653994144D-06
s(48,13)=   -0.63007329286238D-07
s(49,13)=   -0.81663659434407D-07
s(50,13)=    0.92299309105374D-07
s(51,13)=    0.61277108914371D-07
s(52,13)=   -0.31928594863025D-07
s(53,13)=   -0.42066672979787D-07
s(54,13)=    0.20881374939344D-07
s(55,13)=    0.23357346197626D-07
s(56,13)=    0.55376552249704D-08
s(57,13)=    0.26637506077608D-07
s(58,13)=    0.45853451895835D-08
s(59,13)=   -0.37810514693800D-07
s(60,13)=    0.17827259495442D-07
s(61,13)=    0.33606099974048D-07
s(62,13)=    0.11087419305904D-07
s(63,13)=   -0.18858176682939D-07
s(64,13)=   -0.12846138794470D-07
s(65,13)=    0.62418715082964D-08
s(66,13)=    0.17820386364416D-07
s(67,13)=   -0.86672830025847D-08
s(68,13)=    0.42984837759796D-07
s(69,13)=    0.22895749266152D-07
s(70,13)=    0.24048612493333D-07
s(71,13)=    0.35644430659660D-07
s(72,13)=    0.32707523090933D-07
s(73,13)=    0.22954078635023D-07
s(74,13)=    0.23616745111190D-07
s(75,13)=    0.42243201693742D-07
s(76,13)=   -0.41136165002182D-08
s(77,13)=    0.48575954733840D-07
s(78,13)=   -0.10767185336948D-07
s(79,13)=    0.66336704576549D-07
s(80,13)=   -0.44014299099342D-07
s(14,14)=   -0.74240282317235D-06
s(15,14)=   -0.13571091998453D-05
s(16,14)=   -0.93817219005892D-06
s(17,14)=   -0.46937395151752D-06
s(18,14)=    0.19153934057115D-06
s(19,14)=    0.51536173249611D-06
s(20,14)=    0.23375072844855D-06
s(21,14)=   -0.39050845702147D-06
s(22,14)=   -0.51863544921841D-06
s(23,14)=   -0.32332506607369D-06
s(24,14)=    0.39865595042181D-06
s(25,14)=    0.26591743031456D-06
s(26,14)=   -0.81802188971568D-07
s(27,14)=   -0.32291096639111D-07
s(28,14)=    0.12860175423286D-06
s(29,14)=    0.13830153071439D-06
s(30,14)=    0.35042509660984D-07
s(31,14)=   -0.11918938957786D-06
s(32,14)=   -0.18345702668221D-06
s(33,14)=   -0.13524513613963D-07
s(34,14)=    0.77822106322462D-07
s(35,14)=    0.98718013634759D-07
s(36,14)=   -0.52654027292184D-07
s(37,14)=   -0.21310515155309D-06
s(38,14)=   -0.28779444492386D-07
s(39,14)=    0.51388765792751D-07
s(40,14)=    0.51863412740192D-07
s(41,14)=    0.53744070843108D-07
s(42,14)=   -0.91839891177993D-07
s(43,14)=   -0.25400916794614D-07
s(44,14)=    0.54614253806952D-07
s(45,14)=    0.66912389017243D-07
s(46,14)=    0.44942897130642D-07
s(47,14)=   -0.46320048520030D-07
s(48,14)=    0.35565097794370D-07
s(49,14)=    0.12012861482544D-06
s(50,14)=    0.77541596409904D-07
s(51,14)=   -0.35353777578629D-07
s(52,14)=   -0.33571218394550D-07
s(53,14)=    0.30101593756227D-07
s(54,14)=    0.26974668151704D-07
s(55,14)=   -0.15995052196413D-07
s(56,14)=   -0.52692369031442D-08
s(57,14)=   -0.18302276242218D-07
s(58,14)=   -0.21796489423620D-07
s(59,14)=    0.21371899255291D-08
s(60,14)=    0.48769586142551D-07
s(61,14)=   -0.13661773678492D-08
s(62,14)=   -0.79691872794417D-08
s(63,14)=    0.12767619131777D-09
s(64,14)=    0.23277013161202D-08
s(65,14)=   -0.15345742502274D-07
s(66,14)=   -0.10229721399366D-07
s(67,14)=    0.40370458826901D-07
s(68,14)=    0.75088309600281D-08
s(69,14)=    0.69388818240717D-08
s(70,14)=    0.18188183034137D-07
s(71,14)=    0.15492869250359D-07
s(72,14)=   -0.15765666455539D-07
s(73,14)=    0.63309934038770D-08
s(74,14)=   -0.36176635096533D-08
s(75,14)=    0.85598158046522D-08
s(76,14)=   -0.20951869182447D-07
s(77,14)=    0.31145212456406D-07
s(78,14)=   -0.24169629065450D-07
s(79,14)=    0.20848670628370D-07
s(80,14)=   -0.94467363979311D-08
s(15,15)=    0.12600108653866D-06
s(16,15)=    0.32271595745797D-06
s(17,15)=    0.40586600900225D-06
s(18,15)=    0.75111672234084D-06
s(19,15)=    0.50346290834619D-06
s(20,15)=   -0.28657327016529D-06
s(21,15)=   -0.57915871656014D-06
s(22,15)=    0.25969276409394D-06
s(23,15)=    0.70089053193120D-06
s(24,15)=   -0.29883118930912D-07
s(25,15)=   -0.14451769473653D-06
s(26,15)=   -0.48198140002715D-07
s(27,15)=    0.86264490499589D-07
s(28,15)=    0.12886841925463D-06
s(29,15)=    0.13684825500462D-06
s(30,15)=   -0.18797149406364D-07
s(31,15)=    0.15358171628512D-07
s(32,15)=    0.78329095841779D-07
s(33,15)=    0.99712503238671D-07
s(34,15)=    0.67720736038059D-07
s(35,15)=   -0.17093110750896D-06
s(36,15)=   -0.78041032755813D-07
s(37,15)=    0.43826692760960D-07
s(38,15)=   -0.79099615504350D-07
s(39,15)=    0.22676114875732D-07
s(40,15)=   -0.34353729658173D-07
s(41,15)=   -0.54211660036877D-07
s(42,15)=    0.16078659653430D-07
s(43,15)=    0.81221155568093D-07
s(44,15)=    0.26945602783787D-07
s(45,15)=   -0.11550748808496D-06
s(46,15)=    0.64854612383680D-08
s(47,15)=    0.59305343934151D-07
s(48,15)=    0.17334533420205D-07
s(49,15)=   -0.26013901667132D-07
s(50,15)=    0.19676556187785D-07
s(51,15)=    0.43600441789107D-07
s(52,15)=   -0.17613984985822D-07
s(53,15)=    0.13042832953335D-08
s(54,15)=    0.21604516348254D-07
s(55,15)=    0.83649756938879D-08
s(56,15)=   -0.25031837850356D-07
s(57,15)=   -0.10202997939252D-07
s(58,15)=    0.20575074324128D-07
s(59,15)=   -0.12508788110162D-07
s(60,15)=   -0.78480710934675D-08
s(61,15)=    0.26665962449193D-07
s(62,15)=   -0.44229316448830D-07
s(63,15)=   -0.49351437360333D-07
s(64,15)=    0.37661603535803D-07
s(65,15)=    0.36167346808350D-07
s(66,15)=    0.51363730875905D-08
s(67,15)=    0.36600538010786D-08
s(68,15)=    0.20002172436516D-07
s(69,15)=    0.28148380444181D-07
s(70,15)=    0.21821541647166D-07
s(71,15)=    0.17456029122273D-08
s(72,15)=    0.37208080162656D-07
s(73,15)=   -0.15048696636328D-07
s(74,15)=    0.19767509098502D-07
s(75,15)=   -0.47994233940629D-08
s(76,15)=    0.14017048848831D-07
s(77,15)=   -0.11215950504154D-07
s(78,15)=    0.14580230802822D-07
s(79,15)=   -0.20740278485624D-07
s(80,15)=    0.18176997659370D-07
s(16,16)=    0.11999639763360D-06
s(17,16)=    0.56520236837337D-06
s(18,16)=    0.22241299522815D-06
s(19,16)=   -0.13582747558851D-06
s(20,16)=   -0.27521114384639D-06
s(21,16)=   -0.13714024214610D-07
s(22,16)=   -0.42714935026343D-07
s(23,16)=    0.13846890751522D-06
s(24,16)=    0.32611904178744D-06
s(25,16)=    0.24278245987262D-06
s(26,16)=   -0.11672236873981D-06
s(27,16)=   -0.23027802135381D-06
s(28,16)=    0.64961121297676D-07
s(29,16)=    0.23471924216379D-07
s(30,16)=    0.63772434226171D-07
s(31,16)=   -0.73835040421860D-08
s(32,16)=   -0.61206732178322D-07
s(33,16)=   -0.55020494835377D-07
s(34,16)=   -0.42828380877650D-08
s(35,16)=    0.75980026207406D-07
s(36,16)=   -0.16217047689604D-07
s(37,16)=   -0.26913462260087D-08
s(38,16)=    0.18528209966077D-07
s(39,16)=   -0.19324089706529D-07
s(40,16)=   -0.82312570349735D-08
s(41,16)=   -0.54669670503534D-08
s(42,16)=    0.64355045434192D-08
s(43,16)=   -0.91015210205317D-07
s(44,16)=   -0.43879876041803D-08
s(45,16)=    0.46355226427945D-07
s(46,16)=    0.23826788645936D-07
s(47,16)=   -0.26727931830145D-07
s(48,16)=   -0.39825287398680D-07
s(49,16)=    0.27701521836883D-07
s(50,16)=    0.37425392830532D-07
s(51,16)=    0.32479270044008D-08
s(52,16)=   -0.31910867126780D-07
s(53,16)=   -0.87454377854954D-08
s(54,16)=    0.22758936660709D-07
s(55,16)=    0.83734593120842D-08
s(56,16)=    0.84285271391003D-08
s(57,16)=    0.35171485277021D-08
s(58,16)=    0.75832095687950D-09
s(59,16)=    0.43268168776164D-08
s(60,16)=   -0.95942808601991D-08
s(61,16)=    0.95403421274832D-08
s(62,16)=    0.76600691821439D-08
s(63,16)=   -0.14954024697104D-07
s(64,16)=    0.67348966750108D-08
s(65,16)=   -0.98517593621749D-08
s(66,16)=   -0.17651774758307D-07
s(67,16)=   -0.13048185705955D-07
s(68,16)=   -0.10775049205703D-07
s(69,16)=   -0.83608843251429D-08
s(70,16)=   -0.22850830348438D-07
s(71,16)=    0.22617862998312D-07
s(72,16)=   -0.11615177933685D-07
s(73,16)=    0.41249820870472D-08
s(74,16)=   -0.10011928436219D-07
s(75,16)=   -0.94657443577397D-08
s(76,16)=    0.20011823844755D-07
s(77,16)=   -0.33807858449145D-07
s(78,16)=    0.32016564260843D-07
s(79,16)=   -0.29026097371181D-07
s(80,16)=    0.22842416845736D-07
s(17,17)=    0.29248504614490D-06
s(18,17)=   -0.11033089461961D-05
s(19,17)=   -0.43514339734912D-06
s(20,17)=   -0.37811064797666D-06
s(21,17)=   -0.32773850863498D-06
s(22,17)=    0.23530841179301D-06
s(23,17)=    0.78411745594169D-06
s(24,17)=   -0.22796960437638D-07
s(25,17)=   -0.81775824768547D-06
s(26,17)=   -0.54832823201304D-07
s(27,17)=    0.36518698462665D-06
s(28,17)=    0.38124418820090D-06
s(29,17)=   -0.72547858089236D-07
s(30,17)=   -0.18838459479281D-06
s(31,17)=    0.47926572910547D-07
s(32,17)=   -0.20450598643054D-07
s(33,17)=    0.29165593938367D-07
s(34,17)=    0.69125694047977D-07
s(35,17)=   -0.87708918207654D-07
s(36,17)=   -0.11676243249569D-06
s(37,17)=    0.11419120111819D-06
s(38,17)=    0.82080852459118D-07
s(39,17)=    0.23240081183279D-07
s(40,17)=   -0.37879701407989D-07
s(41,17)=   -0.96433041548023D-07
s(42,17)=    0.38916730000796D-07
s(43,17)=    0.10053447126529D-06
s(44,17)=   -0.12019890556586D-07
s(45,17)=   -0.52972547052434D-07
s(46,17)=   -0.68212142440632D-07
s(47,17)=    0.26895107755749D-07
s(48,17)=    0.53081409381668D-07
s(49,17)=   -0.22951427108757D-07
s(50,17)=   -0.76202558407981D-07
s(51,17)=   -0.32349951224846D-08
s(52,17)=    0.11481715329839D-07
s(53,17)=   -0.59524666022700D-08
s(54,17)=    0.18813501639828D-07
s(55,17)=   -0.22940730979658D-07
s(56,17)=   -0.53348901787942D-08
s(57,17)=    0.17072079069305D-07
s(58,17)=    0.12881746245804D-07
s(59,17)=   -0.16563061398827D-07
s(60,17)=    0.19277069773992D-07
s(61,17)=   -0.39978892742045D-08
s(62,17)=   -0.16038939466946D-08
s(63,17)=    0.19128554124721D-07
s(64,17)=    0.11520684868882D-07
s(65,17)=    0.44875099081015D-09
s(66,17)=   -0.19609303908066D-07
s(67,17)=    0.18354623733187D-07
s(68,17)=    0.29301502723368D-07
s(69,17)=    0.28161118621374D-07
s(70,17)=    0.69095268314862D-08
s(71,17)=    0.41547099989174D-08
s(72,17)=    0.30357091060737D-07
s(73,17)=   -0.18609159497309D-07
s(74,17)=    0.43833534507557D-07
s(75,17)=   -0.22021390424856D-07
s(76,17)=    0.21569704607224D-07
s(77,17)=   -0.45260068167049D-08
s(78,17)=    0.30967907611003D-08
s(79,17)=   -0.15022630097898D-07
s(80,17)=    0.19496751617350D-07
s(18,18)=    0.34414061621346D-06
s(19,18)=    0.59127244706872D-06
s(20,18)=    0.72378561403062D-06
s(21,18)=    0.54121452249458D-06
s(22,18)=    0.52011524184896D-07
s(23,18)=   -0.60360508908718D-06
s(24,18)=   -0.54724968646783D-06
s(25,18)=    0.46999608659197D-07
s(26,18)=    0.51981672576762D-06
s(27,18)=    0.14846896174889D-06
s(28,18)=   -0.33779965470001D-06
s(29,18)=   -0.92157585525228D-07
s(30,18)=    0.53492902837147D-08
s(31,18)=    0.15359735482649D-06
s(32,18)=    0.98521227357963D-08
s(33,18)=    0.40659583739091D-07
s(34,18)=   -0.97787967430515D-07
s(35,18)=   -0.57606309931337D-07
s(36,18)=    0.14976973703166D-06
s(37,18)=    0.14384575033233D-06
s(38,18)=   -0.94105038960809D-07
s(39,18)=   -0.17770816278234D-06
s(40,18)=    0.62857306993194D-07
s(41,18)=    0.10564055594076D-06
s(42,18)=    0.26395600000873D-07
s(43,18)=   -0.22586048098759D-07
s(44,18)=   -0.19626803050804D-07
s(45,18)=   -0.26790029669805D-07
s(46,18)=    0.48875275415834D-07
s(47,18)=    0.77336313016548D-07
s(48,18)=   -0.28007132170749D-07
s(49,18)=   -0.56903480123028D-07
s(50,18)=    0.21859930554090D-07
s(51,18)=    0.26056945293627D-07
s(52,18)=   -0.17413777591014D-07
s(53,18)=   -0.65354897869392D-07
s(54,18)=   -0.94184192009345D-08
s(55,18)=    0.51900751533010D-07
s(56,18)=   -0.57408788836675D-08
s(57,18)=   -0.50384189794216D-08
s(58,18)=    0.32982260089339D-07
s(59,18)=    0.24308102652917D-07
s(60,18)=   -0.30314024384954D-07
s(61,18)=    0.46233654698603D-08
s(62,18)=   -0.10034917592270D-07
s(63,18)=   -0.25151953911156D-07
s(64,18)=   -0.27547219768778D-07
s(65,18)=   -0.12825473156738D-07
s(66,18)=   -0.23356444061669D-08
s(67,18)=    0.16196132089985D-07
s(68,18)=    0.11747662784129D-07
s(69,18)=    0.40196515918347D-07
s(70,18)=    0.84991016936327D-08
s(71,18)=   -0.16844142890715D-07
s(72,18)=    0.11605574724861D-07
s(73,18)=   -0.14117429773240D-07
s(74,18)=    0.10821358889875D-07
s(75,18)=   -0.91755950038818D-08
s(76,18)=    0.49454972151494D-08
s(77,18)=    0.88031137082583D-08
s(78,18)=    0.25538539368416D-08
s(79,18)=    0.33067900513529D-08
s(80,18)=   -0.26856994598150D-08
s(19,19)=   -0.83207181083117D-06
s(20,19)=   -0.19243210368148D-07
s(21,19)=   -0.77247799989911D-06
s(22,19)=   -0.41876045076903D-06
s(23,19)=    0.89358853172121D-07
s(24,19)=    0.40011464623022D-06
s(25,19)=   -0.56271728665351D-07
s(26,19)=   -0.16817956463516D-06
s(27,19)=    0.18805463102992D-06
s(28,19)=    0.48807625473716D-07
s(29,19)=   -0.27960321578934D-07
s(30,19)=    0.41704961518042D-07
s(31,19)=    0.10017377545951D-06
s(32,19)=   -0.99040199940945D-07
s(33,19)=   -0.19747524928826D-06
s(34,19)=    0.19999031060789D-06
s(35,19)=    0.73847531274333D-07
s(36,19)=   -0.77229614984630D-07
s(37,19)=   -0.66164882097243D-07
s(38,19)=   -0.57655690293397D-07
s(39,19)=    0.62929739483183D-07
s(40,19)=   -0.64102977423462D-09
s(41,19)=    0.17572199930565D-07
s(42,19)=   -0.42606475303314D-07
s(43,19)=   -0.79925299535828D-07
s(44,19)=   -0.24430081199529D-08
s(45,19)=   -0.40069459859954D-08
s(46,19)=   -0.90487367026726D-08
s(47,19)=   -0.46045184562653D-07
s(48,19)=    0.60516966760385D-08
s(49,19)=   -0.12727570691752D-07
s(50,19)=   -0.61817612034050D-08
s(51,19)=    0.17609953063319D-07
s(52,19)=   -0.43194132518131D-07
s(53,19)=   -0.52548737695946D-08
s(54,19)=   -0.12007188827659D-07
s(55,19)=    0.31582175037706D-08
s(56,19)=   -0.11497396834233D-07
s(57,19)=   -0.32165430513135D-09
s(58,19)=    0.32312688232480D-07
s(59,19)=    0.16256252932751D-08
s(60,19)=   -0.11796624470980D-07
s(61,19)=    0.24179493099147D-07
s(62,19)=   -0.12316628846491D-08
s(63,19)=    0.16770336871055D-07
s(64,19)=   -0.11401084406765D-07
s(65,19)=    0.16191494603284D-08
s(66,19)=    0.19016507804866D-08
s(67,19)=   -0.19853183621290D-08
s(68,19)=   -0.13392201323370D-07
s(69,19)=    0.27863708515734D-07
s(70,19)=    0.37308575408022D-07
s(71,19)=    0.14003822138422D-08
s(72,19)=    0.42743649987506D-07
s(73,19)=    0.24059852936671D-07
s(74,19)=    0.72128040264929D-09
s(75,19)=    0.36498605152585D-07
s(76,19)=   -0.14229114783173D-09
s(77,19)=    0.22176186379688D-07
s(78,19)=    0.16909636249029D-07
s(79,19)=    0.13548014160202D-07
s(80,19)=    0.16098574863556D-07
s(20,20)=    0.68495628517140D-07
s(21,20)=   -0.31271698464022D-06
s(22,20)=    0.57051643647335D-07
s(23,20)=    0.46838729499973D-08
s(24,20)=   -0.14902251350050D-06
s(25,20)=    0.11568289451042D-06
s(26,20)=    0.36424260503754D-06
s(27,20)=   -0.75694899244258D-07
s(28,20)=   -0.42145265720545D-06
s(29,20)=   -0.97179348884585D-07
s(30,20)=    0.18431388850136D-06
s(31,20)=    0.16621600245409D-06
s(32,20)=    0.48973155179554D-07
s(33,20)=   -0.52738506468826D-07
s(34,20)=   -0.89743613430438D-07
s(35,20)=   -0.90198666076864D-07
s(36,20)=    0.91884351911659D-07
s(37,20)=    0.93363810701203D-07
s(38,20)=   -0.33630822298069D-07
s(39,20)=   -0.10605986819334D-06
s(40,20)=   -0.10691444909545D-08
s(41,20)=    0.59276915592852D-07
s(42,20)=    0.15815236396154D-07
s(43,20)=   -0.61838158769288D-07
s(44,20)=   -0.51267588638259D-07
s(45,20)=   -0.28526483573614D-07
s(46,20)=    0.46277117704401D-07
s(47,20)=    0.51681277500150D-07
s(48,20)=   -0.65171325042152D-07
s(49,20)=   -0.21247258261835D-07
s(50,20)=   -0.43771569363668D-08
s(51,20)=    0.34169491656334D-08
s(52,20)=   -0.97534013510140D-08
s(53,20)=   -0.11347312251366D-07
s(54,20)=    0.15467152210774D-07
s(55,20)=    0.40292347301225D-07
s(56,20)=   -0.13079962180250D-07
s(57,20)=   -0.47052245767119D-07
s(58,20)=   -0.82495766579039D-08
s(59,20)=    0.97114458103715D-08
s(60,20)=    0.19913011710781D-07
s(61,20)=    0.36330274198173D-07
s(62,20)=   -0.13478687494604D-07
s(63,20)=   -0.26320561902926D-08
s(64,20)=   -0.16223451913486D-07
s(65,20)=   -0.10574281830612D-07
s(66,20)=    0.13839348502220D-07
s(67,20)=    0.11574284005268D-08
s(68,20)=   -0.11394258845984D-07
s(69,20)=   -0.68757042573032D-08
s(70,20)=    0.26808323992675D-07
s(71,20)=   -0.20990759642604D-07
s(72,20)=    0.17962079718463D-07
s(73,20)=    0.14028330809315D-07
s(74,20)=   -0.24187878478030D-08
s(75,20)=    0.17239582900246D-07
s(76,20)=   -0.18791192207054D-07
s(77,20)=    0.30074987185535D-07
s(78,20)=   -0.30113578588281D-07
s(79,20)=    0.39571391512958D-07
s(80,20)=   -0.34183756475334D-07
s(21,21)=    0.39509696270180D-06
s(22,21)=    0.58412852268006D-06
s(23,21)=    0.65443933539251D-06
s(24,21)=    0.49488306227768D-06
s(25,21)=    0.13005647364853D-06
s(26,21)=   -0.57736507305760D-06
s(27,21)=   -0.28716704730465D-06
s(28,21)=    0.15344083253979D-06
s(29,21)=    0.25221193020450D-06
s(30,21)=    0.77537550128164D-07
s(31,21)=    0.27125424267040D-07
s(32,21)=   -0.22822771898708D-06
s(33,21)=   -0.24125249759175D-06
s(34,21)=    0.14008778097270D-06
s(35,21)=    0.22976364858201D-06
s(36,21)=   -0.33609303386721D-07
s(37,21)=   -0.11099515652416D-06
s(38,21)=   -0.78867820909184D-07
s(39,21)=    0.67538516214823D-07
s(40,21)=    0.80979401970730D-07
s(41,21)=    0.18693961079579D-07
s(42,21)=   -0.76253978031094D-07
s(43,21)=   -0.10584263309521D-06
s(44,21)=    0.59555481450139D-07
s(45,21)=    0.62478998656977D-07
s(46,21)=   -0.10407191719497D-07
s(47,21)=   -0.57007502243515D-07
s(48,21)=   -0.46275770111026D-08
s(49,21)=    0.56834403221335D-08
s(50,21)=    0.39607818527299D-07
s(51,21)=    0.29133654190181D-07
s(52,21)=   -0.21171185313688D-07
s(53,21)=   -0.33129279461783D-07
s(54,21)=    0.11092356803145D-07
s(55,21)=   -0.36856992077205D-08
s(56,21)=   -0.61171426857546D-07
s(57,21)=   -0.44966367677811D-07
s(58,21)=    0.28939146376838D-07
s(59,21)=    0.27997575907820D-07
s(60,21)=   -0.43559182682407D-07
s(61,21)=   -0.18396042855272D-07
s(62,21)=    0.38422180278898D-07
s(63,21)=    0.19545978340560D-07
s(64,21)=    0.11044904958006D-07
s(65,21)=   -0.12556088956423D-07
s(66,21)=   -0.12450193849506D-07
s(67,21)=   -0.90707423538275D-08
s(68,21)=   -0.15648720052535D-07
s(69,21)=   -0.27818304947498D-07
s(70,21)=    0.52571730877252D-08
s(71,21)=   -0.78485105342175D-08
s(72,21)=   -0.24361309893413D-07
s(73,21)=    0.18537092819574D-07
s(74,21)=   -0.31528850699925D-07
s(75,21)=   -0.12295009610784D-07
s(76,21)=   -0.26879033180146D-07
s(77,21)=   -0.48641495118474D-08
s(78,21)=   -0.27262229899390D-07
s(79,21)=   -0.31296070040175D-08
s(80,21)=   -0.18355253625983D-07
s(22,22)=   -0.32472175304397D-06
s(23,22)=   -0.18794503883115D-06
s(24,22)=   -0.70495763436332D-06
s(25,22)=   -0.40985566966945D-06
s(26,22)=    0.82733514635286D-07
s(27,22)=    0.37746620514078D-06
s(28,22)=    0.69614132594899D-07
s(29,22)=   -0.13794766607726D-07
s(30,22)=   -0.22551051537638D-08
s(31,22)=   -0.14602895561313D-06
s(32,22)=   -0.75909083967744D-07
s(33,22)=    0.11406221093765D-06
s(34,22)=    0.12136146886712D-06
s(35,22)=   -0.98589998891183D-07
s(36,22)=   -0.10003304936051D-06
s(37,22)=    0.50298676434516D-07
s(38,22)=    0.60612538803716D-07
s(39,22)=   -0.13103547630008D-07
s(40,22)=   -0.30051950468690D-07
s(41,22)=   -0.26011654021801D-07
s(42,22)=   -0.81220118762205D-07
s(43,22)=    0.31474978731978D-07
s(44,22)=    0.28339826075713D-07
s(45,22)=   -0.38700086509404D-08
s(46,22)=   -0.46011104766624D-07
s(47,22)=   -0.27640275048026D-07
s(48,22)=   -0.10239800122767D-07
s(49,22)=    0.24734785595614D-07
s(50,22)=    0.14992914812251D-07
s(51,22)=    0.25393101822446D-07
s(52,22)=    0.38371746547476D-08
s(53,22)=   -0.19153595362559D-07
s(54,22)=    0.33803865148579D-07
s(55,22)=    0.56285859789713D-07
s(56,22)=   -0.38085937147590D-07
s(57,22)=   -0.63827478261483D-07
s(58,22)=    0.68254911727498D-08
s(59,22)=    0.31916386181269D-07
s(60,22)=   -0.28390268407220D-07
s(61,22)=   -0.51497624546058D-07
s(62,22)=   -0.10102801384017D-07
s(63,22)=    0.21202600642596D-07
s(64,22)=   -0.21960804438917D-07
s(65,22)=   -0.38562372640046D-08
s(66,22)=    0.50321543221399D-08
s(67,22)=   -0.11120594384824D-07
s(68,22)=   -0.34766929179416D-07
s(69,22)=    0.14560299724971D-07
s(70,22)=   -0.17208788676041D-07
s(71,22)=   -0.22412321730335D-07
s(72,22)=    0.34117018314396D-08
s(73,22)=    0.28254216720055D-07
s(74,22)=    0.17709979587694D-07
s(75,22)=    0.13675524334921D-07
s(76,22)=    0.33574552806071D-07
s(77,22)=   -0.49966310691079D-08
s(78,22)=    0.28734719794092D-07
s(79,22)=   -0.15906431570965D-07
s(80,22)=    0.35288814115458D-07
s(23,23)=    0.22719941596897D-06
s(24,23)=   -0.15396408753817D-07
s(25,23)=    0.42480863281236D-06
s(26,23)=    0.20336928244053D-06
s(27,23)=   -0.25760860961128D-06
s(28,23)=    0.13579469239490D-06
s(29,23)=    0.32187304031659D-06
s(30,23)=    0.80394556354472D-07
s(31,23)=   -0.17912333490658D-06
s(32,23)=    0.13220864045900D-06
s(33,23)=    0.11829940030533D-06
s(34,23)=   -0.47860233500407D-07
s(35,23)=    0.54187652460546D-07
s(36,23)=    0.89141708188153D-07
s(37,23)=    0.20676076210812D-07
s(38,23)=   -0.10261855286854D-06
s(39,23)=    0.36368105668424D-08
s(40,23)=    0.80548590305522D-07
s(41,23)=   -0.39402998138252D-07
s(42,23)=   -0.13044776330453D-07
s(43,23)=    0.35230822511696D-07
s(44,23)=   -0.38126285997105D-08
s(45,23)=   -0.41870291407728D-07
s(46,23)=   -0.40542175855106D-07
s(47,23)=    0.52512970896453D-07
s(48,23)=    0.64472339154468D-07
s(49,23)=   -0.13764138097832D-07
s(50,23)=   -0.23262961649640D-07
s(51,23)=   -0.54498090962075D-08
s(52,23)=   -0.40832137320425D-07
s(53,23)=    0.30883788965682D-07
s(54,23)=    0.59554143490002D-08
s(55,23)=    0.15826576402758D-07
s(56,23)=   -0.30546366389131D-08
s(57,23)=   -0.33710759479548D-07
s(58,23)=    0.42432292031115D-07
s(59,23)=    0.19055844739369D-07
s(60,23)=   -0.53493171836413D-07
s(61,23)=   -0.86122918322742D-08
s(62,23)=    0.42435124794369D-07
s(63,23)=    0.24086226827129D-08
s(64,23)=   -0.19299029335542D-07
s(65,23)=    0.27344753084321D-07
s(66,23)=   -0.64811007443060D-08
s(67,23)=    0.40723327050428D-08
s(68,23)=    0.12308992461489D-08
s(69,23)=   -0.13519766695142D-07
s(70,23)=    0.15952942060573D-07
s(71,23)=   -0.19862866658780D-07
s(72,23)=    0.21193483815864D-07
s(73,23)=    0.77746498722503D-08
s(74,23)=    0.46649362127192D-08
s(75,23)=    0.13571741248582D-07
s(76,23)=    0.77914367333071D-09
s(77,23)=    0.12620984254097D-07
s(78,23)=   -0.76843752977740D-08
s(79,23)=    0.14957902201148D-07
s(80,23)=   -0.11704671972042D-08
s(24,24)=    0.13365654513812D-06
s(25,24)=    0.10420440793780D-06
s(26,24)=    0.69892842255922D-07
s(27,24)=    0.20822951271181D-07
s(28,24)=    0.17328224754317D-06
s(29,24)=   -0.31771054193579D-06
s(30,24)=   -0.22658020777197D-06
s(31,24)=    0.11768365987851D-06
s(32,24)=    0.15392478445121D-06
s(33,24)=    0.14970977269238D-06
s(34,24)=   -0.18409701453663D-07
s(35,24)=   -0.13423591418765D-07
s(36,24)=   -0.71647913872669D-07
s(37,24)=   -0.15648959767529D-07
s(38,24)=    0.44608328778969D-07
s(39,24)=   -0.13572564408638D-07
s(40,24)=    0.44809188076787D-08
s(41,24)=   -0.85562804175540D-07
s(42,24)=   -0.17330890954671D-07
s(43,24)=   -0.75641223117804D-08
s(44,24)=    0.33642338348812D-07
s(45,24)=    0.22118487981987D-07
s(46,24)=   -0.76448411311210D-07
s(47,24)=   -0.59672118147122D-07
s(48,24)=    0.36943136668974D-07
s(49,24)=    0.27410232206408D-07
s(50,24)=   -0.36602861967874D-07
s(51,24)=   -0.49738771243184D-07
s(52,24)=    0.24689100919164D-07
s(53,24)=   -0.34229297600612D-08
s(54,24)=   -0.47632654023059D-08
s(55,24)=    0.22289622363861D-08
s(56,24)=    0.44689900647894D-08
s(57,24)=   -0.12928420382710D-07
s(58,24)=   -0.42174766304784D-08
s(59,24)=    0.21273250831510D-07
s(60,24)=   -0.33424298857639D-09
s(61,24)=   -0.12107587705232D-07
s(62,24)=    0.14317042771814D-09
s(63,24)=   -0.20993014525141D-07
s(64,24)=   -0.26778503005350D-07
s(65,24)=   -0.30292522435917D-07
s(66,24)=    0.17601994659046D-07
s(67,24)=    0.33698577624798D-07
s(68,24)=    0.12250204010052D-07
s(69,24)=    0.29217190072430D-07
s(70,24)=   -0.10819728157274D-07
s(71,24)=    0.39233527254074D-07
s(72,24)=    0.48275080073852D-08
s(73,24)=   -0.19264686446292D-07
s(74,24)=    0.17623540189089D-07
s(75,24)=   -0.43220085749331D-07
s(76,24)=    0.10621305199599D-07
s(77,24)=   -0.39080250417911D-07
s(78,24)=    0.21350367060357D-07
s(79,24)=   -0.46835778642193D-07
s(80,24)=    0.28507269098269D-07
s(25,25)=   -0.28343383249738D-06
s(26,25)=   -0.41376619591643D-06
s(27,25)=   -0.58673900423947D-06
s(28,25)=   -0.11103774475890D-06
s(29,25)=    0.28364746018767D-06
s(30,25)=    0.25965497175507D-06
s(31,25)=   -0.94438231555344D-08
s(32,25)=   -0.12000761429814D-06
s(33,25)=   -0.54415022292270D-07
s(34,25)=   -0.10164908575427D-06
s(35,25)=   -0.47046845571941D-07
s(36,25)=    0.11550555976001D-06
s(37,25)=    0.60618145789925D-07
s(38,25)=   -0.41389424022968D-07
s(39,25)=   -0.92966391031709D-07
s(40,25)=   -0.53736514318232D-07
s(41,25)=   -0.54497487858867D-07
s(42,25)=    0.72148313908372D-07
s(43,25)=    0.62334082759228D-07
s(44,25)=   -0.29218597413776D-07
s(45,25)=   -0.61029121259232D-07
s(46,25)=   -0.50983851941039D-07
s(47,25)=    0.18630917759058D-07
s(48,25)=    0.59003820157533D-07
s(49,25)=    0.28403235587668D-07
s(50,25)=   -0.22985892903289D-07
s(51,25)=   -0.33488421991056D-07
s(52,25)=   -0.12563713087659D-07
s(53,25)=    0.28423814410021D-07
s(54,25)=    0.90493803763162D-08
s(55,25)=   -0.34771934325339D-07
s(56,25)=   -0.86197138185889D-08
s(57,25)=   -0.50852900447310D-08
s(58,25)=    0.14205254460726D-07
s(59,25)=   -0.38952178709227D-08
s(60,25)=    0.18170442000334D-07
s(61,25)=    0.15252783788542D-07
s(62,25)=   -0.86902777480379D-08
s(63,25)=   -0.42928136365090D-07
s(64,25)=    0.29468803733683D-08
s(65,25)=   -0.93471377774845D-08
s(66,25)=   -0.10127788250041D-07
s(67,25)=    0.42852726412705D-07
s(68,25)=   -0.18592828646705D-08
s(69,25)=    0.16988875042302D-07
s(70,25)=   -0.18162904280886D-07
s(71,25)=   -0.29915694799725D-08
s(72,25)=   -0.17528084203899D-07
s(73,25)=   -0.25442757363903D-07
s(74,25)=   -0.25468049336641D-08
s(75,25)=   -0.31969511674933D-07
s(76,25)=    0.54675379007774D-09
s(77,25)=   -0.24783534205878D-07
s(78,25)=    0.25779347834228D-08
s(79,25)=   -0.32941556756683D-07
s(80,25)=    0.92235974771148D-08
s(26,26)=    0.43338187234426D-06
s(27,26)=    0.32255747503100D-06
s(28,26)=    0.51542611412935D-06
s(29,26)=    0.26605201857376D-06
s(30,26)=   -0.32663208942466D-06
s(31,26)=   -0.30707617643347D-07
s(32,26)=   -0.65448538173474D-07
s(33,26)=   -0.12366214261228D-06
s(34,26)=   -0.10310360826025D-06
s(35,26)=    0.39829968263619D-07
s(36,26)=    0.76407885525851D-07
s(37,26)=   -0.88826261416065D-07
s(38,26)=   -0.11340067571691D-06
s(39,26)=   -0.62812841030767D-08
s(40,26)=    0.12313744577302D-06
s(41,26)=    0.94051267435214D-07
s(42,26)=   -0.89213215251186D-07
s(43,26)=   -0.85416625541469D-07
s(44,26)=   -0.54529321002864D-07
s(45,26)=    0.17806700371047D-07
s(46,26)=    0.32277381280023D-07
s(47,26)=   -0.13854184904487D-07
s(48,26)=   -0.53769740665614D-07
s(49,26)=   -0.26659909255083D-07
s(50,26)=    0.33559732217126D-08
s(51,26)=    0.12995541672008D-08
s(52,26)=    0.40773157036140D-07
s(53,26)=    0.45692247832734D-07
s(54,26)=   -0.36732719547698D-07
s(55,26)=   -0.92533664791458D-08
s(56,26)=    0.86132244618483D-08
s(57,26)=    0.10846602869552D-07
s(58,26)=   -0.17247917928208D-07
s(59,26)=   -0.28413095500352D-07
s(60,26)=    0.22558628341088D-07
s(61,26)=   -0.41449878368046D-08
s(62,26)=    0.60769760888392D-09
s(63,26)=    0.19623277856483D-07
s(64,26)=   -0.14473895347379D-07
s(65,26)=    0.90417076513148D-08
s(66,26)=    0.16239691540043D-07
s(67,26)=   -0.14184759836242D-07
s(68,26)=   -0.87303978099637D-08
s(69,26)=    0.13089690945973D-07
s(70,26)=    0.62270462211982D-08
s(71,26)=   -0.95019319057695D-08
s(72,26)=    0.13705009796467D-07
s(73,26)=   -0.15879202025637D-07
s(74,26)=    0.52592361908252D-08
s(75,26)=    0.74244869289485D-08
s(76,26)=   -0.12607839461997D-08
s(77,26)=    0.31095151602912D-08
s(78,26)=    0.99793715366390D-09
s(79,26)=    0.51272824362602D-08
s(80,26)=   -0.57770178320492D-08
s(27,27)=   -0.34170443014580D-07
s(28,27)=    0.42909908753145D-07
s(29,27)=    0.70641170336061D-07
s(30,27)=   -0.20849406545547D-06
s(31,27)=   -0.13261660576630D-06
s(32,27)=   -0.10992100095684D-06
s(33,27)=    0.40746677241334D-07
s(34,27)=    0.11510647196811D-06
s(35,27)=   -0.97357586430325D-08
s(36,27)=   -0.61999880756536D-07
s(37,27)=   -0.62936756411612D-07
s(38,27)=    0.10171846329650D-06
s(39,27)=    0.12125582972906D-06
s(40,27)=   -0.16059653639408D-07
s(41,27)=   -0.50701162673004D-07
s(42,27)=   -0.43803028209664D-07
s(43,27)=    0.98345630252460D-07
s(44,27)=    0.37100350641786D-07
s(45,27)=   -0.75542753431233D-08
s(46,27)=   -0.31517137838158D-07
s(47,27)=   -0.12291518438071D-07
s(48,27)=    0.44004759412813D-07
s(49,27)=    0.24694589785261D-07
s(50,27)=   -0.23941178895794D-07
s(51,27)=   -0.22826126011045D-07
s(52,27)=   -0.66949746661043D-08
s(53,27)=    0.18135417630152D-07
s(54,27)=    0.29502084853047D-07
s(55,27)=   -0.41327952897796D-07
s(56,27)=   -0.96129223122568D-08
s(57,27)=    0.29900387474418D-07
s(58,27)=    0.79785180422400D-08
s(59,27)=   -0.26146583614980D-07
s(60,27)=   -0.12639543885775D-07
s(61,27)=    0.34889103817727D-07
s(62,27)=    0.14741697501236D-07
s(63,27)=   -0.14955246973965D-07
s(64,27)=    0.12973524043207D-07
s(65,27)=    0.33743032144455D-07
s(66,27)=    0.10104936959663D-07
s(67,27)=   -0.28782081297409D-07
s(68,27)=    0.13625943373120D-08
s(69,27)=    0.22799370217497D-07
s(70,27)=   -0.91103186516143D-08
s(71,27)=    0.11806369559136D-07
s(72,27)=    0.24960924930024D-07
s(73,27)=    0.13308222061162D-07
s(74,27)=    0.32700211987723D-07
s(75,27)=    0.39456678413889D-07
s(76,27)=    0.28884966557988D-07
s(77,27)=    0.50691069305749D-07
s(78,27)=    0.51318087982712D-08
s(79,27)=    0.56932129094163D-07
s(80,27)=   -0.53130371518947D-08
s(28,28)=   -0.58985837407318D-06
s(29,28)=   -0.22491706106929D-06
s(30,28)=   -0.33657612528292D-06
s(31,28)=   -0.17588007120398D-06
s(32,28)=    0.27198972538494D-06
s(33,28)=    0.18609546807248D-06
s(34,28)=   -0.14344940495241D-07
s(35,28)=    0.78461507290809D-08
s(36,28)=    0.19337648678170D-06
s(37,28)=   -0.72197020817564D-08
s(38,28)=   -0.13766557518626D-06
s(39,28)=   -0.33989003948112D-07
s(40,28)=    0.13716997268007D-06
s(41,28)=    0.43720861794535D-07
s(42,28)=   -0.53775665106866D-07
s(43,28)=   -0.49722185022982D-07
s(44,28)=    0.15179888848508D-08
s(45,28)=    0.49797683624109D-07
s(46,28)=    0.78626078282292D-07
s(47,28)=    0.29271105990369D-08
s(48,28)=   -0.12095485088284D-07
s(49,28)=    0.33592775142139D-08
s(50,28)=    0.14593614729016D-07
s(51,28)=    0.77588273282151D-08
s(52,28)=    0.16854138017509D-07
s(53,28)=    0.51531108406875D-07
s(54,28)=   -0.30953009796593D-07
s(55,28)=   -0.13159235460372D-07
s(56,28)=    0.12573303530113D-07
s(57,28)=    0.53297198998565D-07
s(58,28)=   -0.47285470038410D-07
s(59,28)=   -0.53724988753629D-07
s(60,28)=    0.32119661032479D-07
s(61,28)=    0.18521465315670D-08
s(62,28)=   -0.22863052872106D-07
s(63,28)=   -0.10582187950961D-08
s(64,28)=    0.31079660281112D-07
s(65,28)=    0.11364485225309D-07
s(66,28)=    0.22983251494080D-07
s(67,28)=    0.18454123546409D-07
s(68,28)=   -0.17287640656947D-07
s(69,28)=    0.18760264653712D-07
s(70,28)=    0.75169808875106D-08
s(71,28)=   -0.10219926273799D-08
s(72,28)=   -0.67842615858714D-08
s(73,28)=   -0.14976962407661D-07
s(74,28)=   -0.30469466429983D-07
s(75,28)=   -0.23386278766646D-07
s(76,28)=   -0.33580736734242D-07
s(77,28)=   -0.18394421066540D-07
s(78,28)=   -0.31594991341127D-07
s(79,28)=   -0.86844378634889D-08
s(80,28)=   -0.37054489276426D-07
s(29,29)=    0.15934179139910D-06
s(30,29)=    0.80261183497611D-07
s(31,29)=    0.20452087491675D-06
s(32,29)=    0.19302949858415D-06
s(33,29)=    0.31044324018296D-07
s(34,29)=    0.68372561049159D-07
s(35,29)=   -0.68406151061339D-07
s(36,29)=   -0.18550487603784D-06
s(37,29)=   -0.98042903845356D-07
s(38,29)=    0.22138827313363D-06
s(39,29)=    0.21394124499928D-06
s(40,29)=   -0.15416201266810D-06
s(41,29)=   -0.10576508454304D-06
s(42,29)=   -0.22805079313001D-07
s(43,29)=    0.14827037599014D-06
s(44,29)=    0.75196844827512D-07
s(45,29)=   -0.21389377822180D-07
s(46,29)=   -0.87293909692996D-07
s(47,29)=   -0.10737236474157D-06
s(48,29)=    0.44305518712465D-07
s(49,29)=    0.86990542087037D-07
s(50,29)=   -0.15392293253890D-07
s(51,29)=   -0.27392242860976D-07
s(52,29)=    0.40015795231051D-07
s(53,29)=    0.47976448579178D-08
s(54,29)=    0.43940406414663D-09
s(55,29)=   -0.63088793315440D-08
s(56,29)=    0.17116906947311D-07
s(57,29)=   -0.24362425804244D-07
s(58,29)=   -0.31087444982469D-07
s(59,29)=    0.71810373620219D-08
s(60,29)=    0.44337931570666D-07
s(61,29)=    0.30994772588096D-09
s(62,29)=   -0.15456911802245D-07
s(63,29)=    0.27462379375534D-07
s(64,29)=   -0.18164946905034D-07
s(65,29)=   -0.31507596095268D-07
s(66,29)=   -0.19550489228411D-07
s(67,29)=    0.28930785623828D-07
s(68,29)=   -0.11451115042695D-07
s(69,29)=   -0.58469224805794D-10
s(70,29)=    0.13406669034784D-07
s(71,29)=   -0.14069122767834D-07
s(72,29)=    0.97192394157168D-08
s(73,29)=    0.28957457754385D-08
s(74,29)=   -0.11911867022075D-07
s(75,29)=    0.29834427818979D-07
s(76,29)=   -0.15518302956492D-07
s(77,29)=    0.20382105419020D-07
s(78,29)=   -0.78095369888892D-08
s(79,29)=    0.23614790722157D-07
s(80,29)=   -0.18284196055722D-07
s(30,30)=   -0.73208354445994D-08
s(31,30)=    0.20503197108988D-06
s(32,30)=   -0.60664330159172D-07
s(33,30)=   -0.15149502720144D-06
s(34,30)=   -0.22995321782695D-06
s(35,30)=   -0.11486429179361D-06
s(36,30)=    0.19271900187518D-06
s(37,30)=    0.30515933963324D-06
s(38,30)=   -0.10739721056435D-06
s(39,30)=   -0.23837351416915D-06
s(40,30)=   -0.43570935938292D-07
s(41,30)=    0.10306543986188D-06
s(42,30)=    0.95073214474989D-07
s(43,30)=   -0.31670675362918D-07
s(44,30)=   -0.43072082570833D-07
s(45,30)=   -0.42628687639845D-07
s(46,30)=   -0.46591332486110D-07
s(47,30)=    0.36297872037229D-07
s(48,30)=    0.50066873364918D-08
s(49,30)=   -0.20773077169859D-07
s(50,30)=   -0.56230260521248D-07
s(51,30)=   -0.23225006093591D-07
s(52,30)=    0.13150371355644D-07
s(53,30)=    0.13647901560827D-08
s(54,30)=   -0.34376051509691D-07
s(55,30)=    0.35961991143393D-08
s(56,30)=    0.46380988195390D-07
s(57,30)=    0.97444094524642D-08
s(58,30)=   -0.18212691120108D-07
s(59,30)=   -0.39712009869295D-07
s(60,30)=    0.30781051371393D-07
s(61,30)=    0.51745294237553D-08
s(62,30)=   -0.12859144651056D-07
s(63,30)=   -0.28232237860941D-07
s(64,30)=   -0.58781004111777D-07
s(65,30)=   -0.12682708120815D-07
s(66,30)=    0.15713937458397D-07
s(67,30)=    0.13893147215378D-07
s(68,30)=   -0.14738832789576D-07
s(69,30)=    0.84698454522225D-08
s(70,30)=    0.63319613781522D-08
s(71,30)=   -0.11737291055526D-08
s(72,30)=    0.20959946701376D-07
s(73,30)=    0.17291696084040D-07
s(74,30)=    0.33146482887149D-07
s(75,30)=    0.26747011919294D-07
s(76,30)=    0.32677623283850D-07
s(77,30)=    0.30887204285341D-07
s(78,30)=    0.20456256872414D-07
s(79,30)=    0.22274599354938D-07
s(80,30)=    0.13660597104178D-07
s(31,31)=    0.78953381289128D-07
s(32,31)=    0.55092681530910D-07
s(33,31)=   -0.24066812936350D-07
s(34,31)=    0.83331006094216D-07
s(35,31)=    0.20605195393677D-07
s(36,31)=   -0.16601047334463D-06
s(37,31)=   -0.11457522585646D-06
s(38,31)=    0.55916272890031D-07
s(39,31)=    0.16616746720778D-06
s(40,31)=   -0.37589902744299D-07
s(41,31)=   -0.11581940689458D-06
s(42,31)=    0.44661598930658D-07
s(43,31)=    0.23115976793491D-07
s(44,31)=    0.52789006379119D-08
s(45,31)=   -0.38801867403663D-07
s(46,31)=    0.58426457641535D-07
s(47,31)=   -0.73087469795559D-08
s(48,31)=   -0.15876001289290D-07
s(49,31)=    0.24498735694810D-07
s(50,31)=   -0.12124946374761D-07
s(51,31)=   -0.59637546214590D-07
s(52,31)=    0.15716912069639D-07
s(53,31)=    0.27156878636963D-07
s(54,31)=   -0.20598529734341D-07
s(55,31)=   -0.38168653542138D-07
s(56,31)=    0.18160363777274D-08
s(57,31)=    0.17309169075534D-08
s(58,31)=   -0.57335222802739D-07
s(59,31)=    0.18564145366185D-07
s(60,31)=    0.25766229579466D-07
s(61,31)=   -0.10333355906040D-07
s(62,31)=   -0.18464747507191D-07
s(63,31)=    0.13195542467325D-07
s(64,31)=    0.11997271372907D-08
s(65,31)=    0.30458308482371D-07
s(66,31)=   -0.52514627551796D-08
s(67,31)=   -0.32599509915357D-07
s(68,31)=    0.31415332439130D-07
s(69,31)=   -0.25740231453981D-07
s(70,31)=   -0.16680113918527D-08
s(71,31)=    0.11995163703205D-08
s(72,31)=   -0.63110056965274D-08
s(73,31)=   -0.27153087065670D-07
s(74,31)=   -0.19474317437066D-07
s(75,31)=   -0.20038597783690D-07
s(76,31)=   -0.32275924064217D-07
s(77,31)=   -0.16759991878209D-07
s(78,31)=   -0.35136289991141D-07
s(79,31)=   -0.23546381236919D-07
s(80,31)=   -0.36693928431794D-07
s(32,32)=    0.89374420888902D-07
s(33,32)=   -0.39771492440454D-06
s(34,32)=    0.11469182814266D-06
s(35,32)=    0.25395957524143D-06
s(36,32)=    0.42129618578725D-07
s(37,32)=    0.98208473981999D-07
s(38,32)=    0.99469730680367D-08
s(39,32)=   -0.25745403534494D-06
s(40,32)=   -0.41443252825890D-07
s(41,32)=    0.16820069751870D-06
s(42,32)=    0.10308232626339D-06
s(43,32)=   -0.70038064535491D-07
s(44,32)=   -0.12196361595307D-06
s(45,32)=    0.42272108497804D-07
s(46,32)=    0.59377348437907D-07
s(47,32)=    0.18906628013296D-07
s(48,32)=    0.13649977026942D-07
s(49,32)=   -0.34490511565167D-07
s(50,32)=    0.43027917035399D-09
s(51,32)=    0.50270922076139D-07
s(52,32)=    0.46357304116044D-07
s(53,32)=   -0.24282795615853D-08
s(54,32)=   -0.90643338662239D-08
s(55,32)=   -0.23829540610693D-07
s(56,32)=    0.11711589064585D-07
s(57,32)=    0.35754022645868D-07
s(58,32)=   -0.48182753078852D-08
s(59,32)=   -0.23753633638613D-07
s(60,32)=   -0.22342169368798D-07
s(61,32)=    0.24521309772462D-07
s(62,32)=    0.39161940497231D-07
s(63,32)=   -0.21554336265665D-07
s(64,32)=   -0.28586080845642D-08
s(65,32)=   -0.97613286544951D-08
s(66,32)=   -0.33763392592834D-07
s(67,32)=   -0.38442901236260D-08
s(68,32)=   -0.24250680463403D-08
s(69,32)=   -0.97733825080326D-08
s(70,32)=    0.21796544995104D-07
s(71,32)=   -0.63533297575582D-09
s(72,32)=   -0.24344198631389D-07
s(73,32)=   -0.14144615871724D-08
s(74,32)=   -0.27606198840908D-07
s(75,32)=   -0.26461670063636D-07
s(76,32)=   -0.14595660921976D-07
s(77,32)=   -0.25425669305366D-07
s(78,32)=   -0.81070748800349D-08
s(79,32)=   -0.15184683821167D-07
s(80,32)=   -0.11268838420761D-07
s(33,33)=   -0.94208598006704D-07
s(34,33)=    0.29799617179318D-08
s(35,33)=   -0.22242110984293D-06
s(36,33)=   -0.29066411984174D-06
s(37,33)=   -0.20209988882640D-06
s(38,33)=    0.17897713731673D-07
s(39,33)=    0.13953494947431D-06
s(40,33)=    0.11413168482978D-06
s(41,33)=   -0.48813662280767D-07
s(42,33)=   -0.16926140355608D-06
s(43,33)=   -0.49953698561631D-07
s(44,33)=    0.52300068682533D-07
s(45,33)=    0.53569128202290D-07
s(46,33)=    0.16136373991569D-07
s(47,33)=   -0.14485682922362D-07
s(48,33)=   -0.65624463147833D-07
s(49,33)=   -0.19381554421885D-07
s(50,33)=    0.58309224409575D-07
s(51,33)=    0.98195847838631D-08
s(52,33)=   -0.26550731825424D-08
s(53,33)=    0.97682209709751D-08
s(54,33)=    0.42960391153541D-08
s(55,33)=    0.12797933236365D-07
s(56,33)=    0.38541856808679D-07
s(57,33)=   -0.12234250882913D-07
s(58,33)=    0.15517306753324D-07
s(59,33)=   -0.31018353753256D-07
s(60,33)=    0.22914717662037D-07
s(61,33)=    0.15804975374955D-07
s(62,33)=    0.50757151028676D-08
s(63,33)=   -0.80321210231905D-08
s(64,33)=   -0.12636820680965D-08
s(65,33)=   -0.17081053415114D-07
s(66,33)=   -0.16051416121923D-08
s(67,33)=   -0.28763013588895D-07
s(68,33)=   -0.18131621216844D-07
s(69,33)=   -0.23869942652726D-08
s(70,33)=   -0.61060341390310D-08
s(71,33)=    0.15176146381324D-07
s(72,33)=    0.40800040051723D-08
s(73,33)=   -0.14707421372003D-07
s(74,33)=   -0.88115310774511D-08
s(75,33)=   -0.19825518153810D-07
s(76,33)=    0.64698444424280D-08
s(77,33)=   -0.56801384886818D-08
s(78,33)=   -0.11924496262571D-07
s(79,33)=    0.18518675387124D-07
s(80,33)=   -0.19178282810397D-07
s(34,34)=    0.57891361833442D-07
s(35,34)=    0.35410000719729D-06
s(36,34)=    0.14032205000435D-06
s(37,34)=    0.27999725729822D-07
s(38,34)=    0.18242723452292D-06
s(39,34)=   -0.73527413966796D-07
s(40,34)=   -0.21144239516408D-06
s(41,34)=   -0.56269656280633D-07
s(42,34)=    0.37385595963484D-07
s(43,34)=   -0.61348664789110D-08
s(44,34)=   -0.35479484011716D-07
s(45,34)=    0.25877502023238D-07
s(46,34)=   -0.22933272566050D-07
s(47,34)=   -0.10917202035018D-06
s(48,34)=   -0.45419133217149D-07
s(49,34)=    0.41834613424255D-07
s(50,34)=    0.45564394341709D-07
s(51,34)=    0.29839224484537D-07
s(52,34)=   -0.21684207261625D-07
s(53,34)=   -0.52560461950635D-07
s(54,34)=   -0.49734569633855D-08
s(55,34)=    0.62609322210806D-07
s(56,34)=    0.23037019399123D-08
s(57,34)=    0.28567371343927D-08
s(58,34)=   -0.48283994827938D-08
s(59,34)=    0.85945902364989D-08
s(60,34)=   -0.16766221308052D-07
s(61,34)=    0.20397296405930D-07
s(62,34)=   -0.16348347560272D-07
s(63,34)=   -0.16633471743895D-07
s(64,34)=   -0.20720067762420D-07
s(65,34)=   -0.21684071991585D-07
s(66,34)=    0.43792547302982D-08
s(67,34)=   -0.44138759859463D-08
s(68,34)=    0.70357085654840D-08
s(69,34)=    0.11822517036222D-07
s(70,34)=   -0.12881206313722D-07
s(71,34)=   -0.19482515896062D-07
s(72,34)=   -0.67663252746500D-08
s(73,34)=   -0.27619766132433D-08
s(74,34)=    0.65622048076009D-08
s(75,34)=    0.13198689067665D-07
s(76,34)=   -0.72285366979372D-08
s(77,34)=    0.11673329617847D-07
s(78,34)=    0.58468649983396D-09
s(79,34)=   -0.75417288260744D-08
s(80,34)=    0.13374343775362D-07
s(35,35)=    0.17574893437826D-07
s(36,35)=   -0.18917698132751D-06
s(37,35)=    0.10402270997036D-06
s(38,35)=    0.66196197312449D-07
s(39,35)=    0.10897013392005D-06
s(40,35)=    0.71306369032483D-07
s(41,35)=    0.21379564891250D-07
s(42,35)=   -0.45856563658787D-07
s(43,35)=    0.22318982455819D-07
s(44,35)=    0.95498176986744D-07
s(45,35)=    0.21977868445240D-07
s(46,35)=   -0.72451901992539D-07
s(47,35)=   -0.31058301957535D-07
s(48,35)=    0.16324095077608D-07
s(49,35)=    0.10429426214531D-07
s(50,35)=    0.94221774231592D-08
s(51,35)=    0.12176235287938D-07
s(52,35)=   -0.29172667049278D-07
s(53,35)=   -0.16589494701342D-07
s(54,35)=    0.26089459167721D-07
s(55,35)=    0.78345423397103D-07
s(56,35)=   -0.45583432726474D-08
s(57,35)=   -0.49458722169173D-07
s(58,35)=   -0.24416230031301D-08
s(59,35)=   -0.35624603264824D-08
s(60,35)=   -0.27603732186268D-08
s(61,35)=    0.91364997384810D-08
s(62,35)=    0.50467062635717D-09
s(63,35)=   -0.30528611453612D-07
s(64,35)=   -0.83018345618891D-08
s(65,35)=   -0.81366284734688D-08
s(66,35)=   -0.61185123425041D-08
s(67,35)=   -0.14965236741050D-07
s(68,35)=    0.57111160293270D-08
s(69,35)=   -0.13278372297938D-07
s(70,35)=   -0.24891598456516D-07
s(71,35)=    0.45291250627892D-08
s(72,35)=   -0.25883497180767D-07
s(73,35)=   -0.19594072388370D-07
s(74,35)=   -0.31778035058203D-08
s(75,35)=   -0.26856121714990D-07
s(76,35)=   -0.80253713217046D-08
s(77,35)=   -0.30881564979065D-08
s(78,35)=   -0.46166508756688D-08
s(79,35)=   -0.12906559027396D-07
s(80,35)=    0.13935082715953D-07
s(36,36)=   -0.20219346710156D-06
s(37,36)=   -0.18565703689859D-06
s(38,36)=   -0.18586053945669D-06
s(39,36)=   -0.71406867219384D-07
s(40,36)=   -0.63776123497810D-07
s(41,36)=   -0.47554070004202D-07
s(42,36)=    0.79512097443865D-07
s(43,36)=    0.11649586191650D-06
s(44,36)=   -0.50904431344214D-07
s(45,36)=   -0.75423388300267D-07
s(46,36)=    0.58602766101646D-07
s(47,36)=    0.46850361082441D-07
s(48,36)=    0.32015762470739D-08
s(49,36)=   -0.25335742724641D-07
s(50,36)=    0.26346357568252D-07
s(51,36)=    0.10879802316803D-09
s(52,36)=   -0.66981996635857D-08
s(53,36)=   -0.29306027844671D-08
s(54,36)=    0.16351763245747D-07
s(55,36)=    0.25856474578397D-07
s(56,36)=    0.13700179370321D-08
s(57,36)=   -0.27868568299405D-07
s(58,36)=    0.84458741057917D-08
s(59,36)=    0.37746613907334D-07
s(60,36)=    0.17836526402338D-07
s(61,36)=   -0.29162421924768D-07
s(62,36)=   -0.32836453597881D-07
s(63,36)=   -0.25289019775190D-07
s(64,36)=    0.12202549479010D-07
s(65,36)=    0.19022429206720D-07
s(66,36)=    0.99826883595768D-08
s(67,36)=   -0.12786222890282D-07
s(68,36)=    0.37875948566872D-07
s(69,36)=   -0.29062152850641D-07
s(70,36)=    0.39182214734233D-08
s(71,36)=    0.96862530007850D-09
s(72,36)=    0.90080155735125D-08
s(73,36)=    0.80087161731193D-08
s(74,36)=    0.38858346341065D-08
s(75,36)=    0.11475665208163D-07
s(76,36)=    0.19320641646836D-07
s(77,36)=   -0.95972419672990D-10
s(78,36)=    0.12611076496036D-07
s(79,36)=    0.15892552740599D-07
s(80,36)=   -0.61614213870784D-08
s(37,37)=    0.82499963689137D-07
s(38,37)=    0.32579548514251D-06
s(39,37)=    0.82129925904680D-07
s(40,37)=    0.11123696638987D-06
s(41,37)=    0.10814123499530D-06
s(42,37)=    0.20875184329936D-07
s(43,37)=   -0.53648197571687D-07
s(44,37)=   -0.30090988672809D-07
s(45,37)=    0.32274694930966D-07
s(46,37)=    0.20143534109540D-07
s(47,37)=    0.74917476342432D-08
s(48,37)=    0.40429547752282D-09
s(49,37)=    0.27519728689204D-09
s(50,37)=   -0.19402358295504D-07
s(51,37)=    0.22498893646525D-07
s(52,37)=   -0.10868384024413D-07
s(53,37)=   -0.16461898658393D-07
s(54,37)=   -0.26565798924259D-07
s(55,37)=    0.10062120351788D-07
s(56,37)=   -0.35329493726135D-07
s(57,37)=   -0.34725348468397D-07
s(58,37)=    0.38650836668818D-07
s(59,37)=    0.47928999771503D-07
s(60,37)=   -0.21037700929038D-07
s(61,37)=   -0.23348386689476D-07
s(62,37)=   -0.33175121606318D-08
s(63,37)=    0.10813339572817D-07
s(64,37)=   -0.98037281949979D-08
s(65,37)=    0.16641096036252D-07
s(66,37)=   -0.20239999348763D-08
s(67,37)=   -0.86697474572006D-08
s(68,37)=    0.47546227967430D-07
s(69,37)=   -0.11331337838453D-08
s(70,37)=    0.25441637443214D-07
s(71,37)=    0.82327691620187D-08
s(72,37)=    0.13909970929050D-07
s(73,37)=    0.19755954609375D-07
s(74,37)=    0.12156382486042D-07
s(75,37)=    0.61724690151772D-07
s(76,37)=   -0.72090808068583D-08
s(77,37)=    0.72313923016550D-07
s(78,37)=    0.11118751637188D-07
s(79,37)=    0.55826430259919D-07
s(80,37)=    0.16853621311489D-07
s(38,38)=    0.27223825235871D-07
s(39,38)=   -0.11797781782116D-06
s(40,38)=   -0.27052502411252D-07
s(41,38)=   -0.11205128441237D-06
s(42,38)=   -0.13264118511390D-06
s(43,38)=    0.18137299951394D-07
s(44,38)=    0.25366425674263D-07
s(45,38)=   -0.39955814072124D-07
s(46,38)=    0.40592068927744D-07
s(47,38)=    0.66094637996260D-07
s(48,38)=   -0.59480576006430D-07
s(49,38)=   -0.84877407686127D-07
s(50,38)=    0.27760798622301D-07
s(51,38)=    0.52906318561677D-07
s(52,38)=   -0.28076115668295D-07
s(53,38)=    0.67679253181968D-08
s(54,38)=    0.25023820910351D-07
s(55,38)=   -0.11048359651356D-07
s(56,38)=   -0.76902740526226D-08
s(57,38)=    0.15112418885969D-07
s(58,38)=   -0.49008942050067D-08
s(59,38)=   -0.42835038959172D-07
s(60,38)=   -0.64952785851025D-08
s(61,38)=   -0.45589164370352D-09
s(62,38)=    0.51483704670460D-07
s(63,38)=    0.15181653233649D-07
s(64,38)=   -0.12355377536275D-08
s(65,38)=    0.74811315736456D-08
s(66,38)=    0.13310719431886D-07
s(67,38)=    0.35815003168421D-07
s(68,38)=    0.38129468651773D-07
s(69,38)=   -0.12798508066516D-07
s(70,38)=    0.20274720379769D-07
s(71,38)=   -0.19832787862815D-07
s(72,38)=    0.20011499089307D-07
s(73,38)=    0.18873273793084D-08
s(74,38)=    0.11229187536925D-07
s(75,38)=    0.47505456955470D-07
s(76,38)=    0.73479965358184D-08
s(77,38)=    0.22438227895970D-07
s(78,38)=    0.22393204715969D-07
s(79,38)=    0.33354869887151D-07
s(80,38)=    0.34794547967602D-08
s(39,39)=   -0.23100635007645D-06
s(40,39)=   -0.27132087150901D-07
s(41,39)=   -0.29897783129876D-07
s(42,39)=    0.45265121596468D-09
s(43,39)=    0.42604362762869D-07
s(44,39)=    0.39685524039415D-07
s(45,39)=   -0.23165310740645D-08
s(46,39)=    0.48756411069904D-08
s(47,39)=   -0.37235630952126D-07
s(48,39)=    0.34506784357156D-07
s(49,39)=    0.63774228923100D-07
s(50,39)=    0.37698756216119D-07
s(51,39)=   -0.63429968977591D-07
s(52,39)=   -0.34314284304191D-07
s(53,39)=   -0.14510670384084D-08
s(54,39)=    0.40910022347703D-07
s(55,39)=   -0.10674556987582D-07
s(56,39)=    0.88658932627125D-09
s(57,39)=    0.11712834589807D-07
s(58,39)=    0.27155416083045D-07
s(59,39)=    0.22173455156329D-07
s(60,39)=   -0.47702136132581D-07
s(61,39)=   -0.81528797843478D-08
s(62,39)=    0.31748754541902D-07
s(63,39)=    0.32090501264960D-07
s(64,39)=   -0.25061816495101D-07
s(65,39)=   -0.29457360270732D-08
s(66,39)=   -0.72222793814467D-08
s(67,39)=   -0.26685180842147D-07
s(68,39)=   -0.28893311138220D-07
s(69,39)=    0.49541153990283D-08
s(70,39)=    0.20220364716010D-08
s(71,39)=    0.17994294427792D-08
s(72,39)=    0.38180894115193D-08
s(73,39)=   -0.81797098238474D-08
s(74,39)=    0.16597033668737D-07
s(75,39)=    0.12839057150194D-07
s(76,39)=    0.19423668585256D-08
s(77,39)=    0.23922240782047D-08
s(78,39)=    0.20899374431030D-07
s(79,39)=    0.10939143624707D-07
s(80,39)=    0.33971126423981D-07
s(40,40)=    0.37633854063966D-07
s(41,40)=    0.15948666943187D-06
s(42,40)=    0.13930277589366D-06
s(43,40)=    0.66549248594838D-07
s(44,40)=    0.30909904541082D-07
s(45,40)=   -0.49584762111330D-07
s(46,40)=   -0.44018295076622D-07
s(47,40)=    0.30405240826239D-07
s(48,40)=   -0.39008575624906D-07
s(49,40)=   -0.43461828371225D-07
s(50,40)=   -0.22996589710628D-07
s(51,40)=    0.21230992940219D-08
s(52,40)=    0.25943906999772D-07
s(53,40)=   -0.36429058564067D-08
s(54,40)=   -0.67648234541438D-08
s(55,40)=    0.16197265528847D-08
s(56,40)=    0.15403161449800D-07
s(57,40)=    0.24091423149443D-07
s(58,40)=   -0.25958650400231D-07
s(59,40)=   -0.30090295714212D-07
s(60,40)=    0.22328902928061D-07
s(61,40)=    0.87894198549764D-08
s(62,40)=    0.16722121284547D-07
s(63,40)=    0.24822384736756D-08
s(64,40)=   -0.13224757347678D-07
s(65,40)=   -0.15709124253837D-07
s(66,40)=    0.33822511715929D-08
s(67,40)=   -0.51688211527445D-08
s(68,40)=    0.10228655575234D-08
s(69,40)=   -0.84616272515865D-08
s(70,40)=   -0.14335359703049D-07
s(71,40)=    0.10305188346012D-07
s(72,40)=   -0.17198291224937D-07
s(73,40)=    0.14138217077587D-07
s(74,40)=    0.88432263233331D-08
s(75,40)=   -0.51690588014883D-09
s(76,40)=    0.70026887954992D-08
s(77,40)=   -0.63782994128776D-08
s(78,40)=   -0.72332619059323D-09
s(79,40)=   -0.42733879792314D-08
s(80,40)=   -0.53694026424618D-08
s(41,41)=    0.16361583052024D-06
s(42,41)=   -0.14270106868287D-06
s(43,41)=   -0.21622197626879D-07
s(44,41)=   -0.63359839274286D-07
s(45,41)=   -0.13327859682953D-06
s(46,41)=   -0.69118344660528D-07
s(47,41)=    0.55641732790067D-07
s(48,41)=    0.48222948428843D-08
s(49,41)=    0.94428120484847D-08
s(50,41)=    0.88766067357840D-08
s(51,41)=   -0.42503796506277D-07
s(52,41)=   -0.72736915711942D-07
s(53,41)=    0.48149716729819D-07
s(54,41)=    0.54294075001544D-07
s(55,41)=   -0.29081387916231D-07
s(56,41)=   -0.38764691292697D-07
s(57,41)=    0.22870772994382D-08
s(58,41)=    0.54967585178154D-08
s(59,41)=    0.73958775748365D-08
s(60,41)=    0.17255365571717D-08
s(61,41)=   -0.31948069897762D-07
s(62,41)=    0.97542386001343D-08
s(63,41)=    0.28835943902525D-07
s(64,41)=    0.10094744488281D-07
s(65,41)=   -0.33852368620708D-08
s(66,41)=   -0.15899017687540D-07
s(67,41)=   -0.77570233005493D-08
s(68,41)=   -0.89675282337751D-08
s(69,41)=   -0.10455947048287D-07
s(70,41)=    0.12388769941907D-07
s(71,41)=    0.15698764962763D-07
s(72,41)=   -0.36472860036889D-08
s(73,41)=    0.79706161045882D-08
s(74,41)=    0.28522677490986D-08
s(75,41)=   -0.31802752389756D-08
s(76,41)=   -0.48429000222026D-08
s(77,41)=    0.11117759385866D-07
s(78,41)=    0.57121298774040D-08
s(79,41)=    0.55007918672074D-08
s(80,41)=    0.10741095553646D-07
s(42,42)=   -0.94068189248922D-07
s(43,42)=   -0.78475574428159D-07
s(44,42)=   -0.21291124423595D-07
s(45,42)=    0.72712292796957D-07
s(46,42)=    0.15721274969429D-06
s(47,42)=    0.58508120683539D-07
s(48,42)=    0.50802973859490D-07
s(49,42)=   -0.12163182901950D-07
s(50,42)=   -0.48060046663381D-07
s(51,42)=    0.21016277449898D-07
s(52,42)=    0.70639423628635D-07
s(53,42)=   -0.13442592402474D-07
s(54,42)=   -0.68730543364773D-07
s(55,42)=   -0.26145802044600D-07
s(56,42)=    0.52023033546749D-07
s(57,42)=    0.55216866997477D-07
s(58,42)=   -0.89221508510052D-09
s(59,42)=   -0.28969136741613D-07
s(60,42)=   -0.47262881525008D-08
s(61,42)=    0.49644823278567D-08
s(62,42)=    0.14674788264441D-07
s(63,42)=   -0.13846385454179D-07
s(64,42)=   -0.14302505306491D-07
s(65,42)=    0.24585348266529D-08
s(66,42)=   -0.19036969477211D-07
s(67,42)=    0.19388854163557D-07
s(68,42)=    0.28229706002095D-07
s(69,42)=   -0.47726239636899D-08
s(70,42)=   -0.14550935657971D-07
s(71,42)=    0.49995071977586D-08
s(72,42)=    0.14730620563668D-07
s(73,42)=   -0.19699606456276D-07
s(74,42)=    0.45987930383571D-08
s(75,42)=   -0.12764967951112D-07
s(76,42)=   -0.65963694759358D-08
s(77,42)=    0.13145141507264D-08
s(78,42)=   -0.13092518502489D-07
s(79,42)=    0.22064865780193D-07
s(80,42)=   -0.24244728742283D-07
s(43,43)=   -0.17175615617360D-07
s(44,43)=    0.10258719927514D-06
s(45,43)=   -0.93605270868609D-07
s(46,43)=   -0.29653575732994D-07
s(47,43)=    0.27866959534308D-07
s(48,43)=   -0.67860105464752D-08
s(49,43)=   -0.79717452180601D-07
s(50,43)=    0.17773008164420D-07
s(51,43)=    0.78598952090166D-08
s(52,43)=   -0.57943976379675D-07
s(53,43)=    0.23321268480130D-07
s(54,43)=    0.62760408057170D-07
s(55,43)=   -0.20315736754032D-07
s(56,43)=   -0.50337164678803D-07
s(57,43)=   -0.88098873811898D-08
s(58,43)=   -0.20279041166094D-07
s(59,43)=   -0.18403190307266D-07
s(60,43)=    0.14374202707386D-07
s(61,43)=   -0.39945112984150D-07
s(62,43)=   -0.25125706256732D-07
s(63,43)=    0.11054654858016D-07
s(64,43)=    0.60269555058883D-08
s(65,43)=   -0.13495157784041D-07
s(66,43)=   -0.87315372152658D-08
s(67,43)=    0.12785861073393D-07
s(68,43)=    0.17894615605784D-07
s(69,43)=    0.79559030798355D-08
s(70,43)=   -0.17732617097839D-07
s(71,43)=   -0.25969192658430D-07
s(72,43)=    0.10460427938408D-07
s(73,43)=   -0.25004808857820D-07
s(74,43)=   -0.63896670164526D-08
s(75,43)=    0.51701135999391D-08
s(76,43)=   -0.21221841431201D-08
s(77,43)=    0.11158689146173D-07
s(78,43)=    0.18183513032212D-07
s(79,43)=   -0.22844797597350D-07
s(80,43)=    0.14056355185349D-07
s(44,44)=    0.19649743701830D-06
s(45,44)=   -0.57499258737927D-08
s(46,44)=    0.52681752369751D-07
s(47,44)=   -0.17819888430494D-08
s(48,44)=   -0.14807072952223D-06
s(49,44)=   -0.48971146987035D-07
s(50,44)=    0.50108159887150D-07
s(51,44)=    0.19588526493264D-07
s(52,44)=    0.32170405238601D-07
s(53,44)=    0.65214326301844D-08
s(54,44)=   -0.64208210163437D-07
s(55,44)=   -0.22741537933769D-07
s(56,44)=    0.75815480934254D-07
s(57,44)=    0.23262449193373D-07
s(58,44)=    0.17883947143873D-08
s(59,44)=   -0.36630210562996D-07
s(60,44)=   -0.35269048564885D-07
s(61,44)=    0.46518011697369D-07
s(62,44)=   -0.75934903981466D-09
s(63,44)=   -0.84924522238898D-07
s(64,44)=    0.70285558573001D-07
s(65,44)=    0.38084206122098D-07
s(66,44)=   -0.55269253309188D-07
s(67,44)=    0.47795929245733D-07
s(68,44)=    0.10429536893844D-07
s(69,44)=   -0.39285221715527D-07
s(70,44)=   -0.15979873468228D-07
s(71,44)=    0.70623045381242D-08
s(72,44)=   -0.28306028666131D-07
s(73,44)=   -0.20005201856888D-07
s(74,44)=    0.18800919344200D-07
s(75,44)=   -0.28527082293398D-07
s(76,44)=    0.38386032796612D-07
s(77,44)=   -0.38906574388033D-08
s(78,44)=    0.45375767288144D-08
s(79,44)=    0.56513112459831D-07
s(80,44)=   -0.19255397328033D-07
s(45,45)=   -0.10865082700048D-06
s(46,45)=    0.22449081799764D-07
s(47,45)=    0.11452550163554D-07
s(48,45)=    0.71124305689779D-07
s(49,45)=    0.52740137523616D-07
s(50,45)=    0.13691364917273D-07
s(51,45)=   -0.92681462753584D-08
s(52,45)=   -0.41161019223283D-07
s(53,45)=   -0.17672853921607D-07
s(54,45)=    0.18050146388353D-07
s(55,45)=    0.25155932496392D-07
s(56,45)=   -0.82071719711664D-08
s(57,45)=   -0.36716978291219D-07
s(58,45)=   -0.15962763701682D-07
s(59,45)=    0.20212472566845D-07
s(60,45)=    0.36085656631687D-07
s(61,45)=   -0.62079380929853D-08
s(62,45)=    0.61419783719518D-08
s(63,45)=    0.84355508741461D-08
s(64,45)=   -0.75785809211492D-08
s(65,45)=   -0.31540365747449D-08
s(66,45)=    0.13059243935779D-07
s(67,45)=    0.26074913879975D-07
s(68,45)=    0.28232459150081D-08
s(69,45)=   -0.26158484508234D-08
s(70,45)=   -0.23761634344636D-08
s(71,45)=    0.10813245864093D-07
s(72,45)=    0.33671248775672D-07
s(73,45)=    0.25201712831274D-07
s(74,45)=    0.46237802256802D-08
s(75,45)=    0.13216404078551D-07
s(76,45)=    0.14791030326728D-08
s(77,45)=    0.20800712258151D-07
s(78,45)=   -0.23900170752781D-08
s(79,45)=    0.13299013133179D-07
s(80,45)=    0.43964054988728D-08
s(46,46)=   -0.88356225965378D-07
s(47,46)=   -0.64760424058702D-07
s(48,46)=   -0.39826277930573D-07
s(49,46)=   -0.75906935850517D-07
s(50,46)=    0.98672161286397D-08
s(51,46)=    0.20659809230718D-07
s(52,46)=   -0.58253254250572D-07
s(53,46)=   -0.39453247431162D-08
s(54,46)=   -0.15453009145983D-07
s(55,46)=   -0.10065623045061D-07
s(56,46)=    0.64678697152750D-07
s(57,46)=    0.69102225799176D-07
s(58,46)=    0.12211266501667D-07
s(59,46)=   -0.31701203022193D-07
s(60,46)=   -0.35969777153588D-07
s(61,46)=    0.11654300141702D-07
s(62,46)=   -0.41135319043059D-08
s(63,46)=   -0.74714099779694D-08
s(64,46)=   -0.17431396789403D-08
s(65,46)=   -0.32480237129338D-07
s(66,46)=   -0.22636362469650D-08
s(67,46)=    0.39752977211402D-07
s(68,46)=   -0.41201164173970D-08
s(69,46)=   -0.19071890092420D-08
s(70,46)=   -0.10675912037546D-07
s(71,46)=   -0.16867147025881D-07
s(72,46)=    0.67612347551718D-08
s(73,46)=   -0.90103490948792D-08
s(74,46)=    0.57580821456753D-08
s(75,46)=    0.22528547467380D-07
s(76,46)=    0.19236683341946D-07
s(77,46)=    0.32723809997015D-07
s(78,46)=    0.15739454322480D-07
s(79,46)=    0.20756800713409D-07
s(80,46)=    0.21001965327917D-07
s(47,47)=    0.14713449426786D-06
s(48,47)=   -0.14490673119231D-07
s(49,47)=    0.58241696161682D-07
s(50,47)=   -0.39377256236389D-07
s(51,47)=   -0.66624987073184D-07
s(52,47)=   -0.15390784582860D-07
s(53,47)=    0.23609855345167D-07
s(54,47)=    0.31331887295154D-07
s(55,47)=    0.53148542810258D-07
s(56,47)=    0.56560670984931D-08
s(57,47)=   -0.75416563810170D-07
s(58,47)=   -0.68686622162163D-07
s(59,47)=   -0.61226247817557D-08
s(60,47)=    0.42852299221086D-07
s(61,47)=    0.18568547550282D-07
s(62,47)=   -0.11500339104735D-07
s(63,47)=   -0.43090253077912D-10
s(64,47)=    0.80779319148837D-08
s(65,47)=    0.29102735381094D-07
s(66,47)=    0.73412240002891D-08
s(67,47)=   -0.14108928666898D-07
s(68,47)=   -0.49923372973142D-08
s(69,47)=   -0.19528244915968D-07
s(70,47)=   -0.18079336718715D-07
s(71,47)=   -0.17072216033032D-07
s(72,47)=   -0.19817637709924D-07
s(73,47)=   -0.80644078505365D-08
s(74,47)=    0.15869731761828D-08
s(75,47)=   -0.70347199851673D-08
s(76,47)=    0.26620112804348D-07
s(77,47)=    0.26024658824095D-07
s(78,47)=    0.30006950515015D-07
s(79,47)=    0.41591832853396D-07
s(80,47)=    0.14920507656653D-07
s(48,48)=   -0.85825814090879D-08
s(49,48)=    0.40424701944199D-07
s(50,48)=    0.27571321637424D-07
s(51,48)=    0.58168004029702D-07
s(52,48)=    0.39701694424485D-07
s(53,48)=   -0.47797012917057D-07
s(54,48)=    0.68366852961632D-08
s(55,48)=   -0.29732342377243D-07
s(56,48)=   -0.24298702145464D-07
s(57,48)=    0.15659575360552D-07
s(58,48)=   -0.10742059006886D-08
s(59,48)=    0.33194866195169D-07
s(60,48)=    0.34023828374198D-07
s(61,48)=    0.17927001976014D-08
s(62,48)=   -0.27907183126808D-07
s(63,48)=   -0.13600043604844D-08
s(64,48)=    0.33174845311097D-08
s(65,48)=    0.57892016321730D-08
s(66,48)=    0.42572312647152D-07
s(67,48)=    0.26056621586498D-07
s(68,48)=   -0.59818467585085D-08
s(69,48)=   -0.22532828133696D-07
s(70,48)=   -0.99543454356127D-08
s(71,48)=   -0.69725754731765D-08
s(72,48)=   -0.16967058726840D-07
s(73,48)=   -0.21269195383782D-07
s(74,48)=   -0.43940549399126D-08
s(75,48)=   -0.22513525598338D-07
s(76,48)=   -0.11404204479014D-07
s(77,48)=   -0.14328969595887D-07
s(78,48)=   -0.40363436125631D-07
s(79,48)=    0.36945342616493D-08
s(80,48)=   -0.44130987055631D-07
s(49,49)=   -0.11260526377236D-06
s(50,49)=   -0.10635221628594D-07
s(51,49)=   -0.46731148814368D-07
s(52,49)=   -0.47249561144027D-07
s(53,49)=    0.47589248753034D-07
s(54,49)=    0.37343882195670D-07
s(55,49)=   -0.33983604958279D-08
s(56,49)=    0.26655237671410D-07
s(57,49)=   -0.21709226762108D-07
s(58,49)=   -0.15462558454324D-07
s(59,49)=    0.11967295783703D-07
s(60,49)=    0.68721394585310D-07
s(61,49)=    0.20059662005922D-08
s(62,49)=   -0.58975984561642D-07
s(63,49)=   -0.16660034054480D-09
s(64,49)=    0.91325709578376D-08
s(65,49)=    0.11002498982300D-07
s(66,49)=   -0.22791133905585D-07
s(67,49)=    0.93887282819034D-08
s(68,49)=    0.66094156933658D-08
s(69,49)=   -0.27007489870657D-08
s(70,49)=    0.29936712231289D-07
s(71,49)=    0.17294441517683D-07
s(72,49)=   -0.18010949037563D-07
s(73,49)=    0.18182500570283D-07
s(74,49)=   -0.26102567975112D-08
s(75,49)=   -0.16396642101500D-07
s(76,49)=   -0.83725776177244D-08
s(77,49)=   -0.57245831551524D-08
s(78,49)=   -0.15982186385015D-07
s(79,49)=   -0.19515174299655D-07
s(80,49)=    0.99399969558361D-08
s(50,50)=    0.63310931071963D-07
s(51,50)=   -0.32426966034317D-07
s(52,50)=    0.27758360215090D-07
s(53,50)=   -0.73657968897853D-08
s(54,50)=   -0.83108567351799D-07
s(55,50)=    0.35777826327563D-07
s(56,50)=   -0.33090782260912D-07
s(57,50)=   -0.77057505632743D-09
s(58,50)=    0.32552459967129D-07
s(59,50)=    0.16367778009434D-07
s(60,50)=   -0.50405170622168D-07
s(61,50)=   -0.41804341485482D-07
s(62,50)=    0.56537902950171D-08
s(63,50)=   -0.31947010347959D-07
s(64,50)=   -0.40373010107808D-07
s(65,50)=    0.25315974135676D-07
s(66,50)=    0.36806906177224D-07
s(67,50)=   -0.21929380770970D-07
s(68,50)=   -0.29828341034869D-07
s(69,50)=   -0.48995486788422D-07
s(70,50)=   -0.34696480864670D-07
s(71,50)=   -0.42692532548816D-07
s(72,50)=    0.13317126321926D-07
s(73,50)=   -0.28270557313785D-08
s(74,50)=    0.13636799004785D-07
s(75,50)=    0.14139798253565D-07
s(76,50)=    0.23979956644393D-07
s(77,50)=    0.31689199109142D-08
s(78,50)=    0.38193030778588D-07
s(79,50)=    0.12136339665176D-07
s(80,50)=    0.45078308829265D-07
s(51,51)=    0.21036377950648D-09
s(52,51)=    0.71544713892764D-07
s(53,51)=    0.29826211749575D-07
s(54,51)=    0.55548242411779D-07
s(55,51)=    0.13136323014290D-07
s(56,51)=   -0.21254117753256D-07
s(57,51)=   -0.14524205674212D-07
s(58,51)=   -0.88440457924637D-08
s(59,51)=   -0.29418514317947D-07
s(60,51)=    0.93091893442430D-08
s(61,51)=    0.21319212441324D-07
s(62,51)=    0.23212297956107D-08
s(63,51)=   -0.15081366020778D-07
s(64,51)=    0.26453080769706D-07
s(65,51)=    0.45639925827969D-07
s(66,51)=   -0.54163027325770D-07
s(67,51)=   -0.99129065883353D-08
s(68,51)=    0.10400901583623D-07
s(69,51)=   -0.33509274100876D-08
s(70,51)=   -0.38675690405383D-07
s(71,51)=    0.75392364674283D-08
s(72,51)=    0.26504284412249D-08
s(73,51)=    0.18874294804081D-07
s(74,51)=    0.18106664634767D-07
s(75,51)=    0.35138122344854D-07
s(76,51)=    0.18575295478563D-07
s(77,51)=    0.15983899109469D-07
s(78,51)=    0.24893470078327D-07
s(79,51)=    0.33685374825783D-08
s(80,51)=    0.34671956334330D-07
s(52,52)=   -0.70651930540011D-07
s(53,52)=   -0.41950384236510D-07
s(54,52)=   -0.33194791172937D-07
s(55,52)=   -0.39405376929366D-07
s(56,52)=   -0.25413271840269D-07
s(57,52)=    0.42567249395159D-07
s(58,52)=    0.68249763828673D-07
s(59,52)=    0.39070291708969D-07
s(60,52)=   -0.12569698684176D-07
s(61,52)=   -0.45417003465732D-07
s(62,52)=    0.37861444823120D-07
s(63,52)=   -0.16508195158349D-07
s(64,52)=   -0.52838568901772D-07
s(65,52)=   -0.13031666041770D-07
s(66,52)=    0.88823811106077D-08
s(67,52)=    0.99335228524426D-08
s(68,52)=   -0.38245337922748D-07
s(69,52)=   -0.20037012000485D-07
s(70,52)=   -0.15952667327445D-07
s(71,52)=   -0.30878355546525D-07
s(72,52)=    0.49328447297288D-08
s(73,52)=   -0.25137219351314D-08
s(74,52)=   -0.26317381627957D-07
s(75,52)=   -0.10197093108605D-07
s(76,52)=   -0.28128754088150D-07
s(77,52)=    0.20995695694403D-07
s(78,52)=   -0.26561726441059D-07
s(79,52)=    0.11886220694607D-07
s(80,52)=    0.17753071244395D-07
s(53,53)=    0.97962180638267D-08
s(54,53)=   -0.45183423878581D-07
s(55,53)=   -0.48324882271015D-08
s(56,53)=   -0.34896744810624D-07
s(57,53)=   -0.23690015587009D-07
s(58,53)=    0.38747326759476D-07
s(59,53)=    0.13113418525450D-07
s(60,53)=    0.33962232603379D-08
s(61,53)=   -0.15142561871763D-08
s(62,53)=   -0.38150568764116D-08
s(63,53)=    0.32246974551838D-08
s(64,53)=    0.16024365607515D-07
s(65,53)=   -0.89546432656164D-08
s(66,53)=    0.10225322836894D-08
s(67,53)=   -0.17315652243658D-09
s(68,53)=   -0.27938667603040D-07
s(69,53)=   -0.12245001428012D-08
s(70,53)=    0.25485633978499D-08
s(71,53)=   -0.26784892522196D-07
s(72,53)=   -0.15376969443526D-07
s(73,53)=    0.11406429223802D-07
s(74,53)=    0.53529327504560D-08
s(75,53)=    0.10681455381125D-07
s(76,53)=    0.83390598556217D-08
s(77,53)=    0.19510573348420D-07
s(78,53)=    0.21111791659573D-07
s(79,53)=    0.27664945312980D-07
s(80,53)=    0.26483685256207D-07
s(54,54)=    0.47513077720641D-07
s(55,54)=    0.69657504196109D-07
s(56,54)=    0.65192698677921D-08
s(57,54)=    0.71003456137054D-07
s(58,54)=    0.23737682516615D-07
s(59,54)=   -0.51609046219115D-07
s(60,54)=   -0.34218740125708D-08
s(61,54)=    0.29132938123487D-07
s(62,54)=    0.57302906541253D-08
s(63,54)=   -0.12238570956267D-07
s(64,54)=   -0.19117351656598D-07
s(65,54)=   -0.12367589994689D-07
s(66,54)=   -0.33864393915753D-07
s(67,54)=    0.80448539923812D-08
s(68,54)=   -0.71081539754366D-09
s(69,54)=   -0.40972477236838D-07
s(70,54)=    0.10061847376392D-07
s(71,54)=    0.71702987132682D-08
s(72,54)=    0.15580346865279D-07
s(73,54)=    0.12406339655761D-07
s(74,54)=    0.22252425626581D-07
s(75,54)=    0.22043889058061D-07
s(76,54)=    0.27329153380725D-07
s(77,54)=    0.32570864320265D-07
s(78,54)=    0.32754415245410D-07
s(79,54)=    0.37315320924442D-07
s(80,54)=    0.21477785942637D-07
s(55,55)=   -0.54398801812038D-07
s(56,55)=   -0.66889995585568D-08
s(57,55)=   -0.21058372235604D-07
s(58,55)=    0.37081685359186D-08
s(59,55)=    0.71654115252569D-08
s(60,55)=   -0.17916046801929D-07
s(61,55)=    0.13295320170212D-08
s(62,55)=    0.21185601523102D-07
s(63,55)=    0.37634078455465D-08
s(64,55)=   -0.13822363141095D-07
s(65,55)=    0.47425549068940D-08
s(66,55)=   -0.11252877422222D-07
s(67,55)=   -0.19606476005546D-07
s(68,55)=   -0.11798201001951D-07
s(69,55)=    0.35865636994739D-08
s(70,55)=   -0.10121100056996D-07
s(71,55)=   -0.74298850862301D-08
s(72,55)=    0.13196373038602D-07
s(73,55)=    0.84763716554278D-08
s(74,55)=    0.16066201576666D-07
s(75,55)=    0.39494115243290D-07
s(76,55)=    0.97091141641524D-09
s(77,55)=    0.37145136854224D-07
s(78,55)=    0.18199323664616D-07
s(79,55)=    0.29703556083206D-07
s(80,55)=    0.26213441306470D-07
s(56,56)=    0.62867991403506D-07
s(57,56)=   -0.10499416242556D-06
s(58,56)=   -0.82067734075627D-09
s(59,56)=    0.23715074606582D-08
s(60,56)=   -0.29945665088491D-07
s(61,56)=    0.64656896720423D-07
s(62,56)=    0.27965100961512D-07
s(63,56)=   -0.49862906019927D-07
s(64,56)=   -0.68567957419056D-08
s(65,56)=    0.10491020009097D-07
s(66,56)=    0.83354258331501D-08
s(67,56)=    0.24779337805287D-07
s(68,56)=    0.15281431067374D-07
s(69,56)=   -0.10231918583867D-07
s(70,56)=   -0.41218471396085D-07
s(71,56)=    0.11789407325687D-07
s(72,56)=   -0.14359946183020D-07
s(73,56)=   -0.17711757980606D-07
s(74,56)=   -0.63020580377855D-08
s(75,56)=   -0.31543715348339D-07
s(76,56)=   -0.24275167790506D-07
s(77,56)=   -0.15141392099067D-07
s(78,56)=   -0.81947242971175D-08
s(79,56)=   -0.10894082730477D-07
s(80,56)=    0.12019962018137D-07
s(57,57)=    0.59397751117806D-07
s(58,57)=    0.59131586715798D-07
s(59,57)=   -0.23803196729421D-07
s(60,57)=   -0.85821146174170D-08
s(61,57)=    0.22187588995930D-07
s(62,57)=   -0.81106326561042D-09
s(63,57)=   -0.47795036028833D-07
s(64,57)=    0.29029196654898D-07
s(65,57)=    0.30446553284565D-08
s(66,57)=   -0.69944760388023D-08
s(67,57)=   -0.37134178308185D-07
s(68,57)=   -0.66037981229570D-08
s(69,57)=   -0.82272540097335D-09
s(70,57)=   -0.30806810457634D-07
s(71,57)=   -0.53033059572942D-08
s(72,57)=   -0.97638239605025D-08
s(73,57)=   -0.13427624432836D-08
s(74,57)=   -0.57168157396985D-08
s(75,57)=   -0.48107615558298D-08
s(76,57)=    0.10473602200356D-07
s(77,57)=   -0.17519450355493D-07
s(78,57)=    0.19487549535547D-07
s(79,57)=   -0.19116598372942D-08
s(80,57)=    0.17285407211804D-07
s(58,58)=   -0.66548671848599D-07
s(59,58)=    0.63395712510760D-07
s(60,58)=    0.52290763405760D-07
s(61,58)=   -0.16021621469904D-08
s(62,58)=    0.51950671773504D-08
s(63,58)=   -0.76078998983895D-07
s(64,58)=   -0.57436902457982D-07
s(65,58)=    0.17272421765174D-09
s(66,58)=   -0.77008978846624D-08
s(67,58)=    0.99101754211039D-08
s(68,58)=   -0.92143342790256D-08
s(69,58)=    0.58096749493959D-08
s(70,58)=    0.27298999535042D-08
s(71,58)=    0.11267532846781D-07
s(72,58)=    0.22042424604763D-07
s(73,58)=    0.15132827763681D-07
s(74,58)=    0.24075334526329D-09
s(75,58)=   -0.19292540048521D-07
s(76,58)=   -0.11267257099437D-07
s(77,58)=   -0.11102199313804D-08
s(78,58)=   -0.17385969304236D-07
s(79,58)=   -0.22510357135488D-07
s(80,58)=    0.12383489683129D-07
s(59,59)=    0.11992279747726D-07
s(60,59)=   -0.67555608702795D-07
s(61,59)=    0.35146698133389D-07
s(62,59)=    0.72605128314059D-08
s(63,59)=    0.15792666009731D-07
s(64,59)=    0.46736651429319D-07
s(65,59)=   -0.89443881475131D-08
s(66,59)=   -0.17133916471255D-07
s(67,59)=    0.16914786297516D-07
s(68,59)=   -0.32717129405300D-07
s(69,59)=   -0.78943031116755D-08
s(70,59)=   -0.27967469760782D-07
s(71,59)=    0.60310276155845D-08
s(72,59)=    0.14975836860358D-07
s(73,59)=   -0.11638995483247D-07
s(74,59)=   -0.55780076970927D-08
s(75,59)=    0.17748789793816D-07
s(76,59)=    0.14157337242762D-08
s(77,59)=   -0.25127423122844D-07
s(78,59)=    0.31117970326567D-07
s(79,59)=   -0.71570880714313D-07
s(80,59)=    0.54507766391812D-07
s(60,60)=    0.33806073951768D-07
s(61,60)=    0.10314740165104D-07
s(62,60)=   -0.37478064741543D-07
s(63,60)=   -0.34252775468114D-07
s(64,60)=    0.48642575011430D-08
s(65,60)=    0.77230576340196D-08
s(66,60)=    0.70295432867863D-08
s(67,60)=    0.55957228805018D-08
s(68,60)=   -0.20126527793621D-07
s(69,60)=    0.63914199971802D-08
s(70,60)=    0.34074235087432D-08
s(71,60)=    0.49656381631242D-08
s(72,60)=    0.28421693627136D-07
s(73,60)=    0.98889600397305D-08
s(74,60)=    0.17917471647633D-07
s(75,60)=    0.11360752022602D-07
s(76,60)=    0.96148498384901D-08
s(77,60)=    0.23931814635793D-08
s(78,60)=    0.79182961499720D-08
s(79,60)=   -0.13685622855363D-08
s(80,60)=    0.28073495024369D-07
s(61,61)=   -0.41407423035897D-07
s(62,61)=    0.55979029634436D-07
s(63,61)=    0.25829311090826D-07
s(64,61)=    0.10619620266283D-09
s(65,61)=   -0.17222850065129D-07
s(66,61)=   -0.66553701157538D-07
s(67,61)=    0.13121462219590D-07
s(68,61)=   -0.17207105471159D-07
s(69,61)=   -0.77984903301220D-08
s(70,61)=    0.60514850047506D-07
s(71,61)=    0.20293526490478D-07
s(72,61)=    0.14539943854966D-07
s(73,61)=    0.23092662969080D-09
s(74,61)=    0.30129310624357D-07
s(75,61)=    0.73560537701318D-08
s(76,61)=   -0.61566914554130D-08
s(77,61)=   -0.37151541335308D-07
s(78,61)=   -0.36254812177093D-07
s(79,61)=   -0.14577961521528D-07
s(80,61)=   -0.68865787626503D-07
s(62,62)=   -0.31215048975304D-07
s(63,62)=   -0.16628989592644D-07
s(64,62)=   -0.63131986269691D-09
s(65,62)=    0.24720488679536D-07
s(66,62)=    0.28429277587590D-07
s(67,62)=    0.33906855794909D-07
s(68,62)=    0.17720872432398D-07
s(69,62)=   -0.18494273870061D-07
s(70,62)=    0.12648422962005D-07
s(71,62)=    0.25810442308041D-07
s(72,62)=    0.16440552812412D-07
s(73,62)=    0.47745083076850D-07
s(74,62)=   -0.13048173614386D-07
s(75,62)=    0.96979252154748D-08
s(76,62)=    0.37706894171804D-07
s(77,62)=   -0.17925278771155D-07
s(78,62)=    0.36200148328791D-07
s(79,62)=   -0.15225954817015D-07
s(80,62)=    0.26744454027209D-07
s(63,63)=   -0.22558645110680D-08
s(64,63)=   -0.49151193940631D-08
s(65,63)=   -0.39402932663715D-07
s(66,63)=   -0.46036703667575D-07
s(67,63)=   -0.20977839121650D-07
s(68,63)=    0.16002110773293D-07
s(69,63)=   -0.12683397255976D-07
s(70,63)=    0.17657362300959D-08
s(71,63)=   -0.15143763158202D-07
s(72,63)=   -0.30101354549422D-07
s(73,63)=    0.10328233855826D-08
s(74,63)=   -0.11458505549840D-07
s(75,63)=   -0.60474192455122D-08
s(76,63)=   -0.16121389244904D-08
s(77,63)=   -0.53676492634261D-09
s(78,63)=    0.13707435993822D-07
s(79,63)=   -0.14249148783387D-07
s(80,63)=    0.49464024873434D-07
s(64,64)=   -0.87127007987508D-08
s(65,64)=    0.45115224024110D-08
s(66,64)=    0.37998742533734D-07
s(67,64)=    0.37513599916173D-07
s(68,64)=   -0.56026408905094D-08
s(69,64)=   -0.27550665136848D-07
s(70,64)=   -0.41030973914308D-07
s(71,64)=    0.61569852234094D-08
s(72,64)=    0.16205387444015D-07
s(73,64)=    0.68623221548350D-08
s(74,64)=   -0.19499253018975D-07
s(75,64)=    0.22039916029745D-07
s(76,64)=    0.12512336543455D-07
s(77,64)=    0.10273752680108D-08
s(78,64)=    0.42106957858973D-07
s(79,64)=   -0.53641252498077D-08
s(80,64)=    0.31671170597039D-07
s(65,65)=   -0.86154990439638D-08
s(66,65)=   -0.13856725505752D-07
s(67,65)=    0.61620647865604D-07
s(68,65)=    0.11049275161977D-07
s(69,65)=   -0.48761580928671D-08
s(70,65)=    0.29630825455903D-07
s(71,65)=    0.19096852773588D-08
s(72,65)=   -0.64028481933194D-08
s(73,65)=    0.26046180016634D-07
s(74,65)=    0.41200616073616D-08
s(75,65)=   -0.50220689549103D-08
s(76,65)=    0.40170494053698D-09
s(77,65)=   -0.13684208231148D-07
s(78,65)=    0.17683956031768D-08
s(79,65)=   -0.34641154618707D-07
s(80,65)=    0.45758113500030D-07
s(66,66)=    0.38291404049821D-07
s(67,66)=   -0.53098489196634D-07
s(68,66)=   -0.58400680212276D-07
s(69,66)=   -0.26152182582589D-07
s(70,66)=   -0.14540154170139D-07
s(71,66)=   -0.30976964233253D-07
s(72,66)=   -0.32977526416183D-08
s(73,66)=   -0.14845165931176D-07
s(74,66)=   -0.89271052657336D-08
s(75,66)=   -0.35245788463029D-07
s(76,66)=   -0.50100770484762D-07
s(77,66)=   -0.71556364263464D-08
s(78,66)=   -0.36334712997025D-07
s(79,66)=   -0.58207980725588D-07
s(80,66)=    0.89675259412812D-08
s(67,67)=    0.30399187730913D-07
s(68,67)=   -0.36782955347181D-07
s(69,67)=   -0.13394025889653D-07
s(70,67)=    0.27411801870785D-07
s(71,67)=   -0.30025299874283D-07
s(72,67)=   -0.79953714271451D-08
s(73,67)=    0.46236008509558D-07
s(74,67)=   -0.26169940270101D-07
s(75,67)=    0.20961793125776D-07
s(76,67)=   -0.68893632994121D-07
s(77,67)=   -0.27660409661774D-07
s(78,67)=   -0.73197636664678D-07
s(79,67)=   -0.10205753070715D-06
s(80,67)=    0.14902266951516D-07
s(68,68)=   -0.31480609254161D-07
s(69,68)=    0.13384874521102D-07
s(70,68)=    0.77155956538291D-07
s(71,68)=   -0.66866541422553D-08
s(72,68)=   -0.31147076800038D-07
s(73,68)=    0.13639437365475D-07
s(74,68)=   -0.43860203396215D-07
s(75,68)=   -0.49178247117405D-07
s(76,68)=   -0.26699372673611D-07
s(77,68)=   -0.67382415045807D-08
s(78,68)=    0.14412917377966D-08
s(79,68)=    0.48342457723864D-08
s(80,68)=    0.17390883601169D-07
s(69,69)=   -0.13798990953075D-07
s(70,69)=    0.39375140480278D-07
s(71,69)=   -0.34762758993926D-07
s(72,69)=   -0.27182018774016D-07
s(73,69)=   -0.38370797062333D-08
s(74,69)=   -0.36338066725225D-08
s(75,69)=   -0.15205391846379D-07
s(76,69)=    0.12763287796531D-07
s(77,69)=   -0.11792058947694D-08
s(78,69)=    0.33239509418844D-07
s(79,69)=   -0.10252760607483D-07
s(80,69)=    0.30007375851130D-07
s(70,70)=   -0.42747259376631D-08
s(71,70)=    0.58453055678280D-07
s(72,70)=   -0.35390991644398D-07
s(73,70)=    0.51806667480987D-08
s(74,70)=    0.88448286469802D-08
s(75,70)=   -0.11810715988298D-07
s(76,70)=    0.10938269758611D-07
s(77,70)=    0.23304378325923D-07
s(78,70)=   -0.19234121923198D-07
s(79,70)=    0.21629308852771D-07
s(80,70)=   -0.14769034041885D-07
s(71,71)=    0.19811987374229D-07
s(72,71)=   -0.32001567638459D-07
s(73,71)=    0.16854298985710D-07
s(74,71)=    0.10196641563199D-07
s(75,71)=    0.48703440603701D-08
s(76,71)=    0.11880687702064D-08
s(77,71)=   -0.30652418851886D-07
s(78,71)=   -0.18887828748998D-07
s(79,71)=   -0.30638177018873D-07
s(80,71)=   -0.45786303998232D-07
s(72,72)=    0.14619987101331D-07
s(73,72)=   -0.28385543652245D-08
s(74,72)=    0.20754818105081D-07
s(75,72)=    0.91782741587554D-09
s(76,72)=   -0.11348075853818D-07
s(77,72)=   -0.16743464096756D-07
s(78,72)=   -0.21718438590062D-07
s(79,72)=    0.19282431527263D-07
s(80,72)=   -0.18371376363970D-07
s(73,73)=    0.28162386123265D-08
s(74,73)=   -0.11088863947687D-07
s(75,73)=    0.44883608270606D-08
s(76,73)=   -0.10342902486516D-07
s(77,73)=   -0.20982837873543D-07
s(78,73)=   -0.44623528381893D-08
s(79,73)=   -0.29457122452921D-07
s(80,73)=    0.16518231847476D-07
s(74,74)=   -0.90679068924709D-07
s(75,74)=    0.11943656675702D-07
s(76,74)=   -0.12287476815057D-07
s(77,74)=   -0.86534126471312D-08
s(78,74)=   -0.36699382529634D-07
s(79,74)=    0.39699432381805D-07
s(80,74)=    0.90511319573331D-08
s(75,75)=    0.38068956857789D-08
s(76,75)=   -0.11022066868994D-07
s(77,75)=   -0.66351939341889D-08
s(78,75)=    0.84753383878959D-08
s(79,75)=    0.34003428729587D-07
s(80,75)=    0.45474580203248D-08
s(76,76)=   -0.47755842568499D-07
s(77,76)=   -0.31455773211009D-07
s(78,76)=   -0.12677631609515D-06
s(79,76)=    0.75135393078467D-07
s(80,76)=    0.29637278271697D-07
s(77,77)=   -0.28273200616638D-07
s(78,77)=   -0.30006709911530D-07
s(79,77)=   -0.46078138454365D-07
s(80,77)=    0.47083765379289D-07
s(78,78)=   -0.27390694608716D-07
s(79,78)=    0.10101212031673D-07
s(80,78)=   -0.77640066465978D-07
s(79,79)=    0.64488227584258D-07
s(80,79)=   -0.72740512544236D-08
s(80,80)=   -0.19642737190084D-07

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	r = 3396000.d0 !  MOLA potential surface
	xi = ae/r
	m=0
	x=0.
	sum=0.
	call lgndr(lmax,m,x,plm,root)
	do l=2, lmax
	  sum = sum + c(l,m) *xi**l *plm(l-m+1)
	enddo
	sum = sum+1.
! centrifugal potential	
	v0=(gm*sum/r + 0.5*omega**2 * r**2)
!	write(out,*)'v0=', v0
	return

end subroutine areoid_ini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	subroutine geoid(dlon,dlat,rg)
! marsgeoid.f
! calculate radii of surface of constant gravity potential on rotating body
! given spherical harmonic potential of body at lat,lon,elevation
!
! W = V +Phi = gm/r [1+ (R/r)**n [harmonics]] -1/2 omega**2 r**2

        implicit none

        double precision dlon,dlat
        double precision rg
        
        double precision pi,d2r
	parameter (pi=3.141592653589792D0, d2r=pi/180.d0)
! data structure of gravity field
        integer ndeg,ndeg1,ndeg2,nd2p3
	parameter (ndeg=90,ndeg1=ndeg+1,ndeg2=2*ndeg,nd2p3=ndeg2+3)
        double precision v0,omega,ae,gm
	double precision c(0:ndeg,0:ndeg),s(0:ndeg,0:ndeg)
	common /gmm1/v0,omega,ae,gm,c,s
        integer lmin,lmax
        double precision root(nd2p3)
        double precision r
	common /sqr/ lmin,lmax,root,r
	double precision plm(ndeg1)
	double precision tol
        data tol /0.125d0/  ! tolerence on computed value of rg

        integer i,m,l
        double precision rlon,rlat
        double precision x,cslt,xi,sum,cslm
        double precision diff
! save r
	rg = r

	rlon=dlon*d2r
	rlat=dlat*d2r
	x = sin(rlat)
	cslt= cos(rlat)

 	do i=1,8 ! usually 3 iterations suffice

	  xi = ae/r
	  sum = 0.
	  do m=0, lmax
	   call lgndr(lmax,m,x,plm,root)
	   do l=m, lmax
	     cslm =  c(l,m)*cos(m*rlon) + s(l,m)*sin(m*rlon)
	     sum = sum + cslm*xi**l *plm(l-m+1)
	   enddo
	  enddo
	  sum=sum+1.
! centrifugal potential	
	  rg=(gm*sum+0.5*(omega**2)*(r**3)*(cslt**2))/v0
!          write(out,*) i,r,rg
          diff = r-rg
          r = rg
          if( abs(diff).lt. tol) goto 400
	enddo ! do i=1,8

 400	continue

	end subroutine geoid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE LGNDR(LMAX,M,X,PLM,SQR)
! Return vector of plm's of degree from m to lmax, for order m
      implicit none
      integer LMAX
      integer M
      double precision X
      double precision PLM(*),SQR(*)
      integer IFACT,I,LL
      double precision OM, PMM, PMMP1, PLL, SOMX2
!  needs sqrt of integers from 1 to 2*L+3
!  SIGN MODIF. OF NUM.REC. ALGORITHM, TO MATCH GEOPHYS. CONVENTION.
!  AND NORMALIZATION TO MEAN SQUARE UNITY.

      PMMP1=0 !dummy initialization to get rid of compiler warnings

      OM = 1.
      PMM=SQR(2*M+1)
      IF(M.gt.0) then
        PMM=SQR(2)*PMM
        SOMX2=SQRT((1.-X)*(1.+X))
        IFACT=1
        do I=1,M
          OM=-OM
          PMM=-PMM*(SQR(IFACT)/SQR(IFACT+1))*SOMX2
          IFACT=IFACT+2
        enddo
      ENDIF
      PLM(1)=OM*PMM
      if(LMAX.gt.M) then
        PMMP1=X*SQR(2*M+3)*PMM
        PLM(2)=OM*PMMP1
      endif
      if (LMAX.gt.M+1) then
         DO LL=M+2,LMAX
	   PLL=SQR(2*LL+1)*SQR(2*LL-1)/(SQR(LL+M)*SQR(LL-M)) &
     & * (X*PMMP1-PMM*SQR(LL+M-1)*SQR(LL-M-1)/(SQR(2*LL-1)*SQR(2*LL-3)))
           PLM(LL-M+1)=OM*PLL
           PMM=PMMP1
           PMMP1=PLL
         ENDDO
      endif
      END SUBROUTINE LGNDR

