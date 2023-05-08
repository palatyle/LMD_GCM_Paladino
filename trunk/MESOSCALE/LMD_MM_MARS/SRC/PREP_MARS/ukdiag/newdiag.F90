! ====================================================================
!
! This routine reads in a diagfi.nc file output from the UK-LMD MGCM
! and converts it into a format suitable for producing initial and
! boundary conditions for the mesoscale model. Specifically, the
! arrays are interpolated onto latitudes ranging from +90 --> -90
! and longitudes ranging from -180 --> +180. The surface geopotential
! is also added if it doesn't already exist, and a check is made to
! see if the starting sol is present in the controle array.
!
! Arrays written on regolith layers are now interpolated back on to
! soil layers for use in the mesoscale model (sometimes different
! layers had to be specified for the regolith scheme (compared to
! the soil layers) to ensure numerical stability at high obliquities).
!
! The new file created is called diagfi_new.nc, which can be linked
! to the PREP_MARS directory as input_diagfi.nc
! 
! Liam Steele June 2014, March 2015
!
! ====================================================================

program newdiag
implicit none

include "netcdf.inc"

integer :: i,j,k,ierr,nid,nid2,nid3
integer :: ndims,nvars,ngatts,unlimdimid
integer :: lonid,latid,sigid,soilid,soildim,regid,timeid,contid,phisid
integer :: nlatnew,nlonnew,nsoil,nreg,npos
integer :: vtype,vatts,vdims,solstart
integer, dimension(4) :: dimids
integer, dimension(:), allocatable :: id,dimid,dimarr,loc,corners,edges
real, dimension(:), allocatable :: lat,lon,sig,soil,regol,time,controle
real, dimension(:), allocatable :: newlat,newlon,weight
real, dimension(:,:), allocatable :: geopot,geopot_new,array2,array2_new
real, dimension(:,:,:), allocatable :: array3,array3_new
real, dimension(:,:,:,:), allocatable :: array4,array4_new
real, dimension(:,:,:,:), allocatable :: new_array
character*20, dimension(:), allocatable :: dimname
character*20 :: vname
character*2 :: printid,res
logical :: wphis,skipv

! Open file for reading
ierr = nf_open("diagfi.nc",nf_nowrite,nid)
if (ierr /= nf_noerr) then
  print*, 'Oh dear: cannot find input_diagfi.nc'
  stop
endif

! Check properties of diagfi
ierr = nf_inq(nid, ndims, nvars, ngatts, unlimdimid)
print*, '--------------'
print*, 'File contains:'
write(printid,'(i2)'),ndims
print*, printid//' dimensions'
write(printid,'(i2)'),nvars
print*, printid//' variables'
print*, '--------------'
allocate(dimarr(ndims))
allocate(dimname(ndims))

! Get dimensions
do i = 1, ndims
  ierr = nf_inq_dim(nid,i,dimname(i),dimarr(i))
  if (ierr /= nf_noerr) then
    print*, 'Whoops: problem reading dimensions in diagfi.nc'
    stop
  endif
  if (trim(dimname(i)) == 'latitude' .or. trim(dimname(i)) == 'lat') then
    dimname(i) = 'latitude'
    allocate(lat(dimarr(i)))
    ierr = nf_get_var_real(nid,i,lat)
    nlatnew = dimarr(i)+1
  else if (trim(dimname(i)) == 'longitude' .or. trim(dimname(i)) == 'lon') then
    dimname(i) = 'longitude'
    allocate(lon(dimarr(i)))
    ierr = nf_get_var_real(nid,i,lon)
    nlonnew = dimarr(i)+1
  else if (trim(dimname(i)) == 'sigma') then
    allocate(sig(dimarr(i)))
    ierr = nf_get_var_real(nid,i,sig)
  else if (trim(dimname(i)) == 'soildepth' .or. trim(dimname(i)) == 'soil') then
    dimname(i) = 'soildepth'
    allocate(soil(dimarr(i)))
    ierr = nf_get_var_real(nid,i,soil)
    nsoil = dimarr(i)
    soildim = i
  else if (trim(dimname(i)) == 'regdepth') then
    dimname(i) = 'regdepth'
    allocate(regol(dimarr(i)))
    ierr = nf_get_var_real(nid,i,regol)
    nreg = dimarr(i)
  else if (trim(dimname(i)) == 'Time' .or. trim(dimname(i)) == 'time') then
    dimname(i) = 'Time'
    allocate(time(dimarr(i)))
    ierr = nf_get_var_real(nid,i,time)
  else if (trim(dimname(i)) == 'lentable') then
    allocate(controle(dimarr(i)))
    ierr = nf_get_var_real(nid,i,controle)
  endif
enddo

! Check for starting sol in diagfi
if (controle(4) == 0) then
  print*, 'No starting sol sumber in contole array. Please enter [1-669]:'
  read*, solstart
  do while (solstart > 669 .or. solstart < 0)
    if (solstart > 669) print*, 'Sol > 669. Have another go...'
    if (solstart < 0) print*, 'A negative sol is just silly. Try again...'
    read*, solstart
  enddo
  controle(4) = solstart
else if (controle(4) > 669) then
  controle(4) = mod(controle(4),669.)
endif

! Add soil mid-layer depths if not already present
if (soil(1) == 0 .or. soil(1) == 1) then
  print*, 'Adding the following soil levels (depths in m):'
  do i = 1, nsoil
    soil(i) = 2.e-4*(2.**(i-1.5))
    print*, i, soil(i)
  enddo
endif

! Make new lat/lon arrays
allocate(newlat(nlatnew))
allocate(newlon(nlonnew))
do i = 1, nlatnew
  newlat(i) = 90. + (i-1)*(lat(2)-lat(1))
enddo
do i = 1, nlonnew
  newlon(i) = -180. + (i-1)*(lon(2)-lon(1))
enddo

! Set weights and locations for latitudinal interpolation (assumes linear spacing between lats)
allocate(weight(nlatnew))
allocate(loc(nlatnew))
weight(:) = 0.5
weight(1) = 1.5
weight(nlatnew) = -0.5
do i = 1, nlatnew 
  loc(i) = i-1
enddo
loc(1) = 1
loc(nlatnew) = nlatnew-2

! Create new netcdf file to write data to (doesn't matter about writing descriptions and units)
ierr = nf_create('diagfi_new.nc',nf_clobber,nid2)

! Define dimensions
do i = 1, ndims
  write(printid,'(i2)'),i
  print*, printid//' writing '//trim(dimname(i))//' to file'
  if (trim(dimname(i)) == 'longitude') then
    ierr = nf_def_dim(nid2,'longitude',nlonnew,lonid)
    ierr = nf_def_var(nid2,'longitude',nf_float,1,lonid,lonid)
  else if (trim(dimname(i)) == 'latitude') then
    ierr = nf_def_dim(nid2,'latitude',nlatnew,latid)
    ierr = nf_def_var(nid2,'latitude',nf_float,1,latid,latid)
  else  if (trim(dimname(i)) == 'sigma') then
    ierr = nf_def_dim(nid2,'sigma',size(sig),sigid)
    ierr = nf_def_var(nid2,'sigma',nf_float,1,sigid,sigid)
  else if (trim(dimname(i)) == 'soildepth') then
    ierr = nf_def_dim(nid2,'soildepth',size(soil),soilid)
    ierr = nf_def_var(nid2,'soildepth',nf_float,1,soilid,soilid)
  else if (trim(dimname(i)) == 'regdepth') then
    ierr = nf_def_dim(nid2,'regdepth',size(regol),regid)
    ierr = nf_def_var(nid2,'regdepth',nf_float,1,regid,regid)
  else if (trim(dimname(i)) == 'Time') then
    ierr = nf_def_dim(nid2,'Time',size(time),timeid)
    ierr = nf_def_var(nid2,'Time',nf_float,1,timeid,timeid)
  else if (trim(dimname(i)) == 'lentable') then
    ierr = nf_def_dim(nid2,'lentable',size(controle),contid)
    ierr = nf_def_var(nid2,'controle',nf_float,1,contid,contid)
  endif
  if (ierr /= nf_noerr) then
    print*, "Yikes! Problem defining dimensions in new_diagfi.nc"
    stop
  endif
enddo

ierr = nf_enddef(nid2)

! Write dimension arrays
ierr = nf_put_var_real(nid2,lonid,newlon)
ierr = nf_put_var_real(nid2,latid,newlat)
ierr = nf_put_var_real(nid2,sigid,sig)
ierr = nf_put_var_real(nid2,soilid,soil)
ierr = nf_put_var_real(nid2,regid,regol)
ierr = nf_put_var_real(nid2,timeid,time)
ierr = nf_put_var_real(nid2,contid,controle)

! Read each variable in turn from original file, convert to new lat/lon and write to new file
npos = ndims+1
wphis = .true.
do i = ndims+1, nvars

  ! Read variable dimensions etc.
  ierr = nf_inq_var(nid, i, vname, vtype, vdims, dimids, vatts)
  write(printid,'(i2)'),i
  if (trim(vname) == 'phisinit') wphis = .false.
  allocate(id(vdims))
  allocate(corners(vdims))
  allocate(edges(vdims))
  allocate(dimid(vdims))
  
  ! Require lon/lat to be first two dimensions (everything should satisfy this condition, but you never know!)
  skipv = .true.
  if (vdims >= 2) then
    if (dimname(dimids(1)) == 'longitude' .and. dimname(dimids(2)) == 'latitude')  skipv = .false. 
    if (dimname(dimids(1)) == 'latitude'  .and. dimname(dimids(2)) == 'longitude') skipv = .false. 
    if (skipv) then
      print*, " > "//trim(vname)//": lat and lon aren't first two dimensions. Skipping"
      cycle
    endif
  endif
  
  ! Change certain names to match LMD (more might need added later)
  if (trim(vname) == 'h2ommr')    vname = 'q02'
  if (trim(vname) == 'icemmr')    vname = 'q01'
  if (trim(vname) == 'h2o_ice_s') vname = 'qsurf01'
  
  print*, printid//' writing '//trim(vname)//' to file'
  
  ! Define arrays for specifying location in diagfi
  do k = 1, vdims
    id(k) = dimarr(dimids(k))
    dimid(k) = dimids(k)
    corners(k) = 1
    edges(k) = dimarr(dimids(k))
  enddo
  edges(lonid) = nlonnew
  edges(latid) = nlatnew
  edges(vdims) = 1
  
  ! Allocate and read old arrays
  if (vdims == 2) then
    allocate(array2(id(1),id(2)))
    ierr = nf_get_var_real(nid,i,array2)
  else if (vdims == 3) then
    allocate(array3(id(1),id(2),id(3)))
    ierr = nf_get_var_real(nid,i,array3)
  else if (vdims == 4) then
    allocate(array4(id(1),id(2),id(3),id(4)))
    ierr = nf_get_var_real(nid,i,array4) 
  endif
  
  ! Check if regolith depth is the vertical dimension (if so, will convert fields to soil depths)
  do k = 1, vdims
    if (dimname(dimids(k)).eq.'regdepth') then
    
      if (vdims == 3) then
        allocate(new_array(id(1),id(2),dimarr(soildim),1))
        call interpol_soil(soil,regol,nsoil,nreg,(/id(1:vdims),1/),array3,new_array)
        deallocate(array3)
        allocate(array3(id(1),id(2),dimarr(soildim)))
        array3 = new_array(1:id(1),1:id(2),1:dimarr(soildim),1)
        deallocate(new_array)
      else if (vdims == 4) then
        allocate(new_array(id(1),id(2),dimarr(soildim),id(4)))
        call interpol_soil(soil,regol,nsoil,nreg,id(1:vdims),array4,new_array)
        deallocate(array4)
        allocate(array4(id(1),id(2),dimarr(soildim),id(4)))
        array4 = new_array
        deallocate(new_array)
      endif   
      
      dimids(k) = soildim
      dimid(k) = soildim
      id(k) = dimarr(dimids(k))
      edges(k) = dimarr(dimids(k))

      exit
    endif
  enddo
  
  ! Allocate new arrays
  id(lonid) = nlonnew
  id(latid) = nlatnew
  if (vdims == 2) allocate(array2_new(id(1),id(2)))
  if (vdims == 3) allocate(array3_new(id(1),id(2),id(3)))
  if (vdims == 4) allocate(array4_new(id(1),id(2),id(3),id(4)))
  
  ! Transform onto new lat-lon grid (loops assume array dimensions ordered lon, lat, sig, time)
  if (vdims == 2) then
    do k = 1, nlatnew
      array2_new(1:nlonnew-1,k) = weight(k)*array2(:,loc(k)) + (1-weight(k))*array2(:,loc(k)+1)
    enddo
    array2_new(nlonnew,:) = array2_new(1,:)
    array2_new(:,1) = sum(array2_new(:,2))/nlonnew
    array2_new(:,nlatnew) = sum(array2_new(:,nlatnew-1))/nlonnew
    if (minval(array2) > 0.0) then
      where (array2_new < 0.0) array2_new = 0.0
    endif
  else if (vdims == 3) then
    do k = 1, nlatnew
      array3_new(1:nlonnew-1,k,:) = weight(k)*array3(:,loc(k),:) + (1-weight(k))*array3(:,loc(k)+1,:)
    enddo
    array3_new(nlonnew,:,:) = array3_new(1,:,:)
    do k = 1, id(3)
      array3_new(:,1,k) = sum(array3_new(:,2,k))/nlonnew
      array3_new(:,nlatnew,k) = sum(array3_new(:,nlatnew-1,k))/nlonnew
    enddo
    if (minval(array3) > 0.0) then
      where (array3_new < 0.0) array3_new = 0.0
    endif
  else if (vdims == 4) then
    do k = 1, nlatnew
      array4_new(1:nlonnew-1,k,:,:) = weight(k)*array4(:,loc(k),:,:) + (1-weight(k))*array4(:,loc(k)+1,:,:)
    enddo
    array4_new(nlonnew,:,:,:) = array4_new(1,:,:,:)
    do k = 1, id(4)
      do j = 1, id(3)
        array4_new(:,1,j,k) = sum(array4_new(:,2,j,k))/nlonnew
        array4_new(:,nlatnew,j,k) = sum(array4_new(:,nlatnew-1,j,k))/nlonnew
      enddo
    enddo
    if (minval(array4) > 0.0) then
      where (array4_new < 0.0) array4_new = 0.0
    endif
  else
    print*, ' > '//trim(vname)//': code does not interpolate 1D variables'
  endif
  
  ! Write variable to new file
  ierr = nf_redef(nid2)
  ierr = nf_def_var(nid2,trim(vname),nf_float,vdims,dimid,npos)
  
  ierr = nf_enddef(nid2)
  if (vdims == 2) then
    ierr = nf_put_var_real(nid2,npos,array2_new)
  else
    do k = 1, id(vdims)
      corners(vdims) = k
      if (vdims == 3) ierr = nf_put_vara_real(nid2,npos,corners,edges,array3_new(:,:,k))
      if (vdims == 4) ierr = nf_put_vara_real(nid2,npos,corners,edges,array4_new(:,:,:,k))
    enddo
    if (vdims == 3) deallocate(array3,array3_new)
    if (vdims == 4) deallocate(array4,array4_new)
  endif
  
  if (ierr /= nf_noerr) then
    print*, "Oh no! Problem writing "//trim(vname)
    stop
  endif
  
  deallocate(id,corners,edges,dimid)
  npos = npos + 1
  
enddo

! Add geopotential height to file
if (wphis) then
  print*, 'Adding phisinit to file'
  allocate(id(2))

  if (nlatnew-1 == 8)   res = '5'
  if (nlatnew-1 == 14)  res = '10'
  if (nlatnew-1 == 24)  res = '21'
  if (nlatnew-1 == 36)  res = '31'
  if (nlatnew-1 == 48)  res = '42'
  if (nlatnew-1 == 72)  res = '63'
  if (nlatnew-1 == 96)  res = '85'
  if (nlatnew-1 == 144) res = '127'
  if (nlatnew-1 == 192) res = '170'
  allocate(geopot(nlonnew-1,nlatnew-1))
  allocate(geopot_new(nlonnew,nlatnew))
  ierr = nf_open("phisinit_files/phisinit-T"//trim(res)//".nc",nf_nowrite,nid3)
  ierr = nf_inq_varid(nid3,"phis",phisid)
  if (ierr /= nf_noerr) then
    print*, "Oh no! Can't find geopotential file. Please link from the ukv_phisinit folder:"
    print*, "ln -sf $MMM/your_comp_dir/PREP_MARS/ukv_phisinit/phisinit-T"//trim(res)//".nc ."
    stop
  endif
  ierr = nf_get_var_real(nid3,phisid,geopot)
  ierr = nf_close(nid3)
  
  do k = 1, nlatnew
    geopot_new(1:nlonnew-1,k) = weight(k)*geopot(:,loc(k)) + (1-weight(k))*geopot(:,loc(k)+1)
  enddo
  geopot_new(nlonnew,:) = geopot_new(1,:)
  geopot_new(:,1) = sum(geopot_new(:,2))/nlonnew
  geopot_new(:,nlatnew) = sum(geopot_new(:,nlatnew-1))/nlonnew
  
  id(1) = lonid
  id(2) = latid
  ierr = nf_redef(nid2)
  ierr = nf_def_var(nid2,"phisinit",nf_float,2,id,npos)
  ierr = nf_enddef(nid2)
  ierr = nf_put_var_real(nid2,npos,geopot_new)
  if (ierr /= nf_noerr) then
    print*, "Oh no! Problem writing phisinit"
    stop
  endif
endif

print*, 'Closing files...'
ierr = nf_close(nid)
ierr = nf_close(nid2)

end
