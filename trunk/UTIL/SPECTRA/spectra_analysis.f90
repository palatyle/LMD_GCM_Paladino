!M. Indurain & E. Millour
!March 2015

PROGRAM spectra_analysis

!Compute energy spectrum of a 2d wind field with spherepack routines.
!v = sum(n=1..nlat-1) 0.5*b(0,n)*B(0,n)+0.5*c(0,n)*C(0,n) + sum(m=1..mmax-1)sum(n=m..nlat-1) b(m,n)*B(m,n)+c(m,n)*C(m,n)
!where B(m,n) and C(m,n) are vector spherical harmonics
!Then energy spectra is defined as:
!E(n) = 0.5*|b(0,n)|^2 + 0.5*|c(0,n)|^2 + sum(m=1..n) |b(m,n)|^2 + |c(m,n)|^2

!Wind field read from netcdf file follows lmdz convention:
!  u(nlong+1,nlat) zonal eastward wind, v(nlong+1,nlat) latitudinal wind
!Wind field needed for spherepack routines must be like:
!  u(nlat,nlong) zonal eastward wind, v(nlong,nlat) colatitudinal wind (from North pole)

use netcdf

implicit none

!Physical grid
integer :: nlat,nlong,nlongp1,nalt,ntime
double precision, dimension(:,:), allocatable :: merid_wind,zonal_wind !spherepack convention
double precision, dimension(:,:), allocatable :: merid_wind_lmdz,zonal_wind_lmdz !lmdz convention
!Spectral grid
integer :: mmax !spectra truncation mmax=min(nalt,nlong/2) or min(nalt,(nlong+1)/2)
double precision, dimension (:,:), allocatable :: br,bi,cr,ci !spectra coefficients
double precision, dimension (:,:), allocatable :: norm,norm_b,norm_c !norm of spectra coefficients
double precision, dimension (:), allocatable :: spectra,spectra_b,spectra_c
double precision, dimension (:,:), allocatable :: spectra_config,spectra_b_config,spectra_c_config
double precision, dimension (:,:), allocatable :: spectra_global,spectra_b_global,spectra_c_global
!Spherepack declarations
integer :: ierror,lvhaec,ldwork,lwork,l1,l2,nt=1,ityp=0,idvw,jdvw,mdab,ndab
double precision, dimension (:), allocatable :: wvhaec,work
double precision, dimension (:), allocatable :: dwork
!Netcdf stuff
integer :: idfile,idmerid_wind,idzonal_wind
integer :: idlat,idlong,idalt,idtime
!Input arguments
character (len=100), dimension(:), allocatable :: arg
character (len=100) :: file_netcdf,file_spectra
integer :: n_indt,n_indz,meant,meanz
integer, dimension(:), allocatable :: tab_indt,tab_indz
character (len=100) :: name_u,name_v,name_lat,name_long,name_alt,name_time
character (len=100) :: tab_name_u(4)=(/'u','U','vitu','zonal_wind'/)
character (len=100) :: tab_name_v(4)=(/'v','V','vitv','merid_wind'/)
character (len=100) :: tab_name_lat(2)=(/'latitude','lat'/)
character (len=100) :: tab_name_long(2)=(/'longitude','lon'/)
character (len=100) :: tab_name_alt(3)=(/'altitude','alt','presnivs'/)
character (len=100) :: tab_name_time(3)=(/'time','Time','time_counter'/)
logical :: is_div_rot
integer :: narg
!Other
integer :: i,j,t,z,mt,mz
integer :: number_config	!number of (t,z) configurations wanted by user (from 1 to n_indt*n_indz))
integer :: number_local		!number of (t,z) sub-configurations to mean (from 1 to (meant+1)*(meanz+1))
integer :: number_global	!number of total (t,z) configurations to compute (from 1 to n_indt*n_indz*(meant+1)*(meanz+1))
logical :: is_file,first_reading=.true.
integer, dimension(:), allocatable :: tmp_tab_int
character (len=100) :: tmp_char

!**********
!Initialisation
file_netcdf = '??@@??'
file_spectra = 'spectra'
n_indt = 0
n_indz = 0
allocate(tab_indt(20))
allocate(tab_indz(20))
allocate(tmp_tab_int(20))
tab_indt(:) = 1
tab_indz(:) = 1
meant = 0
meanz = 0
name_u = '??@@??'
name_v = '??@@??'
name_lat = '??@@??'
name_long = '??@@??'
name_alt = '??@@??'
name_time = '??@@??'
is_div_rot=.false.

!**********
!Input reading
narg = command_argument_count()
allocate(arg(narg))
do i = 1, narg
  call get_command_argument(i,arg(i))
end do
i = 1
do while (i .le. narg)
  if (arg(i) == '-o' .or. arg(i) == '-t' .or. arg(i) == '-z' .or. arg(i) == '-mt' .or. arg(i) == '-mz'&
  		.or. arg(i) == '-u' .or. arg(i) == '-v' .or. arg(i) == '-lat' .or. arg(i) == '-lon' &
		.or. arg(i) == '-alt' .or. arg(i) == '-time') then
    select case(arg(i))
    case('-t')
      n_indt = n_indt + 1
      read(arg(i+1),'(i10)' ) tab_indt(n_indt) !temporal indice			
    case('-z')
      n_indz = n_indz + 1
      read(arg(i+1),'(i10)' ) tab_indz(n_indz) !vertical indice
    case('-mt')
      read(arg(i+1),'(i10)' ) meant !extent of temporal mean
    case('-mz')
      read(arg(i+1),'(i10)' ) meanz !extent of vertical mean
    case('-o')
      file_spectra=arg(i+1) !spectra file
    case('-u')
      read(arg(i+1),'(a100)' ) name_u !name of zonal wind
    case('-v')
      read(arg(i+1),'(a100)' ) name_v !name of meridional wind
    case('-lat')
      read(arg(i+1),'(a100)' ) name_lat !name of latitude axis
    case('-lon')
      read(arg(i+1),'(a100)' ) name_long !name of longitude axis
    case('-alt')
      read(arg(i+1),'(a100)' ) name_alt !name of altitude axis
    case('-time')
      read(arg(i+1),'(a100)' ) name_time !name of time axis
    end select
    i = i + 2
    elseif (arg(i) == '-divrot') then
      is_div_rot = .true.
      i = i + 1
    elseif (arg(i) == '-h' .or. arg(i) == '--help') then
      print*,'Usage'
      print*,'spectra_analysis netcdfFile [option]'
      print*,'[-h or --help]	: brief help'
      print*,'[-t int]	: temporal indice (default: 1)'
      print*,'[-z int]	: vertical indice (default: 1)'
      print*,'[-mt int]	: extent of temporal mean (default: 0)'
      print*,'[-mz int]	: extent of vertical mean (default: 0)'
      print*,'[-o str]	: output spectra file name (default: spectra)'
      print*,'[-u str]	: name of zonal wind in netcdf file'
      print*,'[-v str]	: name of meridional wind in netcdf file'
      print*,'[-lat str]	: name of latitude field in netcdf file'
      print*,'[-lon str]	: name of longitude field in netcdf file'
      print*,'[-alt str]	: name of altitude field in netcdf file. ''none'' if no altitude axis.'
      print*,'[-time str]	: name of time field in netcdf file. ''none'' if no time axis.'
      print*,'[-divrot]	: total, divergence and vorticity spectra coefficient in output file'
      print*,''
      print*,'Compute the kinetic energy spectrum of a wind field on a longitude-latitude grid.'
      print*,'The grid must have a redundant point in longitude.'
      print*,'Ideally the analysis works better on a symetric grid (2N points in longitude and N points in latitude).'
      print*,'The output text file (called spectra by default) give the decomposition'
      print*,'of the velocity on the vectorial spherical harmonic basis.'
      stop 'End help'
    else
      file_netcdf = arg(i)
      i = i + 1
  end if
end do

!**********
!Check input/output files
inquire(file=trim(file_netcdf),exist=is_file)
if (.not. is_file) then
  print*,"no netcdf file: ",trim(file_netcdf)
  stop "Stopped"
end if
call system('rm -f '//trim(file_spectra))

!**********
!Netcdf file dimensions
ierror=nf90_open(trim(file_netcdf),NF90_NOWRITE,idfile)

ierror=-99999
i=0
ierror=nf90_inq_dimid(idfile,trim(name_lat),idlat) !try user named latitude
do while (ierror /= 0 .and. i <= size(tab_name_lat)-1)!try automatic named latitude
  i=i+1
  ierror=nf90_inq_dimid(idfile,trim(tab_name_lat(i)),idlat)
end do
if (ierror == 0) then
  ierror=nf90_inquire_dimension(idfile,idlat,tmp_char,nlat)
else
  print*,'latitude dimension not found!'
  stop ' must have this... use ''-lat latitudefieldname'' option.'
end if

ierror=-99999
i=0
ierror=nf90_inq_dimid(idfile,trim(name_long),idlong) !try user named longitude
do while (ierror /= 0 .and. i <= size(tab_name_long)-1)  !try automatic named longitude
  i=i+1
  ierror=nf90_inq_dimid(idfile,trim(tab_name_long(i)),idlong)
end do
if (ierror == 0) then
  ierror=nf90_inquire_dimension(idfile,idlong,tmp_char,nlongp1)
else
  print*,'longitude dimension not found!'
  stop ' must have this... use ''-lon longitudefieldname'' option.'
end if

ierror=-99999
i=0
ierror=nf90_inq_dimid(idfile,trim(name_alt),idalt) !try user named altitude
do while (ierror /= 0 .and. i <= size(tab_name_alt)-1)!try automatic named altitude
  i=i+1
  ierror=nf90_inq_dimid(idfile,trim(tab_name_alt(i)),idalt)
end do
if (ierror == 0) then
  ierror=nf90_inquire_dimension(idfile,idalt,tmp_char,nalt)
else
  if (name_alt == 'none') then
    print*,'netcdf file without altitude axis: all the same.'
  else 
    print*,'altitude dimension not found!'
    stop ' if no altitude axis, use ''-alt none'' option...'
  end if
end if

ierror=-99999
i=0
ierror=nf90_inq_dimid(idfile,trim(name_time),idtime) !try user named time
do while (ierror /= 0 .and. i <= size(tab_name_time)-1)  !try automatic named time
  i=i+1
  ierror=nf90_inq_dimid(idfile,trim(tab_name_time(i)),idtime)
end do
if (ierror == 0) then
  ierror=nf90_inquire_dimension(idfile,idtime,tmp_char,ntime)
else
  if (name_time == 'none') then
    print*,'netcdf file without time axis: all the same.'
  else 
    print*,'time dimension not found!'
    stop ' if no time axis, use ''-time none'' option...'
  end if
end if

!**********
!Check input indices and add average indices
if (n_indt == 0) n_indt = 1
tmp_tab_int(:) = tab_indt(:)
deallocate(tab_indt)
allocate(tab_indt(n_indt*(meant+1)))
do t = 1,n_indt
do mt = 0,meant
  tab_indt((t-1)*(meant+1)+1+mt) = tmp_tab_int(t)+mt
end do
end do  
do t = 1,n_indt*(meant+1)
  if (tab_indt(t) > ntime .or. tab_indt(t) < 1) tab_indt(t) = 1
end do

if (n_indz == 0) n_indz = 1
tmp_tab_int(:) = tab_indz(:)
deallocate(tab_indz)
allocate(tab_indz(n_indz*(meanz+1)))
do z = 1,n_indz
do mz = 0,meanz
  tab_indz((z-1)*(meanz+1)+1+mz) = tmp_tab_int(z)+mz
end do
end do  
do z = 1,n_indz*(meanz+1)
  if (tab_indz(z) > nalt .or. tab_indz(z) < 1) tab_indz(z) = 1
end do

!**********
!Loop over time and altitude indices
allocate(spectra_config(nlat,n_indt*n_indz))
allocate(spectra_b_config(nlat,n_indt*n_indz))
allocate(spectra_c_config(nlat,n_indt*n_indz))
allocate(spectra_global(nlat,n_indt*(meant+1)*n_indz*(meanz+1)))
allocate(spectra_b_global(nlat,n_indt*(meant+1)*n_indz*(meanz+1)))
allocate(spectra_c_global(nlat,n_indt*(meant+1)*n_indz*(meanz+1)))
do t = 1,n_indt
do z = 1,n_indz
number_config = (t-1)*n_indz+z
do mt = 0,meant
do mz = 0,meanz
number_local = mt*(meanz+1)+mz+1
number_global = (number_config-1)*(meant+1)*(meanz+1) + number_local
if (first_reading) then
  print*,'First netcdf file reading...'
  first_reading = .false.
end if

!**********
!Netcdf file reading
ierror=-99999
i=0
ierror=nf90_inq_varid(idfile,trim(name_u),idzonal_wind) !try user named zonal wind
do while (ierror /= 0 .and. i <= size(tab_name_u)-1)        !try automatic named zonal wind
  i=i+1
  ierror=nf90_inq_varid(idfile,trim(tab_name_u(i)),idzonal_wind)
end do
if (ierror == 0) then
  allocate(zonal_wind_lmdz(nlongp1,nlat))
  if (name_alt == 'none') then
    if(name_time == 'none') then
      ierror=nf90_get_var(idfile,idzonal_wind,zonal_wind_lmdz,(/1,1/),(/nlongp1,nlat/))
    else
      ierror=nf90_get_var(idfile,idzonal_wind,zonal_wind_lmdz,(/1,1,tab_indt((t-1)*(meant+1)+mt+1)/),(/nlongp1,nlat,1/))
    end if
  end if
  if (name_alt /= 'none') then
    if(name_time == 'none') then
      ierror=nf90_get_var(idfile,idzonal_wind,zonal_wind_lmdz,(/1,1,tab_indz((z-1)*(meanz+1)+1)/),(/nlongp1,nlat,1/))
    else
      ierror=nf90_get_var(idfile,idzonal_wind,zonal_wind_lmdz,&
      				(/1,1,tab_indz((z-1)*(meanz+1)+1),tab_indt((t-1)*(meant+1)+mt+1)/),(/nlongp1,nlat,1,1/))
    end if
  end if
else
  print*,'zonal wind variable not found!'
  stop ' must have this... use ''-u zonalwindfieldname'' option.'
end if

ierror=-99999
i=0
ierror=nf90_inq_varid(idfile,trim(name_v),idmerid_wind) !try user named meridional wind
do while (ierror /= 0 .and. i <= size(tab_name_v)-1)        !try automatic named meridional wind
  i=i+1
  ierror=nf90_inq_varid(idfile,trim(tab_name_v(i)),idmerid_wind)
end do
if (ierror == 0) then
  allocate(merid_wind_lmdz(nlongp1,nlat))
  if (name_alt == 'none') then
    if(name_time == 'none') then
      ierror=nf90_get_var(idfile,idmerid_wind,merid_wind_lmdz,(/1,1/),(/nlongp1,nlat/))
    else
      ierror=nf90_get_var(idfile,idmerid_wind,merid_wind_lmdz,(/1,1,tab_indt((t-1)*(meant+1)+mt+1)/),(/nlongp1,nlat,1/))
    end if
  end if
  if (name_alt /= 'none') then
    if(name_time == 'none') then
      ierror=nf90_get_var(idfile,idmerid_wind,merid_wind_lmdz,(/1,1,tab_indz((z-1)*(meanz+1)+1)/),(/nlongp1,nlat,1/))
    else
      ierror=nf90_get_var(idfile,idmerid_wind,merid_wind_lmdz,&
      				(/1,1,tab_indz((z-1)*(meanz+1)+1),tab_indt((t-1)*(meant+1)+mt+1)/),(/nlongp1,nlat,1,1/))
    end if
  end if
else
  print*,'meridional wind variable not found!'
  stop ' must have this... use ''-v meridionalwindfieldname'' option.'
end if

!**********From lmdz to spherepack convention
nlong=nlongp1-1
allocate(zonal_wind(nlat,nlong))
allocate(merid_wind(nlat,nlong))
zonal_wind(:,:)=transpose(zonal_wind_lmdz(1:nlong,:))
merid_wind(:,:)=transpose(merid_wind_lmdz(1:nlong,:))
!if (mod(nlat,2) == 0) then
!  merid_wind(1:nlat/2,:)=-merid_wind(1:nlat/2,:)
!else
!  merid_wind(1:(nlat+1)/2,:)=-merid_wind(1:(nlat+1)/2,:)
!end if

!#####
!Spectra computation
!#####
!**********Maximum value for m
if (mod(nlong,2) == 0) then
  mmax = min(nlat,nlong/2)
else
  mmax = min(nlat,(nlong+1)/2)
end if

!**********Vhaeci function (initialisations for Vhaec function)
if (mod(nlong,2) == 0) then
  l1 = min(nlat,nlong/2)
else
  l1 = min(nlat,(nlong+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lvhaec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlong+15
allocate(wvhaec(lvhaec))
wvhaec(:) = 0.
ldwork=2*(nlat+2)
allocate(dwork(ldwork))
dwork(:) = 0.
ierror=3
call vhaeci(nlat,nlong,wvhaec,lvhaec,dwork,ldwork,ierror)
 
!**********Vhaeci function result
select case (ierror)
  case(1) 
    print*,'Vhaeci: ERROR on nlat'
  case(2) 
    print*,'Vhaeci: ERROR on nlong'
  case(3) 
    print*,'Vhaeci: ERROR on lvhaec'
  case(4) 
    print*,'Vhaeci: ERROR on ldwork'
end select

!**********Vhaec function (computes spectra coefficients)
idvw=nlat
jdvw=nlong
mdab=mmax
ndab=nlat
allocate(br(mdab,ndab))
allocate(bi(mdab,ndab))
allocate(cr(mdab,ndab))
allocate(ci(mdab,ndab))
br(:,:) = 0.
bi(:,:) = 0.
cr(:,:) = 0.
ci(:,:) = 0.
if (mod(nlong,2) == 0) then
  l1 = min(nlat,nlong/2)
else
  l1 = min(nlat,(nlong+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lvhaec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlong+15
if (ityp .le. 2) then
  lwork=nlat*(2*nt*nlong+max(6*l2,nlong))
else
  lwork=l2*(2*nt*nlong+max(6*nlat,nlong))
end if
allocate(work(lwork))
work(:) = 0.
ierror=3
call vhaec(nlat,nlong,ityp,nt,merid_wind,zonal_wind,idvw,jdvw,br,bi,cr,ci,mdab,ndab,wvhaec,lvhaec,work,lwork,ierror)

!**********Vhaec function result
select case (ierror)
  case(1) 
    print*,'Vhaec: ERROR on nlat'
  case(2) 
    print*,'Vhaec: ERROR on nlong'
  case(3) 
    print*,'Vhaec: ERROR on ityp'
  case(4) 
    print*,'Vhaec: ERROR on nt'
  case(5) 
    print*,'Vhaec: ERROR on idvw'
  case(6) 
    print*,'Vhaec: ERROR on jdvw'
  case(7) 
    print*,'Vhaec: ERROR on mdab'
  case(8) 
    print*,'Vhaec: ERROR on ndab'
  case(9) 
    print*,'Vhaec: ERROR on lvhaec'
  case(10) 
    print*,'Vhaec: ERROR on lwork'
end select

!**********Energy spectra
allocate(norm(mMax,nLat))
allocate(norm_b(mMax,nLat))
allocate(norm_c(mMax,nLat))
do i = 1,mMax
do j = 1,nLat
  norm(i,j)=(br(i,j)*br(i,j)+bi(i,j)*bi(i,j)+cr(i,j)*cr(i,j)+ci(i,j)*ci(i,j))
  norm_b(i,j)=(br(i,j)*br(i,j)+bi(i,j)*bi(i,j))
  norm_c(i,j)=(cr(i,j)*cr(i,j)+ci(i,j)*ci(i,j))
end do
end do
allocate(spectra(nLat))
allocate(spectra_b(nLat))
allocate(spectra_c(nLat))
do j = 1,nLat
  spectra(j) = 0.5*norm(1,j)
  spectra_b(j) = 0.5*norm_b(1,j)
  spectra_c(j) = 0.5*norm_c(1,j)
  do i = 2,mMax
    spectra(j) = spectra(j) + norm(i,j)
    spectra_b(j) = spectra_b(j) + norm_b(i,j)
    spectra_c(j) = spectra_c(j) + norm_c(i,j)
  end do
end do
!do j = 2,nLat
!  spectra(j) = spectra(j)/((j-1)*j)
!  spectra_b(j) = spectra_b(j)/((j-1)*j)
!  spectra_c(j) = spectra_c(j)/((j-1)*j)
!end do
write(*,'(a16,i5,a3,i5,a1)') 'done spectra (t=',tab_indt((t-1)*(meant+1)+mt+1),',z=',tab_indz((z-1)*(meanz+1)+mz+1),')'

!**********Storage of all spectra
spectra_global(:,number_global) = spectra(:)
spectra_b_global(:,number_global) = spectra_b(:)
spectra_c_global(:,number_global) = spectra_c(:)

!**********Some cleaning
deallocate(br)
deallocate(bi)
deallocate(cr)
deallocate(ci)
deallocate(wvhaec)
deallocate(dwork)
deallocate(work)
deallocate(norm)
deallocate(norm_b)
deallocate(norm_c)
deallocate(spectra)
deallocate(spectra_b)
deallocate(spectra_c)
deallocate(zonal_wind)
deallocate(merid_wind)
deallocate(zonal_wind_lmdz)
deallocate(merid_wind_lmdz)

!**********
!End of loop over time and altitude indices
end do
end do
write(*,'(a24,i5,a3,i5,a8)') 'spectra computed for (t=',tab_indt((t-1)*(meant+1)+1),',z=',tab_indz((z-1)*(meanz+1)+1),') config'
write(*,*) '**********'
end do
end do

!**********
!Spectra mean
do t = 1,n_indt
do z = 1,n_indz
number_config = (t-1)*n_indz+z
spectra_config(:,number_config) = 0.
do mt = 0,meant
do mz = 0,meanz
number_local = mt*(meanz+1)+mz+1
number_global = (number_config-1)*(meant+1)*(meanz+1) + number_local
spectra_config(:,number_config) = spectra_config(:,number_config) + &
		spectra_global(:,number_global)
spectra_b_config(:,number_config) = spectra_b_config(:,number_config) + &
		spectra_b_global(:,number_global)
spectra_c_config(:,number_config) = spectra_c_config(:,number_config) + &
		spectra_c_global(:,number_global)
end do
end do
spectra_config(:,number_config) = spectra_config(:,number_config)/(meanz+1)/(meant+1)
spectra_b_config(:,number_config) = spectra_b_config(:,number_config)/(meanz+1)/(meant+1)
spectra_c_config(:,number_config) = spectra_c_config(:,number_config)/(meanz+1)/(meant+1)
end do
end do

!**********Writing header of file_spectra
open(unit=10,file=file_spectra,action="write",position="rewind")
write(10,'(a10,a2)',advance='no') '#Spherical',' '
do t = 1,n_indt
do z = 1,n_indz
  write(10,'(a50)',advance='no') file_netcdf
end do
end do
write(10,*)
write(10,'(a11,a1)',advance='no') '#wavenumber',' '
do t = 1,n_indt
do z = 1,n_indz
  write(10,'(a2,i5,a3,i5,a35)',advance='no') 't=',tab_indt((t-1)*(meant+1)+1),' z=',tab_indz((z-1)*(meanz+1)+1),' '
end do
end do
write(10,*)
write(10,'(a1,a11)',advance='no') '#',' '
do t = 1,n_indt
do z = 1,n_indz
if (is_div_rot) then
  write(10,'(a8,a8,a8,a8,a8,a10)',advance='no') 'spec_tot',' ','spec_div',' ','spec_rot',' '
else
  write(10,'(a8,a42)',advance='no') 'spec_tot',' '
end if
end do
end do
write(10,*)
!**********Writing file_spectra
do j=1,nLat
  write(10,'(i5,a6)',advance='no') j-1,' '
  do t = 1,n_indt
  do z = 1,n_indz
  number_config = (t-1)*n_indz+z
  if (is_div_rot) then
    write(10,'(e13.6E2,a3,e13.6E2,a3,e13.6E2,a5)',advance='no') spectra_config(j,number_config),' ',&
    			spectra_b_config(j,number_config),' ',spectra_c_config(j,number_config),' '
  else
    write(10,'(e13.6E2,a37)',advance='no') spectra_config(j,number_config),' '
  end if
  end do
  end do
  if (j /= nlat) write(10,*)
end do
close(10) 

print *, "***SUCCESS writing ",trim(file_spectra)

!***********Some cleaning
deallocate(spectra_config)
deallocate(spectra_b_config)
deallocate(spectra_c_config)
deallocate(spectra_global)
deallocate(spectra_b_global)
deallocate(spectra_c_global)
deallocate(tab_indt)
deallocate(tab_indz)
deallocate(tmp_tab_int)
deallocate(arg)

contains
  subroutine check(status)
    integer, intent (in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

END PROGRAM spectra_analysis
