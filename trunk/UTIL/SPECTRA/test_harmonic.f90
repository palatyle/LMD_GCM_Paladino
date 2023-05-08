!M. Indurain
!March 2015

PROGRAM lmdz_harmonic

!Computes wind field following some vectorial harmonics with spherepack routines.
!Wind field follows lmdz convention:
!  u(nlong+1,nlat) zonal eastward wind, v(nlong+1,nlat) latitudinal wind
!Grid resolution is given in the command line:
!./test_harmonic -lat n_latitudinal_bands -lon n_longitudinal_bands. 
!By default, n_latitudinal_bands = 48. n_longitudinal_bands = 64.

use netcdf

implicit none

!List of harmonics we want
 character (len=100), dimension(5) :: harmonic=(/'2','10','22','6_3','15_20'/)
!Physical grid
integer :: nlat,nlatm1,nlong,nlongp1
double precision, dimension(:), allocatable :: latitude,longitude
double precision, dimension(:,:), allocatable :: merid_wind,zonal_wind !spherepack convention
double precision, dimension(:,:), allocatable :: merid_wind_lmdz,zonal_wind_lmdz !lmdz convention
!Spectral grid
integer :: mmax !spectra truncation mmax=min(nalt,nlong/2) or min(nalt,(nlong+1)/2)
double precision, dimension (:,:), allocatable :: br,bi,cr,ci !spectra coefficients
double precision, dimension (:,:), allocatable :: norm,norm_b,norm_c !norm of spectra coefficients
double precision, dimension (:), allocatable :: spectra,spectra_b,spectra_c
!Spherepack declarations
integer :: ierror,lvhsec,ldwork,lwork,l1,l2,nt=1,ityp=0,idvw,jdvw,mdab,ndab
double precision, dimension (:), allocatable :: wvhsec,work
double precision, dimension (:), allocatable :: dwork
!Netcdf stuff
integer :: idfile,idmerid_wind,idzonal_wind
integer :: idlat,idlong,idlat_var,idlong_var
integer :: idmdab,idndab,idbr,idbi,idcr,idci
integer :: idspectra,idspectra_b,idspectra_c
character (len=100) :: file_netcdf
!Other
integer :: i,j,k
double precision, parameter :: pi=2.*ASIN(1.)
character (len=100) :: tmp_char
!Input arguments
character (len=100), dimension(:), allocatable :: arg
integer :: narg

!**********
!Initialisation
nlatm1=48
nlong=64

!**********
!Input reading
narg = command_argument_count()
allocate(arg(narg))
do i = 1, narg
  call get_command_argument(i,arg(i))
end do
i = 1
do while (i .le. narg)
  if (arg(i) == '-lat' .or. arg(i) == '-lon') then
    select case(arg(i))
    case('-lat')
      read(arg(i+1),'(i10)' ) nlatm1 !number of latitudinal bands		
    case('-lon')
      read(arg(i+1),'(i10)' ) nlong !number of longitudinal bands
    end select
    i = i + 2
  elseif (arg(i) == '-h' .or. arg(i) == '--help') then
      print*,'Usage\n test_harmonic [option]\n [-h or --help]\t: brief help'
      print*,'[-lat int]\t: number of latitudinal bands (default: 48)'
      print*,'[-lon int]\t: number of longitudinal bands (default: 64)'
      stop 'End help'
  else
    print*,'No option ',trim(arg(i))
    stop
  end if
end do

!**********
!Physical grid
nlat=nlatm1+1
nlongp1=nlong+1
print*,'harmonics computed on a ',nlatm1,' latitude bands x ',nlong,' longitude bands.'

allocate(latitude(nlat))
allocate(longitude(nlongp1))
do j=1,nlat
	latitude(j)=90.-(j-1)*180./(nlat-1)
end do
do i=1,nlongp1
	longitude(i)=-180.+(i-1)*2*180./nlong
end do
allocate(zonal_wind(nlat,nlong))
allocate(merid_wind(nlat,nlong))
allocate(zonal_wind_lmdz(nlongp1,nlat))
allocate(merid_wind_lmdz(nlongp1,nlat))

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
lvhsec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlong+15
allocate(wvhsec(lvhsec))
wvhsec(:) = 0.
ldwork=2*(nlat+2)
allocate(dwork(ldwork))
dwork(:) = 0.
ierror=3
call vhseci(nlat,nlong,wvhsec,lvhsec,dwork,ldwork,ierror)
 
!**********Vhseci function result
select case (ierror)
  case(1) 
    print*,'Vhseci: ERROR on nlat'
  case(2) 
    print*,'Vhseci: ERROR on nlong'
  case(3) 
    print*,'Vhseci: ERROR on lvhsec'
  case(4) 
    print*,'Vhseci: ERROR on ldwork'
end select

!**********Loop over all harmonics
do k=1,SIZE(harmonic)
file_netcdf='harmonic_'//trim(integer2string(nlongp1-1))//'x'//&
		trim(integer2string(nlat-1))//'_lmdz_'//trim(harmonic(k))//'.nc'

!**********Vhsec function
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
select case(harmonic(k))
  case('2')
    br(3,3)=1.
  case('10')
    br(11,11)=1.
  case('22')
    br(23,23)=1.
  case('6_3')
    br(4,4)=1.
    br(7,7)=1.
  case('15_20')
    br(16,16)=1.
    br(21,21)=1.
end select 
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
lvhsec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlong+15
if (ityp .le. 2) then
  lwork=nlat*(2*nt*nlong+max(6*l2,nlong))
else
  lwork=l2*(2*nt*nlong+max(6*nlat,nlong))
end if
allocate(work(lwork))
work(:) = 0.
ierror=3
call vhsec(nlat,nlong,ityp,nt,merid_wind,zonal_wind,idvw,jdvw,br,bi,cr,ci,mdab,ndab,wvhsec,lvhsec,work,lwork,ierror)

!**********Vhsec function result
select case (ierror)
  case(1) 
    print*,'Vhsec: ERROR on nlat'
  case(2) 
    print*,'Vhsec: ERROR on nlong'
  case(3) 
    print*,'Vhsec: ERROR on ityp'
  case(4) 
    print*,'Vhsec: ERROR on nt'
  case(5) 
    print*,'Vhsec: ERROR on idvw'
  case(6) 
    print*,'Vhsec: ERROR on jdvw'
  case(7) 
    print*,'Vhsec: ERROR on mdab'
  case(8) 
    print*,'Vhsec: ERROR on ndab'
  case(9) 
    print*,'Vhsec: ERROR on lvhsec'
  case(10) 
    print*,'Vhsec: ERROR on lwork'
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

!**********From spherepack to lmdz convention
zonal_wind_lmdz(1:nlong,:)=transpose(zonal_wind(:,:))
merid_wind_lmdz(1:nlong,:)=transpose(merid_wind(:,:))
zonal_wind_lmdz(nlongp1,:)=zonal_wind_lmdz(1,:)
merid_wind_lmdz(nlongp1,:)=merid_wind_lmdz(1,:)

!**********Write netcdf file
  call check(nf90_create(trim(file_netcdf),NF90_CLOBBER,idfile))
 
  call check( nf90_def_dim(idfile,"mdab",mdab,idmdab))
  call check( nf90_def_dim(idfile,"ndab",ndab,idndab))
  call check(nf90_def_var(idfile,"br",NF90_DOUBLE,(/idmdab,idndab/),idbr))
  call check(nf90_def_var(idfile,"bi",NF90_DOUBLE,(/idmdab,idndab/),idbi))
  call check(nf90_def_var(idfile,"cr",NF90_DOUBLE,(/idmdab,idndab/),idcr))
  call check(nf90_def_var(idfile,"ci",NF90_DOUBLE,(/idmdab,idndab/),idci))
  call check(nf90_def_var(idfile,"spectra",NF90_DOUBLE,(/idndab/),idspectra))
  call check(nf90_def_var(idfile,"spectra_div",NF90_DOUBLE,(/idndab/),idspectra_b))
  call check(nf90_def_var(idfile,"spectra_vort",NF90_DOUBLE,(/idndab/),idspectra_c))
  
  call check( nf90_def_dim(idfile,"latitude",nlat,idlat))
  call check( nf90_def_dim(idfile,"longitude",nlongp1,idlong))
  call check(nf90_def_var(idfile,"latitude",NF90_DOUBLE,(/idlat/),idlat_var))
  call check(nf90_def_var(idfile,"longitude",NF90_DOUBLE,(/idlong/),idlong_var))
  call check(nf90_def_var(idfile,"u",NF90_DOUBLE,(/idlong,idlat/),idzonal_wind))
  call check(nf90_def_var(idfile,"v",NF90_DOUBLE,(/idlong,idlat/),idmerid_wind))
  ! End define mode. This tells netCDF we are done defining metadata.
  call check(nf90_enddef(idfile))

  ! Write the pretend data to the file. Although netCDF supports
  ! reading and writing subsets of data, in this case we write all the
  ! data in one operation.
  call check(nf90_put_var(idfile,idbr,br))
  call check(nf90_put_var(idfile,idbi,bi))
  call check(nf90_put_var(idfile,idcr,cr))
  call check(nf90_put_var(idfile,idci,ci))
  call check(nf90_put_var(idfile,idspectra,spectra))
  call check(nf90_put_var(idfile,idspectra_b,spectra_b))
  call check(nf90_put_var(idfile,idspectra_c,spectra_c))
  
  call check(nf90_put_var(idfile,idlat_var,latitude))
  call check(nf90_put_var(idfile,idlong_var,longitude))
  call check(nf90_put_var(idfile,idzonal_wind,zonal_wind_lmdz))
  call check(nf90_put_var(idfile,idmerid_wind,merid_wind_lmdz))
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check(nf90_close(idfile))

  deallocate(br)
  deallocate(bi)
  deallocate(cr)
  deallocate(ci)
  deallocate(work)
  deallocate(norm)
  deallocate(norm_b)
  deallocate(norm_c)
  deallocate(spectra)
  deallocate(spectra_b)
  deallocate(spectra_c)
  
  
  print*,trim(file_netcdf),' computed!'
end do


contains
  subroutine check(status)
    integer, intent (in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
  
  function integer2string(n)
    implicit none
    character(len=10) :: integer2string
    integer, intent(in) :: n
    integer :: nn

    nn=n
    if (n .lt. 0) nn=-n
    if (nn .ge. 0 .and. nn .lt. 10) write(integer2string,'(i1)') nn
    if (nn .ge. 10 .and. nn .lt. 100) write(integer2string,'(i2)') nn
    if (nn .ge. 100 .and. nn .lt. 1000) write(integer2string,'(i3)') nn
    if (nn .ge. 1000 .and. nn .lt. 10000) write(integer2string,'(i4)') nn
    if (nn .ge. 10000 .and. nn .lt. 100000) write(integer2string,'(i4)') nn
    if (nn .ge. 100000 .and. nn .lt. 1000000) write(integer2string,'(i4)') nn
    if (nn .ge. 1000000 .and. nn .lt. 10000000) write(integer2string,'(i4)') nn
    if (nn .ge. 10000000 .and. nn .lt. 100000000) write(integer2string,'(i4)') nn
    
  end function integer2string
  
END PROGRAM lmdz_harmonic
