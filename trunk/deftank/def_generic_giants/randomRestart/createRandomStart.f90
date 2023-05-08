program createRandomStart

USE netcdf

implicit none

!=======================================================================
! Copy a start.nc file and change ucov, vcov, teta variables:
!	teta -> teta +/- dTeta (K)
!	ucov -> (0 +/- duNat)*cu (m/s)
!	vcov -> (0 +/- dvNat)*cv (m/s)
! Perturbations are vertically the same.
!=======================================================================

!Input
 character (len=100) :: fileNetcdfInput ! input start.nc netcdf file
!Output
 character (len=100) :: fileNetcdf ! output start.nc netcdf file
!Usefull stuff
integer :: nLat,nRlatv,nLong,nRlonu,nAlt,nTime ! dimensions  in netcdf file (redundant point in longitude)
real, dimension (:,:,:,:), allocatable :: uCov,vCov,uNat,vNat,uNat2,vNat2,teta
real, dimension (:,:), allocatable :: cu,cv
real :: duNat,dvNat,dTeta
!Netcdf declarations
integer :: idFile,idStuff,idUcov,idVcov,idCu,idCv,idTeta
integer :: ierror
!Temporary variables
 character (len=100) :: format,varTmp
integer :: i,j,k,l,N=1,dimFilter
logical :: isFile
real :: alea
 character (len=200) :: command


print*,'input start.nc?'
read*, fileNetcdfInput
inquire(file=fileNetcdfInput,exist=isFile)
if (isFile .eqv. .false.) then
	stop 'input doesn t exist'
end if

print*,'output start.nc?'
read*, fileNetcdf
inquire(file=fileNetcdf,exist=isFile)
if (fileNetcdfInput == fileNetcdf) then
	!name of the input netcdf file back-up
	varTmp=trim(fileNetcdfInput) // '.old'
	inquire(file=varTmp,exist=isFile)
	do while(isFile .eqv. .true.)
		varTmp=trim(varTmp) // '.old'
		print*,varTmp
		inquire(file=varTmp,exist=isFile)
	end do
	!input netcdf file back-up
	write(command,'(a6,a50,a1,a50)') 'cp -f ',trim(fileNetcdfInput),' ',trim(varTmp)
	call system(command)
else
	inquire(file=fileNetcdf,exist=isFile)
	if (isFile .eqv. .true.) then
		print*,'output name existing. Want to overwrite? (y/n)'
		read*,varTmp
		if (varTmp == 'y') then
			write(command,'(a6,a50)') 'rm -f ',fileNetcdf
			print*,command
 			call system(command)
		else
			stop 'existing output not deleted'
		end if
	end if
end if

print*,'du (m/s)?'
read*, duNat

print*,'dv (m/s)?'
read*, dvNat

print*,'dT (K)?'
read*, dTeta

!print*,'dim average filter? (0 --> no | 1,2,3... --> average)'
!read*, dimFilter
dimFilter=0

!Creation of the output netcdf file if different from input netcdf file
if (fileNetcdfInput /= fileNetcdf) then
	write(command,'(a6,a50,a1,a50)') 'cp -f ',fileNetcdfInput,' ',fileNetcdf
	print*,command
 	call system(command)
end if

!**********
!Netcdf file reading
PRINT*,'**********'
PRINT*,'Netcdf file reading'
!Open
ierror=nf90_open(fileNetcdf,NF90_WRITE,idFile)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped open"
end if
!Read dimension latitude
ierror=nf90_inq_dimid(idFile,'latitude',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nLat)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim lat"
end if
!Read dimension latitude bis
ierror=nf90_inq_dimid(idFile,'rlatv',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nRlatv)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim rlatv"
end if
!Read dimension longitude
ierror=nf90_inq_dimid(idFile,'longitude',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nLong)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim long"
end if
!Read dimension longitude bis
ierror=nf90_inq_dimid(idFile,'rlonu',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nRlonu)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim rlonu"
end if
!Read dimension altitude
ierror=nf90_inq_dimid(idFile,'altitude',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nAlt)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim alt"
end if
!Read dimension time
ierror=nf90_inq_dimid(idFile,'Time',idStuff)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_inquire_dimension(idFile,idStuff,varTmp,nTime)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped dim time"
end if
!Read variables ucov and vcov
allocate(uCov(nRlonu,nLat,nAlt,nTime))
allocate(vCov(nLong,nRlatv,nAlt,nTime))
ierror=nf90_inq_varid(idFile,"ucov",idUcov)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_get_var(idFile,idUcov,uCov)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped get_var ucov"
end if
ierror=nf90_inq_varid(idFile,"vcov",idVcov)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_get_var(idFile,idVcov,vCov)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped get_var vcov"
end if
!Read variables cu and cv
allocate(cu(nRlonu,nLat))
allocate(cv(nLong,nRlatv))
ierror=nf90_inq_varid(idFile,"cu",idCu)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_get_var(idFile,idCu,cu)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped get_var cu"
end if
ierror=nf90_inq_varid(idFile,"cv",idCv)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_get_var(idFile,idCv,cv)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped get_var cv"
end if
!Read variables teta
allocate(teta(nLong,nLat,nAlt,nTime))
ierror=nf90_inq_varid(idFile,"teta",idTeta)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped"
end if
ierror=nf90_get_var(idFile,idTeta,teta)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped get_var teta"
end if

!**********
!uNat
allocate(uNat(nRlonu,nLat,nAlt,nTime))
uNat(:,:,:,:)=0.0
do i=1,nLong
do j=1,nLat
do k=1,nAlt
do l=1,nTime
	call random_number(alea)
	uNat(i,j,k,l) = -duNat+2*duNat*alea
end do
end do
end do
end do
!-/+180 longitude points must be the same
uNat(1,:,:,:)=uNat(nLong,:,:,:)
!**********
!vNat
allocate(vNat(nLong,nRlatv,nAlt,nTime))
vNat(:,:,:,:)=0.0
do i=1,nLong
do j=1,nRlatv
do k=1,nAlt
do l=1,nTime
	call random_number(alea)
	vNat(i,j,k,l) = -dvNat+2*dvNat*alea
end do
end do
end do
end do
!-/+180 longitude points must be the same
vNat(1,:,:,:)=vNat(nLong,:,:,:)
!**********
!mean (if dimFilter=0 mean mean_filter has no effect)
allocate(uNat2(nRlonu,nLat,nAlt,nTime))
allocate(vNat2(nLong,nRlatv,nAlt,nTime))
uNat2(:,:,:,:)=uNat(:,:,:,:)
vNat2(:,:,:,:)=vNat(:,:,:,:)
do k=1,nAlt
do l=1,nTime
	call mean_filter(uNat(:,:,k,l),nRlonu,nLat,dimFilter,uNat2(:,:,k,l))
	call mean_filter(vNat(:,:,k,l),nLong,nRlatv,dimFilter,vNat2(:,:,k,l))
end do
end do
!**********
!new ucov and vcov
do k=1,nAlt
do l=1,nTime
	uCov(:,:,k,l)=uNat2(:,:,k,l)*cu
	vCov(:,:,k,l)=vNat2(:,:,k,l)*cv
end do
end do
!**********
!teta
do i=1,nLong
do j=1,nRlatv
do k=1,nAlt
do l=1,nTime
	call random_number(alea)
	teta(i,j,k,l) = teta(i,j,k,l)-dTeta+2*dTeta*alea
end do
end do
end do
end do
teta(1,:,:,:)=teta(nLong,:,:,:)
!**********
!Write variable ucov and vcov
ierror=nf90_put_var(idFile,idUcov,uCov) !,(1,1,1,1),(nLong,nLat,nAlt,nTime),
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped put ucov"
end if
ierror=nf90_put_var(idFile,idVcov,vCov) !,(1,1,1,1),(nLong,nLat,nAlt,nTime),
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped put vcov"
end if
!**********
!Write variable teta
ierror=nf90_put_var(idFile,idTeta,teta) !,(1,1,1,1),(nLong,nLat,nAlt,nTime),
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped put teta"
end if
!**********
!Close
ierror=nf90_close(idFile)
if(ierror /= nf90_noerr) then 
     print *, trim(nf90_strerror(ierror))
     stop "Stopped close"
end if

contains

subroutine interface_user(fileNetcdf,duNat,dvNat,dTeta)

implicit none

 character (len=*), intent(out) :: fileNetcdf
real, intent(out) :: duNat,dvNat,dTeta






end subroutine interface_user

subroutine mean_filter(matrix,m,n,dimFilter,matrixOut)

implicit none

integer, intent(in) :: m,n,dimFilter
real, dimension(m,n), intent(in) :: matrix
real, dimension(m,n), intent(out) :: matrixOut
real, dimension(m+2*dimFilter,n+2*dimFilter) :: matrixTmp,matrixConv
real, dimension(2*dimFilter+1,2*dimFilter+1) :: kernel
integer :: mm,nn

!matrix containing input matrix plus zero edges
matrixTmp(:,:)=0.0
mm=m+2*dimFilter
nn=n+2*dimFilter
do i=dimFilter+1,mm-dimFilter
	do j=dimFilter+1,nn-dimFilter
		matrixTmp(i,j)=matrix(i-dimFilter,j-dimFilter)
	end do
end do
!print*,matrixTmp(5,:)

!filter creation
kernel(:,:) = 1.0/((2*dimFilter+1)*(2*dimFilter+1))
!print*,kernel

!matrix convolution with the filter
do i=dimFilter+1,mm-dimFilter
	do j=dimFilter+1,nn-dimFilter
		matrixConv(i,j)=sum(matrixTmp(i-dimFilter:i+dimFilter,j-dimFilter:j+dimFilter)*kernel)
	end do
end do
i=5
j=5
!print*,shape(matrixTmp(i-dimFilter:i+dimFilter,j-dimFilter:j+dimFilter))
!print*,matrixConv(5,:)

!extraction of submatrix without zeros edges
matrixOut(:,:)=0.0
do i=1,m
	do j=1,n
		matrixOut(i,j)=matrixConv(i+dimFilter,j+dimFilter)
	end do
end do

end subroutine mean_filter

end program createRandomStart

