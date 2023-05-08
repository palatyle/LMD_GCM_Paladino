










subroutine wstats(ngrid,nom,titre,unite,dim,px)

use statto_mod, only: istats,istime,count
use mod_phys_lmdz_para, only : is_mpi_root, is_master, gather, klon_mpi_begin
use mod_grid_phy_lmdz, only : klon_glo, Grid1Dto2D_glo, &
                              nbp_lon, nbp_lat, nbp_lev
implicit none

include "netcdf.inc"

integer,intent(in) :: ngrid
character (len=*),intent(in) :: nom,titre,unite
integer,intent(in) :: dim
real,intent(in) :: px(ngrid,nbp_lev)
real,allocatable,save :: mean3d(:,:,:),sd3d(:,:,:),dx3(:,:,:)
real,allocatable,save :: mean2d(:,:),sd2d(:,:),dx2(:,:)
character (len=50) :: namebis
character (len=50), save :: firstvar
!$OMP THREADPRIVATE(firstvar)
integer :: ierr,varid,nbdim,nid
integer :: meanid,sdid
integer, dimension(4)  :: id,start,sizes
logical, save :: firstcall=.TRUE.
integer :: l,i,j,ig0
integer,save :: indx

integer, save :: step=0
!$OMP THREADPRIVATE(firstcall,indx,step)

! Added to work in parallel mode
! When not running in parallel mode:
real px3_glop(ngrid,nbp_lev) ! to store a 3D data set on global physics grid
real px3_glo(nbp_lon,nbp_lat,nbp_lev) ! to store a global 3D data set on global lonxlat grid
real px2_glop(ngrid) ! to store a 2D data set on global physics grid
real px2_glo(nbp_lon,nbp_lat) ! to store a 2D data set on global lonxlat grid

! 1. Initialization (creation of stats.nc file)
if (firstcall) then
   firstcall=.false.
   firstvar=trim((nom))
   call inistats(ierr)
   if (klon_glo>1) then ! general case, 3D GCM
     allocate(mean3d(nbp_lon+1,nbp_lat,nbp_lev))
     allocate(sd3d(nbp_lon+1,nbp_lat,nbp_lev))
     allocate(dx3(nbp_lon+1,nbp_lat,nbp_lev))
     allocate(mean2d(nbp_lon+1,nbp_lat))
     allocate(sd2d(nbp_lon+1,nbp_lat))
     allocate(dx2(nbp_lon+1,nbp_lat))
   else ! 1D model
     allocate(mean3d(1,1,nbp_lev))
     allocate(sd3d(1,1,nbp_lev))
     allocate(dx3(1,1,nbp_lev))
     allocate(mean2d(1,1))
     allocate(sd2d(1,1))
     allocate(dx2(1,1))
   endif
endif

if (firstvar==nom) then ! If we're back to the first variable, increment time counter
      step = step + 1
endif

if (mod(step,istats).ne.0) then
  ! if its not time to write to file, exit
   RETURN
endif

! collect fields on a global physics grid
  if (dim.eq.3) then
    px3_glop(:,1:nbp_lev)=px(:,1:nbp_lev)
!  Passage variable physique -->  variable dynamique
    DO l=1,nbp_lev
      DO i=1,nbp_lon
         px3_glo(i,1,l)=px(1,l)
         px3_glo(i,nbp_lat,l)=px(ngrid,l)
      ENDDO
      DO j=2,nbp_lat-1
         ig0= 1+(j-2)*nbp_lon
         DO i=1,nbp_lon
            px3_glo(i,j,l)=px(ig0+i,l)
         ENDDO
      ENDDO
    ENDDO
  else ! dim.eq.2
    px2_glop(:)=px(:,1)
!    Passage variable physique -->  physique dynamique
   DO i=1,nbp_lon
     px2_glo(i,1)=px(1,1)
     px2_glo(i,nbp_lat)=px(ngrid,1)
   ENDDO
   DO j=2,nbp_lat-1
     ig0= 1+(j-2)*nbp_lon
     DO i=1,nbp_lon
        px2_glo(i,j)=px(ig0+i,1)
     ENDDO
   ENDDO
  endif

! 2. Write field to file

if (is_master) then
! only master needs do this

ierr = NF_OPEN("stats.nc",NF_WRITE,nid)

namebis=trim(nom)

! test: check if that variable already exists in file
ierr= NF_INQ_VARID(nid,namebis,meanid)

if (ierr.ne.NF_NOERR) then
  ! variable not in file, create/define it
   if (firstvar==nom) then 
      indx=1
      count(:)=0
   endif


!declaration de la variable

! choix du nom des coordonnees
   ierr= NF_INQ_DIMID(nid,"longitude",id(1))
   ierr= NF_INQ_DIMID(nid,"latitude",id(2))
   if (dim.eq.3) then
      ierr= NF_INQ_DIMID(nid,"altitude",id(3))
      ierr= NF_INQ_DIMID(nid,"Time",id(4))
      nbdim=4
   else if (dim==2) then
      ierr= NF_INQ_DIMID(nid,"Time",id(3))
      nbdim=3
   endif

   write (*,*) "====================="
   write (*,*) "STATS: creation de ",nom
   namebis=trim(nom)
   call def_var_stats(nid,namebis,titre,unite,nbdim,id,meanid,ierr)
   if (dim.eq.3) then
     call inivar(nid,meanid,size(px3_glop,1),dim,indx,px3_glop,ierr)
   else ! dim.eq.2
     call inivar(nid,meanid,size(px2_glop,1),dim,indx,px2_glop,ierr)
   endif
   namebis=trim(nom)//"_sd"
   call def_var_stats(nid,namebis,trim(titre)//" total standard deviation over the season",unite,nbdim,id,sdid,ierr)
   if (dim.eq.3) then
     call inivar(nid,sdid,size(px3_glop,1),dim,indx,px3_glop,ierr)
   else ! dim.eq.2
     call inivar(nid,sdid,size(px2_glop,1),dim,indx,px2_glop,ierr)
   endif

   ierr= NF_CLOSE(nid)
   return

else
   ! variable found in file
   namebis=trim(nom)//"_sd"
   ierr= NF_INQ_VARID(nid,namebis,sdid)

endif

if (firstvar==nom) then 
   count(indx)=count(int(indx))+1
   indx=indx+1
   if (indx>istime) then
      indx=1
   endif
endif

if (count(indx)==0) then
   ! very first time we write the variable to file
   if (dim.eq.3) then
      start=(/1,1,1,indx/)
      if (klon_glo>1) then !general case
        sizes=(/nbp_lon+1,nbp_lat,nbp_lev,1/)
      else
        sizes=(/1,1,nbp_lev,1/)
      endif
      mean3d(:,:,:)=0
      sd3d(:,:,:)=0
   else if (dim.eq.2) then
      start=(/1,1,indx,0/)
      if (klon_glo>1) then !general case
        sizes=(/nbp_lon+1,nbp_lat,1,0/)
      else
        sizes=(/1,1,1,0/)
      endif
      mean2d(:,:)=0
      sd2d(:,:)=0
   endif
else
   ! load values from file
   if (dim.eq.3) then
      start=(/1,1,1,indx/)
      if (klon_glo>1) then !general case
        sizes=(/nbp_lon+1,nbp_lat,nbp_lev,1/)
      else ! 1D model case
        sizes=(/1,1,nbp_lev,1/)
      endif
      ierr = NF_GET_VARA_DOUBLE(nid,meanid,start,sizes,mean3d)
      ierr = NF_GET_VARA_DOUBLE(nid,sdid,start,sizes,sd3d)
      if (ierr.ne.NF_NOERR) then
         write (*,*) "wstats error reading :",trim(nom)
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif

   else if (dim.eq.2) then
      start=(/1,1,indx,0/)
      if (klon_glo>1) then ! general case
        sizes=(/nbp_lon+1,nbp_lat,1,0/)
      else
        sizes=(/1,1,1,0/)
      endif
      ierr = NF_GET_VARA_DOUBLE(nid,meanid,start,sizes,mean2d)
      ierr = NF_GET_VARA_DOUBLE(nid,sdid,start,sizes,sd2d)
      if (ierr.ne.NF_NOERR) then
         write (*,*) "wstats error reading :",trim(nom)
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif
   endif
endif ! of if (count(indx)==0)


! 2.1. Build dx* (data on lon-lat grid, with redundant longitude)

if (dim.eq.3) then
  dx3(1:nbp_lon,:,:)=px3_glo(:,:,:)
  IF (klon_glo>1) THEN ! in 3D, add redundant longitude point
    dx3(nbp_lon+1,:,:)=dx3(1,:,:)
  ENDIF
else ! dim.eq.2
  dx2(1:nbp_lon,:)=px2_glo(:,:)
  IF (klon_glo>1) THEN ! in 3D, add redundant longitude point
    dx2(nbp_lon+1,:)=dx2(1,:)
  ENDIF
endif


! 2.2. Add current values to previously stored sums

if (dim.eq.3) then

   mean3d(:,:,:)=mean3d(:,:,:)+dx3(:,:,:)
   sd3d(:,:,:)=sd3d(:,:,:)+dx3(:,:,:)**2

   ierr = NF_PUT_VARA_DOUBLE(nid,meanid,start,sizes,mean3d)
   ierr = NF_PUT_VARA_DOUBLE(nid,sdid,start,sizes,sd3d)
  if (ierr.ne.NF_NOERR) then
     write (*,*) "wstats error writing :",trim(nom)
     write (*,*) NF_STRERROR(ierr)
     stop ""
  endif

else if (dim.eq.2) then

   mean2d(:,:)= mean2d(:,:)+dx2(:,:)
   sd2d(:,:)=sd2d(:,:)+dx2(:,:)**2

   ierr = NF_PUT_VARA_DOUBLE(nid,meanid,start,sizes,mean2d)
   ierr = NF_PUT_VARA_DOUBLE(nid,sdid,start,sizes,sd2d)
  if (ierr.ne.NF_NOERR) then
     write (*,*) "wstats error writing :",trim(nom)
     write(*,*) "start:",start
     write(*,*) "sizes:",sizes
     write(*,*) "mean2d:",mean2d
     write(*,*) "sd2d:",sd2d
     write (*,*) NF_STRERROR(ierr)
     stop ""
  endif

endif ! of if (dim.eq.3) elseif (dim.eq.2)

  ierr= NF_CLOSE(nid)
endif ! of if (is_master)

end subroutine wstats

!======================================================
subroutine inivar(nid,varid,ngrid,dim,indx,px,ierr)
use mod_grid_phy_lmdz, only : nbp_lon, nbp_lat, nbp_lev, klon_glo

implicit none

include "netcdf.inc"

integer, intent(in) :: nid,varid,dim,indx,ngrid
real, dimension(ngrid,nbp_lev), intent(in) :: px
integer, intent(out) :: ierr

integer :: l,i,j,ig0
integer, dimension(4) :: start,sizes
real, dimension(nbp_lon+1,nbp_lat,nbp_lev) :: dx3
real, dimension(nbp_lon+1,nbp_lat) :: dx2
real :: dx3_1d(nbp_lev) ! for 1D outputs
real :: dx2_1d ! for 1D outputs

if (dim.eq.3) then

   start=(/1,1,1,indx/)
   if (klon_glo>1) then ! general 3D case
     sizes=(/nbp_lon+1,nbp_lat,nbp_lev,1/)
   else
     sizes=(/1,1,nbp_lev,1/)
   endif

!  Passage variable physique -->  variable dynamique

   if (klon_glo>1) then ! general case
    DO l=1,nbp_lev
      DO i=1,nbp_lon+1
         dx3(i,1,l)=px(1,l)
         dx3(i,nbp_lat,l)=px(ngrid,l)
      ENDDO
      DO j=2,nbp_lat-1
         ig0= 1+(j-2)*nbp_lon
         DO i=1,nbp_lon
            dx3(i,j,l)=px(ig0+i,l)
         ENDDO
         dx3(nbp_lon+1,j,l)=dx3(1,j,l)
      ENDDO
    ENDDO
   else ! 1D model case
     dx3_1d(1:nbp_lev)=px(1,1:nbp_lev)
   endif

   if (klon_glo>1) then
     ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,sizes,dx3)
   else
     ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,sizes,dx3_1d)
   endif
  if (ierr.ne.NF_NOERR) then
     write (*,*) "inivar error writing variable"
     write (*,*) NF_STRERROR(ierr)
     stop ""
  endif

else if (dim.eq.2) then

      start=(/1,1,indx,0/)
      if (klon_glo>1) then ! general 3D case
        sizes=(/nbp_lon+1,nbp_lat,1,0/)
      else
        sizes=(/1,1,1,0/)
      endif

!    Passage variable physique -->  physique dynamique

  if (klon_glo>1) then ! general case
  DO i=1,nbp_lon+1
     dx2(i,1)=px(1,1)
     dx2(i,nbp_lat)=px(ngrid,1)
  ENDDO
  DO j=2,nbp_lat-1
     ig0= 1+(j-2)*nbp_lon
     DO i=1,nbp_lon
        dx2(i,j)=px(ig0+i,1)
     ENDDO
     dx2(nbp_lon+1,j)=dx2(1,j)
  ENDDO
  else ! 1D model case
    dx2_1d=px(1,1)
  endif
  
   if (klon_glo>1) then
     ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,sizes,dx2)
   else
     ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,sizes,dx2_1d)
   endif
  if (ierr.ne.NF_NOERR) then
     write (*,*) "inivar error writing variable"
     write (*,*) NF_STRERROR(ierr)
     stop ""
  endif

endif

end subroutine inivar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine def_var_stats(nid,name,title,units,nbdim,dimids,nvarid,ierr)

! This subroutine defines variable 'name' in a (pre-existing and opened)
! NetCDF file (known from its NetCDF ID 'nid').
! The number of dimensions 'nbdim' of the variable, as well as the IDs of
! corresponding dimensions must be set (in array 'dimids').
! Upon successfull definition of the variable, 'nvarid' contains the
! NetCDF ID of the variable.
! The variables' attributes 'title' (Note that 'long_name' would be more
! appropriate) and 'units' are also set. 

implicit none

include "netcdf.inc"

integer,intent(in) :: nid ! NetCDF file ID
character(len=*),intent(in) :: name ! the variable's name
character(len=*),intent(in) :: title ! 'title' attribute of variable
character(len=*),intent(in) :: units ! 'units' attribute of variable
integer,intent(in) :: nbdim ! number of dimensions of the variable
integer,dimension(nbdim),intent(in) :: dimids ! NetCDF IDs of the dimensions
                                              ! the variable is defined along
integer,intent(out) :: nvarid ! NetCDF ID of the variable
integer,intent(out) :: ierr ! returned NetCDF staus code

! 1. Switch to NetCDF define mode 
ierr=NF_REDEF(nid)

! 2. Define the variable
ierr = NF_DEF_VAR (nid,adjustl(name),NF_DOUBLE,nbdim,dimids,nvarid)
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_stats: Failed defining variable "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 3. Write attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",&
                     len_trim(adjustl(title)),adjustl(title))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_stats: Failed writing title attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",&
                     len_trim(adjustl(units)),adjustl(units))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var_stats: Failed writing units attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 4. Switch out of NetCDF define mode
ierr = NF_ENDDEF(nid)

end subroutine def_var_stats

