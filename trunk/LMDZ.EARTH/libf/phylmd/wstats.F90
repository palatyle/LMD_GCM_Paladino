subroutine wstats(ngrid,nom,titre,unite,dim,px)

implicit none

#include "dimensions.h"
#include "statto.h"
#include "netcdf.inc"

integer,intent(in) :: ngrid
character (len=*),intent(in) :: nom,titre,unite
integer,intent(in) :: dim
real, dimension(ngrid,llm),intent(in) :: px
integer,parameter :: iip1=iim+1
integer,parameter :: jjp1=jjm+1
real, dimension(iip1,jjp1,llm) :: mean3d,sd3d,dx3
real, dimension(iip1,jjp1) :: mean2d,sd2d,dx2
character (len=50) :: namebis
character (len=50), save :: firstvar
integer :: ierr,varid,nbdim,nid
integer :: meanid,sdid
integer, dimension(4)  :: id,start,size
logical, save :: firstcall=.TRUE.
integer :: l,i,j,ig0
integer,save :: index

integer, save :: step=0


if (firstcall) then
   firstcall=.false.
   firstvar=trim((nom))
   call inistats(ierr)
endif

if (firstvar==nom) then ! If we're back to the first variable
      step = step + 1
endif

if (mod(step,istats).ne.0) then
   RETURN
endif

ierr = NF_OPEN("stats.nc",NF_WRITE,nid)

namebis=trim(nom)
ierr= NF_INQ_VARID(nid,namebis,meanid)

if (ierr.ne.NF_NOERR) then

   if (firstvar==nom) then 
      index=1
      count=0
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
   call inivar(nid,meanid,ngrid,dim,index,px,ierr)
   namebis=trim(nom)//"_sd"
   call def_var_stats(nid,namebis,trim(titre)//" total standard deviation over the season",unite,nbdim,id,sdid,ierr)
   call inivar(nid,sdid,ngrid,dim,index,px,ierr)

   ierr= NF_CLOSE(nid)
   return

else
   namebis=trim(nom)//"_sd"
   ierr= NF_INQ_VARID(nid,namebis,sdid)

endif

if (firstvar==nom) then 
   count(index)=count(int(index))+1
   index=index+1
   if (index>istime) then
      index=1
   endif
endif

if (count(index)==0) then
   if (dim.eq.3) then
      start=(/1,1,1,index/)
      size=(/iip1,jjp1,llm,1/)
      mean3d=0
      sd3d=0
   else if (dim.eq.2) then
      start=(/1,1,index,0/)
      size=(/iip1,jjp1,1,0/)
      mean2d=0
      sd2d=0
   endif
else
   if (dim.eq.3) then
      start=(/1,1,1,index/)
      size=(/iip1,jjp1,llm,1/)
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid,meanid,start,size,mean3d)
      ierr = NF_GET_VARA_DOUBLE(nid,sdid,start,size,sd3d)
#else
      ierr = NF_GET_VARA_REAL(nid,meanid,start,size,mean3d)
      ierr = NF_GET_VARA_REAL(nid,sdid,start,size,sd3d)
#endif
      if (ierr.ne.NF_NOERR) then
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif

   else if (dim.eq.2) then
      start=(/1,1,index,0/)
      size=(/iip1,jjp1,1,0/)
#ifdef NC_DOUBLE
      ierr = NF_GET_VARA_DOUBLE(nid,meanid,start,size,mean2d)
      ierr = NF_GET_VARA_DOUBLE(nid,sdid,start,size,sd2d)
#else
      ierr = NF_GET_VARA_REAL(nid,meanid,start,size,mean2d)
      ierr = NF_GET_VARA_REAL(nid,sdid,start,size,sd2d)
#endif
      if (ierr.ne.NF_NOERR) then
         write (*,*) NF_STRERROR(ierr)
         stop ""
      endif
   endif
endif

if (dim.eq.3) then

!  Passage variable physique -->  variable dynamique

   DO l=1,llm
      DO i=1,iip1
         dx3(i,1,l)=px(1,l)
         dx3(i,jjp1,l)=px(ngrid,l)
      ENDDO
      DO j=2,jjm
         ig0= 1+(j-2)*iim
         DO i=1,iim
            dx3(i,j,l)=px(ig0+i,l)
         ENDDO
         dx3(iip1,j,l)=dx3(1,j,l)
      ENDDO
   ENDDO

   mean3d(:,:,:)= mean3d(:,:,:)+dx3(:,:,:)
   sd3d(:,:,:)= sd3d(:,:,:)+dx3(:,:,:)**2

#ifdef NC_DOUBLE
   ierr = NF_PUT_VARA_DOUBLE(nid,meanid,start,size,mean3d)
   ierr = NF_PUT_VARA_DOUBLE(nid,sdid,start,size,sd3d)
#else
   ierr = NF_PUT_VARA_REAL(nid,meanid,start,size,mean3d)
   ierr = NF_PUT_VARA_REAL(nid,sdid,start,size,sd3d)
#endif
   if (ierr.ne.NF_NOERR) then
     write (*,*) NF_STRERROR(ierr)
     stop ""
   endif

else if (dim.eq.2) then

!    Passage variable physique -->  physique dynamique

  DO i=1,iip1
     dx2(i,1)=px(1,1)
     dx2(i,jjp1)=px(ngrid,1)
  ENDDO
  DO j=2,jjm
     ig0= 1+(j-2)*iim
     DO i=1,iim
        dx2(i,j)=px(ig0+i,1)
     ENDDO
     dx2(iip1,j)=dx2(1,j)
  ENDDO

   mean2d(:,:)= mean2d(:,:)+dx2(:,:)
   sd2d(:,:)= sd2d(:,:)+dx2(:,:)**2

#ifdef NC_DOUBLE
   ierr = NF_PUT_VARA_DOUBLE(nid,meanid,start,size,mean2d)
   ierr = NF_PUT_VARA_DOUBLE(nid,sdid,start,size,sd2d)
#else
   ierr = NF_PUT_VARA_REAL(nid,meanid,start,size,mean2d)
   ierr = NF_PUT_VARA_REAL(nid,sdid,start,size,sd2d)
#endif
   if (ierr.ne.NF_NOERR) then
     write (*,*) NF_STRERROR(ierr)
     stop ""
   endif

endif

ierr= NF_CLOSE(nid)

end subroutine wstats

!======================================================
subroutine inivar(nid,varid,ngrid,dim,index,px,ierr)

implicit none

include "dimensions.h"
!!!!!!include "dimphys.h"
include "netcdf.inc"

integer, intent(in) :: nid,varid,dim,index,ngrid
real, dimension(ngrid,llm), intent(in) :: px
integer, intent(out) :: ierr

integer,parameter :: iip1=iim+1
integer,parameter :: jjp1=jjm+1

integer :: l,i,j,ig0
integer, dimension(4) :: start,size
real, dimension(iip1,jjp1,llm) :: dx3
real, dimension(iip1,jjp1) :: dx2

if (dim.eq.3) then

   start=(/1,1,1,index/)
   size=(/iip1,jjp1,llm,1/)

!  Passage variable physique -->  variable dynamique

   DO l=1,llm
      DO i=1,iip1
         dx3(i,1,l)=px(1,l)
         dx3(i,jjp1,l)=px(ngrid,l)
      ENDDO
      DO j=2,jjm
         ig0= 1+(j-2)*iim
         DO i=1,iim
            dx3(i,j,l)=px(ig0+i,l)
         ENDDO
         dx3(iip1,j,l)=dx3(1,j,l)
      ENDDO
   ENDDO

#ifdef NC_DOUBLE
   ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,size,dx3)
#else
   ierr = NF_PUT_VARA_REAL(nid,varid,start,size,dx3)
#endif

else if (dim.eq.2) then

      start=(/1,1,index,0/)
      size=(/iip1,jjp1,1,0/)

!    Passage variable physique -->  physique dynamique

  DO i=1,iip1
     dx2(i,1)=px(1,1)
     dx2(i,jjp1)=px(ngrid,1)
  ENDDO
  DO j=2,jjm
     ig0= 1+(j-2)*iim
     DO i=1,iim
        dx2(i,j)=px(ig0+i,1)
     ENDDO
     dx2(iip1,j)=dx2(1,j)
  ENDDO

#ifdef NC_DOUBLE
   ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,size,dx2)
#else
   ierr = NF_PUT_VARA_REAL(nid,varid,start,size,dx2)
#endif

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

#include "netcdf.inc"

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
#ifdef NC_DOUBLE
ierr = NF_DEF_VAR (nid,adjustl(name),NF_DOUBLE,nbdim,dimids,nvarid)
#else
ierr = NF_DEF_VAR (nid,adjustl(name),NF_FLOAT,nbdim,dimids,nvarid)
#endif
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

