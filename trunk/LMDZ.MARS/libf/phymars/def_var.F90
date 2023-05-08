subroutine def_var(nid,name,title,units,nbdim,dimids,nvarid,ierr)

! This subroutine defines variable 'name' in a (pre-existing and opened)
! NetCDF file (known from its NetCDF ID 'nid').
! The number of dimensions 'nbdim' of the variable, as well as the IDs of
! corresponding dimensions must be set (in array 'dimids').
! Upon successfull definition of the variable, 'nvarid' contains the
! NetCDF ID of the variable.
! The variables' attributes 'title' (Note that 'long_name' would be more
! appropriate) and 'units' are also set. 
! Modifs: Aug2010 Ehouarn: enforce outputs to be real*4

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
!#ifdef NC_DOUBLE
!ierr = NF_DEF_VAR (nid,adjustl(name),NF_DOUBLE,nbdim,dimids,nvarid)
!#else
ierr = NF_DEF_VAR (nid,adjustl(name),NF_FLOAT,nbdim,dimids,nvarid)
!#endif
if(ierr/=NF_NOERR) then
   write(*,*) "def_var: Failed defining variable "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 3. Write attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",&
                     len_trim(adjustl(title)),adjustl(title))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var: Failed writing title attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",&
                     len_trim(adjustl(units)),adjustl(units))
if(ierr/=NF_NOERR) then
   write(*,*) "def_var: Failed writing units attribute for "//trim(name)
   write(*,*) NF_STRERROR(ierr)
   stop ""
endif

! 4. Switch out of NetCDF define mode
ierr = NF_ENDDEF(nid)

end
