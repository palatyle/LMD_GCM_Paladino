subroutine def_var(nid,name,title,units,nbdim,dim,nvarid,ierr)

implicit none

include "netcdf.inc"

character (len=*) :: title,units,name
integer :: nid,nbdim,nvarid,ierr
integer, dimension(nbdim) :: dim

ierr=NF_REDEF(nid)
#ifdef NC_DOUBLE
ierr = NF_DEF_VAR (nid,adjustl(name),NF_DOUBLE,nbdim,dim,nvarid)
#else
ierr = NF_DEF_VAR (nid,adjustl(name),NF_FLOAT,nbdim,dim,nvarid)
#endif
if(ierr/=NF_NOERR) then
   write(*,*) NF_STRERROR(ierr)
   stop "in def_var"
endif
ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", len_trim(adjustl(title)),adjustl(title))
if(ierr/=NF_NOERR) then
   write(*,*) NF_STRERROR(ierr)
   stop "in def_var"
endif
ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", len_trim(adjustl(units)),adjustl(units))
if(ierr/=NF_NOERR) then
   write(*,*) NF_STRERROR(ierr)
   stop "in def_var"
endif
ierr = NF_ENDDEF(nid)

end
