subroutine get_iddim(infid,latid,latlength,lonid,lonlength, &
                           altid,altlength,timid,timelength,lmdflag )

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer infid ! NetCDF input file ID
integer lonid,latid,altid,timid
integer lonlength ! # of grid points along longitude
integer latlength ! # of grid points along latitude
integer altlength ! # of grid point along altitude (of input datasets)
integer timelength ! # of points along time
logical lmdflag ! true=LMD, false=CAM

!local
integer tmpdimid ! temporarily store a dimension ID
integer ierr ! NetCDF routines return code

! latitude
ierr=NF_INQ_DIMID(infid,"latitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get latitude dimension ID"
  write(*,*) "  looking for lat dimension instead... "
  ierr=NF_INQ_DIMID(infid,"lat",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get lat dimension ID"
  else
    ierr=NF_INQ_VARID(infid,"lat",latid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get lat ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,latlength)
      if (ierr.ne.NF_NOERR) stop "Error: Failed to get lat length"
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"latitude",latid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get latitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,latlength)
    if (ierr.ne.NF_NOERR) stop "Error: Failed to get latitude length"
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
    ierr=NF_INQ_VARID(infid,"lon",lonid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get lon ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,lonlength)
      if (ierr.ne.NF_NOERR) stop "Error: Failed to get lon length"
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"longitude",lonid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get longitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,lonlength)
    if (ierr.ne.NF_NOERR) stop "Error: Failed to get longitude length"
  endif
endif

lmdflag=.true.
! altitude : pressure levels
ierr=NF_INQ_DIMID(infid,"altitude",tmpdimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get altitude dimension ID"
  write(*,*) "  looking for presnivs dimension instead... "
  ierr=NF_INQ_DIMID(infid,"presnivs",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Could not get presnivs dimension ID"
    write(*,*) "  looking for lev dimension instead... "
    ierr=NF_INQ_DIMID(infid,"lev",tmpdimid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get lev dimension ID"
    else
      ierr=NF_INQ_VARID(infid,"lev",altid)
      if (ierr.ne.NF_NOERR) then
        stop "Error: Failed to get lev ID"
      else
        ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
        if (ierr.ne.NF_NOERR) stop "Error: Failed to get lev length"
        lmdflag=.false.
      endif
    endif
  else
    ierr=NF_INQ_VARID(infid,"presnivs",altid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get presnivs ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
      if (ierr.ne.NF_NOERR) stop "Error: Failed to get presnivs length"
    endif
  endif
else
  ierr=NF_INQ_VARID(infid,"altitude",altid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get altitude ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,altlength)
    if (ierr.ne.NF_NOERR) stop "Error: Failed to get altitude length"
  endif
endif

! time
ierr=NF_INQ_DIMID(infid,"Time",tmpdimid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get Time dimension ID"
  write(*,*) "  looking for time dimension instead... "
  ierr=NF_INQ_DIMID(infid,"time",tmpdimid)
 if (ierr.ne.NF_NOERR) then
  write(*,*) "Could not get time dimension ID"
  write(*,*) "  looking for time_counter dimension instead... "
  ierr=NF_INQ_DIMID(infid,"time_counter",tmpdimid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get time_counter dimension ID"
  else
    ierr=NF_INQ_VARID(infid,"time_counter",timid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get time_counter ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
      if (ierr.ne.NF_NOERR) stop "Error: Failed to get time_counter length"
    endif
  endif
 else
    ierr=NF_INQ_VARID(infid,"time",timid)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to get time ID"
    else
      ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
      if (ierr.ne.NF_NOERR) stop "Error: Failed to get time length"
    endif
 endif
else
  ierr=NF_INQ_VARID(infid,"Time",timid)
  if (ierr.ne.NF_NOERR) then
    stop "Error: Failed to get Time ID"
  else
    ierr=NF_INQ_DIMLEN(infid,tmpdimid,timelength)
    if (ierr.ne.NF_NOERR) stop "Error: Failed to get Time length"
  endif
endif

return
end

!===========================================================================

subroutine get_var2d(infid,dim1,dim2,text,var,ierr1,ierr2)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: infid ! NetCDF input file ID
integer :: dim1,dim2 ! dim length of the 3D variable
character (len=64) :: text ! name of variable to read
real,dimension(dim1,dim2) :: var ! variable to read
integer :: ierr1,ierr2 ! NetCDF routines return code

! local
integer tmpvarid ! temporarily store a variable ID

ierr1=NF_INQ_VARID(infid,trim(text),tmpvarid)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Could not get ID for ",trim(text)
else
  write(*,*) "ID ok for ",trim(text)
  ierr2=NF_GET_VAR_REAL(infid,tmpvarid,var)
endif

return
end

!===========================================================================

subroutine get_var3d(infid,dim1,dim2,dim3,text,var,ierr1,ierr2)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: infid ! NetCDF input file ID
integer :: dim1,dim2,dim3 ! dim length of the 3D variable
character (len=64) :: text ! name of variable to read
real,dimension(dim1,dim2,dim3) :: var ! variable to read
integer :: ierr1,ierr2 ! NetCDF routines return code

! local
integer tmpvarid ! temporarily store a variable ID

ierr1=NF_INQ_VARID(infid,trim(text),tmpvarid)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Could not get ID for ",trim(text)
else
  write(*,*) "ID ok for ",trim(text)
  ierr2=NF_GET_VAR_REAL(infid,tmpvarid,var)
endif

return
end

!===========================================================================

subroutine get_var4d(infid,dim1,dim2,dim3,dim4,text,var,missing,ierr1,ierr2)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: infid ! NetCDF input file ID
integer :: dim1,dim2,dim3,dim4 ! dim length of the 4D variable
character (len=64) :: text ! name of variable to read
real,dimension(dim1,dim2,dim3,dim4) :: var ! variable to read
real :: missing ! missing value
integer :: ierr1,ierr2,miss ! NetCDF routines return code

! local
integer tmpvarid ! temporarily store a variable ID

ierr1=NF_INQ_VARID(infid,trim(text),tmpvarid)
if (ierr1.ne.NF_NOERR) then
  write(*,*) "Could not get ID for ",trim(text)
else
  write(*,*) "ID ok for ",trim(text)
  ierr2=NF_GET_VAR_REAL(infid,tmpvarid,var)
  miss=NF_GET_ATT_REAL(infid,tmpvarid,"missing_value",missing)
endif

return
end

!===========================================================================

subroutine write_dim(outfid,dim1,dim2,dim3,dim4,lon,lat,plev,time,&
                          lon_dimid,lat_dimid,alt_dimid,time_dimid)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: outfid ! NetCDF output file ID
integer :: dim1,dim2,dim3,dim4 ! dim length of the 4D variable
real,dimension(dim1) :: lon ! longitude
real,dimension(dim2) :: lat ! latitude
real,dimension(dim3) :: plev ! Pressure levels (Pa)
real,dimension(dim4) :: time ! time
integer lon_dimid,lat_dimid,alt_dimid,time_dimid ! NetCDF dimension IDs

! local
integer lon_varid,lat_varid,alt_varid,time_varid
character (len=64) :: text ! to store some text
integer ierr ! NetCDF routines return code


!------------
! longitude
!------------

ierr=NF_DEF_DIM(outfid,"longitude",dim1,lon_dimid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define longitude dimension"

ierr=NF_DEF_VAR(outfid,"longitude",NF_REAL,1,lon_dimid,lon_varid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define longitude variable"

! longitude attributes
text='east longitude'
ierr=NF_PUT_ATT_TEXT(outfid,lon_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing long_name for longitude"

text='degrees_east'
ierr=NF_PUT_ATT_TEXT(outfid,lon_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing units for longitude"

!------------
! latitude
!------------

ierr=NF_DEF_DIM(outfid,"latitude",dim2,lat_dimid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define latitude dimension"

ierr=NF_DEF_VAR(outfid,"latitude",NF_REAL,1,lat_dimid,lat_varid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define latitude variable"

! latitude attributes
text='north latitude'
ierr=NF_PUT_ATT_TEXT(outfid,lat_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing long_name for latitude"

text='degrees_north'
ierr=NF_PUT_ATT_TEXT(outfid,lat_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing units for latitude"

!------------
! pressure
!------------

ierr=NF_DEF_DIM(outfid,"presnivs",dim3,alt_dimid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define presnivs dimension"

ierr=NF_DEF_VAR(outfid,"presnivs",NF_REAL,1,alt_dimid,alt_varid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define presnivs variable"

!presnivs attributes
text='Pressure levels'
ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing long_name for presnivs (p levels)"

text='Pa'
ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing units for presnivs (p levels)"

text='down'
ierr=NF_PUT_ATT_TEXT(outfid,alt_varid,'positive',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing positive for presnivs (p levels)"

!------------
! time
!------------

ierr=NF_DEF_DIM(outfid,"Time",dim4,time_dimid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define time dimension"

ierr=NF_DEF_VAR(outfid,"Time",NF_REAL,1,time_dimid,time_varid)
if (ierr.ne.NF_NOERR) stop "Error: Could not define Time variable"

! time attributes
text='Time'
ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'long_name',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing long_name for Time"

text='days since 0000-01-1 00:00:00'
ierr=NF_PUT_ATT_TEXT(outfid,time_varid,'units',len_trim(text),text)
if (ierr.ne.NF_NOERR) stop "Error: Problem writing units for Time"

print*,"End of dim definitions"

!------------
! Switch out of NetCDF define mode
!------------
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch out of define mode"

! Write longitude
ierr=NF_PUT_VAR_REAL(outfid,lon_varid,lon)
if (ierr.ne.NF_NOERR) stop "Error: Could not write longitude data to output file"

! Write latitude
ierr=NF_PUT_VAR_REAL(outfid,lat_varid,lat)
if (ierr.ne.NF_NOERR) stop "Error: Could not write latitude data to output file"

! Write pressure
ierr=NF_PUT_VAR_REAL(outfid,alt_varid,plev)
if (ierr.ne.NF_NOERR) stop "Error: Could not write presnivs data to output file"

! Write time
ierr=NF_PUT_VAR_REAL(outfid,time_varid,time)
if (ierr.ne.NF_NOERR) stop "Error: Could not write Time data to output file"

print*,"Writing dim OK"

return
end


!===========================================================================

subroutine write_var1d(outfid,datashape,dim1,&
                       name,lgname,units,miss_val,var)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: outfid ! NetCDF output file ID
integer :: dim1 ! dim length of the 1D variable
integer :: datashape ! shape of datasets
character (len=10) :: name   ! name of variable
character (len=20) :: lgname ! long name of variable
character (len=10) :: units  ! unit of variable
real :: miss_val ! "missing value" to specify missing data
real,dimension(dim1) :: var ! variable to read

! local
character (len=64) :: text ! to store some text
integer ierr ! NetCDF routines return code
integer varid

!------------
! Define variable
!------------

ierr=NF_REDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch into define mode"

ierr=NF_DEF_VAR(outfid,trim(name),NF_REAL,1,datashape,varid)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not define variable : ",name
      stop
endif

! attributes
ierr=NF_PUT_ATT_TEXT(outfid,varid,'long_name',len_trim(lgname),lgname)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing long_name for ",name
      stop
endif

ierr=NF_PUT_ATT_TEXT(outfid,varid,'units',len_trim(units),units)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing units for ",name
      stop
endif

ierr=NF_PUT_ATT_REAL(outfid,varid,'missing_value',NF_REAL,1,miss_val)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: failed to write missing_value for ",name
      stop
endif

print*,"End of var def : ",name

!------------
! Switch out of NetCDF define mode
!------------
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch out of define mode"

ierr=NF_PUT_VAR_REAL(outfid,varid,var)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write var ",name
      stop
endif

print*,"Writing var OK : ",name

return
end

!===========================================================================

subroutine write_var2d(outfid,datashape,dim1,dim2,&
                       name,lgname,units,miss_val,var)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: outfid ! NetCDF output file ID
integer :: dim1,dim2 ! dim length of the 2D variable
integer,dimension(2) :: datashape ! shape of datasets
character (len=10) :: name   ! name of variable
character (len=20) :: lgname ! long name of variable
character (len=10) :: units  ! unit of variable
real :: miss_val ! "missing value" to specify missing data
real,dimension(dim1,dim2) :: var ! variable to read

! local
character (len=64) :: text ! to store some text
integer ierr ! NetCDF routines return code
integer varid

!------------
! Define variable
!------------

ierr=NF_REDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch into define mode"

ierr=NF_DEF_VAR(outfid,trim(name),NF_REAL,2,datashape,varid)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not define variable : ",name
      stop
endif

! attributes
ierr=NF_PUT_ATT_TEXT(outfid,varid,'long_name',len_trim(lgname),lgname)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing long_name for ",name
      stop
endif

ierr=NF_PUT_ATT_TEXT(outfid,varid,'units',len_trim(units),units)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing units for ",name
      stop
endif

ierr=NF_PUT_ATT_REAL(outfid,varid,'missing_value',NF_REAL,1,miss_val)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: failed to write missing_value for ",name
      stop
endif

print*,"End of var def : ",name

!------------
! Switch out of NetCDF define mode
!------------
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch out of define mode"

ierr=NF_PUT_VAR_REAL(outfid,varid,var)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write var ",name
      stop
endif

print*,"Writing var OK : ",name

return
end

!===========================================================================

subroutine write_var3d(outfid,datashape,dim1,dim2,dim3,&
                       name,lgname,units,miss_val,var)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: outfid ! NetCDF output file ID
integer :: dim1,dim2,dim3 ! dim length of the 3D variable
integer,dimension(3) :: datashape ! shape of datasets
character (len=10) :: name   ! name of variable
character (len=20) :: lgname ! long name of variable
character (len=10) :: units  ! unit of variable
real :: miss_val ! "missing value" to specify missing data
real,dimension(dim1,dim2,dim3) :: var ! variable to read

! local
character (len=64) :: text ! to store some text
integer ierr ! NetCDF routines return code
integer varid

!------------
! Define variable
!------------

ierr=NF_REDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch into define mode"

ierr=NF_DEF_VAR(outfid,trim(name),NF_REAL,3,datashape,varid)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not define variable : ",name
      stop
endif

! attributes
ierr=NF_PUT_ATT_TEXT(outfid,varid,'long_name',len_trim(lgname),lgname)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing long_name for ",name
      stop
endif

ierr=NF_PUT_ATT_TEXT(outfid,varid,'units',len_trim(units),units)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing units for ",name
      stop
endif

ierr=NF_PUT_ATT_REAL(outfid,varid,'missing_value',NF_REAL,1,miss_val)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: failed to write missing_value for ",name
      stop
endif

print*,"End of var def : ",name

!------------
! Switch out of NetCDF define mode
!------------
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch out of define mode"

ierr=NF_PUT_VAR_REAL(outfid,varid,var)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write var ",name
      stop
endif

print*,"Writing var OK : ",name

return
end

!===========================================================================

subroutine write_var4d(outfid,datashape,dim1,dim2,dim3,dim4,&
                       name,lgname,units,miss_val,var)

implicit none

include "netcdf.inc" ! NetCDF definitions

! arguments
integer :: outfid ! NetCDF output file ID
integer :: dim1,dim2,dim3,dim4 ! dim length of the 4D variable
integer,dimension(4) :: datashape ! shape of datasets
character (len=10) :: name   ! name of variable
character (len=20) :: lgname ! long name of variable
character (len=10) :: units  ! unit of variable
real :: miss_val ! "missing value" to specify missing data
real,dimension(dim1,dim2,dim3,dim4) :: var ! variable to read

! local
character (len=64) :: text ! to store some text
integer ierr ! NetCDF routines return code
integer varid

!------------
! Define variable
!------------

ierr=NF_REDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch into define mode"

ierr=NF_DEF_VAR(outfid,trim(name),NF_REAL,4,datashape,varid)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not define variable : ",name
      stop
endif

! attributes
ierr=NF_PUT_ATT_TEXT(outfid,varid,'long_name',len_trim(lgname),lgname)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing long_name for ",name
      stop
endif

ierr=NF_PUT_ATT_TEXT(outfid,varid,'units',len_trim(units),units)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Problem writing units for ",name
      stop
endif

ierr=NF_PUT_ATT_REAL(outfid,varid,'missing_value',NF_REAL,1,miss_val)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: failed to write missing_value for ",name
      stop
endif

print*,"End of var def : ",name

!------------
! Switch out of NetCDF define mode
!------------
ierr=NF_ENDDEF(outfid)
if (ierr.ne.NF_NOERR) stop "Error: Could not switch out of define mode"

ierr=NF_PUT_VAR_REAL(outfid,varid,var)
if (ierr.ne.NF_NOERR) then
      write(*,*) "Error: Could not write var ",name
      stop
endif

print*,"Writing var OK : ",name

return
end


