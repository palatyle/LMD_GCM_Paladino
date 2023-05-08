program concatnc


! ********************************************************
! Program to concatenate data from Netcdf "diagfi" files,
! outputs from the martian GCM
! input : diagfi.nc  / concat.nc / stats.nc kind of files
! author: Y. Wanherdrick
! + aps(), bps() and phisinit() are now also written
!   to output file (E. Millour, september 2006)
!   if available (F. Forget, october 2006)
! + handle 1D data (EM, January 2007)
! + ap(), bp()  (F. Forget, February 2008)
! + handle the possibility that number of GCM layers (aps, bps
!   or sigma) may be different from number of vertical levels
!   of data (which is the case for outputs from 'zrecast')
!   (EM, April 2010)
! + handle absence of ap() and bp() if aps and bps are available 
!    (case of stats file) FF, November 2011
! + read and write controle field, if available. TN, October 2013
! ********************************************************

implicit none

include "netcdf.inc" ! NetCDF definitions

character (len=80), dimension(1000) :: file
! file(): input file(s) names(s)
character (len=30), dimension(16) :: notconcat
! notconcat(): names of the (16) variables that won't be concatenated
character (len=50), dimension(:), allocatable :: var
! var(): name(s) of variable(s) that will be concatenated
character (len=50) :: tmpvar,tmpfile,title,units
! tmpvar(): used to temporarily store a variable name
! tmpfile(): used to temporarily store a file name
! title(): [netcdf] title attribute
! units(): [netcdf] units attribute
character (len=100) :: filename,vartmp
! filename(): output file name
! vartmp(): temporary variable name (read from netcdf input file)
!character (len=1) :: ccopy
! ccpy: 'y' or 'n' answer
character (len=4) :: axis
! axis: "ls" or "sol"
integer :: nid,ierr,miss
! nid: [netcdf] file ID #
! ierr: [netcdf] subroutine returned error code
! miss: [netcdf] subroutine returned error code
integer :: i,j,k,inter
! for various loops
integer :: varid
! varid: [netcdf] variable ID #
integer :: memolatlen=0,memolonlen=0,memoaltlen=0
! memolatlen: # of elements of lat(), read from the first input file
! memolonlen: # of elements of lon(), read from the first input file
! memoaltlen: # of elements of alt(), read from the first input file
real, dimension(:), allocatable:: lat,lon,alt,ctl,time
! lat(): array, stores latitude coordinates
! lon(): array, stores longitude coordinates
! alt(): array, stores altitude coordinates
! ctl(): array, stores controle coordinates
! time(): array, stores time coordinates
integer :: nbvar,nbfile,nbvarfile,ndim
! nbvar: # of variables to concatenate
! nbfile: # number of input file(s)
! nbvarfile: total # of variables in an input file
! ndim: [netcdf] # (3 or 4) of dimensions (for variables)
integer :: latdim,londim,altdim,ctldim,timedim
! latdim: [netcdf] "latitude" dim ID
! londim: [netcdf] "longitude" dim ID
! altdim: [netcdf] "altdim" dim ID
! ctldim: [netcdf] "ctldim" dim ID
! timedim: [netcdf] "timedim" dim ID
integer :: gcmlayerdim ! NetCDF dimension ID for # of layers in GCM
integer :: latvar,lonvar,altvar,ctlvar,timevar
! latvar: [netcdf] ID of "latitude" variable
! lonvar: [netcdf] ID of "longitude" variable
! altvar: [netcdf] ID of "altitude" variable
! ctlvar: [netcdf] ID of "controle" variable
! timevar: [netcdf] ID of "Time" variable
integer :: latlen,lonlen,altlen,ctllen,timelen
! latlen: # of elements of lat() array
! lonlen: # of elements of lon() array
! altlen: # of elements of alt() array
! ctllen: # of elements of ctl() array
! timelen: # of elemnets of time() array
integer :: GCM_layers ! number of GCM atmospheric layers (may not be
! same as altlen if file is an output of zrecast)
integer :: nout,latdimout,londimout,altdimout,timedimout,timevarout
! nout: [netcdf] output file ID
! latdimout: [netcdf] output latitude (dimension) ID
! londimout: [netcdf] output longitude (dimension) ID
! altdimout: [netcdf] output altitude (dimension) ID
! timedimout: [netcdf] output time (dimension) ID
! timevarout: [netcdf] ID of output "Time" variable
integer :: layerdimout ! NetCDF dimension ID for # of layers in GCM
integer :: interlayerdimout ! dimension ID for # of interlayers in GCM
integer :: reptime,rep,varidout
! reptime: total length of concatenated time() arrays
! rep: # number of elements of a time() array to write to the output file
! varidout: [netcdf] variable ID # (of a variable to write to the output file)
integer :: Nnotconcat,var_ok
! Nnotconcat: # of (leading)variables that won't be concatenated
! var_ok: flag (0 or 1)
integer, dimension(4) :: corner,edges,dim
! corner: [netcdf]
! edges: [netcdf]
! dim: [netcdf]
real, dimension(:,:,:,:), allocatable :: var3d
! var3D(,,,): 4D array to store a field
real :: memotime
! memotime: (cumulative) time value, in martian days (sols)
real :: missing
!PARAMETER(missing=1E+20)
! missing: [netcdf] to handle "missing" values when reading/writing files
real, dimension(2) :: valid_range
! valid_range(2): [netcdf] interval in which a value is considered valid

!==============================================================================
! 1.1. Get input file name(s)
!==============================================================================
write(*,*) 
write(*,*) "-- Program written specifically for NetCDF 'diagfi' files from the martian GCM --"
write(*,*) 
write(*,*) "which files do you want to use?"
write(*,*) "<Enter> when list ok"

nbfile=0
read(*,'(a50)') tmpfile
do While (len_trim(tmpfile).ne.0)
   nbfile=nbfile+1
   file(nbfile)=tmpfile
   read(*,'(a50)') tmpfile
enddo

if(nbfile==0) then
   write(*,*) "no file... game over"
   stop ""
endif

!==============================================================================
! 1.2. Ask for starting day value (memotime)
!==============================================================================

write(*,*)
!write(*,*) "Beginning day of the first specified file?"
write(*,*) "Starting day of the run stored in the first input file?"
write(*,*) " (Obsolete if the controle field is present, answer any number)"
write(*,*) "(e.g.: 100 if that run started at time=100 sols)"
read(*,*) memotime

!==============================================================================
! 1.3. Open the first input file
!==============================================================================

ierr = NF_OPEN(file(1),NF_NOWRITE,nid)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR: Pb opening file '//file(1)
   stop ""
endif

ierr=NF_INQ_NVARS(nid,nbvarfile)
! nbvarfile now set to be the (total) number of variables in file

!==============================================================================
! 1.4. Ask for (output) "Time" axis type
!==============================================================================

write(*,*) "Warning: to read the result with grads, choose sol"
write(*,*) "Warning: use ferret to read the non linear scale ls"
write(*,*) "Which time axis should be given in the output? (sol/ls)"
read(*,*) axis
! loop as long as axis is neither "sol" nor "ls"
do while ((axis/="sol").AND.(axis/="ls"))
   read(*,*) axis
enddo

! The following variables don't need to be concatenated
notconcat(1)='Time'
notconcat(2)='controle'
notconcat(3)='rlonu'
notconcat(4)='latitude'
notconcat(5)='longitude'
notconcat(6)='altitude'
notconcat(7)='rlatv'
notconcat(8)='aps'
notconcat(9)='bps'
notconcat(10)='ap'
notconcat(11)='bp'
notconcat(12)='cu'
notconcat(13)='cv'
notconcat(14)='aire'
notconcat(15)='phisinit'
notconcat(16)='soildepth'


!==============================================================================
! 1.5. Get (and check) list of variables to concatenate
!==============================================================================
write(*,*)
   Nnotconcat=0
do i=1,nbvarfile
   ierr=NF_INQ_VARNAME(nid,i,vartmp)
   ! vartmp now contains the "name" of variable of ID # i
   var_ok=0
   do inter=1,size(notconcat)
      if (vartmp.eq.notconcat(inter)) then
         var_ok=1 
         Nnotconcat=Nnotconcat+1
      endif 
   enddo     
   if (var_ok.eq.0)  write(*,*) trim(vartmp)
enddo

! Nnotconcat: # of variables that won't be concatenated
! nbvarfile: total # of variables in file
allocate(var(nbvarfile-Nnotconcat))


write(*,*)
write(*,*) "which variables do you want to concatenate?"
write(*,*) "all / list of <variables> (separated by <Enter>s)"
write(*,*) "(an empty line , i.e: just <Enter>, implies end of list)"
nbvar=0
read(*,'(a50)') tmpvar
do while ((tmpvar/=' ').AND.(trim(tmpvar)/='all'))
   nbvar=nbvar+1
   var(nbvar)=tmpvar
   read(*,'(a50)') tmpvar
enddo

if (tmpvar=="all") then
   if (axis=="ls") then
!      write(*,*) "Do you want to keep the original file? (y/n)"
!      read(*,*) ccopy
!      if ((ccopy=="n").or.(ccopy=="N")) then
!         do i=1,nbfile
!            ierr=NF_CLOSE(nid)
!            ierr = NF_OPEN(file(1),NF_WRITE,nid)
!            call change_time_axis(nid,ierr)
!            ierr=NF_CLOSE(nid)
!            STOP ""
!         enddo 
!      else
         nbvar=nbvarfile-Nnotconcat
         do j=Nnotconcat+1,nbvarfile
            ierr=nf_inq_varname(nid,j,var(j-Nnotconcat))
         enddo
!      endif ! of if ((ccopy=="n").or.(ccopy=="N"))
   endif ! of if (axis=="ls")
! Variables names from the file are catched
   nbvar=nbvarfile-Nnotconcat
   j=1
   do i=1,nbvarfile
      ierr=nf_inq_varname(nid,i,vartmp)
      var_ok=0
      do inter=1,size(notconcat)
        if (vartmp.eq.notconcat(inter)) then
           var_ok=1
        endif
      enddo
      if (var_ok.eq.0) then
         if (j .gt. nbvar) then
           write(*,*) "PROBLEM HERE !", var
           stop
         endif
         var(j) = vartmp 
         write(*,'(a9,1x,i2,1x,a1,1x,a50)') "variable ",j,":",var(j)
         j=j+1
      endif
   enddo
else if(nbvar==0) then
   write(*,*) "no variable... game over"
   stop ""
endif ! of if (tmpvar=="all")

! Name of the new file
!==========================================================
!filename=var(1)
!do i=2, nbvar
!   filename=trim(adjustl(filename))//"_"//var(i)
!enddo
!filename=trim(adjustl(filename))//".nc"

!==============================================================================
! 1.6. Get output file name
!==============================================================================
filename="concat.nc"


!==============================================================================
! 2. Concatenate input file(s) into output file
!==============================================================================

reptime=0

do i=1,nbfile

!==============================================================================
! 2.1. Open input file
!==============================================================================

   if (i/=1) then
      write(*,*) 
      write(*,*) "opening "//trim(file(i))//"..."
      ierr = NF_OPEN(file(i),NF_NOWRITE,nid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Pb opening file '//file(i)
         write(*,*) NF_STRERROR(ierr)
         stop ""
      endif
   endif

!==============================================================================
! 2.2. Read (and check) dimensions of variables from input file
!==============================================================================

   ierr=NF_INQ_DIMID(nid,"latitude",latdim)
   ierr=NF_INQ_VARID(nid,"latitude",latvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <latitude> is missing in file'//file(i)
      stop ""  
   endif
   ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
!  write(*,*) "latlen: ",latlen

   ierr=NF_INQ_DIMID(nid,"longitude",londim)
   ierr=NF_INQ_VARID(nid,"longitude",lonvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <longitude> is missing in file'//file(i)
      stop "" 
   endif
   ierr=NF_INQ_DIMLEN(nid,londim,lonlen)
!  write(*,*) "lonlen: ",lonlen

   ierr=NF_INQ_DIMID(nid,"altitude",altdim)
   ierr=NF_INQ_VARID(nid,"altitude",altvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <altitude> is missing in file'//file(i)
      stop ""
   endif
   ierr=NF_INQ_DIMLEN(nid,altdim,altlen)
!  write(*,*) "altlen: ",altlen

   ierr=NF_INQ_DIMID(nid,"index",ctldim)
   ierr=NF_INQ_VARID(nid,"controle",ctlvar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'Field <controle> is missing in file'//file
      ctllen=0
      !stop ""
   else
      ierr=NF_INQ_DIMLEN(nid,ctldim,ctllen)
   endif
!  write(*,*) "controle: ",controle

! load size of aps() or sigma() (in case it is not altlen)
   ! default is that GCM_layers=altlen
   ! but for outputs of zrecast, it may be a different value
   ierr=NF_INQ_DIMID(nid,"GCM_layers",gcmlayerdim)
   if (ierr.ne.NF_NOERR) then
     ! didn't find a GCM_layers dimension; therefore we have:
     GCM_layers=altlen
   else
     ! load value of GCM_layers
     ierr=NF_INQ_DIMLEN(nid,gcmlayerdim,GCM_layers)
   endif
!   write(*,*) "GCM_layers=",GCM_layers

!==============================================================================
! 2.3. Read (and check compatibility of) dimensions of
!       variables from input file
!==============================================================================

   if (i==1) then ! First call; initialize/allocate
      memolatlen=latlen
      memolonlen=lonlen
      memoaltlen=altlen
      allocate(lat(latlen))
      allocate(lon(lonlen))
      allocate(alt(altlen))
      allocate(ctl(ctllen))
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid,latvar,lat)
      ierr = NF_GET_VAR_DOUBLE(nid,lonvar,lon)
      ierr = NF_GET_VAR_DOUBLE(nid,altvar,alt)
      if (ctllen .ne. -1) ierr = NF_GET_VAR_DOUBLE(nid,ctlvar,ctl)
#else
      ierr = NF_GET_VAR_REAL(nid,latvar,lat)
      ierr = NF_GET_VAR_REAL(nid,lonvar,lon)
      ierr = NF_GET_VAR_REAL(nid,altvar,alt)
      if (ctllen .ne. -1) ierr = NF_GET_VAR_REAL(nid,ctlvar,ctl)
#endif
      if (ctllen .ne. -1) then
         if (modulo(int(memotime),669)/=modulo(int(ctl(4)),669)) then
           write(*,'(2(A,I4),A)') "WARNING: Starting day of the run is not ",&
                                modulo(int(memotime),669)," but ",modulo(int(ctl(4)),669),"!!"
           write(*,*) "Starting day of the run has been corrected."
           memotime=float(modulo(int(ctl(4)),669)) + ctl(27)
           ctl(4) = 0.
           ctl(27) = 0.
         endif
      endif
   ! Initialize output file's lat,lon,alt and time dimensions
      call initiate (filename,lat,lon,alt,ctl,GCM_layers,nout,&
       latdimout,londimout,altdimout,timedimout,&
       layerdimout,interlayerdimout,timevarout)
   ! Initialize output file's aps,bps,ap,bp and phisinit variables
     call init2(nid,lonlen,latlen,altlen,GCM_layers,&
                nout,londimout,latdimout,altdimout,&
                layerdimout,interlayerdimout)
   
   else ! Not a first call,
   ! Check that latitude,longitude and altitude of current input file
   ! are identical to those of the output file
      if (memolatlen/=latlen) then
           write(*,*) "ERROR: Not the same latitude axis"
           stop ""
        else if (memolonlen/=lonlen) then
           write(*,*) "ERROR: Not the same longitude axis"
           stop ""
             else if (memoaltlen/=altlen) then
                write(*,*) "ERROR: Not the same altitude axis"
                stop ""
       endif
   endif ! of if (i==1)

!==============================================================================
! 2.4. Handle "Time" dimension from input file
!==============================================================================

!==============================================================================
! 2.4.1 Read "Time" dimension from input file
!==============================================================================
   ierr=NF_INQ_DIMID(nid,"Time",timedim)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Dimension <Time> is missing in file'//file(i)
      stop ""
   endif
   ierr=NF_INQ_VARID(nid,"Time",timevar)
   if (ierr.NE.NF_NOERR) then
      write(*,*) 'ERROR: Field <Time> is missing in file'//file(i)
      stop ""
   endif
   ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
!  write(*,*) "timelen: ",timelen

   ! allocate time() array and fill it with values from input file
   allocate(time(timelen))

#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid,timevar,time)
#else
   ierr = NF_GET_VAR_REAL(nid,timevar,time)
#endif

!==============================================================================
! 2.4.2 Write/extend "Time" dimension/values in output file
!==============================================================================

   rep=0
   write(*,*)
   write(*,'(a3,1x,f6.1)') "Sol",memotime

   ! Add (memotime) offset and write "concatenated" time values to output file
   do k=reptime+1,reptime+timelen
      rep=rep+1
#ifdef NC_DOUBLE
      ierr= NF_PUT_VARA_DOUBLE(nout,timevarout,k,1,memotime+time(rep))
#else
      ierr= NF_PUT_VARA_REAL(nout,timevarout,k,1,memotime+time(rep))
#endif
   enddo
   ! Compute new time offset (for further concatenations)
   memotime=memotime+time(timelen)

!==============================================================================
! 2.5 Read/write variables
!==============================================================================

   do j=1,nbvar ! loop on variables to read/write

!==============================================================================
! 2.5.1 Check that variable to be read existe in input file
!==============================================================================

      write(*,*) "variable ",var(j)
      ierr=nf_inq_varid(nid,var(j),varid)
      if (ierr.NE.NF_NOERR) then
         write(*,*) 'ERROR: Field <',var(j),'> not found in file'//file(i)
         stop ""
      endif
      ierr=nf_inq_varndims(nid,varid,ndim)

!==============================================================================
! 2.5.2 Prepare things in order to read/write the variable
!==============================================================================

      ! build dim(),corner() and edges() arrays
      ! and allocate var3d() array
      if (ndim==1) then
         allocate(var3d(timelen,1,1,1))
         dim(1)=timedimout

         ! start indexes (where data values will be written)
         corner(1)=reptime+1
         corner(2)=1
         corner(3)=1
         corner(4)=1

	 ! length (along dimensions) of block of data to be written
         edges(1)=timelen 
         edges(2)=1 
         edges(3)=1
         edges(4)=1

      else if (ndim==3) then
         allocate(var3d(lonlen,latlen,timelen,1))
         dim(1)=londimout
         dim(2)=latdimout
         dim(3)=timedimout

         ! start indexes (where data values will be written)
         corner(1)=1
         corner(2)=1
         corner(3)=reptime+1
         corner(4)=1

	 ! length (along dimensions) of block of data to be written
         edges(1)=lonlen 
         edges(2)=latlen 
         edges(3)=timelen
         edges(4)=1
   
      else if (ndim==4) then
         allocate(var3d(lonlen,latlen,altlen,timelen))
         dim(1)=londimout
         dim(2)=latdimout
         dim(3)=altdimout
         dim(4)=timedimout

         ! start indexes (where data values will be written)
         corner(1)=1
         corner(2)=1
         corner(3)=1
         corner(4)=reptime+1

	 ! length (along dimensions) of block of data to be written
         edges(1)=lonlen
         edges(2)=latlen
         edges(3)=altlen
         edges(4)=timelen
      endif

      if (i==1) then ! First call: write some definitions to output file
         units="                                                    "
         title="                                                    "
         ierr=nf_get_att_text(nid,varid,"title",title)
         ierr=nf_get_att_text(nid,varid,"units",units)
         call def_var(nout,var(j),title,units,ndim,dim,varidout,ierr)
      else
         ierr=NF_INQ_VARID(nout,var(j),varidout)
      endif

!==============================================================================
! 2.5.3 Read from input file and write (append) to the output file
!==============================================================================

#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid,varid,var3d)
      ierr= NF_PUT_VARA_DOUBLE(nout,varidout,corner,edges,var3d)
      miss=NF_GET_ATT_DOUBLE(nid,varid,"missing_value",missing)
      miss=NF_GET_ATT_DOUBLE(nid,varid,"valid_range",valid_range)
#else
      ierr = NF_GET_VAR_REAL(nid,varid,var3d)
      ierr= NF_PUT_VARA_REAL(nout,varidout,corner,edges,var3d)
      miss=NF_GET_ATT_REAL(nid,varid,"missing_value",missing)
      miss=NF_GET_ATT_REAL(nid,varid,"valid_range",valid_range)
#endif

      if (ierr.ne.NF_NOERR) then
         write(*,*) 'PUT_VAR ERROR: ',NF_STRERROR(ierr)
         stop ""
      endif

! In case there is a "valid_range" attribute
      ! Write "valid_range" and "missing_value" attributes in output file
      if (miss.eq.NF_NOERR) then
         call missing_value(nout,varidout,valid_range,missing)
      endif

      ! free var3d() array
      deallocate(var3d)

   enddo ! of do j=1,nbvar

   ! Free time() and compute/store array length (for further concatenations)
   deallocate(time)
   reptime=reptime+timelen

   ! Close input file
   ierr=nf_close(nid)

enddo ! of i=1,nbfile

!==============================================================================
! 3. If required, change time axis (from sols to Ls)
!==============================================================================

if (axis=="ls") then
   call change_time_axis(nout,ierr)
endif

! Close output file
ierr=NF_CLOSE(nout)

contains

!******************************************************************************
subroutine initiate (filename,lat,lon,alt,ctl,GCM_layers,nout,&
         latdimout,londimout,altdimout,timedimout,&
         layerdimout,interlayerdimout,timevarout)
!==============================================================================
! Purpose:
! Create and initialize a data file (NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file (created in this subroutine) remains open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
character (len=*), intent(in):: filename
! filename(): the file's name
real, dimension(:), intent(in):: lat
! lat(): latitude
real, dimension(:), intent(in):: lon
! lon(): longitude
real, dimension(:), intent(in):: alt
! alt(): altitude
real, dimension(:), intent(in):: ctl
! ctl(): controle
integer,intent(in) :: GCM_layers ! number of GCM layers
integer, intent(out):: nout
! nout: [netcdf] file ID
integer, intent(out):: latdimout
! latdimout: [netcdf] lat() (i.e.: latitude)  ID
integer, intent(out):: londimout
! londimout: [netcdf] lon()  ID
integer, intent(out):: altdimout
! altdimout: [netcdf] alt()  ID
integer, intent(out):: timedimout
! timedimout: [netcdf] "Time"  ID
integer,intent(out) :: layerdimout
! layerdimout: [netcdf] "GCM_layers" ID
integer,intent(out) :: interlayerdimout
! layerdimout: [netcdf] "GCM_layers+1" ID
integer, intent(out):: timevarout
! timevarout: [netcdf] "Time" (considered as a variable) ID

!==============================================================================
! Local variables:
!==============================================================================
!integer :: latdim,londim,altdim,timedim
integer :: nvarid,ierr
integer :: ctldimout
! nvarid: [netcdf] ID of a variable
! ierr: [netcdf]  return error code (from called subroutines)

!==============================================================================
! 1. Create (and open) output file
!==============================================================================
write(*,*) "creating "//trim(adjustl(filename))//'...'
ierr = NF_CREATE(filename,IOR(NF_CLOBBER,NF_64BIT_OFFSET),nout)
! NB: setting NF_CLOBBER mode means that it's OK to overwrite an existing file
if (ierr.NE.NF_NOERR) then
   WRITE(*,*)'ERROR: Impossible to create the file.'
   stop ""
endif

!==============================================================================
! 2. Define/write "dimensions" and get their IDs
!==============================================================================

ierr = NF_DEF_DIM(nout, "latitude", size(lat), latdimout)
ierr = NF_DEF_DIM(nout, "longitude", size(lon), londimout)
ierr = NF_DEF_DIM(nout, "altitude", size(alt), altdimout)
if (size(ctl).ne.0) ierr = NF_DEF_DIM(nout, "index", size(ctl), ctldimout)
ierr = NF_DEF_DIM(nout, "Time", NF_UNLIMITED, timedimout)
ierr = NF_DEF_DIM(nout, "GCM_layers", GCM_layers, layerdimout)
ierr = NF_DEF_DIM(nout, "GCM_interlayers",GCM_layers+1,interlayerdimout)

! End netcdf define mode
ierr = NF_ENDDEF(nout)

!==============================================================================
! 3. Write "Time" and its attributes
!==============================================================================

call def_var(nout,"Time","Time","years since 0000-00-0 00:00:00",1,&
             (/timedimout/),timevarout,ierr)

!==============================================================================
! 4. Write "latitude" (data and attributes)
!==============================================================================

call def_var(nout,"latitude","latitude","degrees_north",1,&
             (/latdimout/),nvarid,ierr)

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,lat)
#else
ierr = NF_PUT_VAR_REAL (nout,nvarid,lat)
#endif

!==============================================================================
! 4. Write "longitude" (data and attributes)
!==============================================================================

call def_var(nout,"longitude","East longitude","degrees_east",1,&
             (/londimout/),nvarid,ierr)

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,lon)
#else
ierr = NF_PUT_VAR_REAL (nout,nvarid,lon)
#endif

!==============================================================================
! 5. Write "altitude" (data and attributes)
!==============================================================================

! Switch to netcdf define mode
ierr = NF_REDEF (nout)

#ifdef NC_DOUBLE
ierr = NF_DEF_VAR (nout,"altitude",NF_DOUBLE,1,altdimout,nvarid)
#else
ierr = NF_DEF_VAR (nout,"altitude",NF_FLOAT,1,altdimout,nvarid)
#endif

ierr = NF_PUT_ATT_TEXT (nout,nvarid,"long_name",8,"altitude")
ierr = NF_PUT_ATT_TEXT (nout,nvarid,'units',2,"km")
ierr = NF_PUT_ATT_TEXT (nout,nvarid,'positive',2,"up")

! End netcdf define mode
ierr = NF_ENDDEF(nout)

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,alt)
#else
ierr = NF_PUT_VAR_REAL (nout,nvarid,alt)
#endif

!==============================================================================
! 6. Write "controle" (data and attributes)
!==============================================================================

if (size(ctl).ne.0) then
   ! Switch to netcdf define mode
   ierr = NF_REDEF (nout)

#ifdef NC_DOUBLE
   ierr = NF_DEF_VAR (nout,"controle",NF_DOUBLE,1,ctldimout,nvarid)
#else
   ierr = NF_DEF_VAR (nout,"controle",NF_FLOAT,1,ctldimout,nvarid)
#endif

   ierr = NF_PUT_ATT_TEXT (nout,nvarid,"title",18,"Control parameters")

   ! End netcdf define mode
   ierr = NF_ENDDEF(nout)

#ifdef NC_DOUBLE
   ierr = NF_PUT_VAR_DOUBLE (nout,nvarid,ctl)
#else
   ierr = NF_PUT_VAR_REAL (nout,nvarid,ctl)
#endif
endif

end Subroutine initiate
!******************************************************************************
subroutine init2(infid,lonlen,latlen,altlen,GCM_layers, &
                 outfid,londimout,latdimout,altdimout, &
                 layerdimout,interlayerdimout)
!==============================================================================
! Purpose:
! Copy ap() , bp(), aps(), bps(), aire() and phisinit()
! from input file to outpout file
!==============================================================================
! Remarks:
! The NetCDF files must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: infid  ! NetCDF output file ID
integer, intent(in) :: lonlen ! # of grid points along longitude
integer, intent(in) :: latlen ! # of grid points along latitude
integer, intent(in) :: altlen ! # of grid points along latitude
integer, intent(in) :: GCM_layers ! # of GCM atmospheric layers
integer, intent(in) :: outfid ! NetCDF output file ID
integer, intent(in) :: londimout ! longitude dimension ID
integer, intent(in) :: latdimout ! latitude dimension ID
integer, intent(in) :: altdimout ! altitude dimension ID
integer, intent(in) :: layerdimout ! GCM_layers dimension ID
integer, intent(in) :: interlayerdimout ! GCM_layers+1 dimension ID
!==============================================================================
! Local variables:
!==============================================================================
real,dimension(:),allocatable :: aps,bps ! hybrid vertical coordinates
real,dimension(:),allocatable :: ap,bp ! hybrid vertical coordinates
real,dimension(:),allocatable :: sigma ! sigma levels
real,dimension(:,:),allocatable :: aire ! mesh areas
real,dimension(:,:),allocatable :: phisinit ! Ground geopotential
integer :: apsid,bpsid,sigmaid,phisinitid
integer :: apid,bpid
integer :: ierr
integer :: tmpvarid ! temporary variable ID
logical :: area ! is "aire" available ?
logical :: phis ! is "phisinit" available ?
logical :: hybrid ! are "aps" and "bps" available ?
logical :: apbp ! are "ap" and "bp" available ?

!==============================================================================
! 1. Read data from input file
!==============================================================================

! hybrid coordinate aps
!write(*,*) "aps: altlen=",altlen," GCM_layers=",GCM_layers
allocate(aps(GCM_layers),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate aps!"
  stop
endif
ierr=NF_INQ_VARID(infid,"aps",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*) "Ooops. Failed to get aps ID. OK, will look for sigma coord."
  hybrid=.false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,aps)
  hybrid=.true.
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading aps"
  endif

  ! hybrid coordinate bps
!  write(*,*) "bps: altlen=",altlen," GCM_layers=",GCM_layers
  allocate(bps(GCM_layers),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate bps!"
    stop
  endif
  ierr=NF_INQ_VARID(infid,"bps",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to get bps ID."
  endif
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,bps)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading bps"
  endif
endif

! hybrid coordinate ap
allocate(ap(GCM_layers+1),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate ap!"
  stop
else
  ierr=NF_INQ_VARID(infid,"ap",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Ooops. Failed to get ap ID. OK."
    apbp=.false.
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,ap)
    apbp=.true.
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading ap"
    endif
  endif
endif

! hybrid coordinate bp
allocate(bp(GCM_layers+1),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate bp!"
  stop
else
  ierr=NF_INQ_VARID(infid,"bp",tmpvarid)
  if (ierr.ne.NF_NOERR) then
    write(*,*) "Ooops. Failed to get bp ID. OK."
    apbp=.false.
  else
    ierr=NF_GET_VAR_REAL(infid,tmpvarid,bp)
    apbp=.true.
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed reading bp"
    endif
  endif
endif

! sigma levels (if any)
if (.not.hybrid) then
  allocate(sigma(GCM_layers),stat=ierr)
  if (ierr.ne.0) then
    write(*,*) "init2: failed to allocate sigma"
    stop
  endif
  ierr=NF_INQ_VARID(infid,"sigma",tmpvarid)
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,sigma)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading sigma"
  endif
endif ! of if (.not.hybrid)

! mesh area
allocate(aire(lonlen,latlen),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate aire!"
  stop
endif
ierr=NF_INQ_VARID(infid,"aire",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"init2 warning: Failed to get aire ID."
  area = .false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,aire)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading aire"
  endif
  area = .true.
endif

! ground geopotential phisinit
allocate(phisinit(lonlen,latlen),stat=ierr)
if (ierr.ne.0) then
  write(*,*) "init2: failed to allocate phisinit!"
  stop
endif
ierr=NF_INQ_VARID(infid,"phisinit",tmpvarid)
if (ierr.ne.NF_NOERR) then
  write(*,*)"init2 warning: Failed to get phisinit ID."
  phis = .false.
else
  ierr=NF_GET_VAR_REAL(infid,tmpvarid,phisinit)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed reading phisinit"
  endif
  phis = .true.
endif

!==============================================================================
! 2. Write
!==============================================================================

!==============================================================================
! 2.2. Hybrid coordinates ap() , bp(), aps() and bps()
!==============================================================================
if(hybrid) then 
! define aps
  call def_var(nout,"aps","hybrid pressure at midlayers"," ",1,&
             (/layerdimout/),apsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to def_var aps"
  endif

! write aps
#ifdef NC_DOUBLE
  ierr=NF_PUT_VAR_DOUBLE(outfid,apsid,aps)
#else
  ierr=NF_PUT_VAR_REAL(outfid,apsid,aps)
#endif
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to write aps"
  endif

! define bps
  call def_var(nout,"bps","hybrid sigma at midlayers"," ",1,&
             (/layerdimout/),bpsid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to def_var bps"
  endif

! write bps
#ifdef NC_DOUBLE
  ierr=NF_PUT_VAR_DOUBLE(outfid,bpsid,bps)
#else
  ierr=NF_PUT_VAR_REAL(outfid,bpsid,bps)
#endif
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to write bps"
  endif

  if (apbp) then
!   define ap
    call def_var(nout,"ap","hybrid sigma at interlayers"," ",1,&
             (/interlayerdimout/),apid,ierr)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to def_var ap"
    endif

! write ap
#ifdef NC_DOUBLE
    ierr=NF_PUT_VAR_DOUBLE(outfid,apid,ap)
#else
    ierr=NF_PUT_VAR_REAL(outfid,apid,ap)
#endif
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to write ap"
    endif

! define bp
    call def_var(nout,"bp","hybrid sigma at interlayers"," ",1,&
             (/interlayerdimout/),bpid,ierr)
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to def_var bp"
    endif

!   write bp
#ifdef NC_DOUBLE
    ierr=NF_PUT_VAR_DOUBLE(outfid,bpid,bp)
#else
    ierr=NF_PUT_VAR_REAL(outfid,bpid,bp)
#endif
    if (ierr.ne.NF_NOERR) then
      stop "Error: Failed to write bp"
    endif
  endif ! of if (apbp)

else
! define sigma
  call def_var(nout,"sigma","sigma at midlayers"," ",1,&
             (/layerdimout/),sigmaid,ierr)
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to def_var sigma"
  endif
! write sigma
#ifdef NC_DOUBLE
  ierr=NF_PUT_VAR_DOUBLE(outfid,sigmaid,sigma)
#else
  ierr=NF_PUT_VAR_REAL(outfid,sigmaid,sigma)
#endif
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to write sigma"
  endif
endif ! of if (hybrid)

!==============================================================================
! 2.2. aire() and phisinit()
!==============================================================================

if (area) then
  ! define aire
  call def_var(nout,"aire","Mesh area","m2",2,&
           (/londimout,latdimout/),tmpvarid,ierr)
  if (ierr.ne.NF_NOERR) then
     stop "init2 Error: Failed to def_var aire"
  endif
  
  ! write aire
#ifdef NC_DOUBLE
  ierr=NF_PUT_VAR_DOUBLE(outfid,tmpvarid,aire)
#else
  ierr=NF_PUT_VAR_REAL(outfid,tmpvarid,aire)
#endif
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to write aire"
  endif
endif ! of if (area)

IF (phis) THEN

  !define phisinit
   call def_var(nout,"phisinit","Ground level geopotential"," ",2,&
            (/londimout,latdimout/),phisinitid,ierr)
   if (ierr.ne.NF_NOERR) then
     stop "init2 Error: Failed to def_var phisinit"
   endif

  ! write phisinit
#ifdef NC_DOUBLE
  ierr=NF_PUT_VAR_DOUBLE(outfid,phisinitid,phisinit)
#else
  ierr=NF_PUT_VAR_REAL(outfid,phisinitid,phisinit)
#endif
  if (ierr.ne.NF_NOERR) then
    stop "init2 Error: Failed to write phisinit"
  endif

ENDIF ! of IF (phis)


! Cleanup
if (allocated(aps)) deallocate(aps)
if (allocated(bps)) deallocate(bps)
if (allocated(ap)) deallocate(ap)
if (allocated(bp)) deallocate(bp)
if (allocated(sigma)) deallocate(sigma)
if (allocated(phisinit)) deallocate(phisinit)
if (allocated(aire)) deallocate(aire)

end subroutine init2
!******************************************************************************
subroutine def_var(nid,name,title,units,nbdim,dim,nvarid,ierr)
!==============================================================================
! Purpose: Write a variable (i.e: add a variable to a dataset)
! called "name"; along with its attributes "title", "units"...
! to a file (following the NetCDF format)
!==============================================================================
! Remarks:
! The NetCDF file must be open
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nid
! nid: [netcdf] file ID #
character (len=*), intent(in) :: name
! name(): [netcdf] variable's name
character (len=*), intent(in) :: title
! title(): [netcdf] variable's "title" attribute
character (len=*), intent(in) :: units
! unit(): [netcdf] variable's "units" attribute
integer, intent(in) :: nbdim
! nbdim: number of dimensions of the variable
integer, dimension(nbdim), intent(in) :: dim
! dim(nbdim): [netcdf] dimension(s) ID(s)
integer, intent(out) :: nvarid
! nvarid: [netcdf] ID # of the variable
integer, intent(out) :: ierr
! ierr: [netcdf] subroutines returned error code

! Switch to netcdf define mode
ierr=NF_REDEF(nid)

! Insert the definition of the variable
#ifdef NC_DOUBLE
ierr=NF_DEF_VAR(nid,adjustl(name),NF_DOUBLE,nbdim,dim,nvarid)
#else
ierr=NF_DEF_VAR(nid,adjustl(name),NF_FLOAT,nbdim,dim,nvarid)
#endif

! Write the attributes
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"title",len_trim(adjustl(title)),adjustl(title))
ierr=NF_PUT_ATT_TEXT(nid,nvarid,"units",len_trim(adjustl(units)),adjustl(units))

! End netcdf define mode
ierr=NF_ENDDEF(nid)

end subroutine def_var
!******************************************************************************
subroutine change_time_axis(nid,ierr)
!==============================================================================
! Purpose: 
! Read "time" variable from a dataset, convert it from "sol" (martian
! days) to "Ls" (solar longitude) and write the result back in the file
!==============================================================================
! Remarks:
! The NetCDF file must be opened before this subroutine is called
!==============================================================================

implicit none

include "netcdf.inc" ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nid
! nid: [netcdf] file ID
integer, intent(out) :: ierr
! ierr: [netcdf]  return error code

!==============================================================================
! Local variables:
!==============================================================================
integer :: nvarid
! nvarid: ID of the "Time" variable
integer :: timelen
! timelen: size of the arrays
integer :: timedim
! timedim: ID of the "Time" dimension
integer i

real, dimension(:), allocatable :: time,ls
! time(): time, given in sols
! ls(): time, given in Ls (solar longitude)

!==============================================================================
! 1. Read
!==============================================================================

ierr=NF_INQ_DIMID(nid,"Time",timedim)
ierr=NF_INQ_VARID(nid,"Time",nvarid)
if (ierr.NE.NF_NOERR) then
   write(*,*) 'ERROR in change_time_axis: Field <Time> not found'
   print*, NF_STRERROR(ierr)
   stop ""
endif

ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
allocate(time(timelen),ls(timelen))
#ifdef NC_DOUBLE
ierr = NF_GET_VAR_DOUBLE(nid,nvarid,time)
#else
ierr = NF_GET_VAR_REAL(nid,nvarid,time)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*) "ERROR in change_time_axis: Failed to load Time"
  stop
endif

!==============================================================================
! 2. Convert sols to Ls
!==============================================================================

do i=1,timelen
   call sol2ls(time(i),ls(i))
enddo

!==============================================================================
! 2.1 Check if there are not jumps in Ls (rounding problems in sol2ls)
!==============================================================================

do i=1,timelen-1
   if ((ls(i+1)-ls(i)) > 350) then
       write(*,*) "+ 360 deg Ls jump solved:", ls(i), ls(i+1), "at timestep", i
      ls(i+1) = ls(i+1) - 360
       write(*,*) " corrected to now be:", ls(i), ls(i+1)
   else if ((ls(i)-ls(i+1)) > 350) then
       write(*,*) "- 360 deg Ls jump solved:", ls(i), ls(i+1), "at timestep", i
      ls(i+1) = ls(i+1) + 360
       write(*,*) " corrected to now be:", ls(i), ls(i+1)
   endif
enddo

!==============================================================================
! 3. Write
!==============================================================================

#ifdef NC_DOUBLE
ierr = NF_PUT_VAR_DOUBLE(nid,nvarid,ls)
#else
ierr = NF_PUT_VAR_REAL(nid,nvarid,ls)
#endif
if (ierr.ne.NF_NOERR) then
  write(*,*) "ERROR in change_time_axis: Failed to write Ls"
  stop
endif

end subroutine change_time_axis
!******************************************************************************
subroutine sol2ls(sol,Ls)
!==============================================================================
! Purpose: 
! Convert a date/time, given in sol (martian day),
! into solar longitude date/time, in Ls (in degrees),
! where sol=0 is (by definition) the northern hemisphere
!  spring equinox (where Ls=0).
!==============================================================================
! Notes:
! Even though "Ls" is cyclic, if "sol" is greater than N (martian) year,
! "Ls" will be increased by N*360
! Won't work as expected if sol is negative (then again,
! why would that ever happen?)
!==============================================================================

implicit none

!==============================================================================
! Arguments:
!==============================================================================
real,intent(in) :: sol
real,intent(out) :: Ls

!==============================================================================
! Local variables:
!==============================================================================
real year_day,peri_day,timeperi,e_elips,twopi,degrad
data year_day /669./            ! # of sols in a martian year
data peri_day /485.0/           
data timeperi /1.9082314/
data e_elips  /0.093358/
data twopi       /6.2831853/    ! 2.*pi
data degrad   /57.2957795/      ! pi/180

real zanom,xref,zx0,zdx,zteta,zz

integer count_years
integer iter

!==============================================================================
! 1. Compute Ls
!==============================================================================

zz=(sol-peri_day)/year_day
zanom=twopi*(zz-nint(zz))
xref=abs(zanom)

!  The equation zx0 - e * sin (zx0) = xref, solved by Newton
zx0=xref+e_elips*sin(xref)
do iter=1,20 ! typically, 2 or 3 iterations are enough
   zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
   zx0=zx0+zdx
   if(abs(zdx).le.(1.e-7)) then
!      write(*,*)'iter:',iter,'     |zdx|:',abs(zdx)
      exit
   endif
enddo

if(zanom.lt.0.) zx0=-zx0

zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
Ls=zteta-timeperi

if(Ls.lt.0.) then
   Ls=Ls+twopi
else
   if(Ls.gt.twopi) then
      Ls=Ls-twopi
   endif
endif

Ls=degrad*Ls
! Ls is now in degrees

!==============================================================================
! 1. Account for (eventual) years included in input date/time sol
!==============================================================================

count_years=0 ! initialize
zz=sol  ! use "zz" to store (and work on) the value of sol
do while (zz.ge.year_day)
   count_years=count_years+1
   zz=zz-year_day
enddo

! Add 360 degrees to Ls for every year
if (count_years.ne.0) then
   Ls=Ls+360.*count_years
endif


end subroutine sol2ls
!******************************************************************************
subroutine  missing_value(nout,nvarid,valid_range,missing)
!==============================================================================
! Purpose:
! Write "valid_range" and "missing_value" attributes (of a given
! variable) to a netcdf file
!==============================================================================
! Remarks:
! NetCDF file must be open
! Variable (nvarid) ID must be set
!==============================================================================

implicit none

include "netcdf.inc"  ! NetCDF definitions

!==============================================================================
! Arguments:
!==============================================================================
integer, intent(in) :: nout
! nout: [netcdf] file ID #
integer, intent(in) :: nvarid
! varid: [netcdf] variable ID #
real, dimension(2), intent(in) :: valid_range
! valid_range(2): [netcdf] "valid_range" attribute (min and max)
real, intent(in) :: missing
! missing: [netcdf] "missing_value" attribute

!==============================================================================
! Local variables:
!==============================================================================
integer :: ierr
! ierr: [netcdf] subroutine returned error code
!      INTEGER nout,nvarid,ierr
!      REAL missing
!      REAL valid_range(2)

! Switch to netcdf dataset define mode
ierr = NF_REDEF (nout)
if (ierr.ne.NF_NOERR) then
   print*,'missing_value: '
   print*, NF_STRERROR(ierr)
endif

! Write valid_range() attribute
#ifdef NC_DOUBLE
ierr=NF_PUT_ATT_DOUBLE(nout,nvarid,'valid_range',NF_DOUBLE,2,valid_range)
#else
ierr=NF_PUT_ATT_REAL(nout,nvarid,'valid_range',NF_FLOAT,2,valid_range)
#endif

if (ierr.ne.NF_NOERR) then
   print*,'missing_value: valid_range attribution failed'
   print*, NF_STRERROR(ierr)
   !write(*,*) 'NF_NOERR', NF_NOERR
   !CALL abort
   stop ""
endif

! Write "missing_value" attribute
#ifdef NC_DOUBLE
ierr= NF_PUT_ATT_DOUBLE(nout,nvarid,'missing_value',NF_DOUBLE,1,missing)
#else
ierr= NF_PUT_ATT_REAL(nout,nvarid,'missing_value',NF_FLOAT,1,missing)
#endif

if (ierr.NE.NF_NOERR) then
   print*, 'missing_value: missing value attribution failed'
   print*, NF_STRERROR(ierr)
!    WRITE(*,*) 'NF_NOERR', NF_NOERR
!    CALL abort
   stop ""
endif

! End netcdf dataset define mode
ierr = NF_ENDDEF(nout)

end subroutine  missing_value
!******************************************************************************

end program concatnc
