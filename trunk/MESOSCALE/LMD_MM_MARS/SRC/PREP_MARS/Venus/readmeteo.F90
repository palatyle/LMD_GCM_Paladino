program readmeteo

implicit none
include "netcdf.inc"

!------------------------------------------------------------------------!
! Readmeteo prepares an initial state for the WPS pre-processor of WRF   !
!                                                                        !
! Input is a diagfi.nc NETCDF file from a LMD/GCM simulation             !
!                                                                        !
! Outputs are WRFSI intermediate format files ready for metgrid.exe      !
!       http://www.mmm.ucar.edu/wrf/OnLineTutorial/WPS/IM_si.htm         !
!                                                                        !
! **** Please create a WPSFEED folder (or a symbolic link) ****          !       
!                                                                        !
! A. Spiga - 16/04/2007                                                  !
!	     22/05/2007 : need to run zrecast as a preliminary           !
!            07/06/2007 : external parameter file 'readmeteo.def         !
!	     30/07/2007 : manage additional variables (tsoil, emiss,...) !	
!            19/10/2007 : no need to run zrecast at all (eta levels)     !   
!		12/2007 : include co2ice and emissivity			 ! 
!               02/2008 : include water vapor and ice                    !
!               01/2010 : possible use of diagsoil for new physics       ! 
!               02/2011 : add dust tracers, correct surface              !
!									 !
! spiga@lmd.jussieu.fr							 !	
!------------------------------------------------------------------------!



!***************************************************************************
!***************************************************************************
REAL, PARAMETER :: emiss_prescribed=0.95
REAL, PARAMETER :: co2ice_prescribed=0.
CHARACTER (LEN=3), PARAMETER :: ident='LMD'



!***************************************************************************
!***************************************************************************

REAL :: ptop
REAL, PARAMETER :: grav=8.89
LOGICAL, PARAMETER :: blank=.false.



!! Variables to be written in the output file
!! *** NB: conformity with metgrid/src/read_met_module.F90
CHARACTER (LEN=33) :: output
INTEGER :: IFV
CHARACTER (LEN=24) :: HDATE
REAL :: XFCST
CHARACTER (LEN=32) :: SOURCE
CHARACTER (LEN=9) :: FIELD
CHARACTER (LEN=25) :: UNITS
CHARACTER (LEN=46) :: DESC
REAL :: XLVL
INTEGER :: NX,NY,IPROJ
CHARACTER (LEN=8) :: STARTLOC
REAL :: STARTLAT,STARTLON,DELTALAT,DELTALON
REAL, POINTER, DIMENSION(:,:) :: SLAB

!! Variables related to NETCDF format
integer :: nid,nid2,nid3,nvarid,ierr,ierr2
integer :: timedim,londim,latdim,altdim,altdim2
integer :: rlatvdim,rlonudim
integer :: timelen,lonlen,latlen,altlen,altlen2
integer :: rlatvlen,rlonulen
integer :: soillen,soildim
integer :: nphys
integer, dimension(4) :: corner,edges

!! Intermediate data arrays
integer :: k,l,m,n,p,i,j
real, dimension(:), allocatable :: lat,lon,time,alt,aps,bps,levels,vertdsoil
real, dimension(:,:), allocatable :: vide,ones,ghtsfile
real, dimension(:,:), allocatable :: interm
real, dimension(:,:,:), allocatable :: gwparam
real, dimension(:,:,:), allocatable :: psfile,tsfile
real, dimension(:,:,:), allocatable :: emissfile,co2icefile
real, dimension(:,:,:), allocatable :: tnsfile,unsfile,vnsfile
real, dimension(:,:,:,:), allocatable :: tfile,ufile,vfile,rfile,hfile,ghfile
real, dimension(:,:,:,:), allocatable :: pfile
real, dimension(:,:,:,:), allocatable :: eta_gcm
!real, dimension(:,:,:,:), allocatable :: tfileorig,ufileorig,vfileorig
real, dimension(:,:,:,:), allocatable :: tsoilfile, dsoilfile, isoilfile
real, dimension(:,:,:,:), allocatable :: waterfile, watericefile
real, dimension(:,:,:), allocatable :: swatericefile!, swaterfile
real, dimension(:,:,:,:), allocatable :: dustfile,dustnfile
real, dimension(:,:,:,:), allocatable :: co2file
real, dimension(:,:,:,:), allocatable :: ccnfile,ccnnfile

!! Reading the parameter file 
integer :: tmp,increment,FILES
integer :: tmp2,tmp3,tmp4
character*1 :: answer
character*4 :: str
character*2 :: str2, str3, str4
integer, dimension(:), allocatable :: time_out
character*13, dimension(:), allocatable :: date_out
character*19, dimension(:), allocatable :: date_out2

#ifdef PHOTOCHEM
real, dimension(:,:,:,:,:), allocatable :: chemtrac
integer :: nchemtrac,i
CHARACTER*20,DIMENSION(:),ALLOCATABLE :: wtnom
#endif


!***************************************************************************
!***************************************************************************


!!---------------------
!! Open the datafile
!!---------------------
ierr=NF_OPEN ("input_diagfi.nc",NF_NOWRITE,nid)
IF (ierr.NE.NF_NOERR) THEN
   write(*,*)'**** Please create a symbolic link called input_diagfi.nc'
   CALL ABORT
ENDIF


!!--------------------------
!! Ask for data references
!!--------------------------
!write(*,*) "Create n files. What is n ?"
read(*,*) FILES
allocate(time_out(FILES))
allocate(date_out(FILES))
allocate(date_out2(FILES))

!write(*,*) ""
!write(*,*) "INPUT: Diagfi time"
!write(*,*) "list of n subscripts, separated by <Return>s"
increment=1
do while (increment.ne.FILES+1)
    read(*,*) tmp
    time_out(increment)=tmp
    increment=increment+1
enddo

!write(*,*) ""
!write(*,*) "OUTPUT: WRF time"
!write(*,*) "fill 4 lines indicating"
!write(*,*) "-year (4 digits)"
!write(*,*) "-month (2 digits)" 
!write(*,*) "-day (2 digits)"
!write(*,*) "-hour (2 digits)"
increment=1
do while (increment.ne.FILES+1)
    read(*,*) str
    read(*,*) str2
    read(*,*) str3
    read(*,*) str4
    date_out(increment)=str//'-'//str2//'-'//str3//'_'//str4
    date_out2(increment)=str//'-'//str2//'-'//str3//'_'//str4//':00:00'
    !print *,date_out(increment)
    !write(*,*) "ok? (y/n)"
    read(*,*) answer
    if (answer.eq.'n') then
        !write(*,*) "please write again"
    else 
        !write(*,*) "next one, please"    
        increment=increment+1
    endif
enddo


!!-------------------
!! Get array sizes
!!-------------------
SELECT CASE(ident)
CASE('LMD')
ierr=NF_INQ_DIMID(nid,"time_counter",timedim)
CASE('OXF')
ierr=NF_INQ_DIMID(nid,"time",timedim)
END SELECT
ierr=NF_INQ_DIMLEN(nid,timedim,timelen)
  write(*,*) "timelen: ",timelen

SELECT CASE(ident)
CASE('LMD')
ierr=NF_INQ_DIMID(nid,"lat",latdim)
CASE('OXF')
ierr=NF_INQ_DIMID(nid,"lat",latdim)
END SELECT
ierr=NF_INQ_DIMLEN(nid,latdim,latlen)
  write(*,*) "latlen: ",latlen

SELECT CASE(ident)
CASE('LMD')
ierr=NF_INQ_DIMID(nid,"lon",londim)
CASE('OXF')
ierr=NF_INQ_DIMID(nid,"lon",londim)
END SELECT
ierr=NF_INQ_DIMLEN(nid,londim,lonlen)
  write(*,*) "lonlen: ",lonlen

SELECT CASE(ident)
CASE('LMD')
ierr=NF_INQ_DIMID(nid,"presnivs",altdim)
CASE('OXF')
ierr=NF_INQ_DIMID(nid,"sigma",altdim)
END SELECT
ierr=NF_INQ_DIMLEN(nid,altdim,altlen)
  write(*,*) "altlen: ",altlen



!!-------------------------
!! Allocate local arrays
!!-------------------------
allocate(eta_gcm(lonlen,latlen,altlen,timelen))
allocate(tfile(lonlen,latlen,altlen,timelen))
allocate(pfile(lonlen,latlen,altlen,timelen))
allocate(tsoilfile(lonlen,latlen,altlen,timelen))
allocate(dsoilfile(lonlen,latlen,altlen,timelen))
allocate(isoilfile(lonlen,latlen,altlen,timelen))
allocate(vertdsoil(altlen))
!allocate(tfileorig(lonlen,latlen,altlen,timelen))
allocate(ufile(lonlen,latlen,altlen,timelen))
!allocate(ufileorig(lonlen,latlen,altlen,timelen))
allocate(vfile(lonlen,latlen,altlen,timelen))
!allocate(vfileorig(lonlen,latlen,altlen,timelen))
allocate(rfile(lonlen,latlen,altlen,timelen))
allocate(hfile(lonlen,latlen,altlen,timelen))  
allocate(waterfile(lonlen,latlen,altlen,timelen))
allocate(watericefile(lonlen,latlen,altlen,timelen))
!allocate(swaterfile(lonlen,latlen,timelen))
allocate(swatericefile(lonlen,latlen,timelen))
allocate(dustfile(lonlen,latlen,altlen,timelen))
allocate(dustnfile(lonlen,latlen,altlen,timelen))
allocate(co2file(lonlen,latlen,altlen,timelen))
allocate(ccnfile(lonlen,latlen,altlen,timelen))
allocate(ccnnfile(lonlen,latlen,altlen,timelen))
allocate(psfile(lonlen,latlen,timelen)) 
allocate(tsfile(lonlen,latlen,timelen))
allocate(tnsfile(lonlen,latlen,timelen))
allocate(unsfile(lonlen,latlen,timelen))
allocate(vnsfile(lonlen,latlen,timelen))
allocate(emissfile(lonlen,latlen,timelen))
allocate(co2icefile(lonlen,latlen,timelen))
allocate(interm(lonlen,latlen))
allocate(gwparam(lonlen,latlen,5))
allocate(ghtsfile(lonlen,latlen))    !! no scan axis
allocate(ghfile(lonlen,latlen,altlen,timelen))
allocate(vide(lonlen,latlen))
allocate(ones(lonlen,latlen))
allocate(lat(latlen), lon(lonlen), alt(altlen), time(timelen))
allocate(aps(altlen),bps(altlen),levels(altlen))
#ifdef PHOTOCHEM
nchemtrac = 34
allocate(wtnom(nchemtrac))
print*,'PHOTOCHEM2.1'
wtnom(1)  = "co2"
wtnom(2)  = "co"
wtnom(3)  = "h2"
wtnom(4) = "h2o"
wtnom(5)  = "o1d"
wtnom(6)  = "o"
wtnom(7)  = "o2"
wtnom(8)  = "o2dg"
wtnom(9)  = "o3"
wtnom(10)  = "h"
wtnom(11)  = "oh"
wtnom(12)  = "ho2"
wtnom(13) = "h2o2"
wtnom(14)  = "cl"
wtnom(15)  = "clo"
wtnom(16)  = "cl2"
wtnom(17)  = "hcl"
wtnom(18)  = "hocl"
wtnom(19)  = "clco"
wtnom(20)  = "clco3"
wtnom(21)  = "cocl2"
wtnom(22)  = "s"
wtnom(23)  = "so"
wtnom(24)  = "so2"
wtnom(25)  = "so3"
wtnom(26)  = "s2o2"
wtnom(27)  = "ocs"
wtnom(28)  = "hso3"
wtnom(29)  = "h2so4"
wtnom(30)  = "s2"
wtnom(31)  = "clso2"
wtnom(32)  = "oscl"
wtnom(33)  = "h2oliq"
wtnom(34)  = "h2so4liq"
allocate(chemtrac(lonlen,latlen,altlen,timelen,nchemtrac))
chemtrac(:,:,:,:,:)=0
#endif

tfile(:,:,:,:)=0
pfile(:,:,:,:)=0
tsoilfile(:,:,:,:)=0 
isoilfile(:,:,:,:)=0 
dsoilfile(:,:,:,:)=0 
vertdsoil(:)=0.
!tfileorig(:,:,:,:)=0
!ufileorig(:,:,:,:)=0
!vfileorig(:,:,:,:)=0
ufile(:,:,:,:)=0 
vfile(:,:,:,:)=0 
rfile(:,:,:,:)=0 
hfile(:,:,:,:)=0
waterfile(:,:,:,:)=0
watericefile(:,:,:,:)=0
!swaterfile(:,:,:)=0
swatericefile(:,:,:)=0
dustfile(:,:,:,:)=0
dustnfile(:,:,:,:)=0
co2file(:,:,:,:)=0
ccnfile(:,:,:,:)=0
ccnnfile(:,:,:,:)=0
psfile(:,:,:)=0 
tsfile(:,:,:)=0
emissfile(:,:,:)=0
co2icefile(:,:,:)=0
tnsfile(:,:,:)=0
unsfile(:,:,:)=0
vnsfile(:,:,:)=0
interm(:,:)=0
gwparam(:,:,:)=0
ghtsfile(:,:)=0
ghfile(:,:,:,:)=0
vide(:,:)=0
ones(:,:)=0


!!-------------------
!! Read dimensions
!!-------------------

print *,'>>> Read dimensions !'

print *,'Time'
SELECT CASE(ident)
CASE('LMD')
   ierr = NF_INQ_VARID (nid, "time_counter",nvarid)
CASE('OXF')
   ierr = NF_INQ_VARID (nid, "time",nvarid)
END SELECT
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <Time> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, time)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, time)
#endif
   print *,time(1),' ... to ... ',time(timelen)

print *,'Latitude'
SELECT CASE(ident)
CASE('LMD')
   ierr = NF_INQ_VARID (nid, "lat",nvarid)
CASE('OXF')
   ierr = NF_INQ_VARID (nid, "lat",nvarid)
END SELECT
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <latitude> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, lat)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, lat)
#endif
   print *,lat(1),' ... to ... ',lat(latlen),' ... step: ',lat(latlen)-lat(latlen-1)

print *,'Longitude'
SELECT CASE(ident)
CASE('LMD')
   ierr = NF_INQ_VARID (nid, "lon",nvarid)
CASE('OXF')
   ierr = NF_INQ_VARID (nid, "lon",nvarid)
END SELECT
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <longitude> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, lon)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, lon)
#endif
   print *,lon(1),' ... to ... ',lon(lonlen),' ... step: ',lon(lonlen)-lon(lonlen-1)

!SELECT CASE(ident)
!CASE('LMD')
!print *, 'Hybrid coordinates'
!   ierr = NF_INQ_VARID (nid, "aps", nvarid)
!   IF (ierr .NE. NF_NOERR) THEN
!      PRINT *, "Error: Readmeteo <aps> not found"
!      stop
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, aps)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, aps)
!#endif
!   ierr = NF_INQ_VARID (nid, "bps", nvarid)
!   IF (ierr .NE. NF_NOERR) THEN
!      PRINT *, "Error: Readmeteo <bps> not found"
!      stop
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bps)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, bps)
!#endif
!   print *,aps(1),' ... to ... ',aps(altlen)
!   print *,bps(1),' ... to ... ',bps(altlen)
!CASE('OXF')
!   ierr = NF_INQ_VARID (nid, "sigma", nvarid)
!  IF (ierr .NE. NF_NOERR) THEN
!      PRINT *, "Error: Readmeteo <sigma> not found"
!      stop
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, bps)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, bps)
!#endif
!   print *,bps(1),' ... to ... ',bps(altlen)
!END SELECT


!!-------------------------------------------
!! Reading 2D and 3D meteorological fields
!!-------------------------------------------

IF (blank .EQV. .false.) THEN

print *,'>>> Read 2D optional fields !'

print *,'Emissivity'
   ierr = NF_INQ_VARID (nid,"emis",nvarid)
IF (ierr .NE. NF_NOERR) THEN
        PRINT *, '...warning: not found in diagfi !'
        PRINT *, '...will be filled with a prescribed value', emiss_prescribed
  emissfile(:,:,:)=emiss_prescribed
ELSE
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, emissfile)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, emissfile)
#endif
ENDIF   

print *,'CO2 ice'
   ierr = NF_INQ_VARID (nid,"co2ice",nvarid)
IF (ierr .NE. NF_NOERR) THEN
  PRINT *, '...warning: not found in diagfi !'
  PRINT *, '...will be filled with a prescribed value', co2ice_prescribed
  co2icefile(:,:,:)=co2ice_prescribed
ELSE
#ifdef NC_DOUBLE
  ierr = NF_GET_VAR_DOUBLE(nid, nvarid, co2icefile)
#else
  ierr = NF_GET_VAR_REAL(nid, nvarid, co2icefile)
#endif
ENDIF

print *,'>>> Read 2D surface fields !'

print *,'Surface Pressure'
   ierr = NF_INQ_VARID (nid,"psol",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <ps> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, psfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, psfile)
#endif

print *,'Ground Temperature'
   ierr = NF_INQ_VARID (nid,"tsol",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <tsurf> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tsfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, tsfile)
#endif


!!!"atmospheric" surface temperature is taken
!!!from original diagfi.nc first level
!!!... level is ~3-5 meters
!print *,'Near-Surface Temperature'
!   ierr = NF_INQ_VARID (nid,"temp",nvarid)
!   IF (ierr .NE. NF_NOERR) THEN
!      ierr = NF_INQ_VARID (nid,"t",nvarid)
!	IF (ierr .NE. NF_NOERR) THEN
!           PRINT *, "Error: Readmeteo <temp> not found"
!           stop
!	ENDIF
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tfileorig)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, tfileorig)
!#endif
!   tnsfile=tfileorig(:,:,1,:)

!!!"atmospheric" surface u is taken
!!!from original diagfi.nc first level
!!!... level is ~3-5 meters
!print *,'Near-Surface Zonal Wind'
!   ierr = NF_INQ_VARID (nid,"u",nvarid)
!   IF (ierr .NE. NF_NOERR) THEN
!     PRINT *, "Error: Readmeteo <u> not found"
!     stop
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ufileorig)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, ufileorig)
!#endif
!   unsfile=ufileorig(:,:,1,:)

!!!"atmospheric" surface v is taken
!!!from original diagfi.nc first level
!!!... level is ~3-5 meters
!print *,'Near-Surface Meridional Wind'
!   ierr = NF_INQ_VARID (nid,"v",nvarid)
!   IF (ierr .NE. NF_NOERR) THEN
!     PRINT *, "Error: Readmeteo <v> not found"
!     stop
!   ENDIF
!#ifdef NC_DOUBLE
!   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, vfileorig)
!#else
!   ierr = NF_GET_VAR_REAL(nid, nvarid, vfileorig)
!#endif
!   vnsfile=vfileorig(:,:,1,:)

print *,'>>> Read 3D meteorological fields ! - This may take some time ...'

print *,'Geopotential height'
   ierr = NF_INQ_VARID (nid,"geop",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <geop> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ghfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, ghfile)
#endif
ghtsfile=ghfile(:,:,1,1)/grav !surface geop

print *,'Pressure'
   ierr = NF_INQ_VARID (nid,"pres",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      ierr = NF_INQ_VARID (nid,"p",nvarid)
        IF (ierr .NE. NF_NOERR) THEN
          PRINT *, "Error: Readmeteo <p> not found"
          stop
        ENDIF
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, pfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, pfile)
#endif

print *,'Temperature'
   ierr = NF_INQ_VARID (nid,"temp",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      ierr = NF_INQ_VARID (nid,"t",nvarid)
        IF (ierr .NE. NF_NOERR) THEN
          PRINT *, "Error: Readmeteo <t> not found"
          stop
        ENDIF
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, tfile)
#endif
tnsfile=tfile(:,:,1,:)

print *,'Zonal wind'   
   ierr = NF_INQ_VARID (nid,"vitu",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <u> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, ufile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, ufile)
#endif
unsfile=ufile(:,:,1,:)

print *,'Meridional wind'
   ierr = NF_INQ_VARID (nid,"vitv",nvarid)
   IF (ierr .NE. NF_NOERR) THEN
      PRINT *, "Error: Readmeteo <v> not found"
      stop
   ENDIF
#ifdef NC_DOUBLE
   ierr = NF_GET_VAR_DOUBLE(nid, nvarid, vfile)
#else
   ierr = NF_GET_VAR_REAL(nid, nvarid, vfile)
#endif
vnsfile=vfile(:,:,1,:)

!!------------------------
!! special water stuff
!!------------------------
    print *,'Water vapor'
    ierr=NF_INQ_VARID(nid,"q02",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No q02 - Water vapor set to 0"
      waterfile(:,:,:,:)=0.  
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,waterfile)
    endif   

    print *,'Water ice'
    ierr=NF_INQ_VARID(nid,"q01",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No q01 - Water ice set to 0"  
      watericefile(:,:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,watericefile)
    endif
!    print *,'Surface Water vapor'
!    ierr=NF_INQ_VARID(nid,"qsurf02",nvarid)
!    if (ierr.ne.NF_NOERR) then
!      write(*,*) "...No qsurf02 - surface Water vapor set to 0"
!      swaterfile(:,:,:)=0.
!    endif
!    ierr=NF_GET_VAR_REAL(nid,nvarid,swaterfile)

    print *,'Surface Water ice'
!!!!!! ATTENTION ATTENTION 
!!!!!! water ice a la surface est qsurf(ig,nqmx)
    ierr=NF_INQ_VARID(nid,"qsurf02",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No qsurf02 - surface Water ice set to 0"
      swatericefile(:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,swatericefile)
    endif
!!------------------------
!! special water stuff
!!------------------------

    print *,'CO2 mass mixing ratio'
    ierr=NF_INQ_VARID(nid,"co2",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No co2 - co2 mixing ratio set to 0.95"
      co2file(:,:,:,:)=0.95
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,co2file)
    endif

!!------------------------
!! special dust stuff
!!------------------------
    print *,'Dust mass'
    ierr=NF_INQ_VARID(nid,"dustq",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No dustq - Dust mass set to 0"
      dustfile(:,:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,dustfile)
    endif

    print *,'Dust number'
    ierr=NF_INQ_VARID(nid,"dustN",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No dustN - Dust number set to 0"
      dustnfile(:,:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,dustnfile)
    endif

    print *,'CCN mass'
    ierr=NF_INQ_VARID(nid,"ccn",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No ccn - CCN mass set to 0"
      ccnfile(:,:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,ccnfile)
    endif

    print *,'CCN number'
    ierr=NF_INQ_VARID(nid,"ccnN",nvarid)
    if (ierr.ne.NF_NOERR) then
      write(*,*) "...No ccnN - CCN number set to 0"
      ccnnfile(:,:,:,:)=0.
    else
      ierr=NF_GET_VAR_REAL(nid,nvarid,ccnnfile)
    endif
!!------------------------
!! special dust stuff
!!------------------------

!SELECT CASE(ident)
!CASE('LMD')

    print *,'Soil temperatures'
    ierr=NF_INQ_VARID(nid,"tsoil",nvarid)
    if (ierr.ne.NF_NOERR) then
        write(*,*) "...No tsoil - Soil temperatures set isothermal with surface temperature"
        DO l=1,altlen
                tsoilfile(:,:,l,:)=tsfile(:,:,:)
        ENDDO
    else
        ierr=NF_GET_VAR_REAL(nid,nvarid,tsoilfile)
    endif

!!!!!!!!
!!!!!!!! new physics (but still compatible with old physics)
    print *,'Soil depths'
    ierr=NF_INQ_VARID(nid,"soildepth",nvarid)
    if (ierr.ne.NF_NOERR) then
        write(*,*) "...No soildepth - Set to -999"  !!! see soil_settings in LMD physics
        DO l=1,altlen
                vertdsoil(l)=-999.
        ENDDO
    else
        ierr=NF_GET_VAR_REAL(nid,nvarid,vertdsoil)
    endif
    print *, 'wait a minute' !! AS: I know this could be better
    DO m=1,lonlen
     DO n=1,latlen
      DO p=1,timelen
       dsoilfile(m,n,:,p) = vertdsoil(:) 
      ENDDO
     ENDDO
    ENDDO 
    DEALLOCATE(vertdsoil)

    print *,'Soil thermal inertia'
    ierr=NF_INQ_VARID(nid,"inertiedat",nvarid)
    if (ierr.ne.NF_NOERR) then
        write(*,*) "...No soil therm. inert. - Set to -999"
        DO l=1,altlen
                isoilfile(:,:,l,:)=-999.
        ENDDO
    else
        ierr=NF_GET_VAR_REAL(nid,nvarid,isoilfile)
    endif
!!!!!!!!
!!!!!!!! new physics


!!!!!!!!!!!!!!!!!!!!!!!!NEW PHYSICS + PHOTOCHEM
!!!!!!!!!!!!!!!!!!!!!!!!NEW PHYSICS + PHOTOCHEM
#ifdef PHOTOCHEM
    print *,'photochem'
    DO i=1,nchemtrac
     print *,wtnom(i)
     ierr=NF_INQ_VARID(nid,wtnom(i),nvarid)
     if (ierr.ne.NF_NOERR) then
       write(*,*) "...No ",wtnom(i), " - set to 0"
       chemtrac(:,:,:,:,i)=0.
     else
       ierr=NF_GET_VAR_REAL(nid,nvarid,chemtrac(:,:,:,:,i))
     endif
    ENDDO
#endif
!!!!!!!!!!!!!!!!!!!!!!!!NEW PHYSICS + PHOTOCHEM
!!!!!!!!!!!!!!!!!!!!!!!!NEW PHYSICS + PHOTOCHEM



!END SELECT

print*,'VENUS : rotation backward -> inversion of lat and long'
DO i=1,latlen
DO j=1,lonlen
  psfile(j,i,:)=psfile(lonlen+1-j,latlen+1-i,:)
  tsfile(i,j,:)=tsfile(lonlen+1-j,latlen+1-i,:)
  ghtsfile(j,i)=ghtsfile(lonlen+1-j,latlen+1-i)
  ghfile(j,i,:,:)=ghfile(lonlen+1-j,latlen+1-i,:,:)
  pfile(j,i,:,:)=pfile(lonlen+1-j,latlen+1-i,:,:)
  tfile(j,i,:,:)=tfile(lonlen+1-j,latlen+1-i,:,:)
  ufile(j,i,:,:)=ufile(lonlen+1-j,latlen+1-i,:,:)
  vfile(j,i,:,:)=vfile(lonlen+1-j,latlen+1-i,:,:)
  tsoilfile(j,i,:,:)=tsoilfile(lonlen+1-j,latlen+1-i,:,:)
#ifdef PHOTOCHEM
 chemtrac(j,i,:,:,:)=chemtrac(lonlen+1-j,latlen+1-i,:,:,:)
#endif
ENDDO
ENDDO  

ierr=NF_CLOSE(nid)

!!!!lott stuff
!print *,'get lott'
!ierr=NF_OPEN("init_lott.nc",NF_NOWRITE,nid3)
!ierr=NF_INQ_VARID(nid3,"ZMEA",nvarid)
!ierr=NF_GET_VAR_REAL(nid3,nvarid,interm)
!gwparam(:,:,1)=interm(:,:)
!ierr=NF_INQ_VARID(nid3,"ZSTD",nvarid)
!ierr=NF_GET_VAR_REAL(nid3,nvarid,interm)
!gwparam(:,:,2)=interm(:,:)
!ierr=NF_INQ_VARID(nid3,"ZSIG",nvarid)
!ierr=NF_GET_VAR_REAL(nid3,nvarid,interm)
!gwparam(:,:,3)=interm(:,:)
!ierr=NF_INQ_VARID(nid3,"ZGAM",nvarid)
!ierr=NF_GET_VAR_REAL(nid3,nvarid,interm)
!gwparam(:,:,4)=interm(:,:)
!ierr=NF_INQ_VARID(nid3,"ZTHE",nvarid)
!ierr=NF_GET_VAR_REAL(nid3,nvarid,interm)
!gwparam(:,:,5)=interm(:,:)
!ierr=NF_CLOSE(nid3)

ENDIF

!!-----------------------------
!! Loop on the written files
!! >>> what follows is pretty repetitive ... 
!!        gotta do something (one more loop?)
!!-----------------------------

!!! Equivalent eta levels for WRF interpolation in initialize_real
!print *,'Computing equivalent eta levels'
!	
!	!ptop will be passed through RH(surface)
!	ptop=MINVAL(aps(altlen)+bps(altlen)*psfile(:,:,:))
!	print *, 'ptop', ptop
!
!print *, 'sample: equivalent eta levels at i=1,j=1'
!DO k = 1,altlen
!        levels(k)=-k
!        eta_gcm(:,:,k,:)=(aps(k)+bps(k)*psfile(:,:,:)-ptop)/(psfile(:,:,:)-ptop)
!        print *, eta_gcm(1,1,k,1) 
!END DO

!!---better to pass pressure values 
!!---the eta treatment is done in module_initialize_real
DO k = 1,altlen
        levels(k)=-k
        !!dummy decreasing levels
END DO




!! Dummy surface level is XLVL=200100.


DO l=1,FILES


!!---------------------------------------------
!! Write file in intermediate format for WPS
!! 1. Surface data
!!---------------------------------------------
!!
!! a mettre pour tous sinon
!!     WRF_DEBUG: Warning DIM             4 , NAME num_metgrid_levels REDIFINED  by
!!            var DUSTN            26           25  in wrf_io.F90 line         2349
!!

!
! These variables remain the same in the "loop"
!
HDATE=date_out(l)
IFV=4
XFCST=0.
SOURCE=ident
NX=lonlen
NY=latlen
ALLOCATE(SLAB(NX,NY))
IPROJ=0
STARTLOC='SWCORNER'
STARTLAT=lat(1) 
STARTLON=lon(1)
DELTALAT=lat(latlen)-lat(latlen-1)
DELTALON=lon(lonlen)-lon(lonlen-1)
!
! Open the data destination file
!
output="./WPSFEED/"//ident//":"//date_out2(l)       
open(UNIT=1,			&
	FILE=output,		&
	STATUS='REPLACE',	&
	FORM='UNFORMATTED',	&
	ACCESS='SEQUENTIAL')

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='T'
UNITS='K'
DESC='Atmospheric temperature'
XLVL=200100.
SLAB=tnsfile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output
!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='U'
UNITS='m s{-1}'
DESC='Zonal wind'
XLVL=200100.
SLAB=unsfile(:,:,time_out(l))
        ! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON 
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output
!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='V'
UNITS='m s{-1}'
DESC='Meridional wind'
XLVL=200100.
SLAB=-1*vnsfile(:,:,time_out(l)) !VENUS ratoting backwards : v=-1*v
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON 
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output
!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='RH'
UNITS=''
DESC='Customized 2D static field'
XLVL=200100.
SLAB=vide(:,:)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SOILHGT'
UNITS='m' 
DESC='Terrain field of source analysis'  !Ground geopotential height
XLVL=200100.
SLAB=ghtsfile(:,:)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='PSFC'
UNITS='Pa'
DESC='Surface Pressure'
XLVL=200100.
SLAB=psfile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='ST000010'
UNITS=''
DESC='Emissivity'
XLVL=200100.
SLAB=emissfile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON 
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='ST010040'
UNITS=''
DESC='CO2 ice'
XLVL=200100.
SLAB=co2icefile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='ST040100'
UNITS=''
DESC='ZMEA'
XLVL=200100.
SLAB=gwparam(:,:,1)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='ST100200'
UNITS=''
DESC='ZSTD'
XLVL=200100.
SLAB=gwparam(:,:,2)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='LANDSEA'
UNITS='0/1 Flag'
DESC='Land/Sea flag'
XLVL=200100.
SLAB=ones(:,:)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SKINTEMP'
UNITS='K'
DESC='Ground temperature'
XLVL=200100.
SLAB=tsfile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SM000010'
UNITS=''
DESC='ZSIG'
XLVL=200100.
SLAB=gwparam(:,:,3)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SM010040'
UNITS=''
DESC='ZGAM'
XLVL=200100.
SLAB=gwparam(:,:,4)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SM040100'
UNITS=''
DESC='ZTHE'
XLVL=200100.
SLAB=gwparam(:,:,5)
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='SM100200'
UNITS='kg/m2'
DESC='Surf water ice'
XLVL=200100.
SLAB=swatericefile(:,:,time_out(l))
	! And now put everything in the destination file
	! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON 
	! ... Data
	write(1) SLAB
!print *,'The field '//DESC//' was written to '//output


!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='HV'
UNITS='kg/kg'
DESC='Water vapor'
XLVL=200100.
SLAB=waterfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='HI'
UNITS='kg/kg'
DESC='Water ice'
XLVL=200100.
SLAB=watericefile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='TSOIL'
UNITS='K'
DESC='Soil temperature'
XLVL=200100.
SLAB=tsoilfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DSOIL'
UNITS='m'
DESC='Soil depths'
XLVL=200100.
SLAB=dsoilfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='ISOIL'
UNITS='tiu'
DESC='Soil thermal inertia'
XLVL=200100.
SLAB=isoilfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CO2'
UNITS='kg/kg'
DESC='CO2 mixing ratio'
XLVL=200100.
SLAB=co2file(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DUSTQ'
UNITS='kg/kg'
DESC='Dust mixing ratio'
XLVL=200100.
SLAB=dustfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DUSTN'
UNITS='part/kg' 
DESC='Dust number density'
XLVL=200100.
SLAB=dustnfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CCNQ'
UNITS='kg/kg'
DESC='CCN mixing ratio'
XLVL=200100.
SLAB=ccnfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CCNN'
UNITS='part/kg'
DESC='CCN number density'
XLVL=200100.
SLAB=ccnnfile(:,:,1,time_out(l))
        ! And now put everything in the destination file
        ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
        ! ... Data
        write(1) SLAB
!print *,'The field '//DESC//' was written to '//output


!------------------------!
! >>> Write a variable   !
!     PHOTOCHEMISTRY     !
!------------------------!
#ifdef PHOTOCHEM
        FIELD='CO2'
        UNITS='kg/kg'
        DESC='CO2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),1)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='CO'
        UNITS='kg/kg'
        DESC='CO mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),2)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2'
        UNITS='kg/kg'
        DESC='H2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),3)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2O'
        UNITS='kg/kg'
        DESC='H2O mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),4)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='O1D'
        UNITS='kg/kg'
        DESC='O1D mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),5)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='O'
        UNITS='kg/kg'
        DESC='O mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),6)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='O2'
        UNITS='kg/kg'
        DESC='O2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),7)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='O2DG'
        UNITS='kg/kg'
        DESC='O2DG mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),8)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='O3'
        UNITS='kg/kg'
        DESC='O3 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),9)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H'
        UNITS='kg/kg'
        DESC='H mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),10)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='OH'
        UNITS='kg/kg'
        DESC='OH mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),11)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='HO2'
        UNITS='kg/kg'
        DESC='hO2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),12)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2O2'
        UNITS='kg/kg'
        DESC='H2O2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),13)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='Cl'
        UNITS='kg/kg'
        DESC='Cl mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),14)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='ClO'
        UNITS='kg/kg'
        DESC='ClO mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),15)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='Cl2'
        UNITS='kg/kg'
        DESC='Cl2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),16)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='HCl'
        UNITS='kg/kg'
        DESC='HCl mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),17)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='HOCl'
        UNITS='kg/kg'
        DESC='HOCl mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),18)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='ClCO'
        UNITS='kg/kg'
        DESC='ClCO mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),19)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='ClCO3'
        UNITS='kg/kg'
        DESC='ClCO3 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),20)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='COCl2'
        UNITS='kg/kg'
        DESC='COClC2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),21)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='S'
        UNITS='kg/kg'
        DESC='S mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),22)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='SO'
        UNITS='kg/kg'
        DESC='SO mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),23)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='SO2'
        UNITS='kg/kg'
        DESC='SO2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),24)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='SO3'
        UNITS='kg/kg'
        DESC='SO3 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),25)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='S2O2'
        UNITS='kg/kg'
        DESC='S2O2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),26)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='OCS'
        UNITS='kg/kg'
        DESC='OCS mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),27)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='HSO3'
        UNITS='kg/kg'
        DESC='HSO3 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),28)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2SO4'
        UNITS='kg/kg'
        DESC='H2SO4 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),29)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='S2'
        UNITS='kg/kg'
        DESC='S2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),30)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='ClSO2'
        UNITS='kg/kg'
        DESC='ClSO2 mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),31)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='OSCl'
        UNITS='kg/kg'
        DESC='OSCl mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),32)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2Oliq'
        UNITS='kg/kg'
        DESC='H2Oliq mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),33)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB


        FIELD='H2SO4liq'
        UNITS='kg/kg'
        DESC='H2SO4liq mixing ratio'
        XLVL=levels(k)
        SLAB=chemtrac(:,:,1,time_out(l),34)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
        write(1) SLAB
#endif

!!----------------------------------------------------
!! Write file in intermediate format for WPS
!! 2. 3D meteorological data
!! >>> same stuff, but handle vertical profile
!! >>> !! XLVL must be decreasing (pressure levels) 
!!----------------------------------------------------

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='T'
UNITS='K'
DESC='Atmospheric temperature'
DO k = 1,altlen
	XLVL=levels(k)
        SLAB=tfile(:,:,k,time_out(l))
		! And now put everything in the destination file
		! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
		! ... Data
		write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='U'
UNITS='m s{-1}'
DESC='Zonal wind'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=ufile(:,:,k,time_out(l))
		! And now put everything in the destination file
		! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
		! ... Data
		write(1) SLAB
                !write(1) ufile(:,:,k,time_out(l))
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='V'
UNITS='m s{-1}'
DESC='Meridional wind'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=-1*vfile(:,:,k,time_out(l)) ! VENUS rotation backward : v=-1*v
		! And now put everything in the destination file
		! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
		! ... Data
		write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='RH'
UNITS=''
DESC='Customized 2D static field'
DESC='Eta-levels from the GCM'
DESC='Pressure values from the GCM'
DO k = 1,altlen
        XLVL=levels(k)
SELECT CASE(ident)
CASE('LMD')
        SLAB=pfile(:,:,k,time_out(l))
CASE('OXF')
	!SLAB=bps(k)*psfile(:,:,time_out(l))
END SELECT
		! And now put everything in the destination file
		! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
		! ... Data
		write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!    
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------! 
FIELD='HGT'
UNITS='m'
DESC='Height'
DO k = 1,altlen
        XLVL=levels(k)
SELECT CASE(ident)
CASE('LMD')
        SLAB=(ghfile(:,:,k,time_out(l)))/grav
CASE('OXF')
	!SLAB=10000.*alog(610./(bps(k)*psfile(:,:,time_out(l))))
END SELECT
		! And now put everything in the destination file
		! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON  
		! ... Data
		write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='HV'
UNITS='kg/kg'
DESC='Water vapor'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=waterfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='HI'
UNITS='kg/kg'
DESC='Water ice'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=watericefile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='TSOIL'
UNITS='K'
DESC='Soil temperature'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=tsoilfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DSOIL'
UNITS='m'
DESC='Soil depths'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=dsoilfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='ISOIL'
UNITS='tiu'
DESC='Soil thermal inertia'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=isoilfile(:,:,k,time_out(l))
        
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CO2'
UNITS='kg/kg'
DESC='CO2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=co2file(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DUSTQ'
UNITS='kg/kg'
DESC='Dust mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=dustfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='DUSTN'
UNITS='part/kg'
DESC='Dust number density'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=dustnfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CCNQ'
UNITS='kg/kg'
DESC='CCN mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=ccnfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output

!------------------------!
! >>> Write a variable   !
!    ... Copy&Paste part !
!------------------------!
FIELD='CCNN'
UNITS='part/kg'
DESC='CCN number density'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=ccnnfile(:,:,k,time_out(l))
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
!print *,'The field '//DESC//' was written to '//output


!------------------------!
! >>> Write a variable   !
!     PHOTOCHEMISTRY     !
!------------------------!
#ifdef PHOTOCHEM
FIELD='CO2'
UNITS='kg/kg'
DESC='CO2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),1)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='CO'
UNITS='kg/kg'
DESC='CO mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),2)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2'
UNITS='kg/kg'
DESC='H2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),3)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2O'
UNITS='kg/kg'
DESC='H2O mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),4)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='O1D'
UNITS='kg/kg'
DESC='O1D mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),5)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='O'
UNITS='kg/kg'
DESC='O mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),6)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='O2'
UNITS='kg/kg'
DESC='O2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),7)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='O2DG'
UNITS='kg/kg'
DESC='O2DG mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),8)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='O3'
UNITS='kg/kg'
DESC='O3 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),9)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H'
UNITS='kg/kg'
DESC='H mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),10)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='OH'
UNITS='kg/kg'
DESC='OH mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),11)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='HO2'
UNITS='kg/kg'
DESC='hO2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),12)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2O2'
UNITS='kg/kg'
DESC='H2O2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),13)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='Cl'
UNITS='kg/kg'
DESC='Cl mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),14)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='ClO'
UNITS='kg/kg'
DESC='ClO mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),15)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='Cl2'
UNITS='kg/kg'
DESC='Cl2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),16)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='HCl'
UNITS='kg/kg'
DESC='HCl mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),17)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='HOCl'
UNITS='kg/kg'
DESC='HOCl mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),18)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='ClCO'
UNITS='kg/kg'
DESC='ClCO mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),19)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='ClCO3'
UNITS='kg/kg'
DESC='ClCO3 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),20)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='COCl2'
UNITS='kg/kg'
DESC='COClC2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),21)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='S'
UNITS='kg/kg'
DESC='S mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),22)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='SO'
UNITS='kg/kg'
DESC='SO mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),23)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='SO2'
UNITS='kg/kg'
DESC='SO2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),24)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='SO3'
UNITS='kg/kg'
DESC='SO3 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),25)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='S2O2'
UNITS='kg/kg'
DESC='S2O2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),26)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='OCS'
UNITS='kg/kg'
DESC='OCS mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),27)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='HSO3'
UNITS='kg/kg'
DESC='HSO3 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),28)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2SO4'
UNITS='kg/kg'
DESC='H2SO4 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),29)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='S2'
UNITS='kg/kg'
DESC='S2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),30)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='ClSO2'
UNITS='kg/kg'
DESC='ClSO2 mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),31)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='OSCl'
UNITS='kg/kg'
DESC='OSCl mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),32)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2Oliq'
UNITS='kg/kg'
DESC='H2Oliq mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),33)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO

FIELD='H2SO4liq'
UNITS='kg/kg'
DESC='H2SO4liq mixing ratio'
DO k = 1,altlen
        XLVL=levels(k)
        SLAB=chemtrac(:,:,k,time_out(l),34)
                ! And now put everything in the destination file
                ! ... Header
        write(1) IFV
        write(1) HDATE,XFCST,SOURCE,FIELD,UNITS,DESC,XLVL,NX,NY,IPROJ
        write(1) STARTLOC,STARTLAT,STARTLON,DELTALAT,DELTALON
                ! ... Data
                write(1) SLAB
END DO
#endif

print *,'****done file '//output, int(100.*float(l)/float(FILES)), ' % '
close(1)

DEALLOCATE(SLAB)

END DO !! end of the files loop


!!-------------------------
!! DeAllocate local arrays
!!-------------------------
deallocate(eta_gcm)
deallocate(tfile)
deallocate(pfile)
deallocate(tsoilfile)
deallocate(isoilfile)
deallocate(dsoilfile)
!deallocate(tfileorig)
deallocate(ufile)
!deallocate(ufileorig)
deallocate(vfile)
!deallocate(vfileorig)
deallocate(rfile)
deallocate(hfile)  
deallocate(waterfile)
deallocate(watericefile)
!deallocate(swaterfile)
deallocate(swatericefile)
deallocate(dustfile)
deallocate(dustnfile)
deallocate(co2file)
deallocate(ccnfile)
deallocate(ccnnfile)
deallocate(psfile) 
deallocate(tsfile)
deallocate(tnsfile)
deallocate(unsfile)
deallocate(vnsfile)
deallocate(emissfile)
deallocate(co2icefile)
deallocate(ghtsfile)	!! no scan axis
deallocate(ghfile)
deallocate(vide)
deallocate(ones)
deallocate(lat, lon, alt, time)
deallocate(aps,bps,levels)

#ifdef PHOTOCHEM
deallocate(chemtrac)
deallocate(wtnom)
#endif

print *, '------------------------'
print *, 'End of pre-WPS process !'
print *, '------------------------'
print *, 'Here is what you should set in namelist.wps:' 
print *, " start_date = '"//date_out2(1)//"'"
print *, " end_date = '"//date_out2(FILES)//"'"
print *, '------------------------'
print *, 'Any date between those is ok'
print *, 'for example, second available date is ... '//date_out2(2)
print *, '------------------------'

end
