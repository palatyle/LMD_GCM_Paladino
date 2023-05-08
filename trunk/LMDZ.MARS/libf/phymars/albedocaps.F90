subroutine albedocaps(zls,ngrid,piceco2,psolaralb,emisref)

! routine which changes the albedo (and emissivity) of the surface
! depending on the presence of CO2 ice on the surface

! to use the 'getin' routine 
use ioipsl_getincom, only: getin
use geometry_mod, only: latitude ! grid point latitudes (rad)
use surfdat_h, only: TESicealbedo, TESice_Ncoef, TESice_Scoef, &
                     emisice, albedice, watercaptag, albedo_h2o_ice, &
                     emissiv, albedodat
implicit none

include"callkeys.h"

! arguments:
real,intent(in) :: zls ! solar longitude (rad)
integer,intent(in) :: ngrid
real,intent(in) :: piceco2(ngrid) ! amount of CO2 ice on the surface (kg/m2)
real,intent(out) :: psolaralb(ngrid,2) ! albedo of the surface
real,intent(out) :: emisref(ngrid) ! emissivity of the surface


! local variables:
logical,save :: firstcall=.true.
integer :: ig,icap

! 1. Initializations
! AS: OK firstcall absolute
if (firstcall) then
  ! find out if user wants to use TES cap albedoes or not
  TESicealbedo=.false. ! default value 
  write(*,*)" albedocaps: Use TES Cap albedoes ?"
  call getin("TESicealbedo",TESicealbedo)
  write(*,*)" albedocaps: TESicealbedo = ",TESicealbedo
  
  ! if using TES albedoes, load coeffcients
  if (TESicealbedo) then
    write(*,*)" albedocaps: Coefficient for Northern Cap ?"
    TESice_Ncoef=1.0 ! default value
    call getin("TESice_Ncoef",TESice_Ncoef)
    write(*,*)" albedocaps: TESice_Ncoef = ",TESice_Ncoef
    
    write(*,*)" albedocaps: Coefficient for Southern Cap ?"
    TESice_Scoef=1.0 ! default value
    call getin("TESice_Scoef",TESice_Scoef)
    write(*,*)" albedocaps: TESice_Scoef = ",TESice_Scoef
  endif
  
  firstcall=.false.
endif ! of if (firstcall)

do ig=1,ngrid
  if (latitude(ig).lt.0.) then
    icap=2 ! Southern hemisphere
  else
    icap=1 ! Northern hemisphere
  endif

  if (piceco2(ig).gt.0) then
    ! set emissivity of surface to be the ice emissivity
    emisref(ig)=emisice(icap)
    ! set the surface albedo to be the ice albedo
    if (TESicealbedo) then
      call TES_icecap_albedo(zls,ig,psolaralb(ig,1),icap)
      psolaralb(ig,2)=psolaralb(ig,1)
    else
      psolaralb(ig,1)=albedice(icap)
      psolaralb(ig,2)=albedice(icap)
    endif
  else if (watercaptag(ig) .and. water) then
    ! there is a water ice cap: set the surface albedo to the water ice one
    ! to do : emissivity
    emisref(ig) = 1
    psolaralb(ig,1)=albedo_h2o_ice
    psolaralb(ig,2)=albedo_h2o_ice
  else
    ! set emissivity of surface to be bare ground emissivity
    emisref(ig)=emissiv
    ! set the surface albedo to bare ground albedo 
    psolaralb(ig,1)=albedodat(ig)
    psolaralb(ig,2)=albedodat(ig)
  endif ! of if (piceco2(ig).gt.0)
enddo ! of ig=1,ngrid
end subroutine albedocaps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TES_icecap_albedo(zls,ig,alb,icap)

use geometry_mod, only: latitude, longitude ! in radians
use surfdat_h, only: albedice, TESice_Ncoef, TESice_Scoef
use datafile_mod, only: datadir
use netcdf, only: nf90_open, NF90_NOWRITE, NF90_NOERR, &
                  nf90_strerror, nf90_inq_varid, nf90_get_var, nf90_close
                  
implicit none

! arguments:
real,intent(in) :: zls ! solar longitude (rad)
integer,intent(in) :: ig ! grid point index
real,intent(out) :: alb ! (interpolated) TES ice albedo at that grid point
integer :: icap ! =1: Northern hemisphere =2: Southern hemisphere

! local variables:
logical,save :: firstcall=.true.
real,save :: zls_old ! value of zls from a previous call
integer,save :: tinf,tsup ! encompassing time indexes of TES data
real,save :: reltime ! relative position in-between time indexes (in [0;1])
integer :: latinf,latsup ! encompassing latitude indexes of TES data
real :: rellat ! relative position in-between latitude indexes (in[0;1])
integer :: loninf,lonsup ! encompassing longitude indexes of TES data
real :: rellon !relative position in-between longitude indexes (in[0;1])
real,save :: pi,radeg ! to convert radians to degrees
real :: zlsd ! solar longitude, in degrees
real :: latd ! latitude, in degrees
real :: lond ! longitude, in degrees
integer :: i

! TES datasets: (hard coded fixed length/sizes; for now)
integer,parameter :: TESlonsize=72
real,parameter :: TESdeltalon=5.0 ! step in longitude in TES files
! longitudes, in TES files, in degrees, from TESlon(1)=-177.5 to TESlon(72)=177.5
real,save :: TESlon(TESlonsize) 
integer,parameter :: TESlatsize=30
real,parameter :: TESdeltalat=2.0 ! step in latitude in TES files
! latitudes (north hemisphere file), in degrees, from TESlatn(1)=31,
! to TESlatn(30)=89 ; TESlatn(8)=45
real,parameter :: TESlatnmin=45. ! minimum TES latitude (North hemisphere)
real,parameter :: TESlatsmax=-45. ! maximum TES latitude (South hemisphere)
real,save :: TESlatn(TESlatsize)
! latitudes (south hemisphere file), in degrees, from TESlats(1)=-89,
! to TESlats(30)=-31 ; TESlats(23)=-45
real,save :: TESlats(TESlatsize)
integer,parameter :: TESlssize=72
real,parameter :: TESdeltals=5.0 ! step in solar longitude in TES files
! Solar longitude in TES files, TESls(1)=2.5 to TESls(72)=357.5
real,save :: TESls(TESlssize)
! TES North albedo (=-1 for missing values)
real,save :: TESalbn(TESlonsize,TESlatsize,TESlssize)
! TES South albedo (=-1 for missing values)
real,save :: TESalbs(TESlonsize,TESlatsize,TESlssize)
! encompassing nodes arranged as follow : 4 3
real :: val(4)                          ! 1 2 

!NetCDF variables:
integer :: ierr ! NetCDF status
integer :: nid ! NetCDF file ID
integer :: nvarid ! NetCDF variable ID

! 0. Preliminary stuff
if (firstcall) then
! Load TES albedoes for Northern Hemisphere
  ierr=nf90_open(trim(datadir)//"/npsc_albedo.nc",NF90_NOWRITE,nid)
  IF (ierr.NE.NF90_NOERR) THEN
    write(*,*)'Problem opening npsc_albedo.nc (phymars/albedocaps.F90)'
    write(*,*)'It should be in :',trim(datadir),'/'
    write(*,*)'1) You can change this directory address in callfis.def with'
    write(*,*)'   datadir=/path/to/datafiles'
    write(*,*)'2) If necessary, npsc_albedo.nc (and other datafiles)'
    write(*,*)'   can be obtained online on:'
    write(*,*)'   http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
    CALL ABORT
  ELSE
    write(*,*) "albedocaps: using file ",trim(datadir)//"/npsc_albedo.nc"
  ENDIF
  
  ierr=nf90_inq_varid(nid,"longitude",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find longitude in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESlon)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading longitude data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif
  
  ierr=nf90_inq_varid(nid,"latitude",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find latitude in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESlatn)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading latitude data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif
  
  ierr=nf90_inq_varid(nid,"time",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find time in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESls)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading time data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif
  
  ierr=nf90_inq_varid(nid,"albedo",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find albedo in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESalbn)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading albedo data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif

  ierr=nf90_close(nid)

! Load albedoes for Southern Hemisphere
  ierr=nf90_open(trim(datadir)//"/spsc_albedo.nc",NF90_NOWRITE,nid)
  IF (ierr.NE.NF90_NOERR) THEN
    write(*,*)'Problem opening spsc_albedo.nc (phymars/albedocaps.F90)'
    write(*,*)'It should be in :',trim(datadir),'/'
    write(*,*)'1) You can change this directory address in callfis.def with'
    write(*,*)'   datadir=/path/to/datafiles'
    write(*,*)'2) If necessary, spsc_albedo.nc (and other datafiles)'
    write(*,*)'   can be obtained online on:'
    write(*,*)'   http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
    CALL ABORT
  ELSE
    write(*,*) "albedocaps: using file ",trim(datadir)//"/spsc_albedo.nc"
  ENDIF

  ierr=nf90_inq_varid(nid,"latitude",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find latitude in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESlats)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading latitude data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif

  ierr=nf90_inq_varid(nid,"albedo",nvarid)
  if (ierr.ne.NF90_NOERR) then
    write(*,*) "Failed to find albedo in file!"
    write(*,*)trim(nf90_strerror(ierr))
    stop
  else
    ierr=nf90_get_var(nid,nvarid,TESalbs)
    if (ierr.ne.NF90_NOERR) then
      write(*,*) "Failed loading albedo data from file!"
      write(*,*)trim(nf90_strerror(ierr))
      stop
    endif
  endif

  ierr=nf90_close(nid)

! constants:
  pi=acos(-1.)
  radeg=180/pi

  zls_old=-999 ! dummy initialization

  firstcall=.false.
endif ! of if firstcall

! 1. Identify encompassing latitudes

! Check that latitude is such that there is TES data to use
! (ie: latitude 45 deg and poleward) otherwise use 'default' albedoes
latd=latitude(ig)*radeg ! latitude, in degrees
if (icap.eq.1) then
 ! North hemisphere
 if (latd.lt.TESlatnmin) then
   alb=albedice(1)
   ! the job is done; quit this routine
   return
 else
   ! find encompassing latitudes
   if (latd.ge.TESlatn(TESlatsize)) then
     latinf=TESlatsize
     latsup=TESlatsize
     rellat=0.
   else
     do i=1,TESlatsize-1
       if ((latd.ge.TESlatn(i)).and.(latd.lt.TESlatn(i+1))) then
         latinf=i
         latsup=i+1
         rellat=(latd-TESlatn(i))/TESdeltalat
         exit ! found encompassing indexes; quit loop
       endif
     enddo
   endif
 endif ! of if (latd.lt.TESlatnmin)
else ! icap=2
 ! South hemisphere
 if (latd.gt.TESlatsmax) then
   alb=albedice(2)
   ! the job is done; quit this routine
   return
 else
   ! find encompassing latitudes
   if (latd.lt.TESlats(1)) then
     latinf=1
     latsup=1
     rellat=0.
   else
     do i=1,TESlatsize-1
       if ((latd.ge.TESlats(i)).and.(latd.lt.TESlats(i+1))) then
         latinf=i
         latsup=i+1
         rellat=(latd-TESlats(i))/TESdeltalat
         exit ! found encompassing indexes; quit loop
       endif
     enddo
   endif
 endif ! of if (latd.gt.-45.)
endif ! of if (icap.eq.1)

! 2. Identify encompassing time indexes
if (zls.ne.zls_old) then
  zlsd=zls*radeg ! solar longitude, in degrees
  
  if (zlsd.lt.TESls(1)) then
    tinf=TESlssize
    tsup=1
    reltime=0.5+zlsd/TESdeltals
  else
    if (zlsd.ge.TESls(TESlssize)) then
      tinf=TESlssize
      tsup=1
      reltime=(360.-zlsd)/TESdeltals
    else
      ! look for encompassing indexes
      do i=1,TESlssize-1
        if ((zlsd.ge.TESls(i)).and.(zlsd.lt.TESls(i+1))) then
          tinf=i
          tsup=i+1
          reltime=(zlsd-TESls(i))/TESdeltals
          exit ! quit loop, we found the indexes
        endif
      enddo
    endif
  endif ! of if (zlsd.lt.TESls(1))
  
  zls_old=zls ! store current zls
endif ! of if (zls.ne.zls_old)

! 3. Identify encompassing longitudes
lond=longitude(ig)*radeg ! east longitude, in degrees
if (lond.lt.TESlon(1)) then
  loninf=TESlonsize
  lonsup=1
  rellon=0.5+(180.+lond)/TESdeltalon
else
  if (lond.ge.TESlon(TESlonsize)) then
    loninf=TESlonsize
    lonsup=1
    rellon=(180-lond)/TESdeltalon
  else
    do i=1,TESlonsize-1
      if ((lond.ge.TESlon(i)).and.(lond.lt.TESlon(i+1))) then
        loninf=i
        lonsup=i+1
        rellon=(lond-TESlon(i))/TESdeltalon
        exit ! quit loop, we found the indexes
      endif
    enddo
  endif ! of if (lond.ge.TESlon(TESlonsize))
endif ! of if (lond.lt.TESlon(1))

! 4. Use linear interpolation in time to build encompassing nodal values
!    encompassing nodes are arranged as follow : 4 3
!                                                1 2
if (icap.eq.1) then
  ! Northern hemisphere
  val(1)=(1.-reltime)*TESalbn(loninf,latinf,tinf) &
         +reltime*TESalbn(loninf,latinf,tsup)
  val(2)=(1.-reltime)*TESalbn(lonsup,latinf,tinf) &
         +reltime*TESalbn(lonsup,latinf,tsup)
  val(3)=(1.-reltime)*TESalbn(lonsup,latsup,tinf) &
         +reltime*TESalbn(lonsup,latsup,tsup)
  val(4)=(1.-reltime)*TESalbn(loninf,latsup,tinf) &
         +reltime*TESalbn(loninf,latsup,tsup)
else
  ! Southern hemisphere
  val(1)=(1.-reltime)*TESalbs(loninf,latinf,tinf) &
         +reltime*TESalbs(loninf,latinf,tsup)
  val(2)=(1.-reltime)*TESalbs(lonsup,latinf,tinf) &
         +reltime*TESalbs(lonsup,latinf,tsup)
  val(3)=(1.-reltime)*TESalbs(lonsup,latsup,tinf) &
         +reltime*TESalbs(lonsup,latsup,tsup)
  val(4)=(1.-reltime)*TESalbs(loninf,latsup,tinf) &
         +reltime*TESalbs(loninf,latsup,tsup)
endif ! of if (icap.eq.1)

! 5. Use bilinear interpolation to compute albedo
alb=(1.-rellon)*(1.-rellat)*val(1) &
    +rellon*(1.-rellat)*val(2) &
    +rellon*rellat*val(3) &
    +(1.-rellon)*rellat*val(4)

! 6. Apply coefficient to interpolated TES albedo
if (icap.eq.1) then
  alb=alb*TESice_Ncoef
else
  alb=alb*TESice_Scoef
endif ! of if (icap.eq.1)

! Make sure that returned albedo is never greater than 0.90
if (alb.gt.0.90) alb=0.90

end subroutine TES_icecap_albedo

