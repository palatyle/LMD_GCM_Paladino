module phyetat0_mod

implicit none

contains

subroutine phyetat0 (startphy_file, &
                     ngrid,nlayer,fichnom,tab0,Lmodif,nsoil,nq, &
                     day_ini,time,tsurf,tsoil, &
                     emis,q2,qsurf,cloudfrac,totcloudfrac,hice, &
                     rnat,pctsrf_sic,tslab,tsea_ice,sea_ice)

! to use  'getin_p'
      use ioipsl_getin_p_mod, only: getin_p

  use tabfi_mod, only: tabfi
  USE tracer_h, ONLY: noms
  USE surfdat_h, only: phisfi, albedodat, zmea, zstd, zsig, zgam, zthe
  use iostart, only: nid_start, open_startphy, close_startphy, &
                     get_field, get_var, inquire_field, &
                     inquire_dimension, inquire_dimension_length
  use slab_ice_h, only: noceanmx
  use callkeys_mod, only: CLFvarying,surfalbedo,surfemis

  implicit none

!======================================================================
!  Arguments:
!  ---------
!  inputs:
  logical,intent(in) :: startphy_file ! .true. if reading start file
  integer,intent(in) :: ngrid
  integer,intent(in) :: nlayer
  character*(*),intent(in) :: fichnom ! "startfi.nc" file
  integer,intent(in) :: tab0
  integer,intent(in) :: Lmodif
  integer,intent(in) :: nsoil ! # of soil layers
  integer,intent(in) :: nq
  integer,intent(out) :: day_ini
  real,intent(out) :: time

!  outputs:
  real,intent(out) :: tsurf(ngrid) ! surface temperature
  real,intent(out) :: tsoil(ngrid,nsoil) ! soil temperature
  real,intent(out) :: emis(ngrid) ! surface emissivity
  real,intent(out) :: q2(ngrid,nlayer+1) ! 
  real,intent(out) :: qsurf(ngrid,nq) ! tracers on surface
! real co2ice(ngrid) ! co2 ice cover
  real,intent(out) :: cloudfrac(ngrid,nlayer)
  real,intent(out) :: hice(ngrid), totcloudfrac(ngrid)
  real,intent(out) :: pctsrf_sic(ngrid),tslab(ngrid,noceanmx)  
  real,intent(out) :: tsea_ice(ngrid),sea_ice(ngrid)
  real,intent(out) :: rnat(ngrid)

!======================================================================
!  Local variables:

!      INTEGER radpas
!      REAL co2_ppm
!      REAL solaire

      real xmin,xmax ! to display min and max of a field
!
      INTEGER ig,iq,lmax
      INTEGER nid, nvarid
      INTEGER ierr, i, nsrf
!      integer isoil 
!      INTEGER length
!      PARAMETER (length=100)
      CHARACTER*7 str7
      CHARACTER*2 str2
      CHARACTER*1 yes
!
      REAL p_rad,p_omeg,p_g,p_cpp,p_mugaz,p_daysec
      INTEGER nqold

! flag which identifies if 'startfi.nc' file is using old names (qsurf01,...)
!      logical :: oldtracernames=.false.
      integer :: count
      character(len=30) :: txt ! to store some text
      
      INTEGER :: indextime=1 ! index of selected time, default value=1
      logical :: found
      
      character(len=8) :: modname="phyetat0"

!
! ALLOCATE ARRAYS IN surfdat_h
!
IF (.not. ALLOCATED(albedodat)) ALLOCATE(albedodat(ngrid))
IF (.not. ALLOCATED(phisfi)) ALLOCATE(phisfi(ngrid))
IF (.not. ALLOCATED(zmea)) ALLOCATE(zmea(ngrid))
IF (.not. ALLOCATED(zstd)) ALLOCATE(zstd(ngrid))
IF (.not. ALLOCATED(zsig)) ALLOCATE(zsig(ngrid))
IF (.not. ALLOCATED(zgam)) ALLOCATE(zgam(ngrid))
IF (.not. ALLOCATED(zthe)) ALLOCATE(zthe(ngrid))

if (startphy_file) then
  ! open physics initial state file:
  call open_startphy(fichnom)

  ! possibility to modify tab_cntrl in tabfi
  write(*,*)
  write(*,*) 'TABFI in phyeta0: Lmodif=',Lmodif," tab0=",tab0
  call tabfi (ngrid,nid_start,Lmodif,tab0,day_ini,lmax,p_rad, &
                   p_omeg,p_g,p_cpp,p_mugaz,p_daysec,time)

else ! "academic" initialization of planetary parameters via tabfi
  call tabfi (ngrid,0,0,0,day_ini,lmax,p_rad, &
                   p_omeg,p_g,p_cpp,p_mugaz,p_daysec,time)
endif ! of if (startphy_file)

if (startphy_file) then
  ! Load surface geopotential:
  call get_field("phisfi",phisfi,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <phisfi>",1)
  endif
else
  phisfi(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: surface geopotential <phisfi> range:", &
               minval(phisfi), maxval(phisfi)

if (startphy_file) then
  ! Load bare ground albedo:
  call get_field("albedodat",albedodat,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <albedodat>",1)
  endif
else
  ! If no startfi file, use parameter surfalbedo in def file
  surfalbedo=0.5
  call getin_p("surfalbedo",surfalbedo)
  print*,"surfalbedo",surfalbedo
  albedodat(:)=surfalbedo
endif ! of if (startphy_file)
write(*,*) "phyetat0: Bare ground albedo <albedodat> range:", &
             minval(albedodat), maxval(albedodat)

! ZMEA
if (startphy_file) then
  call get_field("ZMEA",zmea,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <ZMEA>",1)
  endif
else
  zmea(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZMEA> range:", &
             minval(zmea), maxval(zmea)

! ZSTD
if (startphy_file) then
  call get_field("ZSTD",zstd,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <ZSTD>",1)
  endif
else
  zstd(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZSTD> range:", &
             minval(zstd), maxval(zstd)

! ZSIG
if (startphy_file) then
  call get_field("ZSIG",zsig,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <ZSIG>",1)
  endif
else
  zsig(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZSIG> range:", &
             minval(zsig), maxval(zsig)

! ZGAM
if (startphy_file) then
  call get_field("ZGAM",zgam,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <ZGAM>",1)
  endif
else
  zgam(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZGAM> range:", &
             minval(zgam), maxval(zgam)

! ZTHE
if (startphy_file) then
  call get_field("ZTHE",zthe,found)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <ZTHE>",1)
  endif
else
  zthe(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZTHE> range:", &
             minval(zthe), maxval(zthe)

! Surface temperature :
if (startphy_file) then
  call get_field("tsurf",tsurf,found,indextime)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <tsurf>",1)
  endif
else
  tsurf(:)=0. ! will be updated afterwards in physiq !
endif ! of if (startphy_file)
write(*,*) "phyetat0: Surface temperature <tsurf> range:", &
             minval(tsurf), maxval(tsurf)

! Surface emissivity
if (startphy_file) then
  call get_field("emis",emis,found,indextime)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <emis>",1)
  endif
else
  ! If no startfi file, use parameter surfemis in def file
  surfemis=1.0
  call getin_p("surfemis",surfemis)
  print*,"surfemis",surfemis
  emis(:)=surfemis
endif ! of if (startphy_file)
write(*,*) "phyetat0: Surface emissivity <emis> range:", &
             minval(emis), maxval(emis)

! Cloud fraction (added by BC 2010)
if (CLFvarying) then
  if (startphy_file) then
    call get_field("cloudfrac",cloudfrac,found,indextime)
    if (.not.found) then
      call abort_physic(modname,"Failed loading <cloudfrac>",1)
    endif
  else
    cloudfrac(:,:)=0.0
  endif ! of if (startphy_file)
  write(*,*) "phyetat0: Cloud fraction <cloudfrac> range:", &
               minval(cloudfrac), maxval(cloudfrac)
else
  cloudfrac(:,:)=0.0
endif ! of if (CLFvarying)

! Total cloud fraction (added by BC 2010)
if (CLFvarying) then
  if (startphy_file) then
    call get_field("totcloudfrac",totcloudfrac,found,indextime)
    if (.not.found) then
      call abort_physic(modname,"Failed loading <totcloudfrac>",1)
    endif
  else
    totcloudfrac(:)=0.0
  endif ! of if (startphy_file)
  write(*,*) "phyetat0: Total cloud fraction <totcloudfrac> range:", &
             minval(totcloudfrac), maxval(totcloudfrac)
else
  totcloudfrac(:)=0.0
endif ! of if (CLFvarying)

! Height of oceanic ice (added by BC 2010)
if (startphy_file) then
  call get_field("hice",hice,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <hice>"
    !  call abort
    hice(:)=0.
  endif
else
  hice(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: Height of oceanic ice <hice> range:", &
             minval(hice), maxval(hice)

! SLAB OCEAN (added by BC 2014)
if (startphy_file) then
  ! nature of the surface
  call get_field("rnat",rnat,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <rnat>"
        rnat(1:ngrid)=1.
  else
      do ig=1,ngrid
        if((nint(rnat(ig)).eq.2).or.(nint(rnat(ig)).eq.0))then
          rnat(ig)=0.
        else
          rnat(ig)=1.
        endif      
      enddo
  endif ! of if (.not.found)
else
  rnat(:)=1.
endif ! of if (startphy_file)
write(*,*) "phyetat0: Nature of surface <rnat> range:", &
             minval(rnat), maxval(rnat)

if (startphy_file) then
  ! Pourcentage of sea ice cover
  call get_field("pctsrf_sic",pctsrf_sic,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <pctsrf_sic>"
    pctsrf_sic(1:ngrid)=0.
  endif
else
  pctsrf_sic(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: Pourcentage of sea ice cover <pctsrf_sic> range:", &
             minval(pctsrf_sic), maxval(pctsrf_sic)

if (startphy_file) then
  ! Slab ocean temperature (2 layers)
  call get_field("tslab",tslab,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <tslab>"
    do iq=1,noceanmx
      tslab(1:ngrid,iq)=tsurf(1:ngrid)
    enddo
  endif
else
  do iq=1,noceanmx
    tslab(1:ngrid,iq)=tsurf(1:ngrid)
  enddo
endif ! of if (startphy_file)
write(*,*) "phyetat0: Slab ocean temperature <tslab> range:", &
             minval(tslab), maxval(tslab)

if (startphy_file) then
  ! Oceanic ice temperature
  call get_field("tsea_ice",tsea_ice,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <tsea_ice>"
    tsea_ice(1:ngrid)=273.15-1.8
  endif
else
  tsea_ice(1:ngrid)=273.15-1.8
endif ! of if (startphy_file)
write(*,*) "phyetat0: Oceanic ice temperature <tsea_ice> range:", &
             minval(tsea_ice), maxval(tsea_ice)

if (startphy_file) then
  !  Oceanic ice quantity (kg/m^2)
  call get_field("sea_ice",sea_ice,found,indextime)
  if (.not.found) then
    write(*,*) "phyetat0: Failed loading <sea_ice>"
    sea_ice(1:ngrid)=0.
  endif
else
  sea_ice(1:ngrid)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: Oceanic ice quantity <sea_ice> range:", &
             minval(sea_ice), maxval(sea_ice)


! pbl wind variance
if (startphy_file) then
  call get_field("q2",q2,found,indextime)
  if (.not.found) then
    call abort_physic(modname,"Failed loading <q2>",1)
  endif
else
  q2(:,:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: PBL wind variance <q2> range:", &
             minval(q2), maxval(q2)

! tracer on surface
if (nq.ge.1) then
  do iq=1,nq
    txt=noms(iq)
    if (startphy_file) then
      call get_field(txt,qsurf(:,iq),found,indextime)
      if (.not.found) then
        write(*,*) "phyetat0: Failed loading <",trim(txt),">"
        write(*,*) "         ",trim(txt)," is set to zero"
        qsurf(:,iq) = 0.
      endif
    else
      qsurf(:,iq)=0.
    endif ! of if (startphy_file)
    write(*,*) "phyetat0: Surface tracer <",trim(txt),"> range:", &
                 minval(qsurf(:,iq)), maxval(qsurf(:,iq))
  enddo! of do iq=1,nq
endif ! of if (nq.ge.1)


if (startphy_file) then
  ! Call to soil_settings, in order to read soil temperatures,
  ! as well as thermal inertia and volumetric heat capacity
  call soil_settings(nid_start,ngrid,nsoil,tsurf,tsoil,indextime)
endif ! of if (startphy_file)
!
! close file:
!
if (startphy_file) call close_startphy

end subroutine phyetat0

end module phyetat0_mod
