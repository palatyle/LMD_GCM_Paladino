module phyetat0_mod

implicit none

contains

subroutine phyetat0 (startphy_file, &
                     ngrid,nlayer,fichnom,tab0,Lmodif,nsoil,nq, &
                     day_ini,time,tsurf,tsoil, &
                     emis,q2,qsurf,tankCH4)

! to use  'getin_p'
      use ioipsl_getin_p_mod, only: getin_p

  use tabfi_mod, only: tabfi
  USE tracer_h, ONLY: noms
  USE surfdat_h, only: phisfi, albedodat, zmea, zstd, zsig, zgam, zthe
  use iostart, only: nid_start, open_startphy, close_startphy, &
                     get_field, get_var, inquire_field, &
                     inquire_dimension, inquire_dimension_length
  use callkeys_mod, only: surfalbedo,surfemis

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
  real,intent(out) :: tankCH4(ngrid)  ! depth of CH4 tank

!======================================================================
!  Local variables:

!      INTEGER radpas
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

! Depth of methane tank (added by JVO 2017)
 if (startphy_file) then
   call get_field("tankCH4",tankCH4,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: Failed loading <tankCH4>"
     !  call abort
     tankCH4(:)=2.
   endif
 else
   tankCH4(:)=2.
 endif ! of if (startphy_file)
 write(*,*) "phyetat0: Depth of methane tank <tankCH4> range:", &
              minval(tankCH4), maxval(tankCH4)

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

! Upper chemistry
if (startphy_file) then
  ! Call to soil_settings, in order to read upper chemistry
  ! pressure grid as well as composition fields
  call chem_settings(nid_start,ngrid,nlayer,indextime)
endif ! of if (startphy_file)

!
! close file:
!
if (startphy_file) call close_startphy

end subroutine phyetat0

end module phyetat0_mod
