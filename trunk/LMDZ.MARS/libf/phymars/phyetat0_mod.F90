module phyetat0_mod

implicit none

contains

subroutine phyetat0 (fichnom,tab0,Lmodif,nsoil,ngrid,nlay,nq, &
                     day_ini,time0,tsurf,tsoil,albedo,emis,q2,qsurf,co2ice, &
                     tauscaling,totcloudfrac,wstar,mem_Mccn_co2,mem_Nccn_co2,&
                     mem_Mh2o_co2,watercap)

  use tracer_mod, only: noms ! tracer names
  use surfdat_h, only: phisfi, albedodat, z0, z0_default,&
                       zmea, zstd, zsig, zgam, zthe, hmons, summit, base
  use iostart, only: nid_start, open_startphy, close_startphy, &
                     get_field, get_var, inquire_field, &
                     inquire_dimension, inquire_dimension_length
  use nonoro_gwd_ran_mod, only: du_nonoro_gwd, dv_nonoro_gwd
  use compute_dtau_mod, only: dtau

  USE ioipsl_getin_p_mod, ONLY : getin_p

  implicit none
  
  include "callkeys.h"
!======================================================================
! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
!  Adaptation à Mars : Yann Wanherdrick 
! Objet: Lecture de l etat initial pour la physique
! Modifs: Aug.2010 EM : use NetCDF90 to load variables (enables using
!                      r4 or r8 restarts independently of having compiled
!                      the GCM in r4 or r8)
!         June 2013 TN : Possibility to read files with a time axis
!         November 2013 EM : Enabeling parallel, using iostart module
!         March 2020 AD: Enabling initialization of physics without startfi
!                        flag: startphy_file
!======================================================================
  INTEGER nbsrf !Mars nbsrf a 1 au lieu de 4
  PARAMETER (nbsrf=1) ! nombre de sous-fractions pour une maille
!======================================================================
!  Arguments:
!  ---------
!  inputs:
!  logical,intent(in) :: startphy_file ! .true. if reading start file
  character*(*),intent(in) :: fichnom ! "startfi.nc" file
  integer,intent(in) :: tab0
  integer,intent(in) :: Lmodif
  integer,intent(in) :: nsoil ! # of soil layers
  integer,intent(in) :: ngrid ! # of atmospheric columns
  integer,intent(in) :: nlay ! # of atmospheric layers
  integer,intent(in) :: nq
  integer :: day_ini
  real :: time0

!  outputs:
  real,intent(out) :: tsurf(ngrid) ! surface temperature
  real,intent(out) :: tsoil(ngrid,nsoil) ! soil temperature
  real,intent(out) :: albedo(ngrid,2) ! surface albedo
  real,intent(out) :: emis(ngrid) ! surface emissivity
  real,intent(out) :: q2(ngrid,nlay+1) ! 
  real,intent(out) :: qsurf(ngrid,nq) ! tracers on surface
  real,intent(out) :: co2ice(ngrid) ! co2 ice cover
  real,intent(out) :: tauscaling(ngrid) ! dust conversion factor
  real,intent(out) :: totcloudfrac(ngrid) ! total cloud fraction
  real,intent(out) :: wstar(ngrid) ! Max vertical velocity in thermals (m/s)
  real,intent(out) :: mem_Mccn_co2(ngrid,nlay) ! Memory of CCN mass of H2O and dust used by CO2
  real,intent(out) :: mem_Nccn_co2(ngrid,nlay) ! Memory of CCN number of H2O and dust used by CO2
  real,intent(out) :: mem_Mh2o_co2(ngrid,nlay) ! Memory of H2O mass integred into CO2 crystal
  real,intent(out) :: watercap(ngrid) ! h2o_ice_cover
!======================================================================
!  Local variables:

      real surffield(ngrid) ! to temporarily store a surface field
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
      REAL p_rad,p_omeg,p_g,p_mugaz,p_daysec
      INTEGER nqold

! flag which identifies if 'startfi.nc' file is using old names (qsurf01,...)
      logical :: oldtracernames=.false.
      integer :: count
      character(len=30) :: txt ! to store some text

! specific for time 
      REAL,ALLOCATABLE :: time(:) ! times stored in start
      INTEGER timelen ! number of times stored in the file
      INTEGER indextime ! index of selected time

      INTEGER :: edges(3),corner(3)
      LOGICAL :: found

      REAL :: timestart ! to pick which initial state to start from
      REAL :: surfemis  ! constant emissivity when no startfi
      REAL :: surfalbedo  ! constant emissivity when no startfi
      CHARACTER(len=5) :: modname="phyetat0"

write(*,*) "phyetat0: startphy_file", startphy_file

if(.not.startphy_file) then
   surfemis = -999
   surfalbedo = -999
   call getin_p("surfemis",surfemis)
   call getin_p("surfalbedo",surfalbedo)
   if (surfemis.eq.-999) then
      call abort_physic(modname, &
               "Missing value for surfemis in def files!",1)
   endif
   if (surfalbedo.eq.-999) then
      call abort_physic(modname, &
               "Missing value for surfalbedo in def files!",1)
   endif
   write(*,*) "phyetat0: surfemis: ", surfemis
   write(*,*) "phyetat0: surfalbedo: ", surfalbedo
endif

if (startphy_file) then
   ! open physics initial state file:
   call open_startphy(fichnom)
   ! possibility to modify tab_cntrl in tabfi
   write(*,*)
   write(*,*) 'TABFI in phyeta0: Lmodif=',Lmodif," tab0=",tab0
   call tabfi (nid_start,Lmodif,tab0,day_ini,lmax,p_rad, &
               p_omeg,p_g,p_mugaz,p_daysec,time0)
else ! "academic" initialization of planetary parameters via tabfi
   call tabfi (0,0,0,day_ini,lmax,p_rad, &
               p_omeg,p_g,p_mugaz,p_daysec,time0)
endif ! of if (startphy_file)

if (startphy_file) then
   ! Load surface geopotential:
   call get_field("phisfi",phisfi,found)
   if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <phisfi>",1)
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <albedodat>",1)
   endif
else ! If no startfi file, use parameter surfalbedo in def file
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <ZMEA>",1)
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <ZSTD>",1)
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <ZSIG>",1)
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <ZGAM>",1)
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
     call abort_physic(modname, &
                "phyetat0: Failed loading <ZTHE>",1)
   endif
else
  zthe(:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: <ZTHE> range:", &
             minval(zthe), maxval(zthe)

! hmons
if (startphy_file) then
   call get_field("hmons",hmons,found)
   if (.not.found) then
     write(*,*) "WARNING: phyetat0: Failed loading <hmons>"
     if (rdstorm) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <hmons>",1)
     else
       write(*,*) "will continue anyway..."
       write(*,*) "because you may not need it."
       hmons(:)=0.
     end if ! if (rdstorm)
   else
     do ig=1,ngrid
       if (hmons(ig) .eq. -999999.)  hmons(ig)=0.
     enddo
   endif ! (.not.found)
else
   hmons(:)=0.
endif ! if (startphy_file)
write(*,*) "phyetat0: <hmons> range:", &
            minval(hmons), maxval(hmons)


! summit
if (startphy_file) then
   call get_field("summit",summit,found)
   if (.not.found) then
     write(*,*) "WARNING: phyetat0: Failed loading <summit>"
     if (rdstorm) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <summit>",1)
     else
       write(*,*) "will continue anyway..."
       write(*,*) "because you may not need it."
       summit(:)=0.
     end if
   else
     do ig=1,ngrid
       if (summit(ig) .eq. -999999.)  summit(ig)=0.
     enddo
   endif ! if (.not.found)
else
   summit(:)=0.  
endif ! if (startphy_file)
write(*,*) "phyetat0: <summit> range:", &
            minval(summit), maxval(summit)


! base 
if (startphy_file) then
   call get_field("base",base,found)
   if (.not.found) then
     write(*,*) "WARNING: phyetat0: Failed loading <base>"
     if (rdstorm) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <base>",1)
     else
       write(*,*) "will continue anyway..."
       write(*,*) "because you may not need it."
       base(:)=0.
     end if
   else
     do ig=1,ngrid
       if (base(ig) .eq. -999999.)  base(ig)=0.
     enddo
   endif ! if(.not.found)
else
   base(:)=0.
endif ! if (startphy_file)
write(*,*) "phyetat0: <base> range:", &
            minval(base), maxval(base)

 
! Time axis
! obtain timestart from run.def
timestart=-9999 ! default value
call getin_p("timestart",timestart)
if (startphy_file) then
   found=inquire_dimension("Time")
   if (.not.found) then
     indextime = 1
     write(*,*) "phyetat0: No time axis found in "//trim(fichnom)
   else
     write(*,*) "phyetat0: Time axis found in "//trim(fichnom)
     timelen=inquire_dimension_length("Time")
     allocate(time(timelen))
     ! load "Time" array:
     call get_var("Time",time,found)
     if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <Time>",1)
     endif
     ! seclect the desired time index
     IF (timestart .lt. 0) THEN  ! default: we use the last time value
       indextime = timelen
     ELSE  ! else we look for the desired value in the time axis
       indextime = 0
       DO i=1,timelen
         IF (abs(time(i) - timestart) .lt. 0.01) THEN
           indextime = i
           EXIT
         ENDIF
       ENDDO
       IF (indextime .eq. 0) THEN
         PRINT*, "Time", timestart," is not in "//trim(fichnom)//"!!"
         PRINT*, "Stored times are:"
         DO i=1,timelen
            PRINT*, time(i)
         ENDDO
         call abort_physic(modname,"phyetat0: Time error",1) 
       ENDIF
     ENDIF ! of IF (timestart .lt. 0)
     ! In startfi the absolute date is day_ini + time0 + time
     ! For now on, in the GCM physics, it is day_ini + time0
     time0 = time(indextime) + time0
     day_ini = day_ini + INT(time0)
     time0 = time0 - INT(time0)    
     PRINT*, "phyetat0: Selected time ",time(indextime), &
             " at index ",indextime
     DEALLOCATE(time)
   endif ! of if Time not found in file
else
  indextime = 1
endif ! if (startphy_file)

! CO2 ice cover
if (startphy_file) then
   call get_field("co2ice",co2ice,found,indextime)
   if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <co2ice>",1)
   endif
else
   co2ice(:)=0.
endif !if (startphy_file)
write(*,*) "phyetat0: CO2 ice cover <co2ice> range:", &
            minval(co2ice), maxval(co2ice)

! Memory of the origin of the co2 particles
if (startphy_file) then
   call get_field("mem_Mccn_co2",mem_Mccn_co2,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <mem_Mccn_co2> not in file"
     mem_Mccn_co2(:,:)=0
   endif
else
   mem_Mccn_co2(:,:)=0
endif !if (startphy_file)
write(*,*) "phyetat0: Memory of CCN mass of H2O and dust used by CO2"
write(*,*) " <mem_Mccn_co2> range:", &
             minval(mem_Mccn_co2), maxval(mem_Mccn_co2)

if (startphy_file) then
   call get_field("mem_Nccn_co2",mem_Nccn_co2,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <mem_Nccn_co2> not in file"
     mem_Nccn_co2(:,:)=0
   endif
else
  mem_Nccn_co2(:,:)=0
endif ! if (startphy_file)
write(*,*) "phyetat0: Memory of CCN number of H2O and dust used by CO2"
write(*,*) " <mem_Nccn_co2> range:", &
             minval(mem_Nccn_co2), maxval(mem_Nccn_co2)

if (startphy_file) then
   call get_field("mem_Mh2o_co2",mem_Mh2o_co2,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <mem_Mh2o_co2> not in file"
     mem_Mh2o_co2(:,:)=0
   endif
else
   mem_Mh2o_co2(:,:)=0
endif ! if (startphy_file)
write(*,*) "phyetat0: Memory of H2O mass integred into CO2 crystal"
write(*,*) " <mem_Mh2o_co2> range:", &
             minval(mem_Mh2o_co2), maxval(mem_Mh2o_co2)

! Dust conversion factor
if (startphy_file) then
   call get_field("tauscaling",tauscaling,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <tauscaling> not in file"
     tauscaling(:) = 1
   endif
else
   tauscaling(:) = 1
endif ! if (startphy_file)
write(*,*) "phyetat0: dust conversion factor <tauscaling> range:", &
            minval(tauscaling), maxval(tauscaling)


! dtau: opacity difference between GCM and dust scenario
if (startphy_file) then
   call get_field("dtau",dtau,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <dtau> not in file; set to zero"
     dtau(:) = 0
   endif
else
   dtau(:)= 0
endif ! if (startphy_file)
write(*,*) "phyetat0: opacity diff wrt scenario <dtau> range:", &
            minval(dtau), maxval(dtau)


! Sub-grid cloud fraction
if (startphy_file) then
   call get_field("totcloudfrac",totcloudfrac,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <totcloudfrac> not in file WARNING put to 1"
     totcloudfrac(:) = 1.0 !valeur par defaut (CLFfixval par defaut)
   endif
else
   totcloudfrac(:)=1.0
endif ! if (startphy_file)
write(*,*) "phyetat0: total cloud fraction <totcloudfrac> range:", &
            minval(totcloudfrac), maxval(totcloudfrac)


! Max vertical velocity in thermals
if (startphy_file) then
   call get_field("wstar",wstar,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: <wstar> not in file! Set to zero"
     wstar(:)=0
   endif
else
   wstar(:)=0
endif ! if (startphy_file)
write(*,*) "phyetat0: Max vertical velocity in thermals <wstar> range:", &
            minval(wstar),maxval(wstar)


! Surface temperature :
if (startphy_file) then !tsurf
   call get_field("tsurf",tsurf,found,indextime)
   if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <tsurf>",1)
   endif
else
  tsurf(:)=0. ! will be updated afterwards in physiq !
endif ! of if (startphy_file)
write(*,*) "phyetat0: Surface temperature <tsurf> range:", &
            minval(tsurf), maxval(tsurf)

! Surface albedo
if (startphy_file) then
   call get_field("albedo",albedo(:,1),found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: Failed loading <albedo>"
     albedo(:,1)=albedodat(:)
   endif
else
   albedo(:,1)=albedodat(:)
endif ! of if (startphy_file)
write(*,*) "phyetat0: Surface albedo <albedo> range:", &
            minval(albedo(:,1)), maxval(albedo(:,1))
albedo(:,2)=albedo(:,1)

! Surface emissivity
if (startphy_file) then
   call get_field("emis",emis,found,indextime)
   if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <emis>",1)
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


! surface roughness length (NB: z0 is a common in surfdat_h)
if (startphy_file) then
   call get_field("z0",z0,found)
   if (.not.found) then
     write(*,*) "phyetat0: Failed loading <z0>"
     write(*,*) 'will use constant value of z0_default:',z0_default
     z0(:)=z0_default
   endif
else
   z0(:)=z0_default
endif ! of if (startphy_file)
write(*,*) "phyetat0: Surface roughness <z0> range:", &
            minval(z0), maxval(z0)


! pbl wind variance
if (startphy_file) then
   call get_field("q2",q2,found,indextime)
   if (.not.found) then
     call abort_physic(modname, &
                "phyetat0: Failed loading <q2>",1)
   endif
else
  q2(:,:)=0.
endif ! of if (startphy_file)
write(*,*) "phyetat0: PBL wind variance <q2> range:", &
            minval(q2), maxval(q2)

! Non-orographic gravity waves
if (startphy_file) then 
   call get_field("du_nonoro_gwd",du_nonoro_gwd,found,indextime)
   if (.not.found) then
      write(*,*) "phyetat0: <du_nonoro_gwd> not in file"
      du_nonoro_gwd(:,:)=0.
   endif
else
du_nonoro_gwd(:,:)=0.
endif ! if (startphy_file)
write(*,*) "phyetat0: Memory of zonal wind tendency due to non-orographic GW"
write(*,*) " <du_nonoro_gwd> range:", &
             minval(du_nonoro_gwd), maxval(du_nonoro_gwd)

if (startphy_file) then 
   call get_field("dv_nonoro_gwd",dv_nonoro_gwd,found,indextime)
   if (.not.found) then
      write(*,*) "phyetat0: <dv_nonoro_gwd> not in file"
      dv_nonoro_gwd(:,:)=0.
   endif
else ! ! if (startphy_file)
dv_nonoro_gwd(:,:)=0.
endif ! if (startphy_file)
write(*,*) "phyetat0: Memory of meridional wind tendency due to non-orographic GW"
write(*,*) " <dv_nonoro_gwd> range:", &
             minval(dv_nonoro_gwd), maxval(dv_nonoro_gwd)

! tracer on surface
if (nq.ge.1) then
  do iq=1,nq
    txt=noms(iq)
    if (txt.eq."h2o_vap") then
      ! There is no surface tracer for h2o_vap;
      ! "h2o_ice" should be loaded instead
      txt="h2o_ice"
      write(*,*) 'phyetat0: loading surface tracer', &
                           ' h2o_ice instead of h2o_vap'
      write(*,*) 'iq = ', iq
    endif
    if (startphy_file) then
       call get_field(txt,qsurf(:,iq),found,indextime)
       if (.not.found) then
         write(*,*) "phyetat0: Failed loading <",trim(txt),">"
         write(*,*) "         ",trim(txt)," is set to zero"
       endif
    else
      qsurf(:,iq)=0.
    endif ! of if (startphy_file)
    write(*,*) "phyetat0: Surface tracer <",trim(txt),"> range:", &
                 minval(qsurf(:,iq)), maxval(qsurf(:,iq))
  enddo ! of do iq=1,nq
endif ! of if (nq.ge.1)

if (startphy_file) then
   call get_field("watercap",watercap,found,indextime)
   if (.not.found) then
     write(*,*) "phyetat0: Failed loading <watercap> : ", &
                          "<watercap> is set to zero"
     watercap(:)=0

     write(*,*) 'Now transfer negative surface water ice to', &
                ' watercap !'
     if (nq.ge.1) then
      do iq=1,nq
       txt=noms(iq)
       if (txt.eq."h2o_ice") then
         do ig=1,ngrid
          if (qsurf(ig,iq).lt.0.0) then
             watercap(ig) = qsurf(ig,iq)
             qsurf(ig,iq) = 0.0
          end if
         end do
       endif
      enddo
     endif ! of if (nq.ge.1)
   endif ! of if (.not.found)
else
   watercap(:)=0
endif ! of if (startphy_file) 
write(*,*) "phyetat0: Surface water ice <watercap> range:", &
            minval(watercap), maxval(watercap)
 


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
