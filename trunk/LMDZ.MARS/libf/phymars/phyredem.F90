module phyredem

implicit none

contains

subroutine physdem0(filename,lonfi,latfi,nsoil,ngrid,nlay,nq, &
                         phystep,day_ini,time,airefi, &
                         alb,ith,pzmea,pzstd,pzsig,pzgam,pzthe, &
                         phmons,psummit,pbase)
! create physics restart file and write time-independent variables
  use comsoil_h, only: inertiedat, volcapa, mlayer
  use geometry_mod, only: cell_area
  use surfdat_h, only: zmea, zstd, zsig, zgam, zthe, &
                       z0_default, albedice, emisice, emissiv, &
                       iceradius, dtemisice, phisfi, z0,   &
                       hmons,summit,base
  use dimradmars_mod, only: tauvis
  use iostart, only : open_restartphy, close_restartphy, & 
                      put_var, put_field, length
  use mod_grid_phy_lmdz, only : klon_glo
  use planete_h, only: aphelie, emin_turb, lmixmin, obliquit, &
                       peri_day, periheli, year_day
  use comcstfi_h, only: g, mugaz, omeg, rad, rcp
  use time_phylmdz_mod, only: daysec
  implicit none
 
  character(len=*), intent(in) :: filename
  real,intent(in) :: lonfi(ngrid)
  real,intent(in) :: latfi(ngrid)
  integer,intent(in) :: nsoil
  integer,intent(in) :: ngrid
  integer,intent(in) :: nlay
  integer,intent(in) :: nq
  real,intent(in) :: phystep
  real,intent(in) :: day_ini
  real,intent(in) :: time
  real,intent(in) :: airefi(ngrid)
  real,intent(in) :: alb(ngrid)
  real,intent(in) :: ith(ngrid,nsoil)
  real,intent(in) :: pzmea(ngrid)
  real,intent(in) :: pzstd(ngrid)
  real,intent(in) :: pzsig(ngrid)
  real,intent(in) :: pzgam(ngrid)
  real,intent(in) :: pzthe(ngrid)
  real,intent(in) :: phmons(ngrid)
  real,intent(in) :: psummit(ngrid)
  real,intent(in) :: pbase(ngrid)

  real :: tab_cntrl(length) ! nb "length=100" defined in iostart module
  
  ! Create physics start file
  call open_restartphy(filename)
  
  ! Build tab_cntrl(:) array
  tab_cntrl(:)=0.0
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Fill control array tab_cntrl(:) with parameters for this run
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Informations on the physics grid
  tab_cntrl(1) = float(klon_glo)  ! total number of nodes on physics grid
  tab_cntrl(2) = float(nlay) ! number of atmospheric layers
  tab_cntrl(3) = day_ini + int(time)      ! initial day 
  tab_cntrl(4) = time -int(time)          ! initial time of day

  ! Informations about Mars, used by dynamics and physics
  tab_cntrl(5) = rad      ! radius of Mars (m) ~3397200
  tab_cntrl(6) = omeg     ! rotation rate (rad.s-1)
  tab_cntrl(7) = g        ! gravity (m.s-2) ~3.72
  tab_cntrl(8) = mugaz    ! Molar mass of the atmosphere (g.mol-1) ~43.49
  tab_cntrl(9) = rcp      !  = r/cp  ~0.256793 (=kappa dans dynamique)
  tab_cntrl(10) = daysec  ! length of a sol (s)  ~88775

  tab_cntrl(11) = phystep  ! time step in the physics
  tab_cntrl(12) = 0.
  tab_cntrl(13) = 0.

  ! Informations about Mars, only for physics
  tab_cntrl(14) = year_day  ! length of year (sols) ~668.6
  tab_cntrl(15) = periheli  ! min. Sun-Mars distance (Mkm) ~206.66
  tab_cntrl(16) = aphelie   ! max. SUn-Mars distance (Mkm) ~249.22
  tab_cntrl(17) = peri_day  ! date of perihelion (sols since N. spring)
  tab_cntrl(18) = obliquit  ! Obliquity of the planet (deg) ~23.98

  ! Boundary layer and turbulence
  tab_cntrl(19) = z0_default   ! default surface roughness (m) ~0.01
  tab_cntrl(20) = lmixmin   ! mixing length ~100
  tab_cntrl(21) = emin_turb ! minimal energy ~1.e-8

  ! Optical properties of polar caps and ground emissivity
  tab_cntrl(22) = albedice(1)  ! Albedo of northern cap ~0.5
  tab_cntrl(23) = albedice(2)  ! Albedo of southern cap ~0.5
  tab_cntrl(24) = emisice(1)   ! Emissivity of northern cap ~0.95
  tab_cntrl(25) = emisice(2)   ! Emissivity of southern cap ~0.95
  tab_cntrl(26) = emissiv      ! Emissivity of martian soil ~.95
  tab_cntrl(31) = iceradius(1) ! mean scat radius of CO2 snow (north)
  tab_cntrl(32) = iceradius(2) ! mean scat radius of CO2 snow (south)
  tab_cntrl(33) = dtemisice(1) ! time scale for snow metamorphism (north)
  tab_cntrl(34) = dtemisice(2) ! time scale for snow metamorphism (south)

  ! dust aerosol properties
  tab_cntrl(27) = tauvis      ! mean visible optical depth

  ! Soil properties:
  tab_cntrl(35) = volcapa ! soil volumetric heat capacity

  ! Write the controle array
  call put_var("controle","Control parameters",tab_cntrl)
  
  ! Write the mid-layer depths
  call put_var("soildepth","Soil mid-layer depth",mlayer)
  
  ! Write longitudes
  call put_field("longitude","Longitudes of physics grid",lonfi)
  
  ! Write latitudes
  call put_field("latitude","Latitudes of physics grid",latfi)
  
  ! Write mesh areas
  call put_field("area","Mesh area",cell_area)
  
  ! Write surface geopotential
  call put_field("phisfi","Geopotential at the surface",phisfi)
  
  ! Write surface albedo
  call put_field("albedodat","Albedo of bare ground",alb)
  
  ! Subgrid topogaphy variables
  call put_field("ZMEA","Relief: mean relief",zmea)
  call put_field("ZSTD","Relief: standard deviation",zstd)
  call put_field("ZSIG","Relief: sigma parameter",zsig)
  call put_field("ZGAM","Relief: gamma parameter",zgam)
  call put_field("ZTHE","Relief: theta parameter",zthe)
  call put_field("hmons","Relief: hmons parameter (summit - base)",hmons)
  call put_field("summit","Relief: altitude from the aeroid of the summit of the highest subgrid topography",summit)
  call put_field("base","Relief: altitude from the aeroid of the base of the highest subgrid topography",base)
     
  ! Soil thermal inertia
  call put_field("inertiedat","Soil thermal inertia",ith)
  
  ! Surface roughness length
  call put_field("z0","Surface roughness length",z0)
  
  ! Close file
  call close_restartphy
  
end subroutine physdem0

subroutine physdem1(filename,nsoil,ngrid,nlay,nq, &
                    phystep,time,tsurf,tsoil,co2ice,albedo,emis,q2,qsurf,&
                    tauscaling,totcloudfrac,wstar, &
                    mem_Mccn_co2,mem_Nccn_co2,mem_Mh2o_co2, watercap)
  ! write time-dependent variable to restart file
  use iostart, only : open_restartphy, close_restartphy, & 
                      put_var, put_field
  use tracer_mod, only: noms ! tracer names
  use nonoro_gwd_ran_mod, only: du_nonoro_gwd, dv_nonoro_gwd
  use compute_dtau_mod, only: dtau

  implicit none
  
  include "callkeys.h"
  
  character(len=*),intent(in) :: filename
  integer,intent(in) :: nsoil
  integer,intent(in) :: ngrid
  integer,intent(in) :: nlay
  integer,intent(in) :: nq
  real,intent(in) :: phystep
  real,intent(in) :: time
  real,intent(in) :: tsurf(ngrid)
  real,intent(in) :: tsoil(ngrid,nsoil)
  real,intent(in) :: co2ice(ngrid)
  real,intent(in) :: albedo(ngrid,2)
  real,intent(in) :: emis(ngrid)
  real,intent(in) :: q2(ngrid,nlay+1)
  real,intent(in) :: qsurf(ngrid,nq)
  real,intent(in) :: tauscaling(ngrid)
  real,intent(in) :: totcloudfrac(ngrid)
  real,intent(in) :: wstar(ngrid)
  real,intent(in) :: mem_Mccn_co2(ngrid,nlay) ! CCN mass of H2O and dust used by CO2
  real,intent(in) :: mem_Nccn_co2(ngrid,nlay) ! CCN number of H2O and dust used by CO2
  real,intent(in) :: mem_Mh2o_co2(ngrid,nlay) ! H2O mass integred into CO2 crystal
  real,intent(in) :: watercap(ngrid)
  
  integer :: iq
  character(len=30) :: txt ! to store some text
  ! indexes of water vapour & water ice tracers (if any):
  integer :: i_h2o_vap=0
  integer :: i_h2o_ice=0

  
  ! Open file
  call open_restartphy(filename)
  
  ! First variable to write must be "Time", in order to correctly
  ! set time counter in file
  call put_var("Time","Temps de simulation",time)
  
  ! CO2 ice layer
  call put_field("co2ice","CO2 ice cover",co2ice,time)

  ! Water ice layer
  call put_field("watercap","H2O ice cover",watercap,time)
  
  ! Surface temperature
  call put_field("tsurf","Surface temperature",tsurf,time)
  
  ! Soil temperature
  call put_field("tsoil","Soil temperature",tsoil,time)
  
  ! Albedo of the surface
  call put_field("albedo","Surface albedo",albedo(:,1),time)
  
  ! Emissivity of the surface
  call put_field("emis","Surface emissivity",emis,time)
  
  ! Planetary Boundary Layer
  call put_field("q2","pbl wind variance",q2,time)

  ! Sub-grid cloud fraction
  call put_field("totcloudfrac","Total cloud fraction",totcloudfrac,time) 
 
  ! Dust conversion factor
  ! Only to be read by newstart to convert to actual dust values
  ! Or by any user who wants to reconstruct dust, opacity from the start files.
  call put_field("tauscaling","dust conversion factor",tauscaling,time)

  if (dustinjection.gt.0) then
    call put_field("dtau","dust opacity difference between GCM and scenario",&
                   dtau,time)
  endif

  if (calltherm) then
    call put_field("wstar","Max vertical velocity in thermals",wstar,time)
  endif

  ! Tracers on the surface
  ! preliminary stuff: look for water vapour & water ice tracers (if any)
  do iq=1,nq
    if (noms(iq).eq."h2o_vap") then
      i_h2o_vap=iq
    endif
    if (noms(iq).eq."h2o_ice") then
      i_h2o_ice=iq
    endif
  enddo
  
  if (nq.gt.0) then
    do iq=1,nq
      txt=noms(iq)
      ! Exception: there is no water vapour surface tracer
      if (txt.eq."h2o_vap") then
        write(*,*)"physdem1: skipping water vapour tracer"
        if (i_h2o_ice.eq.i_h2o_vap) then
          ! then there is no "water ice" tracer; but still
          ! there is some water ice on the surface
          write(*,*)"          writing water ice instead"
          txt="h2o_ice"
        else
          ! there is a "water ice" tracer which has been / will be
          ! delt with in due time
          cycle
        endif ! of if (igcm_h2o_ice.eq.igcm_h2o_vap)
      endif ! of if (txt.eq."h2o_vap")
      call put_field(trim(txt),"tracer on surface",qsurf(:,iq),time)
    enddo
  endif
  ! Memory of the origin of the co2 particles
  if (co2useh2o) then
     call put_field("mem_Mccn_co2","CCN mass of H2O and dust used by CO2",mem_Mccn_co2,time)
     call put_field("mem_Nccn_co2","CCN number of H2O and dust used by CO2",mem_Nccn_co2,time)
     call put_field("mem_Mh2o_co2","H2O mass integred into CO2 crystal",mem_Mh2o_co2,time)
  endif
 
  ! Non-orographic gavity waves
  if (calllott_nonoro) then
     call put_field("du_nonoro_gwd","Zonal wind tendency due to GW",du_nonoro_gwd,time)
     call put_field("dv_nonoro_gwd","Meridional wind tendency due to GW",dv_nonoro_gwd,time)
  endif
  ! Close file
  call close_restartphy
  
end subroutine physdem1

end module phyredem
