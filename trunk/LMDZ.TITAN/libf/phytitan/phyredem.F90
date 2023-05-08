module phyredem

implicit none

contains

subroutine physdem0(filename,lonfi,latfi,nsoil,ngrid,nlay,nq, &
                         phystep,day_ini,time,airefi, &
                         alb,ith,pzmea,pzstd,pzsig,pzgam,pzthe)
! create physics restart file and write time-independent variables
  use comchem_h, only: preskim
  use comsoil_h, only: volcapa, mlayer
  use geometry_mod, only: cell_area
  use surfdat_h, only: zmea, zstd, zsig, zgam, zthe, &
                       emisice, emissiv,             &
                       iceradius, dtemisice, phisfi
  use iostart, only : open_restartphy, close_restartphy, & 
                      put_var, put_field, length
  use mod_grid_phy_lmdz, only : klon_glo
  use planete_mod, only: year_day, periastr, apoastr, peri_day, &
                         obliquit, z0, lmixmin, emin_turb
  use comcstfi_mod, only: rad, omeg, g, mugaz, rcp
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
 
  real :: tab_cntrl(length) ! nb "length=100" defined in iostart module
  
  ! Create physics start file
  call open_restartphy(filename)

! tab_cntrl() contains run parameters
  tab_cntrl(:)=0 ! initialization
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Fill control array tab_cntrl(:) with paramleters for this run
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Informations on the physics grid
  tab_cntrl(1) = float(klon_glo)  ! number of nodes on physics grid
  tab_cntrl(2) = float(nlay) ! number of atmospheric layers
  tab_cntrl(3) = day_ini + int(time)         ! final day 
  tab_cntrl(4) = time -int(time)            ! final time of day

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
  tab_cntrl(15) = periastr  ! min. star-planet distance (AU)
  tab_cntrl(16) = apoastr   ! max. star-planet distance (AU) 
  tab_cntrl(17) = peri_day  ! date of periastron (sols since N. spring)
  tab_cntrl(18) = obliquit  ! Obliquity of the planet (deg) ~23.98

! Boundary layer and turbulence
  tab_cntrl(19) = z0        ! surface roughness (m) ~0.01
  tab_cntrl(20) = lmixmin   ! mixing length ~100
  tab_cntrl(21) = emin_turb ! minimal energy ~1.e-8

! Optical properties of polar caps and ground emissivity
  tab_cntrl(24) = emisice(1)   ! Emissivity of northern cap ~0.95
  tab_cntrl(25) = emisice(2)   ! Emissivity of southern cap ~0.95
  tab_cntrl(26) = emissiv      ! Emissivity of martian soil ~.95
  tab_cntrl(31) = iceradius(1) ! mean scat radius of CO2 snow (north)
  tab_cntrl(32) = iceradius(2) ! mean scat radius of CO2 snow (south)
  tab_cntrl(33) = dtemisice(1) ! time scale for snow metamorphism (north)
  tab_cntrl(34) = dtemisice(2) ! time scale for snow metamorphism (south)

  tab_cntrl(28) = 0. 
  tab_cntrl(29) = 0.
  tab_cntrl(30) = 0.

! Soil properties:
  tab_cntrl(35) = volcapa ! soil volumetric heat capacity

  call put_var("controle","Control parameters",tab_cntrl)
  
  ! Write the mid-layer depths
  call put_var("soildepth","Soil mid-layer depth",mlayer)
  
  ! Write the mid-layer upper chemistry pressure
  call put_var("preskim","Upper chemistry mid-layer pressure",preskim)
  
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
  
  ! Soil thermal inertia
  call put_field("inertiedat","Soil thermal inertia",ith)
  
  ! Close file
  call close_restartphy
  
end subroutine physdem0

subroutine physdem1(filename,nsoil,ngrid,nlay,nq, &
                    phystep,time,tsurf,tsoil,emis,q2,qsurf,tankCH4)
  ! write time-dependent variable to restart file
  use iostart, only : open_restartphy, close_restartphy, & 
                      put_var, put_field
  use comchem_h, only: nkim, cnames, ykim_up
  use tracer_h, only: noms
  use callkeys_mod, only: callchim

  implicit none

  character(len=*),intent(in) :: filename
  integer,intent(in) :: nsoil
  integer,intent(in) :: ngrid
  integer,intent(in) :: nlay
  integer,intent(in) :: nq
  real,intent(in) :: phystep
  real,intent(in) :: time
  real,intent(in) :: tsurf(ngrid)
  real,intent(in) :: tsoil(ngrid,nsoil)
  real,intent(in) :: emis(ngrid)
  real,intent(in) :: q2(ngrid,nlay+1)
  real,intent(in) :: qsurf(ngrid,nq)
  real,intent(in) :: tankCH4(ngrid)

  integer :: iq
  
  ! Open file
  call open_restartphy(filename)

  ! First variable to write must be "Time", in order to correctly
  ! set time counter in file
  !call put_var("Time","Temps de simulation",time)
  
  ! Surface temperature
  call put_field("tsurf","Surface temperature",tsurf)
  
  ! Soil temperature
  call put_field("tsoil","Soil temperature",tsoil)
  
  ! Emissivity of the surface
  call put_field("emis","Surface emissivity",emis)
  
  ! Planetary Boundary Layer
  call put_field("q2","pbl wind variance",q2)

  ! Methane tank depth
  call put_field("tankCH4","Depth of methane tank",tankCH4)
  
  ! Tracers
  if (nq>0) then
    do iq=1,nq
      call put_field(noms(iq),"tracer on surface",qsurf(:,iq))
    enddo
  endif ! of if (nq>0)
  
  ! Upper chemistry
  if (callchim) then
    do iq=1,nkim
      call put_field(trim(cnames(iq))//"_up",trim(cnames(iq))//" in upper atmosphere",ykim_up(iq,:,:))
    enddo
  endif ! of if callchim
  
! close file
      CALL close_restartphy
!$OMP BARRIER

end subroutine physdem1

end module phyredem
