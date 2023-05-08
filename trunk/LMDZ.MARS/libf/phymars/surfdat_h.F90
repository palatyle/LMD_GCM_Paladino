module surfdat_h

  ! arrays are allocated in conf_phys
  real,save,allocatable :: albedodat(:) ! albedo of bare ground
  real,save,allocatable :: phisfi(:) ! geopotential at ground level
  real,save :: albedice(2) ! default albedo for ice (1: North H. 2: South H.)
  real,save :: emisice(2) ! ice emissivity; 1:Northern hemisphere 2:Southern hemisphere
  real,save :: emissiv ! emissivity of bare ground
  logical,save :: TESicealbedo ! use TES ice cap albedoes (if set to .true.)
  logical,save,allocatable :: watercaptag(:) ! flag for water ice surface
  real, save, allocatable :: dryness(:)
      
  logical,save :: temptag !temp tag for water caps
      
  real,save :: albedo_h2o_ice ! water ice albedo
  real,save :: inert_h2o_ice ! water ice thermal inertia
  real,save :: frost_albedo_threshold ! water frost thickness on the ground (kg.m^-2, ie mm)
  real,save :: TESice_Ncoef ! coefficient for TES ice albedo in Northern hemisphere
  real,save :: TESice_Scoef ! coefficient for TES ice albedo in Southern hemisphere
  real,save :: iceradius(2) , dtemisice(2)
  real,save,allocatable :: zmea(:),zstd(:),zsig(:),zgam(:),zthe(:)
  real,save,allocatable :: hmons(:),summit(:),base(:)
  real,save,allocatable :: z0(:) ! surface roughness length (m)
  real,save :: z0_default ! default (constant over planet) surface roughness (m)

  !! variables
  REAL,SAVE,ALLOCATABLE :: tsurf(:)   ! Surface temperature (K)
  REAL,SAVE,ALLOCATABLE :: co2ice(:)  ! co2 ice surface layer (kg.m-2)  
  REAL,SAVE,ALLOCATABLE :: emis(:)    ! Thermal IR surface emissivity
  REAL,SAVE,ALLOCATABLE :: capcal(:) ! surface heat capacity (J m-2 K-1)
  REAL,SAVE,ALLOCATABLE :: fluxgrd(:) ! surface conduction flux (W.m-2)
  REAL,ALLOCATABLE,SAVE :: qsurf(:,:) ! tracer on surface (e.g. kg.m-2)
  REAL,SAVE,ALLOCATABLE :: watercap(:) ! Surface water ice (kg.m-2)

contains

  subroutine ini_surfdat_h(ngrid,nq)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nq ! number of tracers  

    allocate(albedodat(ngrid))
    allocate(phisfi(ngrid))
    allocate(watercaptag(ngrid))
    allocate(dryness(ngrid))
    allocate(zmea(ngrid))
    allocate(zstd(ngrid))
    allocate(zsig(ngrid))
    allocate(zgam(ngrid))
    allocate(zthe(ngrid))
    allocate(z0(ngrid))
    allocate(qsurf(ngrid,nq))
    allocate(tsurf(ngrid))
    allocate(co2ice(ngrid))
    allocate(watercap(ngrid))
    allocate(emis(ngrid))
    allocate(capcal(ngrid))
    allocate(fluxgrd(ngrid))
    allocate(hmons(ngrid))
    allocate(summit(ngrid))
    allocate(base(ngrid))
    
  end subroutine ini_surfdat_h


  subroutine end_surfdat_h

  implicit none

    if (allocated(albedodat))   deallocate(albedodat)
    if (allocated(phisfi))      deallocate(phisfi)
    if (allocated(watercaptag)) deallocate(watercaptag)
    if (allocated(dryness))     deallocate(dryness)
    if (allocated(zmea))        deallocate(zmea)
    if (allocated(zstd))        deallocate(zstd)
    if (allocated(zsig))        deallocate(zsig)
    if (allocated(zgam))        deallocate(zgam)
    if (allocated(zthe))        deallocate(zthe)
    if (allocated(z0))          deallocate(z0)
    if (allocated(qsurf))       deallocate(qsurf)
    if (allocated(tsurf))       deallocate(tsurf)
    if (allocated(co2ice))      deallocate(co2ice)
    if (allocated(watercap))    deallocate(watercap)
    if (allocated(emis))        deallocate(emis)
    if (allocated(capcal))      deallocate(capcal)
    if (allocated(fluxgrd))     deallocate(fluxgrd)
    if (allocated(hmons))       deallocate(hmons)
    if (allocated(summit))      deallocate(summit)
    if (allocated(base))        deallocate(base)
    
  end subroutine end_surfdat_h

end module surfdat_h
