module dimradmars_mod
!   Declaration and settings for radiative transfer calculations
!   Initializations and allocations are done in phys_state_var_init
implicit none
  ! nflev: number of vertical layer
  ! ndlon,ndlo2: number of horizontal points
  ! Splitting of horizontal grid
  ! NDLO2 and ndomainsz for the splitting in the physics call
  ! WARNING:  One must have  1 < ndomainsz =< ngrid
  integer,save :: NFLEV !=nlayer   ! with splitting
  integer,save :: ndomainsz !=(ngrid-1)/20 + 1
  integer,save :: NDLON !=ndomainsz  ! with splitting
  integer,save :: NDLO2 !=NDLON


! Number of kind of tracer radiative properties
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! naerkind is set by reading callphys.def
! -- see conf_phys
! -- value of nsizemax below is comfortably high
!    but beware in case you add a lot of scatterers
  INTEGER, SAVE :: naerkind

  ! AS: previously in aerkind.h
  character*20, SAVE, ALLOCATABLE :: name_iaer(:)  ! name of the scatterers
  integer iaer_dust_conrath ! Typical dust profiles using a
                            ! Conrath type analytical equation
  integer iaer_dust_doubleq ! Dust profile is given by the
                            ! mass mixing ratio of the two-
                            ! moment scheme method (doubleq)
  integer iaer_dust_submicron ! Dust profile is given by a
                              ! submicron population of dust
                              ! particles
  integer iaer_stormdust_doubleq ! Storm dust profile is given by the
                              ! mass mixing ratio of the two moment scheme 
                              ! method (doubleq)
  integer iaer_topdust_doubleq ! top dust profile is given by the
                              ! mass mixing ratio of the two moment scheme 
                              ! method (doubleq)
  integer iaer_h2o_ice ! Water ice particles

  ! AS: was in aeropacity
  INTEGER,SAVE,ALLOCATABLE :: iaerdust(:)

  ! AS: was in suaer
  CHARACTER(LEN=30), SAVE, ALLOCATABLE :: file_id(:,:)

! Reference wavelengths used to compute reference optical depth (m)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  REAL,SAVE,ALLOCATABLE :: longrefir(:),longrefvis(:)
  
! Definition of spectral intervals at thermal infrared wavelengths (LW)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer,parameter :: nir=4 ! Total number of thermal IR bands
  integer,parameter :: nuco2=2 ! number of bands in CO2 bands
  real,parameter :: long1ir=5.E-6 , long2ir=200.E-6
  real,parameter :: long1co2=1.E+0 / 865.E+2 , long2co2=1.E+0 / 500.E+2

!  Warning : the "nir" thermal IR bands are not ordered by wavelength:
!      iir=1 : central 15um CO2 bands     
!      iir=2 : CO2 band wings    [long1co2-long2co2] MINUS central band
!      iir=3 : 9 um band [long1ir - long1co2]
!      iir=4 : Far IR    [long2co2 - long2ir]
    
!  Definition of spectral interval at solar wavelengths (SW)
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer,parameter :: NSUN=2   ! do not change that !
!  Boundaries of spectral intervals (m) : 
  real,parameter :: long1vis=0.1E-6 , long2vis=0.5E-6 , long3vis=5.E-6
!  Fraction of solar energy in solar band #1 [long1vis-long2vis] : 0.274490
!  Fraction of solar energy in solar band #2 [long2vis-long3vis] : 0.725509
  real,save :: sunfr(2) = (/ 0.274490 , 0.725509 /)

! Maximum number of grain size classes
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This parameter has to be set to the maximum number of particle
!   sizes contained in the optical parameter database; For example,
!   if only one grain size is used to describe dust, and 30 are used
!   to describe water-ice crystals in the visible and 15 in the IR,
!   nsizemax has to be set to 30.
! If only one grain size is considered for all the aerosols, set
!   this parameter to 1 and convolution will be turned off during
!   the radiative calculations.

  integer, parameter :: nsizemax = 60
! integer, parameter :: nsizemax = 1

! Various initialisation for LW radiative code
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! npademx : number of Pade coefficients
! nabsmx : ?
! nt_pademx : number of temperature intervals for Pade

  integer,parameter :: npademx=4
  integer,parameter :: nabsmx=2
  integer,parameter :: nt_pademx=19

!!
!! variables
!!
  REAL,SAVE,ALLOCATABLE :: dtrad(:,:) ! Net atm. radiative heating rate (K.s-1)
  REAL,SAVE,ALLOCATABLE :: fluxrad_sky(:) ! rad. flux from sky absorbed by surface (W.m-2)
  REAL,SAVE,ALLOCATABLE :: fluxrad(:) ! Net radiative surface flux (W.m-2)
  REAL,SAVE,ALLOCATABLE :: albedo(:,:) ! Surface albedo in each solar band
  REAL,SAVE,ALLOCATABLE :: tauscaling(:)   ! Convertion factor for qdust and Ndust
  REAL,SAVE,ALLOCATABLE :: totcloudfrac(:) ! total cloud fraction over the column
! aerosol (dust or ice) extinction optical depth  at reference wavelength 
! "longrefvis" set in dimradmars_mod , for one of the "naerkind"  kind of
! aerosol optical properties  :
  REAL,SAVE,ALLOCATABLE :: aerosol(:,:,:)
  REAL,SAVE,ALLOCATABLE :: nueffdust(:,:) ! Dust effective variance

!! ------------------------------------------------------
!! AS: what was previously in yomaer
!   Shortwave
!   ~~~~~~~~~
! 
! tauvis: dust optical depth at reference wavelength  ("longrefvis" set
! in dimradmars_mod : typically longrefvis = 0.67E-6 m, as measured by Viking )

! For the "naerkind" kind of aerosol radiative properties : 
! QVISsQREF  :  Qext / Qext("longrefvis")   <--- For both solar bands
! omegavis   :  sinle scattering albedo     <--- For both solar bands
! gvis       :  assymetry factor            <--- For both solar bands
! 
!   Longwave
!   ~~~~~~~~
! 
! For the "naerkind" kind of aerosol radiative properties : 
! QIRsQREF :  Qext / Qext("longrefvis")     <--- For the nir bandes IR
! omegaIR  :  mean single scattering albedo <--- For the nir bandes IR
! gIR      :  mean assymetry factor         <--- For the nir bandes IR
! 
  real,save :: tauvis
  real,save,allocatable :: QVISsQREF(:,:,:)
  real,save,allocatable :: omegavis(:,:,:)
  real,save,allocatable :: gvis(:,:,:)
  real,save,allocatable :: QIRsQREF(:,:,:)
  real,save,allocatable :: omegaIR(:,:,:)
  real,save,allocatable :: gIR(:,:,:)
! Actual number of grain size classes in each domain for a
!   given aerosol:
  integer,save,allocatable :: nsize(:,:)
! Particle size axis (depend on the kind of aerosol and the
!   radiation domain)
  real,save,allocatable :: radiustab(:,:,:)
! Extinction coefficient at reference wavelengths;
!   These wavelengths are defined in dimradmars_mod, and called
!   longrefvis and longrefir.
  real,save,allocatable :: QREFvis(:,:)
  real,save,allocatable :: QREFir(:,:)
  real,save,allocatable :: omegaREFvis(:,:)
  real,save,allocatable :: omegaREFir(:,:)
!! ------------------------------------------------------

contains

  subroutine ini_dimradmars_mod(ngrid,nlayer)
  
  implicit none
  
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric layers

   nflev=nlayer
!  ndomainsz=ngrid
   ndomainsz=(ngrid-1)/20 + 1
!  ndomainsz=(ngrid-1)/5 + 1
   ndlon=ndomainsz
   ndlo2=ndlon

   allocate(albedo(ngrid,2))
   allocate(dtrad(ngrid,nlayer))
   allocate(fluxrad_sky(ngrid))
   allocate(fluxrad(ngrid))
   allocate(tauscaling(ngrid))
   allocate(nueffdust(ngrid,nlayer))
   allocate(totcloudfrac(ngrid))

  end subroutine ini_dimradmars_mod

  subroutine end_dimradmars_mod

  implicit none

   if (allocated(albedo))      deallocate(albedo)
   if (allocated(dtrad))       deallocate(dtrad)
   if (allocated(fluxrad_sky)) deallocate(fluxrad_sky)
   if (allocated(fluxrad))     deallocate(fluxrad)
   if (allocated(tauscaling))  deallocate(tauscaling)
   if (allocated(nueffdust))   deallocate(nueffdust)
   if (allocated(totcloudfrac))   deallocate(totcloudfrac)

  end subroutine end_dimradmars_mod

 
  subroutine ini_scatterers(ngrid,nlayer)

  implicit none

  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric layers

   !! domain-dependent 
   !! -- only used in physiq_mod & intent(out) in callradite
   if (allocated(aerosol)) deallocate(aerosol)
   allocate(aerosol(ngrid,nlayer,naerkind))

   !! not domain-dependent
   if (.not.allocated(name_iaer)) allocate(name_iaer(naerkind))
   if (.not.allocated(longrefir)) allocate(longrefir(naerkind))
   if (.not.allocated(longrefvis)) allocate(longrefvis(naerkind))
   if (.not.allocated(iaerdust)) allocate(iaerdust(naerkind))
   if (.not.allocated(file_id)) allocate(file_id(naerkind,2))
   if (.not.allocated(QVISsQREF)) allocate(QVISsQREF(nsun,naerkind,nsizemax))
   if (.not.allocated(omegavis)) allocate(omegavis(nsun,naerkind,nsizemax))
   if (.not.allocated(gvis)) allocate(gvis(nsun,naerkind,nsizemax))
   if (.not.allocated(QIRsQREF)) allocate(QIRsQREF(nir,naerkind,nsizemax))
   if (.not.allocated(omegaIR)) allocate(omegaIR(nir,naerkind,nsizemax))
   if (.not.allocated(gIR)) allocate(gIR(nir,naerkind,nsizemax))
   if (.not.allocated(nsize)) allocate(nsize(naerkind,2))
   if (.not.allocated(radiustab)) allocate(radiustab(naerkind,2,nsizemax))
   if (.not.allocated(QREFvis)) allocate(QREFvis(naerkind,nsizemax))
   if (.not.allocated(QREFir)) allocate(QREFir(naerkind,nsizemax))
   if (.not.allocated(omegaREFvis)) allocate(omegaREFvis(naerkind,nsizemax))
   if (.not.allocated(omegaREFir)) allocate(omegaREFir(naerkind,nsizemax))

  end subroutine ini_scatterers

end module dimradmars_mod
