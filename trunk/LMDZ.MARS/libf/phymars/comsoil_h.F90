module comsoil_h

implicit none
! nsoilmx : number of subterranean layers
!EM: old soil routine:      integer, parameter :: nsoilmx = 10
  integer, parameter :: nsoilmx = 18 

  real,save,allocatable,dimension(:) :: layer      ! soil layer depths
  real,save,allocatable,dimension(:) :: mlayer     ! soil mid-layer depths
  real,save,allocatable,dimension(:,:) :: inertiedat ! soil thermal inertia
  real,save :: volcapa    ! soil volumetric heat capacity
       ! NB: volcapa is read fromn control(35) from physicq start file
       !     in physdem (or set via tabfi, or initialized in
       !                 soil_settings.F)

  ! variables (FC: built in firstcall in soil.F)
  REAL,SAVE,ALLOCATABLE :: tsoil(:,:)       ! sub-surface temperatures (K)
  real,save,allocatable :: mthermdiff(:,:)  ! (FC) mid-layer thermal diffusivity
  real,save,allocatable :: thermdiff(:,:)   ! (FC) inter-layer thermal diffusivity
  real,save,allocatable :: coefq(:)         ! (FC) q_{k+1/2} coefficients
  real,save,allocatable :: coefd(:,:)       ! (FC) d_k coefficients
  real,save,allocatable :: alph(:,:)        ! (FC) alpha_k coefficients
  real,save,allocatable :: beta(:,:)        ! beta_k coefficients
  real,save :: mu

contains

  subroutine ini_comsoil_h(ngrid)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  
    allocate(layer(nsoilmx)) !soil layer depths
    allocate(mlayer(0:nsoilmx-1)) ! soil mid-layer depths
    allocate(inertiedat(ngrid,nsoilmx)) ! soil thermal inertia
 
    allocate(tsoil(ngrid,nsoilmx)) ! soil temperatures

    allocate(mthermdiff(ngrid,0:nsoilmx-1))
    allocate(thermdiff(ngrid,nsoilmx-1))
    allocate(coefq(0:nsoilmx-1))
    allocate(coefd(ngrid,nsoilmx-1))
    allocate(alph(ngrid,nsoilmx-1))
    allocate(beta(ngrid,nsoilmx-1))

 
  end subroutine ini_comsoil_h


  subroutine end_comsoil_h

  implicit none

    if (allocated(layer)) deallocate(layer)
    if (allocated(mlayer)) deallocate(mlayer)
    if (allocated(inertiedat)) deallocate(inertiedat)
    if (allocated(tsoil)) deallocate(tsoil)
    if (allocated(mthermdiff)) deallocate(mthermdiff)
    if (allocated(thermdiff)) deallocate(thermdiff)
    if (allocated(coefq)) deallocate(coefq) 
    if (allocated(coefd)) deallocate(coefd)
    if (allocated(alph)) deallocate(alph)
    if (allocated(beta)) deallocate(beta)

  end subroutine end_comsoil_h

end module comsoil_h
