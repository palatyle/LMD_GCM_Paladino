module comsoil_h

implicit none
! nsoilmx : number of subterranean layers
!integer, parameter :: nsoilmx = 18 ! for z1=0.0002 m, depth = 18 m => mars case 
!integer, parameter :: nsoilmx = 13 ! for z1=0.03 m, depth = 104.8 m => earth case
  integer,save :: nsoilmx = 18 ! default, but may be set in callphys.def
! Full soil layer depths are set as: layer(k)=lay1_soil*alpha_soil**(k-1) , k=1,nsoil
! Mid soil layer depths are set as: mlayer(k)=lay1_soil*alpha_soil**(k-1/2) , k=0,nsoil-1 
  real,save :: lay1_soil=2.e-4 ! depth (m) of first full soil layer (may be set in callphys.def)
  real,save :: alpha_soil=2 ! coefficient for soil layer thickness (may be set in callphys.def)
!$OMP THREADPRIVATE(nsoilmx,lay1_soil,alpha_soil)

  real,save,allocatable,dimension(:) :: layer      ! soil layer depths
  real,save,allocatable,dimension(:) :: mlayer     ! soil mid-layer depths
  real,save,allocatable,dimension(:,:) :: inertiedat ! soil thermal inertia
  real,save :: volcapa    ! soil volumetric heat capacity
       ! NB: volcapa is read fromn control(35) from physicq start file
       !     in physdem (or set via tabfi, or initialized in
       !                 soil_settings.F)
!$OMP THREADPRIVATE(layer,mlayer,inertiedat,volcapa)

contains

  subroutine ini_comsoil_h(ngrid)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  
    if (.not.allocated(layer)) allocate(layer(nsoilmx)) !soil layer depths
    if (.not.allocated(mlayer)) allocate(mlayer(0:nsoilmx-1)) ! soil mid-layer depths
    if (.not.allocated(inertiedat)) allocate(inertiedat(ngrid,nsoilmx)) ! soil thermal inertia
  
  end subroutine ini_comsoil_h

end module comsoil_h

