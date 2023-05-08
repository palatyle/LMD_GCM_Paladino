










module turb_mod

  !! variables
  REAL,SAVE,ALLOCATABLE :: q2(:,:)    ! Turbulent Kinetic Energy
  REAL,allocatable,SAVE :: l0(:)
!$OMP THREADPRIVATE(q2,l0)
  REAL,SAVE,ALLOCATABLE :: ustar(:)
  REAL,SAVE,ALLOCATABLE :: wstar(:)
  REAL,SAVE,ALLOCATABLE :: tstar(:)
!$OMP THREADPRIVATE(ustar,wstar,tstar)
  REAL,SAVE,ALLOCATABLE :: hfmax_th(:)
  REAL,SAVE,ALLOCATABLE :: zmax_th(:)
!$OMP THREADPRIVATE(hfmax_th,zmax_th)
  REAL,SAVE,ALLOCATABLE :: sensibFlux(:)
  LOGICAL,SAVE :: turb_resolved = .false.
!$OMP THREADPRIVATE(sensibFlux,turb_resolved)
      ! this is a flag to say 'turbulence is resolved'
      ! mostly for LES use. default is FALSE (for GCM and mesoscale)

!contains

!  subroutine ini_turb_mod(ngrid,nlayer)

!  implicit none
!  integer,intent(in) :: ngrid ! number of atmospheric columns
!  integer,intent(in) :: nlayer ! number of atmospheric layers

!    allocate(q2(ngrid,nlayer+1))
!    allocate(l0(ngrid))
!    allocate(wstar(ngrid))
!   allocate(ustar(ngrid))
!   allocate(tstar(ngrid))
!   allocate(hfmax_th(ngrid))
!   allocate(zmax_th(ngrid))
!   allocate(sensibFlux(ngrid))

!  end subroutine ini_turb_mod

end module turb_mod
