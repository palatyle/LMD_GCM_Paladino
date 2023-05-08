










module conc_mod

implicit none

  real,save,allocatable :: mmean(:,:)  ! mean molecular mass of the atmosphere
  real,save,allocatable :: Akknew(:,:) ! thermal conductivity cofficient
  real,save,allocatable :: cpnew(:,:)  ! specicic heat
  real,save,allocatable :: rnew(:,:)   ! specific gas constant
  
contains

  subroutine ini_conc_mod(ngrid,nlayer)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric levels
  
    allocate(mmean(ngrid,nlayer))
    allocate(Akknew(ngrid,nlayer))
    allocate(cpnew(ngrid,nlayer))
    allocate(rnew(ngrid,nlayer))
    
  end subroutine ini_conc_mod

  subroutine end_conc_mod

  implicit none

    if (allocated(mmean)) deallocate(mmean)
    if (allocated(Akknew)) deallocate(Akknew)
    if (allocated(cpnew)) deallocate(cpnew)
    if (allocated(rnew)) deallocate(rnew)

  end subroutine end_conc_mod


end module conc_mod
