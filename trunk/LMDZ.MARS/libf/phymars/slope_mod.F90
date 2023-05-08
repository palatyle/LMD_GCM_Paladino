module slope_mod

implicit none

  real,save,allocatable :: theta_sl(:) ! slope angle versus horizontal (deg)
  real,save,allocatable :: psi_sl(:)   ! slope orientation (deg)

contains

  subroutine ini_slope_mod(ngrid)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  
  allocate(theta_sl(ngrid))
  allocate(psi_sl(ngrid))
  
  end subroutine ini_slope_mod


  subroutine end_slope_mod

  implicit none

  if (allocated(theta_sl)) deallocate(theta_sl)
  if (allocated(psi_sl)) deallocate(psi_sl)

  end subroutine end_slope_mod
  
end module slope_mod
