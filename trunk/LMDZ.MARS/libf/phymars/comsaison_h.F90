module comsaison_h

implicit none

  logical,save :: callsais
  integer,save :: isaison
  real,save :: dist_sol
  real,save :: declin
  real,save,allocatable :: mu0(:)
  real,save,allocatable :: fract(:)
  real,save,allocatable :: local_time(:) ! local solar time as fraction of day (0,1)

contains

  subroutine ini_comsaison_h(ngrid)
  
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  
    allocate(mu0(ngrid))
    allocate(fract(ngrid))
    allocate(local_time(ngrid))
  
  end subroutine ini_comsaison_h


  subroutine end_comsaison_h

  implicit none

    if (allocated(mu0)) deallocate(mu0)
    if (allocated(fract)) deallocate(fract)
    if (allocated(local_time)) deallocate(local_time)

  end subroutine end_comsaison_h
  
end module comsaison_h
