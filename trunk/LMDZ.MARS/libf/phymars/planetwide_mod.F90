!
! $Id: $
!
module planetwide_mod
! module which contains functions to obtain max/min/... values over the
! entire globe (trivial in serial mode, but not in parallel)

use mod_phys_lmdz_para, only : is_master, gather, bcast
                               
implicit none

interface planetwide_maxval ! maxval() , over the entire planet
  module procedure planetwide_maxval_i1, planetwide_maxval_i2, &
                   planetwide_maxval_r1, planetwide_maxval_r2
end interface

interface planetwide_minval ! minval() , over the entire planet
  module procedure planetwide_minval_i1, planetwide_minval_i2, &
                   planetwide_minval_r1, planetwide_minval_r2
end interface

interface planetwide_sumval ! sum() , over the entire planet
  module procedure planetwide_sumval_i1, planetwide_sumval_i2, &
                   planetwide_sumval_r1, planetwide_sumval_r2
end interface

contains

  subroutine planetwide_maxval_i1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:) ! local grid (klon)
  integer,intent(out) :: values_max
#ifdef CPP_PARA
  integer :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=maxval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=maxval(values)
#endif
  end subroutine planetwide_maxval_i1

  subroutine planetwide_maxval_i2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_max
#ifdef CPP_PARA
  integer :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=maxval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=maxval(values)
#endif
  end subroutine planetwide_maxval_i2

  subroutine planetwide_maxval_r1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_max
#ifdef CPP_PARA
  real :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=maxval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=maxval(values)
#endif
  end subroutine planetwide_maxval_r1

  subroutine planetwide_maxval_r2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_max
#ifdef CPP_PARA
  real :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=maxval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=maxval(values)
#endif
  end subroutine planetwide_maxval_r2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine planetwide_minval_i1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:) ! local grid (klon)
  integer,intent(out) :: values_max
#ifdef CPP_PARA
  integer :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=minval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=minval(values)
#endif
  end subroutine planetwide_minval_i1

  subroutine planetwide_minval_i2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_max
#ifdef CPP_PARA
  integer :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=minval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=minval(values)
#endif
  end subroutine planetwide_minval_i2

  subroutine planetwide_minval_r1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_max
#ifdef CPP_PARA
  real :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=minval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=minval(values)
#endif
  end subroutine planetwide_minval_r1

  subroutine planetwide_minval_r2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_max
#ifdef CPP_PARA
  real :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! extract maximum value
  if (is_master) then
    values_max=minval(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_max)
#else
  values_max=minval(values)
#endif
  end subroutine planetwide_minval_r2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine planetwide_sumval_i1(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:) ! local grid (klon)
  integer,intent(out) :: values_sum
#ifdef CPP_PARA
  integer :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! calculate sum value
  if (is_master) then
    values_sum=SUM(values_glo(:))
  endif
  ! broadcast information to all cores
  call bcast(values_sum)
#else
  values_sum=SUM(values(:))
#endif
  end subroutine planetwide_sumval_i1

  subroutine planetwide_sumval_i2(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_sum
#ifdef CPP_PARA
  integer :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! calculate sum value
  if (is_master) then
    values_sum=SUM(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_sum)
#else
  values_sum=SUM(values)
#endif
  end subroutine planetwide_sumval_i2

  subroutine planetwide_sumval_r1(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_sum
#ifdef CPP_PARA
  real :: values_glo(klon_glo) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! calculate sum value
  if (is_master) then
    values_sum=SUM(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_sum)
#else
  values_sum=SUM(values)
#endif
  end subroutine planetwide_sumval_r1

  subroutine planetwide_sumval_r2(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_sum
#ifdef CPP_PARA
  real :: values_glo(klon_glo,size(values,2)) ! global grid
  
  ! gather field on master:
  call gather(values,values_glo)
  ! calculate sum value
  if (is_master) then
    values_sum=SUM(values_glo)
  endif
  ! broadcast information to all cores
  call bcast(values_sum)
#else
  values_sum=SUM(values)
#endif
  end subroutine planetwide_sumval_r2
  
  
end module planetwide_mod
