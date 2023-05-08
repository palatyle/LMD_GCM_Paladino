










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
  values_max=maxval(values)
  end subroutine planetwide_maxval_i1

  subroutine planetwide_maxval_i2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_max
  values_max=maxval(values)
  end subroutine planetwide_maxval_i2

  subroutine planetwide_maxval_r1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_max
  values_max=maxval(values)
  end subroutine planetwide_maxval_r1

  subroutine planetwide_maxval_r2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_max
  values_max=maxval(values)
  end subroutine planetwide_maxval_r2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine planetwide_minval_i1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:) ! local grid (klon)
  integer,intent(out) :: values_max
  values_max=minval(values)
  end subroutine planetwide_minval_i1

  subroutine planetwide_minval_i2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_max
  values_max=minval(values)
  end subroutine planetwide_minval_i2

  subroutine planetwide_minval_r1(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_max
  values_max=minval(values)
  end subroutine planetwide_minval_r1

  subroutine planetwide_minval_r2(values,values_max)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_max
  values_max=minval(values)
  end subroutine planetwide_minval_r2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine planetwide_sumval_i1(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:) ! local grid (klon)
  integer,intent(out) :: values_sum
  values_sum=SUM(values(:))
  end subroutine planetwide_sumval_i1

  subroutine planetwide_sumval_i2(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  integer,intent(in) :: values(:,:) ! local grid (klon,...)
  integer,intent(out) :: values_sum
  values_sum=SUM(values)
  end subroutine planetwide_sumval_i2

  subroutine planetwide_sumval_r1(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:) ! local grid (klon)
  real,intent(out) :: values_sum
  values_sum=SUM(values)
  end subroutine planetwide_sumval_r1

  subroutine planetwide_sumval_r2(values,values_sum)
  use dimphy, only: klon
  use mod_grid_phy_lmdz, only : klon_glo
  implicit none
  real,intent(in) :: values(:,:) ! local grid (klon,...)
  real,intent(out) :: values_sum
  values_sum=SUM(values)
  end subroutine planetwide_sumval_r2
  
  
end module planetwide_mod
