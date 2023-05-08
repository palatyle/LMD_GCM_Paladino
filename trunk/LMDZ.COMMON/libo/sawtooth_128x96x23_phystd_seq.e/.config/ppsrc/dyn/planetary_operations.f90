










!
! $Id: $
!
module planetary_operations
! module which contains functions to obtain total values over the
! entire globe (trivial in serial mode, but not in parallel)

implicit none

contains

  subroutine planetary_atm_mass_from_ps(ps,value)
  implicit none
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  real,intent(in) :: ps(ip1jmp1)
  real,intent(out) :: value
  integer :: i,j,ij

  ! compute total atmospheric mass over the whole planet
  ! North Pole
  value=ps(1)*airesurg(1)*iim
  do j=2,jjm
    do i=1,iim
      ij=i+(j-1)*iip1
      value=value+ps(ij)*airesurg(ij)
    enddo
  enddo
  ! South pole
  value=value+ps(ip1jmp1-iim)*airesurg(ip1jmp1-iim)*iim
  
  end subroutine planetary_atm_mass_from_ps

  subroutine planetary_atm_mass_from_mass(mass,value)
  implicit none
  include "dimensions.h"
  include "paramet.h"
  real,intent(in) :: mass(ip1jmp1,llm) ! air mass (kg)
  real,intent(out) :: value
  integer :: i,j,ij,l
  real :: column(ip1jmp1) ! columns amount of air (kg)

  ! compute total atmospheric mass over the whole planet
  ! 1. build column amout of tracer (kg)
  column(:)=0
  do ij=1,ip1jmp1
    do l=1,llm
      column(ij)=column(ij)+mass(ij,l)
    enddo
  enddo
  
  ! 2. Compute total tracer mass over the planet
  !North pole
  value=column(1)*iim
  do j=2,jjm
    do i=1,iim
      ij=i+(j-1)*iip1
      value=value+column(ij)
    enddo
  enddo
  ! South pole
  value=value+column(ip1jmp1-iim)*iim
  
  end subroutine planetary_atm_mass_from_mass
  
  subroutine planetary_tracer_amount_from_p(p,q,amount)
  use comconst_mod, ONLY: g
  implicit none
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  real,intent(in) :: p(ip1jmp1,llmp1) ! pressure at vertical mesh interfaces
  real,intent(in) :: q(ip1jmp1,llm) ! 3D tracer (kg/kg_air)
  real,intent(out) :: amount
  integer :: i,j,ij,l
  real :: column(ip1jmp1) ! columns amount of tracer (kg.m-2)
  
  ! 1. build column amout of tracer (kg.m-2)
  column(:)=0
  do ij=1,ip1jmp1
    do l=1,llm
      column(ij)=column(ij)+q(ij,l)*(p(ij,l)-p(ij,l+1))/g
    enddo
  enddo
  
  ! 2. Compute total tracer mass over the planet
  !North pole
  amount=column(1)*aire(1)*iim
  do j=2,jjm
    do i=1,iim
      ij=i+(j-1)*iip1
      amount=amount+column(ij)*aire(ij)
    enddo
  enddo
  ! South pole
  amount=amount+column(ip1jmp1-iim)*aire(ip1jmp1-iim)*iim

  end subroutine planetary_tracer_amount_from_p

  subroutine planetary_tracer_amount_from_mass(mass,q,amount)
  implicit none
  include "dimensions.h"
  include "paramet.h"
  real,intent(in) :: mass(ip1jmp1,llm) ! air mass (kg)
  real,intent(in) :: q(ip1jmp1,llm) ! 3D tracer (kg/kg_air)
  real,intent(out) :: amount
  integer :: i,j,ij,l
  real :: column(ip1jmp1) ! columns amount of tracer (kg)
  
  ! 1. build column amout of tracer (kg)
  column(:)=0
  do ij=1,ip1jmp1
    do l=1,llm
      column(ij)=column(ij)+q(ij,l)*mass(ij,l)
    enddo
  enddo
  
  ! 2. Compute total tracer mass over the planet
  !North pole
  amount=column(1)*iim
  do j=2,jjm
    do i=1,iim
      ij=i+(j-1)*iip1
      amount=amount+column(ij)
    enddo
  enddo
  ! South pole
  amount=amount+column(ip1jmp1-iim)*iim

  end subroutine planetary_tracer_amount_from_mass

end module planetary_operations
