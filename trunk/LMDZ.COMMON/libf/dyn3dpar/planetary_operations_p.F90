!
! $Id: $
!
module planetary_operations_p
! module which contains functions to obtain total values over the
! entire globe (in parallel)

implicit none

contains

  subroutine planetary_atm_mass_from_ps_p(ps,value)
  USE parallel_lmdz, ONLY: Gather_Field, mpi_rank
  USE mod_const_mpi, ONLY: COMM_LMDZ
  implicit none
  include "dimensions.h"
  include "paramet.h"
  include "comgeom.h"
  include 'mpif.h'
  real,intent(in) :: ps(ip1jmp1)
  real,intent(out) :: value
  integer :: i,j,ij
  integer :: ierr

  ! compute total atmospheric mass over the whole planet
  ! 1. Gather ps on master to do the computations
  call Gather_Field(ps,ip1jmp1,1,0)

  ! 2. compute on master
  if (mpi_rank==0) then
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
  endif ! of if (mpi_rank==0)

  ! broadcast the result to all procs
!$OMP CRITICAL (MPI_BCAST)
#ifdef CPP_MPI
  call MPI_BCAST(value,1,MPI_REAL8,0,COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI_BCAST)
  
  end subroutine planetary_atm_mass_from_ps_p

  subroutine planetary_tracer_amount_from_mass_p(mass,q,amount)
  USE parallel_lmdz, ONLY: ij_begin,ij_end, &
                           pole_nord,pole_sud, &
                           Gather_Field, mpi_rank
  USE mod_const_mpi, ONLY: COMM_LMDZ
  implicit none
  include "dimensions.h"
  include "paramet.h"
  include 'mpif.h'
  real,intent(in) :: mass(ip1jmp1,llm) ! air mass (kg)
  real,intent(in) :: q(ip1jmp1,llm) ! 3D tracer (kg/kg_air)
  real,intent(out) :: amount
  integer :: i,j,ij,l
  integer :: ijb,ije
  integer :: ierr
  real :: column(ip1jmp1) ! columns amount of tracer (kg)
  
  ! 1. build column amout of tracer (kg) on each proc
  ijb=ij_begin-iip1
  ije=ij_end+2*iip1
  if (pole_nord) ijb=ij_begin
  if (pole_sud)  ije=ij_end
  
  column(ijb:ije)=0
  do ij=ijb,ije
    do l=1,llm
      column(ij)=column(ij)+q(ij,l)*mass(ij,l)
    enddo
  enddo
  
  ! 2 Gather "column" to do the upcoming computations on master
  call Gather_Field(column,ip1jmp1,1,0)
  
  ! 3. Compute total tracer mass over the planet on master
  if (mpi_rank==0) then
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
  endif ! of if (mpi_rank==0)
  
  ! broadcast the result to all procs
!$OMP CRITICAL (MPI_BCAST)
#ifdef CPP_MPI
  call MPI_BCAST(amount,1,MPI_REAL8,0,COMM_LMDZ,ierr)
#endif
!$OMP END CRITICAL (MPI_BCAST)

  end subroutine planetary_tracer_amount_from_mass_p

end module planetary_operations_p
