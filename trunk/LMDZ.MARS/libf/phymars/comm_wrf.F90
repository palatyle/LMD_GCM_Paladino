module comm_wrf

!!
!! This module is useful to output fields in WRF style
!! -- useful for diagnostics that are not already in modules
!!

  REAL,SAVE,ALLOCATABLE :: comm_HR_SW(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_HR_LW(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_SWDOWNZ(:)
  REAL,SAVE,ALLOCATABLE :: comm_TAU_DUST(:)
  REAL,SAVE,ALLOCATABLE :: comm_RDUST(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_QSURFDUST(:)
  REAL,SAVE,ALLOCATABLE :: comm_MTOT(:)
  REAL,SAVE,ALLOCATABLE :: comm_ICETOT(:)
  REAL,SAVE,ALLOCATABLE :: comm_VMR_ICE(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_TAU_ICE(:)
  REAL,SAVE,ALLOCATABLE :: comm_RICE(:,:)

contains

  subroutine allocate_comm_wrf(ngrid,nlayer)
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric layers
  if (.not.allocated(comm_HR_SW)) allocate(comm_HR_SW(ngrid,nlayer))
  if (.not.allocated(comm_HR_LW)) allocate(comm_HR_LW(ngrid,nlayer))
  if (.not.allocated(comm_SWDOWNZ)) allocate(comm_SWDOWNZ(ngrid))
  if (.not.allocated(comm_TAU_DUST)) allocate(comm_TAU_DUST(ngrid))
  if (.not.allocated(comm_RDUST)) allocate(comm_RDUST(ngrid,nlayer))
  if (.not.allocated(comm_QSURFDUST)) allocate(comm_QSURFDUST(ngrid))
  if (.not.allocated(comm_MTOT)) allocate(comm_MTOT(ngrid))
  if (.not.allocated(comm_ICETOT)) allocate(comm_ICETOT(ngrid))
  if (.not.allocated(comm_VMR_ICE)) allocate(comm_VMR_ICE(ngrid,nlayer))
  if (.not.allocated(comm_TAU_ICE)) allocate(comm_TAU_ICE(ngrid))
  if (.not.allocated(comm_RICE)) allocate(comm_RICE(ngrid,nlayer))
  end subroutine allocate_comm_wrf

  subroutine deallocate_comm_wrf
  implicit none
  deallocate(comm_HR_SW)
  deallocate(comm_HR_LW)
  deallocate(comm_SWDOWNZ)
  deallocate(comm_TAU_DUST)
  deallocate(comm_RDUST)
  deallocate(comm_QSURFDUST)
  deallocate(comm_MTOT)
  deallocate(comm_ICETOT)
  deallocate(comm_VMR_ICE)
  deallocate(comm_TAU_ICE)
  deallocate(comm_RICE)
  end subroutine deallocate_comm_wrf

end module comm_wrf

