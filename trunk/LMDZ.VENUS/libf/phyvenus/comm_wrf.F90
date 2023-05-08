module comm_wrf

!!
!! This module is useful to output fields in WRF style
!!
!USE only for the mescoscale (comm_HR_SW, comm_HR_LW). OUTPUT
  REAL,SAVE,ALLOCATABLE :: comm_HR_SW(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_HR_LW(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_HR_DYN(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_SWDOWNZ(:)
  REAL,SAVE,ALLOCATABLE :: comm_TAU_DUST(:)
  REAL,SAVE,ALLOCATABLE :: comm_RDUST(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_QSURFDUST(:)
  REAL,SAVE,ALLOCATABLE :: comm_MTOT(:)
  REAL,SAVE,ALLOCATABLE :: comm_ICETOT(:)
  REAL,SAVE,ALLOCATABLE :: comm_VMR_ICE(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_TAU_ICE(:)
  REAL,SAVE,ALLOCATABLE :: comm_RICE(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_GEOP(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_DT(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_DT_RAD(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_DT_VDF(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_DT_AJS(:,:)


contains

  subroutine allocate_comm_wrf(ngrid,nlayer)
  implicit none
  integer,intent(in) :: ngrid ! number of atmospheric columns
  integer,intent(in) :: nlayer ! number of atmospheric layers
  allocate(comm_HR_SW(ngrid,nlayer))
  allocate(comm_HR_LW(ngrid,nlayer))
  allocate(comm_HR_DYN(ngrid,nlayer))
  allocate(comm_SWDOWNZ(ngrid))
  allocate(comm_TAU_DUST(ngrid))
  allocate(comm_RDUST(ngrid,nlayer))
  allocate(comm_QSURFDUST(ngrid))
  allocate(comm_MTOT(ngrid))
  allocate(comm_ICETOT(ngrid))
  allocate(comm_VMR_ICE(ngrid,nlayer))
  allocate(comm_TAU_ICE(ngrid))
  allocate(comm_RICE(ngrid,nlayer))
  allocate(comm_GEOP(ngrid,nlayer))
  allocate(comm_DT(ngrid,nlayer))
  allocate(comm_DT_RAD(ngrid,nlayer))
  allocate(comm_DT_VDF(ngrid,nlayer))
  allocate(comm_DT_AJS(ngrid,nlayer))
  end subroutine allocate_comm_wrf

  subroutine deallocate_comm_wrf
  implicit none
  deallocate(comm_HR_SW)
  deallocate(comm_HR_LW)
  deallocate(comm_HR_DYN)
  deallocate(comm_SWDOWNZ)
  deallocate(comm_TAU_DUST)
  deallocate(comm_RDUST)
  deallocate(comm_QSURFDUST)
  deallocate(comm_MTOT)
  deallocate(comm_ICETOT)
  deallocate(comm_VMR_ICE)
  deallocate(comm_TAU_ICE)
  deallocate(comm_RICE)
  deallocate(comm_GEOP)
  deallocate(comm_DT)
  deallocate(comm_DT_RAD)
  deallocate(comm_DT_VDF)
  deallocate(comm_DT_AJS)
  end subroutine deallocate_comm_wrf

end module comm_wrf

