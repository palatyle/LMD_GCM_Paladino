










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
  REAL,SAVE,ALLOCATABLE :: comm_CLOUDFRAC(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_TOTCLOUDFRAC(:)
  REAL,SAVE,ALLOCATABLE :: comm_RAIN(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_SNOW(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_ALBEQ(:) 
  REAL,SAVE,ALLOCATABLE :: comm_FLUXTOP_DN(:)
  REAL,SAVE,ALLOCATABLE :: comm_FLUXABS_SW(:)
  REAL,SAVE,ALLOCATABLE :: comm_FLUXTOP_LW(:)
  REAL,SAVE,ALLOCATABLE :: comm_FLUXSURF_SW(:)
  REAL,SAVE,ALLOCATABLE :: comm_FLUXSURF_LW(:)
  REAL,SAVE,ALLOCATABLE :: comm_FLXGRD(:)
  REAL,SAVE,ALLOCATABLE :: comm_LSCEZ(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_H2OICE_REFF(:,:)
  REAL,SAVE,ALLOCATABLE :: comm_LATENT_HF(:)
  REAL,SAVE,ALLOCATABLE :: comm_REEVAP(:)
  REAL,SAVE,ALLOCATABLE :: comm_SURFRAIN(:)

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
  allocate(comm_CLOUDFRAC(ngrid,nlayer))
  allocate(comm_TOTCLOUDFRAC(ngrid))
  allocate(comm_RAIN(ngrid,nlayer))
  allocate(comm_SNOW(ngrid,nlayer))
  allocate(comm_ALBEQ(ngrid))
  allocate(comm_FLUXTOP_DN(ngrid))
  allocate(comm_FLUXABS_SW(ngrid))
  allocate(comm_FLUXTOP_LW(ngrid))
  allocate(comm_FLUXSURF_SW(ngrid))
  allocate(comm_FLUXSURF_LW(ngrid))
  allocate(comm_FLXGRD(ngrid))
  allocate(comm_LSCEZ(ngrid,nlayer))
  allocate(comm_H2OICE_REFF(ngrid,nlayer))
  allocate(comm_LATENT_HF(ngrid))
  allocate(comm_REEVAP(ngrid))
  allocate(comm_SURFRAIN(ngrid))

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
  deallocate(comm_CLOUDFRAC)
  deallocate(comm_TOTCLOUDFRAC)
  deallocate(comm_RAIN)
  deallocate(comm_SNOW)
  deallocate(comm_ALBEQ)
  deallocate(comm_FLUXTOP_DN)
  deallocate(comm_FLUXABS_SW)
  deallocate(comm_FLUXTOP_LW)
  deallocate(comm_FLUXSURF_SW)
  deallocate(comm_FLUXSURF_LW)
  deallocate(comm_FLXGRD)
  deallocate(comm_LSCEZ)
  deallocate(comm_H2OICE_REFF)
  deallocate(comm_LATENT_HF)
  deallocate(comm_REEVAP)
  deallocate(comm_SURFRAIN)

  end subroutine deallocate_comm_wrf

end module comm_wrf

