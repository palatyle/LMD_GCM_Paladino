!
! $Id $
!
module com_io_dyn_mod

  implicit none 

! Names of various files for outputs (in the dynamics)
  ! to store instantaneous values:
  character(len=18),parameter :: dynhist_file="dyn_hist.nc" ! on scalar grid
  character(len=18),parameter :: dynhistv_file="dyn_histv.nc" ! on v grid
  character(len=18),parameter :: dynhistu_file="dyn_histu.nc" ! on u grid

  ! to store averaged values:
  character(len=18),parameter :: dynhistave_file="dyn_hist_ave.nc"
  character(len=18),parameter :: dynhistvave_file="dyn_histv_ave.nc"
  character(len=18),parameter :: dynhistuave_file="dyn_histu_ave.nc"
  
! Ids of various files for outputs (in the dynamics)

  ! instantaneous (these are set by inithist.F)
  integer :: histid
  integer :: histvid
  integer :: histuid
  
  ! averages (these are set by initdynav.F)
  integer :: histaveid
  integer :: histvaveid
  integer :: histuaveid
  
end module com_io_dyn_mod
