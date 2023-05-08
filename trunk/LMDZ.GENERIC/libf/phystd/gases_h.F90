      module gases_h

      implicit none

!======================================================================C
!
!     gases_h    
!
!     A F90-allocatable version for gases.h -- AS 12/2011
!
!======================================================================C

      ! Set and allocated in su_gases.F90
      integer :: ngasmx
      integer :: vgas
      character*20,allocatable,DIMENSION(:) :: gnom ! name of the gas, read by master
      real,allocatable,DIMENSION(:) :: gfrac

      ! in analogy with tracer.h ...
      integer :: igas_H2
      integer :: igas_He
      integer :: igas_H2O
      integer :: igas_CO2
      integer :: igas_CO
      integer :: igas_N2
      integer :: igas_O2
      integer :: igas_SO2
      integer :: igas_H2S
      integer :: igas_CH4
      integer :: igas_NH3
      integer :: igas_C2H2
      integer :: igas_C2H6
!!$OMP THREADPRIVATE(ngasmx,vgas,gnom,gfrac,&
!	!$OMP igas_H2,igas_He,igas_H2O,igas_CO2,igas_CO,igas_N2,&
!	!$OMP igas_O2,igas_SO2,igas_H2S,igas_CH4,igas_NH3,igas_C2H2,igas_C2H6)

      end module gases_h
