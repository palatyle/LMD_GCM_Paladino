      module gases_h

      implicit none

!============================================================================================C
!
!     gases_h    
!
!     * A F90-allocatable version for gases.h -- AS 12/2011
!
!     * Titan's version : J. Vatant d'Ollone (2017)
!                  * gfrac is now 2D (altitude-dependent for the CIA-calculations) 
!                  * Reference level (nivref) is added for the cpp_mu and rayleigh routines
!===========================================================================================C

      ! Set and allocated in su_gases.F90
      integer :: ngasmx
      integer :: nivref
      character*20,allocatable,DIMENSION(:) :: gnom ! name of the gas, set by master
      real,allocatable,DIMENSION(:,:) :: gfrac

      ! in analogy with tracer.h ...
      
      integer :: igas_N2
      integer :: igas_CH4
      integer :: igas_H2
      
!!$OMP THREADPRIVATE(ngasmx,nivref,gnom,gfrac,&
!	!$OMP igas_H2,igas_N2,igas_CH4)

      end module gases_h
