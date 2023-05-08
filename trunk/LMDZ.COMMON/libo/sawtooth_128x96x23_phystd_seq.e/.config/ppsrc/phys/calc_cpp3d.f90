










      subroutine calc_cpp3d(cppNI,rcpNI,t,p)

!==================================================================
!     Purpose
!     -------
!     Compute the atmospheric specific heat capacity as a 
!     function of pressure and temperature (for CO2 gas only)
!
!     Authors
!     ------- 
!     Robin Wordsworth (2009)
!
!==================================================================
      use comcstfi_mod, only: cpp, r
      implicit none

      !real cp0, dB2dT2
      real cppNI      ! specific heat capacity at const. pressure
      real rcpNI      ! R / specific heat capacity
      real t
      real p

      ! dummy function until I decide what to do!

      cppNI  = cpp
      rcpNI  = R/cppNI

      return
    end subroutine calc_cpp3d
