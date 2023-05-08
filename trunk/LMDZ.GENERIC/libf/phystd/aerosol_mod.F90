!==================================================================
module aerosol_mod
implicit none
save
!==================================================================

!  aerosol indexes: these are initialized to be 0 if the
!                 corresponding aerosol was not activated in callphys.def
!                 -- otherwise a value is given in iniaerosol
      integer :: iaero_co2 = 0 
      integer :: iaero_h2o = 0
      integer :: iaero_dust = 0
      integer :: iaero_h2so4 = 0
      logical :: noaero = .false.

! two-layer simple aerosol model
      integer :: iaero_back2lay = 0
 ! NH3 cloud
      integer :: iaero_nh3 = 0
 ! Auroral aerosols
      integer :: iaero_aurora = 0
!$OMP THREADPRIVATE(iaero_co2,iaero_h2o,iaero_dust,iaero_h2so4,noaero,iaero_back2lay,iaero_nh3,iaero_aurora)
      
!==================================================================
end module aerosol_mod
!==================================================================
