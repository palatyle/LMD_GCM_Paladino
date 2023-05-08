!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE CONSTANTS_MODULE
!
! This module defines constants that are used by other modules 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module constants_module

   real, parameter :: PI = 3.141592653589793
!   real, parameter :: OMEGA_E = 7.292e-5 ! Angular rotation rate of the earth
!*** Mars
   real, parameter :: OMEGA_E = 7.0721E-5 ! Angular rotation rate of Mars
!*** Mars

   real, parameter :: DEG_PER_RAD = 180./PI
   real, parameter :: RAD_PER_DEG = PI/180.
 
   ! Mean Earth Radius in m.  The value below is consistent
   ! with NCEP's routines and grids.
!************************************************************************
! Earth ellipsoid: only used on specific polar stereographic projections
! - equatorial axis
   real, parameter :: A_WGS84 = 6378137.
! - polar axis
   real, parameter :: B_WGS84 = 6356752.314
   real, parameter :: RE_WGS84 = A_WGS84
   real, parameter :: E_WGS84 = 0.081819192
!************************************************************************
!   real, parameter :: EARTH_RADIUS_M = 6370000.   ! same as MM5 system
!*** Mars
!!necessary if WRFSI intermediate file is used
!!instead of WPS intermediate file format
   real, parameter :: EARTH_RADIUS_M = 3397200.   ! same as MM5 system 
!*** Mars
   real, parameter :: EARTH_CIRC_M = 2.*PI*EARTH_RADIUS_M

end module constants_module