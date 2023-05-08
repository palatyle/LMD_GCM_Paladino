      module slab_ice_h

      implicit none

! noceanmx : number of oceanic layers
      integer, parameter :: noceanmx = 2 ! number of oceanic layers

      real, parameter :: ice_cap=2.05e+03       ! J/kg/K
      real, parameter :: ice_den=945.0          ! kg/m3
      real, parameter :: ice_lat=0.334e+05      ! J/kg
      real, parameter :: ice_cond=2.1           ! W/m/K
      real, parameter :: alb_ice_min=0.2       
      real, parameter :: alb_ice_max=0.65 
      real, parameter :: alb_ocean=0.07
      real, parameter :: ice_frac_min=0.0001
      real, parameter :: ice_frac_max=1.!0.9999
      real, parameter :: h_alb_ice=0.5*ice_den   !m->kg/m2
      real, parameter :: h_ice_thin=0.2*ice_den
      real, parameter :: h_ice_thick=2.5*ice_den
      real, parameter :: h_ice_min=0.000001*ice_den
      real, parameter :: h_ice_max=10000.0

      real, parameter :: capcalocean=50.*4.228e+06!121635.0
      real, parameter :: capcalseaice=5.1444e+06*0.15
      real, parameter :: capcalsno=2.3867e+06*0.15

      real, parameter :: epsfra=1.0E-05
      real, parameter :: soil_hdiff=25000.0

!      real, parameter :: calice=1.0/(5.1444e+06*0.15)
!      real, parameter :: tau_gl=86400.*5.
!      real, parameter :: calsno=1./(2.3867e+06*0.15)






      end module slab_ice_h
