      subroutine stokes(p,t,rd,w,rho_aer)

!==================================================================
!     Purpose
!     -------
!     Compute the sedimentation velocity in a low pressure 
!     atmosphere.
!
!     Authors
!     ------- 
!     Francois Forget (1997)
!
!==================================================================

!      use radcommon_h, only : Rgas
      use comcstfi_mod, only: g, pi, avocado

      implicit none

!     input
!     -----
!     pressure (Pa), Temperature (K), particle radius (m), density
      real p, t, rd, rho_aer 

!     output
!     ------
!     sedimentation velocity (m/s, >0)
      real w

!     locally saved variables
!     ---------------------
      real a,b,molrad,visc
      save a,b
!$OMP THREADPRIVATE(a,b)
  
      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall/.true./
!$OMP THREADPRIVATE(firstcall)

      if (firstcall) then

         !print*,'Routine not working: replace Rgas with r'
	 !stop
         !a = 0.707*Rgas/(4*pi*molrad**2 * avocado)

         molrad=2.2e-10   ! CO2 (only used in condense_co2cloud at the moment)
         visc=1.e-5       ! CO2


         a = 0.707*8.31/(4*pi*molrad**2 * avocado)
         b = (2./9.) * rho_aer * g / visc 
	!print*,'molrad=',molrad
	!print*,'visc=',visc
	!print*,'a=',a
	!print*,'b=',b
	!print*,'rho_aer=',rho_aer
	!stop
 
         firstcall=.false.
      end if

!     Sedimentation velocity =
!     Stokes' Law corrected for low pressures by the Cunningham
!     slip-flow correction according to Rossow (Icarus 36, 1-50, 1978)
      w = b * rd*rd * (1 + 1.333* (a*T/P)/rd )  
      return
    end subroutine stokes
