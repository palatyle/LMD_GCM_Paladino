      subroutine lwtt (kdlon,u,up,nu,tr)

c----------------------------------------------------------------------
c     LWTT   computes the longwave transmission functions
c            for all the absorbers in all spectral intervals
c            using pade approximants and horner's algorithm
c----------------------------------------------------------------------

      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "yomlw.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer kdlon            ! part of ngrid
      integer nu               ! 

      real    u (ndlo2,nu)     ! absorber amounts 
      real    up (ndlo2,nu)    ! idem scaled by the pressure

c                                                            outputs:
c                                                            --------
      real    tr (ndlo2,nu)    ! transmission functions

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer ja,jl

      real xn (ndlon)
      real xd (ndlon)
      real ueq (ndlon)

c----------------------------------------------------------------------
c   Transmission by the CO2 15 microns band:
c   ----------------------------------------

      do  ja=1,nu
	    do jl=1,kdlon
c                              equivalent absorber amount (Doppler effect)
c                             --------------------------------------------
	      ueq(jl) = sqrt(up(jl,ja))
     .              +cst_voigt(1,ja)*u(jl,ja)**cst_voigt(2,ja)

c                                                      Horner's algorithm
c                                                      ------------------
            xn(jl) = ga(1,ja) + 
     .      ueq(jl)*(ga(2,ja) + ueq(jl) * ga(3,ja) )
            xd(jl) = gb(1,ja) + ueq(jl)*(gb(2,ja) + 
     .      ueq(jl) * ( gb(3,ja) + ueq(jl)  ))
            tr(jl,ja) = xn(jl) / xd(jl)

	    enddo
      enddo

c----------------------------------------------------------------------
      return
      end
