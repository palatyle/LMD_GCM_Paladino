      subroutine lwb (kdlon,kflev,tlev,tlay,dt0
     .               ,bsurf,btop,blay,blev,dblay,dbsublay)
 
c----------------------------------------------------------------------
c     LWB  computes the planck function and gradient
c          from a polynomial development of planck function
c----------------------------------------------------------------------

      implicit none
 
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "callkeys.h"

#include "yomlw.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c                                                            inputs:
c                                                            -------
      integer kdlon            ! part of ngrid
      integer kflev            ! part of nalyer

      real dt0 (ndlo2)                 ! surface temperature discontinuity
      real tlay (ndlo2,kflev)          ! layer temperature
      real tlev (ndlo2,kflev+1)        ! level temperature

c                                                            outputs:
c                                                            --------
      real bsurf (ndlo2,nir)        ! surface spectral planck function
      real btop (ndlo2,nir)         ! top spectral planck function
      real blev (ndlo2,nir,kflev+1) ! level   spectral planck function
      real blay (ndlo2,nir,kflev)   ! layer   spectral planck function
      real dblay (ndlo2,nir,kflev)  ! layer gradient spectral planck function
      real dbsublay (ndlo2,nir,2*kflev)  ! layer gradient spectral planck
                                         ! function in sub layers

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer jk, jl, jnu, jk1, jk2

      real ztlay (ndlon)
      real ztlev (ndlon)

c----------------------------------------------------------------------
      do jnu=1,nir
c----------------------------------------------------------------------
c         1.1   levels and layers from surface to nlaylte
c               ---------------------------------------

      do jk = 1 , nlaylte
        do jl = 1 , kdlon

c                                                  level planck function 
c                                                  ---------------------
      ztlev(jl)=(tlev(jl,jk)-tstand)/tstand    ! tstand = 200k

      blev(jl,jnu,jk)   = xp(1,jnu)
     .                    +ztlev(jl)*(xp(2,jnu) 
     .                    +ztlev(jl)*(xp(3,jnu) 
     .                    +ztlev(jl)*(xp(4,jnu) 
     .                    +ztlev(jl)*(xp(5,jnu)
     .                    +ztlev(jl)*(xp(6,jnu)    )))))

c                                                  layer planck function 
c                                                  ---------------------
      ztlay(jl)=(tlay(jl,jk)-tstand)/tstand

      blay(jl,jnu,jk) = xp(1,jnu)
     .                  +ztlay(jl)*(xp(2,jnu)
     .                  +ztlay(jl)*(xp(3,jnu)
     .                  +ztlay(jl)*(xp(4,jnu)
     .                  +ztlay(jl)*(xp(5,jnu)
     .                  +ztlay(jl)*(xp(6,jnu)    )))))

c                                               planck function gradient
c                                               ------------------------
      dblay(jl,jnu,jk) = xp(2,jnu)
     .                   +ztlay(jl)*(2*xp(3,jnu)
     .                   +ztlay(jl)*(3*xp(4,jnu)
     .                   +ztlay(jl)*(4*xp(5,jnu)
     .                   +ztlay(jl)*(5*xp(6,jnu)  ))))
      dblay(jl,jnu,jk) = dblay(jl,jnu,jk)/tstand

        enddo
      enddo
     
c----------------------------------------------------------------------
c         1.2   top of the atmosphere and surface
c               --------------------------------

      do jl = 1 , kdlon
c                                                  top of the atmosphere
c                                                  ---------------------
      ztlev(jl) = (tlev(jl,nlaylte+1)-tstand)/tstand

      blev(jl,jnu,nlaylte+1) = xp(1,jnu)
     .                       +ztlev(jl)*(xp(2,jnu)
     .                       +ztlev(jl)*(xp(3,jnu)
     .                       +ztlev(jl)*(xp(4,jnu)
     .                       +ztlev(jl)*(xp(5,jnu)
     .                       +ztlev(jl)*(xp(6,jnu)    )))))
      btop(jl,jnu) = blev(jl,jnu,nlaylte+1)

c                                                                surface
c                                                                -------
      ztlay(jl) = (tlev(jl,1)+dt0(jl)-tstand)/tstand

      bsurf(jl,jnu) = xp(1,jnu)
     .               +ztlay(jl)*(xp(2,jnu)
     .               +ztlay(jl)*(xp(3,jnu)
     .               +ztlay(jl)*(xp(4,jnu)
     .               +ztlay(jl)*(xp(5,jnu)
     .               +ztlay(jl)*(xp(6,jnu)     )))))

      enddo

c----------------------------------------------------------------------
c         1.3   Gradients in sub-layers
c               -----------------------

      do jk=1,nlaylte
        jk2 = 2 * jk
        jk1 = jk2 - 1 
          do jl=1,kdlon
            dbsublay(jl,jnu,jk1)=blay(jl,jnu,jk)-blev(jl,jnu,jk)
            dbsublay(jl,jnu,jk2)=blev(jl,jnu,jk+1)-blay(jl,jnu,jk)
          enddo
      enddo

c----------------------------------------------------------------------
      enddo           ! (do jnu=1,nir)
c----------------------------------------------------------------------

      return
      end
