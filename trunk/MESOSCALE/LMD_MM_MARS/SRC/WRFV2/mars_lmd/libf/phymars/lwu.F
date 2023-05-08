      subroutine lwu (kdlon,kflev
     &                ,dp,plev,tlay,aerosol 
     &                ,aer_t,co2_u,co2_up
     &                ,tautotal,omegtotal,gtotal)

c----------------------------------------------------------------------
c     LWU   computes  - co2: longwave effective absorber amounts including 
c                      pressure and temperature effects 
c                     - aerosols: amounts for every band
c                                transmission for bandes 1 and 2 of co2
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c ATTENTION AUX UNITES:
c le facteur 10*g fait passer des kg m-2 aux g cm-2
c-----------------------------------------------------------------------
c! modif diffusion
c! on ne change rien a la bande CO2 : les quantites d'absorbant CO2
c! sont multipliees par 1.66
c! pview= 1/cos(teta0)=1.66
c-----------------------------------------------------------------------
 
      implicit none

#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "comcstfi.h"

#include "yomaer.h"
#include "yomlw.h"

#include "fisice.h"
#include "callkeys.h"

#include "aerice.h"
 
c----------------------------------------------------------------------
c         0.1   arguments
c               ---------                                    
c                                                            inputs:
c                                                            -------
      integer kdlon        ! part of ngrid
      integer kflev        ! part of nalyer

      real dp (ndlo2,kflev)         ! layer pressure thickness (Pa)
      real plev (ndlo2,kflev+1)     ! level pressure (Pa)
      real tlay (ndlo2,kflev)       ! layer temperature (K)
      real aerosol (ndlo2,kflev,naerkind) ! aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars.h , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.


c                                                            outputs:
c                                                            --------
      real aer_t (ndlo2,nuco2,kflev+1)   ! transmission (aer) 
      real co2_u (ndlo2,nuco2,kflev+1)   ! absorber amounts (co2)
      real co2_up (ndlo2,nuco2,kflev+1)  ! idem scaled by the pressure (co2)

      real tautotal(ndlo2,kflev,nir)  ! \   Total single scattering
      real omegtotal(ndlo2,kflev,nir) !  >  properties (Addition of the
      real gtotal(ndlo2,kflev,nir)    ! /   NAERKIND aerosols properties)

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer jl,jk,jkl,ja,n
 
      real aer_u (ndlon,nir,nflev+1)     ! absorber amounts (aer) extinction
      real co2c           ! co2 concentration (pa/pa)
      real pview          ! cosecant of viewing angle
      real pref           ! reference pressure (1013 mb = 101325 Pa)
      real tx,tx2
      real phi (ndlon,nuco2)
      real psi (ndlon,nuco2)
      real plev2 (ndlon,nflev+1)
      real zzz

      real ray,coefsize,coefsizew,coefsizeg 

c************************************************************************
c----------------------------------------------------------------------
c         0.3  Initialisation  
c               -------------

      pview = 1.66
      co2c = 0.95           
      pref = 101325.

      do jk=1,nlaylte+1
        do jl=1,kdlon
          plev2(jl,jk)=plev(jl,jk)*plev(jl,jk)
        enddo
      enddo

c----------------------------------------------------------------------
c  Computing TOTAL single scattering parameters by adding properties of
c  all the NAERKIND kind of aerosols in each IR band

      call zerophys(ndlon*kflev*nir,tautotal)
      call zerophys(ndlon*kflev*nir,omegtotal)
      call zerophys(ndlon*kflev*nir,gtotal)

      do n=1,naerkind
        do ja=1,nir
          do jk=1,nlaylte
            do jl = 1,kdlon
c TEST : to account for the varying sol/ir optical depth of ice with varying crystal size     
c      : and for varying w,g with varying crystal size     
              if (activice.and.n.eq.naerkind) then 
        ray=min( max(rice(jl,jk)*1.e+6, 1.),10.) 
            if (ja.eq.1.or.ja.eq.2) then 
          coefsize=(-0.00382*ray**3.+0.0503*ray**2.+0.03531)
     &             /QIRsQREF(ja,n)
          coefsizew=(-0.011*ray**2.+0.1824*ray-0.1283)
     &           /omegaIR(ja,n)
          coefsizeg=(-0.00122*ray**3.+0.009161*ray**2.+  
     &          0.1182*ray-0.096877)/gIR(ja,n)
                elseif (ja.eq.3) then
          coefsize=(-0.00324*ray**3.+0.0419*ray**2.+0.0459)
     &             /QIRsQREF(ja,n)
          coefsizew=(-0.01292*ray**2.+0.1963*ray-0.06566)
     &           /omegaIR(ja,n)
          coefsizeg=(0.00271*ray**3.-0.05959*ray**2.+
     &          0.4411*ray-0.2724)/gIR(ja,n)
                elseif (ja.eq.4) then
          coefsize=(-0.0003823*ray**3.+0.0104*ray**2.+0.005966)
     &             /QIRsQREF(ja,n)
              coefsizew=(-0.002*ray**3.+0.02623*ray**2.-0.014465)
     &           /omegaIR(ja,n)
          coefsizeg=(-0.0017192*ray**3.+0.0259*ray**2.-
     &          0.027692*ray+0.016099)/gIR(ja,n)
            endif
        if (coefsize.le.0.or.coefsizew.le.0.or.coefsizeg.le.0) 
     &        stop'pb dans lwu avec prp opt glace'

            tautotal(jl,jk,ja)=tautotal(jl,jk,ja) +
     &          QIRsQREF(ja,n)*aerosol(jl,jk,n)
     &          *coefsize
                omegtotal(jl,jk,ja) =  omegtotal(jl,jk,ja) + 
     &           QIRsQREF(ja,n)*aerosol(jl,jk,n)*omegaIR(ja,n)
     &          *coefsize*coefsizew
                gtotal(jl,jk,ja) =  gtotal(jl,jk,ja) + 
     &           QIRsQREF(ja,n)*aerosol(jl,jk,n)*omegaIR(ja,n)*gIR(ja,n)
     &          *coefsize*coefsizew*coefsizeg     
              else
              tautotal(jl,jk,ja)=tautotal(jl,jk,ja) +
     &           QIRsQREF(ja,n)*aerosol(jl,jk,n)
              omegtotal(jl,jk,ja) =  omegtotal(jl,jk,ja) + 
     &           QIRsQREF(ja,n)*aerosol(jl,jk,n)*omegaIR(ja,n)
              gtotal(jl,jk,ja) =  gtotal(jl,jk,ja) + 
     &           QIRsQREF(ja,n)*aerosol(jl,jk,n)*omegaIR(ja,n)*gIR(ja,n)
              endif
            enddo
          enddo
        enddo
      enddo
      do ja=1,nir
        do jk=1,nlaylte
          do jl = 1,kdlon
             gtotal(jl,jk,ja)=gtotal(jl,jk,ja)/omegtotal(jl,jk,ja)
             omegtotal(jl,jk,ja)=omegtotal(jl,jk,ja)/tautotal(jl,jk,ja)
          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         1.0   cumulative (aerosol) amounts (for every band)
c               ----------------------------

      jk=nlaylte+1
      do ja=1,nir
        do jl = 1 , kdlon
          aer_u(jl,ja,jk)=0.
        enddo
      enddo

      do jk=1,nlaylte
        jkl=nlaylte+1-jk
        do ja=1,nir
          do jl=1,kdlon
             aer_u(jl,ja,jkl)=aer_u(jl,ja,jkl+1)+ tautotal(jl,jkl,ja)
          enddo
        enddo
      enddo                 

c----------------------------------------------------------------------
c         1.0   bands 1 and 2 of co2
c               --------------------

      jk=nlaylte+1
      do ja=1,nuco2
        do jl = 1 , kdlon
          co2_u(jl,ja,jk)=0.
          co2_up(jl,ja,jk)=0.
          aer_t(jl,ja,jk)=1.
        enddo
      enddo

      do jk=1,nlaylte
                       jkl=nlaylte+1-jk
        do ja=1,nuco2
          do jl=1,kdlon

c                 introduces temperature effects on absorber(co2) amounts
c                 -------------------------------------------------------
            tx = sign(min(abs(tlay(jl,jkl)-tref),70.)
     .         ,tlay(jl,jkl)-tref)
            tx2=tx*tx
            phi(jl,ja)=at(1,ja)*tx+bt(1,ja)*tx2
            psi(jl,ja)=at(2,ja)*tx+bt(2,ja)*tx2
            phi(jl,ja)=exp(phi(jl,ja)/cst_voigt(2,ja))
            psi(jl,ja)=exp(2.*psi(jl,ja))

c                                        cumulative absorber(co2) amounts
c                                        --------------------------------
            co2_u(jl,ja,jkl)=co2_u(jl,ja,jkl+1)
     .     +         pview/(10*g)*phi(jl,ja)*dp(jl,jkl)*co2c

            co2_up(jl,ja,jkl)=co2_up(jl,ja,jkl+1)
     .     +         pview/(10*g*2*pref)*psi(jl,ja)
     .     *        (plev2(jl,jkl)-plev2(jl,jkl+1))*co2c

c                                                  (aerosol) transmission
c                                                  ----------------------
c   on calcule directement les transmissions pour les aerosols.
c   on multiplie le Qext  par 1-omega dans la bande du CO2.
c   et pourquoi pas d'abord?  hourdin@lmd.ens.fr

           zzz=pview*(1.-omegtotal(jl,jkl,ja))
           aer_t(jl,ja,jkl)=exp(-zzz*aer_u(jl,ja,jkl))
          
          enddo
        enddo
      enddo                 

c----------------------------------------------------------------------
      return
      end
