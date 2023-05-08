      subroutine growthrate(temp,pmid,psat,rcrystal,res)

      use tracer_mod, only: rho_ice
      USE comcstfi_h
      IMPLICIT NONE

c=======================================================================
c
c     Determination of the water ice crystal growth rate
c
c     Authors: F. Montmessin
c       Adapted for the LMD/GCM by J.-B. Madeleine (October 2011)
c       Use of resistances in the analytical function 
c            instead of growth rate - T. Navarro (2012)
c     
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "microphys.h"

c
c   arguments:
c   ----------

c     Input
      REAL temp     ! temperature in the middle of the layer (K)
      REAL pmid     ! pressure in the middle of the layer (K)
      REAL psat   ! water vapor saturation pressure (Pa) 
      REAL rcrystal ! crystal radius before condensation (m)

c     Output
      REAL res      ! growth resistance (res=Rk+Rd)


c   local:
c   ------

      REAL k,Lv                 
      REAL knudsen           ! Knudsen number (gas mean free path/particle radius)
      REAL afactor,Dv,lambda       ! Intermediate computations for growth rate
      REAL Rk,Rd
      
      

c-----------------------------------------------------------------------
c      Ice particle growth rate by diffusion/impegement of water molecules
c                r.dr/dt = (S-Seq) / (Seq*Rk+Rd)
c        with r the crystal radius, Rk and Rd the resistances due to 
c        latent heat release and to vapor diffusion respectively 
c----------------------------------------------------------------------- 

c     - Equilibrium saturation accounting for KeLvin Effect
c      seq=exp(2*sigh2o*mh2o/(rho_ice*rgp*t*r))
c      (already computed in improvedcloud.F)

c     - Thermal conductibility of CO2
      k  = (0.17913 * temp - 13.9789) * 4.184e-4
c     - Latent heat of h2o (J.kg-1)
      Lv = (2834.3 
     &        - 0.28  * (temp-To) 
     &        - 0.004 * (temp-To) * (temp-To) ) * 1.e+3

c     - Constant to compute gas mean free path
c     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
      afactor = 0.707*rgp/(4 * pi * molco2 * molco2 * nav)

c     - Compute Dv, water vapor diffusion coefficient
c       accounting for both kinetic and continuum regime of diffusion,
c       the nature of which depending on the Knudsen number.

      Dv = 1./3. * sqrt( 8*kbz*temp/(pi*mh2o/nav) )* kbz * temp / 
     &   ( pi * pmid * (molco2+molh2o)*(molco2+molh2o) 
     &        * sqrt(1.+mh2o/mco2) )
      
      knudsen = temp / pmid * afactor / rcrystal
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
      
c      Dv is not corrected. Instead, we use below coefficients coeff1, coeff2
c      Dv      = Dv / (1. + lambda * knudsen)

c     - Compute Rk
      Rk = Lv*Lv* rho_ice * mh2o / (k*rgp*temp*temp)
c     - Compute Rd
      Rd = rgp * temp *rho_ice / (Dv*psat*mh2o)
      
      
      res = Rk + Rd*(1. + lambda * knudsen)
      
      !coeff1 = real(Rk + Rd*(1. + lambda * knudsen))
      !coeff2 = real(Rk + Rd*(1. - lambda * knudsen))
      
c Below are growth rate used for other schemes
c     - Compute growth=rdr/dt, then r(t+1)= sqrt(r(t)**2.+2.*growth*dt)
c      growth = 1. / (seq*Rk+Rd)
c      growth = (ph2o/psat-seq) / (seq*Rk+Rd)
c      rf   = sqrt( max( r**2.+2.*growth*timestep , 0. ) )
c      dr   = rf-r

      RETURN
      END

