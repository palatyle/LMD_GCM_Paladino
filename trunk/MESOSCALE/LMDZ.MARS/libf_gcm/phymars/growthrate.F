      subroutine growthrate(timestep,t,p,ph2o,psat,seq,r,Cste)

      IMPLICIT NONE

c=======================================================================
c
c     Determination of the water ice crystal growth rate
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

c
c   arguments:
c   ----------

      REAL timestep
      REAL t    ! temperature in the middle of the layer (K)
      REAL p    ! pressure in the middle of the layer (K)
      REAL*8 ph2o ! water vapor partial pressure (Pa)
      REAL*8 psat ! water vapor saturation pressure (Pa) 
      REAL r    ! crystal radius before condensation (m)
      REAL seq  ! Equilibrium saturation ratio
      REAL dr   ! crystal radius variation (m)

c   local:
c   ------

      REAL molco2,molh2o
      REAL Mco2,Mh2o,rho_i,sigh2o
      REAL nav,rgp,kbz,pi,To

c     Effective gas molecular radius (m)
      data molco2/2.2e-10/   ! CO2
c     Effective gas molecular radius (m)
      data molh2o/1.2e-10/   ! H2O
c     Molecular weight of CO2
      data Mco2/44.e-3/            ! kg.mol-1
c     Molecular weight of H2O
      data Mh2o/18.e-3/            ! kg.mol-1
c     surface tension of ice/vapor
      data sigh2o/0.12/	     ! N.m
c     Ice density
      data rho_i/917./        ! kg.m-3 also defined in initcld.f
c     Avogadro number
      data nav/6.023e23/
c     Perfect gas constant
      data rgp/8.3143/
c     Boltzman constant
      data kbz/1.381e-23/
c     pi number
      data pi/3.141592654/
c     Reference temperature, T=273,15 K 
      data To/273.15/

      REAL k,Lv                 
      REAL knudsen           ! Knudsen number (gas mean free path/particle radius)
      REAL a,Dv,lambda       ! Intermediate computations for growth rate
      REAL*8 Rk,Rd
      REAL Cste, rf

c-----------------------------------------------------------------------
c      Ice particle growth rate by diffusion/impegement of water molecules
c                r.dr/dt = (S-Seq) / (Seq*Rk+Rd)
c        with r the crystal radius, Rk and Rd the resistances due to 
c        latent heat release and to vapor diffusion respectively 
c----------------------------------------------------------------------- 

c     - Equilibrium saturation accounting for KeLvin Effect
       seq=exp(2*sigh2o*Mh2o/(rho_i*rgp*t*r))

c     - Thermal conductibility of CO2
      k  = (0.17913 * t - 13.9789) * 4.184e-4
c     - Latent heat of h2o (J.kg-1)
      Lv = (2834.3 - 0.28 * (t-To) - 0.004 * (t-To)**2 ) * 1.e+3

c     - Constant to compute gas mean free path
c     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))
      a = 0.707*rgp/(4 * pi* molco2**2  * nav)

c     - Compute Dv, water vapor diffusion coefficient
c       accounting for both kinetic and continuum regime of diffusion,
c       the nature of which depending on the Knudsen number.

      Dv = 1./3. * sqrt( 8*kbz*t/(pi*Mh2o/nav) )* kbz * t / 
     &   ( pi * p * (molco2+molh2o)**2 * sqrt(1.+Mh2o/Mco2) )

      knudsen = t / p * a / r
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
      Dv      = Dv / (1. + lambda * knudsen)

c     - Compute Rk
      Rk = Lv**2 * rho_i * Mh2o / (k*rgp*t**2.)
c     - Compute Rd
      Rd = rgp * t *rho_i / (Dv*psat*Mh2o)

c     - Compute Cste=rdr/dt, then r(t+1)= sqrt(r(t)**2.+2.*Cste*dt)
       Cste = 1. / (seq*Rk+Rd)
c      Cste = (ph2o/psat-seq) / (seq*Rk+Rd)
c      rf   = sqrt( max( r**2.+2.*Cste*timestep , 0. ) )
c      dr   = rf-r

      RETURN
      END

