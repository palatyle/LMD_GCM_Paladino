      subroutine n_methane(ngrid,nq,nbin,
     *                     dt,pl,tl,aerad,
     *                     q,qprime)

      implicit none
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

c  Arguments
c  ---------

      integer ngrid,nq,nbin

      REAL dt              ! physical time step (s)
      REAL pl(ngrid,nz)    ! pressure at each level (mbar)
      REAL tl(ngrid,nz)    ! temperature at each level (K)
      REAL aerad(nbin)     ! Radius array

c    Tracers :
      REAL q(ngrid,nz,nq)         ! tracer (kg/kg)
      REAL qprime(ngrid,nz,nbin)  ! tracer (kg/kg)


c     Local variables
c     ---------------

      integer ntyp
      parameter (ntyp=3)   ! 3 = 1 aerosol + 1 noyau + 1 glace 

      real n_aer(nz,nbin,ntyp)
      real ch4vap(nz)

      integer itrac
      integer ig,i,j,k,l,n   ! Loop integers
      integer ilay,iq

c  Treatment 
c  ---------

      DO ig = 1 , NGRID

c     Set up the aerosol array
        do j = 1, ntyp
          do k = 1, nbin
            itrac = (j-1) * nbin + k
            do l = 1, nz
              n_aer(l,k,j) = max(q(ig,l,itrac),0.)
            enddo 
          enddo
        enddo

c     Set up the methane vapor array
        do l = 1, nz
          ch4vap(l) = q(ig,l,nq)
        enddo



        call nucleacond(ngrid,nbin,dt,ig,pl,tl,aerad,
     &                 n_aer,qprime,ch4vap)


c       Update q arrays
        do j = 1, ntyp
          do k = 1, nbin
            itrac = (j-1) * nbin + k
            do l = 1, nz
              q(ig,l,itrac) = n_aer(l,k,j)
            enddo 
          enddo
        enddo

c     Update methane vapor array
        do l = 1, nz
          q(ig,l,nq) = ch4vap(l)
        enddo

      ENDDO

      return
      END

****************************************************************
      subroutine nucleacond(ngrid,nbin,dt,ig,
     &                      pl,tl,aerad,n_aer,qprime,ch4vap)
*                                                              *
*     This routine updates species concentrations due          *
*     to both nucleation and condensation-induced variations.  *
*     Gain and loss rates associated to each one of these      *
*     processes are computed separately in other routines.     *
*                                                              *
****************************************************************

      implicit none
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

      integer ng,nalt
      parameter(ng=1,nalt=llm)

      real lv
      common/lheat/lv

c  Arguments
c  ---------

      integer ngrid,nbin
      integer ig
      integer ntyp
      parameter (ntyp=3)

      real dt                    ! Global time step
      real pl(ngrid,nz),tl(ngrid,nz)
      real aerad(nbin)
      real ch4vap(nz)            ! Methane vapor mass mixing ratio (kg/m3)
      real ch4vap_old
      real n_aer(nz,nbin,ntyp)  ! number concentrations of particle/each size bin 
      real qprime(ngrid,nz,nbin)  ! tracer (kg/kg)
      REAL total1(nz),total11(nz),total2(nz),total22(nz)
      REAL dmsm,mtot


c  Local
c  -----

      integer i,j,k,l,n,iindice,iselec

      real dQc           ! Amount of condensed methane (kg/m3) during timestep
      real*8 sat_ratio   ! Methane saturation ratio over liquid
      real*8 sat_ratmix   ! Methane saturation ratio over liquid
      real*8 pch4        ! Methane partial pressure (Pa) 
      real qsat          ! Methane mass mixing ratio at saturation (kg/kg of air)
      real qsatmix          ! Methane mass mixing ratio at saturation (kg/kg of air)
      real*8 rate(nbin)    ! Heterogeneous Nucleation rate (s-1)
      real*8 elim          

      real nsav(nbin,ntyp)
      real dn(nbin,ntyp) 
      real rad(nbin)     ! Radius of droplets in each size bin
      real*8 gr(nbin)      ! Growth rate in each bin
      real radius        ! Radius of droplets after growth
      real Qs            ! Mass of condensate required to reach saturation
      real newsat
      real vol(nbin)

      real press
      real sig,temp,seq(nbin)
      real Ctot,up,dwn,newvap,gltot

      real rho_a,cap
      real tempref
      real frac
      real xtime,xtime_prime
      real dQc2

c     Variables for latent heat release
      real lw
      data lw / 510.e+3/
      save lw


c  Treatment
c  ---------

c  Treatment
c  ---------
      do i = 1, nbin
        vol(i) = 4./3. * pi * aerad(i)**3.
      enddo

        do l = 1, nz
        total1(l)=0. !solide
        do k = 1, nbin
        total1(l)=total1(l)+n_aer(l,k,2)*rhoi_ch4
        enddo
        total2(l)=ch4vap(l)
        enddo

c     Start loop over heights
      DO 100 l = 1, nz
   
        iindice=0                ! mettre l'indice à 0

        temp   = tl(ig,l)
        press  = pl(ig,l)
        tempref=temp

c       Save the values of the particle arrays before condensation
        do j = 1, ntyp
          do i = 1, nbin
            nsav(i,j) = n_aer(l,i,j)
          enddo
        enddo


 99     continue   ! DEBUT DE LA BOUCLE CONDITONNELLE

        call ch4sat(temp,press,qsat)
        qsat=qsat*0.85
        qsatmix=qsat*0.85

c       Get the partial presure of methane vapor and its saturation ratio
        pch4       = ch4vap(l) * (Mn2/Mch4) * press
        sat_ratio  = ch4vap(l) / qsat
        sat_ratmix = ch4vap(l) / qsatmix

c       Get the rates of nucleation
        call nuclea(nbin,aerad,pch4,temp,sat_ratio,rate)

c       Get the growth rates of condensation/sublimation
        up   = ch4vap(l)
        dwn  = 1.
        Ctot = ch4vap(l)
        DO i = 1, nbin

        if (n_aer(l,i,3).eq.0) then
          rad(i) = aerad(i)
        else
          rad(i) = ((n_aer(l,i,2)/n_aer(l,i,3) +
     &     qprime(ig,l,i)/n_aer(l,i,3)
     &    +vol(i))*0.75/pi)**(1./3.)
        endif
       

*       Equilibrium saturation ratio (due to curvature effect)
        seq(i) = exp( 2.*sig(temp)*Mch4 / (rhoi_ch4*rgp*temp*rad(i)) )

        call growthrate(dt,temp,press,pch4,
     &                  sat_ratmix,seq(i),rad(i),gr(i))
        up = up + dt * gr(i) * 4. * pi * rhoi_ch4 * rad(i) * seq(i)
     *                 * nsav(i,3) 
        dwn= dwn+ dt * gr(i) * 4. * pi * rhoi_ch4 * rad(i) / qsat
     *                 * nsav(i,3) 
        Ctot= Ctot + rhoi_ch4 * nsav(i,2)

        ENDDO

        newvap = min(up/dwn,Ctot)  
        newvap = max(newvap,0.)  

        gltot=0.
        DO i = 1, nbin
          gr(i)  = gr(i) * ( newvap/qsat - seq(i) )
          if(nsav(i,2).le.0. .and. gr(i).le.0.) then
              n_aer(l,i,2) = 0.
          else
          n_aer(l,i,2) = nsav(i,2) + dt * gr(i) * 4. * pi * rad(i)
     *                               * n_aer(l,i,3)
          if (n_aer(l,i,2).le.0.) then
            n_aer(l,i,1) = n_aer(l,i,1) + n_aer(l,i,3)
            n_aer(l,i,2) = 0.
            n_aer(l,i,3) = 0.
          endif
          gltot=n_aer(l,i,2)*rhoi_ch4+gltot
          endif 

        ENDDO
    
c       Determine the mass of exchanged methane

        dQc = 0.
        DO i = 1, nbin
          dQc = dQc - rhoi_ch4 * ( n_aer(l,i,2) - nsav(i,2) )
        ENDDO


c       Arrays resetted to their initial value before condensation

        do j = 1, ntyp
          do i = 1, nbin
            dn(i,j)      = n_aer(l,i,j) - nsav(i,j)
            n_aer(l,i,j) = nsav(i,j)
          enddo
        enddo


c       Update the c arrays.
c       nucleation & cond/eva tendencies added together.

        do i=1,nbin
          elim         = dt * rate(i)
          n_aer(l,i,1) = n_aer(l,i,1) / (1.+elim)
          n_aer(l,i,3) = n_aer(l,i,3) + elim * n_aer(l,i,1) + dn(i,3)
          n_aer(l,i,1) = n_aer(l,i,1) + dn(i,1)
          n_aer(l,i,2) = n_aer(l,i,2) + dn(i,2)
         if(n_aer(l,i,2).lt.0.) n_aer(l,i,2)=0.

        enddo

        dQc = 0.
        DO i = 1, nbin    ! dQc <0 si glace produite !!
          dQc = dQc - rhoi_ch4 * ( n_aer(l,i,2) - nsav(i,2) )
        ENDDO


       ch4vap(l)    = ch4vap(l) + dQc

100   CONTINUE

        do l = 1, nz
        total11(l)=0.
        do k = 1, nbin
          total11(l)=total11(l)+n_aer(l,k,2)*rhoi_ch4
        enddo
        total22(l)=ch4vap(l) 
        enddo


      return
      end


*******************************************************
* 						      *
      subroutine nuclea(nbin,aerad,pch4,temp,sat,nucrate)
*   This subroutine computes the nucleation rate      *
*   as given in Pruppacher & Klett (1978) in the      *
*   case of water ice forming on a solid substrate.   *
*     Definition refined by Keese (jgr,1989)	      *
*						      *
*******************************************************

      implicit none
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

      integer nbin
      real aerad(nbin)

      real*8 nucrate(nbin)
      real*8 pch4
      real   temp
      real*8 sat

      integer l,i
      real*8 nch4
      real sig            ! Water-ice/air surface tension  (N.m)
      real*8 rstar        ! Radius of the critical germ (m)
      real*8 gstar        ! # of molecules forming a critical embryo
      real*8 x            ! Ratio rstar/radius of the nucleating dust particle
      real fistar         ! Activation energy required to form a critical embryo (J)
      real*8 zeldov       ! Zeldovitch factor (no dim)
      real*8 fshape       ! function defined at the end of the file
      real*8 deltaf

      real nus
      data nus/1.e+13/       ! Jump frequency of a molecule (s-1)
      real m0
      data m0/2.6578e-26/     ! Weight of a methane molecule (kg)
      real vo1
c     data vo1/6.254e-29/     ! Volume of a methane molecule (m3)
      data vo1/3.7649e-5/     ! Volume of a methane molecule (m3)

      real desorp
      data desorp/0.288e-19/ ! Activation energy for desorption of water on a dust-like substrate (J/molecule)
      real surfdif
      data surfdif/0.288e-20/! Estimated activation energy for surface diffusion of water molecules (J/molecule)

      IF (sat .GT. 1.) THEN    ! minimum condition to activate nucleation

        nch4    = pch4 / kbz / temp
        rstar  = 2. * sig(temp) * vo1 / (rgp*temp*log(sat))
        gstar  = 4. * nav * pi * (rstar**3) / (3.*vo1)
c       Loop over size bins
        do i=1,nbin
          x      = aerad(i) / rstar
          x      = aerad(imono) / rstar  ! attention r(5)=monomeres

          fistar = (4./3.*pi) * sig(temp) * (rstar**2.) *
     &             fshape(mtetach4,x)
          deltaf = min( max((2.*desorp-surfdif-fistar)/(kbz*temp)
     &           , -100.), 100.)
          if (deltaf.eq.-100.) then 
            nucrate(i) = 0. 
          else
            zeldov = sqrt ( fistar / (3.*pi*kbz*temp*(gstar**2.)) )
            nucrate(i)  = zeldov * kbz* temp * rstar**2.
     &                  * 4. * pi * ( nch4*aerad(i) )**2.
     &                  / ( fshape(mtetach4,x) * nus * m0 )
     &                  * dexp(deltaf)
 
             if(i.gt.imono)  nucrate(i)= zeldov * kbz* temp * rstar**2.
     &           * 4. * pi * vrat_e**(i-imono)*(nch4*aerad(imono) )**2.
     &            / (fshape(mtetach4,x) * nus * m0 )
     &            * dexp(deltaf)


          endif
        enddo
      ELSE
        do i=1,nbin
          nucrate(i) = 0.
        enddo

      ENDIF

      return
      end

******************************************************************
        subroutine ch4sat(t,p,qsat)
* cette fonction calcule la pression de vapeur saturante du      *
* methane a une altitude donnee z par l'equation de Thek-Stiel   *
******************************************************************

        real rgp
        data rgp/8.3143/

* declaration des variables internes
* ----------------------------------
        real p,tc,pc,tr,tb,tbr,phib,ac,hbr,hvb,avb,a,b
        real qsat

        tc = 190.53
        pc = 45.96 * 0.986923
        tr = t / tc
        tb = 111.63
        tbr= tb / tc
        phib = -35. + 36./tbr + 42.*log(tbr)-tbr**6
        ac   = (0.315*phib + log(pc))/(0.0838*phib-log(tbr))
        hbr  = tbr * log(pc) / (1.-tbr)
        hvb  = rgp*tb*(0.4343*log(pc)-0.68859+0.89584*tbr)/(0.37691-
     s      0.37306*tbr+0.14878/(pc*tbr**2))
        avb  = hvb / (rgp*tc*(1.-tbr)**(0.375))
        a    = 5.2691 + 2.0753*avb - 3.1738*hbr
        b    = 1.042 * ac - 0.46284 * avb

        qsat=pc*exp(avb*(1.14893-1./tr-0.11719*tr-0.03174*tr**2
     s        -0.375*log(tr))+b*((tr**a-1.)/a+0.04*(1./tr-1.)))
        qsat=qsat*1e+5/0.986923
        !  ---> qsat en Pa

        qsat=qsat* 16.0 / (28.0*p)   ! kg/kg
        


        return
        end

c=======================================================================
      subroutine growthrate(timestep,temp,press,pch4,sat,seq,r,Cste)
c
c     Determination of the droplet growth rate
c
c=======================================================================

      IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

c-----------------------------------------------------------------------
c   declarations:
c   -------------

      common/lheat/Lv

c
c   arguments:
c   ----------

      REAL timestep
      REAL temp    ! temperature in the middle of the layer (K)
      REAL press    ! pressure in the middle of the layer (K)
      REAL*8 pch4 ! Methane vapor partial pressure (Pa)
      REAL*8 sat  ! saturation ratio 
      REAL r    ! crystal radius before condensation (m)
      REAL seq  ! Equilibrium saturation ratio

c   local:
c   ------

      REAL psat
      REAL moln2,molch4
      REAL To,wch4,tch4,ftm

c     Effective gas molecular radius (m)
      data moln2/1.75e-10/   ! N2
c     Effective gas molecular radius (m)
      data molch4/2.e-10/   ! CH4

      data tch4/190.4/
      data wch4/1.1e-2/    ! Reid et al (Eq 7-9.4 + Appendix Compound [116])

      REAL k,Lv                 
      REAL knudsen           ! Knudsen number (gas mean free path/particle radius)
      REAL a,Dv,lambda,Rk,Rd ! Intermediate computations for growth rate
      REAL*8 Cste

c-----------------------------------------------------------------------
c      Ice particle growth rate by diffusion/impegement of molecules
c                r.dr/dt = (S-Seq) / (Seq*Rk+Rd)
c        with r the crystal radius, Rk and Rd the resistances due to 
c        latent heat release and to vapor diffusion respectively 
c----------------------------------------------------------------------- 

      psat = pch4 / sat

c     - Thermal conductibility of N2
      k = ( 2.857e-2 * temp - 0.5428  ) * 4.184e-3
c     - Latent heat of ch4 (J.kg-1)
      Lv = 510.e+3
      ftm=1.-temp/tch4
      if (ftm.le.1.e-3)  ftm=1.e-3
      Lv=8.314*tch4*(7.08*ftm**0.354+10.95*wch4*ftm**0.456)/16.e-3
      

c     - Constant to compute gas mean free path
c     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))

      a = 0.707*rgp/(4 * pi* (moln2*1.e10)**2  * (nav*1.e-20))

c     - Compute Dv, methane vapor diffusion coefficient
c       accounting for both kinetic and continuum regime of diffusion,
c       the nature of which depending on the Knudsen number.

      Dv = 1./3. * sqrt( 8*rgp*temp/(pi*Mch4) )* (kbz*1.e20) * temp / 
     &   (pi*press*(moln2*1.e10+molch4*1.e10)**2 * sqrt(1.+Mch4/Mn2) )

      knudsen = temp / press * a / r
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
      Dv      = Dv / (1. + lambda * knudsen)

c     - Compute Rk
      Rk = Lv**2 * rhoi_ch4 * Mch4 / (k*rgp*temp**2.)
c     - Compute Rd
      Rd = rgp * temp *rhoi_ch4 / (Dv*psat*Mch4)

c     - Compute:      rdr/dt = Cste * (S-Seq)
      Cste = 1. / (seq*Rk+Rd)

      RETURN
      END


*********************************************************
      real function sig(t)
*    this function computes the surface tension (N.m)   *
*   between methane and air as a function of temp.    *
*********************************************************

      real t

      sig = ( t/4. + 41.) * 1e-3

      return
      end

*********************************************************
      real*8 function fshape(cost,rap)
*        function computing the f(m,x) factor           *
* related to energy required to form a critical embryo  *
*********************************************************

      implicit none

      real cost
      real*8 rap
      real*8 phi
      real*8 a,b,c

      phi = sqrt( 1. - 2.*cost*rap + rap**2 )

      a = 1. + ( (1.-cost*rap)/phi )**3
      b = (rap**3) * (2.-3.*(rap-cost)/phi+((rap-cost)/phi)**3)
      c = 3. * cost * (rap**2) * ((rap-cost)/phi-1.)

      fshape = 0.5*(a+b+c)

      if (rap.gt.3000.) fshape = ((2.+cost)*(1.-cost)**2)/4.


      return
      end

