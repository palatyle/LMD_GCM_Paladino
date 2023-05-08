      subroutine n_ethane(ngrid,nq,nbin,
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
      parameter (ntyp=3)

      real n_aer(nz,nbin,ntyp)
      real c2h6vap(nz)

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
          c2h6vap(l) = q(ig,l,nq)
        enddo

        call nucleacond2(ngrid,nbin,dt,ig,pl,tl,aerad,
     &                             n_aer,qprime,c2h6vap)

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
      q(ig,l,nq) = c2h6vap(l)

c     if(q(ig,l,nq).le.0.) print*,'in nuc2h6',ig,l, q(ig,l,nq) 

      enddo

      ENDDO

      return
      END

****************************************************************
      subroutine nucleacond2(ngrid,nbin,dt,ig,
     *                      pl,tl,aerad,n_aer,qprime,c2h6vap)
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

* dejà évalué dans n_methane !
      REAL x1(nalt)
      REAL x2(nalt)
      REAL x3(nalt)
      REAL icefrac(nalt)
      REAL pmixch4(nalt)
      REAL pmixc2h6(nalt)
      REAL pmixn2(nalt)

      common/lheat/lv
      common/mixing/x1,x2,x3,icefrac,
     &         pmixch4,pmixc2h6,pmixn2



c  Arguments
c  ---------

      integer ngrid,nbin
      integer ig
      integer ntyp
      parameter (ntyp=3)

      real dt                    ! Global time step
      real pl(ngrid,nz),tl(ngrid,nz)
      real aerad(nbin)
      real c2h6vap(nz)            ! Methane vapor mass mixing ratio (kg/m3)
      real c2h6vap_old
      real n_aer(nz,nbin,ntyp)  ! number concentrations of particle/each size bin 
      real qprime(ngrid,nz,nbin)  ! tracer (kg/kg)
      REAL total1(nz),total11(nz),total2(nz),total22(nz)
      REAL dmsm,mtot





c  Local
c  -----

      integer i,j,k,l,n,iindice,iselec

      real dQc           ! Amount of condensed methane (kg/m3) during timestep
      real*8 sat_ratio     ! Methane saturation ratio over liquid
      real*8 sat_ratmix  ! Methane saturation ratio over liquid
      real*8 pc2h6         ! Methane partial pressure (Pa) 
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
      real sig2,temp,seq(nbin)
      real Ctot,up,dwn,newvap,gltot

      real temp0,temp1,temp2,last_temp
      real qsat1,sat_ratio1,tempf(0:10),sat_ratiof(0:10)
      real rho_a,cap
      real tempref
      real xtime,xtime_prime



c     Variables for latent heat release
      real lw
      data lw / 581.e+3/
      save lw


c  Treatment
c  ---------
      do i = 1, nbin
        vol(i) = 4./3. * pi * aerad(i)**3.
      enddo

      do l = 1, nz
      total1(l)=0. !solide
      do k = 1, nbin
      total1(l)=total1(l)+n_aer(l,k,2)*rhoi_c2h6
      enddo
      total2(l)=c2h6vap(l)
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


 99     continue 


        call c2h6sat(temp,press,qsat)
        qsatmix=qsat
 
c  quantité  pmixc2h6(l) déjà calculé dans cnuages.F et passé dans un common
c       qsat=pmixc2h6(l)*30.0 / (28.0*p)
c       qsatmix=pmixc2h6(l)*30.0 / (28.0*p)
 
c       Get the partial presure of methane vapor and its saturation ratio
        pc2h6      = c2h6vap(l) * (Mn2/Mc2h6) * press
        sat_ratio  = c2h6vap(l) / qsat
        sat_ratmix = c2h6vap(l) / qsatmix

c       Get the rates of nucleation
        call nuclea2(nbin,aerad,pc2h6,temp,sat_ratio,rate)

c       Get the growth rates of condensation/sublimation
        up   = c2h6vap(l)
        dwn  = 1.
        Ctot = c2h6vap(l)
        DO i = 1, nbin

        if (n_aer(l,i,3).eq.0) then
         rad(i) = aerad(i)
        else
         rad(i) = ((n_aer(l,i,2)/n_aer(l,i,3) + 
     &   qprime(ig,l,i)/n_aer(l,i,3)
     &   +vol(i))*0.75/pi)**(1./3.)
        endif


*       Equilibrium saturation ratio (due to curvature effect)
        seq(i) = exp( 2.*sig2(temp)*Mc2h6 /(rhoi_c2h6*rgp*temp*rad(i)))

        call growthrate2(dt,temp,press,pc2h6,
     &                   sat_ratmix,seq(i),rad(i),gr(i))
        up = up + dt * gr(i) * 4. * pi * rhoi_c2h6 * rad(i) * seq(i)
     *                 * nsav(i,3) 
        dwn= dwn+ dt * gr(i) * 4. * pi * rhoi_c2h6 * rad(i) / qsat
     *                 * nsav(i,3) 
        Ctot= Ctot + rhoi_c2h6 * nsav(i,2)

        ENDDO

        newvap = min(up/dwn,Ctot)
        newvap = max(newvap,0.)

        gltot = 0.
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
          gltot=n_aer(l,i,2)*rhoi_c2h6+gltot
          endif

        ENDDO
   
c       Determine the mass of exchanged methane

        dQc = 0.
        DO i = 1, nbin
          dQc = dQc - rhoi_c2h6 * ( n_aer(l,i,2) - nsav(i,2) )
        ENDDO

c       Update the methane vapor mixing ratio implied by 
c       the cond/eva processes.


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
        DO i = 1, nbin
          dQc = dQc - rhoi_c2h6 * ( n_aer(l,i,2) - nsav(i,2) )
        ENDDO


        c2h6vap(l)  = c2h6vap(l) + dQc

100   CONTINUE

      do l = 1, nz
      total11(l)=0. 
      do k = 1, nbin
      total11(l)=total11(l)+n_aer(l,k,2)*rhoi_c2h6
      enddo
      total22(l)=c2h6vap(l)
      enddo

      return
      end


*******************************************************
* 						      *
      subroutine nuclea2(nbin,aerad,pc2h6,temp,sat,nucrate)
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
      real*8 pc2h6
      real   temp
      real*8 sat

      integer l,i
      real*8 nc2h6
      real sig2            ! Water-ice/air surface tension  (N.m)
      real*8 rstar        ! Radius of the critical germ (m)
      real*8 gstar        ! # of molecules forming a critical embryo
      real*8 x            ! Ratio rstar/radius of the nucleating dust particle
      real fistar         ! Activation energy required to form a critical embryo (J)
      real*8 zeldov       ! Zeldovitch factor (no dim)
      real*8 fshape2       ! function defined at the end of the file
      real*8 deltaf

      real nus
      data nus/1.e+13/       ! Jump frequency of a molecule (s-1)
      real m0
      data m0/4.983e-26/     ! Weight of a methane molecule (kg)
      real vo1
c     data vo1/9.094e-29/
      data vo1/5.4746e-5/    ! Volume of 1 mole
      real desorp
      data desorp/0.288e-19/ ! Activation energy for desorption of water on a dust-like substrate (J/molecule)
      real surfdif
      data surfdif/0.288e-20/! Estimated activation energy for surface diffusion of water molecules (J/molecule)

      IF (sat .GT. 1.) then    ! minimum condition to activate nucleation

        nc2h6    = pc2h6 / kbz / temp
        rstar  = 2. * sig2(temp) * vo1 / (rgp*temp*log(sat))
        gstar  = 4. * nav * pi * (rstar**3) / (3.*vo1)
c       Loop over size bins
        do i=1,nbin
          x      = aerad(i) / rstar
          x      = aerad(imono) / rstar  ! r(5)=monomere
          fistar = (4./3.*pi) * sig2(temp) * (rstar**2.) 
     &     *fshape2(mtetac2h6,x)
          deltaf = min( max((2.*desorp-surfdif-fistar)/(kbz*temp)
     &           , -100.), 100.)
          if (deltaf.eq.-100.) then 
            nucrate(i) = 0. 
          else
            zeldov = sqrt ( fistar / (3.*pi*kbz*temp*(gstar**2.)) )
            nucrate(i)  = zeldov * kbz* temp * rstar**2.
     &                  * 4. * pi * ( nc2h6*aerad(i) )**2.
     &                  / ( fshape2(mtetac2h6,x) * nus * m0 )
     &                  * dexp(deltaf)
            if(i.gt.imono)  nucrate(i)= zeldov * kbz* temp * rstar**2.
     &          * 4. * pi * vrat_e**(i-imono)*(nc2h6*aerad(imono) )**2.
     &          / (fshape2(mtetac2h6,x) * nus * m0 )
     &          * dexp(deltaf)

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
        subroutine c2h6sat(t,p,qsat)
*                                                                 *
* cette fonction calcule la pression de vapeur saturante de l'    *
* ethane a une altitude donnee z par Reid et al., p657            *
*                                                                 *
* Compatible avec Barth et al., dans l'intervalle 30-90K          *
*                                                                 *
*                                                                 *
******************************************************************

        real rgp
        data rgp/8.3143/

* declaration des variables internes
* ----------------------------------

        real qsat,t,p
         
c       qsat=10.01-1085.0/(t-0.561)
c       qsat=10.**qsat/760.*1.e5


        a=-6.34307
        b=1.011630
        c=-1.19116
        d=-2.03539
        pc=48.8*1.e5
        tc=305.4

        x=(1.-t/tc)
        if(x.gt. 0.) qsat=(1-x)**(-1)*(a*x+b*x**1.5+c*x**3.+d*x**6.)
        if(x.le. 0.) qsat=a*x/abs(1.-x)     ! approx pour  t > tc
        qsat=pc*exp(qsat)

        qsat=qsat* 30.0 / (28.0*p)  ! kg/kg

        return
        end

c=======================================================================
      subroutine growthrate2(timestep,temp,press,pc2h6,sat,seq,r,Cste)
c
c     Determination of the droplet growth rate
c
c=======================================================================

      IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"

c-----------------------------------------------------------------------
C   DECLARATIONS:
c   -------------

      common/lheat/Lv


c
c   arguments:
c   ----------

      REAL timestep
      REAL temp    ! temperature in the middle of the layer (K)
      REAL press   ! pressure in the middle of the layer (K)
      REAL*8 pc2h6 ! Methane vapor partial pressure (Pa)
      REAL*8 sat  ! saturation ratio 
      REAL r    ! crystal radius before condensation (m)
      REAL seq  ! Equilibrium saturation ratio

c   local:
c   ------

      REAL psat
      REAL moln2,molc2h6
      REAL To,tc2h6,wc2h6       ! Reid et al., (eq 7-9.4 + Appendix compound  [168])
      REAL fte

c     Effective gas molecular radius (m)
      data moln2/1.75e-10/   ! N2
c     Effective gas molecular radius (m)
      data molc2h6/2.22e-10/   ! C2H6
c     Temperature critique  + omega
      data tc2h6/305.4/
      data wc2h6/9.9e-2/

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

      psat = pc2h6 / sat

c     - Thermal conductibility of N2
      
      k = ( 2.857e-2 * temp - 0.5428  ) * 4.184e-3
      
      
c     - Latent heat of c2h6 (J.kg-1)
      Lv =581.e3                      ! eq (7-9.4) Reid et al. 
      fte=(1.-temp/tc2h6)
      if (fte.le.1.e-3)  fte=1.e-3
      Lv=8.314*tc2h6*(7.08*fte**0.354+10.95*wc2h6*fte**0.456)/30.e-3

      

c     - Constant to compute gas mean free path
c     l= (T/P)*a, with a = (  0.707*8.31/(4*pi*molrad**2 * avogadro))

      a = 0.707*rgp/(4 * pi* (moln2*1.e10)**2  * (nav*1.e-20))

c     - Compute Dv, methane vapor diffusion coefficient
c       accounting for both kinetic and continuum regime of diffusion,
c       the nature of which depending on the Knudsen number.

      Dv = 1./3. * sqrt( 8*rgp*temp/(pi*Mc2h6) )* (kbz*1.e20) * temp / 
     & (pi*press*(moln2*1.e10+molc2h6*1.e10)**2 * sqrt(1.+Mc2h6/Mn2) )

      knudsen = temp / press * a / r
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen)
      Dv      = Dv / (1. + lambda * knudsen)

c     - Compute Rk
      Rk = Lv**2 * rhoi_c2h6 * Mc2h6 / (k*rgp*temp**2.)
c     - Compute Rd
      Rd = rgp * temp *rhoi_c2h6 / (Dv*psat*Mc2h6)

c     - Compute:      rdr/dt = Cste * (S-Seq)
      Cste = 1. / (seq*Rk+Rd)

      RETURN
      END


*********************************************************
      real function sig2(t)
*    this function computes the surface tension (N.m)   *
*   between methane and air as a function of temp.    *
*   Eq (12-3-6 et 12-3-7 Reid et al.)
*
*********************************************************

      real t

      pc=48.8*1.01325e5
      tc=305.4
      tb=184.6
      tr=t/tc
      tbr=tb/tc
      if(t.gt.305.0) then 
         tr=305./tc
      endif

      sig2=0.1196*(1.+(tbr*alog(pc/1.01325))/(1.-tbr))-0.279
      sig2=pc**(2./3.)*tc**(1./3.)*sig2*(1.-tr)**(11./9.)
      sig2=sig2*1.e-8




c     sig2 = 21.157*((305.4-t)/(305.4-153.2))**(11./9.)
c     sig2=sig2*1.e-4
c
c       if(t.gt.305.0) then 
c         sig2 = 21.157*((305.4-305.)/(305.4-153.2))**(11./9.)
c         sig2=sig2*1.e-4
c       endif
c
c     print*,t,sig2,sig3

      return
      end

*********************************************************
      real*8 function fshape2(cost,rap)
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

      fshape2 = 0.5*(a+b+c)

      if (rap.gt.3000.) fshape2 = ((2.+cost)*(1.-cost)**2)/4.

      return
      end

