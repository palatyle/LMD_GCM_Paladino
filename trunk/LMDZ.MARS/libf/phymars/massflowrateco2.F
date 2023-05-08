c=======================================================================
      subroutine massflowrateco2(P,T,Sat,Radius,Matm,Ic)
c
c     Determination of the mass transfer rate for CO2 condensation & 
c     sublimation
c
c   inputs: Pressure (P), Temperature (T), saturation ratio (Sat),
c           particle radius (Radius), molecular mass of the atm (Matm)
c   output:  MASS FLUX Ic
c
c   Authors: C. Listowski (2014) then J. Audouard (2016-2017)
c  
c
c   Updates:
c   --------
c   December 2017 - C. Listowski - Simplification of the derivation of
c   massflowrate by using explicit formula for surface temperature,
c   No Newton-Raphson routine anymore- see comment at relevant line
c=======================================================================
      USE comcstfi_h, ONLY: pi

      implicit none

      include "microphys.h"

c   arguments: INPUT
c   ----------
      REAL,INTENT(in) :: T,Matm
      REAL*8,INTENT(in) :: SAT
      REAL,INTENT(in) :: P
      DOUBLE PRECISION,INTENT(in) :: Radius
c   arguments: OUTPUT
c   ----------
      DOUBLE PRECISION,INTENT(out) ::   Ic
c   Local Variables
c   ----------
      DOUBLE PRECISION   Tsurf
      DOUBLE PRECISION   C0,C1,C2
      DOUBLE PRECISION   kmix,Lsub,cond
      DOUBLE PRECISION   Ak

      call  coefffunc(P,T,Sat,Radius,Matm,kmix,Lsub,C0,C1,C2,Ak)
   
      Tsurf = 1./C1*dlog(Sat/Ak) + T 

      !Note by CL - dec 2017 (see also technical note)
      !The above is a simplified version of Tsurf
      !compared to the one used by Listowski et al. 2014 (Ta), where a
      !Newton-Raphson routine must be used. Approximations made by
      !considering the orders of magnitude of the different factors lead to
      !simplification of the equation 5 of Listowski et al. (2014).
      !The error compared to the exact value determined by NR iterations
      !is less than 0.6% for all sizes, pressures, supersaturations
      !relevant to present Mars. Should also be ok for most conditions
      !in ancient Mars (However, needs to be double-cheked, as explained in
      !(Listowski et al. 2013,JGR)

      cond = 4.*pi*Radius*kmix
      Ic = cond*(Tsurf-T)/Lsub
  
      END

c********************************************************************************

      subroutine coefffunc(P,T,S,rc,Matm,kmix,Lsub,C0,C1,C2,Ak)

c********************************************************************************
c defini la fonction eq 6 papier 2014  (Listowski et al., 2014)
      use tracer_mod, only: rho_ice_co2
      USE comcstfi_h, ONLY: pi

      implicit none
      include "microphys.h"

c   arguments: INPUT
c   ----------------
      REAL,INTENT(in) :: P
      REAL,INTENT(in) :: T
      REAL*8,INTENT(in) :: S
      DOUBLE PRECISION,INTENT(in) :: rc
      REAL,INTENT(in) :: Matm !g.mol-1 ( = mmean(ig,l) )
c   arguments: OUTPUT
c   ----------
      DOUBLE PRECISION,INTENT(out) ::  C0,C1,C2
      DOUBLE PRECISION,INTENT(out) ::  kmix,Lsub

c   local:
c   ------
      DOUBLE PRECISION Cpatm,Cpn2,Cpco2
      DOUBLE PRECISION psat, xinf, pco2
      DOUBLE PRECISION Dv     
      DOUBLE PRECISION l0,l1,l2,l3,l4            
      DOUBLE PRECISION knudsen, a, lambda      ! F and S correction
      DOUBLE PRECISION Ak                       ! kelvin factor    
      DOUBLE PRECISION vthatm,lpmt,rhoatm, vthco2 ! for Kn,th


c     DEFINE heat cap. J.kg-1.K-1 and To

      data Cpco2/0.7e3/
      data Cpn2/1e3/

      kmix = 0d0
      Lsub = 0d0

      C0 = 0d0
      C1 = 0d0
      C2 = 0d0

c     Equilibirum pressure over a flat surface
      psat = 1.382 * 1.00e12 * exp(-3182.48/dble(T))  ! (Pa)
c     Compute transport coefficient
      pco2 = psat * dble(S)
c     Latent heat of sublimation if CO2  co2 (J.kg-1)
c     version Azreg_Ainou (J/kg) :
      l0=595594.      
      l1=903.111     
      l2=-11.5959    
      l3=0.0528288 
      l4=-0.000103183
      Lsub = l0 + l1 * dble(T) + l2 * dble(T)**2 + l3 * 
     &     dble(T)**3 + l4 * dble(T)**4 ! J/kg
c     atmospheric density
      rhoatm = dble(P*Matm)/(rgp*dble(T))   ! g.m-3
      rhoatm = rhoatm * 1.00e-3 !kg.m-3
      call  KthMixNEW(kmix,T,pco2/dble(P),rhoatm) ! compute thermal cond of mixture co2/N2
      call  Diffcoeff(P, T, Dv)  

      Dv = Dv * 1.00e-4         !!! cm2.s-1  to m2.s-1

c     ----- FS correction for Diff
      vthco2  = sqrt(8d0*kbz*dble(T)/(dble(pi) * mco2/nav)) ! units OK: m.s-1
      knudsen = 3*Dv / (vthco2 * rc) 
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen) ! pas adaptée, Dahneke 1983? en fait si (Monschick&Black)
      Dv      = Dv / (1. + lambda * knudsen)
c     ----- FS correction for Kth 
      vthatm = sqrt(8d0*kbz*dble(T)/(pi * 1.00e-3*dble(Matm)/nav)) ! Matm/nav = mass of "air molecule" in G , *1e-3 --> kg 
      Cpatm = Cpco2 * pco2/dble(P) + Cpn2 * (1d0 - pco2/dble(P)) !J.kg-1.K-1
      lpmt = 3 * kmix / (rhoatm * vthatm * (Cpatm - 0.5*rgp/
     &     (dble(Matm)*1.00e-3))) ! mean free path related to heat transfer
      knudsen = lpmt / rc
      lambda  = (1.333+0.71/knudsen) / (1.+1./knudsen) ! pas adaptée, Dahneke 1983? en fait si (Monschick&Black)
      kmix    = kmix /  (1. + lambda * knudsen)
c     --------------------- ASSIGN coeff values for FUNCTION
      xinf = dble(S) * psat / dble(P)
      Ak = exp(2d0*sigco2*mco2/(rgp* dble(rho_ice_co2*T* rc) ))
      C0 = mco2*Dv*psat*Lsub/(rgp*dble(T)*kmix)*Ak*exp(-Lsub*mco2/
     &     (rgp*dble(T)))
      C1 = Lsub*mco2/(rgp*dble(T)**2)
      C2 = dble(T) + dble(P)*mco2*Dv*Lsub*xinf/(kmix*rgp*dble(T))

      END


c======================================================================
      subroutine Diffcoeff(P, T, Diff)
c     Compute diffusion coefficient CO2/N2
c     cited in Ilona's lecture - from Reid et al. 1987
c======================================================================
       IMPLICIT NONE

       include "microphys.h"

c      arguments
c     -----------
      
      REAL,INTENT(in) :: P
      REAL,INTENT(in) :: T
      
c     output
c     -----------

      DOUBLE PRECISION,INTENT(out) :: Diff
    
c      local
c     -----------

      REAL Pbar                     !!! has to be in bar for the formula
      DOUBLE PRECISION dva, dvb, Mab  ! Mab has to be in g.mol-1
        
        Pbar = P * 1d-5
    
        Mab = 2. / ( 1./mn2 + 1./mco2 ) * 1000.

  	dva = 26.9        ! diffusion volume of CO2,  Reid et al. 1987 (cited in Ilona's lecture)
  	dvb = 18.5        ! diffusion volume of N2
    
        Diff  = 0.00143 * dble(T)**(1.75) / (dble(Pbar) * sqrt(Mab) 
     &       * (dble(dva)**(1./3.) + dble(dvb)**(1./3.))**2.) !!! in cm2.s-1
  
       RETURN

       END


c======================================================================

         subroutine KthMixNEW(Kthmix,T,x,rho)

c        Compute thermal conductivity of CO2/N2 mixture
c         (***WITHOUT*** USE OF VISCOSITY)

c          (Mason & Saxena, 1958 - Wassiljeva 1904)
c======================================================================

       implicit none
       
       include "microphys.h"
c      arguments
c     -----------
         
         REAL,INTENT(in) :: T
         DOUBLE PRECISION,INTENT(in) :: x
         DOUBLE PRECISION,INTENT(in) :: rho !kg.m-3

c     outputs
c     -----------

         DOUBLE PRECISION,INTENT(out) :: Kthmix

c     local
c    ------------

         DOUBLE PRECISION x1,x2

         DOUBLE PRECISION  Tc1, Tc2, Pc1, Pc2

         DOUBLE PRECISION  A12, A11, A22, A21

         DOUBLE PRECISION  Gamma1, Gamma2, M1, M2
         DOUBLE PRECISION  lambda_trans1, lambda_trans2,epsilon

         DOUBLE PRECISION  kco2, kn2

      x1 = x
      x2 = 1d0 - x

      M1 = mco2
      M2 = mn2

      Tc1 =  304.1282 !(Scalabrin et al. 2006)
      Tc2 =  126.192  ! (Lemmon & Jacobsen 2003)

      Pc1 =  73.773   ! (bars)
      Pc2 =  33.958   ! (bars)
    
      Gamma1 = 210.*(Tc1*M1**(3.)/Pc1**(4.))**(1./6.)
      Gamma2 = 210.*(Tc2*M2**(3.)/Pc2**(4.))**(1./6.)

c Translational conductivities

      lambda_trans1 = ( exp(0.0464 * T/Tc1) - exp(-0.2412 * T/Tc1) )
     &                                                          /Gamma1

      lambda_trans2 = ( exp(0.0464 * T/Tc2) - exp(-0.2412 * T/Tc2) )
     &                                                          /Gamma2
      
c     Coefficient of Mason and Saxena
      epsilon = 1.

      A11 = 1.
	 
      A22 = 1.

      A12 = epsilon * (1. + sqrt(lambda_trans1/lambda_trans2)*
     &                    (M1/M2)**(1./4.))**(2.) / sqrt(8*(1.+ M1/M2))

      A21 = epsilon * (1. + sqrt(lambda_trans2/lambda_trans1)*
     &                    (M2/M1)**(1./4.))**(2.) / sqrt(8*(1.+ M2/M1))

c     INDIVIDUAL COND.

         call KthCO2Scalab(kco2,T,rho)
         call KthN2LemJac(kn2,T,rho)

c     MIXTURE COND.
        Kthmix = kco2*x1 /(x1*A11 + x2*A12) + kn2*x2 /(x1*A21 + x2*A22)
        Kthmix = Kthmix*1e-3   ! from mW.m-1.K-1 to  W.m-1.K-1

        END

c======================================================================
         subroutine KthN2LemJac(kthn2,T,rho)
c        Compute thermal cond of N2 (Lemmon and Jacobsen, 2003)
cWITH viscosity
c======================================================================

       implicit none

        include "microphys.h"
c        include "microphysCO2.h"


c      arguments
c     -----------
         
         REAL,INTENT(in) :: T
         DOUBLE PRECISION,INTENT(in) :: rho !kg.m-3

c     outputs
c     -----------

         DOUBLE PRECISION,INTENT(out) :: kthn2

c     local
c    ------------

        DOUBLE PRECISION g1,g2,g3,g4,g5,g6,g7,g8,g9,g10
        DOUBLE PRECISION h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
        DOUBLE PRECISION n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
        DOUBLE PRECISION d4,d5,d6,d7,d8,d9
        DOUBLE PRECISION l4,l5,l6,l7,l8,l9
        DOUBLE PRECISION t2,t3,t4,t5,t6,t7,t8,t9
        DOUBLE PRECISION gamma4,gamma5,gamma6,gamma7,gamma8,gamma9

        DOUBLE PRECISION Tc,rhoc

        DOUBLE PRECISION tau, delta

        DOUBLE PRECISION visco

        DOUBLE PRECISION k1, k2

         N1 = 1.511d0
         N2 = 2.117d0
         N3 = -3.332d0

         N4 = 8.862
         N5 = 31.11
         N6 = -73.13
         N7 = 20.03
         N8 = -0.7096
         N9 = 0.2672

         t2 = -1.0d0
         t3 = -0.7d0
         t4 = 0.0d0
         t5 = 0.03
         t6 = 0.2
         t7 = 0.8
         t8 = 0.6
         t9 = 1.9
   
         d4 =  1.
         d5 =  2.
         d6 =  3.
         d7 =  4.
         d8 =  8.
         d9 = 10.
   
         l4 = 0. 
         gamma4 = 0.
        
         l5 = 0. 
         gamma5 = 0.
        
         l6 = 1. 
         gamma6 = 1.

         l7 = 2. 
         gamma7 = 1.

         l8 = 2. 
         gamma8 = 1.

         l9 = 2. 
         gamma9 = 1.

c----------------------------------------------------------------------           
         call viscoN2(T,visco)  !! v given in microPa.s
      
         Tc   = 126.192d0
         rhoc = 11.1839  * 1000 * mn2   !!!from mol.dm-3 to kg.m-3

         tau  = Tc / T
         delta = rho/rhoc 

         k1 =  N1 * visco + N2 * tau**t2 + N3 * tau**t3  !!! mW m-1 K-1     
c--------- residual thermal conductivity

         k2 = N4 * tau**t4 * delta**d4 * exp(-gamma4*delta**l4)         
     &  +     N5 * tau**t5 * delta**d5 * exp(-gamma5*delta**l5) 
     &  +     N6 * tau**t6 * delta**d6 * exp(-gamma6*delta**l6) 
     &  +     N7 * tau**t7 * delta**d7 * exp(-gamma7*delta**l7) 
     &  +     N8 * tau**t8 * delta**d8 * exp(-gamma8*delta**l8) 
     &  +     N9 * tau**t9 * delta**d9 * exp(-gamma9*delta**l9) 

         kthn2 = k1 + k2

         END


c======================================================================

         subroutine viscoN2(T,visco)

c        Compute viscosity of N2 (Lemmon and Jacobsen, 2003)

c======================================================================

         implicit none

       include "microphys.h"
c       include "microphysCO2.h"
c      arguments
c     -----------
         
      REAL,INTENT(in) :: T

c     outputs
c     -----------

      DOUBLE PRECISION,INTENT(out) :: visco


c     local
c    ------------

      DOUBLE PRECISION a0,a1,a2,a3,a4 
      DOUBLE PRECISION Tstar,factor,sigma,M2
      DOUBLE PRECISION RGCS
      

c----------------------------------------------------------------------  

  
      factor = 98.94   ! (K)   
  
      sigma  = 0.3656  ! (nm)
  
      a0 =  0.431
      a1 = -0.4623
      a2 =  0.08406
      a3 =  0.005341
      a4 = -0.00331
  
      M2 = mn2 * 1.00e3   !!! to g.mol-1
  
      Tstar = T*1./factor

      RGCS = exp( a0 + a1 * log(Tstar) + a2 * (log(Tstar))**2. + 
     &                a3 * (log(Tstar))**3. + a4 * (log(Tstar))**4. )
  
  
      visco = 0.0266958 * sqrt(M2*T) / ( sigma**2. * RGCS )  !!! microPa.s


      RETURN

      END


c======================================================================

         subroutine KthCO2Scalab(kthco2,T,rho)

c        Compute thermal cond of CO2 (Scalabrin et al. 2006)

c======================================================================

         implicit none



c      arguments
c     -----------

      REAL,INTENT(in) :: T
      DOUBLE PRECISION,INTENT(in) :: rho

c     outputs
c     -----------

      DOUBLE PRECISION,INTENT(out) :: kthco2

c     LOCAL
c     -----------

      DOUBLE PRECISION Tc,Pc,rhoc, Lambdac
      
      DOUBLE PRECISION Tr, rhor, k1, k2

      DOUBLE PRECISION g1,g2,g3,g4,g5,g6,g7,g8,g9,g10
      DOUBLE PRECISION h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
      DOUBLE PRECISION n1,n2,n3,n4,n5,n6,n7,n8,n9,n10

      Tc   = 304.1282   !(K) 
      Pc   = 7.3773e6   !(MPa)
      rhoc = 467.6      !(kg.m-3)
      Lambdac = 4.81384 !(mW.m-1K-1)
  
      g1 = 0.
      g2 = 0.
      g3 = 1.5
      g4 = 0.0
      g5 = 1.0
      g6 = 1.5
      g7 = 1.5
      g8 = 1.5
      g9 = 3.5
      g10 = 5.5

      h1 = 1.
      h2 = 5.
      h3 = 1.
      h4 = 1.
      h5 = 2.
      h6 = 0.
      h7 = 5.0
      h8 = 9.0
      h9 = 0.
      h10 = 0.

      n1 = 7.69857587
      n2 = 0.159885811
      n3 = 1.56918621
      n4 = -6.73400790
      n5 = 16.3890156
      n6 = 3.69415242
      n7 = 22.3205514
      n8 = 66.1420950
      n9 = -0.171779133
      n10 = 0.00433043347

      Tr   = T/Tc
      rhor = rho/rhoc

      k1 = n1*Tr**(g1)*rhor**(h1) + n2*Tr**(g2)*rhor**(h2) 
     &     + n3*Tr**(g3)*rhor**(h3)

      k2 = n4*Tr**(g4)*rhor**(h4) + n5*Tr**(g5)*rhor**(h5)  
     &    + n6*Tr**(g6)*rhor**(h6) + n7*Tr**(g7)*rhor**(h7) 
     &    + n8*Tr**(g8)*rhor**(h8) + n9*Tr**(g9)*rhor**(h9) 
     &    + n10*Tr**(g10)*rhor**(h10)
    
      k2  = exp(-5.*rhor**(2.)) * k2
                
      kthco2 = (k1 + k2) *  Lambdac   ! mW

      END
