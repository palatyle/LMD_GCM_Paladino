SUBROUTINE FLUX(TAIR,PAIR,dt,sig,r,M3g,dM0,dM3)

  !*  Masses flux = Condensation/Evaporation + Thermodynamical equilibrium
  !* 
  !*  Condensation/Evaporation:
  !*    Sulfuric acid drive this process (James 1997): most present in droplets
  !*    more or less 95%
  !*
  !*  Thermodynamical equilibrium:
  !*    Only water flux here (James 1997) => compare to acid vapor, water 
  !*    vapor is most present in the atmosphere
  !*
  !*  ONLY FOR ONE MODE HERE

  use free_param
  use donnees

  IMPLICIT NONE

  real, intent(in), dimension(3) :: M3g ! Third moment of the mode
  real, intent(in) :: TAIR, PAIR, dt    ! Temp, timestep, pressure
  real, intent(in) :: r, sig      ! Mean radius and variance of the mode
  real, intent(out) :: dM3, dM0   ! Tendancy

  real :: RDSA, RCSA               ! Resistance
  real :: A , B, cste, a1, a2, a3  ! Calculus cstes
  real :: alpha_k   ! Function
  real :: MSAD      ! Mass of sulfuric acid in the droplet, in kg
  real :: mk3       ! Tendancy
  real :: gamma


  ! ----- EQUILIBRIUM -----
  CALL WSA_ROSA_NEW(TAIR,PAIR,r,WSAEQ,MSAD) ! Calculation of WSA

  ! ----- CONDENSATION / EVAPORATION -----
  IF (WSAEQ .gt. 0) THEN
     ! Resistance due to the VAPOR diffusion (s/m2)
     ! Here, we supposed a Dirac function for the calculation of D (Kn(r)
     RDSA = (RHOSA*RGAS*TAIR) / (D*MSA*RHOsasat)
     ! Resistance due to the HEAT diffusion (s/m2)
     RCSA = (LSA*RHOSA)/(KAIR*TAIR) * ((LSA*MSA)/(RGAS*TAIR)-1.0D0)

     A = 2.0D0*ST*MSA / (RHOSA*RGAS*TAIR) !m
     B = exp(A/r)  
     cste = 3.0D0/(RCSA+RDSA)

     a1 = SH2SO4-B-A*B/r-(r**2)*B*A/2.0D0*(2.0D0*r+A)/(r**4)
     a2 = A*B/(r**3) * (A+3.0D0*r)                  !m-1
     a3 = (-1.0D0)*A*B * (2.0D0*r+A)/(2.0D0*r**4)   !m-2

     gamma = (a1 * r**(-2) * alpha_k(1,sig)/alpha_k(3,sig) + &
          &  a2 * r**(-1) * alpha_k(2,sig)/alpha_k(3,sig) + &
          &  a3) * cste

     mk3 = (1.D0/dt)*((WSA/WSAEQ) - 1.D0)*dt
     mk3 = mk3 + (gamma/WSAEQ)
     mk3 = 1.D0 - mk3
     mk3 = (1.D0/mk3) * (M3g(1)+M3g(2)+M3g(3))

     ! ----- TOTAL FLUX -----
     !  dm < 0: evaporation and dm > 0: condensation
     dM3 = mk3 - (M3g(1) + M3g(2) + M3g(3)) !m3 s
     dM0 = dM3 / (r**3*alpha_k(3,sig))
  ELSE
     dM3 = 0.D0
     dM0 = 0.D0
  END IF

  WSA = WSAEQ

  RETURN 
END SUBROUTINE FLUX
