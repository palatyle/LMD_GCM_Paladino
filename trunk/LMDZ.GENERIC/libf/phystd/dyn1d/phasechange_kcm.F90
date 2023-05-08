!  ----------------------------------------------------------------
!  Purpose: Thermodynamic data on H2O, NH3
!  Authour: Adapted from various sources by R. Wordsworth (2011)


! ---------------------
! NH3
! ---------------------
subroutine psat_NH3 ( T, p )
!
!  PSAT_NH3 makes a rough estimate of the vapor pressure.
!
!  Interpolated from www.engineeringtoolbox.com data by RDW 21/09/11 for
!  temperatures between 223 and 323 K.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision P, the vapor pressure, in MegaPascals.
!

  implicit none

  double precision p,T

  p=exp(-1.5609d-004*T**2 + 0.1236*T - 9.1530)/1d6

end subroutine psat_NH3

subroutine latheat_NH3 ( T, sc, sv )
!
!  PSAT_NH3 makes a rough estimate of the entropies of condensation /
!  vapourisation.
! 
!  Interpolated from www.engineeringtoolbox.com data by RDW 21/09/11 for
!  temperatures between 223 and 323 K.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision sc, the entropy of condensate, in J kg^-1 K^-1.
!            double precision sv, the entropy of gas, in J kg^-1 K^-1.
!

  implicit none

  double precision T,sc,sv

  sv = 0.0492*T**2 - 40.4199*T + 1.2708e+004
  sc = -0.0215*T**2 + 28.7138*T - 5.5267e+003

  !L = -11.2373*T**2 + 2.5326d+03*T + 1.4099d+06

end subroutine latheat_NH3


! ---------------------
! H2O
! ---------------------

subroutine base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )
!
!*******************************************************************************
!
!! BASE calculates quantities associated with the base Helmholtz function.
!
!
!  Discussion:
!
!    The equation for the base Helmholtz function AB(T,RHO) is:
!
!      AB(T,RHO) = R * T * (
!        - ln ( 1 - y ) 
!        - ( beta - 1 ) / ( 1 - y ) 
!        + ( alpha + beta + 1 ) / ( 2 * ( 1 - y )**2 )
!        + 4 * y * ( ( Bbar / b ) - gamma ) 
!        - 0.5 * ( alpha - beta + 3 ) 
!        + ln ( RHO * R * T / P0 ) )
!                                                      (Equation 2)
!   where 
!
!     y = b * rho / 4, 
!     alpha = 11,
!     beta = 133/3, 
!     gamma = 7/2, 
!     P0 = 0.101325 MegaPascals = 1 atm
!
!   and
!
!     b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!     Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!   where 
!
!     T0=647.073 K and the coefficients b(j) and B(j) are
!    
!     j    b(j)                         B(j)
!    --    -----------                  ----------
!     0    0.7478629                    1.1278334
!     1   -0.3540782                   -0.5944001
!     2    0                           -5.010996
!     3    0.007159876                  0
!     4    0                            0.63684256
!     5   -0.003528426                  0
!
!  For the derived quantities, the following relations are used:
!
!    Pressure:                  PB      = RHO**2 * dAB/dRHO
!    Density derivative:        DPDRB   = 2*PB/RHO + RHO**2 * d2AB/dRHO2
!    Temperature derivative:    DPDTB   = RHO**2 * d2AB/(dRHO dT)
!    Specific entropy:          SB      = ( UB - AB ) / T
!    Specific internal energy:  UB      = AB + T * SB
!    Specific enthalpy:         HB      = UB + PB / RHO
!    Specific heat capacity
!      at constant volume:      CVB     = - T * d2AB/dT2
!    Specific Gibbs function:   GB      = AB + PB / RHO
!
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Modified:
!
!    03 February 2002
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Input, double precision RHO, the density, in G/CM3.
!
!    Output, double precision AB, the base value of the Helmholtz function,
!    in KJ/kg.
!
!    Output, double precision CVB, the base value of the isochoric (constant 
!    volume) heat capacity, in KJ/(kg degrees Kelvin).
!
!    Output, double precision DPDRB, the base value of the partial
!    derivative dP(T,RHO)/dRHO, with T held fixed, in (MegaPascals CM3)/G.
!
!    Output, double precision DPDTB, the base value of the partial
!    derivative dP(T,RHO)/dT, with RHO held fixed, in 
!    MegaPascals/degrees Kelvin.
!
!    Output, double precision GB, the base value of the Gibbs free energy,
!    in KJ/kg.
!
!    Output, double precision HB, the base value of enthalpy, in KJ/kg.
!
!    Output, double precision PB, the base pressure, in MegaPascals.
!
!    Output, double precision SB, the base value of entropy, 
!    in KJ/(kg degrees Kelvin).
!
!    Output, double precision UB, the base value of internal energy, 
!    in KJ/kg.
!
  implicit none
!
  double precision ab
  double precision, parameter :: alpha = 11.0D+00
  double precision b1
  double precision b1t
  double precision b1tt
  double precision b2
  double precision b2t
  double precision b2tt
  double precision, parameter :: beta = 44.333333333333D+00
  double precision cvb
  double precision dpdrb
  double precision dpdtb
  double precision dz
  double precision dz0
  double precision, parameter :: gamma = 3.5D+00
  double precision gascon
  double precision gb
  double precision hb
  double precision, parameter :: p_zero = 0.101325D+00
  double precision pb
  double precision rho
  double precision sb
  double precision t
  double precision ub
  double precision x
  double precision y
  double precision z
  double precision z0
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASE - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASE - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if
!
!  Compute auxilliary quantities for Equation 2.
!
  call bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )

  y = 0.25D+00 * b1 * rho

  x = 1.0D+00 - y
!
!  Evaluate Equation 2.
!
  ab =   - log ( 1.0D+00 - y ) &
         - ( beta - 1.0D+00 ) / ( 1.0D+00 - y ) &
         + ( alpha + beta + 1.0D+00 ) / ( 2.0D+00 * ( 1.0D+00 - y )**2 ) &
         + 4.0D+00 * y * ( ( b2 / b1 ) - gamma ) &
         - 0.5D+00 * ( alpha - beta + 3.0D+00 ) &
         + log ( rho * gascon() * t / p_zero )
!
!  Determine quantities defined in terms of AB.
!
  pb = ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**3 &
    + 4.0D+00 * y * ( b2 / b1 - gamma )

  z0 = ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**3

  z = z0 + 4.0D+00 * y * ( b2 / b1 - gamma )

  dz0 = ( alpha + 2.0D+00 * beta * y ) / ( 1.0D+00 - y )**3 &
    + 3.0D+00 * ( 1.0D+00 + alpha * y + beta * y**2 ) / ( 1.0D+00 - y )**4

  dz = dz0 + 4.0D+00 * ( b2 / b1 - gamma )

  gb = ab + pb

  ub = - t * b1t * ( pb - 1.0D+00 - rho * b2 ) / b1 - rho * t * b2t

  hb = pb + ub
!
!  An incorrect version of this equation began:
!
!    cvb = 2.0D+00 * ub + ( pb - 1.0D+00 ) &
!
!  and caused me no end of trouble.  My fault, JVB, 03 February 2002
!
  cvb = 2.0D+00 * ub + ( z0 - 1.0D+00 ) &
    * ( ( t * b1t / b1 )**2 - t**2 * b1tt / b1 ) &
    - rho * t**2 * ( b2tt - gamma * b1tt ) - ( t * b1t / b1 )**2 * y * dz0

  dpdtb = pb / t + rho * ( 0.25D+00 * ( dz0 + 4.0D+00 * ( b2 / b1 - gamma ) ) &
    * b1t + b2t - b2 / b1 * b1t )

  sb = ub - ab

  dpdrb = pb + y * ( dz0 + 4.0D+00 * ( b2 / b1 - gamma ) )
!
!  Assign dimensions.
!
  ab =    gascon() * t       * ab
  cvb =   gascon()           * cvb
  dpdrb = gascon() * t       * dpdrb
  dpdtb = gascon() * t * rho * dpdtb
  gb =    gascon() * t       * gb
  hb =    gascon() * t       * hb
  pb =    gascon() * t * rho * pb
  sb =    gascon()           * sb
  ub =    gascon() * t       * ub

  return
end subroutine base

subroutine bb ( t, b1, b2, b1t, b2t, b1tt, b2tt )
!
!*******************************************************************************
!
!! BB calculates the B's of equations 3 and 4.  
!
!
!  Discussion:
!
!    Here
!
!      b(T) = b1 * ln(T/T0) + sum(j=0,1,3,5) b(j)*(T0/T)**j  (Equation 3)
!
!      Bbar(T) = sum(j=0,1,2,4) B(j)*(T0/T)**j               (Equation 4).
!
!    where 
!
!      T0 = 647.073 K 
!
!    and the coefficients b(j) and B(j) are
!    
!      j    b(j)                         B(j)
!     --    -----------                  ----------
!      0    0.7478629                    1.1278334
!      1   -0.3540782                   -0.5944001
!      2    0                           -5.010996
!      3    0.007159876                  0
!      4    0                            0.63684256
!      5   -0.003528426                  0
!
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision B1, the coefficient b from equation 3, 
!    in CM3/G.
!
!    Output, double precision B2, the coefficient Bbar from equation 4, 
!    in CM3/G.
!
!    Output, double precision B1T, the derivative dB1/dT, 
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, double precision B2T, the derivative dB2/dT, 
!    in (CM3)/(G Degrees Kelvin).
!
!    Output, double precision B1TT, the second derivative of B1 with 
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
!    Output, double precision B2TT, the second derivative of B2 with 
!    respect to T, in (CM3)/(G (Degrees Kelvin)**2 ).
!
  implicit none
!
  double precision b1
  double precision b1t
  double precision b1tt
  double precision b2
  double precision b2t
  double precision b2tt
  double precision, parameter, dimension ( 10 ) :: bp = (/ &
    0.7478629D+00,   -0.3540782D+00,    0.0D+00,           0.0D+00, &
    0.007159876D+00,  0.0D+00,         -0.003528426D+00,   0.0D+00, &
    0.0D+00,          0.0D+00 /)
  double precision, parameter, dimension ( 10 ) :: bq = (/ &
    1.1278334D+00,    0.0D+00,         -0.5944001D+00,   -5.010996D+00, &
    0.0D+00,          0.63684256D+00,   0.0D+00,          0.0D+00, &
    0.0D+00,          0.0D+00 /)
  integer i
  double precision t
  double precision, parameter :: t_ref = 647.073D+00
  double precision v(10)
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BB - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Set V(I) = ( T_REF / T )**(I-1).
!
  v(1) = 1.0D+00
  do i = 2, 10
    v(i) = v(i-1) * t_ref / t
  end do
!
!  Set B1, B1T, B1TT.
!
  b1 = bp(1) + bp(2) * log ( 1.0D+00 / v(2) )
  b1t = bp(2) * v(2) / t_ref
  b1tt = 0.0D+00
  do i = 3, 10
    b1 = b1 + bp(i) * v(i-1)
    b1t = b1t - dble ( i - 2 ) * bp(i) * v(i-1) / t
    b1tt = b1tt + bp(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b1tt = b1tt -  ( b1t / t )
!
!  Set B2, B2T, B2TT.
!
  b2 = bq(1)
  b2t = 0.0D+00
  b2tt = 0.0D+00
  do i = 3, 10
    b2 = b2 + bq(i) * v(i-1)
    b2t = b2t - dble ( i - 2 ) * bq(i) * v(i-1) / t
    b2tt = b2tt + bq(i) * dble ( i - 2 )**2 * v(i-1) / t**2
  end do

  b2tt = b2tt - ( b2t / t )

  return
end subroutine bb


subroutine ideal ( t, ai, cpi, cvi, gi, hi, si, ui )
!
!*******************************************************************************
!
!! IDEAL computes ideal gas thermodynamic properties of water.
!
!
!  Discussion:
!
!    Values for thermodynamic properties of water in the ideal
!    gas state were reported by Woolley.  The formula for the ideal gas
!    term of the Helmholtz function approximates a term by term summation of 
!    contributions from each of the rotation and vibration states.  
!    The formula, equation #6 in the reference, is:
!   
!    A(ideal)(T) = -R * T * ( 1 + ( C(1)/Tr + C(2) ) * ln(Tr)
!      + Sum ( 3 <= I <= 18) C(I) * Tr**(I-6)
!   
!    where Tr=T/100 K.  The C(i) are tabulated coefficients.  Equation
!    6 can be used for temperatures below 3000 K, and is accurate to
!    within the tolerance of the gas constant for 50<=T<=2000 K.
!     
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision AI, the Helmholtz function, in KJ/kg.
!
!    Output, double precision CPI, the heat capacity at constant pressure,
!    in KJ/(kg degrees Kelvin).
!
!    Output, double precision CVI, the heat capacity at constant volume,
!    in KJ/(kg degrees Kelvin).
!
!    Output, double precision GI, the Gibbs free energy, in KJ/kg.
!
!    Output, double precision HI, the enthalpy, in KJ/kg.
!
!    Output, double precision SI, the entropy, in KJ/(kg degrees Kelvin).
!
!    Output, double precision UI, the internal energy, in KJ/kg.
!
  implicit none
!
  double precision ai
  double precision, parameter, dimension ( 18 ) :: c = (/ &
   19.730271018D+00,      20.9662681977D+00,     -0.483429455355D+00, &
    6.05743189245D+00,    22.56023885D+00,       -9.87532442D+00, &
   -4.3135538513D+00,      0.458155781D+00,      -0.047754901883D+00, &
    0.0041238460633D+00,  -0.00027929052852D+00,  0.14481695261D-04, &
   -0.56473658748D-06,     0.16200446D-07,       -0.3303822796D-09, &
    0.451916067368D-11,   -0.370734122708D-13,    0.137546068238D-15 /)
  double precision cpi
  double precision cvi
  double precision gascon
  double precision gi
  double precision hi
  integer i
  double precision si
  double precision t
  double precision temp
  double precision tt
  double precision ui
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDEAL - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  tt = t / 100.0D+00

  gi = - ( c(1) / tt + c(2) ) * log ( tt )
  do i = 3, 18
    gi = gi - c(i) * tt**(i-6)
  end do

  hi = c(2) + c(1) * ( 1.0D+00 - log ( tt ) ) / tt
  do i = 3, 18
    hi = hi + dble ( i - 6 ) * c(i) * tt**(i-6)
  end do

  cpi = c(2) - c(1) / tt
  do i = 3, 18
    cpi = cpi + dble ( ( i - 6 ) * ( i - 5 ) ) * c(i) * tt**(i-6)
  end do

  ai = gi - 1.0D+00
  ui = hi - 1.0D+00
  cvi = cpi - 1.0D+00
  si = hi - gi
!
!  Assign dimensions.
!
  ai =  gascon() * t * ai
  cpi = gascon()     * cpi
  cvi = gascon()     * cvi
  gi =  gascon() * t * gi
  hi =  gascon() * t * hi
  si =  gascon()     * si
  ui =  gascon() * t * ui

  return
end subroutine ideal


subroutine resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )
!
!*******************************************************************************
!
!! RESID calculates residual contributions to thermodynamic quantities.
!
!
!  Discussion:
!
!    The residual function consists of 40 terms.  The first 36 are
!    used in a global least squares fit to experimental data.
!
!    Three terms were added that contribute only in the immediate
!    neighborhood of the critical point 
!      (tk-5) <= T <= (tk+5) C
!      0.20   <= rho <= 0.44 g/cm3, 
!
!    A single term was added for the region of high pressure and 
!    low temperature: T < 75 C, P > 300 MPa.
!
!    Except in these limited regions, the residual function is
!    given by the first 36 terms.  The equation is
!   
!      A(residual)(rho,T)=
!        sum(i=1 to 36) (g(i)/k(i)) * (T0/T)**(l(i)) (1-exp(-rho))**(k(i))
!      + sum(i=37 to 40) g(i)*delta(i)**(k(i))
!        * exp(-alpha(i)*delta(i)**(k(i)) - beta(i)*tau(i)**2)
!                                                     (Equation 5)
!   
!    where
!
!      g(i) are coefficients determined by fits to data,
!      delta(i) are reduced densities (delta(i)=((rho-rho(i))/rho(i))
!      tau(i) are reduced temperatures (tau(i)=((T-tau(i))/tau(i))
!      rho(i) are specified densities.
!      tau(i) are specified temperatures.
!      The k(i) and l(i) are specified integers.
!   
!  Modified:
!
!    22 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Input, double precision RHO, the density, in G/CM3.
!
!    Output, double precision AR, the residual contribution to the 
!    Helmholtz function, in KJ/kg.
!
!    Output, double precision CVR, the residual contribution to the 
!    isochoric (constant volume) heat capacity, in KJ/(kg degrees Kelvin).
!
!    Output, double precision DPDRR, the residual contribution to 
!    the partial derivative dP(T,RHO)/dRHO, with T held fixed, in 
!    (MegaPascals CM3)/G.
!
!    Output, double precision DPDTR, the residual contribution to 
!    the partial derivative dP(T,RHO)/dT, with RHO held fixed, 
!    in MegaPascals/degrees Kelvin.
!
!    Output, double precision GR, the residual contribution to the Gibbs 
!    function, in KJ/kg.
!
!    Output, double precision HR, the residual contribution to the 
!    enthalpy, in KJ/kg.
!
!    Output, double precision PR, the residual contribution to the pressure, 
!    in MegaPascals.
!
!    Output, double precision SR, the residual contribution to the entropy, 
!    in KJ/(kg degrees Kelvin).
!
!    Output, double precision UR, the residual contribution to the 
!    internal energy, in KJ/kg.
!
  implicit none
!
  double precision, parameter, dimension ( 4 ) :: aad = (/ &
    34.0D+00, 40.0D+00, 30.0D+00, 1050.0D+00 /)
  double precision, parameter, dimension ( 4 ) :: aat = (/ &
    20000.0D+00, 20000.0D+00, 40000.0D+00, 25.0D+00 /)
  double precision, parameter, dimension ( 4 ) :: adz = (/ &
    0.319D+00, 0.319D+00, 0.319D+00, 1.55D+00 /)
  double precision ar
  double precision att
  double precision, parameter, dimension ( 4 ) ::  atz = (/ &
    640.0D+00, 640.0D+00, 641.6D+00, 270.0D+00 /)
  double precision cvr
  double precision dadt
  double precision ddz
  double precision del
  double precision dex
  double precision dfdt
  double precision dpdrr
  double precision dpdtr
  double precision e
  double precision errtol
  double precision ex0
  double precision ex1
  double precision ex2
  double precision fct
  double precision, parameter, dimension ( 40 ) :: g = (/ &
    -530.62968529023D+00,  0.22744901424408D+04, 0.78779333020687D+03, &
    -69.830527374994D+00,  0.17863832875422D+05,-0.39514731563338D+05, &
    0.33803884280753D+05, -0.13855050202703D+05,-0.25637436613260D+06, &
    0.48212575981415D+06, -0.34183016969660D+06, 0.12223156417448D+06, &
    0.11797433655832D+07, -0.21734810110373D+07, 0.10829952168620D+07, &
   -0.25441998064049D+06, -0.31377774947767D+07, 0.52911910757704D+07, &
   -0.13802577177877D+07, -0.25109914369001D+06, 0.46561826115608D+07, &
   -0.72752773275387D+07,  0.41774246148294D+06, 0.14016358244614D+07, &
   -0.31555231392127D+07,  0.47929666384584D+07, 0.40912664781209D+06, &
   -0.13626369388386D+07,  0.69625220862664D+06,-0.10834900096447D+07, &
   -0.22722827401688D+06,  0.38365486000660D+06, 0.68833257944332D+04, &
    0.21757245522644D+05, -0.26627944829770D+04,-0.70730418082074D+05, &
   -0.225D+00, -1.68D+00, 0.055D+00, -93.0D+00 /)
  double precision gascon
  double precision gr
  double precision hr
  integer i
  integer, parameter, dimension ( 40 ) :: ii = (/ &
    0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6, &
    8,8,8,8,2,2,0,4,2,2,2,4 /)
  integer j
  integer, parameter, dimension ( 40 ) :: jj = (/ &
    2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,2,3,5,7,&
    2,3,5,7,1,4,4,4,0,2,0,0 /)
  integer k
  integer l
  integer nc
  double precision pr
  double precision q10
  double precision q20
  double precision q2a
  double precision q5t
  double precision qm
  double precision qp
  double precision qr(11)
  double precision qt(10)
  double precision rho
  double precision sr
  double precision, parameter :: s_ref = 7.6180720166752D+00
  double precision t
  double precision, parameter :: t_ref = 647.073D+00
  double precision tau
  double precision tx
  double precision, parameter :: u_ref = - 4328.4549774261D+00
  double precision ur
  double precision v
!
  errtol = sqrt ( epsilon ( errtol ) )
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESID - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESID - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if

  nc = 36
  dpdrr = 0.0D+00
  pr = 0.0D+00
  ar = 0.0D+00
  dadt = 0.0D+00
  cvr = 0.0D+00
  dpdtr = 0.0D+00

  ex0 = - rho
! ex0 = max ( ex0, - 225.0D+00 )
! ex0 = min ( ex0, 225.0D+00 )
  e = exp ( ex0 )

  q10 = rho * rho * e
  q20 = 1.0D+00 - e

  qr(1) = 0.0D+00
  qr(2) = q10
  do i = 2, 10
    qr(i+1) = qr(i) * q20
  end do

  v = t_ref / t
  qt(1) = t / t_ref
  do i = 2, 10
    qt(i) = qt(i-1) * v
  end do

  do i = 1, nc

    k = ii(i) + 1
    l = jj(i)
    qp = g(i) * qr(k+1) * qt(l+1)
    pr = pr + qp

    dpdrr = dpdrr + ( 2.0D+00 / rho - ( 1.0D+00 - e * dble ( k - 1 ) / &
      ( 1.0D+00 - e ) ) ) * qp

    ar = ar + g(i) * qr(k+2) * qt(l+1) / ( rho**2 * e * dble ( k ) &
      * gascon ( ) * t )

    dfdt = ( 1.0D+00 - e )**k * dble ( 1 - l ) * qt(l+2) / t_ref / dble ( k )

    dadt = dadt + g(i) * dfdt

    dpdtr = dpdtr + g(i) * dfdt * rho**2 * e * dble ( k ) / ( 1.0D+00 - e )

    cvr = cvr + g(i) * dble ( l ) * dfdt / gascon()

  end do

  qp = 0.0D+00
  q2a = 0.0D+00

  do j = 37, 40

    k = ii(j)
    ddz = adz(j-36)
    del = rho / ddz - 1.0D+00

    if ( abs ( del ) < errtol ) then
      del = errtol
    end if

    ex1 = - aad(j-36) * del**k
!   ex1 = max ( ex1, - 225.0D+00 )
!   ex1 = min ( ex1, 225.0D+00 )
    dex = exp ( ex1 ) * del**jj(j)

    att = aat(j-36)
    tx = atz(j-36)
    tau = ( t / tx ) - 1.0D+00

    ex2 = - att * tau**2
!   ex2 = max ( ex2, - 225.0D+00 )
!   ex2 = min ( ex2, 225.0D+00 )
    q10 = dex * exp ( ex2 )

    qm = dble ( jj(j) ) / del - dble ( k ) * aad(j-36) * del**(k-1)
    fct = qm * rho**2 * q10 / ddz

    q5t = fct * ( 2.0D+00 / rho + qm / ddz ) - ( rho / ddz )**2 * q10 * &
      ( dble ( jj(j) ) / del**2 + dble ( k * ( k - 1 ) ) * aad(j-36) * &
      del**(k-2) )

    dpdrr = dpdrr + q5t * g(j)
    qp = qp + g(j) * fct
    dadt = dadt - 2.0D+00 * g(j) * att * tau * q10 / tx
    dpdtr = dpdtr - 2.0D+00 * g(j) * att * tau * fct / tx

    q2a = q2a + t * g(j) * att * ( 4.0D+00 * ex2 + 2.0D+00 ) * q10 / tx**2

    ar = ar + q10 * g(j) / ( gascon() * t )

  end do

  cvr = cvr + q2a / gascon()
  pr = pr + qp
  sr = - dadt / gascon()
  ur = ar + sr
!
!  Assign dimensions.
!
  ar =  gascon() * t *  ar
  cvr = gascon() *     cvr
  sr =  gascon() *      sr
  ur =  gascon() * t *  ur
!
!  Adjust energies.
!
  ar = ar + gascon ( ) * t * s_ref - gascon ( ) * u_ref
  sr = sr - gascon ( ) * s_ref
  ur = ur - gascon ( ) * u_ref

  gr = ar + pr / rho
  hr = ur + pr / rho

  return
end subroutine resid

subroutine psat_H2O ( t, p )
!
!*******************************************************************************
!
!! PSAT_H2O makes a rough estimate of the vapor pressure.
!
!
!  Discussion:
!
!    The calculation agrees with tabulated data to within
!    0.02% for temperature to within a degree or so of the critical
!    temperature.  The approximate vapor pressure can be refined
!    by imposing the condition that the Gibbs functions of the vapor
!    and liquid phases be equal.
!
!  Modified:
!
!    21 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision P, the vapor pressure, in MegaPascals.
!
  implicit none
!
  double precision, parameter, dimension ( 8 ) :: a = (/ &
    -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,  -148.39348D+00, 48.631602D+00 /)
  double precision b
  integer i
  double precision p
  double precision q
  double precision t
  double precision, parameter :: t_ref = 647.25D+00
  double precision v
  double precision w
  double precision z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PSAT_H2O - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  if ( t <= 314.0D+00 ) then

    p = 0.1D+00 * exp ( 6.3573118D+00 - 8858.843D+00 / t &
      + 607.56335D+00 * t**( -0.6D+00 ) )

  else

    v = t / t_ref
    w = abs ( 1.0D+00 - v )
    b = 0.0D+00
    do i = 1, 8
      z = i
      b = b + a(i) * w**( ( z + 1.0D+00 ) / 2.0D+00 )
    end do

    q = b / v
    p = 22.093D+00 * exp ( q )

  end if

  return
end subroutine psat_H2O


subroutine tdpsdt ( t, dp )
!
!*******************************************************************************
!
!! TDPSDT computes the quantity T * dP(Sat)/dT.
!
!
!  Discussion:
!
!    Here T is the temperature and P(Sat) is the vapor pressure.  
!    It is used by TSAT_EST and TSAT.
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Output, double precision DP, the value T*(dP(Sat)/dT), 
!    in MegaPascals.
!
  implicit none
!
  double precision, parameter, dimension ( 8 ) :: a = (/ &
      -7.8889166D+00,   2.5514255D+00,   -6.716169D+00, 33.239495D+00, &
    -105.38479D+00,   174.35319D+00,   -148.39348D+00,  48.631602D+00 /)
  double precision b
  double precision c
  double precision dp
  integer i
  double precision q
  double precision t
  double precision, parameter :: t_ref = 647.25D+00
  double precision v
  double precision w
  double precision y
  double precision z
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TDPSDT - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if

  v = t / t_ref
  w = 1.0D+00 - v
  b = 0.0D+00
  c = 0.0D+00
  do i = 1, 8
    z = dble ( i + 1 ) / 2.0D+00
    y = a(i) * w**z
    c = c + ( y / w ) * ( 0.5D+00 - 0.5D+00 * dble ( i ) - 1.0D+00 / v )
    b = b + y
  end do

  q = b / v
  dp = 22.093D+00 * exp ( q ) * c

  return
end subroutine tdpsdt



subroutine therm ( t, rho, a, cjth, cjtt, cp, cv, dpdr, dpdt, g, h, p, s, u )
!
!*******************************************************************************
!
!! THERM calculates thermodynamic functions given temperature and density.
!
!
!  Discussion:
!
!    Thermodynamic values were calculated from an analytic equation
!    that approximates the Helmholtz function (specific Helmholtz
!    energy) for ordinary water and steam, of the form A=A(rho,T)
!    where A is the Helmholtz function, rho the density, and T
!    the absolute (thermodynamic) temperature.  Any thermodynamic
!    value for any state, liquid, vapor or metastable, may be
!    calculated by differentiation of this equation in accord with
!    the first and second laws of thermodynamics.
! 
!    The International Association for the Properties of Steam
!    has provisionally accepted this formulation for the range
!    273.15 <= T <= 1273.15 degrees Kelvin, where, for 423.15 <= T,
!    the maximum pressure is Pmax = 1500 MPa = 15000 bar, and for
!    273.15 <= T < 423.15, the maximum pressure is
!    Pmax = 100 * (5 + (T-273.15)/15) MPa.
! 
!    Close to the critical point, a small region is excluded:
!    Abs(T-Tk) < 1, abs((rho-rhok)/rhok) < 0.3.
! 
!    The equation has a wider useful range, namely, fluid states
!    of pure, undissociated water and steam defined by
!    260 <= T <= 2500 K and 0 <= P <= 3000 MPa.
! 
!    Thermodynamic property values for specific volume, density,
!    specific internal energy, specific enthalpy, and specific
!    entropy of water and steam were tabulated over the range
!    0 <= t <= 2000 C, 0 <= P <= 3000 MPa.  The reference
!    state is the liquid at the triple point, for which the
!    internal energy and entropy have been assigned the value zero.
!     
!    Thermodynamic quantities are determined from the Helmholtz function
!    A(rho,T), which is computed as the sum of three terms:
!
!      A(rho,T) = A(base)(rho,T) + A(residual)(rho,T) + A(ideal)(T)
!                                                       (Equation 1)
!
!    Because A(rho,T) is everywhere single valued and analytic,
!    we can derive closed form relations for all other properties.
!    In the following, unless otherwise indicated, the independent
!    variables are temperature T and density RHO, and differentiation
!    with respect to one variable is to imply that the other is fixed.
! 
!    Pressure:                  P       = RHO**2 * dA/dRHO
!    Density derivative:        dP/dRHO = 2*P/RHO + RHO**2 * d2A/dRHO2
!    Temperature derivative:    dP/dT   = RHO**2 * d2A/(dRHO dT)
!    Specific entropy:          S       = - dA/dT
!    Specific internal energy:  U       = A + T*S
!    Specific enthalpy:         H       = U + P/RHO
!    Specific heat capacity
!      at constant volume:      Cv      = - T * d2A/dT2
!    Specific Gibbs function:   G       = A + P/RHO
!    Specific heat capacity
!      at constant pressure:    Cp      = Cv + (T*(dP/dT)**2)/(RHO**2*dP/dRHO)
!    Speed of sound:            Omega   = Sqrt ((Cp/Cv) * dP/dRHO)
!    Second virial coefficient: B       = 1/(2*R*T) * (d2P/dRHO2) (at RHO=0)
!    Isothermal Joule-Thomson
!      coefficient:             DeltaT  = (dH/dP) (fixed T) =
!                                         (1/RHO)-(T*dP/dT)/(RHO**2*dP/dRHO)
!    Joule-Thomson coefficient: Mu      = (dT/dP) (fixed H) = DeltaT/Cp
!    Isentropic temperature-
!      pressure coefficient:    BetaS   = (dT/dP) (fixed S) =
!                                         (DeltaT - 1/RHO)/Cp
!   
!  Modified:
!
!    19 November 1998
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Input, double precision T, the temperature, in degrees Kelvin.
!
!    Input, double precision RHO, the fluid density, in G/CM3.
!
!    Output, double precision A, the Helmholtz function, in KJ/kg.
!
!    Output, double precision CJTH, the Joule-Thomson coefficient,
!    in K/MegaPascals.
!
!    Output, double precision CJTT, the isothermal Joule-Thomson coefficient,
!    in CM3/G.
!
!    Output, double precision CP, the isobaric (constant pressure) heat
!    capacity, in KJ/(kg degrees Kelvin).
!
!    Output, double precision CV, the isochoric (constant volume) heat capacity,
!    in KJ/(kg degrees Kelvin).
!
!    Output, double precision DPDR, the partial derivative 
!    dP(T,RHO)/dRHO, with T held fixed, in MegaPascals*CM3/G.
!
!    Output, double precision DPDT, the partial derivative 
!    dP(T,RHO)/dT, with RHO held fixed, in MegaPascals/degrees Kelvin.
!
!    Output, double precision G, the Gibbs free energy, in KJ/kg.
!
!    Output, double precision H, the enthalpy, in KJ/kg.
!
!    Output, double precision P, the pressure, in MegaPascals.
!
!    Output, double precision S, the entropy, in KJ/(kg degrees Kelvin).
!
!    Output, double precision U, the internal energy, in KJ/kg.
!
  implicit none
!
  double precision a
  double precision ab
  double precision ai
  double precision ar
  double precision cjth
  double precision cjtt
  double precision cp
  double precision cpi
  double precision cv
  double precision cvb
  double precision cvi
  double precision cvr
  logical, parameter :: debug = .false.
!  logical, parameter :: debug = .true.
  double precision dpdr
  double precision dpdrb
  double precision dpdrr
  double precision dpdt
  double precision dpdtb
  double precision dpdtr
  double precision g
  double precision gb
  double precision gi
  double precision gr
  double precision h
  double precision hb
  double precision hi
  double precision hr
  double precision p
  double precision pb
  double precision pr
  double precision rho
  double precision s
  double precision sb
  double precision si
  double precision sr
  double precision t
  double precision u
  double precision ub
  double precision ui
  double precision ur
!
!  Refuse to handle zero or negative temperatures.
!
  if ( t <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERM - Fatal error!'
    write ( *, '(a)' ) '  The input temperature T must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was T = ', t
    stop
  end if
!
!  Refuse to handle zero or negative density.
!
  if ( rho <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'THERM - Fatal error!'
    write ( *, '(a)' ) '  The input density RHO must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value was RHO = ', rho
    stop
  end if

  call ideal ( t, ai, cpi, cvi, gi, hi, si, ui )

  call resid ( t, rho, ar, cvr, dpdrr, dpdtr, gr, hr, pr, sr, ur )

  call base ( t, rho, ab, cvb, dpdrb, dpdtb, gb, hb, pb, sb, ub )

  a =       ab +    ar +  ai
  cv =     cvb +   cvr + cvi

  if ( debug ) then
    write ( *, * ) ' '
    write ( *, * ) 'THERM:'
    write ( *, * ) '  CVB = ', cvb
    write ( *, * ) '  CVR = ', cvr
    write ( *, * ) '  CVI = ', cvi
    write ( *, * ) '  CV  = ', cv
  end if


  dpdr = dpdrb + dpdrr
  dpdt = dpdtb + dpdtr
  p =       pb +    pr
  s =       sb +    sr +  si
  u =       ub +    ur +  ui

  if ( debug ) then
    write ( *, * ) ' '
    write ( *, * ) 'THERM:'
    write ( *, * ) '  UB = ', ub
    write ( *, * ) '  UR = ', ur
    write ( *, * ) '  UI = ', ui
  end if

  g = a + p / rho
  h = u + p / rho
  cp = cv + t * dpdt**2 / ( dpdr * rho**2 )
  cjtt = 1.0D+00 / rho - t * dpdt / ( dpdr * rho**2 )
  cjth = - cjtt / cp

  return
end subroutine therm

function gascon ( )
!
!*******************************************************************************
!
!! GASCON returns the value of the specific gas constant.
!
!
!  Note:
!
!    The specific gas constant R is related to the universal gas
!    constant R-bar = 8.31441 J/(mol degrees Kelvin) by the molar mass 
!    M = 18.0152 g/mol:
!
!      R = R-bar / M.
!
!  Reference:
!
!    Lester Haar, John Gallagher and George Kell,
!    NBS/NRC Steam Tables:
!    Thermodynamic and Transport Properties and Computer Programs
!      for Vapor and Liquid States of Water in SI Units,
!    Hemisphere Publishing Corporation, Washington, 1984,
!    TJ270.H3
!
!    C A Meyer, R B McClintock, G J Silvestri, R C Spencer,
!    ASME Steam Tables: Thermodynamic and Transport Properties of Steam,
!    American Society of Mechanical Engineers, 1967.
!
!  Parameters:
!
!    Output, double precision GASCON, the value of the specific gas 
!    constant, in J/(g degrees Kelvin).
!
  implicit none
!
  double precision gascon
!
  gascon = 0.461522D+00
 
  return
end function gascon
