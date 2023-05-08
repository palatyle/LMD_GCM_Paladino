!     Saturated sulfuric acid vapor pressure 
!     over plane sulfuric acid solution.
!
!     Source: Zeleznik
!
!     Input:  TAIR: Temperature  (K)
!             WSA:  Weight fraction of H2SO4  [0;1]
!     Output: H2SO4 vapor pressure 
!             over sulfuric acid solution    (Pa) )
!
!     External functions needed for calculation of activity coefficients
!****************************************************************

FUNCTION PSASAS_ZELE(TAIR,WSA)                                        
  
  IMPLICIT NONE
  
  REAL, INTENT(in) :: TAIR,WSA
  REAL :: ADOT,BDOT,CDOT
  REAL :: MMHGPA, PSASAS_ZELE
  REAL :: pstand,x1,lpar,acidps,act(2)
  
  pstand=1.01325D5 !Pa  1 atm pressure

  !     Mole fraction needed for zeleznik
  x1=(wsa/98.08D0)/(wsa/98.08D0 + ((1.0-wsa)/18.106D0))

  call zeleznik(x1,tair,act)
  
  !     Pure acid satur vapor pressure
  lpar= -11.695D0+LOG(pstand) ! Zeleznik

  acidps=1.0/360.15D0-1.0/tair+0.38D0/545.0D0 &
       & *(1.0+LOG(360.15D0/tair)-360.15D0/tair)
  
  acidps = 10156.0D0*acidps +lpar
  acidps = EXP(acidps)    !Pa

  !     Sat.vap.pres over mixture (flat surface) in Pa:
  psasas_zele=act(2)*acidps



  RETURN 
END FUNCTION PSASAS_ZELE

!****************************************************************
SUBROUTINE Zeleznik(x,T,act)
    
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !     Water and sulfuric acid activities in liquid
  !     aqueous solutions.
  !     Frank J. Zeleznik, Thermodynnamic properties
  !     of the aqueous sulfuric acid system to 220K-350K,
  !     mole fraction 0,...,1
  !     J. Phys. Chem. Ref. Data, Vol. 20, No. 6,pp.1157, 1991
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
  
  IMPLICIT NONE
  
  REAL, INTENT(inout) :: x
  REAL, INTENT(in) :: T
  REAL, INTENT(out), DIMENSION(2):: act
  REAL, external :: activitya, activityw

  if (x .gt. 0.99D0) then
     x=0.99D0
  endif

  act(2)=activitya(x,T) 
  act(1)=activityw(x,T)
  
END SUBROUTINE Zeleznik

!****************************************************************
FUNCTION activitya(xal,T)

  IMPLICIT NONE
  
  REAL, INTENT(in) :: T,xal
  REAL :: activitya 
  REAL, external :: lnAa

  activitya = EXP(lnAa(xal,T)-lnAa(0.999999,T))

END FUNCTION activitya
  
!****************************************************************
FUNCTION activityw(xal,T)
  
  IMPLICIT NONE
  
  REAL, intent(in) :: T,xal 
  REAL :: activityw 
  REAL, external :: lnAw
  
  activityw=EXP(lnAw(xal,T)-lnAw(0.000001,T))
  
END FUNCTION activityw

!****************************************************************
FUNCTION lnAa(x1,T)
  
  IMPLICIT NONE

  REAL :: lnAa
  REAL, INTENT(in) :: T,x1
  REAL, external :: m111,m121,m221,m122 
  REAL, external :: e111,e121,e211,e122,e212,e221

  lnAa=-( &
       &    (2*m111(T)+e111(T)*(2*log(x1)+1))*x1 &
       &    +(2*m121(T)+e211(T)*log(1-x1)+e121(T)*(log(x1)+1))*(1-x1) &
       &    -(m111(T)+e111(T)*(log(x1)+1))*x1*x1 &
       &    -(2*m121(T)+e121(T)*(log(x1)+1)+e211(T)*(log(1-x1)+1) &
       &    -(2*m122(T)+e122(T)*log(x1) &
       &           +e212(T)*log(1-x1))*(1-x1))*x1*(1-x1) &
       &    -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2 &
       &    -x1*(1-x1)*( &
       &                  (6*m122(T)+e122(T)*(3*log(x1)+1) &
       &                          +e212(T)*(3*log(1-x1)+1) &
       &                   )*x1*(1-x1) &
       &                -(2*m122(T)+e122(T)*(log(x1)+1) &
       &                                   +e212(T)*log(1-x1) &
       &                    )*(1-x1)) &
       &     )

END FUNCTION lnAa

!****************************************************************
FUNCTION lnAw(x1,T)
  
  IMPLICIT NONE
  
  REAL, intent(in) :: T,x1
  REAL :: lnAw
  REAL, external :: m111,m121,m221,m122
  REAL, external :: e111,e121,e211,e122,e212,e221
  
  lnAw=-( &
       &  (2*m121(T)+e121(T)*log(x1)+e211(T)*(log(1-x1)+1))*x1 &
       &  +(2*m221(T)+e221(T)*(2*log(1-x1)+1))*(1-x1) &
       &  -(m111(T)+e111(T)*(log(x1)+1))*x1*x1 &
       & -(2*m121(T)+e121(T)*(log(x1)+1) &
       &            +e211(T)*(log(1-x1)+1))*x1*(1-x1) &
       &        -(m221(T)+e221(T)*(log(1-x1)+1))*(1-x1)**2 &
       &   +x1*(2*m122(T)+e122(T)*log(x1)+e212(T)*log(1-x1))*x1*(1-x1) &
       &  +x1*(1-x1)*((2*m122(T)+e122(T)*log(x1) &
       &                        +e212(T)*(log(1-x1)+1))*x1 &
       &               -(6*m122(T)+e122(T)*(3*log(x1)+1) &
       &                        +e212(T)*(3*log(1-x1)+1))*(1-x1)*x1) &
       &     )

END FUNCTION lnAw
!****************************************************************
FUNCTION m111(T)
        
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: m111
    
  m111=-23.524503387D0 +  &
       &    0.0406889449841D0*T -  & 
       &    0.151369362907D-4*T**2+2961.44445015D0/T + &
       &    0.492476973663D0*log(T) 

END FUNCTION m111

!****************************************************************
FUNCTION m121(T)
  
  IMPLICIT NONE
  REAL, INTENT(in) :: T
  REAL :: m121  

  m121=1114.58541077D0-1.1833078936D0*T - &
       &    0.00209946114412D0*T**2-246749.842271D0/T + &
       &    34.1234558134D0*log(T)

END FUNCTION m121

!****************************************************************
FUNCTION m221(T)
  
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: m221
  
  m221=-80.1488100747D0-0.0116246143257D0*T  &
       &    +0.606767928954D-5*T**2+3092.72150882D0/T  &
       &    +12.7601667471D0*log(T)

END FUNCTION m221

!****************************************************************
FUNCTION m122(T)

  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: m122
  
  m122=888.711613784D0-2.50531359687D0*T + &
       &    0.000605638824061D0*T**2-196985.296431D0/T + &
       &    74.550064338D0*log(T)

END FUNCTION m122

!****************************************************************
FUNCTION e111(T)
  
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: e111
  
  e111=2887.31663295D0-3.32602457749D0*T &
       &    -0.2820472833D-2*T**2-528216.112353D0/T & 
       &    +0.68699743564D0*log(T)

END FUNCTION e111

!****************************************************************
FUNCTION e121(T)
    
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: e121
  
  e121=-370.944593249D0-0.690310834523D0*T &
       &    +0.56345508422D-3*T**2-3822.52997064D0/T &
       &    +94.2682037574D0*log(T)

END FUNCTION e121
  
!****************************************************************
FUNCTION e211(T)
 
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: e211
  
  e211=38.3025318809D0-0.0295997878789D0*T &
       &    +0.120999746782D-4*T**2-3246.97498999D0/T &
       &    -3.83566039532D0*log(T)

END FUNCTION e211

!****************************************************************
FUNCTION e221(T)
  
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: e221
  
  e221=2324.76399402D0-0.141626921317D0*T &
       &    -0.00626760562881D0*T**2-450590.687961D0/T &
       &    -61.2339472744D0*log(T)

END FUNCTION e221
  
!****************************************************************
FUNCTION e122(T)
  
  IMPLICIT NONE
  
  REAL, INTENT(in) :: T
  REAL :: e122
  
  e122=-1633.85547832D0-3.35344369968D0*T &
       &    +0.00710978119903D0*T**2+198200.003569D0/T &
       &    +246.693619189D0*log(T)

END FUNCTION e122

!****************************************************************
FUNCTION e212(T)
  
  IMPLICIT NONE
  
  REAL :: e212
  REAL, INTENT(in) :: T
  
  e212=1273.75159848D0+1.03333898148D0*T &
       &    +0.00341400487633D0*T**2+195290.667051D0/T &
       &    -431.737442782D0*log(T)

END FUNCTION e212
