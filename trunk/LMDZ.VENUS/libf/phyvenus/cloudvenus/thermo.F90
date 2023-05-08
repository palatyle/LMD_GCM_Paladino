
!  FUNCTION    FPLAIR     Mean free path of air molecules (m)
!  FUNCTION    VISAIR     Dynamic viscosity of air
!  FUNCTION    DFWVA      Diffusivity of water vapor in air
!  FUNCTION    STSAS      Surface tension of H2SO4 solution/vapor
!  FUNCTION    ROSAS      Density of liquid sulfuric acid solution
!  FUNCTION    waterps    Saturation vapour pressure of pure water
!  FUNCTION    CDTAIR     Thermal conduvtivity of air


!*****************************************************************************
FUNCTION FPLAIR(T,P)

  ! Molecular mean free path of air molecules
  ! Source: Seinfield's book (2006,p.399)

  use free_param
  use donnees

  IMPLICIT NONE

  REAL :: FPLAIR, T, P, VISAIR

  FPLAIR=sqrt((PI*RGAS*T)/(2.0D0*MAIR))*(VISAIR(T)/P)
  
  RETURN

END FUNCTION FPLAIR


!*****************************************************************************
FUNCTION VISAIR(T)

  ! Dynamic viscosity of air.
  ! Source: Jones 1942

  ! Input:  TAIR:   Temperature (K)
  ! Output: VISAIR: Dynamic viscosity of air  (kg/(m s))=(Pa s)
  
  use free_param
  use donnees

  IMPLICIT NONE
  
  REAL :: T, VISAIR
  REAL :: AA, SS, T0

  AA = (5.27D0-3.0D0)/(5.27D0 -1.0D0)
  SS = -0.435D0
  T0 = 200.0D0 

  VISAIR=1015.0D0*((T/T0)**(0.5D0))*(T0**(AA)+SS)/(T**(AA)+SS) 
  VISAIR=VISAIR*1.D-8 
  
  RETURN

END FUNCTION VISAIR


!*****************************************************************************
FUNCTION DFWVA(T,P)

  !     Diffusivity of water vapor in air.
  !     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
  !                               (1980), 13-3, p. 413
  !     The relation D = E0 (T/T0)**n (P0/P); n=1.94 has been used.
  
  !     Input:  TAIR: Temperature (K);  Range:  [180,273]
  !             PAIR: Pressure (Pa)
  !     Output: Diffusivity of water vapor in air  (m**2/sec)
  
  use free_param
  use donnees

  IMPLICIT NONE

  REAL :: E0, D1, P0, T0, T, P, DFWVA

  PARAMETER(E0=0.211D-4,P0=1.01325D+5,T0=273.15D0,D1=E0*P0)
  
  DFWVA=D1*((T/T0)**1.94D0)/P
  
  RETURN

END FUNCTION DFWVA


!*****************************************************************************
FUNCTION STSAS(T,xmass) 
  !     Input:  T: Temperature (K)
  !             xmass: Mass fraction of H2SO4  [0;1]
  !     Output: Surface tension of sulfuric acid solution (N/m)

! about 230-323 K , x=0,...,1
!(valid down to the solid phase limit temp, which depends on molefraction)

  use donnees
      IMPLICIT NONE
      REAL :: STSAS
      REAL, INTENT(IN):: xmass, T
      REAL :: a, b, T1, Tc, xmole

       IF (T .LT. 305.15) THEN

!low temperature surface tension
! Hanna Vehkam‰ki and Markku Kulmala and Ismo Napari
! and Kari E. J. Lehtinen and Claudia Timmreck and Madis Noppel and Ari Laaksonen, 2002,
! An improved parameterization for sulfuric acid/water nucleation rates for tropospheric
!and stratospheric conditions, () J. Geophys. Res., 107, pp. 4622-4631

          a= 0.11864 + xmass* (-0.11651    &   
	             + xmass* ( 0.76852    &
		     + xmass* (-2.40909    &
		     + xmass*  (2.95434    &
		     + xmass* (-1.25852)))))
		     
          b= -0.00015709 + xmass* (0.00040102   &
	                 + xmass*(-0.00239950   &
			 + xmass* (0.007611235  &
			 + xmass*(-0.00937386   &
			 + xmass*0.00389722))))
          STSAS=a+T*b

      ELSE

      xmole = (xmass/MSA)*(1./((xmass/MSA)+(1.-xmass)/MWV))
      
! high temperature surface tension
!H. Vehkam‰ki and M. Kulmala and K.E. J. lehtinen, 2003,
!Modelling binary homogeneous nucleation of water-sulfuric acid vapours:
! parameterisation for high temperature emissions, () Environ. Sci. Technol., 37, 3392-3398

      Tc=    647.15*(1.0-xmole)*(1.0-xmole)   &
         +   900.0 *   xmole   *  xmole       &
         + 3156.186*   xmole   *(1-xmole)          !critical temperature
      T1=1.0-T/Tc

      a= 0.2358 + xmole*(-0.529     &
                + xmole* (4.073     &
		+ xmole*(-12.6707   &
		+ xmole* (15.3552   &
		+ xmole*(-6.3138)))))
		
      b= -0.14738 + xmole* (0.6253   &
                  + xmole*(-5.4808   &
		  + xmole*(17.2366   &
		  + xmole*(-21.0487  &
		  + xmole*(8.719)))))
      STSAS=(a+b*T1)*T1**(1.256)

      END IF
      RETURN

END FUNCTION STSAS


!*****************************************************************************
FUNCTION ROSAS(T,xmass)

!
! calculates the density of the liquid in kg/m^3
! xmass=mass fraction of h2so4, T in kelvins
! Hanna Vehkam‰ki and Markku Kulmala and Ismo Napari
! and Kari E. J. Lehtinen and Claudia Timmreck and Madis Noppel and Ari Laaksonen, 2002,
! An improved parameterization for sulfuric acid/water nucleation rates for tropospheric
!and stratospheric conditions, () J. Geophys. Res., 107, pp. 4622-4631

! about 220-373 K , x=0,...,1
!(valid down to the solid phase limit temp, which depends on molefraction)

      IMPLICIT NONE
      REAL :: ROSAS
      REAL, INTENT(IN) :: T, xmass
      REAL ::  a,b,c


      a= 0.7681724 + xmass* (2.1847140         &
                   + xmass* (7.1630022         &
                   + xmass* (-44.31447         &
                   + xmass* (88.75606          &
                   + xmass*(-75.73729          &
                   + xmass*  23.43228 )))))

      b= 1.808225e-3 + xmass* (-9.294656e-3       &
                     + xmass* (-0.03742148        &
                     + xmass* (0.2565321          &
                     + xmass* (-0.5362872         &
                     + xmass* (0.4857736          &
                     + xmass* (-0.1629592))))))
     
      c= -3.478524e-6 + xmass* (1.335867e-5      &
                      + xmass* (5.195706e-5      &
                      + xmass*(-3.717636e-4      &
                      + xmass* (7.990811e-4      &
                      + xmass*(-7.458060e-4      &
                      + xmass*  2.58139e-4)))))
     
      ROSAS= a+T*(b+c*T) ! g/cm^3
      ROSAS= ROSAS*1.0e3 !kg/m^3
      
      RETURN
END FUNCTION ROSAS


!****************************************************************
FUNCTION waterps(t)

  !     Saturation vapour pressure of pure water in Pa
  !     temperature t in K

  !     for 0 to 100C: Wexler 1976
  !     for <0C (validity range 123-332K): Murphy and Koop 2005

  use free_param
  use donnees

  IMPLICIT NONE

  REAL:: waterps, t,w

  if(t .ge. 273.15D0) then
     waterps=exp(-2991.2729D0*(t**(-2.))-6017.0128D0/t+18.87643854D0 &
          &        -0.028354721D0*t+0.17838301D-4*t**2.-0.84150417D-9*t**3. &
          &        +0.44412543D-12*t**4.+2.858487D0*LOG(t))
  else if(t .lt. 273.15D0) then
     waterps=exp(54.842763D0-6763.22D0/t-4.210D0*LOG(t)+0.000367D0*t &
          &        + tanh(0.0415D0*(t- 218.8D0))*(53.878D0- 1331.22D0/t &
          &        - 9.44523D0*LOG(t) + 0.014025D0*t))
  else
     stop 'no good temperatures in waterps!'
  endif

END FUNCTION waterps


!****************************************************************
FUNCTION CDTAIR(T)

  !     Thermal conduvtivity of air
  !     Source: Prupacher & Klett:Microphysics of clouds and precipitation,
  !                               (1980), p 418, 13-16
  !     Formula used:             CDTAIR=4.381276E-3+7.117560E-5*TAIR
  
  !     Input:  TAIR:     Air temperature (K)
  !     Output:           Thermal conductivity of air (J/(m sec K))
  
  use free_param
  use donnees
  
  IMPLICIT NONE

  REAL :: CDTAIR, T

  CDTAIR=4.381276D-3+7.117560D-5*T

  RETURN

END FUNCTION CDTAIR

