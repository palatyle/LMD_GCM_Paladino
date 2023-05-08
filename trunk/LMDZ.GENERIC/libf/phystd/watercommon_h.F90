module watercommon_h

      implicit none

      real, parameter :: T_coup = 234.0
      real, parameter :: T_h2O_ice_liq = 273.16
      real, parameter :: T_h2O_ice_clouds = T_h2O_ice_liq-15. 
      real, parameter :: mH2O = 18.01528   

      ! benjamin additions
      real, parameter :: RLVTT = 2.257E+6 ! Latent heat of vaporization (J kg-1) 
      real, parameter :: RLSTT = 2.257E+6 ! 2.591E+6 in reality ! Latent heat of sublimation (J kg-1)

      real, parameter :: RLFTT = 3.334E+5 ! Latent heat of fusion (J kg-1) ! entails an energy sink but better description of albedo
      real, parameter :: rhowater = 1.0E+3 ! mass of water (kg/m^3)
      real, parameter :: rhowaterice = 9.2E+2 ! mass of water (kg/m^3)
      real, parameter :: capcal_h2o_liq = 4181.3 ! specific heat capacity of liquid water J/kg/K
      real, parameter :: mx_eau_sol = 150 ! mass of water (kg/m^2)

      real, save :: epsi, RCPD, RCPV, RV, RVTMP2
      real, save :: RETV
      real, save :: RLvCp
!$OMP THREADPRIVATE(epsi,RCPD,RCPV,RV,RVTMP2)
      
      contains

!==================================================================
      subroutine su_watercycle

      use comcstfi_mod, only: r, cpp, mugaz
      implicit none


!==================================================================
!
!     Purpose
!     -------
!     Set up relevant constants and parameters for the water cycle, and water cloud properties
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!     Jeremy Leconte (2012)
!
!==================================================================

      epsi   = mH2O / mugaz
      RCPD   = cpp 

      !RV = 1000.*R/mH2O
      RV = 1000.*8.314/mH2O ! caution! R is R/mugaz already!

      RCPV   = 1.88e3 ! specific heat capacity of water vapor at 350K

      RVTMP2 = RCPV/RCPD-1. ! not currently used...
      
! AB : initializations added for the thermal plume model
      RETV = RV / r - 1.
      RLvCp = RLVTT / RCPD
      
      end subroutine su_watercycle
      


!==================================================================
      subroutine Psat_water(T,p,psat,qsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the saturation vapor pressure and mass mixing ratio at saturation (kg/kg)
!     for a given pressure (Pa) and temperature (K)
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real, intent(in) :: T, p
  
!        output
         real psat,qsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Trefsublimation=7.66
         real,parameter :: Tmin=8.
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

! checked vs. old watersat data 14/05/2012 by JL.

         if (T.gt.T_h2O_ice_liq) then 
            psat = Pref_solid_liquid*Exp(r3vaporization*(T-T_h2O_ice_liq)/(T-Trefvaporization)) ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in psat water"
          !  psat = Pref_solid_liquid*Exp(r3sublimation*(Tmin-T_h2O_ice_liq)/(Tmin-Trefsublimation)) ! min psat  
         ! Ehouarn: gfortran says: Error: Result of EXP underflows its kind,
         !          so set psat to the smallest possible value instead
            psat=tiny(psat)
         else                 
            psat = Pref_solid_liquid*Exp(r3sublimation*(T-T_h2O_ice_liq)/(T-Trefsublimation)) ! solid / vapour
         endif
         if(psat.gt.p) then
            qsat=1.
         else
            qsat=epsi*psat/(p-(1.-epsi)*psat)
         endif
         return
      end subroutine Psat_water




!==================================================================
      subroutine Lcpdqsat_water(T,p,psat,qsat,dqsat,dlnpsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute dqsat=L/cp*d (q_sat)/d T and dlnpsat=L/cp d(ln Psat)/d T
!     for a given temperature (K)! 
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real T, p, psat, qsat
  
!        output
         real dqsat,dlnpsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Tmin=8.
         real,parameter :: Trefsublimation=7.66
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

         real :: dummy

         if (psat.gt.p) then
	    dqsat=0.
	    return
	 endif

         if (T.gt.T_h2O_ice_liq) then
            dummy = r3vaporization*(T_h2O_ice_liq-Trefvaporization)/(T-Trefvaporization)**2  ! liquid / vapour
         else if (T.lt.Tmin) then
	    print*, "careful, T<Tmin in Lcp psat water"
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(Tmin-Trefsublimation)**2  ! solid / vapour
         else               
            dummy = r3sublimation*(T_h2O_ice_liq-Trefsublimation)/(T-Trefsublimation)**2  ! solid / vapour
         endif

         dqsat=RLVTT/RCPD*qsat*(p/(p-(1.-epsi)*psat))*dummy
	 dlnpsat=RLVTT/RCPD*dummy
         return
      end subroutine Lcpdqsat_water




!==================================================================
      subroutine Tsat_water(p,Tsat)

         implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the saturation temperature
!     for a given pressure (Pa) 
!     Based on the Tetens formula from L.Li physical parametrization manual
!
!     Authors
!     -------
!     Jeremy Leconte (2012)
!
!==================================================================

!        input
         real p
  
!        output
         real Tsat

! JL12 variables for tetens formula
         real,parameter :: Pref_solid_liquid=611.14
         real,parameter :: Trefvaporization=35.86
         real,parameter :: Trefsublimation=7.66
         real,parameter :: r3vaporization=17.269
         real,parameter :: r3sublimation=21.875

         if (p.lt.Pref_solid_liquid) then ! solid / vapour
            Tsat =(T_h2O_ice_liq*r3sublimation- Trefsublimation*Log(p/Pref_solid_liquid))/(r3sublimation-Log(p/Pref_solid_liquid))
         else                 ! liquid / vapour
            Tsat =(T_h2O_ice_liq*r3vaporization- Trefvaporization*Log(p/Pref_solid_liquid))/(r3vaporization-Log(p/Pref_solid_liquid))
         endif

         return
      end subroutine Tsat_water
!==================================================================


!==================================================================
subroutine watersat(T,p,qsat)
  implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the water mass mixing ratio at saturation (kg/kg)
!     for a given pressure (Pa) and temperature (K)
!     A replacement for the old watersat.F in the Martian GCM.
!     Based on FCTTRE.h in the LMDTERRE model.
!
!     JL18 watersat was used only in vdifc and thus it was not consistent with other routines (turbdiff, rain, largescale...)
!        which used Psat_water. This is now harmonized
!        we put it here for archival purpose, but it is not used anymore.
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!
!==================================================================

!   input
  real T, p
  
!   output
  real qsat

! checked vs. NIST data 22/06/2010 by RW.
! / by p gives partial pressure
! x by epsi converts to mass mixing ratio

  if (T.lt.T_h2O_ice_liq) then ! solid / vapour
     qsat = 100.0 * 10**(2.07023 - 0.00320991             &
          * T - 2484.896 / T + 3.56654 * alog10(T))
  else                 ! liquid / vapour
     qsat = 100.0 * 10**(23.8319 - 2948.964 / T - 5.028  &
          * alog10(T) - 29810.16 * exp( -0.0699382 * T)  &
          + 25.21935 * exp(-2999.924/T))
  endif
!  qsat=epsi*qsat/p
  if(qsat.gt.p) then
     qsat=1.
  else
     qsat=epsi*qsat/(p-(1.-epsi)*qsat)
  endif
  return
end subroutine watersat

subroutine watersat_grad(T,qsat,dqsat)

  implicit none

!==================================================================
!     Purpose
!     -------
!     Compute the L/cp*d (q_sat)/d T
!     for a given temperature (K)
!
!     JL18: see watersat
!
!     Authors
!     -------
!     Robin Wordsworth (2010)
!
!==================================================================

!   input
  real T,qsat
  
!   output
  real dqsat

!  if (T.lt.T_coup) then ! solid / vapour !why use T_coup?????????? JL12
  if (T.lt.T_h2O_ice_liq) then ! solid / vapour
     dqsat = RLVTT/RCPD*qsat*(3.56654/T             &
          +2484.896*LOG(10.)/T**2                   &
          -0.00320991*LOG(10.))
  else                 ! liquid / vapour
     dqsat = RLVTT/RCPD*qsat*LOG(10.)*              &
          (2948.964/T**2-5.028/LOG(10.)/T           &
          +25.21935*2999.924/T**2*EXP(-2999.924/T)  &
          +29810.16*0.0699382*EXP(-0.0699382*T))
  end if

  return
end subroutine watersat_grad



end module watercommon_h
