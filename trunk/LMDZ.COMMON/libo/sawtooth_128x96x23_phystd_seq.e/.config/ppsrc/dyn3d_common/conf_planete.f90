










!
! $Id$
!
SUBROUTINE conf_planete
!
! if not using IOIPSL, we still need to use (a local version of) getin
USE ioipsl_getincom
USE comconst_mod, ONLY: pi, g, molmass, kappa, cpp, omeg, rad, &
                        year_day, daylen, daysec, ihf
USE comvert_mod, ONLY: preff, pa
IMPLICIT NONE
!
!
!   Declarations :
!   --------------

!
!   local:
!   ------

! ---------------------------------------------
! Initialisations de constantes de la dynamique
! ---------------------------------------------
! Pi
pi=2.*asin(1.)

!Reference surface pressure (Pa)
preff=101325.
CALL getin('preff', preff)
! Reference pressure at which hybrid coord. become purely pressure
! pa=50000.
pa=preff/2.
CALL getin('pa', pa)
! Gravity
g=9.80665
CALL getin('g',g)
! Molar mass of the atmosphere
molmass = 28.9644
CALL getin('molmass',molmass)
! kappa=R/Cp et Cp      
kappa = 2./7.
CALL getin('kappa',kappa)
cpp=8.3145/molmass/kappa*1000.
CALL getin('cpp',cpp)
! Radius of the planet
rad = 6371229. 
CALL getin('radius',rad)
! Length of a standard day (s)
daysec=86400.
CALL getin('daysec',daysec)
! Rotation rate of the planet:
! Length of a solar day, in standard days
daylen = 1.
CALL getin('daylen',daylen)
! Number of days (standard) per year:
year_day = 365.25
CALL getin('year_day',year_day)
! Omega
! omeg=2.*pi/86400.
omeg=2.*pi/daysec*(1./daylen+1./year_day)
CALL getin('omeg',omeg)

! Intrinsic heat flux (default: none) (only used if planet_type="giant")
ihf = 0.
call getin('ihf',ihf)

END SUBROUTINE conf_planete
