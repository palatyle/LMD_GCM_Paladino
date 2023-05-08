SUBROUTINE season_haze(zday,lat,press,fact)

  !     ==============================================================================
  !     Purpose
  !     -------
  !     Compute haze opacity seasonal modulation factor based on Karkoschka 2016
  !     
  !     Authors
  !     -------
  !     J. Vatant d'Ollone (2018)
  !     ==============================================================================

  USE radinc_h

  !-----------------------------------------------------------------------
  !     Declarations:
  !     -------------

  IMPLICIT NONE

  !     Arguments :
  !     -----------
  REAL, INTENT(IN)                        :: zday  ! Time elapsed since Ls=0 (sols)
  REAL, INTENT(IN)                        :: lat   ! latitude of grid point
  REAL, DIMENSION(L_LEVELS), INTENT(IN)   :: press ! layers boundary pressure (mbar)
  REAL, DIMENSION(L_LEVELS), INTENT(OUT)  :: fact  ! Haze opacity seasonal factor

  !     Local variables :
  !     -----------------
  REAL :: M1,M2,D1,D2,E1,E2,B1,B2
  REAL :: pi = 4.0*atan(1.0)
  REAL :: sol2earthyr = 15.945 / 365.25
  REAL :: equinox = 2009.611 ! 11 August 2009 Ls=0
  INTEGER :: i
  !     -------------------------------------------------------------------------

  ! Mi, Ei and Di are fitted from Fig 14 from Karkoschka 2016

  M1 = 0.20543*atan(0.142666*abs(lat)-2.83289)-0.560454
  M2 = 1.09941*atan(-0.0584703*abs(lat)+1.30911)+0.0430241

  D1 = 2.30639*( cos(-0.0311457*(abs(lat)*180.0/pi-90.0)) - cos(0.0311457*90.0) )
  D2 = 68.3087*( cos(-0.0035113*(abs(lat)*180.0/pi-90.0)) - cos(0.0035113*90.0) )

  E1 = 0.644826*atan(0.145421*abs(lat)-4.50363)+2003.35
  E2 = 2001.93

  IF (lat .GE. 0.0 ) THEN
    B1 = M1 + D1 * sin( (equinox+zday*sol2earthyr-E1) * 2*pi / 29.4571 ) ! Eq. 3  from Karkoschka 2016
    B2 = M2 + D2 * sin( (equinox+zday*sol2earthyr-E2) * 2*pi / 29.4571 ) !          """
  ELSE
    B1 = M1 - D1 * sin( (equinox+zday*sol2earthyr-E1) * 2*pi / 29.4571 )
    B2 = M2 - D2 * sin( (equinox+zday*sol2earthyr-E2) * 2*pi / 29.4571 )
  ENDIF

  fact(:) = 1.0 ! default value -> no adjustment of reference haze profile

  DO i=1,L_LEVELS

  IF ( press(i).GT.20.0 ) THEN
    fact(i) = 5.57403*exp(0.0629317*B1)-4.23873 ! Fit from table 5 of Karkoschka 2016
  ELSE IF ( press(i).LT.2.0 ) THEN
    fact(i) = 4.83997*exp(0.0501367*B2)-3.83877 !           """
  ENDIF

    IF ( fact(i) .LT. 1.0E-6 ) fact(i) = 1.0E-6     ! Threshold for extreme values ( fit go neg and no haze at all seems unphysical )

  ENDDO


END SUBROUTINE season_haze
