SUBROUTINE param_slope( &
  !
  ! INPUTS
  !
               &  csza, declin, rho, latitude &
               &  ,taudust, albedo  &
               &  ,theta_s, psi_s &
               &  ,fdir_0, ftot_0  &
  !
  ! OUTPUTS
  !
               &  ,ftot &
 )


!!*****************************************************************************************
!
! SUBROUTINE:
! param_slope
! 
!
! PURPOSE: 
! computes total solar irradiance on a given Martian slope 
!
!
! INPUTS:
! csza		cosine solar zenith angle
! declin	sun declination (rad)
! rho		sun right ascension (rad)
! latitude	latitude (deg)
! taudust	dust optical depth at reference wavelength 0.67 mic. 
! albedo	spectrally integrated surface Lambertian reflection albedo 
! theta_s	slope inclination angle (deg)
!			0 is horizontal, 90 is vertical
! phi_s		slope azimuth (deg) 
!       		0 >> Northward
!       		90 >> Eastward
!       		180 >> Southward
!       		270 >> Westward
! ftot_0	spectrally integrated total irradiance on an horizontal surface (W/m2)
! fdir_0	spectrally integrated direct irradiance on an horizontal surface (W/m2)
!
!
! OUTPUTS:
! ftot		spectrally integrated total irradiance on the slope (W/m2)
!
! REFERENCE: 
! "Fast and accurate estimation of irradiance on Martian slopes"
! A. Spiga & F. Forget
! .....
!
! AUTHOR: 
! A. Spiga (spiga@lmd.jussieu.fr)
! March 2008
!
!!*****************************************************************************************

IMPLICIT NONE

!!
!! INPUT
!!
REAL, INTENT(IN) :: csza, declin, rho, latitude, taudust, theta_s, psi_s, albedo, ftot_0 , fdir_0 

!!
!! LOCAL
!!
REAL :: pi, deg2rad
REAL :: a
REAL :: mu_s, sigma_s
REAL :: fdir, fscat, fscat_0, fref
REAL, DIMENSION(4,2) :: mat_M, mat_N, mat_T
REAL, DIMENSION(2) :: g_vector
REAL, DIMENSION(4) :: s_vector
REAL :: ratio

!!
!! OUTPUT
!!
REAL, INTENT(OUT) :: ftot

!!*****************************************************************************************

!
! Prerequisite
!
  pi = 2.*asin(1.)
  deg2rad = pi/180.
  if ((theta_s > 90.) .or. (theta_s < 0.)) then
        print *, 'please set theta_s between 0 and 90', theta_s
        stop 
  endif

!
! Solar Zenith angle (radian) 
!
  if (csza .lt. 0.01) then
      !print *, 'sun below horizon' 
      !fdir_0=0.
      fdir=0.
      fscat_0=0.
      fscat=0.   
      fref=0.
  else

!!
!! Low incidence fix
!! 
!  if (csza .lt. 0.15) csza = 0.15

!
! 'Slope vs Sun' azimuth (radian)
!
  if ( ( (cos(declin)*sin(rho)) .eq. 0.0 ) &
               .and. &
     ( ( sin(deg2rad*latitude)*cos(declin)*cos(rho) &
              -cos(deg2rad*latitude)*sin(declin) ) .eq. 0.0 ) ) then
    a = deg2rad*psi_s ! some compilator need specfying value for atan2(0,0)  
  else
    a = deg2rad*psi_s + atan2(cos(declin)*sin(rho), &
               sin(deg2rad*latitude)*cos(declin)*cos(rho)-cos(deg2rad*latitude)*sin(declin))
  end if
  
!
! Cosine of slope-sun phase angle 
!
  mu_s = csza*cos(deg2rad*theta_s) - cos(a)*sin(deg2rad*theta_s)*sqrt(1-csza**2)
  if (mu_s .le. 0.) mu_s=0.

!
! Sky-view factor
!
  sigma_s=0.5*(1.+cos(deg2rad*theta_s))

!
! Direct flux on the slope
!
  fdir = fdir_0 * mu_s/csza

!
! Reflected flux on the slope
!
  fref = albedo * (1-sigma_s) * ftot_0

!
! Scattered flux on a flat surface
!
  fscat_0 = ftot_0 - fdir_0

!
! Scattering vector (slope vs sky)
!
  s_vector=(/ 1., exp(-taudust) , sin(deg2rad*theta_s), sin(deg2rad*theta_s)*exp(-taudust) /)

!
! Geometry vector (slope vs sun)
!
  g_vector=(/ mu_s/csza, 1. /)

!
! Coupling matrix
!
  if (csza .ge. 0.5) then
        mat_M(:,1) = (/ -0.264,  1.309,  0.208, -0.828 /)
        mat_M(:,2) = (/  1.291*sigma_s, -1.371*sigma_s, -0.581,  1.641 /)
        mat_N(:,1) = (/  0.911, -0.777, -0.223,  0.623 /)
        mat_N(:,2) = (/ -0.933*sigma_s,  0.822*sigma_s,  0.514, -1.195 /)

  else
        mat_M(:,1) = (/ -0.373,  0.792, -0.095,  0.398 /)
        mat_M(:,2) = (/  1.389*sigma_s, -0.794*sigma_s, -0.325,  0.183 /)
        mat_N(:,1) = (/  1.079,  0.275,  0.419, -1.855 /)
        mat_N(:,2) = (/ -1.076*sigma_s, -0.357*sigma_s, -0.075,  1.844 /)
  endif
  !
  mat_T = mat_M + csza*mat_N


!
! Scattered flux slope ratio
!
  if (deg2rad*theta_s <= 0.0872664626) then
  !
  ! low angles
  !
        s_vector = (/ 1., exp(-taudust) , sin(0.0872664626), sin(0.0872664626)*exp(-taudust) /)
        ratio = DOT_PRODUCT ( MATMUL( s_vector, mat_T), g_vector )
        ratio = 1. + (ratio - 1.)*deg2rad*theta_s/0.0872664626
  else
  !
  ! general case
  !
        ratio= DOT_PRODUCT ( MATMUL( s_vector, mat_T), g_vector )  
                  !
                  ! NB: ratio= DOT_PRODUCT ( s_vector, MATMUL( mat_T, g_vector ) ) is equivalent
  endif

!
! Scattered flux on the slope
!
  fscat = ratio * fscat_0


endif !! if (csza < 0.01)

!
! Total flux on the slope
!
  ftot = fdir + fref + fscat

!!
!! Display results
!!
!  print *, 'sca component 0 ', ftot_0-fdir_0
!  print *, 'dir component 0 ', fdir_0
!  print *, 'scattered component ', fscat
!  print *, 'direct component ', fdir
!  print *, 'reflected component ', fref

END SUBROUTINE param_slope
