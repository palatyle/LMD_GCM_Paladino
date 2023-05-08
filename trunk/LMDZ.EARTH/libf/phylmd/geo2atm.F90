!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/geo2atm.F90,v 1.1 2008-12-05 17:56:40 lsce Exp $
!
SUBROUTINE geo2atm(im, jm, px, py, pz, plon, plat, pu, pv, pr)
  USE dimphy
  USE mod_phys_lmdz_para

  IMPLICIT NONE
  INCLUDE 'dimensions.h'
  INCLUDE 'YOMCST.h'

! Change wind coordinates from cartesian geocentric to local spherical
! NB! Fonctionne probablement uniquement en MPI seul (sans OpenMP)
!
  INTEGER, INTENT (IN)                 :: im, jm
  REAL, DIMENSION (im,jm), INTENT(IN)  :: px, py, pz
  REAL, DIMENSION (im,jm), INTENT(IN)  :: plon, plat
  REAL, DIMENSION (im,jm), INTENT(OUT) :: pu, pv, pr

  REAL :: rad


  rad = rpi / 180.0E0
  
  pu(:,:) = &
       - px(:,:) * SIN(rad * plon(:,:)) &
       + py(:,:) * COS(rad * plon(:,:))

  pv(:,:) = &
       - px(:,:) * SIN(rad * plat(:,:)) * COS(rad * plon(:,:)) &
       - py(:,:) * SIN(rad * plat(:,:)) * SIN(rad * plon(:,:)) &
       + pz(:,:) * COS(rad * plat(:,:))  

  pr(:,:) = &
       + px(:,:) * COS(rad * plat(:,:)) * COS(rad * plon(:,:)) &
       + py(:,:) * COS(rad * plat(:,:)) * SIN(rad * plon(:,:)) &
       + pz(:,:) * SIN(rad * plat(:,:))

  ! Value at North Pole
  IF (is_north_pole) THEN
     pu(:, 1) = -px (1,1)
     pv(:, 1) = -py (1,1)
     pr(:, 1) = 0.0
  ENDIF
  
  ! Value at South Pole     
  IF (is_south_pole) THEN
     pu(:,jm) = -px (1,jm)
     pv(:,jm) = -py (1,jm)
     pr(:,jm) = 0.0
  ENDIF
  
END SUBROUTINE geo2atm
