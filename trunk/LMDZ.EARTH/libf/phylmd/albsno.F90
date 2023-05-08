!
! $Header$
!
SUBROUTINE albsno(klon, knon, dtime, agesno, alb_neig_grid, precip_snow)

  IMPLICIT NONE

! Input arguments
!****************************************************************************************
  INTEGER, INTENT(IN)                  :: klon, knon
  REAL, INTENT(IN)                     :: dtime
  REAL, DIMENSION(klon), INTENT(IN)    :: precip_snow

! In/Output arguments
!****************************************************************************************
  REAL, DIMENSION(klon), INTENT(INOUT) :: agesno

! Output arguments
!****************************************************************************************
  REAL, DIMENSION(klon), INTENT(OUT)   :: alb_neig_grid

! Local variables
!****************************************************************************************
  INTEGER                              :: i, nv
  INTEGER, PARAMETER                   :: nvm = 8 
  REAL                                 :: as
  REAL, DIMENSION(klon,nvm)            :: veget
  REAL, DIMENSION(nvm),SAVE            :: init, decay
  !$OMP THREADPRIVATE(init, decay)

  DATA init /0.55, 0.14, 0.18, 0.29, 0.15, 0.15, 0.14, 0./
  DATA decay/0.30, 0.67, 0.63, 0.45, 0.40, 0.14, 0.06, 1./
!****************************************************************************************

  veget = 0.
  veget(:,1) = 1.     ! desert partout
  DO i = 1, knon
     alb_neig_grid(i) = 0.0
  ENDDO
  DO nv = 1, nvm
     DO i = 1, knon
        as = init(nv)+decay(nv)*EXP(-agesno(i)/5.)
        alb_neig_grid(i) = alb_neig_grid(i) + veget(i,nv)*as
     ENDDO
  ENDDO
  

! modilation en fonction de l'age de la neige
  DO i = 1, knon
     agesno(i)  = (agesno(i) + (1.-agesno(i)/50.)*dtime/86400.)&
          &             * EXP(-1.*MAX(0.0,precip_snow(i))*dtime/0.3)
     agesno(i) =  MAX(agesno(i),0.0)
  ENDDO
  
END SUBROUTINE albsno
