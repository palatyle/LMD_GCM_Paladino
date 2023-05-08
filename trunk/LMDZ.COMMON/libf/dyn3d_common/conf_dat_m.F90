MODULE conf_dat_m
!
!*******************************************************************************

  PRIVATE
  PUBLIC :: conf_dat2d, conf_dat3d


CONTAINS


!-------------------------------------------------------------------------------
!
SUBROUTINE conf_dat2d(title, xd, yd, xf, yf, champd, interbar)
!
!-------------------------------------------------------------------------------
! Author: P. Le Van
!-------------------------------------------------------------------------------
! Purpose: Configure the 2D data field "champd" so that:
!     - Longitudes are in [ -pi     pi   ]
!     - Latitudes  are in [  pi/2. -pi/2 ]
!  * xd / yd are initial lon / lats.
!  * xf / yf are output  lon / lats, possibly modified to satisfy configuration.
!  * interbar is TRUE for barycentric interpolation.
!-------------------------------------------------------------------------------
  USE assert_eq_m, ONLY: assert_eq
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN)    :: title
  REAL,             INTENT(IN)    :: xd(:), yd(:)  ! dim (lons) (lats)
  REAL,             INTENT(INOUT) :: xf(:), yf(:)  ! dim (lons) (lats)
  REAL,             INTENT(INOUT) :: champd(:,:)   ! dim (lons,lats)
  LOGICAL,          INTENT(IN)    :: interbar
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname="conf_dat2d"
  INTEGER :: i, j, ip180, ind, lons, lats
  LOGICAL :: radlon, invlon ,radlat, invlat
  REAL    :: pi, pis2, depi, rlatmin, rlatmax, oldxd1
  REAL, ALLOCATABLE :: xtemp(:) , ytemp(:), champf(:,:)
!-------------------------------------------------------------------------------
  lons=assert_eq(SIZE(xd),SIZE(xf),SIZE(champd,1),TRIM(modname)//" lons")
  lats=assert_eq(SIZE(yd),SIZE(yf),SIZE(champd,2),TRIM(modname)//" lats")
  ALLOCATE(xtemp(lons),ytemp(lats)); xtemp(:)=xd(:); ytemp(:)=yd(:)
  ALLOCATE(champf(lons,lats))
  pi = 2. * ASIN(1.) 
  pis2 = pi/2.
  depi = 2. * pi
  radlon=.FALSE.; invlon=.FALSE.
  IF     (xtemp(1)>=-pi-0.5.AND.xtemp(lons)<=  pi+0.5) THEN
    radlon = .TRUE.;  invlon = .FALSE.
  ELSE IF(xtemp(1)>=   -0.5.AND.xtemp(lons)<=depi+0.5) THEN
    radlon = .TRUE.;  invlon = .TRUE.
  ELSE IF(xtemp(1)>= -180.5.AND.xtemp(lons)<=   180.5) THEN
    radlon = .FALSE.; invlon = .FALSE.
  ELSE IF(xtemp(1)>=   -0.5.AND.xtemp(lons)<=   360.5) THEN
    radlon = .FALSE.; invlon = .TRUE.
  ELSE; WRITE(6,*) 'Problems with data longitudes for file '//TRIM(title)
  END IF
  invlat = ytemp(1)<ytemp(lats)
  rlatmin = MIN( ytemp(1), ytemp(lats) )
  rlatmax = MAX( ytemp(1), ytemp(lats) )
      
  IF     (rlatmin>=-pis2-0.5.AND.rlatmax<=pis2+0.5) THEN; radlat = .TRUE.
  ELSE IF(rlatmin>= -90.-0.5.AND.rlatmax<= 90.+0.5) THEN; radlat = .FALSE.
  ELSE; WRITE(6,*) 'Problems with data latitudes for file '//TRIM(title)
  END IF

  IF(.NOT.radlon) xtemp(:)=xtemp(:)*pi/180.
  IF(.NOT.radlat) ytemp(:)=ytemp(:)*pi/180.

!--- FLIPPED LONGITUDES
  IF(invlon) THEN
    champf(:,:)=champd(:,:); xf(:)=xtemp(:)

  !--- Longitudes rotated to get them in [ -pi pi] interval.
    DO i=1,lons;      IF(xf(i)>pi) EXIT;                   END DO; ip180 = i
    DO i=1,lons;      IF(xf(i)>pi) xf(i)=xf(i)-depi;       END DO
    DO i= ip180,lons; ind=i-ip180+1; xtemp(ind)=xf(i);     END DO
    DO i= ind+1,lons;                xtemp(i  )=xf(i-ind); END DO

  !--- Longitudes rotated in champf
    DO j=1,lats
      DO i=ip180,lons; ind=i-ip180+1; champd(ind,j)=champf(i    ,j); END DO
      DO i=ind+1,lons;                champd(i  ,j)=champf(i-ind,j); END DO
    END DO
  END IF

!--- FLIPPED LATITUDES
  IF(invlat) THEN
    yf(:)=ytemp(:)
    champf(:,:)=champd(:,:)
    ytemp(lats:1:-1)=yf(:)
    DO j=1,lats; champd(:,lats-j+1)=champf(:,j); END DO
  END IF

!--- FOR BARYCENTRIC INTERPOLATION
  IF(interbar) THEN
    oldxd1 = xtemp(1)
    xtemp(1:lons-1)=0.5*(xtemp(1:lons-1)+xtemp(2:lons))
    xtemp(  lons  )=0.5*(xtemp(  lons  )+oldxd1+depi)
    ytemp(1:lats-1)=0.5*(ytemp(1:lats-1)+ytemp(2:lats))
  END IF
  DEALLOCATE(champf)
  xf(:)=xtemp(:); DEALLOCATE(xtemp)
  yf(:)=ytemp(:); DEALLOCATE(ytemp)

END SUBROUTINE conf_dat2d
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE conf_dat3d(title, xd, yd, zd, xf, yf, zf, champd, interbar)
!
!-------------------------------------------------------------------------------
! Author: P. Le Van
!-------------------------------------------------------------------------------
! Purpose: Configure the 3D data field "champd" so that:
!     - Longitudes are in [ -pi     pi   ]
!     - Latitudes  are in [  pi/2. -pi/2 ]
!     - Vertical levels from ground to model top (in Pascals)
!  * xd / yd are initial lon / lats.
!  * xf / yf are output  lon / lats, possibly modified to satisfy configuration.
!  * zf      are output pressures
!  * interbar is TRUE for barycentric interpolation.
!-------------------------------------------------------------------------------
  USE assert_eq_m, ONLY: assert_eq
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  CHARACTER(LEN=*), INTENT(IN)    :: title
  REAL,             INTENT(IN)    :: xd(:), yd(:), zd(:) ! (lons) (lats) (levs)
  REAL,             INTENT(INOUT) :: xf(:), yf(:), zf(:) ! (lons) (lats) (levs)
  REAL,             INTENT(INOUT) :: champd(:,:,:)       ! (lons,lats,levs)
  LOGICAL,          INTENT(IN)    :: interbar
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname="conf_dat3d"
  INTEGER :: i, j, l, ip180, ind, lons, lats, levs
  LOGICAL :: radlon, invlon ,radlat, invlat, invlev
  REAL    :: pi, pis2, depi, presmax, rlatmin, rlatmax, oldxd1
  REAL, ALLOCATABLE :: xtemp(:) , ytemp(:), ztemp(:), champf(:,:,:)
!-------------------------------------------------------------------------------
  lons=assert_eq(SIZE(xd),SIZE(xf),SIZE(champd,1),TRIM(modname)//" lons")
  lats=assert_eq(SIZE(yd),SIZE(yf),SIZE(champd,2),TRIM(modname)//" lats")
  levs=assert_eq(SIZE(zd),SIZE(zf),SIZE(champd,3),TRIM(modname)//" levs")
  ALLOCATE(xtemp(lons),ytemp(lats),ztemp(levs),champf(lons,lats,levs))
  xtemp(:)=xd(:); ytemp(:)=yd(:); ztemp(:)=zd(:)
  pi = 2. * ASIN(1.) 
  pis2 = pi/2.
  depi = 2. * pi
  radlon=.FALSE.; invlon=.FALSE.
  IF     (xtemp(1)>=-pi-0.5.AND.xtemp(lons)<=  pi+0.5) THEN
    radlon = .TRUE.;  invlon = .FALSE.
  ELSE IF(xtemp(1)>=   -0.5.AND.xtemp(lons)<=depi+0.5) THEN
    radlon = .TRUE.;  invlon = .TRUE.
  ELSE IF(xtemp(1)>= -180.5.AND.xtemp(lons)<=   180.5) THEN
    radlon = .FALSE.; invlon = .FALSE.
  ELSE IF(xtemp(1)>=   -0.5.AND.xtemp(lons)<=   360.5) THEN
    radlon = .FALSE.; invlon = .TRUE.
  ELSE; WRITE(6,*) 'Problems with data longitudes for file '//TRIM(title)
  END IF
  invlat = ytemp(1)<ytemp(lats)
  rlatmin = MIN( ytemp(1), ytemp(lats) )
  rlatmax = MAX( ytemp(1), ytemp(lats) )
      
  IF     (rlatmin>=-pis2-0.5.AND.rlatmax<=pis2+0.5) THEN; radlat = .TRUE.
  ELSE IF(rlatmin>= -90.-0.5.AND.rlatmax<= 90.+0.5) THEN; radlat = .FALSE.
  ELSE; WRITE(6,*) 'Problems with data latitudes for file '//TRIM(title)
  END IF

  IF(.NOT.radlon) xtemp(:)=xtemp(:)*pi/180.
  IF(.NOT.radlat) ytemp(:)=ytemp(:)*pi/180.

!--- FLIPPED LONGITUDES
  IF(invlon) THEN
    champf(:,:,:)=champd(:,:,:); xf(:)=xtemp(:)

  !--- Longitudes rotated to get them in [ -pi pi] interval.
    DO i=1,lons;      IF(xf(i)>pi) EXIT;                   END DO; ip180 = i
    DO i=1,lons;      IF(xf(i)>pi) xf(i)=xf(i)-depi;       END DO
    DO i= ip180,lons; ind=i-ip180+1; xtemp(ind)=xf(i);     END DO
    DO i= ind+1,lons;                xtemp(i  )=xf(i-ind); END DO

  !--- Longitudes rotated in champf
    DO l=1,levs
    DO j=1,lats
      DO i=ip180,lons; ind=i-ip180+1; champd(ind,j,l)=champf(i    ,j,l); END DO
      DO i=ind+1,lons;                champd(i  ,j,l)=champf(i-ind,j,l); END DO
    END DO
    END DO
  END IF

!--- FLIPPED LATITUDES
  IF(invlat) THEN
    yf(:)=ytemp(:)
    champf(:,:,:)=champd(:,:,:)
    ytemp(lats:1:-1)=yf(:)
    DO l=1,levs
      DO j=1,lats; champd(:,lats-j+1,l)=champf(:,j,l); END DO
    END DO
  END IF

!--- FOR BARYCENTRIC INTERPOLATION
  IF(interbar) THEN
    oldxd1 = xtemp(1)
    xtemp(1:lons-1)=0.5*(xtemp(1:lons-1)+xtemp(2:lons))
    xtemp(  lons  )=0.5*(xtemp(  lons  )+oldxd1+depi)
    ytemp(1:lats-1)=0.5*(ytemp(1:lats-1)+ytemp(2:lats))
  END IF

!--- FLIPPED LEVELS
  invlev=ztemp(1)<ztemp(levs)
  IF(MAX(ztemp(1),ztemp(levs))<1200.) ztemp(:)=ztemp(:)*100.
  IF(invlev) THEN
    zf(:)=ztemp(:)
    champf(:,:,:)=champd(:,:,:)
    ztemp(levs:1:-1)=zf(:)
    DO l=1,levs; champd(:,:,levs+1-l)=champf(:,:,l); END DO
  END IF

  DEALLOCATE(champf)
  xf(:)=xtemp(:); DEALLOCATE(xtemp)
  yf(:)=ytemp(:); DEALLOCATE(ytemp)
  zf(:)=ztemp(:); DEALLOCATE(ztemp)

END SUBROUTINE conf_dat3d
!
!-------------------------------------------------------------------------------

END MODULE conf_dat_m

