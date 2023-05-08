










!*******************************************************************************
!
MODULE grid_atob_m
!
!*******************************************************************************

  USE assert_eq_m, ONLY: assert_eq

  PRIVATE
  PUBLIC :: grille_m, rugosite, sea_ice, rugsoro

CONTAINS

!-------------------------------------------------------------------------------
!
SUBROUTINE fine2coarse(x_i, y_i, x_o, y_o, d_o1, d_i, msk, d_o2)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL,              INTENT(IN)  :: x_i(:), y_i(:) !-- INPUT  X&Y COOR. (mi)(ni)
  REAL,              INTENT(IN)  :: x_o(:), y_o(:) !-- OUTPUT X&Y COOR. (mo)(no)
  REAL,              INTENT(OUT) :: d_o1(:,:)      !-- OUTPUT FIELD     (mo,no)
  REAL,    OPTIONAL, INTENT(IN)  :: d_i (:,:)      !-- INPUT FIELD      (mi,ni)
  LOGICAL, OPTIONAL, INTENT(IN)  :: msk (:,:)      !-- MASK             (mo,no)
  REAL,    OPTIONAL, INTENT(OUT) :: d_o2(:,:)      !-- OUTPUT FOR d_i^2 (mo,no)
!-------------------------------------------------------------------------------
! Local variables:
  CHARACTER(LEN=256) :: modname="fine2coarse"
  INTEGER :: mi, ni, ii, ji, mo, no, io, jo, nr(2), m1,m2, n1,n2, mx,my, nn, i,j
  LOGICAL :: li, lo, first=.TRUE.
  REAL    :: inc, cpa, spa, crlo(SIZE(x_i))
  REAL, SAVE :: pi, hpi
  INTEGER, DIMENSION(SIZE(x_o),SIZE(y_o)) :: num_tot
  LOGICAL, DIMENSION(SIZE(x_o),SIZE(y_o)) :: found, mask
  REAL,    DIMENSION(SIZE(x_i),SIZE(y_i)) :: dist
  REAL,    DIMENSION(SIZE(x_o))           :: a, b
  REAL,    DIMENSION(SIZE(y_o))           :: c, d
  REAL,    PARAMETER :: thresh=1.E-5
!-------------------------------------------------------------------------------
  IF(first) THEN; pi=4.0*ATAN(1.0); hpi=pi/2.0; first=.FALSE.; END IF
  mi=SIZE(x_i); ni=SIZE(y_i); mo=SIZE(x_o); no=SIZE(y_o)
  m1=m1; m2=mo; mx=mo; IF(PRESENT(msk)) mx=SIZE(msk,1)
  n1=ni; n2=no; my=no; IF(PRESENT(msk)) my=SIZE(msk,2)
  li=PRESENT(d_i ); IF(li) THEN; m1=SIZE(d_i ,1); n1=SIZE(d_i ,2); END IF
  lo=PRESENT(d_o2); IF(lo) THEN; m2=SIZE(d_o2,1); n2=SIZE(d_o2,2); END IF
  mi=assert_eq(mi,m1,TRIM(modname)//" mi")
  ni=assert_eq(ni,n1,TRIM(modname)//" ni")
  mo=assert_eq(mo,m2,mx,SIZE(d_o1,1),TRIM(modname)//" mo")
  no=assert_eq(no,n2,my,SIZE(d_o1,2),TRIM(modname)//" no")
  mask(:,:)=.TRUE.; IF(PRESENT(msk)) mask(:,:)=msk(:,:)

!--- COMPUTE CELLS INTERFACES COORDINATES OF OUTPUT GRID
  b(mo)=x_o(mo)+(x_o(mo)-x_o(mo-1))/2.0; b(1:mo-1)=(x_o(1:mo-1)+x_o(2:mo))/2.0
  d(no)=y_o(no)+(y_o(no)-y_o(no-1))/2.0; d(1:no-1)=(y_o(1:no-1)+y_o(2:no))/2.0
  a(1 )=x_o(1 )-(x_o(2 )-x_o(1   ))/2.0; a(2:mo  )=   b(1:mo-1)
  c(1 )=y_o(1 )-(y_o(2 )-y_o(1   ))/2.0; c(2:no  )=   d(1:no-1)

!--- ACCUMULATE INPUT POINTS ON OUTPUT GRID
  d_o1(:,:)=0.; num_tot(:,:)=0; inc=1.0
  IF(lo) d_o2(:,:)=0.
  DO ji = 1, ni
    DO ii = 1, mi
      IF(li) inc=d_i(ii,ji)
      DO jo = 1, no
        IF((y_i(ji)-c(jo)<thresh.OR.y_i(ji)-d(jo)>thresh).AND.   &
           (y_i(ji)-c(jo)>thresh.OR.y_i(ji)-d(jo)<thresh)) CYCLE
        DO io = 1, mo
          IF((x_i(ii)-a(io)<thresh.OR.x_i(ii)-b(io)>thresh).AND. &
             (x_i(ii)-a(io)>thresh.OR.x_i(ii)-b(io)<thresh)) CYCLE
          num_tot(io,jo)=num_tot(io,jo)+1
          IF(mask(io,jo)) d_o1(io,jo)=d_o1(io,jo)+inc
          IF(.NOT.lo) CYCLE
          IF(mask(io,jo)) d_o2(io,jo)=d_o2(io,jo)+inc*inc
        END DO
      END DO
    END DO
  END DO

!--- CHECK INPUT POINTS HAVE BEEN FOUND IN EACH OUTPUT CELL
  found(:,:)=num_tot(:,:)/=0
  WHERE(found.AND.mask) d_o1(:,:)=d_o1(:,:)/REAL(num_tot(:,:))
  IF(PRESENT(d_o2)) THEN
    WHERE(found.AND.mask) d_o2(:,:)=d_o2(:,:)/REAL(num_tot(:,:))
    RETURN
  END IF
  nn=COUNT(found); IF(nn==0) RETURN

!--- MISSING POINTS ; USE DISTANCE ON THE SPHERE TO FIND NEAREST POINT nr(2)
  DO io = 1, mo
    DO jo = 1, no
      IF(found(io,jo)) CYCLE
!      IF(prt_level>=1) PRINT*, "Problem: point out of domain (i,j)=", io,jo
      crlo(:)=COS(x_o(io)-x_i(:))     !--- COS of points 1 and 2 angle
      cpa=COS(y_o(jo)); spa=SIN(y_o(jo))
      DO j=1,ni; dist(:,j)=ACOS(spa*SIN(y_i(j))+cpa*COS(y_i(j))*crlo(:)); END DO
      nr=MINLOC(dist(:,:))!; IF(prt_level>=1) PRINT*, "Solution: ", nr
      inc=1.0; IF(li) inc=d_i(nr(1),nr(2))
      IF(mask(io,jo)) d_o1(io,jo)=inc
    END DO
  END DO

END SUBROUTINE fine2coarse
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE grille_m(xdata, ydata, entree, x, y, sortie)
!
!-------------------------------------------------------------------------------
! Author:  Z.X. Li (april 1st 1994) (see also A. Harzallah and L. Fairhead)
!-------------------------------------------------------------------------------
! Purpose: Naive method to regrid on a coarser grid.
!   Value at a new point is an average of the old points lcoated in a zone
!   surrounding the considered new point.
!   No ponderation is used (see grille_p)
!
!           (c)
!        ----d-----
!        | . . . .|
!        |        |
!     (b)a . * . .b(a)
!        |        |
!        | . . . .|
!        ----c-----
!           (d)
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: xdata(:), ydata(:)       !--- INPUT  FIELD X AND Y COORD.
  REAL, INTENT(IN)  :: entree(SIZE(xdata),SIZE(ydata)) !--- INPUT  FIELD
  REAL, INTENT(IN)  :: x(:), y(:)               !--- OUTPUT FIELD X AND Y COORD.
  REAL, INTENT(OUT) :: sortie(SIZE(x),SIZE(y))  !--- OUTPUT FIELD 
!-------------------------------------------------------------------------------
  CALL fine2coarse(xdata,ydata,x,y,sortie,entree)

END SUBROUTINE grille_m
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE rugosite(xdata, ydata, entree, x, y, sortie, mask)
!
!-------------------------------------------------------------------------------
! Author:  Z.X. Li (april 1st 1994)
!-------------------------------------------------------------------------------
! Purpose: Remap rugosity length ; constant value (0.001) on oceans.
! Naive method  (see grille_m)
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: xdata(:), ydata(:)      !--- INPUT  FIELD X AND Y COORD.
  REAL, INTENT(IN)  :: entree(SIZE(xdata),SIZE(ydata)) !--- INPUT  FIELD
  REAL, INTENT(IN)  :: x(:), y(:)              !--- OUTPUT FIELD X AND Y COORD.
  REAL, INTENT(OUT) :: sortie(SIZE(x),SIZE(y)) !--- OUTPUT FIELD 
  REAL, INTENT(IN)  :: mask  (SIZE(x),SIZE(y)) !--- MASK
!-------------------------------------------------------------------------------
  CALL fine2coarse(xdata,ydata,x,y,sortie,LOG(entree))
  WHERE(NINT(mask)==1)
    sortie(:,:)=EXP(sortie(:,:))
  ELSE WHERE
    sortie(:,:)=0.001
  END WHERE

END SUBROUTINE rugosite
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE sea_ice(xdata, ydata, glace01, x, y, frac_ice)
!
!-------------------------------------------------------------------------------
! Author:  Z.X. Li (april 1st 1994)
! Purpose: Regrid ice indicator (0 or 1) on a coarser grid to get an ice fract.
! field (between 0 and 1).
! Naive method  (see grille_m)
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: xdata(:), ydata(:)      !--- INPUT  FIELD X AND Y COORD.
  REAL, INTENT(IN)  :: glace01(:,:)            !--- INPUT  FIELD
  REAL, INTENT(IN)  :: x(:), y(:)              !--- OUTPUT FIELD X AND Y COORD.
  REAL, INTENT(OUT) :: frac_ice(SIZE(x),SIZE(y)) !--- OUTPUT FIELD 
!-------------------------------------------------------------------------------
  CALL fine2coarse(xdata,ydata,x,y,frac_ice,msk=NINT(glace01)==1)

END SUBROUTINE sea_ice
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!
SUBROUTINE rugsoro(xrel, yrel, relief, xmod, ymod, rugs)
!
!-------------------------------------------------------------------------------
! Author:  Z.X. Li (april 1st 1994) ; Modifications: august 23rd 1995.
! Purpose: Compute rugosity due to orography by using standard dev in a 1x1 cell
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
  REAL, INTENT(IN)  :: xrel(:), yrel(:)        !--- INPUT  FIELD X AND Y COORD.
  REAL, INTENT(IN)  :: relief(:,:)             !--- INPUT  FIELD
  REAL, INTENT(IN)  :: xmod(:), ymod(:)        !--- OUTPUT FIELD X AND Y COORD.
  REAL, INTENT(OUT) :: rugs(SIZE(xmod),SIZE(ymod)) !--- OUTPUT FIELD 
!-------------------------------------------------------------------------------
! Local variable:
  INTEGER           :: k, nn
  INTEGER, PARAMETER:: itmp=360, jtmp=180
  REAL  :: out(SIZE(xmod),SIZE(xmod)), amin, amax
  REAL  :: cham1tmp(itmp,jtmp), cham2tmp(itmp,jtmp), xtmp(itmp), ytmp(jtmp)
!-------------------------------------------------------------------------------

!--- TEST INPUT FILE FITS FOR ONE DEGREE RESOLUTION
  xtmp(:)=4.0*ATAN(1.0)*[(-1.0+REAL(2*k-1)/REAL(itmp),k=1,itmp)]
  ytmp(:)=2.0*ATAN(1.0)*[(-1.0+REAL(2*k-1)/REAL(jtmp),k=1,jtmp)]
  CALL fine2coarse(xrel,yrel,xtmp,ytmp,cham1tmp,relief,d_o2=cham2tmp)
  cham2tmp(:,:)=cham2tmp(:,:)-cham1tmp(:,:)**2
  nn=COUNT(cham2tmp<=-7.5)
  IF(nn/=0) THEN
    PRINT*,'Problem for rugsoro ; std**2 < -7.5 for several points: ',nn
    CALL ABORT_GCM("", "", 1)
  END IF
  nn=COUNT(cham2tmp<0)
  IF(nn/=0) PRINT*,'Problem for rugsoro ; std**2 < 0. for several points: ',nn
  WHERE(cham2tmp<0.0) cham2tmp=0.0
  cham2tmp(:,:)=SQRT(cham2tmp(:,:))
  amin=MINVAL(cham2tmp); amax=MAXVAL(cham2tmp)
  PRINT*, 'Ecart-type 1x1:', amin, amax

!--- COMPUTE RUGOSITY AT REQUIRED SCALE
  WHERE(cham2tmp<0.001) cham2tmp=0.001
  CALL fine2coarse(xtmp,ytmp,xmod,ymod,out,REAL(LOG(cham2tmp)))
  out=EXP(out)
  amin=MINVAL(out); amax=MAXVAL(out)
  PRINT*, 'Ecart-type du modele:', amin, amax
  out=out/amax*20.0
  amin=MINVAL(out); amax=MAXVAL(out)
  PRINT*, 'Longueur de rugosite du modele:', amin, amax
  rugs=REAL(out)

END SUBROUTINE rugsoro
!
!-------------------------------------------------------------------------------

END MODULE grid_atob_m
!
!*******************************************************************************

