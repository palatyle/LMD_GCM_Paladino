!
MODULE climb_wind_mod
!
! Module to solve the verctical diffusion of the wind components "u" and "v".
!
  USE dimphy

  IMPLICIT NONE

  SAVE
  PRIVATE
  
  REAL, DIMENSION(:),   ALLOCATABLE  :: alf1, alf2
  !$OMP THREADPRIVATE(alf1,alf2)
  REAL, DIMENSION(:,:), ALLOCATABLE  :: Kcoefm
  !$OMP THREADPRIVATE(Kcoefm)
  REAL, DIMENSION(:,:), ALLOCATABLE  :: Ccoef_U, Dcoef_U
  !$OMP THREADPRIVATE(Ccoef_U, Dcoef_U)
  REAL, DIMENSION(:,:), ALLOCATABLE  :: Ccoef_V, Dcoef_V
  !$OMP THREADPRIVATE(Ccoef_V, Dcoef_V)
  REAL, DIMENSION(:), ALLOCATABLE   :: Acoef_U, Bcoef_U
  !$OMP THREADPRIVATE(Acoef_U, Bcoef_U)
  REAL, DIMENSION(:), ALLOCATABLE   :: Acoef_V, Bcoef_V
  !$OMP THREADPRIVATE(Acoef_V, Bcoef_V)
  LOGICAL                            :: firstcall=.TRUE.
  !$OMP THREADPRIVATE(firstcall)

  
  PUBLIC :: climb_wind_down, climb_wind_up

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE climb_wind_init

    INTEGER             :: ierr
    CHARACTER(len = 20) :: modname = 'climb_wind_init'    

!****************************************************************************************
! Allocation of global module variables
!
!****************************************************************************************

    ALLOCATE(alf1(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocate alf2',1)

    ALLOCATE(alf2(klon), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocate alf2',1)

    ALLOCATE(Kcoefm(klon,klev), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocate Kcoefm',1)

    ALLOCATE(Ccoef_U(klon,klev), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocate Ccoef_U',1)

    ALLOCATE(Dcoef_U(klon,klev), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocation Dcoef_U',1)

    ALLOCATE(Ccoef_V(klon,klev), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocation Ccoef_V',1)

    ALLOCATE(Dcoef_V(klon,klev), stat=ierr)
    IF (ierr /= 0) CALL abort_gcm(modname,'Pb in allocation Dcoef_V',1)

    ALLOCATE(Acoef_U(klon), Bcoef_U(klon), Acoef_V(klon), Bcoef_V(klon), STAT=ierr)
    IF ( ierr /= 0 )  PRINT*,' pb in allloc Acoef_U and Bcoef_U, ierr=', ierr

    firstcall=.FALSE.

  END SUBROUTINE climb_wind_init
!
!****************************************************************************************
!
  SUBROUTINE climb_wind_down(knon, dtime, coef_in, pplay, paprs, temp, delp, u_old, v_old, &
       Acoef_U_out, Acoef_V_out, Bcoef_U_out, Bcoef_V_out)
!
! This routine calculates for the wind components u and v,
! recursivly the coefficients C and D in equation 
! X(k) = C(k) + D(k)*X(k-1), X=[u,v], k=[1,klev] is the vertical layer.
!
!
    INCLUDE "YOMCST.h"
! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: knon
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: coef_in
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay ! pres au milieu de couche (Pa)
    REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs ! pression a inter-couche (Pa)
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: temp  ! temperature
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: delp
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: u_old
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: v_old

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef_U_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef_V_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Bcoef_U_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Bcoef_V_out

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon)                    :: u1lay, v1lay
    INTEGER                                  :: k, i


!****************************************************************************************
! Initialize module
    IF (firstcall) CALL climb_wind_init

!****************************************************************************************
! Calculate the coefficients C and D in : u(k) = C(k) + D(k)*u(k-1)
!
!****************************************************************************************
! - Define alpha (alf1 and alf2) 
    alf1(:) = 1.0
    alf2(:) = 1.0 - alf1(:)

! - Calculate the coefficients K
    Kcoefm(:,:) = 0.0
    DO k = 2, klev
       DO i=1,knon
          Kcoefm(i,k) = coef_in(i,k)*RG*RG*dtime/(pplay(i,k-1)-pplay(i,k)) &
               *(paprs(i,k)*2/(temp(i,k)+temp(i,k-1))/RD)**2
       END DO
    END DO

! - Calculate the coefficients C and D, component "u"
    CALL calc_coef(knon, Kcoefm(:,:), delp(:,:), &
         u_old(:,:), alf1(:), alf2(:),  &
         Ccoef_U(:,:), Dcoef_U(:,:), Acoef_U(:), Bcoef_U(:))

! - Calculate the coefficients C and D, component "v"
    CALL calc_coef(knon, Kcoefm(:,:), delp(:,:), &
         v_old(:,:), alf1(:), alf2(:),  &
         Ccoef_V(:,:), Dcoef_V(:,:), Acoef_V(:), Bcoef_V(:))

!****************************************************************************************
! 6)
! Return the first layer in output variables
!
!****************************************************************************************
    Acoef_U_out = Acoef_U
    Bcoef_U_out = Bcoef_U
    Acoef_V_out = Acoef_V
    Bcoef_V_out = Bcoef_V

  END SUBROUTINE climb_wind_down
!
!****************************************************************************************
!
  SUBROUTINE calc_coef(knon, Kcoef, delp, X, alfa1, alfa2, Ccoef, Dcoef, Acoef, Bcoef)
!
! Find the coefficients C and D in fonction of alfa, K and delp
!
! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: knon
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: Kcoef, delp
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: X
    REAL, DIMENSION(klon), INTENT(IN)        :: alfa1, alfa2

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef, Bcoef
    REAL, DIMENSION(klon,klev), INTENT(OUT)  :: Ccoef, Dcoef
  
! local variables
!****************************************************************************************
    INTEGER                                  :: k, i
    REAL                                     :: buf

    INCLUDE "YOMCST.h"
!****************************************************************************************
! 

! Calculate coefficients C and D at top level, k=klev
!
    Ccoef(:,:) = 0.0
    Dcoef(:,:) = 0.0

    DO i = 1, knon
       buf = delp(i,klev) + Kcoef(i,klev)

       Ccoef(i,klev) = X(i,klev)*delp(i,klev)/buf 
       Dcoef(i,klev) = Kcoef(i,klev)/buf
    END DO
    
! 
! Calculate coefficients C and D at top level (klev-1) <= k <= 2
!
    DO k=(klev-1),2,-1
       DO i = 1, knon
          buf = delp(i,k) + Kcoef(i,k) + Kcoef(i,k+1)*(1.-Dcoef(i,k+1))
          
          Ccoef(i,k) = (X(i,k)*delp(i,k) + Kcoef(i,k+1)*Ccoef(i,k+1))/buf
          Dcoef(i,k) = Kcoef(i,k)/buf
       END DO
    END DO

!
! Calculate coeffiecent A and B at surface
!
    DO i = 1, knon
       buf = delp(i,1) + Kcoef(i,2)*(1-Dcoef(i,2))
       Acoef(i) = (X(i,1)*delp(i,1) + Kcoef(i,2)*Ccoef(i,2))/buf
       Bcoef(i) = -RG/buf
    END DO

  END SUBROUTINE calc_coef
!
!****************************************************************************************
!

  SUBROUTINE climb_wind_up(knon, dtime, u_old, v_old, flx_u1, flx_v1,  &
       flx_u_new, flx_v_new, d_u_new, d_v_new)
!
! Diffuse the wind components from the surface layer and up to the top layer. 
! Coefficents A, B, C and D are known from before. Start values for the diffusion are the
! momentum fluxes at surface.
!
! u(k=1) = A + B*flx*dtime
! u(k)   = C(k) + D(k)*u(k-1)  [2 <= k <= klev]
!
!****************************************************************************************
    INCLUDE "YOMCST.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon,klev), INTENT(IN)  :: u_old
    REAL, DIMENSION(klon,klev), INTENT(IN)  :: v_old
    REAL, DIMENSION(klon), INTENT(IN)       :: flx_u1, flx_v1 ! momentum flux

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon,klev), INTENT(OUT) :: flx_u_new, flx_v_new
    REAL, DIMENSION(klon,klev), INTENT(OUT) :: d_u_new, d_v_new

! Local variables
!****************************************************************************************
    REAL, DIMENSION(klon,klev)              :: u_new, v_new
    INTEGER                                 :: k, i
    
!
!****************************************************************************************

! Niveau 1
    DO i = 1, knon
       u_new(i,1) = Acoef_U(i) + Bcoef_U(i)*flx_u1(i)*dtime
       v_new(i,1) = Acoef_V(i) + Bcoef_V(i)*flx_v1(i)*dtime
    END DO

! Niveau 2 jusqu'au sommet klev
    DO k = 2, klev
       DO i=1, knon
          u_new(i,k) = Ccoef_U(i,k) + Dcoef_U(i,k) * u_new(i,k-1)
          v_new(i,k) = Ccoef_V(i,k) + Dcoef_V(i,k) * v_new(i,k-1)
       END DO
    END DO

!****************************************************************************************
! Calcul flux
!
!== flux_u/v est le flux de moment angulaire (positif vers bas)
!== dont l'unite est: (kg m/s)/(m**2 s) 
!
!****************************************************************************************
!
    flx_u_new(:,:) = 0.0
    flx_v_new(:,:) = 0.0

    flx_u_new(1:knon,1)=flx_u1(1:knon)
    flx_v_new(1:knon,1)=flx_v1(1:knon)

! Niveau 2->klev
    DO k = 2, klev
       DO i = 1, knon
          flx_u_new(i,k) = Kcoefm(i,k)/RG/dtime * &
               (u_new(i,k)-u_new(i,k-1))
          
          flx_v_new(i,k) = Kcoefm(i,k)/RG/dtime * &
               (v_new(i,k)-v_new(i,k-1))
       END DO
    END DO

!****************************************************************************************
! Calcul tendances
!
!****************************************************************************************
    d_u_new(:,:) = 0.0
    d_v_new(:,:) = 0.0
    DO k = 1, klev
       DO i = 1, knon
          d_u_new(i,k) = u_new(i,k) - u_old(i,k)
          d_v_new(i,k) = v_new(i,k) - v_old(i,k)
       END DO
    END DO

  END SUBROUTINE climb_wind_up
!
!****************************************************************************************
!
END MODULE climb_wind_mod
