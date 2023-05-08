MODULE climb_hq_mod
!
! Module to solve the verctical diffusion of "q" and "H"; 
! specific humidity and potential energi.
!
  USE dimphy

  IMPLICIT NONE
  SAVE 
  PRIVATE
  PUBLIC :: climb_hq_down, climb_hq_up

  REAL, DIMENSION(:,:), ALLOCATABLE :: gamaq, gamah
  !$OMP THREADPRIVATE(gamaq,gamah)
  REAL, DIMENSION(:,:), ALLOCATABLE :: Ccoef_Q, Dcoef_Q
  !$OMP THREADPRIVATE(Ccoef_Q, Dcoef_Q)
  REAL, DIMENSION(:,:), ALLOCATABLE :: Ccoef_H, Dcoef_H
  !$OMP THREADPRIVATE(Ccoef_H, Dcoef_H)
  REAL, DIMENSION(:), ALLOCATABLE   :: Acoef_Q, Bcoef_Q
  !$OMP THREADPRIVATE(Acoef_Q, Bcoef_Q)
  REAL, DIMENSION(:), ALLOCATABLE   :: Acoef_H, Bcoef_H
  !$OMP THREADPRIVATE(Acoef_H, Bcoef_H)
  REAL, DIMENSION(:,:), ALLOCATABLE :: Kcoefhq
  !$OMP THREADPRIVATE(Kcoefhq)

CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE climb_hq_down(knon, coefhq, paprs, pplay, &
       delp, temp, q, dtime, &
       Acoef_H_out, Acoef_Q_out, Bcoef_H_out, Bcoef_Q_out)

    INCLUDE "YOMCST.h"
! This routine calculates recursivly the coefficients C and D
! for the quantity X=[Q,H] in equation X(k) = C(k) + D(k)*X(k-1), where k is
! the index of the vertical layer.
!
! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: knon
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: coefhq
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay 
    REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs 
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: temp, delp  ! temperature
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: q
    REAL, INTENT(IN)                         :: dtime

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef_H_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef_Q_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Bcoef_H_out
    REAL, DIMENSION(klon), INTENT(OUT)       :: Bcoef_Q_out

! Local variables
!****************************************************************************************
    LOGICAL, SAVE                            :: first=.TRUE.
    !$OMP THREADPRIVATE(first)
    REAL, DIMENSION(klon,klev)               :: local_H
    REAL, DIMENSION(klon)                    :: psref 
    REAL                                     :: delz, pkh
    INTEGER                                  :: k, i, ierr

! Include
!****************************************************************************************
    INCLUDE "compbl.h"    


!****************************************************************************************
! 1)
! Allocation at first time step only
!   
!****************************************************************************************

    IF (first) THEN
       first=.FALSE.
       ALLOCATE(Ccoef_Q(klon,klev), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Ccoef_Q, ierr=', ierr
       
       ALLOCATE(Dcoef_Q(klon,klev), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Dcoef_Q, ierr=', ierr
       
       ALLOCATE(Ccoef_H(klon,klev), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Ccoef_H, ierr=', ierr
       
       ALLOCATE(Dcoef_H(klon,klev), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Dcoef_H, ierr=', ierr
       
       ALLOCATE(Acoef_Q(klon), Bcoef_Q(klon), Acoef_H(klon), Bcoef_H(klon), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Acoef_X and Bcoef_X, ierr=', ierr
       
       ALLOCATE(Kcoefhq(klon,klev), STAT=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in allloc Kcoefhq, ierr=', ierr
       
       ALLOCATE(gamaq(1:klon,2:klev), STAT=ierr)
       IF ( ierr /= 0 ) PRINT*,' pb in allloc gamaq, ierr=', ierr
       
       ALLOCATE(gamah(1:klon,2:klev), STAT=ierr)
       IF ( ierr /= 0 ) PRINT*,' pb in allloc gamah, ierr=', ierr
    END IF

!****************************************************************************************
! 2)
! Definition of the coeficient K 
!
!****************************************************************************************
    Kcoefhq(:,:) = 0.0
    DO k = 2, klev
       DO i = 1, knon
          Kcoefhq(i,k) = &
               coefhq(i,k)*RG*RG*dtime /(pplay(i,k-1)-pplay(i,k)) &
               *(paprs(i,k)*2/(temp(i,k)+temp(i,k-1))/RD)**2
       ENDDO
    ENDDO

!****************************************************************************************
! 3)
! Calculation of gama for "Q" and "H"
!
!****************************************************************************************
!   surface pressure is used as reference
    psref(:) = paprs(:,1) 

!   definition of gama
    IF (iflag_pbl == 1) THEN
       gamaq(:,:) = 0.0
       gamah(:,:) = -1.0e-03
       gamah(:,2) = -2.5e-03
 
! conversion de gama
       DO k = 2, klev
          DO i = 1, knon
             delz = RD * (temp(i,k-1)+temp(i,k)) / & 
                    2.0 / RG / paprs(i,k) * (pplay(i,k-1)-pplay(i,k))
             pkh  = (psref(i)/paprs(i,k))**RKAPPA
          
! convertie gradient verticale d'humidite specifique en difference d'humidite specifique entre centre de couches
             gamaq(i,k) = gamaq(i,k) * delz    
! convertie gradient verticale de temperature en difference de temperature potentielle entre centre de couches 
             gamah(i,k) = gamah(i,k) * delz * RCPD * pkh
          ENDDO
       ENDDO

    ELSE
       gamaq(:,:) = 0.0
       gamah(:,:) = 0.0
    ENDIF
    

!****************************************************************************************    
! 4)
! Calculte the coefficients C and D for specific humidity, q
!
!****************************************************************************************
    
    CALL calc_coef(knon, Kcoefhq(:,:), gamaq(:,:), delp(:,:), q(:,:), &
         Ccoef_Q(:,:), Dcoef_Q(:,:), Acoef_Q, Bcoef_Q)

!****************************************************************************************
! 5)
! Calculte the coefficients C and D for potentiel entalpie, H 
!
!****************************************************************************************
    local_H(:,:) = 0.0

    DO k=1,klev
       DO i = 1, knon
          ! convertie la temperature en entalpie potentielle
          local_H(i,k) = RCPD * temp(i,k) * &
               (psref(i)/pplay(i,k))**RKAPPA
       ENDDO
    ENDDO

    CALL calc_coef(knon, Kcoefhq(:,:), gamah(:,:), delp(:,:), local_H(:,:), &
         Ccoef_H(:,:), Dcoef_H(:,:), Acoef_H, Bcoef_H)
 
!****************************************************************************************
! 6)
! Return the first layer in output variables
!
!****************************************************************************************
    Acoef_H_out = Acoef_H
    Bcoef_H_out = Bcoef_H
    Acoef_Q_out = Acoef_Q
    Bcoef_Q_out = Bcoef_Q

  END SUBROUTINE climb_hq_down
!
!****************************************************************************************
!
  SUBROUTINE calc_coef(knon, Kcoef, gama, delp, X, Ccoef, Dcoef, Acoef, Bcoef)
!
! Calculate the coefficients C and D in : X(k) = C(k) + D(k)*X(k-1)
! where X is H or Q, and k the vertical level k=1,klev
!
    INCLUDE "YOMCST.h"
! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: knon
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: Kcoef, delp
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: X
    REAL, DIMENSION(klon,2:klev), INTENT(IN) :: gama

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)       :: Acoef, Bcoef
    REAL, DIMENSION(klon,klev), INTENT(OUT)  :: Ccoef, Dcoef

! Local variables
!****************************************************************************************
    INTEGER                                  :: k, i
    REAL                                     :: buf

!****************************************************************************************
! Niveau au sommet, k=klev
!
!****************************************************************************************
    Ccoef(:,:) = 0.0
    Dcoef(:,:) = 0.0

    DO i = 1, knon
       buf = delp(i,klev) + Kcoef(i,klev)
       
       Ccoef(i,klev) = (X(i,klev)*delp(i,klev) - Kcoef(i,klev)*gama(i,klev))/buf
       Dcoef(i,klev) = Kcoef(i,klev)/buf
    END DO


!****************************************************************************************
! Niveau  (klev-1) <= k <= 2
!
!****************************************************************************************

    DO k=(klev-1),2,-1
       DO i = 1, knon
          buf = delp(i,k) + Kcoef(i,k) + Kcoef(i,k+1)*(1.-Dcoef(i,k+1))
          Ccoef(i,k) = (X(i,k)*delp(i,k) + Kcoef(i,k+1)*Ccoef(i,k+1) + &
               Kcoef(i,k+1)*gama(i,k+1) - Kcoef(i,k)*gama(i,k))/buf
          Dcoef(i,k) = Kcoef(i,k)/buf
       END DO
    END DO

!****************************************************************************************
! Niveau k=1
!
!****************************************************************************************

    DO i = 1, knon
       buf = delp(i,1) + Kcoef(i,2)*(1.-Dcoef(i,2))
       Acoef(i) = (X(i,1)*delp(i,1) + Kcoef(i,2)*(gama(i,2)+Ccoef(i,2)))/buf
       Bcoef(i) = -1. * RG / buf
    END DO

  END SUBROUTINE calc_coef
!
!****************************************************************************************
!
  SUBROUTINE climb_hq_up(knon, dtime, t_old, q_old, &
       flx_q1, flx_h1, paprs, pplay, &
       flux_q, flux_h, d_q, d_t)
! 
! This routine calculates the flux and tendency of the specific humidity q and 
! the potential engergi H. 
! The quantities q and H are calculated according to 
! X(k) = C(k) + D(k)*X(k-1) for X=[q,H], where the coefficients 
! C and D are known from before and k is index of the vertical layer.
!   
    INCLUDE "YOMCST.h"
! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                      :: knon
    REAL, INTENT(IN)                         :: dtime
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: t_old, q_old
    REAL, DIMENSION(klon), INTENT(IN)        :: flx_q1, flx_h1
    REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay

! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon,klev), INTENT(OUT)  :: flux_q, flux_h, d_q, d_t

! Local variables
!****************************************************************************************
    LOGICAL, SAVE                            :: last=.FALSE.
    REAL, DIMENSION(klon,klev)               :: h_new, q_new
    REAL, DIMENSION(klon)                    :: psref         
    INTEGER                                  :: k, i, ierr

!****************************************************************************************
! 1) 
! Definition of some variables
!
!****************************************************************************************
    flux_q(:,:) = 0.0
    flux_h(:,:) = 0.0
    d_q(:,:)    = 0.0
    d_t(:,:)    = 0.0

    psref(1:knon) = paprs(1:knon,1)  

!****************************************************************************************
! 2)
! Calculation of Q and H
!
!****************************************************************************************

!- First layer
    q_new(1:knon,1) = Acoef_Q(1:knon) + Bcoef_Q(1:knon)*flx_q1(1:knon)*dtime
    h_new(1:knon,1) = Acoef_H(1:knon) + Bcoef_H(1:knon)*flx_h1(1:knon)*dtime
    
!- All the other layers 
    DO k = 2, klev
       DO i = 1, knon
          q_new(i,k) = Ccoef_Q(i,k) + Dcoef_Q(i,k)*q_new(i,k-1)
          h_new(i,k) = Ccoef_H(i,k) + Dcoef_H(i,k)*h_new(i,k-1)
       END DO
    END DO
!****************************************************************************************
! 3)
! Calculation of the flux for Q and H
!
!****************************************************************************************

!- The flux at first layer, k=1
    flux_q(1:knon,1)=flx_q1(1:knon)
    flux_h(1:knon,1)=flx_h1(1:knon)

!- The flux at all layers above surface
    DO k = 2, klev
       DO i = 1, knon
          flux_q(i,k) = (Kcoefhq(i,k)/RG/dtime) * &
               (q_new(i,k)-q_new(i,k-1)+gamaq(i,k))

          flux_h(i,k) = (Kcoefhq(i,k)/RG/dtime) * &
               (h_new(i,k)-h_new(i,k-1)+gamah(i,k)) 
       END DO
    END DO

!****************************************************************************************
! 4)
! Calculation of tendency for Q and H
!
!****************************************************************************************

    DO k = 1, klev
       DO i = 1, knon
          d_t(i,k) = h_new(i,k)/(psref(i)/pplay(i,k))**RKAPPA/RCPD - t_old(i,k)
          d_q(i,k) = q_new(i,k) - q_old(i,k)
       END DO
    END DO

!****************************************************************************************
! Some deallocations
!
!****************************************************************************************
    IF (last) THEN
       DEALLOCATE(Ccoef_Q, Dcoef_Q, Ccoef_H, Dcoef_H,stat=ierr)    
       IF ( ierr /= 0 )  PRINT*,' pb in dealllocate Ccoef_Q, Dcoef_Q, Ccoef_H, Dcoef_H, ierr=', ierr
       DEALLOCATE(Acoef_Q, Bcoef_Q, Acoef_H, Bcoef_H,stat=ierr)    
       IF ( ierr /= 0 )  PRINT*,' pb in dealllocate Acoef_Q, Bcoef_Q, Acoef_H, Bcoef_H, ierr=', ierr
       DEALLOCATE(gamaq, gamah,stat=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in dealllocate gamaq, gamah, ierr=', ierr
       DEALLOCATE(Kcoefhq,stat=ierr)
       IF ( ierr /= 0 )  PRINT*,' pb in dealllocate Kcoefhq, ierr=', ierr
    END IF
  END SUBROUTINE climb_hq_up
!
!****************************************************************************************
!
END MODULE climb_hq_mod

 




