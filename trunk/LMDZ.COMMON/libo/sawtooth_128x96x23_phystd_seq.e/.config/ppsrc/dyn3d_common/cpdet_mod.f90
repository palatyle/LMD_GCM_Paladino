










module cpdet_mod

implicit none

! ADAPTATION OF GCM TO CP(T)
!======================================================================
! S. Lebonnois, 10/2010
!
! Cp must be computed using cpdet(t) to be valid
!
! The Exner function is still pk = RCPD*(play/pref)**RKAPPA
! (RCPD=cpp, RKAPPA=kappa)
!
! One goes from T to teta (potential temperature) using t2tpot(t,teta,pk)
! One goes from teta to T using tpot2t(teta,t,pk)
!
!======================================================================

contains

      SUBROUTINE ini_cpdet
      
      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: nu_venus,t0_venus
      IMPLICIT none
!======================================================================
! Initialization of nu_venus and t0_venus
!======================================================================

      if (cpofT) then
          nu_venus=0.35
          t0_venus=460.
      else
          nu_venus=0.
          t0_venus=0.
      endif

      return
      end subroutine ini_cpdet

!======================================================================
!======================================================================

      FUNCTION cpdet(t)

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,t0_venus,nu_venus
      IMPLICIT none

! for cpp, nu_venus and t0_venus:

      real,intent(in) :: t
      real cpdet

      if (cpofT) then
          cpdet = cpp*(t/t0_venus)**nu_venus
      else
          cpdet = cpp
      endif

      return
      end function cpdet
      
!======================================================================
!======================================================================

      SUBROUTINE t2tpot(npoints, yt, yteta, ypk)
!======================================================================
! Arguments:
!
! yt   --------input-R- Temperature
! yteta-------output-R- Temperature potentielle
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,t0_venus,nu_venus

      IMPLICIT NONE

      integer,intent(in) :: npoints
      REAL,intent(in) :: yt(npoints), ypk(npoints)
      REAL,intent(out) :: yteta(npoints)
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(npoints)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44
      
      DO k = 1, npoints
        ratio_mod(k) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k).gt.ypklim)) then
	   ratio_mod(k) = mmm0 /                                        &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k)/cpp))                    &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k) = max(ratio_mod(k),mmm0/(mmm0+0.56))
	endif
      ENDDO
!----------------------

      if (cpofT) then
          yteta = yt**nu_venus                                          &
     &            - nu_venus * t0_venus**nu_venus * log(ypk/cpp)
          yteta = yteta**(1./nu_venus)
!----------------------
! ATMOSPHERE PROFONDE
          yteta = yteta*ratio_mod
!----------------------

      else
          yteta = yt * cpp/ypk
      endif

      return
      end subroutine t2tpot

!======================================================================
!======================================================================

      SUBROUTINE tpot2t(npoints,yteta, yt, ypk)
!======================================================================
! Arguments:
!
! yteta--------input-R- Temperature potentielle
! yt   -------output-R- Temperature
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      USE control_mod, ONLY: cpofT 
      USE comconst_mod, ONLY: cpp,nu_venus,t0_venus

      IMPLICIT NONE

      integer,intent(in) :: npoints
      REAL,intent(in) :: yteta(npoints), ypk(npoints)
      REAL,intent(out) :: yt(npoints)
      
!----------------------
! ATMOSPHERE PROFONDE
      real ypklim,ypklim2,mmm0
      real ratio_mod(npoints)
      integer k
! theta_mod = theta * ratio_mod => ratio_mod = mmm0/mmm(p)

      ypklim  = cpp*(  6e6/9.2e6)**0.1914
      ypklim2 = cpp*(8.9e6/9.2e6)**0.1914
      mmm0 = 43.44

      DO k = 1, npoints
        ratio_mod(k) = 1.
! ATM PROFONDE DESACTIVEE !
        if (   (1 .EQ. 0)  .AND.(ypk(k).gt.ypklim)) then
	   ratio_mod(k) = mmm0 /                                        &
     &  (mmm0+0.56*(log(ypklim/cpp)-log(ypk(k)/cpp))                    &
     &            /(log(ypklim/cpp)-log(ypklim2/cpp)))
           ratio_mod(k) = max(ratio_mod(k),mmm0/(mmm0+0.56))
	endif
      ENDDO
!----------------------

      if (cpofT) then

!----------------------
! ATMOSPHERE PROFONDE
!----------------------
         yt = (yteta/ratio_mod)**nu_venus                               &
!----------------------
     &       + nu_venus * t0_venus**nu_venus * log(ypk/cpp)
         yt = yt**(1./nu_venus)
      else
          yt = yteta * ypk/cpp
      endif
  
      return
      end subroutine tpot2t

!======================================================================
!======================================================================
! Routines pour les calculs paralleles
!======================================================================
!======================================================================
! Fin routines specifiques parallele
!======================================================================
!======================================================================
!
! ATTENTION
!
! Si un jour on a besoin, il faudra coder les routines 
!    dt2dtpot / dtpto2dt 
!
!======================================================================
!======================================================================
end module cpdet_mod
