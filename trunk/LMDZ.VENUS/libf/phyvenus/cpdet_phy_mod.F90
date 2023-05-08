module cpdet_phy_mod

! Cpdet module to use in physics (same as Cpdet module in dynamics)
! Warning: Changing hard coded formulaes in this module should also be
!          done in dynamics cpdet_mod module  (and vice-versa) !

implicit none

real,save :: cpp ! reference Cp
real,save :: nu_venus
real,save :: t0_venus

contains

      SUBROUTINE init_cpdet_phy(cpp_,nu_venus_,t0_venus_)
      ! initialize module variables
      REAL,INTENT(IN) :: cpp_ ! cpp from dynamics
      REAL,INTENT(IN) :: nu_venus_ ! nu_venus from dynamics
      REAL,INTENT(IN) :: t0_venus_ ! t0_venus from dynamics
      
      cpp=cpp_
      nu_venus=nu_venus_
      t0_venus=t0_venus_
      
      END SUBROUTINE init_cpdet_phy

!======================================================================

      FUNCTION cpdet(t)

      IMPLICIT none

! for cpp, nu_venus and t0_venus:

      real,intent(in) :: t
      real cpdet

      if (t0_venus.ne.0.) then
       cpdet = cpp*(t/t0_venus)**nu_venus
      else
       cpdet = cpp
      endif

      end function cpdet
      
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

      if (t0_venus.ne.0.) then
       yteta = yt**nu_venus                                          &
     &            - nu_venus * t0_venus**nu_venus * log(ypk/cpp)
       yteta = yteta**(1./nu_venus)
      else
       yteta = yt * cpp/ypk
      endif

!----------------------
! ATMOSPHERE PROFONDE
      yteta = yteta*ratio_mod
!----------------------

      end subroutine t2tpot

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

!----------------------
! ATMOSPHERE PROFONDE
!----------------------
      if (t0_venus.ne.0.) then
     yt = (yteta/ratio_mod)**nu_venus                               &
!----------------------
         + nu_venus * t0_venus**nu_venus * log(ypk/cpp)
     yt = yt**(1./nu_venus)
      else
     yt = (yteta/ratio_mod) * ypk/cpp
      endif
  
  
      end subroutine tpot2t

end module cpdet_phy_mod
