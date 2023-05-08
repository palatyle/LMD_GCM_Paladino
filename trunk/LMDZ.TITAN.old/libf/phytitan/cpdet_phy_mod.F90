module cpdet_phy_mod

implicit none

real,save :: cpp ! reference Cp

contains

      SUBROUTINE init_cpdet_phy(cpp_)
      ! initialize module variables
      REAL,INTENT(IN) :: cpp_ ! cpp from dynamics
      
      cpp=cpp_
      
      END SUBROUTINE init_cpdet_phy

!======================================================================

      FUNCTION cpdet(t)

      IMPLICIT none

! for now in Titan Cp does not change with temperature

      real,intent(in) :: t
      real cpdet

      cpdet = cpp
      
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
      
      yteta = yt * cpp/ypk

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
      
      yt = yteta * ypk/cpp

      end subroutine tpot2t

!======================================================================

end module cpdet_phy_mod
