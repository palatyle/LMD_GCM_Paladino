MODULE cpdet_phy_mod
IMPLICIT NONE  

! ADAPTATION GCM POUR CP(T)
!======================================================================
! S. Lebonnois, 10/2007:
!
! VENUS: Cp(T) = cpp*(T/T0)^nu 
! avec T0=460. et nu=0.35
! cpp=RCPD=cp0 = 1000.
! R/RCPD = RKAPPA
!
! La fonction d'Exner reste pk = RCPD*(play/pref)**RKAPPA
! 
! T et teta (temperature potentielle) sont liees par:
! 
!   integrale[teta a T](cp/T dT) = integrale[pref a p](R/p dp)
!
! Dans le cas de l'expression pour Venus, ca donne:
!
!   teta**nu = T**nu - nu * T0**nu * ln[ (p/pref)**RKAPPA ]
! ou
!   teta**nu = T**nu - nu * T0**nu * ln[pk/RCPD]
!
! On passe de T a teta par t2tpot(t,teta,pk)
! On passe de teta a T par tpot2t(teta,t,pk)
!
! Pour DT <-> Dteta, on utilise: dteta = dT *(T/teta)**(nu-1)
! -> routine dt2dtpot(dt,dteta,t,teta) 
! (utilisee seulement pour le contregradient)
!
!======================================================================


      REAL cp0,t0,nu
!      PARAMETER (cp0 = 1000.) !doit etre egal a cpp (dyn) et RCPD (phy)
      PARAMETER (cp0 = 800.) !doit etre egal a cpp (dyn) et RCPD (phy)
!      PARAMETER (t0  = 460.)
      PARAMETER (t0  = 480.)
!      PARAMETER (t0  = 300.)
      PARAMETER (nu  = 0.35)
!      PARAMETER (nu  = 1.0)

CONTAINS

      FUNCTION cpdet(t)
      IMPLICIT NONE

      real cpdet,t

      cpdet = cp0*(t/t0)**nu

      RETURN
      END FUNCTION cpdet
      
!======================================================================
!======================================================================

      SUBROUTINE t2tpot(npoints,yt, yteta, ypk)
      IMPLICIT NONE
!======================================================================
! Arguments:
!
! yt   --------input-R- Temperature
! yteta-------output-R- Temperature potentielle
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      INTEGER npoints
      REAL    yt(npoints), yteta(npoints), ypk(npoints)
      
      yteta = yt**nu - nu * t0**nu * log(ypk/cp0)
      yteta = yteta**(1./nu)
        

      !yteta = yt/ypk


      RETURN
      END SUBROUTINE t2tpot

!======================================================================
!======================================================================

      SUBROUTINE tpot2t(npoints,yteta, yt, ypk)
      IMPLICIT NONE
!======================================================================
! Arguments:
!
! yteta--------input-R- Temperature potentielle
! yt   -------output-R- Temperature
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================

      INTEGER npoints
      REAL yt(npoints), yteta(npoints), ypk(npoints)

      yt = yteta**nu + nu * t0**nu * log(ypk/cp0)
      yt = yt**(1./nu)

      !yt = yteta*ypk

      RETURN
      END SUBROUTINE tpot2t

END MODULE cpdet_phy_mod
