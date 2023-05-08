! ADAPTATION GCM POUR CP(T)
!======================================================================
! S. Lebonnois, 10/2007:
!
! VENUS: Cp(T) = cpp*(T/T0)^nu 
! avec T0=460. et nu=0.35
!c cpp=RCPD=cp0 = 1000.
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

      FUNCTION cpdet(t)
      IMPLICIT none
#include "planet.h"

      real cpdet,t

      if(nu.ne.0.) then
        cpdet = cp0*(t/t0)**nu
      else
        cpdet = cp0
      endif

      return
      end
      
!======================================================================
!======================================================================

      SUBROUTINE t2tpot(npoints,yt, yteta, ypk)
      IMPLICIT none
!======================================================================
! Arguments:
!
! yt   --------input-R- Temperature
! yteta-------output-R- Temperature potentielle
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================
#include "planet.h"

      integer npoints
      REAL    yt(npoints), yteta(npoints), ypk(npoints)
      
      if(nu.ne.0.) then
        yteta = yt**nu - nu * t0**nu * log(ypk/cp0)
        yteta = yteta**(1./nu)
      else
        yteta = yt * cp0/ypk
      endif
        
      end

!======================================================================
!======================================================================

      SUBROUTINE tpot2t(npoints,yteta, yt, ypk)
      IMPLICIT none
!======================================================================
! Arguments:
!
! yteta--------input-R- Temperature potentielle
! yt   -------output-R- Temperature
! ypk  --------input-R- Fonction d'Exner: RCPD*(pplay/pref)**RKAPPA
!
!======================================================================
#include "planet.h"

      integer npoints
      REAL yt(npoints), yteta(npoints), ypk(npoints)

      if(nu.ne.0.) then
        yt = yteta**nu + nu * t0**nu * log(ypk/cp0)
        yt = yt**(1./nu)
      else
        yt = yteta * ypk/cp0
      endif
      
      end

