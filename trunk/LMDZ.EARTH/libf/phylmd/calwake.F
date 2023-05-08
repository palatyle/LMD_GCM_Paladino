!
! $Id: calwake.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE CALWAKE(paprs,pplay,dtime
     :             ,t,q,omgb
     :             ,dt_dwn,dq_dwn,M_dwn,M_up
     :             ,dt_a,dq_a,sigd
     :             ,wdt_PBL,wdq_PBL
     :             ,udt_PBL,udq_PBL
     o             ,wake_deltat,wake_deltaq,wake_dth
     o             ,wake_h,wake_s,wake_dens
     o             ,wake_pe,wake_fip,wake_gfl
     o             ,dt_wake,dq_wake
     o             ,wake_k
     o             ,undi_t,undi_q
     o             ,wake_omgbdth,wake_dp_omgb
     o             ,wake_dtKE,wake_dqKE
     o             ,wake_dtPBL,wake_dqPBL
     o             ,wake_omg,wake_dp_deltomg
     o             ,wake_spread,wake_Cstar,wake_d_deltat_gw
     o             ,wake_ddeltat,wake_ddeltaq)
***************************************************************
*                                                             *
* CALWAKE                                                     *
*           interface avec le schema de calcul de la poche    *
*           froide                                            *
*                                                             *
* written by   : CHERUY Frederique, 13/03/2000, 10.31.05      *
* modified by :  ROEHRIG Romain,    01/30/2007                *
***************************************************************
*
      USE dimphy
      IMPLICIT none
c======================================================================
#include "dimensions.h"
!#include "dimphy.h"
#include "YOMCST.h"

c Arguments
c----------

      INTEGER  i,l,ktopw(klon)
      REAL   dtime

      REAL paprs(klon,klev+1),pplay(klon,klev)
      REAL t(klon,klev), q(klon,klev), omgb(klon,klev)
      REAL dt_dwn(klon,klev), dq_dwn(klon,klev),M_dwn(klon,klev)
      REAL M_up(klon,klev)
      REAL dt_a(klon,klev), dq_a(klon,klev)
      REAL wdt_PBL(klon,klev), wdq_PBL(klon,klev)
      REAL udt_PBL(klon,klev), udq_PBL(klon,klev)
      REAL wake_deltat(klon,klev),wake_deltaq(klon,klev)
      REAL dt_wake(klon,klev),dq_wake(klon,klev)
      REAL wake_d_deltat_gw(klon,klev)
      REAL wake_h(klon),wake_s(klon)
      REAL wake_dth(klon,klev)
      REAL wake_pe(klon),wake_fip(klon),wake_gfl(klon)
      REAL undi_t(klon,klev),undi_q(klon,klev)
      REAL wake_omgbdth(klon,klev),wake_dp_omgb(klon,klev)
      REAL wake_dtKE(klon,klev),wake_dqKE(klon,klev)
      REAL wake_dtPBL(klon,klev),wake_dqPBL(klon,klev)
      REAL wake_omg(klon,klev+1),wake_dp_deltomg(klon,klev)
      REAL wake_spread(klon,klev),wake_Cstar(klon)
      REAL wake_ddeltat(klon,klev),wake_ddeltaq(klon,klev)
      REAL d_deltatw(klon,klev), d_deltaqw(klon,klev)
      INTEGER wake_k(klon)
      REAL sigd(klon)
      REAL wake_dens(klon)

C  Variable internes
C  -----------------

      REAL aire
      REAL p(klon,klev),ph(klon,klev+1),pi(klon,klev)
      REAL te(klon,klev),qe(klon,klev),omgbe(klon,klev+1)
      REAL dtdwn(klon,klev),dqdwn(klon,klev)
      REAL dta(klon,klev),dqa(klon,klev)
      REAL wdtPBL(klon,klev),wdqPBL(klon,klev)
      REAL udtPBL(klon,klev),udqPBL(klon,klev) 
      REAL amdwn(klon,klev),amup(klon,klev)
      REAL dtw(klon,klev),dqw(klon,klev),dth(klon,klev)
      REAL d_deltat_gw(klon,klev)
      REAL dtls(klon,klev),dqls(klon,klev)
      REAL tu(klon,klev),qu(klon,klev)
      REAL hw(klon),sigmaw(klon),wape(klon),fip(klon),gfl(klon)
      REAL omgbdth(klon,klev),dp_omgb(klon,klev)
      REAL dtKE(klon,klev),dqKE(klon,klev)
      REAL dtPBL(klon,klev),dqPBL(klon,klev)
      REAL omg(klon,klev+1),dp_deltomg(klon,klev),spread(klon,klev)
      REAL Cstar(klon)
      REAL sigd0(klon),wdens(klon)

      REAL RDCP

c      print *, '-> calwake, wake_s ', wake_s(1)

      RDCP=1./3.5


c-----------------------------------------------------------
cIM 290108     DO 999 i=1,klon   ! a vectoriser
c----------------------------------------------------------


      DO l=1,klev
      DO i=1,klon 
        p(i,l)= pplay(i,l)
        ph(i,l)= paprs(i,l)
        pi(i,l) = (pplay(i,l)/100000.)**RDCP

        te(i,l) = t(i,l)
        qe(i,l) = q(i,l)
        omgbe(i,l) = omgb(i,l)

        dtdwn(i,l)= dt_dwn(i,l)
        dqdwn(i,l)= dq_dwn(i,l)
        dta(i,l)= dt_a(i,l)
        dqa(i,l)= dq_a(i,l)
        wdtPBL(i,l)= wdt_PBL(i,l)
        wdqPBL(i,l)= wdq_PBL(i,l)
        udtPBL(i,l)= udt_PBL(i,l)
        udqPBL(i,l)= udq_PBL(i,l)
      ENDDO
      ENDDO

      omgbe(:,klev+1) = 0.
      
      DO i=1,klon 
      sigd0(i)=sigd(i)
      ENDDO
c      print*, 'sigd0,sigd', sigd0, sigd(i)
      DO i=1,klon 
      ph(i,klev+1)=0.
      ENDDO

      DO i=1,klon 
      ktopw(i) = wake_k(i)
      ENDDO

      DO l=1,klev
      DO i=1,klon 
        dtw(i,l) = wake_deltat(i,l)
        dqw(i,l) = wake_deltaq(i,l)
      ENDDO
      ENDDO

      DO l=1,klev
      DO i=1,klon 
        dtls(i,l)=dt_wake(i,l)
        dqls(i,l)=dq_wake(i,l)
      ENDDO
      ENDDO

      DO i=1,klon 
      hw(i) = wake_h(i)
      sigmaw(i)= wake_s(i)
      ENDDO

cfkc les flux de masses sont evalues aux niveaux et valent 0 a la surface
cfkc  on veut le flux de masse au milieu des couches

      DO l=1,klev-1
      DO i=1,klon 
        amdwn(i,l)= 0.5*(M_dwn(i,l)+M_dwn(i,l+1))
        amdwn(i,l)= (M_dwn(i,l+1))
      ENDDO
      ENDDO

c au sommet le flux de masse est nul

      DO i=1,klon 
      amdwn(i,klev)=0.5*M_dwn(i,klev)
      ENDDO
c
      DO l = 1,klev
      DO i=1,klon 
        amup(i,l)=M_up(i,l)
      ENDDO
      ENDDO

      call WAKE(p,ph,pi,dtime,sigd0
     $                ,te,qe,omgbe
     $                ,dtdwn,dqdwn,amdwn,amup,dta,dqa
     $                ,wdtPBL,wdqPBL,udtPBL,udqPBL
     $                ,dtw,dqw,dth,hw,sigmaw,wape,fip,gfl
     $                ,dtls,dqls,ktopw
     $                ,omgbdth,dp_omgb,wdens
     $                ,tu,qu
     $                ,dtKE,dqKE
     $                ,dtPBL,dqPBL
     $                ,omg,dp_deltomg,spread
     $                ,Cstar,d_deltat_gw
     $                ,d_deltatw,d_deltaqw)
c
      DO l=1,klev
       DO i=1,klon 
        IF (ktopw(i) .GT. 0) THEN
           wake_deltat(i,l)= dtw(i,l)
           wake_deltaq(i,l)= dqw(i,l)
           wake_d_deltat_gw(i,l)= d_deltat_gw(i,l)
           wake_omgbdth(i,l) = omgbdth(i,l)
           wake_dp_omgb(i,l) = dp_omgb(i,l)
           wake_dtKE(i,l) = dtKE(i,l)
           wake_dqKE(i,l) = dqKE(i,l)
 	   wake_dtPBL(i,l) = dtPBL(i,l)
	   wake_dqPBL(i,l) = dqPBL(i,l)
           wake_omg(i,l) = omg(i,l)
           wake_dp_deltomg(i,l) = dp_deltomg(i,l)
           wake_spread(i,l) = spread(i,l)
           wake_dth(i,l) = dth(i,l)
           dt_wake(i,l) = dtls(i,l)
           dq_wake(i,l) = dqls(i,l)
           undi_t(i,l) = tu(i,l)
           undi_q(i,l) = qu(i,l)
           wake_ddeltat(i,l) = d_deltatw(i,l)
           wake_ddeltaq(i,l) = d_deltaqw(i,l)
        ELSE
           wake_deltat(i,l)= 0.
           wake_deltaq(i,l)= 0.
           wake_d_deltat_gw(i,l)= 0.
           wake_omgbdth(i,l) = 0.
           wake_dp_omgb(i,l) = 0.
           wake_dtKE(i,l) = 0.
           wake_dqKE(i,l) = 0.
 	   wake_dtPBL(i,l) = 0.
	   wake_dqPBL(i,l) = 0.
           wake_omg(i,l) = 0.
           wake_dp_deltomg(i,l) = 0.
           wake_spread(i,l) = 0.
           wake_dth(i,l)=0.
           dt_wake(i,l)=0.
           dq_wake(i,l)=0.
           undi_t(i,l)=te(i,l)
           undi_q(i,l)=qe(i,l)
           wake_ddeltat(i,l) = 0.
           wake_ddeltaq(i,l) = 0.
        ENDIF
       ENDDO
      ENDDO
c
      DO i=1,klon 
       wake_h(i)= hw(i)
       wake_s(i)= sigmaw(i)
       wake_pe(i)= wape(i)
       wake_fip(i)= fip(i)
       wake_gfl(i) = gfl(i)
       wake_k(i) =ktopw(i)
       wake_Cstar(i) = Cstar(i)
       wake_dens(i) = wdens(i)
      ENDDO
c
      RETURN
      END

      SUBROUTINE CALWAKE_scal(paprs,pplay,dtime
     :             ,t,q,omgb
     :             ,dt_dwn,dq_dwn,M_dwn,M_up
     :             ,dt_a,dq_a,sigd
     :             ,wdt_PBL,wdq_PBL
     :             ,udt_PBL,udq_PBL
     o             ,wake_deltat,wake_deltaq,wake_dth
     o             ,wake_h,wake_s,wake_dens
     o             ,wake_pe,wake_fip,wake_gfl
     o             ,dt_wake,dq_wake
     o             ,wake_k
     o             ,undi_t,undi_q
     o             ,wake_omgbdth,wake_dp_omgb
     o             ,wake_dtKE,wake_dqKE
     o             ,wake_dtPBL,wake_dqPBL
     o             ,wake_omg,wake_dp_deltomg
     o             ,wake_spread,wake_Cstar,wake_d_deltat_gw
     o             ,wake_ddeltat,wake_ddeltaq)
***************************************************************
*                                                             *
* CALWAKE                                                     *
*           interface avec le schema de calcul de la poche    *
*           froide                                            *
*                                                             *
* written by   : CHERUY Frederique, 13/03/2000, 10.31.05      *
* modified by :  ROEHRIG Romain,    01/30/2007                *
***************************************************************
*
      USE dimphy
      IMPLICIT none
c======================================================================

#include "dimensions.h"
cccc#include "dimphy.h"
#include "YOMCST.h"

c Arguments
c----------

      INTEGER  i,l,ktopw
      REAL   dtime

      REAL paprs(klon,klev+1),pplay(klon,klev)
      REAL t(klon,klev), q(klon,klev), omgb(klon,klev)
      REAL dt_dwn(klon,klev), dq_dwn(klon,klev),M_dwn(klon,klev)
      REAL M_up(klon,klev)
      REAL dt_a(klon,klev), dq_a(klon,klev)
      REAL wdt_PBL(klon,klev), wdq_PBL(klon,klev)
      REAL udt_PBL(klon,klev), udq_PBL(klon,klev)
      REAL wake_deltat(klon,klev),wake_deltaq(klon,klev)
      REAL dt_wake(klon,klev),dq_wake(klon,klev)
      REAL wake_d_deltat_gw(klon,klev)
      REAL wake_h(klon),wake_s(klon)
      REAL wake_dth(klon,klev)
      REAL wake_pe(klon),wake_fip(klon),wake_gfl(klon)
      REAL undi_t(klon,klev),undi_q(klon,klev)
      REAL wake_omgbdth(klon,klev),wake_dp_omgb(klon,klev)
      REAL wake_dtKE(klon,klev),wake_dqKE(klon,klev)
      REAL wake_dtPBL(klon,klev),wake_dqPBL(klon,klev)
      REAL wake_omg(klon,klev+1),wake_dp_deltomg(klon,klev)
      REAL wake_spread(klon,klev),wake_Cstar(klon)
      REAL wake_ddeltat(klon,klev),wake_ddeltaq(klon,klev)
      REAL d_deltatw(klev), d_deltaqw(klev)
      INTEGER wake_k(klon)
      REAL sigd(klon)
      REAL wake_dens(klon)

C  Variable internes
C  -----------------

      REAL aire
      REAL p(klev),ph(klev+1),pi(klev)
      REAL te(klev),qe(klev),omgbe(klev),dtdwn(klev),dqdwn(klev)
      REAL dta(klev),dqa(klev)
      REAL wdtPBL(klev),wdqPBL(klev)
      REAL udtPBL(klev),udqPBL(klev) 
      REAL amdwn(klev),amup(klev)
      REAL dtw(klev),dqw(klev),dth(klev),d_deltat_gw(klev)
      REAL dtls(klev),dqls(klev)
      REAL tu(klev),qu(klev)
      REAL hw,sigmaw,wape,fip,gfl
      REAL omgbdth(klev),dp_omgb(klev)
      REAL dtKE(klev),dqKE(klev)
      REAL dtPBL(klev),dqPBL(klev)
      REAL omg(klev+1),dp_deltomg(klev),spread(klev),Cstar
      REAL sigd0,wdens

      REAL RDCP

c      print *, '-> calwake, wake_s ', wake_s(1)

      RDCP=1./3.5

c-----------------------------------------------------------
      DO 999 i=1,klon   ! a vectoriser
c----------------------------------------------------------


      DO l=1,klev
        p(l)= pplay(i,l)
        ph(l)= paprs(i,l)
        pi(l) = (pplay(i,l)/100000.)**RDCP

        te(l) = t(i,l)
        qe(l) = q(i,l)
        omgbe(l) = omgb(i,l)

        dtdwn(l)= dt_dwn(i,l)
        dqdwn(l)= dq_dwn(i,l)
        dta(l)= dt_a(i,l)
        dqa(l)= dq_a(i,l)
        wdtPBL(l)= wdt_PBL(i,l)
        wdqPBL(l)= wdq_PBL(i,l)
        udtPBL(l)= udt_PBL(i,l)
        udqPBL(l)= udq_PBL(i,l)
      ENDDO
      
      sigd0=sigd(i)
c      print*, 'sigd0,sigd', sigd0, sigd(i)
      ph(klev+1)=0.

      ktopw = wake_k(i)

      DO l=1,klev
        dtw(l) = wake_deltat(i,l)
        dqw(l) = wake_deltaq(i,l)
      ENDDO

      DO l=1,klev
        dtls(l)=dt_wake(i,l)
        dqls(l)=dq_wake(i,l)
      ENDDO

      hw = wake_h(i)
      sigmaw = wake_s(i)

cfkc les flux de masses sont evalues aux niveaux et valent 0 a la surface
cfkc  on veut le flux de masse au milieu des couches

      DO l=1,klev-1
        amdwn(l)= 0.5*(M_dwn(i,l)+M_dwn(i,l+1))
        amdwn(l)= (M_dwn(i,l+1))
      ENDDO

c au sommet le flux de masse est nul

      amdwn(klev)=0.5*M_dwn(i,klev)
c
      DO l = 1,klev
        amup(l)=M_up(i,l)
      ENDDO

      call WAKE_scal(p,ph,pi,dtime,sigd0
     $                ,te,qe,omgbe
     $                ,dtdwn,dqdwn,amdwn,amup,dta,dqa
     $                ,wdtPBL,wdqPBL,udtPBL,udqPBL
     $                ,dtw,dqw,dth,hw,sigmaw,wape,fip,gfl
     $                ,dtls,dqls,ktopw
     $                ,omgbdth,dp_omgb,wdens
     $                ,tu,qu
     $                ,dtKE,dqKE
     $                ,dtPBL,dqPBL
     $                ,omg,dp_deltomg,spread
     $                ,Cstar,d_deltat_gw
     $                ,d_deltatw,d_deltaqw)

       IF (ktopw .GT. 0) THEN
         DO l=1,klev
           wake_deltat(i,l)= dtw(l)
           wake_deltaq(i,l)= dqw(l)
           wake_d_deltat_gw(i,l)= d_deltat_gw(l)
           wake_omgbdth(i,l) = omgbdth(l)
           wake_dp_omgb(i,l) = dp_omgb(l)
           wake_dtKE(i,l) = dtKE(l)
           wake_dqKE(i,l) = dqKE(l)
 	   wake_dtPBL(i,l) = dtPBL(l)
	   wake_dqPBL(i,l) = dqPBL(l)
           wake_omg(i,l) = omg(l)
           wake_dp_deltomg(i,l) = dp_deltomg(l)
           wake_spread(i,l) = spread(l)
           wake_dth(i,l) = dth(l)
           dt_wake(i,l) = dtls(l)
           dq_wake(i,l) = dqls(l)
           undi_t(i,l) = tu(l)
           undi_q(i,l) = qu(l)
           wake_ddeltat(i,l) = d_deltatw(l)
           wake_ddeltaq(i,l) = d_deltaqw(l)
         ENDDO
       ELSE
         DO l = 1,klev
           wake_deltat(i,l)= 0.
           wake_deltaq(i,l)= 0.
           wake_d_deltat_gw(i,l)= 0.
           wake_omgbdth(i,l) = 0.
           wake_dp_omgb(i,l) = 0.
           wake_dtKE(i,l) = 0.
           wake_dqKE(i,l) = 0.
           wake_omg(i,l) = 0.
           wake_dp_deltomg(i,l) = 0.
           wake_spread(i,l) = 0.
           wake_dth(i,l)=0.
           dt_wake(i,l)=0.
           dq_wake(i,l)=0.
           undi_t(i,l)=te(l)
           undi_q(i,l)=qe(l)
         ENDDO
       ENDIF

       wake_h(i)= hw
       wake_s(i)= sigmaw
       wake_pe(i)= wape
       wake_fip(i)= fip
       wake_gfl(i) = gfl
       wake_k(i) =ktopw
       wake_Cstar(i) = Cstar
       wake_dens(i) = wdens
c
 999  CONTINUE
c
      RETURN
      END
