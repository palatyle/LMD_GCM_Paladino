!
! $Header$
!
c
c
      SUBROUTINE hydrol(dtime,pctsrf,rain_fall,snow_fall,evap,
     .                  agesno, tsol,qsol,snow,runoff)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c date: 19940414
c======================================================================
c
c Traitement de l'hydrologie du sol
c ---------------------------------
c rain_fall: taux de pluie
c snow_fall: taux de neige
c agesno: age de la neige
c evap: taux d'evaporation
c tsol: temperature du sol
c qsol: humidite du sol
c snow: couverture neigeuse
C
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "indicesol.h"
c
      REAL chasno ! epaisseur du sol: 0.15 m
      PARAMETER (chasno=3.334E+05/(2.3867E+06*0.15))
      REAL mx_eau_sol
      PARAMETER (mx_eau_sol=150.0)
c
      REAL dtime
      REAL pctsrf(klon,nbsrf)
      REAL snow(klon,nbsrf), tsol(klon,nbsrf), qsol(klon,nbsrf)
      REAL snow_fall(klon), rain_fall(klon), evap(klon)
      REAL runoff(klon), agesno(klon)
C
      INTEGER i, is
      REAL subli, fsno
C-----------------------------------------------------------------------
      DO 99999 i = 1, klon
c
         runoff(i) = 0.0
c
         is = is_ter
         snow(i,is) = snow(i,is) + snow_fall(i) * dtime * pctsrf(i,is)
         IF (pctsrf(i,is) .GT. epsfra) THEN
            subli = MIN(evap(i)*dtime,snow(i,is))
            snow(i,is) = snow(i,is) - subli
            fsno = MIN(MAX((tsol(i,is)-RTT)/chasno,0.0),snow(i,is))
            snow(i,is) = snow(i,is) - fsno
            tsol(i,is) = tsol(i,is) - fsno*chasno
            qsol(i,is) = qsol(i,is) + (rain_fall(i)-evap(i))*dtime
     .                              + subli + fsno
            qsol(i,is) = MAX(qsol(i,is),0.0)
            runoff(i) = runoff(i) + MAX(qsol(i,is)-mx_eau_sol, 0.0)
     .                            * pctsrf(i,is)
            qsol(i,is) = MIN(qsol(i,is),mx_eau_sol)
ccc         ELSE
ccc            snow(i,is) = 0.0
ccc            qsol(i,is) = 0.0
ccc            tsol(i,is) = 0.0
         ENDIF
c
         is = is_lic
         snow(i,is) = snow(i,is) + snow_fall(i) * dtime * pctsrf(i,is)
         IF (pctsrf(i,is) .GT. epsfra) THEN
            subli = MIN(evap(i)*dtime,snow(i,is))
            snow(i,is) = snow(i,is) - subli
            fsno = MIN(MAX((tsol(i,is)-RTT)/chasno,0.0),snow(i,is))
            snow(i,is) = snow(i,is) - fsno
            tsol(i,is) = tsol(i,is) - fsno*chasno
            qsol(i,is) = qsol(i,is) + (rain_fall(i)-evap(i))*dtime
     .                              + subli + fsno
            qsol(i,is) = MAX(qsol(i,is),0.0)
            runoff(i) = runoff(i) + MAX(qsol(i,is)-mx_eau_sol, 0.0)
     .                            * pctsrf(i,is)
            qsol(i,is) = MIN(qsol(i,is),mx_eau_sol)
c je limite la temperature a RTT-1.8 (il faudrait aussi prendre l'eau de
c la fonte) (Laurent Li, le 14mars98):
cIM cf GK   tsol(i,is) = MIN(tsol(i,is),RTT-1.8)
cIM cf GK : la glace fond a 0C, non pas a -1.8
            tsol(i,is) = MIN(tsol(i,is),RTT)
c
ccc         ELSE
ccc            snow(i,is) = 0.0
ccc            qsol(i,is) = 0.0
ccc            tsol(i,is) = 0.0
         ENDIF
c
         is = is_sic
         qsol(i,is) = 0.0
         snow(i,is) = snow(i,is) + snow_fall(i) * dtime * pctsrf(i,is)
         IF (pctsrf(i,is) .GT. epsfra) THEN
            subli = MIN(evap(i)*dtime,snow(i,is))
            snow(i,is) = snow(i,is) - subli
            fsno = MIN(MAX((tsol(i,is)-RTT)/chasno,0.0),snow(i,is))
            snow(i,is) = snow(i,is) - fsno
            tsol(i,is) = tsol(i,is) - fsno*chasno
c je limite la temperature a RTT-1.8 (il faudrait aussi prendre l'eau de
c la fonte) (Laurent Li, le 14mars98):
cIM cf GK   tsol(i,is) = MIN(tsol(i,is),RTT-1.8)
cIM cf GK : la glace fond a 0C, non pas a -1.8
            tsol(i,is) = MIN(tsol(i,is),RTT)
c
ccc         ELSE
ccc            snow(i,is) = 0.0
ccc            tsol(i,is) = 0.0
         ENDIF
c
         agesno(i) = (agesno(i)+ (1.-agesno(i)/50.)*dtime/86400.)
     .             * EXP(-1.*MAX(0.0,snow_fall(i))*dtime/0.3)
         agesno(i) = MAX(agesno(i),0.0)
c
99999 CONTINUE
c
      RETURN
      END
