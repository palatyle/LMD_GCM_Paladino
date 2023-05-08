!
! $Header$
!
c======================================================================
      SUBROUTINE orbite(xjour,longi,dist)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) (adapte du GCM du LMD) date: 19930818
c Objet: pour un jour donne, calculer la longitude vraie de la terre
c        (par rapport au point vernal-21 mars) dans son orbite solaire
c        calculer aussi la distance terre-soleil (unite astronomique)
c======================================================================
c Arguments:
c xjour--INPUT--R- jour de l'annee a compter du 1er janvier
c longi--OUTPUT-R- longitude vraie en degres par rapport au point
c                  vernal (21 mars) en degres
c dist---OUTPUT-R- distance terre-soleil (par rapport a la moyenne)
      REAL xjour, longi, dist
c======================================================================
#include "YOMCST.h"
C
C  -- Variables dynamiques locales
      REAL pir,xl,xllp,xee,xse,xlam,dlamm,anm,ranm,anv,ranv
C
      pir = 4.0*ATAN(1.0) / 180.0
      xl=R_peri+180.0
      xllp=xl*pir
      xee=R_ecc*R_ecc
      xse=SQRT(1.0-xee)
      xlam = (R_ecc/2.0+R_ecc*xee/8.0)*(1.0+xse)*SIN(xllp)
     .     - xee/4.0*(0.5+xse)*SIN(2.0*xllp)
     .     + R_ecc*xee/8.0*(1.0/3.0+xse)*SIN(3.0*xllp)
      xlam=2.0*xlam/pir
      dlamm=xlam+(xjour-81.0)
      anm=dlamm-xl
      ranm=anm*pir
      xee=xee*R_ecc
      ranv=ranm+(2.0*R_ecc-xee/4.0)*SIN(ranm)
     .         +5.0/4.0*R_ecc*R_ecc*SIN(2.0*ranm)
     .         +13.0/12.0*xee*SIN(3.0*ranm)
c
      anv=ranv/pir
      longi=anv+xl
C
      dist = (1-R_ecc*R_ecc)
     .      /(1+R_ecc*COS(pir*(longi-(R_peri+180.0))))
      RETURN
      END
c======================================================================
      SUBROUTINE angle(longi, lati, frac, muzero)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Calculer la duree d'ensoleillement pour un jour et la hauteur
c        du soleil (cosinus de l'angle zinithal) moyenne sur la journee
c======================================================================
c Arguments:
c longi----INPUT-R- la longitude vraie de la terre dans son plan 
c                   solaire a partir de l'equinoxe de printemps (degre)
c lati-----INPUT-R- la latitude d'un point sur la terre (degre)
c frac-----OUTPUT-R la duree d'ensoleillement dans la journee divisee
c                   par 24 heures (unite en fraction de 0 a 1)
c muzero---OUTPUT-R la moyenne du cosinus de l'angle zinithal sur
c                   la journee (0 a 1)
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
      REAL longi
      REAL lati(klon), frac(klon), muzero(klon)
#include "YOMCST.h"
      REAL lat, omega, lon_sun, lat_sun
      REAL pi_local, incl
      INTEGER i
c
      pi_local = 4.0 * ATAN(1.0)
      incl=R_incl * pi_local / 180.
c
      lon_sun = longi * pi_local / 180.0
      lat_sun = ASIN (sin(lon_sun)*SIN(incl) )
c
      DO i = 1, klon
      lat = lati(i) * pi_local / 180.0
c
      IF ( lat .GE. (pi_local/2.+lat_sun)
     .    .OR. lat.LE.(-pi_local/2.+lat_sun)) THEN
         omega = 0.0   ! nuit polaire
      ELSE IF ( lat.GE.(pi_local/2.-lat_sun)
     .          .OR. lat.LE.(-pi_local/2.-lat_sun)) THEN
         omega = pi_local   ! journee polaire
      ELSE
         omega = -TAN(lat)*TAN(lat_sun)
         omega = ACOS (omega)
      ENDIF
c
      frac(i) = omega / pi_local
c
      IF (omega .GT. 0.0) THEN
         muzero(i) = SIN(lat)*SIN(lat_sun)
     .          + COS(lat)*COS(lat_sun)*SIN(omega) / omega
      ELSE
         muzero(i) = 0.0
      ENDIF
      ENDDO
c
      RETURN
      END
c====================================================================
      SUBROUTINE zenang(longi,gmtime,pdtrad,lat,long,
     s                  pmu0,frac)
      USE dimphy
      IMPLICIT none
c=============================================================
c Auteur : O. Boucher (LMD/CNRS)
c          d'apres les routines zenith et angle de Z.X. Li 
c Objet  : calculer les valeurs moyennes du cos de l'angle zenithal
c          et l'ensoleillement moyen entre gmtime1 et gmtime2 
c          connaissant la declinaison, la latitude et la longitude.
c Rque   : Different de la routine angle en ce sens que zenang 
c          fournit des moyennes de pmu0 et non des valeurs 
c          instantanees, du coup frac prend toutes les valeurs 
c          entre 0 et 1.
c Date   : premiere version le 13 decembre 1994
c          revu pour  GCM  le 30 septembre 1996
c===============================================================
c longi : la longitude vraie de la terre dans son plan
c                  solaire a partir de l'equinoxe de printemps (degre)
c gmtime : temps universel en fraction de jour
c pdtrad : pas de temps du rayonnement (secondes)
c lat------INPUT : latitude en degres
c long-----INPUT : longitude en degres
c pmu0-----OUTPUT: angle zenithal moyen entre gmtime et gmtime+pdtrad
c frac-----OUTPUT: ensoleillement moyen entre gmtime et gmtime+pdtrad
c================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c================================================================
      real, intent(in):: longi, gmtime, pdtrad
      real lat(klon), long(klon), pmu0(klon), frac(klon)
c================================================================
      integer i
      real gmtime1, gmtime2
      real pi_local, deux_pi_local, incl
      real omega1, omega2, omega
c omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
c omega : heure en radian du coucher de soleil 
c -omega est donc l'heure en radian de lever du soleil
      real omegadeb, omegafin
      real zfrac1, zfrac2, z1_mu, z2_mu
      real lat_sun          ! declinaison en radian
      real lon_sun          ! longitude solaire en radian
      real latr             ! latitude du pt de grille en radian
c================================================================
c
      pi_local = 4.0 * ATAN(1.0)
      deux_pi_local = 2.0 * pi_local
      incl=R_incl * pi_local / 180.
c
      lon_sun = longi * pi_local / 180.0
      lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
c
      gmtime1=gmtime*86400.
      gmtime2=gmtime*86400.+pdtrad
c
      DO i = 1, klon
c
      latr = lat(i) * pi_local / 180.
c
c--pose probleme quand lat=+/-90 degres
c
c      omega = -TAN(latr)*TAN(lat_sun)
c      omega = ACOS(omega)
c      IF (latr.GE.(pi_local/2.+lat_sun)
c     .    .OR. latr.LE.(-pi_local/2.+lat_sun)) THEN
c         omega = 0.0       ! nuit polaire
c      ENDIF
c      IF (latr.GE.(pi_local/2.-lat_sun)
c     .          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
c         omega = pi_local  ! journee polaire
c      ENDIF
c
c--remplace par cela (le cas par defaut est different)
c
      omega=0.0  !--nuit polaire
      IF (latr.GE.(pi_local/2.-lat_sun)
     .          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
         omega = pi_local  ! journee polaire
      ENDIF
      IF (latr.LT.(pi_local/2.+lat_sun).AND.
     .    latr.GT.(-pi_local/2.+lat_sun).AND.
     .    latr.LT.(pi_local/2.-lat_sun).AND.
     .    latr.GT.(-pi_local/2.-lat_sun)) THEN
      omega = -TAN(latr)*TAN(lat_sun)
      omega = ACOS(omega)
      ENDIF
c
         omega1 = gmtime1 + long(i)*86400.0/360.0
         omega1 = omega1 / 86400.0*deux_pi_local
         omega1 = MOD (omega1+deux_pi_local, deux_pi_local)
         omega1 = omega1 - pi_local
c
         omega2 = gmtime2 + long(i)*86400.0/360.0
         omega2 = omega2 / 86400.0*deux_pi_local
         omega2 = MOD (omega2+deux_pi_local, deux_pi_local)
         omega2 = omega2 - pi_local
c
      IF (omega1.LE.omega2) THEN  !--on est dans la meme journee locale
c
      IF (omega2.LE.-omega .OR. omega1.GE.omega
     .                     .OR. omega.LT.1e-5) THEN   !--nuit
         frac(i)=0.0
         pmu0(i)=0.0
      ELSE                                              !--jour+nuit/jour
        omegadeb=MAX(-omega,omega1)
        omegafin=MIN(omega,omega2)
        frac(i)=(omegafin-omegadeb)/(omega2-omega1)
        pmu0(i)=SIN(latr)*SIN(lat_sun) + 
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omegafin)-SIN(omegadeb))/
     .          (omegafin-omegadeb)        
      ENDIF
c
      ELSE  !---omega1 GT omega2 -- a cheval sur deux journees
c
c-------------------entre omega1 et pi
      IF (omega1.GE.omega) THEN  !--nuit
         zfrac1=0.0
         z1_mu =0.0
      ELSE                       !--jour+nuit
        omegadeb=MAX(-omega,omega1)
        omegafin=omega
        zfrac1=omegafin-omegadeb
        z1_mu =SIN(latr)*SIN(lat_sun) +
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omegafin)-SIN(omegadeb))/
     .          (omegafin-omegadeb)
      ENDIF 
c---------------------entre -pi et omega2
      IF (omega2.LE.-omega) THEN   !--nuit
         zfrac2=0.0
         z2_mu =0.0
      ELSE                         !--jour+nuit
         omegadeb=-omega
         omegafin=MIN(omega,omega2)
         zfrac2=omegafin-omegadeb
         z2_mu =SIN(latr)*SIN(lat_sun) +
     .           COS(latr)*COS(lat_sun)*
     .           (SIN(omegafin)-SIN(omegadeb))/
     .           (omegafin-omegadeb)
c
      ENDIF
c-----------------------moyenne 
      frac(i)=(zfrac1+zfrac2)/(omega2+deux_pi_local-omega1)
      pmu0(i)=(zfrac1*z1_mu+zfrac2*z2_mu)/MAX(zfrac1+zfrac2,1.E-10)
c
      ENDIF   !---comparaison omega1 et omega2
c
      ENDDO
c
      END
c===================================================================
      SUBROUTINE zenith (longi, gmtime, lat, long,
     s                   pmu0, fract)
      USE dimphy
      IMPLICIT none
c
c Auteur(s): Z.X. Li (LMD/ENS)
c
c Objet: calculer le cosinus de l'angle zenithal du soleil en
c        connaissant la declinaison du soleil, la latitude et la
c        longitude du point sur la terre, et le temps universel
c
c Arguments d'entree:
c     longi  : declinaison du soleil (en degres)
c     gmtime : temps universel en second qui varie entre 0 et 86400
c     lat    : latitude en degres
c     long   : longitude en degres
c Arguments de sortie:
c     pmu0   : cosinus de l'angle zenithal
c
c====================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c====================================================================
      REAL longi, gmtime
      REAL lat(klon), long(klon), pmu0(klon), fract(klon)
c=====================================================================
      INTEGER n
      REAL zpi, zpir, omega, zgmtime
      REAL incl, lat_sun, lon_sun
c----------------------------------------------------------------------
      zpi = 4.0*ATAN(1.0)
      zpir = zpi / 180.0
      zgmtime=gmtime*86400.
c
      incl=R_incl * zpir
c
      lon_sun = longi * zpir
      lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
c
c--initialisation a la nuit
c
      DO n =1, klon
        pmu0(n)=0.
        fract(n)=0.0
      ENDDO
c
c 1 degre en longitude = 240 secondes en temps
c
      DO n = 1, klon
         omega = zgmtime + long(n)*86400.0/360.0
         omega = omega / 86400.0 * 2.0 * zpi
         omega = MOD(omega + 2.0 * zpi, 2.0 * zpi)
         omega = omega - zpi
         pmu0(n) = sin(lat(n)*zpir) * sin(lat_sun)
     .           + cos(lat(n)*zpir) * cos(lat_sun)
     .           * cos(omega)
         pmu0(n) = MAX (pmu0(n), 0.0)
         IF (pmu0(n).GT.1.E-6) fract(n)=1.0
      ENDDO
c
      RETURN
      END
