c====================================================================
      SUBROUTINE zenang(longi,gmtime,pdtrad,lat,long,
     s                  pmu0,frac)
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
c longi----INPUT : la longitude vraie de la terre dans son plan
c                  solaire a partir de l'equinoxe de printemps (degre)
c gmtime---INPUT : temps universel en fraction de jour
c pdtrad---INPUT : pas de temps du rayonnement (secondes)
c lat------INPUT : latitude en degres
c long-----INPUT : longitude en degres
c pmu0-----OUTPUT: angle zenithal moyen entre gmtime et gmtime+pdtrad/RDAY
c frac-----OUTPUT: ensoleillement moyen entre gmtime et gmtime+pdtrad/RDAY
c================================================================
      use dimphy
      IMPLICIT none
#include "YOMCST.h"
c================================================================
      real longi, gmtime, pdtrad
      real lat(klon), long(klon), pmu0(klon), frac(klon)
c================================================================
      integer i
      real gmtime1, gmtime2
      real incl
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
      incl=R_incl * RPI / 180.
c
      lon_sun = longi * RPI / 180.0
      lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
c
      gmtime1=gmtime*RDAY
      gmtime2=gmtime*RDAY+pdtrad
c
      DO i = 1, klon
c
      latr = lat(i) * RPI / 180.
c
c--pose probleme quand lat=+/-90 degres
c
c      omega = -TAN(latr)*TAN(lat_sun)
c      omega = ACOS(omega)
c      IF (latr.GE.(RPI/2.+lat_sun)
c     .    .OR. latr.LE.(-RPI/2.+lat_sun)) THEN
c         omega = 0.0       ! nuit polaire
c      ENDIF
c      IF (latr.GE.(RPI/2.-lat_sun)
c     .          .OR. latr.LE.(-RPI/2.-lat_sun)) THEN
c         omega = RPI  ! journee polaire
c      ENDIF
c
c--remplace par cela (le cas par defaut est different)
c
      omega=0.0  !--nuit polaire
      IF (latr.GE.(RPI/2.-lat_sun)
     .          .OR. latr.LE.(-RPI/2.-lat_sun)) THEN
         omega = RPI  ! journee polaire
      ENDIF
      IF (latr.LT.(RPI/2.+lat_sun).AND.
     .    latr.GT.(-RPI/2.+lat_sun).AND.
     .    latr.LT.(RPI/2.-lat_sun).AND.
     .    latr.GT.(-RPI/2.-lat_sun)) THEN
      omega = -TAN(latr)*TAN(lat_sun)
      omega = ACOS(omega)
      ENDIF
c
         omega1 = gmtime1 + long(i)*RDAY/360.0
         omega1 = omega1*2*RPI / RDAY
         omega1 = MOD (omega1+2*RPI, 2*RPI)
         omega1 = omega1 - RPI
c
         omega2 = gmtime2 + long(i)*RDAY/360.0
         omega2 = omega2*2*RPI / RDAY
         omega2 = MOD (omega2+2*RPI, 2*RPI)
         omega2 = omega2 - RPI
c
      IF (omega1.LE.omega2) THEN  !--on est dans la meme journee locale
c
      IF (omega2.LE.-omega .OR. omega1.GE.omega      !--nuit
     .           .OR. omega.LT.1e-5) THEN            !--nuit polaire
        frac(i)=0.0
        pmu0(i)=SIN(latr)*SIN(lat_sun) + 
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omega2)-SIN(omega1))/
     .          (omega2-omega1)        
      ELSE                                           !--jour+nuit/jour
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
         z1_mu =SIN(latr)*SIN(lat_sun) + 
     .          COS(latr)*COS(lat_sun)*
     .          (-SIN(omega1))/       
     .          (RPI-omega1)
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
         z2_mu =SIN(latr)*SIN(lat_sun) + 
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omega2))/       
     .          (omega2+RPI) 
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
      frac(i)=(zfrac1+zfrac2)/(omega2+2*RPI-omega1)
      if (frac(i).ne.0.) then
       pmu0(i)=(zfrac1*z1_mu+zfrac2*z2_mu)/MAX(zfrac1+zfrac2,1.E-10)
      else
       pmu0(i)=((RPI-omega1)*z1_mu+(omega2+RPI)*z2_mu)/
     .                   (omega2+2*RPI-omega1)
      endif
c
      ENDIF   !---comparaison omega1 et omega2
c
      ENDDO
c
      END
