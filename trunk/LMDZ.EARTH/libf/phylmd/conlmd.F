!
! $Header$
!
      SUBROUTINE conlmd (dtime, paprs, pplay, t, q, conv_q,
     s                   d_t, d_q, rain, snow, ibas, itop)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Schema de convection utilis'e dans le modele du LMD
c        Ajustement humide (Manabe) + Ajustement convectif (Kuo)
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOETHF.h"
c
c Arguments:
c
      REAL dtime              ! pas d'integration (s)
      REAL paprs(klon,klev+1) ! pression inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
      REAL t(klon,klev)       ! temperature (K)
      REAL q(klon,klev)       ! humidite specifique (kg/kg)
      REAL conv_q(klon,klev)  ! taux de convergence humidite (g/g/s)
c
      REAL d_t(klon,klev)     ! incrementation temperature
      REAL d_q(klon,klev)     ! incrementation humidite
      REAL rain(klon)         ! pluies (mm/s)
      REAL snow(klon)         ! neige (mm/s)
      INTEGER ibas(klon)      ! niveau du bas
      INTEGER itop(klon)      ! niveau du haut
c
      LOGICAL usekuo ! utiliser convection profonde (schema Kuo)
      PARAMETER (usekuo=.TRUE.)
c
      REAL d_t_bis(klon,klev)
      REAL d_q_bis(klon,klev)
      REAL rain_bis(klon)
      REAL snow_bis(klon)
      INTEGER ibas_bis(klon)
      INTEGER itop_bis(klon)
      REAL d_ql(klon,klev), d_ql_bis(klon,klev)
      REAL rneb(klon,klev), rneb_bis(klon,klev)
c
      INTEGER i, k
      REAL zlvdcp, zlsdcp, zdelta, zz, za, zb
c
ccc      CALL fiajh ! ancienne version de Convection Manabe
      CALL conman ! nouvelle version de Convection Manabe
     e     (dtime, paprs, pplay, t, q,
     s      d_t, d_q, d_ql, rneb,
     s      rain, snow, ibas, itop)
c
      IF (usekuo) THEN
ccc      CALL fiajc ! ancienne version de Convection Kuo
      CALL conkuo ! nouvelle version de Convection Kuo
     e     (dtime, paprs, pplay, t, q, conv_q,
     s      d_t_bis, d_q_bis, d_ql_bis, rneb_bis, 
     s      rain_bis, snow_bis, ibas_bis, itop_bis)
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = d_t(i,k) + d_t_bis(i,k)
         d_q(i,k) = d_q(i,k) + d_q_bis(i,k)
         d_ql(i,k) = d_ql(i,k) + d_ql_bis(i,k)
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = rain(i) + rain_bis(i)
         snow(i) = snow(i) + snow_bis(i)
         ibas(i) = MIN(ibas(i),ibas_bis(i))
         itop(i) = MAX(itop(i),itop_bis(i))
      ENDDO
      ENDIF
c
c L'eau liquide convective est dispersee dans l'air:
c
      DO k = 1, klev
      DO i = 1, klon
         zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q(i,k))
         zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*q(i,k))
         zdelta = MAX(0.,SIGN(1.,RTT-t(i,k)))
         zz = d_ql(i,k) ! re-evap. de l'eau liquide
         zb = MAX(0.0,zz)
         za = - MAX(0.0,zz) * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
         d_t(i,k) = d_t(i,k) + za
         d_q(i,k) = d_q(i,k) + zb
      ENDDO
      ENDDO
c
      RETURN
      END
      SUBROUTINE conman (dtime, paprs, pplay, t, q,
     s                   d_t, d_q, d_ql, rneb,
     s                   rain, snow, ibas, itop)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19970324
c Objet: ajustement humide convectif avec la possibilite de faire
c        l'ajustement sur une fraction de la maille.
c Methode: On impose une distribution uniforme pour la vapeur d'eau
c au sein d'une maille. On applique la procedure d'ajustement
c successivement a la totalite, 75%, 50%, 25% et 5% de la maille
c jusqu'a ce que l'ajustement a lieu. J'espere que ceci augmente
c les activites convectives et corrige le biais "trop froid et sec"
c du modele.
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c
      REAL dtime              ! pas d'integration (s)
      REAL t(klon,klev)       ! temperature (K)
      REAL q(klon,klev)       ! humidite specifique (kg/kg)
      REAL paprs(klon,klev+1) ! pression inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
c
      REAL d_t(klon,klev)     ! incrementation temperature
      REAL d_q(klon,klev)     ! incrementation humidite
      REAL d_ql(klon,klev)    ! incrementation eau liquide
      REAL rneb(klon,klev)    ! nebulosite
      REAL rain(klon)         ! pluies (mm/s)
      REAL snow(klon)         ! neige (mm/s)
      INTEGER ibas(klon)      ! niveau du bas
      INTEGER itop(klon)      ! niveau du haut
c
      LOGICAL afaire(klon)   ! .TRUE. implique l'ajustement
      LOGICAL accompli(klon) ! .TRUE. si l'ajustement est effectif
c
      INTEGER nb ! nombre de sous-fractions a considere
      PARAMETER (nb=1)
ccc      PARAMETER (nb=3)
c
      REAL ratqs ! largeur de la distribution pour vapeur d'eau
      PARAMETER (ratqs=0.05)
c
      REAL w_q(klon,klev)
      REAL w_d_t(klon,klev), w_d_q(klon,klev), w_d_ql(klon,klev)
      REAL w_rneb(klon,klev)
      REAL w_rain(klon), w_snow(klon)
      INTEGER w_ibas(klon), w_itop(klon)
      REAL zq1, zq2
      INTEGER i, k, n
c
      REAL t_coup
      PARAMETER (t_coup=234.0)
      REAL zdp1, zdp2
      REAL zqs1, zqs2, zdqs1, zdqs2
      REAL zgamdz
      REAL zflo ! flotabilite
      REAL zsat ! sur-saturation
      REAL zdelta, zcor, zcvm5
      LOGICAL imprim
c
      INTEGER ncpt
      SAVE ncpt
c$OMP THREADPRIVATE(ncpt)
      REAL frac(nb) ! valeur de la maille fractionnelle
      SAVE frac
c$OMP THREADPRIVATE(frac)
      INTEGER opt_cld(nb) ! option pour le modele nuageux
      SAVE opt_cld
c$OMP THREADPRIVATE(opt_cld)
      LOGICAL appel1er
      SAVE appel1er
c$OMP THREADPRIVATE(appel1er)
c
c Fonctions thermodynamiques:
c
#include "YOETHF.h"
#include "FCTTRE.h"
c
      DATA frac / 1.0 /
      DATA opt_cld / 4 /
ccc      DATA frac    / 1.0, 0.50, 0.25/
ccc      DATA opt_cld / 4,   4,    4/
c
      DATA appel1er /.TRUE./
      DATA ncpt /0/
c
      IF (appel1er) THEN
         PRINT*, 'conman, nb:', nb
         PRINT*, 'conman, frac:', frac
         PRINT*, 'conman, opt_cld:', opt_cld
         appel1er = .FALSE.
      ENDIF
c
c Initialiser les sorties a zero:
c
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
         d_ql(i,k) = 0.0
         rneb(i,k) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         ibas(i) = klev
         itop(i) = 1
         rain(i) = 0.0
         snow(i) = 0.0
      ENDDO
c
c S'il n'y a pas d'instabilite conditionnelle,
c pas la penne de se fatiguer:
c
      DO i = 1, klon
         afaire(i) = .FALSE.
      ENDDO
      DO k = 1, klev-1
      DO i = 1, klon
         IF (thermcep) THEN
            zdelta=MAX(0.,SIGN(1.,RTT-t(i,k)))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*q(i,k))
            zqs1= R2ES*FOEEW(t(i,k),zdelta)/pplay(i,k)
            zqs1=MIN(0.5,zqs1)
            zcor=1./(1.-RETV*zqs1)
            zqs1=zqs1*zcor
            zdqs1 =FOEDE(t(i,k),zdelta,zcvm5,zqs1,zcor)
c
            zdelta=MAX(0.,SIGN(1.,RTT-t(i,k+1)))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*q(i,k+1))
            zqs2= R2ES*FOEEW(t(i,k+1),zdelta)/pplay(i,k+1)
            zqs2=MIN(0.5,zqs2)
            zcor=1./(1.-RETV*zqs2)
            zqs2=zqs2*zcor
            zdqs2 =FOEDE(t(i,k+1),zdelta,zcvm5,zqs2,zcor)
         ELSE
           IF (t(i,k) .LT. t_coup) THEN
              zqs1= qsats(t(i,k)) / pplay(i,k)
              zdqs1= dqsats(t(i,k),zqs1)
c
              zqs2= qsats(t(i,k+1)) / pplay(i,k+1)
              zdqs2= dqsats(t(i,k+1),zqs2)
           ELSE
              zqs1= qsatl(t(i,k)) / pplay(i,k)
              zdqs1= dqsatl(t(i,k),zqs1)
c
              zqs2= qsatl(t(i,k+1)) / pplay(i,k+1)
              zdqs2= dqsatl(t(i,k+1),zqs2)
           ENDIF
         ENDIF
         zdp1 = paprs(i,k) - paprs(i,k+1)
         zdp2 = paprs(i,k+1) - paprs(i,k+2)
         zgamdz = - (pplay(i,k)-pplay(i,k+1))/paprs(i,k+1)/RCPD
     .          *( RD*(t(i,k)*zdp1+t(i,k+1)*zdp2)/(zdp1+zdp2)
     .            +RLVTT*(zqs1*zdp1+zqs2*zdp2)/(zdp1+zdp2)
     .           ) / (1.0+(zdqs1*zdp1+zdqs2*zdp2)/(zdp1+zdp2) )
         zflo = t(i,k) + zgamdz - t(i,k+1)
         zsat = (q(i,k)-zqs1)*zdp1 + (q(i,k+1)-zqs2)*zdp2
         IF (zflo.GT.0.0) afaire(i) = .TRUE.
c erreur         IF (zflo.GT.0.0 .AND. zsat.GT.0.0) afaire(i) = .TRUE.
      ENDDO
      ENDDO
c
      imprim = MOD(ncpt,48).EQ.0
      DO 99999 n = 1, nb
c
      DO k = 1, klev
      DO i = 1, klon
      IF (afaire(i)) THEN
         zq1 = q(i,k) * (1.0-ratqs)
         zq2 = q(i,k) * (1.0+ratqs)
         w_q(i,k) = zq2 - frac(n)/2.0 * (zq2-zq1)
      ENDIF
      ENDDO
      ENDDO
c
      CALL conmanv (dtime, paprs, pplay, t, w_q,
     e              afaire, opt_cld(n),
     s              w_d_t, w_d_q, w_d_ql, w_rneb,
     s              w_rain, w_snow, w_ibas, w_itop,accompli,imprim)
      DO k = 1, klev
      DO i = 1, klon
      IF (afaire(i) .AND. accompli(i)) THEN
         d_t(i,k) = w_d_t(i,k) * frac(n)
         d_q(i,k) = w_d_q(i,k) * frac(n)
         d_ql(i,k) = w_d_ql(i,k) * frac(n)
         IF (NINT(w_rneb(i,k)).EQ.1) rneb(i,k) = frac(n)
      ENDIF
      ENDDO
      ENDDO
      DO i = 1, klon
      IF (afaire(i) .AND. accompli(i)) THEN
         rain(i) = w_rain(i) * frac(n)
         snow(i) = w_snow(i) * frac(n)
         ibas(i) = MIN(ibas(i),w_ibas(i))
         itop(i) = MAX(itop(i),w_itop(i))
      ENDIF
      ENDDO
      DO i = 1, klon
         IF(afaire(i) .AND. accompli(i)) afaire(i) = .FALSE.
      ENDDO
c
99999 CONTINUE
c
      ncpt = ncpt + 1
c
      RETURN
      END
      SUBROUTINE conmanv (dtime, paprs, pplay, t, q,
     e                    afaire, opt_cld,
     s                    d_t, d_q, d_ql, rneb,
     s                    rain, snow, ibas, itop,accompli,imprim)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: ajustement humide (convection proposee par Manabe).
c        Pour une colonne verticale, il peut avoir plusieurs blocs
c        necessitant l'ajustement. ibas est le bas du plus bas bloc
c        et itop est le haut du plus haut bloc
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c
c Arguments:
c
      REAL dtime              ! pas d'integration (s)
      REAL t(klon,klev)       ! temperature (K)
      REAL q(klon,klev)       ! humidite specifique (kg/kg)
      REAL paprs(klon,klev+1) ! pression inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
      INTEGER opt_cld ! comment traiter l'eau liquide
      LOGICAL afaire(klon) ! .TRUE. si le point est a faire (Input)
      LOGICAL imprim ! .T. pour imprimer quelques diagnostiques
c
      REAL d_t(klon,klev)     ! incrementation temperature
      REAL d_q(klon,klev)     ! incrementation humidite
      REAL d_ql(klon,klev)    ! incrementation eau liquide
      REAL rneb(klon,klev)    ! nebulosite
      REAL rain(klon)         ! pluies (mm/s)
      REAL snow(klon)         ! neige (mm/s)
      INTEGER ibas(klon)      ! niveau du bas
      INTEGER itop(klon)      ! niveau du haut
      LOGICAL accompli(klon) ! .TRUE. si l'ajustement a eu lieu (Output)
c
c Quelques options:
c
      LOGICAL new_top ! re-calculer sommet quand re-ajustement est fait
      PARAMETER (new_top=.FALSE.)
      LOGICAL evap_prec ! evaporation de pluie au-dessous de convection
      PARAMETER (evap_prec=.TRUE.)
      REAL coef_eva
      PARAMETER (coef_eva=1.0E-05)
      REAL t_coup
      PARAMETER (t_coup=234.0)
      REAL seuil_vap
      PARAMETER (seuil_vap=1.0E-10)
      LOGICAL old_tau ! implique precip nulle, si vrai.
      PARAMETER (old_tau=.FALSE.)
      REAL toliq(klon) ! rapport entre l'eau nuageuse et l'eau precipitante
      REAL dpmin, tomax !Epaisseur faible, rapport eau liquide plus grande
      PARAMETER (dpmin=0.15, tomax=0.97)
      REAL dpmax, tomin !Epaisseur grande, rapport eau liquide plus faible
      PARAMETER (dpmax=0.30, tomin=0.05)
      REAL deep_sig, deep_to ! au dela de deep_sig, utiliser deep_to
      PARAMETER (deep_sig=0.50, deep_to=0.05)
      LOGICAL exigent ! implique un calcul supplementaire pour Qs
      PARAMETER (exigent=.FALSE.)
c
      INTEGER kbase
      PARAMETER (kbase=0)
c
c Variables locales:
c
      INTEGER nexpo
      INTEGER i, k, k1min, k1max, k2min, k2max, is
      REAL zgamdz(klon,klev-1)
      REAL zt(klon,klev), zq(klon,klev)
      REAL zqs(klon,klev), zdqs(klon,klev)
      REAL zqmqsdp(klon,klev)
      REAL ztnew(klon,klev), zqnew(klon,klev)
      REAL zcond(klon), zvapo(klon), zrapp(klon)
      REAL zrfl(klon), zrfln, zqev, zqevt
      REAL zsat(klon) ! sur-saturation
      REAL zflo(klon) ! flotabilite
      REAL za(klon), zb(klon), zc(klon)
      INTEGER k1(klon), k2(klon)
      REAL zdelta, zcor, zcvm5
      REAL delp(klon,klev)
      LOGICAL possible(klon), todo(klon), etendre(klon)
      LOGICAL aller(klon), todobis(klon)
      REAL zalfa
      INTEGER nbtodo, nbdone
c
c Fonctions thermodynamiques:
c
#include "YOETHF.h"
#include "FCTTRE.h"
c
      DO k = 1, klev
      DO i = 1, klon
         delp(i,k) = paprs(i,k) - paprs(i,k+1)
      ENDDO
      ENDDO
c
c Initialiser les sorties a zero
c
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
         d_ql(i,k) = 0.0
         rneb(i,k) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         ibas(i) = klev
         itop(i) = 1
         rain(i) = 0.0
         snow(i) = 0.0
         accompli(i) = .FALSE.
      ENDDO
c
c Preparations
c
      DO k = 1, klev
      DO i = 1, klon
      IF (afaire(i)) THEN
         zt(i,k) = t(i,k)
         zq(i,k) = q(i,k) 
c
c        Calculer Qs et L/Cp*dQs/dT
c
         IF (thermcep) THEN
            zdelta=MAX(0.,SIGN(1.,RTT-zt(i,k)))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i,k))
            zqs(i,k)= R2ES*FOEEW(zt(i,k),zdelta)/pplay(i,k)
            zqs(i,k)=MIN(0.5,zqs(i,k))
            zcor=1./(1.-RETV*zqs(i,k))
            zqs(i,k)=zqs(i,k)*zcor
            zdqs(i,k) =FOEDE(zt(i,k),zdelta,zcvm5,zqs(i,k),zcor)
         ELSE
           IF (zt(i,k) .LT. t_coup) THEN
              zqs(i,k)= qsats(zt(i,k)) / pplay(i,k)
              zdqs(i,k)= dqsats(zt(i,k),zqs(i,k))
           ELSE
              zqs(i,k)= qsatl(zt(i,k)) / pplay(i,k)
              zdqs(i,k)= dqsatl(zt(i,k),zqs(i,k))
           ENDIF
         ENDIF
c
c        Calculer (q-qs)*dp
         zqmqsdp(i,k) = (zq(i,k)-zqs(i,k)) * delp(i,k)
      ENDIF
      ENDDO
      ENDDO
c
c-----zgama is the moist convective lapse rate (-dT/dz).
c-----zgamdz(*,k) est la difference minimale autorisee des temperatures
c-----entre deux couches (k et k+1), c.a.d. si T(k+1)-T(k) est inferieur
c-----a zgamdz(*,k), alors ces 2 couches sont instables conditionnellement
c
      DO k = 1, klev-1
      DO i = 1, klon
      IF (afaire(i)) THEN
         zgamdz(i,k) = - (pplay(i,k)-pplay(i,k+1))/paprs(i,k+1)/RCPD
     .          *( RD*(zt(i,k)*delp(i,k)+zt(i,k+1)*delp(i,k+1))
     .               /(delp(i,k)+delp(i,k+1))
     .            +RLVTT*(zqs(i,k)*delp(i,k)+zqs(i,k+1)*delp(i,k+1))
     .                  /(delp(i,k)+delp(i,k+1))
     .           ) / (1.0+(zdqs(i,k)*delp(i,k)+zdqs(i,k+1)*delp(i,k+1))
     .                   /(delp(i,k)+delp(i,k+1)) )
      ENDIF
      ENDDO
      ENDDO
c
c On cherche la presence simultanee d'instabilite conditionnelle
c et de sur-saturation. Sinon, pas la penne de se fatiguer:
c
      DO i = 1, klon
         possible(i) = .FALSE.
      ENDDO
      DO k = 2, klev
      DO i = 1, klon
      IF (afaire(i)) THEN
         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
         zsat(i) = zqmqsdp(i,k) + zqmqsdp(i,k-1)
         IF (zflo(i).GT.0.0 .AND. zsat(i).GT.0.0) possible(i) = .TRUE.
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (possible(i)) THEN
         k1(i) = kbase
         k2(i) = k1(i) + 1
      ENDIF
      ENDDO
c
  810 CONTINUE ! chercher le bas de la colonne a ajuster
c
      k2min = klev
      DO i = 1, klon
         todo(i) = .FALSE.
         aller(i) = .TRUE.
         IF (possible(i)) k2min = MIN(k2min,k2(i))
      ENDDO
      IF (k2min.EQ.klev) GOTO 860
      DO k = k2min, klev-1
      DO i = 1, klon
      IF (possible(i) .AND. k.GE.k2(i) .AND. aller(i)) THEN
         zflo(i) = zt(i,k) + zgamdz(i,k) - zt(i,k+1)
         zsat(i) = zqmqsdp(i,k) + zqmqsdp(i,k+1)
         IF (zflo(i).GT.0.0 .AND. zsat(i).GT.0.0) THEN
            k1(i) = k
            k2(i) = k+1
            todo(i) = .TRUE.
            aller(i) = .FALSE.
         ENDIF
      ENDIF
      ENDDO
      ENDDO
      DO i = 1, klon
      IF (possible(i).AND.aller(i)) THEN
         todo(i) = .FALSE.
         k1(i) = klev
         k2(i) = klev
      ENDIF
      ENDDO
c
CCC      DO i = 1, klon
CCC      IF (possible(i)) THEN
CCC  811    k2(i) = k2(i) + 1
CCC         IF (k2(i) .GT. klev) THEN
CCC            todo(i) = .FALSE.
CCC            GOTO 812
CCC         ENDIF
CCC         k = k2(i)
CCC         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
CCC         zsat(i) = zqmqsdp(i,k) + zqmqsdp(i,k-1)
CCC         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) GOTO 811
CCC         k1(i) = k2(i) - 1
CCC         todo(i) = .TRUE.
CCC      ENDIF
CCC  812 CONTINUE
CCC      ENDDO
c
  820 CONTINUE ! chercher le haut de la colonne
c
      k2min = klev
      DO i = 1, klon
         aller(i) = .TRUE.
         IF (todo(i)) k2min = MIN(k2min,k2(i))
      ENDDO
      IF (k2min.LT.klev) THEN
      DO k = k2min, klev
      DO i = 1, klon
      IF (todo(i) .AND. k.GT.k2(i) .AND. aller(i)) THEN
            zsat(i) = zsat(i) + zqmqsdp(i,k)
            zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
            IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) THEN
               aller(i) = .FALSE.
            ELSE
               k2(i) = k
            ENDIF
      ENDIF
      ENDDO
      ENDDO
c error      is = 0
c error      DO i = 1, klon
c error      IF(todo(i).AND.aller(i)) THEN
c error         is = is + 1
c error         todo(i) = .FALSE.
c error         k2(i) = klev
c error      ENDIF
c error      ENDDO
c error      IF (is.GT.0) THEN
c error         PRINT*, "Bizard. je pourrais continuer mais j arrete"
c error         CALL abort
c error      ENDIF
      ENDIF
c
CCC      DO i = 1, klon
CCC      IF (todo(i)) THEN
CCC  821    CONTINUE
CCC         IF (k2(i) .EQ. klev) GOTO 822
CCC         k = k2(i) + 1
CCC         zsat(i) = zsat(i) + zqmqsdp(i,k)
CCC         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
CCC         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) GOTO 822
CCC         k2(i) = k
CCC         GOTO 821
CCC      ENDIF
CCC  822 CONTINUE
CCC      ENDDO
c
  830 CONTINUE ! faire l'ajustement en sachant k1 et k2
c
      is = 0
      DO i = 1, klon
      IF (todo(i)) THEN
         IF (k2(i).LE.k1(i)) is = is + 1
      ENDIF
      ENDDO
      IF (is.GT.0) THEN
         PRINT*, "Impossible: k1 trop grand ou k2 trop petit"
         PRINT*, "is=", is
         CALL abort
      ENDIF
c
      k1min = klev
      k1max = 1
      k2max = 1
      DO i = 1, klon
      IF (todo(i)) THEN
         k1min = MIN(k1min,k1(i))
         k1max = MAX(k1max,k1(i))
         k2max = MAX(k2max,k2(i))
      ENDIF
      ENDDO
c
      DO i = 1, klon
      IF (todo(i)) THEN
      k = k1(i)
      za(i) = 0.
      zb(i) = ( RCPD*(1.+zdqs(i,k))*(zt(i,k)-za(i))
     .      -RLVTT*(zqs(i,k)-zq(i,k)) )*delp(i,k)
      zc(i) = delp(i,k) * RCPD*(1.+zdqs(i,k))
      ENDIF
      ENDDO
c
      DO k = k1min, k2max
      DO i = 1, klon
      IF (todo(i) .AND. k.GE.(k1(i)+1) .AND. k.LE.k2(i)) THEN
         za(i) = za(i) + zgamdz(i,k-1)
         zb(i) = zb(i)+(RCPD*(1.+zdqs(i,k))*(zt(i,k)-za(i))
     .           -RLVTT*(zqs(i,k)-zq(i,k)) ) * delp(i,k)
         zc(i) = zc(i) + delp(i,k)*RCPD*(1.+zdqs(i,k))
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (todo(i)) THEN
         k = k1(i)
         ztnew(i,k) = zb(i)/zc(i)
         zqnew(i,k) = zqs(i,k) + (ztnew(i,k)-zt(i,k))
     .                          *RCPD/RLVTT*zdqs(i,k)
      ENDIF
      ENDDO
c
      DO k = k1min, k2max
      DO i = 1, klon
      IF (todo(i) .AND. k.GE.(k1(i)+1) .AND. k.LE.k2(i)) THEN
         ztnew(i,k) = ztnew(i,k-1) + zgamdz(i,k-1)
         zqnew(i,k) = zqs(i,k) + (ztnew(i,k)-zt(i,k))
     .                        *RCPD/RLVTT*zdqs(i,k)
      ENDIF
      ENDDO
      ENDDO
c
c Quantite de condensation produite pendant l'ajustement:
c
      DO i = 1, klon
         zcond(i) = 0.0
      ENDDO
      DO k = k1min, k2max
      DO i = 1, klon
      IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i)) THEN
         rneb(i,k) = 1.0
         zcond(i) = zcond(i) + (zq(i,k)-zqnew(i,k)) *delp(i,k)/RG
      ENDIF
      ENDDO
      ENDDO
c
c Si condensation negative, effort completement perdu:
c
      DO i = 1, klon
         IF (todo(i).AND.zcond(i).LE.0.) todo(i) = .FALSE.
      ENDDO
c
c L'ajustement a ete accompli, meme les calculs accessoires
c ne sont pas encore faits:
c
      DO i = 1, klon
         IF (todo(i)) accompli(i) = .TRUE.
      ENDDO
c
c=====
c Une fois que la condensation a lieu, on doit construire un
c "modele nuageux" pour partager la condensation entre l'eau
c liquide nuageuse et la precipitation (leur rapport toliq
c est calcule selon l'epaisseur nuageuse). Je suppose que
c toliq=tomax quand l'epaisseur nuageuse est inferieure a dpmin,
c et que toliq=tomin quand l'epaisseur depasse dpmax (interpolation
c lineaire entre dpmin et dpmax).
c=====
      DO i = 1, klon
      IF (todo(i)) THEN
      toliq(i) = tomax-((paprs(i,k1(i))-paprs(i,k2(i)+1))
     .                 /paprs(i,1)-dpmin)
     .                 *(tomax-tomin)/(dpmax-dpmin)
      toliq(i) = MAX(tomin,MIN(tomax,toliq(i)))
      IF (pplay(i,k2(i))/paprs(i,1) .LE. deep_sig) toliq(i) = deep_to
      IF (old_tau) toliq(i) = 1.0
      ENDIF
      ENDDO
c=====
c On doit aussi determiner la distribution verticale de 
c l'eau nuageuse. Plusieurs options sont proposees:
c
c (0) La condensation precipite integralement (toliq ne sera
c     pas utilise).
c (1) L'eau liquide est distribuee entre k1 et k2 et proportionnelle
c     a la vapeur d'eau locale.
c (2) Elle est distribuee entre k1 et k2 avec une valeur constante.
c (3) Elle est seulement distribuee aux couches ou la vapeur d'eau
c     est effectivement diminuee pendant le processus d'ajustement.
c (4) Elle est en fonction (lineaire ou exponentielle) de la
c     distance (epaisseur en pression) avec le niveau k1 (la couche
c     k1 n'aura donc pas d'eau liquide).
c=====
c
      IF (opt_cld.EQ.0) THEN
c
         DO i = 1, klon
            IF (todo(i)) zrfl(i) = zcond(i) / dtime
         ENDDO
c
      ELSE IF (opt_cld.EQ.1) THEN
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de vapeur d'eau
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
            IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i))
     .         zvapo(i) = zvapo(i) + zqnew(i,k)*delp(i,k)/RG
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
            zrapp(i) = toliq(i) * zcond(i) / zvapo(i)
            zrapp(i) = MAX(0.,MIN(1.,zrapp(i)))
            zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDIF
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i)) THEN
            d_ql(i,k) = d_ql(i,k) + zrapp(i) * zqnew(i,k)
         ENDIF
         ENDDO
         ENDDO
c
      ELSE IF (opt_cld.EQ.2) THEN
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
            IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i))
     .         zvapo(i) = zvapo(i) + delp(i,k)/RG
         ENDDO
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i)) THEN
            d_ql(i,k) = d_ql(i,k) + toliq(i) * zcond(i) / zvapo(i)
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
            IF (todo(i)) zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDDO
c
      ELSE IF (opt_cld.EQ.3) THEN
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite de l'eau strictement condensee
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
            IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i))
     .      zvapo(i) = zvapo(i) + MAX(0.0,zq(i,k)-zqnew(i,k)) 
     .                          * delp(i,k)/RG
         ENDDO
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i) .AND.
     .       zvapo(i).GT.0.0)
     .      d_ql(i,k) = d_ql(i,k) + toliq(i) * zcond(i) / zvapo(i)
     .                            * MAX(0.0,zq(i,k)-zqnew(i,k))
         ENDDO
         ENDDO
         DO i = 1, klon
            IF (todo(i)) zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDDO
c
      ELSE IF (opt_cld.EQ.4) THEN
c
         nexpo = 3
ccc         nexpo = 1 ! distribution lineaire
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse 
         ENDDO                       ! (avec ponderation)
         DO k = k1min, k2max
         DO i = 1, klon
            IF (todo(i) .AND. k.GE.(k1(i)+1) .AND. k.LE.k2(i))
     .         zvapo(i) = zvapo(i) + delp(i,k) / RG
     .                    * (pplay(i,k1(i))-pplay(i,k))**nexpo
         ENDDO
         ENDDO
         DO k = k1min, k2max
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.(k1(i)+1) .AND. k.LE.k2(i))
     .      d_ql(i,k) = d_ql(i,k) + toliq(i) * zcond(i) / zvapo(i)
     .                            * (pplay(i,k1(i))-pplay(i,k))**nexpo
         ENDDO
         ENDDO
         DO i = 1, klon
            IF (todo(i)) zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDDO
c
      ELSE ! valeur non-prevue pour opt_cld
c
         PRINT*, "opt_cld est faux:", opt_cld
         CALL abort
c
      ENDIF ! fin de opt_cld
c
c L'eau precipitante peut etre evaporee:
c
      zalfa = 0.05
      IF (evap_prec .AND. (k1max.GE.2)) THEN
      DO k = k1max-1, 1, -1
      DO i = 1, klon
      IF (todo(i) .AND. k.LT.k1(i) .AND. zrfl(i).GT.0.0) THEN
         zqev = MAX (0.0, (zqs(i,k)-zq(i,k))*zalfa )
         zqevt = coef_eva * (1.0-zq(i,k)/zqs(i,k))*SQRT(zrfl(i))
     .         * delp(i,k)/pplay(i,k)*zt(i,k)*RD/RG
         zqevt = MAX(0.0,MIN(zqevt,zrfl(i))) * RG*dtime/delp(i,k)
         zqev = MIN (zqev, zqevt)
         zrfln = zrfl(i) - zqev*(delp(i,k))/RG/dtime
         zq(i,k) = zq(i,k) - (zrfln-zrfl(i)) 
     .                     * (RG/(delp(i,k)))*dtime
         zt(i,k) = zt(i,k) + (zrfln-zrfl(i))
     .                     * (RG/(delp(i,k)))*dtime
     .                     * RLVTT/RCPD/(1.0+RVTMP2*zq(i,k))
         zrfl(i) = zrfln
      ENDIF
      ENDDO
      ENDDO
      ENDIF
c
c La temperature de la premiere couche determine la pluie ou la neige:
c
      DO i = 1, klon
      IF (todo(i)) THEN
      IF (zt(i,1) .GT. RTT) THEN
         rain(i) = rain(i) + zrfl(i)
      ELSE
         snow(i) = snow(i) + zrfl(i)
      ENDIF
      ENDIF
      ENDDO
c
c Mise a jour de la temperature et de l'humidite
c
      DO k = k1min, k2max
      DO i = 1, klon
      IF (todo(i) .AND. k.GE.k1(i) .AND. k.LE.k2(i)) THEN
         zt(i,k) = ztnew(i,k)
         zq(i,k) = zqnew(i,k)
      ENDIF
      ENDDO
      ENDDO
c
c Re-calculer certaines variables pour etendre et re-ajuster la colonne
c
      IF (exigent) THEN
      DO k = 1, klev
      DO i = 1, klon
      IF (todo(i)) THEN
         IF (thermcep) THEN
            zdelta=MAX(0.,SIGN(1.,RTT-zt(i,k)))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i,k))
            zqs(i,k)= R2ES*FOEEW(zt(i,k),zdelta)/pplay(i,k)
            zqs(i,k)=MIN(0.5,zqs(i,k))
            zcor=1./(1.-RETV*zqs(i,k))
            zqs(i,k)=zqs(i,k)*zcor
            zdqs(i,k) =FOEDE(zt(i,k),zdelta,zcvm5,zqs(i,k),zcor)
         ELSE
           IF (zt(i,k) .LT. t_coup) THEN
              zqs(i,k)= qsats(zt(i,k)) / pplay(i,k)
              zdqs(i,k)= dqsats(zt(i,k),zqs(i,k))
           ELSE
              zqs(i,k)= qsatl(zt(i,k)) / pplay(i,k)
              zdqs(i,k)= dqsatl(zt(i,k),zqs(i,k))
           ENDIF
         ENDIF
      ENDIF
      ENDDO
      ENDDO
      ENDIF
c
      IF (exigent) THEN
      DO k = 1, klev-1
      DO i = 1, klon
      IF (todo(i)) THEN
         zgamdz(i,k) = - (pplay(i,k)-pplay(i,k+1))/paprs(i,k+1)/RCPD
     .          *( RD*(zt(i,k)*delp(i,k)+zt(i,k+1)*delp(i,k+1))
     .               /(delp(i,k)+delp(i,k+1))
     .            +RLVTT*(zqs(i,k)*delp(i,k)+zqs(i,k+1)*delp(i,k+1))
     .                  /(delp(i,k)+delp(i,k+1))
     .           ) / (1.0+(zdqs(i,k)*delp(i,k)+zdqs(i,k+1)*delp(i,k+1))
     .                   /(delp(i,k)+delp(i,k+1)) )
      ENDIF
      ENDDO
      ENDDO
      ENDIF
c
c Puisque l'humidite a ete modifiee, on re-fait (q-qs)*dp
c
      DO k = 1, klev
      DO i = 1, klon
      IF (todo(i)) THEN
         zqmqsdp(i,k) = (zq(i,k)-zqs(i,k))*delp(i,k)
      ENDIF
      ENDDO
      ENDDO
c
c Verifier si l'on peut etendre le bas de la colonne
c
      DO i = 1, klon
         etendre(i) = .FALSE.
      ENDDO
c
      k1max = 1
      DO i = 1, klon
      IF (todo(i) .AND. k1(i).GT.(kbase+1)) THEN
         k = k1(i)
         zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
         zsat(i) = zqmqsdp(i,k) + zqmqsdp(i,k-1)
csc voici l'ancienne ligne:
csc         IF (zflo(i).LE.0.0 .OR. zsat(i).LE.0.0) THEN
csc sylvain: il faut RESPECTER les 2 criteres:
         IF (zflo(i).GT.0.0 .AND. zsat(i).GT.0.0) THEN
            etendre(i) = .TRUE.
            k1(i) = k1(i) - 1
            k1max = MAX(k1max,k1(i))
            aller(i) = .TRUE.
         ENDIF
      ENDIF
      ENDDO
c
      IF (k1max.GT.(kbase+1)) THEN
      DO k = k1max, kbase+1, -1
      DO i = 1, klon
      IF (etendre(i) .AND. k.LT.k1(i) .AND. aller(i)) THEN
         zsat(i) = zsat(i) + zqmqsdp(i,k)
         zflo(i) = zt(i,k) + zgamdz(i,k) - zt(i,k+1)
         IF (zsat(i).LE.0.0 .OR. zflo(i).LE.0.0) THEN
            aller(i) = .FALSE.
         ELSE
            k1(i) = k
         ENDIF
      ENDIF
      ENDDO
      ENDDO
      DO i = 1, klon
         IF (etendre(i).AND.aller(i)) THEN
            k1(i) = 1
         ENDIF
      ENDDO
      ENDIF
c
CCC      DO i = 1, klon
CCC      IF (etendre(i)) THEN
CCC  840    k = k1(i)
CCC         IF (k.GT.1) THEN
CCC            zsat(i) = zsat(i) + zqmqsdp(i,k-1)
CCC            zflo(i) = zt(i,k-1) + zgamdz(i,k-1) - zt(i,k)
CCC            IF (zflo(i).GT.0.0 .AND. zsat(i).GT.0.0) THEN
CCC               k1(i) = k - 1
CCC               GOTO 840
CCC            ENDIF
CCC         ENDIF
CCC      ENDIF
CCC      ENDDO
c
      DO i = 1, klon
         todobis(i) = todo(i)
         todo(i) = .FALSE.
      ENDDO
      is = 0
      DO i = 1, klon
      IF (etendre(i)) THEN
         todo(i) = .TRUE.
         is = is + 1
      ENDIF
      ENDDO
      IF (is.GT.0) THEN
         IF (new_top) THEN
            GOTO 820 ! chercher de nouveau le sommet k2
         ELSE
            GOTO 830 ! supposer que le sommet est celui deja trouve
         ENDIF
      ENDIF
c
      DO i = 1, klon
         possible(i) = .FALSE.
      ENDDO
      is = 0
      DO i = 1, klon
      IF (todobis(i) .AND. k2(i).LT.klev) THEN
         is = is + 1
         possible(i) = .TRUE.
      ENDIF
      ENDDO
      IF (is.GT.0) GOTO 810 !on cherche en haut d'autres blocks 
c     a ajuster a partir du sommet de la colonne precedente
c
  860 CONTINUE ! Calculer les tendances et diagnostiques
ccc      print*, "Apres 860"
c
      DO k = 1, klev
      DO i = 1, klon
      IF (accompli(i)) THEN
         d_t(i,k) = zt(i,k) - t(i,k)
         zq(i,k) = MAX(zq(i,k),seuil_vap)
         d_q(i,k) = zq(i,k) - q(i,k)
      ENDIF
      ENDDO
      ENDDO
c
      DO 888 i = 1, klon
      IF (accompli(i)) THEN
         DO k = 1, klev
         IF (rneb(i,k).GT.0.0) THEN
            ibas(i) = k
            GOTO 807
         ENDIF
         ENDDO
  807    CONTINUE
         DO k = klev, 1, -1
         IF (rneb(i,k).GT.0.0) THEN
            itop(i) = k
            GOTO 808
         ENDIF
         ENDDO
  808    CONTINUE
      ENDIF
  888 CONTINUE
c
      IF (imprim) THEN
         nbtodo = 0
         nbdone = 0
         DO i = 1, klon
            IF (afaire(i)) nbtodo = nbtodo + 1
            IF (accompli(i)) nbdone = nbdone + 1
         ENDDO
         PRINT*, "nbTodo, nbDone=", nbtodo, nbdone
      ENDIF
c
      RETURN
      END
      SUBROUTINE conkuo(dtime, paprs, pplay, t, q, conv_q,
     s                  d_t, d_q, d_ql, rneb, 
     s                  rain, snow, ibas, itop)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Schema de convection de type Kuo (1965).
c        Cette version du code peut calculer le niveau de depart
c N.B. version vectorielle (le 6 oct. 1997)
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c
c Arguments:
c
      REAL dtime               ! intervalle du temps (s)
      REAL paprs(klon,klev+1)  ! pression a inter-couche (Pa)
      REAL pplay(klon,klev)    ! pression au milieu de couche (Pa)
      REAL t(klon,klev)        ! temperature (K)
      REAL q(klon,klev)        ! humidite specifique
      REAL conv_q(klon,klev)   ! taux de convergence humidite (g/g/s)
c
      REAL d_t(klon,klev)      ! incrementation temperature
      REAL d_q(klon,klev)      ! incrementation humidite
      REAL d_ql(klon,klev)     ! incrementation eau liquide
      REAL rneb(klon,klev)     ! nebulosite
      REAL rain(klon)          ! pluies (mm/s)
      REAL snow(klon)          ! neige (mm/s)
      INTEGER itop(klon)       ! niveau du sommet
      INTEGER ibas(klon)       ! niveau du bas
c
      LOGICAL ldcum(klon)      ! convection existe
      LOGICAL todo(klon)
c
c Quelsques options:
c
      LOGICAL calcfcl ! calculer le niveau de convection libre
      PARAMETER (calcfcl=.TRUE.)
      INTEGER ldepar ! niveau fixe de convection libre
      PARAMETER (ldepar=4)
      INTEGER opt_cld ! comment traiter l'eau liquide
      PARAMETER (opt_cld=4) ! valeur possible: 0, 1, 2, 3 ou 4
      LOGICAL evap_prec ! evaporation de pluie au-dessous de convection
      PARAMETER (evap_prec=.TRUE.)
      REAL coef_eva
      PARAMETER (coef_eva=1.0E-05)
      LOGICAL new_deh ! nouvelle facon de calculer dH
      PARAMETER (new_deh=.FALSE.)
      REAL t_coup
      PARAMETER (t_coup=234.0)
      LOGICAL old_tau ! implique precipitation nulle
      PARAMETER (old_tau=.FALSE.)
      REAL toliq(klon) ! rapport entre l'eau nuageuse et l'eau precipitante
      REAL dpmin, tomax !Epaisseur faible, rapport eau liquide plus grande
      PARAMETER (dpmin=0.15, tomax=0.97)
      REAL dpmax, tomin !Epaisseur grande, rapport eau liquide plus faible
      PARAMETER (dpmax=0.30, tomin=0.05)
      REAL deep_sig, deep_to ! au dela de deep_sig, utiliser deep_to
      PARAMETER (deep_sig=0.50, deep_to=0.05)
c
c Variables locales:
c
      INTEGER nexpo
      LOGICAL nuage(klon)
      INTEGER i, k, kbmin, kbmax, khmax
      REAL ztotal(klon,klev), zdeh(klon,klev)
      REAL zgz(klon,klev)
      REAL zqs(klon,klev)
      REAL zdqs(klon,klev)
      REAL ztemp(klon,klev)
      REAL zpres(klon,klev)
      REAL zconv(klon) ! convergence d'humidite
      REAL zvirt(klon) ! convergence virtuelle d'humidite
      REAL zfrac(klon) ! fraction convective
      INTEGER kb(klon), kh(klon)
c
      REAL zcond(klon), zvapo(klon), zrapp(klon)
      REAL zrfl(klon), zrfln, zqev, zqevt
      REAL zdelta, zcvm5, zcor
      REAL zvar
c
      LOGICAL appel1er
      SAVE appel1er
c$OMP THREADPRIVATE(appel1er)
c
c Fonctions thermodynamiques
c
#include "YOETHF.h"
#include "FCTTRE.h"
c
      DATA appel1er /.TRUE./
c
      IF (appel1er) THEN
         PRINT*, 'conkuo, calcfcl:', calcfcl
         IF (.NOT.calcfcl) PRINT*, 'conkuo, ldepar:', ldepar
         PRINT*, 'conkuo, opt_cld:', opt_cld
         PRINT*, 'conkuo, evap_prec:', evap_prec
         PRINT*, 'conkuo, new_deh:', new_deh
         appel1er = .FALSE.
      ENDIF
c
c Initialiser les sorties a zero
c
      DO k = 1, klev
      DO i = 1, klon
         d_q(i,k) = 0.0
         d_t(i,k) = 0.0
         d_ql(i,k) = 0.0
         rneb(i,k) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
         ibas(i) = 0
         itop(i) = 0
      ENDDO
c
c Calculer la vapeur d'eau saturante Qs et sa derive L/Cp * dQs/dT
c
      DO k = 1, klev
      DO i = 1, klon
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-t(i,k)))
           zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*q(i,k))
           zqs(i,k)=R2ES*FOEEW(t(i,k),zdelta)/pplay(i,k)
           zqs(i,k)=MIN(0.5,zqs(i,k))
           zcor=1./(1.-RETV*zqs(i,k))
           zqs(i,k)=zqs(i,k)*zcor
           zdqs(i,k) =FOEDE(t(i,k),zdelta,zcvm5,zqs(i,k),zcor)
         ELSE
           IF (t(i,k).LT.t_coup) THEN
              zqs(i,k) = qsats(t(i,k))/pplay(i,k)
              zdqs(i,k) = dqsats(t(i,k),zqs(i,k))
           ELSE
              zqs(i,k) = qsatl(t(i,k))/pplay(i,k)
              zdqs(i,k) = dqsatl(t(i,k),zqs(i,k))
           ENDIF
         ENDIF
      ENDDO
      ENDDO
c
c Calculer gz (energie potentielle)
c
      DO i = 1, klon
         zgz(i,1) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1)))
     .                   * (paprs(i,1)-pplay(i,1))
      ENDDO
      DO k = 2, klev
      DO i = 1, klon
         zgz(i,k) = zgz(i,k-1)
     .            + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k)
     .                 * (pplay(i,k-1)-pplay(i,k))
      ENDDO
      ENDDO
c
c Calculer l'energie statique humide saturee (Cp*T + gz + L*Qs)
c
      DO k = 1, klev
      DO i = 1, klon
         ztotal(i,k) = RCPD*t(i,k) + RLVTT*zqs(i,k) + zgz(i,k)
      ENDDO
      ENDDO
c
c Determiner le niveau de depart et calculer la difference de
c l'energie statique humide saturee (ztotal) entre la couche
c de depart et chaque couche au-dessus.
c
      IF (calcfcl) THEN
         DO k = 1, klev
         DO i = 1, klon
            zpres(i,k) = pplay(i,k)
            ztemp(i,k) = t(i,k)
         ENDDO
         ENDDO
         CALL kuofcl(ztemp, q, zgz, zpres, ldcum, kb)
         DO i = 1, klon
         IF (ldcum(i)) THEN
            k = kb(i)
            IF (new_deh) THEN
            zdeh(i,k) = ztotal(i,k-1) - ztotal(i,k)
            ELSE
            zdeh(i,k) = RCPD * (t(i,k-1)-t(i,k))
     .                - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                  *(pplay(i,k-1)-pplay(i,k))
     .                + RLVTT*(zqs(i,k-1)-zqs(i,k))
            ENDIF
            zdeh(i,k) = zdeh(i,k) * 0.5
         ENDIF
         ENDDO
         DO k = 1, klev
         DO i = 1, klon
         IF (ldcum(i) .AND. k.GE.(kb(i)+1)) THEN
            IF (new_deh) THEN
               zdeh(i,k) = zdeh(i,k-1) + (ztotal(i,k-1)-ztotal(i,k))
            ELSE
               zdeh(i,k) = zdeh(i,k-1)
     .                   + RCPD * (t(i,k-1)-t(i,k))
     .                   - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                        *(pplay(i,k-1)-pplay(i,k))
     .                   + RLVTT*(zqs(i,k-1)-zqs(i,k))
            ENDIF
         ENDIF
         ENDDO
         ENDDO
      ELSE
         DO i = 1, klon
            k = ldepar
            kb(i) = ldepar
            ldcum(i) = .TRUE.
            IF (new_deh) THEN
            zdeh(i,k) = ztotal(i,k-1) - ztotal(i,k)
            ELSE
            zdeh(i,k) = RCPD * (t(i,k-1)-t(i,k))
     .                - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                  *(pplay(i,k-1)-pplay(i,k))
     .                + RLVTT*(zqs(i,k-1)-zqs(i,k))
            ENDIF
            zdeh(i,k) = zdeh(i,k) * 0.5
         ENDDO
         DO k = ldepar+1, klev
         DO i = 1, klon
         IF (new_deh) THEN
             zdeh(i,k)  = zdeh(i,k-1) + (ztotal(i,k-1)-ztotal(i,k))
         ELSE
             zdeh(i,k) = zdeh(i,k-1)
     .                 + RCPD * (t(i,k-1)-t(i,k))
     .                 - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                      *(pplay(i,k-1)-pplay(i,k))
     .                 + RLVTT*(zqs(i,k-1)-zqs(i,k))
         ENDIF
         ENDDO
         ENDDO
      ENDIF
c
c-----Chercher le sommet du nuage
c-----Calculer la convergence de l'humidite (en kg/m**2 a un facteur
c-----psolpa/RG pres) du bas jusqu'au sommet du nuage.
c-----Calculer la convergence virtuelle pour que toute la maille 
c-----deviennt nuageuse (du bas jusqu'au sommet du nuage)
c
      DO i = 1, klon
         nuage(i) = .TRUE.
         zconv(i) = 0.0
         zvirt(i) = 0.0
         kh(i) = -999
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
      IF (k.GE.kb(i) .AND. ldcum(i)) THEN
         nuage(i) = nuage(i) .AND. zdeh(i,k).GT.0.0
         IF (nuage(i)) THEN
            kh(i) = k
            zconv(i)=zconv(i)+conv_q(i,k)*dtime
     .                       *(paprs(i,k)-paprs(i,k+1))
            zvirt(i)=zvirt(i)+(zdeh(i,k)/RLVTT+zqs(i,k)-q(i,k))
     .                       *(paprs(i,k)-paprs(i,k+1))
         ENDIF
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
         todo(i) = ldcum(i) .AND. kh(i).GT.kb(i) .AND. zconv(i).GT.0.0
      ENDDO
c
      kbmin = klev
      kbmax = 0
      khmax = 0
      DO i = 1, klon
      IF (todo(i)) THEN
         kbmin = MIN(kbmin,kb(i))
         kbmax = MAX(kbmax,kb(i))
         khmax = MAX(khmax,kh(i))
      ENDIF
      ENDDO
c
c-----Calculer la surface couverte par le nuage
c
      DO i = 1, klon
      IF (todo(i)) THEN
      zfrac(i) = MAX(0.0,MIN(zconv(i)/zvirt(i), 1.0))
      ENDIF
      ENDDO
c
c-----Calculs essentiels:
c
      DO i = 1, klon
      IF (todo(i)) THEN
      zcond(i) = 0.0
      ENDIF
      ENDDO
      DO k = kbmin, khmax
      DO i = 1, klon
      IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
         zvar = zdeh(i,k)/(1.+zdqs(i,k))
         d_t(i,k) = zvar * zfrac(i) / RCPD
         d_q(i,k) = (zvar*zdqs(i,k)/RLVTT+zqs(i,k)-q(i,k))*zfrac(i)
     .            - conv_q(i,k)*dtime
         zcond(i) = zcond(i) - d_q(i,k) *(paprs(i,k)-paprs(i,k+1))/RG
         rneb(i,k) = zfrac(i)
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (todo(i) .AND. zcond(i).LT.0.0) THEN
         PRINT*, 'WARNING: cond. negative (Kuo) ',
     .            i,kb(i),kh(i), zcond(i)
         zcond(i) = 0.0
         DO k = kb(i), kh(i)
            d_t(i,k) = 0.0
            d_q(i,k) = 0.0
         ENDDO
         todo(i) = .FALSE. ! effort totalement perdu
      ENDIF
      ENDDO
c
c=====
c Une fois que la condensation a lieu, on doit construire un
c "modele nuageux" pour partager la condensation entre l'eau
c liquide nuageuse et la precipitation (leur rapport toliq
c est calcule selon l'epaisseur nuageuse). Je suppose que
c toliq=tomax quand l'epaisseur nuageuse est inferieure a dpmin,
c et que toliq=tomin quand l'epaisseur depasse dpmax (interpolation
c lineaire entre dpmin et dpmax).
c=====
      DO i = 1, klon
      IF (todo(i)) THEN
      toliq(i) = tomax-((paprs(i,kb(i))-paprs(i,kh(i)+1))
     .               /paprs(i,1)-dpmin)
     .             *(tomax-tomin)/(dpmax-dpmin)
      toliq(i) = MAX(tomin,MIN(tomax,toliq(i)))
      IF (pplay(i,kh(i))/paprs(i,1) .LE. deep_sig) toliq(i) = deep_to
      IF (old_tau) toliq(i) = 1.0
      ENDIF
      ENDDO
c=====
c On doit aussi determiner la distribution verticale de
c l'eau nuageuse. Plusieurs options sont proposees:
c
c (0) La condensation precipite integralement (toliq ne sera
c     pas utilise).
c (1) L'eau liquide est distribuee entre k1 et k2 et proportionnelle
c     a la vapeur d'eau locale.
c (2) Elle est distribuee entre k1 et k2 avec une valeur constante.
c (3) Elle est seulement distribuee aux couches ou la vapeur d'eau
c     est effectivement diminuee pendant le processus d'ajustement.
c (4) Elle est en fonction (lineaire ou exponentielle) de la
c     distance (epaisseur en pression) avec le niveau k1 (la couche
c     k1 n'aura donc pas d'eau liquide).
c=====
c
      IF (opt_cld.EQ.0) THEN
c
         DO i = 1, klon
         IF (todo(i)) zrfl(i) = zcond(i) / dtime
         ENDDO
c
      ELSE IF (opt_cld.EQ.1) THEN
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de vapeur d'eau
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
            zvapo(i) = zvapo(i) + (q(i,k)+d_q(i,k))
     .                     *(paprs(i,k)-paprs(i,k+1))/RG
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
            zrapp(i) = toliq(i) * zcond(i) / zvapo(i)
            zrapp(i) = MAX(0.,MIN(1.,zrapp(i)))
         ENDIF
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
            d_ql(i,k) = zrapp(i) * (q(i,k)+d_q(i,k))
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
         zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDIF
         ENDDO
c
      ELSE IF (opt_cld.EQ.2) THEN
c
         DO i = 1, klon
         IF (todo(i)) zvapo(i) = 0.0 ! quantite integrale de masse
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
            zvapo(i) = zvapo(i) + (paprs(i,k)-paprs(i,k+1))/RG
         ENDIF
         ENDDO
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
            d_ql(i,k) = toliq(i) * zcond(i) / zvapo(i)
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
         zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDIF
         ENDDO
c
      ELSE IF (opt_cld.EQ.3) THEN
c
         DO i = 1, klon
         IF (todo(i)) THEN
         zvapo(i) = 0.0 ! quantite de l'eau strictement condensee
         ENDIF
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i)) THEN
            zvapo(i) = zvapo(i) + MAX(0.0,-d_q(i,k))
     .                    * (paprs(i,k)-paprs(i,k+1))/RG
         ENDIF
         ENDDO
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.kb(i) .AND. k.LE.kh(i) .AND.
     .                     zvapo(i).GT.0.0) THEN
            d_ql(i,k) = d_ql(i,k) + toliq(i) * zcond(i) / zvapo(i)
     .                            * MAX(0.0,-d_q(i,k))
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
         zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDIF
         ENDDO
c
      ELSE IF (opt_cld.EQ.4) THEN
c
         nexpo = 3
ccc         nexpo = 1 ! distribution lineaire
c
         DO i = 1, klon
         IF (todo(i)) THEN
         zvapo(i) = 0.0 ! quantite integrale de masse (avec ponderation)
         ENDIF
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.(kb(i)+1) .AND. k.LE.kh(i)) THEN
            zvapo(i) = zvapo(i) + (paprs(i,k)-paprs(i,k+1)) / RG
     .                    * (pplay(i,kb(i))-pplay(i,k))**nexpo
         ENDIF
         ENDDO
         ENDDO
         DO k = kbmin, khmax
         DO i = 1, klon
         IF (todo(i) .AND. k.GE.(kb(i)+1) .AND. k.LE.kh(i)) THEN
            d_ql(i,k) = d_ql(i,k) + toliq(i) * zcond(i) / zvapo(i)
     .                            * (pplay(i,kb(i))-pplay(i,k))**nexpo
         ENDIF
         ENDDO
         ENDDO
         DO i = 1, klon
         IF (todo(i)) THEN
         zrfl(i) = (1.0-toliq(i)) * zcond(i) / dtime
         ENDIF
         ENDDO
c
      ELSE ! valeur non-prevue pour opt_cld
c
         PRINT*, "opt_cld est faux:", opt_cld
         CALL abort
c
      ENDIF ! fin de opt_cld
c
c L'eau precipitante peut etre re-evaporee:
c
      IF (evap_prec .AND. kbmax.GE.2) THEN
      DO k = kbmax, 1, -1
      DO i = 1, klon
      IF (todo(i) .AND. k.LE.(kb(i)-1) .AND. zrfl(i).GT.0.0) THEN
         zqev = MAX (0.0, (zqs(i,k)-q(i,k))*zfrac(i) )
         zqevt = coef_eva * (1.0-q(i,k)/zqs(i,k))*SQRT(zrfl(i))
     .       * (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*t(i,k)*RD/RG
         zqevt = MAX(0.0,MIN(zqevt,zrfl(i)))
     .         * RG*dtime/(paprs(i,k)-paprs(i,k+1))
         zqev = MIN (zqev, zqevt)
         zrfln = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1))
     .                 /RG/dtime
         d_q(i,k) = - (zrfln-zrfl(i))
     .          * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
         d_t(i,k) = (zrfln-zrfl(i))
     .          * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
     .          * RLVTT/RCPD
         zrfl(i) = zrfln
      ENDIF
      ENDDO
      ENDDO
      ENDIF
c
c La temperature de la premiere couche determine la pluie ou la neige:
c
      DO i = 1, klon
      IF (todo(i)) THEN
      IF (t(i,1) .GT. RTT) THEN
         rain(i) = rain(i) + zrfl(i)
      ELSE
         snow(i) = snow(i) + zrfl(i)
      ENDIF
      ENDIF
      ENDDO
c
      RETURN
      END
      SUBROUTINE kuofcl(pt, pq, pg, pp, LDCUM, kcbot)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19940927
c            adaptation du code de Tiedtke du ECMWF
c Objet: calculer le niveau de convection libre
c        (FCL: Free Convection Level)
c======================================================================
c Arguments:
c pt---input-R- temperature (K)
c pq---input-R- vapeur d'eau (kg/kg)
c pg---input-R- geopotentiel (g*z ou z est en metre)
c pp---input-R- pression (Pa)
c
c LDCUM---output-L- Y-t-il la convection
c kcbot---output-I- Niveau du bas de la convection
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "YOETHF.h"
C
      REAL pt(klon,klev), pq(klon,klev), pg(klon,klev), pp(klon,klev)
      INTEGER  kcbot(klon)
      LOGICAL  LDCUM(klon)
C
      REAL ztu(klon,klev), zqu(klon,klev), zlu(klon,klev)
      REAL zqold(klon), zbuo
      INTEGER is, i, k
c
c klab=1: on est sous le nuage convectif
c klab=2: le bas du nuage convectif
c klab=0: autres couches
      INTEGER klab(klon,klev)
c
c quand lflag=.true., on est sous le nuage, il faut donc appliquer
c le processus d'elevation.
      LOGICAL lflag(klon)
C
      DO k = 1, klev
      DO i = 1, klon
         ztu(i,k) = pt(i,k)
         zqu(i,k) = pq(i,k)
         zlu(i,k) = 0.0
         klab(i,k) = 0
      ENDDO
      ENDDO
C----------------------------------------------------------------------
      DO i = 1, klon
         klab(i,1)=1
         kcbot(i)=2
         LDCUM(i)=.FALSE.
      ENDDO
C
      DO 290 k = 2, klev-1
c
      is=0
      DO i = 1, klon
         if (klab(i,k-1).EQ.1) is = is + 1
         lflag(i) = .FALSE.
         if (klab(i,k-1).EQ.1) lflag(i) = .TRUE.
      ENDDO
      IF (is.EQ.0) GOTO 290
c
c on eleve le parcel d'air selon l'adiabatique sec
c
      DO i = 1, klon
      IF (lflag(i)) THEN
         zqu(i,k) = zqu(i,k-1)
         ztu(i,k) = ztu(i,k-1) + (pg(i,k-1)-pg(i,k))/RCPD
         zbuo = ztu(i,k)*(1.+RETV*zqu(i,k))-
     .          pt(i,k)*(1.+RETV*pq(i,k))+0.5
         IF (zbuo.GT.0.) klab(i,k)=1
         zqold(i) = zqu(i,k)
      ENDIF
      ENDDO
c
c on calcule la condensation eventuelle
c
      CALL adjtq(pp(1,k), ztu(1,k), zqu(1,k), lflag, 1)
c
c s'il y a la condensation et la "buoyancy" force est positive
c c'est bien le bas de la tour de convection
c
      DO i=1, klon
      IF(lflag(i).AND.zqu(i,k).NE.zqold(i)) THEN
         klab(i,k) = 2
         zlu(i,k) = zlu(i,k)+zqold(i)-zqu(i,k)
         zbuo = ztu(i,k)*(1.+RETV*zqu(i,k))-
     .          pt(i,k)*(1.+RETV*pq(i,k))+0.5
         IF (zbuo.GT.0.) THEN
            kcbot(i) = k
            LDCUM(i) = .TRUE.
         ENDIF
      ENDIF
      ENDDO
C
  290 CONTINUE
C
      RETURN
      END
      SUBROUTINE adjtq(pp, pt, pq, LDFLAG, KCALL)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19940927
c            adaptation du code de Tiedtke du ECMWF
c Objet: ajustement entre T et Q
c======================================================================
c Arguments:
c pp---input-R- pression (Pa)
c pt---input/output-R- temperature (K)
c pq---input/output-R- vapeur d'eau (kg/kg)
c======================================================================
C TO PRODUCE T,Q AND L VALUES FOR CLOUD ASCENT
C
C NOTE: INPUT PARAMETER KCALL DEFINES CALCULATION AS
C        KCALL=0    ENV. T AND QS IN*CUINI*
C        KCALL=1  CONDENSATION IN UPDRAFTS  (E.G. CUBASE, CUASC)
C        KCALL=2  EVAPORATION IN DOWNDRAFTS (E.G. CUDLFS,CUDDRAF)
C
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
C
      REAL pt(klon), pq(klon), pp(klon)
      LOGICAL  ldflag(klon)
      INTEGER KCALL
c
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
      REAL zcond(klon), zcond1
      REAL zdelta, zcvm5, zldcp, zqsat, zcor, zdqsat
      INTEGER is, i
#include "YOETHF.h"
#include "FCTTRE.h"
c
      DO i = 1, klon
         zcond(i) = 0.0
      ENDDO
C
      DO 210 i=1, klon
      IF (LDFLAG(i)) THEN
         zdelta=MAX(0.,SIGN(1.,RTT-pt(i)))
         zldcp = RLVTT*(1.-zdelta) + zdelta*RLSTT
         zldcp = zldcp / RCPD/(1.0+RVTMP2*pq(i))
         IF (thermcep) THEN
           zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
           zcvm5 = zcvm5 / RCPD/(1.0+RVTMP2*pq(i))
           zqsat=R2ES*FOEEW (pt(i), zdelta) / pp(i)
           zqsat=MIN(0.5,zqsat)
           zcor=1./(1.-RETV  *zqsat)
           zqsat=zqsat*zcor
           zdqsat = FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor)
         ELSE
           IF (pt(i).LT.t_coup) THEN
              zqsat = qsats(pt(i))/pp(i)
              zdqsat = dqsats(pt(i),zqsat)
           ELSE
              zqsat = qsatl(pt(i))/pp(i)
              zdqsat = dqsatl(pt(i),zqsat)
           ENDIF
         ENDIF
         zcond(i)=(pq(i)-zqsat) / (1. + zdqsat)
         IF(KCALL.EQ.1) zcond(i)=MAX(zcond(i),0.)
         IF(KCALL.EQ.2) zcond(i)=MIN(zcond(i),0.)
         pt(i)=pt(i)+zldcp*zcond(i)
         pq(i)=pq(i)-zcond(i)
      ENDIF
  210 CONTINUE
C
      is = 0
      DO i=1, klon
         if (zcond(i).NE.0.) is = is + 1
      ENDDO
      IF(is.EQ.0) GOTO 230
C
      DO 220 i = 1, klon
      IF(LDFLAG(i).AND.zcond(i).NE.0.) THEN
         zdelta=MAX(0.,SIGN(1.,RTT-pt(i)))
         zldcp = RLVTT*(1.-zdelta) + zdelta*RLSTT
         zldcp = zldcp / RCPD/(1.0+RVTMP2*pq(i))
         IF (thermcep) THEN
           zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
           zcvm5 = zcvm5 / RCPD/(1.0+RVTMP2*pq(i))
           zqsat=R2ES*FOEEW (pt(i), zdelta) / pp(i)
           zqsat=MIN(0.5,zqsat)
           zcor=1./(1.-RETV  *zqsat)
           zqsat=zqsat*zcor
           zdqsat = FOEDE(pt(i), zdelta, zcvm5, zqsat, zcor)
         ELSE
           IF (pt(i).LT.t_coup) THEN
              zqsat = qsats(pt(i))/pp(i)
              zdqsat = dqsats(pt(i),zqsat)
           ELSE
              zqsat = qsatl(pt(i))/pp(i)
              zdqsat = dqsatl(pt(i),zqsat)
           ENDIF
         ENDIF
         zcond1=(pq(i)-zqsat) / (1.+zdqsat)
         pt(i)=pt(i)+zldcp*zcond1
         pq(i)=pq(i)-zcond1
      END IF
  220 CONTINUE
C
  230 CONTINUE
      RETURN
      END
      SUBROUTINE fiajh(dtime, paprs, pplay, t, q,
     .                 d_t, d_q, d_ql, rneb,
     .                 rain, snow, ibas, itop)
      USE dimphy
      IMPLICIT NONE
c
c Ajustement humide (Schema de convection de Manabe)
C.
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c
c Arguments:
c
      REAL dtime        ! intervalle du temps (s)
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (kg/kg)
      REAL paprs(klon,klev+1) ! pression a inter-couche (Pa)
      REAL pplay(klon,klev) ! pression au milieu de couche (Pa)
c
      REAL d_t(klon,klev) ! incrementation pour la temperature
      REAL d_q(klon,klev) ! incrementation pour vapeur d'eau
      REAL d_ql(klon,klev) ! incrementation pour l'eau liquide
      REAL rneb(klon,klev) ! fraction nuageuse
c
      REAL rain(klon)    ! variable non utilisee
      REAL snow(klon)    ! variable non utilisee
      INTEGER ibas(klon) ! variable non utilisee
      INTEGER itop(klon) ! variable non utilisee

      REAL t_coup
      PARAMETER (t_coup=234.0)
      REAL seuil_vap
      PARAMETER (seuil_vap=1.0E-10)
c
c Variables locales:
c 
      INTEGER i, k
      INTEGER k1, k1p, k2, k2p
      LOGICAL itest(klon)
      REAL delta_q(klon, klev)
      REAL cp_new_t(klev)
      REAL cp_delta_t(klev)
      REAL new_qb(klev)
      REAL v_cptj(klev), v_cptjk1, v_ssig
      REAL v_cptt(klon,klev), v_p, v_t
      REAL v_qs(klon,klev), v_qsd(klon,klev)
      REAL zq1(klon), zq2(klon)
      REAL gamcpdz(klon,2:klev)
      REAL zdp, zdpm
c
      REAL zsat ! sur-saturation
      REAL zflo ! flotabilite
c
      REAL local_q(klon,klev),local_t(klon,klev)
c
      REAL zdelta, zcor, zcvm5
C
#include "YOETHF.h"
#include "FCTTRE.h"
C
      DO k = 1, klev
      DO i = 1, klon
         local_q(i,k) = q(i,k)
         local_t(i,k) = t(i,k)
         rneb(i,k) = 0.0
         d_ql(i,k) = 0.0
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
         ibas(i) = 0
         itop(i) = 0
      ENDDO
c
c Calculer v_qs et v_qsd:
c
      DO k = 1, klev
      DO i = 1, klon
         v_cptt(i,k) = RCPD * local_t(i,k)
         v_t = local_t(i,k)
         v_p = pplay(i,k)
c
         IF (thermcep) THEN
            zdelta=MAX(0.,SIGN(1.,RTT-v_t))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*local_q(i,k))
            v_qs(i,k)= R2ES * FOEEW(v_t,zdelta)/v_p
            v_qs(i,k)=MIN(0.5,v_qs(i,k))
            zcor=1./(1.-RETV*v_qs(i,k))
            v_qs(i,k)=v_qs(i,k)*zcor
            v_qsd(i,k) =FOEDE(v_t,zdelta,zcvm5,v_qs(i,k),zcor)
         ELSE
           IF (v_t.LT.t_coup) THEN
              v_qs(i,k) = qsats(v_t) / v_p
              v_qsd(i,k) = dqsats(v_t,v_qs(i,k))
           ELSE
              v_qs(i,k) = qsatl(v_t) / v_p
              v_qsd(i,k) = dqsatl(v_t,v_qs(i,k))
           ENDIF
         ENDIF
      ENDDO
      ENDDO
c
c Calculer Gamma * Cp * dz: (gamm est le gradient critique)
c
      DO k = 2, klev
      DO i = 1, klon
         zdp = paprs(i,k)-paprs(i,k+1)
         zdpm = paprs(i,k-1)-paprs(i,k)
         gamcpdz(i,k) = ( ( RD/RCPD /(zdpm+zdp) * 
     .                      (v_cptt(i,k-1)*zdpm + v_cptt(i,k)*zdp)
     .                     +RLVTT /(zdpm+zdp) * 
     .                      (v_qs(i,k-1)*zdpm + v_qs(i,k)*zdp)
     .                    )* (pplay(i,k-1)-pplay(i,k)) / paprs(i,k) )
     .                / (1.0+(v_qsd(i,k-1)*zdpm+
     .                        v_qsd(i,k)*zdp)/(zdpm+zdp) )
      ENDDO
      ENDDO
C
C------------------------------------ modification des profils instables
      DO 9999 i = 1, klon
      itest(i) = .FALSE.
C
      k1 = 0
      k2 = 1
C
  810 CONTINUE ! chercher k1, le bas de la colonne
      k2 = k2 + 1
      IF (k2 .GT. klev) GOTO 9999
      zflo = v_cptt(i,k2-1) - v_cptt(i,k2) - gamcpdz(i,k2)
      zsat=(local_q(i,k2-1)-v_qs(i,k2-1))*(paprs(i,k2-1)-paprs(i,k2))
     .    +(local_q(i,k2)-v_qs(i,k2))*(paprs(i,k2)-paprs(i,k2+1))
      IF ( zflo.LE.0.0 .OR. zsat.LE.0.0 ) GOTO 810
      k1 = k2 - 1
      itest(i) = .TRUE.
C
  820 CONTINUE ! chercher k2, le haut de la colonne
      IF (k2 .EQ. klev) GOTO 821
      k2p = k2 + 1
      zsat=zsat +(paprs(i,k2p)-paprs(i,k2p+1))
     .          *(local_q(i,k2p)-v_qs(i,k2p))
      zflo = v_cptt(i,k2p-1) - v_cptt(i,k2p) - gamcpdz(i,k2p)
      IF (zflo.LE.0.0 .OR. zsat.LE.0.0) GOTO 821
      k2 = k2p
      GOTO 820
  821 CONTINUE
C
C------------------------------------------------------ ajustement local
  830 CONTINUE ! ajustement proprement dit
      v_cptj(k1) = 0.0
      zdp = paprs(i,k1)-paprs(i,k1+1)
      v_cptjk1 = ( (1.0+v_qsd(i,k1))*(v_cptt(i,k1)+v_cptj(k1))
     .               + RLVTT*(local_q(i,k1)-v_qs(i,k1)) ) * zdp
      v_ssig = zdp * (1.0+v_qsd(i,k1))
C
      k1p = k1 + 1
      DO k = k1p, k2
         zdp = paprs(i,k)-paprs(i,k+1)
         v_cptj(k) = v_cptj(k-1) + gamcpdz(i,k)
         v_cptjk1 = v_cptjk1 + zdp
     .             * ( (1.0+v_qsd(i, k))*(v_cptt(i,k)+v_cptj(k))
     .               + RLVTT*(local_q(i,k)-v_qs(i,k)) )
         v_ssig = v_ssig + zdp *(1.0+v_qsd(i,k))
      ENDDO
C
      DO k = k1, k2
         cp_new_t(k) = v_cptjk1/v_ssig - v_cptj(k)
         cp_delta_t(k) = cp_new_t(k) - v_cptt(i,k)
         new_qb(k) = v_qs(i,k) + v_qsd(i,k)*cp_delta_t(k)/RLVTT
         local_q(i,k) = new_qb(k)
         local_t(i,k) = cp_new_t(k) / RCPD
      ENDDO
C
C--------------------------------------------------- sondage vers le bas
C              -- on redefinit les variables prognostiques dans
C              -- la colonne qui vient d'etre ajustee
C
      DO k = k1, k2
         v_cptt(i,k) = RCPD * local_t(i,k)
         v_t = local_t(i,k)
         v_p = pplay(i,k)
C
         IF (thermcep) THEN
            zdelta=MAX(0.,SIGN(1.,RTT-v_t))
            zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
            zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*local_q(i,k))
            v_qs(i,k)= R2ES * FOEEW(v_t,zdelta)/v_p
            v_qs(i,k)=MIN(0.5,v_qs(i,k))
            zcor=1./(1.-RETV*v_qs(i,k))
            v_qs(i,k)=v_qs(i,k)*zcor
            v_qsd(i,k) =FOEDE(v_t,zdelta,zcvm5,v_qs(i,k),zcor)
         ELSE
           IF (v_t.LT.t_coup) THEN
              v_qs(i,k) = qsats(v_t) / v_p
              v_qsd(i,k) = dqsats(v_t,v_qs(i,k))
           ELSE
              v_qs(i,k) = qsatl(v_t) / v_p
              v_qsd(i,k) = dqsatl(v_t,v_qs(i,k))
           ENDIF
         ENDIF
      ENDDO
      DO k = 2, klev
         zdpm = paprs(i,k-1) - paprs(i,k)
         zdp = paprs(i,k) - paprs(i,k+1)
         gamcpdz(i,k) = ( ( RD/RCPD /(zdpm+zdp) *
     .                      (v_cptt(i,k-1)*zdpm+v_cptt(i,k)*zdp)
     .                     +RLVTT /(zdpm+zdp) *
     .                      (v_qs(i,k-1)*zdpm+v_qs(i,k)*zdp)
     .                    )* (pplay(i,k-1)-pplay(i,k)) / paprs(i,k) )
     .                / (1.0+(v_qsd(i,k-1)*zdpm+v_qsd(i,k)*zdp)
     .                      /(zdpm+zdp) )
      ENDDO
C
C Verifier si l'on peut etendre la colonne vers le bas
C
      IF (k1 .EQ. 1) GOTO 841 ! extension echouee
      zflo = v_cptt(i,k1-1) - v_cptt(i,k1) - gamcpdz(i,k1)
      zsat=(local_q(i,k1-1)-v_qs(i,k1-1))*(paprs(i,k1-1)-paprs(i,k1))
     .   + (local_q(i,k1)-v_qs(i,k1))*(paprs(i,k1)-paprs(i,k1+1))
      IF (zflo.LE.0.0 .OR. zsat.LE.0.0) GOTO 841 ! extension echouee
C
  840 CONTINUE
      k1 = k1 - 1
      IF (k1 .EQ. 1) GOTO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
      zsat = zsat + (local_q(i,k1-1)-v_qs(i,k1-1))
     .             *(paprs(i,k1-1)-paprs(i,k1))
      zflo = v_cptt(i,k1-1) - v_cptt(i,k1) - gamcpdz(i,k1)
      IF (zflo.GT.0.0 .AND. zsat.GT.0.0) THEN
         GOTO 840
      ELSE
         GOTO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
      ENDIF
  841 CONTINUE
C
      GOTO 810 ! chercher d'autres blocks en haut
C
 9999 CONTINUE ! boucle sur tous les points
C-----------------------------------------------------------------------
c
c Determiner la fraction nuageuse (hypothese: la nebulosite a lieu
c a l'endroit ou la vapeur d'eau est diminuee par l'ajustement):
c
      DO k = 1, klev
      DO i = 1, klon
         IF (itest(i)) THEN
         delta_q(i,k) = local_q(i,k) - q(i,k)
         IF (delta_q(i,k).LT.0.) rneb(i,k)  = 1.0
         ENDIF
      ENDDO
      ENDDO
c
c Distribuer l'eau condensee en eau liquide nuageuse (hypothese:
c l'eau liquide est distribuee aux endroits ou la vapeur d'eau
c diminue et d'une maniere proportionnelle a cet diminution):
c
      DO i = 1, klon
         IF (itest(i)) THEN
         zq1(i) = 0.0
         zq2(i) = 0.0
         ENDIF
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         IF (itest(i)) THEN
         zdp = paprs(i,k)-paprs(i,k+1)
         zq1(i) = zq1(i) - delta_q(i,k) * zdp
         zq2(i) = zq2(i) - MIN(0.0, delta_q(i,k)) * zdp
         ENDIF
      ENDDO
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         IF (itest(i)) THEN
         IF (zq2(i).NE.0.0)
     .      d_ql(i,k) = - MIN(0.0,delta_q(i,k))*zq1(i)/zq2(i)
         ENDIF
      ENDDO
      ENDDO
C
      DO k = 1, klev
      DO i = 1, klon
          local_q(i, k) = MAX(local_q(i, k), seuil_vap)
      ENDDO
      ENDDO
C
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = local_t(i,k) - t(i,k)
         d_q(i,k) = local_q(i,k) - q(i,k)
      ENDDO
      ENDDO
c
      RETURN
      END
      SUBROUTINE fiajc(dtime,paprs,pplay,
     .                 t, q,conv_q,
     .                 d_t, d_q, d_ql,rneb,
     .                 rain, snow, ibas, itop)
      USE dimphy
      IMPLICIT NONE
c
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
c
c Options:
c
      INTEGER plb ! niveau de depart pour la convection
      PARAMETER (plb=4)
c
c Mystere: cette option n'est pas innocente pour les resultats !
c Qui peut resoudre ce mystere ? (Z.X.Li mars 1995)
      LOGICAL vector ! calcul vectorise
      PARAMETER (vector=.FALSE.)
c
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
c Arguments:
c
      REAL q(klon,klev) ! humidite specifique (kg/kg)
      REAL t(klon,klev) ! temperature (K)
      REAL paprs(klon,klev+1) ! pression a inter-couche (Pa)
      REAL pplay(klon,klev) ! pression au milieu de couche (Pa)
      REAL dtime ! intervalle du temps (s)
      REAL conv_q(klon,klev) ! taux de convergence de l'humidite
      REAL rneb(klon,klev) ! fraction nuageuse
      REAL d_q(klon,klev) ! incrementaion pour la vapeur d'eau
      REAL d_ql(klon,klev) ! incrementation pour l'eau liquide
      REAL d_t(klon,klev) ! incrementation pour la temperature
      REAL rain(klon) ! variable non-utilisee
      REAL snow(klon) ! variable non-utilisee
      INTEGER itop(klon) ! variable non-utilisee
      INTEGER ibas(klon) ! variable non-utilisee
c
      INTEGER kh(klon), i, k
      LOGICAL nuage(klon), test(klon,klev)
      REAL zconv(klon), zdeh(klon,klev), zvirt(klon)
      REAL zdqs(klon,klev), zqs(klon,klev)
      REAL ztt, zvar, zfrac(klon)
      REAL zq1(klon), zq2(klon)
      REAL zdelta, zcor, zcvm5
C
#include "YOETHF.h"
#include "FCTTRE.h"
c
c Initialiser les sorties:
c
      DO k = 1, klev
      DO i = 1, klon
          rneb(i,k) = 0.0
          d_ql(i,k) = 0.0
          d_t(i,k) = 0.0
          d_q(i,k) = 0.0
      ENDDO
      ENDDO
      DO i = 1, klon
         itop(i) = 0
         ibas(i) = 0
         rain(i) = 0.0
         snow(i) = 0.0
      ENDDO
c
c Calculer Qs et L/Cp * dQs/dT:
c
      DO k = 1, klev
      DO i = 1, klon
         ztt = t(i,k)
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-ztt))
           zcvm5=R5LES*RLVTT*(1.-zdelta) + zdelta*R5IES*RLSTT
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*q(i,k))
           zqs(i,k)= R2ES*FOEEW(ztt,zdelta)/pplay(i,k)
           zqs(i,k)=MIN(0.5,zqs(i,k))
           zcor=1./(1.-RETV*zqs(i,k))
           zqs(i,k)=zqs(i,k)*zcor
           zdqs(i,k) =FOEDE(ztt,zdelta,zcvm5,zqs(i,k),zcor)
         ELSE
           IF (ztt .LT. t_coup) THEN
              zqs(i,k) = qsats(ztt) / pplay(i,k)
              zdqs(i,k) = dqsats(ztt,zqs(i,k))
           ELSE
              zqs(i,k) = qsatl(ztt) / pplay(i,k)
              zdqs(i,k) = dqsatl(ztt,zqs(i,k))
           ENDIF
         ENDIF
      ENDDO
      ENDDO
c
c Determiner la difference de l'energie totale saturee:
c
      DO i = 1, klon
         k = plb
         zdeh(i,k) = RCPD * (t(i,k-1)-t(i,k))
     .             - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                  *(pplay(i,k-1)-pplay(i,k))
     .             + RLVTT*(zqs(i,k-1)-zqs(i,k))
         zdeh(i,k) = zdeh(i,k) * 0.5 ! on prend la moitie
      ENDDO
      DO k = plb+1, klev
      DO i = 1, klon
      zdeh(i,k) = zdeh(i,k-1)
     .             + RCPD * (t(i,k-1)-t(i,k))
     .             - RD *0.5*(t(i,k-1)+t(i,k))/paprs(i,k)
     .                  *(pplay(i,k-1)-pplay(i,k))
     .             + RLVTT*(zqs(i,k-1)-zqs(i,k))
      ENDDO
      ENDDO
c
c Determiner le sommet du nuage selon l'instabilite
c Calculer les convergences d'humidite (reelle et virtuelle)
c
      DO i = 1, klon
         nuage(i) = .TRUE.
         zconv(i) = 0.0
         zvirt(i) = 0.0
         kh(i) = -999
      ENDDO
      DO k = plb, klev
      DO i = 1, klon
         nuage(i) = nuage(i) .AND. zdeh(i,k).GT.0.0
         IF (nuage(i)) THEN
            kh(i)  = k
            zconv(i) = zconv(i)+conv_q(i,k)*dtime
     .                         *(paprs(i,k)-paprs(i,k+1))
            zvirt(i)=zvirt(i)+(zdeh(i,k)/RLVTT+zqs(i,k)-q(i,k))
     .                       *(paprs(i,k)-paprs(i,k+1))
         ENDIF
      ENDDO
      ENDDO
c
      IF (vector) THEN
c
c
      DO k = plb, klev
      DO i = 1, klon
      IF (k.LE.kh(i) .AND. kh(i).GT.plb .AND. zconv(i).GT.0.0) THEN
         test(i,k) = .TRUE.
         zfrac(i) = MAX(0.0,MIN(zconv(i)/zvirt(i),1.0))
      ELSE
         test(i,k) = .FALSE.
      ENDIF
      ENDDO
      ENDDO
c
      DO k = plb, klev
      DO i = 1, klon
      IF (test(i,k)) THEN
         zvar = zdeh(i,k)/(1.0+zdqs(i,k))
         d_q(i,k) = (zvar*zdqs(i,k)/RLVTT+zqs(i,k)-q(i,k))*zfrac(i)
     .            - conv_q(i,k)*dtime
         d_t(i,k) = zvar * zfrac(i) / RCPD
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
         zq1(i) = 0.0
         zq2(i) = 0.0
      ENDDO
      DO k = plb, klev
      DO i = 1, klon
      IF (test(i,k)) THEN
         IF (d_q(i,k).LT.0.0) rneb(i,k) = zfrac(i)
         zq1(i) = zq1(i) - d_q(i,k) * (paprs(i,k)-paprs(i,k+1))
         zq2(i) = zq2(i) - MIN(0.0, d_q(i,k))
     .                   * (paprs(i,k)-paprs(i,k+1))
      ENDIF
      ENDDO
      ENDDO
c
      DO k = plb, klev
      DO i = 1, klon
      IF (test(i,k)) THEN
         IF(zq2(i).NE.0.)d_ql(i,k)=-MIN(0.0,d_q(i,k))*zq1(i)/zq2(i)
      ENDIF
      ENDDO
      ENDDO
c
      ELSE ! (.NOT. vector)
c
      DO 999 i = 1, klon
      IF (kh(i).GT.plb .AND. zconv(i).GT.0.0) THEN
ccc         IF (kh(i).LE.plb) GOTO 999 ! il n'y a pas d'instabilite
ccc         IF (zconv(i).LE.0.0) GOTO 999 ! convergence insuffisante
         zfrac(i)  = MAX(0.0,MIN(zconv(i)/zvirt(i),1.0))
         DO k = plb, kh(i)
            zvar = zdeh(i,k)/(1.0+zdqs(i,k))
            d_q(i,k) = (zvar*zdqs(i,k)/RLVTT+zqs(i,k)-q(i,k))*zfrac(i)
     .               - conv_q(i,k)*dtime
            d_t(i,k) = zvar * zfrac(i) / RCPD
         ENDDO
c
         zq1(i) = 0.0
         zq2(i) = 0.0
         DO k = plb, kh(i)
            IF (d_q(i,k).LT.0.0) rneb(i,k) = zfrac(i)
            zq1(i) = zq1(i) - d_q(i,k) * (paprs(i,k)-paprs(i,k+1))
            zq2(i) = zq2(i) - MIN(0.0, d_q(i,k))
     .                      * (paprs(i,k)-paprs(i,k+1))
         ENDDO
         DO k = plb, kh(i)
            IF(zq2(i).NE.0.)d_ql(i,k)=-MIN(0.0,d_q(i,k))*zq1(i)/zq2(i)
         ENDDO
      ENDIF
  999 CONTINUE
c
      ENDIF ! fin de teste sur vector
c
      RETURN
      END
