!
! $Id: fisrtilp_tr.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
      SUBROUTINE fisrtilp_tr(dtime,paprs,pplay,t,q,ratqs,
     s                   d_t, d_q, d_ql, rneb, radliq, rain, snow,
     s                   pfrac_impa, pfrac_nucl, pfrac_1nucl,
     s                   frac_impa, frac_nucl,
     s                   prfl, psfl,
     s                   RHcl) ! relative humidity in clear sky (needed for aer optical properties; aeropt.F)

c
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c Date: le 20 mars 1995
c Objet: condensation et precipitation stratiforme.
c        schema de nuage
c======================================================================
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
#include "tracstoke.h"
c
c Arguments:
c
      REAL dtime ! intervalle du temps (s)
      REAL paprs(klon,klev+1) ! pression a inter-couche
      REAL pplay(klon,klev) ! pression au milieu de couche
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! humidite specifique (kg/kg)
      REAL d_t(klon,klev) ! incrementation de la temperature (K)
      REAL d_q(klon,klev) ! incrementation de la vapeur d'eau
      REAL d_ql(klon,klev) ! incrementation de l'eau liquide
      REAL rneb(klon,klev) ! fraction nuageuse
      REAL radliq(klon,klev) ! eau liquide utilisee dans rayonnements
      REAL rain(klon) ! pluies (mm/s)
      REAL snow(klon) ! neige (mm/s)
      REAL prfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s)
      REAL psfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s)
      
Cjq   For aerosol opt properties needed (see aeropt.F)
      REAL RHcl(klon,klev)
      
cAA
c Coeffients de fraction lessivee : pour OFF-LINE
c
      REAL pfrac_nucl(klon,klev)
      REAL pfrac_1nucl(klon,klev)
      REAL pfrac_impa(klon,klev)
c
c Fraction d'aerosols lessivee par impaction et par nucleation
c POur ON-LINE
c
      REAL frac_impa(klon,klev)
      REAL frac_nucl(klon,klev)
cAA
c
c Options du programme:
c
      REAL seuil_neb ! un nuage existe vraiment au-dela
      PARAMETER (seuil_neb=0.001)
      REAL ct ! inverse du temps pour qu'un nuage precipite
      PARAMETER (ct=1./1800.)
      REAL cl ! seuil de precipitation
      PARAMETER (cl=2.6e-4)
ccc      PARAMETER (cl=2.3e-4)
ccc      PARAMETER (cl=2.0e-4)
      INTEGER ninter ! sous-intervals pour la precipitation
      PARAMETER (ninter=5)
      LOGICAL evap_prec ! evaporation de la pluie
      PARAMETER (evap_prec=.TRUE.)
      REAL coef_eva
      PARAMETER (coef_eva=2.0E-05)
      LOGICAL calcrat ! calculer ratqs au lieu de fixer sa valeur
      REAL ratqs(klon,klev) ! determine la largeur de distribution de vapeur
      PARAMETER (calcrat=.TRUE.)
      REAL zx_min, rat_max
      PARAMETER (zx_min=1.0, rat_max=0.01)
      REAL zx_max, rat_min
      PARAMETER (zx_max=0.1, rat_min=0.3)
      REAL zx
c
      LOGICAL cpartiel ! condensation partielle
      PARAMETER (cpartiel=.TRUE.)
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
c Variables locales:
c
      INTEGER i, k, n, kk
      REAL zqs(klon), zdqs(klon), zdelta, zcor, zcvm5
      REAL zrfl(klon), zrfln(klon), zqev, zqevt
      REAL zoliq(klon), zcond(klon), zq(klon), zqn(klon), zdelq
      REAL ztglace, zt(klon)
      INTEGER nexpo ! exponentiel pour glace/eau
      REAL zdz(klon),zrho(klon),ztot(klon), zrhol(klon)
      REAL zchau(klon),zfroi(klon),zfice(klon),zneb(klon)
c
      LOGICAL appel1er
      SAVE appel1er
c$OMP THREADPRIVATE(appel1er)
c
c---------------------------------------------------------------
c
cAA Variables traceurs:
cAA  Provisoire !!! Parametres alpha du lessivage
cAA  A priori on a 4 scavenging # possibles
c
      REAL a_tr_sca(4)
      save a_tr_sca
c$OMP THREADPRIVATE(a_tr_sca)
c
c Variables intermediaires
c
      REAL zalpha_tr
      REAL zfrac_lessi
      REAL zprec_cond(klon)
cAA
c---------------------------------------------------------------
c
c Fonctions en ligne:
c
      REAL fallv ! vitesse de chute pour crystaux de glace
      REAL zzz
#include "YOETHF.h"
#include "FCTTRE.h"
      fallv (zzz) = 3.29/2.0 * ((zzz)**0.16)
ccc      fallv (zzz) = 3.29/3.0 * ((zzz)**0.16)
ccc      fallv (zzz) = 3.29 * ((zzz)**0.16)
c
      DATA appel1er /.TRUE./
c
      IF (appel1er) THEN
c
         PRINT*, 'fisrtilp, calcrat:', calcrat
         PRINT*, 'fisrtilp, ninter:', ninter
         PRINT*, 'fisrtilp, evap_prec:', evap_prec
         PRINT*, 'fisrtilp, cpartiel:', cpartiel
         IF (ABS(dtime/REAL(ninter)-360.0).GT.0.001) THEN
          PRINT*, 'fisrtilp: Ce n est pas prevu, voir Z.X.Li', dtime
          PRINT*, 'Je prefere un sous-intervalle de 6 minutes'
          CALL abort
         ENDIF
         appel1er = .FALSE.
c
cAA initialiation provisoire
       a_tr_sca(1) = -0.5
       a_tr_sca(2) = -0.5
       a_tr_sca(3) = -0.5
       a_tr_sca(4) = -0.5
c
cAA Initialisation a 1 des coefs des fractions lessivees 
c
      DO k = 1, klev
       DO i = 1, klon
          pfrac_nucl(i,k)=1.
          pfrac_1nucl(i,k)=1.
          pfrac_impa(i,k)=1.
       ENDDO 
      ENDDO 

      ENDIF          !  test sur appel1er
c
cMAf Initialisation a 0 de zoliq
       DO i = 1, klon
          zoliq(i)=0.
       ENDDO 
c Determiner les nuages froids par leur temperature
c
      ztglace = RTT - 15.0
      nexpo = 6
ccc      nexpo = 1
c
c Initialiser les sorties:
c
      DO k = 1, klev+1
      DO i = 1, klon
         prfl(i,k) = 0.0
         psfl(i,k) = 0.0
      ENDDO
      ENDDO

      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
         d_ql(i,k) = 0.0
         rneb(i,k) = 0.0
         radliq(i,k) = 0.0
         frac_nucl(i,k) = 1. 
         frac_impa(i,k) = 1. 
      ENDDO
      ENDDO
      DO i = 1, klon
         rain(i) = 0.0
         snow(i) = 0.0
      ENDDO
c
c Initialiser le flux de precipitation a zero
c
      DO i = 1, klon
         zrfl(i) = 0.0
         zneb(i) = seuil_neb
      ENDDO
c
c
cAA Pour plus de securite 

      zalpha_tr   = 0.
      zfrac_lessi = 0.

cAA----------------------------------------------------------
c
c Boucle verticale (du haut vers le bas)
c
      DO 9999 k = klev, 1, -1
c
cAA----------------------------------------------------------
c
      DO i = 1, klon
         zt(i)=t(i,k)
         zq(i)=q(i,k)
      ENDDO
c
c Calculer l'evaporation de la precipitation
c
      IF (evap_prec) THEN
      DO i = 1, klon
      IF (zrfl(i) .GT.0.) THEN
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-zt(i)))
           zqs(i)= R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
           zqs(i)=MIN(0.5,zqs(i))
           zcor=1./(1.-RETV*zqs(i))
           zqs(i)=zqs(i)*zcor
         ELSE
           IF (zt(i) .LT. t_coup) THEN
              zqs(i) = qsats(zt(i)) / pplay(i,k)
           ELSE
              zqs(i) = qsatl(zt(i)) / pplay(i,k)
           ENDIF
         ENDIF
         zqev = MAX (0.0, (zqs(i)-zq(i))*zneb(i) )
         zqevt = coef_eva * (1.0-zq(i)/zqs(i)) * SQRT(zrfl(i))
     .         * (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
         zqevt = MAX(0.0,MIN(zqevt,zrfl(i)))
     .         * RG*dtime/(paprs(i,k)-paprs(i,k+1))
         zqev = MIN (zqev, zqevt)
         zrfln(i) = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1))
     .                            /RG/dtime
         zq(i) = zq(i) - (zrfln(i)-zrfl(i))
     .             * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
         zt(i) = zt(i) + (zrfln(i)-zrfl(i))
     .             * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
     .             * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
         zrfl(i) = zrfln(i)
      ENDIF
      ENDDO
      ENDIF
c
c Calculer Qs et L/Cp*dQs/dT:
c
      IF (thermcep) THEN
         DO i = 1, klon
           zdelta = MAX(0.,SIGN(1.,RTT-zt(i)))
           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
           zcvm5 = zcvm5 /RCPD/(1.0+RVTMP2*zq(i))
           zqs(i) = R2ES*FOEEW(zt(i),zdelta)/pplay(i,k)
           zqs(i) = MIN(0.5,zqs(i))
           zcor = 1./(1.-RETV*zqs(i))
           zqs(i) = zqs(i)*zcor
           zdqs(i) = FOEDE(zt(i),zdelta,zcvm5,zqs(i),zcor)
         ENDDO
      ELSE
         DO i = 1, klon
            IF (zt(i).LT.t_coup) THEN
               zqs(i) = qsats(zt(i))/pplay(i,k)
               zdqs(i) = dqsats(zt(i),zqs(i))
            ELSE
               zqs(i) = qsatl(zt(i))/pplay(i,k)
               zdqs(i) = dqsatl(zt(i),zqs(i))
            ENDIF
         ENDDO
      ENDIF
c
c Determiner la condensation partielle et calculer la quantite
c de l'eau condensee:
c
      IF (cpartiel) THEN
         DO i = 1, klon
c
            zdelq = ratqs(i,k) * zq(i)
            rneb(i,k) = (zq(i)+zdelq-zqs(i)) / (2.0*zdelq)
            zqn(i) = (zq(i)+zdelq+zqs(i))/2.0
            IF (rneb(i,k) .LE. 0.0) zqn(i) = 0.0
            IF (rneb(i,k) .GE. 1.0) zqn(i) = zq(i)
            rneb(i,k) = MAX(0.0,MIN(1.0,rneb(i,k)))
            zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)/(1.+zdqs(i))
            
c--Olivier
            RHcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
            IF (rneb(i,k) .LE. 0.0) RHcl(i,k)=zq(i)/zqs(i)
            IF (rneb(i,k) .GE. 1.0) RHcl(i,k)=1.0
c--fin
            
         ENDDO
      ELSE
         DO i = 1, klon
            IF (zq(i).GT.zqs(i)) THEN
               rneb(i,k) = 1.0
            ELSE
               rneb(i,k) = 0.0
            ENDIF
            zcond(i) = MAX(0.0,zq(i)-zqs(i))/(1.+zdqs(i))
         ENDDO
      ENDIF
c
      DO i = 1, klon
         zq(i) = zq(i) - zcond(i)
         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD
      ENDDO
c
c Partager l'eau condensee en precipitation et eau liquide nuageuse
c
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         zoliq(i) = zcond(i)
         zrho(i) = pplay(i,k) / zt(i) / RD
         zdz(i) = (paprs(i,k)-paprs(i,k+1)) / (zrho(i)*RG)
         zfice(i) = 1.0 - (zt(i)-ztglace) / (273.13-ztglace)
         zfice(i) = MIN(MAX(zfice(i),0.0),1.0)
         zfice(i) = zfice(i)**nexpo
         zneb(i) = MAX(rneb(i,k), seuil_neb)
         radliq(i,k) = zoliq(i)/REAL(ninter+1)
      ENDIF
      ENDDO
c
      DO n = 1, ninter
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         zchau(i) = ct*dtime/REAL(ninter) * zoliq(i)
     .          * (1.0-EXP(-(zoliq(i)/zneb(i)/cl)**2)) *(1.-zfice(i))
         zrhol(i) = zrho(i) * zoliq(i) / zneb(i)
         zfroi(i) = dtime/REAL(ninter)/zdz(i)*zoliq(i)
     .              *fallv(zrhol(i)) * zfice(i)
         ztot(i) = zchau(i) + zfroi(i)
         IF (zneb(i).EQ.seuil_neb) ztot(i) = 0.0
         ztot(i) = MIN(MAX(ztot(i),0.0),zoliq(i))
         zoliq(i) = MAX(zoliq(i)-ztot(i), 0.0)
         radliq(i,k) = radliq(i,k) + zoliq(i)/REAL(ninter+1)
      ENDIF
      ENDDO
      ENDDO
c
      DO i = 1, klon
      IF (rneb(i,k).GT.0.0) THEN
         d_ql(i,k) = zoliq(i)
         zrfl(i) = zrfl(i)+ MAX(zcond(i)-zoliq(i),0.0)
     .                    * (paprs(i,k)-paprs(i,k+1))/(RG*dtime)
      ENDIF
      IF (zt(i).LT.RTT) THEN
        psfl(i,k)=zrfl(i)
      ELSE
        prfl(i,k)=zrfl(i)
      ENDIF
      ENDDO
c
c Calculer les tendances de q et de t:
c
      DO i = 1, klon
         d_q(i,k) = zq(i) - q(i,k)
         d_t(i,k) = zt(i) - t(i,k)
      ENDDO
c
cAA--------------- Calcul du lessivage stratiforme  -------------

      DO i = 1,klon
c
         zprec_cond(i) = MAX(zcond(i)-zoliq(i),0.0)
     .                * (paprs(i,k)-paprs(i,k+1))/RG
         IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
cAA lessivage nucleation LMD5 dans la couche elle-meme
            if (t(i,k) .GE. ztglace) THEN
               zalpha_tr = a_tr_sca(3)
            else
               zalpha_tr = a_tr_sca(4)
            endif
            zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
            pfrac_nucl(i,k)=pfrac_nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
            frac_nucl(i,k)= 1.-zneb(i)*zfrac_lessi 
c
c nucleation avec un facteur -1 au lieu de -0.5
            zfrac_lessi = 1. - EXP(-zprec_cond(i)/zneb(i))
            pfrac_1nucl(i,k)=pfrac_1nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
         ENDIF
c
      ENDDO      ! boucle sur i
c
cAA Lessivage par impaction dans les couches en-dessous
      DO kk = k-1, 1, -1
        DO i = 1, klon
          IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
            if (t(i,kk) .GE. ztglace) THEN
              zalpha_tr = a_tr_sca(1)
            else
              zalpha_tr = a_tr_sca(2)
            endif
            zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
            pfrac_impa(i,kk)=pfrac_impa(i,kk)*(1.-zneb(i)*zfrac_lessi)
            frac_impa(i,kk)= 1.-zneb(i)*zfrac_lessi
          ENDIF
        ENDDO
      ENDDO
c
cAA----------------------------------------------------------
c                     FIN DE BOUCLE SUR K   
 9999 CONTINUE
c
cAA-----------------------------------------------------------
c
c Pluie ou neige au sol selon la temperature de la 1ere couche
c
      DO i = 1, klon
      IF ((t(i,1)+d_t(i,1)) .LT. RTT) THEN
         snow(i) = zrfl(i)
      ELSE
         rain(i) = zrfl(i)
      ENDIF
      ENDDO
c
      RETURN
      END
