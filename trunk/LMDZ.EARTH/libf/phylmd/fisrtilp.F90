!
! $Id: fisrtilp.F90 1472 2010-12-23 16:38:42Z lguez $
!
!
SUBROUTINE fisrtilp(dtime,paprs,pplay,t,q,ptconv,ratqs, &
     d_t, d_q, d_ql, rneb, radliq, rain, snow, &
     pfrac_impa, pfrac_nucl, pfrac_1nucl, &
     frac_impa, frac_nucl, &
     prfl, psfl, rhcl, zqta, fraca, &
     ztv, zpspsk, ztla, zthl, iflag_cldcon)

  !
  USE dimphy
  IMPLICIT none
  !======================================================================
  ! Auteur(s): Z.X. Li (LMD/CNRS)
  ! Date: le 20 mars 1995
  ! Objet: condensation et precipitation stratiforme.
  !        schema de nuage
  !======================================================================
  !======================================================================
  !ym include "dimensions.h"
  !ym include "dimphy.h"
  include "YOMCST.h"
  include "tracstoke.h"
  include "fisrtilp.h"
  !
  ! Arguments:
  !
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
  REAL rhcl(klon,klev) ! humidite relative en ciel clair
  REAL rain(klon) ! pluies (mm/s)
  REAL snow(klon) ! neige (mm/s)
  REAL prfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s)
  REAL psfl(klon,klev+1) ! flux d'eau precipitante aux interfaces (kg/m2/s) 
  REAL ztv(klon,klev)
  REAL zqta(klon,klev),fraca(klon,klev) 
  REAL sigma1(klon,klev),sigma2(klon,klev)
  REAL qltot(klon,klev),ctot(klon,klev)
  REAL zpspsk(klon,klev),ztla(klon,klev)
  REAL zthl(klon,klev)

  logical lognormale(klon)

  !AA
  ! Coeffients de fraction lessivee : pour OFF-LINE
  !
  REAL pfrac_nucl(klon,klev)
  REAL pfrac_1nucl(klon,klev)
  REAL pfrac_impa(klon,klev)
  !
  ! Fraction d'aerosols lessivee par impaction et par nucleation
  ! POur ON-LINE
  !
  REAL frac_impa(klon,klev)
  REAL frac_nucl(klon,klev)
  real zct      ,zcl
  !AA
  !
  ! Options du programme:
  !
  REAL seuil_neb ! un nuage existe vraiment au-dela
  PARAMETER (seuil_neb=0.001)

  INTEGER ninter ! sous-intervals pour la precipitation
  INTEGER ncoreczq  
  INTEGER iflag_cldcon
  PARAMETER (ninter=5)
  LOGICAL evap_prec ! evaporation de la pluie
  PARAMETER (evap_prec=.TRUE.)
  REAL ratqs(klon,klev) ! determine la largeur de distribution de vapeur
  logical ptconv(klon,klev) ! determine la largeur de distribution de vapeur

  real zpdf_sig(klon),zpdf_k(klon),zpdf_delta(klon)
  real Zpdf_a(klon),zpdf_b(klon),zpdf_e1(klon),zpdf_e2(klon)
  real erf   
  REAL qcloud(klon)
  !
  LOGICAL cpartiel ! condensation partielle
  PARAMETER (cpartiel=.TRUE.)
  REAL t_coup
  PARAMETER (t_coup=234.0)
  !
  ! Variables locales:
  !
  INTEGER i, k, n, kk
  REAL zqs(klon), zdqs(klon), zdelta, zcor, zcvm5   
  REAL zrfl(klon), zrfln(klon), zqev, zqevt
  REAL zoliq(klon), zcond(klon), zq(klon), zqn(klon), zdelq
  REAL ztglace, zt(klon)
  INTEGER nexpo ! exponentiel pour glace/eau
  REAL zdz(klon),zrho(klon),ztot      , zrhol(klon)
  REAL zchau      ,zfroi      ,zfice(klon),zneb(klon)
  !
  LOGICAL appel1er
  SAVE appel1er
  !$OMP THREADPRIVATE(appel1er)
  !
  !---------------------------------------------------------------
  !
  !AA Variables traceurs:
  !AA  Provisoire !!! Parametres alpha du lessivage
  !AA  A priori on a 4 scavenging # possibles
  !
  REAL a_tr_sca(4)
  save a_tr_sca
  !$OMP THREADPRIVATE(a_tr_sca)
  !
  ! Variables intermediaires
  !
  REAL zalpha_tr
  REAL zfrac_lessi
  REAL zprec_cond(klon)
  !AA
  REAL zmair, zcpair, zcpeau
  !     Pour la conversion eau-neige
  REAL zlh_solid(klon), zm_solid
  !IM 
  !ym      INTEGER klevm1
  !---------------------------------------------------------------
  !
  ! Fonctions en ligne:
  !
  REAL fallvs,fallvc ! vitesse de chute pour crystaux de glace
  REAL zzz
  include "YOETHF.h"
  include "FCTTRE.h"
  fallvc (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_con
  fallvs (zzz) = 3.29/2.0 * ((zzz)**0.16) * ffallv_lsc
  !
  DATA appel1er /.TRUE./
  !ym
  zdelq=0.0

  print*,'NUAGES4 A. JAM'
  IF (appel1er) THEN
     !
     PRINT*, 'fisrtilp, ninter:', ninter
     PRINT*, 'fisrtilp, evap_prec:', evap_prec
     PRINT*, 'fisrtilp, cpartiel:', cpartiel
     IF (ABS(dtime/REAL(ninter)-360.0).GT.0.001) THEN
        PRINT*, 'fisrtilp: Ce n est pas prevu, voir Z.X.Li', dtime
        PRINT*, 'Je prefere un sous-intervalle de 6 minutes'
        !         CALL abort
     ENDIF
     appel1er = .FALSE.
     !
     !AA initialiation provisoire
     a_tr_sca(1) = -0.5
     a_tr_sca(2) = -0.5
     a_tr_sca(3) = -0.5
     a_tr_sca(4) = -0.5
     !
     !AA Initialisation a 1 des coefs des fractions lessivees 
     !
     !cdir collapse
     DO k = 1, klev
        DO i = 1, klon
           pfrac_nucl(i,k)=1.
           pfrac_1nucl(i,k)=1.
           pfrac_impa(i,k)=1.
        ENDDO
     ENDDO

  ENDIF          !  test sur appel1er
  !
  !MAf Initialisation a 0 de zoliq
  !      DO i = 1, klon
  !         zoliq(i)=0.
  !      ENDDO 
  ! Determiner les nuages froids par leur temperature
  !  nexpo regle la raideur de la transition eau liquide / eau glace.
  !
  ztglace = RTT - 15.0
  nexpo = 6
  !cc      nexpo = 1
  !
  ! Initialiser les sorties:
  !
  !cdir collapse
  DO k = 1, klev+1
     DO i = 1, klon
        prfl(i,k) = 0.0
        psfl(i,k) = 0.0
     ENDDO
  ENDDO

  !cdir collapse
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
     zoliq(i)=0.
     !     ENDDO
     !
     ! Initialiser le flux de precipitation a zero
     !
     !     DO i = 1, klon
     zrfl(i) = 0.0
     zneb(i) = seuil_neb
  ENDDO
  !
  !
  !AA Pour plus de securite 

  zalpha_tr   = 0.
  zfrac_lessi = 0.

  !AA----------------------------------------------------------
  !
  ncoreczq=0
  ! Boucle verticale (du haut vers le bas)
  !
  !IM : klevm1
  !ym      klevm1=klev-1
  DO k = klev, 1, -1
     !
     !AA----------------------------------------------------------
     !
     DO i = 1, klon
        zt(i)=t(i,k)
        zq(i)=q(i,k)
     ENDDO
     !
     ! Calculer la varition de temp. de l'air du a la chaleur sensible
     ! transporter par la pluie.
     ! Il resterait a rajouter cet effet de la chaleur sensible sur les
     ! flux de surface, du a la diff. de temp. entre le 1er niveau et la
     ! surface.
     !
     IF(k.LE.klevm1) THEN         
        DO i = 1, klon
           !IM
           zmair=(paprs(i,k)-paprs(i,k+1))/RG
           zcpair=RCPD*(1.0+RVTMP2*zq(i))
           zcpeau=RCPD*RVTMP2
           zt(i) = ( (t(i,k+1)+d_t(i,k+1))*zrfl(i)*dtime*zcpeau &
                + zmair*zcpair*zt(i) ) &
                / (zmair*zcpair + zrfl(i)*dtime*zcpeau)
           !     C        WRITE (6,*) 'cppluie ', zt(i)-(t(i,k+1)+d_t(i,k+1))
        ENDDO
     ENDIF
     !
     !
     ! Calculer l'evaporation de la precipitation
     !


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
              zqevt = coef_eva * (1.0-zq(i)/zqs(i)) * SQRT(zrfl(i)) &
                   * (paprs(i,k)-paprs(i,k+1))/pplay(i,k)*zt(i)*RD/RG
              zqevt = MAX(0.0,MIN(zqevt,zrfl(i))) &
                   * RG*dtime/(paprs(i,k)-paprs(i,k+1))
              zqev = MIN (zqev, zqevt)
              zrfln(i) = zrfl(i) - zqev*(paprs(i,k)-paprs(i,k+1)) &
                   /RG/dtime

              ! pour la glace, on ré-évapore toute la précip dans la
              ! couche du dessous
              ! la glace venant de la couche du dessus est simplement
              ! dans la couche du dessous.

              IF (zt(i) .LT. t_coup.and.reevap_ice) zrfln(i)=0.

              zq(i) = zq(i) - (zrfln(i)-zrfl(i)) &
                   * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime
              zt(i) = zt(i) + (zrfln(i)-zrfl(i)) &
                   * (RG/(paprs(i,k)-paprs(i,k+1)))*dtime &
                   * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
              zrfl(i) = zrfln(i)
           ENDIF
        ENDDO
     ENDIF
     !
     ! Calculer Qs et L/Cp*dQs/dT:
     !
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
     !
     ! Determiner la condensation partielle et calculer la quantite
     ! de l'eau condensee:
     !

     IF (cpartiel) THEN

        !        print*,'Dans partiel k=',k
        !
        !   Calcul de l'eau condensee et de la fraction nuageuse et de l'eau
        !   nuageuse a partir des PDF de Sandrine Bony.
        !   rneb  : fraction nuageuse
        !   zqn   : eau totale dans le nuage
        !   zcond : eau condensee moyenne dans la maille.
        !  on prend en compte le réchauffement qui diminue la partie
        ! condensee
        !
        !   Version avec les raqts

        if (iflag_pdf.eq.0) then

           do i=1,klon
              zdelq = min(ratqs(i,k),0.99) * zq(i)
              rneb(i,k) = (zq(i)+zdelq-zqs(i)) / (2.0*zdelq)
              zqn(i) = (zq(i)+zdelq+zqs(i))/2.0
           enddo

        else
           !
           !   Version avec les nouvelles PDFs.
           do i=1,klon
              if(zq(i).lt.1.e-15) then
                 ncoreczq=ncoreczq+1
                 zq(i)=1.e-15
              endif
           enddo

           if (iflag_cldcon>=5) then

              call cloudth(klon,klev,k,ztv, &
                   zq,zqta,fraca, &
                   qcloud,ctot,zpspsk,paprs,ztla,zthl, &
                   ratqs,zqs,t)

              do i=1,klon
                 rneb(i,k)=ctot(i,k)
                 zqn(i)=qcloud(i)
              enddo

           endif

           if (iflag_cldcon <= 4) then
              lognormale = .true.
           elseif (iflag_cldcon == 6) then
              ! lognormale en l'absence des thermiques
              lognormale = fraca(:,k) < 1e-10
           else
              ! Dans le cas iflag_cldcon=5, on prend systématiquement la
              ! bi-gaussienne
              lognormale = .false.
           end if

           do i=1,klon
              if (lognormale(i)) then
                 zpdf_sig(i)=ratqs(i,k)*zq(i)
                 zpdf_k(i)=-sqrt(log(1.+(zpdf_sig(i)/zq(i))**2))
                 zpdf_delta(i)=log(zq(i)/zqs(i))
                 zpdf_a(i)=zpdf_delta(i)/(zpdf_k(i)*sqrt(2.))
                 zpdf_b(i)=zpdf_k(i)/(2.*sqrt(2.))
                 zpdf_e1(i)=zpdf_a(i)-zpdf_b(i)
                 zpdf_e1(i)=sign(min(abs(zpdf_e1(i)),5.),zpdf_e1(i))
                 zpdf_e1(i)=1.-erf(zpdf_e1(i))
                 zpdf_e2(i)=zpdf_a(i)+zpdf_b(i)
                 zpdf_e2(i)=sign(min(abs(zpdf_e2(i)),5.),zpdf_e2(i))
                 zpdf_e2(i)=1.-erf(zpdf_e2(i))
              endif
           enddo

           do i=1,klon
              if (lognormale(i)) then
                 if (zpdf_e1(i).lt.1.e-10) then
                    rneb(i,k)=0.
                    zqn(i)=zqs(i)
                 else
                    rneb(i,k)=0.5*zpdf_e1(i)
                    zqn(i)=zq(i)*zpdf_e2(i)/zpdf_e1(i)
                 endif
              endif

           enddo


        endif ! iflag_pdf

        DO i=1,klon
           IF (rneb(i,k) .LE. 0.0) THEN
              zqn(i) = 0.0
              rneb(i,k) = 0.0
              zcond(i) = 0.0
              rhcl(i,k)=zq(i)/zqs(i)
           ELSE IF (rneb(i,k) .GE. 1.0) THEN
              zqn(i) = zq(i)
              rneb(i,k) = 1.0                  
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))
              rhcl(i,k)=1.0
           ELSE
              zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)
              rhcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
           ENDIF
        ENDDO
        !         do i=1,klon
        !            IF (rneb(i,k) .LE. 0.0) zqn(i) = 0.0
        !            IF (rneb(i,k) .GE. 1.0) zqn(i) = zq(i)
        !            rneb(i,k) = MAX(0.0,MIN(1.0,rneb(i,k)))
        !c           zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)/(1.+zdqs(i))
        !c  On ne divise pas par 1+zdqs pour forcer a avoir l'eau predite par
        !c  la convection.
        !c  ATTENTION !!! Il va falloir verifier tout ca.
        !            zcond(i) = MAX(0.0,zqn(i)-zqs(i))*rneb(i,k)
        !c           print*,'ZDQS ',zdqs(i)
        !c--Olivier
        !            rhcl(i,k)=(zqs(i)+zq(i)-zdelq)/2./zqs(i)
        !            IF (rneb(i,k) .LE. 0.0) rhcl(i,k)=zq(i)/zqs(i)
        !            IF (rneb(i,k) .GE. 1.0) rhcl(i,k)=1.0
        !c--fin
        !           ENDDO
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
     !
     DO i = 1, klon
        zq(i) = zq(i) - zcond(i)
        !         zt(i) = zt(i) + zcond(i) * RLVTT/RCPD
        zt(i) = zt(i) + zcond(i) * RLVTT/RCPD/(1.0+RVTMP2*zq(i))
     ENDDO
     !
     ! Partager l'eau condensee en precipitation et eau liquide nuageuse
     !
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
     !
     DO n = 1, ninter
        DO i = 1, klon
           IF (rneb(i,k).GT.0.0) THEN
              zrhol(i) = zrho(i) * zoliq(i) / zneb(i)

              IF (zneb(i).EQ.seuil_neb) THEN
                 ztot = 0.0
              ELSE
                 !  quantite d'eau a eliminer: zchau
                 !  meme chose pour la glace: zfroi
                 if (ptconv(i,k)) then
                    zcl   =cld_lc_con
                    zct   =1./cld_tau_con
                    zfroi    = dtime/REAL(ninter)/zdz(i)*zoliq(i) &
                         *fallvc(zrhol(i)) * zfice(i)
                 else
                    zcl   =cld_lc_lsc
                    zct   =1./cld_tau_lsc
                    zfroi    = dtime/REAL(ninter)/zdz(i)*zoliq(i) &
                         *fallvs(zrhol(i)) * zfice(i)
                 endif
                 zchau    = zct   *dtime/REAL(ninter) * zoliq(i) &
                      *(1.0-EXP(-(zoliq(i)/zneb(i)/zcl   )**2)) *(1.-zfice(i))
                 ztot    = zchau    + zfroi
                 ztot    = MAX(ztot   ,0.0)
              ENDIF
              ztot    = MIN(ztot,zoliq(i))
              zoliq(i) = MAX(zoliq(i)-ztot   , 0.0)
              radliq(i,k) = radliq(i,k) + zoliq(i)/REAL(ninter+1)
           ENDIF
        ENDDO
     ENDDO
     !
     DO i = 1, klon
        IF (rneb(i,k).GT.0.0) THEN
           d_ql(i,k) = zoliq(i)
           zrfl(i) = zrfl(i)+ MAX(zcond(i)-zoliq(i),0.0) &
                * (paprs(i,k)-paprs(i,k+1))/(RG*dtime)
        ENDIF
        IF (zt(i).LT.RTT) THEN
           psfl(i,k)=zrfl(i)
        ELSE
           prfl(i,k)=zrfl(i)
        ENDIF
     ENDDO
     !
     ! Calculer les tendances de q et de t:
     !
     DO i = 1, klon
        d_q(i,k) = zq(i) - q(i,k)
        d_t(i,k) = zt(i) - t(i,k)
     ENDDO
     !
     !AA--------------- Calcul du lessivage stratiforme  -------------

     DO i = 1,klon
        !
        zprec_cond(i) = MAX(zcond(i)-zoliq(i),0.0) &
             * (paprs(i,k)-paprs(i,k+1))/RG
        IF (rneb(i,k).GT.0.0.and.zprec_cond(i).gt.0.) THEN
           !AA lessivage nucleation LMD5 dans la couche elle-meme
           if (t(i,k) .GE. ztglace) THEN
              zalpha_tr = a_tr_sca(3)
           else
              zalpha_tr = a_tr_sca(4)
           endif
           zfrac_lessi = 1. - EXP(zalpha_tr*zprec_cond(i)/zneb(i))
           pfrac_nucl(i,k)=pfrac_nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
           frac_nucl(i,k)= 1.-zneb(i)*zfrac_lessi 
           !
           ! nucleation avec un facteur -1 au lieu de -0.5
           zfrac_lessi = 1. - EXP(-zprec_cond(i)/zneb(i))
           pfrac_1nucl(i,k)=pfrac_1nucl(i,k)*(1.-zneb(i)*zfrac_lessi)
        ENDIF
        !
     ENDDO      ! boucle sur i
     !
     !AA Lessivage par impaction dans les couches en-dessous
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
     !
     !AA----------------------------------------------------------
     !                     FIN DE BOUCLE SUR K   
  end DO
  !
  !AA-----------------------------------------------------------
  !
  ! Pluie ou neige au sol selon la temperature de la 1ere couche
  !
  DO i = 1, klon
     IF ((t(i,1)+d_t(i,1)) .LT. RTT) THEN
        snow(i) = zrfl(i)
        zlh_solid(i) = RLSTT-RLVTT
     ELSE
        rain(i) = zrfl(i)
        zlh_solid(i) = 0.
     ENDIF
  ENDDO
  !
  ! For energy conservation : when snow is present, the solification
  ! latent heat is considered.
  DO k = 1, klev
     DO i = 1, klon
        zcpair=RCPD*(1.0+RVTMP2*(q(i,k)+d_q(i,k)))
        zmair=(paprs(i,k)-paprs(i,k+1))/RG
        zm_solid = (prfl(i,k)-prfl(i,k+1)+psfl(i,k)-psfl(i,k+1))*dtime
        d_t(i,k) = d_t(i,k) + zlh_solid(i) *zm_solid / (zcpair*zmair)
     END DO
  END DO
  !

  if (ncoreczq>0) then
     print*,'WARNING : ZQ dans fisrtilp ',ncoreczq,' val < 1.e-15.'
  endif

END SUBROUTINE fisrtilp
