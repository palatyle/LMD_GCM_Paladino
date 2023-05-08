!
MODULE coef_diff_turb_mod
!
! This module contains some procedures for calculation of the coefficients of the
! turbulent diffusion in the atmosphere and coefficients for turbulent diffusion 
! at surface(cdrag)
!
  IMPLICIT NONE
  
CONTAINS
!
!****************************************************************************************
!
  SUBROUTINE coef_diff_turb(dtime, nsrf, knon, ni, &
       ypaprs, ypplay, yu, yv, yq, yt, yts, yrugos, yqsurf, ycdragm, &
       ycoefm, ycoefh ,yq2)
 
    USE dimphy
!
! Calculate coefficients(ycoefm, ycoefh) for turbulent diffusion in the 
! atmosphere 
! NB! No values are calculated between surface and the first model layer. 
!     ycoefm(:,1) and ycoefh(:,1) are not valid !!!
!
!
! Input arguments
!****************************************************************************************
    REAL, INTENT(IN)                           :: dtime
    INTEGER, INTENT(IN)                        :: nsrf, knon
    INTEGER, DIMENSION(klon), INTENT(IN)       :: ni
    REAL, DIMENSION(klon,klev+1), INTENT(IN)   :: ypaprs
    REAL, DIMENSION(klon,klev), INTENT(IN)     :: ypplay
    REAL, DIMENSION(klon,klev), INTENT(IN)     :: yu, yv
    REAL, DIMENSION(klon,klev), INTENT(IN)     :: yq, yt
    REAL, DIMENSION(klon), INTENT(IN)          :: yts, yrugos, yqsurf
    REAL, DIMENSION(klon), INTENT(IN)          :: ycdragm

! InOutput arguments
!****************************************************************************************
    REAL, DIMENSION(klon,klev+1), INTENT(INOUT):: yq2
  
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon,klev), INTENT(OUT)    :: ycoefh
    REAL, DIMENSION(klon,klev), INTENT(OUT)    :: ycoefm

! Other local variables
!****************************************************************************************
    INTEGER                                    :: k, i, j
    REAL, DIMENSION(klon,klev)                 :: ycoefm0, ycoefh0, yzlay, yteta
    REAL, DIMENSION(klon,klev+1)               :: yzlev, q2diag, ykmm, ykmn, ykmq
    REAL, DIMENSION(klon)                      :: yustar

! Include
!****************************************************************************************
    INCLUDE "clesphys.h"
    INCLUDE "indicesol.h"
    INCLUDE "iniprint.h"
    INCLUDE "compbl.h"
    INCLUDE "YOETHF.h"
    INCLUDE "YOMCST.h"


!****************************************************************************************    
! Calcul de coefficients de diffusion turbulent de l'atmosphere : 
! ycoefm(:,2:klev), ycoefh(:,2:klev) 
!
!****************************************************************************************    

    CALL coefkz(nsrf, knon, ypaprs, ypplay, &
         ksta, ksta_ter, &
         yts, yrugos, yu, yv, yt, yq, &
         yqsurf, &
         ycoefm, ycoefh)
  
!****************************************************************************************
! Eventuelle recalcule des coeffeicients de diffusion turbulent de l'atmosphere : 
! ycoefm(:,2:klev), ycoefh(:,2:klev) 
!
!****************************************************************************************

    IF (iflag_pbl.EQ.1) THEN
       CALL coefkz2(nsrf, knon, ypaprs, ypplay, yt, &
            ycoefm0, ycoefh0)

       DO k = 2, klev
          DO i = 1, knon
             ycoefm(i,k) = MAX(ycoefm(i,k),ycoefm0(i,k))
             ycoefh(i,k) = MAX(ycoefh(i,k),ycoefh0(i,k))
          ENDDO
       ENDDO
    ENDIF

  
!****************************************************************************************  
! Calcul d'une diffusion minimale pour les conditions tres stables
!
!****************************************************************************************
    IF (ok_kzmin) THEN
       CALL coefkzmin(knon,ypaprs,ypplay,yu,yv,yt,yq,ycdragm, &
            ycoefm0,ycoefh0)
       
       DO k = 2, klev
          DO i = 1, knon
             ycoefm(i,k) = MAX(ycoefm(i,k),ycoefm0(i,k))
             ycoefh(i,k) = MAX(ycoefh(i,k),ycoefh0(i,k))
          ENDDO
       ENDDO
       
    ENDIF

  
!****************************************************************************************
! MELLOR ET YAMADA adapte a Mars Richard Fournier et Frederic Hourdin
! 
!****************************************************************************************

    IF (iflag_pbl.GE.3) THEN

       yzlay(1:knon,1)= &
            RD*yt(1:knon,1)/(0.5*(ypaprs(1:knon,1)+ypplay(1:knon,1))) &
            *(ypaprs(1:knon,1)-ypplay(1:knon,1))/RG
       DO k=2,klev
          DO i = 1, knon
             yzlay(i,k)= &
                  yzlay(i,k-1)+RD*0.5*(yt(i,k-1)+yt(i,k)) &
                  /ypaprs(i,k)*(ypplay(i,k-1)-ypplay(i,k))/RG
          END DO
       END DO

       DO k=1,klev
          DO i = 1, knon
             yteta(i,k)= &
                  yt(i,k)*(ypaprs(i,1)/ypplay(i,k))**RKAPPA &
                  *(1.+0.61*yq(i,k))
          END DO
       END DO

       yzlev(1:knon,1)=0.
       yzlev(1:knon,klev+1)=2.*yzlay(1:knon,klev)-yzlay(1:knon,klev-1)
       DO k=2,klev
          DO i = 1, knon
             yzlev(i,k)=0.5*(yzlay(i,k)+yzlay(i,k-1))
          END DO
       END DO

!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! Pour memoire, le papier Hourdin et al. 2002 a ete obtenur avec un
!!$! bug sur les coefficients de surface :
!!$!          ycdragh(1:knon) = ycoefm(1:knon,1)
!!$!          ycdragm(1:knon) = ycoefh(1:knon,1)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       CALL ustarhb(knon,yu,yv,ycdragm, yustar)
     
       IF (prt_level > 9) THEN
          WRITE(lunout,*) 'USTAR = ',yustar
       ENDIF
         
!   iflag_pbl peut etre utilise comme longuer de melange
       IF (iflag_pbl.GE.11) THEN
          CALL vdif_kcay(knon,dtime,RG,RD,ypaprs,yt, &
               yzlev,yzlay,yu,yv,yteta, &
               ycdragm,yq2,q2diag,ykmm,ykmn,yustar, &
               iflag_pbl)
       ELSE
          CALL yamada4(knon,dtime,RG,RD,ypaprs,yt, &
               yzlev,yzlay,yu,yv,yteta, &
               ycdragm,yq2,ykmm,ykmn,ykmq,yustar, &
               iflag_pbl)
       ENDIF
       
       ycoefm(1:knon,2:klev)=ykmm(1:knon,2:klev)
       ycoefh(1:knon,2:klev)=ykmn(1:knon,2:klev)
                
    ENDIF !(iflag_pbl.ge.3)

  END SUBROUTINE coef_diff_turb
!
!****************************************************************************************
!
  SUBROUTINE coefkz(nsrf, knon, paprs, pplay, &
       ksta, ksta_ter, &
       ts, rugos, &
       u,v,t,q, &
       qsurf, &
       pcfm, pcfh)
    
    USE dimphy
  
!======================================================================
! Auteur(s) F. Hourdin, M. Forichon, Z.X. Li (LMD/CNRS) date: 19930922
!           (une version strictement identique a l'ancien modele)
! Objet: calculer le coefficient du frottement du sol (Cdrag) et les
!        coefficients d'echange turbulent dans l'atmosphere.
! Arguments:
! nsrf-----input-I- indicateur de la nature du sol
! knon-----input-I- nombre de points a traiter
! paprs----input-R- pregssion a chaque intercouche (en Pa)
! pplay----input-R- pression au milieu de chaque couche (en Pa)
! ts-------input-R- temperature du sol (en Kelvin)
! rugos----input-R- longeur de rugosite (en m)
! u--------input-R- vitesse u
! v--------input-R- vitesse v
! t--------input-R- temperature (K)
! q--------input-R- vapeur d'eau (kg/kg)
!
! pcfm-----output-R- coefficients a calculer (vitesse)
! pcfh-----output-R- coefficients a calculer (chaleur et humidite)
!======================================================================
    INCLUDE "YOETHF.h"
    INCLUDE "FCTTRE.h"
    INCLUDE "iniprint.h"
    INCLUDE "indicesol.h"
    INCLUDE "compbl.h"
    INCLUDE "YOMCST.h"
!
! Arguments:
!
    INTEGER, INTENT(IN)                      :: knon, nsrf
    REAL, INTENT(IN)                         :: ksta, ksta_ter
    REAL, DIMENSION(klon), INTENT(IN)        :: ts
    REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay
    REAL, DIMENSION(klon,klev), INTENT(IN)   :: u, v, t, q
    REAL, DIMENSION(klon), INTENT(IN)        :: rugos
    REAL, DIMENSION(klon), INTENT(IN)        :: qsurf

    REAL, DIMENSION(klon,klev), INTENT(OUT)  :: pcfm, pcfh

!
! Local variables:
!
    INTEGER, DIMENSION(klon)    :: itop ! numero de couche du sommet de la couche limite
!
! Quelques constantes et options:
!
    REAL, PARAMETER :: cepdu2=0.1**2
    REAL, PARAMETER :: CKAP=0.4
    REAL, PARAMETER :: cb=5.0
    REAL, PARAMETER :: cc=5.0
    REAL, PARAMETER :: cd=5.0
    REAL, PARAMETER :: clam=160.0
    REAL, PARAMETER :: ratqs=0.05 ! largeur de distribution de vapeur d'eau
    LOGICAL, PARAMETER :: richum=.TRUE. ! utilise le nombre de Richardson humide
    REAL, PARAMETER :: ric=0.4 ! nombre de Richardson critique
    REAL, PARAMETER :: prandtl=0.4
    REAL kstable ! diffusion minimale (situation stable)
    ! GKtest
    ! PARAMETER (kstable=1.0e-10)
!IM: 261103     REAL kstable_ter, kstable_sinon
!IM: 211003 cf GK   PARAMETER (kstable_ter = 1.0e-6)
!IM: 261103     PARAMETER (kstable_ter = 1.0e-8)
!IM: 261103   PARAMETER (kstable_ter = 1.0e-10)
!IM: 261103   PARAMETER (kstable_sinon = 1.0e-10)
    ! fin GKtest
    REAL, PARAMETER :: mixlen=35.0 ! constante controlant longueur de melange
    INTEGER isommet ! le sommet de la couche limite
    LOGICAL, PARAMETER :: tvirtu=.TRUE. ! calculer Ri d'une maniere plus performante
    LOGICAL, PARAMETER :: opt_ec=.FALSE.! formule du Centre Europeen dans l'atmosphere

!
! Variables locales:
    INTEGER i, k !IM 120704
    REAL zgeop(klon,klev)
    REAL zmgeom(klon)
    REAL zri(klon)
    REAL zl2(klon)
    REAL zdphi, zdu2, ztvd, ztvu, zcdn
    REAL zscf
    REAL zt, zq, zdelta, zcvm5, zcor, zqs, zfr, zdqs
    REAL z2geomf, zalh2, zalm2, zscfh, zscfm
    REAL, PARAMETER :: t_coup=273.15
    LOGICAL, PARAMETER :: check=.FALSE.
!
! contre-gradient pour la chaleur sensible: Kelvin/metre
    REAL gamt(2:klev)

    LOGICAL, SAVE :: appel1er=.TRUE.
    !$OMP THREADPRIVATE(appel1er)
!
! Fonctions thermodynamiques et fonctions d'instabilite
    REAL fsta, fins, x

    fsta(x) = 1.0 / (1.0+10.0*x*(1+8.0*x))
    fins(x) = SQRT(1.0-18.0*x)

    isommet=klev
      
    IF (appel1er) THEN
       IF (prt_level > 9) THEN
          WRITE(lunout,*)'coefkz, opt_ec:', opt_ec
          WRITE(lunout,*)'coefkz, richum:', richum
          IF (richum) WRITE(lunout,*)'coefkz, ratqs:', ratqs
          WRITE(lunout,*)'coefkz, isommet:', isommet
          WRITE(lunout,*)'coefkz, tvirtu:', tvirtu
          appel1er = .FALSE.
       ENDIF
    ENDIF
!
! Initialiser les sorties
!
    DO k = 1, klev
       DO i = 1, knon
          pcfm(i,k) = 0.0
          pcfh(i,k) = 0.0
       ENDDO
    ENDDO
    DO i = 1, knon
       itop(i) = 0
    ENDDO
    
!
! Prescrire la valeur de contre-gradient
!
    IF (iflag_pbl.EQ.1) THEN
       DO k = 3, klev
          gamt(k) = -1.0E-03
       ENDDO
       gamt(2) = -2.5E-03
    ELSE
       DO k = 2, klev
          gamt(k) = 0.0
       ENDDO
    ENDIF
!IM cf JLD/ GKtest
    IF ( nsrf .NE. is_oce ) THEN
!IM 261103     kstable = kstable_ter
       kstable = ksta_ter
    ELSE
!IM 261103     kstable = kstable_sinon
       kstable = ksta
    ENDIF
!IM cf JLD/ GKtest fin

!
! Calculer les geopotentiels de chaque couche
!
    DO i = 1, knon
       zgeop(i,1) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1))) &
            * (paprs(i,1)-pplay(i,1))
    ENDDO
    DO k = 2, klev
       DO i = 1, knon
          zgeop(i,k) = zgeop(i,k-1) &
               + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k) &
               * (pplay(i,k-1)-pplay(i,k))
       ENDDO
    ENDDO

!
! Calculer les coefficients turbulents dans l'atmosphere
!
    DO i = 1, knon
       itop(i) = isommet
    ENDDO


    DO k = 2, isommet
       DO i = 1, knon
          zdu2=MAX(cepdu2,(u(i,k)-u(i,k-1))**2 &
               +(v(i,k)-v(i,k-1))**2)
          zmgeom(i)=zgeop(i,k)-zgeop(i,k-1)
          zdphi =zmgeom(i) / 2.0
          zt = (t(i,k)+t(i,k-1)) * 0.5
          zq = (q(i,k)+q(i,k-1)) * 0.5

!
! Calculer Qs et dQs/dT:
!
          IF (thermcep) THEN
             zdelta = MAX(0.,SIGN(1.,RTT-zt))
             zcvm5 = R5LES*RLVTT/RCPD/(1.0+RVTMP2*zq)*(1.-zdelta) &
                  + R5IES*RLSTT/RCPD/(1.0+RVTMP2*zq)*zdelta
             zqs = R2ES * FOEEW(zt,zdelta) / pplay(i,k)
             zqs = MIN(0.5,zqs)
             zcor = 1./(1.-RETV*zqs)
             zqs = zqs*zcor
             zdqs = FOEDE(zt,zdelta,zcvm5,zqs,zcor)
          ELSE
             IF (zt .LT. t_coup) THEN
                zqs = qsats(zt) / pplay(i,k)
                zdqs = dqsats(zt,zqs)
             ELSE
                zqs = qsatl(zt) / pplay(i,k)
                zdqs = dqsatl(zt,zqs)
             ENDIF
          ENDIF
!
!           calculer la fraction nuageuse (processus humide):
!
          zfr = (zq+ratqs*zq-zqs) / (2.0*ratqs*zq)
          zfr = MAX(0.0,MIN(1.0,zfr))
          IF (.NOT.richum) zfr = 0.0
!
!           calculer le nombre de Richardson:
!
          IF (tvirtu) THEN
             ztvd =( t(i,k) &
                  + zdphi/RCPD/(1.+RVTMP2*zq) &
                  *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                  )*(1.+RETV*q(i,k))
             ztvu =( t(i,k-1) &
                  - zdphi/RCPD/(1.+RVTMP2*zq) &
                  *( (1.-zfr) + zfr*(1.+RLVTT*zqs/RD/zt)/(1.+zdqs) ) &
                  )*(1.+RETV*q(i,k-1))
             zri(i) =zmgeom(i)*(ztvd-ztvu)/(zdu2*0.5*(ztvd+ztvu))
             zri(i) = zri(i) &
                  + zmgeom(i)*zmgeom(i)/RG*gamt(k) &
                  *(paprs(i,k)/101325.0)**RKAPPA &
                  /(zdu2*0.5*(ztvd+ztvu))

          ELSE ! calcul de Ridchardson compatible LMD5

             zri(i) =(RCPD*(t(i,k)-t(i,k-1)) &
                  -RD*0.5*(t(i,k)+t(i,k-1))/paprs(i,k) &
                  *(pplay(i,k)-pplay(i,k-1)) &
                  )*zmgeom(i)/(zdu2*0.5*RCPD*(t(i,k-1)+t(i,k)))
             zri(i) = zri(i) + &
                  zmgeom(i)*zmgeom(i)*gamt(k)/RG &
                  *(paprs(i,k)/101325.0)**RKAPPA &
                  /(zdu2*0.5*(t(i,k-1)+t(i,k)))
          ENDIF
!
!           finalement, les coefficients d'echange sont obtenus:
!
          zcdn=SQRT(zdu2) / zmgeom(i) * RG

          IF (opt_ec) THEN
             z2geomf=zgeop(i,k-1)+zgeop(i,k)
             zalm2=(0.5*ckap/RG*z2geomf &
                  /(1.+0.5*ckap/rg/clam*z2geomf))**2
             zalh2=(0.5*ckap/rg*z2geomf &
                  /(1.+0.5*ckap/RG/(clam*SQRT(1.5*cd))*z2geomf))**2
             IF (zri(i).LT.0.0) THEN  ! situation instable
                zscf = ((zgeop(i,k)/zgeop(i,k-1))**(1./3.)-1.)**3 &
                     / (zmgeom(i)/RG)**3 / (zgeop(i,k-1)/RG)
                zscf = SQRT(-zri(i)*zscf)
                zscfm = 1.0 / (1.0+3.0*cb*cc*zalm2*zscf)
                zscfh = 1.0 / (1.0+3.0*cb*cc*zalh2*zscf)
                pcfm(i,k)=zcdn*zalm2*(1.-2.0*cb*zri(i)*zscfm)
                pcfh(i,k)=zcdn*zalh2*(1.-3.0*cb*zri(i)*zscfh)
             ELSE ! situation stable
                zscf=SQRT(1.+cd*zri(i))
                pcfm(i,k)=zcdn*zalm2/(1.+2.0*cb*zri(i)/zscf)
                pcfh(i,k)=zcdn*zalh2/(1.+3.0*cb*zri(i)*zscf)
             ENDIF
          ELSE
             zl2(i)=(mixlen*MAX(0.0,(paprs(i,k)-paprs(i,itop(i)+1)) &
                  /(paprs(i,2)-paprs(i,itop(i)+1)) ))**2
             pcfm(i,k)=SQRT(MAX(zcdn*zcdn*(ric-zri(i))/ric, kstable))
             pcfm(i,k)= zl2(i)* pcfm(i,k)
             pcfh(i,k) = pcfm(i,k) /prandtl ! h et m different
          ENDIF
       ENDDO
    ENDDO

!
! Au-dela du sommet, pas de diffusion turbulente:
!
    DO i = 1, knon
       IF (itop(i)+1 .LE. klev) THEN
          DO k = itop(i)+1, klev
             pcfh(i,k) = 0.0
             pcfm(i,k) = 0.0
          ENDDO
       ENDIF
    ENDDO
      
  END SUBROUTINE coefkz
!
!****************************************************************************************
!
  SUBROUTINE coefkz2(nsrf, knon, paprs, pplay,t, &
       pcfm, pcfh)

    USE dimphy

!======================================================================
! J'introduit un peu de diffusion sauf dans les endroits
! ou une forte inversion est presente
! On peut dire qu'il represente la convection peu profonde
!
! Arguments:
! nsrf-----input-I- indicateur de la nature du sol
! knon-----input-I- nombre de points a traiter
! paprs----input-R- pression a chaque intercouche (en Pa)
! pplay----input-R- pression au milieu de chaque couche (en Pa)
! t--------input-R- temperature (K)
!
! pcfm-----output-R- coefficients a calculer (vitesse)
! pcfh-----output-R- coefficients a calculer (chaleur et humidite)
!======================================================================
!
! Arguments:
!
    INTEGER, INTENT(IN)                       :: knon, nsrf
    REAL, DIMENSION(klon, klev+1), INTENT(IN) ::  paprs
    REAL, DIMENSION(klon, klev), INTENT(IN)   ::  pplay
    REAL, DIMENSION(klon, klev), INTENT(IN)   :: t(klon,klev)
    
    REAL, DIMENSION(klon, klev), INTENT(OUT)  :: pcfm, pcfh
!
! Quelques constantes et options:
!
    REAL, PARAMETER :: prandtl=0.4
    REAL, PARAMETER :: kstable=0.002
!   REAL, PARAMETER :: kstable=0.001
    REAL, PARAMETER :: mixlen=35.0 ! constante controlant longueur de melange
    REAL, PARAMETER :: seuil=-0.02 ! au-dela l'inversion est consideree trop faible
!    PARAMETER (seuil=-0.04)
!    PARAMETER (seuil=-0.06)
!    PARAMETER (seuil=-0.09)

!
! Variables locales:
!
    INTEGER i, k, invb(knon)
    REAL zl2(knon)
    REAL zdthmin(knon), zdthdp

    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"
!
! Initialiser les sorties
!
    DO k = 1, klev
       DO i = 1, knon
          pcfm(i,k) = 0.0
          pcfh(i,k) = 0.0
       ENDDO
    ENDDO

!
! Chercher la zone d'inversion forte
!
    DO i = 1, knon
       invb(i) = klev
       zdthmin(i)=0.0
    ENDDO
    DO k = 2, klev/2-1
       DO i = 1, knon
          zdthdp = (t(i,k)-t(i,k+1))/(pplay(i,k)-pplay(i,k+1)) &
               - RD * 0.5*(t(i,k)+t(i,k+1))/RCPD/paprs(i,k+1)
          zdthdp = zdthdp * 100.0
          IF (pplay(i,k).GT.0.8*paprs(i,1) .AND. &
               zdthdp.LT.zdthmin(i) ) THEN
             zdthmin(i) = zdthdp
             invb(i) = k
          ENDIF
       ENDDO
    ENDDO

!
! Introduire une diffusion:
!
    IF ( nsrf.EQ.is_oce ) THEN
       DO k = 2, klev
          DO i = 1, knon
!IM cf FH/GK   IF ( (nsrf.NE.is_oce) .OR.  ! si ce n'est pas sur l'ocean
!IM cf FH/GK  .     (invb(i).EQ.klev) .OR. ! s'il n'y a pas d'inversion
      !IM cf JLD/ GKtest TERkz2
      ! IF ( (nsrf.EQ.is_ter) .OR.  ! si on est sur la terre
      ! fin GKtest


! s'il n'y a pas d'inversion ou si l'inversion est trop faible
!          IF ( (nsrf.EQ.is_oce) .AND. &
             IF ( (invb(i).EQ.klev) .OR. (zdthmin(i).GT.seuil) ) THEN 
                zl2(i)=(mixlen*MAX(0.0,(paprs(i,k)-paprs(i,klev+1)) &
                     /(paprs(i,2)-paprs(i,klev+1)) ))**2
                pcfm(i,k)= zl2(i)* kstable
                pcfh(i,k) = pcfm(i,k) /prandtl ! h et m different
             ENDIF
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE coefkz2
!
!****************************************************************************************
!
END MODULE coef_diff_turb_mod
