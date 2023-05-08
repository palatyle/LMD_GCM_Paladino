!
!$Id: clcdrag.F90 1279 2009-12-10 09:02:56Z fairhead $
!
SUBROUTINE clcdrag(knon, nsrf, paprs, pplay,&
     u1, v1, t1, q1, &
     tsurf, qsurf, rugos, &
     pcfm, pcfh)

  USE dimphy
  IMPLICIT NONE
! ================================================================= c
!
! Objet : calcul des cdrags pour le moment (pcfm) et 
!         les flux de chaleur sensible et latente (pcfh).   
!
! ================================================================= c
!
! knon----input-I- nombre de points pour un type de surface
! nsrf----input-I- indice pour le type de surface; voir indicesol.h
! u1-------input-R- vent zonal au 1er niveau du modele
! v1-------input-R- vent meridien au 1er niveau du modele
! t1-------input-R- temperature de l'air au 1er niveau du modele
! q1-------input-R- humidite de l'air au 1er niveau du modele
! tsurf------input-R- temperature de l'air a la surface
! qsurf---input-R- humidite de l'air a la surface
! rugos---input-R- rugosite
!
! pcfm---output-R- cdrag pour le moment 
! pcfh---output-R- cdrag pour les flux de chaleur latente et sensible
!
  INTEGER, INTENT(IN)                      :: knon, nsrf
  REAL, DIMENSION(klon,klev+1), INTENT(IN) :: paprs
  REAL, DIMENSION(klon,klev), INTENT(IN)   :: pplay
  REAL, DIMENSION(klon), INTENT(IN)        :: u1, v1, t1, q1
  REAL, DIMENSION(klon), INTENT(IN)        :: tsurf, qsurf
  REAL, DIMENSION(klon), INTENT(IN)        :: rugos
  REAL, DIMENSION(klon), INTENT(OUT)       :: pcfm, pcfh
!
! ================================================================= c
!
  INCLUDE "YOMCST.h"
  INCLUDE "YOETHF.h"
  INCLUDE "indicesol.h"
  INCLUDE "clesphys.h"
!
! Quelques constantes et options:
!!$PB      REAL, PARAMETER :: ckap=0.35, cb=5.0, cc=5.0, cd=5.0, cepdu2=(0.1)**2
  REAL, PARAMETER :: ckap=0.40, cb=5.0, cc=5.0, cd=5.0, cepdu2=(0.1)**2
!
! Variables locales :
  INTEGER               :: i
  REAL                  :: zdu2, ztsolv
  REAL                  :: ztvd, zscf
  REAL                  :: zucf, zcr
  REAL                  :: friv, frih
  REAL, DIMENSION(klon) :: zcfm1, zcfm2
  REAL, DIMENSION(klon) :: zcfh1, zcfh2
  REAL, DIMENSION(klon) :: zcdn
  REAL, DIMENSION(klon) :: zri
  REAL, DIMENSION(klon) :: zgeop1       ! geopotentiel au 1er niveau du modele
  LOGICAL, PARAMETER    :: zxli=.FALSE. ! calcul des cdrags selon Laurent Li
!
! Fonctions thermodynamiques et fonctions d'instabilite
  REAL                  :: fsta, fins, x
  fsta(x) = 1.0 / (1.0+10.0*x*(1+8.0*x))
  fins(x) = SQRT(1.0-18.0*x)

! ================================================================= c
!
! Calculer le geopotentiel du premier couche de modele
!
  DO i = 1, knon
     zgeop1(i) = RD * t1(i) / (0.5*(paprs(i,1)+pplay(i,1))) &
          * (paprs(i,1)-pplay(i,1))
  END DO
! ================================================================= c
!
! Calculer le frottement au sol (Cdrag)
!
  DO i = 1, knon
     zdu2 = MAX(cepdu2,u1(i)**2+v1(i)**2)
     ztsolv = tsurf(i) * (1.0+RETV*qsurf(i))
     ztvd = (t1(i)+zgeop1(i)/RCPD/(1.+RVTMP2*q1(i))) &
          *(1.+RETV*q1(i))
     zri(i) = zgeop1(i)*(ztvd-ztsolv)/(zdu2*ztvd)
     zcdn(i) = (ckap/LOG(1.+zgeop1(i)/(RG*rugos(i))))**2

!!$        IF (zri(i) .ge. 0.) THEN      ! situation stable
     IF (zri(i) .GT. 0.) THEN      ! situation stable
        zri(i) = MIN(20.,zri(i))
        IF (.NOT.zxli) THEN
           zscf = SQRT(1.+cd*ABS(zri(i)))
           FRIV = AMAX1(1. / (1.+2.*CB*zri(i)/ZSCF), 0.1)
           zcfm1(i) = zcdn(i) * FRIV
           FRIH = AMAX1(1./ (1.+3.*CB*zri(i)*ZSCF), 0.1 )
!!$  PB          zcfh1(i) = zcdn(i) * FRIH
!!$ PB           zcfh1(i) = f_cdrag_stable * zcdn(i) * FRIH
           zcfh1(i) = f_cdrag_ter * zcdn(i) * FRIH
           IF(nsrf.EQ.is_oce) zcfh1(i) = f_cdrag_oce * zcdn(i) * FRIH
!!$ PB
           pcfm(i) = zcfm1(i)
           pcfh(i) = zcfh1(i)
        ELSE
           pcfm(i) = zcdn(i)* fsta(zri(i))
           pcfh(i) = zcdn(i)* fsta(zri(i))
        ENDIF
     ELSE                          ! situation instable
        IF (.NOT.zxli) THEN
           zucf = 1./(1.+3.0*cb*cc*zcdn(i)*SQRT(ABS(zri(i)) &
                *(1.0+zgeop1(i)/(RG*rugos(i)))))
           zcfm2(i) = zcdn(i)*amax1((1.-2.0*cb*zri(i)*zucf),0.1)
!!$PB            zcfh2(i) = zcdn(i)*amax1((1.-3.0*cb*zri(i)*zucf),0.1)
           zcfh2(i) = f_cdrag_ter*zcdn(i)*amax1((1.-3.0*cb*zri(i)*zucf),0.1)
           pcfm(i) = zcfm2(i)
           pcfh(i) = zcfh2(i)
        ELSE
           pcfm(i) = zcdn(i)* fins(zri(i))
           pcfh(i) = zcdn(i)* fins(zri(i))
        ENDIF
        zcr = (0.0016/(zcdn(i)*SQRT(zdu2)))*ABS(ztvd-ztsolv)**(1./3.)
        IF(nsrf.EQ.is_oce) pcfh(i) =f_cdrag_oce* zcdn(i)*(1.0+zcr**1.25)**(1./1.25)
     ENDIF
  END DO

! ================================================================= c
     
  ! IM cf JLD : on seuille cdrag_m et cdrag_h
  IF (nsrf == is_oce) THEN
     DO i=1,knon
        pcfm(i)=MIN(pcfm(i),cdmmax)
        pcfh(i)=MIN(pcfh(i),cdhmax)
     END DO
  END IF

END SUBROUTINE clcdrag
