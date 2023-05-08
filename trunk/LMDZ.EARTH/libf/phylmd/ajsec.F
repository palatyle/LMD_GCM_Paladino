!
! $Header$
!
      SUBROUTINE ajsec(paprs, pplay, t,q,limbas,d_t,d_q)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: ajustement sec (adaptation du GCM du LMD)
c======================================================================
c Arguments:
c t-------input-R- Temperature
c
c d_t-----output-R-Incrementation de la temperature
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev), q(klon,klev)
      REAL d_t(klon,klev), d_q(klon,klev)
c
      INTEGER limbas(klon), limhau ! les couches a ajuster
c
      LOGICAL mixq
ccc      PARAMETER (mixq=.TRUE.)
      PARAMETER (mixq=.FALSE.)
c
      REAL zh(klon,klev)
      REAL zho(klon,klev)
      REAL zq(klon,klev)
      REAL zpk(klon,klev)
      REAL zpkdp(klon,klev)
      REAL hm, sm, qm
      LOGICAL modif(klon), down
      INTEGER i, k, k1, k2
c
c Initialisation:
c
cym
      limhau=klev
  
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
c------------------------------------- detection des profils a modifier
      DO k = 1, limhau
      DO i = 1, klon
         zpk(i,k) = pplay(i,k)**RKAPPA
         zh(i,k) = RCPD * t(i,k)/ zpk(i,k)
         zho(i,k) = zh(i,k)
         zq(i,k) = q(i,k)
      ENDDO
      ENDDO
c
      DO k = 1, limhau
      DO i = 1, klon
         zpkdp(i,k) = zpk(i,k) * (paprs(i,k)-paprs(i,k+1))
      ENDDO
      ENDDO
c
      DO i = 1, klon
         modif(i) = .FALSE.
      ENDDO
      DO k = 2, limhau
      DO i = 1, klon
      IF (.NOT.modif(i).and.k-1>limbas(i)) THEN
         IF ( zh(i,k).LT.zh(i,k-1) ) modif(i) = .TRUE.
      ENDIF
      ENDDO
      ENDDO
c------------------------------------- correction des profils instables
      DO 1080 i = 1, klon
      IF (modif(i)) THEN
          k2 = limbas(i)
 8000     CONTINUE
            k2 = k2 + 1
            IF (k2 .GT. limhau) goto 8001
            IF (zh(i,k2) .LT. zh(i,k2-1)) THEN
              k1 = k2 - 1
              k = k1
              sm = zpkdp(i,k2)
              hm = zh(i,k2)
              qm = zq(i,k2)
 8020         CONTINUE
                sm = sm +zpkdp(i,k)
                hm = hm +zpkdp(i,k) * (zh(i,k)-hm) / sm
                qm = qm +zpkdp(i,k) * (zq(i,k)-qm) / sm
                down = .FALSE.
                IF (k1 .ne. limbas(i)) THEN
                  IF (hm .LT. zh(i,k1-1)) down = .TRUE.
                ENDIF
                IF (down) THEN
                  k1 = k1 - 1
                  k = k1
                ELSE
                  IF ((k2 .EQ. limhau)) GOTO 8021
                  IF ((zh(i,k2+1).GE.hm)) GOTO 8021
                  k2 = k2 + 1
                  k = k2
                ENDIF
              GOTO 8020
 8021         CONTINUE
c------------ nouveau profil : constant (valeur moyenne)
              DO k = k1, k2
                zh(i,k) = hm
                zq(i,k) = qm
              ENDDO
              k2 = k2 + 1
            ENDIF
          GOTO 8000
 8001     CONTINUE
      ENDIF
 1080 CONTINUE
c
      DO k = 1, limhau
      DO i = 1, klon
         d_t(i,k) = (zh(i,k)-zho(i,k))*zpk(i,k)/RCPD
         d_q(i,k) = zq(i,k) - q(i,k)
      ENDDO
      ENDDO
c
! FH : les d_q et d_t sont maintenant calcules de facon a valoir
! effectivement 0. si on ne fait rien.
!
!     IF (limbas.GT.1) THEN
!     DO k = 1, limbas-1
!     DO i = 1, klon
!        d_t(i,k) = 0.0
!        d_q(i,k) = 0.0
!     ENDDO
!     ENDDO
!     ENDIF
c
!     IF (limhau.LT.klev) THEN
!     DO k = limhau+1, klev
!     DO i = 1, klon
!        d_t(i,k) = 0.0
!        d_q(i,k) = 0.0
!     ENDDO
!     ENDDO
!     ENDIF
c
      IF (.NOT.mixq) THEN
      DO k = 1, klev
      DO i = 1, klon
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
      ENDIF
c
      RETURN
      END

      SUBROUTINE ajsec_convV2(paprs, pplay, t,q, d_t,d_q)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: ajustement sec (adaptation du GCM du LMD)
c======================================================================
c Arguments:
c t-------input-R- Temperature
c
c d_t-----output-R-Incrementation de la temperature
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev), q(klon,klev)
      REAL d_t(klon,klev), d_q(klon,klev)
c
      INTEGER limbas, limhau ! les couches a ajuster
ccc      PARAMETER (limbas=klev-3, limhau=klev)
cym      PARAMETER (limbas=1, limhau=klev)
c
      LOGICAL mixq
ccc      PARAMETER (mixq=.TRUE.)
      PARAMETER (mixq=.FALSE.)
c
      REAL zh(klon,klev)
      REAL zq(klon,klev)
      REAL zpk(klon,klev)
      REAL zpkdp(klon,klev)
      REAL hm, sm, qm
      LOGICAL modif(klon), down
      INTEGER i, k, k1, k2
c
c Initialisation:
c
cym
      limbas=1
      limhau=klev
  
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
c------------------------------------- detection des profils a modifier
      DO k = limbas, limhau
      DO i = 1, klon
         zpk(i,k) = pplay(i,k)**RKAPPA
         zh(i,k) = RCPD * t(i,k)/ zpk(i,k)
         zq(i,k) = q(i,k)
      ENDDO
      ENDDO
c
      DO k = limbas, limhau
      DO i = 1, klon
         zpkdp(i,k) = zpk(i,k) * (paprs(i,k)-paprs(i,k+1))
      ENDDO
      ENDDO
c
      DO i = 1, klon
         modif(i) = .FALSE.
      ENDDO
      DO k = limbas+1, limhau
      DO i = 1, klon
      IF (.NOT.modif(i)) THEN
         IF ( zh(i,k).LT.zh(i,k-1) ) modif(i) = .TRUE.
      ENDIF
      ENDDO
      ENDDO
c------------------------------------- correction des profils instables
      DO 1080 i = 1, klon
      IF (modif(i)) THEN
          k2 = limbas
 8000     CONTINUE
            k2 = k2 + 1
            IF (k2 .GT. limhau) goto 8001
            IF (zh(i,k2) .LT. zh(i,k2-1)) THEN
              k1 = k2 - 1
              k = k1
              sm = zpkdp(i,k2)
              hm = zh(i,k2)
              qm = zq(i,k2)
 8020         CONTINUE
                sm = sm +zpkdp(i,k)
                hm = hm +zpkdp(i,k) * (zh(i,k)-hm) / sm
                qm = qm +zpkdp(i,k) * (zq(i,k)-qm) / sm
                down = .FALSE.
                IF (k1 .ne. limbas) THEN
                  IF (hm .LT. zh(i,k1-1)) down = .TRUE.
                ENDIF
                IF (down) THEN
                  k1 = k1 - 1
                  k = k1
                ELSE
                  IF ((k2 .EQ. limhau)) GOTO 8021
                  IF ((zh(i,k2+1).GE.hm)) GOTO 8021
                  k2 = k2 + 1
                  k = k2
                ENDIF
              GOTO 8020
 8021         CONTINUE
c------------ nouveau profil : constant (valeur moyenne)
              DO k = k1, k2
                zh(i,k) = hm
                zq(i,k) = qm
              ENDDO
              k2 = k2 + 1
            ENDIF
          GOTO 8000
 8001     CONTINUE
      ENDIF
 1080 CONTINUE
c
      DO k = limbas, limhau
      DO i = 1, klon
         d_t(i,k) = zh(i,k)*zpk(i,k)/RCPD - t(i,k)
         d_q(i,k) = zq(i,k) - q(i,k)
      ENDDO
      ENDDO
c
      IF (limbas.GT.1) THEN
      DO k = 1, limbas-1
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
      ENDIF
c
      IF (limhau.LT.klev) THEN
      DO k = limhau+1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
      ENDIF
c
      IF (.NOT.mixq) THEN
      DO k = 1, klev
      DO i = 1, klon
         d_q(i,k) = 0.0
      ENDDO
      ENDDO
      ENDIF
c
      RETURN
      END
      SUBROUTINE ajsec_old(paprs, pplay, t, d_t)
      USE dimphy
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: ajustement sec (adaptation du GCM du LMD)
c======================================================================
c Arguments:
c t-------input-R- Temperature
c
c d_t-----output-R-Incrementation de la temperature
c======================================================================
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "YOMCST.h"
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev)
      REAL d_t(klon,klev)
c
      REAL local_h(klon,klev)
      REAL hm, sm
      LOGICAL modif(klon), down
      INTEGER i, l, l1, l2
c------------------------------------- detection des profils a modifier
      DO i = 1, klon
         modif(i)   = .false.
      ENDDO
c
      DO l = 1, klev
      DO i = 1, klon
         local_h(i,l) = RCPD * t(i,l)/ (pplay(i,l)**RKAPPA)
      ENDDO
      ENDDO
c
      DO l = 2, klev
      DO i = 1, klon
         IF ( local_h(i,l).lt.local_h(i,l-1) ) THEN
            modif(i) = .true.
         ELSE
            modif(i) = modif(i)
         ENDIF
      ENDDO
      ENDDO
c------------------------------------- correction des profils instables
      do 1080 i = 1, klon
        if (modif(i)) then
          l2 = 1
 8000     continue
            l2 = l2 + 1
            if (l2 .gt. klev) goto 8001
            if (local_h(i, l2) .lt. local_h(i, l2-1)) then
              l1 = l2 - 1
              l  = l1
              sm = pplay(i,l2)**rkappa * (paprs(i,l2)-paprs(i,l2+1))
              hm = local_h(i, l2)
 8020         continue
                sm = sm +pplay(i,l)**rkappa*(paprs(i,l)-paprs(i,l+1))
                hm = hm +pplay(i,l)**rkappa*(paprs(i,l)-paprs(i,l+1))
     .                         * (local_h(i, l) - hm) / sm
                down = .false.
                if (l1 .ne. 1) then
                  if (hm .lt. local_h(i, l1-1)) then
                    down = .true.
                  end if
                end if
                if (down) then
                  l1 = l1 - 1
                  l  = l1
                else
                  if ((l2 .eq. klev)) GOTO 8021
                  IF ((local_h(i, l2+1).ge.hm)) goto 8021
                  l2 = l2 + 1
                  l  = l2
                end if
              go to 8020
 8021         continue
c------------ nouveau profil : constant (valeur moyenne)
              do 1100 l = l1, l2
                local_h(i, l) = hm
 1100         continue
              l2 = l2 + 1
            end if
          go to 8000
 8001     continue
        end if
 1080 continue
c
      DO l = 1, klev
      DO i = 1, klon
         d_t(i,l) = local_h(i,l)*(pplay(i,l)**rkappa)/RCPD - t(i,l)
      ENDDO
      ENDDO
c
      RETURN
      END
