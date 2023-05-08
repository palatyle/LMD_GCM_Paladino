!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ajsec.F,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
! ADAPTATION GCM POUR CP(T)
      SUBROUTINE ajsec(paprs, pplay, ppk, tfi, ufi, vfi, nq, qfi, 
     .                             d_tfi, d_ufi, d_vfi, d_qfi)

      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lev
      use cpdet_phy_mod, only: t2tpot, tpot2t
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: ajustement sec (adaptation du GCM du LMD)
c S. Lebonnois, 10/2007:
c melange u et v comme dans convadj (MARS)
c======================================================================
c Arguments:
c tfi-------input-R- Temperature
c ufi-------input-R- vent zonal
c vfi-------input-R- vent meridien
c nq--------input-R- nombre de traceurs
c qfi-------input-R- traceurs
c
c d_tfi-----output-R-Incrementation de la temperature
c d_ufi-----output-R-Incrementation du vent zonal
c d_vfi-----output-R-Incrementation du vent meridien
c d_qfi-----output-R-Incrementation des traceurs
c======================================================================
#include "YOMCST.h"
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL ppk(klon,klev)
      INTEGER nq
      REAL tfi(klon,klev), d_tfi(klon,klev)
      REAL ufi(klon,klev), d_ufi(klon,klev)
      REAL vfi(klon,klev), d_vfi(klon,klev)
      REAL qfi(klon,klev,nq), d_qfi(klon,klev,nq)
c
      INTEGER,save :: limbas, limhau ! les couches a ajuster
c
      REAL zh(klon,klev)
      REAL zu(klon,klev),zv(klon,klev)
      REAL zt(klon,klev),zq(klon,klev,nq)
      REAL zdp(klon,klev)
      REAL zpkdp(klon,klev)
      REAL hm,sm,zum,zvm,zalpha,zqm(nq)
      LOGICAL modif(klon), down
      INTEGER i, k, k1, k2, iq
c
c Initialisation:
c
      limbas=1
      limhau=klev

      DO k = 1, klev
      DO i = 1, klon
         d_tfi(i,k) = 0.0
         d_ufi(i,k) = 0.0
         d_vfi(i,k) = 0.0
         d_qfi(i,k,:) = 0.0
         zu(i,k)    = ufi(i,k)
         zv(i,k)    = vfi(i,k)
         zq(i,k,:)  = qfi(i,k,:)
      ENDDO
      ENDDO
c------------------------------------- passage en temperature potentielle
! ADAPTATION GCM POUR CP(T)
      call t2tpot(klon*nbp_lev,tfi,zh,ppk)
c
      DO k = limbas, limhau
      DO i = 1, klon
         zdp(i,k) = paprs(i,k)-paprs(i,k+1)
         zpkdp(i,k) = ppk(i,k) * zdp(i,k)
      ENDDO
      ENDDO
c
c------------------------------------- detection des profils a modifier
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
 8020         CONTINUE
                sm = sm +zpkdp(i,k)
                hm = hm +zpkdp(i,k) * (zh(i,k)-hm) / sm
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
c------------ et melange partiel des vents
              zalpha=0.
              zum=0.
              zvm=0.
              zqm=0.
              DO k = k1, k2
                zalpha=zalpha+ABS(zh(i,k)-hm)*zdp(i,k)
                zh(i,k) = hm
                zum=zum+ufi(i,k)*zdp(i,k)
                zvm=zvm+vfi(i,k)*zdp(i,k)
                do iq=1,nq
                  zqm(iq)=zqm(iq)+qfi(i,k,iq)*zdp(i,k)
                enddo
              ENDDO
              zalpha=zalpha/(hm*(paprs(i,k1)-paprs(i,k2+1)))
              zum=zum/(paprs(i,k1)-paprs(i,k2+1))
              zvm=zvm/(paprs(i,k1)-paprs(i,k2+1))
              do iq=1,nq
                zqm(iq)=zqm(iq)/(paprs(i,k1)-paprs(i,k2+1))
              enddo

              IF(zalpha.GT.1.) THEN
                 PRINT*,'WARNING dans ajsec zalpha=',zalpha
c         STOP
                 zalpha=1.
              ELSE
c                IF(zalpha.LT.0.) STOP
                 IF(zalpha.LT.1.e-5) zalpha=1.e-4
              ENDIF
c ----------------------------
c TEST --- PAS DE MELANGE DE U ni Q
c             zalpha=0.
c ----------------------------

              DO k=k1,k2
                 zu(i,k)=ufi(i,k)+zalpha*(zum-ufi(i,k))
                 zv(i,k)=vfi(i,k)+zalpha*(zvm-vfi(i,k))
                 do iq=1,nq
                   zq(i,k,iq)=qfi(i,k,iq)+zalpha*(zqm(iq)-qfi(i,k,iq))
                 enddo
              ENDDO
              k2 = k2 + 1
            ENDIF
          GOTO 8000
 8001     CONTINUE
      ENDIF
 1080 CONTINUE
c
c------------------------------------- passage en temperature 
c------------------------------------- et calcul du d_t 
! ADAPTATION GCM POUR CP(T)
      call tpot2t(klon*nbp_lev,zh,zt,ppk)

      DO k = limbas, limhau
      DO i = 1, klon
         d_tfi(i,k) = zt(i,k) - tfi(i,k)
         d_ufi(i,k) = zu(i,k) - ufi(i,k)
         d_vfi(i,k) = zv(i,k) - vfi(i,k)
         do iq=1,nq
           d_qfi(i,k,iq) = zq(i,k,iq) - qfi(i,k,iq)
         enddo 
      ENDDO
      ENDDO
c
      IF (limbas.GT.1) THEN
      DO k = 1, limbas-1
      DO i = 1, klon
         d_tfi(i,k) = 0.0
         d_ufi(i,k) = 0.0
         d_vfi(i,k) = 0.0
         do iq=1,nq
           d_qfi(i,k,iq) = 0.0
         enddo
      ENDDO
      ENDDO
      ENDIF
c
      IF (limhau.LT.klev) THEN
      DO k = limhau+1, klev
      DO i = 1, klon
         d_tfi(i,k) = 0.0
         d_ufi(i,k) = 0.0
         d_vfi(i,k) = 0.0
         do iq=1,nq
           d_qfi(i,k,iq) = 0.0
         enddo
      ENDDO
      ENDDO
      ENDIF
c
      RETURN
      END

