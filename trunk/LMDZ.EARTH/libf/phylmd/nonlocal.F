!
! $Header$
!
C======================================================================
      SUBROUTINE nonlocal(knon, paprs, pplay,
     .                    tsol,beta,u,v,t,q,
     .                    cd_h, cd_m, pcfh, pcfm, cgh, cgq)
      USE dimphy
      IMPLICIT none
c======================================================================
c Laurent Li (LMD/CNRS), le 30 septembre 1998
c Couche limite non-locale. Adaptation du code du CCM3.
c Code non teste, donc a ne pas utiliser.
c======================================================================
c Nonlocal scheme that determines eddy diffusivities based on a
c diagnosed boundary layer height and a turbulent velocity scale.
c Also countergradient effects for heat and moisture are included.
c
c For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
c Local versus nonlocal boundary-layer diffusion in a global climate
c model. J. of Climate, vol. 6, 1825-1842.
c======================================================================
#include "YOMCST.h"
#include "iniprint.h"
c
c Arguments:
c
      INTEGER knon ! nombre de points a calculer
      REAL tsol(klon) ! temperature du sol (K)
      REAL beta(klon) ! efficacite d'evaporation (entre 0 et 1)
      REAL paprs(klon,klev+1) ! pression a inter-couche (Pa)
      REAL pplay(klon,klev)   ! pression au milieu de couche (Pa)
      REAL u(klon,klev) ! vitesse U (m/s)
      REAL v(klon,klev) ! vitesse V (m/s)
      REAL t(klon,klev) ! temperature (K)
      REAL q(klon,klev) ! vapeur d'eau (kg/kg)
      REAL cd_h(klon) ! coefficient de friction au sol pour chaleur
      REAL cd_m(klon) ! coefficient de friction au sol pour vitesse
c
      INTEGER isommet
      REAL vk
      PARAMETER (vk=0.40)
      REAL ricr
      PARAMETER (ricr=0.4)
      REAL fak
      PARAMETER (fak=8.5)
      REAL fakn
      PARAMETER (fakn=7.2)
      REAL onet
      PARAMETER (onet=1.0/3.0)
      REAL t_coup
      PARAMETER(t_coup=273.15)
      REAL zkmin
      PARAMETER (zkmin=0.01)
      REAL betam
      PARAMETER (betam=15.0)
      REAL betah
      PARAMETER (betah=15.0)
      REAL betas
      PARAMETER (betas=5.0)
      REAL sffrac
      PARAMETER (sffrac=0.1)
      REAL binm
      PARAMETER (binm=betam*sffrac)
      REAL binh
      PARAMETER (binh=betah*sffrac)
      REAL ccon
      PARAMETER (ccon=fak*sffrac*vk)
c
      REAL z(klon,klev)
      REAL pcfm(klon,klev), pcfh(klon,klev)
c
      INTEGER i, k
      REAL zxt, zxq, zxu, zxv, zxmod, taux, tauy
      REAL zx_alf1, zx_alf2 ! parametres pour extrapolation
      REAL khfs(klon)       ! surface kinematic heat flux [mK/s]
      REAL kqfs(klon)       ! sfc kinematic constituent flux [m/s]
      REAL heatv(klon)      ! surface virtual heat flux
      REAL ustar(klon)
      REAL rino(klon,klev)  ! bulk Richardon no. from level to ref lev
      LOGICAL unstbl(klon)  ! pts w/unstbl pbl (positive virtual ht flx)
      LOGICAL stblev(klon)  ! stable pbl with levels within pbl
      LOGICAL unslev(klon)  ! unstbl pbl with levels within pbl
      LOGICAL unssrf(klon)  ! unstb pbl w/lvls within srf pbl lyr
      LOGICAL unsout(klon)  ! unstb pbl w/lvls in outer pbl lyr
      LOGICAL check(klon)   ! True=>chk if Richardson no.>critcal
      REAL pblh(klon)
      REAL cgh(klon,2:klev) ! counter-gradient term for heat [K/m]
      REAL cgq(klon,2:klev) ! counter-gradient term for constituents
      REAL cgs(klon,2:klev) ! counter-gradient star (cg/flux)
      REAL obklen(klon)
      REAL ztvd, ztvu, zdu2
      REAL therm(klon)      ! thermal virtual temperature excess
      REAL phiminv(klon)    ! inverse phi function for momentum
      REAL phihinv(klon)    ! inverse phi function for heat
      REAL wm(klon)         ! turbulent velocity scale for momentum
      REAL fak1(klon)       ! k*ustar*pblh
      REAL fak2(klon)       ! k*wm*pblh
      REAL fak3(klon)       ! fakn*wstr/wm
      REAL pblk(klon)       ! level eddy diffusivity for momentum
      REAL pr(klon)         ! Prandtl number for eddy diffusivities
      REAL zl(klon)         ! zmzp / Obukhov length
      REAL zh(klon)         ! zmzp / pblh
      REAL zzh(klon)        ! (1-(zmzp/pblh))**2
      REAL wstr(klon)       ! w*, convective velocity scale
      REAL zm(klon)         ! current level height
      REAL zp(klon)         ! current level height + one level up
      REAL zcor, zdelta, zcvm5, zxqs
      REAL fac, pblmin, zmzp, term
c
#include "YOETHF.h"
#include "FCTTRE.h"
c
c Initialisation
c
      isommet=klev

      DO i = 1, klon
         pcfh(i,1) = cd_h(i)
         pcfm(i,1) = cd_m(i)
      ENDDO
      DO k = 2, klev
      DO i = 1, klon
         pcfh(i,k) = zkmin
         pcfm(i,k) = zkmin
         cgs(i,k) = 0.0
         cgh(i,k) = 0.0
         cgq(i,k) = 0.0
      ENDDO
      ENDDO
c
c Calculer les hauteurs de chaque couche
c
      DO i = 1, knon
         z(i,1) = RD * t(i,1) / (0.5*(paprs(i,1)+pplay(i,1)))
     .               * (paprs(i,1)-pplay(i,1)) / RG
      ENDDO
      DO k = 2, klev
      DO i = 1, knon
         z(i,k) = z(i,k-1)
     .              + RD * 0.5*(t(i,k-1)+t(i,k)) / paprs(i,k)
     .                   * (pplay(i,k-1)-pplay(i,k)) / RG
      ENDDO
      ENDDO
c
      DO i = 1, knon
         IF (thermcep) THEN
           zdelta=MAX(0.,SIGN(1.,RTT-tsol(i)))
           zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
           zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q(i,1))
           zxqs= r2es * FOEEW(tsol(i),zdelta)/paprs(i,1)
           zxqs=MIN(0.5,zxqs)
           zcor=1./(1.-retv*zxqs)
           zxqs=zxqs*zcor
         ELSE
           IF (tsol(i).LT.t_coup) THEN
              zxqs = qsats(tsol(i)) / paprs(i,1)
           ELSE
              zxqs = qsatl(tsol(i)) / paprs(i,1)
           ENDIF
         ENDIF
        zx_alf1 = 1.0
        zx_alf2 = 1.0 - zx_alf1
        zxt = (t(i,1)+z(i,1)*RG/RCPD/(1.+RVTMP2*q(i,1)))
     .        *(1.+RETV*q(i,1))*zx_alf1
     .      + (t(i,2)+z(i,2)*RG/RCPD/(1.+RVTMP2*q(i,2)))
     .        *(1.+RETV*q(i,2))*zx_alf2
        zxu = u(i,1)*zx_alf1+u(i,2)*zx_alf2
        zxv = v(i,1)*zx_alf1+v(i,2)*zx_alf2
        zxq = q(i,1)*zx_alf1+q(i,2)*zx_alf2
        zxmod = 1.0+SQRT(zxu**2+zxv**2)
        khfs(i) = (tsol(i)*(1.+RETV*q(i,1))-zxt) *zxmod*cd_h(i)
        kqfs(i) = (zxqs-zxq) *zxmod*cd_h(i) * beta(i)
        heatv(i) = khfs(i) + 0.61*zxt*kqfs(i)
        taux = zxu *zxmod*cd_m(i)
        tauy = zxv *zxmod*cd_m(i)
        ustar(i) = SQRT(taux**2+tauy**2)
        ustar(i) = MAX(SQRT(ustar(i)),0.01)
      ENDDO
c
      DO i = 1, knon
         rino(i,1) = 0.0
         check(i) = .TRUE.
         pblh(i) = z(i,1)
         obklen(i) = -t(i,1)*ustar(i)**3/(RG*vk*heatv(i))
      ENDDO

C
C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C
      fac = 100.0
      DO k = 1, isommet
      DO i = 1, knon
      IF (check(i)) THEN
         zdu2 = (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
         zdu2 = max(zdu2,1.0e-20)
         ztvd =(t(i,k)+z(i,k)*0.5*RG/RCPD/(1.+RVTMP2*q(i,k)))
     .         *(1.+RETV*q(i,k))
         ztvu =(t(i,1)-z(i,k)*0.5*RG/RCPD/(1.+RVTMP2*q(i,1)))
     .         *(1.+RETV*q(i,1))
         rino(i,k) = (z(i,k)-z(i,1))*RG*(ztvd-ztvu)
     .               /(zdu2*0.5*(ztvd+ztvu))
         IF (rino(i,k).GE.ricr) THEN
           pblh(i) = z(i,k-1) + (z(i,k-1)-z(i,k)) *
     .              (ricr-rino(i,k-1))/(rino(i,k-1)-rino(i,k))
           check(i) = .FALSE.
         ENDIF
      ENDIF
      ENDDO
      ENDDO

C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      DO i = 1, knon
        if (check(i)) pblh(i) = z(i,isommet)
      ENDDO
C
C Improve estimate of pbl height for the unstable points.
C Find unstable points (sensible heat flux is upward):
C
      DO i = 1, knon
      IF (heatv(i) .GT. 0.) THEN
        unstbl(i) = .TRUE.
        check(i) = .TRUE.
      ELSE
        unstbl(i) = .FALSE.
        check(i) = .FALSE.
      ENDIF
      ENDDO
C
C For the unstable case, compute velocity scale and the
C convective temperature excess:
C
      DO i = 1, knon
      IF (check(i)) THEN
        phiminv(i) = (1.-binm*pblh(i)/obklen(i))**onet
        wm(i)= ustar(i)*phiminv(i)
        therm(i) = heatv(i)*fak/wm(i)
        rino(i,1) = 0.0
      ENDIF
      ENDDO
C
C Improve pblh estimate for unstable conditions using the
C convective temperature excess:
C
      DO k = 1, isommet
      DO i = 1, knon
      IF (check(i)) THEN
         zdu2 = (u(i,k)-u(i,1))**2+(v(i,k)-v(i,1))**2+fac*ustar(i)**2
         zdu2 = max(zdu2,1.0e-20)
         ztvd =(t(i,k)+z(i,k)*0.5*RG/RCPD/(1.+RVTMP2*q(i,k)))
     .         *(1.+RETV*q(i,k))
         ztvu =(t(i,1)+therm(i)-z(i,k)*0.5*RG/RCPD/(1.+RVTMP2*q(i,1)))
     .         *(1.+RETV*q(i,1))
         rino(i,k) = (z(i,k)-z(i,1))*RG*(ztvd-ztvu)
     .               /(zdu2*0.5*(ztvd+ztvu))
         IF (rino(i,k).GE.ricr) THEN
           pblh(i) = z(i,k-1) + (z(i,k-1)-z(i,k)) *
     .              (ricr-rino(i,k-1))/(rino(i,k-1)-rino(i,k))
           check(i) = .FALSE.
         ENDIF
      ENDIF
      ENDDO
      ENDDO
C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      DO i = 1, knon
        if (check(i)) pblh(i) = z(i,isommet)
      ENDDO
C
C Points for which pblh exceeds number of pbl layers allowed;
C set to maximum
C
      DO i = 1, knon
        IF (check(i)) pblh(i) = z(i,isommet)
      ENDDO
C
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin
C and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
C
      DO i = 1, knon
        pblmin  = 700.0*ustar(i)
        pblh(i) = MAX(pblh(i),pblmin)
      ENDDO
C
C pblh is now available; do preparation for diffusivity calculation:
C
      DO i = 1, knon
        pblk(i) = 0.0
        fak1(i) = ustar(i)*pblh(i)*vk
C
C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.
C
        IF (unstbl(i)) THEN
          zxt=(t(i,1)-z(i,1)*0.5*RG/RCPD/(1.+RVTMP2*q(i,1)))
     .         *(1.+RETV*q(i,1))
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          phihinv(i) = sqrt(1. - binh*pblh(i)/obklen(i))
          wm(i)      = ustar(i)*phiminv(i)
          fak2(i)    = wm(i)*pblh(i)*vk
          wstr(i)    = (heatv(i)*RG*pblh(i)/zxt)**onet
          fak3(i)    = fakn*wstr(i)/wm(i)
        ENDIF
      ENDDO

C Main level loop to compute the diffusivities and
C counter-gradient terms:
C
      DO 1000 k = 2, isommet
C
C Find levels within boundary layer:
C
        DO i = 1, knon
          unslev(i) = .FALSE.
          stblev(i) = .FALSE.
          zm(i) = z(i,k-1)
          zp(i) = z(i,k)
          IF (zkmin.EQ.0.0 .AND. zp(i).GT.pblh(i)) zp(i) = pblh(i)
          IF (zm(i) .LT. pblh(i)) THEN
            zmzp = 0.5*(zm(i) + zp(i))
            zh(i) = zmzp/pblh(i)
            zl(i) = zmzp/obklen(i)
            zzh(i) = 0.
            IF (zh(i).LE.1.0) zzh(i) = (1. - zh(i))**2
C
C stblev for points zm < plbh and stable and neutral
C unslev for points zm < plbh and unstable
C
            IF (unstbl(i)) THEN
              unslev(i) = .TRUE.
            ELSE
              stblev(i) = .TRUE.
            ENDIF
          ENDIF
        ENDDO
C
C Stable and neutral points; set diffusivities; counter-gradient
C terms zero for stable case:
C
        DO i = 1, knon
          IF (stblev(i)) THEN
            IF (zl(i).LE.1.) THEN
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(1. + betas*zl(i))
            ELSE
              pblk(i) = fak1(i)*zh(i)*zzh(i)/(betas + zl(i))
            ENDIF
            pcfm(i,k) = pblk(i)
            pcfh(i,k) = pcfm(i,k)
          ENDIF
        ENDDO
C
C unssrf, unstable within surface layer of pbl
C unsout, unstable within outer   layer of pbl
C
        DO i = 1, knon
          unssrf(i) = .FALSE.
          unsout(i) = .FALSE.
          IF (unslev(i)) THEN
            IF (zh(i).lt.sffrac) THEN
              unssrf(i) = .TRUE.
            ELSE
              unsout(i) = .TRUE.
            ENDIF
          ENDIF
        ENDDO
C
C Unstable for surface layer; counter-gradient terms zero
C
        DO i = 1, knon
          IF (unssrf(i)) THEN
            term = (1. - betam*zl(i))**onet
            pblk(i) = fak1(i)*zh(i)*zzh(i)*term
            pr(i) = term/sqrt(1. - betah*zl(i))
          ENDIF
        ENDDO
C
C Unstable for outer layer; counter-gradient terms non-zero:
C
        DO i = 1, knon
          IF (unsout(i)) THEN
            pblk(i) = fak2(i)*zh(i)*zzh(i)
            cgs(i,k) = fak3(i)/(pblh(i)*wm(i))
            cgh(i,k) = khfs(i)*cgs(i,k)
            pr(i) = phiminv(i)/phihinv(i) + ccon*fak3(i)/fak
            cgq(i,k) = kqfs(i)*cgs(i,k)
          ENDIF
        ENDDO
C
C For all unstable layers, set diffusivities
C
        DO i = 1, knon
        IF (unslev(i)) THEN
            pcfm(i,k) = pblk(i)
            pcfh(i,k) = pblk(i)/pr(i)
        ENDIF
        ENDDO
 1000 continue           ! end of level loop

      RETURN
      END
