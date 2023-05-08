!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/clmain.F,v 1.3 2005/02/07 16:41:35 fairhead Exp $
!
c
c
      SUBROUTINE clmain_ideal(dtime,itap,
     .                  t,u,v,
     .                  rmu0, 
     .                  ts,
     .                  ftsoil,
     .                  paprs,pplay,ppk,radsol,albe,
     .                  solsw, sollw, sollwdown, fder,
     .                  rlon, rlat, cufi, cvfi, 
     .                  debut, lafin, 
     .                  d_t,d_u,d_v,d_ts,
     .                  flux_t,flux_u,flux_v,cdragh,cdragm,
     .                  dflux_t,
     .                  zcoefh,zu1,zv1) 

c---------------------------------------------------------------
c POUR VENUS
c
c Routine pour une Couche Limite ultra-simple: 
c  - Rayleigh friction dans la couche la plus basse, tau=3Ed=2.6e5s
c  - Kedd=0.15 m^2/s

c S Lebonnois, 10/11/08
c---------------------------------------------------------------
      USE ioipsl
      use dimphy
      use mod_grid_phy_lmdz, only: nbp_lev
      IMPLICIT none
c======================================================================
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c Objet: interface de "couche limite" (diffusion verticale)
c Arguments:
c dtime----input-R- interval du temps (secondes)
c itap-----input-I- numero du pas de temps
c t--------input-R- temperature (K)
c u--------input-R- vitesse u
c v--------input-R- vitesse v
c ts-------input-R- temperature du sol (en Kelvin)
c paprs----input-R- pression a intercouche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c radsol---input-R- flux radiatif net (positif vers le sol) en W/m**2
c rlat-----input-R- latitude en degree
c cufi-----input-R- resolution des mailles en x (m)
c cvfi-----input-R- resolution des mailles en y (m)
c
c d_t------output-R- le changement pour "t"
c d_u------output-R- le changement pour "u"
c d_v------output-R- le changement pour "v"
c d_ts-----output-R- le changement pour "ts"
c flux_t---output-R- flux de chaleur sensible (CpT) J/m**2/s (W/m**2)
c                    (orientation positive vers le bas)
c flux_u---output-R- tension du vent X: (kg m/s)/(m**2 s) ou Pascal
c flux_v---output-R- tension du vent Y: (kg m/s)/(m**2 s) ou Pascal
c dflux_t derive du flux sensible
cAA on rajoute en output yu1 et yv1 qui sont les vents dans 
cAA la premiere couche
c======================================================================
c$$$ PB ajout pour soil
#include "dimsoil.h"
#include "iniprint.h"
#include "clesphys.h"
#include "compbl.h"
c
      REAL dtime
      integer itap
      REAL t(klon,klev)
      REAL u(klon,klev), v(klon,klev)
      REAL paprs(klon,klev+1), pplay(klon,klev), radsol(klon)
! ADAPTATION GCM POUR CP(T)
      real ppk(klon,klev)
      REAL rlon(klon), rlat(klon), cufi(klon), cvfi(klon)
      REAL d_t(klon, klev)
      REAL d_u(klon, klev), d_v(klon, klev)
      REAL flux_t(klon,klev)
      REAL dflux_t(klon)

      REAL flux_u(klon,klev), flux_v(klon,klev)
      REAL cdragh(klon), cdragm(klon)
      real rmu0(klon)         ! cosinus de l'angle solaire zenithal
      LOGICAL debut, lafin
c
      REAL ts(klon)
      REAL d_ts(klon)
      REAL albe(klon)
C
      REAL fder(klon)
      REAL sollw(klon), solsw(klon), sollwdown(klon)
cAA
      REAL zcoefh(klon,klev)
      REAL zu1(klon)
      REAL zv1(klon)
cAA
c$$$ PB ajout pour soil
      REAL ftsoil(klon,nsoilmx)
      REAL ytsoil(klon,nsoilmx)
c======================================================================
      EXTERNAL clqh, clvent, coefkz
c======================================================================
      REAL yts(klon)
      REAL yalb(klon)
      REAL yu1(klon), yv1(klon)
      real ysollw(klon), ysolsw(klon), ysollwdown(klon)
      real yfder(klon), ytaux(klon), ytauy(klon)
      REAL yrads(klon)
C
      REAL y_d_ts(klon)
      REAL y_d_t(klon, klev)
      REAL y_d_u(klon, klev), y_d_v(klon, klev)
      REAL y_flux_t(klon,klev)
      REAL y_flux_u(klon,klev), y_flux_v(klon,klev)
      REAL y_dflux_t(klon)
      REAL ycoefh(klon,klev), ycoefm(klon,klev)
      REAL yu(klon,klev), yv(klon,klev)
      REAL yt(klon,klev)
      REAL ypaprs(klon,klev+1), ypplay(klon,klev), ydelp(klon,klev)
c
      REAL ycoefm0(klon,klev), ycoefh0(klon,klev)

      real yzlay(klon,klev),yzlev(klon,klev+1)
      real yteta(klon,klev)
      real ykmm(klon,klev+1),ykmn(klon,klev+1)
      real ykmq(klon,klev+1)
      real yustar(klon),y_cd_m(klon),y_cd_h(klon)
c
#include "YOMCST.h"
      REAL u1lay(klon), v1lay(klon)
      REAL delp(klon,klev)
      INTEGER i, k
      INTEGER ni(klon), knon, j
      
c======================================================================
      REAL zx_alf1, zx_alf2 !valeur ambiante par extrapola.
c======================================================================
c
      LOGICAL zxli ! utiliser un jeu de fonctions simples
      PARAMETER (zxli=.FALSE.)
c
      REAL zt, zdelta, zcor
C
      real taurelax

c=========================================================
c DEBUT
c=========================================================
          
      DO k = 1, klev   ! epaisseur de couche
      DO i = 1, klon
         delp(i,k) = paprs(i,k)-paprs(i,k+1)
      ENDDO
      ENDDO
      DO i = 1, klon  ! vent de la premiere couche
ccc         zx_alf1 = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
         zx_alf1 = 1.0
         zx_alf2 = 1.0 - zx_alf1
         u1lay(i) = u(i,1)*zx_alf1 + u(i,2)*zx_alf2
         v1lay(i) = v(i,1)*zx_alf1 + v(i,2)*zx_alf2
      ENDDO
c
c initialisation:
c
      DO i = 1, klon
         cdragh(i) = 0.0
         cdragm(i) = 0.0
         dflux_t(i) = 0.0
         zu1(i) = 0.0
         zv1(i) = 0.0
      ENDDO
      yts = 0.0
      yalb = 0.0
      yfder = 0.0
      ytaux = 0.0
      ytauy = 0.0
      ysolsw = 0.0
      ysollw = 0.0
      ysollwdown = 0.0
      yu1 = 0.0
      yv1 = 0.0
      yrads = 0.0
      ypaprs = 0.0
      ypplay = 0.0
      ydelp = 0.0
      yu = 0.0
      yv = 0.0
      yt = 0.0
      y_flux_u = 0.0
      y_flux_v = 0.0
      y_d_ts = 0.0
      y_d_t = 0.0
      y_d_u = 0.0 
      y_d_v = 0.0
      y_flux_t = 0.0
C$$ PB
      y_dflux_t = 0.0
      ytsoil = 999999.
      DO i = 1, klon
         d_ts(i) = 0.0
      ENDDO
      flux_t = 0.
      flux_u = 0.
      flux_v = 0.
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
         zcoefh(i,k) = 0.0
      ENDDO
      ENDDO
c
c chercher les indices:
      DO j = 1, klon
         ni(j) = j
      ENDDO
      knon = klon

      DO j = 1, knon
      i = ni(j)
        yts(j) = ts(i)
        yalb(j) = albe(i)
        yfder(j) = fder(i)
        ytaux(j) = flux_u(i,1)
        ytauy(j) = flux_v(i,1)
        ysolsw(j) = solsw(i)
        ysollw(j) = sollw(i)
        ysollwdown(j) = sollwdown(i)
        yu1(j) = u1lay(i)
        yv1(j) = v1lay(i)
        yrads(j) =  ysolsw(j)+ ysollw(j)
        ypaprs(j,klev+1) = paprs(i,klev+1)
      END DO
C
c$$$ PB ajour pour soil
      DO k = 1, nsoilmx
        DO j = 1, knon
          i = ni(j)
          ytsoil(j,k) = ftsoil(i,k)
        END DO  
      END DO 
      DO k = 1, klev
      DO j = 1, knon
      i = ni(j)
        ypaprs(j,k) = paprs(i,k)
        ypplay(j,k) = pplay(i,k)
        ydelp(j,k) = delp(i,k)
        yu(j,k) = u(i,k)
        yv(j,k) = v(i,k)
        yt(j,k) = t(i,k)
      ENDDO
      ENDDO
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c RAYLEIGH FRICTION (implicit scheme) dans 1ere couche
c Ref: thèse de C. Lee Oxford 2006
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      taurelax = 2.6e5
      yu1 = yu1 / (1+dtime/taurelax)
      yv1 = yv1 / (1+dtime/taurelax)
      yu(:,1) = yu(:,1) / (1+dtime/taurelax)
      yv(:,1) = yv(:,1) / (1+dtime/taurelax)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Coefficient de diffusion verticale
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ycoefm = 0.15

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculer la diffusion des vitesses "u" et "v"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CALL clvent(knon,dtime,yu1,yv1,ycoefm,yt,yu,ypaprs,ypplay,ydelp,
     s            y_d_u,y_flux_u)
      CALL clvent(knon,dtime,yu1,yv1,ycoefm,yt,yv,ypaprs,ypplay,ydelp,
     s            y_d_v,y_flux_v)

c pour le couplage
      ytaux = y_flux_u(:,1)
      ytauy = y_flux_v(:,1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c pas de diffusion de "q" et de "h"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ycoefh = 0.

c=========================
c FIN: tendances
c=========================

      DO j = 1, knon
         i = ni(j)
         d_ts(i) = y_d_ts(j)
         albe(i) = yalb(j)
         cdragh(i) = cdragh(i) + ycoefh(j,1)
         cdragm(i) = cdragm(i) + ycoefm(j,1)
         dflux_t(i) = dflux_t(i) + y_dflux_t(j)
         zu1(i) = zu1(i) + yu1(j)
         zv1(i) = zv1(i) + yv1(j)
      END DO

c$$$ PB ajout pour soil
      DO k = 1, nsoilmx
        DO j = 1, knon
         i = ni(j)
         ftsoil(i, k) = ytsoil(j,k)
        ENDDO
      END DO
      
      DO k = 1, klev
        DO j = 1, knon
         i = ni(j)
         flux_t(i,k) = y_flux_t(j,k)
         flux_u(i,k) = y_flux_u(j,k)
         flux_v(i,k) = y_flux_v(j,k)
         d_t(i,k) = d_t(i,k) + y_d_t(j,k)
         d_u(i,k) = d_u(i,k) + y_d_u(j,k)
         d_v(i,k) = d_v(i,k) + y_d_v(j,k)
         zcoefh(i,k) = zcoefh(i,k) + ycoefh(j,k)
        ENDDO
      ENDDO

c --------------------
c TEST!!!!! PAS DE MELANGE PAR TURBULENCE !!!
c       d_u = 0. 
c       d_v = 0.
c       flux_u = 0.
c       flux_v = 0.
c --------------------

c     print*,"y_d_t apres clqh=",y_d_t(klon/2,:)

      RETURN
      END

