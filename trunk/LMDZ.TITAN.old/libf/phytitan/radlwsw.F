      SUBROUTINE radlwsw(dist, rmu0, fract, zzlev,
     .                  paprs, pplay,tsol, pt, nq, nmicro, pq,qaer)
c      
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
c Objet: interface entre le modele et les rayonnements
c Arguments:
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c falbe----input-R- surface albedo
c zzlev----input-R- altitudes des inter-couches (m)
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c tsol-----input-R- temperature du sol (en K)
c t--------input-R- temperature (K)
c nq-------input-R- nombre de traceurs 
c nmicro---input-R- nombre de traceurs microphysiques
c pq-------input-R- traceurs (rapports de melange)
c heat-----output-R- echauffement atmospherique (visible) (K/s)
c cool-----output-R- refroidissement dans l'IR (K/s)
c radsol---output-R- bilan radiatif net au sol (W/m**2) (+ vers le bas)
c topsw----output-R- flux solaire net au sommet de l'atm. (+ vers le bas)
c toplw----output-R- ray. IR net au sommet de l'atmosphere (+ vers le haut)
c solsw----output-R- flux solaire net a la surface (+ vers le bas)
c sollw----output-R- ray. IR net a la surface (+ vers le bas)
c sollwdown-output-R- ray. IR descendant a la surface (+ vers le bas)
c lwnet____output-R- flux IR net (+ vers le haut)
c swnet____output-R- flux solaire net (+ vers le bas)
c
      
c   S. Lebonnois    05/2008
c  VERSION TITAN 

c======================================================================
      use dimphy
      USE phys_state_var_mod, only: falbe,heat,cool,radsol,
     .      topsw,toplw,solsw,sollw,sollwdown,lwnet,swnet,
     .      lwup,lwdn,swup,swdn
      USE write_field_phy
       IMPLICIT none
#include "YOMCST.h"
#include "clesphys.h" 
c
c ARGUMENTS
      INTEGER nq,nmicro
      real rmu0(klon), fract(klon), dist
c
      real zzlev(klon,klev+1),paprs(klon,klev+1), pplay(klon,klev)
      real tsol(klon)
      real pt(klon,klev)
      real pq(klon,klev,nq)
      REAL qaer(klon,klev,nq)
c
c LOCAL VARIABLES
      integer i,k,l,iq
      real zp(klon,klev+1),zt(klon,klev+1),zz(klon,klev+1)
      real zq(klon,klev,nq)
      real zheatc(klon,klev), zcoolc(klon,klev)
      real zheatp(klon,klev), zcoolp(klon,klev)
      REAL zswupc(klon,klev+1),zlwupc(klon,klev+1)
      REAL zswupp(klon,klev+1),zlwupp(klon,klev+1)
      REAL zswdnc(klon,klev+1),zlwdnc(klon,klev+1)
      REAL zswdnp(klon,klev+1),zlwdnp(klon,klev+1)
      REAL zsollwdownc(klon),zsollwdownp(klon)
      INTEGER icld


c =======================================
c INITIALISATIONS
c =======================================

c   passage au pressions en bar avec indice 1 au sommet.
             do l=2,klev+1
                do i=1,klon
                   zp(i,l)=paprs(i,klev+2-l)*1.e-5
                enddo
             enddo
             do i=1,klon
                zp(i,1)=zp(i,2)*.001
             enddo

c     call WriteField_phy('radlwsw_zp',zp,klev+1)

c =======================================
c   altitudes (m) avec indice 1 en haut
             do l=1,klev+1
                do i=1,klon
                   zz(i,l)=zzlev(i,klev+2-l)
                enddo
             enddo

c   temperatures avec indice 1 en haut
             do l=1,klev
                do i=1,klon
                   zt(i,l)=pt(i,klev+1-l)
                enddo
             enddo
             do i=1,klon
                zt(i,klev+1)=tsol(i)
             enddo

c  traceurs avec indice 1 en haut
             do l=1,klev
                do i=1,klon
                 do iq=1,nq
                   zq(i,l,iq)=pq(i,klev+1-l,iq)
                 enddo
                enddo
             enddo

c =======================================
c CALCUL DES TAU V+IR  (dans des common...)
c =======================================

      print*,'On calcule les opacites'

         CALL radtitan(zp,nq,nmicro,zq,qaer)

c CALCUL DU SW 
c =======================================

      print*,'On calcule le rayonnement SW'

       IF (clouds.eq.1) THEN
         ICLD = 1   ! colonne avec nuages
         CALL heating(dist,rmu0,fract,falbe,zheatc,zswupc,zswdnc,icld)
       ELSE
         zheatc  = 0.
         zswupc = 0.
         zswdnc = 0.
       ENDIF 
       ICLD = 0   ! colonne sans nuages
       CALL heating(dist,rmu0,fract,falbe,zheatp,zswupp,zswdnp,icld)

c inversion de l'axe vertical
       do l=1,klev
         do i=1,klon
           heat(i,l)=zheatc(i,klev+1-l)*xnuf +
     &               zheatp(i,klev+1-l)*(1.-xnuf)
         enddo
       enddo
       do l=1,klev+1
         do i=1,klon
           swup(i,l) =zswupc(i,klev+2-l)*xnuf +
     &                zswupp(i,klev+2-l)*(1.-xnuf)
           swdn(i,l) =zswdnc(i,klev+2-l)*xnuf +
     &                zswdnp(i,klev+2-l)*(1.-xnuf)
           swnet(i,l)=swdn(i,l)-swup(i,l)
         enddo
       enddo

      solsw = swnet(:,1)
      topsw = swnet(:,klev+1)

c =======================================
c CALCUL DU LW
c =======================================

      print*,'On calcule le rayonnement LW'

       IF (clouds.eq.1) THEN
         ICLD = 1
         CALL cooling(klon,klev+1,zp,zt,zz,zcoolc,zlwupc,zlwdnc,
     &   zsollwdownc,icld)
       ELSE
         zcoolc      = 0.
         zlwupc      = 0.
         zlwdnc      = 0.
         zsollwdownc = 0.
       ENDIF
       ICLD = 0
       CALL cooling(klon,klev+1,zp,zt,zz,zcoolp,zlwupp,zlwdnp,
     & zsollwdownp,icld)

c inversion de l'axe vertical
       do l=1,klev
         do i=1,klon
           cool(i,l)=zcoolc(i,klev+1-l)*xnuf +
     &               zcoolp(i,klev+1-l)*(1.-xnuf)
         enddo
       enddo
       do l=1,klev+1
         do i=1,klon
           lwup(i,l) =zlwupc(i,klev+2-l)*xnuf +
     &                zlwupp(i,klev+2-l)*(1.-xnuf)
           lwdn(i,l) =zlwdnc(i,klev+2-l)*xnuf +
     &                zlwdnp(i,klev+2-l)*(1.-xnuf)
           lwnet(i,l)=lwup(i,l)-lwdn(i,l)
         enddo
       enddo
   
       do i=1,klon
         sollwdown(i)=zsollwdownc(i)*xnuf +
     &                zsollwdownp(i)*(1.-xnuf)
       enddo

      sollw  = -lwnet(:,1)
      toplw  = lwnet(:,klev+1)
      radsol = solsw+sollw
      
      RETURN
      END
