! $Id: newmicro.F 1337 2010-04-02 11:31:05Z fairhead $
!     
      SUBROUTINE newmicro (paprs, pplay,ok_newmicro,
     .                  t, pqlwp, pclc, pcltau, pclemi,
     .                  pch, pcl, pcm, pct, pctlwp,
     s                  xflwp, xfiwp, xflwc, xfiwc,
     e                  ok_aie, 
     e                  mass_solu_aero, mass_solu_aero_pi, 
     e                  bl95_b0, bl95_b1,
     s                  cldtaupi, re, fl, reliq, reice)

      USE dimphy
      USE phys_local_var_mod, only: scdnc,cldncl,reffclwtop,lcc,
     .                              reffclws,reffclwc,cldnvi,lcc3d,
     .                              lcc3dcon,lcc3dstra
      USE phys_state_var_mod, only: rnebcon,clwcon
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930910
c Objet: Calculer epaisseur optique et emmissivite des nuages
c======================================================================
c Arguments:
c t-------input-R-temperature
c pqlwp---input-R-eau liquide nuageuse dans l'atmosphere (kg/kg)
c pclc----input-R-couverture nuageuse pour le rayonnement (0 a 1)
c 
c ok_aie--input-L-apply aerosol indirect effect or not
c mass_solu_aero-----input-R-total mass concentration for all soluble aerosols[ug/m^3]
c mass_solu_aero_pi--input-R-dito, pre-industrial value
c bl95_b0-input-R-a parameter, may be varied for tests (s-sea, l-land)
c bl95_b1-input-R-a parameter, may be varied for tests (    -"-      )
c      
c cldtaupi-output-R-pre-industrial value of cloud optical thickness, 
c                   needed for the diagnostics of the aerosol indirect 
c                   radiative forcing (see radlwsw)
c re------output-R-Cloud droplet effective radius multiplied by fl [um]
c fl------output-R-Denominator to re, introduced to avoid problems in
c                  the averaging of the output. fl is the fraction of liquid
c                  water clouds within a grid cell           
c pcltau--output-R-epaisseur optique des nuages
c pclemi--output-R-emissivite des nuages (0 a 1)
c======================================================================
C
#include "YOMCST.h"
c
cym#include "dimensions.h"
cym#include "dimphy.h"
#include "nuage.h"
cIM cf. CR: include pour NOVLP et ZEPSEC
#include "radepsi.h"
#include "radopt.h"
c choix de l'hypothese de recouvrememnt nuageuse
      LOGICAL RANDOM,MAXIMUM_RANDOM,MAXIMUM,FIRST
      parameter (RANDOM=.FALSE., MAXIMUM_RANDOM=.TRUE., MAXIMUM=.FALSE.)
c Hypoyhese de recouvrement : MAXIMUM_RANDOM
      INTEGER flag_max
      REAL phase3d(klon, klev),dh(klon, klev),pdel(klon, klev),
     .     zrho(klon, klev)
      REAL tcc(klon), ftmp(klon), lcc_integrat(klon), height(klon)
      REAL thres_tau,thres_neb
      PARAMETER (thres_tau=0.3, thres_neb=0.001)
      REAL t_tmp
      REAL gravit
      PARAMETER (gravit=9.80616)  !m/s2
      REAL pqlwpcon(klon, klev), pqlwpstra(klon, klev) 
c
      REAL paprs(klon,klev+1), pplay(klon,klev)
      REAL t(klon,klev)
c
      REAL pclc(klon,klev)
      REAL pqlwp(klon,klev)
      REAL pcltau(klon,klev), pclemi(klon,klev)
c
      REAL pct(klon), pctlwp(klon), pch(klon), pcl(klon), pcm(klon)
c
      LOGICAL lo
c
      REAL cetahb, cetamb
      PARAMETER (cetahb = 0.45, cetamb = 0.80)
C
      INTEGER i, k
cIM: 091003   REAL zflwp, zradef, zfice, zmsac
      REAL zflwp(klon), zradef, zfice, zmsac
cIM: 091003 rajout
      REAL xflwp(klon), xfiwp(klon)
      REAL xflwc(klon,klev), xfiwc(klon,klev)
c
      REAL radius, rad_chaud
cc      PARAMETER (rad_chau1=13.0, rad_chau2=9.0, rad_froid=35.0)
ccc      PARAMETER (rad_chaud=15.0, rad_froid=35.0)
c sintex initial      PARAMETER (rad_chaud=10.0, rad_froid=30.0)
      REAL coef, coef_froi, coef_chau
      PARAMETER (coef_chau=0.13, coef_froi=0.09)
      REAL seuil_neb
      PARAMETER (seuil_neb=0.001)
      INTEGER nexpo ! exponentiel pour glace/eau
      PARAMETER (nexpo=6)
ccc      PARAMETER (nexpo=1)

c -- sb:
      logical ok_newmicro
c     parameter (ok_newmicro=.FALSE.)
cIM: 091003   real rel, tc, rei, zfiwp
      real rel, tc, rei, zfiwp(klon)
      real k_liq, k_ice0, k_ice, DF
      parameter (k_liq=0.0903, k_ice0=0.005) ! units=m2/g
      parameter (DF=1.66) ! diffusivity factor
c sb --
cjq for the aerosol indirect effect
cjq introduced by Johannes Quaas (quaas@lmd.jussieu.fr), 27/11/2003
cjq      
      LOGICAL ok_aie            ! Apply AIE or not?
      LOGICAL ok_a1lwpdep       ! a1 LWP dependent?
      
      REAL mass_solu_aero(klon, klev)    ! total mass concentration for all soluble aerosols [ug m-3]
      REAL mass_solu_aero_pi(klon, klev) ! - " - (pre-industrial value)
      REAL cdnc(klon, klev)     ! cloud droplet number concentration [m-3]
      REAL re(klon, klev)       ! cloud droplet effective radius [um]
      REAL cdnc_pi(klon, klev)     ! cloud droplet number concentration [m-3] (pi value)
      REAL re_pi(klon, klev)       ! cloud droplet effective radius [um] (pi value)
      
      REAL fl(klon, klev)       ! xliq * rneb (denominator to re; fraction of liquid water clouds within the grid cell)
      
      REAL bl95_b0, bl95_b1     ! Parameter in B&L 95-Formula
      
      REAL cldtaupi(klon, klev) ! pre-industrial cloud opt thickness for diag
cjq-end    
cIM cf. CR:parametres supplementaires
      REAL zclear(klon)
      REAL zcloud(klon) 

c **************************
c *                        *
c * DEBUT PARTIE OPTIMISEE *
c *                        *
c **************************

      REAL diff_paprs(klon, klev), zfice1, zfice2(klon, klev)
      REAL rad_chaud_tab(klon, klev), zflwp_var, zfiwp_var

! Abderrahmane oct 2009
      Real reliq(klon, klev), reice(klon, klev)

c
c Calculer l'epaisseur optique et l'emmissivite des nuages
c
c     IM inversion des DO
      xflwp = 0.d0
      xfiwp = 0.d0
      xflwc = 0.d0
      xfiwc = 0.d0

! Initialisation 
      reliq=0.
      reice=0.

      DO k = 1, klev
         DO i = 1, klon
            diff_paprs(i,k) = (paprs(i,k)-paprs(i,k+1))/RG
         ENDDO
      ENDDO

      IF (ok_newmicro) THEN


         DO k = 1, klev
            DO i = 1, klon
c               zfice2(i,k) = 1.0 - (t(i,k)-t_glace) / (273.13-t_glace)
               zfice2(i,k) = 1.0 - (t(i,k)-t_glace_min) / 
     &                             (t_glace_max-t_glace_min)
               zfice2(i,k) = MIN(MAX(zfice2(i,k),0.0),1.0)
c     IM Total Liquid/Ice water content                                    
               xflwc(i,k) = (1.-zfice2(i,k))*pqlwp(i,k)
               xfiwc(i,k) = zfice2(i,k)*pqlwp(i,k)
c     IM In-Cloud Liquid/Ice water content
c     xflwc(i,k) = xflwc(i,k)+(1.-zfice)*pqlwp(i,k)/pclc(i,k)
c     xfiwc(i,k) = xfiwc(i,k)+zfice*pqlwp(i,k)/pclc(i,k)
            ENDDO
         ENDDO

         IF (ok_aie) THEN
            DO k = 1, klev
               DO i = 1, klon
                                ! Formula "D" of Boucher and Lohmann, Tellus, 1995
                                !             
                  cdnc(i,k) = 10.**(bl95_b0+bl95_b1*
     &               log(MAX(mass_solu_aero(i,k),1.e-4))/log(10.))*1.e6 !-m-3
                                ! Cloud droplet number concentration (CDNC) is restricted
                                ! to be within [20, 1000 cm^3]
                                ! 
                  cdnc(i,k)=MIN(1000.e6,MAX(20.e6,cdnc(i,k)))
                                !
                                !
                  cdnc_pi(i,k) = 10.**(bl95_b0+bl95_b1*
     &               log(MAX(mass_solu_aero_pi(i,k),1.e-4))/log(10.))
     &               *1.e6 !-m-3
                  cdnc_pi(i,k)=MIN(1000.e6,MAX(20.e6,cdnc_pi(i,k)))
               ENDDO
            ENDDO
            DO k = 1, klev
               DO i = 1, klon
!                  rad_chaud_tab(i,k) = 
!     &                 MAX(1.1e6 
!     &                 *((pqlwp(i,k)*pplay(i,k)/(RD * T(i,k)))  
!     &                 /(4./3*RPI*1000.*cdnc(i,k)) )**(1./3.),5.)
                  rad_chaud_tab(i,k) = 
     &                 1.1
     &                 *((pqlwp(i,k)*pplay(i,k)/(RD * T(i,k)))  
     &                 /(4./3*RPI*1000.*cdnc(i,k)) )**(1./3.)
                  rad_chaud_tab(i,k) = MAX(rad_chaud_tab(i,k) * 1e6, 5.) 
               ENDDO            
            ENDDO
         ELSE
            DO k = 1, MIN(3,klev)
               DO i = 1, klon
                  rad_chaud_tab(i,k) = rad_chau2
               ENDDO            
            ENDDO
            DO k = MIN(3,klev)+1, klev
               DO i = 1, klon
                  rad_chaud_tab(i,k) = rad_chau1
               ENDDO            
            ENDDO

         ENDIF
         
         DO k = 1, klev
!            IF(.not.ok_aie) THEN
            rad_chaud = rad_chau1
            IF (k.LE.3) rad_chaud = rad_chau2
!            ENDIF
            DO i = 1, klon
               IF (pclc(i,k) .LE. seuil_neb) THEN
               
c     -- effective cloud droplet radius (microns):
               
c     for liquid water clouds: 
                                ! For output diagnostics
                                !
                                ! Cloud droplet effective radius [um]
                                !
                                ! we multiply here with f * xl (fraction of liquid water
                                ! clouds in the grid cell) to avoid problems in the
                                ! averaging of the output.
                                ! In the output of IOIPSL, derive the real cloud droplet 
                                ! effective radius as re/fl
                                !
                                    
                  fl(i,k) = seuil_neb*(1.-zfice2(i,k))            
                  re(i,k) = rad_chaud_tab(i,k)*fl(i,k)
                  
                  rel = 0.
                  rei = 0.
                  pclc(i,k) = 0.0
                  pcltau(i,k) = 0.0
                  pclemi(i,k) = 0.0
                  cldtaupi(i,k) = 0.0                  
               ELSE

c     -- liquid/ice cloud water paths:
                  
                  zflwp_var= 1000.*(1.-zfice2(i,k))*pqlwp(i,k)/pclc(i,k)
     &                 *diff_paprs(i,k)
                  zfiwp_var= 1000.*zfice2(i,k)*pqlwp(i,k)/pclc(i,k)
     &                 *diff_paprs(i,k)
                  
c     -- effective cloud droplet radius (microns):
               
c     for liquid water clouds: 
                                    
                  IF (ok_aie) THEN
                     radius = 
     &                    1.1
     &                    *((pqlwp(i,k)*pplay(i,k)/(RD * T(i,k)))  
     &                    /(4./3.*RPI*1000.*cdnc_pi(i,k)))**(1./3.)
                     radius = MAX(radius*1e6, 5.) 
                  
                     tc = t(i,k)-273.15
                     rei = 0.71*tc + 61.29 
                     if (tc.le.-81.4) rei = 3.5 
                     if (zflwp_var.eq.0.) radius = 1. 
                     if (zfiwp_var.eq.0. .or. rei.le.0.) rei = 1. 
                     cldtaupi(i,k) = 3.0/2.0 * zflwp_var / radius
     &                    + zfiwp_var * (3.448e-03  + 2.431/rei)

                  ENDIF         ! ok_aie
                                ! For output diagnostics
                                !
                                ! Cloud droplet effective radius [um]
                                !
                                ! we multiply here with f * xl (fraction of liquid water
                                ! clouds in the grid cell) to avoid problems in the
                                ! averaging of the output.
                                ! In the output of IOIPSL, derive the real cloud droplet 
                                ! effective radius as re/fl
                                !
 
                  fl(i,k) = pclc(i,k)*(1.-zfice2(i,k))            
                  re(i,k) = rad_chaud_tab(i,k)*fl(i,k)
                  
                  rel = rad_chaud_tab(i,k)
c     for ice clouds: as a function of the ambiant temperature
c     [formula used by Iacobellis and Somerville (2000), with an 
c     asymptotical value of 3.5 microns at T<-81.4 C added to be 
c     consistent with observations of Heymsfield et al. 1986]:
                  tc = t(i,k)-273.15
                  rei = 0.71*tc + 61.29 
                  if (tc.le.-81.4) rei = 3.5 
c     -- cloud optical thickness :
               
c     [for liquid clouds, traditional formula, 
c     for ice clouds, Ebert & Curry (1992)] 
                  
                 if (zflwp_var.eq.0.) rel = 1. 
                 if (zfiwp_var.eq.0. .or. rei.le.0.) rei = 1. 
                 pcltau(i,k) = 3.0/2.0 * ( zflwp_var/rel )
     &                 + zfiwp_var * (3.448e-03  + 2.431/rei)
c     -- cloud infrared emissivity:
               
c     [the broadband infrared absorption coefficient is parameterized
c     as a function of the effective cld droplet radius]
               
c     Ebert and Curry (1992) formula as used by Kiehl & Zender (1995):
                  k_ice = k_ice0 + 1.0/rei
                  
                  pclemi(i,k) = 1.0
     &                 - EXP( -coef_chau*zflwp_var - DF*k_ice*zfiwp_var)

               ENDIF
              reliq(i,k)=rel
              reice(i,k)=rei 
!              if (i.eq.1) then
!              print*,'Dans newmicro rel, rei :',rel, rei
!              print*,'Dans newmicro reliq, reice :',
!     $             reliq(i,k),reice(i,k)
!              endif

            ENDDO
         ENDDO

         DO k = 1, klev
            DO i = 1, klon
               xflwp(i) = xflwp(i)+ xflwc(i,k) * diff_paprs(i,k)
               xfiwp(i) = xfiwp(i)+ xfiwc(i,k) * diff_paprs(i,k)
            ENDDO
         ENDDO

      ELSE
         DO k = 1, klev
            rad_chaud = rad_chau1
            IF (k.LE.3) rad_chaud = rad_chau2
            DO i = 1, klon
                              
               IF (pclc(i,k) .LE. seuil_neb) THEN

                  pclc(i,k) = 0.0
                  pcltau(i,k) = 0.0
                  pclemi(i,k) = 0.0
                  cldtaupi(i,k) = 0.0

               ELSE

                  zflwp_var = 1000.*pqlwp(i,k)*diff_paprs(i,k)
     &                 /pclc(i,k)
                  
                  zfice1 = MIN(
     &                 MAX( 1.0 - (t(i,k)-t_glace_min) / 
     &                    (t_glace_max-t_glace_min),0.0),1.0)**nexpo
                  
                  radius = rad_chaud * (1.-zfice1) + rad_froid * zfice1
                  coef   = coef_chau * (1.-zfice1) + coef_froi * zfice1

                  pcltau(i,k) = 3.0 * zflwp_var / (2.0 * radius)
                  pclemi(i,k) = 1.0 - EXP( - coef * zflwp_var)

               ENDIF
                              
            ENDDO
         ENDDO
      ENDIF
      
      IF (.NOT.ok_aie) THEN
         DO k = 1, klev
            DO i = 1, klon
               cldtaupi(i,k)=pcltau(i,k)
            ENDDO
         ENDDO               
      ENDIF

ccc   DO k = 1, klev
ccc   DO i = 1, klon
ccc   t(i,k) = t(i,k)
ccc   pclc(i,k) = MAX( 1.e-5 , pclc(i,k) )
ccc   lo = pclc(i,k) .GT. (2.*1.e-5)
ccc   zflwp = pqlwp(i,k)*1000.*(paprs(i,k)-paprs(i,k+1))
ccc   .          /(rg*pclc(i,k))
ccc   zradef = 10.0 + (1.-sigs(k))*45.0
ccc   pcltau(i,k) = 1.5 * zflwp / zradef
ccc   zfice=1.0-MIN(MAX((t(i,k)-263.)/(273.-263.),0.0),1.0)
ccc   zmsac = 0.13*(1.0-zfice) + 0.08*zfice
ccc   pclemi(i,k) = 1.-EXP(-zmsac*zflwp)
ccc   if (.NOT.lo) pclc(i,k) = 0.0
ccc   if (.NOT.lo) pcltau(i,k) = 0.0
ccc   if (.NOT.lo) pclemi(i,k) = 0.0
ccc   ENDDO
ccc   ENDDO
ccccc print*, 'pas de nuage dans le rayonnement'
ccccc DO k = 1, klev
ccccc DO i = 1, klon
ccccc pclc(i,k) = 0.0
ccccc pcltau(i,k) = 0.0
ccccc pclemi(i,k) = 0.0
ccccc ENDDO
ccccc ENDDO
C     
C     COMPUTE CLOUD LIQUID PATH AND TOTAL CLOUDINESS
C     
c     IM cf. CR:test: calcul prenant ou non en compte le recouvrement
c     initialisations
      DO i=1,klon
         zclear(i)=1.
         zcloud(i)=0.
         pch(i)=1.0
         pcm(i) = 1.0
         pcl(i) = 1.0
         pctlwp(i) = 0.0
      ENDDO
C
cIM cf CR DO k=1,klev
      DO k = klev, 1, -1
         DO i = 1, klon
            pctlwp(i) = pctlwp(i) 
     &           + pqlwp(i,k)*diff_paprs(i,k)
         ENDDO
      ENDDO
c     IM cf. CR
      IF (NOVLP.EQ.1) THEN
         DO k = klev, 1, -1
            DO i = 1, klon
               zclear(i)=zclear(i)*(1.-MAX(pclc(i,k),zcloud(i)))
     &              /(1.-MIN(real(zcloud(i), kind=8),1.-ZEPSEC))
               pct(i)=1.-zclear(i) 
               IF (pplay(i,k).LE.cetahb*paprs(i,1)) THEN
                  pch(i) = pch(i)*(1.-MAX(pclc(i,k),zcloud(i)))
     &                 /(1.-MIN(real(zcloud(i), kind=8),1.-ZEPSEC))
               ELSE IF (pplay(i,k).GT.cetahb*paprs(i,1) .AND.
     &                 pplay(i,k).LE.cetamb*paprs(i,1)) THEN
                  pcm(i) = pcm(i)*(1.-MAX(pclc(i,k),zcloud(i)))
     &                 /(1.-MIN(real(zcloud(i), kind=8),1.-ZEPSEC))
               ELSE IF (pplay(i,k).GT.cetamb*paprs(i,1)) THEN 
                  pcl(i) = pcl(i)*(1.-MAX(pclc(i,k),zcloud(i)))
     &                 /(1.-MIN(real(zcloud(i), kind=8),1.-ZEPSEC))
               endif
               zcloud(i)=pclc(i,k)
            ENDDO
         ENDDO
      ELSE IF (NOVLP.EQ.2) THEN
         DO k = klev, 1, -1
            DO i = 1, klon
               zcloud(i)=MAX(pclc(i,k),zcloud(i))
               pct(i)=zcloud(i)
               IF (pplay(i,k).LE.cetahb*paprs(i,1)) THEN
                  pch(i) = MIN(pclc(i,k),pch(i))
               ELSE IF (pplay(i,k).GT.cetahb*paprs(i,1) .AND.
     &                 pplay(i,k).LE.cetamb*paprs(i,1)) THEN
                  pcm(i) = MIN(pclc(i,k),pcm(i))
               ELSE IF (pplay(i,k).GT.cetamb*paprs(i,1)) THEN
                  pcl(i) = MIN(pclc(i,k),pcl(i))
               endif
            ENDDO
         ENDDO
      ELSE IF (NOVLP.EQ.3) THEN
         DO k = klev, 1, -1
            DO i = 1, klon
               zclear(i)=zclear(i)*(1.-pclc(i,k))
               pct(i)=1-zclear(i)
               IF (pplay(i,k).LE.cetahb*paprs(i,1)) THEN
                  pch(i) = pch(i)*(1.0-pclc(i,k))
               ELSE IF (pplay(i,k).GT.cetahb*paprs(i,1) .AND.
     &                 pplay(i,k).LE.cetamb*paprs(i,1)) THEN 
                  pcm(i) = pcm(i)*(1.0-pclc(i,k))
               ELSE IF (pplay(i,k).GT.cetamb*paprs(i,1)) THEN
                  pcl(i) = pcl(i)*(1.0-pclc(i,k))
               endif
            ENDDO
         ENDDO
      ENDIF
      
C     
      DO i = 1, klon
c     IM cf. CR          pct(i)=1.-pct(i)
         pch(i)=1.-pch(i)
         pcm(i)=1.-pcm(i)
         pcl(i)=1.-pcl(i)
      ENDDO

c ========================================================
! DIAGNOSTICS CALCULATION FOR CMIP5 PROTOCOL
c ========================================================
!! change by Nicolas Yan (LSCE) 
!! Cloud Droplet Number Concentration (CDNC) : 3D variable
!! Fractionnal cover by liquid water cloud (LCC3D) : 3D variable
!! Cloud Droplet Number Concentration at top of cloud (CLDNCL) : 2D variable
!! Droplet effective radius at top of cloud (REFFCLWTOP) : 2D variable
!! Fractionnal cover by liquid water at top of clouds (LCC) : 2D variable
      IF (ok_newmicro) THEN
         IF (ok_aie) THEN
            DO k = 1, klev
               DO i = 1, klon
                  phase3d(i,k)=1-zfice2(i,k)
                  IF (pclc(i,k) .LE. seuil_neb) THEN
                     lcc3d(i,k)=seuil_neb*phase3d(i,k)
                  ELSE
                     lcc3d(i,k)=pclc(i,k)*phase3d(i,k)
                  ENDIF
                  scdnc(i,k)=lcc3d(i,k)*cdnc(i,k) ! m-3
               ENDDO
            ENDDO

            DO i=1,klon
               lcc(i)=0.
               reffclwtop(i)=0.
               cldncl(i)=0.
               IF(RANDOM .OR. MAXIMUM_RANDOM) tcc(i) = 1.
               IF(MAXIMUM) tcc(i) = 0.
            ENDDO
     
            FIRST=.TRUE.

            DO i=1,klon
               DO k=klev-1,1,-1 !From TOA down


            ! Test, if the cloud optical depth exceeds the necessary
            ! threshold:

             IF (pcltau(i,k).GT.thres_tau .AND. pclc(i,k).GT.thres_neb)
     .                                                             THEN
               ! To calculate the right Temperature at cloud top,
               ! interpolate it between layers:
                  t_tmp = t(i,k) +
     .              (paprs(i,k+1)-pplay(i,k))/(pplay(i,k+1)-pplay(i,k))
     .              * ( t(i,k+1) - t(i,k) )

                  IF(MAXIMUM) THEN
                    IF(FIRST) THEN
                       write(*,*)'Hypothese de recouvrement: MAXIMUM'
                       FIRST=.FALSE.
                    ENDIF
                    flag_max= -1.
                    ftmp(i) = MAX(tcc(i),pclc(i,k))
                  ENDIF

                  IF(RANDOM) THEN
                    IF(FIRST) THEN
                       write(*,*)'Hypothese de recouvrement: RANDOM'
                       FIRST=.FALSE.
                    ENDIF
                    flag_max= 1.
                    ftmp(i) = tcc(i) * (1-pclc(i,k))
                  ENDIF

                  IF(MAXIMUM_RANDOM) THEN
                    IF(FIRST) THEN
                       write(*,*)'Hypothese de recouvrement: MAXIMUM_
     .                         RANDOM'
                       FIRST=.FALSE.
                    ENDIF
                    flag_max= 1.
                    ftmp(i) = tcc(i) *
     .              (1. - MAX(pclc(i,k),pclc(i,k+1))) /
     .              (1. - MIN(pclc(i,k+1),1.-thres_neb))
                  ENDIF
c Effective radius of cloud droplet at top of cloud (m)
                  reffclwtop(i) = reffclwtop(i) + rad_chaud_tab(i,k) * 
     .           1.0E-06 * phase3d(i,k) * ( tcc(i) - ftmp(i))*flag_max
c CDNC at top of cloud (m-3)
                  cldncl(i) = cldncl(i) + cdnc(i,k) * phase3d(i,k) * 
     .                 (tcc(i) - ftmp(i))*flag_max
c Liquid Cloud Content at top of cloud
                  lcc(i) = lcc(i) + phase3d(i,k) * (tcc(i)-ftmp(i))*
     .                    flag_max
c Total Cloud Content at top of cloud
                  tcc(i)=ftmp(i)
              
          ENDIF ! is there a visible, not-too-small cloud?  
          ENDDO ! loop over k

          IF(RANDOM .OR. MAXIMUM_RANDOM) tcc(i)=1.-tcc(i)
         ENDDO ! loop over i

!! Convective and Stratiform Cloud Droplet Effective Radius (REFFCLWC  REFFCLWS)
            DO i = 1, klon
               DO k = 1, klev
                  pqlwpcon(i,k)=rnebcon(i,k)*clwcon(i,k) ! fraction eau liquide convective
                  pqlwpstra(i,k)=pclc(i,k)*phase3d(i,k)-pqlwpcon(i,k) ! fraction eau liquide stratiforme
                  IF (pqlwpstra(i,k) .LE. 0.0) pqlwpstra(i,k)=0.0
! Convective Cloud Droplet Effective Radius (REFFCLWC) : variable 3D
                  reffclwc(i,k)=1.1
     &                 *((pqlwpcon(i,k)*pplay(i,k)/(RD * T(i,k)))
     &                 /(4./3*RPI*1000.*cdnc(i,k)) )**(1./3.)
                  reffclwc(i,k) = MAX(reffclwc(i,k) * 1e6, 5.)

! Stratiform Cloud Droplet Effective Radius (REFFCLWS) : variable 3D
                  IF ((pclc(i,k)-rnebcon(i,k)) .LE. seuil_neb) THEN ! tout sous la forme convective
                     reffclws(i,k)=0.0
                     lcc3dstra(i,k)= 0.0
                  ELSE
                     reffclws(i,k) = (pclc(i,k)*phase3d(i,k)*
     &                               rad_chaud_tab(i,k)-
     &                            pqlwpcon(i,k)*reffclwc(i,k))
                     IF(reffclws(i,k) .LE. 0.0) reffclws(i,k)=0.0
                     lcc3dstra(i,k)=pqlwpstra(i,k)
                 ENDIF
!Convertion from um to m
                  IF(rnebcon(i,k). LE. seuil_neb) THEN
                    reffclwc(i,k) = reffclwc(i,k)*seuil_neb*clwcon(i,k)
     &                              *1.0E-06
                    lcc3dcon(i,k)= seuil_neb*clwcon(i,k)
                  ELSE
                    reffclwc(i,k) = reffclwc(i,k)*pqlwpcon(i,k)
     &                              *1.0E-06
                    lcc3dcon(i,k) = pqlwpcon(i,k)
                  ENDIF

                  reffclws(i,k) = reffclws(i,k)*1.0E-06

               ENDDO !klev
            ENDDO !klon 

!! Column Integrated Cloud Droplet Number (CLDNVI) : variable 2D
            DO k = 1, klev
               DO i = 1, klon
                   pdel(i,k) = paprs(i,k)-paprs(i,k+1)
                   zrho(i,k)=pplay(i,k)/t(i,k)/RD                  ! kg/m3
                   dh(i,k)=pdel(i,k)/(gravit*zrho(i,k)) ! hauteur de chaque boite (m)
               ENDDO
            ENDDO
c
            DO i = 1, klon
               cldnvi(i)=0.
               lcc_integrat(i)=0.
               height(i)=0.
               DO k = 1, klev
                  cldnvi(i)=cldnvi(i)+cdnc(i,k)*lcc3d(i,k)*dh(i,k)
                  lcc_integrat(i)=lcc_integrat(i)+lcc3d(i,k)*dh(i,k)
                  height(i)=height(i)+dh(i,k)
               ENDDO ! klev
               lcc_integrat(i)=lcc_integrat(i)/height(i)
               IF (lcc_integrat(i) .LE. 1.0E-03) THEN
                  cldnvi(i)=cldnvi(i)*lcc(i)/seuil_neb
               ELSE
                  cldnvi(i)=cldnvi(i)*lcc(i)/lcc_integrat(i)
               ENDIF
            ENDDO ! klon
            
            DO i = 1, klon
               DO k = 1, klev
                IF (scdnc(i,k) .LE. 0.0) scdnc(i,k)=0.0
                IF (reffclws(i,k) .LE. 0.0) reffclws(i,k)=0.0
                IF (reffclwc(i,k) .LE. 0.0) reffclwc(i,k)=0.0
                IF (lcc3d(i,k) .LE. 0.0) lcc3d(i,k)=0.0
                IF (lcc3dcon(i,k) .LE. 0.0) lcc3dcon(i,k)=0.0
                IF (lcc3dstra(i,k) .LE. 0.0) lcc3dstra(i,k)=0.0
               ENDDO
               IF (reffclwtop(i) .LE. 0.0) reffclwtop(i)=0.0
               IF (cldncl(i) .LE. 0.0) cldncl(i)=0.0
               IF (cldnvi(i) .LE. 0.0) cldnvi(i)=0.0
               IF (lcc(i) .LE. 0.0) lcc(i)=0.0
            ENDDO

         ENDIF !ok_aie
      ENDIF !ok newmicro
c
C
      RETURN
      END
