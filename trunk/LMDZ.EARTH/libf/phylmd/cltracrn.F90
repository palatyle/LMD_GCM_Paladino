!$Id $

SUBROUTINE cltracrn( itr, dtime,u1lay, v1lay, &
     cdrag,coef,t,ftsol,pctsrf,               &
     tr,trs,paprs,pplay,delp,                 &
     masktr,fshtr,hsoltr,tautr,vdeptr,        &
     lat,d_tr,d_trs )
  
  USE dimphy
  USE traclmdz_mod, ONLY : id_rn, id_pb
  IMPLICIT NONE
!======================================================================
! Auteur(s): Alex/LMD) date:  fev 99
!            inspire de clqh + clvent
! Objet: diffusion verticale de traceurs avec quantite de traceur ds 
!        le sol ( reservoir de sol de radon ) 
!        
! note : pour l'instant le traceur dans le sol et le flux sont
!        calcules mais ils ne servent que de diagnostiques
!        seule la tendance sur le traceur est sortie (d_tr)
!---------------------------------------------------------------------
! Arguments:
! itr......input-R-  le type de traceur : id_rn(radon), id_pb(plomb)
! dtime....input-R-  intervalle du temps (en secondes) ~ pdtphys
! u1lay....input-R-  vent u de la premiere couche (m/s)
! v1lay....input-R-  vent v de la premiere couche (m/s)
! cdrag....input-R-  cdrag
! coef.....input-R-  le coefficient d'echange (m**2/s) l>1, valable uniquement pour k entre 2 et klev
! t........input-R-  temperature (K)
! paprs....input-R-  pression a inter-couche (Pa)
! pplay....input-R-  pression au milieu de couche (Pa)
! delp.....input-R-  epaisseur de couche (Pa)
! ftsol....input-R-  temperature du sol (en Kelvin)
! tr.......input-R-  traceurs
! trs......input-R-  traceurs dans le sol
! masktr...input-R-  Masque reservoir de sol traceur (1 = reservoir)
! fshtr....input-R-  Flux surfacique de production dans le sol
! tautr....input-R-  Constante de decroissance du traceur
! vdeptr...input-R-  Vitesse de depot sec dans la couche brownienne
! hsoltr...input-R-  Epaisseur equivalente du reservoir de sol
! lat......input-R-  latitude en degree
! d_tr.....output-R- le changement de "tr"
! d_trs....output-R- le changement de "trs"
!======================================================================
  include "YOMCST.h"
  include "indicesol.h"
!
!Entrees
  INTEGER,INTENT(IN)                     :: itr
  REAL,INTENT(IN)                        :: dtime
  REAL,DIMENSION(klon),INTENT(IN)        :: u1lay, v1lay
  REAL,DIMENSION(klon),INTENT(IN)        :: cdrag
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: coef, t
  REAL,DIMENSION(klon,nbsrf),INTENT(IN)  :: ftsol, pctsrf 
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: tr 
  REAL,DIMENSION(klon),INTENT(IN)        :: trs
  REAL,DIMENSION(klon,klev+1),INTENT(IN) :: paprs 
  REAL,DIMENSION(klon,klev),INTENT(IN)   :: pplay, delp
  REAL,DIMENSION(klon),INTENT(IN)        :: masktr 
  REAL,DIMENSION(klon),INTENT(IN)        :: fshtr 
  REAL,INTENT(IN)                        :: hsoltr
  REAL,INTENT(IN)                        :: tautr
  REAL,INTENT(IN)                        :: vdeptr
  REAL,DIMENSION(klon),INTENT(IN)        :: lat   

!Sorties
  REAL,DIMENSION(klon,klev),INTENT(OUT) :: d_tr
  REAL,DIMENSION(klon),INTENT(OUT) :: d_trs  ! (diagnostic) traceur ds le sol

!Locales
  REAL,DIMENSION(klon,klev) :: flux_tr  ! (diagnostic) flux de traceur
  INTEGER                   :: i, k, n, l
  REAL,DIMENSION(klon)      :: rotrhi
  REAL,DIMENSION(klon,klev) :: zx_coef
  REAL,DIMENSION(klon)      :: zx_buf
  REAL,DIMENSION(klon,klev) :: zx_ctr
  REAL,DIMENSION(klon,klev) :: zx_dtr
  REAL,DIMENSION(klon)      :: zx_trs
  REAL                      :: zx_a, zx_b
  
  REAL,DIMENSION(klon,klev) :: local_tr
  REAL,DIMENSION(klon)      :: local_trs
  REAL,DIMENSION(klon)      :: zts      ! champ de temperature du sol
  REAL,DIMENSION(klon)      :: zx_alpha1, zx_alpha2
!======================================================================
!AA Pour l'instant les 4 types de surface ne sont pas pris en compte
!AA On fabrique avec zts un champ de temperature de sol  
!AA que le pondere par la fraction de nature de sol.
 
  DO i = 1,klon
     zts(i) = 0. 
  ENDDO

  DO n=1,nbsrf
     DO i = 1,klon
        zts(i) = zts(i) + ftsol(i,n)*pctsrf(i,n)
     ENDDO
  ENDDO

  DO i = 1,klon
     rotrhi(i) = RD * zts(i) / hsoltr 
  ENDDO

  DO k = 1, klev
     DO i = 1, klon
        local_tr(i,k) = tr(i,k)
     ENDDO
  ENDDO

  DO i = 1, klon
     local_trs(i) = trs(i)
  ENDDO
!======================================================================
!AA   Attention si dans clmain zx_alf1(i) = 1.0 
!AA   Il doit y avoir coherence (dc la meme chose ici)

  DO i = 1, klon
!AA         zx_alpha1(i) = (paprs(i,1)-pplay(i,2))/(pplay(i,1)-pplay(i,2))
     zx_alpha1(i) = 1.0
     zx_alpha2(i) = 1.0 - zx_alpha1(i)
  ENDDO
!======================================================================
  DO i = 1, klon
     zx_coef(i,1) = cdrag(i)*(1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) &
          *pplay(i,1)/(RD*t(i,1))
     zx_coef(i,1) = zx_coef(i,1) * dtime*RG
  ENDDO

  DO k = 2, klev
     DO i = 1, klon
        zx_coef(i,k) = coef(i,k)*RG/(pplay(i,k-1)-pplay(i,k)) &
             *(paprs(i,k)*2/(t(i,k)+t(i,k-1))/RD)**2
        zx_coef(i,k) = zx_coef(i,k) * dtime*RG
     ENDDO
  ENDDO
!======================================================================
  DO i = 1, klon
     zx_buf(i)      = delp(i,klev) + zx_coef(i,klev)
     zx_ctr(i,klev) = local_tr(i,klev)*delp(i,klev)/zx_buf(i)
     zx_dtr(i,klev) = zx_coef(i,klev) / zx_buf(i)
  ENDDO

  DO l = klev-1, 2 , -1
     DO i = 1, klon
        zx_buf(i) = delp(i,l)+zx_coef(i,l)      &
             +zx_coef(i,l+1)*(1.-zx_dtr(i,l+1))
 
        zx_ctr(i,l) = ( local_tr(i,l)*delp(i,l) &
             + zx_coef(i,l+1)*zx_ctr(i,l+1) )/zx_buf(i)
        zx_dtr(i,l) = zx_coef(i,l) / zx_buf(i)
     ENDDO
  ENDDO

  DO i = 1, klon
     zx_buf(i) = delp(i,1) + zx_coef(i,2)*(1.-zx_dtr(i,2))  &
          + masktr(i) * zx_coef(i,1)                        &
          *( zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i,2) )

     zx_ctr(i,1) = ( local_tr(i,1)*delp(i,1)                &
          + zx_ctr(i,2)                                     &
          *(zx_coef(i,2)                                    &
          - masktr(i) * zx_coef(i,1)                        &
          *zx_alpha2(i) ) ) / zx_buf(i)
     zx_dtr(i,1) = masktr(i) * zx_coef(i,1) / zx_buf(i)
  ENDDO
!======================================================================
! Calculer d'abord local_trs nouvelle quantite dans le reservoir
! de sol
!=====================================================================

  DO i = 1, klon
!-------------------------
! Au dessus des continents
!--
! Le pb peut se deposer partout : vdeptr = 10-3 m/s
! Le Rn est traiter commme une couche Brownienne puisque vdeptr = 0.
!-------------------------------------------------------------------
     IF ( NINT(masktr(i)) .EQ. 1 ) THEN
        zx_trs(i) = local_trs(i)
        zx_a = zx_trs(i)                                           &
             +fshtr(i)*dtime*rotrhi(i)                             &
             +rotrhi(i)*masktr(i)*zx_coef(i,1)/RG                  &
             *(zx_ctr(i,1)*(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i,2)) &
             +zx_alpha2(i)*zx_ctr(i,2))
! Pour l'instant, pour aller vite, le depot sec est traite comme une decroissance
        zx_b = 1. + rotrhi(i)*masktr(i)*zx_coef(i,1)/RG            &
             * (1.-zx_dtr(i,1)                                     &
             *(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i,2)))             &
             + dtime / tautr                                       & 
             + dtime * vdeptr / hsoltr 
        zx_trs(i) = zx_a / zx_b
        local_trs(i) = zx_trs(i)
     ENDIF
!--------------------------------------------------------
! Si on est entre 60N et 70N on divise par 2 l'emanation
!--------------------------------------------------------

     IF ( (itr.eq.id_rn.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GE.60..AND.lat(i).LE.70.).OR.      &
          (itr.eq.id_pb.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GE.60..AND.lat(i).LE.70.) ) THEN
        zx_trs(i) = local_trs(i)
        zx_a = zx_trs(i)                                           &
             +(fshtr(i)/2.)*dtime*rotrhi(i)                        & 
             +rotrhi(i)*masktr(i)*zx_coef(i,1)/RG                  &
             *(zx_ctr(i,1)*(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i,2)) &
             +zx_alpha2(i)*zx_ctr(i,2))
        !
        zx_b = 1. + rotrhi(i)*masktr(i)*zx_coef(i,1)/RG  &
             * (1.-zx_dtr(i,1)                           &
             *(zx_alpha1(i)+zx_alpha2(i)*zx_dtr(i,2)))   &
             + dtime / tautr                             &
             + dtime * vdeptr / hsoltr
        ! 
        zx_trs(i) = zx_a / zx_b
        local_trs(i) = zx_trs(i)
     ENDIF

!----------------------------------------------
! Au dessus des oceans et aux hautes latitudes
!--
! au dessous de -60S  pas d'emission de radon au dessus 
! des oceans et des continents
!---------------------------------------------------------------

     IF ( (itr.EQ.id_rn.AND.NINT(masktr(i)).EQ.0).OR.       &
          (itr.EQ.id_rn.AND.NINT(masktr(i)).EQ.1.AND.lat(i).LT.-60.)) THEN
        zx_trs(i) = 0.
        local_trs(i) = 0.
     END IF
!--
! au dessus de 70 N pas d'emission de radon au dessus 
! des oceans et des continents
!--------------------------------------------------------------
     IF ( (itr.EQ.id_rn.AND.NINT(masktr(i)).EQ.0).OR.    &
          (itr.EQ.id_rn.AND.NINT(masktr(i)).EQ.1.AND.lat(i).GT.70.)) THEN
        zx_trs(i) = 0.
        local_trs(i) = 0.
     END IF
!---------------------------------------------
! Au dessus des oceans la source est nulle
!--------------------------------------------

     IF (itr.eq.id_rn.AND.NINT(masktr(i)).EQ.0) THEN
        zx_trs(i) = 0.
        local_trs(i) = 0.
     END IF

  ENDDO    ! sur le i=1,klon
!
!======================================================================
! Une fois on a zx_trs, on peut faire l'iteration 
!====================================================================== 

  DO i = 1, klon
     local_tr(i,1) = zx_ctr(i,1)+zx_dtr(i,1)*zx_trs(i)
  ENDDO
  DO l = 2, klev
     DO i = 1, klon
        local_tr(i,l) = zx_ctr(i,l) + zx_dtr(i,l)*local_tr(i,l-1)
     ENDDO
  ENDDO
!======================================================================
! Calcul du flux de traceur (flux_tr): UA/(m**2 s)
!======================================================================
  DO i = 1, klon
     flux_tr(i,1) = masktr(i)*zx_coef(i,1)/RG                      &
          * (zx_alpha1(i)*local_tr(i,1)+zx_alpha2(i)*local_tr(i,2) &
          -zx_trs(i)) / dtime
  ENDDO
  DO l = 2, klev
     DO i = 1, klon
        flux_tr(i,l) = zx_coef(i,l)/RG                    &
             * (local_tr(i,l)-local_tr(i,l-1)) / dtime
     ENDDO
  ENDDO
!======================================================================
! Calcul des tendances du traceur ds le sol et dans l'atmosphere
!======================================================================
  DO l = 1, klev
     DO i = 1, klon
        d_tr(i,l) = local_tr(i,l) - tr(i,l)
     ENDDO
  ENDDO
  DO i = 1, klon
     d_trs(i) = local_trs(i) - trs(i)
  ENDDO

END SUBROUTINE cltracrn
