!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/interface_surf.F90,v 1.6 2005/02/24 09:58:18 fairhead Exp $
!

  MODULE interface_surf

! Ce module regroupe toutes les routines gerant l'interface entre le modele 
! atmospherique et les modeles de surface (sols continentaux, oceans, glaces)
! Les routines sont les suivantes:
!
!   interfsurf_*: routines d'aiguillage vers les interfaces avec les 
!                 differents modeles de surface
!
! L. Fairhead, LMD, 02/2000

  USE ioipsl

  IMPLICIT none

  PRIVATE
  PUBLIC :: interfsurf,interfsurf_hq 

  INTERFACE interfsurf
    module procedure interfsurf_hq
  END INTERFACE

#include "YOMCST.h"

  CONTAINS
!
!############################################################################
!
! ADAPTATION GCM POUR CP(T)
  SUBROUTINE interfsurf_hq(itime, dtime, rmu0, &
      & klon, iim, jjm, knon, &
      & rlon, rlat, cufi, cvfi, &
      & debut, lafin, soil_model, nsoilmx, tsoil, &
      & zlev,  u1_lay, v1_lay, temp_air, epot_air, & 
      & tq_cdrag, petAcoef, petBcoef, &
      & sollw, sollwdown, swnet, swdown, &
      & fder, taux, tauy, &
      & albedo, &
      & tsurf, pkh1, p1lay, radsol, &
      & fluxsens, dflux_s, &              
      & tsol_rad, tsurf_new, alb_new)

      use write_field_phy
      use cpdet_phy_mod, only: cpdet

      IMPLICIT none

! Cette routine sert d'aiguillage entre l'atmosphere et la surface en general 
! (sols continentaux, oceans, glaces) pour les fluxs de chaleur et d'humidite.
! En pratique l'interface se fait entre la couche limite du modele 
! atmospherique (clmain.F) et les routines de surface (sechiba, oasis, ...)
!
! 
! L.Fairhead 02/2000
!
! input:
!   itime        numero du pas de temps
!   klon         nombre total de points de grille
!   iim, jjm     nbres de pts de grille
!   dtime        pas de temps de la physique (en s)
!   rmu0         cosinus de l'angle solaire zenithal
!   knon         nombre de points de la surface a traiter
!   rlon         longitudes
!   rlat         latitudes
!   cufi,cvfi    resolution des mailles en x et y (m)
!   debut        logical: 1er appel a la physique
!   lafin        logical: dernier appel a la physique
!   zlev         hauteur de la premiere couche
!   u1_lay       vitesse u 1ere couche
!   v1_lay       vitesse v 1ere couche
!   temp_air     temperature de l'air 1ere couche
!   epot_air     temp potentielle de l'air
!   tq_cdrag     cdrag
!   petAcoef     coeff. A de la resolution de la CL pour t
!   petBcoef     coeff. B de la resolution de la CL pour t
!   sollw        flux IR net a la surface
!   sollwdown    flux IR descendant a la surface
!   swnet        flux solaire net
!   swdown       flux solaire entrant a la surface
!   albedo       albedo de la surface
!   tsurf        temperature de surface
!   pkh1         fct Exner à la surface: RCPD*(paprs(1)/preff)**RKAPPA
!   p1lay        pression 1er niveau (milieu de couche)
!   radsol       rayonnement net aus sol (LW + SW)
!   fder         derivee des flux (pour le couplage)
!   taux, tauy   tension de vents
!
! output:
!   fluxsens     flux de chaleur sensible
!   tsol_rad     
!   tsurf_new    temperature au sol
!   alb_new      albedo

#include "iniprint.h"


! Parametres d'entree
  integer, intent(IN) :: itime
  integer, intent(IN) :: iim, jjm
  integer, intent(IN) :: klon
  real, intent(IN) :: dtime
  real, intent(IN)    :: rmu0(klon)
  integer, intent(IN) :: knon
  logical, intent(IN) :: debut, lafin
  real, dimension(klon), intent(IN) :: rlon, rlat
  real, dimension(klon), intent(IN) :: cufi, cvfi
  real, dimension(klon), intent(INOUT) :: tq_cdrag
  real, dimension(klon), intent(IN) :: zlev
  real, dimension(klon), intent(IN) :: u1_lay, v1_lay
  real, dimension(klon), intent(IN) :: temp_air
  real, dimension(klon), intent(IN) :: epot_air
  real, dimension(klon), intent(IN) :: petAcoef
  real, dimension(klon), intent(IN) :: petBcoef
  real, dimension(klon), intent(IN) :: sollw, sollwdown, swnet, swdown
  real, dimension(klon), intent(IN) :: albedo
  real, dimension(klon), intent(IN) :: tsurf, pkh1, p1lay
  REAL, DIMENSION(klon), INTENT(INOUT) :: radsol,fder
  real, dimension(klon), intent(IN) :: taux, tauy
!! PB ajout pour soil
  logical          :: soil_model
  integer          :: nsoilmx
  REAL, DIMENSION(klon, nsoilmx) :: tsoil
  REAL, dimension(klon)          :: soilcap
  REAL, dimension(klon)          :: soilflux
! Parametres de sortie
  real, dimension(klon), intent(OUT):: fluxsens
  real, dimension(klon), intent(OUT):: tsol_rad, tsurf_new, alb_new
  real, dimension(klon), intent(OUT):: dflux_s

! Local
  character (len = 20),save :: modname = 'interfsurf_hq'
  character (len = 80) :: abort_message 
  integer, save        :: error
  integer              :: ii, index
  logical,save              :: check = .false.
  real, dimension(klon):: cal, beta, capsol
  real, dimension(klon):: tsurf_temp, zcp
  INTEGER,dimension(1) :: iloc
  INTEGER                 :: isize
  real, dimension(klon):: fder_prev

  if (check) write(*,*) 'Entree ', modname

! Initialisations diverses
!
  cal = 999999. ; beta = 999999. ; capsol = 999999.
  alb_new = albedo 
  tsurf_new = 999999.

! ADAPTATION GCM POUR CP(T)
       do ii=1,klon
         zcp(ii)=cpdet(tsurf(ii))
       enddo

       IF (soil_model) THEN 
           CALL soil(dtime, knon, tsurf, tsoil,soilcap, soilflux)
           cal(1:knon) = zcp(1:knon) / soilcap(1:knon)
! for tests:
!  call writefield_phy('interfsurf_hq_zcp',zcp,1)
!  call writefield_phy('interfsurf_hq_cal',cal,1)
!  call writefield_phy('interfsurf_hq_soilcap',soilcap,1)
!       print*,"DIAGNOSTIC SOIL"
!       print*,"soilcap=",soilcap
!       print*,"soilflux=",soilflux
!       print*,"radsol=",radsol(knon/2)
           radsol(1:knon) = radsol(1:knon)  + soilflux(1:knon)
       ELSE 
!           abort_message = "PAS DE MODELE DE SOL: CALCUL SOILCAP!!"
!           call abort_gcm(modname,abort_message,1)
! VENUS: Valeur pour inertie = 200:
           soilcap = 14735.
           print*,"PAS DE MODELE DE SOL, soilcap=",soilcap
           cal(1:knon) = zcp(1:knon) / soilcap(1:knon)
       ENDIF
! ADAPTATION GCM POUR CP(T)
       CALL calcul_fluxs( klon, knon, dtime, &
     &   tsurf, zcp, pkh1, p1lay, cal, beta, tq_cdrag, &
     &   radsol, temp_air, u1_lay, v1_lay, &
     &   petAcoef, petBcoef, &
     &   tsurf_new, fluxsens, dflux_s )

  END SUBROUTINE interfsurf_hq

!
!#########################################################################
!
  SUBROUTINE calcul_fluxs( klon, knon, dtime, &
! ADAPTATION GCM POUR CP(T)
     & tsurf, zcp, pkh1, p1lay, cal, beta, coef1lay, &
     & radsol, t1lay, u1lay, v1lay, &
     & petAcoef, petBcoef, &
     & tsurf_new, fluxsens, dflux_s)

  use write_field_phy
  use cpdet_phy_mod, only: t2tpot, tpot2t

  IMPLICIT none

! Cette routine calcule les fluxs en h a l'interface et eventuellement
! une temperature de surface (au cas ou ok_veget = false)
!
! L. Fairhead 4/2000
!
! input:
!   knon         nombre de points a traiter
!   tsurf        temperature de surface
!   zcp          Cp(Tsurf)              
!   pkh1         fct Exner à la surface: RCPD*(paprs(1)/preff)**RKAPPA
!   p1lay        pression 1er niveau (milieu de couche)
!   cal          capacite calorifique du sol
!   beta         evap reelle
!   coef1lay     coefficient d'echange
!   petAcoef     coeff. A de la resolution de la CL pour t
!   petBcoef     coeff. B de la resolution de la CL pour t
!   radsol       rayonnement net aus sol (LW + SW)
!
! output:
!   tsurf_new    temperature au sol
!   fluxsens     flux de chaleur sensible
!   dflux_s      derivee du flux de chaleur sensible / Ts
!

! Parametres d'entree
  integer, intent(IN) :: knon, klon
  real   , intent(IN) :: dtime
  real, dimension(klon), intent(IN) :: petAcoef
  real, dimension(klon), intent(IN) :: petBcoef
! ADAPTATION GCM POUR CP(T)
  real, dimension(klon), intent(IN) :: tsurf,pkh1,zcp
  real, dimension(klon), intent(IN) :: p1lay, cal, beta, coef1lay
  real, dimension(klon), intent(IN) :: radsol
  real, dimension(klon), intent(IN) :: t1lay, u1lay, v1lay

! Parametres sorties
  real, dimension(klon), intent(OUT):: tsurf_new, fluxsens
  real, dimension(klon), intent(OUT):: dflux_s

! Variables locales
  integer :: i
  real, dimension(klon) :: zx_mh, zx_nh, zx_oh
  real, dimension(klon) :: zx_coef
  real, dimension(klon) :: ztetasurf,ztetasurf_new
  real, dimension(klon) :: zx_k1
  real, dimension(klon) :: zx_q_0 , d_ts
  real                  :: zdelta, zcvm5, zcor
!
  logical, save         :: check = .false.
  character (len = 20)  :: modname = 'calcul_fluxs'
  character (len = 80) :: abort_message 
  logical,save         :: first = .true.,second=.false.

  if (check) write(*,*)'Entree ', modname

  IF (check) THEN
      WRITE(*,*)' radsol (min, max)' &
         &     , MINVAL(radsol(1:knon)), MAXVAL(radsol(1:knon))
      CALL flush(6)
  ENDIF

! 
! Initialisation
!
  fluxsens=0.
  dflux_s = 0.
!
  DO i = 1, knon

    zx_coef(i) = coef1lay(i) &
     & * SQRT(u1lay(i)**2+v1lay(i)**2) &
     & * p1lay(i)/(RD*t1lay(i))

  ENDDO


! === Calcul de la temperature de surface ===
! 
! MODIF VENUS:
! Le calcul se fait en temperature potentielle

  call t2tpot(knon,tsurf,ztetasurf,pkh1)

  do i = 1, knon
    zx_k1(i) = zx_coef(i)
  enddo


  do i = 1, knon

! H
    zx_oh(i) = 1. - (zx_k1(i) * petBcoef(i) * dtime)
    zx_mh(i) = zx_k1(i) * petAcoef(i) / zx_oh(i)
! Derives des flux dF/d(teta)s:
    zx_nh(i) = - (zx_k1(i) * zcp(i))/ zx_oh(i)
! Derives des flux dF/dTs (W m-2 K-1):      version Terre
!   zx_nh(i) = - (zx_k1(i) * RCPD * zx_pkh(i))/ zx_oh(i)

! Tsurface  Version Terre
!
!   tsurf_new(i) = (tsurf(i) + cal(i)/(RCPD * zx_pkh(i)) * dtime * &
!    &             (radsol(i) + zx_mh(i)) & 
!    &                 + dif_grnd(i) * t_grnd * dtime)/ &
!    &          ( 1. - dtime * cal(i)/(RCPD * zx_pkh(i)) * &
!    &                       zx_nh(i) &  
!    &                     + dtime * dif_grnd(i))
!
!   d_ts(i) = tsurf_new(i) - tsurf(i)
!   fluxsens(i) = zx_mh(i) + zx_nh(i) * tsurf_new(i)
! Derives des flux dF/dTs (W m-2 K-1):
!   dflux_s(i) = zx_nh(i)

! MODIF VENUS  : on vire dif_grnd (=0) et t_grnd
!                et on travaille en teta

    ztetasurf_new(i) = (ztetasurf(i) + cal(i)/zcp(i) * dtime * &
     &                  (radsol(i) + zx_mh(i)) & 
     &             ) / &
     &             ( 1.      - cal(i)/zcp(i) * dtime * &
     &                      zx_nh(i) )
  ENDDO

    call tpot2t(knon,ztetasurf_new,tsurf_new,pkh1)

  do i = 1, knon
    d_ts(i) = tsurf_new(i) - tsurf(i)
    fluxsens(i) = zx_mh(i) + zx_nh(i) * ztetasurf_new(i)
! Derives des flux dF/dTs (W m-2 K-1):
    dflux_s(i) = zx_nh(i)*ztetasurf(i)/tsurf(i)
  ENDDO

! for tests: write output fields...
!  call writefield_phy('calcul_fluxs_d_ts',d_ts,1)
!  call writefield_phy('calcul_fluxs_fluxsens',fluxsens,1)
!  call writefield_phy('calcul_fluxs_dflux_s',dflux_s,1)

  END SUBROUTINE calcul_fluxs
!
!#########################################################################
!
  END MODULE interface_surf
