!
MODULE calcul_fluxs_mod


CONTAINS
  SUBROUTINE calcul_fluxs( knon, nisurf, dtime, &
       tsurf, p1lay, cal, beta, coef1lay, ps, &
       precip_rain, precip_snow, snow, qsurf, &
       radsol, dif_grnd, t1lay, q1lay, u1lay, v1lay, &
       petAcoef, peqAcoef, petBcoef, peqBcoef, &
       tsurf_new, evap, fluxlat, fluxsens, dflux_s, dflux_l)
    
    USE dimphy, ONLY : klon

! Cette routine calcule les fluxs en h et q a l'interface et eventuellement
! une temperature de surface (au cas ou ok_veget = false)
!
! L. Fairhead 4/2000
!
! input:
!   knon         nombre de points a traiter
!   nisurf       surface a traiter
!   tsurf        temperature de surface
!   p1lay        pression 1er niveau (milieu de couche)
!   cal          capacite calorifique du sol
!   beta         evap reelle
!   coef1lay     coefficient d'echange
!   ps           pression au sol
!   precip_rain  precipitations liquides
!   precip_snow  precipitations solides
!   snow         champs hauteur de neige
!   runoff       runoff en cas de trop plein
!   petAcoef     coeff. A de la resolution de la CL pour t
!   peqAcoef     coeff. A de la resolution de la CL pour q
!   petBcoef     coeff. B de la resolution de la CL pour t
!   peqBcoef     coeff. B de la resolution de la CL pour q
!   radsol       rayonnement net aus sol (LW + SW)
!   dif_grnd     coeff. diffusion vers le sol profond
!
! output:
!   tsurf_new    temperature au sol
!   qsurf        humidite de l'air au dessus du sol
!   fluxsens     flux de chaleur sensible
!   fluxlat      flux de chaleur latente
!   dflux_s      derivee du flux de chaleur sensible / Ts
!   dflux_l      derivee du flux de chaleur latente  / Ts
!

    INCLUDE "YOETHF.h"
    INCLUDE "FCTTRE.h"
    INCLUDE "indicesol.h"
    INCLUDE "YOMCST.h"

! Parametres d'entree
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: knon, nisurf
    REAL   , INTENT(IN)                  :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: petAcoef, peqAcoef
    REAL, DIMENSION(klon), INTENT(IN)    :: petBcoef, peqBcoef
    REAL, DIMENSION(klon), INTENT(IN)    :: ps, q1lay
    REAL, DIMENSION(klon), INTENT(IN)    :: tsurf, p1lay, cal, beta, coef1lay
    REAL, DIMENSION(klon), INTENT(IN)    :: precip_rain, precip_snow ! pas utiles
    REAL, DIMENSION(klon), INTENT(IN)    :: radsol, dif_grnd
    REAL, DIMENSION(klon), INTENT(IN)    :: t1lay, u1lay, v1lay

! Parametres entree-sorties
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(INOUT) :: snow  ! snow pas utile

! Parametres sorties
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)   :: qsurf
    REAL, DIMENSION(klon), INTENT(OUT)   :: tsurf_new, evap, fluxsens, fluxlat
    REAL, DIMENSION(klon), INTENT(OUT)   :: dflux_s, dflux_l

! Variables locales
!****************************************************************************************
    INTEGER                              :: i
    REAL, DIMENSION(klon)                :: zx_mh, zx_nh, zx_oh
    REAL, DIMENSION(klon)                :: zx_mq, zx_nq, zx_oq
    REAL, DIMENSION(klon)                :: zx_pkh, zx_dq_s_dt, zx_qsat, zx_coef
    REAL, DIMENSION(klon)                :: zx_sl, zx_k1
    REAL, DIMENSION(klon)                :: d_ts
    REAL                                 :: zdelta, zcvm5, zx_qs, zcor, zx_dq_s_dh
    REAL                                 :: qsat_new, q1_new
    REAL, PARAMETER                      :: t_grnd = 271.35, t_coup = 273.15
    REAL, PARAMETER                      :: max_eau_sol = 150.0
    CHARACTER (len = 20)                 :: modname = 'calcul_fluxs'
    LOGICAL                              :: fonte_neige
    LOGICAL, SAVE                        :: check = .FALSE.
    !$OMP THREADPRIVATE(check)

! End definition
!****************************************************************************************

    IF (check) WRITE(*,*)'Entree ', modname,' surface = ',nisurf
    
    IF (check) THEN
       WRITE(*,*)' radsol (min, max)', &
            MINVAL(radsol(1:knon)), MAXVAL(radsol(1:knon))
    ENDIF
  
! Traitement neige et humidite du sol
!****************************************************************************************
!
!!$  WRITE(*,*)'test calcul_flux, surface ', nisurf
!!PB test
!!$    if (nisurf == is_oce) then
!!$      snow = 0.
!!$      qsol = max_eau_sol
!!$    else
!!$      where (precip_snow > 0.) snow = snow + (precip_snow * dtime)
!!$      where (snow > epsilon(snow)) snow = max(0.0, snow - (evap * dtime))
!!$!      snow = max(0.0, snow + (precip_snow - evap) * dtime)
!!$      where (precip_rain > 0.) qsol = qsol + (precip_rain - evap) * dtime
!!$    endif 
!!$    IF (nisurf /= is_ter) qsol = max_eau_sol


! 
! Initialisation
!****************************************************************************************
    evap = 0.
    fluxsens=0.
    fluxlat=0.
    dflux_s = 0.
    dflux_l = 0.	
!
! zx_qs = qsat en kg/kg
!****************************************************************************************
    DO i = 1, knon
       zx_pkh(i) = (ps(i)/ps(i))**RKAPPA
       IF (thermcep) THEN
          zdelta=MAX(0.,SIGN(1.,rtt-tsurf(i)))
          zcvm5 = R5LES*RLVTT*(1.-zdelta) + R5IES*RLSTT*zdelta
          zcvm5 = zcvm5 / RCPD / (1.0+RVTMP2*q1lay(i))
          zx_qs= r2es * FOEEW(tsurf(i),zdelta)/ps(i)
          zx_qs=MIN(0.5,zx_qs)
          zcor=1./(1.-retv*zx_qs)
          zx_qs=zx_qs*zcor
          zx_dq_s_dh = FOEDE(tsurf(i),zdelta,zcvm5,zx_qs,zcor) &
               /RLVTT / zx_pkh(i)
       ELSE
          IF (tsurf(i).LT.t_coup) THEN
             zx_qs = qsats(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsats(tsurf(i),zx_qs)/RLVTT &
                  / zx_pkh(i)
          ELSE
             zx_qs = qsatl(tsurf(i)) / ps(i)
             zx_dq_s_dh = dqsatl(tsurf(i),zx_qs)/RLVTT &
                  / zx_pkh(i)
          ENDIF
       ENDIF
       zx_dq_s_dt(i) = RCPD * zx_pkh(i) * zx_dq_s_dh
       zx_qsat(i) = zx_qs
       zx_coef(i) = coef1lay(i) * &
            (1.0+SQRT(u1lay(i)**2+v1lay(i)**2)) * &
            p1lay(i)/(RD*t1lay(i))
       
    ENDDO


! === Calcul de la temperature de surface ===
! zx_sl = chaleur latente d'evaporation ou de sublimation
!****************************************************************************************

    DO i = 1, knon
       zx_sl(i) = RLVTT
       IF (tsurf(i) .LT. RTT) zx_sl(i) = RLSTT
       zx_k1(i) = zx_coef(i)
    ENDDO
    

    DO i = 1, knon
! Q
       zx_oq(i) = 1. - (beta(i) * zx_k1(i) * peqBcoef(i) * dtime)
       zx_mq(i) = beta(i) * zx_k1(i) * &
            (peqAcoef(i) - zx_qsat(i) + &
            zx_dq_s_dt(i) * tsurf(i)) &
            / zx_oq(i)
       zx_nq(i) = beta(i) * zx_k1(i) * (-1. * zx_dq_s_dt(i)) &
            / zx_oq(i)
       
! H
       zx_oh(i) = 1. - (zx_k1(i) * petBcoef(i) * dtime)
       zx_mh(i) = zx_k1(i) * petAcoef(i) / zx_oh(i)
       zx_nh(i) = - (zx_k1(i) * RCPD * zx_pkh(i))/ zx_oh(i)
     
! Tsurface
       tsurf_new(i) = (tsurf(i) + cal(i)/(RCPD * zx_pkh(i)) * dtime * &
            (radsol(i) + zx_mh(i) + zx_sl(i) * zx_mq(i)) & 
            + dif_grnd(i) * t_grnd * dtime)/ &
            ( 1. - dtime * cal(i)/(RCPD * zx_pkh(i)) * ( &
            zx_nh(i) + zx_sl(i) * zx_nq(i)) &  
            + dtime * dif_grnd(i))

!
! Y'a-t-il fonte de neige?
!
!    fonte_neige = (nisurf /= is_oce) .AND. &
!     & (snow(i) > epsfra .OR. nisurf == is_sic .OR. nisurf == is_lic) &
!     & .AND. (tsurf_new(i) >= RTT)
!    if (fonte_neige) tsurf_new(i) = RTT  
       d_ts(i) = tsurf_new(i) - tsurf(i)
!    zx_h_ts(i) = tsurf_new(i) * RCPD * zx_pkh(i)
!    zx_q_0(i) = zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)

!== flux_q est le flux de vapeur d'eau: kg/(m**2 s)  positive vers bas
!== flux_t est le flux de cpt (energie sensible): j/(m**2 s)
       evap(i) = - zx_mq(i) - zx_nq(i) * tsurf_new(i) 
       fluxlat(i) = - evap(i) * zx_sl(i)
       fluxsens(i) = zx_mh(i) + zx_nh(i) * tsurf_new(i)
       
! Derives des flux dF/dTs (W m-2 K-1):
       dflux_s(i) = zx_nh(i)
       dflux_l(i) = (zx_sl(i) * zx_nq(i))

! Nouvelle valeure de l'humidite au dessus du sol
       qsat_new=zx_qsat(i) + zx_dq_s_dt(i) * d_ts(i)
       q1_new = peqAcoef(i) - peqBcoef(i)*evap(i)*dtime
       qsurf(i)=q1_new*(1.-beta(i)) + beta(i)*qsat_new
!
! en cas de fonte de neige
!
!    if (fonte_neige) then
!      bilan_f = radsol(i) + fluxsens(i) - (zx_sl(i) * evap (i)) - &
!     &          dif_grnd(i) * (tsurf_new(i) - t_grnd) - &
!     &          RCPD * (zx_pkh(i))/cal(i)/dtime * (tsurf_new(i) - tsurf(i))
!      bilan_f = max(0., bilan_f)
!      fq_fonte = bilan_f / zx_sl(i)
!      snow(i) = max(0., snow(i) - fq_fonte * dtime)
!      qsol(i) = qsol(i) + (fq_fonte * dtime)
!    endif
!!$    if (nisurf == is_ter)  &
!!$     &  run_off(i) = run_off(i) + max(qsol(i) - max_eau_sol, 0.0)
!!$    qsol(i) = min(qsol(i), max_eau_sol) 
    ENDDO
!
!****************************************************************************************
!
  END SUBROUTINE calcul_fluxs
!
!****************************************************************************************
!
  SUBROUTINE calcul_flux_wind(knon, dtime, &
       u0, v0, u1, v1, cdrag_m, &
       AcoefU, AcoefV, BcoefU, BcoefV, &
       p1lay, t1lay, &
       flux_u1, flux_v1)

    USE dimphy
    INCLUDE "YOMCST.h"

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                  :: knon
    REAL, INTENT(IN)                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)    :: u0, v0  ! u and v at niveau 0
    REAL, DIMENSION(klon), INTENT(IN)    :: u1, v1  ! u and v at niveau 1
    REAL, DIMENSION(klon), INTENT(IN)    :: cdrag_m ! cdrag pour momentum
    REAL, DIMENSION(klon), INTENT(IN)    :: AcoefU, AcoefV, BcoefU, BcoefV
    REAL, DIMENSION(klon), INTENT(IN)    :: p1lay   ! pression 1er niveau (milieu de couche)
    REAL, DIMENSION(klon), INTENT(IN)    :: t1lay   ! temperature 
! Output arguments
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)   :: flux_u1
    REAL, DIMENSION(klon), INTENT(OUT)   :: flux_v1

! Local variables
!****************************************************************************************
    INTEGER                              :: i
    REAL                                 :: mod_wind, buf

!****************************************************************************************
! Calculate the surface flux
!
!****************************************************************************************
    DO i=1,knon
       mod_wind = 1.0 + SQRT((u1(i) - u0(i))**2 + (v1(i)-v0(i))**2)
       buf = cdrag_m(i) * mod_wind * p1lay(i)/(RD*t1lay(i))
       flux_u1(i) = (AcoefU(i) - u0(i)) / (1/buf - BcoefU(i)*dtime )
       flux_v1(i) = (AcoefV(i) - v0(i)) / (1/buf - BcoefV(i)*dtime )
    END DO

  END SUBROUTINE calcul_flux_wind
!
!****************************************************************************************
!
END MODULE calcul_fluxs_mod
