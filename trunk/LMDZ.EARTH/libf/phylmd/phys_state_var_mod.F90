!
! $Id: phys_state_var_mod.F90 1403 2010-07-01 09:02:53Z fairhead $
!
      MODULE phys_state_var_mod
! Variables sauvegardees pour le startphy.nc
!======================================================================
!
!
!======================================================================
! Declaration des variables
      USE dimphy
      INTEGER, PARAMETER :: nlevSTD=17
      INTEGER, PARAMETER :: nlevSTD8=8
      INTEGER, PARAMETER :: nlevSTD3=3
      INTEGER, PARAMETER :: nout=3
      INTEGER, PARAMETER :: napisccp=1
      INTEGER, SAVE :: radpas
!$OMP THREADPRIVATE(radpas)
      REAL, SAVE :: dtime, solaire_etat0
!$OMP THREADPRIVATE(dtime, solaire_etat0)

      REAL, ALLOCATABLE, SAVE :: rlat(:), rlon(:), pctsrf(:,:)
!$OMP THREADPRIVATE(rlat, rlon, pctsrf)
      REAL, ALLOCATABLE, SAVE :: ftsol(:,:)
!$OMP THREADPRIVATE(ftsol)
!      character(len=6), SAVE :: ocean
!!!!!!$OMP THREADPRIVATE(ocean)
!      logical, SAVE :: ok_veget 
!!!!!!$OMP THREADPRIVATE(ok_veget)
      REAL, ALLOCATABLE, SAVE :: falb1(:,:), falb2(:,:)
!$OMP THREADPRIVATE(falb1, falb2)
      REAL, ALLOCATABLE, SAVE :: rain_fall(:), snow_fall(:)
!$OMP THREADPRIVATE( rain_fall, snow_fall)
      REAL, ALLOCATABLE, SAVE :: solsw(:), sollw(:)
!$OMP THREADPRIVATE(solsw, sollw)
      REAL, ALLOCATABLE, SAVE :: radsol(:)
!$OMP THREADPRIVATE(radsol)

!clesphy0 param physiq
!
! Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
!
      REAL, ALLOCATABLE, SAVE :: zmea(:), zstd(:), zsig(:), zgam(:)
!$OMP THREADPRIVATE(zmea, zstd, zsig, zgam)
      REAL, ALLOCATABLE, SAVE :: zthe(:), zpic(:), zval(:)
!$OMP THREADPRIVATE(zthe, zpic, zval)
!     REAL tabcntr0(100)
      REAL, ALLOCATABLE, SAVE :: rugoro(:)
!$OMP THREADPRIVATE(rugoro)
      REAL, ALLOCATABLE, SAVE :: t_ancien(:,:), q_ancien(:,:)
!$OMP THREADPRIVATE(t_ancien, q_ancien)
      REAL, ALLOCATABLE, SAVE :: u_ancien(:,:), v_ancien(:,:)
!$OMP THREADPRIVATE(u_ancien, v_ancien)
      LOGICAL, SAVE :: ancien_ok
!$OMP THREADPRIVATE(ancien_ok)
      REAL, ALLOCATABLE, SAVE :: clwcon(:,:),rnebcon(:,:)
!$OMP THREADPRIVATE(clwcon,rnebcon)
      REAL, ALLOCATABLE, SAVE :: ratqs(:,:)
!$OMP THREADPRIVATE(ratqs)
      REAL, ALLOCATABLE, SAVE :: pbl_tke(:,:,:) ! turb kinetic energy
!$OMP THREADPRIVATE(pbl_tke)
      REAL, ALLOCATABLE, SAVE :: zmax0(:), f0(:) ! 
!$OMP THREADPRIVATE(zmax0,f0)
      REAL, ALLOCATABLE, SAVE :: ema_work1(:,:), ema_work2(:,:)
!$OMP THREADPRIVATE(ema_work1,ema_work2)
      REAL, ALLOCATABLE, SAVE :: entr_therm(:,:), fm_therm(:,:)
!$OMP THREADPRIVATE(entr_therm,fm_therm)
      REAL, ALLOCATABLE, SAVE :: detr_therm(:,:)
!$OMP THREADPRIVATE(detr_therm)
!IM 150408
!     pour phsystoke avec thermiques
      REAL,ALLOCATABLE,SAVE :: clwcon0th(:,:),rnebcon0th(:,:)
!$OMP THREADPRIVATE(clwcon0th,rnebcon0th)
! radiation outputs
      REAL,ALLOCATABLE,SAVE :: swdn0(:,:), swdn(:,:)
!$OMP THREADPRIVATE(swdn0,swdn)
      REAL,ALLOCATABLE,SAVE :: swup0(:,:), swup(:,:)
!$OMP THREADPRIVATE(swup0,swup)
      REAL,ALLOCATABLE,SAVE :: SWdn200clr(:), SWdn200(:)
!$OMP THREADPRIVATE(SWdn200clr,SWdn200)
      REAL,ALLOCATABLE,SAVE :: SWup200clr(:), SWup200(:)
!$OMP THREADPRIVATE(SWup200clr,SWup200)
      REAL,ALLOCATABLE,SAVE :: lwdn0(:,:), lwdn(:,:)
!$OMP THREADPRIVATE(lwdn0,lwdn)
      REAL,ALLOCATABLE,SAVE :: lwup0(:,:), lwup(:,:)
!$OMP THREADPRIVATE(lwup0,lwup)
      REAL,ALLOCATABLE,SAVE :: LWdn200clr(:), LWdn200(:)
!$OMP THREADPRIVATE(LWdn200clr,LWdn200)
      REAL,ALLOCATABLE,SAVE :: LWup200clr(:), LWup200(:)
!$OMP THREADPRIVATE(LWup200clr,LWup200)
      REAL,ALLOCATABLE,SAVE :: LWdnTOA(:), LWdnTOAclr(:)
!$OMP THREADPRIVATE(LWdnTOA,LWdnTOAclr)
! pressure level
      REAL,ALLOCATABLE,SAVE :: tsumSTD(:,:,:)
!$OMP THREADPRIVATE(tsumSTD)
      REAL,ALLOCATABLE,SAVE :: usumSTD(:,:,:), vsumSTD(:,:,:)
!$OMP THREADPRIVATE(usumSTD,vsumSTD)
      REAL,ALLOCATABLE,SAVE :: wsumSTD(:,:,:), phisumSTD(:,:,:)
!$OMP THREADPRIVATE(wsumSTD,phisumSTD)
      REAL,ALLOCATABLE,SAVE :: qsumSTD(:,:,:), rhsumSTD(:,:,:)
!$OMP THREADPRIVATE(qsumSTD,rhsumSTD)
      REAL,ALLOCATABLE,SAVE :: tnondef(:,:,:) 
!$OMP THREADPRIVATE(tnondef)
      REAL,ALLOCATABLE,SAVE :: uvsumSTD(:,:,:)
!$OMP THREADPRIVATE(uvsumSTD)
      REAL,ALLOCATABLE,SAVE :: vqsumSTD(:,:,:)
!$OMP THREADPRIVATE(vqsumSTD)
      REAL,ALLOCATABLE,SAVE :: vTsumSTD(:,:,:)
!$OMP THREADPRIVATE(vTsumSTD)
      REAL,ALLOCATABLE,SAVE :: wqsumSTD(:,:,:)
!$OMP THREADPRIVATE(wqsumSTD)
      REAL,ALLOCATABLE,SAVE :: vphisumSTD(:,:,:)
!$OMP THREADPRIVATE(vphisumSTD)
      REAL,ALLOCATABLE,SAVE :: wTsumSTD(:,:,:)
!$OMP THREADPRIVATE(wTsumSTD)
      REAL,ALLOCATABLE,SAVE :: u2sumSTD(:,:,:)
!$OMP THREADPRIVATE(u2sumSTD)
      REAL,ALLOCATABLE,SAVE :: v2sumSTD(:,:,:)
!$OMP THREADPRIVATE(v2sumSTD)
      REAL,ALLOCATABLE,SAVE :: T2sumSTD(:,:,:)
!$OMP THREADPRIVATE(T2sumSTD)
      REAL,ALLOCATABLE,SAVE :: O3sumSTD(:,:,:), O3daysumSTD(:,:,:)
!$OMP THREADPRIVATE(O3sumSTD,O3daysumSTD) 
!IM begin
      REAL,ALLOCATABLE,SAVE :: wlevSTD(:,:), ulevSTD(:,:), vlevSTD(:,:)
!$OMP THREADPRIVATE(wlevSTD,ulevSTD,vlevSTD)
      REAL,ALLOCATABLE,SAVE :: tlevSTD(:,:), qlevSTD(:,:), rhlevSTD(:,:)
!$OMP THREADPRIVATE(tlevSTD,qlevSTD,rhlevSTD)
      REAL,ALLOCATABLE,SAVE :: philevSTD(:,:)
!$OMP THREADPRIVATE(philevSTD)
      REAL,ALLOCATABLE,SAVE :: uvSTD(:,:)
!$OMP THREADPRIVATE(uvSTD)
      REAL,ALLOCATABLE,SAVE :: vqSTD(:,:)
!$OMP THREADPRIVATE(vqSTD)
      REAL,ALLOCATABLE,SAVE :: vTSTD(:,:)
!$OMP THREADPRIVATE(vTSTD)
      REAL,ALLOCATABLE,SAVE :: wqSTD(:,:)
!$OMP THREADPRIVATE(wqSTD)
      REAL,ALLOCATABLE,SAVE :: vphiSTD(:,:)
!$OMP THREADPRIVATE(vphiSTD)
      REAL,ALLOCATABLE,SAVE :: wTSTD(:,:)
!$OMP THREADPRIVATE(wTSTD)
      REAL,ALLOCATABLE,SAVE :: u2STD(:,:)
!$OMP THREADPRIVATE(u2STD)
      REAL,ALLOCATABLE,SAVE :: v2STD(:,:) 
!$OMP THREADPRIVATE(v2STD)
      REAL,ALLOCATABLE,SAVE :: T2STD(:,:)
!$OMP THREADPRIVATE(T2STD)
      REAL,ALLOCATABLE,SAVE :: O3STD(:,:), O3daySTD(:,:)
!$OMP THREADPRIVATE(O3STD,O3daySTD)
!IM end
      INTEGER,ALLOCATABLE,SAVE :: seed_old(:,:)
!$OMP THREADPRIVATE(seed_old)
      REAL,ALLOCATABLE,SAVE :: zuthe(:),zvthe(:)
!$OMP THREADPRIVATE(zuthe,zvthe)
      REAL,ALLOCATABLE,SAVE :: alb_neig(:)
!$OMP THREADPRIVATE(alb_neig)
!cloud base mass flux
      REAL,ALLOCATABLE,SAVE :: ema_cbmf(:)
!$OMP THREADPRIVATE(ema_cbmf)
!cloud base pressure & cloud top pressure
      REAL,ALLOCATABLE,SAVE :: ema_pcb(:), ema_pct(:)
!$OMP THREADPRIVATE(ema_pcb,ema_pct)
      REAL,ALLOCATABLE,SAVE :: Ma(:,:)        ! undilute upward mass flux
!$OMP THREADPRIVATE(Ma)
      REAL,ALLOCATABLE,SAVE :: qcondc(:,:)    ! in-cld water content from convect
!$OMP THREADPRIVATE(qcondc)
      REAL,ALLOCATABLE,SAVE :: wd(:) ! sb
!$OMP THREADPRIVATE(wd)
      REAL,ALLOCATABLE,SAVE :: sigd(:)
!$OMP THREADPRIVATE(sigd)
!
      REAL,ALLOCATABLE,SAVE :: cin(:)
!$OMP THREADPRIVATE(cin)
! ftd : differential heating between wake and environment
      REAL,ALLOCATABLE,SAVE :: ftd(:,:)
!$OMP THREADPRIVATE(ftd)
! fqd : differential moistening between wake and environment
      REAL,ALLOCATABLE,SAVE :: fqd(:,:)     
!$OMP THREADPRIVATE(fqd)
!34EK
! -- Variables de controle de ALE et ALP
!ALE : Energie disponible pour soulevement : utilisee par la 
!      convection d'Emanuel pour le declenchement et la regulation
      REAL,ALLOCATABLE,SAVE :: ALE(:)
!$OMP THREADPRIVATE(ALE)
!ALP : Puissance  disponible pour soulevement
      REAL,ALLOCATABLE,SAVE :: ALP(:)
!$OMP THREADPRIVATE(ALP)
!
! nouvelles variables pour le couplage convection-couche limite
      REAL,ALLOCATABLE,SAVE :: Ale_bl(:)
!$OMP THREADPRIVATE(Ale_bl)
      REAL,ALLOCATABLE,SAVE :: Alp_bl(:)
!$OMP THREADPRIVATE(Alp_bl)
      INTEGER,ALLOCATABLE,SAVE :: lalim_conv(:)
!$OMP THREADPRIVATE(lalim_conv)
      REAL,ALLOCATABLE,SAVE :: wght_th(:,:)
!$OMP THREADPRIVATE(wght_th)
!
! variables de la wake
! wake_deltat : ecart de temperature avec la zone non perturbee
! wake_deltaq : ecart d'humidite avec la zone non perturbee
! wake_Cstar  : vitesse d'etalement de la poche
! wake_s      : fraction surfacique occupee par la poche froide
! wake_pe     : wake potential energy - WAPE
! wake_fip    : Gust Front Impinging power - ALP
! dt_wake, dq_wake: LS tendencies due to wake
      REAL,ALLOCATABLE,SAVE :: wake_deltat(:,:)
!$OMP THREADPRIVATE(wake_deltat)
      REAL,ALLOCATABLE,SAVE :: wake_deltaq(:,:)
!$OMP THREADPRIVATE(wake_deltaq)
      REAL,ALLOCATABLE,SAVE :: wake_Cstar(:)
!$OMP THREADPRIVATE(wake_Cstar)
      REAL,ALLOCATABLE,SAVE :: wake_s(:)
!$OMP THREADPRIVATE(wake_s)
      REAL,ALLOCATABLE,SAVE :: wake_pe(:)
!$OMP THREADPRIVATE(wake_pe)
      REAL,ALLOCATABLE,SAVE :: wake_fip(:)
!$OMP THREADPRIVATE(wake_fip)
      REAL,ALLOCATABLE,SAVE :: dt_wake(:,:)
!$OMP THREADPRIVATE(dt_wake)
      REAL,ALLOCATABLE,SAVE :: dq_wake(:,:)
!$OMP THREADPRIVATE(dq_wake)
!
! pfrac_impa : Produits des coefs lessivage impaction
! pfrac_nucl : Produits des coefs lessivage nucleation
! pfrac_1nucl: Produits des coefs lessi nucl (alpha = 1) 
      REAL,ALLOCATABLE,SAVE :: pfrac_impa(:,:), pfrac_nucl(:,:)
!$OMP THREADPRIVATE(pfrac_impa,pfrac_nucl)
      REAL,ALLOCATABLE,SAVE :: pfrac_1nucl(:,:)
!$OMP THREADPRIVATE(pfrac_1nucl)
!
      REAL,ALLOCATABLE,SAVE :: total_rain(:), nday_rain(:)  
!$OMP THREADPRIVATE(total_rain,nday_rain)
      REAL,ALLOCATABLE,SAVE :: paire_ter(:)
!$OMP THREADPRIVATE(paire_ter)
! albsol1: albedo du sol total pour SW visible
! albsol2: albedo du sol total pour SW proche IR
      REAL,ALLOCATABLE,SAVE :: albsol1(:), albsol2(:)
!$OMP THREADPRIVATE(albsol1,albsol2)

      REAL, ALLOCATABLE, SAVE:: wo(:, :, :)
      ! column-density of ozone in a layer, in kilo-Dobsons
      ! Third dimension has size 1 or 2.
      ! "wo(:, :, 1)" is for the average day-night field, 
      ! "wo(:, :, 2)" is for daylight time.
      !$OMP THREADPRIVATE(wo)

! heat : chauffage solaire
! heat0: chauffage solaire ciel clair
! cool : refroidissement infrarouge
! cool0 : refroidissement infrarouge ciel clair
! sollwdown : downward LW flux at surface
! sollwdownclr : downward CS LW flux at surface
! toplwdown : downward CS LW flux at TOA
! toplwdownclr : downward CS LW flux at TOA
      REAL,ALLOCATABLE,SAVE :: clwcon0(:,:),rnebcon0(:,:)
!$OMP THREADPRIVATE(clwcon0,rnebcon0)
      REAL,ALLOCATABLE,SAVE :: heat(:,:)   
!$OMP THREADPRIVATE(heat)
      REAL,ALLOCATABLE,SAVE :: heat0(:,:)
!$OMP THREADPRIVATE(heat0)
      REAL,ALLOCATABLE,SAVE :: cool(:,:)
!$OMP THREADPRIVATE(cool)
      REAL,ALLOCATABLE,SAVE :: cool0(:,:)
!$OMP THREADPRIVATE(cool0)
      REAL,ALLOCATABLE,SAVE :: topsw(:), toplw(:)
!$OMP THREADPRIVATE(topsw,toplw)
      REAL,ALLOCATABLE,SAVE :: sollwdown(:)
!$OMP THREADPRIVATE(sollwdown)
      REAL,ALLOCATABLE,SAVE :: sollwdownclr(:)
!$OMP THREADPRIVATE(sollwdownclr)
      REAL,ALLOCATABLE,SAVE :: toplwdown(:)
!$OMP THREADPRIVATE(toplwdown)
      REAL,ALLOCATABLE,SAVE :: toplwdownclr(:)
!$OMP THREADPRIVATE(toplwdownclr)
      REAL,ALLOCATABLE,SAVE :: topsw0(:),toplw0(:),solsw0(:),sollw0(:)
!$OMP THREADPRIVATE(topsw0,toplw0,solsw0,sollw0)
      REAL,ALLOCATABLE,SAVE :: albpla(:)
!$OMP THREADPRIVATE(albpla)
! pbase : cloud base pressure
! bbase : cloud base buoyancy
      REAL,ALLOCATABLE,SAVE :: cape(:)
!$OMP THREADPRIVATE(cape)
      REAL,ALLOCATABLE,SAVE :: pbase(:)
!$OMP THREADPRIVATE(pbase)
      REAL,ALLOCATABLE,SAVE :: bbase(:)
!$OMP THREADPRIVATE(bbase)
!
      REAL,SAVE,ALLOCATABLE :: zqasc(:,:)
!$OMP THREADPRIVATE( zqasc)
      INTEGER,ALLOCATABLE,SAVE :: ibas_con(:), itop_con(:)
!$OMP THREADPRIVATE(ibas_con,itop_con)
      REAL,SAVE,ALLOCATABLE :: rain_con(:)
!$OMP THREADPRIVATE(rain_con)
      REAL,SAVE,ALLOCATABLE :: snow_con(:)
!$OMP THREADPRIVATE(snow_con)
!
      REAL,SAVE,ALLOCATABLE :: rlonPOS(:)
!$OMP THREADPRIVATE(rlonPOS)
      REAL,SAVE,ALLOCATABLE :: newsst(:)
!$OMP THREADPRIVATE(newsst)
      REAL,SAVE,ALLOCATABLE :: u10m(:,:), v10m(:,:)
!$OMP THREADPRIVATE(u10m,v10m)
!
! ok_ade=T -ADE=topswad-topsw
! ok_aie=T ->
!       ok_ade=T -AIE=topswai-topswad
!       ok_ade=F -AIE=topswai-topsw
!
!topswad, solswad : Aerosol direct effect
      REAL,SAVE,ALLOCATABLE :: topswad(:), solswad(:)
!$OMP THREADPRIVATE(topswad,solswad)
!topswai, solswai : Aerosol indirect effect
      REAL,SAVE,ALLOCATABLE :: topswai(:), solswai(:)
!$OMP THREADPRIVATE(topswai,solswai)

      REAL,SAVE,ALLOCATABLE :: tau_aero(:,:,:,:), piz_aero(:,:,:,:), cg_aero(:,:,:,:)
!$OMP THREADPRIVATE(tau_aero, piz_aero, cg_aero)
      REAL,SAVE,ALLOCATABLE :: ccm(:,:,:)
!$OMP THREADPRIVATE(ccm)

CONTAINS

!======================================================================
SUBROUTINE phys_state_var_init(read_climoz)
use dimphy
USE control_mod
use aero_mod
IMPLICIT NONE

integer, intent(in)::  read_climoz
! read ozone climatology
! Allowed values are 0, 1 and 2
! 0: do not read an ozone climatology
! 1: read a single ozone climatology that will be used day and night
! 2: read two ozone climatologies, the average day and night
! climatology and the daylight climatology

#include "indicesol.h"
      ALLOCATE(rlat(klon), rlon(klon))
      ALLOCATE(pctsrf(klon,nbsrf))
      ALLOCATE(ftsol(klon,nbsrf))
      ALLOCATE(falb1(klon,nbsrf))
      ALLOCATE(falb2(klon,nbsrf))
      ALLOCATE(rain_fall(klon))
      ALLOCATE(snow_fall(klon))
      ALLOCATE(solsw(klon), sollw(klon))
      ALLOCATE(radsol(klon))
      ALLOCATE(zmea(klon), zstd(klon), zsig(klon), zgam(klon))
      ALLOCATE(zthe(klon), zpic(klon), zval(klon))

      ALLOCATE(rugoro(klon))
      ALLOCATE(t_ancien(klon,klev), q_ancien(klon,klev))
      ALLOCATE(u_ancien(klon,klev), v_ancien(klon,klev))
      ALLOCATE(clwcon(klon,klev),rnebcon(klon,klev))
      ALLOCATE(ratqs(klon,klev))
      ALLOCATE(pbl_tke(klon,klev+1,nbsrf))
      ALLOCATE(zmax0(klon), f0(klon))
      ALLOCATE(ema_work1(klon,klev), ema_work2(klon,klev))
      ALLOCATE(entr_therm(klon,klev), fm_therm(klon,klev+1))
      ALLOCATE(detr_therm(klon,klev))
!     pour phsystoke avec thermiques
      ALLOCATE(clwcon0th(klon,klev),rnebcon0th(klon,klev))
! radiation outputs
      ALLOCATE(swdn0(klon,klevp1), swdn(klon,klevp1))
      ALLOCATE(swup0(klon,klevp1), swup(klon,klevp1))
      ALLOCATE(lwdn0(klon,klevp1), lwdn(klon,klevp1))
      ALLOCATE(lwup0(klon,klevp1), lwup(klon,klevp1))
      ALLOCATE(SWdn200clr(klon), SWdn200(klon))
      ALLOCATE(SWup200clr(klon), SWup200(klon))
      ALLOCATE(LWdn200clr(klon), LWdn200(klon))
      ALLOCATE(LWup200clr(klon), LWup200(klon))
      ALLOCATE(LWdnTOA(klon), LWdnTOAclr(klon))
! pressure level
      ALLOCATE(tsumSTD(klon,nlevSTD,nout))
      ALLOCATE(usumSTD(klon,nlevSTD,nout), vsumSTD(klon,nlevSTD,nout))
      ALLOCATE(wsumSTD(klon,nlevSTD,nout), phisumSTD(klon,nlevSTD,nout))
      ALLOCATE(qsumSTD(klon,nlevSTD,nout), rhsumSTD(klon,nlevSTD,nout))
      ALLOCATE(tnondef(klon,nlevSTD,nout))
      ALLOCATE(uvsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vqsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vTsumSTD(klon,nlevSTD,nout))
      ALLOCATE(wqsumSTD(klon,nlevSTD,nout))
      ALLOCATE(vphisumSTD(klon,nlevSTD,nout))
      ALLOCATE(wTsumSTD(klon,nlevSTD,nout))
      ALLOCATE(u2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(v2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(T2sumSTD(klon,nlevSTD,nout))
      ALLOCATE(O3sumSTD(klon,nlevSTD,nout))
      ALLOCATE(O3daysumSTD(klon,nlevSTD,nout))
!IM beg
      ALLOCATE(wlevSTD(klon,nlevSTD), ulevSTD(klon,nlevSTD), vlevSTD(klon,nlevSTD))
      ALLOCATE(tlevSTD(klon,nlevSTD), qlevSTD(klon,nlevSTD), rhlevSTD(klon,nlevSTD))
      ALLOCATE(philevSTD(klon,nlevSTD))
      ALLOCATE(uvSTD(klon,nlevSTD),vqSTD(klon,nlevSTD))
      ALLOCATE(vTSTD(klon,nlevSTD),wqSTD(klon,nlevSTD))
      ALLOCATE(vphiSTD(klon,nlevSTD),wTSTD(klon,nlevSTD))
      ALLOCATE(u2STD(klon,nlevSTD),v2STD(klon,nlevSTD))
      ALLOCATE(T2STD(klon,nlevSTD))
      ALLOCATE(O3STD(klon,nlevSTD))
      ALLOCATE(O3daySTD(klon,nlevSTD))
!IM end
      ALLOCATE(seed_old(klon,napisccp))
      ALLOCATE(zuthe(klon),zvthe(klon))
      ALLOCATE(alb_neig(klon))
!cloud base mass flux
      ALLOCATE(ema_cbmf(klon))
!cloud base pressure & cloud top pressure
      ALLOCATE(ema_pcb(klon), ema_pct(klon))
!
      ALLOCATE(Ma(klon,klev))
      ALLOCATE(qcondc(klon,klev))
      ALLOCATE(wd(klon))
      ALLOCATE(sigd(klon))
      ALLOCATE(cin(klon), ALE(klon), ALP(klon))
      ALLOCATE(ftd(klon,klev), fqd(klon,klev))
      ALLOCATE(Ale_bl(klon))
      ALLOCATE(Alp_bl(klon))
      ALLOCATE(lalim_conv(klon))
      ALLOCATE(wght_th(klon,klev))
      ALLOCATE(wake_deltat(klon,klev), wake_deltaq(klon,klev))
      ALLOCATE(wake_Cstar(klon), wake_s(klon))
      ALLOCATE(wake_pe(klon), wake_fip(klon))
      ALLOCATE(dt_wake(klon,klev), dq_wake(klon,klev))
      ALLOCATE(pfrac_impa(klon,klev), pfrac_nucl(klon,klev))
      ALLOCATE(pfrac_1nucl(klon,klev))
      ALLOCATE(total_rain(klon), nday_rain(klon))
      ALLOCATE(paire_ter(klon))
      ALLOCATE(albsol1(klon), albsol2(klon))

      if (read_climoz <= 1) then
         ALLOCATE(wo(klon,klev, 1))
      else
         ! read_climoz == 2
         ALLOCATE(wo(klon,klev, 2))
      end if
      
      ALLOCATE(clwcon0(klon,klev),rnebcon0(klon,klev))
      ALLOCATE(heat(klon,klev), heat0(klon,klev)) 
      ALLOCATE(cool(klon,klev), cool0(klon,klev))
      ALLOCATE(topsw(klon), toplw(klon))
      ALLOCATE(sollwdown(klon), sollwdownclr(klon))
      ALLOCATE(toplwdown(klon), toplwdownclr(klon))
      ALLOCATE(topsw0(klon),toplw0(klon),solsw0(klon),sollw0(klon))
      ALLOCATE(albpla(klon))
      ALLOCATE(cape(klon))
      ALLOCATE(pbase(klon),bbase(klon))
      ALLOCATE(zqasc(klon,klev))
      ALLOCATE(ibas_con(klon), itop_con(klon))
      ALLOCATE(rain_con(klon), snow_con(klon))
      ALLOCATE(rlonPOS(klon))
      ALLOCATE(newsst(klon))
      ALLOCATE(u10m(klon,nbsrf), v10m(klon,nbsrf))
      ALLOCATE(topswad(klon), solswad(klon))
      ALLOCATE(topswai(klon), solswai(klon))
      ALLOCATE(tau_aero(klon,klev,naero_grp,nbands),piz_aero(klon,klev,naero_grp,nbands),cg_aero(klon,klev,naero_grp,nbands))
      ALLOCATE(ccm(klon,klev,nbands))

END SUBROUTINE phys_state_var_init

!======================================================================
SUBROUTINE phys_state_var_end
use dimphy
use control_mod
IMPLICIT NONE
#include "indicesol.h"

      deallocate(rlat, rlon, pctsrf, ftsol, falb1, falb2)
      deallocate(rain_fall, snow_fall, solsw, sollw, radsol)
      deallocate(zmea, zstd, zsig, zgam)
      deallocate(zthe, zpic, zval)
      deallocate(rugoro, t_ancien, q_ancien, clwcon, rnebcon)
      deallocate(        u_ancien, v_ancien                 )
      deallocate(ratqs, pbl_tke)
      deallocate(zmax0, f0)
      deallocate(ema_work1, ema_work2)
      deallocate(entr_therm, fm_therm)
      deallocate(detr_therm)
      deallocate(clwcon0th, rnebcon0th)
! radiation outputs
      deallocate(swdn0, swdn)
      deallocate(swup0, swup)
      deallocate(lwdn0, lwdn)
      deallocate(lwup0, lwup)
      deallocate(SWdn200clr, SWdn200)
      deallocate(SWup200clr, SWup200)
      deallocate(LWdn200clr, LWdn200)
      deallocate(LWup200clr, LWup200)
      deallocate(LWdnTOA, LWdnTOAclr)
! pressure level
      deallocate(tsumSTD)
      deallocate(usumSTD, vsumSTD)
      deallocate(wsumSTD, phisumSTD)
      deallocate(tnondef)
      deallocate(qsumSTD, rhsumSTD)
      deallocate(uvsumSTD)
      deallocate(vqsumSTD)
      deallocate(vTsumSTD)
      deallocate(wqsumSTD)
      deallocate(vphisumSTD)
      deallocate(wTsumSTD)
      deallocate(u2sumSTD)
      deallocate(v2sumSTD)
      deallocate(T2sumSTD)
      deallocate(O3sumSTD)
      deallocate(O3daysumSTD)
!IM beg
      deallocate(wlevSTD,ulevSTD,vlevSTD,tlevSTD,qlevSTD,rhlevSTD,philevSTD)
      deallocate(uvSTD,vqSTD,vTSTD,wqSTD,vphiSTD,wTSTD,u2STD,v2STD,T2STD,O3STD,O3daySTD)
!IM end
      deallocate(seed_old)
      deallocate(zuthe, zvthe)
      deallocate(alb_neig)
      deallocate(ema_cbmf)
      deallocate(ema_pcb, ema_pct)
      deallocate(Ma, qcondc)
      deallocate(wd, sigd)
      deallocate(cin, ALE, ALP)
      deallocate(ftd, fqd)
      deallocate(Ale_bl, Alp_bl)
      deallocate(lalim_conv, wght_th)
      deallocate(wake_deltat, wake_deltaq)
      deallocate(wake_Cstar, wake_s, wake_pe, wake_fip)
      deallocate(dt_wake, dq_wake)
      deallocate(pfrac_impa, pfrac_nucl)
      deallocate(pfrac_1nucl)
      deallocate(total_rain, nday_rain)
      deallocate(paire_ter)
      deallocate(albsol1, albsol2)
      deallocate(wo)
      deallocate(clwcon0,rnebcon0)
      deallocate(heat, heat0) 
      deallocate(cool, cool0)
      deallocate(topsw, toplw)
      deallocate(sollwdown, sollwdownclr)
      deallocate(toplwdown, toplwdownclr)
      deallocate(topsw0,toplw0,solsw0,sollw0)
      deallocate(albpla)
      deallocate(cape)
      deallocate(pbase,bbase)
      deallocate(zqasc)
      deallocate(ibas_con, itop_con)
      deallocate(rain_con, snow_con)
      deallocate(rlonPOS)
      deallocate(newsst)
      deallocate(u10m, v10m)
      deallocate(topswad, solswad)
      deallocate(topswai, solswai)
      deallocate(tau_aero,piz_aero,cg_aero)
      deallocate(ccm)
       
END SUBROUTINE phys_state_var_end

      END MODULE phys_state_var_mod
