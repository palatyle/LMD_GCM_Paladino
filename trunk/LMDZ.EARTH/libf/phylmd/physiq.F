! $Id: physiq.F 1486 2011-02-11 12:07:39Z fairhead $
c#define IO_DEBUG

      SUBROUTINE physiq (nlon,nlev,
     .            debut,lafin,jD_cur, jH_cur,pdtphys,
     .            paprs,pplay,pphi,pphis,presnivs,clesphy0,
     .            u,v,t,qx,
     .            flxmass_w,
     .            d_u, d_v, d_t, d_qx, d_ps
     .            , dudyn
     .            , PVteta)

      USE ioipsl, only: histbeg, histvert, histdef, histend, histsync,
     $     histwrite, ju2ymds, ymds2ju, ioget_year_len
      USE comgeomphy
      USE phys_cal_mod
      USE write_field_phy
      USE dimphy
      USE infotrac
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para
      USE iophy
      USE misc_mod, mydebug=>debug
      USE vampir
      USE pbl_surface_mod, ONLY : pbl_surface
      USE change_srf_frac_mod
      USE surface_data,     ONLY : type_ocean, ok_veget
      USE phys_local_var_mod ! Variables internes non sauvegardees de la physique
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      USE phys_output_var_mod ! Variables pour les ecritures des sorties
      USE fonte_neige_mod, ONLY  : fonte_neige_get_vars
      USE phys_output_mod
      use open_climoz_m, only: open_climoz ! ozone climatology from a file
      use regr_pr_av_m, only: regr_pr_av
      use netcdf95, only: nf95_close
cIM for NMC files
      use netcdf, only: nf90_fill_real
      use mod_phys_lmdz_mpi_data, only: is_mpi_root
      USE aero_mod
      use ozonecm_m, only: ozonecm ! ozone of J.-F. Royer
      use conf_phys_m, only: conf_phys
      use radlwsw_m, only: radlwsw
      USE control_mod


      IMPLICIT none
c======================================================================
c
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c
c Objet: Moniteur general de la physique du modele
cAA      Modifications quant aux traceurs :
cAA                  -  uniformisation des parametrisations ds phytrac
cAA                  -  stockage des moyennes des champs necessaires
cAA                     en mode traceur off-line 
c======================================================================
c   CLEFS CPP POUR LES IO
c   =====================
#define histNMC
c#define histISCCP
c======================================================================
c    modif   ( P. Le Van ,  12/10/98 )
c
c  Arguments:
c
c nlon----input-I-nombre de points horizontaux
c nlev----input-I-nombre de couches verticales, doit etre egale a klev
c debut---input-L-variable logique indiquant le premier passage
c lafin---input-L-variable logique indiquant le dernier passage
c jD_cur       -R-jour courant a l'appel de la physique (jour julien)
c jH_cur       -R-heure courante a l'appel de la physique (jour julien)
c pdtphys-input-R-pas d'integration pour la physique (seconde)
c paprs---input-R-pression pour chaque inter-couche (en Pa)
c pplay---input-R-pression pour le mileu de chaque couche (en Pa)
c pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
c pphis---input-R-geopotentiel du sol
c presnivs-input_R_pressions approximat. des milieux couches ( en PA)
c u-------input-R-vitesse dans la direction X (de O a E) en m/s
c v-------input-R-vitesse Y (de S a N) en m/s
c t-------input-R-temperature (K)
c qx------input-R-humidite specifique (kg/kg) et d'autres traceurs
c d_t_dyn-input-R-tendance dynamique pour "t" (K/s)
c d_q_dyn-input-R-tendance dynamique pour "q" (kg/kg/s)
c flxmass_w -input-R- flux de masse verticale
c d_u-----output-R-tendance physique de "u" (m/s/s)
c d_v-----output-R-tendance physique de "v" (m/s/s)
c d_t-----output-R-tendance physique de "t" (K/s)
c d_qx----output-R-tendance physique de "qx" (kg/kg/s)
c d_ps----output-R-tendance physique de la pression au sol
cIM
c PVteta--output-R-vorticite potentielle a des thetas constantes
c======================================================================
#include "dimensions.h"
      integer jjmp1
      parameter (jjmp1=jjm+1-1/jjm)
      integer iip1
      parameter (iip1=iim+1)

#include "regdim.h"
#include "indicesol.h"
#include "dimsoil.h"
#include "clesphys.h"
#include "temps.h"
#include "iniprint.h"
#include "thermcell.h"
c======================================================================
      LOGICAL ok_cvl  ! pour activer le nouveau driver pour convection KE
      PARAMETER (ok_cvl=.TRUE.)
      LOGICAL ok_gust ! pour activer l'effet des gust sur flux surface
      PARAMETER (ok_gust=.FALSE.)
      integer iflag_radia     ! active ou non le rayonnement (MPL)
      save iflag_radia
c$OMP THREADPRIVATE(iflag_radia)
c======================================================================
      LOGICAL check ! Verifier la conservation du modele en eau
      PARAMETER (check=.FALSE.)
      LOGICAL ok_stratus ! Ajouter artificiellement les stratus
      PARAMETER (ok_stratus=.FALSE.)
c======================================================================
      REAL amn, amx
      INTEGER igout
c======================================================================
c Clef controlant l'activation du cycle diurne:
ccc      LOGICAL cycle_diurne
ccc      PARAMETER (cycle_diurne=.FALSE.)
c======================================================================
c Modele thermique du sol, a activer pour le cycle diurne:
ccc      LOGICAL soil_model
ccc      PARAMETER (soil_model=.FALSE.)
c======================================================================
c Dans les versions precedentes, l'eau liquide nuageuse utilisee dans
c le calcul du rayonnement est celle apres la precipitation des nuages.
c Si cette cle new_oliq est activee, ce sera une valeur moyenne entre
c la condensation et la precipitation. Cette cle augmente les impacts
c radiatifs des nuages.
ccc      LOGICAL new_oliq
ccc      PARAMETER (new_oliq=.FALSE.)
c======================================================================
c Clefs controlant deux parametrisations de l'orographie:
cc      LOGICAL ok_orodr
ccc      PARAMETER (ok_orodr=.FALSE.)
ccc      LOGICAL ok_orolf
ccc      PARAMETER (ok_orolf=.FALSE.)
c======================================================================
      LOGICAL ok_journe ! sortir le fichier journalier
      save ok_journe
c$OMP THREADPRIVATE(ok_journe)
c
      LOGICAL ok_mensuel ! sortir le fichier mensuel
      save ok_mensuel
c$OMP THREADPRIVATE(ok_mensuel)
c
      LOGICAL ok_instan ! sortir le fichier instantane
      save ok_instan
c$OMP THREADPRIVATE(ok_instan)
c
      LOGICAL ok_LES ! sortir le fichier LES 
      save ok_LES                            
c$OMP THREADPRIVATE(ok_LES)                  
c
      LOGICAL callstats ! sortir le fichier stats 
      save callstats                            
c$OMP THREADPRIVATE(callstats)                  
c
      LOGICAL ok_region ! sortir le fichier regional
      PARAMETER (ok_region=.FALSE.)
c======================================================================
      real weak_inversion(klon),dthmin(klon)
      real seuil_inversion
      save seuil_inversion
c$OMP THREADPRIVATE(seuil_inversion)
      integer iflag_ratqs
      save iflag_ratqs
c$OMP THREADPRIVATE(iflag_ratqs)
      REAL lambda_th(klon,klev),zz,znum,zden
      REAL wmax_th(klon)
      REAL zmax_th(klon)
      REAL tau_overturning_th(klon)

      integer lmax_th(klon)
      integer limbas(klon)
      real ratqscth(klon,klev)
      real ratqsdiff(klon,klev)
      real zqsatth(klon,klev)

c======================================================================
c
      INTEGER ivap          ! indice de traceurs pour vapeur d'eau
      PARAMETER (ivap=1)
      INTEGER iliq          ! indice de traceurs pour eau liquide
      PARAMETER (iliq=2)

c
c
c Variables argument:
c
      INTEGER nlon
      INTEGER nlev
      REAL, intent(in):: jD_cur, jH_cur

      REAL pdtphys
      LOGICAL debut, lafin
      REAL paprs(klon,klev+1)
      REAL pplay(klon,klev)
      REAL pphi(klon,klev)
      REAL pphis(klon)
      REAL presnivs(klev)
      REAL znivsig(klev)
      real pir

      REAL u(klon,klev)
      REAL v(klon,klev)
      REAL t(klon,klev),theta(klon,klev)
      REAL qx(klon,klev,nqtot)
      REAL flxmass_w(klon,klev)
      REAL omega(klon,klev) ! vitesse verticale en Pa/s
      REAL d_u(klon,klev)
      REAL d_v(klon,klev)
      REAL d_t(klon,klev)
      REAL d_qx(klon,klev,nqtot)
      REAL d_ps(klon)
      real da(klon,klev),phi(klon,klev,klev),mp(klon,klev)
!IM definition dynamique o_trac dans phys_output_open
!      type(ctrl_out) :: o_trac(nqtot)
c
cIM Amip2 PV a theta constante 
c
      INTEGER nbteta
      PARAMETER(nbteta=3)
      CHARACTER*3 ctetaSTD(nbteta)
      DATA ctetaSTD/'350','380','405'/
      SAVE ctetaSTD
c$OMP THREADPRIVATE(ctetaSTD)
      REAL rtetaSTD(nbteta)
      DATA rtetaSTD/350., 380., 405./
      SAVE rtetaSTD
c$OMP THREADPRIVATE(rtetaSTD)     
c
      REAL PVteta(klon,nbteta)
      REAL zx_tmp_3dte(iim,jjmp1,nbteta)
c
cMI Amip2 PV a theta constante

cym      INTEGER klevp1, klevm1
cym      PARAMETER(klevp1=klev+1,klevm1=klev-1)
cym#include "raddim.h"
c
c
cIM Amip2
c variables a une pression donnee
c
      real rlevSTD(nlevSTD)
      DATA rlevSTD/100000., 92500., 85000., 70000.,
     .60000., 50000., 40000., 30000., 25000., 20000.,
     .15000., 10000., 7000., 5000., 3000., 2000., 1000./
      SAVE rlevstd
c$OMP THREADPRIVATE(rlevstd)
      CHARACTER*4 clevSTD(nlevSTD)
      DATA clevSTD/'1000','925 ','850 ','700 ','600 ',
     .'500 ','400 ','300 ','250 ','200 ','150 ','100 ',
     .'70  ','50  ','30  ','20  ','10  '/
      SAVE clevSTD
c$OMP THREADPRIVATE(clevSTD)
c
      CHARACTER*4 bb2
      CHARACTER*2 bb3

      real twriteSTD(klon,nlevSTD,nfiles)
      real qwriteSTD(klon,nlevSTD,nfiles)
      real rhwriteSTD(klon,nlevSTD,nfiles)
      real phiwriteSTD(klon,nlevSTD,nfiles)
      real uwriteSTD(klon,nlevSTD,nfiles)
      real vwriteSTD(klon,nlevSTD,nfiles)
      real wwriteSTD(klon,nlevSTD,nfiles)
cIM for NMC files
      REAL geo500(klon)
      real :: rlevSTD3(nlevSTD3)
      DATA rlevSTD3/85000., 50000., 25000./
      SAVE rlevSTD3
c$OMP THREADPRIVATE(rlevSTD3)
      real :: rlevSTD8(nlevSTD8)
      DATA rlevSTD8/100000., 85000., 70000., 50000., 25000., 10000.,
     $     5000., 1000./
      SAVE rlevSTD8
c$OMP THREADPRIVATE(rlevSTD8) 
      real twriteSTD3(klon,nlevSTD3)
      real qwriteSTD3(klon,nlevSTD3)
      real rhwriteSTD3(klon,nlevSTD3)
      real phiwriteSTD3(klon,nlevSTD3)
      real uwriteSTD3(klon,nlevSTD3)
      real vwriteSTD3(klon,nlevSTD3)
      real wwriteSTD3(klon,nlevSTD3)
c
      real tnondefSTD8(klon,nlevSTD8)
      real twriteSTD8(klon,nlevSTD8)
      real qwriteSTD8(klon,nlevSTD8)
      real rhwriteSTD8(klon,nlevSTD8)
      real phiwriteSTD8(klon,nlevSTD8)
      real uwriteSTD8(klon,nlevSTD8)
      real vwriteSTD8(klon,nlevSTD8)
      real wwriteSTD8(klon,nlevSTD8)
c
c plevSTD3 END
c
c nout : niveau de output des variables a une pression donnee
      logical oknondef(klon,nlevSTD,nout)
c
c les produits uvSTD, vqSTD, .., T2STD sont calcules
c a partir des valeurs instantannees toutes les 6 h
c qui sont moyennees sur le mois
c
#include "radopt.h"
c
c
c prw: precipitable water
      real prw(klon)

      REAL convliq(klon,klev)  ! eau liquide nuageuse convective
      REAL convfra(klon,klev)  ! fraction nuageuse convective

      REAL cldl_c(klon),cldm_c(klon),cldh_c(klon) !nuages bas, moyen et haut
      REAL cldt_c(klon),cldq_c(klon) !nuage total, eau liquide integree
      REAL cldl_s(klon),cldm_s(klon),cldh_s(klon) !nuages bas, moyen et haut
      REAL cldt_s(klon),cldq_s(klon) !nuage total, eau liquide integree

      INTEGER linv, kp1
c flwp, fiwp = Liquid Water Path & Ice Water Path (kg/m2)
c flwc, fiwc = Liquid Water Content & Ice Water Content (kg/kg)
      REAL flwp(klon), fiwp(klon)
      REAL flwc(klon,klev), fiwc(klon,klev)
      REAL flwp_c(klon), fiwp_c(klon)
      REAL flwc_c(klon,klev), fiwc_c(klon,klev)
      REAL flwp_s(klon), fiwp_s(klon)
      REAL flwc_s(klon,klev), fiwc_s(klon,klev)

cIM ISCCP simulator v3.4
c dans clesphys.h top_height, overlap
cv3.4
      INTEGER debug, debugcol
cym      INTEGER npoints
cym      PARAMETER(npoints=klon) 
c
      INTEGER sunlit(klon) !sunlit=1 if day; sunlit=0 if night
      INTEGER nregISCtot
      PARAMETER(nregISCtot=1) 
c
c imin_debut, nbpti, jmin_debut, nbptj : parametres pour sorties sur 1 region rectangulaire
c y compris pour 1 point
c imin_debut : indice minimum de i; nbpti : nombre de points en direction i (longitude)
c jmin_debut : indice minimum de j; nbptj : nombre de points en direction j (latitude)
      INTEGER imin_debut, nbpti
      INTEGER jmin_debut, nbptj 
cIM parametres ISCCP BEG
      INTEGER nbapp_isccp
!     INTEGER nbapp_isccp,isccppas
!     PARAMETER(isccppas=6) !appel du simulateurs tous les 6pas de temps de la physique
!                           !i.e. toutes les 3 heures 
      INTEGER n
      INTEGER ifreq_isccp(napisccp), freqin_pdt(napisccp)
      DATA ifreq_isccp/3/
      SAVE ifreq_isccp
c$OMP THREADPRIVATE(ifreq_isccp)
      CHARACTER*5 typinout(napisccp)
      DATA typinout/'i3od'/
      SAVE typinout
c$OMP THREADPRIVATE(typinout)
cIM verif boxptop BEG
      CHARACTER*1 verticaxe(napisccp)
      DATA verticaxe/'1'/ 
      SAVE verticaxe
c$OMP THREADPRIVATE(verticaxe)
cIM verif boxptop END
      INTEGER nvlev(napisccp)
c     INTEGER nvlev
      REAL t1, aa
      REAL seed_re(klon,napisccp)
cym !!!! A voir plus tard 
cym      INTEGER iphy(iim,jjmp1)
cIM parametres ISCCP END
c
c ncol = nb. de sous-colonnes pour chaque maille du GCM 
c ncolmx = No. max. de sous-colonnes pour chaque maille du GCM 
c      INTEGER ncol(napisccp), ncolmx, seed(klon,napisccp)
      INTEGER,SAVE :: ncol(napisccp)
c$OMP THREADPRIVATE(ncol)
      INTEGER ncolmx, seed(klon,napisccp)
      REAL nbsunlit(nregISCtot,klon,napisccp)  !nbsunlit : moyenne de sunlit
c     PARAMETER(ncolmx=1500)
      PARAMETER(ncolmx=300)
c
cIM verif boxptop BEG
      REAL vertlev(ncolmx,napisccp)
cIM verif boxptop END
c
      REAL,SAVE :: tautab_omp(0:255),tautab(0:255)
      INTEGER,SAVE :: invtau_omp(-20:45000),invtau(-20:45000)
c$OMP THREADPRIVATE(tautab,invtau)
      REAL emsfc_lw
      PARAMETER(emsfc_lw=0.99)
c     REAL    ran0                      ! type for random number fuction
c
      REAL cldtot(klon,klev)
c variables de haut en bas pour le simulateur ISCCP
      REAL dtau_s(klon,klev) !tau nuages startiformes
      REAL dtau_c(klon,klev) !tau nuages convectifs
      REAL dem_s(klon,klev)  !emissivite nuages startiformes 
      REAL dem_c(klon,klev)  !emissivite nuages convectifs
c
c variables de haut en bas pour le simulateur ISCCP
      REAL pfull(klon,klev)
      REAL phalf(klon,klev+1)
      REAL qv(klon,klev)
      REAL cc(klon,klev)
      REAL conv(klon,klev)
      REAL dtau_sH2B(klon,klev)
      REAL dtau_cH2B(klon,klev)
      REAL at(klon,klev)
      REAL dem_sH2B(klon,klev)
      REAL dem_cH2B(klon,klev)
c
      INTEGER kmax, lmax, lmax3
      PARAMETER(kmax=8, lmax=8, lmax3=3)
      INTEGER kmaxm1, lmaxm1
      PARAMETER(kmaxm1=kmax-1, lmaxm1=lmax-1)
      INTEGER iimx7, jjmx7, jjmp1x7
      PARAMETER(iimx7=iim*kmaxm1, jjmx7=jjm*lmaxm1, 
     .jjmp1x7=jjmp1*lmaxm1)
c
c output from ISCCP simulator
      REAL fq_isccp(klon,kmaxm1,lmaxm1,napisccp)
      REAL fq_is_true(klon,kmaxm1,lmaxm1,napisccp)
      REAL totalcldarea(klon,napisccp) 
      REAL meanptop(klon,napisccp)
      REAL meantaucld(klon,napisccp)
      REAL boxtau(klon,ncolmx,napisccp)
      REAL boxptop(klon,ncolmx,napisccp) 
      REAL zx_tmp_fi3d_bx(klon,ncolmx)
      REAL zx_tmp_3d_bx(iim,jjmp1,ncolmx)
c
      REAL cld_fi3d(klon,lmax3)
      REAL cld_3d(iim,jjmp1,lmax3)
c
      INTEGER iw, iwmax
      REAL wmin, pas_w
c     PARAMETER(wmin=-100.,pas_w=10.,iwmax=30)
cIM 051005     PARAMETER(wmin=-200.,pas_w=10.,iwmax=40)
      PARAMETER(wmin=-100.,pas_w=10.,iwmax=20)
      REAL o500(klon)
c

c sorties ISCCP

      integer nid_isccp
      save nid_isccp        
c$OMP THREADPRIVATE(nid_isccp)

      REAL zx_tau(kmaxm1), zx_pc(lmaxm1), zx_o500(iwmax)
      DATA zx_tau/0.0, 0.3, 1.3, 3.6, 9.4, 23., 60./
      SAVE zx_tau
      DATA zx_pc/180., 310., 440., 560., 680., 800., 1000./
      SAVE zx_pc
c$OMP THREADPRIVATE(zx_tau,zx_pc)
c cldtopres pression au sommet des nuages
      REAL cldtopres(lmaxm1), cldtopres3(lmax3)
      DATA cldtopres/180., 310., 440., 560., 680., 800., 1000./
      DATA cldtopres3/440., 680., 1000./
      SAVE cldtopres,cldtopres3
c$OMP THREADPRIVATE(cldtopres,cldtopres3)
cIM 051005 BEG
      INTEGER komega, nhoriRD 

c taulev: numero du niveau de tau dans les sorties ISCCP
      CHARACTER *4 taulev(kmaxm1)
c     DATA taulev/'tau1','tau2','tau3','tau4','tau5','tau6','tau7'/
      DATA taulev/'tau0','tau1','tau2','tau3','tau4','tau5','tau6'/
      CHARACTER *3 pclev(lmaxm1)
      DATA pclev/'pc1','pc2','pc3','pc4','pc5','pc6','pc7'/
      SAVE taulev,pclev
c$OMP THREADPRIVATE(taulev,pclev)
c
c cnameisccp
      CHARACTER *29 cnameisccp(lmaxm1,kmaxm1)
cIM bad 151205     DATA cnameisccp/'pc< 50hPa, tau< 0.3', 
      DATA cnameisccp/'pc= 50-180hPa, tau< 0.3',
     .                'pc= 180-310hPa, tau< 0.3',
     .                'pc= 310-440hPa, tau< 0.3',
     .                'pc= 440-560hPa, tau< 0.3',
     .                'pc= 560-680hPa, tau< 0.3',
     .                'pc= 680-800hPa, tau< 0.3',
     .                'pc= 800-1000hPa, tau< 0.3',
     .                'pc= 50-180hPa, tau= 0.3-1.3',
     .                'pc= 180-310hPa, tau= 0.3-1.3',
     .                'pc= 310-440hPa, tau= 0.3-1.3',
     .                'pc= 440-560hPa, tau= 0.3-1.3',
     .                'pc= 560-680hPa, tau= 0.3-1.3',
     .                'pc= 680-800hPa, tau= 0.3-1.3',
     .                'pc= 800-1000hPa, tau= 0.3-1.3',
     .                'pc= 50-180hPa, tau= 1.3-3.6',
     .                'pc= 180-310hPa, tau= 1.3-3.6',
     .                'pc= 310-440hPa, tau= 1.3-3.6',
     .                'pc= 440-560hPa, tau= 1.3-3.6',
     .                'pc= 560-680hPa, tau= 1.3-3.6',
     .                'pc= 680-800hPa, tau= 1.3-3.6',
     .                'pc= 800-1000hPa, tau= 1.3-3.6',
     .                'pc= 50-180hPa, tau= 3.6-9.4',
     .                'pc= 180-310hPa, tau= 3.6-9.4',
     .                'pc= 310-440hPa, tau= 3.6-9.4',
     .                'pc= 440-560hPa, tau= 3.6-9.4',
     .                'pc= 560-680hPa, tau= 3.6-9.4',
     .                'pc= 680-800hPa, tau= 3.6-9.4',
     .                'pc= 800-1000hPa, tau= 3.6-9.4',
     .                'pc= 50-180hPa, tau= 9.4-23',
     .                'pc= 180-310hPa, tau= 9.4-23',
     .                'pc= 310-440hPa, tau= 9.4-23',
     .                'pc= 440-560hPa, tau= 9.4-23',
     .                'pc= 560-680hPa, tau= 9.4-23',
     .                'pc= 680-800hPa, tau= 9.4-23',
     .                'pc= 800-1000hPa, tau= 9.4-23',
     .                'pc= 50-180hPa, tau= 23-60',
     .                'pc= 180-310hPa, tau= 23-60',
     .                'pc= 310-440hPa, tau= 23-60',
     .                'pc= 440-560hPa, tau= 23-60',
     .                'pc= 560-680hPa, tau= 23-60',
     .                'pc= 680-800hPa, tau= 23-60',
     .                'pc= 800-1000hPa, tau= 23-60',
     .                'pc= 50-180hPa, tau> 60.',
     .                'pc= 180-310hPa, tau> 60.',
     .                'pc= 310-440hPa, tau> 60.',
     .                'pc= 440-560hPa, tau> 60.',
     .                'pc= 560-680hPa, tau> 60.',
     .                'pc= 680-800hPa, tau> 60.',
     .                'pc= 800-1000hPa, tau> 60.'/
       SAVE cnameisccp
c$OMP THREADPRIVATE(cnameisccp)
c
c     REAL zx_lonx7(iimx7), zx_latx7(jjmp1x7)
c     INTEGER nhorix7
cIM: region='3d' <==> sorties en global
      CHARACTER*3 region
      PARAMETER(region='3d')
c
cIM ISCCP simulator v3.4
c
      logical ok_hf
c
      integer nid_hf, nid_hf3d
      save ok_hf, nid_hf, nid_hf3d
c$OMP THREADPRIVATE(ok_hf, nid_hf, nid_hf3d)
c  QUESTION : noms de variables ?

      INTEGER        longcles
      PARAMETER    ( longcles = 20 )
      REAL clesphy0( longcles      )
c
c Variables propres a la physique
      INTEGER itap
      SAVE itap                   ! compteur pour la physique
c$OMP THREADPRIVATE(itap)
c
      real slp(klon) ! sea level pressure
c
      REAL fevap(klon,nbsrf)
      REAL fluxlat(klon,nbsrf)
c
      REAL qsol(klon)
      REAL,save ::  solarlong0
c$OMP THREADPRIVATE(solarlong0)

c
c  Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
c
cIM 141004     REAL zulow(klon),zvlow(klon),zustr(klon), zvstr(klon)
      REAL zulow(klon),zvlow(klon)
c
      INTEGER igwd,idx(klon),itest(klon)
c
      REAL agesno(klon,nbsrf)
c
c      REAL,allocatable,save :: run_off_lic_0(:)
cc$OMP THREADPRIVATE(run_off_lic_0)
cym      SAVE run_off_lic_0
cKE43
c Variables liees a la convection de K. Emanuel (sb):
c
      REAL bas, top             ! cloud base and top levels
      SAVE bas
      SAVE top
c$OMP THREADPRIVATE(bas, top)

      REAL wdn(klon), tdn(klon), qdn(klon)
c
c=================================================================================================
cCR04.12.07: on ajoute les nouvelles variables du nouveau schema de convection avec poches froides
c Variables liées à la poche froide (jyg)

      REAL mip(klon,klev)  ! mass flux shed by the adiab ascent at each level
      REAL Vprecip(klon,klev+1)   ! precipitation vertical profile
c
      REAL wape_prescr, fip_prescr
      INTEGER it_wape_prescr
      SAVE wape_prescr, fip_prescr, it_wape_prescr
c$OMP THREADPRIVATE(wape_prescr, fip_prescr, it_wape_prescr)
c
c variables supplementaires de concvl
      REAL Tconv(klon,klev)
      REAL ment(klon,klev,klev),sij(klon,klev,klev)
      REAL dd_t(klon,klev),dd_q(klon,klev)

      real, save :: alp_bl_prescr=0.
      real, save :: ale_bl_prescr=0.

      real, save :: ale_max=1000.
      real, save :: alp_max=2.

c$OMP THREADPRIVATE(alp_bl_prescr,ale_bl_prescr)
c$OMP THREADPRIVATE(ale_max,alp_max)

      real ale_wake(klon)
      real alp_wake(klon)
cRC
c Variables liées à la poche froide (jyg et rr)
c Version diagnostique pour l'instant : pas de rétroaction sur la convection

      REAL t_wake(klon,klev),q_wake(klon,klev) ! wake pour la convection

      REAL wake_dth(klon,klev)        ! wake : temp pot difference

      REAL wake_d_deltat_gw(klon,klev)! wake : delta T tendency due to Gravity Wave (/s)
      REAL wake_omgbdth(klon,klev)    ! Wake : flux of Delta_Theta transported by LS omega
      REAL wake_dp_omgb(klon,klev)    ! Wake : vertical gradient of large scale omega
      REAL wake_dtKE(klon,klev)       ! Wake : differential heating (wake - unpertubed) CONV
      REAL wake_dqKE(klon,klev)       ! Wake : differential moistening (wake - unpertubed) CONV
      REAL wake_dtPBL(klon,klev)      ! Wake : differential heating (wake - unpertubed) PBL
      REAL wake_dqPBL(klon,klev)      ! Wake : differential moistening (wake - unpertubed) PBL
      REAL wake_omg(klon,klev)        ! Wake : velocity difference (wake - unpertubed)
      REAL wake_ddeltat(klon,klev),wake_ddeltaq(klon,klev)
      REAL wake_dp_deltomg(klon,klev) ! Wake : gradient vertical de wake_omg
      REAL wake_spread(klon,klev)     ! spreading term in wake_delt
c
cpourquoi y'a pas de save??
      REAL wake_h(klon)               ! Wake : hauteur de la poche froide
c
      INTEGER wake_k(klon)            ! Wake sommet
c
      REAL t_undi(klon,klev)               ! temperature moyenne dans la zone non perturbee
      REAL q_undi(klon,klev)               ! humidite moyenne dans la zone non perturbee
c
cjyg
ccc      REAL wake_pe(klon)              ! Wake potential energy - WAPE 

      REAL wake_gfl(klon)             ! Gust Front Length
      REAL wake_dens(klon)
c
c
      REAL dt_dwn(klon,klev)
      REAL dq_dwn(klon,klev)
      REAL wdt_PBL(klon,klev)
      REAL udt_PBL(klon,klev)
      REAL wdq_PBL(klon,klev)
      REAL udq_PBL(klon,klev)
      REAL M_dwn(klon,klev)
      REAL M_up(klon,klev)
      REAL dt_a(klon,klev)
      REAL dq_a(klon,klev)
      REAL, SAVE :: alp_offset
c$OMP THREADPRIVATE(alp_offset)

c
cRR:fin declarations poches froides
c=======================================================================================================

      REAL zw2(klon,klev+1)
      REAL fraca(klon,klev+1)        
      REAL ztv(klon,klev) 
      REAL zpspsk(klon,klev)
      REAL ztla(klon,klev) 
      REAL zthl(klon,klev)

c Variables locales pour la couche limite (al1):
c
cAl1      REAL pblh(klon)           ! Hauteur de couche limite
cAl1      SAVE pblh
c34EK
c
c Variables locales:
c
      REAL cdragh(klon) ! drag coefficient pour T and Q
      REAL cdragm(klon) ! drag coefficient pour vent
cAA
cAA  Pour phytrac 
cAA
      REAL coefh(klon,klev)     ! coef d'echange pour phytrac, valable pour 2<=k<=klev
      REAL u1(klon)             ! vents dans la premiere couche U
      REAL v1(klon)             ! vents dans la premiere couche V

      REAL zxffonte(klon), zxfqcalving(klon),zxfqfonte(klon)

c@$$      LOGICAL offline           ! Controle du stockage ds "physique"
c@$$      PARAMETER (offline=.false.)
c@$$      INTEGER physid
      REAL frac_impa(klon,klev) ! fractions d'aerosols lessivees (impaction)
      REAL frac_nucl(klon,klev) ! idem (nucleation)
      INTEGER       :: iii
      REAL          :: calday

cIM cf FH pour Tiedtke 080604
      REAL rain_tiedtke(klon),snow_tiedtke(klon)
c
cIM 050204 END
      REAL evap(klon), devap(klon) ! evaporation et sa derivee
      REAL sens(klon), dsens(klon) ! chaleur sensible et sa derivee

      REAL bils(klon) ! bilan de chaleur au sol
      REAL wfbilo(klon,nbsrf) ! bilan d'eau, pour chaque
C                             ! type de sous-surface et pondere par la fraction
      REAL wfbils(klon,nbsrf) ! bilan de chaleur au sol, pour chaque
C                             ! type de sous-surface et pondere par la fraction
      REAL slab_wfbils(klon)  ! bilan de chaleur au sol pour le cas de slab, sur les points d'ocean

      REAL fder(klon)         
      REAL ve(klon) ! integr. verticale du transport meri. de l'energie
      REAL vq(klon) ! integr. verticale du transport meri. de l'eau
      REAL ue(klon) ! integr. verticale du transport zonal de l'energie
      REAL uq(klon) ! integr. verticale du transport zonal de l'eau
c
      REAL frugs(klon,nbsrf)
      REAL zxrugs(klon) ! longueur de rugosite
c
c Conditions aux limites
c
!
      REAL :: day_since_equinox
! Date de l'equinoxe de printemps
      INTEGER, parameter :: mth_eq=3, day_eq=21
      REAL :: jD_eq

      LOGICAL, parameter :: new_orbit = .true.

c
      INTEGER lmt_pas
      SAVE lmt_pas                ! frequence de mise a jour
c$OMP THREADPRIVATE(lmt_pas) 
      real zmasse(klon, llm) 
C     (column-density of mass of air in a cell, in kg m-2)
      real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

cIM sorties
      REAL un_jour
      PARAMETER(un_jour=86400.)
c======================================================================
c
c Declaration des procedures appelees
c
      EXTERNAL angle     ! calculer angle zenithal du soleil
      EXTERNAL alboc     ! calculer l'albedo sur ocean
      EXTERNAL ajsec     ! ajustement sec
      EXTERNAL conlmd    ! convection (schema LMD)
cKE43
      EXTERNAL conema3  ! convect4.3
      EXTERNAL fisrtilp  ! schema de condensation a grande echelle (pluie)
cAA 
      EXTERNAL fisrtilp_tr ! schema de condensation a grande echelle (pluie)
c                          ! stockage des coefficients necessaires au
c                          ! lessivage OFF-LINE et ON-LINE
      EXTERNAL hgardfou  ! verifier les temperatures
      EXTERNAL nuage     ! calculer les proprietes radiatives
CC      EXTERNAL o3cm      ! initialiser l'ozone
      EXTERNAL orbite    ! calculer l'orbite terrestre
      EXTERNAL phyetat0  ! lire l'etat initial de la physique
      EXTERNAL phyredem  ! ecrire l'etat de redemarrage de la physique
      EXTERNAL suphel    ! initialiser certaines constantes
      EXTERNAL transp    ! transport total de l'eau et de l'energie
      EXTERNAL ecribina  ! ecrire le fichier binaire global
      EXTERNAL ecribins  ! ecrire le fichier binaire global
      EXTERNAL ecrirega  ! ecrire le fichier binaire regional
      EXTERNAL ecriregs  ! ecrire le fichier binaire regional
cIM
      EXTERNAL haut2bas  !variables de haut en bas
      INTEGER lnblnk1
      EXTERNAL lnblnk1   !enleve les blancs a la fin d'une variable de type
                         !caracter
      EXTERNAL ini_undefSTD  !initialise a 0 une variable a 1 niveau de pression
      EXTERNAL undefSTD      !somme les valeurs definies d'1 var a 1 niveau de pression
c     EXTERNAL moy_undefSTD  !moyenne d'1 var a 1 niveau de pression
c     EXTERNAL moyglo_aire   !moyenne globale d'1 var ponderee par l'aire de la maille (moyglo_pondaire)
c                            !par la masse/airetot (moyglo_pondaima) et la vraie masse (moyglo_pondmass)
c
c Variables locales
c
      REAL rhcl(klon,klev)    ! humiditi relative ciel clair
      REAL dialiq(klon,klev)  ! eau liquide nuageuse
      REAL diafra(klon,klev)  ! fraction nuageuse
      REAL cldliq(klon,klev)  ! eau liquide nuageuse
      REAL cldfra(klon,klev)  ! fraction nuageuse
      REAL cldtau(klon,klev)  ! epaisseur optique
      REAL cldemi(klon,klev)  ! emissivite infrarouge
c
CXXX PB 
      REAL fluxq(klon,klev, nbsrf)   ! flux turbulent d'humidite
      REAL fluxt(klon,klev, nbsrf)   ! flux turbulent de chaleur
      REAL fluxu(klon,klev, nbsrf)   ! flux turbulent de vitesse u
      REAL fluxv(klon,klev, nbsrf)   ! flux turbulent de vitesse v
c
      REAL zxfluxt(klon, klev)
      REAL zxfluxq(klon, klev)
      REAL zxfluxu(klon, klev)
      REAL zxfluxv(klon, klev)
CXXX
c
      REAL fsollw(klon, nbsrf)   ! bilan flux IR pour chaque sous surface
      REAL fsolsw(klon, nbsrf)   ! flux solaire absorb. pour chaque sous surface
c Le rayonnement n'est pas calcule tous les pas, il faut donc
c                      sauvegarder les sorties du rayonnement
cym      SAVE  heat,cool,albpla,topsw,toplw,solsw,sollw,sollwdown
cym      SAVE  sollwdownclr, toplwdown, toplwdownclr
cym      SAVE  topsw0,toplw0,solsw0,sollw0, heat0, cool0
c
      INTEGER itaprad
      SAVE itaprad
c$OMP THREADPRIVATE(itaprad)
c
      REAL conv_q(klon,klev) ! convergence de l'humidite (kg/kg/s)
      REAL conv_t(klon,klev) ! convergence de la temperature(K/s)
c
      REAL cldl(klon),cldm(klon),cldh(klon) !nuages bas, moyen et haut
      REAL cldt(klon),cldq(klon) !nuage total, eau liquide integree
c
      REAL zxtsol(klon), zxqsurf(klon), zxsnow(klon), zxfluxlat(klon)
      REAL zxsnow_dummy(klon)
c
      REAL dist, rmu0(klon), fract(klon)
      REAL zdtime, zlongi
c
      CHARACTER*2 str2
      CHARACTER*2 iqn
c
      REAL qcheck
      REAL z_avant(klon), z_apres(klon), z_factor(klon)
      LOGICAL zx_ajustq
c
      REAL za, zb
      REAL zx_t, zx_qs, zdelta, zcor, zfra, zlvdcp, zlsdcp
      real zqsat(klon,klev)
      INTEGER i, k, iq, ig, j, nsrf, ll, l, iiq, iff
      REAL t_coup
      PARAMETER (t_coup=234.0)
c
      REAL zphi(klon,klev)
cym A voir plus tard !!
cym      REAL zx_relief(iim,jjmp1)
cym      REAL zx_aire(iim,jjmp1)
c
c Grandeurs de sorties
      REAL s_pblh(klon), s_lcl(klon), s_capCL(klon)
      REAL s_oliqCL(klon), s_cteiCL(klon), s_pblt(klon)
      REAL s_therm(klon), s_trmb1(klon), s_trmb2(klon)
      REAL s_trmb3(klon)
cKE43
c Variables locales pour la convection de K. Emanuel (sb):
c
      REAL upwd(klon,klev)      ! saturated updraft mass flux
      REAL dnwd(klon,klev)      ! saturated downdraft mass flux
      REAL dnwd0(klon,klev)     ! unsaturated downdraft mass flux
      REAL tvp(klon,klev)       ! virtual temp of lifted parcel
      CHARACTER*40 capemaxcels  !max(CAPE)

      REAL rflag(klon)          ! flag fonctionnement de convect
      INTEGER iflagctrl(klon)          ! flag fonctionnement de convect

c -- convect43:
      INTEGER ntra              ! nb traceurs pour convect4.3
      REAL pori_con(klon)    ! pressure at the origin level of lifted parcel
      REAL plcl_con(klon),dtma_con(klon),dtlcl_con(klon)
      REAL dtvpdt1(klon,klev), dtvpdq1(klon,klev)
      REAL dplcldt(klon), dplcldr(klon)
c?     .     condm_con(klon,klev),conda_con(klon,klev),
c?     .     mr_con(klon,klev),ep_con(klon,klev)
c?     .    ,sadiab(klon,klev),wadiab(klon,klev)
c --
c34EK
c
c Variables du changement
c
c con: convection
c lsc: condensation a grande echelle (Large-Scale-Condensation)
c ajs: ajustement sec
c eva: evaporation de l'eau liquide nuageuse
c vdf: couche limite (Vertical DiFfusion)
      REAL rneb(klon,klev)

! tendance nulles
      REAL du0(klon,klev),dv0(klon,klev),dq0(klon,klev),dql0(klon,klev)

c
*********************************************************
*     declarations
      
*********************************************************
cIM 081204 END
c
      REAL pmfu(klon,klev), pmfd(klon,klev)
      REAL pen_u(klon,klev), pen_d(klon,klev)
      REAL pde_u(klon,klev), pde_d(klon,klev)
      INTEGER kcbot(klon), kctop(klon), kdtop(klon)
      REAL pmflxr(klon,klev+1), pmflxs(klon,klev+1)
      REAL prfl(klon,klev+1), psfl(klon,klev+1)
c
      REAL rain_lsc(klon)
      REAL snow_lsc(klon)
c
      REAL ratqss(klon,klev),ratqsc(klon,klev)
      real ratqsbas,ratqshaut,tau_ratqs
      save ratqsbas,ratqshaut,tau_ratqs
c$OMP THREADPRIVATE(ratqsbas,ratqshaut,tau_ratqs)
      real zpt_conv(klon,klev)

c Parametres lies au nouveau schema de nuages (SB, PDF)
      real fact_cldcon
      real facttemps
      logical ok_newmicro
      save ok_newmicro
      real ref_liq(klon,klev), ref_ice(klon,klev)
c$OMP THREADPRIVATE(ok_newmicro)
      save fact_cldcon,facttemps
c$OMP THREADPRIVATE(fact_cldcon,facttemps)
      real facteur

      integer iflag_cldcon
      save iflag_cldcon
c$OMP THREADPRIVATE(iflag_cldcon)
      logical ptconv(klon,klev)
cIM cf. AM 081204 BEG
      logical ptconvth(klon,klev)
cIM cf. AM 081204 END
c
c Variables liees a l'ecriture de la bande histoire physique
c
c======================================================================
c
cIM cf. AM 081204 BEG
c   declarations pour sortir sur une sous-region
      integer imin_ins,imax_ins,jmin_ins,jmax_ins
      save imin_ins,imax_ins,jmin_ins,jmax_ins
c$OMP THREADPRIVATE(imin_ins,imax_ins,jmin_ins,jmax_ins)
c      real lonmin_ins,lonmax_ins,latmin_ins
c     s  ,latmax_ins
c     data lonmin_ins,lonmax_ins,latmin_ins
c    s  ,latmax_ins/
c valeurs initiales     s   -5.,20.,41.,55./   
c    s   100.,130.,-20.,20./
c    s   -180.,180.,-90.,90./
c======================================================================
cIM cf. AM 081204 END

c
      integer itau_w   ! pas de temps ecriture = itap + itau_phy
c
c
c Variables locales pour effectuer les appels en serie
c
      REAL zx_rh(klon,klev)
cIM RH a 2m (la surface)
      REAL rh2m(klon), qsat2m(klon)
      REAL tpot(klon), tpote(klon)
      REAL Lheat

      INTEGER        length
      PARAMETER    ( length = 100 )
      REAL tabcntr0( length       )
c
      INTEGER ndex2d(iim*jjmp1),ndex3d(iim*jjmp1*klev)
cIM
      INTEGER ndex2d1(iwmax)
c
cIM AMIP2 BEG
      REAL moyglo, mountor
cIM 141004 BEG
      REAL zustrdr(klon), zvstrdr(klon)
      REAL zustrli(klon), zvstrli(klon)
      REAL zustrph(klon), zvstrph(klon)
      REAL zustrhi(klon), zvstrhi(klon)
      REAL aam, torsfc
cIM 141004 END
cIM 190504 BEG
      INTEGER ij, imp1jmp1
      PARAMETER(imp1jmp1=(iim+1)*jjmp1)
cym A voir plus tard
      REAL zx_tmp(imp1jmp1), airedyn(iim+1,jjmp1)
      REAL padyn(iim+1,jjmp1,klev+1)
      REAL dudyn(iim+1,jjmp1,klev)
      REAL rlatdyn(iim+1,jjmp1)
cIM 190504 END
      LOGICAL ok_msk
      REAL msk(klon)
cIM 
      REAL airetot, pi
cym A voir plus tard
cym      REAL zm_wo(jjmp1, klev)
cIM AMIP2 END
c
      REAL zx_tmp_fi2d(klon)      ! variable temporaire grille physique
      REAL zx_tmp_fi3d(klon,klev) ! variable temporaire pour champs 3D 
c#ifdef histNMC
cym   A voir plus tard !!!!
cym      REAL zx_tmp_NC(iim,jjmp1,nlevSTD)
      REAL zx_tmp_fiNC(klon,nlevSTD) 
c#endif
      REAL(KIND=8) zx_tmp2_fi3d(klon,klev) ! variable temporaire pour champs 3D 
      REAL zx_tmp_2d(iim,jjmp1), zx_tmp_3d(iim,jjmp1,klev)
      REAL zx_lon(iim,jjmp1), zx_lat(iim,jjmp1)
cIM for NMC files
      REAL missing_val
      REAL, SAVE :: freq_moyNMC(nout)
c$OMP THREADPRIVATE(freq_moyNMC)
c
      INTEGER nid_day, nid_mth, nid_ins, nid_mthnmc, nid_daynmc
      INTEGER nid_hfnmc, nid_day_seri, nid_ctesGCM
      SAVE nid_day, nid_mth, nid_ins, nid_mthnmc, nid_daynmc
      SAVE nid_hfnmc, nid_day_seri, nid_ctesGCM
c$OMP THREADPRIVATE(nid_day, nid_mth, nid_ins)
c$OMP THREADPRIVATE(nid_mthnmc, nid_daynmc, nid_hfnmc)
c$OMP THREADPRIVATE(nid_day_seri,nid_ctesGCM)
c
cIM 280405 BEG
      INTEGER nid_bilKPins, nid_bilKPave
      SAVE nid_bilKPins, nid_bilKPave
c$OMP THREADPRIVATE(nid_bilKPins, nid_bilKPave)
c
      REAL ve_lay(klon,klev) ! transport meri. de l'energie a chaque niveau vert.
      REAL vq_lay(klon,klev) ! transport meri. de l'eau a chaque niveau vert.
      REAL ue_lay(klon,klev) ! transport zonal de l'energie a chaque niveau vert.
      REAL uq_lay(klon,klev) ! transport zonal de l'eau a chaque niveau vert.
c
cIM 280405 END
c
      INTEGER nhori, nvert, nvert1, nvert3
      REAL zsto, zsto1, zsto2
      REAL zstophy, zstorad, zstohf, zstoday, zstomth, zout
      REAL zcals(napisccp), zcalh(napisccp), zoutj(napisccp)
      REAL zout_isccp(napisccp)
      SAVE zcals, zcalh, zoutj, zout_isccp
c$OMP THREADPRIVATE(zcals, zcalh, zoutj, zout_isccp)

      real zjulian
      save zjulian
c$OMP THREADPRIVATE(zjulian)

      character*20 modname
      character*80 abort_message
      logical ok_sync
      real date0
      integer idayref

C essai writephys
      integer fid_day, fid_mth, fid_ins
      parameter (fid_ins = 1, fid_day = 2, fid_mth = 3) 
      integer prof2d_on, prof3d_on, prof2d_av, prof3d_av
      parameter (prof2d_on = 1, prof3d_on = 2,
     .           prof2d_av = 3, prof3d_av = 4)
      character*30 nom_fichier
      character*10 varname
      character*40 vartitle
      character*20 varunits
C     Variables liees au bilan d'energie et d'enthalpi
      REAL ztsol(klon)
      REAL      h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot
     $        , h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot
      SAVE      h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot
     $        , h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot
c$OMP THREADPRIVATE(h_vcol_tot, h_dair_tot, h_qw_tot, h_ql_tot,
c$OMP+              h_qs_tot, qw_tot, ql_tot, qs_tot , ec_tot)
      REAL      d_h_vcol, d_h_dair, d_qt, d_qw, d_ql, d_qs, d_ec
      REAL      d_h_vcol_phy
      REAL      fs_bound, fq_bound
      SAVE      d_h_vcol_phy
c$OMP THREADPRIVATE(d_h_vcol_phy)
      REAL      zero_v(klon)
      CHARACTER*15 ztit
      INTEGER   ip_ebil  ! PRINT level for energy conserv. diag.
      SAVE      ip_ebil
      DATA      ip_ebil/0/
c$OMP THREADPRIVATE(ip_ebil)
      INTEGER   if_ebil ! level for energy conserv. dignostics
      SAVE      if_ebil
c$OMP THREADPRIVATE(if_ebil)
c+jld ec_conser
      REAL ZRCPD
c-jld ec_conser
      REAL t2m(klon,nbsrf)  ! temperature a 2m
      REAL q2m(klon,nbsrf)  ! humidite a 2m

cIM: t2m, q2m, u10m, v10m et t2mincels, t2maxcels
      REAL zt2m(klon), zq2m(klon)             !temp., hum. 2m moyenne s/ 1 maille
      REAL zu10m(klon), zv10m(klon)           !vents a 10m moyennes s/1 maille
      CHARACTER*40 t2mincels, t2maxcels       !t2m min., t2m max
      CHARACTER*40 tinst, tave, typeval
      REAL cldtaupi(klon,klev)  ! Cloud optical thickness for pre-industrial (pi) aerosols

      REAL re(klon, klev)       ! Cloud droplet effective radius
      REAL fl(klon, klev)  ! denominator of re

      REAL re_top(klon), fl_top(klon) ! CDR at top of liquid water clouds

      ! Aerosol optical properties
      CHARACTER*4, DIMENSION(naero_grp) :: rfname 
      REAL, DIMENSION(klon)          :: aerindex     ! POLDER aerosol index
      REAL, DIMENSION(klon,klev)     :: mass_solu_aero    ! total mass concentration for all soluble aerosols[ug/m3]
      REAL, DIMENSION(klon,klev)     :: mass_solu_aero_pi ! - " - (pre-industrial value)
      INTEGER :: naero ! aerosol species 

      ! Parameters
      LOGICAL ok_ade, ok_aie    ! Apply aerosol (in)direct effects or not
      REAL bl95_b0, bl95_b1   ! Parameter in Boucher and Lohmann (1995)
      SAVE ok_ade, ok_aie, bl95_b0, bl95_b1
c$OMP THREADPRIVATE(ok_ade, ok_aie, bl95_b0, bl95_b1)
      LOGICAL, SAVE :: aerosol_couple ! true  : calcul des aerosols dans INCA
                                      ! false : lecture des aerosol dans un fichier
c$OMP THREADPRIVATE(aerosol_couple)    
      INTEGER, SAVE :: flag_aerosol 
c$OMP THREADPRIVATE(flag_aerosol) 
      LOGICAL, SAVE :: new_aod
c$OMP THREADPRIVATE(new_aod) 
   
c
c Declaration des constantes et des fonctions thermodynamiques
c
      LOGICAL,SAVE :: first=.true.
c$OMP THREADPRIVATE(first)

      integer iunit

      integer, save::  read_climoz ! read ozone climatology
C     (let it keep the default OpenMP shared attribute)
C     Allowed values are 0, 1 and 2
C     0: do not read an ozone climatology
C     1: read a single ozone climatology that will be used day and night
C     2: read two ozone climatologies, the average day and night
C     climatology and the daylight climatology

      integer, save:: ncid_climoz ! NetCDF file containing ozone climatologies
C     (let it keep the default OpenMP shared attribute)

      real, pointer, save:: press_climoz(:)
C     (let it keep the default OpenMP shared attribute)
!     edges of pressure intervals for ozone climatologies, in Pa, in strictly
!     ascending order

      integer, save:: co3i = 0
!     time index in NetCDF file of current ozone fields
c$OMP THREADPRIVATE(co3i) 

      integer ro3i
!     required time index in NetCDF file for the ozone fields, between 1
!     and 360

      INTEGER ierr
#include "YOMCST.h"
#include "YOETHF.h"
#include "FCTTRE.h"
cIM 100106 BEG : pouvoir sortir les ctes de la physique
#include "conema3.h"
#include "fisrtilp.h"
#include "nuage.h"
#include "compbl.h"
cIM 100106 END : pouvoir sortir les ctes de la physique
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Declarations pour Simulateur COSP
c============================================================
      real :: mr_ozone(klon,klev)
cIM for NMC files
      missing_val=nf90_fill_real
c======================================================================
! Gestion calendrier : mise a jour du module phys_cal_mod
!
      CALL phys_cal_update(jD_cur,jH_cur)

c======================================================================
! Ecriture eventuelle d'un profil verticale en entree de la physique.
! Utilise notamment en 1D mais peut etre active egalement en 3D
! en imposant la valeur de igout.
c======================================================================

      if (prt_level.ge.1) then
          igout=klon/2+1/klon
         write(lunout,*) 'DEBUT DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
         write(lunout,*)
     s 'nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys'
         write(lunout,*)
     s  nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys 

         write(lunout,*) 'paprs, play, phi, u, v, t'
         do k=1,klev
            write(lunout,*) paprs(igout,k),pplay(igout,k),pphi(igout,k),
     s   u(igout,k),v(igout,k),t(igout,k)
         enddo
         write(lunout,*) 'ovap (g/kg),  oliq (g/kg)'
         do k=1,klev
            write(lunout,*) qx(igout,k,1)*1000,qx(igout,k,2)*1000.
         enddo
      endif

c======================================================================

cym => necessaire pour iflag_con != 2    
      pmfd(:,:) = 0.
      pen_u(:,:) = 0.
      pen_d(:,:) = 0.
      pde_d(:,:) = 0.
      pde_u(:,:) = 0.
      aam=0.

      torsfc=0.
      forall (k=1: llm) zmasse(:, k) = (paprs(:, k)-paprs(:, k+1)) / rg

      if (first) then 
      
cCR:nvelles variables convection/poches froides
      
      print*, '================================================='
      print*, 'Allocation des variables locales et sauvegardees'
      call phys_local_var_init
c
      pasphys=pdtphys
c     appel a la lecture du run.def physique
      call conf_phys(ok_journe, ok_mensuel,
     .     ok_instan, ok_hf,
     .     ok_LES,
     .     callstats,
     .     solarlong0,seuil_inversion,
     .     fact_cldcon, facttemps,ok_newmicro,iflag_radia,
     .     iflag_cldcon,iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs,
     .     ok_ade, ok_aie, aerosol_couple, 
     .     flag_aerosol, new_aod,
     .     bl95_b0, bl95_b1,
     .     iflag_thermals,nsplit_thermals,tau_thermals,
     .     iflag_thermals_ed,iflag_thermals_optflux,
c     nv flags pour la convection et les poches froides
     .     iflag_coupl,iflag_clos,iflag_wake, read_climoz,
     &     alp_offset)
      call phys_state_var_init(read_climoz)
      call phys_output_var_init
      print*, '================================================='
cIM for NMC files
cIM freq_moyNMC = frequences auxquelles on moyenne les champs accumules 
cIM               sur les niveaux de pression standard du NMC
      DO n=1, nout
       freq_moyNMC(n)=freq_outNMC(n)/freq_calNMC(n)
      ENDDO
c
cIM beg
          dnwd0=0.0
          ftd=0.0
          fqd=0.0
          cin=0.
cym Attention pbase pas initialise dans concvl !!!!
          pbase=0
cIM 180608
c         pmflxr=0.
c         pmflxs=0.

        itau_con=0
        first=.false.

      endif  ! first

       modname = 'physiq'
cIM
      IF (ip_ebil_phy.ge.1) THEN
        DO i=1,klon
          zero_v(i)=0.
        END DO 
      END IF 
      ok_sync=.TRUE.

      IF (debut) THEN
         CALL suphel ! initialiser constantes et parametres phys.
      ENDIF

      if(prt_level.ge.1) print*,'CONVERGENCE PHYSIQUE THERM 1 '


c======================================================================
! Gestion calendrier : mise a jour du module phys_cal_mod
!
c     CALL phys_cal_update(jD_cur,jH_cur)

c
c Si c'est le debut, il faut initialiser plusieurs choses
c          ********
c
       IF (debut) THEN
!rv
cCRinitialisation de wght_th et lalim_conv pour la definition de la couche alimentation 
cde la convection a partir des caracteristiques du thermique
         wght_th(:,:)=1.
         lalim_conv(:)=1 
cRC
         u10m(:,:)=0.
         v10m(:,:)=0.
         rain_con(:)=0.
         snow_con(:)=0.
         topswai(:)=0.
         topswad(:)=0.
         solswai(:)=0.
         solswad(:)=0.

         lambda_th(:,:)=0.
         wmax_th(:)=0.
         tau_overturning_th(:)=0.

         IF (config_inca /= 'none') THEN
            ! jg : initialisation jusqu'au ces variables sont dans restart
            ccm(:,:,:) = 0.
            tau_aero(:,:,:,:) = 0.
            piz_aero(:,:,:,:) = 0.
            cg_aero(:,:,:,:) = 0.
         END IF

         rnebcon0(:,:) = 0.0
         clwcon0(:,:) = 0.0
         rnebcon(:,:) = 0.0
         clwcon(:,:) = 0.0

cIM      
         IF (ip_ebil_phy.ge.1) d_h_vcol_phy=0.
c
      print*,'iflag_coupl,iflag_clos,iflag_wake',
     .   iflag_coupl,iflag_clos,iflag_wake
      print*,'CYCLE_DIURNE', cycle_diurne
c
      IF (iflag_con.EQ.2.AND.iflag_cldcon.GT.-1) THEN
         abort_message = 'Tiedtke needs iflag_cldcon=-2 or -1'
         CALL abort_gcm (modname,abort_message,1)
      ENDIF
c
      IF(ok_isccp.AND.iflag_con.LE.2) THEN
         abort_message = 'ISCCP-like outputs may be available for KE
     .(iflag_con >= 3); for Tiedtke (iflag_con=-2) put ok_isccp=n'
         CALL abort_gcm (modname,abort_message,1)
      ENDIF
c
c Initialiser les compteurs:
c
         itap    = 0
         itaprad = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Un petit travail à faire ici.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (iflag_pbl>1) then
             PRINT*, "Using method MELLOR&YAMADA" 
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/02 changement lie a la lecture de nbapp_rad dans phylmd plutot que
! dyn3d
! Attention : la version precedente n'etait pas tres propre.
! Il se peut qu'il faille prendre une valeur differente de nbapp_rad
! pour obtenir le meme resultat.
         dtime=pdtphys
         radpas = NINT( 86400./dtime/nbapp_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CALL phyetat0 ("startphy.nc",clesphy0,tabcntr0)
cIM begin
          print*,'physiq: clwcon rnebcon ratqs',clwcon(1,1),rnebcon(1,1)
     $,ratqs(1,1)
cIM end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
C on remet le calendrier a zero
c
         IF (raz_date .eq. 1) THEN
           itau_phy = 0
         ENDIF

cIM cf. AM 081204 BEG
         PRINT*,'cycle_diurne3 =',cycle_diurne
cIM cf. AM 081204 END
c
         CALL printflag( tabcntr0,radpas,ok_journe,
     ,                    ok_instan, ok_region )
c
         IF (ABS(dtime-pdtphys).GT.0.001) THEN
            WRITE(lunout,*) 'Pas physique n est pas correct',dtime,
     .                        pdtphys
            abort_message='Pas physique n est pas correct '
!           call abort_gcm(modname,abort_message,1)
            dtime=pdtphys
         ENDIF
         IF (nlon .NE. klon) THEN
            WRITE(lunout,*)'nlon et klon ne sont pas coherents', nlon, 
     .                      klon
            abort_message='nlon et klon ne sont pas coherents'
            call abort_gcm(modname,abort_message,1)
         ENDIF
         IF (nlev .NE. klev) THEN
            WRITE(lunout,*)'nlev et klev ne sont pas coherents', nlev,
     .                       klev
            abort_message='nlev et klev ne sont pas coherents'
            call abort_gcm(modname,abort_message,1)
         ENDIF
c
         IF (dtime*REAL(radpas).GT.21600..AND.cycle_diurne) THEN 
           WRITE(lunout,*)'Nbre d appels au rayonnement insuffisant'
           WRITE(lunout,*)"Au minimum 4 appels par jour si cycle diurne"
           abort_message='Nbre d appels au rayonnement insuffisant'
           call abort_gcm(modname,abort_message,1)
         ENDIF
         WRITE(lunout,*)"Clef pour la convection, iflag_con=", iflag_con
         WRITE(lunout,*)"Clef pour le driver de la convection, ok_cvl=",
     .                   ok_cvl
c
cKE43
c Initialisation pour la convection de K.E. (sb):
         IF (iflag_con.GE.3) THEN

         WRITE(lunout,*)"*** Convection de Kerry Emanuel 4.3  "
         WRITE(lunout,*)
     .      "On va utiliser le melange convectif des traceurs qui"
         WRITE(lunout,*)"est calcule dans convect4.3"
         WRITE(lunout,*)" !!! penser aux logical flags de phytrac"

          DO i = 1, klon
           ema_cbmf(i) = 0.
           ema_pcb(i)  = 0.
           ema_pct(i)  = 0.
c          ema_workcbmf(i) = 0.
          ENDDO
cIM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>BEG
          DO i = 1, klon
           ibas_con(i) = 1
           itop_con(i) = 1
          ENDDO
cIM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>END
c===============================================================================
cCR:04.12.07: initialisations poches froides
c Controle de ALE et ALP pour la fermeture convective (jyg)
          if (iflag_wake.eq.1) then
            CALL ini_wake(0.,0.,it_wape_prescr,wape_prescr,fip_prescr
     s                  ,alp_bl_prescr, ale_bl_prescr)
c 11/09/06 rajout initialisation ALE et ALP du wake et PBL(YU)
c        print*,'apres ini_wake iflag_cldcon=', iflag_cldcon
          endif

        do i = 1,klon
         Ale_bl(i)=0.
         Alp_bl(i)=0.
        enddo

c================================================================================

         ENDIF !debut

           DO i=1,klon
             rugoro(i) = f_rugoro * MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
           ENDDO

c34EK
         IF (ok_orodr) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH sans doute a enlever de finitivement ou, si on le garde, l'activer
! justement quand ok_orodr = false.
! ce rugoro est utilise par la couche limite et fait double emploi
! avec les paramétrisations spécifiques de Francois Lott.
!           DO i=1,klon
!             rugoro(i) = MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
!           ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF (ok_strato) THEN
             CALL SUGWD_strato(klon,klev,paprs,pplay)
           ELSE
             CALL SUGWD(klon,klev,paprs,pplay)
           ENDIF
           
           DO i=1,klon
             zuthe(i)=0.
             zvthe(i)=0.
             if(zstd(i).gt.10.)then
               zuthe(i)=(1.-zgam(i))*cos(zthe(i))
               zvthe(i)=(1.-zgam(i))*sin(zthe(i))
             endif
           ENDDO
         ENDIF
c
c
         lmt_pas = NINT(86400./dtime * 1.0)   ! tous les jours
         WRITE(lunout,*)'La frequence de lecture surface est de ', 
     .                   lmt_pas
c
cIM 030306 END

      capemaxcels = 't_max(X)'
      t2mincels = 't_min(X)'
      t2maxcels = 't_max(X)'
      tinst = 'inst(X)'
      tave = 'ave(X)'
cIM cf. AM 081204 BEG
      write(lunout,*)'AVANT HIST IFLAG_CON=',iflag_con
cIM cf. AM 081204 END
c
c=============================================================
c   Initialisation des sorties
c=============================================================

#ifdef CPP_IOIPSL

c$OMP MASTER
       call phys_output_open(jjmp1,nlevSTD,clevSTD,nbteta, 
     &                        ctetaSTD,dtime,ok_veget,
     &                        type_ocean,iflag_pbl,ok_mensuel,ok_journe,
     &                        ok_hf,ok_instan,ok_LES,ok_ade,ok_aie, 
     &                        read_climoz, new_aod, aerosol_couple
     &                        )
c$OMP END MASTER
c$OMP BARRIER

#ifdef histISCCP
#include "ini_histISCCP.h"
#endif

#ifdef histNMC
#include "ini_histhfNMC.h"
#include "ini_histdayNMC.h"
#include "ini_histmthNMC.h"
#endif

#include "ini_histday_seri.h"

#include "ini_paramLMDZ_phy.h"

#endif

cIM 250308bad guide        ecrit_hf2mth = 30*1/ecrit_hf 
         ecrit_hf2mth = ecrit_mth/ecrit_hf

         ecrit_hf = ecrit_hf * un_jour
cIM
         IF(ecrit_day.LE.1.) THEN
          ecrit_day = ecrit_day * un_jour !en secondes
         ENDIF
cIM
         ecrit_mth = ecrit_mth * un_jour
         ecrit_ins = ecrit_ins * un_jour
         ecrit_reg = ecrit_reg * un_jour
         ecrit_tra = ecrit_tra * un_jour
         ecrit_ISCCP = ecrit_ISCCP * un_jour
         ecrit_LES = ecrit_LES * un_jour
c
         PRINT*,'physiq ecrit_ hf day mth reg tra ISCCP hf2mth',
     .   ecrit_hf,ecrit_day,ecrit_mth,ecrit_reg,ecrit_tra,ecrit_ISCCP,
     .   ecrit_hf2mth
cIM 030306 END


cXXXPB Positionner date0 pour initialisation de ORCHIDEE
      date0 = jD_ref 
      WRITE(*,*) 'physiq date0 : ',date0
c
c
c
c Prescrire l'ozone dans l'atmosphere
c
c
cc         DO i = 1, klon
cc         DO k = 1, klev
cc            CALL o3cm (paprs(i,k)/100.,paprs(i,k+1)/100., wo(i,k),20)
cc         ENDDO
cc         ENDDO
c
      IF (config_inca /= 'none') THEN
#ifdef INCA
         CALL VTe(VTphysiq)
         CALL VTb(VTinca)
!         iii = MOD(NINT(xjour),360)
!         calday = REAL(iii) + jH_cur
         calday = REAL(days_elapsed) + jH_cur
         WRITE(lunout,*) 'initial time chemini', days_elapsed, calday

         CALL chemini( 
     $                   rg,
     $                   ra,
     $                   airephy,
     $                   rlat,
     $                   rlon,
     $                   presnivs,
     $                   calday,
     $                   klon,
     $                   nqtot,
     $                   pdtphys,
     $                   annee_ref,
     $                   day_ref, 
     $                   itau_phy)

         CALL VTe(VTinca)
         CALL VTb(VTphysiq)
#endif
      END IF
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nouvelle initialisation pour le rayonnement RRTM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call iniradia(klon,klev,paprs(1,1:klev+1))

C$omp single
      if (read_climoz >= 1) then
         call open_climoz(ncid_climoz, press_climoz)
      END IF
C$omp end single
      ENDIF
!
!   ****************     Fin  de   IF ( debut  )   ***************
!
!
! Incrementer le compteur de la physique
!
      itap   = itap + 1
!
! Update fraction of the sub-surfaces (pctsrf) and 
! initialize, where a new fraction has appeared, all variables depending 
! on the surface fraction.
!
      CALL change_srf_frac(itap, dtime, days_elapsed+1, 
     *     pctsrf, falb1, falb2, ftsol, u10m, v10m, pbl_tke)

! Tendances bidons pour les processus qui n'affectent pas certaines
! variables.
      du0(:,:)=0.
      dv0(:,:)=0.
      dq0(:,:)=0.
      dql0(:,:)=0.
c
c Mettre a zero des variables de sortie (pour securite)
c
      DO i = 1, klon
         d_ps(i) = 0.0
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
         d_t(i,k) = 0.0
         d_u(i,k) = 0.0
         d_v(i,k) = 0.0
      ENDDO
      ENDDO
      DO iq = 1, nqtot
      DO k = 1, klev
      DO i = 1, klon
         d_qx(i,k,iq) = 0.0
      ENDDO
      ENDDO
      ENDDO
      da(:,:)=0.
      mp(:,:)=0.
      phi(:,:,:)=0.
c
c Ne pas affecter les valeurs entrees de u, v, h, et q
c
      DO k = 1, klev
      DO i = 1, klon
         t_seri(i,k)  = t(i,k)
         u_seri(i,k)  = u(i,k)
         v_seri(i,k)  = v(i,k)
         q_seri(i,k)  = qx(i,k,ivap)
         ql_seri(i,k) = qx(i,k,iliq)
         qs_seri(i,k) = 0.
      ENDDO
      ENDDO
      IF (nqtot.GE.3) THEN
      DO iq = 3, nqtot
      DO  k = 1, klev
      DO  i = 1, klon
         tr_seri(i,k,iq-2) = qx(i,k,iq)
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO k = 1, klev
      DO i = 1, klon
         tr_seri(i,k,1) = 0.0
      ENDDO
      ENDDO
      ENDIF
C
      DO i = 1, klon
        ztsol(i) = 0.
      ENDDO
      DO nsrf = 1, nbsrf
        DO i = 1, klon
          ztsol(i) = ztsol(i) + ftsol(i,nsrf)*pctsrf(i,nsrf)
        ENDDO
      ENDDO
cIM
      IF (ip_ebil_phy.ge.1) THEN 
        ztit='after dynamic'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,1,1,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
C     Comme les tendances de la physique sont ajoute dans la dynamique,
C     on devrait avoir que la variation d'entalpie par la dynamique
C     est egale a la variation de la physique au pas de temps precedent.
C     Donc la somme de ces 2 variations devrait etre nulle.
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol+d_h_vcol_phy, d_qt, 0.
     s      , fs_bound, fq_bound )
      END IF 

c Diagnostiquer la tendance dynamique
c
      IF (ancien_ok) THEN
         DO k = 1, klev
         DO i = 1, klon
            d_u_dyn(i,k) = (u_seri(i,k)-u_ancien(i,k))/dtime
            d_v_dyn(i,k) = (v_seri(i,k)-v_ancien(i,k))/dtime
            d_t_dyn(i,k) = (t_seri(i,k)-t_ancien(i,k))/dtime
            d_q_dyn(i,k) = (q_seri(i,k)-q_ancien(i,k))/dtime
         ENDDO
         ENDDO
      ELSE
         DO k = 1, klev
         DO i = 1, klon
            d_u_dyn(i,k) = 0.0
            d_v_dyn(i,k) = 0.0
            d_t_dyn(i,k) = 0.0
            d_q_dyn(i,k) = 0.0
         ENDDO
         ENDDO
         ancien_ok = .TRUE.
      ENDIF
c
c Ajouter le geopotentiel du sol:
c
      DO k = 1, klev
      DO i = 1, klon
         zphi(i,k) = pphi(i,k) + pphis(i)
      ENDDO
      ENDDO
c
c Verifier les temperatures
c
cIM BEG
      IF (check) THEN
       amn=MIN(ftsol(1,is_ter),1000.)
       amx=MAX(ftsol(1,is_ter),-1000.)
       DO i=2, klon
        amn=MIN(ftsol(i,is_ter),amn)
        amx=MAX(ftsol(i,is_ter),amx)
       ENDDO
c
       PRINT*,' debut avant hgardfou min max ftsol',itap,amn,amx
      ENDIF !(check) THEN
cIM END
c
      CALL hgardfou(t_seri,ftsol,'debutphy')
c
cIM BEG
      IF (check) THEN
       amn=MIN(ftsol(1,is_ter),1000.)
       amx=MAX(ftsol(1,is_ter),-1000.)
       DO i=2, klon
        amn=MIN(ftsol(i,is_ter),amn)
        amx=MAX(ftsol(i,is_ter),amx)
       ENDDO
c
       PRINT*,' debut apres hgardfou min max ftsol',itap,amn,amx
      ENDIF !(check) THEN
cIM END
c
c Mettre en action les conditions aux limites (albedo, sst, etc.).
c Prescrire l'ozone et calculer l'albedo sur l'ocean.
c
      if (read_climoz >= 1) then
C        Ozone from a file
!        Update required ozone index:
         ro3i = int((days_elapsed + jh_cur - jh_1jan)
     $        / ioget_year_len(year_cur) * 360.) + 1
         if (ro3i == 361) ro3i = 360
C        (This should never occur, except perhaps because of roundup
C        error. See documentation.)
         if (ro3i /= co3i) then
C           Update ozone field:
            if (read_climoz == 1) then
               call regr_pr_av(ncid_climoz, (/"tro3"/), julien=ro3i,
     $              press_in_edg=press_climoz, paprs=paprs, v3=wo)
            else
C              read_climoz == 2
               call regr_pr_av(ncid_climoz,
     $              (/"tro3         ", "tro3_daylight"/),
     $              julien=ro3i, press_in_edg=press_climoz, paprs=paprs,
     $              v3=wo)
            end if
!           Convert from mole fraction of ozone to column density of ozone in a
!           cell, in kDU:
            forall (l = 1: read_climoz) wo(:, :, l) = wo(:, :, l)
     $           * rmo3 / rmd * zmasse / dobson_u / 1e3
C           (By regridding ozone values for LMDZ only once every 360th of
C           year, we have already neglected the variation of pressure in one
C           360th of year. So do not recompute "wo" at each time step even if
C           "zmasse" changes a little.)
            co3i = ro3i
         end if
      elseif (MOD(itap-1,lmt_pas) == 0) THEN
C        Once per day, update ozone from Royer:
         wo(:, :, 1) = ozonecm(rlat, paprs, rjour=real(days_elapsed+1))
      ENDIF
c
c Re-evaporer l'eau liquide nuageuse
c
      DO k = 1, klev  ! re-evaporation de l'eau liquide nuageuse
      DO i = 1, klon
         zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
c        zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
         zlsdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
         zdelta = MAX(0.,SIGN(1.,RTT-t_seri(i,k)))
         zb = MAX(0.0,ql_seri(i,k))
         za = - MAX(0.0,ql_seri(i,k))
     .                  * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
         t_seri(i,k) = t_seri(i,k) + za
         q_seri(i,k) = q_seri(i,k) + zb
         ql_seri(i,k) = 0.0
         d_t_eva(i,k) = za
         d_q_eva(i,k) = zb
      ENDDO
      ENDDO
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after reevap'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,1,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
         call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
C
      END IF 

c
c=========================================================================
! Calculs de l'orbite.
! Necessaires pour le rayonnement et la surface (calcul de l'albedo).
! doit donc etre placé avant radlwsw et pbl_surface

!!!   jyg 17 Sep 2010 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ymds2ju(year_cur, mth_eq, day_eq,0., jD_eq)
      day_since_equinox = (jD_cur + jH_cur) - jD_eq 
!
!   choix entre calcul de la longitude solaire vraie ou valeur fixee a 
!   solarlong0
      if (solarlong0<-999.) then
       if (new_orbit) then
! calcul selon la routine utilisee pour les planetes
        call solarlong(day_since_equinox, zlongi, dist)
       else    
! calcul selon la routine utilisee pour l'AR4
        CALL orbite(REAL(days_elapsed+1),zlongi,dist)
       endif   
      else    
           zlongi=solarlong0  ! longitude solaire vraie
           dist=1.            ! distance au soleil / moyenne 
      endif   
      if(prt_level.ge.1)                                                &
     &    write(lunout,*)'Longitude solaire ',zlongi,solarlong0,dist

!  Avec ou sans cycle diurne
      IF (cycle_diurne) THEN
        zdtime=dtime*REAL(radpas) ! pas de temps du rayonnement (s)
        CALL zenang(zlongi,jH_cur,zdtime,rlat,rlon,rmu0,fract)
      ELSE
        CALL angle(zlongi, rlat, fract, rmu0)
      ENDIF

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Appel au pbl_surface : Planetary Boudary Layer et Surface
c Cela implique tous les interactions des sous-surfaces et la partie diffusion 
c turbulent du couche limit. 
c 
c Certains varibales de sorties de pbl_surface sont utiliser que pour 
c ecriture des fihiers hist_XXXX.nc, ces sont :
c   qsol,      zq2m,      s_pblh,  s_lcl,
c   s_capCL,   s_oliqCL,  s_cteiCL,s_pblT,
c   s_therm,   s_trmb1,   s_trmb2, s_trmb3,
c   zxrugs,    zu10m,     zv10m,   fder,
c   zxqsurf,   rh2m,      zxfluxu, zxfluxv,
c   frugs,     agesno,    fsollw,  fsolsw,
c   d_ts,      fevap,     fluxlat, t2m,
c   wfbils,    wfbilo,    fluxt,   fluxu, fluxv,
c
c Certains ne sont pas utiliser du tout : 
c   dsens, devap, zxsnow, zxfluxt, zxfluxq, q2m, fluxq
c

      CALL pbl_surface( 
     e     dtime,     date0,     itap,    days_elapsed+1,
     e     debut,     lafin,
     e     rlon,      rlat,      rugoro,  rmu0,     
     e     rain_fall, snow_fall, solsw,   sollw,    
     e     t_seri,    q_seri,    u_seri,  v_seri,   
     e     pplay,     paprs,     pctsrf,            
     +     ftsol,     falb1,     falb2,   u10m,   v10m,
     s     sollwdown, cdragh,    cdragm,  u1,    v1,
     s     albsol1,   albsol2,   sens,    evap,  
     s     zxtsol,    zxfluxlat, zt2m,    qsat2m, 
     s     d_t_vdf,   d_q_vdf,   d_u_vdf, d_v_vdf,
     s     coefh,     slab_wfbils,                
     d     qsol,      zq2m,      s_pblh,  s_lcl,
     d     s_capCL,   s_oliqCL,  s_cteiCL,s_pblT,
     d     s_therm,   s_trmb1,   s_trmb2, s_trmb3,
     d     zxrugs,    zu10m,     zv10m,   fder,
     d     zxqsurf,   rh2m,      zxfluxu, zxfluxv,
     d     frugs,     agesno,    fsollw,  fsolsw,
     d     d_ts,      fevap,     fluxlat, t2m,
     d     wfbils,    wfbilo,    fluxt,   fluxu,  fluxv,
     -     dsens,     devap,     zxsnow,
     -     zxfluxt,   zxfluxq,   q2m,     fluxq, pbl_tke )


!-----------------------------------------------------------------------------------------
! ajout des tendances de la diffusion turbulente
      CALL add_phys_tend(d_u_vdf,d_v_vdf,d_t_vdf,d_q_vdf,dql0,'vdf')
!-----------------------------------------------------------------------------------------

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif


      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after surface_main'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
         call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, sens
     e      , evap  , zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 

c =================================================================== c
c   Calcul de Qsat

      DO k = 1, klev
      DO i = 1, klon
         zx_t = t_seri(i,k)
         IF (thermcep) THEN
            zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
            zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
            zx_qs  = MIN(0.5,zx_qs)
            zcor   = 1./(1.-retv*zx_qs)
            zx_qs  = zx_qs*zcor
         ELSE
           IF (zx_t.LT.t_coup) THEN
              zx_qs = qsats(zx_t)/pplay(i,k)
           ELSE
              zx_qs = qsatl(zx_t)/pplay(i,k)
           ENDIF
         ENDIF
         zqsat(i,k)=zx_qs
      ENDDO
      ENDDO

      if (prt_level.ge.1) then
      write(lunout,*) 'L   qsat (g/kg) avant clouds_gno'
      write(lunout,'(i4,f15.4)') (k,1000.*zqsat(igout,k),k=1,klev)
      endif
c
c Appeler la convection (au choix)
c
      DO k = 1, klev
      DO i = 1, klon
         conv_q(i,k) = d_q_dyn(i,k) 
     .               + d_q_vdf(i,k)/dtime
         conv_t(i,k) = d_t_dyn(i,k) 
     .               + d_t_vdf(i,k)/dtime
      ENDDO
      ENDDO
      IF (check) THEN
         za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
         WRITE(lunout,*) "avantcon=", za
      ENDIF
      zx_ajustq = .FALSE.
      IF (iflag_con.EQ.2) zx_ajustq=.TRUE.
      IF (zx_ajustq) THEN
         DO i = 1, klon
            z_avant(i) = 0.0
         ENDDO
         DO k = 1, klev
         DO i = 1, klon
            z_avant(i) = z_avant(i) + (q_seri(i,k)+ql_seri(i,k))
     .                        *(paprs(i,k)-paprs(i,k+1))/RG
         ENDDO
         ENDDO
      ENDIF

c Calcule de vitesse verticale a partir de flux de masse verticale
      DO k = 1, klev
         DO i = 1, klon
            omega(i,k) = RG*flxmass_w(i,k) / airephy(i)
         END DO
      END DO
      if (prt_level.ge.1) write(lunout,*) 'omega(igout, :) = ',
     $     omega(igout, :)

      IF (iflag_con.EQ.1) THEN
        abort_message ='reactiver le call conlmd dans physiq.F'
        CALL abort_gcm (modname,abort_message,1)
c     CALL conlmd (dtime, paprs, pplay, t_seri, q_seri, conv_q,
c    .             d_t_con, d_q_con,
c    .             rain_con, snow_con, ibas_con, itop_con)
      ELSE IF (iflag_con.EQ.2) THEN
      CALL conflx(dtime, paprs, pplay, t_seri, q_seri,
     e            conv_t, conv_q, -evap, omega,
     s            d_t_con, d_q_con, rain_con, snow_con,
     s            pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     s            kcbot, kctop, kdtop, pmflxr, pmflxs)
      d_u_con = 0.
      d_v_con = 0.

      WHERE (rain_con < 0.) rain_con = 0.
      WHERE (snow_con < 0.) snow_con = 0.
      DO i = 1, klon
         ibas_con(i) = klev+1 - kcbot(i)
         itop_con(i) = klev+1 - kctop(i)
      ENDDO
      ELSE IF (iflag_con.GE.3) THEN
c nb of tracers for the KE convection:
c MAF la partie traceurs est faite dans phytrac
c on met ntra=1 pour limiter les appels mais on peut
c supprimer les calculs / ftra.
              ntra = 1

c=====================================================================================
cajout pour la parametrisation des poches froides: 
ccalcul de t_wake et t_undi: si pas de poches froides, t_wake=t_undi=t_seri 
      do k=1,klev
            do i=1,klon
             if (iflag_wake.eq.1) then
             t_wake(i,k) = t_seri(i,k)
     .           +(1-wake_s(i))*wake_deltat(i,k)
             q_wake(i,k) = q_seri(i,k)
     .           +(1-wake_s(i))*wake_deltaq(i,k)
             t_undi(i,k) = t_seri(i,k)
     .           -wake_s(i)*wake_deltat(i,k)
             q_undi(i,k) = q_seri(i,k)
     .           -wake_s(i)*wake_deltaq(i,k)
             else
             t_wake(i,k) = t_seri(i,k)
             q_wake(i,k) = q_seri(i,k)
             t_undi(i,k) = t_seri(i,k)
             q_undi(i,k) = q_seri(i,k)
             endif
            enddo
         enddo
      
cc--   Calcul de l'energie disponible ALE (J/kg) et de la puissance disponible ALP (W/m2)
cc--    pour le soulevement des particules dans le modele convectif
c
      do i = 1,klon
        ALE(i) = 0.
        ALP(i) = 0.
      enddo
c
ccalcul de ale_wake et alp_wake
       if (iflag_wake.eq.1) then
         if (itap .le. it_wape_prescr) then
          do i = 1,klon
           ale_wake(i) = wape_prescr
           alp_wake(i) = fip_prescr
          enddo
         else
          do i = 1,klon
cjyg  ALE=WAPE au lieu de ALE = 1/2 Cstar**2
ccc           ale_wake(i) = 0.5*wake_cstar(i)**2
           ale_wake(i) = wake_pe(i)
           alp_wake(i) = wake_fip(i)
          enddo
         endif
       else
         do i = 1,klon
           ale_wake(i) = 0.
           alp_wake(i) = 0.
         enddo
       endif
ccombinaison avec ale et alp de couche limite: constantes si pas de couplage, valeurs calculees
cdans le thermique sinon
       if (iflag_coupl.eq.0) then
          if (debut.and.prt_level.gt.9)
     $                     WRITE(lunout,*)'ALE et ALP imposes'
          do i = 1,klon
con ne couple que ale
c           ALE(i) = max(ale_wake(i),Ale_bl(i))
            ALE(i) = max(ale_wake(i),ale_bl_prescr)
con ne couple que alp
c           ALP(i) = alp_wake(i) + Alp_bl(i)
            ALP(i) = alp_wake(i) + alp_bl_prescr
          enddo
       else
         IF(prt_level>9)WRITE(lunout,*)'ALE et ALP couples au thermique'
!         do i = 1,klon
!             ALE(i) = max(ale_wake(i),Ale_bl(i))
! avant        ALP(i) = alp_wake(i) + Alp_bl(i)
!             ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
!         write(20,*)'ALE',ALE(i),Ale_bl(i),ale_wake(i)
!         write(21,*)'ALP',ALP(i),Alp_bl(i),alp_wake(i)
!         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modif FH 2010/04/27. Sans doute temporaire.
! Deux options pour le alp_offset : constant si >Ã‚ 0 ou proportionnel Ãƒa
! w si <0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do i = 1,klon
          ALE(i) = max(ale_wake(i),Ale_bl(i))
          if (alp_offset>=0.) then
            ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
          else
            ALP(i)=alp_wake(i)+Alp_bl(i)+alp_offset*min(omega(i,6),0.)
            if (alp(i)<0.) then
               print*,'ALP ',alp(i),alp_wake(i)
     s         ,Alp_bl(i),alp_offset*min(omega(i,6),0.)
            endif
          endif
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       endif
       do i=1,klon
          if (alp(i)>alp_max) then
             IF(prt_level>9)WRITE(lunout,*)                             &
     &       'WARNING SUPER ALP (seuil=',alp_max,
     ,       '): i, alp, alp_wake,ale',i,alp(i),alp_wake(i),ale(i)
             alp(i)=alp_max
          endif
          if (ale(i)>ale_max) then
             IF(prt_level>9)WRITE(lunout,*)                             &
     &       'WARNING SUPER ALE (seuil=',ale_max,
     ,       '): i, alp, alp_wake,ale',i,ale(i),ale_wake(i),alp(i)
             ale(i)=ale_max
          endif
       enddo

cfin calcul ale et alp
c=================================================================================================


c sb, oct02:
c Schema de convection modularise et vectorise:
c (driver commun aux versions 3 et 4)
c
          IF (ok_cvl) THEN ! new driver for convectL

          CALL concvl (iflag_con,iflag_clos,
     .        dtime,paprs,pplay,t_undi,q_undi,
     .        t_wake,q_wake,wake_s,
     .        u_seri,v_seri,tr_seri,nbtr,
     .        ALE,ALP,
     .        ema_work1,ema_work2,
     .        d_t_con,d_q_con,d_u_con,d_v_con,d_tr,
     .        rain_con, snow_con, ibas_con, itop_con, sigd,
     .        ema_cbmf,upwd,dnwd,dnwd0,
     .        Ma,mip,Vprecip,cape,cin,tvp,Tconv,iflagctrl,
     .        pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr,qcondc,wd,
     .        pmflxr,pmflxs,da,phi,mp,
     .        ftd,fqd,lalim_conv,wght_th)

cIM begin
c       print*,'physiq: cin pbase dnwd0 ftd fqd ',cin(1),pbase(1),
c    .dnwd0(1,1),ftd(1,1),fqd(1,1)
cIM end
cIM cf. FH
              clwcon0=qcondc
              pmfu(:,:)=upwd(:,:)+dnwd(:,:)

              do i = 1, klon
                if (iflagctrl(i).le.1) itau_con(i)=itau_con(i)+1
              enddo

          ELSE ! ok_cvl

c MAF conema3 ne contient pas les traceurs
          CALL conema3 (dtime,
     .        paprs,pplay,t_seri,q_seri,
     .        u_seri,v_seri,tr_seri,ntra,
     .        ema_work1,ema_work2,
     .        d_t_con,d_q_con,d_u_con,d_v_con,d_tr,
     .        rain_con, snow_con, ibas_con, itop_con,
     .        upwd,dnwd,dnwd0,bas,top,
     .        Ma,cape,tvp,rflag,
     .        pbase
     .        ,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr
     .        ,clwcon0)

          ENDIF ! ok_cvl

c
c Correction precip
          rain_con = rain_con * cvl_corr
          snow_con = snow_con * cvl_corr
c

           IF (.NOT. ok_gust) THEN
           do i = 1, klon
            wd(i)=0.0
           enddo
           ENDIF

c =================================================================== c
c Calcul des proprietes des nuages convectifs
c

c   calcul des proprietes des nuages convectifs
             clwcon0(:,:)=fact_cldcon*clwcon0(:,:)
             call clouds_gno
     s       (klon,klev,q_seri,zqsat,clwcon0,ptconv,ratqsc,rnebcon0)

c =================================================================== c

          DO i = 1, klon
           itop_con(i) = min(max(itop_con(i),1),klev)
           ibas_con(i) = min(max(ibas_con(i),1),itop_con(i))
          ENDDO

          DO i = 1, klon
            ema_pcb(i)  = paprs(i,ibas_con(i))
          ENDDO
          DO i = 1, klon
! L'idicage de itop_con peut cacher un pb potentiel
! FH sous la dictee de JYG, CR
            ema_pct(i)  = paprs(i,itop_con(i)+1)

            if (itop_con(i).gt.klev-3) then
              if(prt_level >= 9) then
                write(lunout,*)'La convection monte trop haut '
                write(lunout,*)'itop_con(,',i,',)=',itop_con(i)
              endif
            endif
          ENDDO     
      ELSE IF (iflag_con.eq.0) THEN
          write(lunout,*) 'On n appelle pas la convection'
          clwcon0=0.
          rnebcon0=0.
          d_t_con=0.
          d_q_con=0.
          d_u_con=0.
          d_v_con=0.
          rain_con=0.
          snow_con=0.
          bas=1
          top=1
      ELSE
          WRITE(lunout,*) "iflag_con non-prevu", iflag_con
          CALL abort
      ENDIF

c     CALL homogene(paprs, q_seri, d_q_con, u_seri,v_seri,
c    .              d_u_con, d_v_con)

!-----------------------------------------------------------------------------------------
! ajout des tendances de la diffusion turbulente
      CALL add_phys_tend(d_u_con,d_v_con,d_t_con,d_q_con,dql0,'con')
!-----------------------------------------------------------------------------------------

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif

cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after convect'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
         call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, rain_con, snow_con, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
C
      IF (check) THEN
          za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
          WRITE(lunout,*)"aprescon=", za
          zx_t = 0.0
          za = 0.0
          DO i = 1, klon
            za = za + airephy(i)/REAL(klon)
            zx_t = zx_t + (rain_con(i)+
     .                   snow_con(i))*airephy(i)/REAL(klon)
          ENDDO
          zx_t = zx_t/za*dtime
          WRITE(lunout,*)"Precip=", zx_t
      ENDIF
      IF (zx_ajustq) THEN
          DO i = 1, klon
            z_apres(i) = 0.0
          ENDDO
          DO k = 1, klev
            DO i = 1, klon
              z_apres(i) = z_apres(i) + (q_seri(i,k)+ql_seri(i,k))
     .            *(paprs(i,k)-paprs(i,k+1))/RG
            ENDDO
          ENDDO
          DO i = 1, klon
            z_factor(i) = (z_avant(i)-(rain_con(i)+snow_con(i))*dtime)
     .          /z_apres(i)
          ENDDO
          DO k = 1, klev
            DO i = 1, klon
              IF (z_factor(i).GT.(1.0+1.0E-08) .OR.
     .            z_factor(i).LT.(1.0-1.0E-08)) THEN
                  q_seri(i,k) = q_seri(i,k) * z_factor(i)
              ENDIF
            ENDDO
          ENDDO
      ENDIF
      zx_ajustq=.FALSE.

c
c=============================================================================
cRR:Evolution de la poche froide: on ne fait pas de separation wake/env 
cpour la couche limite diffuse pour l instant
c
      if (iflag_wake.eq.1) then
      DO k=1,klev
        DO i=1,klon
          dt_dwn(i,k)  = ftd(i,k) 
          wdt_PBL(i,k) = 0.
          dq_dwn(i,k)  = fqd(i,k) 
          wdq_PBL(i,k) = 0.
          M_dwn(i,k)   = dnwd0(i,k)
          M_up(i,k)    = upwd(i,k)
          dt_a(i,k)    = d_t_con(i,k)/dtime - ftd(i,k) 
          udt_PBL(i,k) = 0.
          dq_a(i,k)    = d_q_con(i,k)/dtime - fqd(i,k)
          udq_PBL(i,k) = 0.
        ENDDO
      ENDDO
c
ccalcul caracteristiques de la poche froide
      call calWAKE (paprs,pplay,dtime
     :               ,t_seri,q_seri,omega
     :               ,dt_dwn,dq_dwn,M_dwn,M_up
     :               ,dt_a,dq_a,sigd
     :               ,wdt_PBL,wdq_PBL
     :               ,udt_PBL,udq_PBL
     o               ,wake_deltat,wake_deltaq,wake_dth
     o               ,wake_h,wake_s,wake_dens
     o               ,wake_pe,wake_fip,wake_gfl
     o               ,dt_wake,dq_wake
     o               ,wake_k, t_undi,q_undi
     o               ,wake_omgbdth,wake_dp_omgb
     o               ,wake_dtKE,wake_dqKE
     o               ,wake_dtPBL,wake_dqPBL
     o               ,wake_omg,wake_dp_deltomg
     o               ,wake_spread,wake_Cstar,wake_d_deltat_gw
     o               ,wake_ddeltat,wake_ddeltaq)
c
!-----------------------------------------------------------------------------------------
! ajout des tendances des poches froides
! Faire rapidement disparaitre l'ancien dt_wake pour garder un d_t_wake
! coherent avec les autres d_t_...
      d_t_wake(:,:)=dt_wake(:,:)*dtime
      d_q_wake(:,:)=dq_wake(:,:)*dtime
      CALL add_phys_tend(du0,dv0,d_t_wake,d_q_wake,dql0,'wake')
!-----------------------------------------------------------------------------------------

      endif
c
c===================================================================
cJYG
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after wake'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 

c      print*,'apres callwake iflag_cldcon=', iflag_cldcon
c
c===================================================================
c Convection seche (thermiques ou ajustement)
c===================================================================
c
       call stratocu_if(klon,klev,pctsrf,paprs, pplay,t_seri
     s ,seuil_inversion,weak_inversion,dthmin)



      d_t_ajsb(:,:)=0.
      d_q_ajsb(:,:)=0.
      d_t_ajs(:,:)=0.
      d_u_ajs(:,:)=0.
      d_v_ajs(:,:)=0.
      d_q_ajs(:,:)=0.
      clwcon0th(:,:)=0.
c
c      fm_therm(:,:)=0.
c      entr_therm(:,:)=0.
c      detr_therm(:,:)=0.
c
      IF(prt_level>9)WRITE(lunout,*)
     .    'AVANT LA CONVECTION SECHE , iflag_thermals='
     s   ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals
      if(iflag_thermals.lt.0) then
c  Rien
c  ====
         IF(prt_level>9)WRITE(lunout,*)'pas de convection'


      else

c  Thermiques
c  ==========
         IF(prt_level>9)WRITE(lunout,*)'JUSTE AVANT , iflag_thermals='
     s   ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals


         if (iflag_thermals.gt.1) then
         call calltherm(pdtphys
     s      ,pplay,paprs,pphi,weak_inversion
     s      ,u_seri,v_seri,t_seri,q_seri,zqsat,debut
     s      ,d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs
     s      ,fm_therm,entr_therm,detr_therm
     s      ,zqasc,clwcon0th,lmax_th,ratqscth
     s      ,ratqsdiff,zqsatth
con rajoute ale et alp, et les caracteristiques de la couche alim
     s      ,Ale_bl,Alp_bl,lalim_conv,wght_th, zmax0, f0, zw2,fraca
     s      ,ztv,zpspsk,ztla,zthl)

! ----------------------------------------------------------------------
! Transport de la TKE par les panaches thermiques.
! FH : 2010/02/01
      if (iflag_pbl.eq.10) then
      call thermcell_dtke(klon,klev,nbsrf,pdtphys,fm_therm,entr_therm,
     s           rg,paprs,pbl_tke)
      endif
! ----------------------------------------------------------------------

         endif



c  Ajustement sec
c  ==============

! Dans le cas oÃ¹ on active les thermiques, on fait partir l'ajustement
! a partir du sommet des thermiques.
! Dans le cas contraire, on demarre au niveau 1.

         if (iflag_thermals.ge.13.or.iflag_thermals.eq.0) then

         if(iflag_thermals.eq.0) then
            IF(prt_level>9)WRITE(lunout,*)'ajsec'
            limbas(:)=1
         else
            limbas(:)=lmax_th(:)
         endif
  
! Attention : le call ajsec_convV2 n'est maintenu que momentanneement
! pour des test de convergence numerique.
! Le nouveau ajsec est a priori mieux, meme pour le cas 
! iflag_thermals = 0 (l'ancienne version peut faire des tendances
! non nulles numeriquement pour des mailles non concernees.

         if (iflag_thermals.eq.0) then
            CALL ajsec_convV2(paprs, pplay, t_seri,q_seri
     s      , d_t_ajsb, d_q_ajsb)
         else
            CALL ajsec(paprs, pplay, t_seri,q_seri,limbas
     s      , d_t_ajsb, d_q_ajsb)
         endif

!-----------------------------------------------------------------------------------------
! ajout des tendances de l'ajustement sec ou des thermiques
      CALL add_phys_tend(du0,dv0,d_t_ajsb,d_q_ajsb,dql0,'ajsb')
         d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_ajsb(:,:)
         d_q_ajs(:,:)=d_q_ajs(:,:)+d_q_ajsb(:,:)

!-----------------------------------------------------------------------------------------

         endif

      endif
c
c===================================================================
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after dry_adjust'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 


c-------------------------------------------------------------------------
c  Caclul des ratqs
c-------------------------------------------------------------------------

c      print*,'calcul des ratqs'
c   ratqs convectifs a l'ancienne en fonction de q(z=0)-q / q
c   ----------------
c   on ecrase le tableau ratqsc calcule par clouds_gno
      if (iflag_cldcon.eq.1) then
         do k=1,klev
         do i=1,klon
            if(ptconv(i,k)) then
              ratqsc(i,k)=ratqsbas
     s        +fact_cldcon*(q_seri(i,1)-q_seri(i,k))/q_seri(i,k)
            else
               ratqsc(i,k)=0.
            endif
         enddo
         enddo

c-----------------------------------------------------------------------
c  par nversion de la fonction log normale
c-----------------------------------------------------------------------
      else if (iflag_cldcon.eq.4) then
         ptconvth(:,:)=.false.
         ratqsc(:,:)=0.
         if(prt_level.ge.9) print*,'avant clouds_gno thermique'
         call clouds_gno
     s   (klon,klev,q_seri,zqsat,clwcon0th,ptconvth,ratqsc,rnebcon0th)
         if(prt_level.ge.9) print*,' CLOUDS_GNO OK'
       
       endif

c   ratqs stables
c   -------------

      if (iflag_ratqs.eq.0) then

! Le cas iflag_ratqs=0 correspond a la version IPCC 2005 du modele.
         do k=1,klev
            do i=1, klon
               ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)*
     s         min((paprs(i,1)-pplay(i,k))/(paprs(i,1)-30000.),1.) 
            enddo 
         enddo

! Pour iflag_ratqs=1 ou 2, le ratqs est constant au dessus de 
! 300 hPa (ratqshaut), varie lineariement en fonction de la pression
! entre 600 et 300 hPa et est soit constant (ratqsbas) pour iflag_ratqs=1
! soit lineaire (entre 0 a la surface et ratqsbas) pour iflag_ratqs=2
! Il s'agit de differents tests dans la phase de reglage du modele
! avec thermiques.

      else if (iflag_ratqs.eq.1) then

         do k=1,klev
            do i=1, klon
               if (pplay(i,k).ge.60000.) then
                  ratqss(i,k)=ratqsbas
               else if ((pplay(i,k).ge.30000.).and.
     s            (pplay(i,k).lt.60000.)) then
                  ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)*
     s            (60000.-pplay(i,k))/(60000.-30000.)
               else
                  ratqss(i,k)=ratqshaut
               endif
            enddo
         enddo

      else

         do k=1,klev
            do i=1, klon
               if (pplay(i,k).ge.60000.) then
                  ratqss(i,k)=ratqsbas
     s            *(paprs(i,1)-pplay(i,k))/(paprs(i,1)-60000.)
               else if ((pplay(i,k).ge.30000.).and.
     s             (pplay(i,k).lt.60000.)) then
                    ratqss(i,k)=ratqsbas+(ratqshaut-ratqsbas)*
     s              (60000.-pplay(i,k))/(60000.-30000.)
               else
                    ratqss(i,k)=ratqshaut
               endif
            enddo
         enddo
      endif




c  ratqs final
c  -----------

      if (iflag_cldcon.eq.1 .or.iflag_cldcon.eq.2
     s    .or.iflag_cldcon.eq.4) then

! On ajoute une constante au ratqsc*2 pour tenir compte de 
! fluctuations turbulentes de petite echelle

         do k=1,klev
            do i=1,klon
               if ((fm_therm(i,k).gt.1.e-10)) then
                  ratqsc(i,k)=sqrt(ratqsc(i,k)**2+0.05**2)
               endif
            enddo
         enddo

!   les ratqs sont une combinaison de ratqss et ratqsc
       if(prt_level.ge.9)
     $       write(lunout,*)'PHYLMD NOUVEAU TAU_RATQS ',tau_ratqs

         if (tau_ratqs>1.e-10) then
            facteur=exp(-pdtphys/tau_ratqs)
         else
            facteur=0.
         endif
         ratqs(:,:)=ratqsc(:,:)*(1.-facteur)+ratqs(:,:)*facteur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 22/09/2009
! La ligne ci-dessous faisait osciller le modele et donnait une solution
! assymptotique bidon et dÃ©pendant fortement du pas de temps.
!        ratqs(:,:)=sqrt(ratqs(:,:)**2+ratqss(:,:)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ratqs(:,:)=max(ratqs(:,:),ratqss(:,:))
      else
!   on ne prend que le ratqs stable pour fisrtilp
         ratqs(:,:)=ratqss(:,:)
      endif


c
c Appeler le processus de condensation a grande echelle
c et le processus de precipitation
c-------------------------------------------------------------------------
      CALL fisrtilp(dtime,paprs,pplay,
     .           t_seri, q_seri,ptconv,ratqs,
     .           d_t_lsc, d_q_lsc, d_ql_lsc, rneb, cldliq,
     .           rain_lsc, snow_lsc,
     .           pfrac_impa, pfrac_nucl, pfrac_1nucl,
     .           frac_impa, frac_nucl,
     .           prfl, psfl, rhcl, 
     .           zqasc, fraca,ztv,zpspsk,ztla,zthl,iflag_cldcon )

      WHERE (rain_lsc < 0) rain_lsc = 0.
      WHERE (snow_lsc < 0) snow_lsc = 0.
!-----------------------------------------------------------------------------------------
! ajout des tendances de la diffusion turbulente
      CALL add_phys_tend(du0,dv0,d_t_lsc,d_q_lsc,d_ql_lsc,'lsc')
!-----------------------------------------------------------------------------------------
      DO k = 1, klev
      DO i = 1, klon
         cldfra(i,k) = rneb(i,k)
         IF (.NOT.new_oliq) cldliq(i,k) = ql_seri(i,k)
      ENDDO
      ENDDO
      IF (check) THEN
         za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
         WRITE(lunout,*)"apresilp=", za
         zx_t = 0.0
         za = 0.0
         DO i = 1, klon
            za = za + airephy(i)/REAL(klon)
            zx_t = zx_t + (rain_lsc(i)
     .                  + snow_lsc(i))*airephy(i)/REAL(klon)
        ENDDO
         zx_t = zx_t/za*dtime
         WRITE(lunout,*)"Precip=", zx_t
      ENDIF
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after fisrt'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, rain_lsc, snow_lsc, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif

c
c-------------------------------------------------------------------
c  PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT
c-------------------------------------------------------------------

c 1. NUAGES CONVECTIFS
c
cIM cf FH
c     IF (iflag_cldcon.eq.-1) THEN ! seulement pour Tiedtke
      IF (iflag_cldcon.le.-1) THEN ! seulement pour Tiedtke
       snow_tiedtke=0.
c     print*,'avant calcul de la pseudo precip '
c     print*,'iflag_cldcon',iflag_cldcon
       if (iflag_cldcon.eq.-1) then
          rain_tiedtke=rain_con
       else
c       print*,'calcul de la pseudo precip '
          rain_tiedtke=0.
c         print*,'calcul de la pseudo precip 0'
          do k=1,klev
          do i=1,klon
             if (d_q_con(i,k).lt.0.) then
                rain_tiedtke(i)=rain_tiedtke(i)-d_q_con(i,k)/pdtphys
     s         *(paprs(i,k)-paprs(i,k+1))/rg
             endif
          enddo
          enddo
       endif
c
c     call dump2d(iim,jjm,rain_tiedtke(2:klon-1),'PSEUDO PRECIP ')
c

c Nuages diagnostiques pour Tiedtke
      CALL diagcld1(paprs,pplay,
cIM cf FH  .             rain_con,snow_con,ibas_con,itop_con,
     .             rain_tiedtke,snow_tiedtke,ibas_con,itop_con,
     .             diafra,dialiq)
      DO k = 1, klev
      DO i = 1, klon
      IF (diafra(i,k).GT.cldfra(i,k)) THEN
         cldliq(i,k) = dialiq(i,k)
         cldfra(i,k) = diafra(i,k)
      ENDIF
      ENDDO
      ENDDO

      ELSE IF (iflag_cldcon.ge.3) THEN
c  On prend pour les nuages convectifs le max du calcul de la
c  convection et du calcul du pas de temps precedent diminue d'un facteur
c  facttemps
      facteur = pdtphys *facttemps
      do k=1,klev
         do i=1,klon
            rnebcon(i,k)=rnebcon(i,k)*facteur
            if (rnebcon0(i,k)*clwcon0(i,k).gt.rnebcon(i,k)*clwcon(i,k))
     s      then
                rnebcon(i,k)=rnebcon0(i,k)
                clwcon(i,k)=clwcon0(i,k)
            endif
         enddo
      enddo

c
cjq - introduce the aerosol direct and first indirect radiative forcings
cjq - Johannes Quaas, 27/11/2003 (quaas@lmd.jussieu.fr)
      IF (ok_ade.OR.ok_aie) THEN
         IF (.NOT. aerosol_couple)
     &        CALL readaerosol_optic(
     &        debut, new_aod, flag_aerosol, itap, jD_cur-jD_ref,
     &        pdtphys, pplay, paprs, t_seri, rhcl, presnivs, 
     &        mass_solu_aero, mass_solu_aero_pi, 
     &        tau_aero, piz_aero, cg_aero, 
     &        tausum_aero, tau3d_aero)
      ELSE
cIM 170310 BEG
         tausum_aero(:,:,:) = 0.
cIM 170310 END
         tau_aero(:,:,:,:) = 0.
         piz_aero(:,:,:,:) = 0.
         cg_aero(:,:,:,:)  = 0.
      ENDIF

cIM calcul nuages par le simulateur ISCCP
c
#ifdef histISCCP
      IF (ok_isccp) THEN
c
cIM lecture invtau, tautab des fichiers formattes
c
      IF (debut) THEN
c$OMP MASTER
c
      open(99,file='tautab.formatted', FORM='FORMATTED')
      read(99,'(f30.20)') tautab_omp
      close(99)
c
      open(99,file='invtau.formatted',form='FORMATTED')
      read(99,'(i10)') invtau_omp

c     print*,'calcul_simulISCCP invtau_omp',invtau_omp
c     write(6,'(a,8i10)') 'invtau_omp',(invtau_omp(i),i=1,100)

      close(99)
c$OMP END MASTER
c$OMP BARRIER 
      tautab=tautab_omp
      invtau=invtau_omp
c
      ENDIF !debut
c
cIM appel simulateur toutes les  NINT(freq_ISCCP/dtime) heures
       IF (MOD(itap,NINT(freq_ISCCP/dtime)).EQ.0) THEN
#include "calcul_simulISCCP.h"
       ENDIF !(MOD(itap,NINT(freq_ISCCP/dtime))
      ENDIF !ok_isccp
#endif

c   On prend la somme des fractions nuageuses et des contenus en eau

      if (iflag_cldcon>=5) then

! Si on est sur un point touche par la convection profonde et pas
! par les thermiques, on prend la couverture nuageuse et l'eau nuageuse
! de la convection profonde.

         print*,'TEST SCHEMA DE NUAGES '
         do k=1,klev
            do i=1,klon
               if (ptconv(i,k).and. .not. ptconvth(i,k)) then
                   cldfra(i,k)=rnebcon(i,k)
                   cldliq(i,k)=rnebcon(i,k)*clwcon(i,k)
               endif
            enddo
         enddo
      else
! Ancienne version
      cldfra(:,:)=min(max(cldfra(:,:),rnebcon(:,:)),1.)
      cldliq(:,:)=cldliq(:,:)+rnebcon(:,:)*clwcon(:,:)
      endif

      ENDIF
c
c 2. NUAGES STARTIFORMES
c
      IF (ok_stratus) THEN
      CALL diagcld2(paprs,pplay,t_seri,q_seri, diafra,dialiq)
      DO k = 1, klev
      DO i = 1, klon
      IF (diafra(i,k).GT.cldfra(i,k)) THEN
         cldliq(i,k) = dialiq(i,k)
         cldfra(i,k) = diafra(i,k)
      ENDIF
      ENDDO
      ENDDO
      ENDIF
c
c Precipitation totale
c
      DO i = 1, klon
         rain_fall(i) = rain_con(i) + rain_lsc(i)
         snow_fall(i) = snow_con(i) + snow_lsc(i)
      ENDDO
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit="after diagcld"
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
c
c Calculer l'humidite relative pour diagnostique
c
      DO k = 1, klev
      DO i = 1, klon
         zx_t = t_seri(i,k)
         IF (thermcep) THEN
            zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
            zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
            zx_qs  = MIN(0.5,zx_qs)
            zcor   = 1./(1.-retv*zx_qs)
            zx_qs  = zx_qs*zcor
         ELSE
           IF (zx_t.LT.t_coup) THEN
              zx_qs = qsats(zx_t)/pplay(i,k)
           ELSE
              zx_qs = qsatl(zx_t)/pplay(i,k)
           ENDIF
         ENDIF
         zx_rh(i,k) = q_seri(i,k)/zx_qs
         zqsat(i,k)=zx_qs
      ENDDO
      ENDDO

cIM Calcul temp.potentielle a 2m (tpot) et temp. potentielle 
c   equivalente a 2m (tpote) pour diagnostique
c
      DO i = 1, klon
       tpot(i)=zt2m(i)*(100000./paprs(i,1))**RKAPPA
       IF (thermcep) THEN
        IF(zt2m(i).LT.RTT) then
         Lheat=RLSTT
        ELSE
         Lheat=RLVTT
        ENDIF
       ELSE
        IF (zt2m(i).LT.RTT) THEN
         Lheat=RLSTT
        ELSE
         Lheat=RLVTT
        ENDIF
       ENDIF
       tpote(i) = tpot(i)*     
     . EXP((Lheat *qsat2m(i))/(RCPD*zt2m(i)))
      ENDDO

      IF (config_inca /= 'none') THEN
#ifdef INCA
         CALL VTe(VTphysiq)
         CALL VTb(VTinca)
         calday = REAL(days_elapsed + 1) + jH_cur

         call chemtime(itap+itau_phy-1, date0, dtime)
         IF (config_inca == 'aero') THEN
            CALL AEROSOL_METEO_CALC(
     $           calday,pdtphys,pplay,paprs,t,pmflxr,pmflxs,
     $           prfl,psfl,pctsrf,airephy,rlat,rlon,u10m,v10m)
         END IF

         zxsnow_dummy(:) = 0.0

         CALL chemhook_begin (calday,
     $                          days_elapsed+1,
     $                          jH_cur,
     $                          pctsrf(1,1),
     $                          rlat,
     $                          rlon,
     $                          airephy,
     $                          paprs,
     $                          pplay,
     $                          coefh,
     $                          pphi,
     $                          t_seri,
     $                          u,
     $                          v,
     $                          wo(:, :, 1),
     $                          q_seri,
     $                          zxtsol,
     $                          zxsnow_dummy,
     $                          solsw,
     $                          albsol1,
     $                          rain_fall,
     $                          snow_fall,
     $                          itop_con,
     $                          ibas_con,
     $                          cldfra,
     $                          iim,
     $                          jjm,
     $                          tr_seri,
     $                          ftsol,
     $                          paprs,
     $                          cdragh,
     $                          cdragm,
     $                          pctsrf,
     $                          pdtphys,
     $                            itap)

         CALL VTe(VTinca)
         CALL VTb(VTphysiq)
#endif 
      END IF !config_inca /= 'none'
c     
c Calculer les parametres optiques des nuages et quelques
c parametres pour diagnostiques:
c

      IF (aerosol_couple) THEN 
         mass_solu_aero(:,:)    = ccm(:,:,1) 
         mass_solu_aero_pi(:,:) = ccm(:,:,2) 
      END IF

      if (ok_newmicro) then
      CALL newmicro (paprs, pplay,ok_newmicro,
     .            t_seri, cldliq, cldfra, cldtau, cldemi,
     .            cldh, cldl, cldm, cldt, cldq,
     .            flwp, fiwp, flwc, fiwc,
     e            ok_aie,
     e            mass_solu_aero, mass_solu_aero_pi,
     e            bl95_b0, bl95_b1,
     s            cldtaupi, re, fl, ref_liq, ref_ice)
      else
      CALL nuage (paprs, pplay,
     .            t_seri, cldliq, cldfra, cldtau, cldemi,
     .            cldh, cldl, cldm, cldt, cldq,
     e            ok_aie,
     e            mass_solu_aero, mass_solu_aero_pi,
     e            bl95_b0, bl95_b1,
     s            cldtaupi, re, fl)
      
      endif
c
c Appeler le rayonnement mais calculer tout d'abord l'albedo du sol.
c
      IF (MOD(itaprad,radpas).EQ.0) THEN

      DO i = 1, klon
         albsol1(i) = falb1(i,is_oce) * pctsrf(i,is_oce)
     .             + falb1(i,is_lic) * pctsrf(i,is_lic)
     .             + falb1(i,is_ter) * pctsrf(i,is_ter)
     .             + falb1(i,is_sic) * pctsrf(i,is_sic)
         albsol2(i) = falb2(i,is_oce) * pctsrf(i,is_oce)
     .               + falb2(i,is_lic) * pctsrf(i,is_lic)
     .               + falb2(i,is_ter) * pctsrf(i,is_ter)
     .               + falb2(i,is_sic) * pctsrf(i,is_sic)
      ENDDO

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
       call writefield_phy('q_seri',q_seri,llm)
      endif
      
      IF (aerosol_couple) THEN 
#ifdef INCA
         CALL radlwsw_inca 
     e        (kdlon,kflev,dist, rmu0, fract, solaire,
     e        paprs, pplay,zxtsol,albsol1, albsol2, t_seri,q_seri,
     e        wo(:, :, 1),
     e        cldfra, cldemi, cldtau,
     s        heat,heat0,cool,cool0,radsol,albpla,
     s        topsw,toplw,solsw,sollw,
     s        sollwdown,
     s        topsw0,toplw0,solsw0,sollw0,
     s        lwdn0, lwdn, lwup0, lwup, 
     s        swdn0, swdn, swup0, swup,
     e        ok_ade, ok_aie,
     e        tau_aero, piz_aero, cg_aero,
     s        topswad_aero, solswad_aero,
     s        topswad0_aero, solswad0_aero,
     s        topsw_aero, topsw0_aero,
     s        solsw_aero, solsw0_aero,
     e        cldtaupi,
     s        topswai_aero, solswai_aero)
            
#endif
      ELSE

         CALL radlwsw
     e        (dist, rmu0, fract, 
     e        paprs, pplay,zxtsol,albsol1, albsol2, 
     e        t_seri,q_seri,wo,
     e        cldfra, cldemi, cldtau,
     e        ok_ade, ok_aie,
     e        tau_aero, piz_aero, cg_aero,
     e        cldtaupi,new_aod,
     e        zqsat, flwc, fiwc,
     s        heat,heat0,cool,cool0,radsol,albpla,
     s        topsw,toplw,solsw,sollw,
     s        sollwdown,
     s        topsw0,toplw0,solsw0,sollw0,
     s        lwdn0, lwdn, lwup0, lwup, 
     s        swdn0, swdn, swup0, swup,
     s        topswad_aero, solswad_aero,
     s        topswai_aero, solswai_aero,
     o        topswad0_aero, solswad0_aero,
     o        topsw_aero, topsw0_aero,
     o        solsw_aero, solsw0_aero,
     o        topswcf_aero, solswcf_aero)
         

      ENDIF ! aerosol_couple
      itaprad = 0
      ENDIF ! MOD(itaprad,radpas)
      itaprad = itaprad + 1

      IF (iflag_radia.eq.0) THEN
         IF (prt_level.ge.9) THEN
            PRINT *,'--------------------------------------------------'
            PRINT *,'>>>> ATTENTION rayonnement desactive pour ce cas'
            PRINT *,'>>>>           heat et cool mis a zero '
            PRINT *,'--------------------------------------------------'
         END IF
         heat=0.
         cool=0.
      END IF

c
c Ajouter la tendance des rayonnements (tous les pas)
c
      DO k = 1, klev
      DO i = 1, klon
         t_seri(i,k) = t_seri(i,k)
     .               + (heat(i,k)-cool(i,k)) * dtime/RDAY
      ENDDO
      ENDDO
c
      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif
 
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after rad'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy
     e      , topsw, toplw, solsw, sollw, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
c
c
c Calculer l'hydrologie de la surface
c
c      CALL hydrol(dtime,pctsrf,rain_fall, snow_fall, zxevap,
c     .            agesno, ftsol,fqsurf,fsnow, ruis)
c

c
c Calculer le bilan du sol et la derive de temperature (couplage)
c
      DO i = 1, klon
c         bils(i) = radsol(i) - sens(i) - evap(i)*RLVTT
c a la demande de JLD
         bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
      ENDDO
c
cmoddeblott(jan95)
c Appeler le programme de parametrisation de l'orographie
c a l'echelle sous-maille:
c
      IF (ok_orodr) THEN
c
c  selection des points pour lesquels le shema est actif:
        igwd=0
        DO i=1,klon
        itest(i)=0
c        IF ((zstd(i).gt.10.0)) THEN
        IF (((zpic(i)-zmea(i)).GT.100.).AND.(zstd(i).GT.10.0)) THEN
          itest(i)=1
          igwd=igwd+1
          idx(igwd)=i
        ENDIF
        ENDDO
c        igwdim=MAX(1,igwd)
c
        IF (ok_strato) THEN
        
          CALL drag_noro_strato(klon,klev,dtime,paprs,pplay,
     e                   zmea,zstd, zsig, zgam, zthe,zpic,zval,
     e                   igwd,idx,itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrdr, zvstrdr,
     s                   d_t_oro, d_u_oro, d_v_oro)

       ELSE
        CALL drag_noro(klon,klev,dtime,paprs,pplay,
     e                   zmea,zstd, zsig, zgam, zthe,zpic,zval,
     e                   igwd,idx,itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrdr, zvstrdr,
     s                   d_t_oro, d_u_oro, d_v_oro)
       ENDIF
c
c  ajout des tendances
!-----------------------------------------------------------------------------------------
! ajout des tendances de la trainee de l'orographie
      CALL add_phys_tend(d_u_oro,d_v_oro,d_t_oro,dq0,dql0,'oro')
!-----------------------------------------------------------------------------------------
c
      ENDIF ! fin de test sur ok_orodr
c
      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif
      
      IF (ok_orolf) THEN
c
c  selection des points pour lesquels le shema est actif:
        igwd=0
        DO i=1,klon
        itest(i)=0
        IF ((zpic(i)-zmea(i)).GT.100.) THEN
          itest(i)=1
          igwd=igwd+1
          idx(igwd)=i
        ENDIF
        ENDDO
c        igwdim=MAX(1,igwd)
c
        IF (ok_strato) THEN

          CALL lift_noro_strato(klon,klev,dtime,paprs,pplay,
     e                   rlat,zmea,zstd,zpic,zgam,zthe,zpic,zval,
     e                   igwd,idx,itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrli, zvstrli,
     s                   d_t_lif, d_u_lif, d_v_lif               )
        
        ELSE
          CALL lift_noro(klon,klev,dtime,paprs,pplay,
     e                   rlat,zmea,zstd,zpic,
     e                   itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrli, zvstrli,
     s                   d_t_lif, d_u_lif, d_v_lif)
       ENDIF
c   
!-----------------------------------------------------------------------------------------
! ajout des tendances de la portance de l'orographie
      CALL add_phys_tend(d_u_lif,d_v_lif,d_t_lif,dq0,dql0,'lif')
!-----------------------------------------------------------------------------------------
c
      ENDIF ! fin de test sur ok_orolf
C  HINES GWD PARAMETRIZATION

       IF (ok_hines) then

         CALL hines_gwd(klon,klev,dtime,paprs,pplay,
     i                  rlat,t_seri,u_seri,v_seri,
     o                  zustrhi,zvstrhi,
     o                  d_t_hin, d_u_hin, d_v_hin)
c
c  ajout des tendances
        CALL add_phys_tend(d_u_hin,d_v_hin,d_t_hin,dq0,dql0,'lif')

      ENDIF
c

c
cIM cf. FLott BEG
C STRESS NECESSAIRES: TOUTE LA PHYSIQUE

      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif

      DO i = 1, klon
        zustrph(i)=0.
        zvstrph(i)=0.
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
       zustrph(i)=zustrph(i)+(u_seri(i,k)-u(i,k))/dtime*
     c            (paprs(i,k)-paprs(i,k+1))/rg
       zvstrph(i)=zvstrph(i)+(v_seri(i,k)-v(i,k))/dtime*
     c            (paprs(i,k)-paprs(i,k+1))/rg
      ENDDO
      ENDDO
c
cIM calcul composantes axiales du moment angulaire et couple des montagnes
c
      IF (is_sequential) THEN
      
        CALL aaam_bud (27,klon,klev,jD_cur-jD_ref,jH_cur,
     C                 ra,rg,romega,
     C                 rlat,rlon,pphis,
     C                 zustrdr,zustrli,zustrph,
     C                 zvstrdr,zvstrli,zvstrph,
     C                 paprs,u,v,
     C                 aam, torsfc)
       ENDIF
cIM cf. FLott END
cIM
      IF (ip_ebil_phy.ge.2) THEN 
        ztit='after orography'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
         call diagphy(airephy,ztit,ip_ebil_phy
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
c
c
!====================================================================
! Interface Simulateur COSP (Calipso, ISCCP, MISR, ..)
!====================================================================
! Abderrahmane 24.08.09

      IF (ok_cosp) THEN
! adeclarer 
#ifdef CPP_COSP
       IF (MOD(itap,NINT(freq_cosp/dtime)).EQ.0) THEN

       print*,'freq_cosp',freq_cosp
          mr_ozone=wo(:, :, 1) * dobson_u * 1e3 / zmasse
!       print*,'Dans physiq.F avant appel cosp ref_liq,ref_ice=',
!     s        ref_liq,ref_ice
          call phys_cosp(itap,dtime,freq_cosp,
     $                   ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP,
     $                   ecrit_mth,ecrit_day,ecrit_hf,
     $                   klon,klev,rlon,rlat,presnivs,overlap,
     $                   ref_liq,ref_ice,
     $                   pctsrf(:,is_ter)+pctsrf(:,is_lic),
     $                   zu10m,zv10m,
     $                   zphi,paprs(:,1:klev),pplay,zxtsol,t_seri,
     $                   qx(:,:,ivap),zx_rh,cldfra,rnebcon,flwc,fiwc,
     $                   prfl(:,1:klev),psfl(:,1:klev),
     $                   pmflxr(:,1:klev),pmflxs(:,1:klev),
     $                   mr_ozone,cldtau, cldemi)

!     L          calipso2D,calipso3D,cfadlidar,parasolrefl,atb,betamol,
!     L          cfaddbze,clcalipso2,dbze,cltlidarradar,
!     M          clMISR,
!     R          clisccp2,boxtauisccp,boxptopisccp,tclisccp,ctpisccp,
!     I          tauisccp,albisccp,meantbisccp,meantbclrisccp)

         ENDIF

#endif
       ENDIF  !ok_cosp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cAA
cAA Installation de l'interface online-offline pour traceurs
cAA
c====================================================================
c   Calcul  des tendances traceurs
c====================================================================
C

      call phytrac (
     I     itap,     days_elapsed+1,    jH_cur,   debut,
     I     lafin,    dtime,     u, v,     t,
     I     paprs,    pplay,     pmfu,     pmfd, 
     I     pen_u,    pde_u,     pen_d,    pde_d,
     I     cdragh,   coefh,     fm_therm, entr_therm,
     I     u1,       v1,        ftsol,    pctsrf,
     I     rlat,     frac_impa, frac_nucl,rlon,
     I     presnivs, pphis,     pphi,     albsol1,
     I     qx(:,:,ivap),rhcl,   cldfra,   rneb, 
     I     diafra,   cldliq,    itop_con, ibas_con,
     I     pmflxr,   pmflxs,    prfl,     psfl,
     I     da,       phi,       mp,       upwd,     
     I     dnwd,     aerosol_couple,      flxmass_w,
     I     tau_aero, piz_aero,  cg_aero,  ccm,
     I     rfname,
     O     tr_seri)

      IF (offline) THEN

       IF (prt_level.ge.9)
     $    print*,'Attention on met a 0 les thermiques pour phystoke'
	 call phystokenc (
     I                   nlon,klev,pdtphys,rlon,rlat,
     I                   t,pmfu, pmfd, pen_u, pde_u, pen_d, pde_d,
     I                   fm_therm,entr_therm,
     I                   cdragh,coefh,u1,v1,ftsol,pctsrf,
     I                   frac_impa, frac_nucl,
     I                   pphis,airephy,dtime,itap,
     I                   rlon,rlat,qx(:,:,ivap),da,phi,mp,upwd,dnwd)


      ENDIF

c
c Calculer le transport de l'eau et de l'energie (diagnostique)
c
      CALL transp (paprs,zxtsol,
     e                   t_seri, q_seri, u_seri, v_seri, zphi,
     s                   ve, vq, ue, uq)
c
cIM global posePB BEG
      IF(1.EQ.0) THEN
c
      CALL transp_lay (paprs,zxtsol,
     e                   t_seri, q_seri, u_seri, v_seri, zphi,
     s                   ve_lay, vq_lay, ue_lay, uq_lay)
c
      ENDIF !(1.EQ.0) THEN
cIM global posePB END
c Accumuler les variables a stocker dans les fichiers histoire:
c
c+jld ec_conser
      DO k = 1, klev
      DO i = 1, klon
        ZRCPD = RCPD*(1.0+RVTMP2*q_seri(i,k))
        d_t_ec(i,k)=0.5/ZRCPD
     $      *(u(i,k)**2+v(i,k)**2-u_seri(i,k)**2-v_seri(i,k)**2)
      ENDDO
      ENDDO

      DO k = 1, klev
      DO i = 1, klon
        t_seri(i,k)=t_seri(i,k)+d_t_ec(i,k)
        d_t_ec(i,k) = d_t_ec(i,k)/dtime
       END DO 
      END DO 
c-jld ec_conser
cIM
      IF (ip_ebil_phy.ge.1) THEN 
        ztit='after physic'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,1,1,dtime
     e      , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
C     Comme les tendances de la physique sont ajoute dans la dynamique,
C     on devrait avoir que la variation d'entalpie par la dynamique
C     est egale a la variation de la physique au pas de temps precedent.
C     Donc la somme de ces 2 variations devrait etre nulle.

        call diagphy(airephy,ztit,ip_ebil_phy
     e      , topsw, toplw, solsw, sollw, sens
     e      , evap, rain_fall, snow_fall, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
C
      d_h_vcol_phy=d_h_vcol
C
      END IF 
C
c=======================================================================
c   SORTIES
c=======================================================================

cIM Interpolation sur les niveaux de pression du NMC
c   -------------------------------------------------
c
#include "calcul_STDlev.h"
      twriteSTD(:,:,1)=tsumSTD(:,:,1)
      qwriteSTD(:,:,1)=qsumSTD(:,:,1)
      rhwriteSTD(:,:,1)=rhsumSTD(:,:,1)
      phiwriteSTD(:,:,1)=phisumSTD(:,:,1)
      uwriteSTD(:,:,1)=usumSTD(:,:,1)
      vwriteSTD(:,:,1)=vsumSTD(:,:,1)
      wwriteSTD(:,:,1)=wsumSTD(:,:,1)

      twriteSTD(:,:,2)=tsumSTD(:,:,2)
      qwriteSTD(:,:,2)=qsumSTD(:,:,2)
      rhwriteSTD(:,:,2)=rhsumSTD(:,:,2)
      phiwriteSTD(:,:,2)=phisumSTD(:,:,2)
      uwriteSTD(:,:,2)=usumSTD(:,:,2)
      vwriteSTD(:,:,2)=vsumSTD(:,:,2)
      wwriteSTD(:,:,2)=wsumSTD(:,:,2)

      twriteSTD(:,:,3)=tlevSTD(:,:)
      qwriteSTD(:,:,3)=qlevSTD(:,:)
      rhwriteSTD(:,:,3)=rhlevSTD(:,:)
      phiwriteSTD(:,:,3)=philevSTD(:,:)
      uwriteSTD(:,:,3)=ulevSTD(:,:)
      vwriteSTD(:,:,3)=vlevSTD(:,:)
      wwriteSTD(:,:,3)=wlevSTD(:,:)

      twriteSTD(:,:,4)=tlevSTD(:,:)
      qwriteSTD(:,:,4)=qlevSTD(:,:)
      rhwriteSTD(:,:,4)=rhlevSTD(:,:)
      phiwriteSTD(:,:,4)=philevSTD(:,:)
      uwriteSTD(:,:,4)=ulevSTD(:,:)
      vwriteSTD(:,:,4)=vlevSTD(:,:)
      wwriteSTD(:,:,4)=wlevSTD(:,:)
c
cIM initialisation 5eme fichier de sortie 
cIM ajoute 5eme niveau 170310 BEG
      twriteSTD(:,:,5)=tlevSTD(:,:)
      qwriteSTD(:,:,5)=qlevSTD(:,:)
      rhwriteSTD(:,:,5)=rhlevSTD(:,:)
      phiwriteSTD(:,:,5)=philevSTD(:,:)
      uwriteSTD(:,:,5)=ulevSTD(:,:)
      vwriteSTD(:,:,5)=vlevSTD(:,:)
      wwriteSTD(:,:,5)=wlevSTD(:,:)
cIM for NMC files
      DO n=1, nlevSTD3
       DO k=1, nlevSTD
        if(rlevSTD3(n).EQ.rlevSTD(k)) THEN
         twriteSTD3(:,n)=tlevSTD(:,k)
         qwriteSTD3(:,n)=qlevSTD(:,k)
         rhwriteSTD3(:,n)=rhlevSTD(:,k)
         phiwriteSTD3(:,n)=philevSTD(:,k)
         uwriteSTD3(:,n)=ulevSTD(:,k)
         vwriteSTD3(:,n)=vlevSTD(:,k)
         wwriteSTD3(:,n)=wlevSTD(:,k)
        endif !rlevSTD3(n).EQ.rlevSTD(k)
       ENDDO 
      ENDDO 
c
      DO n=1, nlevSTD8
       DO k=1, nlevSTD
        if(rlevSTD8(n).EQ.rlevSTD(k)) THEN
         tnondefSTD8(:,n)=tnondef(:,k,2)
         twriteSTD8(:,n)=tsumSTD(:,k,2)
         qwriteSTD8(:,n)=qsumSTD(:,k,2)
         rhwriteSTD8(:,n)=rhsumSTD(:,k,2)
         phiwriteSTD8(:,n)=phisumSTD(:,k,2)
         uwriteSTD8(:,n)=usumSTD(:,k,2)
         vwriteSTD8(:,n)=vsumSTD(:,k,2)
         wwriteSTD8(:,n)=wsumSTD(:,k,2)
        endif !rlevSTD8(n).EQ.rlevSTD(k)
       ENDDO 
      ENDDO 
c
c slp sea level pressure
      slp(:) = paprs(:,1)*exp(pphis(:)/(RD*t_seri(:,1)))
c
ccc prw = eau precipitable
      DO i = 1, klon
       prw(i) = 0.
       DO k = 1, klev
        prw(i) = prw(i) +
     .           q_seri(i,k)*(paprs(i,k)-paprs(i,k+1))/RG
       ENDDO
      ENDDO
c
cIM initialisation + calculs divers diag AMIP2
c
#include "calcul_divers.h"
c
      IF (config_inca /= 'none') THEN
#ifdef INCA
         CALL VTe(VTphysiq)
         CALL VTb(VTinca)

         CALL chemhook_end (
     $                        dtime,
     $                        pplay,
     $                        t_seri,
     $                        tr_seri,
     $                        nbtr,
     $                        paprs,
     $                        q_seri,
     $                        airephy,
     $                        pphi,
     $                        pphis,
     $                        zx_rh)

         CALL VTe(VTinca)
         CALL VTb(VTphysiq)
#endif
      END IF

c=============================================================
c
c Convertir les incrementations en tendances
c
      if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
      endif

      DO k = 1, klev
      DO i = 1, klon
         d_u(i,k) = ( u_seri(i,k) - u(i,k) ) / dtime
         d_v(i,k) = ( v_seri(i,k) - v(i,k) ) / dtime
         d_t(i,k) = ( t_seri(i,k)-t(i,k) ) / dtime
         d_qx(i,k,ivap) = ( q_seri(i,k) - qx(i,k,ivap) ) / dtime
         d_qx(i,k,iliq) = ( ql_seri(i,k) - qx(i,k,iliq) ) / dtime
      ENDDO
      ENDDO
c
      IF (nqtot.GE.3) THEN
      DO iq = 3, nqtot
      DO  k = 1, klev
      DO  i = 1, klon
         d_qx(i,k,iq) = ( tr_seri(i,k,iq-2) - qx(i,k,iq) ) / dtime
      ENDDO
      ENDDO
      ENDDO
      ENDIF
c
cIM rajout diagnostiques bilan KP pour analyse MJO par Jun-Ichi Yano
cIM global posePB#include "write_bilKP_ins.h"
cIM global posePB#include "write_bilKP_ave.h"
c

c Sauvegarder les valeurs de t et q a la fin de la physique:
c
      DO k = 1, klev
      DO i = 1, klon
         u_ancien(i,k) = u_seri(i,k)
         v_ancien(i,k) = v_seri(i,k)
         t_ancien(i,k) = t_seri(i,k)
         q_ancien(i,k) = q_seri(i,k)
      ENDDO
      ENDDO
c
!==========================================================================
! Sorties des tendances pour un point particulier
! a utiliser en 1D, avec igout=1 ou en 3D sur un point particulier
! pour le debug
! La valeur de igout est attribuee plus haut dans le programme
!==========================================================================

      if (prt_level.ge.1) then
      write(lunout,*) 'FIN DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
      write(lunout,*)
     s 'nlon,klev,nqtot,debut,lafin,jD_cur, jH_cur, pdtphys pct tlos'
      write(lunout,*)
     s  nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur ,pdtphys,
     s  pctsrf(igout,is_ter), pctsrf(igout,is_lic),pctsrf(igout,is_oce),
     s  pctsrf(igout,is_sic)
      write(lunout,*) 'd_t_dyn,d_t_con,d_t_lsc,d_t_ajsb,d_t_ajs,d_t_eva'
      do k=1,klev
         write(lunout,*) d_t_dyn(igout,k),d_t_con(igout,k),
     s   d_t_lsc(igout,k),d_t_ajsb(igout,k),d_t_ajs(igout,k),
     s   d_t_eva(igout,k)
      enddo
      write(lunout,*) 'cool,heat'
      do k=1,klev
         write(lunout,*) cool(igout,k),heat(igout,k)
      enddo

      write(lunout,*) 'd_t_oli,d_t_vdf,d_t_oro,d_t_lif,d_t_ec'
      do k=1,klev
         write(lunout,*) d_t_oli(igout,k),d_t_vdf(igout,k),
     s d_t_oro(igout,k),d_t_lif(igout,k),d_t_ec(igout,k)
      enddo

      write(lunout,*) 'd_ps ',d_ps(igout)
      write(lunout,*) 'd_u, d_v, d_t, d_qx1, d_qx2 '
      do k=1,klev
         write(lunout,*) d_u(igout,k),d_v(igout,k),d_t(igout,k),
     s  d_qx(igout,k,1),d_qx(igout,k,2)
      enddo
      endif

!==========================================================================

c============================================================
c   Calcul de la temperature potentielle
c============================================================
      DO k = 1, klev
      DO i = 1, klon
cJYG/IM theta en debut du pas de temps
cJYG/IM       theta(i,k)=t(i,k)*(100000./pplay(i,k))**(RD/RCPD)
cJYG/IM theta en fin de pas de temps de physique
        theta(i,k)=t_seri(i,k)*(100000./pplay(i,k))**(RD/RCPD)
      ENDDO
      ENDDO
c

c 22.03.04 BEG
c=============================================================
c   Ecriture des sorties
c=============================================================
#ifdef CPP_IOIPSL
 
c Recupere des varibles calcule dans differents modules
c pour ecriture dans histxxx.nc 

      ! Get some variables from module fonte_neige_mod
      CALL fonte_neige_get_vars(pctsrf, 
     .     zxfqcalving, zxfqfonte, zxffonte)


#include "phys_output_write.h"

#ifdef histISCCP
#include "write_histISCCP.h"
#endif

#ifdef histNMC
#include "write_histhfNMC.h"
#include "write_histdayNMC.h"
#include "write_histmthNMC.h"
#endif

#include "write_histday_seri.h"

#include "write_paramLMDZ_phy.h"

#endif

c 22.03.04 END
c
c====================================================================
c Si c'est la fin, il faut conserver l'etat de redemarrage
c====================================================================
c

c        -----------------------------------------------------------------
c        WSTATS: Saving statistics
c        -----------------------------------------------------------------
c        ("stats" stores and accumulates 8 key variables in file "stats.nc"
c        which can later be used to make the statistic files of the run:
c        "stats")          only possible in 3D runs !

         
         IF (callstats) THEN

           call wstats(klon,o_psol%name,"Surface pressure","Pa"
     &                 ,2,paprs(:,1))
           call wstats(klon,o_tsol%name,"Surface temperature","K",
     &                 2,zxtsol)
           zx_tmp_fi2d(:) = rain_fall(:) + snow_fall(:)
           call wstats(klon,o_precip%name,"Precip Totale liq+sol",
     &                 "kg/(s*m2)",2,zx_tmp_fi2d)
           zx_tmp_fi2d(:) = rain_lsc(:) + snow_lsc(:)
           call wstats(klon,o_plul%name,"Large-scale Precip",
     &                 "kg/(s*m2)",2,zx_tmp_fi2d)
           zx_tmp_fi2d(:) = rain_con(:) + snow_con(:)
           call wstats(klon,o_pluc%name,"Convective Precip",
     &                 "kg/(s*m2)",2,zx_tmp_fi2d)
           call wstats(klon,o_sols%name,"Solar rad. at surf.",
     &                 "W/m2",2,solsw)
           call wstats(klon,o_soll%name,"IR rad. at surf.",
     &                 "W/m2",2,sollw)
          zx_tmp_fi2d(:) = topsw(:)-toplw(:)
          call wstats(klon,o_nettop%name,"Net dn radiatif flux at TOA",
     &                 "W/m2",2,zx_tmp_fi2d)



           call wstats(klon,o_temp%name,"Air temperature","K",
     &                 3,t_seri)
           call wstats(klon,o_vitu%name,"Zonal wind","m.s-1",
     &                 3,u_seri)
           call wstats(klon,o_vitv%name,"Meridional wind",
     &                "m.s-1",3,v_seri)
           call wstats(klon,o_vitw%name,"Vertical wind",
     &                "m.s-1",3,omega)
           call wstats(klon,o_ovap%name,"Specific humidity", "kg/kg",
     &                 3,q_seri)
 


           IF(lafin) THEN
             write (*,*) "Writing stats..."
             call mkstats(ierr)
           ENDIF

         ENDIF !if callstats
     

      IF (lafin) THEN
         itau_phy = itau_phy + itap
         CALL phyredem ("restartphy.nc")
!         open(97,form="unformatted",file="finbin")
!         write(97) u_seri,v_seri,t_seri,q_seri
!         close(97)
C$OMP MASTER
         if (read_climoz >= 1) then
            if (is_mpi_root) then
               call nf95_close(ncid_climoz)
            end if
            deallocate(press_climoz) ! pointer
         end if
C$OMP END MASTER
      ENDIF
      
!      first=.false.

      RETURN
      END
      FUNCTION qcheck(klon,klev,paprs,q,ql,aire)
      IMPLICIT none
c
c Calculer et imprimer l'eau totale. A utiliser pour verifier
c la conservation de l'eau
c
#include "YOMCST.h"
      INTEGER klon,klev
      REAL paprs(klon,klev+1), q(klon,klev), ql(klon,klev)
      REAL aire(klon)
      REAL qtotal, zx, qcheck
      INTEGER i, k
c
      zx = 0.0
      DO i = 1, klon
         zx = zx + aire(i)
      ENDDO
      qtotal = 0.0
      DO k = 1, klev
      DO i = 1, klon
         qtotal = qtotal + (q(i,k)+ql(i,k)) * aire(i)
     .                     *(paprs(i,k)-paprs(i,k+1))/RG
      ENDDO
      ENDDO
c
      qcheck = qtotal/zx
c
      RETURN
      END
      SUBROUTINE gr_fi_ecrit(nfield,nlon,iim,jjmp1,fi,ecrit)
      IMPLICIT none
c
c Tranformer une variable de la grille physique a
c la grille d'ecriture
c
      INTEGER nfield,nlon,iim,jjmp1, jjm
      REAL fi(nlon,nfield), ecrit(iim*jjmp1,nfield)
c
      INTEGER i, n, ig
c
      jjm = jjmp1 - 1
      DO n = 1, nfield
         DO i=1,iim
            ecrit(i,n) = fi(1,n)
            ecrit(i+jjm*iim,n) = fi(nlon,n)
         ENDDO
         DO ig = 1, nlon - 2
           ecrit(iim+ig,n) = fi(1+ig,n)
         ENDDO
      ENDDO
      RETURN
      END
