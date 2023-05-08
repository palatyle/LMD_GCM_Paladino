!
! $Id: $
!
      MODULE physiq_mod

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE physiq (nlon,nlev,nqmax,
     .            debut,lafin,rjourvrai,gmtime,pdtphys,
     .            paprs,pplay,ppk,pphi,pphis,presnivs,
     .            u,v,t,qx,
     .            flxmw,
     .            d_u, d_v, d_t, d_qx, d_ps)

c======================================================================
c
c Modifications pour la physique de Venus
c     S. Lebonnois (LMD/CNRS) Septembre 2005
c
c ---------------------------------------------------------------------
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c
c Objet: Moniteur general de la physique du modele
cAA      Modifications quant aux traceurs :
cAA                  -  uniformisation des parametrisations ds phytrac
cAA                  -  stockage des moyennes des champs necessaires
cAA                     en mode traceur off-line 
c    modif   ( P. Le Van ,  12/10/98 )
c
c  Arguments:
c
c nlon----input-I-nombre de points horizontaux
c nlev----input-I-nombre de couches verticales
c nqmax---input-I-nombre de traceurs
c debut---input-L-variable logique indiquant le premier passage
c lafin---input-L-variable logique indiquant le dernier passage
c rjour---input-R-numero du jour de l'experience
c gmtime--input-R-fraction de la journee (0 a 1)
c pdtphys-input-R-pas d'integration pour la physique (seconde)
c paprs---input-R-pression pour chaque inter-couche (en Pa)
c pplay---input-R-pression pour le mileu de chaque couche (en Pa)
c ppk  ---input-R-fonction d'Exner au milieu de couche
c pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
c pphis---input-R-geopotentiel du sol
c presnivs-input_R_pressions approximat. des milieux couches ( en PA)
c u-------input-R-vitesse dans la direction X (de O a E) en m/s
c v-------input-R-vitesse Y (de S a N) en m/s
c t-------input-R-temperature (K)
c qx------input-R-mass mixing ratio traceurs (kg/kg) 
c d_t_dyn-input-R-tendance dynamique pour "t" (K/s)
c flxmw---input-R-flux de masse vertical en kg/s
c
c d_u-----output-R-tendance physique de "u" (m/s/s)
c d_v-----output-R-tendance physique de "v" (m/s/s)
c d_t-----output-R-tendance physique de "t" (K/s)
c d_qx----output-R-tendance physique de "qx" (kg/kg/s)
c d_ps----output-R-tendance physique de la pression au sol
c======================================================================
      USE ioipsl
!      USE histcom ! not needed; histcom is included in ioipsl
      use dimphy
      USE geometry_mod,only: longitude, latitude, ! in radians
     &                       longitude_deg,latitude_deg, ! in degrees
     &                       cell_area,dx,dy
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      USE cpdet_phy_mod, only: cpdet, t2tpot
      USE chemparam_mod
      USE conc
      USE compo_hedin83_mod2
!      use ieee_arithmetic
      use time_phylmdz_mod, only: annee_ref, day_ref, itau_phy
      use mod_grid_phy_lmdz, only: nbp_lon
      use infotrac_phy, only: iflag_trac, tname, ttext
      use vertical_layers_mod, only: pseudoalt
      use turb_mod, only : sens, turb_resolved
      use sed_and_prod_mad, only: aer_sedimentation, drop_sedimentation
#ifdef CPP_XIOS      
      use xios_output_mod, only: initialize_xios_output, 
     &                           update_xios_timestep, 
     &                           send_xios_field
      use wxios, only: wxios_context_init, xios_context_finalize
#endif
#ifdef MESOSCALE
      use comm_wrf
#else 
      use iophy
      use write_field_phy
      use mod_phys_lmdz_omp_data, ONLY: is_omp_master
      USE mod_phys_lmdz_para, only : is_parallel,jj_nb,
     &                               is_north_pole_phy,
     &                               is_south_pole_phy,
     &                               is_master
#endif
      IMPLICIT none
c======================================================================
c   CLEFS CPP POUR LES IO
c   =====================
#ifndef MESOSCALE
c#define histhf
#define histday
#define histmth
#define histins
#endif
c======================================================================
#include "dimsoil.h"
#include "clesphys.h"
#include "iniprint.h"
#include "timerad.h" 
#include "tabcontrol.h"
#include "nirdata.h"
#include "hedin.h"
c======================================================================
      LOGICAL ok_journe ! sortir le fichier journalier
      save ok_journe
c      PARAMETER (ok_journe=.true.)
c
      LOGICAL ok_mensuel ! sortir le fichier mensuel
      save ok_mensuel
c      PARAMETER (ok_mensuel=.true.)
c
      LOGICAL ok_instan ! sortir le fichier instantane
      save ok_instan
c      PARAMETER (ok_instan=.true.)
c
c======================================================================
c
c Variables argument:
c
      INTEGER nlon
      INTEGER nlev
      INTEGER nqmax
      REAL rjourvrai
      REAL gmtime
      REAL pdtphys
      LOGICAL debut, lafin
      REAL paprs(klon,klev+1)
      REAL pplay(klon,klev)
      REAL pphi(klon,klev)
      REAL pphis(klon)
      REAL presnivs(klev)

! ADAPTATION GCM POUR CP(T)
      REAL ppk(klon,klev)

      REAL u(klon,klev)
      REAL v(klon,klev)
      REAL t(klon,klev)
      REAL qx(klon,klev,nqmax)

      REAL d_u_dyn(klon,klev)
      REAL d_t_dyn(klon,klev)

      REAL flxmw(klon,klev)

      REAL d_u(klon,klev)
      REAL d_v(klon,klev)
      REAL d_t(klon,klev)
      REAL d_qx(klon,klev,nqmax)
      REAL d_ps(klon)

      logical ok_hf
      real ecrit_hf
      integer nid_hf
      save ok_hf, ecrit_hf, nid_hf

#ifdef histhf
      data ok_hf,ecrit_hf/.true.,0.25/
#else
      data ok_hf/.false./
#endif

c Variables propres a la physique

      integer,save :: itap        ! physics counter
      REAL delp(klon,klev)        ! epaisseur d'une couche
      REAL omega(klon,klev)       ! vitesse verticale en Pa/s
      
      INTEGER igwd,idx(klon),itest(klon)

c  Diagnostiques 2D de drag_noro, lift_noro et gw_nonoro

      REAL zulow(klon),zvlow(klon)
      REAL zustrdr(klon), zvstrdr(klon)
      REAL zustrli(klon), zvstrli(klon)
      REAL zustrhi(klon), zvstrhi(klon)
      REAL zublstrdr(klon), zvblstrdr(klon)
      REAL znlow(klon), zeff(klon)
      REAL zbl(klon), knu2(klon),kbreak(nlon)
      REAL tau0(klon), ztau(klon,klev)

c Pour calcul GW drag oro et nonoro: CALCUL de N2:
      real zdtlev(klon,klev),zdzlev(klon,klev)
      real ztlev(klon,klev)
      real zn2(klon,klev) ! BV^2 at plev

c Pour les bilans de moment angulaire, 
      integer bilansmc
c Pour le transport de ballons
      integer ballons
c j'ai aussi besoin
c du stress de couche limite a la surface:

      REAL zustrcl(klon),zvstrcl(klon) 

c et du stress total c de la physique:

      REAL zustrph(klon),zvstrph(klon) 

c Variables locales:
c
      REAL cdragh(klon) ! drag coefficient pour T and Q
      REAL cdragm(klon) ! drag coefficient pour vent
c
cAA  Pour  TRACEURS
cAA
      REAL,save,allocatable :: source(:,:)
      REAL ycoefh(klon,klev)    ! coef d'echange pour phytrac
      REAL yu1(klon)            ! vents dans la premiere couche U
      REAL yv1(klon)            ! vents dans la premiere couche V

      REAL dsens(klon) ! derivee chaleur sensible 
      REAL ve(klon) ! integr. verticale du transport meri. de l'energie
      REAL vq(klon) ! integr. verticale du transport meri. de l'eau
      REAL ue(klon) ! integr. verticale du transport zonal de l'energie
      REAL uq(klon) ! integr. verticale du transport zonal de l'eau
c
      REAL Fsedim(klon,klev+1)  ! Flux de sedimentation (kg.m-2)

c======================================================================
c
c Declaration des procedures appelees
c
      EXTERNAL ajsec     ! ajustement sec
      EXTERNAL clmain    ! couche limite 
      EXTERNAL clmain_ideal    ! couche limite simple 
      EXTERNAL hgardfou  ! verifier les temperatures
c     EXTERNAL orbite    ! calculer l'orbite 
      EXTERNAL phyetat0  ! lire l'etat initial de la physique
      EXTERNAL phyredem  ! ecrire l'etat de redemarrage de la physique
      EXTERNAL radlwsw   ! rayonnements solaire et infrarouge
!      EXTERNAL suphec    ! initialiser certaines constantes
      EXTERNAL transp    ! transport total de l'eau et de l'energie
      EXTERNAL printflag
      EXTERNAL zenang
      EXTERNAL diagetpq
      EXTERNAL conf_phys
      EXTERNAL diagphy
      EXTERNAL mucorr
      EXTERNAL nirco2abs
      EXTERNAL nir_leedat
      EXTERNAL nltecool
      EXTERNAL nlte_tcool
      EXTERNAL nlte_setup
      EXTERNAL blendrad
      EXTERNAL nlthermeq
      EXTERNAL euvheat
      EXTERNAL param_read
      EXTERNAL param_read_e107
      EXTERNAL conduction
      EXTERNAL molvis
      EXTERNAL moldiff_red

c
c Variables locales
c
CXXX PB 
      REAL fluxt(klon,klev)   ! flux turbulent de chaleur
      REAL fluxu(klon,klev)   ! flux turbulent de vitesse u
      REAL fluxv(klon,klev)   ! flux turbulent de vitesse v
c
      REAL flux_dyn(klon,klev)  ! flux de chaleur produit par la dynamique
      REAL flux_ajs(klon,klev)  ! flux de chaleur ajustement sec
      REAL flux_ec(klon,klev)   ! flux de chaleur Ec
c
      REAL tmpout(klon,klev)  ! K s-1 
c
      REAL dist, rmu0(klon), fract(klon)
      REAL zdtime, zlongi
c
      INTEGER i, k, iq, ig, j, ll, ilon, ilat, ilev, isoil
c
      REAL zphi(klon,klev)
      REAL zzlev(klon,klev+1),zzlay(klon,klev),z1,z2
      real tsurf(klon)

c va avec nlte_tcool
      INTEGER ierr_nlte
      REAL    varerr

! photochemistry

      integer :: chempas
      real :: nbapp_chem, zctime

! sedimentation

      REAL :: m0_mode1(klev,2),m0_mode2(klev,2)
      REAL :: m3_mode1(klev,3),m3_mode2 (klev,3)
      REAL :: d_drop_sed(klev),d_ccn_sed(klev,2),d_liq_sed(klev,2)
      REAL :: aer_flux(klev)
c
c Variables du changement
c
c ajs: ajustement sec
c vdf: couche limite (Vertical DiFfusion)
      REAL d_t_ajs(klon,klev), d_tr_ajs(klon,klev,nqmax)
      REAL d_u_ajs(klon,klev), d_v_ajs(klon,klev)
c
      REAL d_ts(klon)
c
      REAL d_u_vdf(klon,klev), d_v_vdf(klon,klev)
      REAL d_t_vdf(klon,klev), d_tr_vdf(klon,klev,nqmax)
c
CMOD LOTT: Tendances Orography Sous-maille
      REAL d_u_oro(klon,klev), d_v_oro(klon,klev)
      REAL d_t_oro(klon,klev)
      REAL d_u_lif(klon,klev), d_v_lif(klon,klev)
      REAL d_t_lif(klon,klev)
C          Tendances Ondes de G non oro (runs strato).
      REAL d_u_hin(klon,klev), d_v_hin(klon,klev)
      REAL d_t_hin(klon,klev)

c Tendencies due to radiative scheme   [K/s]
c     d_t_rad,dtsw,dtlw,d_t_nirco2,d_t_nlte,d_t_euv
c are not computed at each physical timestep
c therefore, they are defined and saved in phys_state_var_mod

c Tendencies due to molecular viscosity and conduction
      real d_t_conduc(klon,klev)     ! [K/s]
      real d_u_molvis(klon,klev)     ! (m/s) /s
      real d_v_molvis(klon,klev)     ! (m/s) /s

c Tendencies due to molecular diffusion
      real d_q_moldif(klon,klev,nqmax)

c
c Variables liees a l'ecriture de la bande histoire physique
c
      INTEGER ecrit_mth
      SAVE ecrit_mth   ! frequence d'ecriture (fichier mensuel)
c
      INTEGER ecrit_day
      SAVE ecrit_day   ! frequence d'ecriture (fichier journalier)
c
      INTEGER ecrit_ins
      SAVE ecrit_ins   ! frequence d'ecriture (fichier instantane)
c
      integer itau_w   ! pas de temps ecriture = itap + itau_phy

c Variables locales pour effectuer les appels en serie
c
      REAL t_seri(klon,klev)
      REAL u_seri(klon,klev), v_seri(klon,klev)
c
      REAL :: tr_seri(klon,klev,nqmax)
      REAL :: d_tr(klon,klev,nqmax)

c pour ioipsl
      INTEGER nid_day, nid_mth, nid_ins
      SAVE nid_day, nid_mth, nid_ins
      INTEGER nhori, nvert, idayref
      REAL zsto, zout, zsto1, zsto2, zero
      parameter (zero=0.0e0)
      real zjulian
      save zjulian

      CHARACTER*2  str2
      character*20 modname
      character*80 abort_message
      logical ok_sync

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
      REAL      d_h_vcol, d_h_dair, d_qt, d_qw, d_ql, d_qs, d_ec
      REAL      d_h_vcol_phy
      REAL      fs_bound, fq_bound
      SAVE      d_h_vcol_phy
      REAL      zero_v(klon),zero_v2(klon,klev)
      CHARACTER*15 ztit
      INTEGER   ip_ebil  ! PRINT level for energy conserv. diag.
      SAVE      ip_ebil
      DATA      ip_ebil/2/
      INTEGER   if_ebil ! level for energy conserv. dignostics
      SAVE      if_ebil
c+jld ec_conser
      REAL d_t_ec(klon,klev)    ! tendance du a la conversion Ec -> E thermique
c-jld ec_conser

c TEST VENUS...
      REAL mang(klon,klev)    ! moment cinetique
      REAL mangtot            ! moment cinetique total

c cell_area for outputs in hist*
      REAL cell_area_out(klon)
#ifdef MESOSCALE 
      REAL :: dt_dyn(klev)
#endif 
c Declaration des constantes et des fonctions thermodynamiques
c
#include "YOMCST.h"

c======================================================================
c INITIALISATIONS
c================

      modname = 'physiq'
      ok_sync=.TRUE.

      bilansmc = 0
      ballons  = 0
! NE FONCTIONNENT PAS ENCORE EN PARALLELE !!!
#ifndef MESOSCALE
      if (is_parallel) then
        bilansmc = 0
        ballons  = 0
      endif
#endif
      IF (if_ebil.ge.1) THEN
        DO i=1,klon
          zero_v(i)=0.
        END DO 
        DO i=1,klon
         DO j=1,klev
          zero_v2(i,j)=0.
         END DO 
        END DO 
      END IF 
      
c PREMIER APPEL SEULEMENT
c========================
      IF (debut) THEN
         allocate(source(klon,nqmax))

#ifdef CPP_XIOS
        ! Initialize XIOS context
        write(*,*) "physiq: call wxios_context_init"
        CALL wxios_context_init
#endif

! The call to suphec is now done in iniphysiq_mod (interface)
!         CALL suphec ! initialiser constantes et parametres phys.

         IF (if_ebil.ge.1) d_h_vcol_phy=0.
c
c appel a la lecture du physiq.def 
c
         call conf_phys(ok_journe, ok_mensuel,
     .                  ok_instan,
     .                  if_ebil)

         call phys_state_var_init(nqmax)
c
c Initialising Hedin model for upper atm 
c   (to be revised when coupled to chemistry) :
         call conc_init

! initialise physics counter

      itap    = 0

#ifdef MESOSCALE
      print*,'check pdtphys',pdtphys
      PRINT*,'check phisfi ',pphis(1),pphis(klon)
      PRINT*,'check geop',pphi(1,1),pphi(klon,klev)
      PRINT*,'check radsol',radsol(1),radsol(klon)
      print*,'check ppk',ppk(1,1),ppk(klon,klev)
      print*,'check ftsoil',ftsoil(1,1),ftsoil(klon,nsoilmx)
      print*,'check ftsol',ftsol(1),ftsol(klon)
      print*, "check temp", t(1,1),t(klon,klev)
      print*, "check pres",paprs(1,1),paprs(klon,klev),pplay(1,1),
     .                     pplay(klon,klev)
      print*, "check u", u(1,1),u(klon,klev)
      print*, "check v", v(1,1),v(klon,klev)
      print*,'check falbe',falbe(1),falbe(klon)
      !nqtot=nqmax
      !ALLOCATE(tname(nqtot)) 
      !tname=noms
      zmea=0.
      zstd=0.
      zsig=0.
      zgam=0.
      zthe=0.
      dtime=pdtphys
#else
c         
c Lecture startphy.nc :
c
         CALL phyetat0 ("startphy.nc")
         IF (.not.startphy_file) THEN
           ! Additionnal academic initializations
           ftsol(:)=t(:,1) ! surface temperature as in first atm. layer
           DO isoil=1, nsoilmx
             ! subsurface temperatures equal to surface temperature
             ftsoil(:,isoil)=ftsol(:)
           ENDDO
         ENDIF
#endif

c dtime est defini dans tabcontrol.h et lu dans startphy
c pdtphys est calcule a partir des nouvelles conditions:
c Reinitialisation du pas de temps physique quand changement
         IF (ABS(dtime-pdtphys).GT.0.001) THEN
            WRITE(lunout,*) 'Pas physique a change',dtime,
     .                        pdtphys
c           abort_message='Pas physique n est pas correct '
c           call abort_physic(modname,abort_message,1)
c----------------
c pour initialiser convenablement le time_counter, il faut tenir compte
c du changement de dtime en changeant itau_phy (point de depart)
            itau_phy = NINT(itau_phy*dtime/pdtphys)
c----------------
            dtime=pdtphys
         ENDIF

         radpas = NINT(RDAY/pdtphys/nbapp_rad)

         CALL printflag( ok_journe,ok_instan )

#ifdef CPP_XIOS
         write(*,*) "physiq: call initialize_xios_output"
         call initialize_xios_output(rjourvrai,gmtime,pdtphys,RDAY,
     &                               presnivs,pseudoalt)
#endif

c
c---------
c FLOTT
       IF (ok_orodr) THEN
         DO i=1,klon
         rugoro(i) = MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
         ENDDO
         CALL SUGWD(klon,klev,paprs,pplay)
         DO i=1,klon
         zuthe(i)=0.
         zvthe(i)=0.
         if(zstd(i).gt.10.)then
           zuthe(i)=(1.-zgam(i))*cos(zthe(i))
           zvthe(i)=(1.-zgam(i))*sin(zthe(i))
         endif
         ENDDO
       ENDIF

      if (bilansmc.eq.1) then
C  OUVERTURE D'UN FICHIER FORMATTE POUR STOCKER LES COMPOSANTES
C  DU BILAN DE MOMENT ANGULAIRE.
      open(27,file='aaam_bud.out',form='formatted')
      open(28,file='fields_2d.out',form='formatted')
      write(*,*)'Ouverture de aaam_bud.out (FL Vous parle)'
      write(*,*)'Ouverture de fields_2d.out (FL Vous parle)'
      endif !bilansmc

c--------------SLEBONNOIS
C  OUVERTURE DES FICHIERS FORMATTES CONTENANT LES POSITIONS ET VITESSES
C  DES BALLONS
      if (ballons.eq.1) then
      open(30,file='ballons-lat.out',form='formatted')
      open(31,file='ballons-lon.out',form='formatted')
      open(32,file='ballons-u.out',form='formatted')
      open(33,file='ballons-v.out',form='formatted')
      open(34,file='ballons-alt.out',form='formatted')
      write(*,*)'Ouverture des ballons*.out'
      endif !ballons
c-------------

c---------
C TRACEURS
C source dans couche limite
         source(:,:) = 0.0 ! pas de source, pour l'instant
c---------

c---------
c INITIALIZE THERMOSPHERIC PARAMETERS
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (callthermos) then
            if(solvarmod.eq.0) call param_read
            if(solvarmod.eq.1) call param_read_e107 
         endif

c Initialisation (recomputed in concentration2)
       do ig=1,klon
         do j=1,klev
            rnew(ig,j)=R
            cpnew(ig,j)=cpdet(t(ig,j))
            mmean(ig,j)=RMD
            akknew(ig,j)=1.e-4
            rho(ig,j)=pplay(ig,j)*mmean(ig,j)*1e-3/(rnew(ig,j)*t(ig,j))
          enddo
c        stop

        enddo  
     
      IF(callthermos.or.callnlte.or.callnirco2) THEN  
         call compo_hedin83_init2
      ENDIF
      if (callnlte.and.nltemodel.eq.2) call nlte_setup
      if (callnirco2.and.nircorr.eq.1) call nir_leedat         
c---------
      
c
c Verifications:
c
         IF (nlon .NE. klon) THEN
            WRITE(lunout,*)'nlon et klon ne sont pas coherents', nlon, 
     .                      klon
            abort_message='nlon et klon ne sont pas coherents'
            call abort_physic(modname,abort_message,1)
         ENDIF
         IF (nlev .NE. klev) THEN
            WRITE(lunout,*)'nlev et klev ne sont pas coherents', nlev,
     .                       klev
            abort_message='nlev et klev ne sont pas coherents'
            call abort_physic(modname,abort_message,1)
         ENDIF
c
         IF (dtime*REAL(radpas).GT.(RDAY*0.25).AND.cycle_diurne)
     $    THEN 
           WRITE(lunout,*)'Nbre d appels au rayonnement insuffisant'
           WRITE(lunout,*)"Au minimum 4 appels par jour si cycle diurne"
           abort_message='Nbre d appels au rayonnement insuffisant'
           call abort_physic(modname,abort_message,1)
         ENDIF
c
         WRITE(lunout,*)"Clef pour la convection seche, iflag_ajs=",
     .                   iflag_ajs
c
         ecrit_mth = NINT(RDAY/dtime*ecriphy)  ! tous les ecritphy jours
         IF (ok_mensuel) THEN
         WRITE(lunout,*)'La frequence de sortie mensuelle est de ', 
     .                   ecrit_mth
         ENDIF

         ecrit_day = NINT(RDAY/dtime *1.0)  ! tous les jours
         IF (ok_journe) THEN
         WRITE(lunout,*)'La frequence de sortie journaliere est de ',
     .                   ecrit_day
         ENDIF

         ecrit_ins = NINT(RDAY/dtime*ecriphy)  ! Fraction de jour reglable
         IF (ok_instan) THEN
         WRITE(lunout,*)'La frequence de sortie instant. est de ', 
     .                   ecrit_ins
         ENDIF

c Initialisation des sorties 
c===========================

#ifdef CPP_IOIPSL

#ifdef histhf
#include "ini_histhf.h"
#endif

#ifdef histday
#include "ini_histday.h"
#endif

#ifdef histmth
#include "ini_histmth.h"
#endif

#ifdef histins
#include "ini_histins.h"
#endif

#endif

c
c Initialiser les valeurs de u pour calculs tendances
c (pour T, c'est fait dans phyetat0)
c
      DO k = 1, klev
      DO i = 1, klon
         u_ancien(i,k) = u(i,k)
      ENDDO
      ENDDO

c---------
c       Ecriture fichier initialisation
c       PRINT*,'Ecriture Initial_State.csv'
c       OPEN(88,file='Trac_Point.csv',
c     & form='formatted')
c---------
     
c---------
c       Initialisation des parametres des nuages
c===============================================
     
c MICROPHY SANS CHIMIE: seulement si full microphy (cl_scheme=2)

      if (ok_chem.and..not.ok_cloud) then
        print*,"LA CHIMIE A BESOIN DE LA MICROPHYSIQUE"
        print*,"ok_cloud doit etre = a ok_chem"
      stop
      endif
      if (.not.ok_chem.and.ok_cloud.and.(cl_scheme.eq.1)) then
        print*,"cl_scheme=1 doesnot work without chemistry"
      stop
      endif
       if (.not.ok_chem.and.ok_cloud.and.(cl_scheme.eq.2)) then
        print*,"Full microphysics without chemistry"
c indexation of microphys tracers
        CALL chemparam_ini() 
      endif
    
! number of microphysical tracers

      nmicro = 0
      if (ok_cloud .and. (cl_scheme == 1)) nmicro = 2
      if (ok_cloud .and. (cl_scheme == 2)) nmicro = 12
 
c CAS 1D POUR MICROPHYS Aurelien
      if ((nlon .EQ. 1) .AND. ok_cloud .AND. (cl_scheme.eq.1)) then
        PRINT*,'Open profile_cloud_parameters.csv'
        OPEN(66,file='profile_cloud_parameters.csv',
     &   form='formatted')
      endif

      if ((nlon .EQ. 1) .AND. ok_sedim .AND. (cl_scheme.eq.1)) then
        PRINT*,'Open profile_cloud_sedim.csv'
        OPEN(77,file='profile_cloud_sedim.csv',
     &   form='formatted')
      endif
           
c INIT PHOTOCHEMISTRY ! includes the indexation of microphys tracers
c     if ((nlon .GT. 1) .AND. ok_chem) then
c !!! DONC 3D !!!  POURQUOI ???
      if (ok_chem) then
        CALL chemparam_ini() 
      endif
         
c INIT MICROPHYS SCHEME 1 (AURELIEN)  
      if ((nlon .GT. 1) .AND. ok_cloud .AND. (cl_scheme.eq.1)) then
c !!! DONC 3D !!!
        CALL cloud_ini(nlon,nlev)
      endif

!     initialise mmean

      if(callthermos) then
         call concentrations2(pplay,t,d_t,qx,nqmax)
      endif

      ENDIF ! debut

! ------------------------------------------------------
!   Initializations done at every physical timestep:
! ------------------------------------------------------

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
      DO iq = 1, nqmax
      DO k = 1, klev
      DO i = 1, klon
         d_qx(i,k,iq) = 0.0
      ENDDO
      ENDDO
      ENDDO
c
c Ne pas affecter les valeurs entrees de u, v, h, et q
c
      DO k = 1, klev
      DO i = 1, klon
         t_seri(i,k)  = t(i,k)
         u_seri(i,k)  = u(i,k)
         v_seri(i,k)  = v(i,k)
      ENDDO
      ENDDO
      DO iq = 1, nqmax
      DO  k = 1, klev
      DO  i = 1, klon
         tr_seri(i,k,iq) = qx(i,k,iq)
      ENDDO
      ENDDO
      ENDDO
C
      DO i = 1, klon
          ztsol(i) = ftsol(i)
      ENDDO
C
      IF (if_ebil.ge.1) THEN 
        ztit='after dynamic'
        CALL diagetpq(cell_area,ztit,ip_ebil,1,1,dtime
     e      , t_seri,zero_v2,zero_v2,zero_v2,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
C     Comme les tendances de la physique sont ajoute dans la dynamique,
C     on devrait avoir que la variation d'entalpie par la dynamique
C     est egale a la variation de la physique au pas de temps precedent.
C     Donc la somme de ces 2 variations devrait etre nulle.
        call diagphy(cell_area,ztit,ip_ebil
     e      , zero_v, zero_v, zero_v, zero_v, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol+d_h_vcol_phy, d_qt, 0.
     s      , fs_bound, fq_bound )
      END IF 

c====================================================================
c XIOS outputs

#ifdef CPP_XIOS      
      ! update XIOS time/calendar
      call update_xios_timestep
#endif      

c====================================================================
c Diagnostiquer la tendance dynamique
c
      IF (ancien_ok) THEN
         DO k = 1, klev
         DO i = 1, klon
            d_u_dyn(i,k) = (u_seri(i,k)-u_ancien(i,k))/dtime
            d_t_dyn(i,k) = (t_seri(i,k)-t_ancien(i,k))/dtime
         ENDDO
         ENDDO

! ADAPTATION GCM POUR CP(T)
         do i=1,klon
          flux_dyn(i,1) = 0.0
          do j=2,klev
            flux_dyn(i,j) = flux_dyn(i,j-1)
     . +cpdet(t_seri(i,j-1))/RG*d_t_dyn(i,j-1)*(paprs(i,j-1)-paprs(i,j))
          enddo
         enddo
         
      ELSE
         DO k = 1, klev
         DO i = 1, klon
            d_u_dyn(i,k) = 0.0
            d_t_dyn(i,k) = 0.0
         ENDDO
         ENDDO
         ancien_ok = .TRUE.
      ENDIF
c====================================================================

c Calcule de vitesse verticale a partir de flux de masse verticale
      DO k = 1, klev
       DO i = 1, klon
        omega(i,k) = RG*flxmw(i,k) / cell_area(i)
       END DO
      END DO

c======
c GEOP CORRECTION
c
c Ajouter le geopotentiel du sol:
c
      DO k = 1, klev
      DO i = 1, klon
         zphi(i,k) = pphi(i,k) + pphis(i)
      ENDDO
      ENDDO

c............................
c CETTE CORRECTION VA DE PAIR AVEC DES MODIFS DE LEAPFROG(_p)
c ELLE MARCHE A 50 NIVEAUX (si mmean constante...)
c MAIS PAS A 78 NIVEAUX (quand mmean varie...)
c A ANALYSER PLUS EN DETAIL AVANT D'UTILISER
c............................
c zphi is recomputed (pphi is not ok if mean molecular mass varies) 
c with     dphi = RT/mmean d(ln p) [evaluated at interface]

c     DO i = 1, klon
c       zphi(i,1) = pphis(i) + R*t_seri(i,1)/mmean(i,1)*1000.
c    *                *( log(paprs(i,1)) - log(pplay(i,1)) )      
c       DO k = 2, klev
c        zphi(i,k) = zphi(i,k-1) 
c    *      + R*500.*(t_seri(i,k)/mmean(i,k)+t_seri(i,k-1)/mmean(i,k-1))
c    *          * (log(pplay(i,k-1)) - log(pplay(i,k)))
c       ENDDO
c     ENDDO
c............................
c=====

c   calcul du geopotentiel aux niveaux intercouches
c   ponderation des altitudes au niveau des couches en dp/p

      DO k=1,klev
         DO i=1,klon
            zzlay(i,k)=zphi(i,k)/RG        ! [m]
         ENDDO
      ENDDO
      DO i=1,klon
         zzlev(i,1)=pphis(i)/RG            ! [m]
      ENDDO
      DO k=2,klev
         DO i=1,klon
            z1=(pplay(i,k-1)+paprs(i,k))/(pplay(i,k-1)-paprs(i,k))
            z2=(paprs(i,k)  +pplay(i,k))/(paprs(i,k)  -pplay(i,k))
            zzlev(i,k)=(z1*zzlay(i,k-1)+z2*zzlay(i,k))/(z1+z2)
         ENDDO
      ENDDO
      DO i=1,klon
         zzlev(i,klev+1)=zzlay(i,klev)+(zzlay(i,klev)-zzlev(i,klev))
      ENDDO

c====================================================================
c
c Verifier les temperatures
c
      CALL hgardfou(t_seri,ftsol,'debutphy')

c====================================================================
c Orbite et eclairement
c====================================================================

c Pour VENUS, on fixe l'obliquite a 0 et l'eccentricite a 0.
c donc pas de variations de Ls, ni de dist.
c La seule chose qui compte, c'est la rotation de la planete devant
c le Soleil...
      
      zlongi = 0.0
      dist   = 0.72  ! en UA

c Si on veut remettre l'obliquite a 3 degres et/ou l'eccentricite 
c a sa valeur, et prendre en compte leur evolution, 
c il faudra refaire un orbite.F...
c     CALL orbite(zlongi,dist)

      IF (cycle_diurne) THEN
        zdtime=dtime*REAL(radpas) ! pas de temps du rayonnement (s)
        CALL zenang(zlongi,gmtime,zdtime,latitude_deg,longitude_deg,
     &              rmu0,fract)
      ELSE
        call mucorr(klon,zlongi,latitude_deg,rmu0,fract)
      ENDIF
      
!     print fraction of venus day

      if (is_master) then
         print*, 'gmtime = ', gmtime
      end if

      if (iflag_trac.eq.1) then
!====================================================================
! Case 1: pseudo-chemistry with relaxation toward fixed profile
!====================================================================
       if (tr_scheme.eq.1) then

         call phytrac_relax (debut,lafin,nqmax,
     I                   nlon,nlev,dtime,pplay,
     O                   tr_seri)

       elseif (tr_scheme.eq.2) then
!====================================================================
! Case 2: surface emission
! For the moment, inspired from Mars version
! However, the variable 'source' could be used in physiq
! so the call to phytrac_emiss could be to initialise it.
!====================================================================

         call phytrac_emiss ( (rjourvrai+gmtime)*RDAY,
     I                   debut,lafin,nqmax,
     I                   nlon,nlev,dtime,paprs,
     I                   latitude_deg,longitude_deg,
     O                   tr_seri)

      else if (tr_scheme.eq.3) then
!====================================================================
! Case 3: Full chemistry and/or clouds.
!         routines are called every "chempas" physical timestep.
!
!         if the physics is called 96000 times per venus day:
!
!         nbapp_chem = 24000 => chempas = 4 => zctime = 420 s
!         nbapp_chem = 12000 => chempas = 8 => zctime = 840 s
!====================================================================

         nbapp_chem = 24000
         chempas = nint(rday/pdtphys/nbapp_chem)
         zctime = dtime*real(chempas)             ! chemical timestep

         if (mod(itap,chempas) == 0) then         ! <------- start of chemistry supercycling

!        photochemistry and microphysics

         call phytrac_chimie(debut,
     $                       gmtime,
     $                       nqmax,
     $                       klon,
     $                       latitude_deg,
     $                       longitude_deg,
     $                       nlev,
     $                       zctime,
     $                       t_seri,
     $                       pplay,
     $                       tr_seri,
     $                       d_tr_chem,
     $                       iter)

         if (ok_sedim) then 
            if (cl_scheme == 1) then 

!        sedimentation for simplified microphysics

#ifndef MESOSCALE
               call new_cloud_sedim(klon,
     $                              nlev,
     $                              zctime,
     $                              pplay,
     $                              paprs,
     $                              t_seri,
     $                              tr_seri,
     $                              d_tr_chem,
     $                              d_tr_sed(:,:,1:2),
     $                              nqmax,
     $                              Fsedim)

!        test to avoid nans 

               do k = 1, klev
                  do i = 1, klon
                     if ((d_tr_sed(i,k,1) /= d_tr_sed(i,k,1)) .or.
     $                   (d_tr_sed(i,k,2) /= d_tr_sed(i,k,2))) then
                        print*,'sedim NaN PROBLEM'
                        print*,'d_tr_sed Nan?',d_tr_sed(i,k,:)
                        print*,'Temp',t_seri(i,k)
                        print*,'lat-lon',i,'level',k,'zctime',zctime
                        print*,'F_sed',Fsedim(i,k)
                        d_tr_sed(i,k,:) = 0.
                     end if
                  end do
               end do

!        tendency due to condensation and sedimentation 

               d_tr_sed(:,:,1:2) = d_tr_sed(:,:,1:2)/zctime
               Fsedim(:,1:klev) = Fsedim(:,1:klev)/zctime
               Fsedim(:,klev+1) = 0.

            else if (cl_scheme == 2) then

!        sedimentation for detailed microphysics

               d_tr_sed(:,:,:) = 0.

               do i = 1, klon

!                 mode 1

                  m0_mode1(:,1) = tr_seri(i,:,i_m0_mode1drop)
                  m0_mode1(:,2) = tr_seri(i,:,i_m0_mode1ccn)
                  m3_mode1(:,1) = tr_seri(i,:,i_m3_mode1sa)
                  m3_mode1(:,2) = tr_seri(i,:,i_m3_mode1w)
                  m3_mode1(:,3) = tr_seri(i,:,i_m3_mode1ccn)

                  call drop_sedimentation(zctime,klev,m0_mode1,m3_mode1,
     $                                    paprs(i,:),zzlev(i,:),
     $                                    zzlay(i,:),t_seri(i,:),1,
     $                                    d_ccn_sed(:,1),d_drop_sed,
     $                                    d_ccn_sed(:,2),d_liq_sed)

        d_tr_sed(i,:,i_m0_mode1drop)= d_tr_sed(i,:,i_m0_mode1drop)
     $                              + d_drop_sed
        d_tr_sed(i,:,i_m0_mode1ccn) = d_tr_sed(i,:,i_m0_mode1ccn)
     $                              + d_ccn_sed(:,1)
        d_tr_sed(i,:,i_m3_mode1ccn) = d_tr_sed(i,:,i_m3_mode1ccn)
     $                              + d_ccn_sed(:,2)
        d_tr_sed(i,:,i_m3_mode1sa)  = d_tr_sed(i,:,i_m3_mode1sa)
     $                              + d_liq_sed(:,1)
        d_tr_sed(i,:,i_m3_mode1w)   = d_tr_sed(i,:,i_m3_mode1w)
     $                              + d_liq_sed(:,2)

!                 mode 2

                  m0_mode2(:,1) = tr_seri(i,:,i_m0_mode2drop)
                  m0_mode2(:,2) = tr_seri(i,:,i_m0_mode2ccn)
                  m3_mode2(:,1) = tr_seri(i,:,i_m3_mode2sa)
                  m3_mode2(:,2) = tr_seri(i,:,i_m3_mode2w)
                  m3_mode2(:,3) = tr_seri(i,:,i_m3_mode2ccn)

                  call drop_sedimentation(zctime,klev,m0_mode2,m3_mode2,
     $                                    paprs(i,:),zzlev(i,:),
     $                                    zzlay(i,:),t_seri(i,:),2,
     $                                    d_ccn_sed(:,1),d_drop_sed,
     $                                    d_ccn_sed(:,2),d_liq_sed)

        d_tr_sed(i,:,i_m0_mode2drop)= d_tr_sed(i,:,i_m0_mode2drop)
     $                              + d_drop_sed
        d_tr_sed(i,:,i_m0_mode2ccn) = d_tr_sed(i,:,i_m0_mode2ccn)
     $                              + d_ccn_sed(:,1)
        d_tr_sed(i,:,i_m3_mode2ccn) = d_tr_sed(i,:,i_m3_mode2ccn)
     $                              + d_ccn_sed(:,2)
        d_tr_sed(i,:,i_m3_mode2sa)  = d_tr_sed(i,:,i_m3_mode2sa)
     $                              + d_liq_sed(:,1)
        d_tr_sed(i,:,i_m3_mode2w)   = d_tr_sed(i,:,i_m3_mode2w)
     $                              + d_liq_sed(:,2)

!                 aer

                  call aer_sedimentation(zctime,klev,
     $                                   tr_seri(i,:,i_m0_aer),
     $                                   tr_seri(i,:,i_m3_aer),
     $                                   paprs(i,:),zzlev(i,:),
     $                                   zzlay(i,:),t_seri(i,:),
     $                                   d_tr_sed(i,:,i_m0_aer),
     $                                   d_tr_sed(i,:,i_m3_aer),
     $                                   aer_flux)

               end do
          
!        tendency due to sedimentation 

               do iq = nqmax-nmicro+1,nqmax
                  d_tr_sed(:,:,iq) = d_tr_sed(:,:,iq)/zctime
               end do
#endif
            end if  ! cl_scheme

!        update gaseous tracers (chemistry)

            do iq = 1, nqmax - nmicro
               tr_seri(:,:,iq) = tr_seri(:,:,iq) 
     $                         + d_tr_chem(:,:,iq)*zctime
            end do

!        update condensed tracers (condensation + sedimentation)

            if (cl_scheme == 1) then 
               tr_seri(:,:,i_h2so4liq) = max(tr_seri(:,:,i_h2so4liq) 
     $                                 + d_tr_sed(:,:,1)*zctime, 1.e-30)
               tr_seri(:,:,i_h2oliq)   = max(tr_seri(:,:,i_h2oliq) 
     $                                 + d_tr_sed(:,:,2)*zctime, 1.e-30)
            else if (cl_scheme == 2) then
               do iq = nqmax-nmicro+1,nqmax
                  tr_seri(:,:,iq) = tr_seri(:,:,iq) 
     $                            + d_tr_sed(:,:,iq)*zctime
               end do
            end if  ! cl_scheme

         end if     ! ok_sedim
         end if     ! mod(itap,chempas)  <------- end of chemistry supercycling

!====================================================================
! End Case 3: Full chemistry and/or clouds.
!====================================================================

         end if     ! tr_scheme
      end if        ! iflag_trac

c====================================================================
c Appeler la diffusion verticale (programme de couche limite)
c====================================================================

c-------------------------------
c VENUS TEST: on ne tient pas compte des calculs de clmain mais on force
c l'equilibre radiatif du sol
      if (.not. ok_clmain) then
              if (debut) then
                print*,"ATTENTION, CLMAIN SHUNTEE..."
              endif

      DO i = 1, klon
         sens(i) = 0.0e0 ! flux de chaleur sensible au sol
         fder(i) = 0.0e0
         dlw(i)  = 0.0e0
      ENDDO

c Incrementer la temperature du sol
c
      DO i = 1, klon
         d_ts(i)  = dtime * radsol(i)/22000. !valeur calculee par GCM pour I=200
         ftsol(i) = ftsol(i) + d_ts(i)
         do j=1,nsoilmx
           ftsoil(i,j)=ftsol(i)
         enddo
      ENDDO

c-------------------------------
      else
c-------------------------------

      fder = dlw

! ADAPTATION GCM POUR CP(T)

      if (physideal) then
       CALL clmain_ideal(dtime,itap,
     e            t_seri,u_seri,v_seri,
     e            rmu0, 
     e            ftsol,
     $            ftsoil,
     $            paprs,pplay,ppk,radsol,falbe,
     e            solsw, sollw, sollwdown, fder,
     e            longitude_deg, latitude_deg, dx, dy,   
     e            debut, lafin,
     s            d_t_vdf,d_u_vdf,d_v_vdf,d_ts,
     s            fluxt,fluxu,fluxv,cdragh,cdragm,
     s            dsens,
     s            ycoefh,yu1,yv1) 
      else
       CALL clmain(dtime,itap,
     e            t_seri,u_seri,v_seri,
     e            rmu0, 
     e            ftsol,
     $            ftsoil,
     $            paprs,pplay,ppk,radsol,falbe,
     e            solsw, sollw, sollwdown, fder,
     e            longitude_deg, latitude_deg, dx, dy,   
     e            debut, lafin,
     s            d_t_vdf,d_u_vdf,d_v_vdf,d_ts,
     s            fluxt,fluxu,fluxv,cdragh,cdragm,
     s            dsens,
     s            ycoefh,yu1,yv1) 
      endif

CXXX Incrementation des flux
      DO i = 1, klon
         sens(i) = - fluxt(i,1) ! flux de chaleur sensible au sol
         fder(i) = dlw(i) + dsens(i) 
      ENDDO
CXXX
      IF (.not. turb_resolved) then !True only for LES
        DO k = 1, klev
        DO i = 1, klon
           t_seri(i,k) = t_seri(i,k) + d_t_vdf(i,k)
           d_t_vdf(i,k)= d_t_vdf(i,k)/dtime          ! K/s
           u_seri(i,k) = u_seri(i,k) + d_u_vdf(i,k)
           d_u_vdf(i,k)= d_u_vdf(i,k)/dtime          ! (m/s)/s
           v_seri(i,k) = v_seri(i,k) + d_v_vdf(i,k)
           d_v_vdf(i,k)= d_v_vdf(i,k)/dtime          ! (m/s)/s
        ENDDO
        ENDDO
      ENDIF
C TRACEURS

      if (iflag_trac.eq.1) then
         DO k = 1, klev
         DO i = 1, klon
            delp(i,k) = paprs(i,k)-paprs(i,k+1)
         ENDDO
         ENDDO
   
         DO iq=1, nqmax
     
             CALL cltrac(dtime,ycoefh,t_seri,
     s               tr_seri(:,:,iq),source(:,iq),
     e               paprs, pplay,delp,
     s               d_tr_vdf(:,:,iq))
     
             tr_seri(:,:,iq) = tr_seri(:,:,iq) + d_tr_vdf(:,:,iq)
             d_tr_vdf(:,:,iq)= d_tr_vdf(:,:,iq)/dtime          ! /s

         ENDDO !nqmax

       endif

      IF (if_ebil.ge.2) THEN 
        ztit='after clmain'
        CALL diagetpq(cell_area,ztit,ip_ebil,2,1,dtime
     e      , t_seri,zero_v2,zero_v2,zero_v2,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
         call diagphy(cell_area,ztit,ip_ebil
     e      , zero_v, zero_v, zero_v, zero_v, sens
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
C
c
c Incrementer la temperature du sol
c
      DO i = 1, klon
         ftsol(i) = ftsol(i) + d_ts(i)
      ENDDO

c Calculer la derive du flux infrarouge
c
      DO i = 1, klon
            dlw(i) = - 4.0*RSIGMA*ftsol(i)**3 
      ENDDO

c-------------------------------
      endif  ! fin du VENUS TEST

      ! tests: output tendencies
!      call writefield_phy('physiq_d_t_vdf',d_t_vdf,klev)
!      call writefield_phy('physiq_d_u_vdf',d_u_vdf,klev)
!      call writefield_phy('physiq_d_v_vdf',d_v_vdf,klev)
!      call writefield_phy('physiq_d_ts',d_ts,1)

c
c Appeler l'ajustement sec
c
c===================================================================
c Convection seche 
c===================================================================
c
      d_t_ajs(:,:)=0.
      d_u_ajs(:,:)=0.
      d_v_ajs(:,:)=0.
      d_tr_ajs(:,:,:)=0.
c
      IF(prt_level>9)WRITE(lunout,*)
     .    'AVANT LA CONVECTION SECHE , iflag_ajs='
     s   ,iflag_ajs

      if(iflag_ajs.eq.0) then
c  Rien
c  ====
         IF(prt_level>9)WRITE(lunout,*)'pas de convection'

      else if(iflag_ajs.eq.1) then

c  Ajustement sec
c  ==============
         IF(prt_level>9)WRITE(lunout,*)'ajsec'

! ADAPTATION GCM POUR CP(T)
         CALL ajsec(paprs, pplay, ppk, t_seri, u_seri, v_seri, nqmax,
     .              tr_seri, d_t_ajs, d_u_ajs, d_v_ajs, d_tr_ajs)

! ADAPTATION GCM POUR CP(T)
         do i=1,klon
          flux_ajs(i,1) = 0.0
          do j=2,klev
            flux_ajs(i,j) = flux_ajs(i,j-1)
     .        + cpdet(t_seri(i,j-1))/RG*d_t_ajs(i,j-1)/dtime
     .                                 *(paprs(i,j-1)-paprs(i,j))
          enddo
         enddo
         
         t_seri(:,:) = t_seri(:,:) + d_t_ajs(:,:)
         d_t_ajs(:,:)= d_t_ajs(:,:)/dtime          ! K/s
         u_seri(:,:) = u_seri(:,:) + d_u_ajs(:,:)
         d_u_ajs(:,:)= d_u_ajs(:,:)/dtime          ! (m/s)/s
         v_seri(:,:) = v_seri(:,:) + d_v_ajs(:,:)
         d_v_ajs(:,:)= d_v_ajs(:,:)/dtime          ! (m/s)/s

         if (iflag_trac.eq.1) then
           tr_seri(:,:,:) = tr_seri(:,:,:) + d_tr_ajs(:,:,:)
           d_tr_ajs(:,:,:)= d_tr_ajs(:,:,:)/dtime  ! /s
         endif
      endif

      ! tests: output tendencies
!      call writefield_phy('physiq_d_t_ajs',d_t_ajs,klev)
!      call writefield_phy('physiq_d_u_ajs',d_u_ajs,klev)
!      call writefield_phy('physiq_d_v_ajs',d_v_ajs,klev)
c
      IF (if_ebil.ge.2) THEN 
        ztit='after dry_adjust'
        CALL diagetpq(cell_area,ztit,ip_ebil,2,2,dtime
     e      , t_seri,zero_v2,zero_v2,zero_v2,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(cell_area,ztit,ip_ebil
     e      , zero_v, zero_v, zero_v, zero_v, sens
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 

c------------------------------------
c    Compute mean mass, cp and R :
c------------------------------------

      if(callthermos) then
         call concentrations2(pplay,t_seri,d_t,tr_seri, nqmax)
      endif

c====================================================================
c RAYONNEMENT
c====================================================================
      if (mod(itap,radpas) == 0) then

      dtimerad = dtime*REAL(radpas)  ! pas de temps du rayonnement (s)
c====================================================================

      if (callnlte .or. callthermos) then
         if (ok_chem) then

!           nlte : use computed chemical species
 
            co2vmr_gcm(:,:) = tr_seri(:,:,i_co2)*mmean(:,:)/m_tr(i_co2)
            covmr_gcm(:,:)  = tr_seri(:,:,i_co)*mmean(:,:)/m_tr(i_co)
            ovmr_gcm(:,:)   = tr_seri(:,:,i_o)*mmean(:,:)/m_tr(i_o)
            n2vmr_gcm(:,:)  = tr_seri(:,:,i_n2)*mmean(:,:)/m_tr(i_n2)

         else

!           nlte : use hedin climatology

            call compo_hedin83_mod(pplay,rmu0,   
     $                 co2vmr_gcm,covmr_gcm,ovmr_gcm,n2vmr_gcm,nvmr_gcm)
         end if
      end if

c   NLTE cooling from CO2 emission
c   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF(callnlte) THEN                                 
            if(nltemodel.eq.0.or.nltemodel.eq.1) then
                CALL nltecool(klon, klev, nqmax, pplay*9.869e-6, t_seri,
     $                    tr_seri, d_t_nlte)
            else if(nltemodel.eq.2) then                                
               CALL nlte_tcool(klon,klev,pplay*9.869e-6,              
     $               t_seri,zzlay,co2vmr_gcm, n2vmr_gcm, covmr_gcm, 
     $               ovmr_gcm,d_t_nlte,ierr_nlte,varerr )
                  if(ierr_nlte.gt.0) then
                     write(*,*)
     $                'WARNING: nlte_tcool output with error message',
     $                'ierr_nlte=',ierr_nlte,'varerr=',varerr
                     write(*,*)'I will continue anyway'
                  endif
             endif
             
        ELSE
 
          d_t_nlte(:,:)=0.

        ENDIF        

c      Find number of layers for LTE radiation calculations

      IF(callnlte .or. callnirco2) 
     $        CALL nlthermeq(klon, klev, paprs, pplay) 

c
c       LTE radiative transfert / solar / IR matrix
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (physideal) then
       CALL radlwsw_newtoncool
     e            (dist, rmu0, fract, zzlev,
     e             paprs, pplay,ftsol, t_seri)
      else
       CALL radlwsw
     e            (dist, rmu0, fract, zzlev,
     e             paprs, pplay,ftsol, t_seri)
      endif

c albedo variations: test for Yeon Joo Lee
c  increment to increase it for 20 Vd => +80% 
c       heat(:,:)=heat(:,:)*(1.+0.80*((rjourvrai-356)+gmtime)/20.)
c  or to decrease it for 20 Vd => 1/1.8 
c       heat(:,:)=heat(:,:)/(1.+0.80*((rjourvrai-356)+gmtime)/20.)

c       CO2 near infrared absorption
c      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        d_t_nirco2(:,:)=0.
        if (callnirco2) then
           call nirco2abs (klon, klev, pplay, dist, nqmax, tr_seri,
     .                 rmu0, fract, d_t_nirco2)
        endif


c          Net atmospheric radiative heating rate (K.s-1)
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        IF(callnlte.or.callnirco2) THEN 
           CALL blendrad(klon, klev, pplay,heat,
     &          cool, d_t_nirco2,d_t_nlte, dtsw, dtlw)
        ELSE
           dtsw(:,:)=heat(:,:)
           dtlw(:,:)=-1*cool(:,:)
        ENDIF

         DO k=1,klev
            DO i=1,klon
               d_t_rad(i,k) = dtsw(i,k) + dtlw(i,k)   ! K/s
            ENDDO
         ENDDO


cc---------------------------------------------

c          EUV heating rate (K.s-1)
c          ~~~~~~~~~~~~~~~~~~~~~~~~

        d_t_euv(:,:)=0.

        IF (callthermos) THEN

c           call euvheat(klon, klev,t_seri,paprs,pplay,zzlay,
c     $          rmu0,pdtphys,gmtime,rjourvrai, co2vmr_gcm, n2vmr_gcm, 
c     $          covmr_gcm, ovmr_gcm,d_t_euv )
           call euvheat(klon, klev, nqmax, t_seri,paprs,pplay,zzlay,
     $         rmu0,dtimerad,gmtime,rjourvrai,
     $         tr_seri, d_tr, d_t_euv )
                 
           DO k=1,klev
              DO ig=1,klon
                 d_t_rad(ig,k)=d_t_rad(ig,k)+d_t_euv(ig,k)
                
              ENDDO
           ENDDO

        ENDIF  ! callthermos

c====================================================================
      ENDIF    ! radpas
c====================================================================
c
c Ajouter la tendance des rayonnements (tous les pas)
c
      DO k = 1, klev
      DO i = 1, klon
         t_seri(i,k) = t_seri(i,k) + d_t_rad(i,k) * dtime
      ENDDO
      ENDDO

! increment physics counter

      itap   = itap + 1

! CONDUCTION  and  MOLECULAR VISCOSITY
c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        d_t_conduc(:,:)=0.
        d_u_molvis(:,:)=0.
        d_v_molvis(:,:)=0.

        IF (callthermos) THEN

           tsurf(:)=t_seri(:,1)
           call conduction(klon, klev,pdtphys,
     $            pplay,paprs,t_seri,
     $            tsurf,zzlev,zzlay,d_t_conduc)

            call molvis(klon, klev,pdtphys,
     $            pplay,paprs,t_seri,
     $            u,tsurf,zzlev,zzlay,d_u_molvis)

            call molvis(klon, klev, pdtphys,
     $            pplay,paprs,t_seri,
     $            v,tsurf,zzlev,zzlay,d_v_molvis)

            DO k=1,klev
               DO ig=1,klon
                  t_seri(ig,k)= t_seri(ig,k)+ d_t_conduc(ig,k)*dtime ! [K]
                  u_seri(ig,k)= u_seri(ig,k)+ d_u_molvis(ig,k)*dtime ! m/s
                  v_seri(ig,k)= v_seri(ig,k)+ d_v_molvis(ig,k)*dtime ! m/s
               ENDDO
            ENDDO
        ENDIF


!  --  MOLECULAR DIFFUSION ---

          d_q_moldif(:,:,:)=0

         IF (callthermos .and. ok_chem) THEN

             call moldiff_red(klon, klev, nqmax,
     &                   pplay,paprs,t_seri, tr_seri, pdtphys,
     &                   zzlay,d_t_euv,d_t_conduc,d_q_moldif)


! --- update tendencies tracers ---

          DO iq = 1, nqmax
           DO k=1,klev
              DO ig=1,klon
                tr_seri(ig,k,iq)= tr_seri(ig,k,iq)+ 
     &                           d_q_moldif(ig,k,iq)*dtime ! [Kg/kg]?
              ENDDO
            ENDDO
           ENDDO
           

         ENDIF  ! callthermos & ok_chem

c====================================================================
      ! tests: output tendencies
!      call writefield_phy('physiq_dtrad',dtrad,klev)
 
      IF (if_ebil.ge.2) THEN 
        ztit='after rad'
        CALL diagetpq(cell_area,ztit,ip_ebil,2,2,dtime
     e      , t_seri,zero_v2,zero_v2,zero_v2,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(cell_area,ztit,ip_ebil
     e      , topsw, toplw, solsw, sollw, zero_v
     e      , zero_v, zero_v, zero_v, ztsol
     e      , d_h_vcol, d_qt, d_ec
     s      , fs_bound, fq_bound )
      END IF 
c

c====================================================================
c   Calcul  des gravity waves  FLOTT
c====================================================================
c
c     if (ok_orodr.or.ok_gw_nonoro) then

c  CALCUL DE N2    
c   UTILISE LA RELATION ENTRE N2 ET STABILITE
c   N2 = RG/T (dT/dz+RG/cp(T))
c   ET DONC EN N'UTILISE QUE LA TEMPERATURE, PAS teta.

       do i=1,klon
        do k=2,klev
          ztlev(i,k)  = (t_seri(i,k)+t_seri(i,k-1))/2.
        enddo
       enddo
       do i=1,klon
        do k=2,klev
          ztlev(i,k)  = (t_seri(i,k)+t_seri(i,k-1))/2.
          zdtlev(i,k) =  t_seri(i,k)-t_seri(i,k-1)
          zdzlev(i,k) = (zphi(i,k)-zphi(i,k-1))/RG
          zn2(i,k) = RG/ztlev(i,k) * ( zdtlev(i,k)/zdzlev(i,k)
     .                                  + RG/cpdet(ztlev(i,k)) )
          zn2(i,k) = max(zn2(i,k),1.e-12)  ! securite
        enddo
        zn2(i,1) = 1.e-12  ! securite
       enddo

c     endif
      
c ----------------------------ORODRAG
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
c A ADAPTER POUR VENUS!!!  [ TN: c'est fait ! ]
        CALL drag_noro(klon,klev,dtime,paprs,pplay,pphi,zn2,
     e                   zmea,zstd, zsig, zgam, zthe,zpic,zval,
     e                   igwd,idx,itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrdr, zvstrdr,
     s                   d_t_oro, d_u_oro, d_v_oro,
     s                   zublstrdr,zvblstrdr,znlow,zeff,zbl,
     s                   ztau,tau0,knu2,kbreak)

c       print*,"d_u_oro=",d_u_oro(klon/2,:)
c  ajout des tendances
           t_seri(:,:) = t_seri(:,:) + d_t_oro(:,:)
           d_t_oro(:,:)= d_t_oro(:,:)/dtime          ! K/s
           u_seri(:,:) = u_seri(:,:) + d_u_oro(:,:)
           d_u_oro(:,:)= d_u_oro(:,:)/dtime          ! (m/s)/s
           v_seri(:,:) = v_seri(:,:) + d_v_oro(:,:)
           d_v_oro(:,:)= d_v_oro(:,:)/dtime          ! (m/s)/s
c    
      ELSE
         d_t_oro = 0.
         d_u_oro = 0.
         d_v_oro = 0.
         zustrdr = 0.
         zvstrdr = 0.
         zublstrdr = 0.
         zvblstrdr = 0.
         znlow = 0.
         zeff = 0.
         zbl = 0
         knu2 = 0
         kbreak = 0
         ztau = 0
         tau0 = 0.
c
      ENDIF ! fin de test sur ok_orodr
c
      ! tests: output tendencies
!      call writefield_phy('physiq_d_t_oro',d_t_oro,klev)
!      call writefield_phy('physiq_d_u_oro',d_u_oro,klev)
!      call writefield_phy('physiq_d_v_oro',d_v_oro,klev)

c ----------------------------OROLIFT
      IF (ok_orolf) THEN
       print*,"ok_orolf NOT IMPLEMENTED !"
       stop
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
c A ADAPTER POUR VENUS!!!
c            CALL lift_noro(klon,klev,dtime,paprs,pplay,
c     e                   latitude_deg,zmea,zstd,zpic,zgam,zthe,zpic,zval,
c     e                   igwd,idx,itest,
c     e                   t_seri, u_seri, v_seri,
c     s                   zulow, zvlow, zustrli, zvstrli,
c     s                   d_t_lif, d_u_lif, d_v_lif               )

c
c  ajout des tendances
           t_seri(:,:) = t_seri(:,:) + d_t_lif(:,:)
           d_t_lif(:,:)= d_t_lif(:,:)/dtime          ! K/s
           u_seri(:,:) = u_seri(:,:) + d_u_lif(:,:)
           d_u_lif(:,:)= d_u_lif(:,:)/dtime          ! (m/s)/s
           v_seri(:,:) = v_seri(:,:) + d_v_lif(:,:)
           d_v_lif(:,:)= d_v_lif(:,:)/dtime          ! (m/s)/s
c
      ELSE
         d_t_lif = 0.
         d_u_lif = 0.
         d_v_lif = 0.
         zustrli = 0.
         zvstrli = 0.
c
      ENDIF ! fin de test sur ok_orolf

c ---------------------------- NON-ORO GRAVITY WAVES
       IF(ok_gw_nonoro) then

      call flott_gwd_ran(klon,klev,dtime,pplay,zn2,
     e               t_seri, u_seri, v_seri, paprs(klon/2+1,:),
     o               zustrhi,zvstrhi,
     o               d_t_hin, d_u_hin, d_v_hin)

c  ajout des tendances

         t_seri(:,:) = t_seri(:,:) + d_t_hin(:,:)
         d_t_hin(:,:)= d_t_hin(:,:)/dtime          ! K/s
         u_seri(:,:) = u_seri(:,:) + d_u_hin(:,:)
         d_u_hin(:,:)= d_u_hin(:,:)/dtime          ! (m/s)/s
         v_seri(:,:) = v_seri(:,:) + d_v_hin(:,:)
         d_v_hin(:,:)= d_v_hin(:,:)/dtime          ! (m/s)/s

      ELSE
         d_t_hin = 0.
         d_u_hin = 0.
         d_v_hin = 0.
         zustrhi = 0.
         zvstrhi = 0.

      ENDIF ! fin de test sur ok_gw_nonoro

      ! tests: output tendencies
!      call writefield_phy('physiq_d_t_hin',d_t_hin,klev)
!      call writefield_phy('physiq_d_u_hin',d_u_hin,klev)
!      call writefield_phy('physiq_d_v_hin',d_v_hin,klev)

c====================================================================
c Transport de ballons 
c====================================================================
      if (ballons.eq.1) then
        CALL ballon(30,pdtphys,rjourvrai,gmtime*RDAY,
     &              latitude_deg,longitude_deg,
c    C              t,pplay,u,v,pphi)   ! alt above surface (smoothed for GCM)
     C              t,pplay,u,v,zphi)   ! alt above planet average radius
      endif !ballons

c====================================================================
c Bilan de mmt angulaire
c====================================================================
      if (bilansmc.eq.1) then
CMODDEB FLOTT
C  CALCULER LE BILAN DE MOMENT ANGULAIRE (DIAGNOSTIQUE)
C STRESS NECESSAIRES: COUCHE LIMITE ET TOUTE LA PHYSIQUE

      DO i = 1, klon
        zustrph(i)=0.
        zvstrph(i)=0.
        zustrcl(i)=0.
        zvstrcl(i)=0.
      ENDDO
      DO k = 1, klev
      DO i = 1, klon
       zustrph(i)=zustrph(i)+(u_seri(i,k)-u(i,k))/dtime*
     c            (paprs(i,k)-paprs(i,k+1))/rg
       zvstrph(i)=zvstrph(i)+(v_seri(i,k)-v(i,k))/dtime*
     c            (paprs(i,k)-paprs(i,k+1))/rg
       zustrcl(i)=zustrcl(i)+d_u_vdf(i,k)*
     c            (paprs(i,k)-paprs(i,k+1))/rg
       zvstrcl(i)=zvstrcl(i)+d_v_vdf(i,k)*
     c            (paprs(i,k)-paprs(i,k+1))/rg
      ENDDO
      ENDDO

      CALL aaam_bud (27,klon,klev,rjourvrai,gmtime*RDAY,
     C               ra,rg,romega,
     C               latitude_deg,longitude_deg,pphis,
     C               zustrdr,zustrli,zustrcl,
     C               zvstrdr,zvstrli,zvstrcl,
     C               paprs,u,v)
                     
CCMODFIN FLOTT
      endif !bilansmc

c====================================================================
c====================================================================
c Calculer le transport de l'eau et de l'energie (diagnostique)
c
c  A REVOIR POUR VENUS...
c
c     CALL transp (paprs,ftsol,
c    e                   t_seri, q_seri, u_seri, v_seri, zphi,
c    s                   ve, vq, ue, uq)
c
c====================================================================
c+jld ec_conser
      DO k = 1, klev
      DO i = 1, klon
        d_t_ec(i,k)=0.5/cpdet(t_seri(i,k))
     $      *(u(i,k)**2+v(i,k)**2-u_seri(i,k)**2-v_seri(i,k)**2)
        t_seri(i,k)=t_seri(i,k)+d_t_ec(i,k)
        d_t_ec(i,k) = d_t_ec(i,k)/dtime
       END DO 
      END DO 
         do i=1,klon
          flux_ec(i,1) = 0.0
          do j=2,klev
            flux_ec(i,j) = flux_ec(i,j-1)
     . +cpdet(t_seri(i,j-1))/RG*d_t_ec(i,j-1)*(paprs(i,j-1)-paprs(i,j))
          enddo
         enddo
         
c-jld ec_conser
c====================================================================
      IF (if_ebil.ge.1) THEN 
        ztit='after physic'
        CALL diagetpq(cell_area,ztit,ip_ebil,1,1,dtime
     e      , t_seri,zero_v2,zero_v2,zero_v2,u_seri,v_seri,paprs,pplay
     s      , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
C     Comme les tendances de la physique sont ajoute dans la dynamique,
C     on devrait avoir que la variation d'entalpie par la dynamique
C     est egale a la variation de la physique au pas de temps precedent.
C     Donc la somme de ces 2 variations devrait etre nulle.
        call diagphy(cell_area,ztit,ip_ebil
     e      , topsw, toplw, solsw, sollw, sens
     e      , zero_v, zero_v, zero_v, ztsol
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

c Convertir les incrementations en tendances
c
      DO k = 1, klev
      DO i = 1, klon
         d_u(i,k) = ( u_seri(i,k) - u(i,k) ) / dtime
         d_v(i,k) = ( v_seri(i,k) - v(i,k) ) / dtime
         d_t(i,k) = ( t_seri(i,k) - t(i,k) ) / dtime
      ENDDO
      ENDDO
c
      DO iq = 1, nqmax
      DO  k = 1, klev
      DO  i = 1, klon
         d_qx(i,k,iq) = ( tr_seri(i,k,iq) - qx(i,k,iq) ) / dtime
      ENDDO
      ENDDO
      ENDDO
      
c------------------------
c Calcul moment cinetique
c------------------------
c TEST VENUS...
c     mangtot = 0.0
c     DO k = 1, klev
c     DO i = 1, klon
c       mang(i,k) = RA*cos(latitude(i))
c    .     *(u_seri(i,k)+RA*cos(latitude(i))*ROMEGA)
c    .     *cell_area(i)*(paprs(i,k)-paprs(i,k+1))/RG
c       mangtot=mangtot+mang(i,k)
c     ENDDO
c     ENDDO
c     print*,"Moment cinetique total = ",mangtot
c
c------------------------
c
c Sauvegarder les valeurs de t et u a la fin de la physique:
c
      DO k = 1, klev
      DO i = 1, klon
         u_ancien(i,k) = u_seri(i,k)
         t_ancien(i,k) = t_seri(i,k)
      ENDDO
      ENDDO
c
c=============================================================
c   Ecriture des sorties
c=============================================================
#ifndef MESOSCALE       
#ifdef CPP_IOIPSL

#ifdef histhf
#include "write_histhf.h"
#endif

#ifdef histday
#include "write_histday.h"
#endif

#ifdef histmth
#include "write_histmth.h"
#endif

#ifdef histins
#include "write_histins.h"
#endif

#endif

! XIOS outputs
! This can be done ANYWHERE in the physics routines !

#ifdef CPP_XIOS      
! Send fields to XIOS: (NB these fields must also be defined as
! <field id="..." /> in context_lmdz_physics.xml to be correctly used)
      
! 2D fields

      CALL send_xios_field("phis",pphis)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      CALL send_xios_field("aire",cell_area_out)
      CALL send_xios_field("tsol",ftsol)
      CALL send_xios_field("psol",paprs(:,1))
      CALL send_xios_field("cdragh",cdragh)
      CALL send_xios_field("cdragm",cdragm)

      CALL send_xios_field("tops",topsw)
      CALL send_xios_field("topl",toplw)
      CALL send_xios_field("sols",solsw)
      CALL send_xios_field("soll",sollw)

! 3D fields

      CALL send_xios_field("temp",t_seri)
      CALL send_xios_field("pres",pplay)
      CALL send_xios_field("geop",zphi)
      CALL send_xios_field("vitu",u_seri)
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
      CALL send_xios_field("vitv",-1.*v_seri)
      CALL send_xios_field("vitw",omega)
      CALL send_xios_field("Kz",ycoefh)
      CALL send_xios_field("mmean",mmean)
      CALL send_xios_field("rho",rho)
      CALL send_xios_field("BV2",zn2)

      CALL send_xios_field("dudyn",d_u_dyn)
      CALL send_xios_field("duvdf",d_u_vdf)
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
      CALL send_xios_field("dvvdf",-1.*d_v_vdf)
      CALL send_xios_field("duajs",d_u_ajs)
      CALL send_xios_field("dugwo",d_u_oro)
      CALL send_xios_field("dugwno",d_u_hin)
      CALL send_xios_field("dumolvis",d_u_molvis)
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
      CALL send_xios_field("dvmolvis",-1.*d_v_molvis)
      CALL send_xios_field("dtdyn",d_t_dyn)
      CALL send_xios_field("dtphy",d_t)
      CALL send_xios_field("dtvdf",d_t_vdf)
      CALL send_xios_field("dtajs",d_t_ajs)
      CALL send_xios_field("dtswr",dtsw)
      CALL send_xios_field("dtswrNLTE",d_t_nirco2)
      CALL send_xios_field("dtswrLTE",heat)
      CALL send_xios_field("dtlwr",dtlw)
      CALL send_xios_field("dtlwrNLTE",d_t_nlte)
      CALL send_xios_field("dtlwrLTE",-1.*cool)
      CALL send_xios_field("dteuv",d_t_euv)
      CALL send_xios_field("dtcond",d_t_conduc)
      CALL send_xios_field("dtec",d_t_ec)

      CALL send_xios_field("SWnet",swnet(:,1:klev))
      CALL send_xios_field("LWnet",lwnet(:,1:klev))
      CALL send_xios_field("fluxvdf",fluxt)
      CALL send_xios_field("fluxdyn",flux_dyn)
      CALL send_xios_field("fluxajs",flux_ajs)
      CALL send_xios_field("fluxec",flux_ec)

! when using tracers

      if (iflag_trac == 1) then

! tracers in gas phase, volume mixing ratio

         do iq = 1,nqmax - nmicro
            call send_xios_field(tname(iq),
     $                           qx(:,:,iq)*mmean(:,:)/m_tr(iq))
         end do

! tracers in liquid phase, volume mixing ratio

         if ((tr_scheme == 3) .and. (cl_scheme == 1)) THEN
            call send_xios_field(tname(i_h2oliq),
     $             qx(:,:,i_h2oliq)*mmean(:,:)/m_tr(i_h2oliq))
            call send_xios_field(tname(i_h2so4liq),
     $             qx(:,:,i_h2so4liq)*mmean(:,:)/m_tr(i_h2so4liq))
            if (ok_sedim) then
               call send_xios_field("Fsedim",fsedim(:,1:klev))
            end if
         end if

! chemical iterations

         call send_xios_field("iter",real(iter))

      end if

      IF (callthermos .and. ok_chem) THEN
       CALL send_xios_field("d_qmoldifCO2",d_q_moldif(:,:,i_co2))
       CALL send_xios_field("d_qmoldifO3p",d_q_moldif(:,:,i_o))
       CALL send_xios_field("d_qmoldifN2",d_q_moldif(:,:,i_n2))
      ENDIF

      if (lafin.and.is_omp_master) then
        write(*,*) "physiq: call xios_context_finalize"
        call xios_context_finalize
      endif

#endif
#else
! Outputs MESOSCALE
      CALL allocate_comm_wrf(klon,klev)
      comm_HR_SW(1:klon,1:klev) = dtsw(1:klon,1:klev)
      comm_HR_LW(1:klon,1:klev) = dtlw(1:klon,1:klev)
      comm_DT_RAD(1:klon,1:klev) = d_t_rad(1:klon,1:klev)
      IF (turb_resolved) THEN
        open(17,file='hrdyn.txt',form='formatted',status='old')
        rewind(17)
        DO k=1,klev
          read(17,*) dt_dyn(k)
        ENDDO
        close(17)

        do i=1,klon
          d_t(i,:)=d_t(i,:)+dt_dyn(:)
          comm_HR_DYN(i,:) = dt_dyn(:)
        enddo
       ELSE
         comm_HR_DYN(1:klon,1:klev) = d_t_dyn(1:klon,1:klev)
         comm_DT_VDF(1:klon,1:klev) = d_t_vdf(1:klon,1:klev)
         comm_DT_AJS(1:klon,1:klev) = d_t_ajs(1:klon,1:klev)
       ENDIF
      comm_DT(1:klon,1:klev)=d_t(1:klon,1:klev)
#endif


c====================================================================
c Si c'est la fin, il faut conserver l'etat de redemarrage
c====================================================================
c
      IF (lafin) THEN
         itau_phy = itau_phy + itap
         CALL phyredem ("restartphy.nc")
     
c--------------FLOTT
CMODEB LOTT
C  FERMETURE DU FICHIER FORMATTE CONTENANT LES COMPOSANTES
C  DU BILAN DE MOMENT ANGULAIRE.
      if (bilansmc.eq.1) then
        write(*,*)'Fermeture de aaam_bud.out (FL Vous parle)'
        close(27)                                     
        close(28)                                     
      endif !bilansmc
CMODFIN
c-------------
c--------------SLEBONNOIS
C  FERMETURE DES FICHIERS FORMATTES CONTENANT LES POSITIONS ET VITESSES
C  DES BALLONS
      if (ballons.eq.1) then
        write(*,*)'Fermeture des ballons*.out'
        close(30)                                     
        close(31)                                     
        close(32)                                     
        close(33)                                     
        close(34)                                     
      endif !ballons
c-------------
      ENDIF
      
      END SUBROUTINE physiq

      END MODULE physiq_mod

