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
c Modifications pour la physique de Titan
c     S. Lebonnois (LMD/CNRS) Juin 2013: Parallelisation
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
c gmtime--input-R-temps universel dans la journee (0 a RDAY s)
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
      USE infotrac_phy, ONLY: iflag_trac, tname, ttext
      use dimphy
      USE geometry_mod, ONLY: longitude, latitude, ! in radians
     &                        longitude_deg, latitude_deg, ! in degrees
     &                        cell_area, dx, dy
      use cpdet_phy_mod, only: cpdet, t2tpot
      USE mod_phys_lmdz_para, only : is_parallel,jj_nb,
     &                               is_north_pole_phy,
     &                               is_south_pole_phy
      USE phys_state_var_mod ! Variables sauvegardees de la physique
      USE iophy
      USE common_mod, only: rmcbar,xfbar,ncount,
     &      flxesp_i,tau_drop,tau_aer,solesp,precip,
     &      evapch4,occcld_m,occcld,satch4,satc2h6,satc2h2,rmcloud,
     &      TauHID,TauHVD,TauGID,TauGVD,TauCID,TauCVD,NSPECV,NSPECI,
     &      common_init

      USE moyzon_mod
      USE write_field_phy
      USE time_phylmdz_mod, only: itau_phy,day_ref,annee_ref,nday
      USE logic_mod, only: moyzon_ch,moyzon_mu
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev
      IMPLICIT none
c======================================================================
c   CLEFS CPP POUR LES IO
c   =====================
#define histday
#define histmth
#define histins
c======================================================================
#include "dimensions.h"
      integer jjmp1
      parameter (jjmp1=jjm+1-1/jjm)
#include "dimsoil.h"
#include "clesphys.h"
#include "iniprint.h"
#include "tabcontrol.h"
#include "comorbit.h"
#include "microtab.h"
#include "itemps.h"
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

c Variables propres a la physique
c
      REAL,save,allocatable :: rlev(:,:) ! altitude a chaque niveau (interface inferieure de la couche) 
      INTEGER,save :: itap        ! compteur pour la physique
      REAL delp(klon,klev)        ! epaisseur d'une couche
      REAL omega(klon,klev)
      
      INTEGER igwd,idx(klon),itest(klon)
c
c  Diagnostiques 2D de drag_noro, lift_noro et gw_nonoro

      REAL zulow(klon),zvlow(klon)
      REAL zustrdr(klon), zvstrdr(klon)
      REAL zustrli(klon), zvstrli(klon)
      REAL zustrhi(klon), zvstrhi(klon)

c Pour calcul GW drag oro et nonoro: CALCUL de N2:
      real zdzlev(klon,klev)
      real ztlev(klon,klev),zpklev(klon,klev)
      real ztetalay(klon,klev),ztetalev(klon,klev)
      real zdtetalev(klon,klev)
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
      integer nmicro
      save    nmicro
      character*8 nom
      REAL qaer(klon,klev,nqmax)

      REAL ycoefh(klon,klev)    ! coef d'echange pour phytrac
      REAL yu1(klon)            ! vents dans la premiere couche U
      REAL yv1(klon)            ! vents dans la premiere couche V

      REAL sens(klon), dsens(klon) ! chaleur sensible et sa derivee
      REAL ve(klon) ! integr. verticale du transport meri. de l'energie
      REAL vq(klon) ! integr. verticale du transport meri. de l'eau
      REAL ue(klon) ! integr. verticale du transport zonal de l'energie
      REAL uq(klon) ! integr. verticale du transport zonal de l'eau
c

c======================================================================
c
c Declaration des procedures appelees
c
      EXTERNAL ajsec     ! ajustement sec
      EXTERNAL clmain    ! couche limite 
      EXTERNAL hgardfou  ! verifier les temperatures
      EXTERNAL orbite    ! calculer l'orbite 
      EXTERNAL phyetat0  ! lire l'etat initial de la physique
      EXTERNAL phyredem  ! ecrire l'etat de redemarrage de la physique
      EXTERNAL radlwsw   ! rayonnements solaire et infrarouge
      EXTERNAL suphec    ! initialiser certaines constantes
c     EXTERNAL transp    ! transport total de l'eau et de l'energie
      EXTERNAL abort_gcm
      EXTERNAL printflag
      EXTERNAL zenang
      EXTERNAL diagetpq
      EXTERNAL conf_phys
      EXTERNAL diagphy
      EXTERNAL mucorr
      EXTERNAL phytrac
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
      REAL    dtimerad
      INTEGER itaprad
      SAVE itaprad,dtimerad
      REAL zdtime
c
c CHIMIE

      REAL    dtimechim
      INTEGER itapchim,appel_chim
      SAVE itapchim,dtimechim

c ORBITE

      REAL dist, rmu0(klon), fract(klon), pdecli
      REAL rmu0bar(klon), fractbar(klon)
      REAL zday
      REAL zls,zlsdeg,zlsm1
      save zlsm1
c
      INTEGER i, k, iq, ig, j, ll, l
c
      REAL zphi(klon,klev)
      REAL zzlev(klon,klev+1),zzlay(klon,klev),z1,z2
c
c Variables du changement
c
c ajs: ajustement sec
c vdf: couche limite (Vertical DiFfusion)
c mph: microphysique
c kim: chimie
      REAL d_t_ajs(klon,klev), d_tr_ajs(klon,klev,nqmax)
      REAL d_u_ajs(klon,klev), d_v_ajs(klon,klev)
c
      REAL d_ts(klon)
c
      REAL d_u_vdf(klon,klev), d_v_vdf(klon,klev)
      REAL d_t_vdf(klon,klev), d_tr_vdf(klon,klev,nqmax)
c
      REAL d_tr_mph(klon,klev,nqmax),d_tr_kim(klon,klev,nqmax)

CMOD LOTT: Tendances Orography Sous-maille
      REAL d_u_oro(klon,klev), d_v_oro(klon,klev)
      REAL d_t_oro(klon,klev)
      REAL d_u_lif(klon,klev), d_v_lif(klon,klev)
      REAL d_t_lif(klon,klev)
C          Tendances Ondes de G non oro (runs strato).
      REAL d_u_hin(klon,klev), d_v_hin(klon,klev)
      REAL d_t_hin(klon,klev)

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
      REAL tr_seri(klon,klev,nqmax)
      REAL d_tr(klon,klev,nqmax)
c
c pour ioipsl
      INTEGER nid_day, nid_mth, nid_ins
      SAVE nid_day, nid_mth, nid_ins
      INTEGER nhori, nvert, idayref
      REAL zsto, zout, zsto1, zsto2, zero
      parameter (zero=0.0e0)
      real zjulian
      save zjulian
      REAL tmpout(klon,klev)  ! pour sorties

      CHARACTER*1  str1
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

c Temporaire avant de trouver mieux :
c Recuperation des TAU du TR
      REAL t_tauhvd(klon,klev),t_khvd(klon,klev)
      REAL t_tcld(klon,klev),t_kcld(klon,klev)
      REAL t_kcvd(klon,klev)

       REAL ch4(klon,jjm+1),dch4(jjm+1)
       INTEGER ig0
       integer ich4
       common/ch4ind/ich4

c     flux de chaleur latente d'evaporation CH4
      REAL fclat(klon)
c     reservoir de surface
      REAL,save,allocatable :: reservoir(:)

c cell_area for outputs in hist*
      REAL cell_area_out(klon)

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
      if (is_parallel) then
        bilansmc = 0
        ballons  = 0
      endif

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
         allocate(rlev(klon,klevp1))
         allocate(source(klon,nqmax))
         allocate(reservoir(klon))

         CALL suphec ! initialiser constantes et parametres phys.

         IF (if_ebil.ge.1) d_h_vcol_phy=0.
c
c appel a la lecture du physiq.def 
c
         call conf_phys(ok_mensuel,ok_journe,
     .                  ok_instan,
     .                  if_ebil)

         call phys_state_var_init
         call common_init
c
c Initialiser les compteurs:
c
         itap    = 0
         itaprad = 0
         itapchim    = 1

c init rnuabar
         ncount(:,:) = 0
         rmcbar  = 0.
         xfbar   = 0.
         
c         
c Lecture startphy.nc :
c
         CALL phyetat0 ("startphy.nc")

c dtime est defini dans tabcontrol.h et lu dans startphy
c pdtphys est calcule a partir des nouvelles conditions:
c Reinitialisation du pas de temps physique quand changement
         IF (ABS(dtime-pdtphys).GT.0.001) THEN
            WRITE(lunout,*) 'Pas physique a change',dtime,
     .                        pdtphys
c           abort_message='Pas physique n est pas correct '
c           call abort_gcm(modname,abort_message,1)
c----------------
c pour initialiser convenablement le time_counter, il faut tenir compte
c du changement de dtime en changeant itau_phy (point de depart)
            itau_phy = NINT(itau_phy*dtime/pdtphys)
c----------------
            dtime=pdtphys
         ENDIF

         radpas  = NINT( RDAY/pdtphys/nbapp_rad)
         chimpas =   radpas*nbapp_rad/nbapp_chim

         CALL printflag( ok_mensuel,ok_journe,ok_instan )
c
c Initialiser les pas de temps:
c
      dtimerad = dtime*REAL(radpas)  ! pas de temps du rayonnement (s)
c      PRINT*,'dtimerad,dtime,radpas',dtimerad,dtime,radpas
            
      dtimechim = dtime*REAL(chimpas)  ! pas de temps de la chimie (s)
c      PRINT*,'dtimechim,dtime,chimpas',dtimechim,dtime,chimpas


c INITIALISATION ORBITE

         CALL iniorbit(aphelie,periheli,year_day,peri_day,obliquit)

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
         source = 0.0 ! pas de source, pour l'instant
C
c Si microphysique offline, pas besoin d'avoir de traceurs microphysiques
c car on lit les profils verticaux des qaer dans une look-up table pour 
c le rayonnement. 

c  calcul de nmicro
c !!!! Les traceurs microphysiques doivent etre toujours en premiers!!

      nmicro = 0
      do iq=1,nqmax
         nom = tname(iq)
c        print*,iq,"nom=",nom,"tname=",tname(iq)
         print*,iq,"nom=",nom
         if (nom(1:1).eq."q") then
	   nmicro = nmicro+1
	 endif
      enddo
      print*,"nmicro=",nmicro

c --------------
c Verifications:
c --------------
         IF ((nmicro.eq.0).and.(microfi.eq.1)) THEN
           abort_message="MICROPHYSIQUE ONLINE, MAIS NMICRO=0..."
           call abort_gcm(modname,abort_message,1) 
         ENDIF
         IF (microfi.lt.1.and.clouds.eq.1) THEN
          write(lunout,*)"microfi.lt.1.and.clouds.eq.1"
          abort_message = 
     &    "Impossible de faire des nuages sans microphysique..."
          call abort_gcm(modname,abort_message,1)
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

         IF (((moyzon_mu).and.(microfi.ne.1)).or.
     .       ((.not.moyzon_mu).and.(microfi.eq.1))) THEN
           abort_message="Microphysic 2D and moyzon_mu not compatible"
           write(lunout,*) "moyzon_mu=",moyzon_mu
           write(lunout,*) "microfi=",microfi
           call abort_gcm(modname,abort_message,1) 
         ENDIF
         IF (((moyzon_ch).and.(.not.chimi)).or.
     .       ((.not.moyzon_ch).and.(chimi))) THEN
           abort_message="Chemistry and moyzon_ch not compatible"
           write(lunout,*) "moyzon_ch=",moyzon_ch
           write(lunout,*) "chimi=",chimi
           call abort_gcm(modname,abort_message,1) 
         ENDIF

         IF (dtime*REAL(radpas).GT.(RDAY*0.25).AND.cycle_diurne)
     $    THEN 
           WRITE(lunout,*)'Nbre d appels au rayonnement insuffisant'
           WRITE(lunout,*)"Au minimum 4 appels par jour si cycle diurne"
           abort_message='Nbre d appels au rayonnement insuffisant'
           call abort_gcm(modname,abort_message,1)
         ENDIF
c
         WRITE(lunout,*)"Clef pour la convection seche, iflag_ajs=",
     .                   iflag_ajs
c
         ecrit_mth = NINT(RDAY/dtime) *nday  ! tous les nday jours
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

      ENDIF ! debut
c====================================================================
c======================================================================

c   Creer un reservoir de surface infini 
c
      reservoir(:) = 2.

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
c
c Calcule de vitesse verticale a partir de flux de masse verticale
      DO k = 1, klev
       DO i = 1, klon
        omega(i,k) = RG*flxmw(i,k) / cell_area(i)
       END DO
      END DO

c Ajouter le geopotentiel du sol:
c
      DO k = 1, klev
      DO i = 1, klon
         zphi(i,k) = pphi(i,k) + pphis(i)
      ENDDO
      ENDDO

c     call WriteField_phy('physiq_pphi',pphi,klev)
c     call WriteField_phy('physiq_pphis',pphis,1)

c   calcul du geopotentiel aux niveaux intercouches
c   ponderation des altitudes au niveau des couches en dp/p

      DO l=1,klev
         DO i=1,klon
c           zzlay(i,l)=zphi(i,l)/RG
c SI ON TIENT COMPTE DE LA VARIATION DE G AVEC L'ALTITUDE:
            zzlay(i,l)=RG*RA*RA/(RG*RA-zphi(i,l))-RA
         ENDDO
      ENDDO
      DO i=1,klon
c        zzlev(i,1)=0.
c CORRECTION 13/01/2011  
c (correspond a la position de la surface en ce point vs RA)
         zzlev(i,1)=pphis(i)/RG
      ENDDO
      DO l=2,klev
         DO i=1,klon
            z1=(pplay(i,l-1)+paprs(i,l))/(pplay(i,l-1)-paprs(i,l))
            z2=(paprs(i,l)  +pplay(i,l))/(paprs(i,l)  -pplay(i,l))
            zzlev(i,l)=(z1*zzlay(i,l-1)+z2*zzlay(i,l))/(z1+z2)
         ENDDO
      ENDDO
      DO i=1,klon
         zzlev(i,klev+1)=zzlay(i,klev)+(zzlay(i,klev)-zzlev(i,klev))
      ENDDO

! zonal averages needed
      if (moyzon_ch.or.moyzon_mu) then

c      zzlaybar(1,:)=(zphibar(1,:)+zphisbar(1))/RG
c SI ON TIENT COMPTE DE LA VARIATION DE G AVEC L'ALTITUDE:
       zzlaybar(1,:)=RG*RA*RA/(RG*RA-(zphibar(1,:)+zphisbar(1)))-RA
       zzlevbar(1,1)=zphisbar(1)/RG
       DO l=2,klev
            z1=(zplaybar(1,l-1)+zplevbar(1,l))/
     .            (zplevbar(1,l-1)-zplevbar(1,l))
            z2=(zplevbar(1,l)  +zplaybar(1,l))/
     .            (zplevbar(1,l)  -zplaybar(1,l))
            zzlevbar(1,l)=(z1*zzlaybar(1,l-1)+z2*zzlaybar(1,l))/(z1+z2)
       ENDDO
       zzlevbar(1,klev+1)=zzlaybar(1,klev)+
     .            (zzlaybar(1,klev)-zzlevbar(1,klev))

       DO i=2,klon
        if (latitude(i).ne.latitude(i-1)) then
         DO l=1,klev
c         zzlaybar(i,l)=(zphibar(i,l)+zphisbar(i))/RG
c SI ON TIENT COMPTE DE LA VARIATION DE G AVEC L'ALTITUDE:
          zzlaybar(i,l)=RG*RA*RA/(RG*RA-(zphibar(i,l)+zphisbar(i)))-RA
         ENDDO
         zzlevbar(i,1)=zphisbar(i)/RG
         DO l=2,klev
            z1=(zplaybar(i,l-1)+zplevbar(i,l))/
     .            (zplevbar(i,l-1)-zplevbar(i,l))
            z2=(zplevbar(i,l)  +zplaybar(i,l))/
     .            (zplevbar(i,l)  -zplaybar(i,l))
            zzlevbar(i,l)=(z1*zzlaybar(i,l-1)+z2*zzlaybar(i,l))/(z1+z2)
         ENDDO
         zzlevbar(i,klev+1)=zzlaybar(i,klev)+
     .              (zzlaybar(i,klev)-zzlevbar(i,klev))
        else
         zzlaybar(i,:)=zzlaybar(i-1,:)
         zzlevbar(i,:)=zzlevbar(i-1,:)
        endif
       ENDDO

      endif  ! moyzon

c     call WriteField_phy('physiq_zphi',zphi,klev)
c     call WriteField_phy('physiq_zzlay',zzlay,klev)
c     call WriteField_phy('physiq_zzlev',zzlev,klev+1)
c- - - - - - - - - - - - - - - -
c DIAGNOSTIQUE GRILLE VERTICALE
c- - - - - - - - - - - - - - - -
c     print*,"DIAGNOSTIQUE GRILLE VERTICALE"
c     i=klon/2
c     print*,"Niveau  Pression  Altitude    (lev puis lay)"
c     do l=1,klev
c      print*,l,paprs(i,l),zzlev(i,l)
c      print*,l,pplay(i,l),zzlay(i,l)
c     enddo
c     print*,klev+1,paprs(i,klev+1),zzlev(i,klev+1)
c     stop

c====================================================================
c
c Verifier les temperatures
c
      CALL hgardfou(t_seri,ftsol,'debutphy')
c====================================================================
c
c Incrementer le compteur de la physique
c
      itap   = itap + 1

c====================================================================
c
c Epaisseurs couches

      DO k = 1, klev
      DO i = 1, klon
         delp(i,k) = paprs(i,k)-paprs(i,k+1)
      ENDDO
      ENDDO

c====================================================================
c Orbite et eclairement
c====================================================================

c Pour TITAN:
c  calcul de la longitude solaire
          CALL solarlong(rjourvrai+gmtime,zls)
          zlsdeg = zls*180./RPI      ! zls est en radians !!
          print*,'Ls',zlsdeg

      CALL orbite(zls,dist,pdecli) 
      IF (debut) zlsm1=zls

c dans zenang, Ls en degres ; dans mucorr, Ls en radians
      call mucorr(klon,zls,latitude_deg,rmu0bar,fractbar)
      IF (cycle_diurne) THEN
        zdtime=dtime*REAL(radpas) ! pas de temps du rayonnement (s)
        CALL zenang(zlsdeg,gmtime,zdtime,latitude_deg,longitude_deg,
     &              rmu0,fract)
      ELSE
        rmu0  = rmu0bar
        fract = fractbar
      ENDIF
      
c====================================================================
c Appeler la diffusion verticale (programme de couche limite)
c====================================================================

c-------------------------------
c TEST: on ne tient pas compte des calculs de clmain mais on force
c l'equilibre radiatif du sol
      if (1.eq.0) then
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

c     print*,"radsol avant clmain=",radsol(klon/2)
c     print*,"solsw avant clmain=",solsw(klon/2)
c     print*,"sollw avant clmain=",sollw(klon/2)

! ADAPTATION GCM POUR CP(T)

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

c     print*,"radsol apres clmain=",radsol(klon/2)
c     print*,"solsw apres clmain=",solsw(klon/2)
c     print*,"sollw apres clmain=",sollw(klon/2)

CXXX Incrementation des flux
      DO i = 1, klon
         sens(i) = - fluxt(i,1) ! flux de chaleur sensible au sol
         fder(i) = dlw(i) + dsens(i) 
      ENDDO
CXXX

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

c     call WriteField_phy('physiq_dtvdf',d_t_vdf,klev)
c     call WriteField_phy('physiq_duvdf',d_u_vdf,klev)
c     call WriteField_phy('physiq_dvvdf',d_v_vdf,klev)

C TRACEURS

      d_tr_vdf = 0.
      if (iflag_trac.eq.1) then
         DO iq=1, nqmax
             CALL cltrac(dtime,ycoefh,t_seri,
     s               tr_seri(1,1,iq),source,
     e               paprs, pplay,delp,
     s               d_tr_vdf(1,1,iq))
             tr_seri(:,:,iq) = tr_seri(:,:,iq) + d_tr_vdf(:,:,iq)
             d_tr_vdf(:,:,iq)= d_tr_vdf(:,:,iq)/dtime          ! /s
         ENDDO
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
c     print*,'Tsol avant clmain:',ftsol(1)
      DO i = 1, klon
         ftsol(i) = ftsol(i) + d_ts(i)
      ENDDO
c     print*,'DTsol apres clmain:',d_ts(klon/2)
c     print*,'Tsol apres clmain:',ftsol(1)

c Calculer la derive du flux infrarouge
c
      DO i = 1, klon
            dlw(i) = - 4.0*emis*RSIGMA*ftsol(i)**3 
      ENDDO

c-------------------------------
      endif  ! fin du TEST

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
     .                                 *delp(i,j-1)
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

c     call WriteField_phy('physiq_dtajs',d_t_ajs,klev)
c     call WriteField_phy('physiq_duajs',d_u_ajs,klev)
c     call WriteField_phy('physiq_dvajs',d_v_ajs,klev)

      endif
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

c====================================================================
c   MICROPHYSIQUE ET CHIMIE
c====================================================================

      d_tr_mph(:,:,:)=0.
      d_tr_kim(:,:,:)=0.
      
c on recupere tr_seri inchange, d_tr_micro, d_tr_chim, tous les trois sur nqmax
c on recupere aussi qaer pour le mettre dans les sorties
c  si microfi=1, sortie de qaer(1:nmicro)
c  si nmicro != nqmax et si chimi, sortie de tr_seri(nmicro+1:nqmax) 

c faire un test comme pour rayonnement, avec chimi en plus comme flag, 
c pour voir si chimie appelee -> bouleen, qui passe dans phytrac.
c faut aussi le pas de temps chimique: dtimechim, a passer..

      appel_chim = 0
      IF (MOD(itapchim,chimpas).EQ.0) THEN
c             print*,'CHIMIE ', 
c    $             ' (itapchim=',itapchim,'/chimpas=',chimpas,')'
       appel_chim = 1
       itapchim = 0
      ENDIF
      itapchim = itapchim + 1

      if (iflag_trac.eq.1) then
c        call WriteField_phy('physiq_qaer01',
c    .                          qaer(:,:,1),klev)
c        call WriteField_phy('physiq_qaer10',
c    .                          qaer(:,:,10),klev)
c        call WriteField_phy('physiq_tr_seri01',
c    .                          tr_seri(:,:,1),klev)
c        call WriteField_phy('physiq_tr_seri10',
c    .                          tr_seri(:,:,10),klev)

c         call begintime(tt0)
c in phytrac call, mu0 and fract are only used in brume
c so we need to pass either rmu0 ou rmu0bar depending on
c moyzon_mu
       if (moyzon_mu) then
         call phytrac (debut,lafin,
     .                 nqmax,nmicro,dtime,appel_chim,dtimechim,
     .                 paprs,pplay,delp,t,rmu0bar,fractbar,pdecli,zls,
     .                 yu1,yv1,zzlev,zzlay,ftsol,
     .                 tr_seri,qaer,d_tr_mph,d_tr_kim,
     .                 fclat,reservoir)
       else
         call phytrac (debut,lafin,
     .                 nqmax,nmicro,dtime,appel_chim,dtimechim,
     .                 paprs,pplay,delp,t,rmu0,fract,pdecli,zls,
     .                 yu1,yv1,zzlev,zzlay,ftsol,
     .                 tr_seri,qaer,d_tr_mph,d_tr_kim,
     .                 fclat,reservoir)
       endif

c         call endtime(tt0,tt1)
c         ttphytra=ttphytra+tt1

c ----- ICI on ajuste radsol en tenant compte du flux de chaleur latente 
c       d'evaporation du reservoir.
c       NOTE : c'est pas tres elegant mais ca permet d'eviter d'aller
c              toucher a clmain.
        if (clouds.eq.1) then
          radsol(:) = radsol(:)+fclat(:)    !test pas de flx de chaleur latente
        endif

        if (microfi.ge.1) then
         tr_seri(:,:,1:nmicro) = tr_seri(:,:,1:nmicro)
     .                        + d_tr_mph(:,:,1:nmicro)*dtime
c        call WriteField_phy('physiq_d_tr_mph01',
c    .                          d_tr_mph(:,:,1),klev)
c        call WriteField_phy('physiq_d_tr_mph10',
c    .                          d_tr_mph(:,:,10),klev)
	endif
c       PAS ELEGANT mais je n'ai pas trouve d'autres solutions :
c       Il semblerait qu'il y ait un probleme lorsque les tendances de traceurs
c       retourne des traceurs nuls et il y a parfois des valeurs negatives qui trainent.
c       Pour ne diffuser le probleme, on force les valeurs negatives a ZERO.
        DO iq=1,nmicro
          DO i=1,klon
            DO l=1,klev
              if (tr_seri(i,l,iq).lt.0.) then
                 tr_seri(i,l,iq) = 0.
              endif
            ENDDO
          ENDDO
        ENDDO

c condensation:
c       NE PAS OUBLIER LA CONDENSATION DES NUAGES !!!!
        if ((clouds.eq.1.or.(chimi)).and.nqmax.gt.nmicro) then
          tr_seri(:,:,nmicro+1:nqmax) = tr_seri(:,:,nmicro+1:nqmax)
     .                         + d_tr_mph(:,:,nmicro+1:nqmax)*dtime
        endif

c chimie:
        if ((chimi).and.(nqmax.gt.nmicro)) then
         tr_seri(:,:,:) = tr_seri(:,:,:) + d_tr_kim(:,:,:)*dtime
	endif

      endif  !iflag_trac=1

c       ch4=0.
c       do l=1,llm
c         ch4(1,l) = tr_seri(1,l,ich4)
c         do j=2,jjm
c           ig0=1+(j-2)*iim
c           do i=1,iim
c             ch4(j,l)= ch4(j,l)  + tr_seri(ig0+i,l,ich4)/iim
c           enddo
c         enddo
c         ch4(jjm+1,l) = tr_seri(klon,l,ich4)
c       enddo
c       do j=1,jjm+1
c         write(501,*) j,ch4(j,1)
c       enddo
c       do l=1,llm
c         write(502,'(I3,49(ES24.17,1X))') l, 
c     &   (ch4(j,l),j=1,jjm+1)
c       enddo
c       write(501,*) ""
c       write(502,*) ""

c------------------
c test condensation
c     do i=1,nqmax 
c       if(tname(i).eq."HCN") then
c          print*,"HCN="
c          do k=1,klev
c           print*,k,tr_seri(klon/2,k,i),d_tr_mph(klon/2,k,i)*dtime
c    v      ,d_tr_kim(klon/2,k,i)*dtime
c          enddo
c          stop
c       endif
c     enddo
c------------------

c====================================================================
c RAYONNEMENT
c====================================================================

      IF (MOD(itaprad,radpas).EQ.0) THEN
c             print*,'RAYONNEMENT ', 
c    $             ' (itaprad=',itaprad,'/radpas=',radpas,')'

c ATTENTION, (klon/2) ne marche pas toujours............
c     print*,"radsol avant radlwsw=",radsol(klon/2)
c     print*,"solsw avant radlwsw=",solsw(klon/2)
c     print*,"sollw avant radlwsw=",sollw(klon/2)
c     print*,"avant radlwsw"

c   ----------------
c   Calcul du rayon moyen des gouttes et des fractions volumique pour le TR
c  ----------------
      IF (clouds.eq.1) THEN
        DO i=1,klon
          DO j=1,klev
            rmcbar(i,j)=rmcbar(i,j)/MAX(REAL(ncount(i,j)),1.)
            xfbar(i,j,:)=xfbar(i,j,:)/MAX(REAL(ncount(i,j)),1.)
          ENDDO 
        ENDDO
      ENDIF
     
c      call begintime(tt0)
      CALL radlwsw
     e            (dist, rmu0, fract, zzlev,
     e             paprs, pplay,ftsol, t_seri, nqmax, nmicro,
     c             tr_seri, qaer)
c     print*,"apres radlwsw"

c      call endtime(tt0,tt1)
c      ttrad=ttrad+tt1

c     print*,"apres radlwsw"
c     mise a zero du rayon moyen des gouttes et des fractions volumique 
      IF (clouds.eq.1) THEN
        rmcbar(:,:)  = 0.
        xfbar(:,:,:) = 0.
        ncount(:,:)  = 0
      ENDIF

      itaprad = 0
      DO k = 1, klev
       DO i = 1, klon
         dtrad(i,k) = heat(i,k)-cool(i,k)     !K/s
       ENDDO
      ENDDO

c     call WriteField_phy('physiq_heat',heat,klev)
c     call WriteField_phy('physiq_cool',cool,klev)

      ENDIF
      itaprad = itaprad + 1
c====================================================================
c
c Ajouter la tendance des rayonnements (tous les pas)
c
      DO k = 1, klev
      DO i = 1, klon
         t_seri(i,k) = t_seri(i,k) + dtrad(i,k) * dtime
      ENDDO
      ENDDO
 
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
      if (ok_orodr.or.ok_gw_nonoro) then
c  CALCUL DE N2
       do i=1,klon
        do k=2,klev
	  ztlev(i,k)  = (t_seri(i,k)+t_seri(i,k-1))/2.
	  zpklev(i,k) = sqrt(ppk(i,k)*ppk(i,k-1))
	enddo
       enddo
       call t2tpot(klon*klev,ztlev, ztetalev,zpklev)
       call t2tpot(klon*klev,t_seri,ztetalay,ppk)
       do i=1,klon
        do k=2,klev
	  zdtetalev(i,k) = ztetalay(i,k)-ztetalay(i,k-1)
	  zdzlev(i,k)    = (zphi(i,k)-zphi(i,k-1))/RG
          zn2(i,k) = RG*zdtetalev(i,k)/(ztetalev(i,k)*zdzlev(i,k))
          zn2(i,k) = max(zn2(i,k),1.e-12)  ! securite
	enddo
        zn2(i,1) = 1.e-12  ! securite
       enddo

      endif
      
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
c A ADAPTER POUR VENUS!!!
        CALL drag_noro(klon,klev,dtime,paprs,pplay,pphi,zn2,
     e                   zmea,zstd, zsig, zgam, zthe,zpic,zval,
     e                   igwd,idx,itest,
     e                   t_seri, u_seri, v_seri,
     s                   zulow, zvlow, zustrdr, zvstrdr,
     s                   d_t_oro, d_u_oro, d_v_oro)

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
c
      ENDIF ! fin de test sur ok_orodr
c
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
c A ADAPTER POUR VENUS ET TITAN!!!
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

        abort_message="Option non developpee pour Titan"
        call abort_gcm(modname,abort_message,1)
c A FAIRE POUR TITAN
c      call flott_gwd_ran(klon,klev,dtime,pplay,zn2,
c     e               t_seri, u_seri, v_seri,
c     o               zustrhi,zvstrhi,
c     o               d_t_hin, d_u_hin, d_v_hin)

c  ajout des tendances

c         t_seri(:,:) = t_seri(:,:) + d_t_hin(:,:)
c         d_t_hin(:,:)= d_t_hin(:,:)/dtime          ! K/s
c         u_seri(:,:) = u_seri(:,:) + d_u_hin(:,:)
c         d_u_hin(:,:)= d_u_hin(:,:)/dtime          ! (m/s)/s
c         v_seri(:,:) = v_seri(:,:) + d_v_hin(:,:)
c         d_v_hin(:,:)= d_v_hin(:,:)/dtime          ! (m/s)/s

      ELSE
         d_t_hin = 0.
         d_u_hin = 0.
         d_v_hin = 0.
         zustrhi = 0.
         zvstrhi = 0.

      ENDIF ! fin de test sur ok_gw_nonoro

c====================================================================
c Transport de ballons 
c====================================================================
      if (ballons.eq.1) then
         CALL ballon(30,pdtphys,rjourvrai,gmtime,
     &               latitude_deg,longitude_deg,
c    C               t,pplay,u,v,pphi)   ! alt above surface (smoothed for GCM)
     C               t,pplay,u,v,zphi)   ! alt above planet average radius
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

      CALL aaam_bud (27,klon,klev,rjourvrai,gmtime,
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
c  A REVOIR POUR VENUS ET TITAN...
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
     . +cpdet(t_seri(i,j-1))/RG*d_t_ec(i,j-1)*delp(i,j-1)
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
c TEST...
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
      
#ifdef CPP_IOIPSL

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

c====================================================================
c Si c'est la fin, il faut conserver l'etat de redemarrage
c====================================================================
c
      IF (lafin) THEN
         itau_phy = itau_phy + itap
         lsinit   = zlsdeg
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

