
!
! $Id: conf_phys.F90 1486 2011-02-11 12:07:39Z fairhead $
!
!
!
module conf_phys_m

   implicit none

contains

  subroutine conf_phys(ok_journe, ok_mensuel, ok_instan, ok_hf, &
                       ok_LES,&
                       callstats,&
                       solarlong0,seuil_inversion, &
                       fact_cldcon, facttemps,ok_newmicro,iflag_radia,&
                       iflag_cldcon, &
                       iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs, &
  		       ok_ade, ok_aie, aerosol_couple, &
                       flag_aerosol, new_aod, &
                       bl95_b0, bl95_b1,&
                       iflag_thermals,nsplit_thermals,tau_thermals, &
                       iflag_thermals_ed,iflag_thermals_optflux, &
                       iflag_coupl,iflag_clos,iflag_wake, read_climoz, &
                       alp_offset)

   use IOIPSL
   USE surface_data
   USE phys_cal_mod
   USE carbon_cycle_mod, ONLY : carbon_cycle_tr, carbon_cycle_cpl
   use control_mod

 include "conema3.h"
 include "fisrtilp.h"
 include "nuage.h"
 include "YOMCST.h"
 include "YOMCST2.h"
!IM : on inclut/initialise les taux de CH4, N2O, CFC11 et CFC12
include "clesphys.h"
include "compbl.h"
include "comsoil.h"
!
! Configuration de la "physique" de LMDZ a l'aide de la fonction
! GETIN de IOIPSL
!
! LF 05/2001
!

!
! type_ocean:      type d'ocean (force, slab, couple)
! version_ocean:   version d'ocean (opa8/nemo pour type_ocean=couple ou 
!                                   sicOBS pour type_ocean=slab)
! ok_veget:   type de modele de vegetation
! ok_journe:  sorties journalieres
! ok_hf:  sorties haute frequence
! ok_mensuel: sorties mensuelles
! ok_instan:  sorties instantanees
! ok_ade, ok_aie: apply or not aerosol direct and indirect effects
! bl95_b*: parameters in the formula to link CDNC to aerosol mass conc 
!


! Sortie:
  logical              :: ok_newmicro
  integer              :: iflag_radia
  logical              :: ok_journe, ok_mensuel, ok_instan, ok_hf
  logical              :: ok_LES
  LOGICAL              :: callstats
  LOGICAL              :: ok_ade, ok_aie, aerosol_couple
  INTEGER              :: flag_aerosol
  LOGICAL              :: new_aod
  REAL                 :: bl95_b0, bl95_b1
  real                 :: fact_cldcon, facttemps,ratqsbas,ratqshaut,tau_ratqs
  integer              :: iflag_cldcon
  integer              :: iflag_ratqs

  character (len = 6),SAVE  :: type_ocean_omp, version_ocean_omp, ocean_omp
  CHARACTER(len = 8),SAVE   :: aer_type_omp
  logical,SAVE              :: ok_veget_omp, ok_newmicro_omp
  logical,SAVE        :: ok_journe_omp, ok_mensuel_omp, ok_instan_omp, ok_hf_omp        
  logical,SAVE        :: ok_LES_omp   
  LOGICAL,SAVE        :: callstats_omp
  LOGICAL,SAVE        :: ok_ade_omp, ok_aie_omp, aerosol_couple_omp
  INTEGER, SAVE       :: flag_aerosol_omp
  LOGICAL, SAVE       :: new_aod_omp
  REAL,SAVE           :: bl95_b0_omp, bl95_b1_omp
  REAL,SAVE           :: freq_ISCCP_omp, ecrit_ISCCP_omp
  REAL,SAVE           :: freq_COSP_omp
  real,SAVE           :: fact_cldcon_omp, facttemps_omp,ratqsbas_omp
  real,SAVE           :: ratqshaut_omp
  real,SAVE           :: tau_ratqs_omp
  integer,SAVE        :: iflag_radia_omp
  integer,SAVE        :: iflag_rrtm_omp
  integer,SAVE        :: iflag_cldcon_omp, ip_ebil_phy_omp
  integer,SAVE        :: iflag_ratqs_omp

  Real,SAVE           :: f_cdrag_ter_omp,f_cdrag_oce_omp
  Real,SAVE           :: f_rugoro_omp   

! Local
  integer              :: numout = 6
  real                 :: zzz

  real :: seuil_inversion
  real,save :: seuil_inversion_omp

  integer :: iflag_thermals,nsplit_thermals
  integer,SAVE :: iflag_thermals_ed_omp,iflag_thermals_optflux_omp
  integer :: iflag_thermals_ed,iflag_thermals_optflux
  integer,SAVE :: iflag_thermals_omp,nsplit_thermals_omp
  real :: tau_thermals
  real,save :: tau_thermals_omp
  integer :: iflag_coupl
  integer :: iflag_clos
  integer :: iflag_wake
  real :: alp_offset
  REAL, SAVE :: alp_offset_omp
  integer,SAVE :: iflag_coupl_omp,iflag_clos_omp,iflag_wake_omp
  integer,SAVE :: iflag_cvl_sigd_omp
  REAL, SAVE :: supcrit1_omp, supcrit2_omp
  INTEGER, SAVE :: iflag_mix_omp
  real, save :: scut_omp, qqa1_omp, qqa2_omp, gammas_omp, Fmax_omp, alphas_omp

  REAL,SAVE :: R_ecc_omp,R_peri_omp,R_incl_omp,solaire_omp,co2_ppm_omp
  REAL,SAVE :: RCO2_omp,CH4_ppb_omp,RCH4_omp,N2O_ppb_omp,RN2O_omp,CFC11_ppt_omp
  REAL,SAVE :: RCFC11_omp,CFC12_ppt_omp,RCFC12_omp,epmax_omp
  LOGICAL,SAVE :: ok_adj_ema_omp
  INTEGER,SAVE :: iflag_clw_omp
  REAL,SAVE :: cld_lc_lsc_omp,cld_lc_con_omp,cld_tau_lsc_omp,cld_tau_con_omp
  REAL,SAVE :: ffallv_lsc_omp, ffallv_con_omp,coef_eva_omp
  LOGICAL,SAVE :: reevap_ice_omp
  INTEGER,SAVE :: iflag_pdf_omp
  REAL,SAVE :: rad_froid_omp, rad_chau1_omp, rad_chau2_omp
  REAL,SAVE :: t_glace_min_omp, t_glace_max_omp
  REAL,SAVE :: inertie_sol_omp,inertie_sno_omp,inertie_ice_omp
  REAL,SAVE :: qsol0_omp
  REAL      :: solarlong0
  REAL,SAVE :: solarlong0_omp
  INTEGER,SAVE :: top_height_omp,overlap_omp
  REAL,SAVE :: cdmmax_omp,cdhmax_omp,ksta_omp,ksta_ter_omp
  LOGICAL,SAVE :: ok_kzmin_omp
  REAL, SAVE ::  fmagic_omp, pmagic_omp
  INTEGER,SAVE :: iflag_pbl_omp,lev_histhf_omp,lev_histday_omp,lev_histmth_omp
  Integer, save :: lev_histins_omp, lev_histLES_omp 
  INTEGER, SAVE :: lev_histdayNMC_omp
  LOGICAL, SAVE :: ok_histNMC_omp(3)
  REAL, SAVE :: freq_outNMC_omp(3), freq_calNMC_omp(3)
  CHARACTER*4, SAVE :: type_run_omp
  LOGICAL,SAVE :: ok_isccp_omp
  LOGICAL,SAVE :: ok_cosp_omp
  LOGICAL,SAVE :: ok_mensuelCOSP_omp,ok_journeCOSP_omp,ok_hfCOSP_omp
  REAL,SAVE :: lonmin_ins_omp, lonmax_ins_omp, latmin_ins_omp, latmax_ins_omp
  REAL,SAVE :: ecrit_hf_omp, ecrit_day_omp, ecrit_mth_omp, ecrit_reg_omp
  REAL,SAVE :: ecrit_ins_omp
  REAL,SAVE :: ecrit_LES_omp
  REAL,SAVE :: ecrit_tra_omp
  REAL,SAVE :: cvl_corr_omp
  LOGICAL,SAVE :: ok_lic_melt_omp
!
  LOGICAL,SAVE :: cycle_diurne_omp,soil_model_omp,new_oliq_omp
  LOGICAL,SAVE :: ok_orodr_omp, ok_orolf_omp, ok_limitvrai_omp
  INTEGER, SAVE :: nbapp_rad_omp, iflag_con_omp
  LOGICAL,SAVE :: ok_strato_omp
  LOGICAL,SAVE :: ok_hines_omp
  LOGICAL,SAVE      :: carbon_cycle_tr_omp
  LOGICAL,SAVE      :: carbon_cycle_cpl_omp

  integer, intent(out):: read_climoz ! read ozone climatology, OpenMP shared
  ! Allowed values are 0, 1 and 2
  ! 0: do not read an ozone climatology
  ! 1: read a single ozone climatology that will be used day and night
  ! 2: read two ozone climatologies, the average day and night
  ! climatology and the daylight climatology

!$OMP MASTER 
!Config Key  = type_ocean 
!Config Desc = Type d'ocean
!Config Def  = force
!Config Help = Type d'ocean utilise: force, slab,couple
!
  type_ocean_omp = 'force '
  call getin('type_ocean', type_ocean_omp)
!
!Config Key  = version_ocean 
!Config Desc = Version d'ocean
!Config Def  = xxxxxx
!Config Help = Version d'ocean utilise: opa8/nemo/sicOBS/xxxxxx
!
  version_ocean_omp = 'xxxxxx'
  call getin('version_ocean', version_ocean_omp)

!Config Key  = OCEAN
!Config Desc = Old parameter name for type_ocean
!Config Def  = yyyyyy
!Config Help = This is only for testing purpose
!
  ocean_omp = 'yyyyyy'
  call getin('OCEAN', ocean_omp)
  IF (ocean_omp /= 'yyyyyy') THEN
     WRITE(numout,*)'ERROR!! Old variable name OCEAN used in parmeter file.'
     WRITE(numout,*)'Variable OCEAN has been replaced by the variable type_ocean.'
     WRITE(numout,*)'You have to update your parameter file physiq.def to succed running'
     CALL abort_gcm('conf_phys','Variable OCEAN no longer existing, use variable name type_ocean',1)
  END IF

!
!Config Key  = VEGET 
!Config Desc = Type de modele de vegetation
!Config Def  = .false.
!Config Help = Type de modele de vegetation utilise
!
  ok_veget_omp = .false.
  call getin('VEGET', ok_veget_omp)
!
!Config Key  = OK_journe
!Config Desc = Pour des sorties journalieres 
!Config Def  = .false.
!Config Help = Pour creer le fichier histday contenant les sorties
!              journalieres 
!
  ok_journe_omp = .false.
  call getin('OK_journe', ok_journe_omp)
!
!Config Key  = ok_hf
!Config Desc = Pour des sorties haute frequence
!Config Def  = .false.
!Config Help = Pour creer le fichier histhf contenant les sorties
!              haute frequence ( 3h ou 6h)
!
  ok_hf_omp = .false.
  call getin('ok_hf', ok_hf_omp)
!
!Config Key  = OK_mensuel
!Config Desc = Pour des sorties mensuelles 
!Config Def  = .true.
!Config Help = Pour creer le fichier histmth contenant les sorties
!              mensuelles 
!
  ok_mensuel_omp = .true.
  call getin('OK_mensuel', ok_mensuel_omp)
!
!Config Key  = OK_instan
!Config Desc = Pour des sorties instantanees 
!Config Def  = .false.
!Config Help = Pour creer le fichier histins contenant les sorties
!              instantanees 
!
  ok_instan_omp = .false.
  call getin('OK_instan', ok_instan_omp)
!
!Config Key  = ok_ade
!Config Desc = Aerosol direct effect or not?
!Config Def  = .false.
!Config Help = Used in radlwsw.F
!
  ok_ade_omp = .false.
  call getin('ok_ade', ok_ade_omp)

!
!Config Key  = ok_aie
!Config Desc = Aerosol indirect effect or not?
!Config Def  = .false.
!Config Help = Used in nuage.F and radlwsw.F
!
  ok_aie_omp = .false.
  call getin('ok_aie', ok_aie_omp)

!
!Config Key  = aerosol_couple
!Config Desc = read aerosol in file or calcul by inca
!Config Def  = .false.
!Config Help = Used in physiq.F
!
  aerosol_couple_omp = .false.
  CALL getin('aerosol_couple',aerosol_couple_omp)

!
!Config Key  = flag_aerosol
!Config Desc = which aerosol is use for coupled model
!Config Def  = 1
!Config Help = Used in physiq.F
!
! - flag_aerosol=1 => so4 only (defaut) 
! - flag_aerosol=2 => bc  only 
! - flag_aerosol=3 => pom only
! - flag_aerosol=4 => seasalt only 
! - flag_aerosol=5 => dust only
! - flag_aerosol=6 => all aerosol

  flag_aerosol_omp = 1
  CALL getin('flag_aerosol',flag_aerosol_omp)

! Temporary variable for testing purpose!!
!Config Key  = new_aod
!Config Desc = which calcul of aeropt
!Config Def  = false
!Config Help = Used in physiq.F
!
  new_aod_omp = .true.
  CALL getin('new_aod',new_aod_omp)

! 
!Config Key  = aer_type 
!Config Desc = Use a constant field for the aerosols 
!Config Def  = scenario 
!Config Help = Used in readaerosol.F90 
! 
  aer_type_omp = 'scenario' 
  call getin('aer_type', aer_type_omp) 

!
!Config Key  = bl95_b0
!Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
!Config Def  = .false.
!Config Help = Used in nuage.F
!
  bl95_b0_omp = 2.
  call getin('bl95_b0', bl95_b0_omp)

!Config Key  = bl95_b1
!Config Desc = Parameter in CDNC-maer link (Boucher&Lohmann 1995)
!Config Def  = .false.
!Config Help = Used in nuage.F
!
  bl95_b1_omp = 0.2
  call getin('bl95_b1', bl95_b1_omp)

!Config Key  = freq_ISCCP
!Config Desc = Frequence d'appel du simulateur ISCCP en secondes;
!              par defaut 10800, i.e. 3 heures 
!Config Def  = 10800.
!Config Help = Used in ini_histISCCP.h
!
  freq_ISCCP_omp = 10800.
  call getin('freq_ISCCP', freq_ISCCP_omp)
!
!Config Key  = ecrit_ISCCP
!Config Desc = Frequence d'ecriture des resultats du simulateur ISCCP en nombre de jours;
!              par defaut 1., i.e. 1 jour
!Config Def  = 1.
!Config Help = Used in ini_histISCCP.h
!
!
  ecrit_ISCCP_omp = 1.
  call getin('ecrit_ISCCP', ecrit_ISCCP_omp)

!Config Key  = freq_COSP
!Config Desc = Frequence d'appel du simulateur COSP en secondes;
!              par defaut 10800, i.e. 3 heures
!Config Def  = 10800.
!Config Help = Used in ini_histdayCOSP.h
!
  freq_COSP_omp = 10800.
  call getin('freq_COSP', freq_COSP_omp)

!
!Config Key  = ip_ebil_phy
!Config Desc = Niveau de sortie pour les diags bilan d'energie 
!Config Def  = 0
!Config Help = 
!               
  ip_ebil_phy_omp = 0
  call getin('ip_ebil_phy', ip_ebil_phy_omp)
!
!Config Key  = seuil_inversion
!Config Desc = Seuil ur dTh pour le choix entre les schemas de CL
!Config Def  = -0.1
!Config Help = 
!               
  seuil_inversion_omp = -0.1
  call getin('seuil_inversion', seuil_inversion_omp)

!!
!! Constante solaire & Parametres orbitaux & taux gaz effet de serre BEG
!!
!Config Key  = R_ecc
!Config Desc = Excentricite
!Config Def  = 0.016715
!Config Help = 
!               
!valeur AMIP II
  R_ecc_omp = 0.016715
  call getin('R_ecc', R_ecc_omp)
!!
!Config Key  = R_peri
!Config Desc = Equinoxe
!Config Def  = 
!Config Help = 
!               
!
!valeur AMIP II
  R_peri_omp = 102.7
  call getin('R_peri', R_peri_omp)
!!
!Config Key  = R_incl
!Config Desc = Inclinaison
!Config Def  = 
!Config Help = 
!               
!
!valeur AMIP II
  R_incl_omp = 23.441
  call getin('R_incl', R_incl_omp)
!!
!Config Key  = solaire
!Config Desc = Constante solaire en W/m2
!Config Def  = 1365.
!Config Help = 
!               
!
!valeur AMIP II
  solaire_omp = 1365.
  call getin('solaire', solaire_omp)
!!
!Config Key  = co2_ppm
!Config Desc = concentration du gaz carbonique en ppmv
!Config Def  = 348.
!Config Help = 
!               
!
!valeur AMIP II
  co2_ppm_omp = 348.
  call getin('co2_ppm', co2_ppm_omp)
!!
!Config Key  = RCO2
!Config Desc = Concentration du CO2
!Config Def  = co2_ppm * 1.0e-06  * 44.011/28.97
!Config Def  = 348. * 1.0e-06  * 44.011/28.97
!Config Help = 
!               
! RCO2 = 5.286789092164308E-04
!ancienne valeur
  RCO2_omp = co2_ppm_omp * 1.0e-06  * 44.011/28.97 ! pour co2_ppm=348.

!!  call getin('RCO2', RCO2)
!!
!Config Key  = RCH4
!Config Desc = Concentration du CH4
!Config Def  = 1.65E-06* 16.043/28.97
!Config Help = 
!               
!
!valeur AMIP II
!OK  RCH4 = 1.65E-06* 16.043/28.97
! RCH4 = 9.137366240938903E-07
!
!ancienne valeur
! RCH4 = 1.72E-06* 16.043/28.97
!OK call getin('RCH4', RCH4)
  zzz = 1650.
  call getin('CH4_ppb', zzz)
  CH4_ppb_omp = zzz
  RCH4_omp = CH4_ppb_omp * 1.0E-09 * 16.043/28.97
!!
!Config Key  = RN2O
!Config Desc = Concentration du N2O
!Config Def  = 306.E-09* 44.013/28.97
!Config Help = 
!               
!
!valeur AMIP II
!OK  RN2O = 306.E-09* 44.013/28.97
! RN2O = 4.648939592682085E-07
!
!ancienne valeur
! RN2O = 310.E-09* 44.013/28.97
!OK  call getin('RN2O', RN2O)
  zzz=306.
  call getin('N2O_ppb', zzz)
  N2O_ppb_omp = zzz
  RN2O_omp = N2O_ppb_omp * 1.0E-09 * 44.013/28.97
!!
!Config Key  = RCFC11
!Config Desc = Concentration du CFC11
!Config Def  = 280.E-12* 137.3686/28.97
!Config Help = 
!               
!
!OK RCFC11 = 280.E-12* 137.3686/28.97
  zzz = 280.
  call getin('CFC11_ppt',zzz)
  CFC11_ppt_omp = zzz
  RCFC11_omp=CFC11_ppt_omp* 1.0E-12 * 137.3686/28.97
! RCFC11 = 1.327690990680013E-09
!OK call getin('RCFC11', RCFC11)
!!
!Config Key  = RCFC12
!Config Desc = Concentration du CFC12
!Config Def  = 484.E-12* 120.9140/28.97
!Config Help = 
!               
!
!OK RCFC12 = 484.E-12* 120.9140/28.97
  zzz = 484.
  call getin('CFC12_ppt',zzz)
  CFC12_ppt_omp = zzz
  RCFC12_omp = CFC12_ppt_omp * 1.0E-12 * 120.9140/28.97
! RCFC12 = 2.020102726958923E-09
!OK call getin('RCFC12', RCFC12)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! Constantes precedemment dans dyn3d/conf_gcm

!Config  Key  = cycle_diurne
!Config  Desc = Cycle ddiurne
!Config  Def  = y
!Config  Help = Cette option permet d'eteidre le cycle diurne.
!Config         Peut etre util pour accelerer le code !
       cycle_diurne_omp = .TRUE.
       CALL getin('cycle_diurne',cycle_diurne_omp)

!Config  Key  = soil_model
!Config  Desc = Modele de sol
!Config  Def  = y
!Config  Help = Choix du modele de sol (Thermique ?)
!Config         Option qui pourait un string afin de pouvoir
!Config         plus de choix ! Ou meme une liste d'options !
       soil_model_omp = .TRUE.
       CALL getin('soil_model',soil_model_omp)

!Config  Key  = new_oliq
!Config  Desc = Nouvelle eau liquide
!Config  Def  = y
!Config  Help = Permet de mettre en route la
!Config         nouvelle parametrisation de l'eau liquide !
       new_oliq_omp = .TRUE.
       CALL getin('new_oliq',new_oliq_omp)

!Config  Key  = ok_orodr
!Config  Desc = Orodr ???
!Config  Def  = y
!Config  Help = Y en a pas comprendre !
!Config         
       ok_orodr_omp = .TRUE.
       CALL getin('ok_orodr',ok_orodr_omp)

!Config  Key  =  ok_orolf
!Config  Desc = Orolf ??
!Config  Def  = y
!Config  Help = Connais pas !
       ok_orolf_omp = .TRUE.
       CALL getin('ok_orolf', ok_orolf_omp)

!Config  Key  = ok_limitvrai
!Config  Desc = Force la lecture de la bonne annee
!Config  Def  = n
!Config  Help = On peut forcer le modele a lire le
!Config         fichier SST de la bonne annee. C'est une tres bonne
!Config         idee, pourquoi ne pas mettre toujours a y ???
       ok_limitvrai_omp = .FALSE.
       CALL getin('ok_limitvrai',ok_limitvrai_omp)

!Config  Key  = nbapp_rad
!Config  Desc = Frequence d'appel au rayonnement
!Config  Def  = 12
!Config  Help = Nombre  d'appels des routines de rayonnements
!Config         par jour.
       nbapp_rad_omp = 12
       CALL getin('nbapp_rad',nbapp_rad_omp)

!Config  Key  = iflag_con
!Config  Desc = Flag de convection
!Config  Def  = 2
!Config  Help = Flag  pour la convection les options suivantes existent :
!Config         1 pour LMD,
!Config         2 pour Tiedtke,
!Config         3 pour CCM(NCAR)  
       iflag_con_omp = 2
       CALL getin('iflag_con',iflag_con_omp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Constante solaire & Parametres orbitaux & taux gaz effet de serre END
!!
!! KE
!

!Config key  = cvl_corr
!Config Desc = Facteur multiplication des precip convectives dans KE
!Config Def  = 1.00
!Config Help = 1.02 pour un moderne ou un pre-ind. A ajuster pour un glaciaire
  cvl_corr_omp = 1.00
  CALL getin('cvl_corr', cvl_corr_omp)


!Config Key  = epmax
!Config Desc = Efficacite precip
!Config Def  = 0.993
!Config Help = 
!
  epmax_omp = .993
  call getin('epmax', epmax_omp)
!
!Config Key  = ok_adj_ema
!Config Desc =  
!Config Def  = false
!Config Help = 
!
  ok_adj_ema_omp = .false.
  call getin('ok_adj_ema',ok_adj_ema_omp)
!
!Config Key  = iflag_clw
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_clw_omp = 0
  call getin('iflag_clw',iflag_clw_omp)
!
!Config Key  = cld_lc_lsc 
!Config Desc =  
!Config Def  = 2.6e-4
!Config Help = 
!
  cld_lc_lsc_omp = 2.6e-4
  call getin('cld_lc_lsc',cld_lc_lsc_omp)
!
!Config Key  = cld_lc_con
!Config Desc =  
!Config Def  = 2.6e-4
!Config Help = 
!
  cld_lc_con_omp = 2.6e-4
  call getin('cld_lc_con',cld_lc_con_omp)
!
!Config Key  = cld_tau_lsc
!Config Desc =  
!Config Def  = 3600.
!Config Help = 
!
  cld_tau_lsc_omp = 3600.
  call getin('cld_tau_lsc',cld_tau_lsc_omp)
!
!Config Key  = cld_tau_con
!Config Desc =  
!Config Def  = 3600.
!Config Help = 
!
  cld_tau_con_omp = 3600.
  call getin('cld_tau_con',cld_tau_con_omp)
!
!Config Key  = ffallv_lsc
!Config Desc =  
!Config Def  = 1.
!Config Help = 
!
  ffallv_lsc_omp = 1.
  call getin('ffallv_lsc',ffallv_lsc_omp)
!
!Config Key  = ffallv_con
!Config Desc =  
!Config Def  = 1.
!Config Help = 
!
  ffallv_con_omp = 1.
  call getin('ffallv_con',ffallv_con_omp)
!
!Config Key  = coef_eva
!Config Desc =  
!Config Def  = 2.e-5
!Config Help = 
!
  coef_eva_omp = 2.e-5
  call getin('coef_eva',coef_eva_omp)
!
!Config Key  = reevap_ice
!Config Desc =  
!Config Def  = .false.
!Config Help = 
!
  reevap_ice_omp = .false.
  call getin('reevap_ice',reevap_ice_omp)

!Config Key  = iflag_ratqs
!Config Desc =
!Config Def  = 1
!Config Help =
!
  iflag_ratqs_omp = 1
  call getin('iflag_ratqs',iflag_ratqs_omp)

!
!Config Key  = iflag_radia 
!Config Desc =  
!Config Def  = 1
!Config Help = 
!
  iflag_radia_omp = 1
  call getin('iflag_radia',iflag_radia_omp)

!
!Config Key  = iflag_rrtm 
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_rrtm_omp = 0
  call getin('iflag_rrtm',iflag_rrtm_omp)

!
!Config Key  = iflag_cldcon 
!Config Desc =  
!Config Def  = 1
!Config Help = 
!
  iflag_cldcon_omp = 1
  call getin('iflag_cldcon',iflag_cldcon_omp)

!
!Config Key  = iflag_pdf 
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_pdf_omp = 0
  call getin('iflag_pdf',iflag_pdf_omp)
!
!Config Key  = fact_cldcon
!Config Desc =  
!Config Def  = 0.375
!Config Help = 
!
  fact_cldcon_omp = 0.375
  call getin('fact_cldcon',fact_cldcon_omp)

!
!Config Key  = facttemps
!Config Desc =  
!Config Def  = 1.e-4
!Config Help = 
!
  facttemps_omp = 1.e-4
  call getin('facttemps',facttemps_omp)

!
!Config Key  = ok_newmicro
!Config Desc =  
!Config Def  = .true.
!Config Help = 
!
  ok_newmicro_omp = .true.
  call getin('ok_newmicro',ok_newmicro_omp)
!
!Config Key  = ratqsbas
!Config Desc =  
!Config Def  = 0.01
!Config Help = 
!
  ratqsbas_omp = 0.01
  call getin('ratqsbas',ratqsbas_omp)
!
!Config Key  = ratqshaut
!Config Desc =  
!Config Def  = 0.3
!Config Help = 
!
  ratqshaut_omp = 0.3
  call getin('ratqshaut',ratqshaut_omp)

!Config Key  = tau_ratqs
!Config Desc =  
!Config Def  = 1800.
!Config Help = 
!
  tau_ratqs_omp = 1800.
  call getin('tau_ratqs',tau_ratqs_omp)

!
!-----------------------------------------------------------------------
! Longitude solaire pour le calcul de l'ensoleillement en degre
! si on veut imposer la saison. Sinon, solarlong0=-999.999
!Config Key  = solarlong0
!Config Desc =  
!Config Def  = -999.999 
!Config Help = 
!
  solarlong0_omp = -999.999
  call getin('solarlong0',solarlong0_omp)
!
!-----------------------------------------------------------------------
!  Valeur imposee de l'humidite du sol pour le modele bucket.
!Config Key  = qsol0
!Config Desc =  
!Config Def  = -1.
!Config Help = 
!
  qsol0_omp = -1.
  call getin('qsol0',qsol0_omp)
!
!-----------------------------------------------------------------------
!
!Config Key  = inertie_ice
!Config Desc =  
!Config Def  = 2000.
!Config Help = 
!
  inertie_ice_omp = 2000.
  call getin('inertie_ice',inertie_ice_omp)
!
!Config Key  = inertie_sno
!Config Desc =  
!Config Def  = 2000.
!Config Help = 
!
  inertie_sno_omp = 2000.
  call getin('inertie_sno',inertie_sno_omp)
!
!Config Key  = inertie_sol
!Config Desc =  
!Config Def  = 2000.
!Config Help = 
!
  inertie_sol_omp = 2000.
  call getin('inertie_sol',inertie_sol_omp)

!
!Config Key  = rad_froid
!Config Desc =  
!Config Def  = 35.0
!Config Help = 
!
  rad_froid_omp = 35.0
  call getin('rad_froid',rad_froid_omp)

!
!Config Key  = rad_chau1
!Config Desc =  
!Config Def  = 13.0
!Config Help = 
!
  rad_chau1_omp = 13.0
  call getin('rad_chau1',rad_chau1_omp)

!
!Config Key  = rad_chau2
!Config Desc =  
!Config Def  = 9.0
!Config Help = 
!
  rad_chau2_omp = 9.0
  call getin('rad_chau2',rad_chau2_omp)

!
!Config Key  = t_glace_min
!Config Desc =  
!Config Def  = 258.
!Config Help = 
!
  t_glace_min_omp = 258.
  call getin('t_glace_min',t_glace_min_omp)

!
!Config Key  = t_glace_max
!Config Desc =  
!Config Def  = 273.13
!Config Help = 
!
  t_glace_max_omp = 273.13
  call getin('t_glace_max',t_glace_max_omp)

!
!Config Key  = top_height
!Config Desc =
!Config Def  = 3
!Config Help =
!
  top_height_omp = 3
  call getin('top_height',top_height_omp)

!
!Config Key  = overlap
!Config Desc =
!Config Def  = 3
!Config Help =
!
  overlap_omp = 3
  call getin('overlap',overlap_omp)


!
!
!Config Key  = cdmmax
!Config Desc =
!Config Def  = 1.3E-3
!Config Help =
!
  cdmmax_omp = 1.3E-3
  call getin('cdmmax',cdmmax_omp)

!
!Config Key  = cdhmax
!Config Desc =
!Config Def  = 1.1E-3
!Config Help =
!
  cdhmax_omp = 1.1E-3
  call getin('cdhmax',cdhmax_omp)

!261103
!
!Config Key  = ksta
!Config Desc =
!Config Def  = 1.0e-10
!Config Help =
!
  ksta_omp = 1.0e-10
  call getin('ksta',ksta_omp)

!
!Config Key  = ksta_ter
!Config Desc =
!Config Def  = 1.0e-10
!Config Help =
!
  ksta_ter_omp = 1.0e-10
  call getin('ksta_ter',ksta_ter_omp)

!
!Config Key  = ok_kzmin
!Config Desc =
!Config Def  = .true.
!Config Help =
!
  ok_kzmin_omp = .true.
  call getin('ok_kzmin',ok_kzmin_omp)

!
!Config Key  = fmagic
!Config Desc = additionnal multiplicator factor used for albedo
!Config Def  = 1.
!Config Help = additionnal multiplicator factor used in albedo.F
!
  fmagic_omp = 1.
  call getin('fmagic',fmagic_omp)

!
!Config Key  = pmagic
!Config Desc = additional factor used for albedo
!Config Def  = 0.
!Config Help = additional factor used in albedo.F
!
  pmagic_omp = 0.
  call getin('pmagic',pmagic_omp)


!Config Key = ok_lic_melt
!Config Desc = Prise en compte de la fonte de la calotte dans le bilan d'eau
!Config Def  = .false.
!Config Help = mettre a .false. pour assurer la conservation en eau
  ok_lic_melt_omp = .false.
  call getin('ok_lic_melt', ok_lic_melt_omp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR THE PLANETARY BOUNDARY LAYER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Config Key  = iflag_pbl
!Config Desc =
!Config Def  = 1
!Config Help =
!
  iflag_pbl_omp = 1
  call getin('iflag_pbl',iflag_pbl_omp)
!
!Config Key  = iflag_thermals
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_thermals_omp = 0
  call getin('iflag_thermals',iflag_thermals_omp)
!
!
!Config Key  = iflag_thermals_ed
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_thermals_ed_omp = 0
  call getin('iflag_thermals_ed',iflag_thermals_ed_omp)
!
!
!Config Key  = iflag_thermals_optflux
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_thermals_optflux_omp = 0
  call getin('iflag_thermals_optflux',iflag_thermals_optflux_omp)
!
!
!Config Key  = nsplit_thermals
!Config Desc =
!Config Def  = 1
!Config Help =
!
  nsplit_thermals_omp = 1
  call getin('nsplit_thermals',nsplit_thermals_omp)

!Config Key  = tau_thermals
!Config Desc =
!Config Def  = 0.
!Config Help =
!
  tau_thermals_omp = 0.
  call getin('tau_thermals',tau_thermals_omp)

!
!Config Key  = iflag_coupl
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_coupl_omp = 0
  call getin('iflag_coupl',iflag_coupl_omp)

!
!Config Key  = iflag_clos
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_clos_omp = 1
  call getin('iflag_clos',iflag_clos_omp)
!
!Config Key  = iflag_cvl_sigd
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_cvl_sigd_omp = 0
  call getin('iflag_cvl_sigd',iflag_cvl_sigd_omp)

!Config Key  = iflag_wake
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  iflag_wake_omp = 0
  call getin('iflag_wake',iflag_wake_omp)

!Config Key  = alp_offset
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  alp_offset_omp = 0.
  call getin('alp_offset',alp_offset_omp)

!
!Config Key  = lev_histhf
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histhf_omp = 1
  call getin('lev_histhf',lev_histhf_omp)

!
!Config Key  = lev_histday
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histday_omp = 1
  call getin('lev_histday',lev_histday_omp)

!
!Config Key  = lev_histmth
!Config Desc =
!Config Def  = 2
!Config Help =
!
  lev_histmth_omp = 2
  call getin('lev_histmth',lev_histmth_omp)
!
!Config Key  = lev_histins
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histins_omp = 1
  call getin('lev_histins',lev_histins_omp)
  !
!Config Key  = lev_histLES
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histLES_omp = 1
  call getin('lev_histLES',lev_histLES_omp)
!
!Config Key  = lev_histdayNMC
!Config Desc =
!Config Def  = 8
!Config Help =
!
  lev_histdayNMC_omp = 8
  call getin('lev_histdayNMC',lev_histdayNMC_omp)
!
!histNMC BEG
!Config Key  = ok_histNMC
!Config Desc = ok_histNMC(1) = frequence de sortie fichiers histmthNMC
!Config Desc = ok_histNMC(2) = frequence de sortie fichiers histdayNMC
!Config Desc = ok_histNMC(3) = frequence de sortie fichiers histhfNMC
!Config Def  = n, n, n
!Config Help =
!
  ok_histNMC_omp(1) = .false.
  ok_histNMC_omp(2) = .false.
  ok_histNMC_omp(3) = .false.
  call getin('ok_histNMC',ok_histNMC_omp)
!
!Config Key  = freq_outNMC
!Config Desc = freq_outNMC(1) = frequence de sortie fichiers histmthNMC
!Config Desc = freq_outNMC(2) = frequence de sortie fichiers histdayNMC
!Config Desc = freq_outNMC(3) = frequence de sortie fichiers histhfNMC
!Config Def  = 2592000., 86400., 21600.
!Config Help =
!
! freq_outNMC_omp(1) = 2592000.
  freq_outNMC_omp(1) = mth_len*86400.
  freq_outNMC_omp(2) = 86400.
  freq_outNMC_omp(3) = 21600.
  call getin('freq_outNMC',freq_outNMC_omp)
!
!Config Key  = freq_calNMC
!Config Desc = freq_calNMC(1) = frequence de calcul fichiers histmthNMC
!Config Desc = freq_calNMC(2) = frequence de calcul fichiers histdayNMC
!Config Desc = freq_calNMC(3) = frequence de calcul fichiers histhfNMC
!Config Def  = pasphys
!Config Help =
!
  freq_calNMC_omp(1) = pasphys
  freq_calNMC_omp(2) = pasphys
  freq_calNMC_omp(3) = pasphys
  call getin('freq_calNMC',freq_calNMC_omp)
!
!Config Key  = type_run
!Config Desc =
!Config Def  = 'AMIP'/'CFMIP'  ou 'CLIM'/'ENSP'
!Config Help =
!
  type_run_omp = 'AMIP'
  call getin('type_run',type_run_omp)

!
!Config Key  = ok_isccp
!Config Desc =
!Config Def  = .true.
!Config Help =
!
! ok_isccp = .true.
  ok_isccp_omp = .false.
  call getin('ok_isccp',ok_isccp_omp)

!
!Config Key  = ok_cosp
!Config Desc =
!Config Def  = .false.
!Config Help =
!
  ok_cosp_omp = .false.
  call getin('ok_cosp',ok_cosp_omp)

!
!Config Key  = ok_mensuelCOSP
!Config Desc =
!Config Def  = .true.
!Config Help =
!
  ok_mensuelCOSP_omp = .true.
  call getin('ok_mensuelCOSP',ok_mensuelCOSP_omp)

!
!Config Key  = ok_journeCOSP
!Config Desc =
!Config Def  = .true.
!Config Help = 
!
  ok_journeCOSP_omp = .true.
  call getin('ok_journeCOSP',ok_journeCOSP_omp)

!
!Config Key  = ok_hfCOSP
!Config Desc =
!Config Def  = .false.
!Config Help =
!
  ok_hfCOSP_omp = .false.
  call getin('ok_hfCOSP',ok_hfCOSP_omp)

!
! coordonnees (lonmin_ins, lonmax_ins, latmin_ins, latmax_ins) pour la zone 
! avec sorties instantannees tous les pas de temps de la physique => "histbilKP_ins.nc"
!
!Config Key  = lonmin_ins
!Config Desc = 100.  
!Config Def  = longitude minimale sorties "bilKP_ins"
!Config Help = 
!
   lonmin_ins_omp = 100.
   call getin('lonmin_ins',lonmin_ins_omp)
!
!Config Key  = lonmax_ins
!Config Desc = 130. 
!Config Def  = longitude maximale sorties "bilKP_ins"
!Config Help =
!
   lonmax_ins_omp = 130.
   call getin('lonmax_ins',lonmax_ins_omp)
!
!Config Key  = latmin_ins
!Config Desc = -20.  
!Config Def  = latitude minimale sorties "bilKP_ins"
!Config Help = 
!
   latmin_ins_omp = -20.
   call getin('latmin_ins',latmin_ins_omp)
!
!Config Key  = latmax_ins
!Config Desc = 20. 
!Config Def  = latitude maximale sorties "bilKP_ins"
!Config Help =
!
   latmax_ins_omp = 20.
   call getin('latmax_ins',latmax_ins_omp)
!
!Config Key  = ecrit_hf
!Config Desc =
!Config Def  = 1./8. !toutes les 3h
!Config Help =
!
  ecrit_hf_omp = 1./8.
  call getin('ecrit_hf',ecrit_hf_omp)
!
!Config Key  = ecrit_ins
!Config Desc =
!Config Def  = 1./48. ! toutes les 1/2 h
!Config Help =
!
  ecrit_ins_omp = 1./48.
  call getin('ecrit_ins',ecrit_ins_omp)
!
!Config Key  = ecrit_day
!Config Desc =
!Config Def  = 1.0 !tous les jours
!Config Help = nombre de jours pour ecriture fichier histday.nc
!
  ecrit_day_omp = 1.0
  call getin('ecrit_day',ecrit_day_omp)
!
!Config Key  = ecrit_mth
!Config Desc =
!Config Def  = 30. !tous les 30jours (1 fois par mois)
!Config Help =
!
  ecrit_mth_omp = 30.
  call getin('ecrit_mth',ecrit_mth_omp)
!
!Config Key  = ecrit_tra
!Config Desc =
!Config Def  = 30. !tous les 30jours (1 fois par mois)
!Config Help =
!
  ecrit_tra_omp = 30.
  call getin('ecrit_tra',ecrit_tra_omp)
!
!Config Key  = ecrit_reg
!Config Desc =
!Config Def  = 0.25  !4 fois par jour
!Config Help =
!
  ecrit_reg_omp = 0.25   !4 fois par jour
  call getin('ecrit_reg',ecrit_reg_omp)
!
!
!
! PARAMETRES CDRAG
!
!Config Key  = f_cdrag_ter
!Config Desc =
!Config Def  = 0.8
!Config Help =
!
  f_cdrag_ter_omp = 0.8
  call getin('f_cdrag_ter',f_cdrag_ter_omp)
!
!Config Key  = f_cdrag_oce
!Config Desc =
!Config Def  = 0.8
!Config Help =
!
  f_cdrag_oce_omp = 0.8
  call getin('f_cdrag_oce',f_cdrag_oce_omp)
!
! RUGORO
!Config Key  = f_rugoro
!Config Desc =
!Config Def  = 0.
!Config Help =
!
  f_rugoro_omp = 0.
  call getin('f_rugoro',f_rugoro_omp)

! PARAMETERS FOR CONVECTIVE INHIBITION BY TROPOS. DRYNESS
!
!Config Key  = supcrit1
!Config Desc =
!Config Def  = .540
!Config Help =
!
  supcrit1_omp = .540
  call getin('supcrit1',supcrit1_omp)

!
!Config Key  = supcrit2
!Config Desc =
!Config Def  = .600
!Config Help =
!
  supcrit2_omp = .600
  call getin('supcrit2',supcrit2_omp)

!
! PARAMETERS FOR THE MIXING DISTRIBUTION
! iflag_mix: 0=OLD, 
!            1=NEW (JYG),            
!            2=NEW + conv. depth inhib. by tropos. dryness
! '2' is NOT operationnal and should not be used.
!
!Config Key  = iflag_mix
!Config Desc =
!Config Def  = 1
!Config Help =
!
  iflag_mix_omp = 1
  call getin('iflag_mix',iflag_mix_omp)

!
!Config Key  = scut
!Config Desc =
!Config Def  = 0.95
!Config Help =
!
  scut_omp = 0.95
  call getin('scut',scut_omp)

!
!Config Key  = qqa1
!Config Desc =
!Config Def  = 1.0
!Config Help =
!
  qqa1_omp = 1.0
  call getin('qqa1',qqa1_omp)

!
!Config Key  = qqa2
!Config Desc =
!Config Def  = 0.0
!Config Help =
!
  qqa2_omp = 0.0
  call getin('qqa2',qqa2_omp)

!
!Config Key  = gammas
!Config Desc =
!Config Def  = 0.05
!Config Help =
!
  gammas_omp = 0.05
  call getin('gammas',gammas_omp)

!
!Config Key  = Fmax
!Config Desc =
!Config Def  = 0.65
!Config Help =
!
  Fmax_omp = 0.65
  call getin('Fmax',Fmax_omp)

!
!Config Key  = alphas  
!Config Desc =
!Config Def  = -5.
!Config Help =
!
  alphas_omp = -5.
  call getin('alphas',alphas_omp)

!Config key = ok_strato
!Config  Desc = activation de la version strato
!Config  Def  = .FALSE.
!Config  Help = active la version stratosphérique de LMDZ de F. Lott

  ok_strato_omp=.FALSE.
  CALL getin('ok_strato',ok_strato_omp)
      
!Config  key = ok_hines
!Config  Desc = activation de la parametrisation de hines
!Config  Def  = .FALSE.
!Config  Help = Clefs controlant la parametrization de Hines
!               Et la sponge layer (Runs Stratospheriques)

  ok_hines_omp=.FALSE.
  CALL getin('ok_hines',ok_hines_omp)

!Config Key  = OK_LES                                               
!Config Desc = Pour des sorties LES                                 
!Config Def  = .false.                                              
!Config Help = Pour creer le fichier histLES contenant les sorties  
!              LES                                                  
!                                                                   
  ok_LES_omp = .false.                                              
  call getin('OK_LES', ok_LES_omp)                                  

!Config Key  = callstats                                               
!Config Desc = Pour des sorties callstats                                 
!Config Def  = .false.                                              
!Config Help = Pour creer le fichier stats contenant les sorties  
!              stats                                                  
!                                                                   
  callstats_omp = .false.                                              
  call getin('callstats', callstats_omp)                                  
!
!Config Key  = ecrit_LES
!Config Desc = Frequence d'ecriture des resultats du LES en nombre de jours;
!              par defaut 1., i.e. 1 jour
!Config Def  = 1./8.
!Config Help = ... 
!
!
  ecrit_LES_omp = 1./8.
  call getin('ecrit_LES', ecrit_LES_omp)
!
  read_climoz = 0 ! default value
  call getin('read_climoz', read_climoz)

  carbon_cycle_tr_omp=.FALSE.
  CALL getin('carbon_cycle_tr',carbon_cycle_tr_omp)

  carbon_cycle_cpl_omp=.FALSE.
  CALL getin('carbon_cycle_cpl',carbon_cycle_cpl_omp)

!$OMP END MASTER
!$OMP BARRIER

    R_ecc = R_ecc_omp
    R_peri = R_peri_omp
    R_incl = R_incl_omp
    solaire = solaire_omp
    co2_ppm = co2_ppm_omp
    RCO2 = RCO2_omp
    CH4_ppb = CH4_ppb_omp
    RCH4 = RCH4_omp
    N2O_ppb = N2O_ppb_omp
    RN2O = RN2O_omp
    CFC11_ppt = CFC11_ppt_omp
    RCFC11 = RCFC11_omp
    CFC12_ppt = CFC12_ppt_omp
    RCFC12 = RCFC12_omp

    cycle_diurne = cycle_diurne_omp
    soil_model = soil_model_omp
    new_oliq = new_oliq_omp
    ok_orodr = ok_orodr_omp
    ok_orolf = ok_orolf_omp
    ok_limitvrai = ok_limitvrai_omp
    nbapp_rad = nbapp_rad_omp
    iflag_con = iflag_con_omp

    epmax = epmax_omp
    ok_adj_ema = ok_adj_ema_omp
    iflag_clw = iflag_clw_omp
    cld_lc_lsc = cld_lc_lsc_omp
    cld_lc_con = cld_lc_con_omp
    cld_tau_lsc = cld_tau_lsc_omp
    cld_tau_con = cld_tau_con_omp
    ffallv_lsc = ffallv_lsc_omp
    ffallv_con = ffallv_con_omp
    coef_eva = coef_eva_omp
    reevap_ice = reevap_ice_omp
    iflag_pdf = iflag_pdf_omp
    solarlong0 = solarlong0_omp
    qsol0 = qsol0_omp
    inertie_sol = inertie_sol_omp
    inertie_ice = inertie_ice_omp
    inertie_sno = inertie_sno_omp
    rad_froid = rad_froid_omp
    rad_chau1 = rad_chau1_omp
    rad_chau2 = rad_chau2_omp
    t_glace_min = t_glace_min_omp
    t_glace_max = t_glace_max_omp
    top_height = top_height_omp
    overlap = overlap_omp
    cdmmax = cdmmax_omp
    cdhmax = cdhmax_omp
    ksta = ksta_omp
    ksta_ter = ksta_ter_omp
    ok_kzmin = ok_kzmin_omp
    fmagic = fmagic_omp
    pmagic = pmagic_omp
    iflag_pbl = iflag_pbl_omp
    lev_histhf = lev_histhf_omp
    lev_histday = lev_histday_omp
    lev_histmth = lev_histmth_omp
    lev_histins = lev_histins_omp
    lev_histLES = lev_histLES_omp
    lev_histdayNMC = lev_histdayNMC_omp
    ok_histNMC(:) = ok_histNMC_omp(:)
    freq_outNMC(:) = freq_outNMC_omp(:)
    freq_calNMC(:) = freq_calNMC_omp(:)

    type_ocean = type_ocean_omp
    version_ocean = version_ocean_omp
    ok_veget = ok_veget_omp
    ok_newmicro = ok_newmicro_omp
    ok_journe = ok_journe_omp
    ok_hf = ok_hf_omp
    ok_mensuel = ok_mensuel_omp
    ok_instan = ok_instan_omp
    freq_ISCCP = freq_ISCCP_omp
    ecrit_ISCCP = ecrit_ISCCP_omp
    freq_COSP = freq_COSP_omp
    ok_ade = ok_ade_omp
    ok_aie = ok_aie_omp
    aerosol_couple = aerosol_couple_omp
    flag_aerosol=flag_aerosol_omp
    new_aod=new_aod_omp
    aer_type = aer_type_omp
    bl95_b0 = bl95_b0_omp
    bl95_b1 = bl95_b1_omp
    fact_cldcon = fact_cldcon_omp
    facttemps = facttemps_omp
    ratqsbas = ratqsbas_omp
    ratqshaut = ratqshaut_omp
    tau_ratqs = tau_ratqs_omp

    iflag_radia = iflag_radia_omp
    iflag_rrtm = iflag_rrtm_omp
    iflag_cldcon = iflag_cldcon_omp
    iflag_ratqs = iflag_ratqs_omp
    ip_ebil_phy = ip_ebil_phy_omp
    iflag_thermals = iflag_thermals_omp
    iflag_thermals_ed = iflag_thermals_ed_omp
    iflag_thermals_optflux = iflag_thermals_optflux_omp
    nsplit_thermals = nsplit_thermals_omp
    tau_thermals = tau_thermals_omp
    iflag_coupl = iflag_coupl_omp
    iflag_clos = iflag_clos_omp
    iflag_wake = iflag_wake_omp
    alp_offset = alp_offset_omp
    iflag_cvl_sigd = iflag_cvl_sigd_omp
    type_run = type_run_omp
    ok_isccp = ok_isccp_omp
    ok_cosp = ok_cosp_omp
    ok_mensuelCOSP = ok_mensuelCOSP_omp
    ok_journeCOSP = ok_journeCOSP_omp
    ok_hfCOSP = ok_hfCOSP_omp
    seuil_inversion=seuil_inversion_omp
    lonmin_ins = lonmin_ins_omp
    lonmax_ins = lonmax_ins_omp
    latmin_ins = latmin_ins_omp
    latmax_ins = latmax_ins_omp
    ecrit_hf   = ecrit_hf_omp
    ecrit_ins   = ecrit_ins_omp
    ecrit_day = ecrit_day_omp
    ecrit_mth = ecrit_mth_omp
    ecrit_tra = ecrit_tra_omp
    ecrit_reg = ecrit_reg_omp
    cvl_corr = cvl_corr_omp
    ok_lic_melt = ok_lic_melt_omp
    f_cdrag_ter=f_cdrag_ter_omp
    f_cdrag_oce=f_cdrag_oce_omp
    f_rugoro=f_rugoro_omp
    supcrit1 = supcrit1_omp
    supcrit2 = supcrit2_omp
    iflag_mix = iflag_mix_omp
    scut = scut_omp
    qqa1 = qqa1_omp
    qqa2 = qqa2_omp
    gammas = gammas_omp
    Fmax = Fmax_omp
    alphas = alphas_omp
    ok_strato = ok_strato_omp
    ok_hines = ok_hines_omp
    ok_LES = ok_LES_omp
    callstats = callstats_omp
    ecrit_LES = ecrit_LES_omp
    carbon_cycle_tr = carbon_cycle_tr_omp
    carbon_cycle_cpl = carbon_cycle_cpl_omp

! Test of coherence between type_ocean and version_ocean
    IF (type_ocean=='couple' .AND. (version_ocean/='opa8' .AND. version_ocean/='nemo') ) THEN
       WRITE(numout,*)' ERROR version_ocean=',version_ocean,' not valid in coupled configuration'
       CALL abort_gcm('conf_phys','version_ocean not valid',1)
    END IF

    IF (type_ocean=='slab' .AND. version_ocean=='xxxxxx') THEN
       version_ocean='sicOBS'
    ELSE IF (type_ocean=='slab' .AND. version_ocean/='sicOBS') THEN
       WRITE(numout,*)' ERROR version_ocean=',version_ocean,' not valid with slab ocean'
       CALL abort_gcm('conf_phys','version_ocean not valid',1)
    END IF

! Test sur new_aod. Ce flag permet de retrouver les resultats de l'AR4
! il n'est utilisable que lors du couplage avec le SO4 seul 
    IF (ok_ade .OR. ok_aie) THEN 
       IF ( .NOT. new_aod .AND.  flag_aerosol .NE. 1) THEN
          CALL abort_gcm('conf_phys','new_aod=.FALSE. not compatible avec flag_aerosol=1',1)
       END IF
    END IF

!$OMP MASTER

  write(numout,*)' ##############################################'
  write(numout,*)' Configuration des parametres de la physique: '
  write(numout,*)' Type ocean = ', type_ocean
  write(numout,*)' Version ocean = ', version_ocean
  write(numout,*)' Config veget = ', ok_veget
  write(numout,*)' Sortie journaliere = ', ok_journe
  write(numout,*)' Sortie haute frequence = ', ok_hf
  write(numout,*)' Sortie mensuelle = ', ok_mensuel
  write(numout,*)' Sortie instantanee = ', ok_instan
  write(numout,*)' Frequence appel simulateur ISCCP, freq_ISCCP =', freq_ISCCP
  write(numout,*)' Frequence appel simulateur ISCCP, ecrit_ISCCP =', ecrit_ISCCP
  write(numout,*)' Frequence appel simulateur COSP, freq_COSP =', freq_COSP
  write(numout,*)' Sortie bilan d''energie, ip_ebil_phy =', ip_ebil_phy
  write(numout,*)' Excentricite = ',R_ecc
  write(numout,*)' Equinoxe = ',R_peri
  write(numout,*)' Inclinaison =',R_incl
  write(numout,*)' Constante solaire =',solaire
  write(numout,*)' co2_ppm =',co2_ppm
  write(numout,*)' RCO2 = ',RCO2
  write(numout,*)' CH4_ppb =',CH4_ppb,' RCH4 = ',RCH4
  write(numout,*)' N2O_ppb =',N2O_ppb,' RN2O =  ',RN2O
  write(numout,*)' CFC11_ppt=',CFC11_ppt,' RCFC11 =  ',RCFC11
  write(numout,*)' CFC12_ppt=',CFC12_ppt,' RCFC12 =  ',RCFC12
  write(numout,*)' cvl_corr=', cvl_corr
  write(numout,*)'ok_lic_melt=', ok_lic_melt
  write(numout,*)'cycle_diurne=',cycle_diurne
  write(numout,*)'soil_model=',soil_model
  write(numout,*)'new_oliq=',new_oliq
  write(numout,*)'ok_orodr=',ok_orodr
  write(numout,*)'ok_orolf=',ok_orolf
  write(numout,*)'ok_limitvrai=',ok_limitvrai
  write(numout,*)'nbapp_rad=',nbapp_rad
  write(numout,*)'iflag_con=',iflag_con
  write(numout,*)' epmax = ', epmax
  write(numout,*)' ok_adj_ema = ', ok_adj_ema
  write(numout,*)' iflag_clw = ', iflag_clw
  write(numout,*)' cld_lc_lsc = ', cld_lc_lsc
  write(numout,*)' cld_lc_con = ', cld_lc_con
  write(numout,*)' cld_tau_lsc = ', cld_tau_lsc
  write(numout,*)' cld_tau_con = ', cld_tau_con
  write(numout,*)' ffallv_lsc = ', ffallv_lsc
  write(numout,*)' ffallv_con = ', ffallv_con
  write(numout,*)' coef_eva = ', coef_eva
  write(numout,*)' reevap_ice = ', reevap_ice
  write(numout,*)' iflag_pdf = ', iflag_pdf
  write(numout,*)' iflag_cldcon = ', iflag_cldcon
  write(numout,*)' iflag_radia = ', iflag_radia
  write(numout,*)' iflag_rrtm = ', iflag_rrtm
  write(numout,*)' iflag_ratqs = ', iflag_ratqs
  write(numout,*)' seuil_inversion = ', seuil_inversion
  write(numout,*)' fact_cldcon = ', fact_cldcon
  write(numout,*)' facttemps = ', facttemps
  write(numout,*)' ok_newmicro = ',ok_newmicro 
  write(numout,*)' ratqsbas = ',ratqsbas 
  write(numout,*)' ratqshaut = ',ratqshaut 
  write(numout,*)' tau_ratqs = ',tau_ratqs 
  write(numout,*)' top_height = ',top_height 
  write(numout,*)' rad_froid = ',rad_froid
  write(numout,*)' rad_chau1 = ',rad_chau1
  write(numout,*)' rad_chau2 = ',rad_chau2
  write(numout,*)' t_glace_min = ',t_glace_min
  write(numout,*)' t_glace_max = ',t_glace_max
  write(numout,*)' overlap = ',overlap 
  write(numout,*)' cdmmax = ',cdmmax 
  write(numout,*)' cdhmax = ',cdhmax 
  write(numout,*)' ksta = ',ksta 
  write(numout,*)' ksta_ter = ',ksta_ter 
  write(numout,*)' ok_kzmin = ',ok_kzmin 
  write(numout,*)' fmagic = ',fmagic
  write(numout,*)' pmagic = ',pmagic
  write(numout,*)' ok_ade = ',ok_ade
  write(numout,*)' ok_aie = ',ok_aie
  write(numout,*)' aerosol_couple = ', aerosol_couple
  write(numout,*)' flag_aerosol = ', flag_aerosol
  write(numout,*)' new_aod = ', new_aod
  write(numout,*)' aer_type = ',aer_type
  write(numout,*)' bl95_b0 = ',bl95_b0
  write(numout,*)' bl95_b1 = ',bl95_b1
  write(numout,*)' lev_histhf = ',lev_histhf 
  write(numout,*)' lev_histday = ',lev_histday 
  write(numout,*)' lev_histmth = ',lev_histmth 
  write(numout,*)' lev_histins = ',lev_histins
  write(numout,*)' lev_histLES = ',lev_histLES
  write(numout,*)' lev_histdayNMC = ',lev_histdayNMC
  write(numout,*)' ok_histNMC = ',ok_histNMC
  write(numout,*)' freq_outNMC = ',freq_outNMC
  write(numout,*)' freq_calNMC = ',freq_calNMC
  write(numout,*)' iflag_pbl = ', iflag_pbl
  write(numout,*)' iflag_thermals = ', iflag_thermals
  write(numout,*)' iflag_thermals_ed = ', iflag_thermals_ed
  write(numout,*)' iflag_thermals_optflux = ', iflag_thermals_optflux
  write(numout,*)' iflag_clos = ', iflag_clos
  write(numout,*)' type_run = ',type_run 
  write(numout,*)' ok_isccp = ',ok_isccp 
  write(numout,*)' ok_cosp = ',ok_cosp
  write(numout,*)' ok_mensuelCOSP = ',ok_mensuelCOSP
  write(numout,*)' ok_journeCOSP = ',ok_journeCOSP
  write(numout,*)' ok_hfCOSP =',ok_hfCOSP
  write(numout,*)' solarlong0 = ', solarlong0
  write(numout,*)' qsol0 = ', qsol0
  write(numout,*)' inertie_sol = ', inertie_sol
  write(numout,*)' inertie_ice = ', inertie_ice
  write(numout,*)' inertie_sno = ', inertie_sno
  write(numout,*)' f_cdrag_ter = ',f_cdrag_ter
  write(numout,*)' f_cdrag_oce = ',f_cdrag_oce
  write(numout,*)' f_rugoro = ',f_rugoro
  write(numout,*)' supcrit1 = ', supcrit1
  write(numout,*)' supcrit2 = ', supcrit2
  write(numout,*)' iflag_mix = ', iflag_mix
  write(numout,*)' scut = ', scut
  write(numout,*)' qqa1 = ', qqa1
  write(numout,*)' qqa2 = ', qqa2
  write(numout,*)' gammas = ', gammas
  write(numout,*)' Fmax = ', Fmax
  write(numout,*)' alphas = ', alphas
  write(numout,*)' iflag_wake = ', iflag_wake
  write(numout,*)' alp_offset = ', alp_offset

  write(numout,*)' lonmin lonmax latmin latmax bilKP_ins =',&
   lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
  write(numout,*)' ecrit_ hf, ins, day, mth, reg, tra, ISCCP, LES',&
   ecrit_hf, ecrit_ins, ecrit_day, ecrit_mth, ecrit_reg, ecrit_tra, ecrit_ISCCP, ecrit_LES

  write(numout,*) 'ok_strato = ', ok_strato
  write(numout,*) 'ok_hines = ',  ok_hines
  write(numout,*) 'read_climoz = ', read_climoz
  write(numout,*) 'carbon_cycle_tr = ', carbon_cycle_tr
  write(numout,*) 'carbon_cycle_cpl = ', carbon_cycle_cpl
  
!$OMP END MASTER

  return
  
  end subroutine conf_phys

end module conf_phys_m
!
!#################################################################
!

   subroutine conf_interface(tau_calv)

   use IOIPSL
   implicit none

! Configuration de l'interace atm/surf
!
! tau_calv:    temps de relaxation pour la fonte des glaciers

  REAL          :: tau_calv
  REAL,SAVE     :: tau_calv_omp

! Local
  integer              :: numout = 6
!
!Config Key  = tau_calv
!Config Desc = temps de relaxation pour fonte des glaciers en jours
!Config Def  = 1 an 
!Config Help = 
!
  tau_calv_omp = 360.*10.
!$OMP MASTER
  call getin('tau_calv',tau_calv_omp)
!$OMP END MASTER
!$OMP BARRIER

  tau_calv=tau_calv_omp
  
!$OMP MASTER
  write(numout,*)' ##############################################'
  WRITE(numout,*)' Configuration de l''interface atm/surfaces  : '
  WRITE(numout,*)' tau_calv = ',tau_calv
!$OMP END MASTER

  return

  end subroutine conf_interface
