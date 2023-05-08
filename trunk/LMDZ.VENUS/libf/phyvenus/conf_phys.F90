!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/conf_phys.F90,v 1.3 2005/02/07 15:15:31 fairhead Exp $
!
!
!

  subroutine conf_phys(ok_journe, ok_mensuel, ok_instan, &
 &                     if_ebil)
   use init_print_control_mod, only: init_print_control
   use print_control_mod, only: lunout
   use IOIPSL, only: getin

   implicit none

   include "YOMCST.h"
   include "clesphys.h"
   include "compbl.h"

! ok_journe:  sorties journalieres
! ok_mensuel: sorties mensuelles
! ok_instan:  sorties instantanees


! Sortie:
  logical,intent(out)  :: ok_journe, ok_mensuel, ok_instan        
  integer,intent(out)  :: if_ebil



!
! Configuration de la "physique" de LMDZ a l'aide de la fonction
! GETIN de IOIPSL
!
! LF 05/2001
!
!--- Ca lit le physiq.def ---

       ! initialize print_control module variables
       call init_print_control 

!******************* parametres anciennement lus dans gcm.def

       ! do we read a startphy.nc file? (default: .true.)
       startphy_file=.true.
       CALL getin("startphy_file",startphy_file)

!Config  Key  = cycle_diurne
!Config  Desc = Cycle diurne
!Config  Def  = y
!Config  Help = Cette option permet d'eteidre le cycle diurne.
!Config         Peut etre util pour accelerer le code !
       cycle_diurne = .TRUE.       
       call getin('cycle_diurne',cycle_diurne)

!Config  Key  = soil_model
!Config  Desc = Modele de sol
!Config  Def  = y
!Config  Help = Choix du modele de sol (Thermique ?)
!Config         Option qui pourait un string afin de pouvoir
!Config         plus de choix ! Ou meme une liste d'options !
       soil_model = .true.
       call getin('soil_model',soil_model)

!Config  Key  = ok_orodr
!Config  Desc = Oro drag
!Config  Def  = y
!Config  Help = GW drag orographie
!Config         
       ok_orodr = .false.
       call getin('ok_orodr',ok_orodr)

!Config  Key  =  ok_orolf
!Config  Desc = Oro lift
!Config  Def  = n
!Config  Help = GW lift orographie (pas utilise)
       ok_orolf = .false.
       call getin('ok_orolf', ok_orolf)

!Config  Key  = ok_gw_nonoro
!Config  Desc = Gravity waves parameterization
!Config  Def  = n
!Config  Help = GW drag non-orographique
       ok_gw_nonoro = .false.
       call getin('ok_gw_nonoro',ok_gw_nonoro)

!Config  Key  = nbapp_rad
!Config  Desc = Frequence d'appel au rayonnement
!Config  Def  = 12
!Config  Help = Nombre  d'appels des routines de rayonnements
!Config         par jour.
       nbapp_rad = 12
       call getin('nbapp_rad',nbapp_rad)
       print*,"nbapp_rad",nbapp_rad
!Config  Key  = nbapp_chim
!Config  Desc = Frequence d'appel a la chimie
!Config  Def  = 1
!Config  Help = Nombre  d'appels des routines de chimie
!Config         par jour.
       nbapp_chim = 1
       call getin('nbapp_chim',nbapp_chim)

!Config  Key  = iflag_con
!Config  Desc = Flag de convection
!Config  Def  = 0
!Config  Help = Flag  pour la convection les options suivantes existent :
!Config         0 : ajsec simple (VENUS, TITAN)
!Config         1 pour LMD,
!Config         2 pour Tiedtke,
!Config         3 pour CCM(NCAR)  
       iflag_con = 0
       call getin('iflag_con',iflag_con)

!******************* fin parametres anciennement lus dans gcm.def

!Config Key  = OK_journe
!Config Desc = Pour des sorties journalieres 
!Config Def  = .false.
!Config Help = Pour creer le fichier histday contenant les sorties
!              journalieres 
!
  ok_journe = .false.
  call getin('OK_journe', ok_journe)
!
!Config Key  = OK_mensuel
!Config Desc = Pour des sorties mensuelles 
!Config Def  = .false.
!Config Help = Pour creer le fichier histmth contenant les sorties
!              mensuelles 
!
  ok_mensuel = .false.
  call getin('OK_mensuel', ok_mensuel)
!
!Config Key  = OK_instan
!Config Desc = Pour des sorties instantanees 
!Config Def  = .false.
!Config Help = Pour creer le fichier histins contenant les sorties
!              instantanees 
!
  ok_instan = .false.
  call getin('OK_instan', ok_instan)
!
!Config  Key  = ecritphy
!Config  Desc = Frequence d'ecriture dans histmth et histins
!Config  Def  = 1
!Config  Help = frequence de l'ecriture du fichier histmth et histins
!Config         en jours.
!
       ecriphy = 1.
       call getin('ecritphy', ecriphy)
!
!
!Config Key  = if_ebil
!Config Desc = Niveau de sortie pour les diags bilan d'energie 
!Config Def  = 0
!Config Help = 
!               
!
  if_ebil = 0
  call getin('if_ebil', if_ebil)
!!
!! Constante solaire & Parametres orbitaux & taux gaz effet de serre BEG
!!
!Config Key  = R_ecc
!Config Desc = Excentricite
!Config Def  = 0.006787
!Config Help = 
!               
! VENUS
! R_ecc = 0.006787
  R_ecc   = 0.0
  call getin('R_ecc', R_ecc)
!!
!Config Key  = R_peri
!Config Desc = Equinoxe
!Config Def  = 
!Config Help = 
!               
! VENUS
  R_peri = 0.
  call getin('R_peri', R_peri)
!!
!Config Key  = R_incl
!Config Desc = Inclinaison
!Config Def  = 
!Config Help = 
!               
! VENUS
  R_incl = 0.0
  call getin('R_incl', R_incl)
!
!Config Key  = solaire
!Config Desc = Constante solaire en W/m2
! VENUS
!Config Def  = 2620.
!Config Help = 
!
  solaire = 2620.
    call getin('solaire', solaire)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR THE PLANETARY BOUNDARY LAYER AND SOIL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Config Key  = iflag_pbl
!Config Desc =
!Config Def  = 1
!Config Help =
!
! 2   = calculs Cd et K simples pour VENUS :
!       parametres = z0, lmixmin, ksta (en dur: umin2,ric,cepdu2,karman)
! 1   = calculs Cd et K issus LMDZ Terre
!       parametres = ksta, ok_kzmin (et plein d'autres en dur...)
! 6-9 = schema des thermiques Fred
  iflag_pbl = 1
  call getin('iflag_pbl',iflag_pbl)

!
!Config Key  = ksta
!Config Desc =
!Config Def  = 1.0e-7
!Config Help =
!
  ksta = 1.0e-7
  call getin('ksta',ksta)

!
!Config Key  = z0
!Config Desc =
!Config Def  = 1.0e-2
!Config Help =
!
  z0 = 1.0e-2
  call getin('z0',z0)

!
!Config Key  = lmixmin
!Config Desc =
!Config Def  = 35.
!Config Help =
!
  lmixmin = 35.
  call getin('lmixmin',lmixmin)

!
!Config Key  = ok_kzmin
!Config Desc =
!Config Def  = .false.
!Config Help =
!
  ok_kzmin = .false.
  call getin('ok_kzmin',ok_kzmin)

  ok_clmain = .true.
  call getin('ok_clmain',ok_clmain)

  physideal = .false.
  call getin('physideal',physideal)

!Config Key  = iflag_ajs
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_ajs = 0
  call getin('iflag_ajs',iflag_ajs)

!
!Config Key  = inertie
!Config Desc =
!Config Def  = 2000.
!Config Help =
!
  inertie = 2000.
  call getin('inertie',inertie)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR THE OUTPUT LEVELS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Config Key  = lev_histins
!Config Desc =
!Config Def  = 0
!Config Help =
!
  lev_histins = 0
  call getin('lev_histins',lev_histins)

!
!Config Key  = lev_histday
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histday = 1
  call getin('lev_histday',lev_histday)

!
!Config Key  = lev_histmth
!Config Desc =
!Config Def  = 2
!Config Help =
!
  lev_histmth = 2
  call getin('lev_histmth',lev_histmth)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR THE TRACERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Config Key  = tr_scheme
!Config Desc =
!Config Def  = 0
!Config Help =
!
! 0   = Nothing is done (passive tracers)
! 1   = pseudo-chemistry with relaxation toward fixed profile
!       See Marcq&Lebonnois 2013
! 2   = surface emission
!       For the moment, inspired from Mars version
!       However, the variable 'source' could be used in physiq
!       so the call to phytrac_emiss could be to initialise it.
! 3   = Full chemistry and/or clouds => phytrac_chimie
!       Need ok_chem or ok_cloud
  tr_scheme = 0
  call getin('tr_scheme',tr_scheme)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   PARAMETRES DE LA CHIMIE/NUAGE dans physiq.def
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!Config Key  = reinit_trac
!Config Desc =  
!Config Def  = .FALSE.
!Config Help = 
!
  reinit_trac = .FALSE.
  call getin('reinit_trac',reinit_trac)
  
!
!Config Key  = ok_cloud
!Config Desc =  
!Config Def  = .FALSE.
!Config Help = 
!
  ok_cloud = .false.
  call getin('ok_cloud',ok_cloud)

!
!Config Key  = cl_scheme
!Config Desc =
!Config Def  = 2
!Config Help =
!
! 1   = Simple microphysics (Aurelien Stolzenbach's PhD)
! 2   = Full microphysics (momentum scheme, Sabrina Guilbon's PhD)

  cl_scheme = 2
  call getin('cl_scheme',cl_scheme)

!
!Config Key  = ok_chem
!Config Desc =  
!Config Def  = .FALSE.
!Config Help = 
!
  ok_chem = .false.
  call getin('ok_chem',ok_chem)

  if (((tr_scheme.ne.3).and.(ok_chem.or.ok_cloud)).or. &
      ((tr_scheme.eq.3).and.(.not.ok_chem.and..not.ok_cloud))) then
    write(*,*) "Attention, incoherence :" 
    write(*,*) "tr_scheme=",tr_scheme," / ok_chem=",ok_chem, &
                                     " / ok_cloud=",ok_cloud
    write(*,*) "Verifier votre physiq.def"
    stop
  endif

!
!Config Key  = ok_sedim
!Config Desc =  
!Config Def  = .FALSE.
!Config Help = 
!
  ok_sedim = .false.
  call getin('ok_sedim',ok_sedim)

!
!Config Key  = nb_mode
!Config Desc =  
!Config Def  = 0
!Config Help = 
!
  nb_mode = 0
  call getin('nb_mode',nb_mode)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR NLTE PHYSICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!Config Key  = callnlte
!Config Desc = 
!Config Def  = .false.
!Config Help = 
!
  callnlte = .false.
  call getin('callnlte',callnlte)

!
!Config Key  = callnirco2
!Config Desc = 
!Config Def  = .false.
!
  callnirco2 = .false.
  call getin('callnirco2',callnirco2)

!
!Config Key  = nircorr
!Config Desc = 
!Config Def  = 0
!Config Help =
!
  nircorr = 0
  call getin('nircorr',nircorr)

!
!Config Key  = callthermos
!Config Desc = 
!Config Def  = .false.
!Config Help =
!
  callthermos = .false.
  call getin('callthermos',callthermos)

!
!Config Key  = nltemodel
!Config Desc = 
!Config Def  = 0
!Config Help =
!
  nltemodel = 0
  call getin('nltemodel',nltemodel)

!
!Config Key  = solvarmod
!Config Desc =
!Config Def  = 1
!Config Help =
!
  solvarmod = 1
  call getin('solvarmod',solvarmod)

!
!Config Key  = solarcondate
!Config Desc = 
!Config Def  = 1993.4  ## Average solar cycle condition
!Config Help =
!
  solarcondate = 1993.4
  call getin('solarcondate',solarcondate)

!
!Config Key  = euveff
!Config Desc =
!Config Def  = 0.21
!Config Help =
!
  euveff = 0.21
  call getin('euveff',euveff)

!
!
!Config Key  = 
!Config Desc =  
!Config Def  =
!Config Help = 
!
!   =
!  call getin('',)
!
!
!
!

  write(lunout,*)' ##############################################'
  write(lunout,*)' Configuration des parametres de la physique: '
  write(lunout,*)' cycle_diurne = ', cycle_diurne
  write(lunout,*)' soil_model = ', soil_model
  write(lunout,*)' ok_orodr = ', ok_orodr
  write(lunout,*)' ok_orolf = ', ok_orolf
  write(lunout,*)' ok_gw_nonoro = ', ok_gw_nonoro
  write(lunout,*)' nbapp_rad = ', nbapp_rad
  write(lunout,*)' nbapp_chim = ', nbapp_chim
  write(lunout,*)' iflag_con = ', iflag_con
  write(lunout,*)' Sortie journaliere = ', ok_journe
  write(lunout,*)' Sortie mensuelle = ', ok_mensuel
  write(lunout,*)' Sortie instantanee = ', ok_instan
  write(lunout,*)' frequence sorties = ', ecriphy  
  write(lunout,*)' Sortie bilan d''energie, if_ebil =', if_ebil
  write(lunout,*)' Excentricite = ',R_ecc
  write(lunout,*)' Equinoxe = ',R_peri
  write(lunout,*)' Inclinaison =',R_incl
  write(lunout,*)' tr_scheme = ', tr_scheme
  write(lunout,*)' iflag_pbl = ', iflag_pbl
  write(lunout,*)' z0 = ',z0 
  write(lunout,*)' lmixmin = ',lmixmin 
  write(lunout,*)' ksta = ',ksta 
  write(lunout,*)' ok_kzmin = ',ok_kzmin 
  write(lunout,*)' inertie = ', inertie 
  write(lunout,*)' ok_clmain = ',ok_clmain
  write(lunout,*)' physideal = ',physideal
  write(lunout,*)' iflag_ajs = ', iflag_ajs
  write(lunout,*)' lev_histins = ',lev_histins 
  write(lunout,*)' lev_histday = ',lev_histday 
  write(lunout,*)' lev_histmth = ',lev_histmth 
  write(lunout,*)' reinit_trac = ',reinit_trac
  write(lunout,*)' ok_cloud = ',ok_cloud
  write(lunout,*)' cl_scheme = ',cl_scheme
  write(lunout,*)' ok_chem = ',ok_chem
  write(lunout,*)' ok_sedim = ',ok_sedim
  write(lunout,*)' nb_mode = ',nb_mode
  write(lunout,*)' callnlte = ',callnlte
  write(lunout,*)' nltemodel = ',nltemodel
  write(lunout,*)' callnirco2 = ',callnirco2
  write(lunout,*)' nircorr = ',nircorr
  write(lunout,*)' callthermos = ',callthermos
  write(lunout,*)' solvarmod = ',solvarmod
  write(lunout,*)' solarcondate = ',solarcondate
  write(lunout,*)' euveff = ',euveff

  end subroutine conf_phys

