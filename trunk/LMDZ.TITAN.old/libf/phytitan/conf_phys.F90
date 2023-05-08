!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/conf_phys.F90,v 1.3 2005/02/07 15:15:31 fairhead Exp $
!
!
!

  subroutine conf_phys(ok_mensuel,ok_journe,ok_instan,if_ebil)

   use IOIPSL
   implicit none

#include "YOMCST.h"
#include "clesphys.h"
#include "compbl.h"
#include "comorbit.h"

! ok_journe:  sorties journalieres
! ok_mensuel: sorties mensuelles
! ok_instan:  sorties instantanees


! Sortie:
  logical              :: ok_journe, ok_mensuel, ok_instan        
  integer              :: if_ebil

! Local
  integer              :: numout = 6

!
! Configuration de la "physique" de LMDZ a l'aide de la fonction
! GETIN de IOIPSL
!
! LF 05/2001
!
!--- Ca lit le physiq.def ---

!******************* parametres anciennement lus dans gcm.def

!Config  Key  = cycle_diurne
!Config  Desc = Cycle ddiurne
!Config  Def  = y
!Config  Help = Cette option permet d'eteidre le cycle diurne.
!Config         Peut etre util pour accelerer le code !
       cycle_diurne = .TRUE.
       CALL getin('cycle_diurne',cycle_diurne)

!Config  Key  = soil_model
!Config  Desc = Modele de sol
!Config  Def  = y
!Config  Help = Choix du modele de sol (Thermique ?)
!Config         Option qui pourait un string afin de pouvoir
!Config         plus de choix ! Ou meme une liste d'options !
       soil_model = .TRUE.
       CALL getin('soil_model',soil_model)

!Config  Key  = ok_orodr
!Config  Desc = Oro drag
!Config  Def  = y
!Config  Help = GW drag orographie
!Config         
       ok_orodr = .TRUE.
       CALL getin('ok_orodr',ok_orodr)

!Config  Key  =  ok_orolf
!Config  Desc = Oro lift
!Config  Def  = n
!Config  Help = GW lift orographie (pas utilise)
       ok_orolf = .TRUE.
       CALL getin('ok_orolf', ok_orolf)

!Config  Key  = ok_gw_nonoro
!Config  Desc = Gravity waves parameterization
!Config  Def  = n
!Config  Help = GW drag non-orographique
       ok_gw_nonoro = .FALSE.
       CALL getin('ok_gw_nonoro',ok_gw_nonoro)

!Config  Key  = nbapp_rad
!Config  Desc = Frequence d'appel au rayonnement
!Config  Def  = 12
!Config  Help = Nombre  d'appels des routines de rayonnements
!Config         par jour.
       nbapp_rad = 12
       CALL getin('nbapp_rad',nbapp_rad)

!Config  Key  = nbapp_chim
!Config  Desc = Frequence d'appel a la chimie
!Config  Def  = 1
!Config  Help = Nombre  d'appels des routines de chimie
!Config         par jour.
       nbapp_chim = 1
       CALL getin('nbapp_chim',nbapp_chim)

!Config  Key  = iflag_con
!Config  Desc = Flag de convection
!Config  Def  = 0
!Config  Help = Flag  pour la convection les options suivantes existent :
!Config         0 : ajsec simple (VENUS, TITAN)
!Config         1 pour LMD,
!Config         2 pour Tiedtke,
!Config         3 pour CCM(NCAR)  
       iflag_con = 0
       CALL getin('iflag_con',iflag_con)

!******************* fin parametres anciennement lus dans gcm.def

!Config Key  = OK_mensuel
!Config Desc = Pour des sorties mensuelles 
!Config Def  = .true.
!Config Help = Pour creer le fichier histmth contenant les sorties
!              mensuelles 
!
  ok_mensuel = .true.
  call getin('OK_mensuel', ok_mensuel)
!
!Config Key  = OK_journe
!Config Desc = Pour des sorties journalieres 
!Config Def  = .false.
!Config Help = Pour creer le fichier histday contenant les sorties
!              journalieres 
!
  ok_journe = .false.
  call getin('OK_journe', ok_journe)
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
!
!Config  Key  = ecritphy
!Config  Desc = Frequence d'ecriture dans histins
!Config  Def  = 1
!Config  Help = frequence de l'ecriture du fichier histins
!Config         en jours.
!
       ecriphy = 1.
       CALL getin('ecritphy', ecriphy)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Constante solaire & Parametres orbitaux 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
! TITAN         ! Valeurs par defaut d'apres Fig Tokano.
!
!!
!Config Key  = year_day
!Config Desc = Duree de l'annee en jour
!Config Def  = 
!Config Help = 
!               
  year_day = 673.
  call getin('year_day', year_day)
!
!Config Key  = peri_day
!Config Desc = position du perihelie en jour
!Config Def  = 
!Config Help = 
!               
  peri_day = 533.
  call getin('peri_day', peri_day)
!
!Config Key  = periheli
!Config Desc = Distance au soleil au perihelie
!Config Def  = 
!Config Help = 
!               
  periheli = 1354.5
  call getin('periheli', periheli)
!!
!Config Key  = aphelie
!Config Desc = Distance au soleil a l'aphelie
!Config Def  = 
!Config Help = 
!               
  aphelie = 1506.0
  call getin('aphelie', aphelie)
!!
!Config Key  = obliquit
!Config Desc = Obliquite
!Config Def  = 
!Config Help = 
!               
  obliquit = 26.7
  call getin('obliquit', obliquit)
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


!Config Key  = iflag_ajs
!Config Desc =
!Config Def  = 0
!Config Help =
!
  iflag_ajs = 1
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
!Config Key  = emis
!Config Desc =
!Config Def  = 0.95
!Config Help =
!
  emis = 0.95
  call getin('emis',emis)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parametres CHIMIE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Config Key  = chimi
!Config Desc =
!Config Def  = .false.
!Config Help =
!
  chimi = .false.
  call getin('chimi',chimi)

!
!Config Key  = vchim
!Config Desc =
!Config Def  = 1
!Config Help =
!
  vchim = 1
  call getin('vchim',vchim)

!
!Config Key  = aerprod
!Config Desc =
!Config Def  = 0
!Config Help =
!
  aerprod = 0
  call getin('aerprod',aerprod)

!
!Config Key  = htoh2
!Config Desc =
!Config Def  = 1
!Config Help =
!
  htoh2 = 1
  call getin('htoh2',htoh2)

!
!Config Key  = ylellouch
!Config Desc =
!Config Def  = .true.
!Config Help =
!
  ylellouch = .true.
  call getin('ylellouch',ylellouch)

!
!Config Key  = hcnrad
!Config Desc =
!Config Def  = .false.
!Config Help =
!
  hcnrad = .false.
  call getin('hcnrad',hcnrad)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parametres MICROPHYSIQUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Config Key  = microfi
!Config Desc =
!Config Def  = 1
!Config Help =
!
  microfi = 1
  call getin('microfi',microfi)

!
!Config Key  = tx
!Config Desc =
!Config Def  = 3.5
!Config Help =
!
  tx = 3.5
  call getin('tx',tx)

!
!Config Key  = tcorrect
!Config Desc =
!Config Def  = 1.
!Config Help =
!
  tcorrect = 1.
  call getin('tcorrect',tcorrect)

!
!Config Key  = xvis
!Config Desc = Facteur d ajustement des proprietes vis des aerosols
!Config Def  = 1.5
!Config Help =
!
  xvis = 1.0
  call getin('xvis',xvis)
!
!Config Key  = xir
!Config Desc = Facteur d ajustement des proprietes IR des aerosols
!Config Def  = 0.5
!Config Help =
!
  xir = 1.0
  call getin('xir',xir)
!
!Config Key  = p_prodaer
!Config Desc = pressure level for aerosol production (in Pa)
!Config Def  = 1.
!Config Help =
!
  p_prodaer = 1.
  call getin('p_prodaer',p_prodaer)
!
!Config Key  = cutoff
!Config Desc =
!Config Def  = 2
!Config Help =
!
  cutoff = 2
  call getin('cutoff',cutoff)

!
!Config Key  = clouds
!Config Desc = activation des nuages
!Config Def  = 1
!Config Help =
!
  clouds = 1
  call getin('clouds',clouds)
  if (microfi.lt.1) clouds = 0      ! On  ne fait pas de nuages sans microphysique !
  if (clouds.eq.1) cutoff = 0       ! si nuages, il faut mettre ca

!
!Config Key  = xnuf
!Config Desc = fraction nuageuse
!Config Def  = 0.5
!Config Help =
!
  xnuf = 0.5
  call getin('xnuf',xnuf)
  xnuf = amax1(xnuf,0.1)          ! On garde au minimum 10% de nuages.
  if (clouds.eq.0) xnuf = 0.      ! Si il n'y pas de nuages, on ne met pas de fraction
                                  ! nuageuse -> permet de retomber sur le TR habituel.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARAMETER FOR THE OUTPUT LEVELS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!Config Key  = lev_histmth
!Config Desc =
!Config Def  = 2
!Config Help =
! 
  lev_histmth = 2
  call getin('lev_histmth',lev_histmth)

!
!Config Key  = lev_histday
!Config Desc =
!Config Def  = 1
!Config Help =
!
  lev_histday = 1
  call getin('lev_histday',lev_histday)

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

  write(numout,*)' ##############################################'
  write(numout,*)' Configuration des parametres de la physique: '
  write(numout,*)' cycle_diurne = ', cycle_diurne
  write(numout,*)' soil_model = ', soil_model
  write(numout,*)' ok_orodr = ', ok_orodr
  write(numout,*)' ok_orolf = ', ok_orolf
  write(numout,*)' ok_gw_nonoro = ', ok_gw_nonoro
  write(numout,*)' nbapp_rad = ', nbapp_rad
  write(numout,*)' nbapp_chim = ', nbapp_chim
  write(numout,*)' iflag_con = ', iflag_con
  write(numout,*)' Sortie mensuelle = ', ok_mensuel
  write(numout,*)' Sortie journaliere = ', ok_journe
  write(numout,*)' Sortie instantanee = ', ok_instan
  write(numout,*)' frequence sorties = ', ecriphy  
  write(numout,*)' Sortie bilan d energie, if_ebil =', if_ebil
  write(numout,*)' Duree de l annee = ',year_day
  write(numout,*)' Position du perihelie = ',peri_day
  write(numout,*)' Perihelie = ',periheli
  write(numout,*)' Aphelie = ',aphelie
  write(numout,*)' Obliquite =',obliquit
  write(numout,*)' iflag_pbl = ', iflag_pbl
  write(numout,*)' z0 = ',z0 
  write(numout,*)' lmixmin = ',lmixmin 
  write(numout,*)' ksta = ',ksta 
  write(numout,*)' ok_kzmin = ',ok_kzmin 
  write(numout,*)' inertie = ', inertie 
  write(numout,*)' emis = ', emis 
  write(numout,*)' iflag_ajs = ', iflag_ajs
  write(numout,*)' chimi = ', chimi
  write(numout,*)' vchim = ', vchim
  write(numout,*)' aerprod = ', aerprod
  write(numout,*)' htoh2 = ', htoh2
  write(numout,*)' ylellouch = ', ylellouch
  write(numout,*)' hcnrad = ', hcnrad
  write(numout,*)' microfi = ', microfi
  write(numout,*)' tx = ', tx
  write(numout,*)' tcorrect = ', tcorrect
  write(numout,*)' xvis = ', xvis
  write(numout,*)' xir = ', xir
  write(numout,*)' p_prodaer = ', p_prodaer
  write(numout,*)' cutoff = ', cutoff
  write(numout,*)' clouds = ', clouds
  write(numout,*)' xnuf = ', xnuf
  write(numout,*)' lev_histmth = ',lev_histmth
  write(numout,*)' lev_histday = ',lev_histday 

  return

  end subroutine conf_phys

