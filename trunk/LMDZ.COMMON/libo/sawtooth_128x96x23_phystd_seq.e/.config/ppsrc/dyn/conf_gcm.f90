










!
! $Id: conf_gcm.F 1418 2010-07-19 15:11:24Z jghattas $

SUBROUTINE conf_gcm( tapedef, etatinit )

  ! if not using IOIPSL, we still need to use (a local version of) getin
  use ioipsl_getincom
  USE control_mod, ONLY: anneeref, config_inca, day_step, dayref, &
                         dissip_period, fractday, iapp_tracvl, &
                         iconser, iecri, ip_ebil_dyn, iperiod, &
                         iphysiq, less1day, nday, ndynstep, nsplit_phys, &
                         offline, ok_dyn_ave, ok_dyn_ins, ok_dynzon, &
                         output_grads_dyn, periodav, planet_type, &
                         raz_date, resetvarc, starttime, timestart, &
                         ecritstart, cpofT,force_conserv_tracer
  USE infotrac, ONLY : type_trac
  use assert_m, only: assert
  use sponge_mod, only: callsponge,mode_sponge,nsponge,tetasponge
  USE comconst_mod, ONLY: dissip_factz,dissip_deltaz,dissip_zref, 		&
		dissip_fac_mid,dissip_fac_up,dissip_hdelta,dissip_pupstart, 	&
		mode_top_bound,tau_top_bound,iflag_top_bound,ngroup
  USE logic_mod, ONLY: tidal,purmats,ok_guide,read_start,iflag_phys,		&
		iflag_trac,ok_strato,ok_gradsfile,ok_limit,ok_etat0,		&
		moyzon_mu,moyzon_ch,ok_strato,fxyhypb,ysinus,read_orop
  USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy,		&
		alphax,alphay,taux,tauy
  USE temps_mod, ONLY: calend

  IMPLICIT NONE
!-----------------------------------------------------------------------
!     Auteurs :   L. Fairhead , P. Le Van  .
!
!     Arguments :
!
!     tapedef   :
!     etatinit  :     = TRUE   , on ne  compare pas les valeurs des para- 
!     -metres  du zoom  avec  celles lues sur le fichier start .
!
  LOGICAL,INTENT(IN) :: etatinit
  INTEGER,INTENT(IN) :: tapedef

!   Declarations :
!   --------------
  include "dimensions.h"
  include "paramet.h"
  include "comdissnew.h"
  include "iniprint.h"

!   local:
!   ------

  REAL clonn,clatt,grossismxx,grossismyy
  REAL dzoomxx,dzoomyy, tauxx,tauyy
  LOGICAL  fxyhypbb, ysinuss
  LOGICAL use_filtre_fft
!
!  -------------------------------------------------------------------
!
!       .........     Version  du 29/04/97       ..........
!
!   Nouveaux parametres nitergdiv,nitergrot,niterh,tetagdiv,tetagrot,
!      tetatemp   ajoutes  pour la dissipation   .
!
!   Autre parametre ajoute en fin de liste de tapedef : ** fxyhypb ** 
!
!  Si fxyhypb = .TRUE. , choix de la fonction a derivee tangente hyperb.
!    Sinon , choix de fxynew  , a derivee sinusoidale  ..
!
!   ......  etatinit = . TRUE. si defrun  est appele dans ETAT0_LMD  ou
!         LIMIT_LMD  pour l'initialisation de start.dat (dic) et
!                de limit.dat ( dic)                        ...........
!           Sinon  etatinit = . FALSE .
!
!   Donc etatinit = .F.  si on veut comparer les valeurs de  grossismx ,
!    grossismy,clon,clat, fxyhypb  lues sur  le fichier  start  avec
!   celles passees  par run.def ,  au debut du gcm, apres l'appel a 
!    lectba .  
!   Ces parmetres definissant entre autres la grille et doivent etre
!   pareils et coherents , sinon il y aura  divergence du gcm .
!
!-----------------------------------------------------------------------
!   initialisations:
!   ----------------

!Config  Key  = lunout
!Config  Desc = unite de fichier pour les impressions
!Config  Def  = 6
!Config  Help = unite de fichier pour les impressions 
!Config         (defaut sortie standard = 6)
  lunout=6
  CALL getin('lunout', lunout)
  IF (lunout /= 5 .and. lunout /= 6) THEN
    OPEN(UNIT=lunout,FILE='lmdz.out',ACTION='write',                      &
          STATUS='unknown',FORM='formatted')
  ENDIF

!Config  Key  = prt_level
!Config  Desc = niveau d'impressions de débogage
!Config  Def  = 0
!Config  Help = Niveau d'impression pour le débogage
!Config         (0 = minimum d'impression)
  prt_level = 0
  CALL getin('prt_level',prt_level)

!-----------------------------------------------------------------------
!  Parametres de controle du run:
!-----------------------------------------------------------------------
!Config  Key  = planet_type
!Config  Desc = planet type ("earth", "mars", "venus", ...)
!Config  Def  = earth
!Config  Help = this flag sets the type of atymosphere that is considered
  planet_type="earth"
  CALL getin('planet_type',planet_type)

!Config  Key  = calend
!Config  Desc = type de calendrier utilise
!Config  Def  = earth_360d
!Config  Help = valeur possible: earth_360d, earth_365d, earth_366d
!Config         
  calend = 'earth_360d'
  CALL getin('calend', calend)

!Config  Key  = dayref
!Config  Desc = Jour de l'etat initial
!Config  Def  = 1
!Config  Help = Jour de l'etat initial ( = 350  si 20 Decembre ,
!Config         par expl. ,comme ici ) ... A completer
  dayref=1
  CALL getin('dayref', dayref)

!Config  Key  = anneeref
!Config  Desc = Annee de l'etat initial
!Config  Def  = 1998
!Config  Help = Annee de l'etat  initial 
!Config         (   avec  4  chiffres   ) ... A completer
  anneeref = 1998
  CALL getin('anneeref',anneeref)

!Config  Key  = raz_date
!Config  Desc = Remise a zero de la date initiale
!Config  Def  = 0 (pas de remise a zero)
!Config  Help = Remise a zero de la date initiale 
!Config         0 pas de remise a zero, on garde la date du fichier restart
!Config         1 prise en compte de la date de gcm.def avec remise a zero
!Config         des compteurs de pas de temps
  raz_date = 0
  CALL getin('raz_date', raz_date)

!Config  Key  = resetvarc
!Config  Desc = Reinit des variables de controle
!Config  Def  = n
!Config  Help = Reinit des variables de controle
  resetvarc = .false.
  CALL getin('resetvarc',resetvarc)

!Config  Key  = nday
!Config  Desc = Nombre de jours d'integration
!Config  Def  = 10
!Config  Help = Nombre de jours d'integration
!Config         ... On pourait aussi permettre des mois ou des annees !
  nday = 10
  CALL getin('nday',nday)

  ! alternative to specifying nday (see also 'less1day' and 'fractday'
  ! options below: sopecify numbre of dynamic steps to run:
  ndynstep = -9999 ! default value ; if ndynstep <0 then option not used.
  call getin('ndynstep',ndynstep)
      
!Config  Key  = starttime
!Config  Desc = Heure de depart de la simulation
!Config  Def  = 0
!Config  Help = Heure de depart de la simulation
!Config         en jour
  starttime = 0
  CALL getin('starttime',starttime)
      
  ! Mars: time of start for run in "start.nc" (when there are multiple time
  !       steps stored in the file)
  timestart=-9999 ! default value; if <0, use last stored time
  call getin("timestart",timestart)
      
!Config  Key  = ecritstart
!Config  Desc = Mars - frequency of restart.nc output (dyn timesteps)
!Config  Def  = 0
!Config  Help = Mars - frequency of restart.nc output (dyn timesteps)
  ecritstart = 0
  CALL getin('ecritstart',ecritstart)

!Config  Key  = less1day
!Config  Desc = Possibilite d'integrer moins d'un jour
!Config  Def  = n
!Config  Help = Possibilite d'integrer moins d'un jour
  less1day = .false.
  CALL getin('less1day',less1day)

!Config  Key  = fractday
!Config  Desc = integration sur une fraction de jour
!Config  Def  = 0.01
!Config  Help = integration sur une fraction de jour
  fractday = 0.01
  CALL getin('fractday',fractday)

!Config  Key  = day_step
!Config  Desc = nombre de pas par jour
!Config  Def  = 240 
!Config  Help = nombre de pas par jour (multiple de iperiod) (
!Config          ici pour  dt = 1 min ) 
  day_step = 240 
  CALL getin('day_step',day_step)

!Config  Key  = nsplit_phys
!Config  Desc = nombre de subdivisions par pas physique
!Config  Def  = 1 
!Config  Help = nombre de subdivisions par pas physique
  nsplit_phys = 1 
  CALL getin('nsplit_phys',nsplit_phys)

!Config  Key  = iperiod
!Config  Desc = periode pour le pas Matsuno
!Config  Def  = 5
!Config  Help = periode pour le pas Matsuno (en pas de temps)
  iperiod = 5
  CALL getin('iperiod',iperiod)

!Config  Key  = iapp_tracvl
!Config  Desc = frequence du groupement des flux 
!Config  Def  = iperiod
!Config  Help = frequence du groupement des flux (en pas de temps) 
  iapp_tracvl = iperiod
  CALL getin('iapp_tracvl',iapp_tracvl)

! Enforce tracer conservation in integrd & caladvtrac ?
  force_conserv_tracer = .false. 
  CALL getin('force_conserv_tracer',force_conserv_tracer)

!Config  Key  = iconser
!Config  Desc = periode de sortie des variables de controle
!Config  Def  = 240  
!Config  Help = periode de sortie des variables de controle
!Config         (En pas de temps)
  iconser = 240  
  CALL getin('iconser', iconser)

!Config  Key  = iecri
!Config  Desc = periode d'ecriture du fichier histoire
!Config  Def  = 1
!Config  Help = periode d'ecriture du fichier histoire (en jour) 
  iecri = 1
  CALL getin('iecri',iecri)


!Config  Key  = periodav
!Config  Desc = periode de stockage fichier histmoy
!Config  Def  = 1
!Config  Help = periode de stockage fichier histmoy (en jour) 
  periodav = 1.
  CALL getin('periodav',periodav)

!Config  Key  = output_grads_dyn
!Config  Desc = output dynamics diagnostics in 'dyn.dat' file
!Config  Def  = n
!Config  Help = output dynamics diagnostics in Grads-readable 'dyn.dat' file
  output_grads_dyn=.false.
  CALL getin('output_grads_dyn',output_grads_dyn)

!Config  Key  = dissip_period
!Config  Desc = periode de la dissipation 
!Config  Def  = 0
!Config  Help = periode de la dissipation 
!Config  dissip_period=0 => la valeur sera calcule dans inidissip       
!Config  dissip_period>0 => on prend cette valeur
  dissip_period = 0
  call getin('idissip',dissip_period) ! old Mars/Genreic model parameter
  ! if there is a "dissip_period" in run.def, it overrides "idissip"
  CALL getin('dissip_period',dissip_period)

!cc  ....   P. Le Van , modif le 29/04/97 .pour la dissipation  ...
!cc

!Config  Key  = lstardis
!Config  Desc = choix de l'operateur de dissipation
!Config  Def  = y
!Config  Help = choix de l'operateur de dissipation
!Config         'y' si on veut star et 'n' si on veut non-start !
!Config         Moi y en a pas comprendre ! 
  lstardis = .TRUE.
  CALL getin('lstardis',lstardis)


!Config  Key  = nitergdiv
!Config  Desc = Nombre d'iteration de gradiv
!Config  Def  = 1
!Config  Help = nombre d'iterations de l'operateur de dissipation 
!Config         gradiv
  nitergdiv = 1
  CALL getin('nitergdiv',nitergdiv)

!Config  Key  = nitergrot
!Config  Desc = nombre d'iterations de nxgradrot
!Config  Def  = 2
!Config  Help = nombre d'iterations de l'operateur de dissipation  
!Config         nxgradrot
  nitergrot = 2
  CALL getin('nitergrot',nitergrot)


!Config  Key  = niterh
!Config  Desc = nombre d'iterations de divgrad
!Config  Def  = 2
!Config  Help = nombre d'iterations de l'operateur de dissipation
!Config         divgrad
  niterh = 2
  CALL getin('niterh',niterh)


!Config  Key  = tetagdiv
!Config  Desc = temps de dissipation pour div
!Config  Def  = 7200
!Config  Help = temps de dissipation des plus petites longeur 
!Config         d'ondes pour u,v (gradiv)
  tetagdiv = 7200.
  CALL getin('tetagdiv',tetagdiv)

!Config  Key  = tetagrot
!Config  Desc = temps de dissipation pour grad
!Config  Def  = 7200
!Config  Help = temps de dissipation des plus petites longeur 
!Config         d'ondes pour u,v (nxgradrot)
  tetagrot = 7200.
  CALL getin('tetagrot',tetagrot)

!Config  Key  = tetatemp 
!Config  Desc = temps de dissipation pour h
!Config  Def  = 7200
!Config  Help =  temps de dissipation des plus petites longeur 
!Config         d'ondes pour h (divgrad)   
  tetatemp  = 7200.
  CALL getin('tetatemp',tetatemp )

! For Earth model only:
! Parametres controlant la variation sur la verticale des constantes de
! dissipation.
! Pour le moment actifs uniquement dans la version a 39 niveaux
! avec ok_strato=y

  dissip_factz=4.
  dissip_deltaz=10.
  dissip_zref=30.
  CALL getin('dissip_factz',dissip_factz )
  CALL getin('dissip_deltaz',dissip_deltaz )
  CALL getin('dissip_zref',dissip_zref )

! For other planets:
! Parametres controlant la variation sur la verticale des constantes de
! dissipation.
! Actifs uniquement avec ok_strato=y

  dissip_fac_mid=2.
  dissip_fac_up=10.
  dissip_deltaz=10.! Intervalle (km) pour le changement mid / up
  dissip_hdelta=5. ! scale height (km) dans la zone de la transition(m)
  dissip_pupstart=1.e3  ! pression (Pa) au bas la transition mid / up
  CALL getin('dissip_fac_mid',dissip_fac_mid )
  CALL getin('dissip_fac_up',dissip_fac_up )
  CALL getin('dissip_deltaz',dissip_deltaz )
  CALL getin('dissip_hdelta',dissip_hdelta )
  CALL getin('dissip_pupstart',dissip_pupstart )

! top_bound sponge: only active if iflag_top_bound!=0
!                   iflag_top_bound=0 for no sponge
!                   iflag_top_bound=1 for sponge over 4 topmost layers
!                   iflag_top_bound=2 for sponge from top to ~1% of top layer pressure
  iflag_top_bound=0
  CALL getin('iflag_top_bound',iflag_top_bound)

! mode_top_bound : fields towards which sponge relaxation will be done:
!                  mode_top_bound=0: no relaxation
!                  mode_top_bound=1: u and v relax towards 0
!                  mode_top_bound=2: u and v relax towards their zonal mean
!                  mode_top_bound=3: u,v and pot. temp. relax towards their zonal mean
  mode_top_bound=3
  CALL getin('mode_top_bound',mode_top_bound)

! top_bound sponge : inverse of charactericstic relaxation time scale for sponge
  tau_top_bound=1.e-5
  CALL getin('tau_top_bound',tau_top_bound)

! the other possible sponge layer (sponge_mod)
  callsponge=.false. ! default value; don't use the sponge
  call getin("callsponge",callsponge)
  ! check that user is not trying to use both sponge models
  if ((iflag_top_bound.ge.1).and.callsponge) then
    write(lunout,*)'Bad choice of options:'
    write(lunout,*)' iflag_top_bound=',iflag_top_bound
    write(lunout,*)' and callsponge=.true.'
    write(lunout,*)'But both sponge models should not be', &
                   ' used simultaneously!'
    stop
  endif
       
! nsponge: number of atmospheric layers over which the sponge extends
  nsponge=3 ! default value
  call getin("nsponge",nsponge)

! mode_sponge: (quenching is towards ... over the upper nsponge layers)
!      0: (h=hmean,u=v=0)
!      1: (h=hmean,u=umean,v=0)
!      2: (h=hmean,u=umean,v=vmean)"
  mode_sponge=2 ! default value
  call getin("mode_sponge",mode_sponge)

! tetasponge: characteristic time scale (seconds) at topmost layer
!            (time scale then doubles with decreasing layer index)."
  tetasponge=50000.0
  call getin("tetasponge",tetasponge)

! ngroup: to group longitudinaly near the pole (groupe/groupeun routines)
!        (implies that iim has to be a multiple of 2**ngroup)
  ngroup=3
  CALL getin('ngroup',ngroup)

! FOR TITAN: tidal forces
  if (planet_type=="titan") then
    tidal=.TRUE.
    CALL getin('tidal',tidal)
  else
    tidal=.false.
  endif

!Config  Key  = coefdis
!Config  Desc = coefficient pour gamdissip
!Config  Def  = 0
!Config  Help = coefficient pour gamdissip  
  coefdis = 0.
  CALL getin('coefdis',coefdis)

!Config  Key  = purmats
!Config  Desc = Schema d'integration
!Config  Def  = n
!Config  Help = Choix du schema d'integration temporel.
!Config         y = pure Matsuno sinon c'est du Matsuno-leapfrog
  purmats = .FALSE.
  CALL getin('purmats',purmats)

!Config  Key  = ok_guide
!Config  Desc = Guidage
!Config  Def  = n
!Config  Help = Guidage
  ok_guide = .FALSE.
  CALL getin('ok_guide',ok_guide)

!Config  Key  =  read_start
!Config  Desc = Initialize model using a 'start.nc' file
!Config  Def  = y
!Config  Help = y: intialize dynamical fields using a 'start.nc' file
!               n: fields are initialized by 'iniacademic' routine
  read_start= .true.
  CALL getin('read_start',read_start)

!Config  Key  = iflag_phys
!Config  Desc = Avec ls physique 
!Config  Def  = 1
!Config  Help = Permet de faire tourner le modele sans 
!Config         physique.
  iflag_phys = 1
  CALL getin('iflag_phys',iflag_phys)

!Config  Key  =  iphysiq
!Config  Desc = Periode de la physique
!Config  Def  = 5
!Config  Help = Periode de la physique en pas de temps de la dynamique.
  iphysiq = 5
  CALL getin('iphysiq', iphysiq)

!Config  Key  =  cpofT
!Config  Desc = dependence of Cp on T
!Config  Def  = False
!Config  Help = dependence of Cp on T (true or false)
  cpofT = .False.
  if (planet_type.eq."venus") then
   cpofT = .True.
  endif
  CALL getin('cpofT', cpofT)

!Config  Key  = iflag_trac
!Config  Desc = traceurs presents ou non
!Config  Def  = 1
!Config  Help = Permet de faire tourner le modele sans traceurs
!Config         
  iflag_trac = 1 
  CALL getin('iflag_trac',iflag_trac)

!Config  Key  = ip_ebil_dyn
!Config  Desc = PRINT level for energy conserv. diag.
!Config  Def  = 0
!Config  Help = PRINT level for energy conservation diag. ;
!               les options suivantes existent :
!Config         0 pas de print
!Config         1 pas de print
!Config         2 print,
  ip_ebil_dyn = 0
  CALL getin('ip_ebil_dyn',ip_ebil_dyn)

!Config  Key  = offline
!Config  Desc = Nouvelle eau liquide
!Config  Def  = n
!Config  Help = Permet de mettre en route la
!Config         nouvelle parametrisation de l'eau liquide !
  offline = .FALSE.
  CALL getin('offline',offline)

!Config  Key  = type_trac
!Config  Desc = Choix de couplage avec model de chimie INCA ou REPROBUS
!Config  Def  = lmdz
!Config  Help = 
!Config         'lmdz' = pas de couplage, pur LMDZ
!Config         'inca' = model de chime INCA 
!Config         'repr' = model de chime REPROBUS
  type_trac = 'lmdz'
  CALL getin('type_trac',type_trac)

!Config  Key  = config_inca
!Config  Desc = Choix de configuration de INCA
!Config  Def  = none
!Config  Help = Choix de configuration de INCA :
!Config         'none' = sans INCA
!Config         'chem' = INCA avec calcul de chemie
!Config         'aero' = INCA avec calcul des aerosols 
  config_inca = 'none'
  CALL getin('config_inca',config_inca)

!Config  Key  = ok_dynzon 
!Config  Desc = calcul et sortie des transports
!Config  Def  = n 
!Config  Help = Permet de mettre en route le calcul des transports
!Config          
  ok_dynzon = .FALSE.
  CALL getin('ok_dynzon',ok_dynzon) 

!Config  Key  = ok_dyn_ins
!Config  Desc = sorties instantanees dans la dynamique
!Config  Def  = n 
!Config  Help = 
!Config          
  ok_dyn_ins = .FALSE. 
  CALL getin('ok_dyn_ins',ok_dyn_ins) 

!Config  Key  = ok_dyn_ave
!Config  Desc = sorties moyennes dans la dynamique
!Config  Def  = n 
!Config  Help = 
!Config          
  ok_dyn_ave = .FALSE. 
  CALL getin('ok_dyn_ave',ok_dyn_ave) 

!Config  Key  = use_filtre_fft
!Config  Desc = flag d'activation des FFT pour le filtre
!Config  Def  = false
!Config  Help = permet d'activer l'utilisation des FFT pour effectuer
!Config         le filtrage aux poles.
! Le filtre fft n'est pas implemente dans dyn3d 
  use_filtre_fft=.FALSE.
  CALL getin('use_filtre_fft',use_filtre_fft)

  IF (use_filtre_fft) THEN
    write(lunout,*)'STOP !!!'
    write(lunout,*)'use_filtre_fft not implemented in dyn3d'
    STOP 1
  ENDIF
      
!Config key = ok_strato
!Config  Desc = activation de la version strato
!Config  Def  = .FALSE.
!Config  Help = active la version stratosphérique de LMDZ de F. Lott

  ok_strato=.TRUE.
  CALL getin('ok_strato',ok_strato)

! NB: vert_prof_dissip is Earth-specific; should not impact other models
  if (planet_type=="earth") then
    vert_prof_dissip = merge(1, 0, ok_strato .and. llm==39)
    CALL getin('vert_prof_dissip', vert_prof_dissip)
    call assert(vert_prof_dissip == 0 .or. vert_prof_dissip ==  1, &
         "bad value for vert_prof_dissip")
  else
    vert_prof_dissip=0 ! default for planets !
    if (planet_type=="mars") then
      vert_prof_dissip=1 ! use fac_mid & fac_up & startalt & delta
    endif
  endif

!Config  Key  = ok_gradsfile
!Config  Desc = activation des sorties grads du guidage
!Config  Def  = n
!Config  Help = active les sorties grads du guidage

  ok_gradsfile = .FALSE.
  CALL getin('ok_gradsfile',ok_gradsfile)

!Config  Key  = ok_limit
!Config  Desc = creation des fichiers limit dans create_etat0_limit
!Config  Def  = y
!Config  Help = production du fichier limit.nc requise

  ok_limit = .TRUE.
  CALL getin('ok_limit',ok_limit)

!Config  Key  = ok_etat0
!Config  Desc = creation des fichiers etat0 dans create_etat0_limit
!Config  Def  = y
!Config  Help = production des fichiers start.nc, startphy.nc requise

  ok_etat0 = .TRUE.
  CALL getin('ok_etat0',ok_etat0)

!Config  Key  = read_orop
!Config  Desc = lecture du fichier de params orographiques sous maille
!Config  Def  = f
!Config  Help = lecture fichier plutot que grid_noro

  read_orop = .FALSE.
  CALL getin('read_orop',read_orop)

!----------------------------------------
! Parameters for zonal averages in the case of Titan
  moyzon_mu = .false.
  moyzon_ch = .false.
  if (planet_type=="titan") then
    CALL getin('moyzon_mu', moyzon_mu)
    CALL getin('moyzon_ch', moyzon_ch)
  endif
!----------------------------------------

!----------------------------------------
!cc  ....   P. Le Van , ajout  le 7/03/95 .pour le zoom ...
!     .........   (  modif  le 17/04/96 )   .........
!
! ZOOM PARAMETERS ... the ones read in start.nc prevail anyway ! (SL, 2012)
!
!----------------------------------------
  test_etatinit: IF (.not. etatinit) then
     !Config  Key  = clon
     !Config  Desc = centre du zoom, longitude
     !Config  Def  = 0
     !Config  Help = longitude en degres du centre 
     !Config         du zoom
     clonn = 0.
     CALL getin('clon',clonn)

     !Config  Key  = clat
     !Config  Desc = centre du zoom, latitude
     !Config  Def  = 0
     !Config  Help = latitude en degres du centre du zoom
     !Config         
     clatt = 0.
     CALL getin('clat',clatt)

     IF( ABS(clat - clatt).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de clat passee par run.def', &
             ' est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     !Config  Key  = grossismx 
     !Config  Desc = zoom en longitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la longitude
     grossismxx = 1.0
     CALL getin('grossismx',grossismxx)

     IF( ABS(grossismx - grossismxx).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de grossismx passee par ', &
             'run.def est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     !Config  Key  = grossismy
     !Config  Desc = zoom en latitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la latitude
     grossismyy = 1.0
     CALL getin('grossismy',grossismyy)

     IF( ABS(grossismy - grossismyy).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de grossismy passee par ', &
             'run.def est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     IF( grossismx.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismx < 1 .   *** '
        STOP
     ELSE
        alphax = 1. - 1./ grossismx
     ENDIF

     IF( grossismy.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismy < 1 .   *** '
        STOP
     ELSE
        alphay = 1. - 1./ grossismy
     ENDIF

     write(lunout,*)'conf_gcm: alphax alphay',alphax,alphay

     !    alphax et alphay sont les anciennes formulat. des grossissements

     !Config  Key  = fxyhypb
     !Config  Desc = Fonction  hyperbolique
     !Config  Def  = y
     !Config  Help = Fonction  f(y)  hyperbolique  si = .true.  
     !Config         sinon  sinusoidale
     fxyhypbb = .TRUE.
     CALL getin('fxyhypb',fxyhypbb)

     IF( .NOT.fxyhypb )  THEN
        IF( fxyhypbb )     THEN
           write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
           write(lunout,*)' *** fxyhypb lu sur le fichier start est ', &
                'F alors  qu il est  T  sur  run.def  ***'
           STOP
        ENDIF
     ELSE
        IF( .NOT.fxyhypbb )   THEN
           write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
           write(lunout,*)' ***  fxyhypb lu sur le fichier start est ', &
                'T alors  qu il est  F  sur  run.def  ****  '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = dzoomx
     !Config  Desc = extension en longitude
     !Config  Def  = 0
     !Config  Help = extension en longitude  de la zone du zoom  
     !Config         ( fraction de la zone totale)
     dzoomxx = 0.0
     CALL getin('dzoomx',dzoomxx)

     IF( fxyhypb )  THEN
        IF( ABS(dzoomx - dzoomxx).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de dzoomx passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = dzoomy
     !Config  Desc = extension en latitude
     !Config  Def  = 0
     !Config  Help = extension en latitude de la zone  du zoom  
     !Config         ( fraction de la zone totale)
     dzoomyy = 0.0
     CALL getin('dzoomy',dzoomyy)

     IF( fxyhypb )  THEN
        IF( ABS(dzoomy - dzoomyy).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de dzoomy passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = taux
     !Config  Desc = raideur du zoom en  X
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  X
     tauxx = 3.0
     CALL getin('taux',tauxx)

     IF( fxyhypb )  THEN
        IF( ABS(taux - tauxx).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de taux passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = tauyy
     !Config  Desc = raideur du zoom en  Y
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  Y
     tauyy = 3.0
     CALL getin('tauy',tauyy)

     IF( fxyhypb )  THEN
        IF( ABS(tauy - tauyy).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de tauy passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !c
     IF( .NOT.fxyhypb  )  THEN

        !Config  Key  = ysinus
        !Config  IF   = !fxyhypb
        !Config  Desc = Fonction en Sinus
        !Config  Def  = y
        !Config  Help = Fonction  f(y) avec y = Sin(latit.) si = .true. 
        !Config         sinon y = latit.
        ysinuss = .TRUE.
        CALL getin('ysinus',ysinuss)

        IF( .NOT.ysinus )  THEN
           IF( ysinuss )     THEN
              write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
              write(lunout,*)' *** ysinus lu sur le fichier start est F', &
                   ' alors  qu il est  T  sur  run.def  ***'
              STOP
           ENDIF
        ELSE
           IF( .NOT.ysinuss )   THEN
              write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
              write(lunout,*)' *** ysinus lu sur le fichier start est T', &
                   ' alors  qu il est  F  sur  run.def  ****  '
              STOP
           ENDIF
        ENDIF
     ENDIF ! of IF( .NOT.fxyhypb  )
  else
     !Config  Key  = clon
     !Config  Desc = centre du zoom, longitude
     !Config  Def  = 0
     !Config  Help = longitude en degres du centre 
     !Config         du zoom
     clon = 0.
     CALL getin('clon',clon)

     !Config  Key  = clat
     !Config  Desc = centre du zoom, latitude
     !Config  Def  = 0
     !Config  Help = latitude en degres du centre du zoom
     !Config         
     clat = 0.
     CALL getin('clat',clat)

     !Config  Key  = grossismx 
     !Config  Desc = zoom en longitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la longitude
     grossismx = 1.0
     CALL getin('grossismx',grossismx)

     !Config  Key  = grossismy
     !Config  Desc = zoom en latitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la latitude
     grossismy = 1.0
     CALL getin('grossismy',grossismy)

     IF( grossismx.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismx < 1 .   *** '
        STOP
     ELSE
        alphax = 1. - 1./ grossismx
     ENDIF

     IF( grossismy.LT.1. )  THEN
        write(lunout,*) 'conf_gcm: ***ATTENTION !! grossismy < 1 . *** '
        STOP
     ELSE
        alphay = 1. - 1./ grossismy
     ENDIF

     write(lunout,*)'conf_gcm: alphax alphay ',alphax,alphay

     !    alphax et alphay sont les anciennes formulat. des grossissements

     !Config  Key  = fxyhypb
     !Config  Desc = Fonction  hyperbolique
     !Config  Def  = y
     !Config  Help = Fonction  f(y)  hyperbolique  si = .true.  
     !Config         sinon  sinusoidale
     fxyhypb = .TRUE.
     CALL getin('fxyhypb',fxyhypb)

     !Config  Key  = dzoomx
     !Config  Desc = extension en longitude
     !Config  Def  = 0
     !Config  Help = extension en longitude  de la zone du zoom  
     !Config         ( fraction de la zone totale)
     dzoomx = 0.2
     CALL getin('dzoomx',dzoomx)
     call assert(dzoomx < 1, "conf_gcm: dzoomx must be < 1")

     !Config  Key  = dzoomy
     !Config  Desc = extension en latitude
     !Config  Def  = 0
     !Config  Help = extension en latitude de la zone  du zoom  
     !Config         ( fraction de la zone totale)
     dzoomy = 0.2
     CALL getin('dzoomy',dzoomy)
     call assert(dzoomy < 1, "conf_gcm: dzoomy must be < 1")

     !Config  Key  = taux
     !Config  Desc = raideur du zoom en  X
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  X
     taux = 3.0
     CALL getin('taux',taux)

     !Config  Key  = tauy
     !Config  Desc = raideur du zoom en  Y
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  Y
     tauy = 3.0
     CALL getin('tauy',tauy)

     !Config  Key  = ysinus
     !Config  IF   = !fxyhypb
     !Config  Desc = Fonction en Sinus
     !Config  Def  = y
     !Config  Help = Fonction  f(y) avec y = Sin(latit.) si = .true. 
     !Config         sinon y = latit.
     ysinus = .TRUE.
     CALL getin('ysinus',ysinus)

  end IF test_etatinit
!----------------------------------------


 write(lunout,*)' #########################################'
 write(lunout,*)' Configuration des parametres lus via run.def '
 write(lunout,*)' planet_type = ', planet_type
 write(lunout,*)' calend = ', calend
 write(lunout,*)' dayref = ', dayref
 write(lunout,*)' anneeref = ', anneeref
 write(lunout,*)' nday = ', nday
 if (ndynstep.ne.-9999) write(lunout,*)' ndynstep = ', ndynstep
 if (less1day) write(lunout,*)' fractday = ', fractday
 write(lunout,*)' day_step = ', day_step
 write(lunout,*)' iperiod = ', iperiod
 write(lunout,*)' nsplit_phys = ', nsplit_phys
 write(lunout,*)' iconser = ', iconser
 write(lunout,*)' iecri = ', iecri
 write(lunout,*)' periodav = ', periodav 
 write(lunout,*)' output_grads_dyn = ', output_grads_dyn
 write(lunout,*)' dissip_period = ', dissip_period
 write(lunout,*)' lstardis = ', lstardis
 write(lunout,*)' nitergdiv = ', nitergdiv
 write(lunout,*)' nitergrot = ', nitergrot
 write(lunout,*)' niterh = ', niterh
 write(lunout,*)' tetagdiv = ', tetagdiv
 write(lunout,*)' tetagrot = ', tetagrot
 write(lunout,*)' tetatemp = ', tetatemp
 write(lunout,*)' coefdis = ', coefdis
 write(lunout,*)' purmats = ', purmats
 write(lunout,*)' read_start = ', read_start
 write(lunout,*)' iflag_phys = ', iflag_phys
 write(lunout,*)' iphysiq = ', iphysiq
 write(lunout,*)' cpofT = ', cpofT
 write(lunout,*)' iflag_trac = ', iflag_trac
 write(lunout,*)' iapp_tracvl = ', iapp_tracvl
 write(lunout,*)' force_conserv_tracer = ', force_conserv_tracer
 write(lunout,*)' clon = ', clon
 write(lunout,*)' clat = ', clat
 write(lunout,*)' grossismx = ', grossismx
 write(lunout,*)' grossismy = ', grossismy
 write(lunout,*)' fxyhypb = ', fxyhypb
 write(lunout,*)' dzoomx = ', dzoomx
 write(lunout,*)' dzoomy = ', dzoomy
 write(lunout,*)' taux = ', taux
 write(lunout,*)' tauy = ', tauy
 write(lunout,*)' offline = ', offline
 write(lunout,*)' type_trac = ', type_trac
 write(lunout,*)' config_inca = ', config_inca
 write(lunout,*)' ok_dynzon = ', ok_dynzon
 write(lunout,*)' ok_dyn_ins = ', ok_dyn_ins 
 write(lunout,*)' ok_dyn_ave = ', ok_dyn_ave 
 write(lunout,*)' ok_strato = ', ok_strato
 write(lunout,*)' ok_gradsfile = ', ok_gradsfile
 write(lunout,*)' ok_limit = ', ok_limit
 write(lunout,*)' ok_etat0 = ', ok_etat0
 write(lunout,*)' read_orop = ', read_orop
 if (planet_type=="titan") then
   write(lunout,*)' moyzon_mu = ', moyzon_mu
   write(lunout,*)' moyzon_ch = ', moyzon_ch
 endif

END SUBROUTINE conf_gcm
