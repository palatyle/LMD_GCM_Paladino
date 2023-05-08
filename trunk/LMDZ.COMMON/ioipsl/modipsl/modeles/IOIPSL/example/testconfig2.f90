PROGRAM testconfig2
!-
!$Id: testconfig2.f90 386 2008-09-04 08:38:48Z bellier $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
  !
  USE getincom
  !
  !
  !      This program will do some basic tests on the getin module
  !
  !
  IMPLICIT NONE
  !
  INTEGER :: dayref=181
  INTEGER :: anneeref = 1998
  INTEGER :: nday = 10
  INTEGER :: day_step = 240 
  INTEGER :: iperiod = 5
  INTEGER :: iconser = 240  
  INTEGER :: iecri = 1
  INTEGER :: idissip = 10
  INTEGER :: nitergdiv = 1
  INTEGER :: nitergrot = 2
  INTEGER :: niterh = 2
  INTEGER :: iphysiq = 5
  INTEGER :: ecritphy = 1
  INTEGER :: nbapp_rad = 12
  INTEGER :: iflag_con = 2

  REAL :: periodav = 1.
  REAL :: tetagdiv = 7200.
  REAL :: tetagrot = 7200.
  REAL :: tetatemp  = 7200.
  REAL :: coefdis = 0.
  REAL :: clonn = 0.
  REAL :: clatt = 0.
  REAL :: grossismxx = 1.0
  REAL :: grossismyy = 1.0
  REAL :: dzoomxx = 0.0
  REAL :: dzoomyy = 0.0
  REAL :: tauxx = 3.0
  REAL :: tauyy = 3.0
  REAL :: clon = 0.
  REAL :: clat = 0.
  REAL :: grossismx = 1.0
  REAL :: grossismy = 1.0
  REAL :: dzoomx = 0.0
  REAL :: dzoomy = 0.0
  REAL :: taux = 3.0
  REAL :: tauy = 3.0

  LOGICAL :: lstardis = .TRUE.
  LOGICAL :: purmats = .FALSE.
  LOGICAL :: physic = .TRUE.
  LOGICAL :: cycle_diurne = .TRUE.
  LOGICAL :: soil_model = .TRUE.
  LOGICAL :: new_oliq = .TRUE.
  LOGICAL :: ok_orodr = .TRUE.
  LOGICAL :: ok_orolf = .TRUE.
  LOGICAL :: ok_limitvrai = .FALSE.
  LOGICAL :: fxyhypbb = .TRUE.
  LOGICAL :: ysinuss = .TRUE.
  LOGICAL :: fxyhypb = .TRUE.
  LOGICAL :: ysinus = .TRUE.
!
! For time keeping
!
  INTEGER :: date_time(8), date_start, date_end
  CHARACTER(LEN=10) big_ben(3)
!
!
!
  CALL DATE_AND_TIME(big_ben(1), big_ben(2), big_ben(3), date_time)
  date_start = date_time(8) + 1000*date_time(7) + 60*1000*date_time(7)
!
!Config  Key  = dayref
!Config  Desc = Jour de l'etat initial
!Config  Def  = 181
!Config  Help = Jour de l'etat initial ( = 350  si 20 Decembre ,
!Config         par expl. ,comme ici ) ... A completer
  dayref=181
  CALL getin('dayref', dayref)

!Config  Key  = anneeref
!Config  Desc = Annee de l'etat initial
!Config  Def  = 1998
!Config  Help = Annee de l'etat  initial 
!Config         (   avec  4  chiffres   ) ... A completer
  anneeref = 1998
  CALL getin('anneeref',anneeref)

!Config  Key  = nday
!Config  Desc = Nombre de jours d'integration
!Config  Def  = 10
!Config  Help = Nombre de jours d'integration
!Config         ... On pourait aussi permettre des mois ou des annees !
  nday = 10
  CALL getin('nday',nday)

!Config  Key  = day_step
!Config  Desc = nombre de pas par jour
!Config  Def  = 240 
!Config  Help = nombre de pas par jour (multiple de iperiod) (
!Config          ici pour  dt = 1 min ) 
  day_step = 240 
  CALL getin('day_step',day_step)

!Config  Key  = iperiod
!Config  Desc = periode pour le pas Matsuno 
!Config  Def  = 5
!Config  Help = periode pour le pas Matsuno (en pas de temps) 
  iperiod = 5
  CALL getin('iperiod',iperiod)


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

!Config  Key  = idissip
!Config  Desc = periode de la dissipation 
!Config  Def  = 10
!Config  Help = periode de la dissipation 
!Config         (en pas) ... a completer !
  idissip = 10
  CALL getin('idissip',idissip)

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

!Config  Key  = physic
!Config  Desc = Avec ls physique 
!Config  Def  = y
!Config  Help = Permet de faire tourner le modele sans 
!Config         physique.
  physic = .TRUE.
  CALL getin('physic',physic)


!Config  Key  =  iphysiq
!Config  Desc = Periode de la physique
!Config  Def  = 5
!Config  Help = Periode de la physique en pas de temps de la dynamique.
  iphysiq = 5
  CALL getin('iphysiq', iphysiq)

!Config  Key  = ecritphy
!Config  Desc = Frequence d'ecriture de la physique
!Config  Def  = 1
!Config  Help = frequence  de l'ecriture du fichier histphy
!Config         en jours.
  ecritphy = 1
  CALL getin('ecritphy',ecritphy)

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

!Config  Key  = new_oliq
!Config  Desc = Nouvelle eau liquide
!Config  Def  = y
!Config  Help = Permet de mettre en route la
!Config         nouvelle parametrisation de l'eau liquide !
  new_oliq = .TRUE.
  CALL getin('new_oliq',new_oliq)

!Config  Key  = ok_orodr
!Config  Desc = Orodr ???
!Config  Def  = y
!Config  Help = Y en a pas comprendre !
!Config         
  ok_orodr = .TRUE.
  CALL getin('ok_orodr',ok_orodr)

!Config  Key  =  ok_orolf
!Config  Desc = Orolf ??
!Config  Def  = y
!Config  Help = Connais pas !
  ok_orolf = .TRUE.
  CALL getin('ok_orolf', ok_orolf)

!Config  Key  = ok_limitvrai
!Config  Desc = Force la lecture de la bonne annee
!Config  Def  = n
!Config  Help = On peut forcer le modele a lire le
!Config         fichier SST de la bonne annee. C'est une tres bonne
!Config         idee, pourquoi ne pas mettre toujours a y ???
  ok_limitvrai = .FALSE.
  CALL getin('ok_limitvrai',ok_limitvrai)

!Config  Key  = nbapp_rad
!Config  Desc = Frequence d'appel au rayonnement
!Config  Def  = 12
!Config  Help = Nombre  d'appels des routines de rayonnements
!Config         par jour.
  nbapp_rad = 12
  CALL getin('nbapp_rad',nbapp_rad)

!Config  Key  = iflag_con
!Config  Desc = Flag de convection
!Config  Def  = 2
!Config  Help = Flag  pour la convection les options suivantes existent :
!Config         1 pour LMD,
!Config         2 pour Tiedtke,
!Config         3 pour CCM(NCAR)  
  iflag_con = 2
  CALL getin('iflag_con',iflag_con)

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

!Config  Key  = grossismx 
!Config  Desc = zoom en longitude
!Config  Def  = 1.0
!Config  Help = facteur de grossissement du zoom,
!Config         selon la longitude
  grossismxx = 1.0
  CALL getin('grossismx',grossismxx)

!Config  Key  = grossismy
!Config  Desc = zoom en latitude
!Config  Def  = 1.0
!Config  Help = facteur de grossissement du zoom,
!Config         selon la latitude
  grossismyy = 1.0
  CALL getin('grossismy',grossismyy)

!Config  Key  = fxyhypb
!Config  Desc = Fonction  hyperbolique
!Config  Def  = y
!Config  Help = Fonction  f(y)  hyperbolique  si = .true.  
!Config         sinon  sinusoidale
  fxyhypbb = .TRUE.
  CALL getin('fxyhypb',fxyhypbb)

!Config  Key  = dzoomx
!Config  Desc = extension en longitude
!Config  Def  = 0
!Config  Help = extension en longitude  de la zone du zoom  
!Config         ( fraction de la zone totale)
  dzoomxx = 0.0
  CALL getin('dzoomx',dzoomxx)

!Config  Key  = dzoomy
!Config  Desc = extension en latitude
!Config  Def  = 0
!Config  Help = extension en latitude de la zone  du zoom  
!Config         ( fraction de la zone totale)
  dzoomyy = 0.0
  CALL getin('dzoomy',dzoomyy)
      
!Config  Key  = taux
!Config  Desc = raideur du zoom en  X
!Config  Def  = 3
!Config  Help = raideur du zoom en  X
  tauxx = 3.0
  CALL getin('taux',tauxx)

!Config  Key  = tauyy
!Config  Desc = raideur du zoom en  Y
!Config  Def  = 3
!Config  Help = raideur du zoom en  Y
  tauyy = 3.0
  CALL getin('tauy',tauyy)

!Config  Key  = ysinus
!Config  IF   = !fxyhypb
!Config  Desc = Fonction en Sinus
!Config  Def  = y
!Config  Help = Fonction  f(y) avec y = Sin(latit.) si = .true. 
!Config         sinon y = latit.
  ysinuss = .TRUE.
  CALL getin('ysinus',ysinuss)

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
  dzoomx = 0.0
  CALL getin('dzoomx',dzoomx)

!Config  Key  = dzoomy
!Config  Desc = extension en latitude
!Config  Def  = 0
!Config  Help = extension en latitude de la zone  du zoom  
!Config         ( fraction de la zone totale)
  dzoomy = 0.0
  CALL getin('dzoomy',dzoomy)

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

!!
!!
  CALL DATE_AND_TIME(big_ben(1), big_ben(2), big_ben(3), date_time)
  date_end = date_time(8) + 1000*date_time(7) + 60*1000*date_time(7)
  WRITE(*,*) 'Time before dump in milli-sec :', date_end - date_start
!!
  CALL getin_dump()

!!
  CALL DATE_AND_TIME(big_ben(1), big_ben(2), big_ben(3), date_time)
  date_end = date_time(8) + 1000*date_time(7) + 60*1000*date_time(7)
  WRITE(*,*) 'Total time in milli-sec :', date_end - date_start

END PROGRAM testconfig2
