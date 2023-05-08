










!
! $Id: gcm.F 1446 2010-10-22 09:27:25Z emillour $
!
!
!
PROGRAM gcm

  ! if not using IOIPSL, we still need to use (a local version of) getin
  USE ioipsl_getincom



  USE filtreg_mod
  USE infotrac
  USE control_mod, only: planet_type,nday,day_step,iperiod,iphysiq, &
                             raz_date,anneeref,starttime,dayref,    &
                             ok_dyn_ins,ok_dyn_ave,iecri,periodav,  &
                             less1day,fractday,ndynstep,nsplit_phys
  USE mod_const_mpi, ONLY: COMM_LMDZ
  use cpdet_mod, only: ini_cpdet
  USE temps_mod, ONLY: calend,start_time,annee_ref,day_ref, &
     		itau_dyn,itau_phy,day_ini,jD_ref,jH_ref,day_end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! A nettoyer. On ne veut qu'une ou deux routines d'interface 
! dynamique -> physique pour l'initialisation
! Ehouarn: the following are needed with (parallel) physics:
  USE iniphysiq_mod, ONLY: iniphysiq

  USE comconst_mod, ONLY: daysec,dtvr,dtphys,rad,g,r,cpp
  USE logic_mod, ONLY: read_start,iflag_phys,ok_guide,ecripar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !      ......   Version  du 10/01/98    ..........

  !             avec  coordonnees  verticales hybrides 
  !   avec nouveaux operat. dissipation * ( gradiv2,divgrad2,nxgraro2 )

  !=======================================================================
  !
  !   Auteur:  P. Le Van /L. Fairhead/F.Hourdin
  !   -------
  !
  !   Objet:
  !   ------
  !
  !   GCM LMD nouvelle grille
  !
  !=======================================================================
  !
  !  ... Dans inigeom , nouveaux calculs pour les elongations  cu , cv
  !      et possibilite d'appeler une fonction f(y)  a derivee tangente
  !      hyperbolique a la  place de la fonction a derivee sinusoidale.
  !  ... Possibilite de choisir le schema pour l'advection de
  !        q  , en modifiant iadv dans traceur.def  (MAF,10/02) .
  !
  !      Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron,10/99)
  !      Pour Van-Leer iadv=10
  !
  !-----------------------------------------------------------------------
  !   Declarations:
  !   -------------

  include "dimensions.h"
  include "paramet.h"
  include "comdissnew.h"
  include "comgeom.h"
!!!!!!!!!!!#include "control.h"
!#include "com_io_dyn.h"
  include "iniprint.h"
  include "tracstoke.h"

  REAL zdtvr

  !   variables dynamiques
  REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
  REAL teta(ip1jmp1,llm)                 ! temperature potentielle 
  REAL, ALLOCATABLE, DIMENSION(:,:,:):: q! champs advectes
  REAL ps(ip1jmp1)                       ! pression  au sol
  REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
  REAL masse(ip1jmp1,llm)                ! masse d'air
  REAL phis(ip1jmp1)                     ! geopotentiel au sol
  REAL phi(ip1jmp1,llm)                  ! geopotentiel
  REAL w(ip1jmp1,llm)                    ! vitesse verticale

  ! variables dynamiques intermediaire pour le transport

  !   variables pour le fichier histoire
  REAL dtav      ! intervalle de temps elementaire

  REAL time_0

  LOGICAL lafin
  INTEGER ij,iq,l,i,j


  real time_step, t_wrt, t_ops

  LOGICAL first

  !      LOGICAL call_iniphys
  !      data call_iniphys/.true./

  !+jld variables test conservation energie
  !      REAL ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
  !     Tendance de la temp. potentiel d (theta)/ d t due a la 
  !     tansformation d'energie cinetique en energie thermique
  !     cree par la dissipation
  REAL dhecdt(ip1jmp1,llm)
  !      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
  !      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
  CHARACTER (len=15) :: ztit
  !-jld 


  character (len=80) :: dynhist_file, dynhistave_file
  character (len=20) :: modname
  character (len=80) :: abort_message
  ! locales pour gestion du temps
  INTEGER :: an, mois, jour
  REAL :: heure
  logical use_filtre_fft

!-----------------------------------------------------------------------
!   Initialisations:
!   ----------------

  abort_message = 'last timestep reached'
  modname = 'gcm'
  lafin    = .FALSE.
  dynhist_file = 'dyn_hist.nc'
  dynhistave_file = 'dyn_hist_ave.nc'



!----------------------------------------------------------------------
!  lecture des fichiers gcm.def ou run.def
!  ---------------------------------------
!
  CALL conf_gcm( 99, .TRUE. )
  if (mod(iphysiq, iperiod) /= 0) call abort_gcm("conf_gcm", &
       "iphysiq must be a multiple of iperiod", 1)

  use_filtre_fft=.FALSE.
  CALL getin('use_filtre_fft',use_filtre_fft)
  IF (use_filtre_fft) call abort_gcm("gcm",'FFT filter is not available in the ' &
          // 'sequential version of the dynamics.', 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation de XIOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/02
! A nettoyer. On ne veut qu'une ou deux routines d'interface 
! dynamique -> physique pour l'initialisation
!#ifdef 1
!  CALL init_phys_lmdz(iim,jjp1,llm,1,(/(jjm-1)*iim+2/))
!!      call initcomgeomphy ! now done in iniphysiq
!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Initialisations pour Cp(T) Venus
  call ini_cpdet
!
!-----------------------------------------------------------------------
!   Choix du calendrier
!   -------------------

!      calend = 'earth_365d'

!-----------------------------------------------------------------------
  !
  !
  !------------------------------------
  !   Initialisation partie parallele
  !------------------------------------

  !
  !
  !-----------------------------------------------------------------------
  !   Initialisation des traceurs
  !   ---------------------------
  !  Choix du nombre de traceurs et du schema pour l'advection
  !  dans fichier traceur.def, par default ou via INCA
  call infotrac_init

  ! Allocation de la tableau q : champs advectes   
  allocate(q(ip1jmp1,llm,nqtot))

  !-----------------------------------------------------------------------
  !   Lecture de l'etat initial :
  !   ---------------------------

  !  lecture du fichier start.nc
  if (read_start) then
     ! we still need to run iniacademic to initialize some
     ! constants & fields, if we run the 'newtonian' or 'SW' cases:
     if (iflag_phys.ne.1) then
        CALL iniacademic(vcov,ucov,teta,q,masse,ps,phis,time_0)
     endif

     CALL dynetat0("start.nc",vcov,ucov, &
                    teta,q,masse,ps,phis, time_0)
       
     ! Load relaxation fields (simple nudging). AS 09/2013
     ! ---------------------------------------------------
     if (planet_type.eq."generic") then
       if (ok_guide) then
         CALL relaxetat0("relax.nc")
       endif
     endif
 
     !       write(73,*) 'ucov',ucov
     !       write(74,*) 'vcov',vcov
     !       write(75,*) 'teta',teta
     !       write(76,*) 'ps',ps
     !       write(77,*) 'q',q

  endif ! of if (read_start)


  ! le cas echeant, creation d un etat initial
  IF (prt_level > 9) WRITE(lunout,*) &
       'GCM: AVANT iniacademic AVANT AVANT AVANT AVANT'
  if (.not.read_start) then
     CALL iniacademic(vcov,ucov,teta,q,masse,ps,phis,time_0)
  endif


  !-----------------------------------------------------------------------
  !   Lecture des parametres de controle pour la simulation :
  !   -------------------------------------------------------
  !  on recalcule eventuellement le pas de temps

  IF(MOD(day_step,iperiod).NE.0) THEN
     abort_message = &
       'Il faut choisir un nb de pas par jour multiple de iperiod'
     call abort_gcm(modname,abort_message,1)
  ENDIF

  IF(MOD(day_step,iphysiq).NE.0) THEN
     abort_message = &
       'Il faut choisir un nb de pas par jour multiple de iphysiq'
     call abort_gcm(modname,abort_message,1)
  ENDIF

  zdtvr    = daysec/REAL(day_step)
  IF(dtvr.NE.zdtvr) THEN
     WRITE(lunout,*) &
          'WARNING!!! changement de pas de temps',dtvr,'>',zdtvr
  ENDIF

  !
  ! on remet le calendrier � zero si demande
  !
  IF (start_time /= starttime) then
     WRITE(lunout,*)' GCM: Attention l''heure de depart lue dans le' &
     ,' fichier restart ne correspond pas � celle lue dans le run.def'
     IF (raz_date == 1) then
        WRITE(lunout,*)'Je prends l''heure lue dans run.def'
        start_time = starttime
     ELSE
        call abort_gcm("gcm", "'Je m''arrete'", 1)
     ENDIF
  ENDIF
  IF (raz_date == 1) THEN
     annee_ref = anneeref
     day_ref = dayref
     day_ini = dayref
     itau_dyn = 0
     itau_phy = 0
     time_0 = 0.
     write(lunout,*) &
         'GCM: On reinitialise a la date lue dans gcm.def'
  ELSE IF (annee_ref .ne. anneeref .or. day_ref .ne. dayref) THEN
     write(lunout,*) &
        'GCM: Attention les dates initiales lues dans le fichier'
     write(lunout,*) &
        ' restart ne correspondent pas a celles lues dans '
     write(lunout,*)' gcm.def'
     write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
     write(lunout,*)' day_ref=',day_ref," dayref=",dayref
     write(lunout,*)' Pas de remise a zero'
  ENDIF

!      if (annee_ref .ne. anneeref .or. day_ref .ne. dayref) then
!        write(lunout,*)
!     .  'GCM: Attention les dates initiales lues dans le fichier'
!        write(lunout,*)
!     .  ' restart ne correspondent pas a celles lues dans '
!        write(lunout,*)' gcm.def'
!        write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
!        write(lunout,*)' day_ref=',day_ref," dayref=",dayref
!        if (raz_date .ne. 1) then
!          write(lunout,*)
!     .    'GCM: On garde les dates du fichier restart'
!        else
!          annee_ref = anneeref
!          day_ref = dayref
!          day_ini = dayref
!          itau_dyn = 0
!          itau_phy = 0
!          time_0 = 0.
!          write(lunout,*)
!     .   'GCM: On reinitialise a la date lue dans gcm.def'
!        endif
!      ELSE
!        raz_date = 0
!      endif

  ! Ehouarn: we still need to define JD_ref and JH_ref
  ! and since we don't know how many days there are in a year
  ! we set JD_ref to 0 (this should be improved ...)
  jD_ref=0
  jH_ref=0

  if (iflag_phys.eq.1) then
     ! these initialisations have already been done (via iniacademic)
     ! if running in SW or Newtonian mode
     !-----------------------------------------------------------------------
     !   Initialisation des constantes dynamiques :
     !   ------------------------------------------
     dtvr = zdtvr
     CALL iniconst

     !-----------------------------------------------------------------------
     !   Initialisation de la geometrie :
     !   --------------------------------
     CALL inigeom

     !-----------------------------------------------------------------------
     !   Initialisation du filtre :
     !   --------------------------
     CALL inifilr
  endif ! of if (iflag_phys.eq.1)
  !
  !-----------------------------------------------------------------------
  !   Initialisation de la dissipation :
  !   ----------------------------------

  CALL inidissip( lstardis, nitergdiv, nitergrot, niterh   , &
                  tetagdiv, tetagrot , tetatemp, vert_prof_dissip)

  !-----------------------------------------------------------------------
  !   Initialisation des I/O :
  !   ------------------------

  if (nday>=0) then ! standard case
     day_end=day_ini+nday
  else ! special case when nday <0, run -nday dynamical steps
     day_end=day_ini-nday/day_step
  endif
  if (less1day) then
     day_end=day_ini+floor(time_0+fractday)
  endif
  if (ndynstep.gt.0) then
     day_end=day_ini+floor(time_0+float(ndynstep)/float(day_step))
  endif
      
  WRITE(lunout,'(a,i7,a,i7)') &
               "run from day ",day_ini,"  to day",day_end



  !-----------------------------------------------------------------------
  !   Initialisation de la physique :
  !   -------------------------------

  IF ((iflag_phys==1).or.(iflag_phys>=100)) THEN
     ! Physics:
     CALL iniphysiq(iim,jjm,llm, &
          (jjm-1)*iim+2,comm_lmdz, &
          daysec,day_ini,dtphys/nsplit_phys, &
          rlatu,rlatv,rlonu,rlonv,aire,cu,cv,rad,g,r,cpp, &
          iflag_phys)
  ENDIF ! of IF ((iflag_phys==1).or.(iflag_phys>=100))


  if (planet_type=="mars") then
    ! For Mars we transmit day_ini
    CALL dynredem0("restart.nc", day_ini, phis)
  else
    CALL dynredem0("restart.nc", day_end, phis)
  endif
  ecripar = .TRUE.

! #endif of #ifdef CPP_IOIPSL

  !  Choix des frequences de stokage pour le offline
  !      istdyn=day_step/4     ! stockage toutes les 6h=1jour/4
  !      istdyn=day_step/12     ! stockage toutes les 2h=1jour/12
  istdyn=day_step/4     ! stockage toutes les 6h=1jour/12
  istphy=istdyn/iphysiq     


  !
  !-----------------------------------------------------------------------
  !   Integration temporelle du modele :
  !   ----------------------------------

  !       write(78,*) 'ucov',ucov
  !       write(78,*) 'vcov',vcov
  !       write(78,*) 'teta',teta
  !       write(78,*) 'ps',ps
  !       write(78,*) 'q',q


  CALL leapfrog(ucov,vcov,teta,ps,masse,phis,q,time_0)

END PROGRAM gcm

