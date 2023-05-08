! $Id$
!
MODULE infotrac

! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
! CR: add number of tracers for water (for Earth model only!!)
  INTEGER, SAVE :: nqo

! nbtr : number of tracers not including higher order of moment or water vapor or liquid
!        number of tracers used in the physics
  INTEGER, SAVE :: nbtr

! CRisi: nb of father tracers (i.e. directly advected by air)
  INTEGER, SAVE :: nqperes

! Name variables
  CHARACTER(len=30), ALLOCATABLE, DIMENSION(:), SAVE :: tname ! tracer short name for restart and diagnostics
  CHARACTER(len=33), ALLOCATABLE, DIMENSION(:), SAVE :: ttext ! tracer long name for diagnostics

! iadv  : index of trasport schema for each tracer
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: iadv

! niadv : vector keeping the coorspondance between all tracers(nqtot) treated in the 
!         dynamic part of the code and the tracers (nbtr+2) used in the physics part of the code. 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: niadv ! equivalent dyn / physique

! CRisi: arrays for sons
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: nqfils
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: nqdesc ! number of sons + all gran-sons over all generations
  INTEGER, SAVE :: nqdesc_tot
  INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE    :: iqfils
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: iqpere

! conv_flg(it)=0 : convection desactivated for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: conv_flg
! pbl_flg(it)=0  : boundary layer diffusion desactivaded for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: pbl_flg

  CHARACTER(len=4),SAVE :: type_trac
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE, SAVE :: solsym

    ! CRisi: specific stuff for isotopes
    LOGICAL,SAVE :: ok_isotopes,ok_iso_verif,ok_isotrac,ok_init_iso
    INTEGER :: niso_possibles   
    PARAMETER ( niso_possibles=5) ! 5 possible water isotopes
    REAL, DIMENSION (niso_possibles),SAVE :: tnat,alpha_ideal
    LOGICAL, DIMENSION(niso_possibles),SAVE ::  use_iso
    INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE ::  iqiso ! donne indice iq en fn de (ixt,phase) 
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  iso_num ! donne numéro iso entre 1 et niso_possibles en fn de nqtot
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  iso_indnum ! donne numéro iso entre 1 et niso effectif en fn de nqtot
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  zone_num ! donne numéro de la zone de tracage en fn de nqtot
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE ::  phase_num ! donne numéro de la zone de tracage en fn de nqtot
    INTEGER, DIMENSION(niso_possibles), SAVE :: indnum_fn_num ! donne indice entre entre 1 et niso en fonction du numéro d isotope entre 1 et niso_possibles
    INTEGER, ALLOCATABLE, DIMENSION(:,:), SAVE ::  index_trac ! numéro ixt en fn izone, indnum entre 1 et niso
    INTEGER,SAVE :: niso,ntraceurs_zone,ntraciso

#ifdef CPP_StratAer
!--CK/OB for stratospheric aerosols
  INTEGER, SAVE :: nbtr_bin
  INTEGER, SAVE :: nbtr_sulgas
  INTEGER, SAVE :: id_OCS_strat
  INTEGER, SAVE :: id_SO2_strat
  INTEGER, SAVE :: id_H2SO4_strat
  INTEGER, SAVE :: id_BIN01_strat
  INTEGER, SAVE :: id_TEST_strat
#endif

CONTAINS

  SUBROUTINE infotrac_init
    USE control_mod, ONLY: planet_type, config_inca
#ifdef REPROBUS
    USE CHEM_REP, ONLY : Init_chem_rep_trac
#endif
    IMPLICIT NONE
!=======================================================================
!
!   Auteur:  P. Le Van /L. Fairhead/F.Hourdin
!   -------
!   Modif special traceur F.Forget 05/94
!   Modif M-A Filiberti 02/02 lecture de traceur.def
!
!   Objet:
!   ------
!   GCM LMD nouvelle grille
!
!=======================================================================
!   ... modification de l'integration de q ( 26/04/94 ) ....
!-----------------------------------------------------------------------
! Declarations

    INCLUDE "dimensions.h"
    INCLUDE "iniprint.h"

! Local variables
    INTEGER, ALLOCATABLE, DIMENSION(:) :: hadv  ! index of horizontal trasport schema
    INTEGER, ALLOCATABLE, DIMENSION(:) :: vadv  ! index of vertical trasport schema

    INTEGER, ALLOCATABLE, DIMENSION(:) :: hadv_inca  ! index of horizontal trasport schema
    INTEGER, ALLOCATABLE, DIMENSION(:) :: vadv_inca  ! index of vertical trasport schema

    CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: tnom_0  ! tracer short name
    CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: tnom_transp ! transporting fluid short name: CRisi
    CHARACTER(len=3), DIMENSION(30) :: descrq
    CHARACTER(len=1), DIMENSION(3)  :: txts
    CHARACTER(len=2), DIMENSION(9)  :: txtp
    CHARACTER(len=30)               :: str1,str2
  
    INTEGER :: nqtrue  ! number of tracers read from tracer.def, without higer order of moment
    INTEGER :: iq, new_iq, iiq, jq, ierr, ierr2, ierr3
    INTEGER :: ifils,ipere,generation ! CRisi
    LOGICAL :: continu,nouveau_traceurdef
    INTEGER :: IOstatus ! gestion de la retrocompatibilite de traceur.def
    CHARACTER(len=30) :: tchaine    
    
    character(len=80) :: line ! to store a line of text
 
    character(len=*),parameter :: modname="infotrac_init"
!-----------------------------------------------------------------------
! Initialization :
!
    txts=(/'x','y','z'/)
    txtp=(/'x ','y ','z ','xx','xy','xz','yy','yz','zz'/)

    descrq(14)='VLH'
    descrq(10)='VL1'
    descrq(11)='VLP'
    descrq(12)='FH1'
    descrq(13)='FH2'
    descrq(16)='PPM'
    descrq(17)='PPS'
    descrq(18)='PPP'
    descrq(20)='SLP'
    descrq(30)='PRA'
    
    IF (planet_type=='earth') THEN
     ! Coherence test between parameter type_trac, config_inca and preprocessing keys
     IF (type_trac=='inca') THEN
       WRITE(lunout,*) 'You have choosen to couple with INCA chemestry model : type_trac=', &
            type_trac,' config_inca=',config_inca
       IF (config_inca/='aero' .AND. config_inca/='aeNP' .AND. config_inca/='chem') THEN
          WRITE(lunout,*) 'Incoherence between type_trac and config_inca. Model stops. Modify run.def'
          CALL abort_gcm('infotrac_init','Incoherence between type_trac and config_inca',1)
       END IF
#ifndef INCA
       WRITE(lunout,*) 'To run this option you must add cpp key INCA and compile with INCA code'
       CALL abort_gcm('infotrac_init','You must compile with cpp key INCA',1)
#endif
     ELSE IF (type_trac=='repr') THEN
       WRITE(lunout,*) 'You have choosen to couple with REPROBUS chemestry model : type_trac=', type_trac
#ifndef REPROBUS
       WRITE(lunout,*) 'To run this option you must add cpp key REPROBUS and compile with REPRPBUS code'
       CALL abort_gcm('infotrac_init','You must compile with cpp key REPROBUS',1)
#endif
    ELSE IF (type_trac == 'coag') THEN
       WRITE(lunout,*) 'Tracers are treated for COAGULATION tests : type_trac=', type_trac
#ifndef CPP_StratAer
       WRITE(lunout,*) 'To run this option you must add cpp key StratAer and compile with StratAer code'
       CALL abort_gcm('infotrac_init','You must compile with cpp key StratAer',1)
#endif
     ELSE IF (type_trac == 'lmdz') THEN
       WRITE(lunout,*) 'Tracers are treated in LMDZ only : type_trac=', type_trac
     ELSE
       WRITE(lunout,*) 'type_trac=',type_trac,' not possible. Model stops'
       CALL abort_gcm('infotrac_init','bad parameter',1)
     END IF

     ! Test if config_inca is other then none for run without INCA
     IF (type_trac/='inca' .AND. config_inca/='none') THEN
       WRITE(lunout,*) 'config_inca will now be changed to none as you do not couple with INCA model'
       config_inca='none'
     END IF
    ELSE
     type_trac='plnt'  ! planets... May want to dissociate between each later.
    ENDIF ! of IF (planet_type=='earth')

!-----------------------------------------------------------------------
!
! 1) Get the true number of tracers + water vapor/liquid
!    Here true tracers (nqtrue) means declared tracers (only first order)
!
!-----------------------------------------------------------------------
    IF (planet_type=='earth') THEN
     IF (type_trac == 'lmdz' .OR. type_trac == 'repr' .OR. type_trac == 'coag') THEN
       OPEN(90,file='traceur.def',form='formatted',status='old', iostat=ierr)
       IF(ierr.EQ.0) THEN
          WRITE(lunout,*) trim(modname),': Open traceur.def : ok'
          READ(90,*) nqtrue
          WRITE(lunout,*) trim(modname),' nqtrue=',nqtrue
       ELSE 
          WRITE(lunout,*) trim(modname),': Problem in opening traceur.def'
          WRITE(lunout,*) trim(modname),': WARNING using defaut values'
          nqtrue=4 ! Defaut value
       END IF
       ! For Earth, water vapour & liquid tracers are not in the physics
       nbtr=nqtrue-2
     ELSE ! type_trac=inca
       ! The traceur.def file is used to define the number "nqo" of water phases 
       ! present in the simulation. Default : nqo = 2.
       OPEN(90,file='traceur.def',form='formatted',status='old', iostat=ierr)
       IF(ierr.EQ.0) THEN
          WRITE(lunout,*) trim(modname),': Open traceur.def : ok'
          READ(90,*) nqo
       ELSE 
          WRITE(lunout,*) trim(modname),': Using default value for nqo'
          nqo=2
       ENDIF
       IF (nqo /= 2 .AND. nqo /= 3 ) THEN
          WRITE(lunout,*) trim(modname),': nqo=',nqo, ' is not allowded. Only 2 or 3 water phases allowed'
          CALL abort_gcm('infotrac_init','Bad number of water phases',1)
       END IF
       ! nbtr has been read from INCA by init_const_lmdz() in gcm.F 
#ifdef INCA
       CALL Init_chem_inca_trac(nbtr)
#endif
       nqtrue=nbtr+nqo

       ALLOCATE(hadv_inca(nbtr), vadv_inca(nbtr)) 

     END IF   ! type_trac

     IF (nqtrue < 2) THEN
       WRITE(lunout,*) trim(modname),': nqtrue=',nqtrue, ' is not allowded. 2 tracers is the minimum'
       CALL abort_gcm('infotrac_init','Not enough tracers',1)
     END IF

!jyg<
! Transfert number of tracers to Reprobus
!!    IF (type_trac == 'repr') THEN
!!#ifdef REPROBUS
!!       CALL Init_chem_rep_trac(nbtr)
!!#endif
!!    END IF
!>jyg

    ELSE  ! not Earth
       OPEN(90,file='traceur.def',form='formatted',status='old', iostat=ierr)
       IF(ierr.EQ.0) THEN
          WRITE(lunout,*) 'Open traceur.def : ok'
          READ(90,*) nqtrue
       ELSE 
          WRITE(lunout,*) 'Problem in opening traceur.def'
          WRITE(lunout,*) 'ATTENTION using defaut values: nqtrue=1'
          nqtrue=1 ! Defaut value
       END IF
       ! Other planets (for now); we have the same number of tracers
       ! in the dynamics than in the physics
       nbtr=nqtrue
       nqo=0
     
    ENDIF  ! planet_type
!
! Allocate variables depending on nqtrue and nbtr
!
    ALLOCATE(tnom_0(nqtrue), hadv(nqtrue), vadv(nqtrue),tnom_transp(nqtrue))
!
!jyg<
!!    ALLOCATE(conv_flg(nbtr), pbl_flg(nbtr), solsym(nbtr))
!!    conv_flg(:) = 1 ! convection activated for all tracers
!!    pbl_flg(:)  = 1 ! boundary layer activated for all tracers
!>jyg

!-----------------------------------------------------------------------
! 2)     Choix  des schemas d'advection pour l'eau et les traceurs
!
!     iadv = 1    schema  transport type "humidite specifique LMD"
!     iadv = 2    schema   amont
!     iadv = 14   schema  Van-leer + humidite specifique 
!                            Modif F.Codron
!     iadv = 10   schema  Van-leer (retenu pour l'eau vapeur et liquide)
!     iadv = 11   schema  Van-Leer pour hadv et version PPM (Monotone) pour vadv
!     iadv = 12   schema  Frederic Hourdin I
!     iadv = 13   schema  Frederic Hourdin II
!     iadv = 16   schema  PPM Monotone(Collela & Woodward 1984)
!     iadv = 17   schema  PPM Semi Monotone (overshoots autorisés)
!     iadv = 18   schema  PPM Positif Defini (overshoots undershoots autorisés)
!     iadv = 20   schema  Slopes
!     iadv = 30   schema  Prather
!
!        Dans le tableau q(ij,l,iq) : iq = 1  pour l'eau vapeur
!                                     iq = 2  pour l'eau liquide
!       Et eventuellement             iq = 3,nqtot pour les autres traceurs
!
!        iadv(1): choix pour l'eau vap. et  iadv(2) : choix pour l'eau liq.
!------------------------------------------------------------------------
!
!    Get choice of advection schema from file tracer.def or from INCA
!---------------------------------------------------------------------
    IF (planet_type=='earth') THEN
     IF (type_trac == 'lmdz' .OR. type_trac == 'repr' .OR. type_trac == 'coag') THEN
       IF(ierr.EQ.0) THEN
          ! Continue to read tracer.def
          DO iq=1,nqtrue

             write(*,*) 'infotrac 237: iq=',iq
             ! CRisi: ajout du nom du fluide transporteur
             ! mais rester retro compatible
             READ(90,'(I2,X,I2,X,A)',IOSTAT=IOstatus) hadv(iq),vadv(iq),tchaine
             write(lunout,*) 'iq,hadv(iq),vadv(iq)=',iq,hadv(iq),vadv(iq)
             write(lunout,*) 'tchaine=',trim(tchaine)
!             write(*,*) 'infotrac 238: IOstatus=',IOstatus
             if (IOstatus.ne.0) then
                CALL abort_gcm('infotrac_init','Pb dans la lecture de traceur.def',1)
             endif
             ! Y-a-t-il 1 ou 2 noms de traceurs? -> On regarde s'il y a un
             ! espace ou pas au milieu de la chaine.
             continu=.true.
             nouveau_traceurdef=.false.
             iiq=1
             do while (continu)
                if (tchaine(iiq:iiq).eq.' ') then
                  nouveau_traceurdef=.true.
                  continu=.false.
                else if (iiq.lt.LEN_TRIM(tchaine)) then
                  iiq=iiq+1
                else
                  continu=.false.     
                endif
             enddo
             write(*,*) 'iiq,nouveau_traceurdef=',iiq,nouveau_traceurdef
             if (nouveau_traceurdef) then
                write(lunout,*) 'C''est la nouvelle version de traceur.def'
                tnom_0(iq)=tchaine(1:iiq-1)
                tnom_transp(iq)=tchaine(iiq+1:15)
             else
                write(lunout,*) 'C''est l''ancienne version de traceur.def'
                write(lunout,*) 'On suppose que les traceurs sont tous d''air'
                tnom_0(iq)=tchaine
                tnom_transp(iq) = 'air'
             endif
             write(lunout,*) 'tnom_0(iq)=<',trim(tnom_0(iq)),'>'
             write(lunout,*) 'tnom_transp(iq)=<',trim(tnom_transp(iq)),'>'

          END DO ! DO iq=1,nqtrue
          CLOSE(90)  

       ELSE ! Without tracer.def, set default values (for Earth!)
         if ((nqtrue==4).and.(planet_type=="earth")) then
          hadv(1) = 14
          vadv(1) = 14
          tnom_0(1) = 'H2Ov'
          tnom_transp(1) = 'air'
          hadv(2) = 10
          vadv(2) = 10
          tnom_0(2) = 'H2Ol'
          tnom_transp(2) = 'air'
          hadv(3) = 10
          vadv(3) = 10
          tnom_0(3) = 'RN'
          tnom_transp(3) = 'air'
          hadv(4) = 10
          vadv(4) = 10
          tnom_0(4) = 'PB'
          tnom_transp(4) = 'air'
         else
           ! Error message, we need a traceur.def file
           write(lunout,*) trim(modname),&
           ': Cannot set default tracer names!'
           write(lunout,*) trim(modname),' Make a traceur.def file!!!'
           CALL abort_gcm('infotrac_init','Need a traceur.def file!',1)
         endif ! of if (nqtrue==4)
       END IF ! of IF(ierr.EQ.0) 
       
       WRITE(lunout,*) trim(modname),': Valeur de traceur.def :'
       WRITE(lunout,*) trim(modname),': nombre de traceurs ',nqtrue
       DO iq=1,nqtrue
          WRITE(lunout,*) hadv(iq),vadv(iq),tnom_0(iq),tnom_transp(iq)
       END DO

!       IF ( planet_type=='earth') THEN
         !CR: nombre de traceurs de l eau
         IF (tnom_0(3) == 'H2Oi') then
            nqo=3
         ELSE
            nqo=2
         ENDIF
         ! For Earth, water vapour & liquid tracers are not in the physics
         nbtr=nqtrue-nqo
!       ELSE
!         ! Other planets (for now); we have the same number of tracers
!         ! in the dynamics than in the physics
!         nbtr=nqtrue
!       ENDIF

#ifdef CPP_StratAer
       IF (type_trac == 'coag') THEN
         nbtr_bin=0
         nbtr_sulgas=0
         DO iq=1,nqtrue
           IF (tnom_0(iq)(1:3)=='BIN') THEN !check if tracer name contains 'BIN'
             nbtr_bin=nbtr_bin+1
           ENDIF
           IF (tnom_0(iq)(1:3)=='GAS') THEN !check if tracer name contains 'GAS'
             nbtr_sulgas=nbtr_sulgas+1
           ENDIF
         ENDDO
         print*,'nbtr_bin=',nbtr_bin
         print*,'nbtr_sulgas=',nbtr_sulgas
         DO iq=1,nqtrue
           IF (tnom_0(iq)=='GASOCS') THEN
             id_OCS_strat=iq-nqo
           ENDIF
           IF (tnom_0(iq)=='GASSO2') THEN
             id_SO2_strat=iq-nqo
           ENDIF
           IF (tnom_0(iq)=='GASH2SO4') THEN
             id_H2SO4_strat=iq-nqo
           ENDIF
           IF (tnom_0(iq)=='BIN01') THEN
             id_BIN01_strat=iq-nqo
           ENDIF
           IF (tnom_0(iq)=='GASTEST') THEN
             id_TEST_strat=iq-nqo
           ENDIF
         ENDDO
         print*,'id_OCS_strat  =',id_OCS_strat
         print*,'id_SO2_strat  =',id_SO2_strat
         print*,'id_H2SO4_strat=',id_H2SO4_strat
         print*,'id_BIN01_strat=',id_BIN01_strat
       ENDIF
#endif

     ENDIF  ! (type_trac == 'lmdz' .OR. type_trac == 'repr' .OR. type_trac = 'coag')
!jyg<
!
! Transfert number of tracers to Reprobus
    IF (type_trac == 'repr') THEN
#ifdef REPROBUS
       CALL Init_chem_rep_trac(nbtr)
#endif
    END IF
!
! Allocate variables depending on nbtr
!
    ALLOCATE(conv_flg(nbtr), pbl_flg(nbtr), solsym(nbtr))
    conv_flg(:) = 1 ! convection activated for all tracers
    pbl_flg(:)  = 1 ! boundary layer activated for all tracers
!
!!    ELSE  ! type_trac=inca : config_inca='aero' ou 'chem'
!
    IF (type_trac == 'inca') THEN   ! config_inca='aero' ou 'chem'
!>jyg
! le module de chimie fournit les noms des traceurs
! et les schemas d'advection associes. excepte pour ceux lus 
! dans traceur.def 
       IF (ierr .eq. 0) then 
          DO iq=1,nqo

             write(*,*) 'infotrac 237: iq=',iq
             ! CRisi: ajout du nom du fluide transporteur
             ! mais rester retro compatible
             READ(90,'(I2,X,I2,X,A)',IOSTAT=IOstatus) hadv(iq),vadv(iq),tchaine
             write(lunout,*) 'iq,hadv(iq),vadv(iq)=',iq,hadv(iq),vadv(iq)
             write(lunout,*) 'tchaine=',trim(tchaine)
             write(*,*) 'infotrac 238: IOstatus=',IOstatus
             if (IOstatus.ne.0) then
                CALL abort_gcm('infotrac_init','Pb dans la lecture de traceur.def',1)
             endif
             ! Y-a-t-il 1 ou 2 noms de traceurs? -> On regarde s'il y a un
             ! espace ou pas au milieu de la chaine.
             continu=.true.
             nouveau_traceurdef=.false.
             iiq=1
             do while (continu)
                if (tchaine(iiq:iiq).eq.' ') then
                  nouveau_traceurdef=.true.
                  continu=.false.
                else if (iiq.lt.LEN_TRIM(tchaine)) then
                  iiq=iiq+1
                else
                  continu=.false.
                endif
             enddo
             write(*,*) 'iiq,nouveau_traceurdef=',iiq,nouveau_traceurdef
             if (nouveau_traceurdef) then
                write(lunout,*) 'C''est la nouvelle version de traceur.def'
                tnom_0(iq)=tchaine(1:iiq-1)
                tnom_transp(iq)=tchaine(iiq+1:15)
             else
                write(lunout,*) 'C''est l''ancienne version de traceur.def'
                write(lunout,*) 'On suppose que les traceurs sont tous d''air'
                tnom_0(iq)=tchaine
                tnom_transp(iq) = 'air'
             endif
             write(lunout,*) 'tnom_0(iq)=<',trim(tnom_0(iq)),'>'
             write(lunout,*) 'tnom_transp(iq)=<',trim(tnom_transp(iq)),'>'

          END DO !DO iq=1,nqtrue
          CLOSE(90)  
       ELSE  !! if traceur.def doesn't exist 
          tnom_0(1)='H2Ov'
          tnom_transp(1) = 'air'
          tnom_0(2)='H2Ol'
          tnom_transp(2) = 'air'
          hadv(1) = 10
          hadv(2) = 10 
          vadv(1) = 10 
          vadv(2) = 10 
       ENDIF
     
#ifdef INCA
       CALL init_transport( &
            hadv_inca, &
            vadv_inca, &
            conv_flg, &
            pbl_flg,  &
            tracnam)
#endif

!jyg<
       DO iq = nqo+1, nqtrue
          tnom_0(iq)=solsym(iq-nqo)
          tnom_transp(iq) = 'air'
       END DO
!!       DO iq =3,nqtrue
!!          tnom_0(iq)=solsym(iq-2)
!!       END DO
!!       nqo = 2
!>jyg

     END IF ! (type_trac == 'inca')

    ELSE  ! not Earth
       ! Other planets (for now); we have the same number of tracers
       ! in the dynamics than in the physics
       nbtr=nqtrue
       ! NB: Reading a traceur.def with isotopes remains to be done...
       IF(ierr.EQ.0) THEN
          ! Continue to read tracer.def
          DO iq=1,nqtrue
             !READ(90,*) hadv(iq),vadv(iq),tnom_0(iq)
            ! try to be smart when reading traceur.def
            read(90,'(80a)') line ! store the line from traceur.def
            ! assume format is hadv,vadv,tnom_0
            read(line,*,iostat=ierr2) hadv(iq),vadv(iq),tnom_0(iq)
            if (ierr2.ne.0) then
              ! maybe format is tnom0,hadv,vadv
              read(line,*,iostat=ierr3) tnom_0(iq),hadv(iq),vadv(iq)
              if (ierr3.ne.0) then
                ! assume only tnom0 is provided (havd and vad default to 10)
                read(line,*) tnom_0(iq)
                hadv(iq)=10
                vadv(iq)=10
              endif
            endif ! of if(ierr2.ne.0)
            tnom_transp(iq)='air' ! no isotopes... for now...
          END DO ! of DO iq=1,nqtrue
          CLOSE(90)  
       ELSE ! Without tracer.def
          hadv(1) = 10
          vadv(1) = 10
          tnom_0(1) = 'dummy'
          tnom_transp(1)='air'
       END IF
       
       WRITE(lunout,*) trim(modname),': Valeur de traceur.def :'
       WRITE(lunout,*) trim(modname),': nombre de traceurs ',nqtrue
       DO iq=1,nqtrue
          WRITE(lunout,*) hadv(iq),vadv(iq),tnom_0(iq),tnom_transp(iq)
       END DO

    ENDIF  ! planet_type

!-----------------------------------------------------------------------
!
! 3) Verify if advection schema 20 or 30 choosen
!    Calculate total number of tracers needed: nqtot
!    Allocate variables depending on total number of tracers
!-----------------------------------------------------------------------
    new_iq=0
    DO iq=1,nqtrue
       ! Add tracers for certain advection schema
       IF (hadv(iq)<20 .AND. vadv(iq)<20 ) THEN
          new_iq=new_iq+1  ! no tracers added
       ELSE IF (hadv(iq)==20 .AND. vadv(iq)==20 ) THEN
          new_iq=new_iq+4  ! 3 tracers added
       ELSE IF (hadv(iq)==30 .AND. vadv(iq)==30 ) THEN
          new_iq=new_iq+10 ! 9 tracers added
       ELSE
          WRITE(lunout,*) trim(modname),': This choice of advection schema is not available',iq,hadv(iq),vadv(iq)
          CALL abort_gcm('infotrac_init','Bad choice of advection schema - 1',1)
       END IF
    END DO
    
    IF (new_iq /= nqtrue) THEN
       ! The choice of advection schema imposes more tracers
       ! Assigne total number of tracers
       nqtot = new_iq

       WRITE(lunout,*) trim(modname),': The choice of advection schema for one or more tracers'
       WRITE(lunout,*) 'makes it necessary to add tracers'
       WRITE(lunout,*) trim(modname)//': ',nqtrue,' is the number of true tracers'
       WRITE(lunout,*) trim(modname)//': ',nqtot, ' is the total number of tracers needed'

    ELSE
       ! The true number of tracers is also the total number
       nqtot = nqtrue
    END IF

!
! Allocate variables with total number of tracers, nqtot
!
    ALLOCATE(tname(nqtot), ttext(nqtot))
    ALLOCATE(iadv(nqtot), niadv(nqtot))

!-----------------------------------------------------------------------
!
! 4) Determine iadv, long and short name
!
!-----------------------------------------------------------------------
    new_iq=0
    DO iq=1,nqtrue
       new_iq=new_iq+1

       ! Verify choice of advection schema
       IF (hadv(iq)==vadv(iq)) THEN
          iadv(new_iq)=hadv(iq)
       ELSE IF (hadv(iq)==10 .AND. vadv(iq)==16) THEN
          iadv(new_iq)=11
       ELSE
          WRITE(lunout,*)trim(modname),': This choice of advection schema is not available',iq,hadv(iq),vadv(iq)

          CALL abort_gcm('infotrac_init','Bad choice of advection schema - 2',1)
       END IF
      
       str1=tnom_0(iq)
       tname(new_iq)= tnom_0(iq)
       IF (iadv(new_iq)==0) THEN
          ttext(new_iq)=trim(str1)
       ELSE
          ttext(new_iq)=trim(tnom_0(iq))//descrq(iadv(new_iq))
       END IF

       ! schemas tenant compte des moments d'ordre superieur
       str2=ttext(new_iq)
       IF (iadv(new_iq)==20) THEN
          DO jq=1,3
             new_iq=new_iq+1
             iadv(new_iq)=-20
             ttext(new_iq)=trim(str2)//txts(jq)
             tname(new_iq)=trim(str1)//txts(jq)
          END DO
       ELSE IF (iadv(new_iq)==30) THEN
          DO jq=1,9
             new_iq=new_iq+1
             iadv(new_iq)=-30
             ttext(new_iq)=trim(str2)//txtp(jq)
             tname(new_iq)=trim(str1)//txtp(jq)
          END DO
       END IF
    END DO

!
! Find vector keeping the correspodence between true and total tracers
!
    niadv(:)=0
    iiq=0
    DO iq=1,nqtot
       IF(iadv(iq).GE.0) THEN
          ! True tracer
          iiq=iiq+1
          niadv(iiq)=iq
       ENDIF
    END DO


    WRITE(lunout,*) trim(modname),': Information stored in infotrac :'
    WRITE(lunout,*) trim(modname),': iadv  niadv tname  ttext :'
    DO iq=1,nqtot
       WRITE(lunout,*) iadv(iq),niadv(iq),&
       ' ',trim(tname(iq)),' ',trim(ttext(iq))
    END DO

!
! Test for advection schema. 
! This version of LMDZ only garantees iadv=10 and iadv=14 (14 only for water vapour) .
!
    DO iq=1,nqtot
       IF (iadv(iq)/=10 .AND. iadv(iq)/=14 .AND. iadv(iq)/=0) THEN
          WRITE(lunout,*)trim(modname),' STOP : The option iadv=',iadv(iq),' is not tested in this version of LMDZ'
          CALL abort_gcm('infotrac_init','In this version only iadv=10 and iadv=14 is tested!',1)
       ELSE IF (iadv(iq)==14 .AND. iq/=1) THEN
          WRITE(lunout,*)trim(modname),'STOP : The option iadv=',iadv(iq),' is not tested in this version of LMDZ'
          CALL abort_gcm('infotrac_init','In this version iadv=14 is only permitted for water vapour!',1)
       END IF
    END DO

!-----------------------------------------------------------------------
!
! 5) Determine father/son relations for isotopes and carrying fluid
!
!-----------------------------------------------------------------------

! CRisi: quels sont les traceurs fils et les traceurs pères.
! initialiser tous les tableaux d'indices liés aux traceurs familiaux
! + vérifier que tous les pères sont écrits en premières positions
    ALLOCATE(nqfils(nqtot),nqdesc(nqtot))    
    ALLOCATE(iqfils(nqtot,nqtot))   
    ALLOCATE(iqpere(nqtot))
    nqperes=0
    nqfils(:)=0
    nqdesc(:)=0
    iqfils(:,:)=0
    iqpere(:)=0
    nqdesc_tot=0    
    DO iq=1,nqtot
      if (tnom_transp(iq) == 'air') then
        ! ceci est un traceur père
        WRITE(lunout,*) 'Le traceur',iq,', appele ',trim(tnom_0(iq)),', est un pere'
        nqperes=nqperes+1
        iqpere(iq)=0
      else !if (tnom_transp(iq) == 'air') then
        ! ceci est un fils. Qui est son père?
        WRITE(lunout,*) 'Le traceur',iq,', appele ',trim(tnom_0(iq)),', est un fils'
        continu=.true.
        ipere=1
        do while (continu)           
          if (tnom_transp(iq) == tnom_0(ipere)) then
            ! Son père est ipere
            WRITE(lunout,*) 'Le traceur',iq,'appele ', &
      &          trim(tnom_0(iq)),' est le fils de ',ipere,'appele ',trim(tnom_0(ipere))
            nqfils(ipere)=nqfils(ipere)+1  
            iqfils(nqfils(ipere),ipere)=iq
            iqpere(iq)=ipere         
            continu=.false.
          else !if (tnom_transp(iq) == tnom_0(ipere)) then
            ipere=ipere+1 
            if (ipere.gt.nqtot) then
                WRITE(lunout,*) 'Le traceur',iq,'appele ', &
      &          trim(tnom_0(iq)),', est orpelin.'
                CALL abort_gcm('infotrac_init','Un traceur est orphelin',1)
            endif !if (ipere.gt.nqtot) then
          endif !if (tnom_transp(iq) == tnom_0(ipere)) then
        enddo !do while (continu)
      endif !if (tnom_transp(iq) == 'air') then
    enddo !DO iq=1,nqtot
    WRITE(lunout,*) 'infotrac: nqperes=',nqperes    
    WRITE(lunout,*) 'nqfils=',nqfils
    WRITE(lunout,*) 'iqpere=',iqpere
    WRITE(lunout,*) 'iqfils=',iqfils

! Calculer le nombre de descendants à partir de iqfils et de nbfils
    DO iq=1,nqtot    
      generation=0
      continu=.true.
      ifils=iq
      do while (continu)
        ipere=iqpere(ifils)
        if (ipere.gt.0) then
         nqdesc(ipere)=nqdesc(ipere)+1   
         nqdesc_tot=nqdesc_tot+1      
         iqfils(nqdesc(ipere),ipere)=iq
         ifils=ipere
         generation=generation+1
        else !if (ipere.gt.0) then
         continu=.false.
        endif !if (ipere.gt.0) then
      enddo !do while (continu)    
      WRITE(lunout,*) 'Le traceur ',iq,', appele ',trim(tnom_0(iq)),' est un traceur de generation: ',generation
    enddo !DO iq=1,nqtot
    WRITE(lunout,*) 'infotrac: nqdesc=',nqdesc
    WRITE(lunout,*) 'iqfils=',iqfils
    WRITE(lunout,*) 'nqdesc_tot=',nqdesc_tot

! Interdire autres schémas que 10 pour les traceurs fils, et autres schémas
! que 10 et 14 si des pères ont des fils
    do iq=1,nqtot
      if (iqpere(iq).gt.0) then
        ! ce traceur a un père qui n'est pas l'air
        ! Seul le schéma 10 est autorisé
        if (iadv(iq)/=10) then
           WRITE(lunout,*)trim(modname),' STOP : The option iadv=',iadv(iq),' is not implemented for sons'
          CALL abort_gcm('infotrac_init','Sons should be advected by scheme 10',1)
        endif
        ! Le traceur père ne peut être advecté que par schéma 10 ou 14:
        IF (iadv(iqpere(iq))/=10 .AND. iadv(iqpere(iq))/=14) THEN
          WRITE(lunout,*)trim(modname),' STOP : The option iadv=',iadv(iq),' is not implemented for fathers'
          CALL abort_gcm('infotrac_init','Fathers should be advected by scheme 10 ou 14',1)
        endif !IF (iadv(iqpere(iq))/=10 .AND. iadv(iqpere(iq))/=14) THEN
     endif !if (iqpere(iq).gt.0) the
    enddo !do iq=1,nqtot


! detecter quels sont les traceurs isotopiques parmi des traceurs
    call infotrac_isoinit(tnom_0,nqtrue)
        
!-----------------------------------------------------------------------
! Finalize :
!
    DEALLOCATE(tnom_0, hadv, vadv,tnom_transp)


  END SUBROUTINE infotrac_init

!-----------------------------------------------------------------------

  SUBROUTINE infotrac_isoinit(tnom_0,nqtrue)

#ifdef CPP_IOIPSL
  use IOIPSL
#else
  ! if not using IOIPSL, we still need to use (a local version of) getin
  use ioipsl_getincom
#endif
  implicit none
 
    ! inputs
    INTEGER nqtrue
    CHARACTER(len=*) tnom_0(nqtrue)
    
    ! locals    
    CHARACTER(len=3), DIMENSION(niso_possibles) :: tnom_iso
    INTEGER, ALLOCATABLE,DIMENSION(:,:) :: nb_iso,nb_traciso
    INTEGER, ALLOCATABLE,DIMENSION(:) :: nb_isoind
    INTEGER :: ntraceurs_zone_prec,iq,phase,ixt,iiso,izone
    CHARACTER(len=30) :: tnom_trac
    INCLUDE "iniprint.h"

    tnom_iso=(/'eau','HDO','O18','O17','HTO'/)

    ALLOCATE(nb_iso(niso_possibles,nqo))
    ALLOCATE(nb_isoind(nqo))
    ALLOCATE(nb_traciso(niso_possibles,nqo))
    ALLOCATE(iso_num(nqtot))
    ALLOCATE(iso_indnum(nqtot))
    ALLOCATE(zone_num(nqtot))
    ALLOCATE(phase_num(nqtot))
     
    iso_num(:)=0
    iso_indnum(:)=0
    zone_num(:)=0
    phase_num(:)=0
    indnum_fn_num(:)=0
    use_iso(:)=.false.  
    nb_iso(:,:)=0  
    nb_isoind(:)=0     
    nb_traciso(:,:)=0
    niso=0
    ntraceurs_zone=0  
    ntraceurs_zone_prec=0
    ntraciso=0

    do iq=nqo+1,nqtot
!       write(lunout,*) 'infotrac 569: iq,tnom_0(iq)=',iq,tnom_0(iq)
       do phase=1,nqo    
        do ixt= 1,niso_possibles   
         tnom_trac=trim(tnom_0(phase))//'_'
         tnom_trac=trim(tnom_trac)//trim(tnom_iso(ixt))
!         write(*,*) 'phase,ixt,tnom_trac=',phase,ixt,tnom_trac      
         IF (tnom_0(iq) == tnom_trac) then
!          write(lunout,*) 'Ce traceur est un isotope'
          nb_iso(ixt,phase)=nb_iso(ixt,phase)+1   
          nb_isoind(phase)=nb_isoind(phase)+1   
          iso_num(iq)=ixt
          iso_indnum(iq)=nb_isoind(phase)
          indnum_fn_num(ixt)=iso_indnum(iq)
          phase_num(iq)=phase
!          write(lunout,*) 'iso_num(iq)=',iso_num(iq)
!          write(lunout,*) 'iso_indnum(iq)=',iso_indnum(iq)
!          write(lunout,*) 'indnum_fn_num(ixt)=',indnum_fn_num(ixt)
!          write(lunout,*) 'phase_num(iq)=',phase_num(iq)
          goto 20
         else if (iqpere(iq).gt.0) then          
          if (tnom_0(iqpere(iq)) == tnom_trac) then
!           write(lunout,*) 'Ce traceur est le fils d''un isotope'
           ! c'est un traceur d'isotope
           nb_traciso(ixt,phase)=nb_traciso(ixt,phase)+1
           iso_num(iq)=ixt
           iso_indnum(iq)=indnum_fn_num(ixt)
           zone_num(iq)=nb_traciso(ixt,phase)
           phase_num(iq)=phase
!           write(lunout,*) 'iso_num(iq)=',iso_num(iq)
!           write(lunout,*) 'phase_num(iq)=',phase_num(iq)
!           write(lunout,*) 'zone_num(iq)=',zone_num(iq)
           goto 20
          endif !if (tnom_0(iqpere(iq)) == trim(tnom_0(phase))//trim(tnom_iso(ixt))) then
         endif !IF (tnom_0(iq) == trim(tnom_0(phase))//trim(tnom_iso(ixt))) then
        enddo !do ixt= niso_possibles
       enddo !do phase=1,nqo
  20   continue
      enddo !do iq=1,nqtot

!      write(lunout,*) 'iso_num=',iso_num
!      write(lunout,*) 'iso_indnum=',iso_indnum
!      write(lunout,*) 'zone_num=',zone_num  
!      write(lunout,*) 'phase_num=',phase_num
!      write(lunout,*) 'indnum_fn_num=',indnum_fn_num

      do ixt= 1,niso_possibles  
       
       if (nqo.gt.0) then ! Ehouarn: because tests below only valid if nqo>=1,

        if (nb_iso(ixt,1).eq.1) then
          ! on vérifie que toutes les phases ont le même nombre de
          ! traceurs
          do phase=2,nqo
            if (nb_iso(ixt,phase).ne.nb_iso(ixt,1)) then
!              write(lunout,*) 'ixt,phase,nb_iso=',ixt,phase,nb_iso(ixt,phase)
              CALL abort_gcm('infotrac_init','Phases must have same number of isotopes',1)
            endif
          enddo !do phase=2,nqo

          niso=niso+1
          use_iso(ixt)=.true.
          ntraceurs_zone=nb_traciso(ixt,1)

          ! on vérifie que toutes les phases ont le même nombre de
          ! traceurs
          do phase=2,nqo
            if (nb_traciso(ixt,phase).ne.ntraceurs_zone) then
              write(lunout,*) 'ixt,phase,nb_traciso=',ixt,phase,nb_traciso(ixt,phase)
              write(lunout,*) 'ntraceurs_zone=',ntraceurs_zone
              CALL abort_gcm('infotrac_init','Phases must have same number of tracers',1)
            endif  
          enddo  !do phase=2,nqo
          ! on vérifie que tous les isotopes ont le même nombre de
          ! traceurs
          if (ntraceurs_zone_prec.gt.0) then               
            if (ntraceurs_zone.eq.ntraceurs_zone_prec) then
              ntraceurs_zone_prec=ntraceurs_zone
            else !if (ntraceurs_zone.eq.ntraceurs_zone_prec) then
              write(*,*) 'ntraceurs_zone_prec,ntraceurs_zone=',ntraceurs_zone_prec,ntraceurs_zone   
              CALL abort_gcm('infotrac_init', &
               &'Isotope tracers are not well defined in traceur.def',1)           
            endif !if (ntraceurs_zone.eq.ntraceurs_zone_prec) then
           endif !if (ntraceurs_zone_prec.gt.0) then

        else if (nb_iso(ixt,1).ne.0) then
           WRITE(lunout,*) 'nqo,ixt=',nqo,ixt
           WRITE(lunout,*) 'nb_iso(ixt,1)=',nb_iso(ixt,1)    
           CALL abort_gcm('infotrac_init','Isotopes are not well defined in traceur.def',1)     
        endif   !if (nb_iso(ixt,1).eq.1) then       

     endif ! of if (nqo.gt.0)

    enddo ! do ixt= niso_possibles

    ! dimensions isotopique:
    ntraciso=niso*(ntraceurs_zone+1)
!    WRITE(lunout,*) 'niso=',niso
!    WRITE(lunout,*) 'ntraceurs_zone,ntraciso=',ntraceurs_zone,ntraciso    
 
    ! flags isotopiques:
    if (niso.gt.0) then
        ok_isotopes=.true.
    else
        ok_isotopes=.false.
    endif 
!    WRITE(lunout,*) 'ok_isotopes=',ok_isotopes
 
    if (ok_isotopes) then
        ok_iso_verif=.false.
        call getin('ok_iso_verif',ok_iso_verif)
        ok_init_iso=.false.
        call getin('ok_init_iso',ok_init_iso)
        tnat=(/1.0,155.76e-6,2005.2e-6,0.004/100.,0.0/)
        alpha_ideal=(/1.0,1.01,1.006,1.003,1.0/)
    endif !if (ok_isotopes) then  
!    WRITE(lunout,*) 'ok_iso_verif=',ok_iso_verif
!    WRITE(lunout,*) 'ok_init_iso=',ok_init_iso

    if (ntraceurs_zone.gt.0) then
        ok_isotrac=.true.
    else
        ok_isotrac=.false.
    endif   
!    WRITE(lunout,*) 'ok_isotrac=',ok_isotrac

    ! remplissage du tableau iqiso(ntraciso,phase)
    ALLOCATE(iqiso(ntraciso,nqo))   
    iqiso(:,:)=0     
    do iq=1,nqtot
        if (iso_num(iq).gt.0) then
          ixt=iso_indnum(iq)+zone_num(iq)*niso
          iqiso(ixt,phase_num(iq))=iq
        endif
    enddo
!    WRITE(lunout,*) 'iqiso=',iqiso

    ! replissage du tableau index_trac(ntraceurs_zone,niso)
    ALLOCATE(index_trac(ntraceurs_zone,niso))  
    if (ok_isotrac) then
        do iiso=1,niso
          do izone=1,ntraceurs_zone 
             index_trac(izone,iiso)=iiso+izone*niso
          enddo
        enddo
    else !if (ok_isotrac) then      
        index_trac(:,:)=0.0
    endif !if (ok_isotrac) then
!    write(lunout,*) 'index_trac=',index_trac    

! Finalize :
    DEALLOCATE(nb_iso)

  END SUBROUTINE infotrac_isoinit

!-----------------------------------------------------------------------

! Ehouarn: routine iniadvtrac => from Mars/generic; does essentially the
!          same job as infotrac_init. To clean up and merge at some point...
      subroutine iniadvtrac(nq,numvanle)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine which initializes tracer names and advection schemes
! reads these infos from file 'traceur.def' but uses default values
! if that file is not found.
! Ehouarn Millour. Oct. 2008  (made this LMDZ4-like) for future compatibility 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

!#include "dimensions.h"
!#include "advtrac.h"
!#include "control.h"

! routine arguments:
      INTEGER,INTENT(out) :: nq ! number of tracers
      INTEGER,INTENT(out) :: numvanle

! local variables:
      LOGICAL :: first
      INTEGER :: iq
      INTEGER :: ierr
      CHARACTER(len=3) :: qname

! Look for file traceur.def
      OPEN(90,file='traceur.def',form='formatted',status='old', &
              iostat=ierr)
      IF (ierr.eq.0) THEN
        write(*,*) "iniadvtrac: Reading file traceur.def"
        ! read number of tracers:
        read(90,*,iostat=ierr) nq
        if (ierr.ne.0) then
          write(*,*) "iniadvtrac: error reading number of tracers"
          write(*,*) "   (first line of traceur.def) "
          stop
        endif
        
        ! allocate arrays:
        allocate(iadv(nq))
        allocate(tname(nq))
        
        ! initialize advection schemes to Van-Leer for all tracers
        do iq=1,nq
          iadv(iq)=3 ! Van-Leer
        enddo
        
        do iq=1,nq
        ! minimal version, just read in the tracer names, 1 per line
          read(90,*,iostat=ierr) tname(iq)
          if (ierr.ne.0) then
            write(*,*) 'iniadvtrac: error reading tracer names...'
            stop
          endif
        enddo !of do iq=1,nq
        close(90) ! done reading tracer names, close file
      ENDIF ! of IF (ierr.eq.0)

!  ....  Choix  des shemas d'advection pour l'eau et les traceurs  ...
!  ...................................................................
!
!     iadv = 1    shema  transport type "humidite specifique LMD"  
!     iadv = 2    shema   amont
!     iadv = 3    shema  Van-leer
!     iadv = 4    schema  Van-leer + humidite specifique
!                        Modif F.Codron
! 
!
      DO  iq = 1, nq-1
       IF( iadv(iq).EQ.1 ) PRINT *,' Choix du shema humidite specifique'&
       ,' pour le traceur no ', iq
       IF( iadv(iq).EQ.2 ) PRINT *,' Choix du shema  amont',' pour le'  &
       ,' traceur no ', iq
       IF( iadv(iq).EQ.3 ) PRINT *,' Choix du shema  Van-Leer ',' pour' &
       ,'le traceur no ', iq

       IF( iadv(iq).EQ.4 )  THEN
         PRINT *,' Le shema  Van-Leer + humidite specifique ',          &
       ' est  uniquement pour la vapeur d eau .'
         PRINT *,' Corriger iadv( ',iq, ')  et repasser ! '
         CALL ABORT
       ENDIF

       IF( iadv(iq).LE.0.OR.iadv(iq).GT.4 )   THEN
        PRINT *,' Erreur dans le choix de iadv (nqtot).Corriger et '    &
       ,' repasser car  iadv(iq) = ', iadv(iq)
         CALL ABORT
       ENDIF
      ENDDO

       IF( iadv(nq).EQ.1 ) PRINT *,' Choix du shema humidite '          &
       ,'specifique pour la vapeur d''eau'
       IF( iadv(nq).EQ.2 ) PRINT *,' Choix du shema  amont',' pour la'  &
       ,' vapeur d''eau '
       IF( iadv(nq).EQ.3 ) PRINT *,' Choix du shema  Van-Leer '         &
       ,' pour la vapeur d''eau'
       IF( iadv(nq).EQ.4 ) PRINT *,' Choix du shema  Van-Leer + '       &
       ,' humidite specifique pour la vapeur d''eau'
!
       IF( (iadv(nq).LE.0).OR.(iadv(nq).GT.4) )   THEN
        PRINT *,' Erreur dans le choix de iadv (nqtot).Corriger et '    &
       ,' repasser car  iadv(nqtot) = ', iadv(nqtot)
         CALL ABORT
       ENDIF

      first = .TRUE.
      numvanle = nq + 1
      DO  iq = 1, nq
        IF(((iadv(iq).EQ.3).OR.(iadv(iq).EQ.4)).AND.first ) THEN
          numvanle = iq
          first    = .FALSE. 
        ENDIF
      ENDDO
!
      DO  iq = 1, nq

      IF( (iadv(iq).NE.3.AND.iadv(iq).NE.4).AND.iq.GT.numvanle )  THEN
          PRINT *,' Il y a discontinuite dans le choix du shema de ',   &
          'Van-leer pour les traceurs . Corriger et repasser . '
           CALL ABORT
      ENDIF

      ENDDO
!
      end subroutine iniadvtrac


END MODULE infotrac
