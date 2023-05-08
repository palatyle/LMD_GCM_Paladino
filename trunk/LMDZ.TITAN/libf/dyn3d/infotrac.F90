MODULE infotrac

IMPLICIT NONE
! nqtot : total number of tracers and higher order of moment, water vapor and liquid included
  INTEGER, SAVE :: nqtot
! CR: add number of tracers for water (for Earth model only!!)
  INTEGER, SAVE :: nqo

! nbtr : number of tracers not including higher order of moment or water vapor or liquid
!        number of tracers used in the physics
  INTEGER, SAVE :: nbtr

! Name variables
  CHARACTER(len=30), ALLOCATABLE, DIMENSION(:), SAVE :: tname ! tracer short name for restart and diagnostics
  CHARACTER(len=33), ALLOCATABLE, DIMENSION(:), SAVE :: ttext ! tracer long name for diagnostics

! iadv  : index of trasport schema for each tracer
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: iadv

! niadv : vector keeping the coorspondance between all tracers(nqtot) treated in the 
!         dynamic part of the code and the tracers (nbtr+2) used in the physics part of the code. 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE    :: niadv ! equivalent dyn / physique

! conv_flg(it)=0 : convection desactivated for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: conv_flg
! pbl_flg(it)=0  : boundary layer diffusion desactivaded for tracer number it 
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: pbl_flg

  CHARACTER(len=4),SAVE :: type_trac
  CHARACTER(len=8),DIMENSION(:),ALLOCATABLE, SAVE :: solsym
 
CONTAINS

  SUBROUTINE infotrac_init
    USE control_mod
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

    CHARACTER(len=30), ALLOCATABLE, DIMENSION(:) :: tnom_0  ! tracer short name
    CHARACTER(len=3), DIMENSION(30) :: descrq
    CHARACTER(len=1), DIMENSION(3)  :: txts
    CHARACTER(len=2), DIMENSION(9)  :: txtp
    CHARACTER(len=30)               :: str1,str2
  
    INTEGER :: nqtrue  ! number of tracers read from tracer.def, without higer order of moment
    INTEGER :: iq, new_iq, iiq, jq, ierr, ierr2, ierr3
    
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
     IF (type_trac == 'lmdz' .OR. type_trac == 'repr') THEN
       OPEN(90,file='traceur.def',form='formatted',status='old', iostat=ierr)
       IF(ierr.EQ.0) THEN
          WRITE(lunout,*) trim(modname),': Open traceur.def : ok'
          READ(90,*) nqtrue
       ELSE 
          WRITE(lunout,*) trim(modname),': Problem in opening traceur.def'
          WRITE(lunout,*) trim(modname),': WARNING using defaut values'
          nqtrue=4 ! Defaut value
       END IF
       ! For Earth, water vapour & liquid tracers are not in the physics
       nbtr=nqtrue-2
     ELSE ! type_trac=inca
       ! nbtr has been read from INCA by init_cont_lmdz() in gcm.F 
       nqtrue=nbtr+2
     END IF

     IF (nqtrue < 2) THEN
       WRITE(lunout,*) trim(modname),': nqtrue=',nqtrue, ' is not allowded. 2 tracers is the minimum'
       CALL abort_gcm('infotrac_init','Not enough tracers',1)
     END IF

! Transfert number of tracers to Reprobus
     IF (type_trac == 'repr') THEN
#ifdef REPROBUS
       CALL Init_chem_rep_trac(nbtr)
#endif
     END IF

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
     
    ENDIF  ! planet_type
!
! Allocate variables depending on nqtrue and nbtr
!
    ALLOCATE(tnom_0(nqtrue), hadv(nqtrue), vadv(nqtrue))
    ALLOCATE(conv_flg(nbtr), pbl_flg(nbtr), solsym(nbtr))
    conv_flg(:) = 1 ! convection activated for all tracers
    pbl_flg(:)  = 1 ! boundary layer activated for all tracers

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
     IF (type_trac == 'lmdz' .OR. type_trac == 'repr') THEN
       IF(ierr.EQ.0) THEN
          ! Continue to read tracer.def
          DO iq=1,nqtrue
             READ(90,*) hadv(iq),vadv(iq),tnom_0(iq)
          END DO
          CLOSE(90)  
       ELSE ! Without tracer.def, set default values (for Earth!)
         if ((nqtrue==4).and.(planet_type=="earth")) then
          hadv(1) = 14
          vadv(1) = 14
          tnom_0(1) = 'H2Ov'
          hadv(2) = 10
          vadv(2) = 10
          tnom_0(2) = 'H2Ol'
          hadv(3) = 10
          vadv(3) = 10
          tnom_0(3) = 'RN'
          hadv(4) = 10
          vadv(4) = 10
          tnom_0(4) = 'PB'
         else
           ! Error message, we need a traceur.def file
           write(lunout,*) trim(modname),&
           ': Cannot set default tracer names!'
           write(lunout,*) trim(modname),' Make a traceur.def file!!!'
           CALL abort_gcm('infotrac_init','Need a traceur.def file!',1)
         endif ! of if (nqtrue==4)
       END IF
       
!CR: nombre de traceurs de l eau
       if (tnom_0(3) == 'H2Oi') then
          nqo=3
       else
          nqo=2
       endif

       WRITE(lunout,*) trim(modname),': Valeur de traceur.def :'
       WRITE(lunout,*) trim(modname),': nombre de traceurs ',nqtrue
       DO iq=1,nqtrue
          WRITE(lunout,*) hadv(iq),vadv(iq),tnom_0(iq)
       END DO

     ELSE  ! type_trac=inca : config_inca='aero' ou 'chem'
! le module de chimie fournit les noms des traceurs
! et les schemas d'advection associes.
     
#ifdef INCA
       CALL init_transport( &
            hadv, &
            vadv, &
            conv_flg, &
            pbl_flg,  &
            tracnam)
#endif
       tnom_0(1)='H2Ov'
       tnom_0(2)='H2Ol'

       DO iq =3,nqtrue
          tnom_0(iq)=solsym(iq-2)
       END DO
       nqo = 2

     END IF ! type_trac

    ELSE  ! not Earth
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
          END DO ! of DO iq=1,nqtrue
          CLOSE(90)  
       ELSE ! Without tracer.def
          hadv(1) = 10
          vadv(1) = 10
          tnom_0(1) = 'dummy'
       END IF
       
       WRITE(lunout,*) trim(modname),': Valeur de traceur.def :'
       WRITE(lunout,*) trim(modname),': nombre de traceurs ',nqtrue
       DO iq=1,nqtrue
          WRITE(lunout,*) hadv(iq),vadv(iq),tnom_0(iq)
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
! Finalize :
!
    DEALLOCATE(tnom_0, hadv, vadv)


  END SUBROUTINE infotrac_init

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

!       IF( iadv(nq).EQ.1 ) PRINT *,' Choix du shema humidite '          &
!       ,'specifique pour la vapeur d''eau'
!       IF( iadv(nq).EQ.2 ) PRINT *,' Choix du shema  amont',' pour la'  &
!       ,' vapeur d''eau '
!       IF( iadv(nq).EQ.3 ) PRINT *,' Choix du shema  Van-Leer '         &
!       ,' pour la vapeur d''eau'
!       IF( iadv(nq).EQ.4 ) PRINT *,' Choix du shema  Van-Leer + '       &
!       ,' humidite specifique pour la vapeur d''eau'
!
!       IF( (iadv(nq).LE.0).OR.(iadv(nq).GT.4) )   THEN
!        PRINT *,' Erreur dans le choix de iadv (nqtot).Corriger et '    &
!       ,' repasser car  iadv(nqtot) = ', iadv(nqtot)
!         CALL ABORT
!       ENDIF

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
