!
MODULE oasis
!
! This module contains subroutines for initialization, sending and receiving 
! towards the coupler OASIS3. It also contains some parameters for the coupling.
!
! This module should always be compiled. With the coupler OASIS3 available the cpp key
! CPP_COUPLE should be set and the entier of this file will then be compiled. 
! In a forced mode CPP_COUPLE should not be defined and the compilation ends before 
! the CONTAINS, without compiling the subroutines.
!
  USE dimphy 
  USE mod_phys_lmdz_para
  USE write_field_phy

#ifdef CPP_COUPLE
  USE mod_prism_proto
  USE mod_prism_def_partition_proto
  USE mod_prism_get_proto
  USE mod_prism_put_proto
#endif
  
  IMPLICIT NONE
  
  ! Id for fields sent to ocean
  INTEGER, PARAMETER :: ids_tauxxu = 1
  INTEGER, PARAMETER :: ids_tauyyu = 2
  INTEGER, PARAMETER :: ids_tauzzu = 3
  INTEGER, PARAMETER :: ids_tauxxv = 4
  INTEGER, PARAMETER :: ids_tauyyv = 5
  INTEGER, PARAMETER :: ids_tauzzv = 6
  INTEGER, PARAMETER :: ids_windsp = 7
  INTEGER, PARAMETER :: ids_shfice = 8
  INTEGER, PARAMETER :: ids_shfoce = 9
  INTEGER, PARAMETER :: ids_shftot = 10
  INTEGER, PARAMETER :: ids_nsfice = 11
  INTEGER, PARAMETER :: ids_nsfoce = 12
  INTEGER, PARAMETER :: ids_nsftot = 13
  INTEGER, PARAMETER :: ids_dflxdt = 14
  INTEGER, PARAMETER :: ids_totrai = 15
  INTEGER, PARAMETER :: ids_totsno = 16
  INTEGER, PARAMETER :: ids_toteva = 17
  INTEGER, PARAMETER :: ids_icevap = 18
  INTEGER, PARAMETER :: ids_ocevap = 19
  INTEGER, PARAMETER :: ids_calvin = 20
  INTEGER, PARAMETER :: ids_liqrun = 21
  INTEGER, PARAMETER :: ids_runcoa = 22
  INTEGER, PARAMETER :: ids_rivflu = 23
  INTEGER, PARAMETER :: ids_atmco2 = 24
  INTEGER, PARAMETER :: ids_taumod = 25
  INTEGER, PARAMETER :: maxsend    = 25  ! Maximum number of fields to send
  
  ! Id for fields received from ocean
  INTEGER, PARAMETER :: idr_sisutw = 1
  INTEGER, PARAMETER :: idr_icecov = 2
  INTEGER, PARAMETER :: idr_icealw = 3
  INTEGER, PARAMETER :: idr_icetem = 4
  INTEGER, PARAMETER :: idr_curenx = 5
  INTEGER, PARAMETER :: idr_cureny = 6
  INTEGER, PARAMETER :: idr_curenz = 7
  INTEGER, PARAMETER :: idr_oceco2 = 8
  INTEGER, PARAMETER :: maxrecv    = 8  ! Maximum number of fields to receive
  

  TYPE, PUBLIC ::   FLD_CPL            ! Type for coupling field information
     CHARACTER(len = 8) ::   name      ! Name of the coupling field   
     LOGICAL            ::   action    ! To be exchanged or not
     INTEGER            ::   nid       ! Id of the field
  END TYPE FLD_CPL

  TYPE(FLD_CPL), DIMENSION(maxsend), SAVE, PUBLIC :: infosend   ! Information for sending coupling fields
  TYPE(FLD_CPL), DIMENSION(maxrecv), SAVE, PUBLIC :: inforecv   ! Information for receiving coupling fields
  
  LOGICAL,SAVE :: cpl_current
!$OMP THREADPRIVATE(cpl_current)

#ifdef CPP_COUPLE

CONTAINS

  SUBROUTINE inicma
!************************************************************************************
!**** *INICMA*  - Initialize coupled mode communication for atmosphere
!                 and exchange some initial information with Oasis
!
!     Rewrite to take the PRISM/psmile library into account
!     LF 09/2003
!
    USE IOIPSL
    USE surface_data, ONLY : version_ocean
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl

    INCLUDE "dimensions.h"
    INCLUDE "iniprint.h"

! Local variables
!************************************************************************************
    INTEGER                            :: comp_id
    INTEGER                            :: ierror, il_commlocal
    INTEGER                            :: il_part_id
    INTEGER, DIMENSION(3)              :: ig_paral
    INTEGER, DIMENSION(2)              :: il_var_nodims
    INTEGER, DIMENSION(4)              :: il_var_actual_shape
    INTEGER                            :: il_var_type
    INTEGER                            :: jf
    CHARACTER (len = 6)                :: clmodnam
    CHARACTER (len = 20)               :: modname = 'inicma'
    CHARACTER (len = 80)               :: abort_message 
    LOGICAL                            :: cpl_current_omp

!*    1. Initializations
!        ---------------
!************************************************************************************
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' ROUTINE INICMA'
    WRITE(lunout,*) ' **************'
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' '

!
! Define the model name
!
    clmodnam = 'lmdz.x'       ! as in $NBMODEL in Cpl/Nam/namcouple.tmp


!************************************************************************************
! Define if coupling ocean currents or not
!************************************************************************************
!$OMP MASTER
    cpl_current_omp = .FALSE.
    CALL getin('cpl_current', cpl_current_omp)
!$OMP END MASTER
!$OMP BARRIER
    cpl_current = cpl_current_omp
    WRITE(lunout,*) 'Couple ocean currents, cpl_current = ',cpl_current 

!************************************************************************************
! Define coupling variables
!************************************************************************************

! Atmospheric variables to send

!$OMP MASTER
    infosend(:)%action = .FALSE.

    infosend(ids_tauxxu)%action = .TRUE. ; infosend(ids_tauxxu)%name = 'COTAUXXU'
    infosend(ids_tauyyu)%action = .TRUE. ; infosend(ids_tauyyu)%name = 'COTAUYYU'
    infosend(ids_tauzzu)%action = .TRUE. ; infosend(ids_tauzzu)%name = 'COTAUZZU'
    infosend(ids_tauxxv)%action = .TRUE. ; infosend(ids_tauxxv)%name = 'COTAUXXV'
    infosend(ids_tauyyv)%action = .TRUE. ; infosend(ids_tauyyv)%name = 'COTAUYYV'
    infosend(ids_tauzzv)%action = .TRUE. ; infosend(ids_tauzzv)%name = 'COTAUZZV'
    infosend(ids_windsp)%action = .TRUE. ; infosend(ids_windsp)%name = 'COWINDSP'
    infosend(ids_shfice)%action = .TRUE. ; infosend(ids_shfice)%name = 'COSHFICE'
    infosend(ids_nsfice)%action = .TRUE. ; infosend(ids_nsfice)%name = 'CONSFICE'
    infosend(ids_dflxdt)%action = .TRUE. ; infosend(ids_dflxdt)%name = 'CODFLXDT'
    infosend(ids_calvin)%action = .TRUE. ; infosend(ids_calvin)%name = 'COCALVIN'
    
    IF (version_ocean=='nemo') THEN
        infosend(ids_shftot)%action = .TRUE. ; infosend(ids_shftot)%name = 'COQSRMIX'
        infosend(ids_nsftot)%action = .TRUE. ; infosend(ids_nsftot)%name = 'COQNSMIX'
        infosend(ids_totrai)%action = .TRUE. ; infosend(ids_totrai)%name = 'COTOTRAI'
        infosend(ids_totsno)%action = .TRUE. ; infosend(ids_totsno)%name = 'COTOTSNO'
        infosend(ids_toteva)%action = .TRUE. ; infosend(ids_toteva)%name = 'COTOTEVA'
        infosend(ids_icevap)%action = .TRUE. ; infosend(ids_icevap)%name = 'COICEVAP'
        infosend(ids_liqrun)%action = .TRUE. ; infosend(ids_liqrun)%name = 'COLIQRUN'
        infosend(ids_taumod)%action = .TRUE. ; infosend(ids_taumod)%name = 'COTAUMOD'
        IF (carbon_cycle_cpl) THEN
            infosend(ids_atmco2)%action = .TRUE. ; infosend(ids_atmco2)%name = 'COATMCO2'
        ENDIF
        
    ELSE IF (version_ocean=='opa8') THEN
        infosend(ids_shfoce)%action = .TRUE. ; infosend(ids_shfoce)%name = 'COSHFOCE'
        infosend(ids_nsfoce)%action = .TRUE. ; infosend(ids_nsfoce)%name = 'CONSFOCE'
        infosend(ids_icevap)%action = .TRUE. ; infosend(ids_icevap)%name = 'COTFSICE'
        infosend(ids_ocevap)%action = .TRUE. ; infosend(ids_ocevap)%name = 'COTFSOCE'
        infosend(ids_totrai)%action = .TRUE. ; infosend(ids_totrai)%name = 'COTOLPSU'
        infosend(ids_totsno)%action = .TRUE. ; infosend(ids_totsno)%name = 'COTOSPSU'
        infosend(ids_runcoa)%action = .TRUE. ; infosend(ids_runcoa)%name = 'CORUNCOA'
        infosend(ids_rivflu)%action = .TRUE. ; infosend(ids_rivflu)%name = 'CORIVFLU'
   ENDIF
        
! Oceanic variables to receive

   inforecv(:)%action = .FALSE.

   inforecv(idr_sisutw)%action = .TRUE. ; inforecv(idr_sisutw)%name = 'SISUTESW'
   inforecv(idr_icecov)%action = .TRUE. ; inforecv(idr_icecov)%name = 'SIICECOV'
   inforecv(idr_icealw)%action = .TRUE. ; inforecv(idr_icealw)%name = 'SIICEALW'
   inforecv(idr_icetem)%action = .TRUE. ; inforecv(idr_icetem)%name = 'SIICTEMW'
   
   IF (cpl_current ) THEN
       inforecv(idr_curenx)%action = .TRUE. ; inforecv(idr_curenx)%name = 'CURRENTX'
       inforecv(idr_cureny)%action = .TRUE. ; inforecv(idr_cureny)%name = 'CURRENTY'
       inforecv(idr_curenz)%action = .TRUE. ; inforecv(idr_curenz)%name = 'CURRENTZ'
   ENDIF

   IF (carbon_cycle_cpl ) THEN
       inforecv(idr_oceco2)%action = .TRUE. ; inforecv(idr_oceco2)%name = 'SICO2FLX'
   ENDIF

!************************************************************************************
! Here we go: psmile initialisation
!************************************************************************************
    IF (is_sequential) THEN
       CALL prism_init_comp_proto (comp_id, clmodnam, ierror)
       
       IF (ierror .NE. PRISM_Ok) THEN
          abort_message=' Probleme init dans prism_init_comp '
          CALL abort_gcm(modname,abort_message,1)
       ELSE
          WRITE(lunout,*) 'inicma : init psmile ok '
       ENDIF
    ENDIF

    CALL prism_get_localcomm_proto (il_commlocal, ierror)
!************************************************************************************
! Domain decomposition
!************************************************************************************
    ig_paral(1) = 1                            ! apple partition for //
    ig_paral(2) = (jj_begin-1)*iim+ii_begin-1  ! offset
    ig_paral(3) = (jj_end*iim+ii_end) - (jj_begin*iim+ii_begin) + 1

    IF (mpi_rank==mpi_size-1) ig_paral(3)=ig_paral(3)+iim-1
    WRITE(lunout,*) mpi_rank,'ig_paral--->',ig_paral(2),ig_paral(3)
    
    ierror=PRISM_Ok
    CALL prism_def_partition_proto (il_part_id, ig_paral, ierror)

    IF (ierror .NE. PRISM_Ok) THEN
       abort_message=' Probleme dans prism_def_partition '
       CALL abort_gcm(modname,abort_message,1)
    ELSE
       WRITE(lunout,*) 'inicma : decomposition domaine psmile ok '
    ENDIF

    il_var_nodims(1) = 2
    il_var_nodims(2) = 1

    il_var_actual_shape(1) = 1
    il_var_actual_shape(2) = iim
    il_var_actual_shape(3) = 1
    il_var_actual_shape(4) = jjm+1
   
    il_var_type = PRISM_Real

!************************************************************************************
! Oceanic Fields to receive
! Loop over all possible variables
!************************************************************************************
    DO jf=1, maxrecv
       IF (inforecv(jf)%action) THEN
          CALL prism_def_var_proto(inforecv(jf)%nid, inforecv(jf)%name, il_part_id, &
               il_var_nodims, PRISM_In, il_var_actual_shape, il_var_type, &
               ierror)
          IF (ierror .NE. PRISM_Ok) THEN
             WRITE(lunout,*) 'inicma : Problem with prism_def_var_proto for field : ',&
                  inforecv(jf)%name
             abort_message=' Problem in call to prism_def_var_proto for fields to receive'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
    END DO
    
!************************************************************************************
! Atmospheric Fields to send
! Loop over all possible variables
!************************************************************************************
    DO jf=1,maxsend
       IF (infosend(jf)%action) THEN
          CALL prism_def_var_proto(infosend(jf)%nid, infosend(jf)%name, il_part_id, &
               il_var_nodims, PRISM_Out, il_var_actual_shape, il_var_type, &
               ierror)
          IF (ierror .NE. PRISM_Ok) THEN
             WRITE(lunout,*) 'inicma : Problem with prism_def_var_proto for field : ',&
                  infosend(jf)%name
             abort_message=' Problem in call to prism_def_var_proto for fields to send'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
    END DO
    
!************************************************************************************
! End definition
!************************************************************************************
    CALL prism_enddef_proto(ierror)
    IF (ierror .NE. PRISM_Ok) THEN
       abort_message=' Problem in call to prism_endef_proto'
       CALL abort_gcm(modname,abort_message,1)
    ELSE
       WRITE(lunout,*) 'inicma : endef psmile ok '
    ENDIF

!$OMP END MASTER
    
  END SUBROUTINE inicma

!
!************************************************************************************
!

  SUBROUTINE fromcpl(ktime, tab_get)
! ======================================================================
! L. Fairhead (09/2003) adapted From L.Z.X Li: this subroutine reads the SST 
! and Sea-Ice provided by the coupler. Adaptation to psmile library
!======================================================================
!
    INCLUDE "dimensions.h"
    INCLUDE "iniprint.h"
! Input arguments
!************************************************************************************
    INTEGER, INTENT(IN)                               ::  ktime

! Output arguments
!************************************************************************************
    REAL, DIMENSION(iim, jj_nb,maxrecv), INTENT(OUT) :: tab_get

! Local variables
!************************************************************************************
    INTEGER                       :: ierror, i
    INTEGER                       :: istart,iend
    CHARACTER (len = 20)          :: modname = 'fromcpl'
    CHARACTER (len = 80)          :: abort_message 
    REAL, DIMENSION(iim*jj_nb)    :: field

!************************************************************************************
    WRITE (lunout,*) ' '
    WRITE (lunout,*) 'Fromcpl: Reading fields from CPL, ktime=',ktime
    WRITE (lunout,*) ' '
    
    istart=ii_begin
    IF (is_south_pole) THEN
       iend=(jj_end-jj_begin)*iim+iim
    ELSE
       iend=(jj_end-jj_begin)*iim+ii_end
    ENDIF
    
    DO i = 1, maxrecv
      IF (inforecv(i)%action) THEN
          field(:) = -99999.
          CALL prism_get_proto(inforecv(i)%nid, ktime, field(istart:iend), ierror)
          tab_get(:,:,i) = RESHAPE(field(:),(/iim,jj_nb/))
       
          IF (ierror .NE. PRISM_Ok .AND. ierror.NE.PRISM_Recvd .AND. &
             ierror.NE.PRISM_FromRest &
             .AND. ierror.NE.PRISM_Input .AND. ierror.NE.PRISM_RecvOut &
             .AND. ierror.NE.PRISM_FromRestOut) THEN
              WRITE (lunout,*)  'Error with receiving filed : ', inforecv(i)%name, ktime   
              abort_message=' Problem in prism_get_proto '
              CALL abort_gcm(modname,abort_message,1)
          ENDIF
      ENDIF
    END DO
    
    
  END SUBROUTINE fromcpl

!
!************************************************************************************
! 

  SUBROUTINE intocpl(ktime, last, tab_put) 
! ======================================================================
! L. Fairhead (09/2003) adapted From L.Z.X Li: this subroutine provides the 
! atmospheric coupling fields to the coupler with the psmile library.
! IF last time step, writes output fields to binary files.
! ======================================================================
!
! 
    INCLUDE "dimensions.h"
    INCLUDE "iniprint.h"
! Input arguments
!************************************************************************************
    INTEGER, INTENT(IN)                              :: ktime
    LOGICAL, INTENT(IN)                              :: last
    REAL, DIMENSION(iim, jj_nb, maxsend), INTENT(IN) :: tab_put

! Local variables
!************************************************************************************
    LOGICAL                          :: checkout
    INTEGER                          :: istart,iend
    INTEGER                          :: wstart,wend
    INTEGER                          :: ierror, i
    REAL, DIMENSION(iim*jj_nb)       :: field
    CHARACTER (len = 20),PARAMETER   :: modname = 'intocpl'
    CHARACTER (len = 80)             :: abort_message 

!************************************************************************************
    checkout=.FALSE.

    WRITE(lunout,*) ' '
    WRITE(lunout,*) 'Intocpl: sending fields to CPL, ktime= ', ktime
    WRITE(lunout,*) 'last = ', last
    WRITE(lunout,*)


    istart=ii_begin
    IF (is_south_pole) THEN
       iend=(jj_end-jj_begin)*iim+iim
    ELSE
       iend=(jj_end-jj_begin)*iim+ii_end
    ENDIF
    
    IF (checkout) THEN   
       wstart=istart
       wend=iend
       IF (is_north_pole) wstart=istart+iim-1
       IF (is_south_pole) wend=iend-iim+1
       
       DO i = 1, maxsend
          IF (infosend(i)%action) THEN
             field = RESHAPE(tab_put(:,:,i),(/iim*jj_nb/))
             CALL writefield_phy(infosend(i)%name,field(wstart:wend),1)
          END IF
       END DO
    END IF

!************************************************************************************
! PRISM_PUT
!************************************************************************************

    DO i = 1, maxsend
      IF (infosend(i)%action) THEN
          field = RESHAPE(tab_put(:,:,i),(/iim*jj_nb/))
          CALL prism_put_proto(infosend(i)%nid, ktime, field(istart:iend), ierror)
          
          IF (ierror .NE. PRISM_Ok .AND. ierror.NE.PRISM_Sent .AND. ierror.NE.PRISM_ToRest &
             .AND. ierror.NE.PRISM_LocTrans .AND. ierror.NE.PRISM_Output .AND. &
             ierror.NE.PRISM_SentOut .AND. ierror.NE.PRISM_ToRestOut) THEN
              WRITE (lunout,*) 'Error with sending field :', infosend(i)%name, ktime   
              abort_message=' Problem in prism_put_proto '
              CALL abort_gcm(modname,abort_message,1)
          ENDIF
      ENDIF
    END DO
   
!************************************************************************************
! Finalize PSMILE for the case is_sequential, if parallel finalization is done 
! from Finalize_parallel in dyn3dpar/parallel.F90
!************************************************************************************

    IF (last) THEN
       IF (is_sequential) THEN 
          CALL prism_terminate_proto(ierror)
          IF (ierror .NE. PRISM_Ok) THEN
             abort_message=' Problem in prism_terminate_proto '
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
    ENDIF
    
    
  END SUBROUTINE intocpl

#endif
  
END MODULE oasis
