subroutine inimufi(ptimestep)

  use mmp_gcm
  use callkeys_mod, only : callclouds, p_prod, tx_prod, rc_prod, air_rad, eff_gz
  use tracer_h
  use comcstfi_mod, only : g, rad, mugaz
  use datafile_mod

  implicit none


  !============================================================================
  !
  !     Purpose
  !     -------
  !     This routine call mm_initialize which perform the global initialization
  !     for the microphysical YAMMS model.
  !     It also performs some sanity checks on microphysical tracer names. 
  !
  !     Authors
  !     -------
  !     J. Vatant d'Ollone, J.Burgalat, S.Lebonnois (09/2017)
  !
  !============================================================================


  !--------------------------
  ! 0. Declarations
  ! -------------------------

  real, intent(in)                   :: ptimestep    ! Timestep (s)

  integer :: i,idx
  character(len=30), dimension(4), parameter :: aernames = &
     (/"mu_m0as              ", "mu_m3as              ",   &
       "mu_m0af              ", "mu_m3af              "/)
  CHARACTER(len=30), DIMENSION(2), PARAMETER :: ccnnames = &
     (/"mu_m0n               ", "mu_m3n               "/)
  logical :: err

  ! PATCH : YAMMS now allows to enable/disable effective g computations:
  mm_use_effg = eff_gz

  !----------------------------------------------------
  ! 1. Call microphysical model global initialization
  ! ---------------------------------------------------

  ! enable log for what it's worth...
  ! mm_log = .true.
  
  call mmp_initialize(ptimestep,p_prod,tx_prod,rc_prod, &
        rad,g,air_rad,mugaz,callclouds,config_mufi)

   ! -------------------------
   ! 2. Check names of tracers
   ! -------------------------
   err = .false.
   DO i=1,size(aernames)
     idx = indexOfTracer(TRIM(aernames(i)),.false.)
     IF (idx <= 0) THEN
       WRITE(*,*) "inimufi: '"//TRIM(aernames(i))//"' not found in tracers table."
       err = .true.
     ENDIF
   ENDDO
   IF (callclouds) THEN
     DO i=1,size(ccnnames)
       idx = indexOfTracer(TRIM(ccnnames(i)),.false.)
       IF (idx <= 0) THEN
         WRITE(*,*) "inimufi: '"//TRIM(ccnnames(i))//"' not found in tracers table."
         err = .true.
       ENDIF
     ENDDO
     ALLOCATE(ices_indx(mm_nesp))
     ices_indx(:) = -1
     DO i=1,mm_nesp 
       idx = indexOfTracer("mu_m3"//TRIM(mm_spcname(i)),.false.)
       IF (idx <= 0) THEN
         WRITE(*,*) "inimufi: 'mu_m3"//TRIM(mm_spcname(i))//"' not found in tracers table."
         err = .true.
       ELSE
         ices_indx(i) = idx
       ENDIF
     ENDDO
     IF (err) THEN
       WRITE(*,*) "You loose.... Try again"
       STOP
     ENDIF
   ENDIF

end subroutine inimufi
