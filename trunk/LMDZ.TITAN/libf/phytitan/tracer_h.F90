! WARNING (JB 27/03/2018)
!
! OpenMP directives in this module are inconsistent:
!   - sizes (nmicro...) are thread private \___ BUT BOTH ARE INITIALIZED IN THE SAME ROUTINE 
!   - indexes array are common             /
!
!  Tracers sizes do not need to be private.
!  In such case, OMP THREADPRIVATE should be removed and initracer2 should be called within an OMP SINGLE statement.
!
MODULE tracer_h
  !! Stores data related to physics tracers.
  !!
  !! The module stores public global variables related to the number of tracers
  !! available in the physics and their kind:
  !!
  !! Currently, tracers can be used either for chemistry process (nchimi) or
  !! microphysics (nmicro).
  !!
  !! The subroutine "initracer2" initializes and performs sanity check of
  !! the tracer definitions given in traceur.def and the required tracers in physics
  !! (based on the run parameters).
  !!
  !! The module provides additional methods: 
  !! 
  !!   - indexOfTracer : search for the index of a tracer in the global table (tracers_h:noms) by name.
  !!   - nameOfTracer  : get the name of tracer from a given index (of the global table).
  !!   - dumpTracers   : print the names of all tracers indexes given in argument.
  !!
  IMPLICIT NONE

  INTEGER, SAVE :: nqtot_p  = 0 !! Total number of physical tracers
  INTEGER, SAVE :: nmicro = 0 !! Number of microphysics tracers.
  INTEGER, SAVE :: nice   = 0 !! Number of microphysics ice tracers (subset of nmicro).
  INTEGER, SAVE :: nchimi = 0 !! Number of chemical (gaz species) tracers.
  !$OMP THREADPRIVATE(nqtot_p,nmicro,nice,nchimi)

  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: chimi_indx  !! Indexes of all chemical species tracers 
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: micro_indx  !! Indexes of all microphysical tracers
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ices_indx   !! Indexes of all ice microphysical tracers

  CHARACTER(len=30), DIMENSION(:), ALLOCATABLE, SAVE :: noms      !! name of the tracer
  REAL, DIMENSION(:), ALLOCATABLE, SAVE              :: mmol      !! mole mass of tracer(g/mol-1)
  REAL, DIMENSION(:), ALLOCATABLE, SAVE              :: rat_mmol  !! molar mass ratio 
  REAL, DIMENSION(:), ALLOCATABLE, SAVE              :: rho_q     !! tracer densities (kg.m-3)
  !$OMP THREADPRIVATE(noms,mmol,rat_mmol,rho_q)



  ! tracer indexes: these are initialized in initracer and should be 0 if the
  !                 corresponding tracer does not exist


  CONTAINS

  SUBROUTINE initracer2(nq,nametrac,talk_to_me)
    !! Initialize tracer names and attributes.
    !!
    !! The method initializes the list of tracer names used in the physics from their
    !! dynamics counterpart.
    !!
    !! In addition, it initializes arrays of indexes for the different sub-processes of the physics:
    !!
    !!   - tracers_h:micro_indxs, the array of tracers indexes used for the microphysics.
    !!   - tracers_h:chimi_indxs, the array of tracers indexes used for the chemistry.
    !!
    !! The method also initializes the molar mass array (tracers_h:mmol) for the chemistry and the
    !! molar mass ratio (tracers_h:rat_mmol).
    !!
    !! @note
    !! Strict checking of chemical species name is performed here if the chemistry is activated
    !! (see callchim variable). All the values of 'cnames' must be found in the tracers names
    !! related to chemistry.
    !! @note
    !! Tests are more permissive for the microphysics and is only based on the mimimum number of
    !! tracers expected. Strict name checking is performed in inimufi.
    USE callkeys_mod
    USE comcstfi_mod, only: mugaz
    USE comchem_h, only: nkim, cnames, cmmol
    IMPLICIT NONE

    INTEGER, INTENT(in)                          :: nq         !! Total number of tracers (fixed at compile time)
    character(len=30), DIMENSION(nq), INTENT(in) :: nametrac   !! name of the tracer from dynamics (from 'traceurs.def')
    LOGICAL, INTENT(in), OPTIONAL                :: talk_to_me !! Enable verbose mode.

    LOGICAL                                      :: verb,found
    CHARACTER(len=30)                            :: str

    INTEGER :: i,j,n

    verb = .true. ; IF (PRESENT(talk_to_me)) verb = talk_to_me

    ! nqtot_p could be used everywhere in the physic :)
    nqtot_p=nq

    IF (.NOT.ALLOCATED(noms)) ALLOCATE(noms(nq))
    noms(:)=nametrac(:)

    IF (.NOT.ALLOCATED(rho_q)) ALLOCATE(rho_q(nq)) ! Defined for all tracers, currently initialized to 0.0
    rho_q(:) = 0.0

    ! Defined for all tracers, (actually) initialized only for chemical tracers
    IF (.NOT.ALLOCATED(mmol)) ALLOCATE(mmol(nq))
    IF (.NOT.ALLOCATED(rat_mmol)) ALLOCATE(rat_mmol(nq))
    mmol(:)  = 0.0
    rat_mmol(:) = 1.0

    ! Compute number of microphysics tracers:
    ! By convention they all have the prefix "mu_" (case sensitive !)
    nmicro = 0
    IF (callmufi) THEN
      DO i=1,nq
        str = noms(i)
        IF (str(1:3) == "mu_") nmicro = nmicro+1
      ENDDO
      ! Checking the expected number of tracers:
      !   no cloud:  4 ; w cloud :  4 + 2 + (1+)
      ! Note that we do not make assumptions on the number of chemical species for clouds, this
      ! will be checked in inimufi.
      IF (callclouds) THEN
        IF (nmicro < 7) THEN
          WRITE(*,'((a),I3,(a))') "initracer2:error: Inconsistent number of microphysical tracers &
            &(expected at least 7 tracers,",nmicro," given)"
          CALL abort_gcm("initracer2", "inconsistent number of tracers", 42)
          STOP
        ENDIF
      ELSE IF (nmicro < 4) THEN
          WRITE(*,'((a),I3,(a))') "initracer2:error: Inconsistent number of microphysical tracers &
            &(expected at least 4 tracers,",nmicro," given)"
          CALL abort_gcm("initracer2", "inconsistent number of tracers", 42)
      ELSE IF  (nmicro > 4) THEN
        WRITE(*,'(a)') "initracer2:info: I was expecting only four tracers, you gave me &
          &more. I'll just pretend nothing happen !"
      ENDIF
      ! microphysics indexes share the same values than original tracname.
      IF (.NOT.ALLOCATED(micro_indx)) ALLOCATE(micro_indx(nmicro))
      j = 1
      DO i=1,nq
        str = noms(i)
        IF (str(1:3) == "mu_") THEN
          micro_indx(j) = i 
          j=j+1
        ENDIF
      ENDDO
    ELSE
      IF (.NOT.ALLOCATED(micro_indx)) ALLOCATE(micro_indx(nmicro))
    ENDIF

    ! Compute number of chemical species:
    !   simply assume that all other tracers ARE chemical species
    nchimi = nqtot_p - nmicro

    ! Titan chemistry requires exactly 44 tracers:
    ! Test should be in callchim condition

    IF (callchim) THEN
      IF (nchimi .NE. nkim) THEN
        WRITE(*,*) "initracer2:error: Inconsistent number of chemical species given (",nkim," expected)"
        CALL abort_gcm("initracer2", "inconsistent number of tracers", 42)
      ENDIF
      IF (.NOT.ALLOCATED(chimi_indx)) ALLOCATE(chimi_indx(nchimi))
      n = 0 ! counter on chimi_indx 
      DO j=1,nkim
        found = .false.
        DO i=1,nq
          IF (TRIM(cnames(j)) == TRIM(noms(i))) THEN
            n = n + 1
            chimi_indx(n) = i
            mmol(i) = cmmol(j)
            rat_mmol(i) = cmmol(j)/mugaz
            found = .true.
            EXIT
          ENDIF
        ENDDO
        IF (.NOT.found) THEN
          WRITE(*,*) "initracer2:error: "//TRIM(cnames(j))//" is missing from tracers definition file."
        ENDIF
      ENDDO
      IF (n .NE. nkim) THEN
        WRITE(*,*) "initracer2:error: Inconsistent number of chemical species given (",nkim," expected)"
        CALL abort_gcm("initracer2", "inconsistent number of tracers", 42)
      ENDIF
    ELSE
      IF (.NOT.ALLOCATED(chimi_indx)) ALLOCATE(chimi_indx(0))
    ENDIF
    IF (verb) THEN
      IF (callmufi.OR.callchim) WRITE(*,*) "===== INITRACER2 SPEAKING ====="
      IF (callmufi) THEN
         WRITE(*,*) "Found ",nmicro, "microphysical tracers"
        call dumpTracers(micro_indx)
        WRITE(*,*) "-------------------------------"
      ENDIF
      IF (callchim) THEN
        WRITE(*,*) "Found ",nchimi, "chemical tracers"
        call dumpTracers(chimi_indx)
        WRITE(*,*) "-------------------------------"
      ENDIF
    ENDIF

  END SUBROUTINE initracer2

  FUNCTION indexOfTracer(name, sensitivity) RESULT(idx)
    !! Get the index of a tracer by name.
    !!
    !! The function searches in the global tracer table (tracer_h:noms)
    !! for the given name and returns the first index matching "name".
    !! 
    !! If no name in the table matches the given one, -1 is returned !
    !!
    !! @warning
    !! initracers must be called before any use of this function.
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)  :: name         !! Name of the tracer to search.
    LOGICAL, OPTIONAL, INTENT(in) :: sensitivity  !! Case sensitivity (true by default).
    INTEGER :: idx                                !! Index of the first tracer matching name or -1 if not found.
    LOGICAL                  :: zsens
    INTEGER                  :: j
    CHARACTER(len=LEN(name)) :: zname
    zsens = .true. ; IF(PRESENT(sensitivity)) zsens = sensitivity
    idx = -1
    IF (.NOT.ALLOCATED(noms)) RETURN
    IF (zsens) THEN
      DO j=1,SIZE(noms) 
        IF (TRIM(noms(j)) == TRIM(name)) THEN
          idx = j ; RETURN 
        ENDIF 
      ENDDO
    ELSE
      zname = to_lower(name)
      DO j=1,SIZE(noms) 
        IF (TRIM(to_lower(noms(j))) == TRIM(zname)) THEN
          idx = j ; RETURN 
        ENDIF 
      ENDDO
    ENDIF

    CONTAINS

    FUNCTION to_lower(istr) RESULT(ostr)
      !! Lower case conversion function.
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: istr
      CHARACTER(len=LEN(istr)) :: ostr
      INTEGER :: i, ic
      ostr = istr
      DO i = 1, LEN_TRIM(istr)
        ic = ICHAR(istr(i:i))
        IF (ic >= 65 .AND. ic < 90) ostr(i:i) = char(ic + 32)
      ENDDO 
    END FUNCTION to_lower
  END FUNCTION indexOfTracer

  FUNCTION nameOfTracer(indx) RESULT(name)
    !! Get the name of a tracer by index.
    !!
    !! The function searches in the global tracer table (tracer_h:noms)
    !! and returns the name of the tracer at given index.
    !!
    !! If the index is out of range an empty string is returned.
    !!
    !! @warning
    !! initracers must be called before any use of this function.
    IMPLICIT NONE
    INTEGER, INTENT(in) :: indx   !! Index of the tracer name to retrieve.
    CHARACTER(len=30)   :: name   !! Name of the tracer at given index.
    name = '' 
    IF (.NOT.ALLOCATED(noms)) RETURN
    IF (indx <= 0 .OR. indx > SIZE(noms)) RETURN
    name = noms(indx)
  END FUNCTION nameOfTracer

  SUBROUTINE dumpTracers(indexes)
    !! Print the names of the given list of tracers indexes.
    INTEGER, DIMENSION(:), INTENT(in) :: indexes
    INTEGER :: i,idx,nt
    CHARACTER(len=:), ALLOCATABLE :: suffix
    
    IF (.NOT.ALLOCATED(noms)) THEN
      WRITE(*,'(a)') "[tracers_h:dump_tracers] warning: 'noms' is not allocated, tracers_h:initracer2 has not be called yet"
      RETURN
    ENDIF
    nt = size(noms)
    WRITE(*,"(a)") "local -> global : name"
    DO i=1,size(indexes)
      idx = indexes(i)
      IF (idx < 1 .OR. idx > nt) THEN
        ! WRITE(*,'((a),I3,(a),I3,(a))') "index out of range (",idx,"/",nt,")"
        CYCLE
      ENDIF
      IF (ANY(chimi_indx == idx)) THEN
        suffix = ' (chimi)'
      ELSE IF (ANY(micro_indx == idx)) THEN
        suffix = ' (micro'
        IF (ALLOCATED(ices_indx)) THEN
          IF (ANY(ices_indx == idx)) suffix=suffix//", ice"
        ENDIF
        suffix=suffix//")"
      ELSE 
        suffix=" ()"
      ENDIF
      WRITE(*,'(I5,(a),I6,(a))') i," -> ",idx ," : "//TRIM(noms(i))//suffix
    ENDDO
  END SUBROUTINE dumpTracers


END MODULE tracer_h

