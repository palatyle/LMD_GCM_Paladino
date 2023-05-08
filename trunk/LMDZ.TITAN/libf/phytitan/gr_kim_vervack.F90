SUBROUTINE gr_kim_vervack

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Purpose : * Calculates the number of layers needed for upper chemistry
  ! -------   based on the GCM vertical grid size, depending on Ptop=ap.
  !           * Generates also the pressure grid at mid-layer for upper levels.
  !           * We use an upper atmosphere profile based on Voyager 1 data 
  !           (Vervack et al, 2004) to have dz=10km between Ptop and 1300km.
  !
  ! Author : Jan Vatant d'Ollone (2017)
  ! ------
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE comchem_h, ONLY: nlaykim_up, preskim, grkim_dz
  USE comvert_mod, ONLY: ap, aps, bp, bps, preff 
  USE datafile_mod, ONLY: datadir

  IMPLICIT NONE

  INCLUDE "dimensions.h"

  REAL*8  :: ptop, zptop_ver, zkim
  REAL*8  :: factz

  REAL*8  :: rver(131),tver(131),ctver(131),pver(131) 

  INTEGER :: i, j, itop, jlay, jold
  
  LOGICAL :: foundver

  LOGICAL :: crop

  ! ------------------------
  ! 1. Read Vervack profile
  ! ------------------------

  ! Check the file exists
  INQUIRE(file=TRIM(datadir)//'/tcp.ver',exist=foundver)
  IF(.NOT.foundver) THEN
    WRITE(*,*)'The file ',TRIM(datadir)//'/tcp.ver'
    WRITE(*,*)'was not found by gr_kim_vervack.F90, exiting'
    WRITE(*,*)'Check that your path to datagcm:',trim(datadir)
    WRITE(*,*)'is correct. You can change it in callphys.def with:'
    WRITE(*,*)'datadir = /absolute/path/to/datagcm'
    WRITE(*,*)'Also check that file tcp.ver exists there.'
  ENDIF

  ! Read file
  OPEN(11,file=TRIM(datadir)//'/tcp.ver',status='old')
  READ(11,*)
  DO i=1,131
     READ(11,*) rver(i),tver(i),ctver(i),pver(i)
     pver(i) = pver(i)*100.0 ! mbar->Pa
  ENDDO
  CLOSE(11)
  
  ! --------------------------------------------------------------
  ! 2. Define ptop as the value aps should have if it wasn't zero
  ! assuming ap(llm)-aps(llm) half pressure thickness of top-layer
  ! --------------------------------------------------------------
  
  ! NB : At the top of the model we are in pure pressure coord. -> ap
  ! ( except for 1D where we have only bp )
  
  IF (jjm.GT.1) THEN
    ptop = 2.0*aps(llm) - ap(llm)
  ELSE
    ptop = preff*(2.0*bps(llm) - bp(llm))
  ENDIF
  
  ! --------------------------------------------
  ! 3. Interpolate Ptop and equivalent altitude 
  ! --------------------------------------------
  
  itop=1
  
  DO i=2,131
    IF ( ptop .LT. pver(i) ) THEN
      itop=i
    ENDIF
  ENDDO
  
  ! Crazy case in a far far away future where GCM top will reach 1300km
  IF ( itop .eq. 131 ) THEN
    WRITE(*,*) " You've reach the bounds of Vervack profile ... "
    WRITE(*,*) " Congrats but it wasn't predicted in 2017 !"
    WRITE(*,*) " I'll stop here ... "
    CALL abort
  ENDIF
    
  factz = ( ptop - pver(itop) ) / ( pver(itop+1) - pver(itop) )
  
  zptop_ver = rver(itop)*(1.-factz) + rver(itop+1)*factz
  
  ! ---------------------------------------------------------
  ! 4. Find zkim max assuming dz=10km and hence nlaykim_up
  ! ---------------------------------------------------------
  
  zkim = zptop_ver
  i=0
  
  DO WHILE ( zkim.lt.rver(131) )

    zkim = zkim + grkim_dz
    i=i+1
  
  ENDDO
  
  ! We want the ceiling at 1300km sharp, so we will either crop
  ! the last layer or remove it and enlarge the "n-1"th
  IF ( zkim .GT. rver(131)+0.5*grkim_dz ) THEN 
    nlaykim_up = i-1
    crop = .FALSE.
  ELSE
    nlaykim_up = i
    crop = .TRUE.
  ENDIF
  
  ! -----------------------------------------------------------------
  ! 5. Calculates preskim grid interpolating back on Vervack profile
  ! -----------------------------------------------------------------

  ALLOCATE(preskim(nlaykim_up))
  
  jold=2

  zkim = zptop_ver + 0.5*grkim_dz ! We want values at mid-layer here !
  
  DO i=1,nlaykim_up
  
    DO j=jold,131
      IF ( zkim .GT. rver(j) ) THEN
        jlay = j
      ENDIF
    ENDDO

    jold = jlay ! keep in memory where we were in the above loop

    ! We want the ceiling at 1300km sharp, we readjust the size of last layer
    IF (i.eq.nlaykim_up) THEN
      zkim = 0.5* ( zkim - 0.5*grkim_dz +  rver(131) )
    ENDIF
    
    factz = ( zkim - rver(jlay) ) / ( rver(jlay+1) - rver(jlay) )
    preskim(i) = pver(jlay)**(1.-factz) * pver(jlay+1)**factz 
    
    zkim = zkim + grkim_dz 

  ENDDO

  
END SUBROUTINE gr_kim_vervack
