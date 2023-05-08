MODULE init_print_control_mod

! init_print_control to initialize print_control_mod variables
! not included there because of circular dependecy issues

CONTAINS

  SUBROUTINE init_print_control
  USE print_control_mod, ONLY : set_print_control, &
                                prt_level, lunout, debug 
  USE ioipsl_getin_p_mod, ONLY : getin_p
  USE mod_phys_lmdz_para, ONLY: is_omp_root, is_master
  IMPLICIT NONE

!    INTEGER :: lunout ! default output file identifier (6==screen)
!    INTEGER :: prt_level ! Output level
!    LOGICAL :: debug ! flag to specify if in "debug mode"
    LOGICAL :: opened
    INTEGER :: number
    
    !Config  Key  = prt_level
    !Config  Desc = niveau d'impressions de débogage
    !Config  Def  = 0
    !Config  Help = Niveau d'impression pour le débogage
    !Config         (0 = minimum d'impression)
!    prt_level = 0 ! default set in  print_control_mod)
    CALL getin_p('prt_level',prt_level)

    !Config  Key  = lunout
    !Config  Desc = unite de fichier pour les impressions
    !Config  Def  = 6
    !Config  Help = unite de fichier pour les impressions 
    !Config         (defaut sortie standard = 6)
!    lunout=6 ! default set in  print_control_mod)
    CALL getin_p('lunout', lunout)

    IF (is_omp_root) THEN
      IF (lunout /= 5 .and. lunout /= 6) THEN
         INQUIRE(FILE='lmdz.out_0000',OPENED=opened,NUMBER=number)
         IF (opened) THEN
           lunout=number
         ELSE
           OPEN(UNIT=lunout,FILE='lmdz.out_0000',ACTION='write',  &
                STATUS='unknown',FORM='formatted')
         ENDIF
      ENDIF
    ENDIF

    !Config  Key  = debug
    !Config  Desc = mode debogage
    !Config  Def  = false
    !Config  Help = positionne le mode debogage

!    debug = .FALSE. ! default set in  print_control_mod)
    CALL getin_p('debug',debug)
    
    IF (is_master) THEN
      WRITE(lunout,*)"init_print_control: prt_level=",prt_level
      WRITE(lunout,*)"init_print_control: lunout=",lunout
      WRITE(lunout,*)"init_print_control: debug=",debug      
    ENDIF
    
    CALL set_print_control(lunout,prt_level,debug)

  END SUBROUTINE init_print_control  

END MODULE init_print_control_mod

