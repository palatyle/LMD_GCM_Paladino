










! $Id: $
MODULE print_control_mod

  INTEGER,SAVE :: lunout=6 ! default output file identifier (6==screen)
  INTEGER,SAVE :: prt_level=0 ! debug output level
  LOGICAL,SAVE :: debug=.FALSE. ! flag to specify if in "debug mode"
!$OMP THREADPRIVATE(lunout,prt_level,debug)

  ! NB: Module variable Initializations done by set_print_control
  !     routine from init_print_control_mod to avoid circular
  !     module dependencies

CONTAINS

  SUBROUTINE set_print_control(lunout_,prt_level_,debug_)
  IMPLICIT NONE
    INTEGER :: lunout_
    INTEGER :: prt_level_
    LOGICAL :: debug_
      
    lunout = lunout_
    prt_level = prt_level_
    debug = debug_
    
  END SUBROUTINE set_print_control

END MODULE print_control_mod
