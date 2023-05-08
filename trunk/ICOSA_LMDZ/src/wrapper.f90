SUBROUTINE initialize_external_physics
  USE interface_icosa_lmdz_mod, ONLY: initialize_physics
  IMPLICIT NONE
  
    CALL initialize_physics  

END SUBROUTINE initialize_external_physics


SUBROUTINE external_physics
  USE interface_icosa_lmdz_mod, ONLY: physics
  IMPLICIT NONE
  
  CALL physics

END SUBROUTINE external_physics
