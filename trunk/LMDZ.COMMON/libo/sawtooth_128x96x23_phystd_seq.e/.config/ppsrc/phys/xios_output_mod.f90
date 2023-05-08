










MODULE xios_output_mod

 IMPLICIT NONE
 
 INTEGER,PRIVATE,SAVE :: time_it=0 ! store number of iterations with calls to XIOS since start
! does not need to be threadprivate; managed by omp master

 CHARACTER(LEN=*), PARAMETER :: context_id= "LMDZ" ! same as in context_lmdz_physics.xml
 

END MODULE xios_output_mod
