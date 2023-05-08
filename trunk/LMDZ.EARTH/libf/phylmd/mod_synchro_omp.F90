MODULE mod_synchro_omp

  LOGICAL,SAVE,ALLOCATABLE :: flag_omp(:)
  
CONTAINS

  SUBROUTINE Init_synchro_omp
  USE mod_phys_lmdz_para 
  IMPLICIT NONE
    
    IF (is_omp_root) THEN
      ALLOCATE(flag_omp(0:omp_size-1))
      flag_omp(:)=.FALSE.
    ENDIF
!$OMP BARRIER

  END SUBROUTINE Init_Synchro_omp
  
  SUBROUTINE Synchro_omp
  USE mod_phys_lmdz_para
  IMPLICIT NONE
  
    flag_omp(omp_rank)=.TRUE.
!$OMP BARRIER
    DO WHILE (.NOT. ALL(flag_omp))
!$OMP BARRIER
    ENDDO
!$OMP BARRIER        
    flag_omp(omp_rank)=.FALSE.
!$OMP BARRIER

   END SUBROUTINE Synchro_omp

END MODULE mod_synchro_omp
