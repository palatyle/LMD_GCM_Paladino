!
!$Id: mod_phys_lmdz_omp_data.F90 1403 2010-07-01 09:02:53Z fairhead $
!
MODULE mod_phys_lmdz_omp_data

  INTEGER,SAVE :: omp_size
  INTEGER,SAVE :: omp_rank
  LOGICAL,SAVE :: is_omp_root
  LOGICAL,SAVE :: is_using_omp
  
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_nb
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_begin
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_end    
  
  INTEGER,SAVE :: klon_omp
  INTEGER,SAVE :: klon_omp_begin
  INTEGER,SAVE :: klon_omp_end
!$OMP  THREADPRIVATE(omp_rank,klon_omp,is_omp_root,klon_omp_begin,klon_omp_end)

CONTAINS
  
  SUBROUTINE Init_phys_lmdz_omp_data(klon_mpi)
    USE dimphy
    IMPLICIT NONE
    INTEGER, INTENT(in) :: klon_mpi

    INTEGER :: i

    CHARACTER (LEN=20) :: modname='Init_phys_lmdz_omp_data'
    CHARACTER (LEN=80) :: abort_message


#ifdef CPP_OMP    
    INTEGER :: OMP_GET_NUM_THREADS
    EXTERNAL OMP_GET_NUM_THREADS
    INTEGER :: OMP_GET_THREAD_NUM
    EXTERNAL OMP_GET_THREAD_NUM
#endif  

#ifdef CPP_OMP
!$OMP MASTER
        is_using_omp=.TRUE.
        omp_size=OMP_GET_NUM_THREADS()
!$OMP END MASTER
        omp_rank=OMP_GET_THREAD_NUM()    
#else    
    is_using_omp=.FALSE.
    omp_size=1
    omp_rank=0
#endif

   is_omp_root=.FALSE.
!$OMP MASTER
   IF (omp_rank==0) THEN
     is_omp_root=.TRUE.
   ELSE
     abort_message = 'ANORMAL : OMP_MASTER /= 0'
     CALL abort_gcm (modname,abort_message,1)
   ENDIF
!$OMP END MASTER


!$OMP MASTER 
    ALLOCATE(klon_omp_para_nb(0:omp_size-1))
    ALLOCATE(klon_omp_para_begin(0:omp_size-1))
    ALLOCATE(klon_omp_para_end(0:omp_size-1))
    
    DO i=0,omp_size-1
      klon_omp_para_nb(i)=klon_mpi/omp_size
      IF (i<MOD(klon_mpi,omp_size)) klon_omp_para_nb(i)=klon_omp_para_nb(i)+1
    ENDDO
    
    klon_omp_para_begin(0) = 1
    klon_omp_para_end(0) = klon_omp_para_nb(0)
    
    DO i=1,omp_size-1
      klon_omp_para_begin(i)=klon_omp_para_end(i-1)+1
      klon_omp_para_end(i)=klon_omp_para_begin(i)+klon_omp_para_nb(i)-1
    ENDDO
!$OMP END MASTER
!$OMP BARRIER
   
    klon_omp=klon_omp_para_nb(omp_rank)
    klon_omp_begin=klon_omp_para_begin(omp_rank)
    klon_omp_end=klon_omp_para_end(omp_rank)
    
    CALL Print_module_data
    
  END SUBROUTINE Init_phys_lmdz_omp_data

  SUBROUTINE Print_module_data
  IMPLICIT NONE

!$OMP CRITICAL  
  PRINT *,'--------> TASK ',omp_rank
  PRINT *,'omp_size =',omp_size
  PRINT *,'omp_rank =',omp_rank
  PRINT *,'is_omp_root =',is_omp_root
  PRINT *,'klon_omp_para_nb =',klon_omp_para_nb
  PRINT *,'klon_omp_para_begin =',klon_omp_para_begin
  PRINT *,'klon_omp_para_end =',klon_omp_para_end    
  PRINT *,'klon_omp =',klon_omp
  PRINT *,'klon_omp_begin =',klon_omp_begin
  PRINT *,'klon_omp_end =',klon_omp_end    
!$OMP END CRITICAL

  END SUBROUTINE Print_module_data
END MODULE mod_phys_lmdz_omp_data
