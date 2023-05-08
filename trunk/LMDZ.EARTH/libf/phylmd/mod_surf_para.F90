MODULE mod_surf_para
  IMPLICIT NONE
  
  INTERFACE gather_surf
    MODULE PROCEDURE gather_surf_i,gather_surf_r
  END INTERFACE gather_surf
  
  INTERFACE gather_surf_omp
    MODULE PROCEDURE gather_surf_omp_i,gather_surf_omp_r
  END INTERFACE gather_surf_omp

  INTERFACE gather_surf_mpi
    MODULE PROCEDURE gather_surf_mpi_i,gather_surf_mpi_r
  END INTERFACE gather_surf_mpi

  INTERFACE scatter_surf
    MODULE PROCEDURE scatter_surf_i,scatter_surf_r
  END INTERFACE scatter_surf
  
  INTERFACE scatter_surf_omp
    MODULE PROCEDURE scatter_surf_omp_i,scatter_surf_omp_r
  END INTERFACE scatter_surf_omp

  INTERFACE scatter_surf_mpi
    MODULE PROCEDURE scatter_surf_mpi_i,scatter_surf_mpi_r
  END INTERFACE scatter_surf_mpi
  
  
  INTEGER,SAVE             :: knon_omp
  INTEGER,SAVE             :: knon_omp_begin
  INTEGER,SAVE             :: knon_omp_end
!$OMP THREADPRIVATE(knon_omp,knon_omp_begin,knon_omp_end)
  INTEGER,ALLOCATABLE,SAVE :: knon_omp_para(:)
  INTEGER,ALLOCATABLE,SAVE :: knon_omp_begin_para(:)
  INTEGER,ALLOCATABLE,SAVE :: knon_omp_end_para(:)
  
  INTEGER,SAVE             :: knon_mpi
  INTEGER,ALLOCATABLE,SAVE :: knon_mpi_para(:)
  INTEGER,ALLOCATABLE,SAVE :: knon_mpi_begin_para(:)
  INTEGER,ALLOCATABLE,SAVE :: knon_mpi_end_para(:)
  
  INTEGER,SAVE             :: knon_glo
  INTEGER,SAVE,ALLOCATABLE :: knon_glo_para(:)  
  INTEGER,ALLOCATABLE,SAVE :: knon_glo_begin_para(:)
  INTEGER,ALLOCATABLE,SAVE :: knon_glo_end_para(:)
  
  
CONTAINS

  SUBROUTINE Init_surf_para(knon)
  USE mod_phys_lmdz_para, mpi_rank_root=>mpi_root
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
    INTEGER :: knon
    INTEGER :: i,ierr
    
    knon_omp=knon
    IF (is_omp_root) THEN
      ALLOCATE(knon_omp_para(0:omp_size-1))
      ALLOCATE(knon_omp_begin_para(0:omp_size-1))
      ALLOCATE(knon_omp_end_para(0:omp_size-1))
    ENDIF
!$OMP BARRIER
    knon_omp_para(omp_rank)=knon
!$OMP BARRIER
    IF (is_omp_root) THEN 
      knon_omp_begin_para(0)=1
      knon_omp_end_para(0)=knon_omp_para(0)
      DO i=1,omp_size-1
        knon_omp_begin_para(i)=knon_omp_end_para(i-1)+1
        knon_omp_end_para(i)=knon_omp_begin_para(i)+knon_omp_para(i)-1
      ENDDO
    ENDIF 
!$OMP BARRIER
    knon_omp_begin=knon_omp_begin_para(omp_rank)
    knon_omp_end=knon_omp_end_para(omp_rank)
!$OMP BARRIER    
    IF (is_omp_root) THEN
      knon_mpi=sum(knon_omp_para)
      ALLOCATE(knon_mpi_para(0:mpi_size-1))
      ALLOCATE(knon_mpi_begin_para(0:mpi_size-1))
      ALLOCATE(knon_mpi_end_para(0:mpi_size-1))
      
      ALLOCATE(knon_glo_para(0:mpi_size*omp_size-1))
      ALLOCATE(knon_glo_begin_para(0:mpi_size*omp_size-1))
      ALLOCATE(knon_glo_end_para(0:mpi_size*omp_size-1))
      
      IF (is_using_mpi) THEN
#ifdef CPP_MPI
        CALL MPI_ALLGather(knon_mpi,1,MPI_INTEGER,knon_mpi_para,1,MPI_INTEGER,COMM_LMDZ_PHY,ierr)
        CALL MPI_ALLGather(knon_omp_para,omp_size,MPI_INTEGER,knon_glo_para,omp_size,MPI_INTEGER,COMM_LMDZ_PHY,ierr)
#endif
      ELSE
        knon_mpi_para(:)=knon_mpi
        knon_glo_para(:)=knon_omp_para(:)
      ENDIF     
      
      knon_glo=sum(knon_mpi_para(:))
      
      knon_mpi_begin_para(0)=1
      knon_mpi_end_para(0)=knon_mpi_para(0)
      DO i=1,mpi_size-1
        knon_mpi_begin_para(i)=knon_mpi_end_para(i-1)+1
        knon_mpi_end_para(i)=knon_mpi_begin_para(i)+knon_mpi_para(i)-1
      ENDDO
      
      knon_glo_begin_para(0)=1
      knon_glo_end_para(0)=knon_glo_para(0)
      DO i=1,mpi_size*omp_size-1
        knon_glo_begin_para(i)=knon_glo_end_para(i-1)+1
        knon_glo_end_para(i)= knon_glo_begin_para(i)+knon_glo_para(i)-1
      ENDDO
   ENDIF
!$OMP BARRIER

  END SUBROUTINE Init_surf_para

 
  SUBROUTINE Finalize_surf_para
  USE mod_phys_lmdz_para

!$OMP BARRIER   
   IF (is_omp_root) THEN
      DEALLOCATE(knon_omp_para)
      DEALLOCATE(knon_omp_begin_para)
      DEALLOCATE(knon_omp_end_para)
      DEALLOCATE(knon_mpi_para)
      DEALLOCATE(knon_mpi_begin_para)
      DEALLOCATE(knon_mpi_end_para)
      DEALLOCATE(knon_glo_para)  
      DEALLOCATE(knon_glo_begin_para)
      DEALLOCATE(knon_glo_end_para)
    ENDIF
    
  END SUBROUTINE Finalize_surf_para
  
  
  SUBROUTINE gather_surf_i(FieldIn, FieldOut)
  USE mod_phys_lmdz_para
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
    INTEGER :: FieldTmp(knon_mpi)
    
    CALL gather_surf_omp_i(FieldIn,FieldTmp)
    IF (is_omp_root) CALL gather_surf_mpi_i(FieldTmp,FieldOut)
    
  END SUBROUTINE gather_surf_i


  SUBROUTINE gather_surf_omp_i(FieldIn,FieldOut)
  USE mod_phys_lmdz_para
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
  
    INTEGER,SAVE,ALLOCATABLE :: Field_tmp(:)
    
    IF (is_omp_root) ALLOCATE(Field_tmp(knon_mpi))
!$OMP BARRIER
    Field_tmp(knon_omp_begin:knon_omp_end)=FieldIn(:)
!$OMP BARRIER        
    IF (is_omp_root) FieldOut(:)=Field_tmp(:)
!$OMP BARRIER
    IF (is_omp_root) DEALLOCATE(Field_tmp)
    
  END SUBROUTINE  gather_surf_omp_i
  
     
  SUBROUTINE gather_surf_mpi_i(FieldIn,FieldOut)
  USE mod_phys_lmdz_para, mpi_rank_root => mpi_root
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
    INTEGER :: ierr
    
    IF (is_using_mpi) THEN
#ifdef CPP_MPI
      CALL MPI_Gatherv(FieldIn,knon_mpi,MPI_INTEGER,                                &
                       FieldOut,knon_mpi_para,knon_mpi_begin_para(:)-1,MPI_INTEGER, &
                       mpi_rank_root,COMM_LMDZ_PHY,ierr)
#endif
    ELSE
      FieldOut(:)=FieldIn(:)
    ENDIF
  
  END SUBROUTINE gather_surf_mpi_i
  




  SUBROUTINE gather_surf_r(FieldIn, FieldOut)
  USE mod_phys_lmdz_para
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
    REAL :: FieldTmp(knon_mpi)
    
    CALL gather_surf_omp_r(FieldIn,FieldTmp)
    IF (is_omp_root) CALL gather_surf_mpi_r(FieldTmp,FieldOut)
    
  END SUBROUTINE gather_surf_r


  SUBROUTINE gather_surf_omp_r(FieldIn,FieldOut)
  USE mod_phys_lmdz_para
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
  
    REAL,SAVE,ALLOCATABLE :: Field_tmp(:)
    
    IF (is_omp_root) ALLOCATE(Field_tmp(knon_mpi))
!$OMP BARRIER
    Field_tmp(knon_omp_begin:knon_omp_end)=FieldIn(:)
!$OMP BARRIER        
    IF (is_omp_root) FieldOut(:)=Field_tmp(:)
!$OMP BARRIER
    IF (is_omp_root) DEALLOCATE(Field_tmp)
    
  END SUBROUTINE  gather_surf_omp_r
  
     
  SUBROUTINE gather_surf_mpi_r(FieldIn,FieldOut)
  USE mod_phys_lmdz_para, mpi_rank_root => mpi_root
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
    REAL :: ierr
    
    IF (is_using_mpi) THEN
#ifdef CPP_MPI
      CALL MPI_Gatherv(FieldIn,knon_mpi,MPI_REAL_LMDZ,                                 &
                       FieldOut,knon_mpi_para,knon_mpi_begin_para(:)-1,MPI_REAL_LMDZ,  &
                       mpi_rank_root,COMM_LMDZ_PHY,ierr)            
#endif
    ELSE
      FieldOut(:)=FieldIn(:)
    ENDIF
  
  END SUBROUTINE gather_surf_mpi_r




  SUBROUTINE scatter_surf_i(FieldIn, FieldOut)
  USE mod_phys_lmdz_para
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
    INTEGER :: FieldTmp(knon_mpi)
    
    IF (is_omp_root) CALL scatter_surf_mpi_i(FieldIn,FieldTmp)
    CALL scatter_surf_omp_i(FieldTmp,FieldOut)
    
  END SUBROUTINE scatter_surf_i


  SUBROUTINE scatter_surf_omp_i(FieldIn,FieldOut)
  USE mod_phys_lmdz_para
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
  
    INTEGER,SAVE,ALLOCATABLE :: Field_tmp(:)
    
    IF (is_omp_root) ALLOCATE(Field_tmp(knon_mpi))
    IF (is_omp_root) Field_tmp(:)=FieldIn(:)
!$OMP BARRIER        
    FieldOut(:)=Field_tmp(knon_omp_begin:knon_omp_end)
!$OMP BARRIER
    IF (is_omp_root) DEALLOCATE(Field_tmp)
    
  END SUBROUTINE  scatter_surf_omp_i
  
     
  SUBROUTINE scatter_surf_mpi_i(FieldIn,FieldOut)
  USE mod_phys_lmdz_para, mpi_rank_root => mpi_root
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
    INTEGER :: FieldIn(:)
    INTEGER :: FieldOut(:)
    INTEGER :: ierr
    
    IF (is_using_mpi) THEN
#ifdef CPP_MPI
      CALL MPI_Scatterv(FieldIn,knon_mpi_para,knon_mpi_begin_para(:)-1,MPI_INTEGER,   &
                        FieldOut,knon_mpi,MPI_INTEGER,                                &
                        mpi_rank_root,COMM_LMDZ_PHY,ierr)
#endif
    ELSE
      FieldOut(:)=FieldIn(:)
    ENDIF
  
  END SUBROUTINE scatter_surf_mpi_i



  SUBROUTINE scatter_surf_r(FieldIn, FieldOut)
  USE mod_phys_lmdz_para
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
    REAL :: FieldTmp(knon_mpi)
    
    IF (is_omp_root) CALL scatter_surf_mpi_r(FieldIn,FieldTmp)
    CALL scatter_surf_omp_r(FieldTmp,FieldOut)
    
  END SUBROUTINE scatter_surf_r


  SUBROUTINE scatter_surf_omp_r(FieldIn,FieldOut)
  USE mod_phys_lmdz_para
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
  
    INTEGER,SAVE,ALLOCATABLE :: Field_tmp(:)
    
    IF (is_omp_root) ALLOCATE(Field_tmp(knon_mpi))
    IF (is_omp_root) Field_tmp(:)=FieldIn(:)
!$OMP BARRIER        
    FieldOut(:)=Field_tmp(knon_omp_begin:knon_omp_end)
!$OMP BARRIER
    IF (is_omp_root) DEALLOCATE(Field_tmp)
    
  END SUBROUTINE  scatter_surf_omp_r
  
     
  SUBROUTINE scatter_surf_mpi_r(FieldIn,FieldOut)
  USE mod_phys_lmdz_para, mpi_rank_root => mpi_root
#ifdef CPP_MPI
  INCLUDE 'mpif.h'
#endif
    REAL :: FieldIn(:)
    REAL :: FieldOut(:)
    INTEGER :: ierr
    
    IF (is_using_mpi) THEN
#ifdef CPP_MPI
      CALL MPI_Scatterv(FieldIn,knon_mpi_para,knon_mpi_begin_para(:)-1,MPI_INTEGER,   &
                        FieldOut,knon_mpi,MPI_INTEGER,                                &
                        mpi_rank_root,COMM_LMDZ_PHY,ierr)
#endif
    ELSE
      FieldOut(:)=FieldIn(:)
    ENDIF
  
  END SUBROUTINE scatter_surf_mpi_r

END MODULE mod_surf_para
