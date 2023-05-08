! 
! $Id: mod_interface_dyn_phys.F90 2351 2015-08-25 15:14:59Z emillour $
!
MODULE mod_interface_dyn_phys
  INTEGER,SAVE,dimension(:),allocatable :: index_i
  INTEGER,SAVE,dimension(:),allocatable :: index_j
  
  
CONTAINS
  
#ifdef CPP_PARA
! Interface with parallel physics,
  SUBROUTINE Init_interface_dyn_phys
    USE mod_phys_lmdz_mpi_data
    IMPLICIT NONE
    include 'dimensions.h'    
    
    INTEGER :: i,j,k
    
    ALLOCATE(index_i(klon_mpi))
    ALLOCATE(index_j(klon_mpi))
    
    k=1
    IF (is_north_pole_dyn) THEN
      index_i(k)=1
      index_j(k)=1
      k=2
    ELSE
      DO i=ii_begin,iim
	index_i(k)=i
	index_j(k)=jj_begin
	k=k+1
       ENDDO
    ENDIF
    
    DO j=jj_begin+1,jj_end-1
      DO i=1,iim
	index_i(k)=i
	index_j(k)=j
	k=k+1
      ENDDO
    ENDDO
    
    IF (is_south_pole_dyn) THEN
      index_i(k)=1
      index_j(k)=jj_end
    ELSE
      DO i=1,ii_end
	index_i(k)=i
	index_j(k)=jj_end
	k=k+1
       ENDDO
    ENDIF
  
  END SUBROUTINE Init_interface_dyn_phys 
#else
  SUBROUTINE Init_interface_dyn_phys
  ! dummy routine for seq case
  END SUBROUTINE Init_interface_dyn_phys 
#endif
! of #ifdef CPP_PARA
END MODULE mod_interface_dyn_phys
