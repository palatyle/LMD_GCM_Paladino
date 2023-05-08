










!
!$Id: mod_phys_lmdz_omp_data.F90 2429 2016-01-27 12:43:09Z fairhead $
!
MODULE mod_phys_lmdz_omp_data

  INTEGER,SAVE :: omp_size
  INTEGER,SAVE :: omp_rank
  LOGICAL,SAVE :: is_omp_root
  LOGICAL,SAVE :: is_omp_master  ! alias of is_omp_root
  LOGICAL,SAVE :: is_using_omp
  LOGICAL,SAVE :: is_north_pole_phy, is_south_pole_phy
  
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_nb
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_begin
  INTEGER,SAVE,DIMENSION(:),ALLOCATABLE :: klon_omp_para_end    
  
  INTEGER,SAVE :: klon_omp
  INTEGER,SAVE :: klon_omp_begin
  INTEGER,SAVE :: klon_omp_end
!$OMP  THREADPRIVATE(omp_rank,klon_omp,is_omp_root,is_omp_master,klon_omp_begin,klon_omp_end)
!$OMP  THREADPRIVATE(is_north_pole_phy, is_south_pole_phy)

CONTAINS
  
  SUBROUTINE Init_phys_lmdz_omp_data(klon_mpi)
    USE dimphy 
    USE mod_phys_lmdz_mpi_data, ONLY : is_north_pole_dyn, is_south_pole_dyn
    IMPLICIT NONE
    INTEGER, INTENT(in) :: klon_mpi

    INTEGER :: i

    CHARACTER (LEN=20) :: modname='Init_phys_lmdz_omp_data'
    CHARACTER (LEN=80) :: abort_message



    is_using_omp=.FALSE.
    omp_size=1
    omp_rank=0

   is_omp_root=.FALSE.
!$OMP MASTER
   IF (omp_rank==0) THEN
     is_omp_root=.TRUE.
   ELSE
     abort_message = 'ANORMAL : OMP_MASTER /= 0'
     CALL abort_physic (modname,abort_message,1)
   ENDIF
!$OMP END MASTER
   is_omp_master=is_omp_root

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

   if ((is_north_pole_dyn) .AND. (omp_rank == 0 )) then
      is_north_pole_phy = .TRUE.
    else
      is_north_pole_phy = .FALSE.
    endif
    if ((is_south_pole_dyn) .AND. (omp_rank == omp_size-1)) then
      is_south_pole_phy = .TRUE.
    else
      is_south_pole_phy = .FALSE.
    endif
   
    klon_omp=klon_omp_para_nb(omp_rank)
    klon_omp_begin=klon_omp_para_begin(omp_rank)
    klon_omp_end=klon_omp_para_end(omp_rank)
    
    CALL Print_module_data
    
  END SUBROUTINE Init_phys_lmdz_omp_data

  SUBROUTINE Print_module_data
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE

!$OMP CRITICAL  
  WRITE(lunout,*)'--------> TASK ',omp_rank
  WRITE(lunout,*)'omp_size =',omp_size
  WRITE(lunout,*)'omp_rank =',omp_rank
  WRITE(lunout,*)'is_omp_root =',is_omp_root
  WRITE(lunout,*)'klon_omp_para_nb =',klon_omp_para_nb
  WRITE(lunout,*)'klon_omp_para_begin =',klon_omp_para_begin
  WRITE(lunout,*)'klon_omp_para_end =',klon_omp_para_end    
  WRITE(lunout,*)'klon_omp =',klon_omp
  WRITE(lunout,*)'klon_omp_begin =',klon_omp_begin
  WRITE(lunout,*)'klon_omp_end =',klon_omp_end    
!$OMP END CRITICAL

  END SUBROUTINE Print_module_data
END MODULE mod_phys_lmdz_omp_data
