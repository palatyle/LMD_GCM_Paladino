!
!$Header$
!
MODULE mod_phys_lmdz_mpi_data
  USE mod_const_mpi
  
  INTEGER,SAVE :: ii_begin
  INTEGER,SAVE :: ii_end
  INTEGER,SAVE :: jj_begin
  INTEGER,SAVE :: jj_end
  INTEGER,SAVE :: jj_nb
  INTEGER,SAVE :: ij_begin
  INTEGER,SAVE :: ij_end
  INTEGER,SAVE :: ij_nb
  INTEGER,SAVE :: klon_mpi_begin
  INTEGER,SAVE :: klon_mpi_end
  INTEGER,SAVE :: klon_mpi
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: jj_para_nb
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: jj_para_begin
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: jj_para_end

  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ii_para_begin
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ii_para_end

  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ij_para_nb
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ij_para_begin
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ij_para_end

  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: klon_mpi_para_nb
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: klon_mpi_para_begin
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: klon_mpi_para_end 

  
  INTEGER,SAVE :: mpi_rank
  INTEGER,SAVE :: mpi_size
  INTEGER,SAVE :: mpi_root
  LOGICAL,SAVE :: is_mpi_root
  LOGICAL,SAVE :: is_using_mpi
  
  
  LOGICAL,SAVE :: is_north_pole
  LOGICAL,SAVE :: is_south_pole
  INTEGER,SAVE :: COMM_LMDZ_PHY

CONTAINS
  
  SUBROUTINE Init_phys_lmdz_mpi_data(iim,jjp1,nb_proc,distrib)
  USE mod_const_mpi, ONLY : COMM_LMDZ
  IMPLICIT NONE
    INTEGER,INTENT(in) :: iim
    INTEGER,INTENT(in) :: jjp1
    INTEGER,INTENT(in) :: nb_proc
    INTEGER,INTENT(in) :: distrib(0:nb_proc-1)
    
    INTEGER :: ierr
    INTEGER :: klon_glo
    INTEGER :: i
    
#ifdef CPP_MPI
    is_using_mpi=.TRUE.
#else
    is_using_mpi=.FALSE.
#endif
    
    if (iim.eq.1) then
       klon_glo=1
    else
       klon_glo=iim*(jjp1-2)+2
    endif
    
    COMM_LMDZ_PHY=COMM_LMDZ

    IF (is_using_mpi) THEN    
#ifdef CPP_MPI
      CALL MPI_COMM_SIZE(COMM_LMDZ_PHY,mpi_size,ierr)    
      CALL MPI_COMM_RANK(COMM_LMDZ_PHY,mpi_rank,ierr)
#endif
    ELSE
      mpi_size=1
      mpi_rank=0
    ENDIF
    
    IF (mpi_rank == 0) THEN
      mpi_root = 0
      is_mpi_root = .true.
    ENDIF
    
    IF (mpi_rank == 0) THEN 
      is_north_pole = .TRUE.
    ELSE
      is_north_pole = .FALSE.
    ENDIF
    
    IF (mpi_rank == mpi_size-1) THEN
      is_south_pole = .TRUE.
    ELSE
      is_south_pole = .FALSE.
    ENDIF
    
    ALLOCATE(jj_para_nb(0:mpi_size-1))
    ALLOCATE(jj_para_begin(0:mpi_size-1))
    ALLOCATE(jj_para_end(0:mpi_size-1))
    
    ALLOCATE(ij_para_nb(0:mpi_size-1))
    ALLOCATE(ij_para_begin(0:mpi_size-1))
    ALLOCATE(ij_para_end(0:mpi_size-1))
    
    ALLOCATE(ii_para_begin(0:mpi_size-1))
    ALLOCATE(ii_para_end(0:mpi_size-1))

    ALLOCATE(klon_mpi_para_nb(0:mpi_size-1))
    ALLOCATE(klon_mpi_para_begin(0:mpi_size-1))
    ALLOCATE(klon_mpi_para_end(0:mpi_size-1))
  
      
    klon_mpi_para_nb(0:mpi_size-1)=distrib(0:nb_proc-1)

    DO i=0,mpi_size-1
      IF (i==0) THEN 
        klon_mpi_para_begin(i)=1
      ELSE 
        klon_mpi_para_begin(i)=klon_mpi_para_end(i-1)+1
      ENDIF
        klon_mpi_para_end(i)=klon_mpi_para_begin(i)+klon_mpi_para_nb(i)-1
    ENDDO


    DO i=0,mpi_size-1
      
      IF (i==0) THEN
        ij_para_begin(i) = 1
      ELSE
        ij_para_begin(i) = klon_mpi_para_begin(i)+iim-1
      ENDIF

      jj_para_begin(i) = (ij_para_begin(i)-1)/iim + 1
      ii_para_begin(i) = MOD(ij_para_begin(i)-1,iim) + 1

      
      ij_para_end(i) = klon_mpi_para_end(i)+iim-1
      jj_para_end(i) = (ij_para_end(i)-1)/iim + 1
      ii_para_end(i) = MOD(ij_para_end(i)-1,iim) + 1


      ij_para_nb(i) = ij_para_end(i)-ij_para_begin(i)+1
      jj_para_nb(i) = jj_para_end(i)-jj_para_begin(i)+1
         
    ENDDO
  
    ii_begin = ii_para_begin(mpi_rank)
    ii_end   = ii_para_end(mpi_rank)
    jj_begin = jj_para_begin(mpi_rank)
    jj_end   = jj_para_end(mpi_rank)
    jj_nb    = jj_para_nb(mpi_rank)
    ij_begin = ij_para_begin(mpi_rank)
    ij_end   = ij_para_end(mpi_rank)
    ij_nb    = ij_para_nb(mpi_rank)
    klon_mpi_begin = klon_mpi_para_begin(mpi_rank)
    klon_mpi_end   = klon_mpi_para_end(mpi_rank)
    klon_mpi       = klon_mpi_para_nb(mpi_rank)
   
    CALL Print_module_data
    
  END SUBROUTINE Init_phys_lmdz_mpi_data

  SUBROUTINE print_module_data
  IMPLICIT NONE
  
  
    PRINT *, 'ii_begin =', ii_begin
    PRINT *, 'ii_end =', ii_end
    PRINT *, 'jj_begin =',jj_begin
    PRINT *, 'jj_end =', jj_end
    PRINT *, 'jj_nb =', jj_nb
    PRINT *, 'ij_begin =', ij_begin
    PRINT *, 'ij_end =', ij_end
    PRINT *, 'ij_nb =', ij_nb
    PRINT *, 'klon_mpi_begin =', klon_mpi_begin
    PRINT *, 'klon_mpi_end =', klon_mpi_end
    PRINT *, 'klon_mpi =', klon_mpi
    PRINT *, 'jj_para_nb =', jj_para_nb
    PRINT *, 'jj_para_begin =', jj_para_begin
    PRINT *, 'jj_para_end =', jj_para_end
    PRINT *, 'ii_para_begin =', ii_para_begin
    PRINT *, 'ii_para_end =', ii_para_end
    PRINT *, 'ij_para_nb =', ij_para_nb
    PRINT *, 'ij_para_begin =', ij_para_begin
    PRINT *, 'ij_para_end =', ij_para_end
    PRINT *, 'klon_mpi_para_nb =', klon_mpi_para_nb
    PRINT *, 'klon_mpi_para_begin =', klon_mpi_para_begin
    PRINT *, 'klon_mpi_para_end  =', klon_mpi_para_end 
    PRINT *, 'mpi_rank =', mpi_rank
    PRINT *, 'mpi_size =', mpi_size
    PRINT *, 'mpi_root =', mpi_root
    PRINT *, 'is_mpi_root =', is_mpi_root
    PRINT *, 'is_north_pole =', is_north_pole
    PRINT *, 'is_south_pole =', is_south_pole
    PRINT *, 'COMM_LMDZ_PHY =', COMM_LMDZ_PHY
  
  END SUBROUTINE print_module_data
  
END MODULE mod_phys_lmdz_mpi_data
