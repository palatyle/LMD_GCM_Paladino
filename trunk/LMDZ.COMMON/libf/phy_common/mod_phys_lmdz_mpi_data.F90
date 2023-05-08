!
!$Id$
!
MODULE mod_phys_lmdz_mpi_data
  
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
  INTEGER,SAVE :: mpi_master
  LOGICAL,SAVE :: is_mpi_root
  LOGICAL,SAVE :: is_using_mpi
  
  
  LOGICAL,SAVE :: is_north_pole_dyn
  LOGICAL,SAVE :: is_south_pole_dyn
  INTEGER,SAVE :: COMM_LMDZ_PHY
  INTEGER,SAVE :: MPI_REAL_LMDZ   ! MPI_REAL8

CONTAINS
  
  SUBROUTINE init_phys_lmdz_mpi_data(nbp, nbp_lon, nbp_lat, communicator)
  IMPLICIT NONE
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    INTEGER,INTENT(IN) :: nbp
    INTEGER,INTENT(IN) :: nbp_lon
    INTEGER,INTENT(IN) :: nbp_lat
    INTEGER,INTENT(IN) :: communicator
    
    INTEGER,ALLOCATABLE :: distrib(:)
    INTEGER :: ierr
    INTEGER :: klon_glo
    INTEGER :: i
    
#ifdef CPP_MPI
    is_using_mpi=.TRUE.
#else
    is_using_mpi=.FALSE.
#endif
    
    if ((nbp_lon.eq.1).and.(nbp_lat.eq.1)) then ! running 1D column model
       klon_glo=1
    else
    ! The usual global physics grid: 1 point for each pole and nbp_lon points
    ! for all other latitudes
       klon_glo=nbp_lon*(nbp_lat-2)+2
    endif
    
    COMM_LMDZ_PHY=communicator

    IF (is_using_mpi) THEN    
#ifdef CPP_MPI
      MPI_REAL_LMDZ=MPI_REAL8
      CALL MPI_COMM_SIZE(COMM_LMDZ_PHY,mpi_size,ierr)    
      CALL MPI_COMM_RANK(COMM_LMDZ_PHY,mpi_rank,ierr)
#endif
    ELSE
      mpi_size=1
      mpi_rank=0
    ENDIF
    
    ALLOCATE(distrib(0:mpi_size-1))

    IF (is_using_mpi) THEN    
#ifdef CPP_MPI
    CALL MPI_ALLGATHER(nbp,1,MPI_INTEGER,distrib,1,MPI_INTEGER,COMM_LMDZ_PHY,ierr)
#endif
    ELSE
     distrib(:)=nbp
    ENDIF


    IF (mpi_rank == 0) THEN
      mpi_master = 0
      is_mpi_root = .true.
    ENDIF
    
    IF (mpi_rank == 0) THEN 
      is_north_pole_dyn = .TRUE.
    ELSE
      is_north_pole_dyn = .FALSE.
    ENDIF
    
    IF (mpi_rank == mpi_size-1) THEN
      is_south_pole_dyn = .TRUE.
    ELSE
      is_south_pole_dyn = .FALSE.
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
  
      
    klon_mpi_para_nb(0:mpi_size-1)=distrib(0:mpi_size-1)

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
        ij_para_begin(i) = klon_mpi_para_begin(i)+nbp_lon-1
      ENDIF

      jj_para_begin(i) = (ij_para_begin(i)-1)/nbp_lon + 1
      ii_para_begin(i) = MOD(ij_para_begin(i)-1,nbp_lon) + 1

      
      ij_para_end(i) = klon_mpi_para_end(i)+nbp_lon-1
      jj_para_end(i) = (ij_para_end(i)-1)/nbp_lon + 1
      ii_para_end(i) = MOD(ij_para_end(i)-1,nbp_lon) + 1


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
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
  
    WRITE(lunout,*) 'ii_begin =', ii_begin
    WRITE(lunout,*) 'ii_end =', ii_end
    WRITE(lunout,*) 'jj_begin =',jj_begin
    WRITE(lunout,*) 'jj_end =', jj_end
    WRITE(lunout,*) 'jj_nb =', jj_nb
    WRITE(lunout,*) 'ij_begin =', ij_begin
    WRITE(lunout,*) 'ij_end =', ij_end
    WRITE(lunout,*) 'ij_nb =', ij_nb
    WRITE(lunout,*) 'klon_mpi_begin =', klon_mpi_begin
    WRITE(lunout,*) 'klon_mpi_end =', klon_mpi_end
    WRITE(lunout,*) 'klon_mpi =', klon_mpi
    WRITE(lunout,*) 'jj_para_nb =', jj_para_nb
    WRITE(lunout,*) 'jj_para_begin =', jj_para_begin
    WRITE(lunout,*) 'jj_para_end =', jj_para_end
    WRITE(lunout,*) 'ii_para_begin =', ii_para_begin
    WRITE(lunout,*) 'ii_para_end =', ii_para_end
    WRITE(lunout,*) 'ij_para_nb =', ij_para_nb
    WRITE(lunout,*) 'ij_para_begin =', ij_para_begin
    WRITE(lunout,*) 'ij_para_end =', ij_para_end
    WRITE(lunout,*) 'klon_mpi_para_nb =', klon_mpi_para_nb
    WRITE(lunout,*) 'klon_mpi_para_begin =', klon_mpi_para_begin
    WRITE(lunout,*) 'klon_mpi_para_end  =', klon_mpi_para_end 
    WRITE(lunout,*) 'mpi_rank =', mpi_rank
    WRITE(lunout,*) 'mpi_size =', mpi_size
    WRITE(lunout,*) 'mpi_master =', mpi_master
    WRITE(lunout,*) 'is_mpi_root =', is_mpi_root
    WRITE(lunout,*) 'is_north_pole_dyn =', is_north_pole_dyn
    WRITE(lunout,*) 'is_south_pole_dyn =', is_south_pole_dyn
    WRITE(lunout,*) 'COMM_LMDZ_PHY =', COMM_LMDZ_PHY
  
  END SUBROUTINE print_module_data
  
END MODULE mod_phys_lmdz_mpi_data
