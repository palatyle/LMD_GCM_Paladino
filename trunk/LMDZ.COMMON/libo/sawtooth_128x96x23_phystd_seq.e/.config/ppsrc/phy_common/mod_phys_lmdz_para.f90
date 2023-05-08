










!
! $Id: mod_phys_lmdz_para.F90 2429 2016-01-27 12:43:09Z fairhead $
!
MODULE mod_phys_lmdz_para
  USE mod_phys_lmdz_transfert_para
  USE mod_phys_lmdz_mpi_data
  USE mod_phys_lmdz_omp_data
    
  INTEGER,SAVE :: klon_loc
  LOGICAL,SAVE :: is_sequential
  LOGICAL,SAVE :: is_parallel
  LOGICAL,SAVE :: is_master

  
!$OMP THREADPRIVATE(klon_loc,is_master)
  
CONTAINS

  SUBROUTINE Init_phys_lmdz_para(nbp,nbp_lon,nbp_lat,communicator)
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbp
    INTEGER,INTENT(IN) :: nbp_lon
    INTEGER,INTENT(IN) :: nbp_lat
    INTEGER,INTENT(IN) :: communicator

    CALL Init_phys_lmdz_mpi_data(nbp,nbp_lon,nbp_lat,communicator)
!$OMP PARALLEL
    CALL Init_phys_lmdz_omp_data(klon_mpi)
    klon_loc=klon_omp
    IF (is_mpi_root .AND. is_omp_root) THEN 
       is_master=.TRUE.
     ELSE
       is_master=.FALSE.
     ENDIF
     CALL Test_transfert
!$OMP END PARALLEL    
     IF (is_using_mpi .OR. is_using_omp) THEN
       is_sequential=.FALSE.
       is_parallel=.TRUE.
     ELSE
       is_sequential=.TRUE.
       is_parallel=.FALSE.
     ENDIF


      
  END SUBROUTINE Init_phys_lmdz_para

  SUBROUTINE Test_transfert
  USE mod_grid_phy_lmdz
  USE print_control_mod, ONLY: lunout
  IMPLICIT NONE
 
    REAL :: Test_Field1d_glo(klon_glo,nbp_lev)
    REAL :: tmp1d_glo(klon_glo,nbp_lev)
    REAL :: Test_Field2d_glo(nbp_lon,nbp_lat,nbp_lev)
    REAL :: tmp2d_glo(nbp_lon,nbp_lat,nbp_lev)
    REAL :: Test_Field1d_loc(klon_loc,nbp_lev)
    REAL :: Test_Field2d_loc(nbp_lon,jj_nb,nbp_lev)
    REAL :: CheckSum
    
    INTEGER :: i,l
  
    Test_Field1d_glo = 0.
    Test_Field2d_glo = 0.
    Test_Field1d_loc = 0.
    Test_Field2d_loc = 0.
  
    IF (is_mpi_root) THEN
!$OMP MASTER
      DO l=1,nbp_lev
        DO i=1,klon_glo
!          Test_Field1d_glo(i,l)=MOD(i,10)+10*(l-1)
           Test_Field1d_glo(i,l)=1
        ENDDO
      ENDDO
!$OMP END MASTER  
    ENDIF
  
    CALL Scatter(Test_Field1d_glo,Test_Field1d_loc)
    CALL Gather(Test_Field1d_loc,tmp1d_glo)
  
    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo-tmp1d_glo)
      WRITE(lunout,*) "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF
    
    CALL grid1dTo2d_glo(Test_Field1d_glo,Test_Field2d_glo)
    CALL scatter2D(Test_Field2d_glo,Test_Field1d_loc)
    CALL gather2d(Test_Field1d_loc,Test_Field2d_glo)
    CALL grid2dTo1d_glo(Test_Field2d_glo,tmp1d_glo)

    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo-tmp1d_glo)
      WRITE(lunout,*) "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF

    CALL bcast(Test_Field1d_glo)
    CALL reduce_sum(Test_Field1d_glo,tmp1d_glo)

    IF (is_mpi_root) THEN
!$OMP MASTER  
      Checksum=sum(Test_Field1d_glo*omp_size*mpi_size-tmp1d_glo)
      WRITE(lunout,*) "------> Checksum =",Checksum," MUST BE 0"
!$OMP END MASTER
    ENDIF
    
     
   END SUBROUTINE Test_transfert
  
END MODULE mod_phys_lmdz_para
    
