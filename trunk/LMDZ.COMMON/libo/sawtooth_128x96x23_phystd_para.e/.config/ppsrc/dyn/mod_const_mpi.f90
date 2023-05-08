










! 
! $Id: mod_const_mpi.F90 2055 2014-06-04 12:33:27Z acaubel $
!
MODULE mod_const_mpi
  IMPLICIT NONE
  INTEGER,SAVE :: COMM_LMDZ
  INTEGER,SAVE :: MPI_REAL_LMDZ
 

CONTAINS 

  SUBROUTINE Init_const_mpi
! if not using IOIPSL, we still need to use (a local version of) getin
    USE ioipsl_getincom, only: getin
! Use of Oasis-MCT coupler 
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER             :: ierr
    INTEGER             :: comp_id
    INTEGER             :: thread_required
    INTEGER             :: thread_provided
    CHARACTER(len = 6)  :: type_ocean

!$OMP MASTER
    type_ocean = 'force '
    CALL getin('type_ocean', type_ocean)
!$OMP END MASTER
!$OMP BARRIER

    IF (type_ocean=='couple') THEN
      MPI_REAL_LMDZ=MPI_REAL8
    ELSE
      CALL init_mpi
    ENDIF

  END SUBROUTINE Init_const_mpi
  
  SUBROUTINE Init_mpi
  IMPLICIT NONE
     INCLUDE 'mpif.h'
    INTEGER             :: ierr
    INTEGER             :: thread_required
    INTEGER             :: thread_provided

!$OMP MASTER
      thread_required=MPI_THREAD_SERIALIZED

      CALL MPI_INIT_THREAD(thread_required,thread_provided,ierr)
      IF (thread_provided < thread_required) THEN
        PRINT *,'Warning : The multithreaded level of MPI librairy do not provide the requiered level',  &
                ' in mod_const_mpi::Init_const_mpi'
      ENDIF
      COMM_LMDZ=MPI_COMM_WORLD
      MPI_REAL_LMDZ=MPI_REAL8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation de XIOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP END MASTER

   END SUBROUTINE Init_mpi
    
END MODULE mod_const_mpi

