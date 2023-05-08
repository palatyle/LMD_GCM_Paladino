! 
! $Id: mod_const_mpi.F90 2055 2014-06-04 12:33:27Z acaubel $
!
MODULE mod_const_mpi
  IMPLICIT NONE
  INTEGER,SAVE :: COMM_LMDZ
  INTEGER,SAVE :: MPI_REAL_LMDZ
 

CONTAINS 

  SUBROUTINE Init_const_mpi
#ifdef CPP_IOIPSL
    USE IOIPSL, ONLY: getin
#else
! if not using IOIPSL, we still need to use (a local version of) getin
    USE ioipsl_getincom, only: getin
#endif
! Use of Oasis-MCT coupler 
#ifdef CPP_OMCT
    USE mod_prism
#endif
#ifdef CPP_XIOS
    USE wxios, only: wxios_init
#endif
    IMPLICIT NONE
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif

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
#ifdef CPP_COUPLE
!$OMP MASTER
#ifdef CPP_XIOS
      CALL wxios_init("LMDZ", outcom=COMM_LMDZ, type_ocean=type_ocean)
#else
       CALL prism_init_comp_proto (comp_id, 'LMDZ', ierr)
       CALL prism_get_localcomm_proto(COMM_LMDZ,ierr)
#endif
!$OMP END MASTER
#endif
#ifdef CPP_MPI
      MPI_REAL_LMDZ=MPI_REAL8
#endif
    ELSE
      CALL init_mpi
    ENDIF

  END SUBROUTINE Init_const_mpi
  
  SUBROUTINE Init_mpi
#ifdef CPP_XIOS
    USE wxios, only: wxios_init
#endif
  IMPLICIT NONE
#ifdef CPP_MPI
     INCLUDE 'mpif.h'
#endif
    INTEGER             :: ierr
    INTEGER             :: thread_required
    INTEGER             :: thread_provided

#ifdef CPP_MPI
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
#ifdef CPP_XIOS
      CALL wxios_init("LMDZ", outcom=COMM_LMDZ)
#endif
!$OMP END MASTER
#else
#ifdef CPP_XIOS
!$OMP MASTER
      CALL wxios_init("LMDZ")
!$OMP END MASTER
#endif
#endif

   END SUBROUTINE Init_mpi
    
END MODULE mod_const_mpi

