!
! $Id: getparam.F90 2094 2014-07-16 16:55:47Z lguez $
!
MODULE getparam
#ifdef CPP_IOIPSL
   USE IOIPSL
#else
! if not using IOIPSL, we still need to use (a local version of) getin
   USE ioipsl_getincom
#endif

   INTERFACE getpar
     MODULE PROCEDURE getparamr,getparami,getparaml
   END INTERFACE
   private getparamr,getparami,getparaml

   INTEGER, PARAMETER :: out_eff=99

CONTAINS
  SUBROUTINE ini_getparam(fichier)
  USE parallel_lmdz
    !
    IMPLICIT NONE
    !
    CHARACTER*(*) :: fichier
    IF (mpi_rank==0) OPEN(out_eff,file=fichier,status='unknown',form='formatted')
    
  END SUBROUTINE ini_getparam

  SUBROUTINE fin_getparam
  USE parallel_lmdz
    !
    IMPLICIT NONE
    !
      IF (mpi_rank==0) CLOSE(out_eff)

  END SUBROUTINE fin_getparam

  SUBROUTINE getparamr(TARGET,def_val,ret_val,comment)
  USE parallel_lmdz
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    REAL :: def_val
    REAL :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET,ret_val)

    IF (mpi_rank==0) THEN
      write(out_eff,*) '######################################'
      write(out_eff,*) '#### ',comment,' #####'
      write(out_eff,*) TARGET,'=',ret_val
    ENDIF
    
  END SUBROUTINE getparamr

  SUBROUTINE getparami(TARGET,def_val,ret_val,comment)
  USE parallel_lmdz
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    INTEGER :: def_val
    INTEGER :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET,ret_val)

    IF (mpi_rank==0) THEN
      write(out_eff,*) '######################################'
      write(out_eff,*) '#### ',comment,' #####'
      write(out_eff,*) comment
      write(out_eff,*) TARGET,'=',ret_val
    ENDIF
    
  END SUBROUTINE getparami

  SUBROUTINE getparaml(TARGET,def_val,ret_val,comment)
  USE parallel_lmdz
    !
    IMPLICIT NONE
    !
    !   Get a real scalar. We first check if we find it
    !   in the database and if not we get it from the run.def
    !
    !   getinr1d and getinr2d are written on the same pattern
    !
    CHARACTER*(*) :: TARGET
    LOGICAL :: def_val
    LOGICAL :: ret_val
    CHARACTER*(*) :: comment

    ret_val=def_val
    call getin(TARGET,ret_val)

    IF (mpi_rank==0) THEN
      write(out_eff,*) '######################################'
      write(out_eff,*) '#### ',comment,' #####'
      write(out_eff,*) TARGET,'=',ret_val
    ENDIF
       
  END SUBROUTINE getparaml


END MODULE getparam
