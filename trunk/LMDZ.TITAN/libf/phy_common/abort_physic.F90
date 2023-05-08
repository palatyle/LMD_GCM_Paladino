! $Id: $
      SUBROUTINE abort_physic(modname, message, ierr)
     
#ifdef CPP_IOIPSL
      USE IOIPSL
#else
! if not using IOIPSL, we still need to use (a local version of) getin_dump
      USE ioipsl_getincom
#endif
      USE mod_phys_lmdz_para
      USE print_control_mod, ONLY: lunout
      IMPLICIT NONE
!
! Stops the simulation cleanly, closing files and printing various
! comments
!
!  Input: modname = name of calling program
!         message = stuff to print
!         ierr    = severity of situation ( = 0 normal )

      character(len=*), intent(in):: modname
      integer ierr, ierror_mpi
      character(len=*), intent(in):: message

      write(lunout,*) 'in abort_physic'
#ifdef CPP_IOIPSL
!$OMP MASTER
      call histclo
      call restclo
      if (mpi_rank .eq. 0) then
         call getin_dump
      endif
!$OMP END MASTER
#endif

      write(lunout,*) 'Stopping in ', modname
      write(lunout,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(lunout,*) 'Everything is cool'
      else
        write(lunout,*) 'Houston, we have a problem, ierr = ', ierr
#ifdef CPP_MPI
!$OMP CRITICAL (MPI_ABORT_PHYSIC)
        call MPI_ABORT(COMM_LMDZ_PHY, 1, ierror_mpi)
!$OMP END CRITICAL (MPI_ABORT_PHYSIC)
#else
        stop 1
#endif          
      endif
      END
