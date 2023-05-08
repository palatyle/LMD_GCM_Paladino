










! $Id: $
      SUBROUTINE abort_physic(modname, message, ierr)
     
! if not using IOIPSL, we still need to use (a local version of) getin_dump
      USE ioipsl_getincom
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

      write(lunout,*) 'Stopping in ', modname
      write(lunout,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(lunout,*) 'Everything is cool'
      else
        write(lunout,*) 'Houston, we have a problem, ierr = ', ierr
        stop 1
      endif
      END
