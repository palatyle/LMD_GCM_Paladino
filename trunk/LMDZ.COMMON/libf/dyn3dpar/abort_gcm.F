!
! $Id: abort_gcm.F 1907 2013-11-26 13:10:46Z lguez $
!
c
c
      SUBROUTINE abort_gcm(modname, message, ierr)
     
#ifdef CPP_IOIPSL
      USE IOIPSL
#else
! if not using IOIPSL, we still need to use (a local version of) getin_dump
      USE ioipsl_getincom
#endif
      USE parallel_lmdz




#include "iniprint.h"
 
C
C Stops the simulation cleanly, closing files and printing various
C comments
C
C  Input: modname = name of calling program
C         message = stuff to print
C         ierr    = severity of situation ( = 0 normal )

      character(len=*), intent(in):: modname
      integer ierr, ierror_mpi
      character(len=*), intent(in):: message

      write(lunout,*) 'in abort_gcm'
#ifdef CPP_IOIPSL
c$OMP MASTER
      call histclo
      call restclo
      if (MPI_rank .eq. 0) then
         call getin_dump
      endif
c$OMP END MASTER
#endif



c     call histclo(2)
c     call histclo(3)
c     call histclo(4)
c     call histclo(5)
      write(lunout,*) 'Stopping in ', modname
      write(lunout,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(lunout,*) 'Everything is cool'
      else
        write(lunout,*) 'Houston, we have a problem, ierr = ', ierr
#ifdef CPP_MPI
C$OMP CRITICAL (MPI_ABORT_GCM)
        call MPI_ABORT(COMM_LMDZ, 1, ierror_mpi)
C$OMP END CRITICAL (MPI_ABORT_GCM)
#else
        stop 1
#endif          
      endif
      END

