










!
! $Id: abort_gcm.F 1907 2013-11-26 13:10:46Z lguez $
!
c
c
      SUBROUTINE abort_gcm(modname, message, ierr)
     
! if not using IOIPSL, we still need to use (a local version of) getin_dump
      USE ioipsl_getincom
      USE parallel_lmdz




!
! $Header$
!
!
! gestion des impressions de sorties et de d�bogage
! lunout:    unit� du fichier dans lequel se font les sorties 
!                           (par defaut 6, la sortie standard)
! prt_level: niveau d'impression souhait� (0 = minimum)
!
      INTEGER lunout, prt_level
      COMMON /comprint/ lunout, prt_level
 
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
C$OMP CRITICAL (MPI_ABORT_GCM)
        call MPI_ABORT(COMM_LMDZ, 1, ierror_mpi)
C$OMP END CRITICAL (MPI_ABORT_GCM)
      endif
      END

