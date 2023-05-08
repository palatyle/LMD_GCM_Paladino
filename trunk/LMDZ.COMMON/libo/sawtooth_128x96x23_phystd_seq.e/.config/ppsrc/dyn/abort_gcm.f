










!
! $Id: abort_gcm.F 1425 2010-09-02 13:45:23Z lguez $
!
c
c
      SUBROUTINE abort_gcm(modname, message, ierr)
     
! if not using IOIPSL, we still need to use (a local version of) getin_dump
      USE ioipsl_getincom


!
! $Header$
!
!
! gestion des impressions de sorties et de débogage
! lunout:    unité du fichier dans lequel se font les sorties 
!                           (par defaut 6, la sortie standard)
! prt_level: niveau d'impression souhaité (0 = minimum)
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
      integer, intent(in):: ierr
      character(len=*), intent(in):: message

      write(lunout,*) 'in abort_gcm'


      call getin_dump
c     call histclo(2)
c     call histclo(3)
c     call histclo(4)
c     call histclo(5)
      write(lunout,*) 'Stopping in ', modname
      write(lunout,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(lunout,*) 'Everything is cool'
        stop
      else
        write(lunout,*) 'Houston, we have a problem ', ierr
        stop 1
      endif
      END
