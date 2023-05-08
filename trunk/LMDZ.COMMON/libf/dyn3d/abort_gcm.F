!
! $Id: abort_gcm.F 1425 2010-09-02 13:45:23Z lguez $
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

#ifdef CPP_XIOS
    ! ug Pour les sorties XIOS
      USE wxios
#endif

#include "iniprint.h"
 
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

#ifdef CPP_XIOS
    !Fermeture propre de XIOS
      CALL wxios_close()
#endif

#ifdef CPP_IOIPSL
      call histclo
      call restclo
#endif
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
