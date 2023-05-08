      SUBROUTINE abort_gcm(modname, message, ierr)
     
C
C Stops the simulation cleanly, closing files and printing various
C comments
C
C  Input: modname = name of calling program
C         message = stuff to print
C         ierr    = severity of situation ( = 0 normal )

      character*20 modname
      integer ierr
      character*80 message

      write(*,*) 'in abort_gcm'
      write(*,*) 'Stopping in ', modname
      write(*,*) 'Reason = ',message
      if (ierr .eq. 0) then
        write(*,*) 'Everything is cool'
      else
        write(*,*) 'Houston, we have a problem ', ierr
      endif
      STOP
      END
