!***********************************************************************

      subroutine read_phototable

!***********************************************************************
!
!   subject:
!   --------
!
!   read photolysis lookup table
!
!   VERSION: 8/10/2014
!
!   Author:   Franck Lefevre
!
!   Arguments:
!   ----------
!
!   The output variable is jphot and is put in common chimiedata.
!
!***********************************************************************

      use ioipsl_getincom 
      use ioipsl_getin_p_mod, only: getin_p
      use datafile_mod
      implicit none

#include "chimiedata.h"

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     local:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer :: fic, ij, iozo, isza, itemp, iz, itau, ierr
      real    :: xsza

      character(len = 128) :: phototable ! photolysis table file name

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! set photolysis table input file name
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



! look for a " phototable= ..." option in def files
     write(*,*) "Directory where external input files are:"
      phototable = "jmars.20140930" ! default
     call getin_p("phototable",phototable) ! default path
     write(*,*) " phototable = ",trim(phototable)


      fic = 81

      open(fic, form = 'formatted', status = 'old',                &
           file =trim(datadir)//"/"//trim(phototable),iostat=ierr)

      if (ierr /= 0) THEN
        write(*,*)'Error : cannot open photolysis lookup table ', trim(phototable)
        write(*,*)'It should be in :',trim(datadir),'/'
        write(*,*)'1) You can change this directory in callphys.def'
        write(*,*)'   with:'
        write(*,*)'   datadir=/path/to/the/directory'
        write(*,*)'2) You can change the input phototable file name in'
        write(*,*)'   callphys.def with:'
        write(*,*)'   phototable=filename'
        stop
      end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! read photolys table
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*, 'read photolysis lookup table ',trim(phototable)

      do itau = 1,ntau
         do itemp = 1,ntemp
            do iozo = 1,nozo
               do isza = 1,nsza
                  do iz = nz,1,-1
                     read(fic,*) colairtab(iz), xsza, table_ozo(iozo)
                     read(fic,'(7e11.4)') (jphot(itemp,isza,iz,iozo,itau,ij), ij= 1,nd)
                     do ij = 1,nd
                        if (jphot(itemp,isza,iz,iozo,itau,ij) == 1.e-30) then
                           jphot(itemp,isza,iz,iozo,itau,ij) = 0.
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do

      print*, 'lookup table...ok'
      close(fic)

      return
      end
