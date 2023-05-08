      subroutine init_chimie_B
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Scheme B
c
c   subject:
c   --------
c
c   Initialization of the photochemical scheme
c
c   VERSION: 30/04/2003
c
c   This is dependent of the chemistry version:
c   for scheme B its only purpose is reading the
c   photodissociation rates lookup table
c
c   Author:   Franck Lefevre for scheme B
c
c   Arguments:
c   ----------
c
c   The output variable is jphot and is put in common chimiedata.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
c
#include "dimensions.h"
#include "chimiedata.h"
#include "datafile.h"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     local:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer fic, ij, iozo, isza, itemp, iz, ierr
c
      real xsza
      integer lnblnk
      external lnblnk 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      fic = 81

      open(fic, form = 'formatted', status = 'old',
     $file =datafile(1:lnblnk(datafile))//'/jmars.20030707',iostat=ierr)

      IF (ierr.NE.0) THEN
        write(*,*)'Error : cannot open file jmars.20030707 '
        write(*,*)'(in aeronomars/init_chimie_B.F)'
        write(*,*)'It should be in :',datafile(1:lnblnk(datafile)),'/'
        write(*,*)'1) You can change this directory address in '
        write(*,*)'   file phymars/datafile.h'
        write(*,*)'2) If necessary jmars.20030707 (and other datafiles)'
        write(*,*)'   can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~forget/datagcm/datafile'
        STOP
      ENDIF

c
      print*, 'lecture table des J...'
c
      do itemp = 1,ntemp
         do iozo = 1,nozo
            do isza = 1,nsza
               do iz = 200,0,-1
                  read(fic,*) colairtab(iz), xsza, table_ozo(iozo)
                  read(fic,*) (jphot(itemp,isza,iz,iozo,ij), ij=1,7)
                  read(fic,*) (jphot(itemp,isza,iz,iozo,ij), ij=8,11)
                  do ij = 1,nd
                     if (jphot(itemp,isza,iz,iozo,ij) .eq. 1.e-30) then
                        jphot(itemp,isza,iz,iozo,ij) = 0.
                     end if
                  end do
               end do
            end do
         end do
      end do
c
      print*, 'lecture table des J...ok'
      close(fic)
c
      return
      end
