










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

!--------------------------------------------
!     data for photochemistry
!--------------------------------------------

!--------------------------------------------
!     dimensions of photolysis lookup table
!--------------------------------------------

      integer, parameter :: nd    = 13  ! species
      integer, parameter :: nz    = 143 ! altitude
      integer, parameter :: nozo  = 7   ! ozone
      integer, parameter :: nsza  = 27  ! solar zenith angle
      integer, parameter :: ntemp = 4   ! temperature
      integer, parameter :: ntau  = 8   ! dust

!--------------------------------------------

      real, parameter :: kb = 1.3806e-23

      common/chimiedata/jphot,colairtab,table_ozo

      real jphot(ntemp,nsza,nz,nozo,ntau,nd)
      real colairtab(nz)
      real szatab(nsza)
      real table_ozo(nozo)
      real tautab(ntau)

      data szatab/0.,  5., 10., 15., 20., 25.,                          &
     &            30., 35., 40., 45., 50., 55.,                         &
     &            60., 65., 70., 75., 80., 82.,                         &
     &            84., 86., 88., 90., 91., 92.,                         &
     &            93., 94., 95./

      data tautab/0., 0.2, 0.4, 0.6, 0.8, 1., 2., 4./

!--------------------------------------------
!     number of reactions in ASIS solver
!--------------------------------------------

      integer, parameter :: nb_phot_max       = 18
      integer, parameter :: nb_reaction_3_max = 6
      integer, parameter :: nb_reaction_4_max = 30

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
