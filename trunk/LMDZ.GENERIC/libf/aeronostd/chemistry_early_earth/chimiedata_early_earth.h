!--------------------------------------------
!     data for photochemistry
!--------------------------------------------

!--------------------------------------------
!     dimensions of photolysis lookup table
!--------------------------------------------

      integer, parameter :: nd    = 29  ! species
      integer, parameter :: nz    = 143 ! altitude
      integer, parameter :: nozo  = 4   ! ozone
      integer, parameter :: nsza  = 21  ! solar zenith angle
      integer, parameter :: ntemp = 1   ! temperature
      integer, parameter :: ntau  = 1   ! dust
      integer, parameter :: nch4  = 9   ! ch4

!--------------------------------------------

      real, parameter :: kb = 1.3806e-23

      common/chimiedata_early_earth/jphot,colairtab,table_ozo

      real jphot(ntemp,nsza,nz,nozo,ntau,nch4,nd)
      real colairtab(nz)
      real szatab(nsza)
      real table_ozo(nozo)
      real tautab(ntau)
      real table_ch4(nch4)


!      data szatab/0.,  5., 10., 15., 20., 25.,                          &
!     &            30., 35., 40., 45., 50., 55.,                         &
!     &            60., 65., 70., 75., 80., 82.,                         &
!     &            84., 86., 88., 90., 91., 92.,                         &
!     &            93., 94., 95./

      data szatab/0., 10., 20., 30., 40., 50.,                          &
     &            60., 65., 70., 75., 80., 82.,                         &
     &            84., 86., 88., 90., 91., 92.,                         &
     &            93., 94., 95./                                        

!      data szatab/40./

!      data tautab/0., 0.2, 0.4, 0.6, 0.8, 1., 2., 4./
      data tautab/0./

!      data table_ch4/0.e-3, 1.e-3/
      data table_ch4/0., 1.e-6, 3.e-6, 1.e-5, 3.e-5, 1.e-4, 3.e-4,      &
     &            1.e-3, 3.e-3/
!--------------------------------------------
!     number of reactions in ASIS solver
!--------------------------------------------

!      integer, parameter :: nb_phot_max       = 36/full
      integer, parameter :: nb_phot_max       = 29

!      integer, parameter :: nb_reaction_3_max = 6
!     integer, parameter :: nb_reaction_3_max = 10/full
      integer, parameter :: nb_reaction_3_max = 10

!      integer, parameter :: nb_reaction_4_max = 30
!      integer, parameter :: nb_reaction_4_max = 110/full
!      integer, parameter :: nb_reaction_4_max = 101/no nitrogen
      integer, parameter :: nb_reaction_4_max = 101


