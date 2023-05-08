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
