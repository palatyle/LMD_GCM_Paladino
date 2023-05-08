      subroutine surfacearea(ngrid, nlay, naerkind, ptimestep,
     $                       pplay, pzlay, 
     $                       pt, pq, pdq, nq, rdust, rice, tau, 
     $                       tauscaling,
     $                       surfdust, surfice)

      use tracer_mod, only: nuice_sed, igcm_dust_number,
     &                      igcm_ccn_number, varian, ccn_factor
      use conc_mod, only: rnew
      USE comcstfi_h
      implicit none

!==========================================================================
!     calculation of the ice and dust surface area (m2/m3)
!     available for heterogeneous reactions
!
!     Franck Lefevre
!     version 1.2 april 2012
!==========================================================================

      include "callkeys.h"

! input

      integer,intent(in) :: ngrid, nlay, naerkind
      integer,intent(in) :: nq               ! number of tracers
      real,intent(in) :: ptimestep           ! physics time step (s)
      real,intent(in) :: pplay(ngrid,nlay)   ! pressure at mid-layers (Pa)
      real,intent(in) :: pzlay(ngrid,nlay)   ! altitude at mid-layers (m)
      real,intent(in) :: pt(ngrid,nlay)      ! temperature at mid-layers (K)
      real,intent(in) :: pq(ngrid,nlay,nq)   ! tracers (kg/kg)
      real,intent(in) :: pdq(ngrid,nlay,nq)  ! physical tendency (kg/kg.s-1)
      real,intent(in) :: rdust(ngrid,nlay)   ! dust geometric mean radius (m)
      real,intent(in) :: rice(ngrid,nlay)    ! ice mass mean radius (m)
      real,intent(in) :: tau(ngrid,naerkind) ! column dust optical depth at each point
      real,intent(in) :: tauscaling(ngrid)   ! conversion factor for dust amount

! output

      real,intent(out) :: surfdust(ngrid,nlay) ! dust surface area (m2/m3)
      real,intent(out) :: surfice(ngrid,nlay)  ! water-ice surface area (m2/m3)

! local

      integer    :: l, ig
      real       :: rho                     ! density (kg/m3)
      real       :: dustnd, icend           ! uodated dust and ice number densities (kg/kg)
      real, save :: factor_ice, factor_dust ! multiplying factor to compute total surface area
                                            ! from the mass-mean radius
      real       :: sigma_ice, sigma_dust   ! variance of the ice and dust distributions
      real       :: ccntyp                  ! typical dust number density (#/kg)
                                            ! (microphys = false)
      real       :: rdusttyp                ! typical dust radius (m)
                                            ! (microphys = false)

      logical, save :: firstcall = .true.

!==========================================================================

      if (firstcall) then ! compute the multiplying factors
         sigma_dust  = varian
         sigma_ice   = sqrt(log(nuice_sed + 1.))
         factor_dust = exp(0.5*(log(sigma_dust))**2)
         factor_ice  = exp(0.5*(log(sigma_ice))**2)
         write(*,*) 'surfacearea : factor_dust = ', factor_dust 
         write(*,*) 'surfacearea : factor_ice  = ', factor_ice 
         firstcall = .false.
      end if

      if (microphys) then ! improvedclouds
         do l = 1,nlay
            do ig = 1,ngrid
!              atmospheric density
               rho = pplay(ig,l)/(rnew(ig,l)*pt(ig,l))
!              updated dust number density
               dustnd = pq(ig,l,igcm_dust_number)
     $                + pdq(ig,l,igcm_dust_number)*ptimestep
!              updated ice number density
               icend  = pq(ig,l,igcm_ccn_number)
     $                + pdq(ig,l,igcm_ccn_number)*ptimestep
!              dust surface area
               surfdust(ig,l) = factor_dust*dustnd*rho*tauscaling(ig)
     $                          *4.*pi*rdust(ig,l)**2
!              ice surface area
               surfice(ig,l)  = factor_ice*icend*rho*tauscaling(ig)
     $                          *4.*pi*rice(ig,l)**2
            end do
         end do
      else               ! simpleclouds
         do l = 1,nlay
            do ig = 1,ngrid
!              atmospheric density
               rho = pplay(ig,l)/(rnew(ig,l)*pt(ig,l))
!              typical dust radius
               rdusttyp = max(.8e-6*exp(-pzlay(ig,l)/18000.),1.e-9)
!              typical dust number density
               ccntyp = 1.3e+8*max(tau(ig,1),0.001)/0.1
     $                  *exp(-pzlay(ig,l)/10000.)
               ccntyp = ccntyp/ccn_factor
               if (rice(ig,l) .gt. rdust(ig,l)) then
                  surfdust(ig,l) = factor_dust*ccntyp*(ccn_factor - 1.)
     $                             *rho*4.*pi*rdusttyp**2
                  surfice(ig,l)  = factor_ice*ccntyp*4.*pi*rice(ig,l)**2
               else
                  surfdust(ig,l) = factor_dust*ccntyp*ccn_factor
     $                             *rho*4.*pi*rdusttyp**2
                  surfice(ig,l)  = 0.
               end if
            end do
         end do
      end if         ! of microphys

! write diagnostics in micron2/cm3
      
      if (callstats) then
        call wstats(ngrid,"surfdust", "Dust surface area",
     $            "micron2 cm-3",3,surfdust*1.e6)
        call wstats(ngrid,"surfice", "Ice cloud surface area",
     $            "micron2 cm-3",3,surfice*1.e6)
      endif
      call writediagfi(ngrid,"surfdust", "Dust surface area",
     $            "micron2 cm-3",3,surfdust*1.e6)
      call writediagfi(ngrid,"surfice", "Ice cloud surface area",
     $            "micron2 cm-3",3,surfice*1.e6)

      return
      end
