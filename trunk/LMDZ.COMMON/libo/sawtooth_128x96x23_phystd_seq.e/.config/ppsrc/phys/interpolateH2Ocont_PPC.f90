










     subroutine interpolateH2Ocont_PPC(wn,temp,presS,presF,abcoef,firstcall)

!==================================================================
!     
!     Purpose
!     -------
!     Calculates the H2O continuum opacity, using the formulae
!     provided in Pierrehumbert, PPC (2010). As this is based on
!     the CKD continuum, it provides a useful check for the 
!     implementation of the more general interpolateH2Ocont_CKD.F90.
!
!     Authors
!     -------
!     R. Wordsworth (2012)
!     
!==================================================================

      use watercommon_h, only: mH2O
      use datafile_mod, only: datadir
      implicit none

      ! input
      double precision wn                 ! wavenumber             (cm^-1)
      double precision temp               ! temperature            (Kelvin)
      double precision presS              ! self-pressure          (Pascals)
      double precision presF              ! foreign (air) pressure (Pascals)

      ! parameters
      double precision, parameter :: T0 = 296.0
      double precision, parameter :: p0 = 1.D+4

      ! variables
      double precision rho_w, x

      ! output
      double precision abcoef             ! absorption coefficient (m^-1)

      logical firstcall

      x = wn - 2500.

      if(firstcall)then ! called by sugas_corrk only
         print*,'----------------------------------------------------'
         print*,'Testing H2O continuum...'

         print*,'interpolateH2Ocont: At wavenumber ',wn,' cm^-1'
         print*,'   temperature ',temp,' K'
         print*,'   H2O pressure ',presS,' Pa'

         rho_w = presS/((8.31446/(mH2O/1000.))*temp)

         if(wn.gt.500 .and. wn.lt.1400)then
            abcoef = exp(12.167 - 0.050898*wn + 8.3207e-5*wn**2 - 7.0748e-8*wn**3 + 2.3261e-11*wn**4)*(T0/temp)**4.25*(presS/p0)
         elseif(wn.gt.2100 .and. wn.lt.3000)then
            abcoef = exp(-6.0055 - 0.0021363*x + 6.4723e-7*x**2 - 1.493e-8*x**3 + 2.5621e-11*x**4 + 7.328e-14*x**5)*(T0/temp)**4.25*(presS/p0)
         else
            abcoef = 0.0
         endif
         abcoef = abcoef*rho_w

         print*,'The self absorption is ',abcoef,' m^-1'
         print*,'And optical depth / km : ',1000.0*abcoef

      else

         rho_w = presS/((8.31446/(mH2O/1000.))*temp)

         if(wn.gt.500 .and. wn.lt.1400)then
            abcoef = exp(12.167 - 0.050898*wn + 8.3207e-5*wn**2 - 7.0748e-8*wn**3 + 2.3261e-11*wn**4)*(T0/temp)**4.25*(presS/p0)
         elseif(wn.gt.2100 .and. wn.lt.3000)then
            abcoef = exp(-6.0055 - 0.0021363*x + 6.4723e-7*x**2 - 1.493e-8*x**3 + 2.5621e-11*x**4 + 7.328e-14*x**5)*(T0/temp)**4.25*(presS/p0)
         else
            abcoef = 0.0
         endif
         abcoef = abcoef*rho_w

         ! unlike for Rayleigh scattering, we do not currently weight by the BB function
         ! however our bands are normally thin, so this is no big deal.

      endif

      return
    end subroutine interpolateH2Ocont_PPC

