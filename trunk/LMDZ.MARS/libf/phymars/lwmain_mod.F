       MODULE lwmain_mod

       IMPLICIT NONE

       CONTAINS
       
       subroutine lwmain (ig0,icount,kdlon,kflev
     .                   ,dp,dt0,emis
     .                   ,plev,tlev,tlay,aerosol,coolrate
     .                   ,fluxground,fluxtop
     .                   ,netrad
     &                   ,QIRsQREF3d,omegaIR3d,gIR3d
     &                   ,co2ice)

c----------------------------------------------------------------------
c     LWMAIN     organizes the LTE longwave calculations
c     for layer 1 to layer "nlaylte" (stored in  "yomlw_h")
c----------------------------------------------------------------------

      use dimradmars_mod, only: ndlo2, nflev, nir, ndlon, nuco2
      use dimradmars_mod, only: naerkind
      use yomlw_h, only: nlaylte, xi
      implicit none
      
c     declarations
c     -------------   
      include "callkeys.h"
      include "comg1d.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------
c     inputs/outputs:
c     -------
      integer, intent(in) :: ig0
      integer, intent(in) :: icount
      integer, intent(in) :: kdlon            ! part of ngrid
      integer, intent(in) :: kflev            ! part of nlayer

      real, intent(in) :: dp (ndlo2,kflev)         ! layer pressure thickness (Pa)
      real, intent(in) :: dt0 (ndlo2)              ! surface temperature discontinuity (K)
      real, intent(in) :: emis (ndlo2)             ! surface emissivity
      real, intent(in) :: plev (ndlo2,kflev+1)     ! level pressure (Pa)
      real, intent(in) :: tlev (ndlo2,kflev+1)     ! level temperature (K)
      real, intent(in) :: tlay (ndlo2,kflev)       ! layer temperature (K)
      real, intent(in) :: aerosol(ndlo2,kflev,naerkind)      !  aerosol extinction optical
c                         depth at reference wavelength "longrefvis" set
c                         in dimradmars_mod , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.


c                                                            outputs:
c                                                            --------
      real, intent(out) :: coolrate(ndlo2,kflev)      ! cooling rate (K/s)
      real, intent(out) :: fluxground(ndlo2)          ! downward ground flux (W/m2)
      real, intent(out) :: fluxtop(ndlo2)             ! outgoing upward flux (W/m2) ("OLR")
      real, intent(out) :: netrad (ndlo2,kflev)       ! radiative budget (W/m2)
c     Aerosol optical properties
      real, intent(in) :: QIRsQREF3d(ndlo2,kflev,nir,naerkind)
      real, intent(in) :: omegaIR3d(ndlo2,kflev,nir,naerkind)
      real, intent(in) :: gIR3d(ndlo2,kflev,nir,naerkind)
      real, intent(in) :: co2ice(ndlo2)           ! co2 ice surface layer (kg.m-2)
c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      real aer_t (ndlon,nuco2,nflev+1)  ! transmission (aer)
      real co2_u (ndlon,nuco2,nflev+1)  ! absorber amounts (co2)
      real co2_up (ndlon,nuco2,nflev+1) ! idem scaled by the pressure (co2)

      real bsurf (ndlon,nir)            ! surface spectral planck function
      real btop (ndlon,nir)             ! top spectral planck function
      real blev (ndlon,nir,nflev+1)     ! level   spectral planck function
      real blay (ndlon,nir,nflev)       ! layer   spectral planck function
      real dblay (ndlon,nir,nflev)      ! layer gradient spectral planck function
      real dbsublay (ndlon,nir,2*nflev) ! layer gradient spectral planck function
                                        ! in sub layers

      real tautotal(ndlon,nflev,nir)  ! \   Total single scattering
      real omegtotal(ndlon,nflev,nir) !  >  properties (Addition of the
      real gtotal(ndlon,nflev,nir)    ! /   NAERKIND aerosols prop.)

      real newcoolrate(ndlon,nflev) ! cooling rate (K/s) / with implicite scheme

      real emis_gaz(ndlo2)         ! emissivity for gaz computations

      integer jk,jkk,ja,jl


c----------------------------------------------------------------------
c         0.3   Initialisation
c               --------------

      DO jl=1 , kdlon
         IF(co2ice(jl) .GT. 20.e-3) THEN
             emis_gaz(jl)=1. 
         ELSE
             emis_gaz(jl)=emis(jl)
         ENDIF
      ENDDO

c----------------------------------------------------------------------
c         1.0   planck function
c               ---------------

      call lwb ( kdlon, kflev, tlev, tlay, dt0
     .         , bsurf, btop, blay, blev, dblay, dbsublay)

c----------------------------------------------------------------------
c         2.0   absorber amounts
c               ----------------

      call lwu ( kdlon, kflev
     .         , dp, plev, tlay, aerosol
     &         , QIRsQREF3d,omegaIR3d,gIR3d
     .         , aer_t, co2_u, co2_up
     .         , tautotal,omegtotal,gtotal)

c----------------------------------------------------------------------
c         3.0   transmission functions / exchange coefficiants
c               ----------------------------------------------

c                                                                distants
c                                                                --------
                    if( mod(icount-1,ilwd).eq.0) then

c     print*, 'CALL of DISTANTS'
      call lwxd ( ig0, kdlon, kflev, emis_gaz
     .          , aer_t, co2_u, co2_up)

                    endif
c                                                              neighbours
c                                                              ----------
                    if( mod(icount-1,ilwn).eq.0) then

c     print*, 'CALL of NEIGHBOURS'
      call lwxn ( ig0, kdlon, kflev
     .          , dp
     .          , aer_t, co2_u, co2_up)

                    endif
c                                                              boundaries
c                                                              ----------
                    if( mod(icount-1,ilwb).eq.0) then

c     print*, 'CALL of BOUNDARIES'
      call lwxb ( ig0, kdlon, kflev, emis_gaz
     .          , aer_t, co2_u, co2_up)

                    endif

c----------------------------------------------------------------------
c         4.0   cooling rate
c               ------------

      call lwflux ( ig0, kdlon, kflev, dp
     .            , bsurf, btop, blev, blay, dbsublay
     .            , tlay, tlev, dt0      ! pour sortie dans g2d uniquement
     .            , emis
     .            , tautotal,omegtotal,gtotal
     .            , coolrate, fluxground, fluxtop
     .            , netrad)

c     do jk = 1, nlaylte
c       print*,coolrate(1,jk)
c     enddo
      
c       do jkk = 0 , nlaylte+1
c         do jk = 0 , nlaylte+1
c           do ja = 1 , nuco2
c             do jl = 1 , ngrid
c      if (xi (jl,ja,jk,jkk) .LT. 0
c    .       .OR. xi (jl,ja,jk,jkk) .GT. 1 ) then
c                 print*,'xi bande',ja,jk,jkk,xi (jl,ja,jk,jkk)
c      endif
c             enddo
c           enddo
c         enddo
c       enddo

c----------------------------------------------------------------------
c
c          5.    shema semi-implicite  (lwi) 
c                ---------------------------
c
c
      call lwi (ig0,kdlon,kflev,netrad,dblay,dp
     .          , newcoolrate)
c
c  Verif que   (X sol,space) + somme(X i,sol) = 1
c
      do jkk = 1 , nlaylte
        do jl = 1 , kdlon
c     print*,'NEW et OLD coolrate :',jkk,newcoolrate(jl,jkk) 
c    .  ,coolrate(jl,jkk)
      coolrate(jl,jkk) = newcoolrate(jl,jkk)
        enddo
      enddo
c
c----------------------------------------------------------------------

      END SUBROUTINE lwmain

      END MODULE lwmain_mod 

