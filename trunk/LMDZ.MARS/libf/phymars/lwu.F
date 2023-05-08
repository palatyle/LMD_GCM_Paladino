      subroutine lwu (kdlon,kflev
     &                ,dp,plev,tlay,aerosol 
     &                ,QIRsQREF3d,omegaIR3d,gIR3d
     &                ,aer_t,co2_u,co2_up
     &                ,tautotal,omegtotal,gtotal)

c----------------------------------------------------------------------
c     LWU   computes  - co2: longwave effective absorber amounts including 
c                      pressure and temperature effects 
c                     - aerosols: amounts for every band
c                                transmission for bandes 1 and 2 of co2
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c ATTENTION AUX UNITES:
c le facteur 10*g fait passer des kg m-2 aux g cm-2
c-----------------------------------------------------------------------
c! modif diffusion
c! on ne change rien a la bande CO2 : les quantites d'absorbant CO2
c! sont multipliees par 1.66
c! pview= 1/cos(teta0)=1.66
c
c Modif J.-B. Madeleine: Computing optical properties of dust and
c   water-ice crystals in each gridbox. Optical parameters of
c   water-ice clouds are convolved to crystal sizes predicted by
c   the microphysical scheme.
c
c MODIF : FF : removing the monster bug on water ice clouds 11/2010
c
c MODIF : TN : corrected bug if very big water ice clouds 04/2012
c-----------------------------------------------------------------------
 
      use dimradmars_mod, only: ndlo2, nir, nuco2, ndlon, nflev
      use dimradmars_mod, only: naerkind
      use yomlw_h, only: nlaylte, tref, at, bt, cst_voigt
      use comcstfi_h, only: g
      implicit none

      include "callkeys.h"

c----------------------------------------------------------------------
c         0.1   arguments
c               ---------                                    
c                                                            inputs:
c                                                            -------
      integer,intent(in) :: kdlon ! part of ngrid
      integer,intent(in) :: kflev ! part of nalyer

      real,intent(in) :: dp(ndlo2,kflev) ! layer pressure thickness (Pa)
      real,intent(in) :: plev(ndlo2,kflev+1) ! level pressure (Pa)
      real,intent(in) :: tlay(ndlo2,kflev) ! layer temperature (K)
      real,intent(in) :: aerosol(ndlo2,kflev,naerkind) ! aerosol extinction optical depth
c                         at reference wavelength "longrefvis" set
c                         in dimradmars_mod , in each layer, for one of
c                         the "naerkind" kind of aerosol optical properties.
      real,intent(in) :: QIRsQREF3d(ndlo2,kflev,nir,naerkind) ! 3d ext. coef.
      real,intent(in) :: omegaIR3d(ndlo2,kflev,nir,naerkind) ! 3d ssa
      real,intent(in) :: gIR3d(ndlo2,kflev,nir,naerkind) ! 3d assym. param.

c                                                            outputs:
c                                                            --------
      real,intent(out) :: aer_t(ndlo2,nuco2,kflev+1)   ! transmission (aer) 
      real,intent(out) :: co2_u(ndlo2,nuco2,kflev+1)   ! absorber amounts (co2)
      real,intent(out) :: co2_up(ndlo2,nuco2,kflev+1)  ! idem scaled by the pressure (co2)

      real,intent(out) :: tautotal(ndlo2,kflev,nir)  ! \   Total single scattering
      real,intent(out) :: omegtotal(ndlo2,kflev,nir) !  >  properties (Addition of the
      real,intent(out) :: gtotal(ndlo2,kflev,nir)    ! /   NAERKIND aerosols properties)

c----------------------------------------------------------------------
c         0.2   local arrays
c               ------------

      integer jl,jk,jkl,ja,n
 
      real aer_a (ndlon,nir,nflev+1) ! absorber amounts (aer) ABSORPTION
      real co2c           ! co2 concentration (pa/pa)
      real pview          ! cosecant of viewing angle
      real pref           ! reference pressure (1013 mb = 101325 Pa)
      real tx,tx2
      real phi (ndlon,nuco2)
      real psi (ndlon,nuco2)
      real plev2 (ndlon,nflev+1)
      real zzz

      real ray,coefsize,coefsizew,coefsizeg 

c************************************************************************
c----------------------------------------------------------------------
c         0.3  Initialisation  
c               -------------

      pview = 1.66
      co2c = 0.95           
      pref = 101325.

      do jk=1,nlaylte+1
        do jl=1,kdlon
          plev2(jl,jk)=plev(jl,jk)*plev(jl,jk)
        enddo
      enddo

c----------------------------------------------------------------------
c  Computing TOTAL single scattering parameters by adding properties of
c  all the NAERKIND kind of aerosols in each IR band

      tautotal(:,:,:)=0
      omegtotal(:,:,:)=0
      gtotal(:,:,:)=0

      do n=1,naerkind
        do ja=1,nir
          do jk=1,nlaylte
            do jl = 1,kdlon
              tautotal(jl,jk,ja)=tautotal(jl,jk,ja) +
     &           QIRsQREF3d(jl,jk,ja,n)*aerosol(jl,jk,n)
              omegtotal(jl,jk,ja) =  omegtotal(jl,jk,ja) + 
     &           QIRsQREF3d(jl,jk,ja,n)*aerosol(jl,jk,n)*
     &           omegaIR3d(jl,jk,ja,n)
              gtotal(jl,jk,ja) =  gtotal(jl,jk,ja) + 
     &           QIRsQREF3d(jl,jk,ja,n)*aerosol(jl,jk,n)*
     &           omegaIR3d(jl,jk,ja,n)*gIR3d(jl,jk,ja,n)
            enddo
          enddo
        enddo
      enddo
      do ja=1,nir
        do jk=1,nlaylte
          do jl = 1,kdlon
             gtotal(jl,jk,ja)=gtotal(jl,jk,ja)/omegtotal(jl,jk,ja)
             omegtotal(jl,jk,ja)=omegtotal(jl,jk,ja)/tautotal(jl,jk,ja)
          enddo
        enddo
      enddo

c----------------------------------------------------------------------
c         1.0   cumulative (aerosol) amounts (for every band)
c               ----------------------------

      jk=nlaylte+1
      do ja=1,nir
        do jl = 1 , kdlon
          aer_a(jl,ja,jk)=0.
        enddo
      enddo

      do jk=1,nlaylte
        jkl=nlaylte+1-jk
        do ja=1,nir
          do jl=1,kdlon
             aer_a(jl,ja,jkl)=aer_a(jl,ja,jkl+1)+
     &           tautotal(jl,jkl,ja)*(1.-omegtotal(jl,jkl,ja))
          enddo
        enddo
      enddo                 

c----------------------------------------------------------------------
c         1.0   bands 1 and 2 of co2
c               --------------------

      jk=nlaylte+1
      do ja=1,nuco2
        do jl = 1 , kdlon
          co2_u(jl,ja,jk)=0.
          co2_up(jl,ja,jk)=0.
          aer_t(jl,ja,jk)=1.
        enddo
      enddo

      do jk=1,nlaylte
                       jkl=nlaylte+1-jk
        do ja=1,nuco2
          do jl=1,kdlon

c                 introduces temperature effects on absorber(co2) amounts
c                 -------------------------------------------------------
            tx = sign(min(abs(tlay(jl,jkl)-tref),70.)
     .         ,tlay(jl,jkl)-tref)
            tx2=tx*tx
            phi(jl,ja)=at(1,ja)*tx+bt(1,ja)*tx2
            psi(jl,ja)=at(2,ja)*tx+bt(2,ja)*tx2
            phi(jl,ja)=exp(phi(jl,ja)/cst_voigt(2,ja))
            psi(jl,ja)=exp(2.*psi(jl,ja))

c                                        cumulative absorber(co2) amounts
c                                        --------------------------------
            co2_u(jl,ja,jkl)=co2_u(jl,ja,jkl+1)
     .     +         pview/(10*g)*phi(jl,ja)*dp(jl,jkl)*co2c

            co2_up(jl,ja,jkl)=co2_up(jl,ja,jkl+1)
     .     +         pview/(10*g*2*pref)*psi(jl,ja)
     .     *        (plev2(jl,jkl)-plev2(jl,jkl+1))*co2c


c                                                  (aerosol) transmission
c                                                  ----------------------
c   on calcule directement les transmissions pour les aerosols.
c   on multiplie le Qext  par 1-omega dans la bande du CO2.
c   et pourquoi pas d'abord?  hourdin@lmd.ens.fr

c TN 04/12 : if very big water ice clouds, aer_t is strictly rounded 
c to zero in lower levels, which is a source of NaN
           !aer_t(jl,ja,jkl)=exp(-pview*aer_a(jl,ja,jkl))
           aer_t(jl,ja,jkl)=max(exp(-pview*aer_a(jl,ja,jkl)),1e-30) 
           
          
          enddo
        enddo
      enddo      
      
c----------------------------------------------------------------------
      end
