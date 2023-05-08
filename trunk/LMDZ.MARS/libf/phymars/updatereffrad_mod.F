      MODULE updatereffrad_mod
      
      IMPLICIT NONE
      
      CONTAINS
      
      SUBROUTINE updatereffrad(ngrid,nlayer,
     &                rdust,rstormdust,rtopdust,rice,nuice,
     &                reffrad,nueffrad,
     &                pq,tauscaling,tau,pplay)
       USE updaterad, ONLY: updaterdust, updaterice_micro,
     &                      updaterice_typ
       use tracer_mod, only: nqmx, igcm_dust_mass, igcm_dust_number,
     &                       igcm_h2o_ice, igcm_ccn_mass, radius,
     &                       igcm_ccn_number, nuice_ref, varian,
     &                       ref_r0, igcm_dust_submicron,
     &                       igcm_stormdust_mass,igcm_stormdust_number,
     &                       igcm_topdust_mass,igcm_topdust_number
       USE dimradmars_mod, only: nueffdust,naerkind,
     &            name_iaer,
     &            iaer_dust_conrath,iaer_dust_doubleq,
     &            iaer_dust_submicron,iaer_h2o_ice,
     &            iaer_stormdust_doubleq,iaer_topdust_doubleq 

       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Subroutine designed to update the aerosol size distribution used by
c     the radiative transfer scheme. This size distribution is assumed
c     to be a log-normal distribution, with effective radius "reffrad" and
c     variance "nueffrad".
c   At firstcall, "rice" and "nuice" are not known, because
c     the H2O ice microphysical scheme is called after the radiative
c     transfer in physiq.F. That's why we assess the size of the 
c     water-ice particles at firstcall (see part 1.2 below).
c
c   author:   
c   ------
c   J.-B. Madeleine (2009-2010)
c
c=======================================================================
c
c    Declarations :
c    -------------
c
      include "callkeys.h"

c-----------------------------------------------------------------------
c     Inputs/outputs:
c     ------

      INTEGER, INTENT(in) :: ngrid,nlayer
c     Ice geometric mean radius (m)
      REAL, INTENT(out) :: rice(ngrid,nlayer)
c     Estimated effective variance of the size distribution (n.u.)
      REAL, INTENT(out) :: nuice(ngrid,nlayer)
c     Tracer mass mixing ratio (kg/kg)
      REAL, INTENT(in) :: pq(ngrid,nlayer,nqmx)
      REAL, INTENT(out) :: rdust(ngrid,nlayer) ! Dust geometric mean radius (m)
      REAL, INTENT(out) :: rstormdust(ngrid,nlayer) ! Dust geometric mean radius (m)   
      REAL, INTENT(out) :: rtopdust(ngrid,nlayer) ! Dust geometric mean radius (m) 
      REAL, INTENT(in) :: pplay(ngrid,nlayer) ! altitude at the middle of the layers
      REAL, INTENT(in) :: tau(ngrid,naerkind)
c     Aerosol effective radius used for radiative transfer (meter)
      REAL, INTENT(out) :: reffrad(ngrid,nlayer,naerkind)
c     Aerosol effective variance used for radiative transfer (n.u.)
      REAL, INTENT(out) :: nueffrad(ngrid,nlayer,naerkind)
      REAL, INTENT(in) :: tauscaling(ngrid)         ! Convertion factor for qccn and Nccn
      
c     Local variables:
c     ---------------

      INTEGER :: ig,l          ! 3D grid indices
      INTEGER :: iaer          ! Aerosol index

c     Number of cloud condensation nuclei near the surface
c     (only used at firstcall). This value is taken from 
c     Montmessin et al. 2004 JGR 109 E10004 p5 (2E6 part m-3), and
c     converted to part kg-1 using a typical atmospheric density.

      REAL, PARAMETER :: ccn0 = 1.3E8
      
c     For microphysics only:      
      REAL Mo,No                       ! Mass and number of ccn
      REAL rhocloud(ngrid,nlayer)  ! Cloud density (kg.m-3)

      LOGICAL,SAVE :: firstcall=.true.

      REAL CBRT
      EXTERNAL CBRT

c==================================================================
c 1. Update radius from fields from dynamics or initial state
c==================================================================

c       1.1 Dust particles
c       ------------------
        IF (doubleq.AND.active) THEN
          DO l=1,nlayer
            DO ig=1, ngrid
              call updaterdust(pq(ig,l,igcm_dust_mass),
     &                         pq(ig,l,igcm_dust_number),rdust(ig,l))
              nueffdust(ig,l) = exp(varian**2.)-1.
             ENDDO
           ENDDO
        ELSE
          DO l=1,nlayer
            DO ig=1, ngrid
              rdust(ig,l) = 0.8E-6
              nueffdust(ig,l) = 0.3
            ENDDO
          ENDDO
        ENDIF

        ! updating radius of stormdust particles
        IF (rdstorm.AND.active) THEN
          DO l=1,nlayer
            DO ig=1, ngrid
              call updaterdust(pq(ig,l,igcm_stormdust_mass),
     &                 pq(ig,l,igcm_stormdust_number),rstormdust(ig,l))
              nueffdust(ig,l) = exp(varian**2.)-1.
             ENDDO
           ENDDO
        ENDIF

        ! updating radius of topdust particles
        IF (slpwind.AND.active) THEN
          DO l=1,nlayer
            DO ig=1, ngrid
              call updaterdust(pq(ig,l,igcm_topdust_mass),
     &                 pq(ig,l,igcm_topdust_number),rtopdust(ig,l))
              nueffdust(ig,l) = exp(varian**2.)-1.
             ENDDO
           ENDDO
        ENDIF
        
c       1.2 Water-ice particles
c       -----------------------

        IF (water.AND.activice) THEN 
         IF (microphys) THEN
        
c    At firstcall, the true number and true mass of cloud condensation nuclei are not known.
c    Indeed it is scaled on the prescribed dust opacity via a 'tauscaling' coefficient
c    computed after radiative transfer. If tauscaling is not in startfi, we make an assumption for its value.

          IF (firstcall) THEN
            !IF (minval(tauscaling).lt.0) tauscaling(:) = 1.e-3 ! default value when non-read in startfi is -1
            !IF (freedust)                tauscaling(:) = 1.    ! if freedust, enforce no rescaling at all
            firstcall = .false.
          ENDIF
 
          DO l=1,nlayer
            DO ig=1,ngrid
              call updaterice_micro(pq(ig,l,igcm_h2o_ice),
     &                              pq(ig,l,igcm_ccn_mass),
     &                              pq(ig,l,igcm_ccn_number),
     &                              tauscaling(ig),rice(ig,l),
     &                              rhocloud(ig,l))
              nuice(ig,l) = nuice_ref
            ENDDO
          ENDDO
          
        ELSE ! if not microphys
         
          DO l=1,nlayer
            DO ig=1,ngrid    
              call updaterice_typ(pq(ig,l,igcm_h2o_ice),
     &                          tau(ig,1),pplay(ig,l),rice(ig,l)) 
              nuice(ig,l) = nuice_ref
            ENDDO
          ENDDO
 
        ENDIF ! of if microphys
       ENDIF ! of if (water.AND.activice)

c==================================================================
c 2. Radius used in the radiative transfer code (reffrad)
c==================================================================

      DO iaer = 1, naerkind ! Loop on aerosol kind
        aerkind: SELECT CASE (name_iaer(iaer))
c==================================================================
        CASE("dust_conrath") aerkind         ! Typical dust profile
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rdust(ig,l) * 
     &          (1.e0 + nueffdust(ig,l))**2.5
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        CASE("dust_doubleq") aerkind! Two-moment scheme for dust
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rdust(ig,l) * ref_r0
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        CASE("dust_submicron") aerkind   ! Small dust population
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer)=radius(igcm_dust_submicron)
              nueffrad(ig,l,iaer)=0.03
            ENDDO
          ENDDO     
c==================================================================
        CASE("h2o_ice") aerkind             ! Water ice crystals
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
c             About reffice, do not confuse the mass mean radius
c             (rayon moyen massique) and the number median radius
c             (or geometric mean radius, rayon moyen géométrique).
c             rice is a mass mean radius, whereas rdust
c             is a geometric mean radius:
c             number median rad = mass mean rad x exp(-1.5 sigma0^2)
c             (Montmessin et al. 2004 paragraph 30). Therefore:
              reffrad(ig,l,iaer)=rice(ig,l)*(1.+nuice_ref)
              nueffrad(ig,l,iaer)=nuice_ref
            ENDDO
          ENDDO
c==================================================================
        CASE("stormdust_doubleq") aerkind! Two-moment scheme for
c       stormdust; same distribution than normal dust
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rstormdust(ig,l) * ref_r0
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        CASE("topdust_doubleq") aerkind! MV18: Two-moment scheme for
c       topdust; same distribution than normal dust
c==================================================================
          DO l=1,nlayer
            DO ig=1,ngrid
              reffrad(ig,l,iaer) = rtopdust(ig,l) * ref_r0
              nueffrad(ig,l,iaer) = nueffdust(ig,l)
            ENDDO
          ENDDO
c==================================================================
        END SELECT aerkind
      ENDDO ! iaer (loop on aerosol kind)

      END SUBROUTINE updatereffrad
      
      END MODULE updatereffrad_mod
