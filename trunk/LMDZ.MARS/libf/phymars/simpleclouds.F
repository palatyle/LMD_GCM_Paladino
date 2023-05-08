      subroutine simpleclouds(ngrid,nlay,ptimestep,
     &             pplay,pzlay,pt,pdt,
     &             pq,pdq,pdqcloud,pdtcloud,
     &             nq,tau,rice)
      USE updaterad
      USE watersat_mod, ONLY: watersat
      use tracer_mod, only: igcm_h2o_vap, igcm_h2o_ice
      USE comcstfi_h
      use dimradmars_mod, only: naerkind
      implicit none
c------------------------------------------------------------------
c  This routine is used to form clouds when a parcel of the GCM is
c    saturated. It is a simplified scheme, and there is almost no
c    microphysics involved. When the air is saturated, water-ice
c    clouds form on a fraction of the dust particles, specified by
c    the constant called "ccn_factor". There is no supersaturation,
c    and no nucleation rates computed. A more accurate scheme can
c    be found in the routine called "improvedclouds.F".

c  Modif de zq si saturation dans l'atmosphere
c  si zq(ig,l)> zqsat(ig,l) ->    zq(ig,l)=zqsat(ig,l)
c  Le test est effectue de bas en haut. L'eau condensee
c    (si saturation) est remise dans la couche en dessous.
c  L'eau condensee dans la couche du bas est deposee a la surface

c  Authors: Franck Montmessin (water ice scheme)
c           Francois Forget (changed nuclei density & outputs)
c           Ehouarn Millour (sept.2008, tracers are now handled
c                                   by name and not fixed index)
c           J.-B. Madeleine (developed a single routine called
c                            simpleclouds.F, and corrected calculations
c                            of the typical CCN profile, Oct. 2011)
c------------------------------------------------------------------
#include "callkeys.h"

c------------------------------------------------------------------
c     Arguments:
c     ---------
c     Inputs:
      INTEGER ngrid,nlay
      integer nq                 ! nombre de traceurs
      REAL ptimestep             ! pas de temps physique (s)
      REAL pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)
      REAL pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL pt(ngrid,nlay)        ! temperature at the middle of the
                                 !   layers (K)
      REAL pdt(ngrid,nlay)       ! tendance temperature des autres
                                 !   param.
      real pq(ngrid,nlay,nq)     ! traceur (kg/kg)
      real pdq(ngrid,nlay,nq)    ! tendance avant condensation
                                 !   (kg/kg.s-1)
      REAL tau(ngrid,naerkind)   ! Column dust optical depth at each point

c     Output:
      REAL rice(ngrid,nlay)      ! Ice mass mean radius (m)
                                 ! (r_c in montmessin_2004)
      real pdqcloud(ngrid,nlay,nq) ! tendance de la condensation
                                   !   H2O(kg/kg.s-1)
      REAL pdtcloud(ngrid,nlay)    ! tendance temperature due
                                   !   a la chaleur latente

c------------------------------------------------------------------
c     Local variables:

      REAL rhocloud(ngrid,nlay)  ! Cloud density (kg.m-3)

      INTEGER ig,l

      REAL zq(ngrid,nlay,nq)    ! local value of tracers
      REAL zq0(ngrid,nlay,nq)   ! local initial value of tracers
      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL zqsat(ngrid,nlay)    ! saturation
      REAL*8 dzq                      ! masse de glace echangee (kg/kg)
      REAL lw                         !Latent heat of sublimation (J.kg-1) 
      REAL,PARAMETER :: To=273.15     ! reference temperature, T=273.15 K
      real rdusttyp(ngrid,nlay) ! Typical dust geom. mean radius (m)
      REAL ccntyp(ngrid,nlay)
                                      ! Typical dust number density (#/kg)
c     CCN reduction factor
c      REAL, PARAMETER :: ccn_factor = 4.5  !! comme TESTS_JB // 1. avant
      

c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c    On "update" la valeur de q(nq) (water vapor) et temperature.
c    On effectue qqes calculs preliminaires sur les couches : 

      do l=1,nlay
        do ig=1,ngrid
          zq(ig,l,igcm_h2o_vap)=
     &      pq(ig,l,igcm_h2o_vap)+pdq(ig,l,igcm_h2o_vap)*ptimestep
          zq(ig,l,igcm_h2o_vap)=max(zq(ig,l,igcm_h2o_vap),1.E-30) ! FF 12/2004 
          zq0(ig,l,igcm_h2o_vap)=zq(ig,l,igcm_h2o_vap)
          zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep

          zq(ig,l,igcm_h2o_ice)=
     &      pq(ig,l,igcm_h2o_ice)+pdq(ig,l,igcm_h2o_ice)*ptimestep
          zq(ig,l,igcm_h2o_ice)=max(zq(ig,l,igcm_h2o_ice),0.) ! FF 12/2004 
          zq0(ig,l,igcm_h2o_ice)=zq(ig,l,igcm_h2o_ice)
        enddo
      enddo


      pdqcloud(1:ngrid,1:nlay,1:nq)=0
      pdtcloud(1:ngrid,1:nlay)=0

c     ----------------------------------------------
c
c
c     Rapport de melange a saturation dans la couche l : -------
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call watersat(ngrid*nlay,zt,pplay,zqsat)

c     taux de condensation (kg/kg/s-1) dans les differentes couches
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do l=1,nlay
        do ig=1,ngrid

          if (zq(ig,l,igcm_h2o_vap).ge.zqsat(ig,l))then  !  Condensation
            dzq=zq(ig,l,igcm_h2o_vap)-zqsat(ig,l)                
          elseif(zq(ig,l,igcm_h2o_vap).lt.zqsat(ig,l))then  ! Sublimation
            dzq=-min(zqsat(ig,l)-zq(ig,l,igcm_h2o_vap),
     &               zq(ig,l,igcm_h2o_ice))
          endif

c         Water Mass change
c         ~~~~~~~~~~~~~~~~~
          zq(ig,l,igcm_h2o_ice)=zq(ig,l,igcm_h2o_ice)+dzq
          zq(ig,l,igcm_h2o_vap)=zq(ig,l,igcm_h2o_vap)-dzq
          

        enddo ! of do ig=1,ngrid
      enddo ! of do l=1,nlay

c     Tendance finale
c     ~~~~~~~~~~~~~~~
      do l=1, nlay
        do ig=1,ngrid
          pdqcloud(ig,l,igcm_h2o_vap)=(zq(ig,l,igcm_h2o_vap)
     &                            -zq0(ig,l,igcm_h2o_vap))/ptimestep
          pdqcloud(ig,l,igcm_h2o_ice) =
     &      (zq(ig,l,igcm_h2o_ice) - zq0(ig,l,igcm_h2o_ice))/ptimestep
          lw=(2834.3-0.28*(zt(ig,l)-To)-0.004*(zt(ig,l)-To)**2)*1.e+3
          pdtcloud(ig,l)=-pdqcloud(ig,l,igcm_h2o_vap)*lw/cpp
        end do
      end do

c     ice crystal radius
      do l=1, nlay
        do ig=1,ngrid
          call updaterice_typ(zq(ig,l,igcm_h2o_ice),
     &       tau(ig,1),pzlay(ig,l),rice(ig,l))
        end do
      end do

c------------------------------------------------------------------
      return
      end
