      MODULE rocketduststorm_mod

      IMPLICIT NONE

      REAL, SAVE, ALLOCATABLE :: dustliftday(:) ! dust lifting rate (s-1)
      
      CONTAINS

!=======================================================================
! ROCKET DUST STORM - vertical transport and detrainment 
!=======================================================================
! calculation of the vertical flux
! call of van_leer : Van Leer transport scheme of the dust tracers
! detrainement of stormdust in the background dust
! -----------------------------------------------------------------------
! Authors: M. Vals; C. Wang; F. Forget; T. Bertrand 
! Institution: Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------

      SUBROUTINE rocketduststorm(ngrid,nlayer,nq,ptime,ptimestep,      &
                                 pq,pdqfi,pt,pdtfi,pplev,pplay,pzlev,  &
                                 pzlay,pdtsw,pdtlw,                    &
!             input for radiative transfer
                                 clearatm,icount,zday,zls,             &
                                 tsurf,igout,totstormfract,            &
                                 tauscaling,                           &
!             input sub-grid scale cloud
                                 clearsky,totcloudfrac,                &
!             input sub-grid scale topography
                                 nohmons,alpha_hmons,                  &
!             output
                                 pdqrds,wrad,dsodust,dsords,dsotop,    &
                                 tauref)

      USE tracer_mod, only: igcm_stormdust_mass,igcm_stormdust_number &
                            ,igcm_dust_mass,igcm_dust_number          &
                            ,rho_dust 
      USE comcstfi_h, only: r,g,cpp,rcp
      USE dimradmars_mod, only: albedo,naerkind
      USE comsaison_h, only: dist_sol,mu0,fract
      USE surfdat_h, only: emis,co2ice,zmea, zstd, zsig, hmons
      USE callradite_mod, only: callradite
      IMPLICIT NONE

      include "callkeys.h"

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------

      INTEGER, INTENT(IN) :: ngrid ! number of horizontal grid points
      INTEGER, INTENT(IN) :: nlayer ! number of vertical grid points
      INTEGER, INTENT(IN) :: nq ! number of tracer species
      REAL, INTENT(IN) :: ptime
      REAL, INTENT(IN) :: ptimestep

      REAL, INTENT(IN) :: pq(ngrid,nlayer,nq) ! advected field nq
      REAL, INTENT(IN) :: pdqfi(ngrid,nlayer,nq)! tendancy field pq
      REAL, INTENT(IN) :: pt(ngrid,nlayer)      ! temperature at mid-layer (K)
      REAL, INTENT(IN) :: pdtfi(ngrid,nlayer)   ! tendancy temperature at mid-layer (K)

      REAL, INTENT(IN) :: pplay(ngrid,nlayer)     ! pressure at middle of the layers
      REAL, INTENT(IN) :: pplev(ngrid,nlayer+1) ! pressure at intermediate levels
      REAL, INTENT(IN) :: pzlay(ngrid,nlayer)     ! altitude at the middle of the layers
      REAL, INTENT(IN) :: pzlev(ngrid,nlayer+1)   ! altitude at layer boundaries

      REAL, INTENT(IN) :: pdtsw(ngrid,nlayer)     ! (K/s) env
      REAL, INTENT(IN) :: pdtlw(ngrid,nlayer)     ! (K/s) env

!     input for second radiative transfer 
      LOGICAL, INTENT(IN) :: clearatm
      INTEGER, INTENT(INOUT) :: icount
      REAL, INTENT(IN) :: zday
      REAL, INTENT(IN) :: zls
      REAL, INTENT(IN) :: tsurf(ngrid)
      INTEGER, INTENT(IN) :: igout
      REAL, INTENT(IN) :: totstormfract(ngrid)
      REAL, INTENT(INOUT) :: tauscaling(ngrid)
      
!     sbgrid scale water ice clouds
      logical, intent(in) :: clearsky
      real, intent(in) :: totcloudfrac(ngrid)

!     sbgrid scale topography
      LOGICAL, INTENT(IN) :: nohmons
      REAL, INTENT(IN) :: alpha_hmons(ngrid)   
 
!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------

      REAL, INTENT(OUT) :: pdqrds(ngrid,nlayer,nq) ! tendancy field for dust when detraining
      REAL, INTENT(OUT) :: wrad(ngrid,nlayer+1)   ! vertical speed within the rocket dust storm
      REAL, INTENT(INOUT) :: dsodust(ngrid,nlayer) ! density scaled opacity of env. dust
      REAL, INTENT(INOUT) :: dsords(ngrid,nlayer) ! density scaled opacity of storm dust
      REAL, INTENT(INOUT) :: dsotop(ngrid,nlayer) ! density scaled opacity of topmons dust
      REAL, INTENT(OUT) :: tauref(ngrid)

!--------------------------------------------------------
! Local variables 
!--------------------------------------------------------
      INTEGER l,ig,iq,ll
!     local variables from callradite.F       
      REAL zdtlw1(ngrid,nlayer)    ! (K/s) storm
      REAL zdtsw1(ngrid,nlayer)    ! (K/s) storm
      REAL zt(ngrid,nlayer)        ! actual temperature at mid-layer (K)
      REAL zdtvert(ngrid,nlayer)   ! dT/dz , lapse rate
      REAL ztlev(ngrid,nlayer)     ! temperature at intermediate levels l+1/2 without last level

      REAL zdtlw1_lev(nlayer),zdtsw1_lev(nlayer) ! rad. heating rate at intermediate levels l+1/2 for stormdust
      REAL zdtlw_lev(nlayer),zdtsw_lev(nlayer)   ! rad. heating rate at intermediate levels l+1/2 for background dust

      REAL zq_stormdust_mass(ngrid,nlayer) ! intermediate tracer stormdust mass 
      REAL zq_stormdust_number(ngrid,nlayer) ! intermediate tracer stormdust number
      REAL zq_dust_mass(ngrid,nlayer)           ! intermediate tracer dust mass 
      REAL zq_dust_number(ngrid,nlayer)         ! intermediate tracer dust number

      REAL mr_stormdust_mass(ngrid,nlayer) ! intermediate mixing ratio to calculate van leer transport with the "real" concentration (stormdust mass)
      REAL mr_stormdust_number(ngrid,nlayer) ! intermediate mixing ratio to calculate van leer transport with the "real" concentration (stormdust number)
      REAL mr_dust_mass(ngrid,nlayer) ! intermediate mixing ratio to calculate van leer transport with the "real" concentration (dust mass)
      REAL mr_dust_number(ngrid,nlayer) ! intermediate mixing ratio to calculate van leer transport with the "real" concentration (sdust number)
                   
      REAL dqvl_stormdust_mass(ngrid,nlayer)    ! tendancy of vertical transport (stormdust mass)
      REAL dqvl_stormdust_number(ngrid,nlayer)  ! tendancy of vertical transport (stormdust number)
      REAL dqvl_dust_mass(ngrid,nlayer)    ! tendancy of vertical transport (dust mass)
      REAL dqvl_dust_number(ngrid,nlayer)  ! tendancy of vertical transport (dust number)
      REAL dqdet_stormdust_mass(ngrid,nlayer)   ! tendancy of detrainement (stormdust mass)
      REAL dqdet_stormdust_number(ngrid,nlayer) ! tendancy of detrainement (stormdust number)

      REAL masse_col(nlayer)     ! mass of atmosphere (kg/m2)
      REAL zq(ngrid,nlayer,nq)   ! updated tracers
      
      REAL w(ngrid,nlayer)          ! air mass flux (calculated with the vertical wind velocity profile) used as input in Van Leer (kgair/m2)
      REAL wqmass(ngrid,nlayer+1)   ! tracer (dust_mass) mass flux in Van Leer (kg/m2)
      REAL wqnumber(ngrid,nlayer+1) ! tracer (dust_number) mass flux in Van Leer (kg/m2)

      LOGICAL storm(ngrid)    ! true when there is a dust storm (if the opacity is high): trigger the rocket dust storm scheme
      REAL detrain(ngrid,nlayer)  ! coefficient for detrainment : % of stormdust detrained
      INTEGER scheme(ngrid)   ! triggered scheme
            
      REAL,PARAMETER:: coefmin =0.025 ! 0<coefmin<1 Minimum fraction of stormdust detrained  
      REAL,PARAMETER:: wmin =0.25 ! stormdust detrainment if wrad < wmin  
      REAL,PARAMETER:: wmax =10.   ! maximum vertical velocity of the rocket dust storms (m/s)

!     subtimestep
      INTEGER tsub
      INTEGER nsubtimestep    !number of subtimestep when calling van_leer 
      REAL subtimestep        !ptimestep/nsubtimestep
      REAL dtmax              !considered time needed for dust to cross one layer
      REAL,PARAMETER:: secu=3.!3.      !coefficient on wspeed to avoid dust crossing many layers during subtimestep

!     diagnostics 
      REAL lapserate(ngrid,nlayer) 
      REAL deltahr(ngrid,nlayer+1)
    
      LOGICAL,SAVE :: firstcall=.true.

!     variables for the radiative transfer
      REAL  fluxsurf_lw1(ngrid)
      REAL  fluxsurf_sw1(ngrid,2)
      REAL  fluxtop_lw1(ngrid)
      REAL  fluxtop_sw1(ngrid,2)
      REAL  tau(ngrid,naerkind)
      REAL  aerosol(ngrid,nlayer,naerkind)
      REAL  taucloudtes(ngrid)
      REAL  rdust(ngrid,nlayer)
      REAL  rstormdust(ngrid,nlayer)
      REAL  rtopdust(ngrid,nlayer)
      REAL  rice(ngrid,nlayer)
      REAL  nuice(ngrid,nlayer)


      ! **********************************************************************
      ! **********************************************************************
      ! Rocket dust storm parametrization to reproduce the detached dust layers
      ! during the dust storm season:
      !     The radiative warming due to the presence of storm dust is 
      !     balanced by the adiabatic cooling. The tracer "storm dust"  
      !     is transported by the upward/downward flow. 
      ! **********************************************************************
      ! **********************************************************************      
      !!    1. Radiative transfer in storm dust
      !!    2. Compute vertical velocity for storm dust
      !!      case 1 storm = false: nothing to do
      !!      case 2 rocket dust storm (storm=true)
      !!    3. Vertical transport (Van Leer)
      !!    4. Detrainment
      ! **********************************************************************
      ! **********************************************************************

      
      ! **********************************************************************
      ! Initializations
      ! **********************************************************************
      storm(:)=.false.
      pdqrds(:,:,:) = 0.
      mr_dust_mass(:,:)=0.
      mr_dust_number(:,:)=0.
      mr_stormdust_mass(:,:)=0.
      mr_stormdust_number(:,:)=0.
      dqvl_dust_mass(:,:)=0.
      dqvl_dust_number(:,:)=0.
      dqvl_stormdust_mass(:,:)=0.
      dqvl_stormdust_number(:,:)=0.
      dqdet_stormdust_mass(:,:)=0.
      dqdet_stormdust_number(:,:)=0.
      wrad(:,:)=0.
      w(:,:)=0.
      wqmass(:,:)=0.
      wqnumber(:,:)=0.
      zdtvert(:,:)=0.
      lapserate(:,:)=0.
      deltahr(:,:)=0.
      scheme(:)=0
      detrain(:,:)=1.

      !! no update for the stormdust tracer and temperature fields 
      !! because previous callradite was for background dust
      zq(1:ngrid,1:nlayer,1:nq)=pq(1:ngrid,1:nlayer,1:nq)
      zt(1:ngrid,1:nlayer)=pt(1:ngrid,1:nlayer)


      zq_dust_mass(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_dust_mass)
      zq_dust_number(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_dust_number)
      zq_stormdust_mass(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_stormdust_mass)
      zq_stormdust_number(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_stormdust_number)

      ! *********************************************************************
      ! 0. Check if there is a storm
      ! *********************************************************************
      DO ig=1,ngrid
        storm(ig)=.false.
        DO l=1,nlayer
          IF (zq(ig,l,igcm_stormdust_mass) &
          .gt. zq(ig,l,igcm_dust_mass)*(1.E-4)) THEN
            storm(ig)=.true.
            EXIT
          ENDIF
        ENDDO     
      ENDDO
      
      ! *********************************************************************
      ! 1. Call the second radiative transfer for stormdust, obtain the extra heating
      ! *********************************************************************
      CALL callradite(icount,ngrid,nlayer,nq,zday,zls,pq,albedo,          &
                 emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,      &
                 zdtlw1,zdtsw1,fluxsurf_lw1,fluxsurf_sw1,fluxtop_lw1,     &
                 fluxtop_sw1,tauref,tau,aerosol,dsodust,tauscaling,       &
                 taucloudtes,rdust,rice,nuice,co2ice,rstormdust,rtopdust, &
                 totstormfract,clearatm,dsords,dsotop,alpha_hmons,nohmons,&
                 clearsky,totcloudfrac)

      ! **********************************************************************
      ! 2. Compute vertical velocity for storm dust
      ! **********************************************************************
        !! **********************************************************************
        !! 2.1 Nothing to do when no storm
        !!             no storm
        DO ig=1,ngrid        
          IF (.NOT.(storm(ig))) then
            scheme(ig)=1
            cycle
          ENDIF ! IF (storm(ig))
        ENDDO ! DO ig=1,ngrid                 
        
        !! **********************************************************************
        !! 2.2 Calculation of the extra heating : computing heating rates
        !! gradient at boundaries of each layer, start from surface
        DO ig=1,ngrid
          IF (storm(ig)) THEN

            scheme(ig)=2
      
            !! computing heating rates gradient at boundraies of each layer
            !! start from surface
            zdtlw1_lev(1)=0.
            zdtsw1_lev(1)=0.
            zdtlw_lev(1)=0.
            zdtsw_lev(1)=0.
            ztlev(ig,1)=zt(ig,1)

            DO l=1,nlayer-1
              !! Calculation for the dust storm fraction
              zdtlw1_lev(l+1)=(zdtlw1(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+ &
                           zdtlw1(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /  &
                              (pzlay(ig,l+1)-pzlay(ig,l))
            
              zdtsw1_lev(l+1)=(zdtsw1(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+ &
                           zdtsw1(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /  &
                              (pzlay(ig,l+1)-pzlay(ig,l))
              !! Calculation for the background dust fraction
              zdtlw_lev(l+1)=(pdtlw(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+   &
                           pdtlw(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /   &
                              (pzlay(ig,l+1)-pzlay(ig,l))
           
              zdtsw_lev(l+1)=(pdtsw(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+   &
                           pdtsw(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /   &
                              (pzlay(ig,l+1)-pzlay(ig,l))
           
              ztlev(ig,l+1)=(zt(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+          &
                           zt(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /      &
                              (pzlay(ig,l+1)-pzlay(ig,l))
            ENDDO ! DO l=1,nlayer-1

            !! This is the env. lapse rate 
            zdtvert(ig,1)=0.
            DO l=1,nlayer-1
              zdtvert(ig,l+1)=(ztlev(ig,l+1)-ztlev(ig,l))/(pzlay(ig,l+1)-pzlay(ig,l))
              lapserate(ig,l+1)=zdtvert(ig,l+1)
            ENDDO

            !! **********************************************************************
            !! 2.3 Calculation of the vertical velocity : extra heating
            !! balanced by adiabatic cooling 
            
            DO l=1,nlayer
              deltahr(ig,l)=(zdtlw1_lev(l)+zdtsw1_lev(l))  &
                                          -(zdtlw_lev(l)+zdtsw_lev(l))
              wrad(ig,l)=-deltahr(ig,l)/(g/cpp+   &
                                         max(zdtvert(ig,l),-0.99*g/cpp))       
              !! Limit vertical wind in case of lapse rate close to adiabatic
              wrad(ig,l)=max(wrad(ig,l),-wmax) 
              wrad(ig,l)=min(wrad(ig,l),wmax)
            ENDDO

          ENDIF ! IF (storm(ig))
        ENDDO ! DO ig=1,ngrid

      ! **********************************************************************
      ! 3. Vertical transport
      ! **********************************************************************
        !! **********************************************************************
        !! 3.1 Transport of background dust + storm dust (concentrated)
        DO ig=1,ngrid
          IF (storm(ig)) THEN 
            DO l=1,nlayer
              mr_dust_mass(ig,l) = zq_dust_mass(ig,l)
              mr_dust_number(ig,l) = zq_dust_number(ig,l)
              mr_stormdust_mass(ig,l) = zq_dust_mass(ig,l)+ &
                                      zq_stormdust_mass(ig,l)/totstormfract(ig)
              mr_stormdust_number(ig,l) = zq_dust_number(ig,l)+ &
                                      zq_stormdust_number(ig,l)/totstormfract(ig)
            ENDDO ! DO l=1,nlayer
          ENDIF ! IF (storm(ig))
        ENDDO ! DO ig=1,ngrid

        DO ig=1,ngrid
          IF (storm(ig)) THEN
        !! **********************************************************************
        !! 3.2 Compute the subtimestep to conserve the mass in the Van Leer transport
            dtmax=ptimestep
            DO l=2,nlayer
              IF (wrad(ig,l).lt.0.) THEN
                 dtmax=min(dtmax,(pzlev(ig,l)-pzlev(ig,l-1))/  &
                                   (secu*abs(wrad(ig,l))))
              ELSE IF (wrad(ig,l).gt.0.) then
                 dtmax=min(dtmax,(pzlev(ig,l+1)-pzlev(ig,l))/  &
                                   (secu*abs(wrad(ig,l))))
              ENDIF
            ENDDO
            nsubtimestep= int(ptimestep/dtmax) 
            subtimestep=ptimestep/float(nsubtimestep)
            !! Mass flux generated by wup in kg/m2
            DO l=1,nlayer
               w(ig,l)=wrad(ig,l)*pplev(ig,l)/(r*ztlev(ig,l)) &
                        *subtimestep
            ENDDO ! l=1,nlayer

        !! **********************************************************************
        !! 3.3 Vertical transport by a Van Leer scheme
            !! Mass of atmosphere in the layer
            DO l=1,nlayer
               masse_col(l)=(pplev(ig,l)-pplev(ig,l+1))/g
            ENDDO
            !! Mass flux in kg/m2 if you are not using the subtimestep
            !DO l=1,nlayer
            !   w(ig,l)=wrad(ig,l)*(pplev(ig,l)/(r*ztlev(ig,l)))*ptimestep
            !ENDDO
            !! Loop over the subtimestep
            DO tsub=1,nsubtimestep
              !! Van Leer scheme
              wqmass(ig,:)=0.
              wqnumber(ig,:)=0.
              CALL van_leer(nlayer,mr_stormdust_mass(ig,:),2.,   &
                    masse_col,w(ig,:),wqmass(ig,:))
              CALL van_leer(nlayer,mr_stormdust_number(ig,:),2.,  &
                    masse_col,w(ig,:),wqnumber(ig,:))
            ENDDO !tsub=...           
             
          ENDIF ! IF storm(ig)
        ENDDO ! DO ig=1,ngrid  

        !! **********************************************************************
        !! 3.4 Re-calculation of the mixing ratios after vertical transport
        DO ig=1,ngrid
         IF (storm(ig)) THEN
           DO l=1,nlayer
           
             !! General and "healthy" case
             IF (mr_stormdust_mass(ig,l).ge.mr_dust_mass(ig,l)) THEN
               zq_dust_mass(ig,l) = mr_dust_mass(ig,l)
               zq_dust_number(ig,l) = mr_dust_number(ig,l)
               zq_stormdust_mass(ig,l) = totstormfract(ig)*(mr_stormdust_mass(ig,l)-mr_dust_mass(ig,l))
               zq_stormdust_number(ig,l) = totstormfract(ig)*(mr_stormdust_number(ig,l)-mr_dust_number(ig,l))
             !! Particular case: there is less than initial dust mixing ratio after the vertical transport
             ELSE
               zq_dust_mass(ig,l) = (1.-totstormfract(ig))*mr_dust_mass(ig,l)+totstormfract(ig)*mr_stormdust_mass(ig,l)
               zq_dust_number(ig,l) = (1.-totstormfract(ig))*mr_dust_number(ig,l)+totstormfract(ig)*mr_stormdust_number(ig,l)
               zq_stormdust_mass(ig,l) = 0.
               zq_stormdust_number(ig,l) = 0.
             ENDIF
             
           ENDDO ! DO l=1,nlayer            
         ENDIF ! IF storm(ig)
        ENDDO ! DO ig=1,ngrid

        !! **********************************************************************
        !! 3.5 Calculation of the tendencies of the vertical transport
        DO ig=1,ngrid
         IF (storm(ig)) THEN
           DO l=1,nlayer
             dqvl_stormdust_mass(ig,l) = (zq_stormdust_mass(ig,l)- &
                                   zq(ig,l,igcm_stormdust_mass)) /ptimestep
             dqvl_stormdust_number(ig,l) = (zq_stormdust_number(ig,l)- &
                                     zq(ig,l,igcm_stormdust_number)) /ptimestep
             dqvl_dust_mass(ig,l) = (zq_dust_mass(ig,l)-zq(ig,l,igcm_dust_mass)) /ptimestep
             dqvl_dust_number(ig,l) = (zq_dust_number(ig,l)-zq(ig,l,igcm_dust_number)) /ptimestep
           ENDDO
         ENDIF ! IF storm(ig)
        ENDDO ! DO ig=1,ngrid           

      ! **********************************************************************
      ! 4. Detrainment: stormdust is converted to background dust
      ! **********************************************************************
        !! **********************************************************************
        !! 4.1 Compute the coefficient of detrainmen
        DO ig=1,ngrid
          DO l=1,nlayer
            IF ((max(abs(wrad(ig,l)),abs(wrad(ig,l+1))) .lt.  & 
                          wmin) .or. (zq_dust_mass(ig,l) .gt.  &
                                     10000.*zq_stormdust_mass(ig,l))) THEN
               detrain(ig,l)=1.
            ELSE IF (max(abs(wrad(ig,l)),abs(wrad(ig,l+1)))   &
                                                       .le. wmax) THEN
               detrain(ig,l)=coeff_detrainment*                 &
                             (((1-coefmin)/(wmin-wmax)**2)*     &
                             (max(abs(wrad(ig,l)),abs(wrad(ig,l+1)))-wmax)**2 &
                              +coefmin)
            ELSE IF (max(abs(wrad(ig,l)),abs(wrad(ig,l+1))).gt. wmax ) THEN
               detrain(ig,l)=coefmin
            ELSE
               detrain(ig,l)=coefmin
            ENDIF
          ENDDO ! DO l=1,nlayer
        ENDDO ! DO ig=1,ngrid
        
        !! **********************************************************************
        !! 4.2 Calculation of the tendencies of the detrainment
        DO ig=1,ngrid
          DO l=1,nlayer
            dqdet_stormdust_mass(ig,l)=-detrain(ig,l)*zq_stormdust_mass(ig,l) &
                                                        /ptimestep
            dqdet_stormdust_number(ig,l)=-detrain(ig,l)*zq_stormdust_number(ig,l) &
                                                        /ptimestep
          ENDDO ! DO l=1,nlayer
        ENDDO ! DO ig=1,ngrid
           
      ! **********************************************************************
      ! 5. Final tendencies: vertical transport + detrainment
      ! **********************************************************************
        DO ig=1,ngrid
          DO l=1,nlayer      
          pdqrds(ig,l,igcm_stormdust_mass)=dqdet_stormdust_mass(ig,l) &
                                                 +dqvl_stormdust_mass(ig,l)
          pdqrds(ig,l,igcm_stormdust_number)=dqdet_stormdust_number(ig,l) &
                                                 +dqvl_stormdust_number(ig,l)
          pdqrds(ig,l,igcm_dust_mass)= -dqdet_stormdust_mass(ig,l) &
                                       +dqvl_dust_mass(ig,l) 
          pdqrds(ig,l,igcm_dust_number)= -dqdet_stormdust_number(ig,l) &
                                       +dqvl_dust_number(ig,l) 
          ENDDO ! DO l=1,nlayer
        ENDDO ! DO ig=1,ngrid

!      ! **********************************************************************
!      ! 6. To prevent from negative values
!      ! **********************************************************************
!        DO ig=1,ngrid
!          DO l=1,nlayer
!            IF ((pq(ig,l,igcm_stormdust_mass) &
!               +pdqrds(ig,l,igcm_stormdust_mass)*ptimestep .le. 0.) .or. &
!              (pq(ig,l,igcm_stormdust_number) &
!               +pdqrds(ig,l,igcm_stormdust_number)*ptimestep .le. 0.)) THEN
!               pdqrds(ig,l,igcm_stormdust_mass)=-pq(ig,l,igcm_stormdust_mass)/ptimestep
!               pdqrds(ig,l,igcm_stormdust_number)=-pq(ig,l,igcm_stormdust_number)/ptimestep
!            ENDIF
!          ENDDO ! nlayer
!        ENDDO ! DO ig=1,ngrid
!
!        DO ig=1,ngrid
!          DO l=1,nlayer           
!            IF ((pq(ig,l,igcm_dust_mass) &
!               +pdqrds(ig,l,igcm_dust_mass)*ptimestep .le. 0.) .or. &
!              (pq(ig,l,igcm_dust_number) &
!               +pdqrds(ig,l,igcm_dust_number)*ptimestep .le. 0.)) THEN
!               pdqrds(ig,l,igcm_dust_mass)=-pq(ig,l,igcm_dust_mass)/ptimestep
!               pdqrds(ig,l,igcm_dust_number)=-pq(ig,l,igcm_dust_number)/ptimestep
!            ENDIF
!          ENDDO ! nlayer
!        ENDDO ! DO ig=1,ngrid
        
!=======================================================================
! WRITEDIAGFI

        call WRITEDIAGFI(ngrid,'lapserate','lapse rate in the storm', &
           &                        'k/m',3,lapserate)
        call WRITEDIAGFI(ngrid,'deltahr','extra heating rates', &
           &                        'K/s',3,deltahr)
        call writediagfi(ngrid,'scheme','which scheme',&
                                                   ' ',2,real(scheme))

        END SUBROUTINE rocketduststorm

!=======================================================================  
! **********************************************************************
! Subroutine to determine the vertical transport with 
! a Van Leer advection scheme (copied from the sedimentation scheme --> see vlz_fi.F)
!***********************************************************************
      SUBROUTINE van_leer(nlay,q,pente_max,masse,w,wq)

      IMPLICIT NONE

!--------------------------------------------------------
! Input/Output Variables
!--------------------------------------------------------
      INTEGER,INTENT(IN) :: nlay       ! number of atmospheric layers
      REAL,INTENT(IN) ::  masse(nlay)  ! mass of the layer Dp/g
      REAL,INTENT(IN) :: pente_max     != 2 conseillee
      REAL,INTENT(INOUT) :: q(nlay)    ! mixing ratio (kg/kg)
      REAL,INTENT(INOUT) :: w(nlay)    ! atmospheric mass "transferred" at each timestep (kg.m-2)
      REAL,INTENT(INOUT) :: wq(nlay+1)

!--------------------------------------------------------
! Local Variables
!--------------------------------------------------------

      INTEGER i,l,j,ii
      REAL dzq(nlay),dzqw(nlay),adzqw(nlay),dzqmax
      REAL newmasse
      REAL sigw, Mtot, MQtot
      INTEGER m

! **********************************************************************
!  Mixing ratio vertical gradient at the levels
! **********************************************************************
      do l=2,nlay
            dzqw(l)=q(l-1)-q(l)
            adzqw(l)=abs(dzqw(l))
      enddo

      do l=2,nlay-1
            if(dzqw(l)*dzqw(l+1).gt.0.) then
                dzq(l)=0.5*(dzqw(l)+dzqw(l+1))
            else
                dzq(l)=0.
            endif
            dzqmax=pente_max*min(adzqw(l),adzqw(l+1))
            dzq(l)=sign(min(abs(dzq(l)),dzqmax),dzq(l))
      enddo

      dzq(1)=0.
      dzq(nlay)=0.

! **********************************************************************
!  Vertical advection
! ********************************************************************** 

       !! No flux at the model top:
       wq(nlay+1)=0.

       !! Surface flux up:
       if(w(1).lt.0.) wq(1)=0. ! warning : not always valid

       do l = 1,nlay-1

!      1) Compute wq where w < 0 (up) (UPWARD TRANSPORT)
!      ===============================

         if(w(l+1).le.0)then
!         Regular scheme (transfered mass < 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(-w(l+1).le.masse(l))then
            sigw=w(l+1)/masse(l)
            wq(l+1)=w(l+1)*(q(l)-0.5*(1.+sigw)*dzq(l))
!!-------------------------------------------------------
!          The following part should not be needed in the
!          case of an integration over subtimesteps
!!-------------------------------------------------------
!         Extended scheme (transfered mass > 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else 
             m = l-1
             Mtot = masse(m+1)
             MQtot = masse(m+1)*q(m+1)
             if (m.le.0)goto 77
             do while(-w(l+1).gt.(Mtot+masse(m)))
!             do while(-w(l+1).gt.Mtot)
                m=m-1
                Mtot = Mtot + masse(m+1)
                MQtot = MQtot + masse(m+1)*q(m+1)
                if (m.le.0)goto 77
             end do
 77          continue

             if (m.gt.0) then
                sigw=(w(l+1)+Mtot)/masse(m)
                wq(l+1)= -(MQtot + (-w(l+1)-Mtot)*         &
                (q(m)-0.5*(1.+sigw)*dzq(m))  )
             else
                w(l+1) = -Mtot
                wq(l+1) = -MQtot 
             end if
          endif ! (-w(l+1).le.masse(l))
     
!      2) Compute wq where w > 0 (down) (DOWNWARD TRANSPORT)     
!      ===============================

         else if(w(l).gt.0.)then

!         Regular scheme (transfered mass < 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(l).le.masse(l))then
            sigw=w(l)/masse(l)
            wq(l)=w(l)*(q(l)+0.5*(1.-sigw)*dzq(l))            
!!-------------------------------------------------------
!          The following part should not be needed in the
!          case of an integration over subtimesteps
!!-------------------------------------------------------
!         Extended scheme (transfered mass > 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else 
            m=l
            Mtot = masse(m)
            MQtot = masse(m)*q(m)
            if(m.ge.nlay)goto 88
            do while(w(l).gt.(Mtot+masse(m+1)))
                m=m+1
                Mtot = Mtot + masse(m)
                MQtot = MQtot + masse(m)*q(m)
                if(m.ge.nlay)goto 88
            end do
 88         continue
            if (m.lt.nlay) then
                sigw=(w(l)-Mtot)/masse(m+1)
                wq(l)=(MQtot + (w(l)-Mtot)* &
                          (q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
            else
                w(l) = Mtot
                wq(l) = MQtot 
            end if
          end if

         end if ! w<0 (up)

       enddo ! l = 1,nlay-1
       
       do l = 1,nlay

!         it cannot entrain more than available mass !
          if ( (wq(l+1)-wq(l)) .lt. -(masse(l)*q(l)) ) then
            wq(l+1) = wq(l)-masse(l)*q(l)
          end if

          q(l)=q(l) +  (wq(l+1)-wq(l))/masse(l)

       enddo
       
      END SUBROUTINE van_leer

!=======================================================================
! Initialization of the module variables

       subroutine ini_rocketduststorm_mod(ngrid)
       
       implicit none
       
       integer, intent(in) :: ngrid
       
       allocate(dustliftday(ngrid))
       
       end subroutine ini_rocketduststorm_mod
       
       subroutine end_rocketduststorm_mod
       
       implicit none
       
       if (allocated(dustliftday)) deallocate(dustliftday)

       end subroutine end_rocketduststorm_mod       
      
      END MODULE rocketduststorm_mod
