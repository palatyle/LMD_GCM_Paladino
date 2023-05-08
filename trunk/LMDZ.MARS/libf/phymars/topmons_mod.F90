      MODULE topmons_mod

      IMPLICIT NONE

!     sub-grid scale mountain mesh fraction
      REAL, SAVE, ALLOCATABLE :: alpha_hmons(:)

      CONTAINS

!=======================================================================
! ENTRAINMENT OF DUST ABOVE THE SUB-GRID SCALE TOPOGRAPHY
!=======================================================================
! calculation of the vertical wind velocity of topdust tracers
! transport scheme of topdust tracers
! detrainement of topdust into background dust
! -----------------------------------------------------------------------
! Authors: M. Vals; C. Wang; T. Bertrand; F. Forget 
! Institution: Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------

      SUBROUTINE topmons(ngrid,nlayer,nq,ptime,ptimestep,              &
                                 pq,pdq,pt,pdt,pplev,pplay,pzlev,      &
                                 pzlay,pdtsw,pdtlw,                    &
!             input for radiative transfer
                                 icount,zday,zls,tsurf,igout,aerosol,  &
                                 tauscaling,                           &
!             input sub-grid scale rocket dust storm
                                 totstormfract,clearatm,               &
!             input sub-grid scale cloud
                                 clearsky,totcloudfrac,                &
!             input sub-grid scale mountain
                                 nohmons,hsummit,                      &
!             output
                                 pdqtop,wfin,dsodust,dsords,dsotop,    &
                                 tauref)

      USE tracer_mod, only: igcm_topdust_mass,igcm_topdust_number &
                            ,igcm_dust_mass,igcm_dust_number      &
                            ,rho_dust 
      USE comcstfi_h, only: r,g,cpp,rcp
      USE dimradmars_mod, only: albedo,naerkind
      USE comsaison_h, only: dist_sol,mu0,fract
      USE surfdat_h, only: emis,co2ice,hmons,summit
      USE callradite_mod, only: callradite

      IMPLICIT NONE

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------

      INTEGER, INTENT(IN) :: ngrid ! number of horizontal grid points
      INTEGER, INTENT(IN) :: nlayer ! number of vertical grid points
      INTEGER, INTENT(IN) :: nq ! number of tracer species
      REAL, INTENT(IN) :: ptime
      REAL, INTENT(IN) :: ptimestep

      REAL, INTENT(IN) :: pq(ngrid,nlayer,nq) ! advected field nq
      REAL, INTENT(IN) :: pdq(ngrid,nlayer,nq)! tendancy field pq
      REAL, INTENT(IN) :: pt(ngrid,nlayer)    ! temperature at mid-layer (K)
      REAL, INTENT(IN) :: pdt(ngrid,nlayer)   ! tendancy temperature at mid-layer (K)

      REAL, INTENT(IN) :: pplay(ngrid,nlayer)     ! pressure at middle of the layers
      REAL, INTENT(IN) :: pplev(ngrid,nlayer+1)   ! pressure at intermediate levels
      REAL, INTENT(IN) :: pzlay(ngrid,nlayer)     ! altitude at the middle of the layers
      REAL, INTENT(IN) :: pzlev(ngrid,nlayer+1)   ! altitude at layer boundaries

      REAL, INTENT(IN) :: pdtsw(ngrid,nlayer)     ! (K/s) env
      REAL, INTENT(IN) :: pdtlw(ngrid,nlayer)     ! (K/s) env

!     input for second radiative transfer 
      INTEGER, INTENT(INOUT) :: icount
      REAL, INTENT(IN) :: zday
      REAL, INTENT(IN) :: zls
      REAL, INTENT(IN) :: tsurf(ngrid)
      INTEGER, INTENT(IN) :: igout
      REAL, INTENT(INOUT) :: tauscaling(ngrid)
!     input sub-grid scale rocket dust storm
      LOGICAL, INTENT(IN) :: clearatm
      REAL, INTENT(IN) :: totstormfract(ngrid)
!     input sub-grid scale water ice clouds
      LOGICAL, INTENT(IN) :: clearsky
      REAL, INTENT(IN) :: totcloudfrac(ngrid)
!     input sub-grid scale mountain
      LOGICAL, INTENT(IN) :: nohmons
      REAL, INTENT(IN) :: hsummit(ngrid)

!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------

      REAL, INTENT(OUT) :: pdqtop(ngrid,nlayer,nq)  ! tendancy field for dust
      REAL, INTENT(OUT) :: wfin(ngrid,nlayer+1)     ! final vertical wind velocity: combination of both wup and wrad
                                                    ! wfin < 0 means upward wind (in pressure levels/s)

!     output for second radiative transfer 
      REAL, INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind)
      REAL, INTENT(INOUT) :: dsodust(ngrid,nlayer)
      REAL, INTENT(INOUT) :: dsords(ngrid,nlayer)
      REAL, INTENT(INOUT) :: dsotop(ngrid,nlayer)
      REAL, INTENT(OUT) :: tauref(ngrid)

!--------------------------------------------------------
! Local variables 
!--------------------------------------------------------
      LOGICAL,SAVE :: firstcall=.true.
      INTEGER l,ig,tsub,iq,ll
      REAL zq0(ngrid,nlayer,nq)     ! initial tracers
      REAL zq(ngrid,nlayer,nq)      ! updated tracers
      REAL,PARAMETER:: mu0lim=0.001 ! solar zenith angle limit to determine the daytime

!     Variables for radiative transfer
      REAL  fluxsurf_lw1(ngrid)
      REAL  fluxsurf_sw1(ngrid,2)
      REAL  fluxtop_lw1(ngrid)
      REAL  fluxtop_sw1(ngrid,2)
      REAL  tau(ngrid,naerkind)
      REAL  taucloudtes(ngrid)
      REAL  rdust(ngrid,nlayer)
      REAL  rstormdust(ngrid,nlayer)
      REAL  rtopdust(ngrid,nlayer)
      REAL  rice(ngrid,nlayer)
      REAL  nuice(ngrid,nlayer)

!     Temperature profile
      REAL zt(ngrid,nlayer)                                  ! actual temperature at mid-layer (K)
      REAL ztlev(ngrid,nlayer)                               ! temperature at intermediate levels l+1/2 without last level
      REAL zdtlw1(ngrid,nlayer)                              ! (K/s) topdust
      REAL zdtsw1(ngrid,nlayer)                              ! (K/s) topdust
      REAL zdtlw1_lev(ngrid,nlayer),zdtsw1_lev(ngrid,nlayer) ! rad. heating rate at intermediate levels l+1/2
      REAL zdtlw_lev(ngrid,nlayer),zdtsw_lev(ngrid,nlayer)   ! rad. heating rate at intermediate levels l+1/2
      REAL zdtvert(ngrid,nlayer)                             ! dT/dz , lapse rate
      REAL deltahr(ngrid,nlayer+1)                           ! extra heating between the concentrated and non-concentrated dust

      REAL newzt(ngrid,nlayer)    ! recalculated temperature with extra radiative heating
      REAL t_top(ngrid,nlayer)    ! estimated temperature above the mountain
      REAL t_top_lev(ngrid,nlayer)! estimated temperature above the mountain at the levels
      REAL dt_top(ngrid)          ! temperature difference between the summit and the environment
      INTEGER lmons(ngrid)        ! layer of the mountain's slope
      INTEGER lsummit(ngrid)      ! layer of the mountain's summit
      REAL t_env(ngrid,nlayer)    ! recalculated temperature of the environment air "next" to the mountain

!     Vertical velocity profile
      REAL wrad(ngrid,nlayer+1)  ! vertical wind velocity resulting from radiative heating
      REAL wup(ngrid,nlayer+1)   ! vertical wind velocity calculated above the sub-grid scale mountain
      REAL wup2(ngrid,nlayer+1)  ! square of the vertical velocity calculated above the sub-grid scale mountain
      REAL w0(ngrid)             ! vertical wind velocity convergence at the top of the sub-grid scale mountain 
      INTEGER lwmax(ngrid)       ! level of the maximum vertical velocity

      REAL,PARAMETER:: wmax=10.  ! maximum vertical velocity resulting from radiative heating (m/s)
      REAL,PARAMETER:: secu=3.   ! coefficient on wfin to avoid dust crossing many layers during subtimestep
      REAL,PARAMETER:: k0=0.25   !max 10m/s: k0=1.      !max 3m/s: k0=0.25
      REAL,PARAMETER:: k1=5.e-4  !max 10m/s: k1=5.e-4   !max 3m/s: k1=5.e-4
      REAL,PARAMETER:: k2=5.e-3  !max 10m/s: k2=30.e-3  !max 3m/s: k2=5.e-3

!     Entrainment
      REAL entr(ngrid)                                             ! entrainment flux
      REAL rho(ngrid,nlayer),rhobarz(ngrid,nlayer+1)                 ! density, density at the levels
      REAL masse(ngrid,nlayer)                                     ! air mass
      REAL masse_pbl(ngrid)                                        ! total air mass within the PBL
      REAL masse_pbl_dust_mass(ngrid),masse_pbl_dust_number(ngrid) ! total air mass + dust mass/number within the PBL (kg/m2)
      REAL dqm_mass_pbl(ngrid),dqm_number_pbl(ngrid)               ! total mass of the entrained dust mass/number from the PBL (kg/m2)

!     Van Leer
      REAL zq_dust_mass(ngrid,nlayer),zq_dust_number(ngrid,nlayer)       ! mixing ratio background dust (kg/kg)
      REAL zq_topdust_mass(ngrid,nlayer),zq_topdust_number(ngrid,nlayer) ! mixing ratio topdust (kg/kg)
      REAL mr_dust_mass(ngrid,nlayer),mr_dust_number(ngrid,nlayer)       ! mixing ratio background dust (kg/kg)
      REAL mr_topdust_mass(ngrid,nlayer),mr_topdust_number(ngrid,nlayer) ! mixing ratio background dust + concentrated topdust (kg/kg)

      REAL w(ngrid,nlayer)          ! vertical flux above the sub-grid scale mountain
      REAL wqmass(ngrid,nlayer+1)   ! flux mass
      REAL wqnumber(ngrid,nlayer+1) ! flux number
      REAL qbar_mass_pbl(ngrid)     ! total dust mass mixing ratio within the PBL for Van Leer
      REAL qbar_number_pbl(ngrid)   ! total dust number mixing ratio within the PBL for Van Leer
      REAL masse_col(nlayer)        ! masse of the atmosphere in one column

      REAL dqvl_dust_mass(ngrid,nlayer)       ! tendancy pdq dust mass after vertical transport
      REAL dqvl_dust_number(ngrid,nlayer)     ! tendancy pdq dust number after vertical transport
      REAL dqvl_topdust_mass(ngrid,nlayer)    ! tendancy pdq topdust mass after vertical transport
      REAL dqvl_topdust_number(ngrid,nlayer)  ! tendancy pdq topdust number after vertical transport

      INTEGER nsubtimestep(ngrid)    ! number of subtimestep when calling Van Leer scheme 
      REAL subtimestep(ngrid)        ! ptimestep/nsubtimestep
      REAL dtmax                     ! considered minimum time interval needed for the dust to cross one vertical layer

!     Detrainment      
      REAL coefdetrain(ngrid,nlayer)          ! coefficient for detrainment : % of stormdust detrained
      REAL dqdet_topdust_mass(ngrid,nlayer)   ! tendancy pdq topdust mass after detrainment only
      REAL dqdet_topdust_number(ngrid,nlayer) ! tendancy pdq topdust number after detrainment only

!     Diagnostics
      REAL    zlaywmax(ngrid)
      REAL    detrain_rate(ngrid,nlayer)

      ! **********************************************************************
      ! **********************************************************************
      ! Parametrization of the entrainment by slope wind above the sub-grid
      ! scale topography
      ! **********************************************************************
      ! **********************************************************************      
      !!    1. Radiative transfer in topdust
      !!    2. Compute vertical velocity for topdust
      !!    3. Vertical transport
      !!    4. Detrainment
      ! **********************************************************************
      ! **********************************************************************
     
      ! **********************************************************************
      ! 0. Initializations
      ! **********************************************************************
      aerosol(:,:,:)=0.
      pdqtop(:,:,:) = 0.
      dqvl_topdust_mass(:,:)=0.
      dqvl_topdust_number(:,:)=0.
      dqvl_dust_mass(:,:)=0.
      dqvl_dust_number(:,:)=0.
      dqdet_topdust_mass(:,:)=0.
      dqdet_topdust_number(:,:)=0.
      wup(:,:)=0.
      wup2(:,:)=0.
      wrad(:,:)=0.
      w0(:)=0.
      wfin(:,:)=0.      
      w(:,:)=0.
      zq(:,:,:) = 0.
      zq0(:,:,:) = 0.
      entr(:) = 0.
      mr_dust_mass(:,:) = 0.
      mr_dust_number(:,:) = 0.
      mr_topdust_mass(:,:) = 0.
      mr_topdust_number(:,:) = 0.
      t_top(:,:) = 0.
      dt_top(:) = 0.
      newzt(:,:) = 0.
      lmons(:) = 1
      lsummit(:) =1
      lwmax(:) = 1
      deltahr(:,:) = 0.
      zdtvert(:,:) = 0.
      wqmass(:,:) = 0.
      wqnumber(:,:) = 0.       
      coefdetrain(:,:) = 1.
      dqm_mass_pbl(:)=0.
      dqm_number_pbl(:)=0.
      nsubtimestep(:)=1
      masse_pbl(:)=0.
      detrain_rate(:,:) = 0.

      !! Update the initial temperature
      zt(1:ngrid,1:nlayer)=pt(1:ngrid,1:nlayer)+pdt(1:ngrid,1:nlayer)*ptimestep
      newzt(1:ngrid,1:nlayer)=zt(1:ngrid,1:nlayer)
      t_env(1:ngrid,1:nlayer)=zt(1:ngrid,1:nlayer)
      t_top(1:ngrid,1:nlayer)=zt(1:ngrid,1:nlayer)

      !! Update the initial mixing ratios
      zq0(1:ngrid,1:nlayer,1:nq)=pq(1:ngrid,1:nlayer,1:nq)+pdq(1:ngrid,1:nlayer,1:nq)*ptimestep ! update of the background dust after rocket dust storm scheme
      zq(1:ngrid,1:nlayer,1:nq)=zq0(1:ngrid,1:nlayer,1:nq)
      zq_dust_mass(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_dust_mass)
      zq_dust_number(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_dust_number)
      zq_topdust_mass(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_topdust_mass)
      zq_topdust_number(1:ngrid,1:nlayer)=zq(1:ngrid,1:nlayer,igcm_topdust_number)
         
! ----------------------------------------------------------------------
! 1. Radiative heating
! ----------------------------------------------------------------------

       ! *********************************************************************
       ! 1.1. Call the second radiative transfer for topdust, obtain the extra heating
       ! *********************************************************************
        CALL callradite(icount,ngrid,nlayer,nq,zday,zls,zq,albedo,     &
                 emis,mu0,pplev,pplay,pt,tsurf,fract,dist_sol,igout,   &
                 zdtlw1,zdtsw1,fluxsurf_lw1,fluxsurf_sw1,fluxtop_lw1,  &
                 fluxtop_sw1,tauref,tau,aerosol,dsodust,tauscaling,    &
                 taucloudtes,rdust,rice,nuice,co2ice,rstormdust,rtopdust, &
                 totstormfract,clearatm,dsords,dsotop,alpha_hmons,nohmons,&
                 clearsky,totcloudfrac)
       ! **********************************************************************
       ! 1.2. Compute radiative vertical velocity for entrained topdust
       ! **********************************************************************
        DO ig=1,ngrid
          IF ( (mu0(ig) .gt. mu0lim) .and. (alpha_hmons(ig) .gt. 0.) ) THEN
            !! **********************************************************************
            !! Temperature profile above the mountain and in the close environment
            call t_topmons(nlayer,summit(ig),hsummit(ig),hmons(ig),zt(ig,:),pzlay(ig,:), &
                               t_top(ig,:),dt_top(ig),t_env(ig,:),lsummit(ig),lmons(ig))
            !! **********************************************************************
            !! Calculation of the extra heating : computing heating rates
            !! gradient at boundaries of each layer
            zdtlw1_lev(ig,1)=0.
            zdtsw1_lev(ig,1)=0.
            zdtlw_lev(ig,1)=0.
            zdtsw_lev(ig,1)=0.
            ztlev(ig,1)=zt(ig,1)
            t_top_lev(ig,1)=tsurf(ig)
            
            DO l=1,nlayer-1
              !! Calculation for the background dust AND the concentrated topdust
              zdtlw1_lev(ig,l+1)=(zdtlw1(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+ &
                           zdtlw1(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /  &
                              (pzlay(ig,l+1)-pzlay(ig,l))
           
              zdtsw1_lev(ig,l+1)=(zdtsw1(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+ &
                           zdtsw1(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /  &
                              (pzlay(ig,l+1)-pzlay(ig,l))
              !! Calculation for the background dust only
              zdtlw_lev(ig,l+1)=(pdtlw(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+   &
                           pdtlw(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /   &
                              (pzlay(ig,l+1)-pzlay(ig,l))
           
              zdtsw_lev(ig,l+1)=(pdtsw(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+   &
                           pdtsw(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /   &
                              (pzlay(ig,l+1)-pzlay(ig,l))
              !! Temperature at the levels
              ztlev(ig,l+1)=(zt(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+          &
                           zt(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /      &
                              (pzlay(ig,l+1)-pzlay(ig,l))
              !! Temperature above the mountain at the levels
              t_top_lev(ig,l+1)=(t_top(ig,l)*(pzlay(ig,l+1)-pzlev(ig,l+1))+          &
                           t_top(ig,l+1)*(pzlev(ig,l+1)-pzlay(ig,l)))  /      &
                              (pzlay(ig,l+1)-pzlay(ig,l))
            ENDDO ! DO l=1,nlayer-1
            !! **********************************************************************
            !! Vertical velocity profile for extra heating balanced by adiabatic cooling
            !! Gradient at boundaries of each layer, start from surface
            zdtvert(ig,1)=0.
            DO l=1,nlayer-1
              zdtvert(ig,l+1)=(t_top_lev(ig,l+1)-t_top_lev(ig,l))/(pzlay(ig,l+1)-pzlay(ig,l))
            ENDDO
            !! Extra heating balanced by adiabatic cooling
            DO l=1,nlayer
              deltahr(ig,l)=(zdtlw1_lev(ig,l)+zdtsw1_lev(ig,l))  &
                                          -(zdtlw_lev(ig,l)+zdtsw_lev(ig,l))
              wrad(ig,l)=-deltahr(ig,l)/(g/cpp+   &
                                         max(zdtvert(ig,l),-0.99*g/cpp))
              !! Limit of the vertical wind velocity in case of lapse rate close to adiabatic
              wrad(ig,l)=max(wrad(ig,l),-wmax) 
              wrad(ig,l)=min(wrad(ig,l),wmax)
              !! ATTENTION !! to simplify here, no downward velocity calculated
              if (wrad(ig,l).gt.0.) then
                wrad(ig,l)=0.  
              endif
            ENDDO
          ENDIF ! IF ((mu0(ig) .gt. mu0lim) .and. alpha_hmons(ig) .gt. 0.)
        ENDDO ! DO ig=1,ngrid

! ----------------------------------------------------------------------
! 2. Vertical velocity profile
! ----------------------------------------------------------------------

       ! **********************************************************************
       ! 2.1 Compute total vertical velocity for entrained topdust: buoyancy + radiative
       ! **********************************************************************
       DO ig=1,ngrid
         IF ( (mu0(ig) .gt. mu0lim) .and. (alpha_hmons(ig) .gt. 0.) ) THEN
           !! Positive buoyancy: negative vertical velocity entrains UP
           IF (dt_top(ig) .gt. 0.) THEN
             !! Convergence of vertical winds at the top of the mountain
             w0(ig)=-k0*g*sqrt((dt_top(ig)/t_env(ig,lsummit(ig))))
             !! Case the vertical velocity at the top of the mountain > wrad
             IF (w0(ig).lt.wrad(ig,lsummit(ig)+1)) then
               wup(ig,lsummit(ig)+1)=w0(ig)
               wfin(ig,lsummit(ig)+1)=w0(ig)
               newzt(ig,lsummit(ig))= t_top(ig,lsummit(ig))
               !! Temperature profile from the top of the summit to the top of the atmosphere
               DO l=lsummit(ig)+1,nlayer-1
                 !ATTENTION: test without wrad, T is only decreasing because of the adiabatic gradient
                 !newzt(ig,l)=newzt(ig,l-1)-(g/cpp)*        &
                 !               (pzlay(ig,l)-pzlay(ig,l-1))
                 !! Case of superadiabatic layer
                 IF (zdtvert(ig,l) .lt. -g/cpp) then
                    newzt(ig,l)=t_top(ig,l)
                 !! New temperature due to radiative heating balanced by adiabatic cooling
                 ELSE
                    newzt(ig,l)=newzt(ig,l-1)+ &
                               ( (deltahr(ig,l)/(-wfin(ig,l)))-g/cpp )*  &
                               ( pzlay(ig,l)-pzlay(ig,l-1) )
                 ENDIF ! (zdtvert(ig,l) .lt. -g/cpp)
                 !! Vertical velocity profile due to the presence of the mountain with the new temperature
                 !! Square of the vertical velocity profile with the new temperature
                 wup2(ig,l+1)=(1.-2*k1*(pzlev(ig,l+1)-pzlev(ig,l)))*&
                              (wup(ig,l)**2)+2.*k2*g*(pzlev(ig,l+1) & 
                             -pzlev(ig,l))*((newzt(ig,l)       &
                             -t_env(ig,l))/t_env(ig,l))
                 !! Vertical velocity profile if W2 > 0
                 IF (wup2(ig,l+1) .gt. 0.) THEN
                   wup(ig,l+1)=-sqrt(wup2(ig,l+1))
                   wfin(ig,l+1)=wup(ig,l+1)
                   IF (wup(ig,l+1) .gt. wrad(ig,l+1)) then
                     DO ll=l,nlayer-1
                       newzt(ig,ll)=zt(ig,ll)
                       wfin(ig,ll+1)=wrad(ig,ll+1)
                     ENDDO 
                     EXIT
                   ENDIF
                 !! Vertical velocity profile if W2 < 0
                 ELSE
                   DO ll=l,nlayer-1
                     newzt(ig,ll)=zt(ig,ll)
                     wfin(ig,ll+1)=wrad(ig,ll+1)
                   ENDDO
                     EXIT
                 ENDIF ! IF (wup2(ig,l+1) .gt. 0.)   
               ENDDO ! DO l=lsummit(ig)+1,nlayer-1
             !! Case the vertical velocity at the top of the mountain < wrad
             ELSE
               DO ll=lsummit(ig),nlayer-1
                  newzt(ig,ll)=zt(ig,ll)
                  wfin(ig,ll+1)=wrad(ig,ll+1)
               ENDDO
             ENDIF !(w0(ig).lt.wrad(ig,lsummit(ig)+1))
           !! Positive buoyancy: positive vertical velocity entrains DOWN
           ELSE
             DO l=lsummit(ig)+1,nlayer
                wfin(ig,l)=wrad(ig,l)
             ENDDO
           ENDIF ! IF (dt_top(ig) .gt. 0.)            
       ! **********************************************************************
       ! 2.2 Find the level lwmax of the maximum vertical velocity caused by the
       ! mountain presence (wup), where to entrain the dust from the PBL
       ! **********************************************************************
           IF (minloc(wup(ig,:),dim=1).lt.nlayer) THEN
             lwmax(ig)=minloc(wup(ig,:),dim=1)
           ENDIF
           zlaywmax(ig)=pzlay(ig,lwmax(ig)) ! wup(lwmax) acts on level lwmax, going to layer lwmax
       ! **********************************************************************
       ! 2.3 Compute the subtimestep to conserve the mass in the Van Leer transport
       ! **********************************************************************
           !! Calculation of the subtimesteps
           dtmax=ptimestep
           DO l=lwmax(ig),nlayer
             IF (wfin(ig,l).lt.0.) THEN
               dtmax=min(dtmax,(pzlev(ig,l)-pzlev(ig,l-1))/  &
                                 (secu*abs(wfin(ig,l))))
             ELSE IF (wfin(ig,l).gt.0.) then
               dtmax=min(dtmax,(pzlev(ig,l+1)-pzlev(ig,l))/  &
                                 (secu*abs(wfin(ig,l))))
             ENDIF
           ENDDO
           nsubtimestep(ig)= int(ptimestep/dtmax) 
           subtimestep(ig)=ptimestep/float(nsubtimestep(ig))
           !! Mass flux generated by wup in kg/m2
           DO l=1,nlayer!+1
             w(ig,l)=wfin(ig,l)*pplev(ig,l)/(r*ztlev(ig,l)) &
                      *subtimestep(ig)
           ENDDO ! l=1,nlayer
                        
          ENDIF ! ((mu0(ig) .gt. mu0lim)...
        ENDDO ! DO ig=1,ngrid 

! ----------------------------------------------------------------------
! 3. Dust entrainment by slope winds above the mountain
! ----------------------------------------------------------------------
        !! rho on the layer
        rho(:,:)=pplay(:,:)/(r*pt(:,:))
        !! rhobarz(l) = rho(l+1/2): rho on the levels
        rhobarz(:,1)=rho(:,1)
        DO l=2,nlayer
          rhobarz(:,l)=pplev(:,l)/(r*0.5*(pt(:,l)+pt(:,l-1)))
        ENDDO
        rhobarz(:,nlayer+1)=0. !top of the model
        !! Mass computation
        DO l=1,nlayer
          masse(:,l)=(pplev(:,l)-pplev(:,l+1))/g
        ENDDO

        DO ig=1,ngrid
          IF ( (mu0(ig) .gt. mu0lim) .and. (alpha_hmons(ig) .gt. 0.) ) THEN
            !! Total air mass within the PBL before entrainment (=> by PBL we mean between the surface and the layer where the vertical wind is maximum)
            masse_pbl(ig)=0.
            DO l=1,lwmax(ig)-1
              masse_pbl(ig)=masse_pbl(ig)+(pplev(ig,l)-pplev(ig,l+1))/g
            ENDDO
          ENDIF ! ((mu0(ig) .gt. mu0lim)...
        ENDDO ! ig=1,ngrid
       ! **********************************************************************
       ! 3.1. Entrainment
       ! **********************************************************************
        DO ig=1,ngrid
          IF ( (mu0(ig) .gt. mu0lim) .and. (alpha_hmons(ig) .gt. 0.) .and. (masse_pbl(ig) .gt. 0.) ) THEN
            !! Transport of background dust + concentrated topdust above lwmax
            DO l=lwmax(ig),nlayer
              mr_dust_mass(ig,l) = zq_dust_mass(ig,l)
              mr_dust_number(ig,l) = zq_dust_number(ig,l)
              mr_topdust_mass(ig,l) = zq_dust_mass(ig,l) &
                                     +zq_topdust_mass(ig,l)/alpha_hmons(ig)
              mr_topdust_number(ig,l) = zq_dust_number(ig,l) &
                                     +zq_topdust_number(ig,l)/alpha_hmons(ig)
            ENDDO ! l=lwmax(ig),nlayer
            !! Entrainment flux from the PBL
            entr(ig) = alpha_hmons(ig)*wup(ig,lwmax(ig))*rhobarz(ig,lwmax(ig))/masse_pbl(ig)
            !! If entrains more than available masse_pbl (abs(entr) x mass_pbl x ptimestep > masse_pbl)
            if (abs(entr(ig)) .gt. 1./ptimestep) then
              entr(ig) = -(1./ptimestep)
            end if
       ! **********************************************************************
       ! 3.2. Vertical transport: Van Leer
       ! **********************************************************************
            !! Loop over the subtimesteps
            DO tsub=1,nsubtimestep(ig)
              !! qbar_pbl: total mixing ratio within the PBL
              masse_pbl_dust_mass(ig)=0.
              masse_pbl_dust_number(ig)=0.
              DO l=1,lwmax(ig)-1
                masse_pbl_dust_mass(ig)=masse_pbl_dust_mass(ig)+zq_dust_mass(ig,l)*(pplev(ig,l)-pplev(ig,l+1))/g
                masse_pbl_dust_number(ig)=masse_pbl_dust_number(ig)+zq_dust_number(ig,l)*(pplev(ig,l)-pplev(ig,l+1))/g
              ENDDO
              qbar_mass_pbl(ig)=masse_pbl_dust_mass(ig)/masse_pbl(ig)
              qbar_number_pbl(ig)=masse_pbl_dust_mass(ig)/masse_pbl(ig)
              !! integration of the equation dq/dt = -(entr)q, with entr constant --> w0 < 0 when up so the equation is dq/dt = (entr)q
              dqm_mass_pbl(:)=0.
              dqm_number_pbl(:)=0.
              DO l=1,lwmax(ig)-1 ! this time we take the dust from the surface up to the level, where the velocity is maximum
                !! Total dust entrained from the PBL:
                dqm_mass_pbl(ig)=dqm_mass_pbl(ig)+(1.-exp(entr(ig)*subtimestep(ig)))*zq_dust_mass(ig,l)*(pplev(ig,l)-pplev(ig,l+1))/g
                dqm_number_pbl(ig)=dqm_number_pbl(ig)+(1.-exp(entr(ig)*subtimestep(ig)))*zq_dust_number(ig,l)*(pplev(ig,l)-pplev(ig,l+1))/g
                !! Background dust after entrainment (explicit: qdust(ig,l)=qdust(ig,l)-entr(ig)*qdust(ig,l)*ptimestep)
                zq_dust_mass(ig,l)=zq_dust_mass(ig,l)-(1.-exp(entr(ig)*subtimestep(ig)))*zq_dust_mass(ig,l)
                zq_dust_number(ig,l)=zq_dust_number(ig,l)-(1.-exp(entr(ig)*subtimestep(ig)))*zq_dust_number(ig,l)
              ENDDO
              !! Van Leer scheme (from lwmax to nlayer)
              wqmass(ig,:)=0.
              wqnumber(ig,:)=0.
              CALL van_leer(nlayer,mr_topdust_mass(ig,:),2.,lwmax(ig),  &
                  masse(ig,:),w(ig,:),wqmass(ig,:),masse_pbl(ig),dqm_mass_pbl(ig),alpha_hmons(ig),qbar_mass_pbl(ig))
              CALL van_leer(nlayer,mr_topdust_number(ig,:),2.,lwmax(ig),  &
                  masse(ig,:),w(ig,:),wqnumber(ig,:),masse_pbl(ig),dqm_number_pbl(ig),alpha_hmons(ig),qbar_number_pbl(ig))
            ENDDO !tsub=...
       ! **********************************************************************
       ! 3.3. Re-calculation of the mixing ratios after vertical transport
       ! **********************************************************************
            !! Redistribution from lwmax to nlayer
            DO l=lwmax(ig),nlayer
              !! General and "healthy" case
              IF (mr_topdust_mass(ig,l).ge.mr_dust_mass(ig,l)) THEN
                zq_dust_mass(ig,l) = mr_dust_mass(ig,l)
                zq_dust_number(ig,l) = mr_dust_number(ig,l)
                zq_topdust_mass(ig,l) = alpha_hmons(ig)*(mr_topdust_mass(ig,l)-mr_dust_mass(ig,l))
                zq_topdust_number(ig,l) = alpha_hmons(ig)*(mr_topdust_number(ig,l)-mr_dust_number(ig,l))
              !! Particular case: there is less than initial dust mixing ratio after the vertical transport
              ELSE
                zq_dust_mass(ig,l) = (1.-alpha_hmons(ig))*mr_dust_mass(ig,l)+alpha_hmons(ig)*mr_topdust_mass(ig,l)
                zq_dust_number(ig,l) = (1.-alpha_hmons(ig))*mr_dust_number(ig,l)+alpha_hmons(ig)*mr_topdust_number(ig,l)
                zq_topdust_mass(ig,l) = 0.
                zq_topdust_number(ig,l) = 0.
              ENDIF
            ENDDO ! DO l=lwmax(ig),nlayer
       ! **********************************************************************
       ! 3.4. Calculation of the tendencies of the vertical transport from the surface up to the top of the atmosphere
       ! **********************************************************************
            DO l=1,nlayer
              dqvl_topdust_mass(ig,l) = (zq_topdust_mass(ig,l)- &
                                    zq0(ig,l,igcm_topdust_mass)) /ptimestep
              dqvl_topdust_number(ig,l) = (zq_topdust_number(ig,l)- &
                                      zq0(ig,l,igcm_topdust_number)) /ptimestep
              dqvl_dust_mass(ig,l) = (zq_dust_mass(ig,l)-zq0(ig,l,igcm_dust_mass)) /ptimestep
              dqvl_dust_number(ig,l) = (zq_dust_number(ig,l)-zq0(ig,l,igcm_dust_number)) /ptimestep
           ENDDO

          ENDIF ! ((mu0(ig) .gt. mu0lim)...
        ENDDO ! ig=1,ngrid

! ----------------------------------------------------------------------
! 4. Detrainment: topdust is converted to background dust
! ----------------------------------------------------------------------

        !! **********************************************************************
        !! 4.1 Compute the coefficient of detrainment
        !! **********************************************************************
        DO ig=1,ngrid
          !DO l=lwmax(ig),nlayer-1
          DO l=1,nlayer!-1
            !! Detrainment during the day
            IF ( (mu0(ig) .gt. mu0lim) .and. (zq_topdust_mass(ig,l) .gt. zq_dust_mass(ig,l)*0.01)) THEN
               coefdetrain(ig,l)=1.*( rhobarz(ig,l+1)*abs(wfin(ig,l+1)) - rhobarz(ig,l)*abs(wfin(ig,l)) ) / masse(ig,l)
               !! Detrainment when abs(w(l)) > abs(w(l+1)), i.e. coefdetrain < 0
               if ( coefdetrain(ig,l).lt.0.) then
                dqdet_topdust_mass(ig,l)=-(1.-exp(coefdetrain(ig,l)*ptimestep))*zq_topdust_mass(ig,l)/ptimestep
                dqdet_topdust_number(ig,l)=-(1.-exp(coefdetrain(ig,l)*ptimestep))*zq_topdust_number(ig,l)/ptimestep
                ! for diagnostics
                detrain_rate(ig,l)=(1.-exp(coefdetrain(ig,l)*ptimestep))
                !! One cannot detrain more than available topdust
                if ( (abs(dqdet_topdust_mass(ig,l)).gt.zq_topdust_mass(ig,l)/ptimestep) .or. (abs(dqdet_topdust_number(ig,l)).gt.zq_topdust_number(ig,l)/ptimestep) ) then
                  dqdet_topdust_mass(ig,l)=-zq_topdust_mass(ig,l)/ptimestep
                  dqdet_topdust_number(ig,l)=-zq_topdust_number(ig,l)/ptimestep
                  detrain_rate(ig,l)=1.
                endif
               !else!entrainment
               !  dqdet_topdust_mass(ig,l)=-(1.-exp(coefdetrain(ig,l)*ptimestep))*zq_dust_mass(ig,l)/ptimestep
               !  dqdet_topdust_number(ig,l)=-(1.-exp(coefdetrain(ig,l)*ptimestep))*zq_dust_number(ig,l)/ptimestep
               endif
            !! Full detrainment during the night imposed
            ELSE
               dqdet_topdust_mass(ig,l)=-zq_topdust_mass(ig,l)/ptimestep
               dqdet_topdust_number(ig,l)=-zq_topdust_number(ig,l)/ptimestep
               detrain_rate(ig,l)=1.
            ENDIF
          ENDDO ! DO l=1,nlayer
        ENDDO ! DO ig=1,ngrid

! ----------------------------------------------------------------------
! 5. Final tendencies: vertical transport + detrainment
! ----------------------------------------------------------------------      

        DO ig=1,ngrid
          DO l=1,nlayer      
          pdqtop(ig,l,igcm_topdust_mass)=dqdet_topdust_mass(ig,l) &
                                                 +dqvl_topdust_mass(ig,l)
          pdqtop(ig,l,igcm_topdust_number)=dqdet_topdust_number(ig,l) &
                                                 +dqvl_topdust_number(ig,l)
          pdqtop(ig,l,igcm_dust_mass)= -dqdet_topdust_mass(ig,l) &
                                       +dqvl_dust_mass(ig,l) 
          pdqtop(ig,l,igcm_dust_number)= -dqdet_topdust_number(ig,l) &
                                       +dqvl_dust_number(ig,l) 
          ENDDO ! DO l=1,nlayer
        ENDDO ! DO ig=1,ngrid


! **********************************************************************
! WRITEDIAGFI
! **********************************************************************
        IF (ngrid.eq.1) THEN
               CALL WRITEDIAGFI(ngrid,'wup_top', &
                'wup_top','',1,wup(:,:))
               CALL WRITEDIAGFI(ngrid,'wrad_top', &
                'wrad_top','',1,wrad(:,:))
               CALL WRITEDIAGFI(ngrid,'wfin_top', &
                'wfin_top','',1,wfin(:,:))
               CALL WRITEDIAGFI(ngrid,'wlwmax_top', &
                'wlwmax_top','',0,wup(:,lwmax(1)))
               CALL WRITEDIAGFI(ngrid,'w0_top', &
                'w0_top','',0,wup(:,lsummit(1)+1))
               CALL WRITEDIAGFI(ngrid,'alpha_hmons', &
                'alpha_hmons','',0,alpha_hmons)
               CALL WRITEDIAGFI(ngrid,'zt_top', &
                'zt_top','',1,t_top(:,:))
               CALL WRITEDIAGFI(ngrid,'dt_top', &
                'dt_top','',0,dt_top(:))
               CALL WRITEDIAGFI(ngrid,'envtemp', &
                'envtemp','',1,zt(:,:))
               CALL WRITEDIAGFI(ngrid,'t_env', &
                't_env','',1,t_env(:,:))
               CALL WRITEDIAGFI(ngrid,'newzt_top', &
                'newzt_top','',1,newzt(:,:))
               CALL WRITEDIAGFI(ngrid,'deltahr_top', &
                'deltahr_top','',1,deltahr(:,:))
               CALL WRITEDIAGFI(ngrid,'rholev', &
                'rholev','',1,rho(:,:))
               CALL WRITEDIAGFI(ngrid,'rhobarz', &
                'rhobarz','',1,rhobarz(:,:))
               CALL WRITEDIAGFI(ngrid,'zlsummit', &
                'zlsummit','',0,pzlay(:,lsummit(1)))
               CALL WRITEDIAGFI(ngrid,'zlwmax', &
                'zlwmax','',0,pzlay(:,lwmax(1)))
               CALL WRITEDIAGFI(ngrid,'pzlev_top', &
                'pzlev_top','',1,pzlev(:,:))
               CALL WRITEDIAGFI(ngrid,'coefdetrain', &
                'coefdetrain','',1,coefdetrain(:,:))
               CALL WRITEDIAGFI(ngrid,'zdtvert', &
                'zdtvert','',1,zdtvert(:,:))
               CALL WRITEDIAGFI(ngrid,'entr', &
                'entr','',0,entr(:))
        ELSE
               CALL WRITEDIAGFI(ngrid,'wup_top', &
                'wup_top','',3,wup(:,:))
               CALL WRITEDIAGFI(ngrid,'wup2_top', &
                'wup2_top','',3,wup2(:,:))
               CALL WRITEDIAGFI(ngrid,'wrad_top', &
                'wrad_top','',3,wrad(:,:))
               CALL WRITEDIAGFI(ngrid,'wfin_top', &
                'wfin_top','',3,wfin(:,:))
               CALL WRITEDIAGFI(ngrid,'alpha_hmons', &
                'alpha_hmons','',2,alpha_hmons)
               CALL WRITEDIAGFI(ngrid,'zt_top', &
                'zt_top','',3,t_top(:,:))
               CALL WRITEDIAGFI(ngrid,'ztemp', &
                'envtemp','',3,zt(:,:))
               CALL WRITEDIAGFI(ngrid,'zt_env', &
                't_env','',3,t_env(:,:))
               CALL WRITEDIAGFI(ngrid,'zdtvert_top', &
                'zdtvert_top','',3,zdtvert(:,:))
               CALL WRITEDIAGFI(ngrid,'newzt_top', &
                'newzt_top','',3,newzt(:,:))
               CALL WRITEDIAGFI(ngrid,'deltahr_top', &
                'deltahr_top','',3,deltahr(:,:))
               CALL WRITEDIAGFI(ngrid,'rhobarz', &
                'rhobarz','',3,rhobarz(:,:))
               CALL WRITEDIAGFI(ngrid,'zlwmax', &
                'zlwmax','',2,zlaywmax(:))
               CALL WRITEDIAGFI(ngrid,'coefdetrain', &
                'coefdetrain','',3,coefdetrain(:,:))
               CALL WRITEDIAGFI(ngrid,'detrain_rate', &
                              'detrain_rate', &
                              '',3,detrain_rate(:,:))

        ENDIF
        
        END SUBROUTINE topmons

!=======================================================================     

! **********************************************************************
! Subroutine to determine both temperature profiles above and near 
! a sub-grid scale mountain
!***********************************************************************
       SUBROUTINE t_topmons(nlayer,summit,hsummit,hmons,zt,zlay,t_top,dt_top,t_env,lsummit,lmons)

       USE comcstfi_h, only: g,cpp

       IMPLICIT NONE

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------
       integer,intent(in) :: nlayer
       real,intent(in) :: summit          ! Height of the mountain summit
       real,intent(in) :: hmons           ! Height of the mountain slope
       real,intent(in) :: hsummit         ! Height of the mountain summit from the GCM surface
       real,intent(in) :: zt(nlayer)      ! large scale temperature profile
       real,intent(in) :: zlay(nlayer)    ! height at the middle of each layer
!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------
       real,intent(out) :: t_top(nlayer) ! temperature  above the mountain
       real,intent(out) :: dt_top        ! temperature difference between the slope and
                                         ! the environment
       real,intent(out) :: t_env(nlayer)   ! temperature of the environment next to the mountain
       integer,intent(out) :: lsummit      ! layer reached by the summit height in the GCM
       integer,intent(out) :: lmons        ! layer of the real temperature to be considered near to the mountain top

!--------------------------------------------------------
! Local Variables
!--------------------------------------------------------
       integer l,ll
       real zlmons

! **********************************************************************
       t_env(:)=zt(:)
       t_top(:)=zt(:)
       dt_top=0.
       lsummit=1
       lmons=1

       !! Layer where the dust is injected
       do while (hsummit .ge. zlay(lsummit))
           !! highest layer reached by the mountain (we don't interpolate and define zlay(lsummit) as the top of the mountain)
           lsummit=lsummit+1
       enddo
       !! temperature above the mountain is set to the surface temperature 
       t_top(lsummit)=zt(1)
       do l=lsummit,nlayer-1
          !! temperature above the mountain follows the adiabatic
          t_top(l+1)=t_top(l)-(zlay(l+1)-zlay(l))*g/cpp
       enddo

       !! Layer where to take the temperature above the mountain
       do while (hmons .ge. zlay(lmons))
           lmons=lmons+1
       enddo
       !! Real temperature of the environment near to the mountain is t_env(lsummit)=zt(lmons)
       !! (More simple and makes no real difference is to impose t_env(:)=zt(:))
       if (lmons.gt.lsummit) then
         t_env(lsummit)=zt(lmons)
         zlmons=zlay(lmons)
         do l=0,nlayer-lsummit-1!2
           zlmons=zlmons+(zlay(lsummit+l+1)-zlay(lsummit+l))
           ll=lmons
           do while (zlmons .ge. zlay(ll) .and. ll.lt.nlayer )
             ll=ll+1
           enddo
           ll=ll-1
           !! Interpolation above lsummit
           t_env(lsummit+l+1)= zt(ll) + ((zlmons-zlay(ll))*(zt(ll+1)-zt(ll))/(zlay(ll+1)-zlay(ll)))
         enddo
         t_env(nlayer)=zt(nlayer)
       endif ! (lmons.gt.lsummit)
       do l=1,nlayer
          !! temperature above the mountain follows the adiabatic up to reach the environment temperature 
          if (t_top(l).le.t_env(l)) then
            t_top(l)=t_env(l)
          endif
       enddo
       !!  Temperature delta created at the top of the mountain
       dt_top = t_top(lsummit)-t_env(lsummit)

       END SUBROUTINE t_topmons

!=======================================================================  
! **********************************************************************
! Subroutine to determine the vertical transport with 
! a Van Leer advection scheme (copied from the sedimentation scheme --> see vlz_fi.F)
!***********************************************************************
      SUBROUTINE van_leer(nlay,q,pente_max,lwmax,masse,w,wq,masse_pbl,dqm_pbl,alpha_hmons,qbar_pbl) 

      IMPLICIT NONE

!--------------------------------------------------------
! Input/Output Variables
!--------------------------------------------------------
      INTEGER,INTENT(IN) :: nlay       ! number of atmospheric layers
      REAL,INTENT(IN) :: masse(nlay)   ! mass of the layer Dp/g
      REAL,INTENT(IN) :: pente_max     != 2 conseillee
      INTEGER,INTENT(IN) :: lwmax      ! layer of dust injection above the mountain, where the vertical velocity is maximum
      REAL,INTENT(INOUT) :: w(nlay)    ! atmospheric mass "transferred" at each timestep (kg.m-2)
      REAL,INTENT(INOUT) :: wq(nlay+1) 
      REAL,INTENT(INOUT) :: q(nlay)    ! mixing ratio (kg/kg)
      REAL,INTENT(IN) :: masse_pbl
      REAL,INTENT(IN) :: dqm_pbl
      REAL,INTENT(IN) :: alpha_hmons
      REAL,INTENT(IN) :: qbar_pbl

!--------------------------------------------------------
! Local Variables
!--------------------------------------------------------

      INTEGER i,l,j,ii
      REAL dzq(nlay),dzqw(nlay),adzqw(nlay),dzqmax
      REAL sigw, Mtot, MQtot
      INTEGER m

! **********************************************************************
!  Mixing ratio vertical gradient at the levels
! **********************************************************************
      do l=lwmax+1,nlay !l=2,nlay
            dzqw(l)=q(l-1)-q(l)
            adzqw(l)=abs(dzqw(l))
      enddo

      do l=lwmax+1,nlay-1 !l=2,nlay-1
            if(dzqw(l)*dzqw(l+1).gt.0.) then
                dzq(l)=0.5*(dzqw(l)+dzqw(l+1))
            else
                dzq(l)=0.
            endif
            dzqmax=pente_max*min(adzqw(l),adzqw(l+1))
            dzq(l)=sign(min(abs(dzq(l)),dzqmax),dzq(l))
      enddo

      dzq(lwmax)=0. !! dzq(1)=0. => Attention: in the case of van leer only above the mountain it is very important to initialize the transport from lwmax !!
      dzq(nlay)=0.
! **********************************************************************
!  Vertical advection
! ********************************************************************** 

       !! No flux at the model top:
       wq(nlay+1)=0.

       !! Mass flux injected at lwmax
       wq(lwmax) = -dqm_pbl/alpha_hmons ! dqm_pbl = alpha_hmons x wup(wmax) x rho(lwmax) x ptimestep x qbar_pbl
                                        ! /alpha_hmons because the variable considered is mr_topdust
       do l = lwmax,nlay-1

!      1) Compute wq where w < 0 (up) (UPWARD TRANSPORT)     
!      ===============================

         if (w(l+1).le.0) then
!         Regular scheme (transfered mass < 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(-w(l+1).le.masse(l))then
            sigw=w(l+1)/masse(l)
            wq(l+1)=w(l+1)*(q(l)-0.5*(1.+sigw)*dzq(l))
!!-------------------------------------------------------
!          The following part should not be needed in the
!          case of an integration over subtimesteps
!!-------------------------------------------------------
!!         Extended scheme (transfered mass > 1 layer)
!!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!          else
!!             m = l-1
!!             Mtot = masse(m+1)
!!             MQtot = masse(m+1)*q(m+1)
!!             if (m.lt.lwmax)goto 77!if (m.le.0)goto 77
!!             do while(-w(l+1).gt.(Mtot+masse(m)))
!!                m=m-1
!!                Mtot = Mtot + masse(m+1)
!!                MQtot = MQtot + masse(m+1)*q(m+1)
!!                if (m.lt.lwmax)goto 77!if (m.le.0)goto 77
!!             end do
!! 77          continue
!!
!!              if (m.gt.lwmax) then!if (m.gt.0) then
!!                sigw=(w(l+1)+Mtot)/masse(m)
!!                wq(l+1)= -(MQtot + (-w(l+1)-Mtot)*         &
!!                (q(m)-0.5*(1.+sigw)*dzq(m)))
!!              else if ((-w(l+1).gt.(Mtot))) then
!!                Mtot=Mtot+masse_pbl!*alpha_hmons
!!                MQtot=MQtot+masse_pbl*qbar_pbl!*alpha_hmons
!!                sigw=(w(l+1)+Mtot)/masse(m)
!!                wq(l+1)= -(MQtot + (-w(l+1)-Mtot)*qbar_pbl)
!!              else
!!                w(l+1) = -Mtot
!!                wq(l+1) = -MQtot
!!              end if
!!
          endif

!      2) Compute wq where w > 0 (down) (DOWNWARD TRANSPORT)     
!      ===============================

         else if (w(l).gt.0.) then
!         Regular scheme (transfered mass < 1 layer)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(l).le.masse(l))then
            sigw=w(l)/masse(l)
            wq(l)=w(l)*(q(l)+0.5*(1.-sigw)*dzq(l))
!!-------------------------------------------------------
!          The following part should not be needed in the
!          case of an integration over subtimesteps
!!-------------------------------------------------------
!!         Extended scheme (transfered mass > 1 layer)
!!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!          else
!!            m=l
!!            Mtot = masse(m)
!!            MQtot = masse(m)*q(m)
!!            if(m.ge.nlay)goto 88
!!            do while(w(l).gt.(Mtot+masse(m+1)))
!!               m=m+1
!!                Mtot = Mtot + masse(m)
!!                MQtot = MQtot + masse(m)*q(m)
!!                if(m.ge.nlay)goto 88
!!            end do
!! 88         continue
!!            if (m.lt.nlay) then
!!                sigw=(w(l)-Mtot)/masse(m+1)
!!                wq(l)=(MQtot + (w(l)-Mtot)* &
!!                          (q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
!!            else
!!                w(l) = Mtot
!!                wq(l) = MQtot 
!!            end if
          end if ! (w(l).le.masse(l))

         end if ! (w(l+1).le.0)

       enddo ! l = lwmax+1,nlay-1

! **********************************************************************
!  Mixing ratio update after the vertical transport
! ********************************************************************** 
       do l = lwmax,nlay-1 !l = 1,nlay-1  ! loop different than when w>0

         if ( (wq(l+1)-wq(l)) .lt. -(masse(l)*q(l)) .and. (abs(wq(l+1)).gt. 0.) .and.  (abs(wq(l)).gt. 0.) ) then
           wq(l+1) = wq(l)-masse(l)*q(l)
         end if

         q(l)=q(l) +  (wq(l+1)-wq(l))/masse(l)

       enddo


      end subroutine van_leer

!********************************************************************************
       ! initialization module variables
       subroutine ini_topmons_mod(ngrid,nlayer)
       
       implicit none
       
       integer, intent(in) :: ngrid
       integer, intent(in) :: nlayer
       
       allocate(alpha_hmons(ngrid))

       end subroutine ini_topmons_mod
       
       subroutine end_topmons_mod
       
       implicit none
       
       if (allocated(alpha_hmons)) deallocate(alpha_hmons)

       end subroutine end_topmons_mod       
      
      END MODULE topmons_mod
