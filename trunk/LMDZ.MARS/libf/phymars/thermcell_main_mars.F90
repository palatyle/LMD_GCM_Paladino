!=======================================================================
! THERMCELL_MAIN_MARS
!=======================================================================
! This routine is called by calltherm_interface and is inside a sub-timestep
! loop. It computes thermals properties from parametrized entrainment and
! detrainment rate as well as the source profile.
! Mass flux are then computed and temperature and CO2 MMR are transported.
!=======================================================================
! Author : A. Colaitis 2011-01-05 (with updates 2011-2013)
!          after C. Rio and F. Hourdin
! Institution : Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------
! Corresponding author : A. Spiga aymeric.spiga_AT_upmc.fr
! -----------------------------------------------------------------------
! ASSOCIATED FILES
! --> calltherm_interface.F90
! --> thermcell_dqup.F90
! --> comtherm_h.F90
!=======================================================================
! Reference paper:
! A. Colaïtis, A. Spiga, F. Hourdin, C. Rio, F. Forget, and E. Millour. 
! A thermal plume model for the Martian convective boundary layer. 
! Journal of Geophysical Research (Planets), 118:1468-1487, July 2013.
! http://dx.doi.org/10.1002/jgre.20104
! http://arxiv.org/abs/1306.6215
! -----------------------------------------------------------------------
! Reference paper for terrestrial plume model:
! C. Rio and F. Hourdin.
! A thermal plume model for the convective boundary layer : Representation of cumulus clouds.
! Journal of the Atmospheric Sciences, 65:407-425, 2008.
! -----------------------------------------------------------------------

      SUBROUTINE thermcell_main_mars(ngrid,nlayer,nq &
     &                  ,tracer,igcm_co2 &
     &                  ,ptimestep  &
     &                  ,pplay,pplev,pphi,zlev,zlay  &
     &                  ,pu,pv,pt,pq,pq2  &
     &                  ,pdtadj,pdqadj  &
     &                  ,fm,entr,detr,lmax,zmax,limz  &
     &                  ,zw2,fraca &
     &                  ,zpopsk,heatFlux,heatFlux_down &
     &                  ,buoyancyOut, buoyancyEst)

      USE comtherm_h
#ifndef MESOSCALE
      use planetwide_mod, only: planetwide_maxval
#endif
      ! SHARED VARIABLES. This needs adaptations in another climate model.
      ! contains physical constant values such as
      ! "g" : gravitational acceleration (m.s-2)
      ! "r" : recuced gas constant (J.K-1.mol-1)
      USE comcstfi_h

      IMPLICIT NONE

!=======================================================================

! ============== INPUTS ==============

      INTEGER, INTENT(IN) :: ngrid ! number of horizontal grid points
      INTEGER, INTENT(IN) :: nlayer ! number of vertical grid points
      INTEGER, INTENT(IN) :: nq ! number of tracer species
      LOGICAL, INTENT(IN) :: tracer ! =.true. if tracers are present and to be transported
      INTEGER, INTENT(IN) :: igcm_co2 ! index of the CO2 tracer in mixing ratio array
                                      ! --> 0 if no tracer is CO2 (or no tracer at all)
                                      ! --> this prepares special treatment for polar night mixing
      REAL, INTENT(IN) :: ptimestep !subtimestep (s)
      REAL, INTENT(IN) :: pt(ngrid,nlayer) !temperature (K)
      REAL, INTENT(IN) :: pu(ngrid,nlayer) !u component of the wind (ms-1)
      REAL, INTENT(IN) :: pv(ngrid,nlayer) !v component of the wind (ms-1)
      REAL, INTENT(IN) :: pq(ngrid,nlayer,nq) !tracer concentration (kg/kg)
      REAL, INTENT(IN) :: pq2(ngrid,nlayer) ! Turbulent Kinetic Energy
      REAL, INTENT(IN) :: pplay(ngrid,nlayer) !Pressure at the middle of the layers (Pa)
      REAL, INTENT(IN) :: pplev(ngrid,nlayer+1) !intermediate pressure levels (Pa)
      REAL, INTENT(IN) :: pphi(ngrid,nlayer) !Geopotential at the middle of the layers (m2s-2)
      REAL, INTENT(IN) :: zlay(ngrid,nlayer) ! altitude at the middle of the layers
      REAL, INTENT(IN) :: zlev(ngrid,nlayer+1) ! altitude at layer boundaries

! ============== OUTPUTS ==============

! TEMPERATURE
      REAL, INTENT(OUT) :: pdtadj(ngrid,nlayer) !temperature change from thermals dT/dt (K/s)

! DIAGNOSTICS
      REAL, INTENT(OUT) :: zw2(ngrid,nlayer+1) ! vertical velocity (m/s)
      REAL, INTENT(OUT) :: heatFlux(ngrid,nlayer)   ! interface heatflux
      REAL, INTENT(OUT) :: heatFlux_down(ngrid,nlayer) ! interface heat flux from downdraft

      INTEGER, INTENT(OUT) :: limz ! limit vertical index for integration

! ============== LOCAL ================
      REAL :: pdqadj(ngrid,nlayer,nq) !tracer change from thermals dq/dt, only for CO2 (the rest can be advected outside of the loop)

! dummy variables when output not needed :

      REAL :: buoyancyOut(ngrid,nlayer)  ! interlayer buoyancy term
      REAL :: buoyancyEst(ngrid,nlayer)  ! interlayer estimated buoyancy term

! ============== LOCAL ================

      INTEGER ig,k,l,ll,iq
      INTEGER lmax(ngrid),lmin(ngrid),lalim(ngrid)
      REAL zmax(ngrid)
      REAL ztva(ngrid,nlayer),zw_est(ngrid,nlayer+1),ztva_est(ngrid,nlayer)
      REAL zh(ngrid,nlayer)
      REAL zdthladj(ngrid,nlayer)
      REAL zdthladj_down(ngrid,nlayer)
      REAL ztvd(ngrid,nlayer)
      REAL ztv(ngrid,nlayer)
      REAL zu(ngrid,nlayer),zv(ngrid,nlayer),zo(ngrid,nlayer)
      REAL zva(ngrid,nlayer)
      REAL zua(ngrid,nlayer)

      REAL zta(ngrid,nlayer)
      REAL fraca(ngrid,nlayer+1)
      REAL q2(ngrid,nlayer)
      REAL rho(ngrid,nlayer),rhobarz(ngrid,nlayer),masse(ngrid,nlayer)
      REAL zpopsk(ngrid,nlayer)

      REAL wmax(ngrid)
      REAL fm(ngrid,nlayer+1),entr(ngrid,nlayer),detr(ngrid,nlayer)

      REAL fm_down(ngrid,nlayer+1)

      REAL ztla(ngrid,nlayer)

      REAL f_star(ngrid,nlayer+1),entr_star(ngrid,nlayer)
      REAL detr_star(ngrid,nlayer)
      REAL alim_star_tot(ngrid)
      REAL alim_star(ngrid,nlayer)
      REAL alim_star_clos(ngrid,nlayer)
      REAL f(ngrid)

      REAL detrmod(ngrid,nlayer)

      REAL teta_th_int(ngrid,nlayer)
      REAL teta_env_int(ngrid,nlayer)
      REAL teta_down_int(ngrid,nlayer)

      CHARACTER (LEN=80) :: abort_message
      INTEGER ndt

! ============= PLUME VARIABLES ============

      REAL w_est(ngrid,nlayer+1)
      REAL wa_moy(ngrid,nlayer+1)
      REAL wmaxa(ngrid)
      REAL zdz,zbuoy(ngrid,nlayer),zw2m
      LOGICAL activecell(ngrid),activetmp(ngrid)
      INTEGER tic

! ==========================================

! ============= HEIGHT VARIABLES ===========

      REAL num(ngrid)
      REAL denom(ngrid)
      REAL zlevinter(ngrid)

! =========================================

! ============= CLOSURE VARIABLES =========

      REAL zdenom(ngrid)
      REAL alim_star2(ngrid)
      REAL alim_star_tot_clos(ngrid)
      INTEGER llmax

! =========================================

! ============= FLUX2 VARIABLES ===========

      INTEGER ncorecfm1,ncorecfm2,ncorecfm3,ncorecalpha
      INTEGER ncorecfm4,ncorecfm5,ncorecfm6,ncorecfm7,ncorecfm8
      REAL zfm
      REAL f_old,ddd0,eee0,ddd,eee,zzz
      REAL fomass_max,alphamax

! =========================================

! ============== Theta_M Variables ========

      REAL m_co2, m_noco2, A , B
      SAVE A, B
      REAL zhc(ngrid,nlayer)
      REAL ratiom(ngrid,nlayer)

! =========================================

!-----------------------------------------------------------------------
!   initialization:
!   ---------------

      entr(:,:)=0. ! entrainment mass flux
      detr(:,:)=0. ! detrainment mass flux
      fm(:,:)=0. ! upward mass flux
      zhc(:,:)=pt(:,:)/zpopsk(:,:) ! potential temperature
      ndt=1

!.......................................................................
!  Special treatment for co2:
!.......................................................................
! **********************************************************************
! In order to take into account the effect of vertical molar mass
! gradient on convection, we define modified theta that depends
! on the mass mixing ratio of Co2 in the cell.
! See for details:
!
! Forget, F. and Millour, E. et al. "Non condensable gas enrichment and depletion
! in the martian polar regions", third international workshop on the Mars Atmosphere:
! Modeling and Observations, 1447, 9106. year: 2008
!
! This is especially important for modelling polar convection.
! **********************************************************************
       if (igcm_co2.ne.0) then

         m_co2 = 44.01E-3  ! CO2 molecular mass (kg/mol)
         m_noco2 = 33.37E-3  ! Non condensible mol mass (kg/mol)
         ! Compute A and B coefficient use to compute
         ! mean molecular mass Mair defined by
         ! 1/Mair = q(igcm_co2)/m_co2 + (1-q(igcm_co2))/m_noco2
         ! 1/Mair = A*q(igcm_co2) + B
         A =(1/m_co2 - 1/m_noco2)
         B=1/m_noco2

!     Special case if one of the tracers is CO2 gas
         DO l=1,nlayer
           DO ig=1,ngrid
            ztv(ig,l) = zhc(ig,l)*(A*pq(ig,l,igcm_co2)+B)
           ENDDO
         ENDDO
       else
          ztv(:,:)=zhc(:,:)
       end if

!------------------------------------------------------------------------
! where are the different quantities defined ?
!------------------------------------------------------------------------
!                       --------------------
!
!
!                       + + + + + + + + + + +
!
!
!  wa, fraca, wd, fracd --------------------   zlev(2), rhobarz
!  wh,wt,wo ...
!
!                       + + + + + + + + + + +  zh,zu,zv,zo,rho
!
!
!                       --------------------   zlev(1)
!                       \\\\\\\\\\\\\\\\\\\\
!
!

!-----------------------------------------------------------------------
!   Densities at layer and layer interface (see above), mass:
!-----------------------------------------------------------------------

      rho(:,:)=pplay(:,:)/(r*pt(:,:))

      rhobarz(:,1)=rho(:,1)

      do l=2,nlayer
          rhobarz(:,l)=pplev(:,l)/(r*0.5*(pt(:,l)+pt(:,l-1)))
      enddo

! mass computation
      do l=1,nlayer
         masse(:,l)=(pplev(:,l)-pplev(:,l+1))/g
      enddo


!-----------------------------------------------------------------
!   Schematic representation of an updraft:
!------------------------------------------------------------------
!
!             /|\
!    --------  |  F_k+1 -------   
!                              ----> D_k
!             /|\              <---- E_k , A_k
!    --------  |  F_k --------- 
!                              ----> D_k-1
!                              <---- E_k-1 , A_k-1
!
!
!    ---------------------------
!
!    ----- F_lmax+1=0 ----------         \
!            lmax     (zmax)              |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------          |
!                                         |  E
!    ---------------------------          |  D
!                                         |
!    ---------------------------          |
!                                         |
!    ---------------------------  \       |
!            lalim                 |      |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!                                  | A    |
!    ---------------------------   |      |
!                                  |      |
!    ---------------------------   |      |
!    lmin  (=1 pour le moment)     |      |
!    ----- F_lmin=0 ------------  /      /
!
!    ---------------------------
!    //////////////////////////
!

!=============================================================================
! Mars version: no phase change is considered, we use a "dry" definition
! for the potential temperature.
!=============================================================================

!------------------------------------------------------------------
!  1. alim_star is the source layer vertical profile in the lowest layers
! of the thermal plume. Computed from the air buoyancy
!  2. lmin and lalim are the indices of begining and end of source profile
!------------------------------------------------------------------
!
      entr_star(:,:)=0. ; detr_star(:,:)=0. 
      alim_star(:,:)=0. ; alim_star_tot(:)=0.
      lmin(:)=1

!-----------------------------------------------------------------------------
!  3. wmax and zmax are maximum vertical velocity and altitude of a 
!     conservative plume (entrainment = detrainment = 0) using only
!     the source layer. This is a CAPE computation used for determining 
!     the closure mass flux. 
!-----------------------------------------------------------------------------

! ===========================================================================
! ===================== PLUME ===============================================
! ===========================================================================

! Initialization
      ztva(:,:)=ztv(:,:) ! temperature in the updraft = temperature of the env.
      ztva_est(:,:)=ztva(:,:) ! estimated temp. in the updraft
      ztla(:,:)=0. !intermediary variable
      zdz=0. !layer thickness
      zbuoy(:,:)=0. !buoyancy
      w_est(:,:)=0. !estimated vertical velocity
      f_star(:,:)=0. !non-dimensional upward mass flux f*
      wa_moy(:,:)=0. !vertical velocity

! Some more initializations
      wmaxa(:)=0.
      lalim(:)=1

!-------------------------------------------------------------------------
! We consider as an activecell columns where the two first layers are
! convectively unstable
! When it is the case, we compute the source layer profile (alim_star)
! see paper appendix 4.1 for details on the source layer
!-------------------------------------------------------------------------

      activecell(:)=ztv(:,1)>ztv(:,2)
          do ig=1,ngrid
            if (ztv(ig,1)>=(ztv(ig,2))) then
               alim_star(ig,1)=MAX((ztv(ig,1)-ztv(ig,2)),0.)  &
     &                       *sqrt(zlev(ig,2))
               lalim(ig)=2
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,1)
            endif
         enddo

       do l=2,nlayer-1
         do ig=1,ngrid
           if (ztv(ig,l)>(ztv(ig,l+1)) .and. ztv(ig,1)>=ztv(ig,l) &
     &             .and. (alim_star(ig,l-1).ne. 0.)) then
               alim_star(ig,l)=MAX((ztv(ig,l)-ztv(ig,l+1)),0.)  &
     &                       *sqrt(zlev(ig,l+1))
                lalim(ig)=l+1
               alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
           endif
        enddo
      enddo
      do l=1,nlayer
         do ig=1,ngrid
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo

      alim_star_tot(:)=1.

! We compute the initial squared velocity (zw2) and non-dimensional upward mass flux
! (f_star) in the first and second layer from the source profile.

      do ig=1,ngrid
          if (activecell(ig)) then
          ztla(ig,1)=ztv(ig,1)
          f_star(ig,1)=0.
          f_star(ig,2)=alim_star(ig,1)
          zw2(ig,2)=2.*g*(ztv(ig,1)-ztv(ig,2))/ztv(ig,2)  &
      &                     *(zlev(ig,2)-zlev(ig,1))  &
      &                    *0.4*pphi(ig,1)/(pphi(ig,2)-pphi(ig,1)) !0.4=von Karman constant
          w_est(ig,2)=zw2(ig,2)
          endif
      enddo

!==============================================================================
!==============================================================================
!==============================================================================
! LOOP ON VERTICAL LEVELS
!==============================================================================
      do l=2,nlayer-1
!==============================================================================
!==============================================================================
!==============================================================================


! is the thermal plume still active ?
        do ig=1,ngrid
             activecell(ig)=activecell(ig) &
      &                 .and. zw2(ig,l)>1.e-10 &
      &                 .and. f_star(ig,l)+alim_star(ig,l)>1.e-10
        enddo

!---------------------------------------------------------------------------
!
! .I. INITIALIZATION
!
! Computations of the temperature and buoyancy properties in layer l, 
! without accounting for entrainment and detrainment. We are therefore
! assuming constant temperature in the updraft
!
! This computation yields an estimation of the buoyancy (zbuoy) and thereforce
! an estimation of the velocity squared (w_est)
!---------------------------------------------------------------------------

        do ig=1,ngrid
          if(activecell(ig)) then
            ztva_est(ig,l)=ztla(ig,l-1)

            zdz=zlev(ig,l+1)-zlev(ig,l)
            zbuoy(ig,l)=g*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)

            ! Estimated vertical velocity squared 
            ! (discretized version of equation 12 in paragraph 40 of paper)

            if (((a1*zbuoy(ig,l)/w_est(ig,l)-b1) .gt. 0.) .and. (w_est(ig,l) .ne. 0.)) then
              w_est(ig,l+1)=Max(0.0001,w_est(ig,l)+2.*zdz*a1*zbuoy(ig,l)-2.*zdz*w_est(ig,l)*b1 &
     & -2.*(1.-omega)*zdz*w_est(ig,l)*ae*(a1*zbuoy(ig,l)/w_est(ig,l)-b1)**be)
            else
              w_est(ig,l+1)=Max(0.0001,w_est(ig,l)+2.*zdz*a1inv*zbuoy(ig,l)-2.*zdz*w_est(ig,l)*b1inv)
            endif
            if (w_est(ig,l+1).lt.0.) then
              w_est(ig,l+1)=zw2(ig,l)
            endif
          endif ! of if(activecell(ig))
        enddo ! of do ig=1,ngrid

!-------------------------------------------------
! Compute corresponding non-dimensional (ND) entrainment and detrainment rates
!-------------------------------------------------

        do ig=1,ngrid
         if (activecell(ig)) then

          zw2m=w_est(ig,l+1)
          zdz=zlev(ig,l+1)-zlev(ig,l)

          if((a1*(zbuoy(ig,l)/zw2m)-b1).gt.0.) then 

         ! ND entrainment rate, see equation 16 of paper (paragraph 43)

            entr_star(ig,l)=f_star(ig,l)*zdz*  &
        &   MAX(0.,ae*(a1*(zbuoy(ig,l)/zw2m)-b1)**be)
          
          else
            entr_star(ig,l)=0.
          endif

          if(zbuoy(ig,l) .gt. 0.) then
             if(l .lt. lalim(ig)) then

                detr_star(ig,l)=0.
             else

         ! ND detrainment rate, see paragraph 44 of paper

                detr_star(ig,l) = f_star(ig,l)*zdz*ad

             endif
          else
            detr_star(ig,l)=f_star(ig,l)*zdz*                        &
                &     MAX(ad,bd*zbuoy(ig,l)/zw2m)

          endif

! If we are still in the source layer, we define the source layer entr. rate  (alim_star) as the
! maximum between the source entrainment rate and the estimated entrainment rate.

          if (l.lt.lalim(ig)) then
            alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
            entr_star(ig,l)=0. 
          endif

! Compute the non-dimensional upward mass flux at layer l+1
! using equation 11 of appendix 4.2 in paper

            f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

          endif ! of if (activecell(ig))
        enddo ! of do ig=1,ngrid

! -----------------------------------------------------------------------------------
!
! .II. CONVERGENCE LOOP
!
! We have estimated a vertical velocity profile and refined the source layer profile 
! We now conduct iterations to compute:
!
! - the temperature inside the updraft from the estimated entrainment/source, detrainment,
! and upward mass flux.
! - the buoyancy from the new temperature inside the updraft
! - the vertical velocity from the new buoyancy
! - the entr., detr. and upward mass flux from the new buoyancy and vertical velocity
!
! This loop (tic) converges quickly. We have hardcoded 6 iterations from empirical observations.
! Convergence occurs in 1 or 2 iterations in most cases.
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------
      DO tic=0,5 ! internal convergence loop
! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------

! Is the cell still active ?
      activetmp(:)=activecell(:) .and. f_star(:,l+1)>1.e-10

! If the cell is active, compute temperature inside updraft
      do ig=1,ngrid
       if (activetmp(ig)) then

           ztla(ig,l)=(f_star(ig,l)*ztla(ig,l-1)+  &
     &            (alim_star(ig,l)+entr_star(ig,l))*ztv(ig,l))  &
     &            /(f_star(ig,l+1)+detr_star(ig,l))
       endif
      enddo

! Is the cell still active with respect to temperature variations ?
      activetmp(:)=activetmp(:).and.(abs(ztla(:,l)-ztva(:,l)).gt.0.01)

! Compute new buoyancy and vertical velocity
      do ig=1,ngrid
        zdz=zlev(ig,l+1)-zlev(ig,l)
        if (activetmp(ig)) then 
           ztva(ig,l) = ztla(ig,l)
           zbuoy(ig,l)=g*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)

          ! (discretized version of equation 12 in paragraph 40 of paper)
           if (((a1*zbuoy(ig,l)/zw2(ig,l)-b1) .gt. 0.) .and. &
                (zw2(ig,l) .ne. 0.) ) then
             zw2(ig,l+1)=Max(0.,zw2(ig,l)+2.*zdz*a1*zbuoy(ig,l)-       &
                             2.*zdz*zw2(ig,l)*b1-2.*(1.-omega)*zdz*zw2(ig,l)* &
                             ae*(a1*zbuoy(ig,l)/zw2(ig,l)-b1)**be)
           else
             zw2(ig,l+1)=Max(0.,zw2(ig,l)+2.*zdz*a1inv*zbuoy(ig,l) &
                             -2.*zdz*zw2(ig,l)*b1inv)
           endif
        endif
      enddo

! ================ RECOMPUTE ENTR, DETR, and F FROM NEW W2 ===================
         ! ND entrainment rate, see equation 16 of paper (paragraph 43)
         ! ND detrainment rate, see paragraph 44 of paper

      do ig=1,ngrid
        if (activetmp(ig)) then

          zw2m=zw2(ig,l+1)
          zdz=zlev(ig,l+1)-zlev(ig,l)
          if(zw2m .gt. 0) then
            if((a1*(zbuoy(ig,l)/zw2m)-b1) .gt. 0.) then
              entr_star(ig,l)=f_star(ig,l)*zdz*  &
        &        MAX(0.,ae*(a1*(zbuoy(ig,l)/zw2m)-b1)**be)
            else
              entr_star(ig,l)=0.
            endif

            if(zbuoy(ig,l) .gt. 0.) then
              if(l .lt. lalim(ig)) then

                detr_star(ig,l)=0.

              else
                 detr_star(ig,l) = f_star(ig,l)*zdz*ad

              endif
            else
              detr_star(ig,l)=f_star(ig,l)*zdz*                   &
                &     MAX(ad,bd*zbuoy(ig,l)/zw2m)

            endif
          else
            entr_star(ig,l)=0.
            detr_star(ig,l)=0.
          endif ! of if(zw2m .gt. 0)

! If we are still in the source layer, we define the source layer entr. rate  (alim_star) as the
! maximum between the source entrainment rate and the estimated entrainment rate.

        if (l.lt.lalim(ig)) then
          alim_star(ig,l)=max(alim_star(ig,l),entr_star(ig,l))
          entr_star(ig,l)=0.
        endif

! Compute the non-dimensional upward mass flux at layer l+1
! using equation 11 of appendix 4.2 in paper

        f_star(ig,l+1)=f_star(ig,l)+alim_star(ig,l)+entr_star(ig,l)  &
     &              -detr_star(ig,l)

        endif ! of if (activetmp(ig))
      enddo ! of do ig=1,ngrid
! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------
      ENDDO   ! of internal convergence loop DO tic=0,5
! -----------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------

!---------------------------------------------------------------------------
! Miscellaneous computations for height 
!---------------------------------------------------------------------------

      do ig=1,ngrid
        if (zw2(ig,l+1)>0. .and. zw2(ig,l+1).lt.1.e-10) then
          IF (thermverbose) THEN
           print*,'thermcell_plume, particular case in velocity profile'
          ENDIF
          zw2(ig,l+1)=0.
        endif

        if (zw2(ig,l+1).lt.0.) then
           zw2(ig,l+1)=0.
        endif
        wa_moy(ig,l+1)=sqrt(zw2(ig,l+1))

        if (wa_moy(ig,l+1).gt.wmaxa(ig)) then
            wmaxa(ig)=wa_moy(ig,l+1)
        endif
      enddo

!=========================================================================
!=========================================================================
!=========================================================================
! END OF THE LOOP ON VERTICAL LEVELS
      enddo ! of do l=2,nlayer-1
!=========================================================================
!=========================================================================
!=========================================================================

! Recompute the source layer total entrainment alim_star_tot
! as alim_star may have been modified in the above loop. Renormalization of
! alim_star.

       do ig=1,ngrid
          alim_star_tot(ig)=0.
       enddo
       do ig=1,ngrid
          do l=1,lalim(ig)-1
          alim_star_tot(ig)=alim_star_tot(ig)+alim_star(ig,l)
          enddo
       enddo

      do l=1,nlayer
         do ig=1,ngrid
            if (alim_star_tot(ig) > 1.e-10 ) then
               alim_star(ig,l)=alim_star(ig,l)/alim_star_tot(ig)
            endif
         enddo
      enddo

! ===========================================================================
! ================= FIN PLUME ===============================================
! ===========================================================================

! ===========================================================================
! ================= HEIGHT ==================================================
! ===========================================================================

! WARNING, W2 (squared velocity) IS TRANSFORMED IN ITS SQUARE ROOT HERE

!-------------------------------------------------------------------------------
! Computations of the thermal height zmax and maximum vertical velocity wmax
!-------------------------------------------------------------------------------

! Index of the thermal plume height
      do ig=1,ngrid
         lmax(ig)=lalim(ig)
      enddo
      do ig=1,ngrid
         do l=nlayer,lalim(ig)+1,-1
            if (zw2(ig,l).le.1.e-10) then
               lmax(ig)=l-1 
            endif
         enddo
      enddo

! Particular case when the thermal reached the model top, which is not a good sign
      do ig=1,ngrid
      if ( zw2(ig,nlayer) > 1.e-10 ) then
          print*,'thermcell_main_mars: WARNING !!!!! W2 non-zero in last layer for ig=',ig
          lmax(ig)=nlayer
      endif
      enddo

! Maximum vertical velocity zw2
      do ig=1,ngrid
         wmax(ig)=0.
      enddo

      do l=1,nlayer
         do ig=1,ngrid
            if (l.le.lmax(ig)) then
                if (zw2(ig,l).lt.0.)then
!                  print*,'pb2 zw2<0',zw2(ig,l)
                  zw2(ig,l)=0.
                endif
                zw2(ig,l)=sqrt(zw2(ig,l))
                wmax(ig)=max(wmax(ig),zw2(ig,l))
            else
                 zw2(ig,l)=0.
            endif
          enddo
      enddo

! Height of the thermal plume, defined as the following:
! zmax=Integral[z*w(z)*dz]/Integral[w(z)*dz]
!
      do  ig=1,ngrid
         zmax(ig)=0.
         zlevinter(ig)=zlev(ig,1)
      enddo

         num(:)=0.
         denom(:)=0.
         do ig=1,ngrid
          do l=1,nlayer
             num(ig)=num(ig)+zw2(ig,l)*zlev(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
             denom(ig)=denom(ig)+zw2(ig,l)*(zlev(ig,l+1)-zlev(ig,l))
          enddo
       enddo
       do ig=1,ngrid
       if (denom(ig).gt.1.e-10) then
          zmax(ig)=2.*num(ig)/denom(ig)
       endif
       enddo

! ===========================================================================
! ================= FIN HEIGHT ==============================================
! ===========================================================================

#ifdef MESOSCALE
      limz= nlayer-5 ! the most important is limz > max(PBLheight)+2
                     ! nlayer-5 is more than enough!
#else
      call planetwide_maxval(lmax,limz)
      limz=limz+2
#endif
      
      if (limz .ge. nlayer) then
        print*,'thermals have reached last layer of the model'
        print*,'this is not good !'
        limz=nlayer
      endif
! alim_star_clos is the source profile used for closure. It consists of the
! modified source profile in the source layers, and the entrainment profile 
! above it.

      alim_star_clos(:,:)=entr_star(:,:)+alim_star(:,:)

! ===========================================================================
! ============= CLOSURE =====================================================
! ===========================================================================

!-------------------------------------------------------------------------------
! Closure, determination of the upward mass flux
!-------------------------------------------------------------------------------
! Init.

      alim_star2(:)=0.
      alim_star_tot_clos(:)=0.
      f(:)=0.

! llmax is the index of the heighest thermal in the simulation domain
#ifdef MESOSCALE
      !! AS: THIS IS PARALLEL SENSITIVE!!!!! to be corrected?
      llmax=1
      do ig=1,ngrid
         if (lalim(ig)>llmax) llmax=lalim(ig)
      enddo
#else
      call planetwide_maxval(lalim,llmax)
#endif

! Integral of a**2/(rho* Delta z), see equation 13 of appendix 4.2 in paper

      do k=1,llmax-1
         do ig=1,ngrid
            if (k<lalim(ig)) then
         alim_star2(ig)=alim_star2(ig)+alim_star_clos(ig,k)*alim_star_clos(ig,k)  &
      &                    /(rho(ig,k)*(zlev(ig,k+1)-zlev(ig,k)))
         alim_star_tot_clos(ig)=alim_star_tot_clos(ig)+alim_star_clos(ig,k)
      endif
         enddo
      enddo
 
! Closure mass flux, equation 13 of appendix 4.2 in paper

      do ig=1,ngrid
         if (alim_star2(ig)>1.e-10) then
             f(ig)=wmax(ig)*alim_star_tot_clos(ig)/  &
      &     (max(500.,zmax(ig))*r_aspect_thermals*alim_star2(ig))

         endif
      enddo

! ===========================================================================
! ============= FIN CLOSURE =================================================
! ===========================================================================


! ===========================================================================
! ============= FLUX2 =======================================================
! ===========================================================================

!-------------------------------------------------------------------------------
! With the closure mass flux, we can compute the entrainment, detrainment and
! upward mass flux from the non-dimensional ones.
!-------------------------------------------------------------------------------

      fomass_max=0.8 !maximum mass fraction of a cell that can go upward in an
                     ! updraft
      alphamax=0.5 !maximum updraft coverage in a cell


!    these variables allow to follow corrections made to the mass flux when thermverbose=.true.
      ncorecfm1=0
      ncorecfm2=0
      ncorecfm3=0
      ncorecfm4=0
      ncorecfm5=0
      ncorecfm6=0
      ncorecfm7=0
      ncorecfm8=0
      ncorecalpha=0

!-------------------------------------------------------------------------
! Multiply by the closure mass flux
!-------------------------------------------------------------------------

      do l=1,limz
         entr(:,l)=f(:)*(entr_star(:,l)+alim_star(:,l))
         detr(:,l)=f(:)*detr_star(:,l)
      enddo

! Reconstruct the updraft mass flux everywhere

      do l=1,limz
         do ig=1,ngrid
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now we will reconstruct once again the upward
! mass flux, but we will apply corrections
! in some cases. We can compare to the
! previously computed mass flux (above)
!
! This verification is done level by level
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do l=1,limz !loop on the levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Upward mass flux at level l+1

         do ig=1,ngrid
            if (l.lt.lmax(ig)) then
               fm(ig,l+1)=fm(ig,l)+entr(ig,l)-detr(ig,l)
            elseif(l.eq.lmax(ig)) then
               fm(ig,l+1)=0.
               detr(ig,l)=fm(ig,l)+entr(ig,l)
            else
               fm(ig,l+1)=0.
            endif
         enddo


!-------------------------------------------------------------------------
! Upward mass flux should be positive
!-------------------------------------------------------------------------

         do ig=1,ngrid

            if (fm(ig,l+1).lt.0.) then
               if((l+1) .eq. lmax(ig)) then
               detr(ig,l)=detr(ig,l)+fm(ig,l+1)
               fm(ig,l+1)=0.
               entr(ig,l+1)=0.
               ncorecfm2=ncorecfm2+1
               else
               IF (thermverbose) THEN
          print*,'fm(l+1)<0 : ig, l+1,lmax :',ig,l+1,lmax(ig),fm(ig,l+1)
               ENDIF
               ncorecfm1=ncorecfm1+1
               fm(ig,l+1)=fm(ig,l)
               detr(ig,l)=entr(ig,l)
               endif
            endif

         enddo

!-------------------------------------------------------------------------
! Detrainment should be lower than upward mass flux
!-------------------------------------------------------------------------

         do ig=1,ngrid
            if (detr(ig,l).gt.fm(ig,l)) then
               ncorecfm6=ncorecfm6+1
               detr(ig,l)=fm(ig,l)
               entr(ig,l)=fm(ig,l+1)

! When detrainment is stronger than upward mass flux, and we are above the
! thermal last level, the plume is stopped

            if(l.gt.lmax(ig)) then
               detr(ig,l)=0.
               fm(ig,l+1)=0.
               entr(ig,l)=0.
            endif
            
            endif

         enddo

!-------------------------------------------------------------------------
! Check again for mass flux positivity
!-------------------------------------------------------------------------

         do ig=1,ngrid
            if (fm(ig,l+1).lt.0.) then
               detr(ig,l)=detr(ig,l)+fm(ig,l+1)
               fm(ig,l+1)=0.
               ncorecfm2=ncorecfm2+1
            endif
         enddo

!-----------------------------------------------------------------------
! Fractional coverage should be less than 1
!-----------------------------------------------------------------------

        do ig=1,ngrid
           if (zw2(ig,l+1).gt.1.e-10) then
           zfm=rhobarz(ig,l+1)*zw2(ig,l+1)*alphamax
           if ( fm(ig,l+1) .gt. zfm) then
              f_old=fm(ig,l+1)
              fm(ig,l+1)=zfm
              detr(ig,l)=detr(ig,l)+f_old-fm(ig,l+1)
              ncorecalpha=ncorecalpha+1
           endif
           endif

        enddo

      enddo ! on vertical levels

!-----------------------------------------------------------------------
!
! We limit the total mass going from one level to the next, compared to the
! initial total mass fo the cell
!
!-----------------------------------------------------------------------

      do l=1,limz
         do ig=1,ngrid
            eee0=entr(ig,l)
            ddd0=detr(ig,l)
            eee=entr(ig,l)-masse(ig,l)*fomass_max/ptimestep
            ddd=detr(ig,l)-eee
            if (eee.gt.0.) then
                ncorecfm3=ncorecfm3+1
                entr(ig,l)=entr(ig,l)-eee
                if ( ddd.gt.0.) then
! The entrainment is too strong but we can compensate the excess by a detrainment decrease
                   detr(ig,l)=ddd
                else
! The entrainment is too strong and we compensate the excess by a stronger entrainment
! in the layer above
                   if(l.eq.lmax(ig)) then
                      detr(ig,l)=fm(ig,l)+entr(ig,l)
                   else
                      entr(ig,l+1)=entr(ig,l+1)-ddd
                      detr(ig,l)=0.
                      fm(ig,l+1)=fm(ig,l)+entr(ig,l)
                      detr(ig,l)=0.
                   endif
                endif
            endif
         enddo
      enddo

! Check again that everything cancels at zmax
      do ig=1,ngrid
         fm(ig,lmax(ig)+1)=0.
         entr(ig,lmax(ig))=0.
         detr(ig,lmax(ig))=fm(ig,lmax(ig))+entr(ig,lmax(ig))
      enddo

!-----------------------------------------------------------------------
! Summary of the number of modifications that were necessary (if thermverbose=.true.
! and only if there were a lot of them)
!-----------------------------------------------------------------------

!IM 090508 beg
      IF (thermverbose) THEN
      if (ncorecfm1+ncorecfm2+ncorecfm3+ncorecfm4+ncorecfm5+ncorecalpha > ngrid/4. ) then
         print*,'thermcell warning : large number of corrections'
         print*,'PB thermcell : on a du coriger ',ncorecfm1,'x fm1',&
     &     ncorecfm2,'x fm2',ncorecfm3,'x fm3 et', &
     &     ncorecfm4,'x fm4',ncorecfm5,'x fm5 et', &
     &     ncorecfm6,'x fm6', &
     &     ncorecfm7,'x fm7', &
     &     ncorecfm8,'x fm8', &
     &     ncorecalpha,'x alpha'
      endif
      ENDIF

! ===========================================================================
! ============= FIN FLUX2 ===================================================
! ===========================================================================


! ===========================================================================
! ============= TRANSPORT ===================================================
! ===========================================================================

!------------------------------------------------------------------
! vertical transport computation
!------------------------------------------------------------------

! ------------------------------------------------------------------
! IN THE UPDRAFT
! ------------------------------------------------------------------

      zdthladj(:,:)=0.
! Based on equation 14 in appendix 4.2

      do ig=1,ngrid
         if(lmax(ig) .gt. 1) then
         do k=1,lmax(ig)
            zdthladj(ig,k)=(1./masse(ig,k))*(fm(ig,k+1)*ztv(ig,k+1)-    &
     &   fm(ig,k)*ztv(ig,k)+fm(ig,k)*ztva(ig,k)-fm(ig,k+1)*ztva(ig,k+1))
            if (ztv(ig,k) + ptimestep*zdthladj(ig,k) .le. 0.) then
      IF (thermverbose) THEN
              print*,'Teta<0 in thermcell_dTeta up: qenv .. dq : ', ztv(ig,k),ptimestep*zdthladj(ig,k)
      ENDIF
              if(ztv(ig,k) .gt. 0.) then
              zdthladj(ig,k)=0.
              endif
            endif
         enddo
         endif
      enddo

! ------------------------------------------------------------------
! DOWNDRAFT PARAMETERIZATION
! ------------------------------------------------------------------

      ztvd(:,:)=ztv(:,:)
      fm_down(:,:)=0.
      do ig=1,ngrid
         if (lmax(ig) .gt. 1) then
         do l=1,lmax(ig)
              if(zlay(ig,l) .le. zmax(ig)) then

! see equation 18 of paragraph 48 in paper
              fm_down(ig,l) =fm(ig,l)* &
     &      max(fdfu,-4*max(0.,(zlay(ig,l)/zmax(ig)))-0.6)
              endif

            if(zlay(ig,l) .le. zmax(ig)) then            
! see equation 19 of paragraph 49 in paper
          ztvd(ig,l)=min(ztv(ig,l),ztv(ig,l)*((zlay(ig,l)/zmax(ig))/400. + 0.997832))
             else
          ztvd(ig,l)=ztv(ig,l)
            endif

         enddo
         endif
      enddo

! ------------------------------------------------------------------
! TRANSPORT IN DOWNDRAFT
! ------------------------------------------------------------------

       zdthladj_down(:,:)=0.

      do ig=1,ngrid
         if(lmax(ig) .gt. 1) then
! No downdraft in the very-near surface layer, we begin at k=3
! Based on equation 14 in appendix 4.2
 
         do k=3,lmax(ig)
            zdthladj_down(ig,k)=(1./masse(ig,k))*(fm_down(ig,k+1)*ztv(ig,k+1)- &
     & fm_down(ig,k)*ztv(ig,k)+fm_down(ig,k)*ztvd(ig,k)-fm_down(ig,k+1)*ztvd(ig,k+1))
            if (ztv(ig,k) + ptimestep*zdthladj_down(ig,k) .le. 0.) then
      IF (thermverbose) THEN
              print*,'q<0 in thermcell_dTeta down: qenv .. dq : ', ztv(ig,k),ptimestep*zdthladj_down(ig,k)
      ENDIF
              if(ztv(ig,k) .gt. 0.) then
              zdthladj(ig,k)=0.
              endif
            endif
         enddo
         endif
      enddo

!------------------------------------------------------------------
! Final fraction coverage with the clean upward mass flux, computed at interfaces
!------------------------------------------------------------------
      fraca(:,:)=0.
      do l=2,limz
         do ig=1,ngrid
            if (zw2(ig,l).gt.1.e-10) then
            fraca(ig,l)=fm(ig,l)/(rhobarz(ig,l)*zw2(ig,l))
            else
            fraca(ig,l)=0.
            endif
         enddo
      enddo

!------------------------------------------------------------------
! Transport of C02 Tracer
!------------------------------------------------------------------

! We only transport co2 tracer because it is coupled to the scheme through theta_m
! The rest is transported outside the sub-timestep loop

      ratiom(:,:)=1.

      if (igcm_co2.ne.0) then
      detrmod(:,:)=0.
      do k=1,limz
         do ig=1,ngrid
            detrmod(ig,k)=fm(ig,k)-fm(ig,k+1) &
     &      +entr(ig,k)
            if (detrmod(ig,k).lt.0.) then
               entr(ig,k)=entr(ig,k)-detrmod(ig,k)
               detrmod(ig,k)=0.
            endif
         enddo
      enddo

      call thermcell_dqup(ngrid,nlayer,ptimestep     &
     &     ,fm,entr,detrmod,  &
     &    masse,pq(:,:,igcm_co2),pdqadj(:,:,igcm_co2),ndt,limz)

! Compute the ratio between theta and theta_m
      
       do l=1,limz
          do ig=1,ngrid
             ratiom(ig,l)=1./(A*(pq(ig,l,igcm_co2)+pdqadj(ig,l,igcm_co2)*ptimestep)+B)
          enddo
       enddo

       endif

!------------------------------------------------------------------
!  incrementation dt
!------------------------------------------------------------------

      pdtadj(:,:)=0.
      do l=1,limz
         do ig=1,ngrid
         pdtadj(ig,l)=(zdthladj(ig,l)+zdthladj_down(ig,l))*zpopsk(ig,l)*ratiom(ig,l)
         enddo
      enddo

! ===========================================================================
! ============= FIN TRANSPORT ===============================================
! ===========================================================================


!------------------------------------------------------------------
!   Diagnostics for outputs
!------------------------------------------------------------------
! We compute interface values for teta env and th. The last interface
! value does not matter, as the mass flux is 0 there.

    
      do l=1,nlayer-1
       do ig=1,ngrid
        teta_th_int(ig,l)=0.5*(ztva(ig,l+1)+ztva(ig,l))*ratiom(ig,l)
        teta_down_int(ig,l) = 0.5*(ztvd(ig,l+1)+ztvd(ig,l))*ratiom(ig,l)
        teta_env_int(ig,l)=0.5*(ztv(ig,l+1)+ztv(ig,l))*ratiom(ig,l)
       enddo
      enddo
      do ig=1,ngrid
       teta_th_int(ig,nlayer)=teta_th_int(ig,nlayer-1)
       teta_env_int(ig,nlayer)=teta_env_int(ig,nlayer-1)
       teta_down_int(ig,nlayer)=teta_down_int(ig,nlayer-1)
      enddo
        heatFlux(:,:)=0.
        buoyancyOut(:,:)=0.
        buoyancyEst(:,:)=0.
        heatFlux_down(:,:)=0.
      do l=1,limz
       do ig=1,ngrid
        heatFlux(ig,l)=fm(ig,l)*(teta_th_int(ig,l)-teta_env_int(ig,l))/(rhobarz(ig,l))
        buoyancyOut(ig,l)=g*(ztva(ig,l)-ztv(ig,l))/ztv(ig,l)
        buoyancyEst(ig,l)=g*(ztva_est(ig,l)-ztv(ig,l))/ztv(ig,l)
        heatFlux_down(ig,l)=fm_down(ig,l)*(teta_down_int(ig,l)-teta_env_int(ig,l))/rhobarz(ig,l)
       enddo
      enddo

      return
      end
