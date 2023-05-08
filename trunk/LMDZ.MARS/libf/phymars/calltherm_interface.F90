!=======================================================================
! CALLTHERM_INTERFACE
!=======================================================================
! Main interface to the Martian thermal plume model
! This interface handles sub-timesteps for this model
! A call to this interface must be inserted in the main 'physics' routine
!   NB: for information:
!   In the Mars LMD-GCM, the thermal plume model is called after 
!   the vertical turbulent mixing scheme (Mellor and Yamada)  
!   and the surface layer scheme (Richardson-based surface layer + subgrid gustiness)
!   Other routines called before the thermals model are 
!   radiative transfer and (orographic) gravity wave drag.
! -----------------------------------------------------------------------
! Author : A. Colaitis 2011-01-05 (with updates 2011-2013)
!          after C. Rio and F. Hourdin
! Institution : Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------
! Corresponding author : A. Spiga aymeric.spiga_AT_upmc.fr
! -----------------------------------------------------------------------
! ASSOCIATED FILES
! --> thermcell_main_mars.F90
! --> thermcell_dqup.F90
! --> comtherm_h.F90
! -----------------------------------------------------------------------
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

      SUBROUTINE calltherm_interface (ngrid,nlayer,nq, &
     & tracer,igcm_co2, &
     & zzlev,zzlay, &
     & ptimestep,pu,pv,pt,pq,pdu,pdv,pdt,pdq,q2, &
     & pplay,pplev,pphi,zpopsk, &
     & pdu_th,pdv_th,pdt_th,pdq_th,lmax,zmaxth,pbl_dtke, &
     & pdhdif,hfmax,wstar,sensibFlux)

      use comtherm_h
      use tracer_mod, only: nqmx,noms

      ! SHARED VARIABLES. This needs adaptations in another climate model.
      ! contains physical constant values such as
      ! "g" : gravitational acceleration (m.s-2)
      ! "r" : recuced gas constant (J.K-1.mol-1)
      ! "cpp" : specific heat of the atmosphere (J.kg-1.K-1)
      USE comcstfi_h

      implicit none

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------

      INTEGER, INTENT(IN) :: ngrid ! number of horizontal grid points
      INTEGER, INTENT(IN) :: nlayer ! number of vertical grid points
      INTEGER, INTENT(IN) :: nq ! number of tracer species
      REAL, INTENT(IN) :: ptimestep !timestep (s)
      REAL, INTENT(IN) :: pplev(ngrid,nlayer+1) !intermediate pressure levels (Pa)
      REAL, INTENT(IN) :: pplay(ngrid,nlayer) !Pressure at the middle of the layers (Pa)
      REAL, INTENT(IN) :: pphi(ngrid,nlayer) !Geopotential at the middle of the layers (m2s-2)
      REAL, INTENT(IN) :: pu(ngrid,nlayer),pv(ngrid,nlayer) !u,v components of the wind (ms-1)
      REAL, INTENT(IN) :: pt(ngrid,nlayer),pq(ngrid,nlayer,nq)!temperature (K) and tracer concentration (kg/kg)
      REAL, INTENT(IN) :: zzlay(ngrid,nlayer) ! altitude at the middle of the layers
      REAL, INTENT(IN) :: zzlev(ngrid,nlayer+1) ! altitude at layer boundaries
      LOGICAL, INTENT(IN) :: tracer ! =.true. if tracers are present and to be transported
      INTEGER, INTENT(IN) :: igcm_co2 ! index of the CO2 tracer in mixing ratio array
                                      ! --> 0 if no tracer is CO2 (or no tracer at all)
                                      ! --> this prepares special treatment for polar night mixing
                                      !       (see thermcell_main_mars)
      REAL, INTENT(IN) :: pdu(ngrid,nlayer),pdv(ngrid,nlayer) ! wind velocity change from routines called 
                                                                      ! before thermals du/dt (m/s/s)
      REAL, INTENT(IN) :: pdq(ngrid,nlayer,nq) ! tracer concentration change from routines called 
                                                     ! before thermals dq/dt (kg/kg/s)
      REAL, INTENT(IN) :: pdt(ngrid,nlayer) ! temperature change from routines called before thermals dT/dt (K/s)
      REAL, INTENT(IN) :: q2(ngrid,nlayer+1) ! turbulent kinetic energy
      REAL, INTENT(IN) :: zpopsk(ngrid,nlayer) ! ratio of pressure at middle of layer to surface pressure, 
                                                   ! to the power r/cp, i.e. zpopsk=(pplay(ig,l)/pplev(ig,1))**rcp
      REAL, INTENT(IN) :: pdhdif(ngrid,nlayer) ! potential temperature change from turbulent diffusion scheme dT/dt (K/s)
      REAL, INTENT(IN) :: sensibFlux(ngrid) ! sensible heat flux computed from surface layer scheme

!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------

      REAL, INTENT(OUT) :: pdu_th(ngrid,nlayer) ! wind velocity change from thermals du/dt (m/s/s)
      REAL, INTENT(OUT) :: pdv_th(ngrid,nlayer) ! wind velocity change from thermals dv/dt (m/s/s)
      REAL, INTENT(OUT) :: pdt_th(ngrid,nlayer) ! temperature change from thermals dT/dt (K/s)
      REAL, INTENT(OUT) :: pdq_th(ngrid,nlayer,nq) ! tracer change from thermals dq/dt (kg/kg/s)
      INTEGER, INTENT(OUT) :: lmax(ngrid) ! layer number reacher by thermals in grid point
      REAL, INTENT(OUT) :: zmaxth(ngrid) ! equivalent to lmax, but in (m), interpolated
      REAL, INTENT(OUT) :: pbl_dtke(ngrid,nlayer+1) ! turbulent kinetic energy change from thermals dtke/dt
      REAL, INTENT(OUT) :: wstar(ngrid) ! free convection velocity (m/s)

!--------------------------------------------------------
! Thermals local variables
!--------------------------------------------------------
      REAL zu(ngrid,nlayer), zv(ngrid,nlayer)
      REAL zt(ngrid,nlayer)
      REAL d_t_ajs(ngrid,nlayer)
      REAL d_u_ajs(ngrid,nlayer), d_q_ajs(ngrid,nlayer,nq)
      REAL d_v_ajs(ngrid,nlayer) 
      REAL fm_therm(ngrid,nlayer+1), entr_therm(ngrid,nlayer)
      REAL detr_therm(ngrid,nlayer),detrmod(ngrid,nlayer)
      REAL zw2(ngrid,nlayer+1)
      REAL fraca(ngrid,nlayer+1),zfraca(ngrid,nlayer+1)
      REAL q_therm(ngrid,nlayer), pq_therm(ngrid,nlayer,nq)
      REAL q2_therm(ngrid,nlayer), dq2_therm(ngrid,nlayer)
      REAL lmax_real(ngrid)
      REAL masse(ngrid,nlayer)

      INTEGER l,ig,iq,ii(1),k
      CHARACTER (LEN=20) modname

!--------------------------------------------------------
! Local variables for sub-timestep
!--------------------------------------------------------

      REAL d_t_the(ngrid,nlayer), d_q_the(ngrid,nlayer,nq)
      INTEGER isplit
      REAL fact
      REAL zfm_therm(ngrid,nlayer+1),zdt
      REAL zentr_therm(ngrid,nlayer),zdetr_therm(ngrid,nlayer)
      REAL zheatFlux(ngrid,nlayer)
      REAL zheatFlux_down(ngrid,nlayer)
      REAL zbuoyancyOut(ngrid,nlayer)
      REAL zbuoyancyEst(ngrid,nlayer)
      REAL zzw2(ngrid,nlayer+1)
      REAL zmax(ngrid)
      INTEGER ndt,limz

!--------------------------------------------------------
! Diagnostics
!--------------------------------------------------------

      REAL heatFlux(ngrid,nlayer)
      REAL heatFlux_down(ngrid,nlayer)
      REAL buoyancyOut(ngrid,nlayer)
      REAL buoyancyEst(ngrid,nlayer)
      REAL hfmax(ngrid),wmax(ngrid)
      REAL pbl_teta(ngrid),dteta(ngrid,nlayer)
      REAL rpdhd(ngrid,nlayer)
      REAL wtdif(ngrid,nlayer),rho(ngrid,nlayer)
      REAL wtth(ngrid,nlayer)

! **********************************************************************
! Initializations
! **********************************************************************

      lmax(:)=0
      pdu_th(:,:)=0.
      pdv_th(:,:)=0.
      pdt_th(:,:)=0.
      entr_therm(:,:)=0.
      detr_therm(:,:)=0.
      q2_therm(:,:)=0.
      dq2_therm(:,:)=0.
      pbl_dtke(:,:)=0.
      fm_therm(:,:)=0.
      zw2(:,:)=0.
      fraca(:,:)=0.
      zfraca(:,:)=0.
      if (tracer) then
         pdq_th(:,:,:)=0.
      end if
      d_t_ajs(:,:)=0.
      d_u_ajs(:,:)=0.
      d_v_ajs(:,:)=0.
      d_q_ajs(:,:,:)=0.
      heatFlux(:,:)=0.
      heatFlux_down(:,:)=0.
      buoyancyOut(:,:)=0.
      buoyancyEst(:,:)=0.
      zmaxth(:)=0.
      lmax_real(:)=0.


! **********************************************************************
! Preparing inputs for the thermals: increment tendancies 
! from other subroutines called before the thermals model
! **********************************************************************

       zu(:,:)=pu(:,:)+pdu(:,:)*ptimestep ! u-component of wind velocity
       zv(:,:)=pv(:,:)+pdv(:,:)*ptimestep ! v-component of wind velocity
       zt(:,:)=pt(:,:)+pdt(:,:)*ptimestep ! temperature

       pq_therm(:,:,:)=0. ! tracer concentration

       if(qtransport_thermals) then
          if(tracer) then ! tracer is a logical that is true if tracer must be transported in the GCM physics
                pq_therm(:,:,:)=pq(:,:,:)+pdq(:,:,:)*ptimestep ! tracer concentration
          endif
       endif


       IF(dtke_thermals) THEN
          DO l=1,nlayer
              q2_therm(:,l)=0.5*(q2(:,l)+q2(:,l+1))
          ENDDO
       ENDIF

! **********************************************************************
! --> CALLTHERM
! SUB-TIMESTEP LOOP
! **********************************************************************

         zdt=ptimestep/REAL(nsplit_thermals) !subtimestep

         DO isplit=1,nsplit_thermals !temporal loop on the subtimestep

! Initialization of intermediary variables

         zzw2(:,:)=0.
         zmax(:)=0.
         lmax(:)=0

         if (nq .ne. 0 .and. igcm_co2 .ne. 0) then !initialize co2 tracer tendancy
            d_q_the(:,:,igcm_co2)=0.
         endif

! CALL to main thermal routine
             CALL thermcell_main_mars(ngrid,nlayer,nq &
     &      ,tracer,igcm_co2 &
     &      ,zdt  &
     &      ,pplay,pplev,pphi,zzlev,zzlay  &
     &      ,zu,zv,zt,pq_therm,q2_therm  &
     &      ,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm,zdetr_therm,lmax,zmax,limz  &
     &      ,zzw2,fraca,zpopsk &
     &      ,zheatFlux,zheatFlux_down &
     &      ,zbuoyancyOut,zbuoyancyEst)

      fact=1./REAL(nsplit_thermals)

! Update thermals tendancies
            d_t_the(:,:)=d_t_the(:,:)*ptimestep*fact !temperature
            if (igcm_co2 .ne. 0) then
               d_q_the(:,:,igcm_co2)=d_q_the(:,:,igcm_co2)*ptimestep*fact !co2 mass mixing ratio
            endif
            zmaxth(:)=zmaxth(:)+zmax(:)*fact !thermals height
            lmax_real(:)=lmax_real(:)+float(lmax(:))*fact !thermals height index
            fm_therm(:,:)=fm_therm(:,:)  & !upward mass flux
     &      +zfm_therm(:,:)*fact
            entr_therm(:,:)=entr_therm(:,:)  & !entrainment mass flux
     &       +zentr_therm(:,:)*fact
            detr_therm(:,:)=detr_therm(:,:)  & !detrainment mass flux
     &       +zdetr_therm(:,:)*fact
            zfraca(:,:)=zfraca(:,:) + fraca(:,:)*fact !updraft fractional coverage
            heatFlux(:,:)=heatFlux(:,:) & !upward heat flux
     &       +zheatFlux(:,:)*fact
            heatFlux_down(:,:)=heatFlux_down(:,:) & !downward heat flux
     &       +zheatFlux_down(:,:)*fact
            buoyancyOut(:,:)=buoyancyOut(:,:) & !plume final buoyancy
     &       +zbuoyancyOut(:,:)*fact
            buoyancyEst(:,:)=buoyancyEst(:,:) & !plume estimated buoyancy used for vertical velocity computation
     &       +zbuoyancyEst(:,:)*fact
            zw2(:,:)=zw2(:,:) + zzw2(:,:)*fact !vertical velocity

! Save tendancies
           d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_the(:,:) !temperature tendancy (delta T)
            if (igcm_co2 .ne. 0) then
               d_q_ajs(:,:,igcm_co2)=d_q_ajs(:,:,igcm_co2)+d_q_the(:,:,igcm_co2) !tracer tendancy (delta q)
            endif

!  Increment temperature and co2 concentration for next pass in subtimestep loop
            zt(:,:) = zt(:,:) + d_t_the(:,:) !temperature
            if (igcm_co2 .ne. 0) then
             pq_therm(:,:,igcm_co2) = &
     &          pq_therm(:,:,igcm_co2) + d_q_the(:,:,igcm_co2) !co2 tracer
            endif


         ENDDO ! isplit
!****************************************************************

! Now that we have computed total entrainment and detrainment, we can
! advect u, v, and q in thermals. (potential temperature and co2 MMR 
! have already been advected in thermcell_main because they are coupled
! to the determination of the thermals caracteristics). This is done 
! separatly because u,v, and q are not used in thermcell_main for
! any thermals-related computation : they are purely passive.

! mass of cells
      do l=1,nlayer
         masse(:,l)=(pplev(:,l)-pplev(:,l+1))/g
      enddo

! recompute detrainment mass flux from entrainment and upward mass flux
! this ensure mass flux conservation
      detrmod(:,:)=0.
      do l=1,limz
         do ig=1,ngrid
            detrmod(ig,l)=fm_therm(ig,l)-fm_therm(ig,l+1) &
     &      +entr_therm(ig,l)
            if (detrmod(ig,l).lt.0.) then
               entr_therm(ig,l)=entr_therm(ig,l)-detrmod(ig,l)
               detrmod(ig,l)=0.
            endif
         enddo
      enddo

      ! u component of wind velocity advection in thermals
      ! result is a derivative (d_u_ajs in m/s/s)
      ndt=10
      call thermcell_dqup(ngrid,nlayer,ptimestep                &
     &      ,fm_therm,entr_therm,detrmod,  &
     &     masse,zu,d_u_ajs,ndt,limz)

      ! v component of wind velocity advection in thermals
      ! result is a derivative (d_v_ajs in m/s/s)
      call thermcell_dqup(ngrid,nlayer,ptimestep    &
     &       ,fm_therm,entr_therm,detrmod,  &
     &     masse,zv,d_v_ajs,ndt,limz)

      ! non co2 tracers advection in thermals
      ! result is a derivative (d_q_ajs in kg/kg/s)

      if (nq .ne. 0.) then
      DO iq=1,nq
      if (iq .ne. igcm_co2) then
      call thermcell_dqup(ngrid,nlayer,ptimestep     &
     &     ,fm_therm,entr_therm,detrmod,  &
     &    masse,pq_therm(:,:,iq),d_q_ajs(:,:,iq),ndt,limz)
      endif
      ENDDO
      endif

      ! tke advection in thermals
      ! result is a tendancy (d_u_ajs in J)
      if (dtke_thermals) then
      call thermcell_dqup(ngrid,nlayer,ptimestep     &
     &     ,fm_therm,entr_therm,detrmod,  &
     &    masse,q2_therm,dq2_therm,ndt,limz)
      endif

      ! compute wmax for diagnostics
      DO ig=1,ngrid
         wmax(ig)=MAXVAL(zw2(ig,:))
      ENDDO

! **********************************************************************
! **********************************************************************
! **********************************************************************
! CALLTHERM END
! **********************************************************************
! **********************************************************************
! **********************************************************************


! **********************************************************************
! Preparing outputs
! **********************************************************************

      do l=1,limz
        pdu_th(:,l)=d_u_ajs(:,l)
        pdv_th(:,l)=d_v_ajs(:,l)
      enddo

           ! if tracers are transported in thermals, update output variables, else these are 0.
           if(qtransport_thermals) then
              if(tracer) then
               do iq=1,nq
                if (iq .ne. igcm_co2) then
                  do l=1,limz
                     pdq_th(:,l,iq)=d_q_ajs(:,l,iq) !non-co2 tracers d_q_ajs are dq/dt (kg/kg/s)
                  enddo
                else
                  do l=1,limz
                     pdq_th(:,l,iq)=d_q_ajs(:,l,iq)/ptimestep !co2 tracer d_q_ajs is dq (kg/kg)
                  enddo
                endif
               enddo
              endif
           endif

           ! if tke is transported in thermals, update output variable, else this is 0.
           IF(dtke_thermals) THEN
              DO l=2,nlayer
                 pbl_dtke(:,l)=0.5*(dq2_therm(:,l-1)+dq2_therm(:,l))
              ENDDO
  
              pbl_dtke(:,1)=0.5*dq2_therm(:,1)
              pbl_dtke(:,nlayer+1)=0.
           ENDIF

           ! update output variable for temperature. d_t_ajs is delta T in (K), pdt_th is dT/dt in (K/s)
           do l=1,limz
              pdt_th(:,l)=d_t_ajs(:,l)/ptimestep
           enddo


! **********************************************************************
! SURFACE LAYER INTERFACE
! Compute the free convection velocity w* scale for surface layer gustiness
! speed parameterization. The computed value of w* will be used at the next
! timestep to modify surface-atmosphere exchange fluxes, because the surface
! layer scheme and diffusion are called BEFORE the thermals. (outside of these
! routines)
! **********************************************************************

! Potential temperature gradient

      dteta(:,nlayer)=0.
      DO l=1,nlayer-1
         DO ig=1, ngrid
            dteta(ig,l) = ((zt(ig,l+1)-zt(ig,l))/zpopsk(ig,l))          &
     &              /(zzlay(ig,l+1)-zzlay(ig,l))
         ENDDO
      ENDDO

! Computation of the PBL mixed layer temperature

      DO ig=1, ngrid
         ii=MINLOC(abs(dteta(ig,1:lmax(ig))))
         pbl_teta(ig) = zt(ig,ii(1))/zpopsk(ig,ii(1))
      ENDDO

! In order to have an accurate w*, we must add the heat flux from the
! diffusion scheme to the computation of the maximum heat flux hfmax
! Here pdhdif is the potential temperature change from the diffusion
! scheme (Mellor and Yamada, see paper section 6, paragraph 57) 

! compute rho as it is after the diffusion

      rho(:,:)=pplay(:,:)                                               &
     & /(r*(pt(:,:)+pdhdif(:,:)*zpopsk(:,:)*ptimestep))

! integrate -rho*pdhdif

      rpdhd(:,:)=0.

      DO ig=1,ngrid
       DO l=1,lmax(ig)
        rpdhd(ig,l)=0.
        DO k=1,l
         rpdhd(ig,l)=rpdhd(ig,l)-rho(ig,k)*pdhdif(ig,k)*                &
     & (zzlev(ig,k+1)-zzlev(ig,k))
        ENDDO
        rpdhd(ig,l)=rpdhd(ig,l)-sensibFlux(ig)/cpp
       ENDDO
      ENDDO

! compute w'theta' (vertical turbulent flux of temperature) from 
! the diffusion scheme

      wtdif(:,:)=rpdhd(:,:)/rho(:,:)

! Now we compute the contribution of the thermals to w'theta':
! compute rho as it is after the thermals

      rho(:,:)=pplay(:,:)/(r*(zt(:,:)))

! integrate -rho*pdhdif

      DO ig=1,ngrid
       DO l=1,lmax(ig)
        rpdhd(ig,l)=0.
        DO k=1,l
         rpdhd(ig,l)=rpdhd(ig,l)-rho(ig,k)*(pdt_th(ig,k)/zpopsk(ig,k))* &
     & (zzlev(ig,k+1)-zzlev(ig,k))
        ENDDO
        rpdhd(ig,l)=rpdhd(ig,l)+                                        &
     &    rho(ig,1)*(heatFlux(ig,1)+heatFlux_down(ig,1))
       ENDDO
      ENDDO
      rpdhd(:,nlayer)=0.

! compute w'teta' from thermals

      wtth(:,:)=rpdhd(:,:)/rho(:,:)

! Add vertical turbulent heat fluxes from the thermals and the diffusion scheme
! and compute the maximum

      DO ig=1,ngrid
        hfmax(ig)=MAXVAL(wtth(ig,:)+wtdif(ig,:))
      ENDDO

! Finally we can compute the free convection velocity scale
! We follow Spiga et. al 2010 (QJRMS)
! ------------

      DO ig=1, ngrid
         IF (zmax(ig) .gt. 0.) THEN
            wstar(ig)=(g*zmaxth(ig)*hfmax(ig)/pbl_teta(ig))**(1./3.)
         ELSE
            wstar(ig)=0.
         ENDIF
      ENDDO

       END
