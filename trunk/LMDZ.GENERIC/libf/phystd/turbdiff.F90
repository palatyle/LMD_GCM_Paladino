      subroutine turbdiff(ngrid,nlay,nq,rnat,          &
          ptimestep,pcapcal,lecrit,                    &   
          pplay,pplev,pzlay,pzlev,pz0,                 &
          pu,pv,pt,ppopsk,pq,ptsrf,pemis,pqsurf,       &
          pdtfi,pdqfi,pfluxsrf,            &
          Pdudif,pdvdif,pdtdif,pdtsrf,sensibFlux,pq2,  &
          pdqdif,pdqevap,pdqsdif,flux_u,flux_v,lastcall)

      use watercommon_h, only : RLVTT, T_h2O_ice_liq, RCPD, mx_eau_sol,Psat_water, Lcpdqsat_water
      use radcommon_h, only : sigma, glat
      use surfdat_h, only: dryness
      use tracer_h, only: igcm_h2o_vap, igcm_h2o_ice
      use comcstfi_mod, only: rcp, g, r, cpp
      use callkeys_mod, only: water,tracer,nosurf
      use turb_mod, only : ustar
#ifdef MESOSCALE
      use comm_wrf, only : comm_LATENT_HF
#endif
      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Turbulent diffusion (mixing) for pot. T, U, V and tracers
!     
!     Implicit scheme
!     We start by adding to variables x the physical tendencies
!     already computed. We resolve the equation:
!
!     x(t+1) =  x(t) + dt * (dx/dt)phys(t)  +  dt * (dx/dt)difv(t+1)
!     
!     Authors
!     ------- 
!     F. Hourdin, F. Forget, R. Fournier (199X)
!     R. Wordsworth, B. Charnay (2010)
!     J. Leconte (2012): To f90
!         - Rewritten the diffusion scheme to conserve total enthalpy
!               by accounting for dissipation of turbulent kinetic energy.
!         - Accounting for lost mean flow kinetic energy should come soon.
!     
!==================================================================

!-----------------------------------------------------------------------
!     declarations
!     ------------

!     arguments
!     ---------
      INTEGER,INTENT(IN) :: ngrid
      INTEGER,INTENT(IN) :: nlay
      REAL,INTENT(IN) :: ptimestep
      REAL,INTENT(IN) :: pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      REAL,INTENT(IN) :: pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
      REAL,INTENT(IN) :: pu(ngrid,nlay),pv(ngrid,nlay)
      REAL,INTENT(IN) :: pt(ngrid,nlay),ppopsk(ngrid,nlay)
      REAL,INTENT(IN) :: ptsrf(ngrid) ! surface temperature (K)
      REAL,INTENT(IN) :: pemis(ngrid)
      REAL,INTENT(IN) :: pdtfi(ngrid,nlay)
      REAL,INTENT(IN) :: pfluxsrf(ngrid)
      REAL,INTENT(OUT) :: pdudif(ngrid,nlay),pdvdif(ngrid,nlay)
      REAL,INTENT(OUT) :: pdtdif(ngrid,nlay)
      REAL,INTENT(OUT) :: pdtsrf(ngrid) ! tendency (K/s) on surface temperature
      REAL,INTENT(OUT) :: sensibFlux(ngrid)
      REAL,INTENT(IN) :: pcapcal(ngrid)
      REAL,INTENT(INOUT) :: pq2(ngrid,nlay+1)
      REAL,INTENT(OUT) :: flux_u(ngrid),flux_v(ngrid)
      REAL,INTENT(IN) :: rnat(ngrid)      
      LOGICAL,INTENT(IN) :: lastcall ! not used

!     Arguments added for condensation
      logical,intent(in) :: lecrit ! not used.
      REAL,INTENT(IN) :: pz0

!     Tracers
!     --------
      integer,intent(in) :: nq 
      real,intent(in) :: pqsurf(ngrid,nq)
      real,intent(in) :: pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq) 
      real,intent(out) :: pdqdif(ngrid,nlay,nq) 
      real,intent(out) :: pdqsdif(ngrid,nq) 
      
!     local
!     -----
      integer ilev,ig,ilay,nlev

      REAL z4st,zdplanck(ngrid)
      REAL zkv(ngrid,nlay+1),zkh(ngrid,nlay+1)
      REAL zcdv(ngrid),zcdh(ngrid)
      REAL zcdv_true(ngrid),zcdh_true(ngrid)
      REAL zu(ngrid,nlay),zv(ngrid,nlay)
      REAL zh(ngrid,nlay),zt(ngrid,nlay)
      REAL ztsrf(ngrid)
      REAL z1(ngrid),z2(ngrid)
      REAL zmass(ngrid,nlay)
      REAL zfluxv(ngrid,nlay),zfluxt(ngrid,nlay),zfluxq(ngrid,nlay)
      REAL zb0(ngrid,nlay)
      REAL zExner(ngrid,nlay),zovExner(ngrid,nlay)
      REAL zcv(ngrid,nlay),zdv(ngrid,nlay)  !inversion coefficient for winds
      REAL zct(ngrid,nlay),zdt(ngrid,nlay)  !inversion coefficient for temperature
      REAL zcq(ngrid,nlay),zdq(ngrid,nlay)  !inversion coefficient for tracers
      REAL zcst1
      REAL zu2!, a
      REAL zcq0(ngrid),zdq0(ngrid)
      REAL zx_alf1(ngrid),zx_alf2(ngrid)

      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)
      
!     Tracers
!     -------
      INTEGER iq
      REAL zq(ngrid,nlay,nq)
      REAL zqnoevap(ngrid,nlay) !special for water case to compute where evaporated water goes.
      REAL pdqevap(ngrid,nlay) !special for water case to compute where evaporated water goes.
      REAL zdmassevap(ngrid)
      REAL rho(ngrid)         ! near-surface air density
      REAL kmixmin

!     Variables added for implicit latent heat inclusion
!     --------------------------------------------------
      real dqsat(ngrid),psat_temp,qsat(ngrid),psat(ngrid)

      integer, save :: ivap, iliq, iliq_surf,iice_surf ! also make liq for clarity on surface...
!$OMP THREADPRIVATE(ivap,iliq,iliq_surf,iice_surf)

      real, parameter :: karman=0.4
      real cd0, roughratio

      real dqsdif_total(ngrid) 
      real zq0(ngrid) 


!     Coherence test
!     --------------

      IF (firstcall) THEN

         if(water)then
             ivap=igcm_h2o_vap
             iliq=igcm_h2o_ice
             iliq_surf=igcm_h2o_vap
	     iice_surf=igcm_h2o_ice ! simply to make the code legible               
                                  ! to be generalised
	 else
             ivap=0
             iliq=0
             iliq_surf=0
	     iice_surf=0 ! simply to make the code legible               	  
         endif
         sensibFlux(:)=0.

         firstcall=.false.
      ENDIF

!-----------------------------------------------------------------------
!     1. Initialisation
!     -----------------

      nlev=nlay+1

!     Calculate rho*dz, (P/Ps)**(R/cp) and dt*rho/dz=dt*rho**2 g/dp
!     with rho=p/RT=p/ (R Theta) (p/ps)**kappa
!     ---------------------------------------------

      DO ilay=1,nlay
         DO ig=1,ngrid
            zmass(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/glat(ig)
	    zExner(ig,ilay)=(pplev(ig,ilay)/pplev(ig,1))**rcp
	    zovExner(ig,ilay)=1./ppopsk(ig,ilay)
         ENDDO
      ENDDO

      zcst1=4.*g*ptimestep/(R*R)
      DO ilev=2,nlev-1
         DO ig=1,ngrid
            zb0(ig,ilev)=pplev(ig,ilev)/(pt(ig,ilev-1)+pt(ig,ilev))
            zb0(ig,ilev)=zcst1*zb0(ig,ilev)*zb0(ig,ilev)/(pplay(ig,ilev-1)-pplay(ig,ilev))
         ENDDO
      ENDDO
      DO ig=1,ngrid
         zb0(ig,1)=ptimestep*pplev(ig,1)/(R*ptsrf(ig))
      ENDDO
      dqsdif_total(:)=0.0

!-----------------------------------------------------------------------
!     2. Add the physical tendencies computed so far
!     ----------------------------------------------

      DO ilev=1,nlay
         DO ig=1,ngrid
            zu(ig,ilev)=pu(ig,ilev)
            zv(ig,ilev)=pv(ig,ilev)
            zt(ig,ilev)=pt(ig,ilev)+pdtfi(ig,ilev)*ptimestep
            zh(ig,ilev)=pt(ig,ilev)*zovExner(ig,ilev) !for call vdif_kc, but could be moved and computed there
        ENDDO
      ENDDO
      if(tracer) then
         DO iq =1, nq
            DO ilev=1,nlay
               DO ig=1,ngrid
                  zq(ig,ilev,iq)=pq(ig,ilev,iq) + pdqfi(ig,ilev,iq)*ptimestep
               ENDDO
            ENDDO
         ENDDO
         if (water) then
            DO ilev=1,nlay
               DO ig=1,ngrid
                  zqnoevap(ig,ilev)=pq(ig,ilev,ivap) + pdqfi(ig,ilev,ivap)*ptimestep
               ENDDO
            ENDDO
         Endif
      end if

!-----------------------------------------------------------------------
!     3. Turbulence scheme
!     --------------------
!
!     Source of turbulent kinetic energy at the surface
!     ------------------------------------------------- 
!     Formula is Cd_0 = (karman / log[1+z1/z0])^2

      DO ig=1,ngrid
         roughratio = 1. + pzlay(ig,1)/pz0
         cd0 = karman/log(roughratio)
         cd0 = cd0*cd0
         zcdv_true(ig) = cd0
         zcdh_true(ig) = cd0
       if(nosurf)then
         zcdv_true(ig)=0.D+0 !JL12 disable atm/surface momentum flux
         zcdh_true(ig)=0.D+0 !JL12 disable sensible heat flux 
       endif
      ENDDO

      DO ig=1,ngrid
         zu2=pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1)
         zcdv(ig)=zcdv_true(ig)*sqrt(zu2)
         zcdh(ig)=zcdh_true(ig)*sqrt(zu2)
      ENDDO

!     Turbulent diffusion coefficients in the boundary layer
!     ------------------------------------------------------ 

      call vdif_kc(ngrid,nlay,ptimestep,g,pzlev,pzlay,pu,pv,zh,zcdv_true,pq2,zkv,zkh) !JL12 why not call vdif_kc with updated winds and temperature 
      
!     Adding eddy mixing to mimic 3D general circulation in 1D
!     R. Wordsworth & F. Forget (2010)
      if ((ngrid.eq.1)) then
         kmixmin = 1.0e-2       ! minimum eddy mix coeff in 1D
         do ilev=1,nlay
            do ig=1,ngrid
               zkh(ig,ilev) = max(kmixmin,zkh(ig,ilev))
               zkv(ig,ilev) = max(kmixmin,zkv(ig,ilev)) 
            end do
         end do
      end if

!JL12 change zkv at the surface by zcdv to calculate the surface momentum flux properly
      DO ig=1,ngrid
         zkv(ig,1)=zcdv(ig)
      ENDDO
!we treat only winds, energy and tracers coefficients will be computed with upadted winds
!JL12 calculate the flux coefficients (tables multiplied elements by elements)
      zfluxv(1:ngrid,1:nlay)=zkv(1:ngrid,1:nlay)*zb0(1:ngrid,1:nlay)
      
!-----------------------------------------------------------------------
!     4. Implicit inversion of u
!     --------------------------

!     u(t+1) =  u(t) + dt * {(du/dt)phys}(t)  +  dt * {(du/dt)difv}(t+1)
!     avec
!     /zu/ = u(t) + dt * {(du/dt)phys}(t)   (voir paragraphe 2.)
!     et
!     dt * {(du/dt)difv}(t+1) = dt * {(d/dz)[ Ku (du/dz) ]}(t+1)
!     donc les entrees sont /zcdv/ pour la condition a la limite sol
!     et /zkv/ = Ku

      DO ig=1,ngrid
         z1(ig)=1./(zmass(ig,nlay)+zfluxv(ig,nlay))
         zcv(ig,nlay)=zmass(ig,nlay)*zu(ig,nlay)*z1(ig)
         zdv(ig,nlay)=zfluxv(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(zmass(ig,ilay)+zfluxv(ig,ilay) + zfluxv(ig,ilay+1)*(1.-zdv(ig,ilay+1)))
            zcv(ig,ilay)=(zmass(ig,ilay)*zu(ig,ilay)+zfluxv(ig,ilay+1)*zcv(ig,ilay+1))*z1(ig)
            zdv(ig,ilay)=zfluxv(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zu(ig,1)=zcv(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zu(ig,ilay)=zcv(ig,ilay)+zdv(ig,ilay)*zu(ig,ilay-1)
         ENDDO
      ENDDO

!-----------------------------------------------------------------------
!     5. Implicit inversion of v
!     --------------------------

!     v(t+1) =  v(t) + dt * {(dv/dt)phys}(t)  +  dt * {(dv/dt)difv}(t+1)
!     avec
!     /zv/ = v(t) + dt * {(dv/dt)phys}(t)   (voir paragraphe 2.)
!     et
!     dt * {(dv/dt)difv}(t+1) = dt * {(d/dz)[ Kv (dv/dz) ]}(t+1)
!     donc les entrees sont /zcdv/ pour la condition a la limite sol
!     et /zkv/ = Kv

      DO ig=1,ngrid
         z1(ig)=1./(zmass(ig,nlay)+zfluxv(ig,nlay)) 
         zcv(ig,nlay)=zmass(ig,nlay)*zv(ig,nlay)*z1(ig)
         zdv(ig,nlay)=zfluxv(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(zmass(ig,ilay)+zfluxv(ig,ilay)+zfluxv(ig,ilay+1)*(1.-zdv(ig,ilay+1)))
            zcv(ig,ilay)=(zmass(ig,ilay)*zv(ig,ilay)+zfluxv(ig,ilay+1)*zcv(ig,ilay+1))*z1(ig)
            zdv(ig,ilay)=zfluxv(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zv(ig,1)=zcv(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zv(ig,ilay)=zcv(ig,ilay)+zdv(ig,ilay)*zv(ig,ilay-1)
         ENDDO
      ENDDO

!     Calcul of wind stress

      DO ig=1,ngrid
         flux_u(ig) = zfluxv(ig,1)/ptimestep*zu(ig,1)
         flux_v(ig) = zfluxv(ig,1)/ptimestep*zv(ig,1)
      ENDDO


!----------------------------------------------------------------------------
!     6. Implicit inversion of h, not forgetting the coupling with the ground

!     h(t+1) =  h(t) + dt * {(dh/dt)phys}(t)  +  dt * {(dh/dt)difv}(t+1)
!     avec
!     /zh/ = h(t) + dt * {(dh/dt)phys}(t)   (voir paragraphe 2.)
!     et
!     dt * {(dh/dt)difv}(t+1) = dt * {(d/dz)[ Kh (dh/dz) ]}(t+1)
!     donc les entrees sont /zcdh/ pour la condition de raccord au sol
!     et /zkh/ = Kh

!     Using the wind modified by friction for lifting and sublimation
!     ---------------------------------------------------------------
      DO ig=1,ngrid
         zu2      = zu(ig,1)*zu(ig,1)+zv(ig,1)*zv(ig,1)
         zcdv(ig) = zcdv_true(ig)*sqrt(zu2)
         zcdh(ig) = zcdh_true(ig)*sqrt(zu2)
         zkh(ig,1)= zcdh(ig)
         ustar(ig)=sqrt(zcdv_true(ig))*sqrt(zu2)
      ENDDO


!     JL12 calculate the flux coefficients (tables multiplied elements by elements)
!     ---------------------------------------------------------------
      zfluxq(1:ngrid,1:nlay)=zkh(1:ngrid,1:nlay)*zb0(1:ngrid,1:nlay) !JL12 we save zfluxq which doesn't need the Exner factor
      zfluxt(1:ngrid,1:nlay)=zfluxq(1:ngrid,1:nlay)*zExner(1:ngrid,1:nlay)

      DO ig=1,ngrid
         z1(ig)=1./(zmass(ig,nlay)+zfluxt(ig,nlay)*zovExner(ig,nlay))
         zct(ig,nlay)=zmass(ig,nlay)*zt(ig,nlay)*z1(ig)
         zdt(ig,nlay)=zfluxt(ig,nlay)*zovExner(ig,nlay-1)*z1(ig)
      ENDDO

      DO ilay=nlay-1,2,-1
         DO ig=1,ngrid
            z1(ig)=1./(zmass(ig,ilay)+zfluxt(ig,ilay)*zovExner(ig,ilay)+   &
            zfluxt(ig,ilay+1)*(zovExner(ig,ilay)-zdt(ig,ilay+1)*zovExner(ig,ilay+1)))
            zct(ig,ilay)=(zmass(ig,ilay)*zt(ig,ilay)+zfluxt(ig,ilay+1)*zct(ig,ilay+1)*zovExner(ig,ilay+1))*z1(ig)
            zdt(ig,ilay)=zfluxt(ig,ilay)*z1(ig)*zovExner(ig,ilay-1)
         ENDDO
      ENDDO

!JL12 we treat last point afterward because zovExner(ig,ilay-1) does not exist there
      DO ig=1,ngrid
         z1(ig)=1./(zmass(ig,1)+zfluxt(ig,1)*zovExner(ig,1)+  &
             zfluxt(ig,2)*(zovExner(ig,1)-zdt(ig,2)*zovExner(ig,2)))
         zct(ig,1)=(zmass(ig,1)*zt(ig,1)+zfluxt(ig,2)*zct(ig,2)*zovExner(ig,2))*z1(ig)
         zdt(ig,1)=zfluxt(ig,1)*z1(ig)
      ENDDO


!     Calculate (d Planck / dT) at the interface temperature
!     ------------------------------------------------------

      z4st=4.0*sigma*ptimestep
      DO ig=1,ngrid
         zdplanck(ig)=z4st*pemis(ig)*ptsrf(ig)*ptsrf(ig)*ptsrf(ig)
      ENDDO

!     Calculate temperature tendency at the interface (dry case)
!     ----------------------------------------------------------
!     Sum of fluxes at interface at time t + \delta t gives change in T:
!       radiative fluxes
!       turbulent convective (sensible) heat flux
!       flux (if any) from subsurface

      if(.not.water) then

         DO ig=1,ngrid
            z1(ig)=pcapcal(ig)*ptsrf(ig)+cpp*zfluxt(ig,1)*zct(ig,1)*zovExner(ig,1) &
                + pfluxsrf(ig)*ptimestep + zdplanck(ig)*ptsrf(ig) 
            z2(ig) = pcapcal(ig)+zdplanck(ig)+cpp*zfluxt(ig,1)*(1.-zovExner(ig,1)*zdt(ig,1)) 
            ztsrf(ig) = z1(ig) / z2(ig)
            pdtsrf(ig) = (ztsrf(ig) - ptsrf(ig))/ptimestep
            zt(ig,1)   = zct(ig,1) + zdt(ig,1)*ztsrf(ig)
         ENDDO
! JL12 note that the black body radiative flux emitted by the surface has been updated by the implicit scheme 


!     Recalculate temperature to top of atmosphere, starting from ground
!     ------------------------------------------------------------------

         DO ilay=2,nlay
            DO ig=1,ngrid
               zt(ig,ilay)=zct(ig,ilay)+zdt(ig,ilay)*zt(ig,ilay-1)
            ENDDO
         ENDDO

      endif                     ! not water

!-----------------------------------------------------------------------
!     TRACERS (no vapour)
!     -------

      if(tracer) then

!     Calculate vertical flux from the bottom to the first layer (dust)
!     -----------------------------------------------------------------
         do ig=1,ngrid
            rho(ig) = zb0(ig,1) /ptimestep
         end do

         pdqsdif(:,:)=0.

!     Implicit inversion of q
!     -----------------------
         do iq=1,nq 

            if (iq.ne.ivap) then

               DO ig=1,ngrid
                  z1(ig)=1./(zmass(ig,nlay)+zfluxq(ig,nlay))
                  zcq(ig,nlay)=zmass(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
                  zdq(ig,nlay)=zfluxq(ig,nlay)*z1(ig)
               ENDDO 
            
               DO ilay=nlay-1,2,-1
                  DO ig=1,ngrid
                     z1(ig)=1./(zmass(ig,ilay)+zfluxq(ig,ilay)+zfluxq(ig,ilay+1)*(1.-zdq(ig,ilay+1)))
                     zcq(ig,ilay)=(zmass(ig,ilay)*zq(ig,ilay,iq)+zfluxq(ig,ilay+1)*zcq(ig,ilay+1))*z1(ig)
                     zdq(ig,ilay)=zfluxq(ig,ilay)*z1(ig)
                  ENDDO
               ENDDO

               if ((water).and.(iq.eq.iliq)) then
                  ! special case for condensed water tracer: do not include
                  ! h2o ice tracer from surface (which is set when handling
                  ! h2o vapour case (see further down).
                  ! zb(ig,1)=0 if iq ne ivap
                  DO ig=1,ngrid
                     z1(ig)=1./(zmass(ig,1)+zfluxq(ig,2)*(1.-zdq(ig,2)))
                     zcq(ig,1)=(zmass(ig,1)*zq(ig,1,iq)+zfluxq(ig,2)*zcq(ig,2))*z1(ig)
                  ENDDO
               else             ! general case
                  do ig=1,ngrid
                     z1(ig)=1./(zmass(ig,1)+zfluxq(ig,2)*(1.-zdq(ig,2)))
                     zcq(ig,1)=(zmass(ig,1)*zq(ig,1,iq)+zfluxq(ig,2)*zcq(ig,2)+(-pdqsdif(ig,iq))*ptimestep)*z1(ig)
                          ! tracer flux from surface
                          ! currently pdqsdif always zero here,
                          ! so last line is superfluous
                  enddo
               endif            ! of if (water.and.(iq.eq.igcm_h2o_ice))


!     Starting upward calculations for simple tracer mixing (e.g., dust)
               do ig=1,ngrid
                  zq(ig,1,iq)=zcq(ig,1)
               end do

               do ilay=2,nlay
                  do ig=1,ngrid
                     zq(ig,ilay,iq)=zcq(ig,ilay)+zdq(ig,ilay)*zq(ig,ilay-1,iq)
                  end do
               end do

            endif               ! if (iq.ne.ivap)

!     Calculate temperature tendency including latent heat term
!     and assuming an infinite source of water on the ground
!     ------------------------------------------------------------------

            if (water.and.(iq.eq.ivap)) then 
            
               ! compute evaporation efficiency
               do ig=1,ngrid
                  if(nint(rnat(ig)).eq.1)then
                     dryness(ig)=pqsurf(ig,iliq_surf)+pqsurf(ig,iice_surf) 
                     dryness(ig)=MIN(1.,2*dryness(ig)/mx_eau_sol)
                     dryness(ig)=MAX(0.,dryness(ig))
                  endif
               enddo

               do ig=1,ngrid
                ! Calculate the value of qsat at the surface (water)
                call Psat_water(ptsrf(ig),pplev(ig,1),psat(ig),qsat(ig))
                call Lcpdqsat_water(ptsrf(ig),pplev(ig,1),psat(ig),qsat(ig),dqsat(ig),psat_temp)
                dqsat(ig)=dqsat(ig)*RCPD/RLVTT
	       enddo

! coefficients for q

               do ig=1,ngrid
                  z1(ig)=1./(zmass(ig,nlay)+zfluxq(ig,nlay))
                  zcq(ig,nlay)=zmass(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
                  zdq(ig,nlay)=zfluxq(ig,nlay)*z1(ig)
               enddo 
          
               do ilay=nlay-1,2,-1
                  do ig=1,ngrid
                     z1(ig)=1./(zmass(ig,ilay)+zfluxq(ig,ilay)+zfluxq(ig,ilay+1)*(1.-zdq(ig,ilay+1)))
                     zcq(ig,ilay)=(zmass(ig,ilay)*zq(ig,ilay,iq)+zfluxq(ig,ilay+1)*zcq(ig,ilay+1))*z1(ig)
                     zdq(ig,ilay)=zfluxq(ig,ilay)*z1(ig)
                  enddo
               enddo

               do ig=1,ngrid
                  z1(ig)=1./(zmass(ig,1)+zfluxq(ig,1)*dryness(ig)+zfluxq(ig,2)*(1.-zdq(ig,2)))
                  zcq(ig,1)=(zmass(ig,1)*zq(ig,1,iq)+zfluxq(ig,2)*zcq(ig,2))*z1(ig)
                  zdq(ig,1)=dryness(ig)*zfluxq(ig,1)*z1(ig)
               enddo

              do ig=1,ngrid
!calculation of surface temperature
                  zdq0(ig) = dqsat(ig)
                  zcq0(ig) = qsat(ig)-dqsat(ig)*ptsrf(ig)

                  z1(ig) = pcapcal(ig)*ptsrf(ig) +cpp*zfluxt(ig,1)*zct(ig,1)*zovExner(ig,1)   &
		      + zdplanck(ig)*ptsrf(ig) + pfluxsrf(ig)*ptimestep                       &
                      + zfluxq(ig,1)*dryness(ig)*RLVTT*((zdq(ig,1)-1.0)*zcq0(ig)+zcq(ig,1))

                  z2(ig) = pcapcal(ig) + cpp*zfluxt(ig,1)*(1.-zovExner(ig,1)*zdt(ig,1))       &
                      + zdplanck(ig)+zfluxq(ig,1)*dryness(ig)*RLVTT*zdq0(ig)*(1.0-zdq(ig,1))

                  ztsrf(ig) = z1(ig) / z2(ig)

! calculation of qs and q1
                  zq0(ig)     = zcq0(ig)+zdq0(ig)*ztsrf(ig)
                  zq(ig,1,iq) = zcq(ig,1)+zdq(ig,1)*zq0(ig)

! calculation of evaporation              
                  dqsdif_total(ig)=zfluxq(ig,1)*dryness(ig)*(zq(ig,1,ivap)-zq0(ig))

!     --------------------------------------------------------
!     Now check if we've taken too much water from the surface
!     This can only occur on the continent 
!     If we do, we recompute Tsurf, T1 and q1 accordingly
                  if((-dqsdif_total(ig).gt.(pqsurf(ig,iice_surf)+pqsurf(ig,iliq_surf))).and.rnat(ig).eq.1)then
                      !water flux * ptimestep
	              dqsdif_total(ig)=-(pqsurf(ig,iice_surf)+pqsurf(ig,iliq_surf))

                      !recompute surface temperature  
                      z1(ig) = pcapcal(ig)*ptsrf(ig) +cpp*zfluxq(ig,1)*zct(ig,1)*zovExner(ig,1)   &
		        + zdplanck(ig)*ptsrf(ig) + pfluxsrf(ig)*ptimestep                       &
                        + RLVTT*dqsdif_total(ig)
                      z2(ig) = pcapcal(ig) + cpp*zfluxq(ig,1)*(1.-zovExner(ig,1)*zdt(ig,1))       &
                        + zdplanck(ig)
                      ztsrf(ig) = z1(ig) / z2(ig)

                      !recompute q1 with new water flux from surface  
                      zq(ig,1,iq) = (zmass(ig,1)*(pq(ig,1,iq)+ptimestep*pdqfi(ig,1,iq))  &
		                            +zfluxq(ig,2)*zcq(ig,2)-dqsdif_total(ig))     &
                                 / (zmass(ig,1)+(1.-zdq(ig,2))*zfluxq(ig,2))		     
                  end if
		  
! calculation surface T tendency  and T(1)           
		  pdtsrf(ig) = (ztsrf(ig) - ptsrf(ig))/ptimestep
                  zt(ig,1)   = zct(ig,1) + zdt(ig,1)*ztsrf(ig)		      
               enddo


! recalculate temperature and q(vap) to top of atmosphere, starting from ground
               do ilay=2,nlay
                  do ig=1,ngrid
                     zq(ig,ilay,iq)=zcq(ig,ilay)+zdq(ig,ilay)*zq(ig,ilay-1,iq)
                     zt(ig,ilay)=zct(ig,ilay)+zdt(ig,ilay)*zt(ig,ilay-1)
                  end do
               end do


               do ig=1,ngrid
!     --------------------------------------------------------------------------
!     On the ocean, if T > 0 C then the vapour tendency must replace the ice one
!     The surface vapour tracer is actually liquid. To make things difficult.

                  if (nint(rnat(ig)).eq.0) then ! unfrozen ocean
                     
                     pdqsdif(ig,iliq_surf)=dqsdif_total(ig)/ptimestep
                     pdqsdif(ig,iice_surf)=0.0

                  elseif (nint(rnat(ig)).eq.1) then ! (continent)
!     If water is evaporating / subliming, we take it from ice before liquid
!     -- is this valid??
                     if(dqsdif_total(ig).lt.0)then
                        if (-dqsdif_total(ig).gt.pqsurf(ig,iice_surf))then
                           pdqsdif(ig,iice_surf) = -pqsurf(ig,iice_surf)/ptimestep ! removes all the ice!
                           pdqsdif(ig,iliq_surf) = dqsdif_total(ig)/ptimestep- pdqsdif(ig,iice_surf) ! take the remainder from the liquid instead
			else               
                           pdqsdif(ig,iice_surf)=dqsdif_total(ig)/ptimestep
			   pdqsdif(ig,iliq_surf)=0.
			end if
	             else !dqsdif_total(ig).ge.0
                        !If water vapour is condensing, we must decide whether it forms ice or liquid.
                        if(ztsrf(ig).gt.T_h2O_ice_liq)then
                           pdqsdif(ig,iice_surf)=0.0
                           pdqsdif(ig,iliq_surf)=dqsdif_total(ig)/ptimestep
                        else
                           pdqsdif(ig,iice_surf)=dqsdif_total(ig)/ptimestep
                           pdqsdif(ig,iliq_surf)=0.0
                        endif               
                     endif

                  elseif (nint(rnat(ig)).eq.2) then ! (continental glaciers)
		     pdqsdif(ig,iliq_surf)=0.0
                     pdqsdif(ig,iice_surf)=dqsdif_total(ig)/ptimestep

		  endif !rnat
               end do            ! of DO ig=1,ngrid

           endif                ! if (water et iq=ivap)
        end do                  ! of do iq=1,nq

        if (water) then  ! special case where we recompute water mixing without any evaporation.
	                 !    The difference with the first calculation then tells us where evaporated water has gone

            DO ig=1,ngrid
               z1(ig)=1./(zmass(ig,nlay)+zfluxq(ig,nlay))
               zcq(ig,nlay)=zmass(ig,nlay)*zqnoevap(ig,nlay)*z1(ig)
               zdq(ig,nlay)=zfluxq(ig,nlay)*z1(ig)
            ENDDO 
            
            DO ilay=nlay-1,2,-1
               DO ig=1,ngrid
                  z1(ig)=1./(zmass(ig,ilay)+zfluxq(ig,ilay)+zfluxq(ig,ilay+1)*(1.-zdq(ig,ilay+1)))
                  zcq(ig,ilay)=(zmass(ig,ilay)*zqnoevap(ig,ilay)+zfluxq(ig,ilay+1)*zcq(ig,ilay+1))*z1(ig)
                  zdq(ig,ilay)=zfluxq(ig,ilay)*z1(ig)
               ENDDO
            ENDDO

            do ig=1,ngrid
               z1(ig)=1./(zmass(ig,1)+zfluxq(ig,2)*(1.-zdq(ig,2)))
               zcq(ig,1)=(zmass(ig,1)*zqnoevap(ig,1)+zfluxq(ig,2)*zcq(ig,2))*z1(ig)
            enddo

!     Starting upward calculations for simple tracer mixing (e.g., dust)
            do ig=1,ngrid
               zqnoevap(ig,1)=zcq(ig,1)
            end do

            do ilay=2,nlay
               do ig=1,ngrid
                  zqnoevap(ig,ilay)=zcq(ig,ilay)+zdq(ig,ilay)*zqnoevap(ig,ilay-1)
               end do
            end do

         endif               ! if water
	
	
      endif                     ! tracer


!-----------------------------------------------------------------------
!     8. Final calculation of the vertical diffusion tendencies
!     -----------------------------------------------------------------

      do ilev = 1, nlay
         do ig=1,ngrid
            pdudif(ig,ilev)=(zu(ig,ilev)-(pu(ig,ilev)))/ptimestep
            pdvdif(ig,ilev)=(zv(ig,ilev)-(pv(ig,ilev)))/ptimestep
            pdtdif(ig,ilev)=( zt(ig,ilev)- pt(ig,ilev))/ptimestep-pdtfi(ig,ilev)
         enddo
      enddo
      
      DO ig=1,ngrid ! computing sensible heat flux (atm => surface)
	 sensibFlux(ig)=cpp*zfluxt(ig,1)/ptimestep*(zt(ig,1)*zovExner(ig,1)-ztsrf(ig))
      ENDDO

      if (tracer) then
         do iq = 1, nq
            do ilev = 1, nlay
               do ig=1,ngrid
                  pdqdif(ig,ilev,iq)=(zq(ig,ilev,iq)-(pq(ig,ilev,iq)+pdqfi(ig,ilev,iq)*ptimestep))/ptimestep
               enddo
            enddo
         enddo
	 if (water) then
            do ilev = 1, nlay
               do ig=1,ngrid
                  pdqevap(ig,ilev)=(zq(ig,ilev,ivap)-zqnoevap(ig,ilev))/ptimestep
               enddo
            enddo
            do ig=1,ngrid
	       zdmassevap(ig)=SUM(pdqevap(ig,:)*zmass(ig,:))*ptimestep
	    end do	    
	 endif
      endif 

      if(water)then
#ifndef MESOSCALE
         call writediagfi(ngrid,'beta','Dryness coefficient',' ',2,dryness)
#endif
         if (tracer) then
#ifndef MESOSCALE
            call writediagfi(ngrid,'evap_surf_flux','surface latent heat flux','W.m-2',2,RLVTT*dqsdif_total/ptimestep)
            call writediagfi(ngrid,'fluxsurf_rad','total IR and VIS surface flux','W.m-2',2,pfluxsrf)
            call writediagfi(ngrid,'dqevap','evaporated water vapor specific concentration','s-1',3,pdqevap)
#else
            comm_LATENT_HF(:)=0.0
            comm_LATENT_HF(1:ngrid)=RLVTT*dqsdif_total(1:ngrid)/ptimestep
#endif
	 endif
      endif

!      if(lastcall)then
!        if(ngrid.eq.1)then
!           print*,'Saving k.out...'
!           OPEN(12,file='k.out',form='formatted')
!           DO ilay=1,nlay
!              write(12,*) zkh(1,ilay), pplay(1,ilay)
!           ENDDO
!           CLOSE(12)
!         endif
!      endif

      end
