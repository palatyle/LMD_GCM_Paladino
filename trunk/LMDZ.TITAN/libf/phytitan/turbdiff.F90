      subroutine turbdiff(ngrid,nlay,nq,               &
          ptimestep,pcapcal,lecrit,                    &   
          pplay,pplev,pzlay,pzlev,pz0,                 &
          pu,pv,pt,ppopsk,pq,ptsrf,pemis,pqsurf,       &
          pdtfi,pdqfi,pfluxsrf,            &
          Pdudif,pdvdif,pdtdif,pdtsrf,sensibFlux,pq2,  &
          pdqdif,pdqsdif,flux_u,flux_v,lastcall)

      use radcommon_h, only : sigma, gzlat
      use comcstfi_mod, only: rcp, g, r, cpp
      use callkeys_mod, only: tracer,nosurf

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
      REAL zdmassevap(ngrid)
      REAL rho(ngrid)         ! near-surface air density
      REAL kmixmin


      real, parameter :: karman=0.4
      real cd0, roughratio

      real dqsdif_total(ngrid) 
      real zq0(ngrid) 


!     Coherence test
!     --------------

      IF (firstcall) THEN      	  

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
            zmass(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/gzlat(ig,ilay)
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
               
               do ig=1,ngrid
                  z1(ig)=1./(zmass(ig,1)+zfluxq(ig,2)*(1.-zdq(ig,2)))
                  zcq(ig,1)=(zmass(ig,1)*zq(ig,1,iq)+zfluxq(ig,2)*zcq(ig,2)+(-pdqsdif(ig,iq))*ptimestep)*z1(ig)
                  ! tracer flux from surface
                  ! currently pdqsdif always zero here,
                  ! so last line is superfluous
               enddo

!     Starting upward calculations for simple tracer mixing (e.g., dust)
               do ig=1,ngrid
                  zq(ig,1,iq)=zcq(ig,1)
               end do

               do ilay=2,nlay
                  do ig=1,ngrid
                     zq(ig,ilay,iq)=zcq(ig,ilay)+zdq(ig,ilay)*zq(ig,ilay-1,iq)
                  end do
               end do
        end do                  ! of do iq=1,nq
                       
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
