










      subroutine vdifc(ngrid,nlay,nq,rnat,ppopsk,         
     &     ptimestep,pcapcal,lecrit,                        
     &     pplay,pplev,pzlay,pzlev,pz0,
     &     pu,pv,ph,pq,ptsrf,pemis,pqsurf,
     &     pdhfi,pdqfi,pfluxsrf,
     &     pdudif,pdvdif,pdhdif,pdtsrf,sensibFlux,pq2,
     &     pdqdif,pdqsdif,lastcall)

      use watercommon_h, only : RLVTT, T_h2O_ice_liq, RCPD, mx_eau_sol 
     &    ,Psat_water, Lcpdqsat_water
      use radcommon_h, only : sigma
      USE surfdat_h
      USE tracer_h
      use comcstfi_mod, only: g, r, cpp, rcp
      use callkeys_mod, only: water,tracer,nosurf

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
!     
!==================================================================

!-----------------------------------------------------------------------
!     declarations
!     ------------


!     arguments
!     ---------
      INTEGER ngrid,nlay
      REAL ptimestep
      REAL pplay(ngrid,nlay),pplev(ngrid,nlay+1)
      REAL pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
      REAL pu(ngrid,nlay),pv(ngrid,nlay),ph(ngrid,nlay)
      REAL ptsrf(ngrid),pemis(ngrid)
      REAL pdhfi(ngrid,nlay)
      REAL pfluxsrf(ngrid)
      REAL pdudif(ngrid,nlay),pdvdif(ngrid,nlay),pdhdif(ngrid,nlay)
      REAL pdtsrf(ngrid),sensibFlux(ngrid),pcapcal(ngrid)
      REAL pq2(ngrid,nlay+1)
      
      real rnat(ngrid)      

!     Arguments added for condensation
      REAL ppopsk(ngrid,nlay)
      logical lecrit
      REAL pz0

!     Tracers
!     --------
      integer nq 
      real pqsurf(ngrid,nq)
      real pq(ngrid,nlay,nq), pdqfi(ngrid,nlay,nq) 
      real pdqdif(ngrid,nlay,nq) 
      real pdqsdif(ngrid,nq) 
      
!     local
!     -----
      integer ilev,ig,ilay,nlev

      REAL z4st,zdplanck(ngrid)
      REAL zkv(ngrid,nlay+1),zkh(ngrid,nlay+1)
      REAL zcdv(ngrid),zcdh(ngrid)
      REAL zcdv_true(ngrid),zcdh_true(ngrid)
      REAL zu(ngrid,nlay),zv(ngrid,nlay)
      REAL zh(ngrid,nlay)
      REAL ztsrf2(ngrid)
      REAL z1(ngrid),z2(ngrid)
      REAL za(ngrid,nlay),zb(ngrid,nlay)
      REAL zb0(ngrid,nlay)
      REAL zc(ngrid,nlay),zd(ngrid,nlay)
      REAL zcst1
      REAL zu2!, a
      REAL zcq(ngrid,nlay),zdq(ngrid,nlay)
      REAL evap(ngrid)
      REAL zcq0(ngrid),zdq0(ngrid)
      REAL zx_alf1(ngrid),zx_alf2(ngrid)

      LOGICAL firstcall
      SAVE firstcall
!$OMP THREADPRIVATE(firstcall)
      
      LOGICAL lastcall

!     variables added for CO2 condensation
!     ------------------------------------
      REAL hh                   !, zhcond(ngrid,nlay)
!     REAL latcond,tcond1mb
!     REAL acond,bcond
!     SAVE acond,bcond
!!$OMP THREADPRIVATE(acond,bcond)
!     DATA latcond,tcond1mb/5.9e5,136.27/

!     Tracers
!     -------
      INTEGER iq
      REAL zq(ngrid,nlay,nq)
      REAL zq1temp(ngrid)
      REAL rho(ngrid)         ! near-surface air density
      DATA firstcall/.true./
      REAL kmixmin

!     Variables added for implicit latent heat inclusion
!     --------------------------------------------------
      real dqsat(ngrid),psat_temp,qsat(ngrid),psat(ngrid)

      integer ivap, iice ! also make liq for clarity on surface...
      save ivap, iice
!$OMP THREADPRIVATE(ivap,iice)

      real, parameter :: karman=0.4
      real cd0, roughratio

      logical forceWC
      real masse, Wtot, Wdiff

      real dqsdif_total(ngrid) 
      real zq0(ngrid) 

      forceWC=.true.
!      forceWC=.false.


!     Coherence test
!     --------------

      IF (firstcall) THEN
!     To compute: Tcond= 1./(bcond-acond*log(.0095*p)) (p in pascal)
!     bcond=1./tcond1mb
!     acond=r/latcond
!     PRINT*,'In vdifc: Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
!     PRINT*,'          acond,bcond',acond,bcond

         if(water)then
!                iliq=igcm_h2o_vap
                ivap=igcm_h2o_vap
                iice=igcm_h2o_ice ! simply to make the code legible               
                                  ! to be generalised later
         endif

         firstcall=.false.
      ENDIF

!-----------------------------------------------------------------------
!     1. Initialisation
!     -----------------

      nlev=nlay+1

!     Calculate rho*dz and dt*rho/dz=dt*rho**2 g/dp
!     with rho=p/RT=p/ (R Theta) (p/ps)**kappa
!     ---------------------------------------------

      DO ilay=1,nlay
         DO ig=1,ngrid
            za(ig,ilay)=(pplev(ig,ilay)-pplev(ig,ilay+1))/g
         ENDDO
      ENDDO

      zcst1=4.*g*ptimestep/(R*R)
      DO ilev=2,nlev-1
         DO ig=1,ngrid
            zb0(ig,ilev)=pplev(ig,ilev)*
     s           (pplev(ig,1)/pplev(ig,ilev))**rcp /
     s           (ph(ig,ilev-1)+ph(ig,ilev))
            zb0(ig,ilev)=zcst1*zb0(ig,ilev)*zb0(ig,ilev)/
     s           (pplay(ig,ilev-1)-pplay(ig,ilev))
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
            zh(ig,ilev)=ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep
         ENDDO
      ENDDO
      if(tracer) then
         DO iq =1, nq
            DO ilev=1,nlay
               DO ig=1,ngrid
                  zq(ig,ilev,iq)=pq(ig,ilev,iq) + 
     &                 pdqfi(ig,ilev,iq)*ptimestep
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
         roughratio = 1.E+0 + pzlay(ig,1)/pz0
         cd0 = karman/log(roughratio)
         cd0 = cd0*cd0
         zcdv_true(ig) = cd0
         zcdh_true(ig) = cd0
         if (nosurf) then
             zcdv_true(ig) = 0.   !! disable sensible momentum flux
             zcdh_true(ig) = 0.   !! disable sensible heat flux
         endif
      ENDDO

      DO ig=1,ngrid
         zu2=pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1)
         zcdv(ig)=zcdv_true(ig)*sqrt(zu2)
         zcdh(ig)=zcdh_true(ig)*sqrt(zu2)
      ENDDO

!     Turbulent diffusion coefficients in the boundary layer
!     ------------------------------------------------------ 

      call vdif_kc(ngrid,nlay,ptimestep,g,pzlev,pzlay
     &     ,pu,pv,ph,zcdv_true
     &     ,pq2,zkv,zkh)

!     Adding eddy mixing to mimic 3D general circulation in 1D
!     R. Wordsworth & F. Forget (2010)
      if ((ngrid.eq.1)) then
         kmixmin = 1.0e-2       ! minimum eddy mix coeff in 1D
         do ilev=1,nlay
            do ig=1,ngrid
               !zkh(ig,ilev) = 1.0
               zkh(ig,ilev) = max(kmixmin,zkh(ig,ilev))
               zkv(ig,ilev) = max(kmixmin,zkv(ig,ilev))
            end do
         end do
      end if

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
      
      CALL multipl((nlay-1)*ngrid,zkv(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdv,zb0,zb)

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zu(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $           zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zu(ig,ilay)+
     $           zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zu(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zu(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zu(ig,ilay-1)
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
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zv(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,1,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $           zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zv(ig,ilay)+
     $           zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         zv(ig,1)=zc(ig,1)
      ENDDO
      DO ilay=2,nlay
         DO ig=1,ngrid
            zv(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zv(ig,ilay-1)
         ENDDO
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
         ENDDO

      CALL multipl((nlay-1)*ngrid,zkh(1,2),zb0(1,2),zb(1,2))
      CALL multipl(ngrid,zcdh,zb0,zb)

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
         zc(ig,nlay)=za(ig,nlay)*zh(ig,nlay)*z1(ig)
         zd(ig,nlay)=zb(ig,nlay)*z1(ig)
      ENDDO

      DO ilay=nlay-1,2,-1
         DO ig=1,ngrid
            z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     &           zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
            zc(ig,ilay)=(za(ig,ilay)*zh(ig,ilay)+
     &           zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
            zd(ig,ilay)=zb(ig,ilay)*z1(ig)
         ENDDO
      ENDDO

      DO ig=1,ngrid
         z1(ig)=1./(za(ig,1)+zb(ig,1)+
     &        zb(ig,2)*(1.-zd(ig,2)))
         zc(ig,1)=(za(ig,1)*zh(ig,1)+
     &        zb(ig,2)*zc(ig,2))*z1(ig)
         zd(ig,1)=zb(ig,1)*z1(ig)
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

            z1(ig) = pcapcal(ig)*ptsrf(ig) + cpp*zb(ig,1)*zc(ig,1)
     &           + zdplanck(ig)*ptsrf(ig) + pfluxsrf(ig)*ptimestep
            z2(ig) = pcapcal(ig) + cpp*zb(ig,1)*(1.-zd(ig,1)) 
     &           +zdplanck(ig)
            ztsrf2(ig) = z1(ig) / z2(ig)
            pdtsrf(ig) = (ztsrf2(ig) - ptsrf(ig))/ptimestep
            zh(ig,1)   = zc(ig,1) + zd(ig,1)*ztsrf2(ig)
         ENDDO

!     Recalculate temperature to top of atmosphere, starting from ground
!     ------------------------------------------------------------------

         DO ilay=2,nlay
            DO ig=1,ngrid
               hh = zh(ig,ilay-1)
               zh(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*hh
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

         pdqsdif(1:ngrid,1:nq)=0

!     Implicit inversion of q
!     -----------------------
         do iq=1,nq 

            if (iq.ne.igcm_h2o_vap) then

               DO ig=1,ngrid
                  z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
                  zcq(ig,nlay)=za(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
                  zdq(ig,nlay)=zb(ig,nlay)*z1(ig)
               ENDDO 
            
               DO ilay=nlay-1,2,-1
                  DO ig=1,ngrid
                     z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     &                    zb(ig,ilay+1)*(1.-zdq(ig,ilay+1)))
                     zcq(ig,ilay)=(za(ig,ilay)*zq(ig,ilay,iq)+
     &                    zb(ig,ilay+1)*zcq(ig,ilay+1))*z1(ig)
                     zdq(ig,ilay)=zb(ig,ilay)*z1(ig)
                  ENDDO
               ENDDO

               if ((water).and.(iq.eq.iice)) then
                  ! special case for water ice tracer: do not include
                  ! h2o ice tracer from surface (which is set when handling
                  ! h2o vapour case (see further down).
                  ! zb(ig,1)=0 if iq ne ivap
                  DO ig=1,ngrid
                     z1(ig)=1./(za(ig,1)+
     &                    zb(ig,2)*(1.-zdq(ig,2)))
                     zcq(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     &                    zb(ig,2)*zcq(ig,2))*z1(ig)
                  ENDDO
               else             ! general case
                  DO ig=1,ngrid
                     z1(ig)=1./(za(ig,1)+
     &                    zb(ig,2)*(1.-zdq(ig,2)))
                     zcq(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     &                    zb(ig,2)*zcq(ig,2)
     &                    +(-pdqsdif(ig,iq))*ptimestep)*z1(ig)
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
                     zq(ig,ilay,iq)=zcq(ig,ilay)+
     $                    zdq(ig,ilay)*zq(ig,ilay-1,iq)
                  end do
               end do

            endif               ! if (iq.ne.igcm_h2o_vap)

!     Calculate temperature tendency including latent heat term
!     and assuming an infinite source of water on the ground
!     ------------------------------------------------------------------

            if (water.and.(iq.eq.igcm_h2o_vap)) then 
            
               ! compute evaporation efficiency
               do ig = 1, ngrid
                  if(nint(rnat(ig)).eq.1)then
                     dryness(ig)=pqsurf(ig,ivap)+pqsurf(ig,iice) 
                     dryness(ig)=MIN(1.,2*dryness(ig)/mx_eau_sol)
                     dryness(ig)=MAX(0.,dryness(ig))
                  endif
               enddo

               do ig=1,ngrid
                ! Calculate the value of qsat at the surface (water)
                call Psat_water(ptsrf(ig),pplev(ig,1),psat(ig),qsat(ig))
                call Lcpdqsat_water(ptsrf(ig),pplev(ig,1),psat(ig),
     &             qsat(ig),dqsat(ig),psat_temp)
                dqsat(ig)=dqsat(ig)*RCPD/RLVTT
               enddo

! coefficients for q

               do ig=1,ngrid
                  z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
                  zcq(ig,nlay)=za(ig,nlay)*zq(ig,nlay,iq)*z1(ig)
                  zdq(ig,nlay)=zb(ig,nlay)*z1(ig)
               enddo 
            
               do ilay=nlay-1,2,-1
                  do ig=1,ngrid
                     z1(ig)=1./(za(ig,ilay)+zb(ig,ilay)+
     $                    zb(ig,ilay+1)*(1.-zdq(ig,ilay+1)))
                     zcq(ig,ilay)=(za(ig,ilay)*zq(ig,ilay,iq)+
     $                    zb(ig,ilay+1)*zcq(ig,ilay+1))*z1(ig)
                     zdq(ig,ilay)=zb(ig,ilay)*z1(ig)
                  enddo
               enddo

               do ig=1,ngrid
                  z1(ig)=1./(za(ig,1)+zb(ig,1)*dryness(ig)+
     $                 zb(ig,2)*(1.-zdq(ig,2)))
                  zcq(ig,1)=(za(ig,1)*zq(ig,1,iq)+
     $                 zb(ig,2)*zcq(ig,2))*z1(ig)
                  zdq(ig,1)=dryness(ig)*zb(ig,1)*z1(ig)
               enddo

! calculation of h0 and h1
               do ig=1,ngrid
                  zdq0(ig) = dqsat(ig)
                  zcq0(ig) = qsat(ig)-dqsat(ig)*ptsrf(ig)

                  z1(ig) = pcapcal(ig)*ptsrf(ig) +cpp*zb(ig,1)*zc(ig,1)
     &                 + zdplanck(ig)*ptsrf(ig) + pfluxsrf(ig)*ptimestep 
     &                 +  zb(ig,1)*dryness(ig)*RLVTT*
     &                 ((zdq(ig,1)-1.0)*zcq0(ig)+zcq(ig,1))

                  z2(ig) = pcapcal(ig) + cpp*zb(ig,1)*(1.-zd(ig,1))
     &                 +zdplanck(ig)
     &                 +zb(ig,1)*dryness(ig)*RLVTT*zdq0(ig)*
     &                 (1.0-zdq(ig,1))

                  ztsrf2(ig) = z1(ig) / z2(ig)
                  pdtsrf(ig) = (ztsrf2(ig) - ptsrf(ig))/ptimestep
                  zh(ig,1)   = zc(ig,1) + zd(ig,1)*ztsrf2(ig)
               enddo

! calculation of qs and q1
               do ig=1,ngrid
                  zq0(ig)     = zcq0(ig)+zdq0(ig)*ztsrf2(ig)
                  zq(ig,1,iq) = zcq(ig,1)+zdq(ig,1)*zq0(ig)
               enddo

! calculation of evaporation              
               do ig=1,ngrid 
                  evap(ig)= zb(ig,1)*dryness(ig)*(zq(ig,1,ivap)-zq0(ig))
                  dqsdif_total(ig)=evap(ig)
               enddo

! recalculate temperature and q(vap) to top of atmosphere, starting from ground
               do ilay=2,nlay
                  do ig=1,ngrid
                     zq(ig,ilay,iq)=zcq(ig,ilay)
     &                    +zdq(ig,ilay)*zq(ig,ilay-1,iq)
                     zh(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zh(ig,ilay-1)
                  end do
               end do

               do ig=1,ngrid

!     --------------------------------------------------------------------------
!     On the ocean, if T > 0 C then the vapour tendency must replace the ice one
!     The surface vapour tracer is actually liquid. To make things difficult.

                  if (nint(rnat(ig)).eq.0) then ! unfrozen ocean
                     
                     pdqsdif(ig,ivap)=dqsdif_total(ig)/ptimestep
                     pdqsdif(ig,iice)=0.0


                  elseif (nint(rnat(ig)).eq.1) then ! (continent)

!     --------------------------------------------------------
!     Now check if we've taken too much water from the surface
!     This can only occur on the continent 

!     If water is evaporating / subliming, we take it from ice before liquid
!     -- is this valid??
                     if(dqsdif_total(ig).lt.0)then
                        pdqsdif(ig,iice)=dqsdif_total(ig)/ptimestep
                        pdqsdif(ig,iice)=max(-pqsurf(ig,iice)/ptimestep
     &                       ,pdqsdif(ig,iice))
                     endif
                     ! sublimation only greater than qsurf(ice)
                     ! ----------------------------------------
                     ! we just convert some liquid to vapour too
                     ! if latent heats are the same, no big deal
                     if (-dqsdif_total(ig).gt.pqsurf(ig,iice))then            
                       pdqsdif(ig,iice) = -pqsurf(ig,iice)/ptimestep ! removes all the ice!
                       pdqsdif(ig,ivap) = dqsdif_total(ig)/ptimestep
     &                       - pdqsdif(ig,iice) ! take the remainder from the liquid instead
                       pdqsdif(ig,ivap) = max(-pqsurf(ig,ivap)/ptimestep
     &                       ,pdqsdif(ig,ivap))
                    endif

                 endif          ! if (rnat.ne.1)

!     If water vapour is condensing, we must decide whether it forms ice or liquid.
                 if(dqsdif_total(ig).gt.0)then ! a bug was here!
                    if(ztsrf2(ig).gt.T_h2O_ice_liq)then
                       pdqsdif(ig,iice)=0.0
                       pdqsdif(ig,ivap)=dqsdif_total(ig)/ptimestep
                    else
                       pdqsdif(ig,iice)=dqsdif_total(ig)/ptimestep
                       pdqsdif(ig,ivap)=0.0
                    endif
                 endif

              end do            ! of DO ig=1,ngrid
           endif                ! if (water et iq=ivap)
        end do                  ! of do iq=1,nq
      endif                     ! traceur


!-----------------------------------------------------------------------
!     8. Final calculation of the vertical diffusion tendencies
!     -----------------------------------------------------------------

      do ilev = 1, nlay
         do ig=1,ngrid
            pdudif(ig,ilev)=(zu(ig,ilev)-
     &           (pu(ig,ilev)))/ptimestep
            pdvdif(ig,ilev)=(zv(ig,ilev)-
     &           (pv(ig,ilev)))/ptimestep
            hh = ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep 

            pdhdif(ig,ilev)=( zh(ig,ilev)- hh )/ptimestep
         enddo
      enddo

      DO ig=1,ngrid  ! computing sensible heat flux (atm => surface)
	 sensibFlux(ig)=cpp*zb(ig,1)/ptimestep*(zh(ig,1)-ztsrf2(ig))
      ENDDO      

      if (tracer) then
         do iq = 1, nq
            do ilev = 1, nlay
               do ig=1,ngrid
                  pdqdif(ig,ilev,iq)=(zq(ig,ilev,iq)-
     &           (pq(ig,ilev,iq)+pdqfi(ig,ilev,iq)*ptimestep))/
     &           ptimestep
               enddo
            enddo
         enddo

         if(water.and.forceWC)then ! force water conservation in model
                                   ! we calculate the difference and add it to the ground
                                   ! this is ugly and should be improved in the future
            do ig=1,ngrid
               Wtot=0.0
               do ilay=1,nlay
                  masse = (pplev(ig,ilay) - pplev(ig,ilay+1))/g
!                  Wtot=Wtot+masse*(zq(ig,ilay,iice)-
!     &                 (pq(ig,ilay,iice)+pdqfi(ig,ilay,iice)*ptimestep))
                  Wtot=Wtot+masse*(zq(ig,ilay,ivap)-
     &                 (pq(ig,ilay,ivap)+pdqfi(ig,ilay,ivap)*ptimestep))
               enddo
               Wdiff=Wtot/ptimestep+pdqsdif(ig,ivap)+pdqsdif(ig,iice)

               if(ztsrf2(ig).gt.T_h2O_ice_liq)then
                  pdqsdif(ig,ivap)=pdqsdif(ig,ivap)-Wdiff
               else
                  pdqsdif(ig,iice)=pdqsdif(ig,iice)-Wdiff
               endif
            enddo

         endif

      endif 

      if(water)then
      call writediagfi(ngrid,'beta','Dryness coefficient',' ',2,dryness)
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


      return
      end
