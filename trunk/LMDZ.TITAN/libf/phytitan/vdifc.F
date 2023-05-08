      subroutine vdifc(ngrid,nlay,nq,ppopsk,         
     &     ptimestep,pcapcal,lecrit,                        
     &     pplay,pplev,pzlay,pzlev,pz0,
     &     pu,pv,ph,pq,ptsrf,pemis,pqsurf,
     &     pdhfi,pdqfi,pfluxsrf,
     &     pdudif,pdvdif,pdhdif,pdtsrf,sensibFlux,pq2,
     &     pdqdif,pdqsdif,lastcall)

      use radcommon_h, only : sigma
      USE surfdat_h
      USE tracer_h
      use comcstfi_mod, only: g, r, cpp, rcp
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
      REAL hh

!     Tracers
!     -------
      INTEGER iq
      REAL zq(ngrid,nlay,nq)
      REAL zq1temp(ngrid)
      REAL rho(ngrid)         ! near-surface air density
      REAL qsat(ngrid)
      DATA firstcall/.true./
      REAL kmixmin

      real, parameter :: karman=0.4
      real cd0, roughratio

      real masse, Wtot, Wdiff

      real dqsdif_total(ngrid) 
      real zq0(ngrid) 


!     Coherence test
!     --------------

      IF (firstcall) THEN
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


!-----------------------------------------------------------------------
!     TRACERS (no vapour)
!     -------

      if(tracer) then

!     Calculate vertical flux from the bottom to the first layer (dust)
!     -----------------------------------------------------------------
         do ig=1,ngrid  
            rho(ig) = zb0(ig,1) /ptimestep
         end do

         call zerophys(ngrid*nq,pdqsdif)

!     Implicit inversion of q
!     -----------------------
         do iq=1,nq 

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
