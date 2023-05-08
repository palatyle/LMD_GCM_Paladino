subroutine hydrol(ngrid,nq,ptimestep,rnat,tsurf,  &
     qsurf,dqsurf,dqs_hyd,pcapcal,                &
     albedo,albedo_bareground,                    &
     albedo_snow_SPECTV,albedo_co2_ice_SPECTV,    &
     mu0,pdtsurf,pdtsurf_hyd,hice,                &
     pctsrf_sic,sea_ice)

  use ioipsl_getin_p_mod, only: getin_p
  use watercommon_h, only: T_h2O_ice_liq, RLFTT, rhowater, mx_eau_sol
  USE surfdat_h
  use comdiurn_h
  USE geometry_mod, only: cell_area
  USE tracer_h
  use slab_ice_h
  use callkeys_mod, only: albedosnow,albedoco2ice,ok_slab_ocean,Tsaldiff,maxicethick,co2cond
  use radinc_h, only : L_NSPECTV

  implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Calculate the surface hydrology and albedo changes.
!     
!     Authors
!     ------- 
!     Adapted from LMDTERRE by B. Charnay (2010). Further 
!     Modifications by R. Wordsworth (2010).
!     Spectral albedo by M. Turbet (2015).
!     
!     Called by
!     ---------
!     physiq.F
!     
!     Calls
!     -----
!     none
!
!     Notes
!     -----
!     rnat is terrain type: 0-ocean; 1-continent
!     
!==================================================================

      integer ngrid,nq

!     Inputs
!     ------
      real snowlayer
      parameter (snowlayer=33.0)        ! 33 kg/m^2 of snow, equal to a layer of 3.3 cm 
      real oceantime
      parameter (oceantime=10*24*3600)

      logical,save :: oceanbulkavg ! relax ocean temperatures to a GLOBAL mean value?
      logical,save :: activerunoff ! enable simple runoff scheme?
      logical,save :: oceanalbvary ! ocean albedo varies with the diurnal cycle?
!$OMP THREADPRIVATE(oceanbulkavg,activerunoff,oceanalbvary)

!     Arguments
!     ---------
      real rnat(ngrid) ! I changed this to integer (RW)
      real,dimension(:),allocatable,save :: runoff
      real totalrunoff, tsea, oceanarea
      save oceanarea
!$OMP THREADPRIVATE(runoff,oceanarea)

      real ptimestep
      real mu0(ngrid)
      real qsurf(ngrid,nq), tsurf(ngrid)
      real dqsurf(ngrid,nq), pdtsurf(ngrid)
      real hice(ngrid)
      real albedo(ngrid,L_NSPECTV)
      real albedo_bareground(ngrid)
      real albedo_snow_SPECTV(L_NSPECTV)
      real albedo_co2_ice_SPECTV(L_NSPECTV)
      real pctsrf_sic(ngrid), sea_ice(ngrid)

      real oceanarea2

!     Output
!     ------
      real dqs_hyd(ngrid,nq)
      real pdtsurf_hyd(ngrid)

!     Local
!     -----
      real a,b,E
      integer ig,iq, nw
      real fsnoi, subli, fauxo
      real twater(ngrid)
      real pcapcal(ngrid)
      real hicebis(ngrid)
      real zqsurf(ngrid,nq)
      real ztsurf(ngrid)
      real albedo_sic, alb_ice
      real zfra

      integer, save :: ivap, iliq, iice
!$OMP THREADPRIVATE(ivap,iliq,iice)

      logical, save :: firstcall
!$OMP THREADPRIVATE(firstcall)

      data firstcall /.true./


      if(firstcall)then

         oceanbulkavg=.false.
         oceanalbvary=.false.
         write(*,*)"Activate runnoff into oceans?"
         activerunoff=.false.
         call getin_p("activerunoff",activerunoff)
         write(*,*)" activerunoff = ",activerunoff
	 
	 
	 
         if (activerunoff) then
           ALLOCATE(runoff(ngrid))
           runoff(1:ngrid)=0
         endif

         ivap=igcm_h2o_vap
         iliq=igcm_h2o_vap
         iice=igcm_h2o_ice
        
         write(*,*) "hydrol: ivap=",ivap
         write(*,*) "        iliq=",iliq
         write(*,*) "        iice=",iice

!     Here's the deal: iice is used in place of igcm_h2o_ice both on the 
!                      surface and in the atmosphere. ivap is used in
!                      place of igcm_h2o_vap ONLY in the atmosphere, while
!                      iliq is used in place of igcm_h2o_vap ONLY on the 
!                      surface.
!                      Soon to be extended to the entire water cycle...

!     Total ocean surface area
         oceanarea=0.
         do ig=1,ngrid
            if(nint(rnat(ig)).eq.0)then
               oceanarea=oceanarea+cell_area(ig)
            endif
         enddo

         if(oceanbulkavg.and.(oceanarea.le.0.))then
            print*,'How are we supposed to average the ocean'
            print*,'temperature, when there are no oceans?'
            call abort
         endif

         if(activerunoff.and.(oceanarea.le.0.))then
            print*,'You have enabled runoff, but you have no oceans.'
            print*,'Where did you think the water was going to go?'
            call abort
         endif
         
         firstcall = .false.
      endif

!     add physical tendencies already calculated
!     ------------------------------------------

      do ig=1,ngrid
         ztsurf(ig) = tsurf(ig) + ptimestep*pdtsurf(ig)
         pdtsurf_hyd(ig)=0.0
         do iq=1,nq
            zqsurf(ig,iq) = qsurf(ig,iq) + ptimestep*dqsurf(ig,iq)
         enddo    
      enddo
 
      do ig=1,ngrid
         do iq=1,nq
            dqs_hyd(ig,iq) = 0.0
         enddo
      enddo

      do ig = 1, ngrid

!     Ocean
!     -----
         if(nint(rnat(ig)).eq.0)then

!     re-calculate oceanic albedo
!            if(diurnal.and.oceanalbvary)then
!               fauxo      = ( 1.47 - ACOS( mu0(ig) ) )/0.15 ! where does this come from (Benjamin)?
!               albedo(ig) = 1.1*( .03 + .630/( 1. + fauxo*fauxo))
!               albedo(ig) = MAX(MIN(albedo(ig),0.60),0.04)
!            else
               do nw=1,L_NSPECTV
                  albedo(ig,nw) = alb_ocean ! For now, alb_ocean is defined in slab_ice_h.F90. Later we could introduce spectral dependency for alb_ocean.
	       enddo
!            end if


            if(ok_slab_ocean) then
          
               zfra = MAX(0.0,MIN(1.0,zqsurf(ig,iice)/45.0))     ! Snow Fraction (Critical height 45kg/m2~15cm)
               alb_ice=alb_ice_max-(alb_ice_max-alb_ice_min) &   ! Ice Albedo
               *exp(-sea_ice(ig)/h_alb_ice)
               ! Albedo final calculation :
               do nw=1,L_NSPECTV
                  albedo(ig,nw) = pctsrf_sic(ig)*                                        &
                                 (albedo_snow_SPECTV(nw)*zfra + alb_ice*(1.0-zfra))      &
                               + (1.-pctsrf_sic(ig))*alb_ocean
	       enddo

               ! Oceanic ice height, just for diagnostics
               hice(ig)    = MIN(10.,sea_ice(ig)/rhowater)
            else !ok_slab_ocean


!     calculate oceanic ice height including the latent heat of ice formation
!     hice is the height of oceanic ice with a maximum of maxicethick.
               hice(ig)    = zqsurf(ig,iice)/rhowater ! update hice to include recent snowfall

!              twater(ig)  = tsurf(ig) + ptimestep*zdtsurf(ig) &
               twater(ig)  = ztsurf(ig) - hice(ig)*RLFTT*rhowater/pcapcal(ig)
               ! this is temperature water would have if we melted entire ocean ice layer
               hicebis(ig) = hice(ig)
               hice(ig)    = 0.

               if(twater(ig) .lt. T_h2O_ice_liq)then
                  E=min((T_h2O_ice_liq+Tsaldiff-twater(ig))*pcapcal(ig),RLFTT*rhowater*maxicethick)
                  hice(ig)        = E/(RLFTT*rhowater)
                  hice(ig)        = max(hice(ig),0.0)
                  hice(ig)        = min(hice(ig),maxicethick)
                  pdtsurf_hyd(ig) = (hice(ig) - hicebis(ig))*RLFTT*rhowater/pcapcal(ig)/ptimestep              
                  do nw=1,L_NSPECTV
		     albedo(ig,nw) = albedo_snow_SPECTV(nw) ! Albedo of ice has been replaced by albedo of snow here. MT2015.
	          enddo  

!                 if (zqsurf(ig,iice).ge.snowlayer) then
!                    albedo(ig) = albedoice
!                 else
!                    albedo(ig) = albedoocean &
!                    + (albedosnow - albedoocean)*zqsurf(ig,iice)/snowlayer
!                 endif

               else

                  pdtsurf_hyd(ig) = -hicebis(ig)*RLFTT*rhowater/pcapcal(ig)/ptimestep
                  DO nw=1,L_NSPECTV
                     albedo(ig,nw) = alb_ocean
	          ENDDO               

               endif

               zqsurf(ig,iliq) = zqsurf(ig,iliq)-(hice(ig)*rhowater-zqsurf(ig,iice))
               zqsurf(ig,iice) = hice(ig)*rhowater

            endif!(ok_slab_ocean)


!     Continent
!     ---------
         elseif (nint(rnat(ig)).eq.1) then

!     melt the snow
            if(ztsurf(ig).gt.T_h2O_ice_liq)then
               if(zqsurf(ig,iice).gt.1.0e-8)then

                  a     = (ztsurf(ig)-T_h2O_ice_liq)*pcapcal(ig)/RLFTT
                  b     = zqsurf(ig,iice)
                  fsnoi = min(a,b)

                  zqsurf(ig,iice) = zqsurf(ig,iice) - fsnoi
                  zqsurf(ig,iliq) = zqsurf(ig,iliq) + fsnoi

!                 thermal effects
                  pdtsurf_hyd(ig) = -fsnoi*RLFTT/pcapcal(ig)/ptimestep  

               endif
            else

!     freeze the water
               if(zqsurf(ig,iliq).gt.1.0e-8)then

                  a     = -(ztsurf(ig)-T_h2O_ice_liq)*pcapcal(ig)/RLFTT
                  b     = zqsurf(ig,iliq)
                  
                  fsnoi = min(a,b)

                  zqsurf(ig,iice) = zqsurf(ig,iice) + fsnoi
                  zqsurf(ig,iliq) = zqsurf(ig,iliq) - fsnoi

!                 thermal effects
                  pdtsurf_hyd(ig) = +fsnoi*RLFTT/pcapcal(ig)/ptimestep  

               endif
            endif
            
!     deal with runoff
            if(activerunoff)then

               runoff(ig) = max(zqsurf(ig,iliq) - mx_eau_sol, 0.0)
               if(ngrid.gt.1)then ! runoff only exists in 3D
                  if(runoff(ig).ne.0.0)then
                     zqsurf(ig,iliq) = mx_eau_sol
!                    runoff is added to ocean at end
                  endif
               end if

            endif

!     re-calculate continental albedo
            DO nw=1,L_NSPECTV
	       albedo(ig,nw) = albedo_bareground(ig)
	    ENDDO
            if (zqsurf(ig,iice).ge.snowlayer) then
               DO nw=1,L_NSPECTV
                  albedo(ig,nw) = albedo_snow_SPECTV(nw)
	       ENDDO
            else
               DO nw=1,L_NSPECTV
                  albedo(ig,nw) = albedo_bareground(ig)                            & 
                               + (albedo_snow_SPECTV(nw) - albedo_bareground(ig))  &
                                 *zqsurf(ig,iice)/snowlayer
	       ENDDO
            endif

         else

            print*,'Surface type not recognised in hydrol.F!'
            print*,'Exiting...'
            call abort

         endif

      end do ! ig=1,ngrid

!     perform crude bulk averaging of temperature in ocean
!     ----------------------------------------------------
      if(oceanbulkavg)then

         oceanarea2=0.       
         DO ig=1,ngrid
            if((nint(rnat(ig)).eq.0).and.(hice(ig).eq.0.))then
               oceanarea2=oceanarea2+cell_area(ig)*pcapcal(ig)
            end if
         END DO
       
         tsea=0.
         DO ig=1,ngrid
            if((nint(rnat(ig)).eq.0).and.(hice(ig).eq.0.))then       
               tsea=tsea+ztsurf(ig)*cell_area(ig)*pcapcal(ig)/oceanarea2
            end if
         END DO

         DO ig=1,ngrid
            if((nint(rnat(ig)).eq.0).and.(hice(ig).eq.0))then
               pdtsurf_hyd(ig) = pdtsurf_hyd(ig) + (tsea-ztsurf(ig))/oceantime
            end if
         END DO

         print*,'Mean ocean temperature               = ',tsea,' K'

      endif

!     shove all the runoff water into the ocean
!     -----------------------------------------
      if(activerunoff)then

         totalrunoff=0.
         do ig=1,ngrid
            if (nint(rnat(ig)).eq.1) then
               totalrunoff = totalrunoff + cell_area(ig)*runoff(ig)
            endif
         enddo
         
         do ig=1,ngrid
            if (nint(rnat(ig)).eq.0) then
               zqsurf(ig,iliq) = zqsurf(ig,iliq) + &
                    totalrunoff/oceanarea
            endif
         enddo

      endif         


!     Re-add the albedo effects of CO2 ice if necessary
!     -------------------------------------------------
      if(co2cond)then

         do ig=1,ngrid
            if (qsurf(ig,igcm_co2_ice).gt.1.) then ! Condition changed - Need now ~1 mm CO2 ice coverage. MT2015
	       DO nw=1,L_NSPECTV
                  albedo(ig,nw) = albedo_co2_ice_SPECTV(nw)
	       ENDDO
            endif
         enddo ! ngrid
         
      endif ! co2cond


      do ig=1,ngrid ! We calculate here the tracer tendencies. Don't forget that we have to retrieve the dqsurf tendencies we added at the beginning of the routine !
         dqs_hyd(ig,iliq)=(zqsurf(ig,iliq) - qsurf(ig,iliq))/ptimestep - dqsurf(ig,iliq)
         dqs_hyd(ig,iice)=(zqsurf(ig,iice) - qsurf(ig,iice))/ptimestep - dqsurf(ig,iice)
      enddo

      if (activerunoff) then
        call writediagfi(ngrid,'runoff','Runoff amount',' ',2,runoff)
      endif

      return
    end subroutine hydrol
