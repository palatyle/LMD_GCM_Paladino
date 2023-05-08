      subroutine callcorrk(ngrid,nlayer,pq,nq,qsurf,zday,      &
          albedo,albedo_equivalent,emis,mu0,pplev,pplay,pt,    & 
          tsurf,fract,dist_star,                               &
          dtlw,dtsw,fluxsurf_lw,                               &
          fluxsurf_sw,fluxsurfabs_sw,fluxtop_lw,               &
          fluxabs_sw,fluxtop_dn,                               &
          OLR_nu,OSR_nu,                                       &
          int_dtaui,int_dtauv,                                 &
          lastcall)

      use mod_phys_lmdz_para, only : is_master
      use radinc_h
      use radcommon_h
      use gases_h
      USE tracer_h
      use callkeys_mod, only: global1d, szangle
      use comcstfi_mod, only: pi, mugaz, cpp
      use callkeys_mod, only: diurnal,tracer,seashaze,corrk_recombin,   &
                              strictboundcorrk,specOLR,diagdtau
      use geometry_mod, only: latitude

      implicit none

!==================================================================
!
!     Purpose
!     -------
!     Solve the radiative transfer using the correlated-k method for
!     the gaseous absorption and the Toon et al. (1989) method for
!     scatttering due to aerosols.
!
!     Authors
!     ------- 
!     Emmanuel 01/2001, Forget 09/2001
!     Robin Wordsworth (2009)
!     Jan Vatant d'Ollone (2018) -> corrk recombining case
!
!==================================================================

!-----------------------------------------------------------------------
!     Declaration of the arguments (INPUT - OUTPUT) on the LMD GCM grid
!     Layer #1 is the layer near the ground. 
!     Layer #nlayer is the layer at the top.
!-----------------------------------------------------------------------


      ! INPUT
      INTEGER,INTENT(IN) :: ngrid                  ! Number of atmospheric columns.
      INTEGER,INTENT(IN) :: nlayer                 ! Number of atmospheric layers.
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq)       ! Tracers (X/kg).
      INTEGER,INTENT(IN) :: nq                     ! Number of tracers.
      REAL,INTENT(IN) :: qsurf(ngrid,nq)           ! Tracers on surface (kg.m-2).
      REAL,INTENT(IN) :: zday		           ! Time elapsed since Ls=0 (sols).
      REAL,INTENT(IN) :: albedo(ngrid,L_NSPECTV)   ! Spectral Short Wavelengths Albedo. By MT2015
      REAL,INTENT(IN) :: emis(ngrid)               ! Long Wave emissivity.
      REAL,INTENT(IN) :: mu0(ngrid)                ! Cosine of sun incident angle.
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1)     ! Inter-layer pressure (Pa).
      REAL,INTENT(IN) :: pplay(ngrid,nlayer)       ! Mid-layer pressure (Pa).
      REAL,INTENT(IN) :: pt(ngrid,nlayer)          ! Air temperature (K).
      REAL,INTENT(IN) :: tsurf(ngrid)              ! Surface temperature (K).
      REAL,INTENT(IN) :: fract(ngrid)              ! Fraction of day.
      REAL,INTENT(IN) :: dist_star                 ! Distance star-planet (AU).
      logical,intent(in) :: lastcall               ! Signals last call to physics.
      
      ! OUTPUT
      REAL,INTENT(OUT) :: dtlw(ngrid,nlayer)             ! Heating rate (K/s) due to LW radiation.
      REAL,INTENT(OUT) :: dtsw(ngrid,nlayer)             ! Heating rate (K/s) due to SW radiation.
      REAL,INTENT(OUT) :: fluxsurf_lw(ngrid)             ! Incident LW flux to surf (W/m2).
      REAL,INTENT(OUT) :: fluxsurf_sw(ngrid)             ! Incident SW flux to surf (W/m2)
      REAL,INTENT(OUT) :: fluxsurfabs_sw(ngrid)          ! Absorbed SW flux by the surface (W/m2). By MT2015.
      REAL,INTENT(OUT) :: fluxtop_lw(ngrid)              ! Outgoing LW flux to space (W/m2).
      REAL,INTENT(OUT) :: fluxabs_sw(ngrid)              ! SW flux absorbed by the planet (W/m2).
      REAL,INTENT(OUT) :: fluxtop_dn(ngrid)              ! Incident top of atmosphere SW flux (W/m2).
      REAL,INTENT(OUT) :: OLR_nu(ngrid,L_NSPECTI)        ! Outgoing LW radition in each band (Normalized to the band width (W/m2/cm-1).
      REAL,INTENT(OUT) :: OSR_nu(ngrid,L_NSPECTV)        ! Outgoing SW radition in each band (Normalized to the band width (W/m2/cm-1).
      REAL,INTENT(OUT) :: albedo_equivalent(ngrid)       ! Spectrally Integrated Albedo. For Diagnostic. By MT2015
      REAL,INTENT(OUT) :: int_dtaui(ngrid,nlayer,L_NSPECTI) ! VI optical thickness of layers within narrowbands for diags ().
      REAL,INTENT(OUT) :: int_dtauv(ngrid,nlayer,L_NSPECTV) ! IR optical thickness of layers within narrowbands for diags ().
      
      
!-----------------------------------------------------------------------
!     Declaration of the variables required by correlated-k subroutines
!     Numbered from top to bottom (unlike in the GCM)
!-----------------------------------------------------------------------

      REAL*8 tmid(L_LEVELS),pmid(L_LEVELS)
      REAL*8 tlevrad(L_LEVELS),plevrad(L_LEVELS)

      ! Optical values for the optci/cv subroutines
      REAL*8 stel(L_NSPECTV),stel_fract(L_NSPECTV)
      REAL*8 dtaui(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 dtauv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 cosbv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 cosbi(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 wbari(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      REAL*8 wbarv(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 tauv(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      REAL*8 taucumv(L_LEVELS,L_NSPECTV,L_NGAUSS)
      REAL*8 taucumi(L_LEVELS,L_NSPECTI,L_NGAUSS)

      REAL*8 nfluxtopv,nfluxtopi,nfluxtop,fluxtopvdn
      REAL*8 nfluxoutv_nu(L_NSPECTV)                 ! Outgoing band-resolved VI flux at TOA (W/m2).
      REAL*8 nfluxtopi_nu(L_NSPECTI)                 ! Net band-resolved IR flux at TOA (W/m2).
      REAL*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI)         ! For 1D diagnostic.
      REAL*8 fmneti(L_NLAYRAD),fmnetv(L_NLAYRAD)
      REAL*8 fluxupv(L_NLAYRAD),fluxupi(L_NLAYRAD)
      REAL*8 fluxdnv(L_NLAYRAD),fluxdni(L_NLAYRAD)
      REAL*8 albi,acosz
      REAL*8 albv(L_NSPECTV)                         ! Spectral Visible Albedo.

      INTEGER ig,l,k,nw,iq,ip,ilay,it
      
      LOGICAL found

      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1)
      real*8 taugsurfi(L_NSPECTI,L_NGAUSS-1)

      logical OLRz
      real*8 NFLUXGNDV_nu(L_NSPECTV)
   
      ! Included by MT for albedo calculations.      
      REAL*8 albedo_temp(L_NSPECTV) ! For equivalent albedo calculation.
      REAL*8 surface_stellar_flux   ! Stellar flux reaching the surface. Useful for equivalent albedo calculation.
      
      ! For variable haze
      REAL*8 seashazefact(L_LEVELS)

      ! For muphys optics
      REAL*8 pqmo(ngrid,nlayer,nmicro)  ! Tracers for microphysics optics (X/m2).
      REAL*8 i2e(ngrid,nlayer)          ! int 2 ext factor ( X.kg-1 -> X.m-2 for optics )
      
      ! For corr-k recombining
      REAL*8 pqr(ngrid,L_PINT,L_REFVAR)  ! Tracers for corr-k recombining (mol/mol).
      REAL*8 fact, tmin, tmax
      
      LOGICAL usept(L_PINT,L_NTREF)  ! mask if pfref grid point will be used
      INTEGER inflay(L_PINT)         ! nearest inferior GCM layer for pfgasref grid points
      
!=======================================================================
!             I.  Initialization on every call   
!=======================================================================
 

      ! How much light do we get ?
      do nw=1,L_NSPECTV
         stel(nw)=stellarf(nw)/(dist_star**2)
      end do

      ! Convert (microphysical) tracers for optics: X.kg-1 --> X.m-2
      ! NOTE: it should be moved somewhere else: calmufi performs the same kind of
      ! computations... waste of time...
      i2e(:,1:nlayer) = ( pplev(:,1:nlayer)-pplev(:,2:nlayer+1) ) / gzlat(:,1:nlayer)
      pqmo(:,:,:) = 0.0
      DO iq=1,nmicro
        pqmo(:,:,iq) = pq(:,:,iq)*i2e(:,:)
      ENDDO

      ! Default value for fixed species for whom vmr has been taken
      ! into account while computing high-resolution spectra
      if (corrk_recombin) pqr(:,:,:) = 1.0

!-----------------------------------------------------------------------    
      do ig=1,ngrid ! Starting Big Loop over every GCM column
!-----------------------------------------------------------------------
         
         ! Recombine reference corr-k if needed
         if (corrk_recombin) then
         
         ! NB : To have decent CPU time recombining is not done on all gridpoints and wavelenghts but we
         ! calculate a gasi/v_recomb variable on the reference corrk-k T,P grid (only for T,P values
         ! who match the atmospheric conditions ) which is then processed as a standard pre-mix in
         ! optci/v routines, but updated every time tracers on the ref P grid have varied > 1%.

            ! Extract tracers for variable radiative species
            ! Also find the nearest GCM layer under each ref pressure
            do ip=1,L_PINT
               
               ilay=0
               found = .false.
               do l=1,nlayer
                  if ( pplay(ig,l) .gt. 10.0**(pfgasref(ip)+2.0) ) then ! pfgasref=log(p[mbar])
                     found=.true.
                     ilay=l
                  endif
               enddo       

               if (.not. found ) then ! set to min
                  do iq=1,L_REFVAR
                     if ( radvar_mask(iq) ) then
                        pqr(ig,ip,iq) = pq(ig,1,radvar_indx(iq)) / rat_mmol(radvar_indx(iq)-nmicro) ! mol/mol
                     endif
                  enddo
               else 
                  if (ilay==nlayer) then ! set to max
                     do iq=1,L_REFVAR
                        if ( radvar_mask(iq) ) then
                           pqr(ig,ip,iq) = pq(ig,nlayer,radvar_indx(iq)) / rat_mmol(radvar_indx(iq)-nmicro) ! mol/mol
                        endif
                     enddo
                  else ! standard
                     fact = ( 10.0**(pfgasref(ip)+2.0) - pplay(1,ilay+1) ) / ( pplay(1,ilay) - pplay(1,ilay+1) ) ! pfgasref=log(p[mbar])
                     do iq=1,L_REFVAR
                        if ( radvar_mask(iq) ) then
                           pqr(ig,ip,iq) = pq(ig,ilay,radvar_indx(iq))**fact * pq(ig,ilay+1,radvar_indx(iq))**(1.0-fact)
                           pqr(ig,ip,iq) = pqr(ig,ip,iq) / rat_mmol(radvar_indx(iq)-nmicro) ! mol/mol
                        endif
                     enddo
                  endif ! if ilay==nlayer
               endif ! if not found

               inflay(ip) = ilay

            enddo ! ip=1,L_PINT

            ! NB : The following usept is a trick to call recombine only for the reference T-P
            ! grid points that are useful given the temperature range at this altitude
            ! It saves a looot of time - JVO 18
            usept(:,:) = .true.
            do ip=1,L_PINT-1
              if ( inflay(ip+1)==nlayer ) then
                usept(ip,:) = .false.
              endif
              if ( inflay(ip) == 0 ) then
                usept(ip+1:,:) = .false.
              endif
              if ( usept(ip,1) ) then ! if not all false 
                tmin = minval(pt(ig,min(inflay(ip+1)+1,nlayer):max(inflay(ip),1)))
                tmax = maxval(pt(ig,min(inflay(ip+1)+1,nlayer):max(inflay(ip),1)))
                do it=1,L_NTREF-1
                  if ( tgasref(it+1) .lt. tmin ) then
                    usept(ip,it) = .false.
                  endif
                enddo
                do it=2,L_NTREF
                  if ( tgasref(it-1) .gt. tmax ) then
                    usept(ip,it) = .false.
                  endif
                enddo
                ! in case of out-of-bounds
                if ( tgasref(1)         .lt. tmin ) usept(ip,1) = .true.
                if ( tgasref(L_NTREF)   .gt. tmax ) usept(ip,L_NTREF) = .true.
              endif
            enddo ! ip=1,L_PINT-1
            ! deal with last bound
            if ( inflay(L_PINT-1).ne.0 ) usept(L_PINT,:) = usept(L_PINT-1,:)


            do ip=1,L_PINT

               ! Recombine k at (useful only!) reference T-P values if tracers or T have enough varied
               do it=1,L_NTREF

                 if ( usept(ip,it) .eqv. .false. ) cycle

                 do l=1,L_REFVAR
                    if ( abs( (pqr(ig,ip,l) - pqrold(ip,l)) / max(1.0e-30,pqrold(ip,l))) .GT. 0.01  & ! +- 1%
                         .or. ( useptold(ip,it) .eqv. .false.  ) ) then ! in case T change but not the tracers
                       call recombin_corrk( pqr(ig,ip,:),ip,it )
                       exit ! one is enough to trigger the update
                    endif
                 enddo
                 
               enddo

            enddo ! ip=1,L_PINT

            useptold(:,:)=usept(:,:)

       endif ! if corrk_recombin

!=======================================================================
!              II.  Transformation of the GCM variables
!=======================================================================


         ! Albedo and Emissivity.
         albi=1-emis(ig)   ! Long Wave.
         DO nw=1,L_NSPECTV ! Short Wave loop.
            albv(nw)=albedo(ig,nw)
         ENDDO

      if ((ngrid.eq.1).and.(global1d)) then ! Fixed zenith angle 'szangle' in 1D simulations w/ globally-averaged sunlight.
         acosz = cos(pi*szangle/180.0)
         print*,'acosz=',acosz,', szangle=',szangle
      else
         acosz=mu0(ig) ! Cosine of sun incident angle : 3D simulations or local 1D simulations using latitude.
      endif

!-----------------------------------------------------------------------
!     Pressure and temperature
!-----------------------------------------------------------------------

      DO l=1,nlayer
         plevrad(2*l)   = pplay(ig,nlayer+1-l)/scalep
         plevrad(2*l+1) = pplev(ig,nlayer+1-l)/scalep
         tlevrad(2*l)   = pt(ig,nlayer+1-l)
         tlevrad(2*l+1) = (pt(ig,nlayer+1-l)+pt(ig,max(nlayer-l,1)))/2
      END DO
      
      plevrad(1) = 0.
      plevrad(2) = 0.   !! Trick to have correct calculations of fluxes in gflux(i/v).F, but the pmid levels are not impacted by this change.

      tlevrad(1) = tlevrad(2)
      tlevrad(2*nlayer+1)=tsurf(ig)
      
      pmid(1) = max(pgasmin,0.0001*plevrad(3))
      pmid(2) =  pmid(1)

      tmid(1) = tlevrad(2)
      tmid(2) = tmid(1)
    
      DO l=1,L_NLAYRAD-1
         tmid(2*l+1) = tlevrad(2*l+1)
         tmid(2*l+2) = tlevrad(2*l+1)
         pmid(2*l+1) = plevrad(2*l+1)
         pmid(2*l+2) = plevrad(2*l+1)
      END DO
      pmid(L_LEVELS) = plevrad(L_LEVELS)
      tmid(L_LEVELS) = tlevrad(L_LEVELS)

!!Alternative interpolation:
!         pmid(3) = pmid(1)
!         pmid(4) = pmid(1) 
!         tmid(3) = tmid(1)
!         tmid(4) = tmid(1)
!      DO l=2,L_NLAYRAD-1
!         tmid(2*l+1) = tlevrad(2*l)
!         tmid(2*l+2) = tlevrad(2*l)
!         pmid(2*l+1) = plevrad(2*l)
!         pmid(2*l+2) = plevrad(2*l)
!      END DO
!      pmid(L_LEVELS) = plevrad(L_LEVELS-1)
!      tmid(L_LEVELS) = tlevrad(L_LEVELS-1)

      ! Test for out-of-bounds pressure.
      if(plevrad(3).lt.pgasmin)then
         print*,'Minimum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         call abort
      elseif(plevrad(L_LEVELS).gt.pgasmax)then
         print*,'Maximum pressure is outside the radiative'
         print*,'transfer kmatrix bounds, exiting.'
         call abort
      endif

      ! Test for out-of-bounds temperature.
      do k=1,L_LEVELS
         if(tlevrad(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds'
            print*,"k=",k," tlevrad(k)=",tlevrad(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tlevrad=tgasmin' 
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              !tlevrad(k)=tgasmin
            endif
         elseif(tlevrad(k).gt.tgasmax)then
!            print*,'Maximum temperature is outside the radiative'
!            print*,'transfer kmatrix bounds, exiting.'
!            print*,"k=",k," tlevrad(k)=",tlevrad(k)
!            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              call abort
            else
!              print*,'***********************************************'
!              print*,'we allow model to continue with tlevrad=tgasmax'  
!              print*,'  ... we assume we know what you are doing ... '
!              print*,'  ... but do not let this happen too often ... '
!              print*,'***********************************************'
              !tlevrad(k)=tgasmax
            endif
         endif
      enddo
      do k=1,L_NLAYRAD+1
         if(tmid(k).lt.tgasmin)then
            print*,'Minimum temperature is outside the radiative'
            print*,'transfer kmatrix bounds, exiting.'
            print*,"k=",k," tmid(k)=",tmid(k)
            print*,"tgasmin=",tgasmin
            if (strictboundcorrk) then
              call abort
            else
              print*,'***********************************************'
              print*,'we allow model to continue with tmid=tgasmin'
              print*,'  ... we assume we know what you are doing ... '
              print*,'  ... but do not let this happen too often ... '
              print*,'***********************************************'
              tmid(k)=tgasmin
            endif
         elseif(tmid(k).gt.tgasmax)then
!            print*,'Maximum temperature is outside the radiative'
!            print*,'transfer kmatrix bounds, exiting.'
!            print*,"k=",k," tmid(k)=",tmid(k)
!            print*,"tgasmax=",tgasmax
            if (strictboundcorrk) then
              call abort
            else
!              print*,'***********************************************'
!              print*,'we allow model to continue with tmid=tgasmin'
!              print*,'  ... we assume we know what you are doing ... '
!              print*,'  ... but do not let this happen too often ... '
!              print*,'***********************************************'
              tmid(k)=tgasmax
            endif
         endif
      enddo

!=======================================================================
!          III. Calling the main radiative transfer subroutines
!=======================================================================

         Cmk(:)      = 0.01 * 1.0 / (gzlat(ig,:) * mugaz * 1.672621e-27) ! q_main=1.0 assumed.
         gzlat_ig(:) = gzlat(ig,:)
         
         ! Compute the haze seasonal modulation factor
         if (seashaze) call season_haze(zday,latitude(ig),plevrad,seashazefact)

!-----------------------------------------------------------------------
!        Short Wave Part
!-----------------------------------------------------------------------

         if(fract(ig) .ge. 1.0e-4) then ! Only during daylight.
            if((ngrid.eq.1).and.(global1d))then
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw)* 0.25 / acosz ! globally averaged = divide by 4, and we correct for solar zenith angle
               end do
            else
               do nw=1,L_NSPECTV
                  stel_fract(nw)= stel(nw) * fract(ig)
               end do
            endif
            
            call optcv(pqmo(ig,:,:),nlayer,plevrad,tmid,pmid,   &
                 dtauv,tauv,taucumv,wbarv,cosbv,tauray,taugsurf,seashazefact)

            call sfluxv(dtauv,tauv,taucumv,albv,dwnv,wbarv,cosbv,  &
                 acosz,stel_fract,                                 &
                 nfluxtopv,fluxtopvdn,nfluxoutv_nu,nfluxgndv_nu,   &
                 fmnetv,fluxupv,fluxdnv,fzerov,taugsurf)

         else ! During the night, fluxes = 0.
            nfluxtopv       = 0.0d0
            fluxtopvdn      = 0.0d0
            nfluxoutv_nu(:) = 0.0d0
            nfluxgndv_nu(:) = 0.0d0
            do l=1,L_NLAYRAD
               fmnetv(l)=0.0d0
               fluxupv(l)=0.0d0
               fluxdnv(l)=0.0d0
            end do
         end if


         ! Equivalent Albedo Calculation (for OUTPUT). MT2015
         if(fract(ig) .ge. 1.0e-4) then ! equivalent albedo makes sense only during daylight.       
            surface_stellar_flux=sum(nfluxgndv_nu(1:L_NSPECTV))      
            if(surface_stellar_flux .gt. 1.0e-3) then ! equivalent albedo makes sense only if the stellar flux received by the surface is positive.
               DO nw=1,L_NSPECTV                  
                  albedo_temp(nw)=albedo(ig,nw)*nfluxgndv_nu(nw)
               ENDDO
               albedo_temp(1:L_NSPECTV)=albedo_temp(1:L_NSPECTV)/surface_stellar_flux
               albedo_equivalent(ig)=sum(albedo_temp(1:L_NSPECTV))
            else
               albedo_equivalent(ig)=0.0 ! Spectrally Integrated Albedo not defined for non-irradiated grid points. So we arbitrary set the equivalent albedo to 0.
            endif
         else
            albedo_equivalent(ig)=0.0 ! Spectrally Integrated Albedo not defined for non-irradiated grid points. So we arbitrary set the equivalent albedo to 0.
         endif


!-----------------------------------------------------------------------
!        Long Wave Part
!-----------------------------------------------------------------------

         call optci(pqmo(ig,:,:),nlayer,plevrad,tlevrad,tmid,pmid,   &
              dtaui,taucumi,cosbi,wbari,taugsurfi,seashazefact)

         call sfluxi(plevrad,tlevrad,dtaui,taucumi,ubari,albi,      &
              wnoi,dwni,cosbi,wbari,nfluxtopi,nfluxtopi_nu,         & 
              fmneti,fluxupi,fluxdni,fluxupi_nu,fzeroi,taugsurfi)

!-----------------------------------------------------------------------
!     Transformation of the correlated-k code outputs
!     (into dtlw, dtsw, fluxsurf_lw, fluxsurf_sw, fluxtop_lw, fluxtop_sw)

!     Flux incident at the top of the atmosphere
         fluxtop_dn(ig)=fluxtopvdn 

         fluxtop_lw(ig)  = real(nfluxtopi)
         fluxabs_sw(ig)  = real(-nfluxtopv)
         fluxsurf_lw(ig) = real(fluxdni(L_NLAYRAD))
         fluxsurf_sw(ig) = real(fluxdnv(L_NLAYRAD))
         
!        Flux absorbed by the surface. By MT2015.         
         fluxsurfabs_sw(ig) = fluxsurf_sw(ig)*(1.-albedo_equivalent(ig))

         if(fluxtop_dn(ig).lt.0.0)then
            print*,'Achtung! fluxtop_dn has lost the plot!'
            print*,'fluxtop_dn=',fluxtop_dn(ig)
            print*,'acosz=',acosz
            print*,'temp=   ',pt(ig,:)
            print*,'pplay=  ',pplay(ig,:)
            call abort
         endif

!     Spectral output, for exoplanet observational comparison
         if(specOLR)then
            do nw=1,L_NSPECTI 
               OLR_nu(ig,nw)=nfluxtopi_nu(nw)/DWNI(nw) !JL Normalize to the bandwidth
            end do
            do nw=1,L_NSPECTV 
               !GSR_nu(ig,nw)=nfluxgndv_nu(nw)
               OSR_nu(ig,nw)=nfluxoutv_nu(nw)/DWNV(nw) !JL Normalize to the bandwidth
            end do
         endif

!     Finally, the heating rates

         DO l=2,L_NLAYRAD
            dtsw(ig,L_NLAYRAD+1-l)=(fmnetv(l)-fmnetv(l-1))  &
                *gzlat(ig,L_NLAYRAD+1-l)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
            dtlw(ig,L_NLAYRAD+1-l)=(fmneti(l)-fmneti(l-1))  &
                *gzlat(ig,L_NLAYRAD+1-l)/(cpp*scalep*(plevrad(2*l+1)-plevrad(2*l-1)))
         END DO      

!     These are values at top of atmosphere
         dtsw(ig,L_NLAYRAD)=(fmnetv(1)-nfluxtopv)           &
             *gzlat(ig,L_NLAYRAD)/(cpp*scalep*(plevrad(3)-plevrad(1)))
         dtlw(ig,L_NLAYRAD)=(fmneti(1)-nfluxtopi)           &
             *gzlat(ig,L_NLAYRAD)/(cpp*scalep*(plevrad(3)-plevrad(1)))


      !  Optical thickness diagnostics (added by JVO)
      if (diagdtau) then
        do l=1,L_NLAYRAD
          do nw=1,L_NSPECTV
            int_dtauv(ig,l,nw) = 0.0d0
             DO k=1,L_NGAUSS
              ! Output exp(-tau) because gweight ponderates exp and not tau itself
              int_dtauv(ig,l,nw)= int_dtauv(ig,l,nw) + exp(-dtauv(l,nw,k))*gweight(k)
             ENDDO
          enddo
          do nw=1,L_NSPECTI
           int_dtaui(ig,l,nw) = 0.0d0
             DO k=1,L_NGAUSS
              ! Output exp(-tau) because gweight ponderates exp and not tau itself
              int_dtaui(ig,l,nw)= int_dtaui(ig,l,nw) + exp(-dtaui(l,nw,k))*gweight(k)
             ENDDO
          enddo
        enddo
      endif        


!-----------------------------------------------------------------------    
      end do ! End of big loop over every GCM column.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     Additional diagnostics
!-----------------------------------------------------------------------

      ! IR spectral output, for exoplanet observational comparison
      if(lastcall.and.(ngrid.eq.1))then  ! could disable the 1D output, they are in the diagfi and diagspec... JL12

         print*,'Saving scalar quantities in surf_vals.out...'
         print*,'psurf = ', pplev(1,1),' Pa'
         open(116,file='surf_vals.out')
         write(116,*) tsurf(1),pplev(1,1),fluxtop_dn(1),         &
                      real(-nfluxtopv),real(nfluxtopi) 
         close(116)


!          USEFUL COMMENT - Do Not Remove.
!
!           if(specOLR)then
!               open(117,file='OLRnu.out')
!               do nw=1,L_NSPECTI
!                  write(117,*) OLR_nu(1,nw)
!               enddo
!               close(117)
!
!               open(127,file='OSRnu.out')
!               do nw=1,L_NSPECTV
!                  write(127,*) OSR_nu(1,nw)
!               enddo
!               close(127)
!           endif

           ! OLR vs altitude: do it as a .txt file.
         OLRz=.false.
         if(OLRz)then
            print*,'saving IR vertical flux for OLRz...'
            open(118,file='OLRz_plevs.out')
            open(119,file='OLRz.out')
            do l=1,L_NLAYRAD
               write(118,*) plevrad(2*l)
               do nw=1,L_NSPECTI
                  write(119,*) fluxupi_nu(l,nw) 
               enddo
            enddo 
            close(118)
            close(119)
         endif

      endif

      if (lastcall) then
        IF( ALLOCATED( gasi ) ) DEALLOCATE( gasi )
        IF( ALLOCATED( gasv ) ) DEALLOCATE( gasv )
        IF( ALLOCATED( gasi_recomb ) ) DEALLOCATE( gasi_recomb )
        IF( ALLOCATED( gasv_recomb ) ) DEALLOCATE( gasv_recomb )
        IF( ALLOCATED( pqrold ) ) DEALLOCATE( pqrold )
        IF( ALLOCATED( useptold ) ) DEALLOCATE( useptold )
!$OMP BARRIER
!$OMP MASTER
        IF( ALLOCATED( pgasref ) ) DEALLOCATE( pgasref )
        IF( ALLOCATED( tgasref ) ) DEALLOCATE( tgasref )
        IF( ALLOCATED( pfgasref ) ) DEALLOCATE( pfgasref )
        IF( ALLOCATED( gweight ) ) DEALLOCATE( gweight )
!$OMP END MASTER
!$OMP BARRIER
      endif


    end subroutine callcorrk
