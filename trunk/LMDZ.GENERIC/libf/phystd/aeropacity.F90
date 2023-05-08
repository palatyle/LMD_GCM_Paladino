      Subroutine aeropacity(ngrid,nlayer,nq,pplay,pplev,pq, &
         aerosol,reffrad,QREFvis3d,QREFir3d,tau_col,cloudfrac,totcloudfrac,clearsky)

       use radinc_h, only : L_TAUMAX,naerkind
       use aerosol_mod
       USE tracer_h, only: noms,rho_co2,rho_ice
       use comcstfi_mod, only: g, pi
       use geometry_mod, only: latitude
       use callkeys_mod, only: aerofixco2,aerofixh2o,kastprof,cloudlvl,	&
		CLFvarying,CLFfixval,dusttau,			 	&
		pres_bottom_tropo,pres_top_tropo,obs_tau_col_tropo,	&
		pres_bottom_strato,pres_top_strato,obs_tau_col_strato,  &
                tau_nh3_cloud, pres_nh3_cloud 
                  
       implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Compute aerosol optical depth in each gridbox.
!     
!     Authors
!     ------- 
!     F. Forget
!     F. Montmessin (water ice scheme) 
!     update J.-B. Madeleine (2008)
!     dust removal, simplification by Robin Wordsworth (2009)
!
!     Input
!     ----- 
!     ngrid             Number of horizontal gridpoints
!     nlayer            Number of layers
!     nq                Number of tracers
!     pplev             Pressure (Pa) at each layer boundary
!     pq                Aerosol mixing ratio
!     reffrad(ngrid,nlayer,naerkind)         Aerosol effective radius
!     QREFvis3d(ngrid,nlayer,naerkind) \ 3d extinction coefficients
!     QREFir3d(ngrid,nlayer,naerkind)  / at reference wavelengths
!
!     Output
!     ------
!     aerosol            Aerosol optical depth in layer l, grid point ig
!     tau_col            Total column optical depth at grid point ig
!
!=======================================================================

      INTEGER,INTENT(IN) :: ngrid  ! number of atmospheric columns
      INTEGER,INTENT(IN) :: nlayer ! number of atmospheric layers
      INTEGER,INTENT(IN) :: nq     ! number of tracers
      REAL,INTENT(IN) :: pplay(ngrid,nlayer) ! mid-layer pressure (Pa)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq) ! tracers (.../kg_of_air)
      REAL,INTENT(OUT) :: aerosol(ngrid,nlayer,naerkind) ! aerosol optical depth
      REAL,INTENT(IN) :: reffrad(ngrid,nlayer,naerkind) ! aerosol effective radius
      REAL,INTENT(IN) :: QREFvis3d(ngrid,nlayer,naerkind) ! extinction coefficient in the visible
      REAL,INTENT(IN) :: QREFir3d(ngrid,nlayer,naerkind)
      REAL,INTENT(OUT):: tau_col(ngrid) !column integrated visible optical depth
      ! BENJAMIN MODIFS
      real,intent(in) :: cloudfrac(ngrid,nlayer) ! cloud fraction
      real,intent(out) :: totcloudfrac(ngrid) ! total cloud fraction
      logical,intent(in) :: clearsky

      real aerosol0, obs_tau_col_aurora, pm, pente_cloud
      
      real dp_strato(ngrid)
      real dp_tropo(ngrid)

      INTEGER l,ig,iq,iaer

      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)
      REAL CBRT
      EXTERNAL CBRT

      INTEGER,SAVE :: i_co2ice=0      ! co2 ice
      INTEGER,SAVE :: i_h2oice=0      ! water ice
!$OMP THREADPRIVATE(i_co2ice,i_h2oice)
      CHARACTER(LEN=20) :: tracername ! to temporarily store text

      ! for fixed dust profiles
      real topdust, expfactor, zp
      REAL taudusttmp(ngrid) ! Temporary dust opacity used before scaling
      REAL tauh2so4tmp(ngrid) ! Temporary h2so4 opacity used before scaling

      real CLFtot

      ! identify tracers
      IF (firstcall) THEN

        write(*,*) "Tracers found in aeropacity:"
        do iq=1,nq
          tracername=noms(iq)
          if (tracername.eq."co2_ice") then
            i_co2ice=iq
          write(*,*) "i_co2ice=",i_co2ice

          endif
          if (tracername.eq."h2o_ice") then
            i_h2oice=iq
            write(*,*) "i_h2oice=",i_h2oice
          endif
        enddo

        if (noaero) then
          print*, "No active aerosols found in aeropacity"
        else
          print*, "If you would like to use aerosols, make sure any old"
          print*, "start files are updated in newstart using the option"
          print*, "q=0"
          write(*,*) "Active aerosols found in aeropacity:"
        endif

        if ((iaero_co2.ne.0).and.(.not.noaero)) then
          print*, 'iaero_co2=  ',iaero_co2
        endif
        if (iaero_h2o.ne.0) then
          print*,'iaero_h2o=  ',iaero_h2o    
        endif
        if (iaero_dust.ne.0) then
          print*,'iaero_dust= ',iaero_dust
        endif
        if (iaero_h2so4.ne.0) then
          print*,'iaero_h2so4= ',iaero_h2so4
        endif
        if (iaero_back2lay.ne.0) then
          print*,'iaero_back2lay= ',iaero_back2lay
        endif
        if (iaero_nh3.ne.0) then
          print*,'iaero_nh3= ',iaero_nh3
        endif
        if (iaero_aurora.ne.0) then
          print*,'iaero_aurora= ',iaero_aurora
        endif

        firstcall=.false.
      ENDIF ! of IF (firstcall)


!     ---------------------------------------------------------
!==================================================================
!    CO2 ice aerosols
!==================================================================

      if (iaero_co2.ne.0) then
           iaer=iaero_co2
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation
            if (noaero) then ! aerosol set to zero
             aerosol(1:ngrid,1:nlayer,iaer)=0.0
            elseif (aerofixco2.or.(i_co2ice.eq.0)) then !  CO2 ice cloud prescribed
               aerosol(1:ngrid,1:nlayer,iaer)=1.e-9
               !aerosol(1:ngrid,12,iaer)=4.0 ! single cloud layer option
            else
               DO ig=1, ngrid
                  DO l=1,nlayer-1 ! to stop the rad tran bug

                     aerosol0 =                         &
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_co2 * reffrad(ig,l,iaer) )  ) *   &
                          ( pq(ig,l,i_co2ice) + 1.E-9 ) *         &
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g
                     aerosol0           = max(aerosol0,1.e-9)
                     aerosol0           = min(aerosol0,L_TAUMAX)
                     aerosol(ig,l,iaer) = aerosol0
!                     aerosol(ig,l,iaer) = 0.0
!                     print*, aerosol(ig,l,iaer)
!        using cloud fraction
!                     aerosol(ig,l,iaer) = -log(1 - CLF + CLF*exp(-aerosol0/CLF))
!                     aerosol(ig,l,iaer) = min(aerosol(ig,l,iaer),L_TAUMAX)


                  ENDDO
               ENDDO
            end if ! if fixed or varying
      end if ! if CO2 aerosols   
!==================================================================
!     Water ice / liquid 
!==================================================================

      if (iaero_h2o.ne.0) then 
           iaer=iaero_h2o
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation
            if (aerofixh2o.or.(i_h2oice.eq.0).or.clearsky) then
               aerosol(1:ngrid,1:nlayer,iaer) =1.e-9

               ! put cloud at cloudlvl
               if(kastprof.and.(cloudlvl.ne.0.0))then
                  ig=1
                  do l=1,nlayer
                     if(int(cloudlvl).eq.l)then
                     !if(cloudlvl.gt.(pplay(ig,l)/pplev(ig,1)))then
                        print*,'Inserting cloud at level ',l
                        !aerosol(ig,l,iaer)=10.0

                        rho_ice=920.0

                        ! the Kasting approximation
                        aerosol(ig,l,iaer) =                      &
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_ice * reffrad(ig,l,iaer) )  ) *   &
                          !( pq(ig,l,i_h2oice) + 1.E-9 ) *         &
                          ( 4.0e-4 + 1.E-9 ) *         &
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g


                        open(115,file='clouds.out',form='formatted')
                        write(115,*) l,aerosol(ig,l,iaer)
                        close(115)

                        return
                     endif
                  end do

                  call abort
               endif

            else

               do ig=1, ngrid
                  !do l=1,nlayer-1 ! to stop the rad tran bug
                  do l=1,nlayer !JL18 if aerosols are present in the last layer we must account for them. Provides better upper boundary condition in the IR. They must however be put to zero in the sw (see optcv)
                                ! same correction should b-probably be done for other aerosol types.
                     aerosol(ig,l,iaer) =                                    & !modification by BC
                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
                          ( rho_ice * reffrad(ig,l,iaer) )  ) *   &
                          !  pq(ig,l,i_h2oice) *                   & !JL I dropped the +1e-9 here to have the same
                          !( pplev(ig,l) - pplev(ig,l+1) ) / g       !   opacity in the clearsky=true and the 
                                                                     !   clear=false/pq=0 case
                          ( pq(ig,l,i_h2oice) + 1.E-9 ) *         & ! Doing this makes the code unstable, so I have restored it (RW)
                          ( pplev(ig,l) - pplev(ig,l+1) ) / g 

                  enddo
               enddo

               if(CLFvarying)then
                  call totalcloudfrac(ngrid,nlayer,nq,cloudfrac,totcloudfrac,pplev,pq,aerosol(1,1,iaer))
                  do ig=1, ngrid
                     !do l=1,nlayer-1 ! to stop the rad tran bug
                     do l=1,nlayer !JL18 if aerosols are present in the last layer we must account for them. Provides better upper boundary condition in the IR. They must however be put to zero in the sw (see optcv)
                        CLFtot  = max(totcloudfrac(ig),0.01)
                        aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/CLFtot
                        aerosol(ig,l,iaer) = max(aerosol(ig,l,iaer),1.e-9)
                     enddo
                  enddo
               else
                  do ig=1, ngrid
                     !do l=1,nlayer-1 ! to stop the rad tran bug
                     do l=1,nlayer !JL18 if aerosols are present in the last layer we must account for them. Provides better upper boundary condition in the IR. They must however be put to zero in the sw (see optcv)
                        CLFtot  = CLFfixval
                        aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/CLFtot
                        aerosol(ig,l,iaer) = max(aerosol(ig,l,iaer),1.e-9)
                     enddo
                  enddo
              end if!(CLFvarying)
            endif !(aerofixed.or.(i_h2oice.eq.0).or.clearsky)
	      
      end if ! End if h2o aerosol

!==================================================================
!             Dust 
!==================================================================
      if (iaero_dust.ne.0) then
          iaer=iaero_dust
!         1. Initialization 
          aerosol(1:ngrid,1:nlayer,iaer)=0.0
          
          topdust=30.0 ! km  (used to be 10.0 km) LK

!       2. Opacity calculation

!           expfactor=0.
           DO l=1,nlayer-1
             DO ig=1,ngrid
!             Typical mixing ratio profile

                 zp=(pplev(ig,1)/pplay(ig,l))**(70./topdust)
                 expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)

!             Vertical scaling function
              aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1)) &
               *expfactor


             ENDDO
           ENDDO

!          Rescaling each layer to reproduce the choosen (or assimilated)
!          dust extinction opacity at visible reference wavelength, which
!          is scaled to the surface pressure pplev(ig,1)

            taudusttmp(1:ngrid)=0.
              DO l=1,nlayer
                DO ig=1,ngrid
                   taudusttmp(ig) = taudusttmp(ig) &
                          +  aerosol(ig,l,iaer)
                ENDDO
              ENDDO
            DO l=1,nlayer-1
               DO ig=1,ngrid
                  aerosol(ig,l,iaer) = max(1E-20, &
                          dusttau &
                       *  pplev(ig,1) / pplev(ig,1) &
                       *  aerosol(ig,l,iaer) &
                       /  taudusttmp(ig))

              ENDDO
            ENDDO
      end if ! If dust aerosol   

!==================================================================
!           H2SO4 
!==================================================================
! added by LK
      if (iaero_h2so4.ne.0) then
         iaer=iaero_h2so4

!       1. Initialization
         aerosol(1:ngrid,1:nlayer,iaer)=0.0


!       2. Opacity calculation

!           expfactor=0.
         DO l=1,nlayer-1
            DO ig=1,ngrid
!              Typical mixing ratio profile

               zp=(pplev(ig,1)/pplay(ig,l))**(70./30) !emulating topdust
               expfactor=max(exp(0.007*(1.-max(zp,1.))),1.e-3)

!             Vertical scaling function
               aerosol(ig,l,iaer)= (pplev(ig,l)-pplev(ig,l+1))*expfactor

            ENDDO
         ENDDO
         tauh2so4tmp(1:ngrid)=0.
         DO l=1,nlayer
            DO ig=1,ngrid
               tauh2so4tmp(ig) = tauh2so4tmp(ig) + aerosol(ig,l,iaer)
            ENDDO
         ENDDO
         DO l=1,nlayer-1
            DO ig=1,ngrid
               aerosol(ig,l,iaer) = max(1E-20, &
                          1 &
                       *  pplev(ig,1) / pplev(ig,1) &
                       *  aerosol(ig,l,iaer) &
                       /  tauh2so4tmp(ig))

            ENDDO
         ENDDO

! 1/700. is assuming a "sulfurtau" of 1
! Sulfur aerosol routine to be improved.
!                     aerosol0 =                         &
!                          (  0.75 * QREFvis3d(ig,l,iaer) /        &
!                          ( rho_h2so4 * reffrad(ig,l,iaer) )  ) *   &
!                          ( pq(ig,l,i_h2so4) + 1.E-9 ) *         &
!                          ( pplev(ig,l) - pplev(ig,l+1) ) / g
!                     aerosol0           = max(aerosol0,1.e-9)
!                     aerosol0           = min(aerosol0,L_TAUMAX)
!                     aerosol(ig,l,iaer) = aerosol0

!                  ENDDO
!               ENDDO
      end if
 
           
!     ---------------------------------------------------------
!==================================================================
!    Two-layer aerosols (unknown composition)
!    S. Guerlet (2013) - Modif by J. Vatant d'Ollone (2020)
!==================================================================

      if (iaero_back2lay .ne.0) then
           iaer=iaero_back2lay
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.0
!       2. Opacity calculation


!       JVO 20 : Modif to have each of the layers (strato and tropo) correctly normalized
!                Otherwise we previously had the total optical depth correct but for each
!                separately, so  it didn't match the input values + what's more normalizing
!                to the sum was making them non-independent : eg changing tau_tropo was
!                affecting stratopsheric values of optical depth ...
!
!                Note that the main consequence of the former version bug was (in most cases)
!                to strongly underestimate the stratospheric optical depths compared to the
!                required values, eg, with tau_tropo=10 and tau_strato=0.1, you actually ended
!                with an actual tau_strato of 1E-4 ... !
!
!                NB : Because of the extra transition opacity if the layers are non contiguous,
!                be aware that at the the bottom we have tau > tau_strato + tau_tropo

         DO ig=1,ngrid
          dp_tropo(ig)  = 0.D0
          dp_strato(ig) = 0.D0
          DO l=1,nlayer-1
             aerosol(ig,l,iaer) = ( pplev(ig,l) - pplev(ig,l+1) )
             !! 1. below tropospheric layer: no aerosols
             IF (pplev(ig,l) .gt. pres_bottom_tropo) THEN
               aerosol(ig,l,iaer) = 0.D0
             !! 2. tropo layer
             ELSEIF (pplev(ig,l) .le. pres_bottom_tropo .and. pplev(ig,l) .ge. pres_top_tropo) THEN
               dp_tropo(ig) = dp_tropo(ig) + aerosol(ig,l,iaer)
             !! 3. linear transition 
             ! JVO 20 : This interpolation needs to be done AFTER we set strato and tropo (see below)
             !! 4. strato layer 
             ELSEIF (pplev(ig,l) .le. pres_bottom_strato .and. pplev(ig,l) .ge. pres_top_strato) THEN
               dp_strato(ig) = dp_strato(ig) + aerosol(ig,l,iaer)
             !! 5. above strato layer: no aerosols
             ELSEIF (pplev(ig,l) .lt. pres_top_strato) THEN
               aerosol(ig,l,iaer) = 0.D0
             ENDIF
	  ENDDO
         ENDDO

!       3. Re-normalize to the (input) observed (total) column (for each of the layers)

         DO ig=1,ngrid
          DO l=1,nlayer-1
               IF (pplev(ig,l) .le. pres_bottom_tropo .and. pplev(ig,l) .ge. pres_top_tropo) THEN
                 aerosol(ig,l,iaer) = obs_tau_col_tropo*aerosol(ig,l,iaer)/dp_tropo(ig)
               ELSEIF (pplev(ig,l) .lt. pres_top_tropo .and. pplev(ig,l) .gt. pres_bottom_strato) THEN
                 expfactor=log(pplev(ig,l)/pres_top_tropo)/log(pres_bottom_strato/pres_top_tropo)
                 aerosol(ig,l,iaer) = (obs_tau_col_strato/dp_strato(ig))**expfactor     &
                                    * (obs_tau_col_tropo/dp_tropo(ig))**(1.0-expfactor) &
                                    * aerosol(ig,l,iaer)
               ELSEIF (pplev(ig,l) .le. pres_bottom_strato .and. pplev(ig,l) .ge. pres_top_strato) THEN
                 aerosol(ig,l,iaer) = obs_tau_col_strato*aerosol(ig,l,iaer)/dp_strato(ig)
               ENDIF
               write(*,*), aerosol(ig,l,iaer)
            ENDDO
         ENDDO


      end if ! if Two-layer aerosols  

!==================================================================
!    Saturn/Jupiter ammonia cloud = thin cloud (scale height 0.2 hard coded...)
!    S. Guerlet (2013)
!==================================================================

      if (iaero_nh3 .ne.0) then
           iaer=iaero_nh3
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.D0
!       2. Opacity calculation
          DO ig=1,ngrid

           DO l=1,nlayer-1
            !! 1. below cloud layer: no opacity
	    
            IF (pplev(ig,l) .gt. pres_nh3_cloud ) THEN
            aerosol(ig,l,iaer) = 0.D0            

             ELSEIF (pplev(ig,l) .le. pres_nh3_cloud ) THEN 
	     pente_cloud=5. !!(hard-coded, correspond to scale height 0.2)
             aerosol(ig,l,iaer) = ((pplev(ig,l)/pres_nh3_cloud)**(pente_cloud))*tau_nh3_cloud 

             ENDIF
            ENDDO

          END DO
	  
!       3. Re-normalize to observed total column
         tau_col(:)=0.0
         DO l=1,nlayer
          DO ig=1,ngrid
               tau_col(ig) = tau_col(ig) &
                     + aerosol(ig,l,iaer)/tau_nh3_cloud
            ENDDO
         ENDDO

         DO ig=1,ngrid
           DO l=1,nlayer-1
                aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/tau_col(ig)
           ENDDO
         ENDDO

     end if ! if NH3 cloud  
!==================================================================
!    Jovian auroral aerosols (unknown composition) NON-GENERIC: vertical and meridional profile tuned to observations
!    S. Guerlet (2015)
!==================================================================


      if (iaero_aurora .ne.0) then
           iaer=iaero_aurora
!       1. Initialization
            aerosol(1:ngrid,1:nlayer,iaer)=0.D0 
	 pm = 2000. !!case study: maxi aerosols at 20 hPa
!       2. Opacity calculation
          DO ig=1,ngrid

	  !! Test Jupiter (based on Zhang et al 2013 observations, but a bit different), decembre 2015
              DO l=1,nlayer
       		aerosol(ig,l,iaer) = (pplev(ig,l)/pm)**2 * exp(-(pplev(ig,l)/pm)**2)
              ENDDO
          ENDDO
	 
 !       3. Meridional distribution, and re-normalize to observed total column
         tau_col(:)=0.D0
         DO ig=1,ngrid
	  !!Jupiter
	  !!Hem sud:
          IF (latitude(ig)*180.D0/pi .lt. -45.D0 .and. latitude(ig)*180.D0/pi .gt. -70.) THEN
 	  obs_tau_col_aurora= 10.D0**(-0.06D0*latitude(ig)*180.D0/pi-3.4D0) 
          ELSEIF (latitude(ig)*180.D0/pi .lt. -37.D0 .and. latitude(ig)*180.D0/pi .ge. -45.) THEN
 	  obs_tau_col_aurora= 10.D0**(-0.3D0*latitude(ig)*180.D0/pi-14.3D0) 
           ELSEIF (latitude(ig)*180./pi .le. -70. ) THEN
 	  obs_tau_col_aurora= 10**(0.06*70.-3.4D0) 
	  !!Hem Nord:  
          ELSEIF (latitude(ig)*180.D0/pi .gt. 30.D0 .and. latitude(ig)*180.D0/pi .lt. 70.) THEN
	  obs_tau_col_aurora= 10.D0**(0.03D0*latitude(ig)*180.D0/pi-1.17D0)  
          ELSEIF (latitude(ig)*180.D0/pi .gt. 22.D0 .and. latitude(ig)*180.D0/pi .le. 30.) THEN
	  obs_tau_col_aurora= 10.D0**(0.3D0*latitude(ig)*180.D0/pi-9.4D0)  
          ELSEIF (latitude(ig)*180.D0/pi .ge. 70.) THEN
	  obs_tau_col_aurora= 10**(0.03*70.-1.17D0)  
          ELSEIF (latitude(ig)*180.D0/pi .ge. -37. .and. latitude(ig)*180.D0/pi .le. 22.) THEN
	 obs_tau_col_aurora = 0.001D0    !!Jupiter: mini pas a zero
	  ENDIF

 	  DO l=1,nlayer  
               tau_col(ig) = tau_col(ig) + aerosol(ig,l,iaer)/obs_tau_col_aurora
          ENDDO
         ENDDO

         DO ig=1,ngrid
           DO l=1,nlayer-1
                aerosol(ig,l,iaer)=aerosol(ig,l,iaer)/tau_col(ig)
           ENDDO
         ENDDO


      end if ! if Auroral aerosols  


! --------------------------------------------------------------------------
! Column integrated visible optical depth in each point (used for diagnostic)

      tau_col(:)=0.0
      do iaer = 1, naerkind
         do l=1,nlayer
            do ig=1,ngrid
               tau_col(ig) = tau_col(ig) + aerosol(ig,l,iaer)
            end do
         end do
      end do

      do ig=1,ngrid
         do l=1,nlayer
            do iaer = 1, naerkind
               if(aerosol(ig,l,iaer).gt.1.e3)then
                  print*,'WARNING: aerosol=',aerosol(ig,l,iaer)
                  print*,'at ig=',ig,',  l=',l,', iaer=',iaer
                  print*,'QREFvis3d=',QREFvis3d(ig,l,iaer)
                  print*,'reffrad=',reffrad(ig,l,iaer)
               endif
            end do
         end do
      end do

      do ig=1,ngrid
         if(tau_col(ig).gt.1.e3)then
            print*,'WARNING: tau_col=',tau_col(ig)
            print*,'at ig=',ig
            print*,'aerosol=',aerosol(ig,:,:)
            print*,'QREFvis3d=',QREFvis3d(ig,:,:)
            print*,'reffrad=',reffrad(ig,:,:)
         endif
      end do
      return
    end subroutine aeropacity
      
