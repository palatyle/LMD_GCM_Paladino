










      subroutine condense_co2(ngrid,nlayer,nq,ptimestep, &
          pcapcal,pplay,pplev,ptsrf,pt,                  &
          pdt,pdtsrf,pq,pdq,                             &
          pqsurf,pdqsurfc,albedo,pemisurf,               &
          albedo_bareground,albedo_co2_ice_SPECTV,       &
          pdtc,pdtsrfc,pdpsrfc,pdqc)

      use radinc_h, only : L_NSPECTV, naerkind
      use gases_h, only: gfrac, igas_co2
      use radii_mod, only : co2_reffrad
      use aerosol_mod, only : iaero_co2
      USE surfdat_h, only: emisice, emissiv
      USE geometry_mod, only: latitude ! in radians
      USE tracer_h, only: noms, rho_co2
      use comcstfi_mod, only: g, r, cpp

      implicit none

!==================================================================
!     Purpose
!     -------
!     Condense and/or sublime CO2 ice on the ground and in the atmosphere, and sediment the ice.
!  
!     Inputs
!     ------
!     ngrid                 Number of vertical columns.
!     nlayer                Number of vertical layers.
!     nq                    Number of tracers.
!     ptimestep             Duration of the physical timestep (s).
!     pplay(ngrid,nlayer)   Pressure layers (Pa).
!     pplev(ngrid,nlayer+1) Pressure levels (Pa).
!     pt(ngrid,nlayer)      Atmospheric Temperatures (K).
!     ptsrf(ngrid)          Surface temperatures (K).
!     pq(ngrid,nlayer,nq)   Atmospheric tracers mixing ratios (kg/kg of air).
!     pqsurf(ngrid,nq)      Surface tracers (kg/m2).
!     
!     pdt(ngrid,nlayer)     Time derivative before condensation/sublimation of pt.
!     pdtsrf(ngrid)         Time derivative before condensation/sublimation of ptsrf.
!     pdq(ngrid,nlayer,nq)  Time derivative before condensation/sublimation of
!
!     albedo_bareground(ngrid)           Albedo of the bare ground.
!     albedo_co2_ice_SPECTV(L_NSPECTV)   Spectral albedo of CO2 ice.
!     
!     Outputs
!     -------
!     pdpsrfc(ngrid)          \       Contribution of condensation/sublimation 
!     pdtc(ngrid,nlayer)       \            to the time derivatives of 
!     pdtsrfc(ngrid)           /     Surface Pressure, Atmospheric Temperatures,
!     pdqsurfc(ngrid)         /         Surface Temperatures, Surface Tracers,
!     pdqc(ngrid,nlayer,nq)  /                and Atmospheric Tracers.*
!
!     pemisurf(ngrid)              Emissivity of the surface.
!     
!     Both
!     ----
!     albedo(ngrid,L_NSPECTV)      Spectral albedo of the surface.
!
!     Authors
!     ------- 
!     Francois Forget (1996)
!     Converted to Fortran 90 and slightly modified by R. Wordsworth (2009)
!     Includes simplifed nucleation by J. Leconte (2011)
!     
!==================================================================

!--------------------------
!        Arguments
!--------------------------


      INTEGER,INTENT(IN) :: ngrid
      INTEGER,INTENT(IN) :: nlayer
      INTEGER,INTENT(IN) :: nq
      REAL,INTENT(IN) :: ptimestep 
      REAL,INTENT(IN) :: pcapcal(ngrid)
      REAL,INTENT(IN) :: pplay(ngrid,nlayer)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1)
      REAL,INTENT(IN) :: ptsrf(ngrid)
      REAL,INTENT(IN) :: pt(ngrid,nlayer)
      REAL,INTENT(IN) :: pdt(ngrid,nlayer)
      REAL,INTENT(IN) :: pdtsrf(ngrid)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq)
      REAL,INTENT(IN) :: pqsurf(ngrid,nq)
      REAL,INTENT(IN) :: pdq(ngrid,nlayer,nq)
      REAL,INTENT(IN) :: albedo_bareground(ngrid)
      REAL,INTENT(IN) :: albedo_co2_ice_SPECTV(L_NSPECTV)
      REAL,INTENT(INOUT) :: albedo(ngrid,L_NSPECTV)
      REAL,INTENT(OUT) :: pemisurf(ngrid)
      REAL,INTENT(OUT) :: pdtc(ngrid,nlayer)
      REAL,INTENT(OUT) :: pdtsrfc(ngrid)
      REAL,INTENT(OUT) :: pdpsrfc(ngrid)
      REAL,INTENT(OUT) :: pdqc(ngrid,nlayer,nq)
      REAL,INTENT(OUT) :: pdqsurfc(ngrid)

!------------------------------
!       Local variables
!------------------------------

      INTEGER l,ig,icap,ilay,iq,nw,igas,it

      REAL reffrad(ngrid,nlayer) ! Radius (m) of the CO2 ice particles.
      REAL*8 zt(ngrid,nlayer)    ! Updated Atmospheric Temperatures (K).
      REAL ztsrf(ngrid)          ! Updated Surface Temperatures (K).
      REAL zq(ngrid,nlayer,nq)   ! Updated Atmospheric tracers mixing ratios (kg/kg of air).
      REAL piceco2(ngrid)        ! Updated Surface Tracer (kg/m2).
      REAL ztcond (ngrid,nlayer) ! Atmospheric Temperatures of condensation of CO2.
      REAL ztnuc (ngrid,nlayer)  ! Atmospheric Nucleation Temperatures.
      REAL ztcondsol(ngrid)      ! Temperatures of condensation of CO2 at the surface.
      REAL zcondices(ngrid)      ! Condensation rate on the ground (kg/m2/s).
      REAL zfallice(ngrid)       ! Flux of ice falling on the surface (kg/m2/s).
      REAL Mfallice(ngrid)       ! Total amount of ice fallen to the ground during the timestep (kg/m2).
      REAL wq(ngrid,nlayer+1)    ! Total amount of ice fallen to the ground during the timestep (kg/m2).
      REAL subptimestep          ! Duration of the subtimestep (s) for the sedimentation.
      Integer Ntime              ! Number of subtimesteps.
      REAL masse (ngrid,nlayer)  ! Mass of atmospheric layers (kg/m2)
      REAL w(ngrid,nlayer,nq)    !
      REAL vstokes,reff          !
      REAL ppco2                 !


!------------------------------------------
!         Saved local variables
!------------------------------------------


      REAL,SAVE :: latcond=5.9e5
      REAL,SAVE :: ccond
      REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: emisref
!$OMP THREADPRIVATE(latcond,ccond,emisref)
      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)
      INTEGER,SAVE :: i_co2ice=0      ! co2 ice
!$OMP THREADPRIVATE(i_co2ice)
      CHARACTER(LEN=20) :: tracername ! to temporarily store text


!------------------------------------------------
!       Initialization at the first call
!------------------------------------------------


      IF (firstcall) THEN

         ALLOCATE(emisref(ngrid))
         ! Find CO2 ice tracer.
         do iq=1,nq
            tracername=noms(iq)
            if (tracername.eq."co2_ice") then
               i_co2ice=iq
            endif
         enddo
         
         write(*,*) "condense_co2: i_co2ice=",i_co2ice       

         if((i_co2ice.lt.1))then
            print*,'In condens_cloud but no CO2 ice tracer, exiting.'
            print*,'Still need generalisation to arbitrary species!'
            stop
         endif

         ccond=cpp/(g*latcond)
         print*,'In condens_cloud: ccond=',ccond,' latcond=',latcond

         ! Prepare special treatment if gas is not pure CO2
         ! if (addn2) then
         !    m_co2   = 44.01E-3 ! CO2 molecular mass (kg/mol)   
         !    m_noco2 = 28.02E-3 ! N2 molecular mass (kg/mol)  
         ! Compute A and B coefficient use to compute
         ! mean molecular mass Mair defined by
         ! 1/Mair = q(ico2)/m_co2 + (1-q(ico2))/m_noco2
         ! 1/Mair = A*q(ico2) + B
         ! A = (1/m_co2 - 1/m_noco2)
         ! B = 1/m_noco2
         ! endif

         ! Minimum CO2 mixing ratio below which mixing occurs with layer above : qco2min =0.75  

           firstcall=.false.
      ENDIF


!------------------------------------------------
!        Tendencies initially set to 0
!------------------------------------------------


      pdqc(1:ngrid,1:nlayer,1:nq)     = 0.
      pdtc(1:ngrid,1:nlayer)          = 0.
      zq(1:ngrid,1:nlayer,1:nq)       = 0.
      zt(1:ngrid,1:nlayer)            = 0.
      Mfallice(1:ngrid)               = 0.      
      zfallice(1:ngrid)               = 0.
      zcondices(1:ngrid)              = 0.
      pdtsrfc(1:ngrid)                = 0.
      pdpsrfc(1:ngrid)                = 0.
      pdqsurfc(1:ngrid)               = 0.


!----------------------------------
!     Atmospheric condensation
!----------------------------------


!     Compute CO2 Volume mixing ratio
!     -------------------------------
!      if (addn2) then
!         DO l=1,nlayer
!            DO ig=1,ngrid
!              qco2=pq(ig,l,ico2)+pdq(ig,l,ico2)*ptimestep
!              Mean air molecular mass = 1/(q(ico2)/m_co2 + (1-q(ico2))/m_noco2)
!              mmean=1/(A*qco2 +B)
!              vmr_co2(ig,l) = qco2*mmean/m_co2 
!            ENDDO
!         ENDDO
!      else
!         DO l=1,nlayer
!            DO ig=1,ngrid
!              vmr_co2(ig,l)=0.5
!            ENDDO
!         ENDDO
!      end if


      ! Forecast the atmospheric frost temperature 'ztcond' and nucleation temperature 'ztnuc'.
      DO l=1,nlayer
         DO ig=1,ngrid
            ppco2=gfrac(igas_CO2)*pplay(ig,l)
            call get_tcond_co2(ppco2,ztcond(ig,l))
            call get_tnuc_co2(ppco2,ztnuc(ig,l))
         ENDDO
      ENDDO
      
      ! Initialize zq and zt at the beginning of the sub-timestep loop and qsurf.
      DO ig=1,ngrid
         piceco2(ig)=pqsurf(ig,i_co2ice)
         DO l=1,nlayer
            zt(ig,l)=pt(ig,l)
            zq(ig,l,i_co2ice)=pq(ig,l,i_co2ice)
            IF( zq(ig,l,i_co2ice).lt.-1.e-6 ) THEN
               print*,'Uh-oh, zq = ',zq(ig,l,i_co2ice),'at ig,l=',ig,l
               if(l.eq.1)then
                  print*,'Perhaps the atmosphere is collapsing on surface...?'
               endif
            END IF
         ENDDO
      ENDDO

      ! Calculate the mass of each atmospheric layer (kg.m-2)
      do  ilay=1,nlayer
         DO ig=1,ngrid
            masse(ig,ilay)=(pplev(ig,ilay) - pplev(ig,ilay+1)) /g
         end do
      end do


!-----------------------------------------------------------
!     START CONDENSATION/SEDIMENTATION SUB-TIME LOOP
!-----------------------------------------------------------


      Ntime =  20  ! number of sub-timestep 
      subptimestep = ptimestep/float(Ntime)           

      ! Add the tendencies from other physical processes at each subtimstep.
      DO it=1,Ntime
         DO l=1,nlayer
            DO ig=1,ngrid
               zt(ig,l)   = zt(ig,l)   + pdt(ig,l)   * subptimestep
               zq(ig,l,i_co2ice) = zq(ig,l,i_co2ice) + pdq(ig,l,i_co2ice) * subptimestep
            END DO
         END DO

         ! Gravitational sedimentation starts.
            
         ! Sedimentation computed from radius computed from q in module radii_mod.
	 call co2_reffrad(ngrid,nlayer,nq,zq,reffrad)
	 
         DO  ilay=1,nlayer
            DO ig=1,ngrid

               reff = reffrad(ig,ilay)

               call stokes                      &
                   (pplev(ig,ilay),pt(ig,ilay), &
                    reff,vstokes,rho_co2)

               !w(ig,ilay,i_co2ice) = 0.0
               w(ig,ilay,i_co2ice) = vstokes *  subptimestep * &
                   pplev(ig,ilay)/(r*pt(ig,ilay))

            END DO
         END DO

         ! Computing q after sedimentation
         call vlz_fi(ngrid,nlayer,zq(1,1,i_co2ice),2.,masse,w(1,1,i_co2ice),wq)


         ! Progressively accumulating the flux to the ground.
         ! Mfallice is the total amount of ice fallen to the ground.
         DO ig=1,ngrid
            Mfallice(ig) =  Mfallice(ig) + wq(ig,i_co2ice)
         END DO

!----------------------------------------------------------
!       Condensation / sublimation in the atmosphere
!----------------------------------------------------------
!     (MODIFICATIONS FOR EARLY MARS: falling heat neglected, condensation of CO2 into tracer i_co2ice)


         DO l=nlayer , 1, -1
            DO ig=1,ngrid
               pdtc(ig,l)=0.

               ! ztcond-> ztnuc in test beneath to nucleate only when super saturation occurs(JL 2011)
               IF ((zt(ig,l).LT.ztnuc(ig,l)).or.(zq(ig,l,i_co2ice).gt.1.E-10)) THEN 
                  pdtc(ig,l)   = (ztcond(ig,l) - zt(ig,l))/subptimestep
                  pdqc(ig,l,i_co2ice) = pdtc(ig,l)*ccond*g

                  ! Case when the ice from above sublimes entirely
                  IF ((zq(ig,l,i_co2ice).lt.-pdqc(ig,l,i_co2ice)*subptimestep) &
                      .AND. (zq(ig,l,i_co2ice).gt.0)) THEN

                     pdqc(ig,l,i_co2ice) = -zq(ig,l,i_co2ice)/subptimestep
                     pdtc(ig,l)   =-zq(ig,l,i_co2ice)/(ccond*g*subptimestep)

                  END IF

                  ! Temperature and q after condensation
                  zt(ig,l)   = zt(ig,l)   + pdtc(ig,l)   * subptimestep
                  zq(ig,l,i_co2ice) = zq(ig,l,i_co2ice) + pdqc(ig,l,i_co2ice) * subptimestep
               END IF

            ENDDO
         ENDDO
         
      ENDDO! end of subtimestep loop.

      ! Computing global tendencies after the subtimestep.
      DO l=1,nlayer
         DO ig=1,ngrid
            pdtc(ig,l) = &
              (zt(ig,l) - (pt(ig,l) + pdt(ig,l)*ptimestep))/ptimestep
            pdqc(ig,l,i_co2ice) = &
              (zq(ig,l,i_co2ice)-(pq(ig,l,i_co2ice)+pdq(ig,l,i_co2ice)*ptimestep))/ptimestep
         END DO
      END DO
      DO ig=1,ngrid
         zfallice(ig) = Mfallice(ig)/ptimestep
      END DO


!-----------------------------------------------------------------------
!              Condensation/sublimation on the ground
!-----------------------------------------------------------------------


      ! Forecast of ground temperature ztsrf and frost temperature ztcondsol.
      DO ig=1,ngrid
         ppco2=gfrac(igas_CO2)*pplay(ig,1)
         call get_tcond_co2(ppco2,ztcondsol(ig))
         
         ztsrf(ig) = ptsrf(ig)

         if((ztsrf(ig).le.ztcondsol(ig)+2.0).and.(ngrid.eq.1))then
            print*,'CO2 is condensing on the surface in 1D. This atmosphere is doomed.'
            print*,'T_surf = ',ztsrf,'K'
            print*,'T_cond = ',ztcondsol,'K'
            open(116,file='surf_vals.out')
            write(116,*) 0.0, pplev(1,1), 0.0, 0.0
            close(116)
            call abort
         endif

         ztsrf(ig) = ptsrf(ig) + pdtsrf(ig)*ptimestep

      ENDDO
     
      DO ig=1,ngrid
      
         IF(ig.GT.ngrid/2+1) THEN
            icap=2
         ELSE
            icap=1
         ENDIF

         ! Loop over where we have condensation / sublimation
         IF ((ztsrf(ig) .LT. ztcondsol(ig)) .OR.           &        ! ground condensation
                    (zfallice(ig).NE.0.) .OR.              &        ! falling snow
                    ((ztsrf(ig) .GT. ztcondsol(ig)) .AND.  &        ! ground sublimation
                    ((piceco2(ig)+zfallice(ig)*ptimestep) .NE. 0.))) THEN


            ! Condensation or partial sublimation of CO2 ice
            zcondices(ig)=pcapcal(ig)*(ztcondsol(ig)-ztsrf(ig)) &
                /(latcond*ptimestep)         
            pdtsrfc(ig) = (ztcondsol(ig) - ztsrf(ig))/ptimestep

            ! If the entire CO_2 ice layer sublimes
            ! (including what has just condensed in the atmosphere)
            IF((piceco2(ig)/ptimestep+zfallice(ig)).LE. &
                -zcondices(ig))THEN
               zcondices(ig) = -piceco2(ig)/ptimestep - zfallice(ig)
               pdtsrfc(ig)=(latcond/pcapcal(ig))*       &
                   (zcondices(ig))
            END IF

            ! Changing CO2 ice amount and pressure
            piceco2(ig)  =  piceco2(ig)   + pdqsurfc(ig)*ptimestep
            pdqsurfc(ig) =  zcondices(ig) + zfallice(ig)
            pdpsrfc(ig)  = -pdqsurfc(ig)*g

            IF(ABS(pdpsrfc(ig)*ptimestep).GT.pplev(ig,1)) THEN
               PRINT*,'STOP in condens in condense_co2'
               PRINT*,'condensing more than total mass'
               PRINT*,'Grid point ',ig
               PRINT*,'Ps = ',pplev(ig,1)
               PRINT*,'d Ps = ',pdpsrfc(ig)
               STOP
            ENDIF
         END IF
         
      ENDDO ! end of ngrid loop.


!---------------------------------------------------------------------------------------------
!     Surface albedo and emissivity of the ground below the snow (emisref)
!---------------------------------------------------------------------------------------------


      DO ig=1,ngrid
      
         IF(latitude(ig).LT.0.) THEN
            icap=2 ! Southern Hemisphere
         ELSE
            icap=1 ! Nortnern hemisphere
         ENDIF

         if(.not.piceco2(ig).ge.0.) THEN
            if(piceco2(ig).le.-1.e-8) print*,   &
              'WARNING : in condense_co2cloud: piceco2(',ig,')=', piceco2(ig)
            piceco2(ig)=0.
         endif
         if (piceco2(ig) .gt. 1.) then  ! CO2 Albedo condition changed to ~1 mm coverage. Change by MT2015.
	    DO nw=1,L_NSPECTV
               albedo(ig,nw) = albedo_co2_ice_SPECTV(nw)
	    ENDDO
            emisref(ig)   = emisice(icap)
         else
            DO nw=1,L_NSPECTV
               albedo(ig,nw) = albedo_bareground(ig) ! Note : If you have some water, it will be taken into account in the "hydrol" routine.
	    ENDDO
            emisref(ig)   = emissiv
            pemisurf(ig)  = emissiv
         end if
         
      END DO

      return
      
      end subroutine condense_co2



!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------



      subroutine get_tcond_co2(p,tcond)  ! Calculates the condensation temperature for CO2
     

      implicit none

      real p, peff, tcond
      real, parameter :: ptriple=518000.0

      peff=p

      if(peff.lt.ptriple) then
         tcond = (-3167.8)/(log(.01*peff)-23.23) ! Fanale's formula.
      else
         tcond = 684.2-92.3*log(peff)+4.32*log(peff)**2 ! liquid-vapour transition (based on CRC handbook 2003 data)        
      endif
      return

      end subroutine get_tcond_co2



!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------



      subroutine get_tnuc_co2(p,tnuc)
      ! Calculates the nucleation temperature for CO2, based on a simple super saturation criterion. JL 2011.

      use callkeys_mod, only: co2supsat

      implicit none

      real p, peff, tnuc
      real, parameter :: ptriple=518000.0

      peff=p/co2supsat

      if(peff.lt.ptriple) then
         tnuc = (-3167.8)/(log(.01*peff)-23.23) ! Fanale's formula
      else
         tnuc = 684.2-92.3*log(peff)+4.32*log(peff)**2 
         ! liquid-vapour transition (based on CRC handbook 2003 data)
      endif
      
      return

      end subroutine get_tnuc_co2
