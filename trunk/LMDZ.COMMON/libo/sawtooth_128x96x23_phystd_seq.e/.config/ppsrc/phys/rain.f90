










subroutine rain(ngrid,nlayer,nq,ptimestep,pplev,pplay,t,pdt,pq,pdq,d_t,dqrain,dqsrain,dqssnow,reevap_precip,rneb)


  use ioipsl_getin_p_mod, only: getin_p
  use watercommon_h, only: T_h2O_ice_liq,T_h2O_ice_clouds, RLVTT, RCPD, RCPV, RV, RVTMP2,Psat_water,Tsat_water,rhowater
  use radii_mod, only: h2o_cloudrad
  USE tracer_h, only: igcm_h2o_vap, igcm_h2o_ice
  use comcstfi_mod, only: g, r
  implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Calculates H2O precipitation using simplified microphysics.
!     
!     Authors
!     -------
!     Adapted from the LMDTERRE code by R. Wordsworth (2009)
!     Added rain vaporization in case of T>Tsat
!     Original author Z. X. Li (1993)
!     
!==================================================================

!     Arguments
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers
      real,intent(in) :: ptimestep    ! time interval
      real,intent(in) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      real,intent(in) :: pplay(ngrid,nlayer)   ! mid-layer pressure (Pa)
      real,intent(in) :: t(ngrid,nlayer) ! input temperature (K)
      real,intent(in) :: pdt(ngrid,nlayer) ! input tendency on temperature (K/s)      
      real,intent(in) :: pq(ngrid,nlayer,nq)  ! tracers (kg/kg)
      real,intent(in) :: pdq(ngrid,nlayer,nq) ! input tendency on tracers
      real,intent(out) :: d_t(ngrid,nlayer) ! temperature tendency (K/s)
      real,intent(out) :: dqrain(ngrid,nlayer,nq) ! tendency of H2O precipitation (kg/kg.s-1)
      real,intent(out) :: dqsrain(ngrid)  ! rain flux at the surface (kg.m-2.s-1)
      real,intent(out) :: dqssnow(ngrid)  ! snow flux at the surface (kg.m-2.s-1)
      real,intent(out) :: reevap_precip(ngrid)  ! re-evaporation flux of precipitation integrated over the atmospheric column (kg.m-2.s-1)
      real,intent(in) :: rneb(ngrid,nlayer) ! cloud fraction

      REAL zt(ngrid,nlayer)         ! working temperature (K)
      REAL ql(ngrid,nlayer)         ! liquid water (Kg/Kg)
      REAL q(ngrid,nlayer)          ! specific humidity (Kg/Kg)
      REAL d_q(ngrid,nlayer)        ! water vapor increment
      REAL d_ql(ngrid,nlayer)       ! liquid water / ice increment

!     Subroutine options
      REAL,PARAMETER :: seuil_neb=0.001  ! Nebulosity threshold

      INTEGER,save :: precip_scheme      ! id number for precipitaion scheme
!     for simple scheme  (precip_scheme=1)
      REAL,SAVE :: rainthreshold                ! Precipitation threshold in simple scheme
!     for sundquist scheme  (precip_scheme=2-3)
      REAL,SAVE :: cloud_sat                    ! Precipitation threshold in non simple scheme
      REAL,SAVE :: precip_timescale             ! Precipitation timescale
!     for Boucher scheme  (precip_scheme=4)
      REAL,SAVE :: Cboucher ! Precipitation constant in Boucher 95 scheme
      REAL,PARAMETER :: Kboucher=1.19E8
      REAL,SAVE :: c1
!$OMP THREADPRIVATE(precip_scheme,rainthreshold,cloud_sat,precip_timescale,Cboucher,c1)

      INTEGER,PARAMETER :: ninter=5

      logical,save :: evap_prec ! Does the rain evaporate?
!$OMP THREADPRIVATE(evap_prec)

!     for simple scheme
      real,parameter :: t_crit=218.0
      real lconvert

!     Local variables
      INTEGER i, k, n
      REAL zqs(ngrid,nlayer),Tsat(ngrid,nlayer), zdelta, zcor
      REAL zrfl(ngrid), zrfln(ngrid), zqev, zqevt

      REAL zoliq(ngrid)
      REAL zdz(ngrid),zrho(ngrid),ztot(ngrid), zrhol(ngrid)
      REAL zchau(ngrid),zfroi(ngrid),zfrac(ngrid),zneb(ngrid)

      real reffh2oliq(ngrid,nlayer),reffh2oice(ngrid,nlayer)
  
      real ttemp, ptemp, psat_tmp
      real tnext(ngrid,nlayer)

      real l2c(ngrid,nlayer)
      real dWtot


!     Indices of water vapour and water ice tracers
      INTEGER, SAVE :: i_vap=0  ! water vapour
      INTEGER, SAVE :: i_ice=0  ! water ice
!$OMP THREADPRIVATE(i_vap,i_ice)

      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)

!     Online functions
      REAL fallv, fall2v, zzz ! falling speed of ice crystals
      fallv (zzz) = 3.29 * ((zzz)**0.16)
      fall2v (zzz) =10.6 * ((zzz)**0.31)  !for use with radii


      IF (firstcall) THEN

         i_vap=igcm_h2o_vap
         i_ice=igcm_h2o_ice
        
         write(*,*) "rain: i_ice=",i_ice
         write(*,*) "      i_vap=",i_vap

         PRINT*, 'in rain.F, ninter=', ninter
         PRINT*, 'in rain.F, evap_prec=', evap_prec

         write(*,*) "Precipitation scheme to use?"
         precip_scheme=1 ! default value
         call getin_p("precip_scheme",precip_scheme)
         write(*,*) " precip_scheme = ",precip_scheme

         if (precip_scheme.eq.1) then
	    write(*,*) "rainthreshold in simple scheme?"
            rainthreshold=0. ! default value
            call getin_p("rainthreshold",rainthreshold)
            write(*,*) " rainthreshold = ",rainthreshold

         else if (precip_scheme.eq.2.or.precip_scheme.eq.3) then
	    write(*,*) "cloud water saturation level in non simple scheme?"
            cloud_sat=2.6e-4   ! default value
            call getin_p("cloud_sat",cloud_sat)
            write(*,*) " cloud_sat = ",cloud_sat
	    write(*,*) "precipitation timescale in non simple scheme?"
            precip_timescale=3600.  ! default value
            call getin_p("precip_timescale",precip_timescale)
            write(*,*) " precip_timescale = ",precip_timescale

         else if (precip_scheme.eq.4) then
	    write(*,*) "multiplicative constant in Boucher 95 precip scheme"
            Cboucher=1.   ! default value
            call getin_p("Cboucher",Cboucher)
            write(*,*) " Cboucher = ",Cboucher	
	    c1=1.00*1.097/rhowater*Cboucher*Kboucher  

	 endif

         write(*,*) "re-evaporate precipitations?"
         evap_prec=.true. ! default value
         call getin_p("evap_prec",evap_prec)
         write(*,*) " evap_prec = ",evap_prec

         firstcall = .false.
      ENDIF ! of IF (firstcall)

!     GCM -----> subroutine variables
      DO k = 1, nlayer
      DO i = 1, ngrid

	 zt(i,k)   = t(i,k)+pdt(i,k)*ptimestep ! a big fat bug was here
         q(i,k)    = pq(i,k,i_vap)+pdq(i,k,i_vap)*ptimestep
         ql(i,k)   = pq(i,k,i_ice)+pdq(i,k,i_ice)*ptimestep

         !q(i,k)    = pq(i,k,i_vap)!+pdq(i,k,i_vap)
         !ql(i,k)   = pq(i,k,i_ice)!+pdq(i,k,i_ice)

         if(q(i,k).lt.0.)then ! if this is not done, we don't conserve water
            q(i,k)=0.
         endif
         if(ql(i,k).lt.0.)then
            ql(i,k)=0.
         endif

      ENDDO
      ENDDO

!     Initialise the outputs
      d_t(1:ngrid,1:nlayer) = 0.0
      d_q(1:ngrid,1:nlayer) = 0.0
      d_ql(1:ngrid,1:nlayer) = 0.0
      zrfl(1:ngrid) = 0.0
      zrfln(1:ngrid) = 0.0

      ! calculate saturation mixing ratio
      DO k = 1, nlayer
         DO i = 1, ngrid
            ptemp = pplay(i,k)
	    call Psat_water(zt(i,k) ,ptemp,psat_tmp,zqs(i,k))
	    call Tsat_water(ptemp,Tsat(i,k))
         ENDDO
      ENDDO

      ! get column / layer conversion factor
      DO k = 1, nlayer
         DO i = 1, ngrid
            l2c(i,k)=(pplev(i,k)-pplev(i,k+1))/g
         ENDDO
      ENDDO

      ! Vertical loop (from top to bottom)
      ! We carry the rain with us and calculate that added by warm/cold precipitation
      ! processes and that subtracted by evaporation at each level.
      DO k = nlayer, 1, -1

         IF (evap_prec) THEN ! note no rneb dependence!
            DO i = 1, ngrid
               IF (zrfl(i) .GT.0.) THEN

                  if(zt(i,k).gt.Tsat(i,k))then
!!		     treat the case where all liquid water should boil
		     zqev=MIN((zt(i,k)-Tsat(i,k))*RCPD*l2c(i,k)/RLVTT/ptimestep,zrfl(i))
		     zrfl(i)=MAX(zrfl(i)-zqev,0.)
		     d_q(i,k)=zqev/l2c(i,k)*ptimestep
		     d_t(i,k) = - d_q(i,k) * RLVTT/RCPD
		  else
                     zqev = MAX (0.0, (zqs(i,k)-q(i,k)))*l2c(i,k)/ptimestep !there was a bug here
                     zqevt= 2.0e-5*(1.0-q(i,k)/zqs(i,k))    & !default was 2.e-5
                        *sqrt(zrfl(i))*l2c(i,k)/pplay(i,k)*zt(i,k)*R ! BC modif here
                     zqevt = MAX (zqevt, 0.0)
                     zqev  = MIN (zqev, zqevt)
                     zqev  = MAX (zqev, 0.0)
                     zrfln(i)= zrfl(i) - zqev
		     zrfln(i)= max(zrfln(i),0.0)

                     d_q(i,k) = - (zrfln(i)-zrfl(i))/l2c(i,k)*ptimestep
                     !d_t(i,k) = d_q(i,k) * RLVTT/RCPD!/(1.0+RVTMP2*q(i,k)) ! double BC modif here
                     d_t(i,k) = - d_q(i,k) * RLVTT/RCPD ! was bugged!
                     zrfl(i)  = zrfln(i)
		  end if
		     

               ENDIF ! of IF (zrfl(i) .GT.0.)
            ENDDO
         ENDIF ! of IF (evap_prec)

         zoliq(1:ngrid) = 0.0


         if(precip_scheme.eq.1)then

            DO i = 1, ngrid
               ttemp = zt(i,k)
               IF (ttemp .ge. T_h2O_ice_liq) THEN
                  lconvert=rainthreshold
               ELSEIF (ttemp .gt. t_crit) THEN
                  lconvert=rainthreshold*(1.- t_crit/ttemp)
                  lconvert=MAX(0.0,lconvert)              
               ELSE
                  lconvert=0.
               ENDIF


               IF (ql(i,k).gt.1.e-9) then
                  zneb(i)  = MAX(rneb(i,k), seuil_neb)
                  IF ((ql(i,k)/zneb(i)).gt.lconvert)THEN ! precipitate!
                     d_ql(i,k) = -MAX((ql(i,k)-lconvert*zneb(i)),0.0)
                     zrfl(i)   = zrfl(i) - d_ql(i,k)*l2c(i,k)/ptimestep
                  ENDIF
               ENDIF
            ENDDO

         elseif (precip_scheme.ge.2) then
          
           DO i = 1, ngrid
               IF (rneb(i,k).GT.0.0) THEN
                  zoliq(i) = ql(i,k)
                  zrho(i)  = pplay(i,k) / ( zt(i,k) * R )
                  zdz(i)   = (pplev(i,k)-pplev(i,k+1)) / (zrho(i)*g)
                  zfrac(i) = (zt(i,k)-T_h2O_ice_clouds) / (T_h2O_ice_liq-T_h2O_ice_clouds)
                  zfrac(i) = MAX(zfrac(i), 0.0)
                  zfrac(i) = MIN(zfrac(i), 1.0)
                  zneb(i)  = MAX(rneb(i,k), seuil_neb)
               ENDIF
           ENDDO

 !recalculate liquid water particle radii
	   call h2o_cloudrad(ngrid,nlayer,ql,reffh2oliq,reffh2oice)

           SELECT CASE(precip_scheme)
 !precip scheme from Sundquist 78
 	   CASE(2)

            DO n = 1, ninter
               DO i = 1, ngrid
                  IF (rneb(i,k).GT.0.0) THEN
                     ! this is the ONLY place where zneb, precip_timescale and cloud_sat are used

                     zchau(i) = (ptimestep/(FLOAT(ninter)*precip_timescale)) * zoliq(i)      &
                          * (1.0-EXP(-(zoliq(i)/zneb(i)/cloud_sat)**2)) * zfrac(i)
                     zrhol(i) = zrho(i) * zoliq(i) / zneb(i)
                     zfroi(i) = ptimestep/FLOAT(ninter)/zdz(i)*zoliq(i)    &
                          *fall2v(reffh2oice(i,k)) * (1.0-zfrac(i)) ! zfroi behaves oddly...
                     ztot(i)  = zchau(i) + zfroi(i)

                     IF (zneb(i).EQ.seuil_neb) ztot(i) = 0.0
                     ztot(i)  = MIN(MAX(ztot(i),0.0),zoliq(i))
                     zoliq(i) = MAX(zoliq(i)-ztot(i), 0.0)

                  ENDIF
               ENDDO
            ENDDO            

 !precip scheme modified from Sundquist 78 (in q**3)
	   CASE(3)	    
	    
            DO n = 1, ninter
               DO i = 1, ngrid
                  IF (rneb(i,k).GT.0.0) THEN
                     ! this is the ONLY place where zneb, precip_timescale and cloud_sat are used

                     zchau(i) = (ptimestep/(FLOAT(ninter)*precip_timescale*cloud_sat**2)) * (zoliq(i)/zneb(i))**3  
                     zrhol(i) = zrho(i) * zoliq(i) / zneb(i)
                     zfroi(i) = ptimestep/FLOAT(ninter)/zdz(i)*zoliq(i)    &
                          *fall2v(reffh2oice(i,k)) * (1.0-zfrac(i)) ! zfroi behaves oddly...
                     ztot(i)  = zchau(i) + zfroi(i)

                     IF (zneb(i).EQ.seuil_neb) ztot(i) = 0.0
                     ztot(i)  = MIN(MAX(ztot(i),0.0),zoliq(i))
                     zoliq(i) = MAX(zoliq(i)-ztot(i), 0.0)

                  ENDIF
               ENDDO
            ENDDO            

 !precip scheme modified from Boucher 95
	   CASE(4)

            DO n = 1, ninter
               DO i = 1, ngrid
                  IF (rneb(i,k).GT.0.0) THEN
                     ! this is the ONLY place where zneb and c1 are used

                     zchau(i) = ptimestep/FLOAT(ninter) *c1* zrho(i) &
                                    *(zoliq(i)/zneb(i))**2*reffh2oliq(i,k)*zneb(i)* zfrac(i) 
                     zrhol(i) = zrho(i) * zoliq(i) / zneb(i)
                     zfroi(i) = ptimestep/FLOAT(ninter)/zdz(i)*zoliq(i)    &
                          *fall2v(reffh2oice(i,k)) * (1.0-zfrac(i)) ! zfroi behaves oddly...
                     ztot(i)  = zchau(i) + zfroi(i)

                     IF (zneb(i).EQ.seuil_neb) ztot(i) = 0.0
                     ztot(i)  = MIN(MAX(ztot(i),0.0),zoliq(i))
                     zoliq(i) = MAX(zoliq(i)-ztot(i), 0.0)

                  ENDIF
               ENDDO
            ENDDO            

           END SELECT ! precip_scheme 

            ! Change in cloud density and surface H2O values
            DO i = 1, ngrid
               IF (rneb(i,k).GT.0.0) THEN
                  d_ql(i,k) = (zoliq(i) - ql(i,k))!/ptimestep
                  zrfl(i)   = zrfl(i)+ MAX(ql(i,k)-zoliq(i),0.0)*l2c(i,k)/ptimestep
               ENDIF
            ENDDO


         endif ! if precip_scheme=1

      ENDDO ! of DO k = nlayer, 1, -1

!     Rain or snow on the ground
      DO i = 1, ngrid
         if(zrfl(i).lt.0.0)then
            print*,'Droplets of negative rain are falling...'
            call abort
         endif
         IF (t(i,1) .LT. T_h2O_ice_liq) THEN
            dqssnow(i) = zrfl(i)
            dqsrain(i) = 0.0
         ELSE
            dqssnow(i) = 0.0
            dqsrain(i) = zrfl(i) ! liquid water = ice for now
         ENDIF
      ENDDO

!     now subroutine -----> GCM variables
      if (evap_prec) then
        dqrain(1:ngrid,1:nlayer,i_vap)=d_q(1:ngrid,1:nlayer)/ptimestep
        d_t(1:ngrid,1:nlayer)=d_t(1:ngrid,1:nlayer)/ptimestep
        do i=1,ngrid
           reevap_precip(i)=0.
           do k=1,nlayer
              reevap_precip(i)=reevap_precip(i)+dqrain(i,k,i_vap)*l2c(i,k)
           enddo
        enddo
      else
        dqrain(1:ngrid,1:nlayer,i_vap)=0.0
        d_t(1:ngrid,1:nlayer)=0.0
      endif
      dqrain(1:ngrid,1:nlayer,i_ice) = d_ql(1:ngrid,1:nlayer)/ptimestep

    end subroutine rain
