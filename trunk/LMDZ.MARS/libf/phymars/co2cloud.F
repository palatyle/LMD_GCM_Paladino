      MODULE co2cloud_mod

      IMPLICIT NONE 

      DOUBLE PRECISION,allocatable,save :: mem_Mccn_co2(:,:) ! Memory of CCN mass of H2O and dust used by CO2
      DOUBLE PRECISION,allocatable,save :: mem_Mh2o_co2(:,:) ! Memory of H2O mass integred into CO2 crystal
      DOUBLE PRECISION,allocatable,save :: mem_Nccn_co2(:,:) ! Memory of CCN number of H2O and dust used by CO2

      CONTAINS

      SUBROUTINE co2cloud(ngrid,nlay,ptimestep, 
     &                pplev,pplay,pdpsrf,pzlay,pt,pdt,
     &                pq,pdq,pdqcloudco2,pdtcloudco2,
     &                nq,tau,tauscaling,rdust,rice,riceco2,nuice,
     &                rsedcloudco2,rhocloudco2,
     &                rsedcloud,rhocloud,pzlev,pdqs_sedco2,
     &                pdu,pu,pcondicea)
      USE ioipsl_getincom, only: getin
      use dimradmars_mod, only: naerkind
      USE comcstfi_h, only: pi, g, cpp
      USE updaterad, only: updaterice_microco2, updaterice_micro,
     &                     updaterdust
      use conc_mod, only: mmean,rnew
      use tracer_mod, only: nqmx, igcm_co2, igcm_co2_ice,
     &     igcm_dust_mass, igcm_dust_number,igcm_h2o_ice,
     &     igcm_ccn_mass,igcm_ccn_number,
     &     igcm_ccnco2_mass, igcm_ccnco2_number,
     &     rho_dust, nuiceco2_sed, nuiceco2_ref,
     &     rho_ice_co2,r3n_q,rho_ice,nuice_sed
      USE newsedim_mod, ONLY: newsedim
      USE datafile_mod, ONLY: datadir
      USE improvedco2clouds_mod, ONLY: improvedco2clouds

      IMPLICIT NONE

      include "callkeys.h"
      include "microphys.h"

c=======================================================================
c CO2 clouds formation
c
c  There is a time loop specific to cloud formation 
c  due to timescales smaller than the GCM integration timestep.
c  microphysics subroutine is improvedCO2clouds.F
c  the microphysics time step is a fraction of the physical one
c  the integer imicroco2 must be set in callphys.def  
c
c  The co2 clouds tracers (co2_ice, ccn mass and concentration) are 
c  sedimented at each microtimestep. pdqs_sedco2 keeps track of the 
c  CO2 flux at the surface 
c
c  Authors: 09/2016 Joachim Audouard & Constantino Listowski 
c  Adaptation of the water ice clouds scheme (with specific microphysics)
c  of Montmessin, Navarro & al.
c
c 07/2017 J.Audouard
c Several logicals and integer must be set to .true. in callphys.def 
c if not, default values are .false.
c co2clouds=.true. call this routine
c co2useh2o=.true. allow the use of water ice particles as CCN for CO2
c meteo_flux=.true. supply meteoritic particles 
c CLFvaryingCO2=.true. allows a subgrid temperature distribution
c of amplitude spantCO2(=integer in callphys.def, typically 3)
c satindexco2=.true. allows the filtering out of the sub-grid T distribution
c                    if the GW saturates in the column. Based on Spiga et al
c                    2012
c                    An index is computed for the column, and the sub-grid T 
c                    distribution is applied if the index remains < 0.1
c                    setting to .false. applies the sub-grid T everywhere.
c                    default value is .true., only applies if 
c                    CLFvaryingCO2=.true. anyway.
c imicroco2=50 
c
c The subgrid Temperature distribution is modulated (0 or 1) by Spiga et
c al. (GRL 2012) Saturation Index to account for GW propagation or 
c dissipation upwards.
c
c 4D and column opacities are computed using Qext values at 1µm. 
c=======================================================================

c-----------------------------------------------------------------------
c   arguments:
c   -------------

      INTEGER, INTENT(IN) :: ngrid,nlay
      REAL, INTENT(IN) :: ptimestep            ! pas de temps physique (s)
      REAL, INTENT(IN) :: pplev(ngrid,nlay+1)   ! Inter-layer pressures (Pa)
      REAL, INTENT(IN) :: pplay(ngrid,nlay)     ! mid-layer pressures (Pa)
      REAL, INTENT(IN) :: pdpsrf(ngrid)         ! tendency on surface pressure
      REAL, INTENT(IN) :: pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL, INTENT(IN) :: pt(ngrid,nlay)        ! temperature at the middle of the layers (K)
      REAL, INTENT(IN) :: pdt(ngrid,nlay)       ! tendency on temperature from other parametrizations
      real, INTENT(IN) :: pq(ngrid,nlay,nq)     ! tracers (kg/kg)
      real, INTENT(IN) :: pdq(ngrid,nlay,nq)    ! tendencies before condensation  (kg/kg.s-1)
      real, intent(OUT) :: pdqcloudco2(ngrid,nlay,nq) ! tendency due to CO2 condensation (kg/kg.s-1)
      real, intent(OUT) :: pcondicea(ngrid,nlay)
      real, intent(OUT) :: pdtcloudco2(ngrid,nlay)    ! tendency on temperature due to latent heat
      INTEGER, INTENT(IN) :: nq                 ! number of tracers
      REAL, INTENT(IN) :: tau(ngrid,naerkind) ! Column dust optical depth at each point
      REAL, INTENT(IN) :: tauscaling(ngrid)   ! Convertion factor for dust amount
      REAL, INTENT(OUT) :: rdust(ngrid,nlay)   ! Dust geometric mean radius (m)
      real, intent(OUT) :: rice(ngrid,nlay)    ! Water Ice mass mean radius (m)
                                ! used for nucleation of CO2 on ice-coated ccns
      DOUBLE PRECISION, INTENT(out) :: riceco2(ngrid,nlay)    ! Ice mass mean radius (m)
                               ! (r_c in montmessin_2004)
      REAL, INTENT(IN) :: nuice(ngrid,nlay)   ! Estimated effective variance
                               !   of the size distribution
      real, intent(OUT) :: rsedcloudco2(ngrid,nlay) ! Cloud sedimentation radius
      real, intent(OUT) :: rhocloudco2(ngrid,nlay)  ! Cloud density (kg.m-3)
      real, intent(OUT) :: rsedcloud(ngrid,nlay) ! Water Cloud sedimentation radius
      real, intent(OUT) :: rhocloud(ngrid,nlay)  ! Water Cloud density (kg.m-3)
      real, intent(IN) :: pzlev(ngrid,nlay+1) ! altitude at the boundaries of the layers
      real, intent(OUT) :: pdqs_sedco2(ngrid) ! CO2 flux at the surface
      REAL, INTENT(IN) :: pdu(ngrid,nlay),pu(ngrid,nlay) !Zonal Wind: zu=pu+pdu*ptimestep

c   local:
c   ------
        
      ! for time loop
      INTEGER microstep  ! time subsampling step variable
      INTEGER, SAVE :: imicroco2     ! time subsampling for coupled water microphysics & sedimentation
      REAL, SAVE :: microtimestep ! integration timestep for coupled water microphysics & sedimentation
      
      ! tendency given by clouds (inside the micro loop)
      REAL subpdqcloudco2(ngrid,nlay,nq) ! cf. pdqcloud
      REAL subpdtcloudco2(ngrid,nlay)    ! cf. pdtcloud

      ! global tendency (clouds+physics)
      REAL sum_subpdq(ngrid,nlay,nq)      ! cf. pdqcloud
      REAL sum_subpdt(ngrid,nlay)         ! cf. pdtcloud
      real wq(ngrid,nlay+1)  !  ! displaced tracer mass (kg.m-2) during microtimestep because sedim (?/m-2)

      REAL satuco2(ngrid,nlay)  ! co2 satu ratio for output
      REAL zqsatco2(ngrid,nlay) ! saturation co2

      DOUBLE PRECISION rho_ice_co2T(ngrid,nlay) !T-dependant CO2 ice density
      DOUBLE PRECISION :: myT   ! temperature scalar for co2 density computation

      INTEGER iq,ig,l,i
      LOGICAL,SAVE :: firstcall=.true.
      DOUBLE PRECISION Nccnco2, Niceco2,Nco2,Qccnco2
      real :: beta ! for sedimentation

      real epaisseur (ngrid,nlay) ! Layer thickness (m)
      real masse (ngrid,nlay) ! Layer mass (kg.m-2)
      real ztsed(ngrid,nlay) ! tracers with real-time value in microtimeloop
      real zqsed(ngrid,nlay,nq)
      real zqsed0(ngrid,nlay,nq) !For sedimentation tendancy
      real subpdqsed(ngrid,nlay,nq)
      real sum_subpdqs_sedco2(ngrid) ! CO2 flux at the surface
     
!     What we need for Qext reading and tau computation : size distribution
      DOUBLE PRECISION vrat_cld ! Volume ratio
      DOUBLE PRECISION, SAVE :: rb_cldco2(nbinco2_cld+1) ! boundary values of each rad_cldco2 bin (m)
      DOUBLE PRECISION, PARAMETER :: rmin_cld = 1.e-9 ! Minimum radius (m)
      DOUBLE PRECISION, PARAMETER :: rmax_cld = 5.e-6 ! Maximum radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmin_cld =1.e-10! Minimum boundary radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmax_cld = 2.e-4 ! Maximum boundary radius (m)
      DOUBLE PRECISION dr_cld(nbinco2_cld) ! width of each rad_cldco2 bin (m)
      DOUBLE PRECISION vol_cld(nbinco2_cld) ! particle volume for each bin (m3)
      REAL, SAVE :: sigma_iceco2 ! Variance of the ice and CCN distributions
      logical :: file_ok !Qext file reading
      double precision :: radv(10000),Qextv1mic(10000)
      double precision, save :: Qext1bins(100)
      double precision :: Qtemp
      double precision :: ltemp1(10000),ltemp2(10000)
      integer :: nelem,lebon1,lebon2
      integer,parameter :: uQext=555
      DOUBLE PRECISION n_aer(nbinco2_cld),Rn,No,n_derf,dev2
      DOUBLE PRECISION Qext1bins2(ngrid,nlay)   
      DOUBLE PRECISION tau1mic(ngrid) !co2 ice column opacity at 1µm
     
        ! For sub grid T distribution

      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL :: zq(ngrid, nlay,nq)

      real :: rhocloudco2t(ngrid,nlay)  ! Cloud density (kg.m-3)

      DOUBLE PRECISION :: tcond(ngrid,nlay) !CO2 condensation temperature
      REAL ::  zqvap(ngrid,nlay)
      REAL ::  zqice(ngrid,nlay)
      REAL ::  spant,zdelt ! delta T for the temperature distribution
      REAL ::  pteff(ngrid, nlay)! effective temperature in the cloud,neb
      REAL ::  pqeff(ngrid, nlay, nq)! effective tracers quantities in the cloud
      REAL ::  co2cloudfrac(ngrid,nlay) ! cloud fraction
      REAL ::  mincloud ! min cloud frac
      DOUBLE PRECISION:: rho,zu,NN,gradT !For Saturation Index computation
      DOUBLE PRECISION :: SatIndex(ngrid,nlay),SatIndexmap(ngrid)

c D. BARDET : sensibility test
      REAL :: No_dust(ngrid,nlay)
      REAL :: Mo_dust(ngrid,nlay)
      
c      logical :: CLFvaryingCO2
c     ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
        if (nq.gt.nqmx) then
           write(*,*) 'stop in co2cloud (nq.gt.nqmx)!'
           write(*,*) 'nq=',nq,' nqmx=',nqmx
           stop
        endif
        write(*,*) "co2cloud.F: rho_ice_co2 = ",rho_ice_co2
        write(*,*) "co2cloud: igcm_co2=",igcm_co2
        write(*,*) "            igcm_co2_ice=",igcm_co2_ice
                
        write(*,*) "time subsampling for microphysic ?"
#ifdef MESOSCALE
        imicroco2 = 2
#else
        imicroco2 = 30
#endif
        call getin("imicroco2",imicroco2)
        write(*,*)"imicroco2 = ",imicroco2
        
        microtimestep = ptimestep/real(imicroco2)
        write(*,*)"Physical timestep is",ptimestep 
        write(*,*)"CO2 Microphysics timestep is",microtimestep 

c     Compute the size bins of the distribution of CO2 ice particles
c --> used for opacity calculations

c       rad_cldco2 is the primary radius grid used for microphysics computation.
c       The grid spacing is computed assuming a constant volume ratio
c       between two consecutive bins; i.e. vrat_cld.
c       vrat_cld is determined from the boundary values of the size grid: 
c       rmin_cld and rmax_cld.
c       The rb_cldco2 array contains the boundary values of each rad_cldco2 bin.
c       dr_cld is the width of each rad_cldco2 bin.
        sigma_iceco2 = sqrt(log(1.+nuiceco2_sed))
c       Volume ratio between two adjacent bins
   !     vrat_cld 
        vrat_cld = log(rmax_cld/rmin_cld) / float(nbinco2_cld-1) *3.
        vrat_cld = exp(vrat_cld)
        rb_cldco2(1)  = rbmin_cld
        rad_cldco2(1) = rmin_cld
        vol_cld(1) = 4./3. * dble(pi) * rmin_cld*rmin_cld*rmin_cld
        do i=1,nbinco2_cld-1
          rad_cldco2(i+1)  = rad_cldco2(i) * vrat_cld**(1./3.)
          vol_cld(i+1)  = vol_cld(i) * vrat_cld
        enddo
        do i=1,nbinco2_cld
          rb_cldco2(i+1)= ( (2.*vrat_cld) / (vrat_cld+1.) )**(1./3.) *
     &      rad_cldco2(i)
          dr_cld(i)  = rb_cldco2(i+1) - rb_cldco2(i)
        enddo
        rb_cldco2(nbinco2_cld+1) = rbmax_cld
        dr_cld(nbinco2_cld)   = rb_cldco2(nbinco2_cld+1) -
     &       rb_cldco2(nbinco2_cld)

c   read the Qext values 
        INQUIRE(FILE=TRIM(datadir)//
     &       '/optprop_co2ice_1mic.dat', EXIST=file_ok)
        IF (.not. file_ok) THEN 
           write(*,*) 'file optprop_co2ice_1mic.dat should be in '
     &          ,trim(datadir)
           STOP
        endif
!        open(newunit=uQext,file=trim(datadir)//
        open(unit=uQext,file=trim(datadir)//
     &       '/optprop_co2ice_1mic.dat'
     &       ,FORM='formatted')
        read(uQext,*) !skip 1 line
        do i=1,10000 
           read(uQext,'(E11.5)') radv(i)
        enddo
        read(uQext,*) !skip 1 line
        do i=1,10000 
           read(uQext,'(E11.5)') Qextv1mic(i)
        enddo
        close(uQext)
c   innterpol the Qext values 
        !rice_out=rad_cldco2
        do i=1,nbinco2_cld
           ltemp1=abs(radv(:)-rb_cldco2(i))
           ltemp2=abs(radv(:)-rb_cldco2(i+1))
           lebon1=minloc(ltemp1,DIM=1)
           lebon2=min(minloc(ltemp2,DIM=1),10000)
           nelem=lebon2-lebon1+1.
           Qtemp=0d0
           do l=0,nelem
              Qtemp=Qtemp+Qextv1mic(min(lebon1+l,10000)) !mean value in the interval
           enddo
           Qtemp=Qtemp/nelem
           Qext1bins(i)=Qtemp
        enddo
        Qext1bins(:)=Qext1bins(:)*rad_cldco2(:)*rad_cldco2(:)*pi
!     The actuall tau computation and output is performed in co2cloud.F

        print*,'--------------------------------------------'
        print*,'Microphysics co2: size bin-Qext information:'
        print*,'   i, rad_cldco2(i), Qext1bins(i)'
        do i=1,nbinco2_cld
          write(*,'(i3,3x,3(e13.6,4x))') i, rad_cldco2(i),
     &    Qext1bins(i)
        enddo
        print*,'--------------------------------------------'


        do i=1,nbinco2_cld+1
            rb_cldco2(i) = log(rb_cldco2(i))  !! we save that so that it is not computed  at each timestep and gridpoint
        enddo
        if (CLFvaryingCO2) then
          write(*,*) 
          write(*,*) "CLFvaryingCO2 is set to true is callphys.def"
          write(*,*) "The temperature field is enlarged to +/-",spantCO2 
          write(*,*) "for the CO2 microphysics "
        endif

        firstcall=.false.
      ENDIF                     ! of IF (firstcall)
c ===========================================================================   
c     Initialization
c ===========================================================================
      dev2 = 1. / ( sqrt(2.) * sigma_iceco2 )
      beta=0.85
      sum_subpdq(1:ngrid,1:nlay,1:nq) = 0
      sum_subpdt(1:ngrid,1:nlay)      = 0
      subpdqcloudco2(1:ngrid,1:nlay,1:nq) = 0
      subpdtcloudco2(1:ngrid,1:nlay)      = 0

      wq(:,:)=0
      ! default value if no ice
      rhocloudco2(1:ngrid,1:nlay) = rho_dust
      rhocloudco2t(1:ngrid,1:nlay) = rho_dust
      epaisseur(1:ngrid,1:nlay)=0
      masse(1:ngrid,1:nlay)=0

      zqsed0(1:ngrid,1:nlay,1:nq)=0
      sum_subpdqs_sedco2(1:ngrid)=0
      subpdqsed(1:ngrid,1:nlay,1:nq)=0
      
      do  l=1,nlay
        do ig=1, ngrid
          masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
          epaisseur(ig,l)= pzlev(ig,l+1) - pzlev(ig,l)
       enddo
      enddo

c ==========================================================================
c   0.  Representation of sub-grid water ice clouds
c ==========================================================================
      IF (CLFvaryingCO2) THEN

         spant=spantCO2         ! delta T for the temprature distribution
         mincloud=0.1           ! min co2cloudfrac when there is ice  
         pteff(:,:)=pt(:,:) 
         co2cloudfrac(:,:)=mincloud
         
c Tendencies
         DO l=1,nlay
            DO ig=1,ngrid
               zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep
            ENDDO
         ENDDO
         DO l=1,nlay
            DO ig=1,ngrid
               DO iq=1,nq
                  zq(ig,l,iq)=pq(ig,l,iq)+pdq(ig,l,iq)*ptimestep 
               ENDDO
            ENDDO
         ENDDO
         zqvap=zq(:,:,igcm_co2)
         zqice=zq(:,:,igcm_co2_ice)


         call WRITEDIAGFI(ngrid,"co2cloud_pzlev","pzlev","km",3,
     &        pzlev)
         call WRITEDIAGFI(ngrid,"co2cloud_pzlay","pzlay","km",3,
     &        pzlay)
         call WRITEDIAGFI(ngrid,"co2cloud_pplay","pplay","Pa",3,
     &        pplay)

         if (satindexco2) then !logical in callphys.def
           DO l=12,26
             ! layers 12 --> 26 ~ 12->85 km
             DO ig=1,ngrid
             ! compute N^2 static stability
             gradT=(zt(ig,l+1)-zt(ig,l))/(pzlev(ig,l+1)-pzlev(ig,l))
             NN=sqrt(g/zt(iq,l)*(max(gradT,-g/cpp)+g/cpp))
             ! compute absolute value of zonal wind field 
             zu=abs(pu(ig,l)  + pdu(ig,l)*ptimestep)
             ! compute background density
             rho=pplay(ig,l)/(rnew(ig,l)*zt(ig,l))
             !saturation index
             SatIndex(ig,l)=sqrt(7.5e-7*150.e3/(2.*pi)*NN/
     &                                               (rho*zu*zu*zu))
             ENDDO
           ENDDO
            !Then compute Satindex map 
            ! layers 12 --> 26 ~ 12->85 km
           DO ig=1,ngrid
             SatIndexmap(ig)=maxval(SatIndex(ig,12:26))
           ENDDO

           call WRITEDIAGFI(ngrid,"SatIndexmap","SatIndexmap","km",2,
     &         SatIndexmap)
         else
           do ig=1,ngrid
             SatIndexmap(ig)=0.05 !maxval(SatIndex(ig,12:26))
           enddo
         endif ! of if (satindexco2)

!Modulate the DeltaT by GW propagation index : 
         ! Saturation index S in Spiga 2012 paper 
         !Assuming like in the paper, 
         !GW phase speed (stationary waves) c=0 m.s-1
         !lambdaH =150 km
         !Fo=7.5e-7 J.m-3
        
         CALL tcondco2(ngrid,nlay,pplay,zqvap,tcond)
         zdelt=spant  
         DO ig=1,ngrid
              
           IF (SatIndexmap(ig) .le. 0.1) THEN
             DO l=1,nlay-1
         
               IF (tcond(ig,l) .ge. (zt(ig,l)+zdelt) 
     &             .or. tcond(ig,l) .le. 0 ) THEN !The entire fraction is saturated
                  pteff(ig,l)=zt(ig,l)
                  co2cloudfrac(ig,l)=1.
               ELSE IF (tcond(ig,l) .le. (zt(ig,l)-zdelt)) THEN ! No saturation at all
                  pteff(ig,l)=zt(ig,l)-zdelt
                  co2cloudfrac(ig,l)=mincloud
               ELSE
                  co2cloudfrac(ig,l)=(tcond(ig,l)-zt(ig,l)+zdelt)/
     &                 (2.0*zdelt)
                  pteff(ig,l)=(tcond(ig,l)+zt(ig,l)-zdelt)/2. !Mean temperature of the cloud fraction
               END IF           !ig if (tcond(ig,l) ...
               pteff(ig,l)=pteff(ig,l)-pdt(ig,l)*ptimestep
               IF (co2cloudfrac(ig,l).le. mincloud) THEN
                  co2cloudfrac(ig,l)=mincloud
               ELSE IF (co2cloudfrac(ig,l).gt. 1) THEN
                  co2cloudfrac(ig,l)=1.
               END IF
             ENDDO
           ELSE
! SatIndex not favorable for GW : leave pt untouched
             pteff(ig,l)=pt(ig,l)
             co2cloudfrac(ig,l)=mincloud
           ENDIF                 ! of if(SatIndexmap...
         ENDDO ! of DO ig=1,ngrid
! Totalcloud frac of the column missing here
c
c No sub-grid cloud representation (CLFvarying=false)
      ELSE
         DO l=1,nlay
            DO ig=1,ngrid
               pteff(ig,l)=pt(ig,l)
            END DO
         END DO
      END IF                    ! end if (CLFvaryingco2)
c =============================================================================
c microtimestep timeloop for microphysics:
c 0.Stepped entry for tendancies
c 1.Compute sedimentation and update tendancies
c 2.Call co2clouds microphysics
c 3.Update tendancies 
c =============================================================================
      DO microstep=1,imicroco2
c Temperature tendency subpdt
        ! If imicro=1 subpdt is the same as pdt
        DO l=1,nlay
          DO ig=1,ngrid
               sum_subpdt(ig,l) = sum_subpdt(ig,l)
     &              + pdt(ig,l) ! At each micro timestep we add pdt in order to have a stepped entry
               sum_subpdq(ig,l,igcm_dust_mass) = 
     &              sum_subpdq(ig,l,igcm_dust_mass)
     &              + pdq(ig,l,igcm_dust_mass)
               sum_subpdq(ig,l,igcm_dust_number) = 
     &              sum_subpdq(ig,l,igcm_dust_number)
     &              + pdq(ig,l,igcm_dust_number)

               sum_subpdq(ig,l,igcm_ccnco2_mass) = 
     &              sum_subpdq(ig,l,igcm_ccnco2_mass)
     &              + pdq(ig,l,igcm_ccnco2_mass)
               sum_subpdq(ig,l,igcm_ccnco2_number) = 
     &              sum_subpdq(ig,l,igcm_ccnco2_number)
     &              + pdq(ig,l,igcm_ccnco2_number)

               sum_subpdq(ig,l,igcm_co2_ice) = 
     &              sum_subpdq(ig,l,igcm_co2_ice)
     &              + pdq(ig,l,igcm_co2_ice)
               sum_subpdq(ig,l,igcm_co2) = 
     &              sum_subpdq(ig,l,igcm_co2)
     &              + pdq(ig,l,igcm_co2)
c D.BARDET :
               if (co2useh2o) then
               sum_subpdq(ig,l,igcm_h2o_ice) = 
     &              sum_subpdq(ig,l,igcm_h2o_ice)
     &              + pdq(ig,l,igcm_h2o_ice)

               sum_subpdq(ig,l,igcm_ccn_mass) = 
     &              sum_subpdq(ig,l,igcm_ccn_mass)
     &              + pdq(ig,l,igcm_ccn_mass)
               sum_subpdq(ig,l,igcm_ccn_number) = 
     &              sum_subpdq(ig,l,igcm_ccn_number)
     &              + pdq(ig,l,igcm_ccn_number)
                endif
          ENDDO
        ENDDO
c Effective tracers quantities in the cloud fraction
        IF (CLFvaryingCO2) THEN     
            pqeff(:,:,:)=pq(:,:,:) ! prevent from buggs (A. Pottier)
            pqeff(:,:,igcm_ccnco2_mass) =pq(:,:,igcm_ccnco2_mass)/
     &           co2cloudfrac(:,:)
            pqeff(:,:,igcm_ccnco2_number)= 
     &           pq(:,:,igcm_ccnco2_number)/co2cloudfrac(:,:)
            pqeff(:,:,igcm_co2_ice)= pq(:,:,igcm_co2_ice)/
     &           co2cloudfrac(:,:)
        ELSE
            pqeff(:,:,:)=pq(:,:,:)
        END IF  
       
c ========================================================================
c 1.SEDIMENTATION : update tracers, compute parameters,
c   call to sedimentation routine, update tendancies
c ========================================================================
        IF (sedimentation) THEN
        
        DO l=1, nlay
          DO ig=1,ngrid             
             ztsed(ig,l)=pteff(ig,l)
     &            +sum_subpdt(ig,l)*microtimestep
             zqsed(ig,l,:)=pqeff(ig,l,:)
     &            +sum_subpdq(ig,l,:)*microtimestep
             rho_ice_co2T(ig,l)=1000.*(1.72391-2.53e-4*
     &            ztsed(ig,l)-2.87e-6*
     &            ztsed(ig,l)*ztsed(ig,l))
             
             rho_ice_co2=rho_ice_co2T(ig,l)
             Niceco2=max(zqsed(ig,l,igcm_co2_ice),1.e-30)
             Nccnco2=max(zqsed(ig,l,igcm_ccnco2_number),
     &            1.e-30)
             Qccnco2=max(zqsed(ig,l,igcm_ccnco2_mass),
     &            1.e-30)
             call updaterice_microco2(Niceco2,
     &            Qccnco2,Nccnco2,
     &            tauscaling(ig),riceco2(ig,l),rhocloudco2t(ig,l))
             if (Niceco2 .le. 1.e-25 
     &            .or. Nccnco2*tauscaling(ig) .le. 1) THEN
                riceco2(ig,l)=1.e-9
             endif
             rhocloudco2t(ig,l)=min(max(rhocloudco2t(ig,l)
     &            ,rho_ice_co2),rho_dust)
             rsedcloudco2(ig,l)=max(riceco2(ig,l)*
     &            (1.+nuiceco2_sed)*(1.+nuiceco2_sed)*(1.+nuiceco2_sed),
     &            riceco2(ig,l))
          ENDDO
        ENDDO
! Gravitational sedimentation       
        zqsed0(:,:,igcm_co2_ice)=zqsed(:,:,igcm_co2_ice)
        zqsed0(:,:,igcm_ccnco2_mass)=zqsed(:,:,igcm_ccnco2_mass)
        zqsed0(:,:,igcm_ccnco2_number)=zqsed(:,:,igcm_ccnco2_number) 
! We save actualized tracer values to compute sedimentation tendancies
        call newsedim(ngrid,nlay,ngrid*nlay,ngrid*nlay,
     &     microtimestep,pplev,masse,epaisseur,ztsed,
     &     rsedcloudco2,rhocloudco2t,
     &     zqsed(:,:,igcm_co2_ice),wq,beta) !  3 traceurs

! sedim at the surface of co2 ice : keep track of it for physiq_mod
        do ig=1,ngrid 
          sum_subpdqs_sedco2(ig)=
     &         sum_subpdqs_sedco2(ig)+ wq(ig,1)/microtimestep
        end do

        call newsedim(ngrid,nlay,ngrid*nlay,ngrid*nlay,
     &     microtimestep,pplev,masse,epaisseur,ztsed,
     &     rsedcloudco2,rhocloudco2t,
     &     zqsed(:,:,igcm_ccnco2_mass),wq,beta) 

        call newsedim(ngrid,nlay,ngrid*nlay,ngrid*nlay,
     &     microtimestep,pplev,masse,epaisseur,ztsed,
     &     rsedcloudco2,rhocloudco2t,
     &     zqsed(:,:,igcm_ccnco2_number),wq,beta) 

        DO l = 1, nlay            !Compute tendencies
          DO ig=1,ngrid
            subpdqsed(ig,l,igcm_ccnco2_mass)=
     &           (zqsed(ig,l,igcm_ccnco2_mass)-
     &           zqsed0(ig,l,igcm_ccnco2_mass))/microtimestep
            subpdqsed(ig,l,igcm_ccnco2_number)=
     &           (zqsed(ig,l,igcm_ccnco2_number)-
     &           zqsed0(ig,l,igcm_ccnco2_number))/microtimestep
            subpdqsed(ig,l,igcm_co2_ice)=
     &           (zqsed(ig,l,igcm_co2_ice)-
     &           zqsed0(ig,l,igcm_co2_ice))/microtimestep
          ENDDO
        ENDDO
! update subtimestep tendencies with sedimentation input
        DO l=1,nlay
         DO ig=1,ngrid
            sum_subpdq(ig,l,igcm_ccnco2_mass) =
     &           sum_subpdq(ig,l,igcm_ccnco2_mass)
     &           +subpdqsed(ig,l,igcm_ccnco2_mass)
            sum_subpdq(ig,l,igcm_ccnco2_number) =
     &           sum_subpdq(ig,l,igcm_ccnco2_number)
     &           +subpdqsed(ig,l,igcm_ccnco2_number)
            sum_subpdq(ig,l,igcm_co2_ice) =
     &           sum_subpdq(ig,l,igcm_co2_ice)
     &           +subpdqsed(ig,l,igcm_co2_ice)
         ENDDO
        ENDDO
        
        END IF !(end if sedimentation)
        
c ==============================================================================
c      2.  Main call to the cloud schemes:
c ==============================================================================
        CALL improvedco2clouds(ngrid,nlay,microtimestep,
     &     pplay,pplev,pteff,sum_subpdt,
     &     pqeff,sum_subpdq,subpdqcloudco2,subpdtcloudco2,
     &     nq,tauscaling,mem_Mccn_co2,mem_Mh2o_co2,mem_Nccn_co2,
     &     No_dust,Mo_dust)
c D. BARDET: sensibility test
c      call WRITEDIAGFI(ngrid,"No_dust","Nombre particules de poussiere"
c     &        ,"part/kg",3,No_dust)
c      call WRITEDIAGFI(ngrid,"Mo_dust","Masse particules de poussiere"
c     &        ,"kg/kg ",3,Mo_dust)
c ==============================================================================
c      3.  Updating tendencies after cloud scheme:
c ==============================================================================
        DO l=1,nlay
          DO ig=1,ngrid
               sum_subpdt(ig,l) =
     &              sum_subpdt(ig,l) + subpdtcloudco2(ig,l)

               sum_subpdq(ig,l,igcm_dust_mass) =
     &              sum_subpdq(ig,l,igcm_dust_mass)
     &              + subpdqcloudco2(ig,l,igcm_dust_mass)
               sum_subpdq(ig,l,igcm_dust_number) =
     &              sum_subpdq(ig,l,igcm_dust_number)
     &              + subpdqcloudco2(ig,l,igcm_dust_number)

               sum_subpdq(ig,l,igcm_ccnco2_mass) =
     &              sum_subpdq(ig,l,igcm_ccnco2_mass)
     &              + subpdqcloudco2(ig,l,igcm_ccnco2_mass)
               sum_subpdq(ig,l,igcm_ccnco2_number) =
     &              sum_subpdq(ig,l,igcm_ccnco2_number)
     &              + subpdqcloudco2(ig,l,igcm_ccnco2_number)

               sum_subpdq(ig,l,igcm_co2_ice) =
     &              sum_subpdq(ig,l,igcm_co2_ice)
     &              + subpdqcloudco2(ig,l,igcm_co2_ice)
               sum_subpdq(ig,l,igcm_co2) =
     &              sum_subpdq(ig,l,igcm_co2)
     &              + subpdqcloudco2(ig,l,igcm_co2)
c D.BARDET :
               if (co2useh2o) then
               sum_subpdq(ig,l,igcm_h2o_ice) =
     &              sum_subpdq(ig,l,igcm_h2o_ice)
     &              + subpdqcloudco2(ig,l,igcm_h2o_ice)

               sum_subpdq(ig,l,igcm_ccn_mass) =
     &              sum_subpdq(ig,l,igcm_ccn_mass)
     &              + subpdqcloudco2(ig,l,igcm_ccn_mass)
               sum_subpdq(ig,l,igcm_ccn_number) =
     &              sum_subpdq(ig,l,igcm_ccn_number)
     &              + subpdqcloudco2(ig,l,igcm_ccn_number)
                endif
          ENDDO
        ENDDO
      ENDDO                     ! of DO microstep=1,imicro
      
c------------------------------------------------
c   Compute final tendencies after time loop:
c------------------------------------------------
c Condensation/sublimation tendency after clouds scheme (to replace
c zcondicea in newcondens.F

      DO l=nlay, 1, -1
        DO ig = 1, ngrid
           pcondicea(ig,l) = sum_subpdq(ig,l,igcm_co2_ice)
     &                       /real(imicroco2)
        ENDDO
      ENDDO

c CO2 flux at surface (kg.m-2.s-1)
      do ig=1,ngrid 
         pdqs_sedco2(ig)=sum_subpdqs_sedco2(ig)/real(imicroco2)
      enddo
c Temperature tendency pdtcloud
      DO l=1,nlay
        DO ig=1,ngrid
             pdtcloudco2(ig,l) =
     &         sum_subpdt(ig,l)/real(imicroco2)-pdt(ig,l)
        ENDDO
      ENDDO
c Tracers tendencies pdqcloud
      DO l=1,nlay
        DO ig=1,ngrid        
             pdqcloudco2(ig,l,igcm_co2_ice) = 
     &            sum_subpdq(ig,l,igcm_co2_ice)/real(imicroco2)
     &            - pdq(ig,l,igcm_co2_ice)
             pdqcloudco2(ig,l,igcm_co2) = 
     &            sum_subpdq(ig,l,igcm_co2)/real(imicroco2)
     &            - pdq(ig,l,igcm_co2)
c D.BARDET :
             if (co2useh2o) then
             pdqcloudco2(ig,l,igcm_h2o_ice) = 
     &            sum_subpdq(ig,l,igcm_h2o_ice)/real(imicroco2)
     &            - pdq(ig,l,igcm_h2o_ice)

             pdqcloudco2(ig,l,igcm_ccn_mass) = 
     &            sum_subpdq(ig,l,igcm_ccn_mass)/real(imicroco2)
     &            - pdq(ig,l,igcm_ccn_mass)

             pdqcloudco2(ig,l,igcm_ccn_number) = 
     &            sum_subpdq(ig,l,igcm_ccn_number)/real(imicroco2)
     &            - pdq(ig,l,igcm_ccn_number)
             endif
        
             pdqcloudco2(ig,l,igcm_ccnco2_mass) = 
     &            sum_subpdq(ig,l,igcm_ccnco2_mass)/real(imicroco2)
     &            - pdq(ig,l,igcm_ccnco2_mass)
        
             pdqcloudco2(ig,l,igcm_ccnco2_number) = 
     &            sum_subpdq(ig,l,igcm_ccnco2_number)/real(imicroco2)
     &            - pdq(ig,l,igcm_ccnco2_number)

             pdqcloudco2(ig,l,igcm_dust_mass) = 
     &            sum_subpdq(ig,l,igcm_dust_mass)/real(imicroco2)
     &            - pdq(ig,l,igcm_dust_mass)

             pdqcloudco2(ig,l,igcm_dust_number) = 
     &            sum_subpdq(ig,l,igcm_dust_number)/real(imicroco2)
     &            - pdq(ig,l,igcm_dust_number)
        ENDDO
      ENDDO
c Due to stepped entry, other processes tendencies can add up to negative values
c Therefore, enforce positive values and conserve mass
      DO l=1,nlay
        DO ig=1,ngrid
             IF ((pqeff(ig,l,igcm_ccnco2_number) + 
     &            ptimestep* (pdq(ig,l,igcm_ccnco2_number) + 
     &            pdqcloudco2(ig,l,igcm_ccnco2_number))
     &            .lt. 1.)
     &            .or. (pqeff(ig,l,igcm_ccnco2_mass) + 
     &            ptimestep* (pdq(ig,l,igcm_ccnco2_mass) + 
     &            pdqcloudco2(ig,l,igcm_ccnco2_mass))
     &            .lt. 1.e-20)) THEN
                pdqcloudco2(ig,l,igcm_ccnco2_number) =
     &               - pqeff(ig,l,igcm_ccnco2_number)/ptimestep 
     &               - pdq(ig,l,igcm_ccnco2_number)+1.
                
                pdqcloudco2(ig,l,igcm_dust_number) =  
     &               -pdqcloudco2(ig,l,igcm_ccnco2_number)

                pdqcloudco2(ig,l,igcm_ccnco2_mass) =
     &               - pqeff(ig,l,igcm_ccnco2_mass)/ptimestep
     &               - pdq(ig,l,igcm_ccnco2_mass)+1.e-20

                pdqcloudco2(ig,l,igcm_dust_mass) = 
     &               -pdqcloudco2(ig,l,igcm_ccnco2_mass)
             ENDIF
        ENDDO
      ENDDO
      DO l=1,nlay
        DO ig=1,ngrid
             IF ( (pqeff(ig,l,igcm_dust_number) + 
     &            ptimestep* (pdq(ig,l,igcm_dust_number) + 
     &            pdqcloudco2(ig,l,igcm_dust_number)) .le. 1.)
     &            .or. (pqeff(ig,l,igcm_dust_mass)+ 
     &            ptimestep* (pdq(ig,l,igcm_dust_mass) + 
     &            pdqcloudco2(ig,l,igcm_dust_mass))
     &            .le. 1.e-20)) then                  
                pdqcloudco2(ig,l,igcm_dust_number) =
     &               - pqeff(ig,l,igcm_dust_number)/ptimestep 
     &               - pdq(ig,l,igcm_dust_number)+1.

                pdqcloudco2(ig,l,igcm_ccnco2_number) =  
     &               -pdqcloudco2(ig,l,igcm_dust_number)

                pdqcloudco2(ig,l,igcm_dust_mass) =
     &               - pqeff(ig,l,igcm_dust_mass)/ptimestep
     &               - pdq(ig,l,igcm_dust_mass) +1.e-20

                pdqcloudco2(ig,l,igcm_ccnco2_mass) = 
     &               -pdqcloudco2(ig,l,igcm_dust_mass)
             ENDIF
        ENDDO
      ENDDO
! pq+ptime*(pdq+pdqc)=1 ! pdqc=1-pq/ptime-pdq      
      DO l=1,nlay
        DO ig=1,ngrid
             IF (pqeff(ig,l,igcm_co2_ice) + ptimestep*
     &       (pdq(ig,l,igcm_co2_ice) + pdqcloudco2(ig,l,igcm_co2_ice)) 
     &       .lt. 1.e-15) THEN
           pdqcloudco2(ig,l,igcm_co2_ice) = 
     &     - pqeff(ig,l,igcm_co2_ice)/ptimestep-pdq(ig,l,igcm_co2_ice)
           pdqcloudco2(ig,l,igcm_co2) = -pdqcloudco2(ig,l,igcm_co2_ice)
          ENDIF    
          IF (pqeff(ig,l,igcm_co2) + ptimestep*
     &       (pdq(ig,l,igcm_co2) + pdqcloudco2(ig,l,igcm_co2)) 
     &       .lt. 0.1) THEN
           pdqcloudco2(ig,l,igcm_co2) = 
     &     - pqeff(ig,l,igcm_co2)/ptimestep - pdq(ig,l,igcm_co2)
           pdqcloudco2(ig,l,igcm_co2_ice)= -pdqcloudco2(ig,l,igcm_co2)
          ENDIF
        ENDDO
      ENDDO

c Update clouds parameters values in the cloud fraction (for output)
      DO l=1, nlay
        DO ig=1,ngrid

              Niceco2=pqeff(ig,l,igcm_co2_ice) +                   
     &             (pdq(ig,l,igcm_co2_ice) + 
     &             pdqcloudco2(ig,l,igcm_co2_ice))*ptimestep
              Nco2=pqeff(ig,l,igcm_co2) +                   
     &             (pdq(ig,l,igcm_co2) + 
     &             pdqcloudco2(ig,l,igcm_co2))*ptimestep
              Nccnco2=max((pqeff(ig,l,igcm_ccnco2_number) +                 
     &             (pdq(ig,l,igcm_ccnco2_number) +               
     &             pdqcloudco2(ig,l,igcm_ccnco2_number))*ptimestep)
     &             ,1.e-30)
              Qccnco2=max((pqeff(ig,l,igcm_ccnco2_mass) +                 
     &             (pdq(ig,l,igcm_ccnco2_mass) +               
     &             pdqcloudco2(ig,l,igcm_ccnco2_mass))*ptimestep)
     &             ,1.e-30)
              
              myT=pteff(ig,l)+(pdt(ig,l)+pdtcloudco2(ig,l))*ptimestep
              rho_ice_co2T(ig,l)=1000.*(1.72391-2.53e-4* 
     &             myT-2.87e-6* myT* myT)
              rho_ice_co2=rho_ice_co2T(ig,l)
c             rho_ice_co2 is shared by tracer_mod and used in updaterice
c     Compute particle size
              call updaterice_microco2(Niceco2,
     &             Qccnco2,Nccnco2,
     &             tauscaling(ig),riceco2(ig,l),rhocloudco2(ig,l))
              
              if ( (Niceco2 .le. 1.e-25 .or. 
     &             Nccnco2*tauscaling(ig) .le. 1.) )THEN
                 riceco2(ig,l)=0.
                 Qext1bins2(ig,l)=0.
              else
c     Compute opacities
                No=Nccnco2*tauscaling(ig)
                Rn=-dlog(riceco2(ig,l))
                n_derf = erf( (rb_cldco2(1)+Rn) *dev2)
                Qext1bins2(ig,l)=0.
                do i = 1, nbinco2_cld 
                 n_aer(i) = -0.5 * No * n_derf !! this ith previously computed
                 n_derf = erf((rb_cldco2(i+1)+Rn) *dev2)
                 n_aer(i) = n_aer(i) + 0.5 * No * n_derf
                 Qext1bins2(ig,l)=Qext1bins2(ig,l)+Qext1bins(i)*n_aer(i)
                enddo
              endif
     
c D.BARDET : update rice water only if co2 use h2o ice as CCN
          if (co2useh2o) then
             call updaterice_micro(
     &       pqeff(ig,l,igcm_h2o_ice) +                    ! ice mass
     &       (pdq(ig,l,igcm_h2o_ice) +                     ! ice mass
     &       pdqcloudco2(ig,l,igcm_h2o_ice))*ptimestep,    ! ice mass
     &       pqeff(ig,l,igcm_ccn_mass) +                   ! ccn mass
     &       (pdq(ig,l,igcm_ccn_mass) +                    ! ccn mass
     &       pdqcloudco2(ig,l,igcm_ccn_mass))*ptimestep,   ! ccn mass
     &       pqeff(ig,l,igcm_ccn_number) +                 ! ccn number
     &       (pdq(ig,l,igcm_ccn_number) +                  ! ccn number
     &       pdqcloudco2(ig,l,igcm_ccn_number))*ptimestep, ! ccn number
     &       tauscaling(ig),rice(ig,l),rhocloud(ig,l))
          endif

          call updaterdust(
     &    pqeff(ig,l,igcm_dust_mass) +                   ! dust mass
     &   (pdq(ig,l,igcm_dust_mass) +                     ! dust mass
     &    pdqcloudco2(ig,l,igcm_dust_mass))*ptimestep,   ! dust mass
     &    pqeff(ig,l,igcm_dust_number) +                 ! dust number
     &   (pdq(ig,l,igcm_dust_number) +                   ! dust number
     &    pdqcloudco2(ig,l,igcm_dust_number))*ptimestep, ! dust number
     &    rdust(ig,l))

        ENDDO
      ENDDO 
     
c ======================================================================     
c A correction if a lot of subliming CO2 fills the 1st layer FF04/2005
c Then that should not affect the ice particle radius
c ======================================================================         
      do ig=1,ngrid
        if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,2)))then
             if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,3)))
     &       riceco2(ig,2)=riceco2(ig,3) 
             riceco2(ig,1)=riceco2(ig,2)
        endif
      end do
       
      DO l=1,nlay
         DO ig=1,ngrid
           rsedcloud(ig,l)=max(rice(ig,l)*
     &                 (1.+nuice_sed)*(1.+nuice_sed)*(1.+nuice_sed),
     &                    rdust(ig,l))
!          rsedcloud(ig,l)=min(rsedcloud(ig,l),1.e-4)
         ENDDO
      ENDDO
       
      DO l=1,nlay
         DO ig=1,ngrid
           rsedcloudco2(ig,l)=max(riceco2(ig,l)*
     &           (1.+nuiceco2_sed)*(1.+nuiceco2_sed)*(1.+nuiceco2_sed),
     &                    rdust(ig,l))
c          rsedcloudco2(ig,l)=min(rsedcloudco2(ig,l),1.e-5)
         ENDDO
      ENDDO
       
      call co2sat(ngrid*nlay,pteff+(pdt+pdtcloudco2)*ptimestep
     &      ,pplay,zqsatco2)
      do l=1,nlay
        do ig=1,ngrid 
             satuco2(ig,l) = (pqeff(ig,l,igcm_co2)  + 
     &            (pdq(ig,l,igcm_co2) +
     &            pdqcloudco2(ig,l,igcm_co2))*ptimestep)*
     &            (mmean(ig,l)/44.01)*pplay(ig,l)/zqsatco2(ig,l)
        enddo
      enddo
! Everything modified by CO2 microphysics must be wrt co2cloudfrac
      IF (CLFvaryingCO2) THEN
        DO l=1,nlay
          DO ig=1,ngrid

            pdqcloudco2(ig,l,igcm_ccnco2_mass)=
     &        pdqcloudco2(ig,l,igcm_ccnco2_mass)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_ccnco2_number)=
     &        pdqcloudco2(ig,l,igcm_ccnco2_number)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_dust_mass)=
     &        pdqcloudco2(ig,l,igcm_dust_mass)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_dust_number)=
     &        pdqcloudco2(ig,l,igcm_dust_number)*co2cloudfrac(ig,l)
c D.BARDET
            if (co2useh2o) then
            pdqcloudco2(ig,l,igcm_h2o_ice)=
     &        pdqcloudco2(ig,l,igcm_h2o_ice)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_ccn_mass)=
     &        pdqcloudco2(ig,l,igcm_ccn_mass)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_ccn_number)=
     &        pdqcloudco2(ig,l,igcm_ccn_number)*co2cloudfrac(ig,l)            
            endif

            pdqcloudco2(ig,l,igcm_co2_ice)=
     &        pdqcloudco2(ig,l,igcm_co2_ice)*co2cloudfrac(ig,l)

            pdqcloudco2(ig,l,igcm_co2)=
     &        pdqcloudco2(ig,l,igcm_co2)*co2cloudfrac(ig,l)

            pdtcloudco2(ig,l)=pdtcloudco2(ig,l)*co2cloudfrac(ig,l)

            Qext1bins2(ig,l)=Qext1bins2(ig,l)*co2cloudfrac(ig,l)
          ENDDO
        ENDDO   
      ENDIF
! opacity in mesh ig is the sum over l of Qext1bins2: Is this true ?
      tau1mic(:)=0.
      do l=1,nlay
        do ig=1,ngrid 
          tau1mic(ig)=tau1mic(ig)+Qext1bins2(ig,l)
        enddo
      enddo
!Outputs:
      call WRITEDIAGFI(ngrid,"SatIndex","SatIndex"," ",3,
     &        SatIndex)
      call WRITEDIAGFI(ngrid,"satuco2","vap in satu","kg/kg",3,
     &        satuco2)
      call WRITEdiagfi(ngrid,"riceco2","ice radius","m"
     &        ,3,riceco2)         
      call WRITEdiagfi(ngrid,"co2cloudfrac","co2 cloud fraction"
     &        ," ",3,co2cloudfrac)
      call WRITEdiagfi(ngrid,"rsedcloudco2","rsed co2"
     &        ,"m",3,rsedcloudco2)
      call WRITEdiagfi(ngrid,"Tau3D1mic"," co2 ice opacities"
     &        ," ",3,Qext1bins2)
      call WRITEdiagfi(ngrid,"tau1mic","co2 ice opacity 1 micron"
     &        ," ",2,tau1mic)
      call WRITEDIAGFI(ngrid,"mem_Nccn_co2","CCN number used by CO2"
     &        ,"kg/kg ",3,mem_Nccn_co2) 
      call WRITEDIAGFI(ngrid,"mem_Mccn_co2","CCN mass used by CO2"
     &        ,"kg/kg ",3,mem_Mccn_co2) 
      call WRITEDIAGFI(ngrid,"mem_Mh2o_co2","H2O mass in CO2 crystal"
     &        ,"kg/kg ",3,mem_Mh2o_co2)         
c D.BARDET: sensibility test 
c      call WRITEDIAGFI(ngrid,"No_dust","Nombre particules de poussiere"
c     &        ,"part/kg",3,No_dust)
c      call WRITEDIAGFI(ngrid,"Mo_dust","Masse particules de poussiere"
c     &        ,"kg/kg ",3,Mo_dust) 
      END SUBROUTINE co2cloud

c ===================================================================
c Subroutines used to write variables of memory in start files       
c ===================================================================

      SUBROUTINE ini_co2cloud(ngrid,nlayer)
  
      IMPLICIT NONE

      INTEGER, INTENT (in) :: ngrid  ! number of atmospheric columns
      INTEGER, INTENT (in) :: nlayer ! number of atmospheric layers

         allocate(mem_Nccn_co2(ngrid,nlayer))
         allocate(mem_Mccn_co2(ngrid,nlayer))
         allocate(mem_Mh2o_co2(ngrid,nlayer))

      END SUBROUTINE ini_co2cloud
c ----------------------------------
      SUBROUTINE end_co2cloud

      IMPLICIT NONE

         if (allocated(mem_Nccn_co2)) deallocate(mem_Nccn_co2)
         if (allocated(mem_Mccn_co2)) deallocate(mem_Mccn_co2)
         if (allocated(mem_Mh2o_co2)) deallocate(mem_Mh2o_co2)

      END SUBROUTINE end_co2cloud

      END MODULE co2cloud_mod
