      MODULE improvedco2clouds_mod

      IMPLICIT NONE

      CONTAINS

      subroutine improvedco2clouds(ngrid,nlay,microtimestep,
     &             pplay,pplev,pteff,sum_subpdt,
     &             pqeff,sum_subpdq,subpdqcloudco2,subpdtcloudco2,
     &             nq,tauscaling,
     &             mem_Mccn_co2,mem_Mh2o_co2,mem_Nccn_co2,
     &             No_dust,Mo_dust)
      USE comcstfi_h, only: pi, g, cpp
      USE updaterad, only: updaterice_micro, updaterice_microco2,
     &                     updaterccnCO2
      use tracer_mod, only: igcm_dust_mass, igcm_dust_number, rho_dust,
     &                      igcm_h2o_ice, igcm_ccn_mass,
     &                      igcm_ccn_number, nuice_sed,
     &                      igcm_co2, igcm_co2_ice, igcm_ccnco2_mass,
     &                      igcm_ccnco2_number, nuiceco2_sed,
     &                      rho_ice_co2,nuiceco2_ref
      use conc_mod, only: mmean
      use datafile_mod, only: datadir

      implicit none
      
c------------------------------------------------------------------
c  This routine is used to form CO2 clouds when a parcel of the GCM is
c    saturated. It includes the ability to have supersaturation, a
c    computation of the nucleation rates, growthrates and the
c    scavenging of dust particles by clouds.
c  It is worth noting that the amount of dust is computed using the
c    dust optical depth computed in aeropacity.F. That's why
c    the variable called "tauscaling" is used to convert
c    pq(dust_mass) and pq(dust_number), which are relative
c    quantities, to absolute and realistic quantities stored in zq.
c    This has to be done to convert the inputs into absolute
c    values, but also to convert the outputs back into relative
c    values which are then used by the sedimentation and advection
c    schemes.
c CO2 ice particles can nucleate on both dust and on water ice particles
c When CO2 ice is deposited onto a water ice particles, the particle is
c removed from the water tracers.
c Memory of the origin of the co2 particles is kept and thus the 
c water cycle shouldn't be modified by this.
c WARNING: no sedimentation of the water ice origin is performed
c in the microphysical timestep in co2cloud.F. 

c  Authors of the water ice clouds microphysics
c J.-B. Madeleine, based on the work by Franck Montmessin
c           (October 2011)
c           T. Navarro, debug,correction, new scheme (October-April 2011)
c           A. Spiga, optimization (February 2012)
c Adaptation for CO2 clouds by Joachim Audouard (09/16), based on the work
c of Constantino Listowski 
c There is an energy limit to how much co2 can sublimate/condensate. It is 
c defined by the difference of the GCM temperature with the co2 condensation 
c temperature. 
c Warning:
c If meteoritic particles are activated and turn into co2 ice particles,
c then they will be reversed in the dust tracers if the cloud sublimates
 
c------------------------------------------------------------------
      include "callkeys.h"
      include "microphys.h"
c------------------------------------------------------------------
c     Arguments:

      INTEGER,INTENT(in) :: ngrid,nlay
      integer,intent(in) :: nq         ! number of tracers
      REAL,INTENT(in) :: microtimestep     ! physics time step (s)
      REAL,INTENT(in) :: pplay(ngrid,nlay)     ! mid-layer pressure (Pa)
      REAL,INTENT(in) :: pplev(ngrid,nlay+1)   ! inter-layer pressure (Pa)
      REAL,INTENT(in) :: pteff(ngrid,nlay) ! temperature at the middle of the
                                 !   layers (K)
      REAL,INTENT(in) :: sum_subpdt(ngrid,nlay) ! tendency on temperature from
                                 !  previous physical parametrizations
      REAL,INTENT(in) :: pqeff(ngrid,nlay,nq) ! tracers (kg/kg)
      REAL,INTENT(in) :: sum_subpdq(ngrid,nlay,nq) ! tendencies on tracers 
                                 !  before condensation (kg/kg.s-1)
      REAL,INTENT(in) :: tauscaling(ngrid) ! Convertion factor for qdust and Ndust
c     Outputs:
      REAL,INTENT(out) :: subpdqcloudco2(ngrid,nlay,nq) ! tendency on tracers
                                   ! due to CO2 condensation (kg/kg.s-1)
      ! condensation si igcm_co2_ice 
      REAL,INTENT(out) :: subpdtcloudco2(ngrid,nlay)  ! tendency on temperature due
                                   ! to latent heat
      REAL,INTENT(out) :: No_dust(ngrid,nlay)
      REAL,INTENT(out) :: Mo_dust(ngrid,nlay)

c------------------------------------------------------------------
c     Local variables:
      LOGICAL,SAVE :: firstcall=.true.
      REAL*8   derf ! Error function
      INTEGER ig,l,i

      real masse (ngrid,nlay) ! Layer mass (kg.m-2)
      REAL rice(ngrid,nlay)    ! Water Ice mass mean radius (m)
                                ! used for nucleation of CO2 on ice-coated ccns
      REAL rccnh2o(ngrid,nlay)    ! Water Ice mass mean radius (m)
      REAL zq(ngrid,nlay,nq)  ! local value of tracers
      REAL zq0(ngrid,nlay,nq) ! local initial value of tracers
      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL zqsat(ngrid,nlay)    ! saturation vapor pressure for CO2
      real tcond(ngrid,nlay)
      real zqco2(ngrid,nlay)
      REAL lw                         !Latent heat of sublimation (J.kg-1) 
      REAL,save :: l0,l1,l2,l3,l4
      DOUBLE PRECISION dMice           ! mass of condensed ice
      DOUBLE PRECISION sumcheck
      DOUBLE PRECISION facteurmax!for energy limit on mass growth
      DOUBLE PRECISION pco2,psat  ! Co2 vapor partial pressure (Pa)
      DOUBLE PRECISION satu ! Co2 vapor saturation ratio over ice
      DOUBLE PRECISION Mo ,No
      
c D.BARDET: sensibility test
c      REAL, SAVE :: No ! when sensibility test
c      DOUBLE PRECISION No_dust(ngrid,nlay) ! when sensibility test
c      DOUBLE PRECISION Mo_dust(ngrid,nlay) ! when sensibility test 
     
      DOUBLE PRECISION  Rn, Rm, dev2,dev3, n_derf, m_derf
      DOUBLE PRECISION mem_Mccn_co2(ngrid,nlay) ! Memory of CCN mass of H2O and dust used by CO2
      DOUBLE PRECISION mem_Mh2o_co2(ngrid,nlay) ! Memory of H2O mass integred into CO2 crystal
      DOUBLE PRECISION mem_Nccn_co2(ngrid,nlay) ! Memory of CCN number of H2O and dust used by CO2
      DOUBLE PRECISION interm1,interm2,interm3     
 
!     Radius used by the microphysical scheme (m)
      DOUBLE PRECISION n_aer(nbinco2_cld) ! number concentration volume-1 of particle/each size bin
      DOUBLE PRECISION m_aer(nbinco2_cld) ! mass mixing ratio of particle/each size bin
      DOUBLE PRECISION m_aer_h2oice2(nbinco2_cld) ! mass mixing ratio of particle/each size bin

      DOUBLE PRECISION n_aer_h2oice(nbinco2_cld) ! Same - for CO2 nucleation
      DOUBLE PRECISION m_aer_h2oice(nbinco2_cld) ! Same - for CO2 nucleation
      DOUBLE PRECISION rad_h2oice(nbinco2_cld) 

c      REAL*8 sigco2      ! Co2-ice/air surface tension  (N.m)
c      EXTERNAL sigco2

      DOUBLE PRECISION dN,dM, dNh2o, dMh2o, dNN,dMM,dNNh2o,dMMh2o
      DOUBLE PRECISION dMh2o_ice,dMh2o_ccn

      DOUBLE PRECISION rate(nbinco2_cld)  ! nucleation rate
      DOUBLE PRECISION rateh2o(nbinco2_cld)  ! nucleation rate
      REAL seq
      DOUBLE PRECISION rho_ice_co2T(ngrid,nlay)
      DOUBLE PRECISION riceco2(ngrid,nlay)      ! CO2Ice mean radius (m)
      REAL rhocloud(ngrid,nlay) ! Cloud density (kg.m-3)
                  
      REAL rhocloudco2(ngrid,nlay)  ! Cloud density (kg.m-3)
      REAL rdust(ngrid,nlay) ! Dust geometric mean radius (m)

c      REAL res      ! Resistance growth
      DOUBLE PRECISION Ic_rice      ! Mass transfer rate CO2 ice crystal
      DOUBLE PRECISION ratioh2o_ccn
      DOUBLE PRECISION vo2co2

c     Parameters of the size discretization used by the microphysical scheme
      DOUBLE PRECISION, PARAMETER :: rmin_cld = 1.e-9   ! Minimum radius (m)
      DOUBLE PRECISION, PARAMETER :: rmax_cld = 5.e-6   ! Maximum radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmin_cld =1.e-10  ! Minimum bounary radius (m)
      DOUBLE PRECISION, PARAMETER :: rbmax_cld = 2.e-4  ! Maximum boundary radius (m)
      DOUBLE PRECISION vrat_cld                         ! Volume ratio
      DOUBLE PRECISION rb_cldco2(nbinco2_cld+1)         ! boundary values of each rad_cldco2 bin (m)
      SAVE rb_cldco2
      DOUBLE PRECISION dr_cld(nbinco2_cld)              ! width of each rad_cldco2 bin (m)
      DOUBLE PRECISION vol_cld(nbinco2_cld)             ! particle volume for each bin (m3)

      DOUBLE PRECISION Proba,Masse_atm,drsurdt,reff,Probah2o
      REAL sigma_iceco2   ! Variance of the co2 ice and CCN distributions
      SAVE sigma_iceco2 
      REAL sigma_ice      ! Variance of the h2o ice and CCN distributions
      SAVE sigma_ice
      DOUBLE PRECISION Niceco2,Qccnco2,Nccnco2

!     Variables for the meteoritic flux:
      integer,parameter :: nbin_meteor=100
      integer,parameter :: nlev_meteor=130
      double precision meteor_ccn(ngrid,nlay,100) !100=nbinco2_cld !!! 
      double precision,save :: meteor(130,100)
      double precision mtemp(100),pression_meteor(130)
      logical file_ok
      integer read_ok
      integer nelem,lebon1,lebon2
      double precision :: ltemp1(130),ltemp2(130)
      integer ibin,j
      integer,parameter :: uMeteor=666

      IF (firstcall) THEN
!=============================================================
! 0. Definition of the size grid
!=============================================================
c       rad_cldco2 is the primary radius grid used for microphysics computation.
c       The grid spacing is computed assuming a constant volume ratio
c       between two consecutive bins; i.e. vrat_cld.
c       vrat_cld is determined from the boundary values of the size grid: 
c       rmin_cld and rmax_cld.
c       The rb_cldco2 array contains the boundary values of each rad_cldco2 bin.
c       dr_cld is the width of each rad_cldco2 bin.

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
        print*, ' '
        print*,'Microphysics co2: size bin information:'
        print*,'i,rb_cldco2(i), rad_cldco2(i),dr_cld(i)'
        print*,'-----------------------------------'
        do i=1,nbinco2_cld
          write(*,'(i3,3x,3(e13.6,4x))') i,rb_cldco2(i), rad_cldco2(i),
     &      dr_cld(i)
        enddo
        write(*,'(i3,3x,e13.6)') nbinco2_cld+1,rb_cldco2(nbinco2_cld+1)
        print*,'-----------------------------------'
        do i=1,nbinco2_cld+1
            rb_cldco2(i) = log(rb_cldco2(i))  !! we save that so that it is not computed
                                         !! at each timestep and gridpoint
        enddo
c       Contact parameter of co2 ice on dst ( m=cos(theta) )
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        mteta  = 0.952
c        mtetaco2 = 0.952
c        write(*,*) 'co2_param contact parameter:', mtetaco2

c       Volume of a co2 molecule (m3)
        vo1co2 = m0co2 / dble(rho_ice_co2) ! m0co2 et non mco2

c       Variance of the ice and CCN distributions
        sigma_iceco2 = sqrt(log(1.+nuiceco2_sed))
        sigma_ice = sqrt(log(1.+nuice_sed))

  
        write(*,*) 'Variance of ice & CCN CO2 distribs :', sigma_iceco2
        write(*,*) 'nuice for co2 ice sedimentation:', nuiceco2_sed
        write(*,*) 'Volume of a co2 molecule:', vo1co2
        
        if (co2useh2o) then
           write(*,*)
           write(*,*) "co2useh2o=.true. in callphys.def"
           write(*,*) "This means water ice particles can "
           write(*,*) "serve as CCN for CO2 microphysics"
        endif

        if (meteo_flux) then
           write(*,*)
           write(*,*) "meteo_flux=.true. in callphys.def"
           write(*,*) "meteoritic dust particles are available"
           write(*,*) "for co2 ice nucleation! "
           write(*,*) "Flux given by J. Plane (pressions,size bins)"
           ! Initialisation of the flux: it is constant and is it saved
           !We must interpolate the table to the GCM pressures
           INQUIRE(FILE=TRIM(datadir)//
     &       '/Meteo_flux_Plane.dat', EXIST=file_ok)
           IF (.not. file_ok) THEN 
              write(*,*) 'file Meteo_flux_Plane.dat should be in '
     &             ,trim(datadir)
              STOP
           endif
!used Variables
           open(unit=uMeteor,file=trim(datadir)//
     &          '/Meteo_flux_Plane.dat'
     &          ,FORM='formatted')
!13000 records (130 pressions x 100 bin sizes)
           read(uMeteor,*) !skip 1 line
           do i=1,130 
              read(uMeteor,*,iostat=read_ok) pression_meteor(i)
              if (read_ok==0) then
                write(*,*) pression_meteor(i)
              else
                write(*,*) 'Error reading Meteo_flux_Plane.dat'
                call abort_physic("CO2clouds",
     &                    "Error reading Meteo_flux_Plane.dat",1)
              endif
           enddo
           read(uMeteor,*)       !skip 1 line
           do i=1,130 
              do j=1,100        ! les mêmes 100 bins size que la distri nuclea: on touche pas
                 read(uMeteor,'(F12.6)',iostat=read_ok) meteor(i,j)
                 if (read_ok/=0) then
                   write(*,*) 'Error reading Meteo_flux_Plane.dat'
                   call abort_physic("CO2clouds",
     &                    "Error reading Meteo_flux_Plane.dat",1)
                 endif
              enddo
!On doit maintenant réinterpoler le tableau(130,100) sur les pressions du GCM (nlay,100)
           enddo
           close(uMeteor)
        write(*,*) "File meteo_flux read, end of firstcall in co2 micro"
        endif                     !of if meteo_flux
      
        ! Parameter values for Latent heat computation
        l0=595594d0      
        l1=903.111d0    
        l2=-11.5959d0   
        l3=0.0528288d0
        l4=-0.000103183d0
        
c D.BARDET:
c        No = 1.e10        
        firstcall=.false.
      END IF
!=============================================================
! 1. Initialisation
!=============================================================

      meteor_ccn(:,:,:)=0.
      rice(:,:) = 1.e-8
      riceco2(:,:) = 1.e-11

c     Initialize the tendencies
      subpdqcloudco2(1:ngrid,1:nlay,1:nq)=0.
      subpdtcloudco2(1:ngrid,1:nlay)=0.
      
c pteff temperature layer; sum_subpdt dT.s-1
c pqeff traceur kg/kg; sum_subpdq tendance idem .s-1
      zt(1:ngrid,1:nlay) = 
     &      pteff(1:ngrid,1:nlay) + 
     &      sum_subpdt(1:ngrid,1:nlay) * microtimestep 
      zq(1:ngrid,1:nlay,1:nq) = 
     &      pqeff(1:ngrid,1:nlay,1:nq) + 
     &      sum_subpdq(1:ngrid,1:nlay,1:nq) * microtimestep
      WHERE( zq(1:ngrid,1:nlay,1:nq) < 1.e-30 )
     &     zq(1:ngrid,1:nlay,1:nq) = 1.e-30
         
         zq0(1:ngrid,1:nlay,1:nq) = zq(1:ngrid,1:nlay,1:nq)
!=============================================================
! 2. Compute saturation
!=============================================================
      dev2 = 1. / ( sqrt(2.) * sigma_iceco2 )
      dev3 = 1. / ( sqrt(2.) * sigma_ice )
      call co2sat(ngrid*nlay,zt,pplay,zqsat) !zqsat is psat(co2)
      zqco2=zq(:,:,igcm_co2)+zq(:,:,igcm_co2_ice)
      CALL tcondco2(ngrid,nlay,pplay,zqco2,tcond)
!=============================================================
! 3. Bonus: additional meteoritic particles for nucleation
!=============================================================
      if (meteo_flux) then
         !pression_meteo(130)
         !pplev(ngrid,nlay+1)
         !meteo(130,100)
         !resultat: meteo_ccn(ngrid,nlay,100)
         do l=1,nlay
            do ig=1,ngrid
               masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g 
               ltemp1=abs(pression_meteor(:)-pplev(ig,l))
               ltemp2=abs(pression_meteor(:)-pplev(ig,l+1))
               lebon1=minloc(ltemp1,DIM=1) 
               lebon2=minloc(ltemp2,DIM=1)
               nelem=lebon2-lebon1+1.
               mtemp(:)=0d0     !mtemp(100) : valeurs pour les 100bins
               do ibin=1,100
                  mtemp(ibin)=sum(meteor(lebon1:lebon2,ibin))
               enddo
               meteor_ccn(ig,l,:)=mtemp(:)/nelem/masse(ig,l) !Par kg air
csi par m carre, x epaisseur/masse pour par kg/air 
               !write(*,*) "masse air ig l=",masse(ig,l)
               !check original unit with J. Plane
            enddo
         enddo
      endif
c ------------------------------------------------------------------------
c ---------  Actual microphysics : Main loop over the GCM's grid ---------
c ------------------------------------------------------------------------
       DO l=1,nlay
         DO ig=1,ngrid
c       Get the partial pressure of co2 vapor and its saturation ratio
           pco2 = zq(ig,l,igcm_co2) * (mmean(ig,l)/44.01) * pplay(ig,l)
           satu = pco2 / zqsat(ig,l)

           rho_ice_co2T(ig,l)=1000.*(1.72391-2.53e-4*zt(ig,l)
     &          -2.87e-6*zt(ig,l)*zt(ig,l)) !T-dependant CO2 ice density
           vo2co2 = m0co2 / dble(rho_ice_co2T(ig,l))
           rho_ice_co2=rho_ice_co2T(ig,l)

!=============================================================
!4. Nucleation
!=============================================================
           IF ( satu .ge. 1 ) THEN ! if there is condensation

          call updaterccnCO2(zq(ig,l,igcm_dust_mass),
     &          zq(ig,l,igcm_dust_number),rdust(ig,l),tauscaling(ig))

c D.BARDET: sensibility test
c              rdust=2.e-6
      
c Expand the dust moments into a binned distribution
              
              n_aer(:)=0.
              m_aer(:)=0.

              Mo =4.*pi*rho_dust*No*rdust(ig,l)**(3.)
     &           *exp(9.*nuiceco2_ref/2.)/3. ! in Madeleine et al 2011
     
              No = zq(ig,l,igcm_dust_number)* tauscaling(ig)+1.e-30

              No_dust=No
              Mo_dust=Mo

              Rn = rdust(ig,l)
              Rn = -log(Rn) 
              Rm = Rn - 3. * sigma_iceco2*sigma_iceco2  
              n_derf = derf( (rb_cldco2(1)+Rn) *dev2)
              m_derf = derf( (rb_cldco2(1)+Rm) *dev2)

              do i = 1, nbinco2_cld
                 n_aer(i) = -0.5 * No * n_derf !! this ith previously computed
                 m_aer(i) = -0.5 * Mo * m_derf !! this ith previously computed
                 n_derf = derf((rb_cldco2(i+1)+Rn) *dev2)
                 m_derf = derf((rb_cldco2(i+1)+Rm) *dev2)
                 n_aer(i) = n_aer(i) + 0.5 * No * n_derf
                 m_aer(i) = m_aer(i) + 0.5 * Mo * m_derf
              enddo

c Ajout meteor_ccn particles aux particules de poussière background
              if (meteo_flux) then
                 do i = 1, nbinco2_cld
                    n_aer(i) = n_aer(i) + meteor_ccn(ig,l,i) 
                    m_aer(i) = m_aer(i) + 4./3.*pi*rho_dust
     &                *meteor_ccn(ig,l,i)*rad_cldco2(i)*rad_cldco2(i)
     &                *rad_cldco2(i)
                 enddo
              endif
              
c Same but with h2o particles as CCN only if co2useh2o=.true.
              
              n_aer_h2oice(:)=0.
              m_aer_h2oice(:)=0.

              if (co2useh2o) then
                call updaterice_micro(zq(ig,l,igcm_h2o_ice),
     &               zq(ig,l,igcm_ccn_mass),zq(ig,l,igcm_ccn_number),
     &                 tauscaling(ig),rice(ig,l),rhocloud(ig,l))
                Mo = zq(ig,l,igcm_h2o_ice) +
     &               zq(ig,l,igcm_ccn_mass)*tauscaling(ig)+1.e-30 
                     ! Total mass of H20 crystals,CCN included
                No = zq(ig,l,igcm_ccn_number)* tauscaling(ig) + 1.e-30
                Rn = rice(ig,l)
                Rn = -log(Rn) 
                Rm = Rn - 3. * sigma_ice*sigma_ice  
                n_derf = derf( (rb_cldco2(1)+Rn) *dev3)
                m_derf = derf( (rb_cldco2(1)+Rm) *dev3)
                do i = 1, nbinco2_cld
                  n_aer_h2oice(i) = -0.5 * No * n_derf 
                  m_aer_h2oice(i) = -0.5 * Mo * m_derf 
                  n_derf = derf( (rb_cldco2(i+1)+Rn) *dev3)
                  m_derf = derf( (rb_cldco2(i+1)+Rm) *dev3)
                  n_aer_h2oice(i) = n_aer_h2oice(i) + 0.5 * No * n_derf 
                  m_aer_h2oice(i) = m_aer_h2oice(i) + 0.5 * Mo * m_derf 
                  rad_h2oice(i) = rad_cldco2(i)
                enddo
              endif
             

! Call to nucleation routine
              call nucleaco2(dble(pco2),zt(ig,l),dble(satu)
     &             ,n_aer,rate,n_aer_h2oice
     &             ,rad_h2oice,rateh2o,vo2co2)
              dN = 0.
              dM = 0.
              dNh2o = 0.
              dMh2o = 0.
              do i = 1, nbinco2_cld
                 Proba    =1.0-exp(-1.*microtimestep*rate(i))
                 dN       = dN + n_aer(i) * Proba
                 dM       = dM + m_aer(i) * Proba             
              enddo
              if (co2useh2o) then
                 do i = 1, nbinco2_cld
                    Probah2o = 1.0-exp(-1.*microtimestep*rateh2o(i))
                    dNh2o    = dNh2o + n_aer_h2oice(i) * Probah2o
                    dMh2o    = dMh2o + m_aer_h2oice(i) * Probah2o
                 enddo
              endif

! dM  masse activée (kg) et dN nb particules par  kg d'air
! Now increment CCN tracers and update dust tracers
              dNN= dN/tauscaling(ig)
              dMM= dM/tauscaling(ig)
              dNN=min(dNN,zq(ig,l,igcm_dust_number))
              dMM=min(dMM,zq(ig,l,igcm_dust_mass))
              zq(ig,l,igcm_ccnco2_mass)   = 
     &             zq(ig,l,igcm_ccnco2_mass)   + dMM
              zq(ig,l,igcm_ccnco2_number) =
     &             zq(ig,l,igcm_ccnco2_number) + dNN
              zq(ig,l,igcm_dust_mass)= zq(ig,l,igcm_dust_mass)-dMM 
              zq(ig,l,igcm_dust_number)=zq(ig,l,igcm_dust_number)-dNN


c Update CCN for CO2 nucleating on H2O CCN :
              ! Warning: must keep memory of it 
              if (co2useh2o) then
                 dNNh2o=dNh2o/tauscaling(ig)
                 dNNh2o=min(dNNh2o,zq(ig,l,igcm_ccn_number))
                 ratioh2o_ccn=1./(zq(ig,l,igcm_h2o_ice)
     &                +zq(ig,l,igcm_ccn_mass)*tauscaling(ig))  
                 dMh2o_ice=dMh2o*zq(ig,l,igcm_h2o_ice)*ratioh2o_ccn
                 dMh2o_ccn=dMh2o*zq(ig,l,igcm_ccn_mass)*
     &                tauscaling(ig)*ratioh2o_ccn
                 dMh2o_ccn=dMh2o_ccn/tauscaling(ig)
                 dMh2o_ccn=min(dMh2o_ccn,zq(ig,l,igcm_ccn_mass))
                 dMh2o_ice=min(dMh2o_ice,zq(ig,l,igcm_h2o_ice))
                 zq(ig,l,igcm_ccnco2_mass)   = 
     &                zq(ig,l,igcm_ccnco2_mass)  + dMh2o_ice+dMh2o_ccn
                 zq(ig,l,igcm_ccnco2_number) = 
     &                zq(ig,l,igcm_ccnco2_number) + dNNh2o
                zq(ig,l,igcm_ccn_number)=zq(ig,l,igcm_ccn_number)-dNNh2o 
                zq(ig,l,igcm_h2o_ice) = zq(ig,l,igcm_h2o_ice)-dMh2o_ice
                zq(ig,l,igcm_ccn_mass)= zq(ig,l,igcm_ccn_mass)-dMh2o_ccn
                mem_Mh2o_co2(ig,l)=mem_Mh2o_co2(ig,l)+dMh2o_ice
                mem_Mccn_co2(ig,l)=mem_Mccn_co2(ig,l)+dMh2o_ccn
                mem_Nccn_co2(ig,l)=mem_Nccn_co2(ig,l)+dNNh2o
             endif ! of if co2useh2o
           ENDIF   ! of is satu >1

!=============================================================
! 5. Ice growth: scheme for radius evolution
!=============================================================

c We trigger crystal growth if and only if there is at least one nuclei (N>1).
c Indeed, if we are supersaturated and still don't have at least one nuclei, we should better wait
c to avoid unrealistic value for nuclei radius and so on for cases that remain negligible.
             IF (zq(ig,l,igcm_ccnco2_number)
     &         * tauscaling(ig)+1.e-30.ge. 1)THEN   ! we trigger crystal growth
              
              call updaterice_microco2(dble(zq(ig,l,igcm_co2_ice)),
     &            dble(zq(ig,l,igcm_ccnco2_mass)),
     &            dble(zq(ig,l,igcm_ccnco2_number)),
     &            tauscaling(ig),riceco2(ig,l),rhocloudco2(ig,l))

              Ic_rice=0.
              lw = l0 + l1 * zt(ig,l) + l2 * zt(ig,l)**2 + 
     &             l3 * zt(ig,l)**3 + l4 * zt(ig,l)**4 !J.kg-1
              facteurmax=abs(Tcond(ig,l)-zt(ig,l))*cpp/lw
              !specific heat of co2 ice = 1000 J.kg-1.K-1 
              !specific heat of atm cpp = 744.5 J.kg-1.K-1

c call scheme of microphys. mass growth for CO2
              call massflowrateco2(pplay(ig,l),zt(ig,l),
     &             satu,riceco2(ig,l),mmean(ig,l),Ic_rice) 
c Ic_rice Mass transfer rate (kg/s) for a rice particle >0 si croissance ! 
             
              if ( (Ic_rice.ne.Ic_rice) ! will be true if it is Nan  
     &               .or. (Ic_rice.eq.0.) ) then
                 Ic_rice=0.
                 subpdtcloudco2(ig,l)=-sum_subpdt(ig,l)
                 dMice=0
                 
              else
                 dMice=zq(ig,l,igcm_ccnco2_number)*Ic_rice*microtimestep
     &                *tauscaling(ig) ! Kg par kg d'air, >0 si croissance !
                 !kg.s-1 par particule * nb particule par kg air*s
                 ! = kg par kg air 
              
              dMice = max(dMice,max(-facteurmax,-zq(ig,l,igcm_co2_ice)))
              dMice = min(dMice,min(facteurmax,zq(ig,l,igcm_co2)))
! facteurmax maximum quantity of CO2 that can sublime/condense according to available thermal energy
! latent heat release       >0 if growth i.e. if dMice >0
              subpdtcloudco2(ig,l)=dMice*lw/cpp/microtimestep
! kgco2/kgair* J/kgco2 * 1/(J.kgair-1.K-1)/s= K par seconde
              !Now update tracers 
              zq(ig,l,igcm_co2_ice) = zq(ig,l,igcm_co2_ice)+dMice
              zq(ig,l,igcm_co2) = zq(ig,l,igcm_co2)-dMice
              endif

!=============================================================
! 6. Dust cores releasing if no more co2 ice :
!=============================================================

              if (zq(ig,l,igcm_co2_ice).le. 1.e-25)THEN
! On sublime tout
                 if (co2useh2o) then
                   if (mem_Mccn_co2(ig,l) .gt. 0) then
                    zq(ig,l,igcm_ccn_mass)=zq(ig,l,igcm_ccn_mass)
     &                   +mem_Mccn_co2(ig,l)
                   endif
                   if (mem_Mh2o_co2(ig,l) .gt. 0) then
                    zq(ig,l,igcm_h2o_ice)=zq(ig,l,igcm_h2o_ice)
     &                   +mem_Mh2o_co2(ig,l)
                   endif
                 
                   if (mem_Nccn_co2(ig,l) .gt. 0) then
                    zq(ig,l,igcm_ccn_number)=zq(ig,l,igcm_ccn_number)
     &                   +mem_Nccn_co2(ig,l)
                   endif
                 endif
                    zq(ig,l,igcm_dust_mass) = 
     &                   zq(ig,l,igcm_dust_mass)
     &                   + zq(ig,l,igcm_ccnco2_mass)-
     &                   (mem_Mh2o_co2(ig,l)+mem_Mccn_co2(ig,l))
                    zq(ig,l,igcm_dust_number) = 
     &                   zq(ig,l,igcm_dust_number)
     &                   + zq(ig,l,igcm_ccnco2_number)
     &                   -mem_Nccn_co2(ig,l)
                 
                    zq(ig,l,igcm_co2) = zq(ig,l,igcm_co2) 
     &                   + zq(ig,l,igcm_co2_ice)
                 
                 zq(ig,l,igcm_ccnco2_mass)=0.
                 zq(ig,l,igcm_co2_ice)=0.
                 zq(ig,l,igcm_ccnco2_number)=0.
                 mem_Nccn_co2(ig,l)=0.
                 mem_Mh2o_co2(ig,l)=0.
                 mem_Mccn_co2(ig,l)=0.
                 riceco2(ig,l)=0.

              endif !of if co2_ice <1e-25

              ENDIF            ! of if NCCN > 1
          ENDDO                ! of ig loop
        ENDDO                  ! of nlayer loop  

!=============================================================
! 7. END: get cloud tendencies 
!=============================================================

          ! Get cloud tendencies
        subpdqcloudco2(1:ngrid,1:nlay,igcm_co2) =
     &       (zq(1:ngrid,1:nlay,igcm_co2) - 
     &       zq0(1:ngrid,1:nlay,igcm_co2))/microtimestep
        subpdqcloudco2(1:ngrid,1:nlay,igcm_co2_ice) =
     &       (zq(1:ngrid,1:nlay,igcm_co2_ice) -
     &       zq0(1:ngrid,1:nlay,igcm_co2_ice))/microtimestep
        
        if (co2useh2o) then
           subpdqcloudco2(1:ngrid,1:nlay,igcm_h2o_ice) =
     &         (zq(1:ngrid,1:nlay,igcm_h2o_ice) - 
     &         zq0(1:ngrid,1:nlay,igcm_h2o_ice))/microtimestep
           subpdqcloudco2(1:ngrid,1:nlay,igcm_ccn_mass) =
     &         (zq(1:ngrid,1:nlay,igcm_ccn_mass) -
     &         zq0(1:ngrid,1:nlay,igcm_ccn_mass))/microtimestep
           subpdqcloudco2(1:ngrid,1:nlay,igcm_ccn_number) =
     &         (zq(1:ngrid,1:nlay,igcm_ccn_number) -
     &         zq0(1:ngrid,1:nlay,igcm_ccn_number))/microtimestep
        endif

        subpdqcloudco2(1:ngrid,1:nlay,igcm_ccnco2_mass) =
     &       (zq(1:ngrid,1:nlay,igcm_ccnco2_mass) -
     &       zq0(1:ngrid,1:nlay,igcm_ccnco2_mass))/microtimestep
        subpdqcloudco2(1:ngrid,1:nlay,igcm_ccnco2_number) =
     &       (zq(1:ngrid,1:nlay,igcm_ccnco2_number) -
     &       zq0(1:ngrid,1:nlay,igcm_ccnco2_number))/microtimestep
        subpdqcloudco2(1:ngrid,1:nlay,igcm_dust_mass) =
     &       (zq(1:ngrid,1:nlay,igcm_dust_mass) -
     &       zq0(1:ngrid,1:nlay,igcm_dust_mass))/microtimestep 
        subpdqcloudco2(1:ngrid,1:nlay,igcm_dust_number) =
     &       (zq(1:ngrid,1:nlay,igcm_dust_number) -
     &       zq0(1:ngrid,1:nlay,igcm_dust_number))/microtimestep


c     TEST D.BARDET 
      call WRITEDIAGFI(ngrid,"No_dust","Nombre particules de poussière"
     &        ,"part/kg",3,No_dust)
      call WRITEDIAGFI(ngrid,"Mo_dust","Masse particules de poussière"
     &        ,"kg/kg ",3,Mo_dust)      

        END SUBROUTINE improvedco2clouds

        END MODULE improvedco2clouds_mod

