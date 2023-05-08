MODULE calchim_mod

  IMPLICIT NONE

  INTEGER, SAVE :: ichemistry ! compute chemistry every ichemistry physics step
  REAL,SAVE,ALLOCATABLE :: zdqchim(:,:,:) ! Tendancy on pq due to photochemistry
  REAL,SAVE,ALLOCATABLE :: zdqschim(:,:) ! Tendancy on qsurf due to photochemistry

  CONTAINS

      subroutine calchim(ngrid,nlayer,nq,                           &
                         ptimestep,pplay,pplev,pt,pdt,dist_sol,mu0,         &
                         zzlev,zzlay,zday,pq,pdq,dqchim,dqschim,dqcloud,    &
                         dqscloud,tau,co2ice,                               &
                         pu,pdu,pv,pdv,surfdust,surfice)

      use tracer_mod, only: igcm_co2, igcm_co, igcm_o, igcm_o1d, igcm_o2, &
                            igcm_o3, igcm_h, igcm_h2, igcm_oh, igcm_ho2,  &
                            igcm_h2o2, igcm_ch4, igcm_n2, igcm_h2o_vap,   &
                            igcm_no, igcm_n, igcm_no2, igcm_n2d,          &
                            igcm_o2plus, igcm_co2plus, igcm_oplus,        &
                            igcm_coplus, igcm_cplus, igcm_nplus,          &
                            igcm_noplus, igcm_n2plus, igcm_hplus,         &
                            igcm_hco2plus, igcm_hcoplus, igcm_elec, mmol

      use conc_mod, only: mmean ! mean molecular mass of the atmosphere
      use comcstfi_h, only: pi
      use photolysis_mod, only: init_photolysis, nphot
      use iono_h, only: temp_elect

      implicit none

!=======================================================================
!
!   subject:
!   --------
!
!  Prepare the call for the photochemical module, and send back the
!  tendencies from photochemistry in the chemical species mass mixing ratios
!
!   Arguments:
!   ----------
!
!  Input:
!
!    ptimestep              timestep (s)
!    pplay(ngrid,nlayer)    Pressure at the middle of the layers (Pa)
!    pplev(ngrid,nlayer+1)  Intermediate pressure levels (Pa)
!    pt(ngrid,nlayer)       Temperature (K)
!    pdt(ngrid,nlayer)      Temperature tendency (K)
!    pu(ngrid,nlayer)       u component of the wind (ms-1)
!    pdu(ngrid,nlayer)      u component tendency (K)
!    pv(ngrid,nlayer)       v component of the wind (ms-1)
!    pdv(ngrid,nlayer)      v component tendency (K)
!    dist_sol               distance of the sun (AU)
!    mu0(ngrid)             cos of solar zenith angle (=1 when sun at zenith)
!    pq(ngrid,nlayer,nq)    advected fields, ie chemical species here
!    pdq(ngrid,nlayer,nq)   previous tendencies on pq
!    tau(ngrid)             dust optical depth
!    co2ice(ngrid)          co2 ice surface layer (kg.m-2)
!    surfdust(ngrid,nlayer) dust surface area (m2/m3)
!    surfice(ngrid,nlayer)  ice surface area (m2/m3)
!
!  Output:
!
!    dqchim(ngrid,nlayer,nq) tendencies on pq due to chemistry
!    dqschim(ngrid,nq)       tendencies on qsurf 
!
!=======================================================================

include "callkeys.h"

!     input:

      integer,intent(in) :: ngrid    ! number of atmospheric columns
      integer,intent(in) :: nlayer   ! number of atmospheric layers
      integer,intent(in) :: nq       ! number of tracers
      real :: ptimestep
      real :: pplay(ngrid,nlayer)    ! pressure at the middle of the layers
      real :: zzlay(ngrid,nlayer)    ! pressure at the middle of the layers
      real :: pplev(ngrid,nlayer+1)  ! intermediate pressure levels
      real :: zzlev(ngrid,nlayer+1)  ! altitude at layer boundaries
      real :: pt(ngrid,nlayer)       ! temperature
      real :: pdt(ngrid,nlayer)      ! temperature tendency
      real :: pu(ngrid,nlayer)       ! u component of the wind (m.s-1)
      real :: pdu(ngrid,nlayer)      ! u component tendency
      real :: pv(ngrid,nlayer)       ! v component of the wind (m.s-1)
      real :: pdv(ngrid,nlayer)      ! v component tendency
      real :: dist_sol               ! distance of the sun (AU)
      real :: mu0(ngrid)             ! cos of solar zenith angle (=1 when sun at zenith)
      real :: pq(ngrid,nlayer,nq)    ! tracers mass mixing ratio
      real :: pdq(ngrid,nlayer,nq)   ! previous tendencies
      real :: zday                   ! date (time since Ls=0, in martian days)
      real :: tau(ngrid)             ! dust optical depth
      real :: co2ice(ngrid)          ! co2 ice surface layer (kg.m-2)
      real :: surfdust(ngrid,nlayer) ! dust surface area (m2/m3)
      real :: surfice(ngrid,nlayer)  ! ice surface area (m2/m3)

!     output:

      real :: dqchim(ngrid,nlayer,nq)   ! tendencies on pq due to chemistry
      real :: dqschim(ngrid,nq)         ! tendencies on qsurf 
      real :: dqcloud(ngrid,nlayer,nq)  ! tendencies on pq due to condensation
      real :: dqscloud(ngrid,nq)        ! tendencies on qsurf 

!     local variables:

      integer,save :: nbq                   ! number of tracers used in the chemistry
      integer,allocatable,save :: niq(:)    ! array storing the indexes of the tracers
      integer :: iloc(1)                    ! index of major species
      integer :: ig,l,i,iq,iqmax
      integer :: foundswitch, lswitch
      integer,save :: chemthermod

      integer,save :: i_co2  = 0
      integer,save :: i_co   = 0
      integer,save :: i_o    = 0
      integer,save :: i_o1d  = 0
      integer,save :: i_o2   = 0
      integer,save :: i_o3   = 0
      integer,save :: i_h    = 0
      integer,save :: i_h2   = 0
      integer,save :: i_oh   = 0
      integer,save :: i_ho2  = 0
      integer,save :: i_h2o2 = 0
      integer,save :: i_ch4  = 0
      integer,save :: i_n2   = 0
      integer,save :: i_h2o  = 0
      integer,save :: i_n    = 0
      integer,save :: i_no   = 0
      integer,save :: i_no2  = 0
      integer,save :: i_n2d  = 0
      integer,save :: i_co2plus=0
      integer,save :: i_oplus=0
      integer,save :: i_o2plus=0
      integer,save :: i_coplus=0
      integer,save :: i_cplus=0
      integer,save :: i_nplus=0
      integer,save :: i_noplus=0
      integer,save :: i_n2plus=0
      integer,save :: i_hplus=0
      integer,save :: i_hco2plus=0
      integer,save :: i_hcoplus=0
      integer,save :: i_elec=0

      integer :: ig_vl1

      integer :: nb_reaction_3_max     ! number of quadratic reactions
      integer :: nb_reaction_4_max     ! number of bimolecular reactions
      integer :: nquench               ! number of quenching + heterogeneous reactions
      integer :: nphotion              ! number of photoionizations
      integer :: nb_phot_max           ! total number of photolysis+photoionizations+quenching reactions


      real    :: latvl1, lonvl1
      real    :: zq(ngrid,nlayer,nq)   ! pq+pdq*ptimestep before chemistry
                                       ! new mole fraction after
      real    :: zt(ngrid,nlayer)      ! temperature
      real    :: zu(ngrid,nlayer)      ! u component of the wind
      real    :: zv(ngrid,nlayer)      ! v component of the wind
      real    :: taucol                ! dust optical depth at the surface
      real    :: kb                    ! boltzmann constant

      logical, save :: firstcall = .true.
      logical, save :: depos       ! switch for dry deposition
      logical, save :: ionchem     ! switch for ionospheric chemistry
      logical, save :: jonline     ! switch for online photodissociation rates or lookup table
      logical, save :: unichim     ! only one unified chemical scheme at all 
                                   ! layers (default), or upper atmospheric 
                                   ! scheme in the thermosphere

!     for each column of atmosphere:

      real :: zpress(nlayer)       ! Pressure (mbar)
      real :: zdens(nlayer)        ! Density  (cm-3)
      real :: ztemp(nlayer)        ! Temperature (K)
      real :: ztelec(nlayer)       ! Electronic temperature (K)
      real :: zlocal(nlayer)       ! Altitude (km)
      real :: zycol(nlayer,nq)     ! Composition (mole fractions)
      real :: zmmean(nlayer)       ! Mean molecular mass (g.mole-1)
      real :: szacol               ! Solar zenith angle
      real :: surfice1d(nlayer)    ! Ice surface area (cm2/cm3)
      real :: surfdust1d(nlayer)   ! Dust surface area (cm2/cm3)
      real :: jo3(nlayer)          ! Photodissociation rate O3->O1D (s-1)
      real :: jh2o(nlayer)         ! Photodissociation rate H2O->H+OH (s-1)
      real :: em_no(nlayer)        ! NO nightglow emission rate
      real :: em_o2(nlayer)        ! O2 nightglow emission rate      

      integer :: iter(nlayer)      ! Number of chemical iterations
                                   ! within one physical timestep 

!     for output:

      logical,save :: output        ! to issue calls to writediagfi and stats
      real :: jo3_3d(ngrid,nlayer)  ! Photodissociation rate O3->O1D (s-1)
      real :: jh2o_3d(ngrid,nlayer)  ! Photodissociation rate H2O->H+OH (s-1)
      real :: emission_no(ngrid,nlayer) !NO emission rate
      real :: emission_o2(ngrid,nlayer) !O2 emission rate
      real :: iter_3d(ngrid,nlayer) ! Number of chemical iterations
                                    !  within one physical timestep

!=======================================================================
!     initialization of the chemistry (first call only)
!=======================================================================

      if (firstcall) then

!=======================================================================
!     main dashboard for the chemistry
!=======================================================================

        unichim = .true.     ! true : unified chemistry ! false : separate models in lower and upper atmosphere
        jonline = .true.     ! true : on-line calculation of photodissociation rates ! false : lookup table
        ionchem = .false.    ! switch for ionospheric chemistry
        depos   = .false.    ! switch for dry deposition
        output  = .false.    ! true : triggers writing of specific outputs (photolysis rates, emission rates, etc)

         if (photochem) then
            if (jonline) then
               print*,'calchim: Read UV absorption cross-sections'
               call init_photolysis     ! on-line photolysis
            else
               print*,'calchim: Read photolysis lookup table'
               call read_phototable     ! off-line photolysis
            end if
         end if

         if(.not.unichim) then
            !Read reaction rates from external file if the upper atmospheric
            !chemistry is called
            call chemthermos_readini
         endif

         ! find index of chemical tracers to use
         allocate(niq(nq))
         ! Listed here are all tracers that can go into photochemistry
         nbq = 0 ! to count number of tracers
         ! Species ALWAYS present if photochem=.T.
         i_co2 = igcm_co2
         if (i_co2 == 0) then
            write(*,*) "calchim: Error; no CO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co2
         end if
         i_co = igcm_co
         if (i_co == 0) then
            write(*,*) "calchim: Error; no CO tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_co
         end if
         i_o = igcm_o
         if (i_o == 0) then
            write(*,*) "calchim: Error; no O tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o
         end if
         i_o1d = igcm_o1d
         if (i_o1d == 0) then
            write(*,*) "calchim: Error; no O1D tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o1d
         end if
         i_o2 = igcm_o2
         if (i_o2 == 0) then
            write(*,*) "calchim: Error; no O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o2
         end if
         i_o3 = igcm_o3
         if (i_o3 == 0) then
            write(*,*) "calchim: Error; no O3 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_o3
         end if
         i_h = igcm_h
         if (i_h == 0) then
            write(*,*) "calchim: Error; no H tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h
         end if
         i_h2 = igcm_h2
         if (i_h2 == 0) then
            write(*,*) "calchim: Error; no H2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2
         end if
         i_oh = igcm_oh
         if (i_oh == 0) then
            write(*,*) "calchim: Error; no OH tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_oh
         end if
         i_ho2 = igcm_ho2
         if (i_ho2 == 0) then
            write(*,*) "calchim: Error; no HO2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_ho2
         end if
         i_h2o2 = igcm_h2o2
         if (i_h2o2 == 0) then
            write(*,*) "calchim: Error; no H2O2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o2
         end if
         i_ch4 = igcm_ch4
         if (i_ch4 == 0) then
            write(*,*) "calchim: Error; no CH4 tracer !!!"
            write(*,*) "CH4 will be ignored in the chemistry"
         else
            nbq = nbq + 1
            niq(nbq) = i_ch4
         end if
         i_n2 = igcm_n2
         if (i_n2 == 0) then
            write(*,*) "calchim: Error; no N2 tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_n2
         end if
         i_n = igcm_n
         if (i_n == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no N tracer !!!"
               stop
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_n
         end if
         i_n2d = igcm_n2d
         if (i_n2d == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no N2D tracer !!!"
               stop
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_n2d
         end if
         i_no = igcm_no
         if (i_no == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no NO tracer !!!"
               stop
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_no
         end if
         i_no2 = igcm_no2
         if (i_no2 == 0) then
            if (photochem) then
               write(*,*) "calchim: Error; no NO2 tracer !!!"
               stop
            end if
         else
            nbq = nbq + 1
            niq(nbq) = i_no2
         end if
         i_h2o = igcm_h2o_vap
         if (i_h2o == 0) then
            write(*,*) "calchim: Error; no water vapor tracer !!!"
            stop
         else
            nbq = nbq + 1
            niq(nbq) = i_h2o
         end if
         i_o2plus = igcm_o2plus
         if (i_o2plus == 0) then
            write(*,*) "calchim: no O2+ tracer "
            write(*,*) "Only neutral chemistry"
         else
            nbq = nbq + 1
            niq(nbq) = i_o2plus
            ionchem = .true.
            write(*,*) "calchim: O2+ tracer found in traceur.def"
            write(*,*) "Ion chemistry included"
         end if
         i_co2plus = igcm_co2plus
         if(ionchem) then
            if (i_co2plus == 0) then
               write(*,*) "calchim: Error, no CO2+ tracer !!!"
               write(*,*) "CO2+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_co2plus
            end if
         else
            if (i_co2plus /= 0) then
               write(*,*) "calchim: Error: CO2+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
               stop
            endif
         endif
         i_oplus=igcm_oplus
         if(ionchem) then
            if (i_oplus == 0) then
               write(*,*) "calchim: Error, no O+ tracer !!!"
               write(*,*) "O+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_oplus
            end if
         else
            if (i_oplus /= 0) then
               write(*,*) "calchim: Error: O+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
               stop
            endif
         endif
         i_noplus=igcm_noplus
         if(ionchem) then
            if (i_noplus == 0) then
               write(*,*) "calchim: Error, no NO+ tracer !!!"
               write(*,*) "NO+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_noplus
            end if
         else
            if (i_noplus /= 0) then
               write(*,*) "calchim: Error: NO+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_coplus=igcm_coplus
         if(ionchem) then
            if (i_coplus == 0) then
               write(*,*) "calchim: Error, no CO+ tracer !!!"
               write(*,*) "CO+ is needed if O2+ is in traceur.def"
               stop               
            else
               nbq = nbq + 1
               niq(nbq) = i_coplus
            end if
         else
            if (i_coplus /= 0) then
               write(*,*) "calchim: Error: CO+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_cplus=igcm_cplus
         if(ionchem) then
            if (i_cplus == 0) then
               write(*,*) "calchim: Error, no C+ tracer !!!"
               write(*,*) "C+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_cplus
            end if
         else
            if (i_cplus /= 0) then
               write(*,*) "calchim: Error: C+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_n2plus=igcm_n2plus
         if(ionchem) then
            if (i_n2plus == 0) then
               write(*,*) "calchim: Error, no N2+ tracer !!!"
               write(*,*) "N2+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_n2plus
            end if
         else
            if (i_n2plus /= 0) then
               write(*,*) "calchim: Error: N2+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_nplus=igcm_nplus
         if(ionchem) then
            if (i_nplus == 0) then
               write(*,*) "calchim: Error, no N+ tracer !!!"
               write(*,*) "N+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_nplus
            end if
         else
            if (i_nplus /= 0) then
               write(*,*) "calchim: Error: N+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_hplus=igcm_hplus
         if(ionchem) then
            if (i_hplus == 0) then
               write(*,*) "calchim: Error, no H+ tracer !!!"
               write(*,*) "H+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_hplus
            end if
         else
            if (i_hplus /= 0) then
               write(*,*) "calchim: Error: H+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_hco2plus=igcm_hco2plus
         if(ionchem) then
            if (i_hco2plus == 0) then
               write(*,*) "calchim: Error, no HCO2+ tracer !!!"
               write(*,*) "HCO2+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_hco2plus
            end if
         else
            if (i_hco2plus /= 0) then
               write(*,*) "calchim: Error: HCO2+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_hcoplus=igcm_hcoplus
         if(ionchem) then
            if (i_hcoplus == 0) then
               write(*,*) "calchim: Error, no HCO+ tracer !!!"
               write(*,*) "HCO+ is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_hcoplus
            end if
         else
            if (i_hcoplus /= 0) then
               write(*,*) "calchim: Error: HCO+ is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         i_elec = igcm_elec
         if(ionchem) then
            if (i_elec == 0) then
               write(*,*) "calchim: Error, no e- tracer !!!"
               write(*,*) "e- is needed if O2+ is in traceur.def"
               stop
            else
               nbq = nbq + 1
               niq(nbq) = i_elec
            end if
         else
            if (i_elec /= 0) then
               write(*,*) "calchim: Error: e- is present, but O2+ is not!!!"
               write(*,*) "Both must be in traceur.def if ion chemistry wanted"
            endif
         endif
         write(*,*) 'calchim: found nbq    = ',nbq,' tracers'
               
         firstcall = .false.
      end if ! if (firstcall)

! Initializations

      zycol(:,:)    = 0.
      dqchim(:,:,:) = 0.
      dqschim(:,:)  = 0.

      kb = 1.3806e-23

!     latvl1= 22.27
!     lonvl1= -47.94
!     ig_vl1= 1+ int( (1.5-(latvl1-90.)*48./180.)  -2 )*64. +    &
!             int(1.5+(lonvl1+180)*64./360.)

!=======================================================================
!     loop over grid
!=======================================================================

      do ig = 1,ngrid
         
         foundswitch = 0
         do l = 1,nlayer
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               zq(ig,l,iq) = pq(ig,l,iq) + pdq(ig,l,iq)*ptimestep
               zycol(l,iq) = zq(ig,l,iq)*mmean(ig,l)/mmol(iq)
            end do
            zt(ig,l)  = pt(ig,l) + pdt(ig,l)*ptimestep
            zu(ig,l)  = pu(ig,l) + pdu(ig,l)*ptimestep
            zv(ig,l)  = pv(ig,l) + pdv(ig,l)*ptimestep
            zpress(l) = pplay(ig,l)/100.
            ztemp(l)  = zt(ig,l)
            zdens(l)  = zpress(l)/(kb*1.e4*ztemp(l))
            zlocal(l) = zzlay(ig,l)/1000.
            zmmean(l) = mmean(ig,l)
            ztelec(l) = temp_elect(zlocal(l),ztemp(l),1)
            !Electronic temperature. Index 1 -> Viking; Index 2-> MAVEN

!           surfdust1d and surfice1d: conversion from m2/m3 to cm2/cm3

            surfdust1d(l) = surfdust(ig,l)*1.e-2
            surfice1d(l)  = surfice(ig,l)*1.e-2

!           search for switch index between regions  

            if (unichim) then
               lswitch = nlayer + 1
            else
               if (foundswitch == 0 .and. pplay(ig,l) < 1.e-1) then
                  lswitch = l
                  foundswitch = 1
               end if
            endif
         end do ! of do l=1,nlayer

         szacol = acos(mu0(ig))*180./pi
         taucol = tau(ig)

!=======================================================================
!     call chemical subroutines
!=======================================================================

         if (photochem) then
            ! set number of reactions, depending on ion chemistry or not
            if (ionchem) then
               nb_reaction_4_max = 67   ! set number of bimolecular reactions
               nb_reaction_3_max = 6    ! set number of quadratic reactions
               nquench           = 9    ! set number of quenching + heterogeneous reactions
               nphotion          = 18   ! set number of photoionizations
            else
               nb_reaction_4_max = 31   ! set number of bimolecular reactions
               nb_reaction_3_max = 6    ! set number of quadratic reactions
               nquench           = 9    ! set number of quenching + heterogeneous reactions
               nphotion          = 0    ! set number of photoionizations
            end if

!        nb_phot_max is the total number of processes that are treated
!        numerically as a photolysis:

            nb_phot_max = nphot + nphotion + nquench

!        call main photochemistry routine

            call photochemistry(nlayer,nq,nbq,ionchem,nb_reaction_3_max,  &
                           nb_reaction_4_max, nb_phot_max, nphotion,      &
                           jonline, ig,lswitch,zycol,szacol,ptimestep,    &
                           zpress,zlocal,ztemp,ztelec,zdens,zmmean,       &
                           dist_sol,zday,surfdust1d,surfice1d,            &
                           jo3,jh2o,taucol,iter)

!        ozone photolysis, for output

            do l = 1,nlayer
               jo3_3d(ig,l) = jo3(l)
               jh2o_3d(ig,l) = jh2o(l)
               iter_3d(ig,l) = iter(l)
            end do

!        condensation of h2o2

            call perosat(ngrid, nlayer, nq,                        &
                         ig,ptimestep,pplev,pplay,                 &
                         ztemp,zycol,dqcloud,dqscloud)

!        case of separate photochemical model in the thermosphere

            if (.not.unichim) then
               chemthermod = 3   !C/O/H/N/ions
               call chemthermos(ig,nlayer,lswitch,chemthermod,zycol,ztemp,&
                    zdens,zpress,zlocal,szacol,ptimestep,zday,&
                    em_no,em_o2)
               do l = 1,nlayer
                  emission_no(ig,l) = em_no(l)
                  emission_o2(ig,l) = em_o2(l)
               end do
            end if

         end if  ! photochem

!        dry deposition

         if (depos) then
            call deposition(ngrid, nlayer, nq,                      &
                            ig, ig_vl1, pplay, pplev, zzlay, zzlev, & 
                            zu, zv, zt, zycol, ptimestep, co2ice)
         end if

!=======================================================================
!     tendencies
!=======================================================================

!     index of the most abundant species at each level

!         major(:) = maxloc(zycol, dim = 2)

!     tendency for the most abundant species = - sum of others

         do l = 1,nlayer
            iloc=maxloc(zycol(l,:))
            iqmax=iloc(1)
            do i = 1,nbq
               iq = niq(i) ! get tracer index
               if (iq /= iqmax) then
                  dqchim(ig,l,iq) = (zycol(l,iq)*mmol(iq)/mmean(ig,l)  &
                                   - zq(ig,l,iq))/ptimestep
                  dqchim(ig,l,iqmax) = dqchim(ig,l,iqmax)              &
                                     - dqchim(ig,l,iq) 
               end if
            end do
         end do ! of do l = 1,nlayer

!=======================================================================
!     end of loop over grid
!=======================================================================

      end do ! of do ig=1,ngrid

!=======================================================================
!     write outputs
!=======================================================================

! value of parameter 'output' to trigger writting of outputs
! is set above at the declaration of the variable.

      if (photochem .and. output) then
         if (ngrid > 1) then
            call writediagfi(ngrid,'jo3','j o3->o1d',    &
                             's-1',3,jo3_3d(1,1))
            call writediagfi(ngrid,'jh2o','jh2o',    &
                             's-1',3,jh2o_3d(1,1))
            call writediagfi(ngrid,'iter','iterations',  &
                             ' ',3,iter_3d(1,1))

            if (.not. unichim) then
               call writediagfi(ngrid,'emission_no',        &
                    'NO nightglow emission rate','cm-3 s-1',3,emission_no)
               call writediagfi(ngrid,'emission_o2',        &
                    'O2 nightglow emission rate','cm-3 s-1',3,emission_o2)
            endif
            
            if (callstats) then
               call wstats(ngrid,'jo3','j o3->o1d',       &
                           's-1',3,jo3_3d(1,1))
               call wstats(ngrid,'emission_no',           &
                   'NO nightglow emission rate','cm-3 s-1',3,emission_no)
               call wstats(ngrid,'emission_o2',           &
                   'O2 nightglow emission rate','cm-3 s-1',3,emission_o2)
               call wstats(ngrid,'mmean','mean molecular mass',       &
                           'g.mole-1',3,mmean(1,1))
            endif
         end if ! of if (ngrid.gt.1)
      end if ! of if (output)

      end subroutine calchim


    subroutine ini_calchim_mod(ngrid,nlayer,nq)
  
      implicit none
  
      integer,intent(in) :: ngrid ! number of atmospheric columns
      integer,intent(in) :: nlayer ! number of atmospheric layers
      integer,intent(in) :: nq ! number of tracers

      allocate(zdqchim(ngrid,nlayer,nq))
      zdqchim(:,:,:)=0
      allocate(zdqschim(ngrid,nq))
      zdqschim(:,:)=0

    end subroutine ini_calchim_mod


    subroutine end_calchim_mod

      implicit none

      if (allocated(zdqchim))      deallocate(zdqchim)
      if (allocated(zdqschim))      deallocate(zdqschim)

    end subroutine end_calchim_mod

END MODULE calchim_mod 

