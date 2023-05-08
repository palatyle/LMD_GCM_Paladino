











       module tracer_h

       implicit none

       integer, save :: nqtot ! total number of tracers
!$OMP THREADPRIVATE(nqtot)

       character*30, save, allocatable :: noms(:)   ! name of the tracer
       real, save, allocatable :: mmol(:)     ! mole mass of tracer (g/mol-1) 
       real, save, allocatable :: radius(:)   ! dust and ice particle radius (m)
       real, save, allocatable :: rho_q(:)    ! tracer densities (kg.m-3)
       real, save, allocatable :: qext(:)     ! Single Scat. Extinction coeff at 0.67 um
       real, save, allocatable :: alpha_lift(:)  ! saltation vertical flux/horiz flux ratio (m-1)
       real, save, allocatable :: alpha_devil(:) ! lifting coeeficient by dust devil
       real, save, allocatable :: qextrhor(:) ! Intermediate for computing opt. depth from q

       real,save :: varian      ! Characteristic variance of log-normal distribution
       real,save :: r3n_q     ! used to compute r0 from number and mass mixing ratio
       real,save :: rho_dust     ! Mars dust density (kg.m-3)
       real,save :: rho_ice     ! Water ice density (kg.m-3)
       real,save :: rho_co2     ! CO2 ice density (kg.m-3)
       ! real,save :: rho_volc    ! Volcanic ash density (kg.m-3)
       real,save :: rho_h2so4   ! Liquid sulfate aerosol density (kg.m-3)
       real,save :: ref_r0        ! for computing reff=ref_r0*r0 (in log.n. distribution)
!$OMP THREADPRIVATE(noms,mmol,radius,rho_q,qext,alpha_lift,alpha_devil,qextrhor, &
	!$OMP varian,r3n_q,rho_dust,rho_ice,rho_co2,ref_r0)

! tracer indexes: these are initialized in initracer and should be 0 if the
!                 corresponding tracer does not exist
       ! dust
       integer,save,allocatable :: igcm_dustbin(:) ! for dustbin 'dust' tracers
       ! dust, special doubleq case
       integer,save :: igcm_dust_mass   ! dust mass mixing ratio (for transported dust)
       integer,save :: igcm_dust_number ! dust number mixing ratio (transported dust)
       ! water
       integer,save :: igcm_h2o_vap ! water vapour
       integer,save :: igcm_h2o_ice ! water ice
       ! chemistry:
       integer,save :: igcm_co2
       integer,save :: igcm_co
       integer,save :: igcm_o
       integer,save :: igcm_o1d
       integer,save :: igcm_o2
       integer,save :: igcm_o3
       integer,save :: igcm_h
       integer,save :: igcm_h2
       integer,save :: igcm_oh
       integer,save :: igcm_ho2
       integer,save :: igcm_h2o2
       integer,save :: igcm_n2
       integer,save :: igcm_ar
       integer,save :: igcm_n
       integer,save :: igcm_no
       integer,save :: igcm_no2
       integer,save :: igcm_n2d
       integer,save :: igcm_ch4

       integer,save :: igcm_ch3
       integer,save :: igcm_ch
       integer,save :: igcm_3ch2
       integer,save :: igcm_1ch2
       integer,save :: igcm_cho
       integer,save :: igcm_ch2o
       integer,save :: igcm_ch3o
       integer,save :: igcm_c
       integer,save :: igcm_c2
       integer,save :: igcm_c2h
       integer,save :: igcm_c2h2
       integer,save :: igcm_c2h3
       integer,save :: igcm_c2h4
       integer,save :: igcm_c2h6
       integer,save :: igcm_ch2co
       integer,save :: igcm_ch3co
       integer,save :: igcm_hcaer
       ! Volcano
       integer,save :: igcm_volc_1=0
       integer,save :: igcm_volc_2=0
       integer,save :: igcm_volc_3=0
       integer,save :: igcm_volc_4=0
       integer,save :: igcm_volc_5=0
       integer,save :: igcm_volc_6=0
       integer,save :: igcm_h2so4=0

       ! other tracers
       integer,save :: igcm_ar_n2 ! for simulations using co2 +neutral gaz
       integer,save :: igcm_co2_ice ! CO2 ice 
!$OMP THREADPRIVATE(igcm_dustbin,igcm_dust_mass,igcm_dust_number,igcm_h2o_vap,igcm_h2o_ice, &
	!$OMP igcm_co2,igcm_co,igcm_o,igcm_o1d,igcm_o2,igcm_o3,igcm_h,igcm_h2,igcm_oh,	    &
	!$OMP igcm_ho2,igcm_h2o2,igcm_n2,igcm_ar,igcm_ar_n2,igcm_co2_ice,                   &
        !$OMP igcm_n,igcm_no,igcm_no2,igcm_n2d,igcm_ch4,igcm_ch3,igcm_ch,igcm_3ch2,         &
        !$OMP igcm_1ch2,igcm_cho,igcm_ch2o,igcm_ch3o,igcm_c,igcm_c2,igcm_c2h,igcm_c2h2,     &
        !$OMP igcm_c2h3,igcm_c2h4,igcm_c2h6,igcm_ch2co,igcm_ch3co,igcm_hcaer)

       end module tracer_h

