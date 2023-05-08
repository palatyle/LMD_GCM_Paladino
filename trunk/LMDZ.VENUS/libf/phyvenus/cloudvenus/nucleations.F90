
!  SUBROUTINE    HONUMO        HOmogeneous NUcleation with MOments
!  SUBROUTINE    HETNUMO       HETerogeneous NUcleation with MOments
!  SUBROUTINE    newnuklefit   Calculates jnuc and rc
!  SUBROUTINE    heteronuk     Calculates the sum of activated aerosols


!****************************************************************************
! HOmogeneous NUcleation with MOments
!
! This subroutine uses the newnuclefit function. 
! Hanna Vehkamäki et al., 2002 : An improved parameterization for 
! sulfuric acid/water nucleation rates for tropospheric 
! and stratospheric conditions, () J. Geophys. Res., 107, pp. 4622-4631
!
! we determine here the tendency of M0
! the created paticles will be in mode 1 (the smallest particles).
! 
! The inputs mk0 and mk3 are not use in this function.
!
! An implicit method are used (easy)
!
! by S.Guilbon (LATMOS), 2016.
!
SUBROUTINE HONUMO(TAIR,dt,mk0,mk3,dM0_homo,dM3_homo)

  use free_param
  use donnees

  IMPLICIT NONE

  ! INputs
  real, intent(in), dimension(2) :: mk0   ! first moment
  real, intent(in), dimension(3) :: mk3   ! third moment
  real, intent(in)  :: TAIR, dt           ! temperature, timestep
  ! Outputs 
  real, intent(out) :: dM0_homo, dM3_homo ! Variations of moments

  !Local variables
  real :: jnuc, rc        ! Nucleation rate and critical radius
  real :: m_H2SO4, m_H2O  ! Total condensed masses
  real :: CWV, ntotal

  ! Calculate jnuc and rc
  CALL newnuklefit(TAIR,h2so4_m3,jnuc,rc,ntotal)     

  ! initialization
  dM0_homo = 0.D0
  dM3_homo = 0.D0

  ! Conditions to have nucleation process
  IF (jnuc .ge. 1.D-1 .and. jnuc .le. 1.D16 .and. ntotal .gt. 4.0) then
     CWV = MAIR/MWV
     rc = 1.D-9  !m  Because of the calculation of WSA

     ! Modal tendancies: Total nb density and vol. of created particles
     ! Same thing with implicit scheme (cf thesis)
     dM0_homo = dt * jnuc          !#.m(air)-3 
     dM3_homo = dt * jnuc * rc**3  !m3.m(air)-3
  ENDIF

  RETURN 

END SUBROUTINE HONUMO


!****************************************************************************
! HETerogeneous NUcleation with MOment
!
! Here, we use a simple representation of nucleation rate 
! and heterogeneous nucleation, like in VenLA (see Maattanen Anni).
! 
! The number of activated ccn will become the number of created droplets
! ra is the mean radius of the aerosol distribution.
!
! for mk0 calculation, method 1 and method 2 give same results :)
!  mk0 = Jhetsum ! Number of activated CCN, in # - (method 1)
! Here, we use implicit scheme (cf thesis)
!
SUBROUTINE HETNUMO(TAIR,dt,ra,mk0aer,mk3aer,mk0drop, mk3drop, &
     & dM0_aer,dM3_aer,dM0_het,dM3_het)

  use free_param
  use donnees
  
  IMPLICIT NONE
  
  ! Inputs and Outputs 
  real, intent(in)  :: dt,TAIR,ra           ! time step, temperature, mean radius
  real, intent(in)  :: mk0aer, mk3aer       ! original moments M(t)
  real, intent(in), dimension(2) :: mk0drop ! original moment M(t)
  real, intent(in), dimension(3) :: mk3drop ! original moment M(t)
  real, intent(out) :: dM0_het, dM3_het     ! tendencies of liquid droplets
  real, intent(out) :: dM0_aer, dM3_aer     ! tendencies of aerosols
  ! Local variables
  real :: Jhetsum, alpha_loc1, alpha_loc2
  real :: mk0, mk3, alpha_k

  call build_radius_grid(nbin,rmin,rmax,vratio)
  call heteronuk(TAIR,ra,mk0aer,Jhetsum)

  write(*,*)'mk0aer',mk0aer

  alpha_loc1 = alpha_k(2,sig_aer) 
  mk0 = 1.D0 + (4.D0*PI*Jhetsum*r_aer*r_aer*alpha_loc1)*dt
  mk0 = (1.D0/mk0)   
  mk0 = mk0 * mk0aer ! - (method 2)
  mk0 = mk0 - mk0aer ! Delta betwen 2 time steps = tendancy

  alpha_loc2 = alpha_k(5,sig_aer)/alpha_k(3,sig_aer)
  mk3 = 1.D0 + (4.D0*PI*Jhetsum*r_aer*r_aer*alpha_loc2)*dt
  mk3 = 1.D0/mk3
  mk3 = mk3 * mk3aer
  mk3 = mk3 - mk3aer ! Delta betwen 2 time steps = tendancy
  write(*,*)'HET - mk3 - methode 2', mk3

  ! Modal tendancies
  ! the lost quantity in aerosol will go to droplets.
  dM0_aer = mk0
  dM3_aer = mk3
  dM0_het = mk0  * (-1)      
  dM3_het = mk3  * (-1)

END SUBROUTINE HETNUMO


!*****************************************************************************
SUBROUTINE newnuklefit(t,rhoa,jnuc,rc,ntotal)

  !    Fortran 90 subroutine binapara
  !
  !    Calculates parametrized values of particle formation rate,
  !    mole fraction of sulphuric acid,
  !    total number of particles, and the radius of the critical cluster
  !    in H2O-H2SO4 system if temperature, saturatio ratio of water and 
  !    sulfuric acid concentration  are given. 
  !
  !    Calculates also the kinetic limit and the particle formation rate
  !    above this limit (in which case we set ntotal=1 and x=1)
  !
  !    NB: produces nucleation rate   1/(cm^3 s) as a function of rh and t
  !    of the kinetic limit as a function of rh and t
  !
  !    Copyright (C)2016 Määttänen et al. 2016
  !
  !    anni.maattanen@latmos.ipsl.fr
  !    joonas.merikanto@fmi.fi
  !    hanna.vehkamaki@helsinki.fi
  !
  !    modified by S.Guilbon for the Venus GCM (2016)
  !
  ! SG: change the unity of jnuc here
  ! SG: change inputs and outputs

  use free_param
  use donnees

  ! Inputs
  real,intent(inout) :: t        ! temperature in K (165-400 K)
  real,intent(inout) :: rhoa     ! sulfuric acid concentration in 1/m3 (10^10 - 10^19) 

  ! Outputs
  real,intent(out) :: jnuc       ! nucleation rate in 1/m3s (J>10^-1 1/m3s)
  real,intent(out) :: rc         ! radius of the critical cluster in m 
  real,intent(out) :: ntotal       ! total number of molecules in the critical cluster

  ! Local variables
  real :: x          ! mole fraction of H2SO4 in the critical cluster     
  real :: nac        ! local variable for number of acid molecules

  ! Initialization
  jnuc = 0.D0
  ntotal = 0.D0
  rc = 0.D0
  x = 0.D0   

  if(t < 165.D0) then
     print *,'Warning: temperature < 165.0 K, using 165.0 K'
     t=165.D0
  end if

  if(t > 400.D0) then
     print *,'Warning: temperature > 400. K, using 400. K'
     t=400.D0
  end if

  if(rh < 1.D-5) then
     print *,'Warning: saturation ratio of water < 1.e-5, using 1.e-5'
     rh=1.D-5
  end if

  if(rh > 1.D0) then
     print *,'Warning: saturation ratio of water > 1 using 1'
     rh=1.D0
  end if

  if(rhoa < 1.D10) then
     print *,'Warning: sulfuric acid concentration < 1e4 1/cm3, using 1e4 1/cm3'
     rhoa=1.D10
  end if

  if(rhoa > 1.D19) then
     print *,'Warning: sulfuric acid concentration > 1e13 1/cm3, using 1e13 1/cm3'
     rhoa=1.D19
  end if

  x=  7.9036365428891719D-1 - 2.8414059650092153D-3*t + 1.4976802556584141D-2*LOG(rh) &
       & - 2.4511581740839115D-4*t*LOG(rh) + 3.4319869471066424D-3 *LOG(rh)**2     &  
       & - 2.8799393617748428D-5*t*LOG(rh)**2 + 3.0174314126331765D-4*LOG(rh)**3 & 
       & - 2.2673492408841294D-6*t*LOG(rh)**3 - 4.3948464567032377D-3*LOG(rhoa/1.D6)&
       & + 5.3305314722492146D-5*t*LOG(rhoa/1.D6)


  jnuc= 2.1361182605986115D-1 + 3.3827029855551838D0 *t -3.2423555796175563D-2*t**2 +  &
       &  7.0120069477221989D-5*t**3 + 8.0286874752695141D0/x +  &
       &  -2.6939840579762231D-1*LOG(rh) + 1.6079879299099518D0*t*LOG(rh) +  &
       &  -1.9667486968141933D-2*t**2*LOG(rh) +  &
       &  5.5244755979770844D-5*t**3*LOG(rh) + (7.8884704837892468D0*LOG(rh))/x +  &
       &  4.6374659198909596D0*LOG(rh)**2 - 8.2002809894792153D-2*t*LOG(rh)**2 +  &
       &  8.5077424451172196D-4*t**2*LOG(rh)**2 +  &
       &  -2.6518510168987462D-6*t**3*LOG(rh)**2 +  &
       &  (-1.4625482500575278D0*LOG(rh)**2)/x - 5.2413002989192037D-1*LOG(rh)**3 +  &
       &  5.2755117653715865D-3*t*LOG(rh)**3 +  &
       &  -2.9491061332113830D-6*t**2*LOG(rh)**3 +  &
       &  -2.4815454194486752D-8*t**3*LOG(rh)**3 +  &
       &  (-5.2663760117394626D-2*LOG(rh)**3)/x +  &
       &  1.6496664658266762D0*LOG(rhoa/1.D6) +  &
       &  -8.0809397859218401D-1*t*LOG(rhoa/1.D6) +  &
       &  8.9302927091946642D-3*t**2*LOG(rhoa/1.D6) +  &
       &  -1.9583649496497497D-5*t**3*LOG(rhoa/1.D6) +  &
       &  (-8.9505572676891685D0*LOG(rhoa/1.D6))/x +  &
       &  -3.0025283601622881D+1*LOG(rh)*LOG(rhoa/1.D6) +  &
       &  3.0783365644763633D-1*t*LOG(rh)*LOG(rhoa/1.D6) +  &
       &  -7.4521756337984706D-4*t**2*LOG(rh)*LOG(rhoa/1.D6) +  &
       &  -5.7651433870681853D-7*t**3*LOG(rh)*LOG(rhoa/1.D6) +  &
       &  (1.2872868529673207D0*LOG(rh)*LOG(rhoa/1.D6))/x +  &
       &  -6.1739867501526535D-1*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &  7.2347385705333975D-3*t*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &  -3.0640494530822439D-5*t**2*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &  6.5944609194346214D-8*t**3*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &  (-2.8681650332461055D-2*LOG(rh)**2*LOG(rhoa/1.D6))/x +  &
       &  6.5213802375160306D0*LOG(rhoa/1.D6)**2 +  &
       &  -4.7907162004793016D-2*t*LOG(rhoa/1.D6)**2 +  &
       &  -1.0727890114215117D-4*t**2*LOG(rhoa/1.D6)**2 +  &
       &  5.6401818280534507D-7*t**3*LOG(rhoa/1.D6)**2 +  &
       &  (5.4113070888923009D-1*LOG(rhoa/1.D6)**2)/x +  &
       &  5.2062808476476330D-1*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &  -6.0696882500824584D-3*t*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &  2.3851383302608477D-5*t**2*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &  -1.5243837103067096D-8*t**3*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &  (-5.6543192378015687D-2*LOG(rh)*LOG(rhoa/1.D6)**2)/x +  &
       &  -1.1630806410696815D-1*LOG(rhoa/1.D6)**3 +  &
       &  1.3806404273119610D-3*t*LOG(rhoa/1.D6)**3 +  &
       &  -2.0199865087650833D-6*t**2*LOG(rhoa/1.D6)**3 +  &
       &  -3.0200284885763192D-9*t**3*LOG(rhoa/1.D6)**3 +  &
       &  (-6.9425267104126316D-3*LOG(rhoa/1.D6)**3)/x

  ntotal =-3.5863435141979573D-3 - 1.0098670235841110D-1 *t + 8.9741268319259721D-4 *t**2 - 1.4855098605195757D-6*t**3  &
       &   - 1.2080330016937095D-1 /x + 1.1902674923928015D-3*LOG(rh) - 1.9211358507172177D-2*t*LOG(rh) +  &
       &   2.4648094311204255D-4*t**2*LOG(rh) - 7.5641448594711666D-7*t**3*LOG(rh) +  &
       &   (-2.0668639384228818D-02*LOG(rh))/x - 3.7593072011595188D-2*LOG(rh)**2 + 9.0993182774415718D-4 *t*LOG(rh)**2 +  &
       &   -9.5698412164297149D-6*t**2*LOG(rh)**2 + 3.7163166416110421D-8*t**3*LOG(rh)**2 +  &
       &   (1.1026579525210847D-2*LOG(rh)**2)/x + 1.1530844115561925D-2 *LOG(rh)**3 +  &
       &   - 1.8083253906466668D-4 *t*LOG(rh)**3 + 8.0213604053330654D-7*t**2*LOG(rh)**3 +  &
       &   -8.5797885383051337D-10*t**3*LOG(rh)**3 + (1.0243693899717402D-3*LOG(rh)**3)/x +  &
       &   -1.7248695296299649D-2*LOG(rhoa/1.D6) + 1.1294004162437157D-2*t*LOG(rhoa/1.D6) +  &
       &   -1.2283640163189278D-4*t**2*LOG(rhoa/1.D6) + 2.7391732258259009D-7*t**3*LOG(rhoa/1.D6) +  &
       &   (6.8505583974029602D-2*LOG(rhoa/1.D6))/x +2.9750968179523635D-1*LOG(rh)*LOG(rhoa/1.D6) +  &
       &   -3.6681154503992296D-3 *t*LOG(rh)*LOG(rhoa/1.D6) + 1.0636473034653114D-5*t**2*LOG(rh)*LOG(rhoa/1.D6) +  &
       &   5.8687098466515866D-9*t**3*LOG(rh)*LOG(rhoa/1.D6) + (-5.2028866094191509D-3*LOG(rh)*LOG(rhoa/1.D6))/x +  &
       &   7.6971988880587231D-4*LOG(rh)**2*LOG(rhoa/1.D6) - 2.4605575820433763D-5*t*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &   2.3818484400893008D-7*t**2*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &   -8.8474102392445200D-10*t**3*LOG(rh)**2*LOG(rhoa/1.D6) +  &
       &   (-1.6640566678168968D-4*LOG(rh)**2*LOG(rhoa/1.D6))/x - 7.7390093776705471D-2*LOG(rhoa/1.D6)**2 +  &
       &   5.8220163188828482D-4*t*LOG(rhoa/1.D6)**2 + 1.2291679321523287D-6*t**2*LOG(rhoa/1.D6)**2 +  &
       &   -7.4690997508075749D-9*t**3*LOG(rhoa/1.D6)**2 + (-5.6357941220497648D-3*LOG(rhoa/1.D6)**2)/x +  &
       &   -4.7170109625089768D-3*LOG(rh)*LOG(rhoa/1.D6)**2 + 6.9828868534370193D-5*t*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &   -3.1738912157036403D-7*t**2*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &   2.3975538706787416D-10*t**3*LOG(rh)*LOG(rhoa/1.D6)**2 +  &
       &   (4.2304213386288567D-4*LOG(rh)*LOG(rhoa/1.D6)**2)/x + 1.3696520973423231D-3*LOG(rhoa/1.D6)**3 +  &
       &   -1.6863387574788199D-5*t*LOG(rhoa/1.D6)**3 + 2.7959499278844516D-8*t**2*LOG(rhoa/1.D6)**3 +  &
       &   3.9423927013227455D-11*t**3*LOG(rhoa/1.D6)**3 + (8.6136359966337272D-5*LOG(rhoa/1.D6)**3)/x

  ntotal=EXP(ntotal)

  rc=EXP(-22.378268374023630D0 + 0.44462953606125100D0*x + 0.33499495707849131D0*LOG(ntotal))

  jnuc=EXP(jnuc)

  if(jnuc < 1.D-7) then
     print *,'Warning: nucleation rate < 1e-7/cm3s, using 0.0/cm3s,'
     jnuc = 0.0D0
  endif

  jnuc = jnuc*1.D6 !1/(m^3s)

  nac =x*ntotal

  if (nac .lt. 1.D0) then
     print *, 'Warning: number of acid molecules < 1 in nucleation regime, setting x=1'
     x=1.D0
  endif

  RETURN

END SUBROUTINE newnuklefit


!*****************************************************************************
SUBROUTINE heteronuk(TAIR,ra,mk0i,ACTOT)

  use donnees
  use free_param
  IMPLICIT NONE 

  ! inputs/outputs
  REAL, intent(in) :: TAIR, mk0i, ra
  REAL, intent(out) :: ACTOT

  ! Local variables
  INTEGER :: i
  INTEGER :: AER
  REAL :: NDnuk(nbin,2)
  REAL :: PSASAS

  ! Initialization
  AER = 2
  ntot_aer = mk0i
  ACTOT=0.D0 

  call logdist(sig_aer,ntot_aer,ra,rad_cld,NDnuk(:,AER)) ! NDnuk in #

  do i=1,nbin
     ! Saturation ratio with curvature effect
     PSASAS = RHOsasat * EXP((2.D0*MSA*ST)/(RGAS*TAIR*RHOSA*rad_cld(i))) 
     SH2SO4 = PPSA / PSASAS
   
     if (SH2SO4 .ge. 1.D0 .and. NDnuk(i,AER) .gt. 0.) then
        ! ACTOT = ACTOT + NDnuk(i,AER)/(vratio * rad_cld(i))
        ACTOT = ACTOT + NDnuk(i,AER) ! Total AER number, #
        NDnuk(i,AER) = 0.D0          ! Now, it's empty !
     endif
  enddo
END SUBROUTINE heteronuk
