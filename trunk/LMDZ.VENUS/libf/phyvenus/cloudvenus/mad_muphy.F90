!-----------------------------------------------------------------------
! S.GUILBON - LATMOS - 2016.
!-----------------------------------------------------------------------
! PROGRAM MAD contains
!
!  SUBROUTINE   MAD_MUPHY        Microphysical processes with moments
!  SUBROUTINE   change_layer     Need if layer changes
!  SUBROUTINE   change_wsa       Need if wsa changes
!  SUBROUTINE   change_vapor     Need if number concentration of vapor changes
!  FUNCTION     logn
!  FUNCTION     moment

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
SUBROUTINE MAD_MUPHY(dt,TAIR,PAIR,MRWV_loc,MRSA_loc,M0_aeros,M3_aeros,  &
                           M0_m1drop,M0_m1ccn,M3_m1sa,M3_m1w,M3_m1ccn,  &
	                   M0_m2drop,M0_m2ccn,M3_m2sa,M3_m2w,M3_m2ccn)
  ! Loop on lat and long outside this subroutine => in phytrac_chimie
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  use free_param
  use donnees

  IMPLICIT NONE

  integer :: i, mode

  ! LAYER
  real, intent(in) :: TAIR, PAIR, dt 
  real, intent(inout) :: MRWV_loc, MRSA_loc

  ! MOMENTS
  real, intent(inout):: M0_aeros, M3_aeros
  real, intent(inout):: M0_m1drop,M0_m1ccn,M0_m2drop,M0_m2ccn
  real, intent(inout):: M3_m1sa,M3_m1w,M3_m1ccn,M3_m2sa,M3_m2w,M3_m2ccn

  ! MODAL TENDENCIES
  real, dimension(2) :: M0_m1, M0_m2 
  real, dimension(3) :: M3_m1, M3_m2
  real :: p, k                        ! Modal order
  real :: sum_M0, sum_M3              ! Sum of moment
  real :: dM0_merge_m1, dM3_merge_m1  ! Mode merging mode 1
  real :: dM0_merge_m2, dM3_merge_m2  ! Mode merging mode 2
  real :: dM0_flux_m1,  dM0_flux_m2   ! Masse flux, only M0, mode 1 and 2
  real :: dM3_flux_m1,  dM3_flux_m2   ! Masse flux, only M3, mode 1 and 2
  real :: dM0_aeros, dM3_aeros, ra    ! Heterogeneous nucleation, only aerosol
  real :: dM0_hom, dM3_hom            ! Homogeneous nucleation, only mode 1
  real :: dM0_het, dM3_het            ! Heterogeneous nucleation, only mode 1
  real :: dM0_m1_drop, dM0_m1_ccn, dM3_m1_SA, dM3_m1_WV, dM3_m1_ccn, dM0_aeros_m2
  real :: dM0_m2_drop, dM0_m2_ccn, dM3_m2_SA, dM3_m2_WV, dM3_m2_ccn, dM3_aeros_m2
  real :: dM3_m1, dM3_m2, dM3, dM3_SA1, dM3_SA2, dM3_WV1, dM3_WV2, dM3_SA, dM3_WV
  real :: frac_ccn_m1, frac_ccn_m2

  ! CONSERVATION Variables
  real :: var_mode1_M0, var_mode1_M3, var_mode2_M0, var_mode2_M3
  real :: var_MRSA_m1, var_MRWV_m1
  real :: var_MRSA_m2, var_MRWV_m2
  real :: var_MR_m1, var_MR_m2, var_MRSA, var_MRWV
  real :: factor, diff, check, MRTOT

  ! FUNCTIONS
  real :: moment, D0, alpha_k

  ! OUTPUTS and GRAPHES
  real, dimension(nbin) :: n_g1
  integer :: cpt

  !-------------

  ! Initialization of vectors
  M0_m1(1)=M0_m1drop
  M0_m1(2)=M0_m1ccn
  M0_m2(1)=M0_m2drop
  M0_m2(2)=M0_m2ccn
  M3_m1(1)=M3_m1sa
  M3_m1(2)=M3_m1w
  M3_m1(3)=M3_m1ccn
  M3_m2(1)=M3_m2sa
  M3_m2(2)=M3_m2w
  M3_m2(3)=M3_m2ccn
  
  !-------------

  ! Initialization of tendencies   
  dM0_hom = 0.D0      ;  dM3_hom = 0.D0
  dM0_het = 0.D0      ;  dM3_het = 0.D0 
  dM0_aeros = 0.D0    ;  dM3_aeros = 0.D0  
  dM0_flux_m1 = 0.D0  ;  dM0_flux_m2 = 0.D0
  dM3_flux_m1 = 0.D0  ;  dM3_flux_m2 = 0.D0
  dM0_merge_m1 = 0.D0 ;  dM3_merge_m1 = 0.D0 
  dM0_merge_m2 = 0.D0 ;  dM3_merge_m2 = 0.D0 

  !-------------

  ! To controle the conservation !
  ! Total mass (vapor + cond) of species (h2o + h2so4)
  sum_M3 = M3_m1(1) + M3_m1(2) + M3_m1(3) ! mode 1
  MRTOT = sum_M3/(exp(4.5D0*log(sigma1)))
  sum_M3 = M3_m2(1) + M3_m2(2) + M3_m2(3) ! mode 2
  MRTOT = MRTOT + sum_M3/(exp(4.5D0*log(sigma2)))
  MRTOT = MRTOT*(4.0D0/3.0D0)*PI*RHOSA    ! total condensed mass of SA
  MRTOT = MRTOT*WSA/(MSA*E) + MRTOT*(1.D0-WSA)/(MWV*E) + MRSA_loc + MRWV_loc ! in vapor quantity
  check = MRTOT

!-----------------------------------------------------------------------
!                    INITIALIZATION - THERMO
!-----------------------------------------------------------------------
  WSA = 0.01D0    ! Init. of the fraction of sulfuric acid [0.1;1]
  WSAEQ = 0.0D0   ! WSA at equilibrium

  call change_wsa(TAIR,PAIR,MRSA_loc)            ! Need if wsa changes
  call change_layer(TAIR,PAIR)                   ! Need if layer changes
  call change_vapor(TAIR,PAIR,MRSA_loc,MRWV_loc) ! Need if nb vapor concentration changes

  ! Number concentration of H2SO4 (m-3)
  h2so4_m3 = MRSA_loc*PAIR/(KBZ*TAIR) 

!-----------------------------------------------------------------------
!                           MODE MERGING
!----------------------------------------------------------------------- 
  IF (SEDIM .and. MERGE) THEN
     redge = (r1*r2)**0.5D0 ! Edge radius of mode 1
!     write(*,*)'The intersection size or virtual bonduary', redge

     IF (r1 .EQ. redge) THEN
!        write(*,*)'Mode 1 merge to mode 2', redge, r1
        CALL MERGING(r1,N1,sigma1,sigma2,dM0_merge_m1,dM3_merge_m1, &
             &   dM0_merge_m2,dM3_merge_m2)
        ! mode 1
        M0_m1(1) = dM0_merge_m1
        M0_m1(2) = dM0_merge_m1*frac_ccn_m1
        M3_m1(1) = dM3_merge_m1*(1.D0-frac_ccn_m1)*WSA
        M3_m1(2) = dM3_merge_m1*(1.D0-frac_ccn_m1)*(1.D0-WSA)
        M3_m1(3) = dM3_merge_m1*frac_ccn_m1
        ! mode 2
        M0_m2(1) = M0_m2(1) + dM0_merge_m2
        M0_m2(2) = M0_m2(2) + dM0_merge_m2*(frac_ccn_m2+frac_ccn_m1)
        M3_m2(1) = M3_m2(1) + dM3_merge_m2*(1.D0-(frac_ccn_m2+frac_ccn_m1))*WSA
        M3_m2(2) = M3_m2(2) + dM3_merge_m2*(1.D0-(frac_ccn_m2+frac_ccn_m1))*(1.D0-WSA)
        M3_m2(3) = M3_m2(3) + dM3_merge_m2*(frac_ccn_m2+frac_ccn_m1)
     ENDIF

     IF (r2 .EQ. redge) THEN
!        write(*,*)'Mode 2 merge to mode 1',redge,r2
        CALL MERGING(r2,N2,sigma2,sigma1,dM0_merge_m2,dM3_merge_m2, &
             &   dM0_merge_m1,dM3_merge_m1)
        ! mode 1
        M0_m2(1) = dM0_merge_m2
        M0_m2(2) = dM0_merge_m2*frac_ccn_m2
        M3_m2(1) = dM3_merge_m2*(1.D0-frac_ccn_m2)*WSA
        M3_m2(2) = dM3_merge_m2*(1.D0-frac_ccn_m2)*(1.D0-WSA)
        M3_m2(3) = dM3_merge_m2*frac_ccn_m2
        ! mode 2
        M0_m1(1) = M0_m1(1) + dM0_merge_m1
        M0_m1(2) = M0_m1(2) + dM0_merge_m1*(frac_ccn_m2+frac_ccn_m1)
        M3_m1(1) = M3_m1(1) + dM3_merge_m1*(1.D0-(frac_ccn_m2+frac_ccn_m1))*WSA
        M3_m1(2) = M3_m1(2) + dM3_merge_m1*(1.D0-(frac_ccn_m2+frac_ccn_m1))*(1.D0-WSA)
        M3_m1(3) = M3_m1(3) + dM3_merge_m1*(frac_ccn_m2+frac_ccn_m1)
     ENDIF
  ENDIF


!-----------------------------------------------------------------------
!                           NUCLEATION                       Only mode 1
!-----------------------------------------------------------------------
  IF(NUCLEA) then ! HOmogeneous NUcleation with MOment
     CALL HONUMO(TAIR,dt,M0_m1,M3_m1,dM0_hom,dM3_hom) !#.m(air)-3 and m3.m(air)-3
  ENDIF

  IF(HETNUCLEA) then ! HETerogeneous NUcleation with MOment
     ! Calculation of ra (mean radius of aerosols) before heterogeneous nucleation
     ! Needed ra to calculate the bin idstrib to determine Jhet.
     IF (M0_aeros .le. ZERO) THEN ! Mean radii MODE AER
        ra = 0.D0
     ELSE
        ra = (1.D0/alpha_k(3,sig_aer) * M3_aeros/M0_aeros)**(1.D0/3.D0)     
     ENDIF

     CALL HETNUMO(TAIR,dt,ra,M0_aeros,M3_aeros,M0_m1,M3_m1, &
          & dM0_aeros,dM3_aeros,dM0_het,dM3_het)
  ENDIF


!-----------------------------------------------------------------------
!                           MASSES FLUX                       Mode 1 & 2
!-----------------------------------------------------------------------
  ! We can have condensation/evaporation if we have some particles :)
  ! Conditions : flag + mode + gas
  IF (MASSFLUX .and. (M0_m1(1).gt.ZERO)) THEN ! Mode 1
     CALL FLUX(TAIR,PAIR,dt,sigma1,r1,M3_m1,dM0_flux_m1,dM3_flux_m1)
  ENDIF

  IF (MASSFLUX .and. (M0_m2(1).gt.ZERO))  THEN ! Mode 2
     CALL FLUX(TAIR,PAIR,dt,sigma2,r2,M3_m2,dM0_flux_m2,dM3_flux_m2)
  ENDIF


!-----------------------------------------------------------------------
!                           CONSERVATION                      Mode 1 & 2
!----------------------------------------------------------------------- 
  ! M3(1) sulfuric acid
  ! M3(2) water
  ! M3(3) condensed cloud nuclei
  !-------------------------------------!


  !--------AEROSOLS / HETERONUC---------!
  ! We just modify the tendencies. They will be add to moments.
  IF ((M0_aeros + dM0_aeros) .le. ZERO) THEN
     dM0_aeros = (-1.D0) * M0_aeros
     dM3_aeros = (-1.D0) * M3_aeros
     dM0_het = M0_aeros
     dM3_het = M3_aeros
  ENDIF
  !-------------------------------------!


  !----------SUM OF TENDENCIES----------!
  ! +5% of ccn vol to condensed mass to active the mass flux in next time step.
  ! Mode 1
  dM0_m1_drop = dM0_hom + dM0_het + dM0_flux_m1
  dM0_m1_ccn = dM0_het
  dM3_m1_SA = (dM3_hom + dM3_het + dM3_flux_m1) * WSA
  dM3_m1_WV = (dM3_hom + dM3_het + dM3_flux_m1) * (1.D0-WSA)
  dM3_m1_ccn = dM3_het
  ! Mode 2
  dM0_m2_drop = dM0_flux_m2
  dM0_m2_ccn = 0.D0
  dM3_m2_SA = dM3_flux_m2 * WSA
  dM3_m2_WV = dM3_flux_m2 * (1.D0-WSA)
  dM3_m2_ccn = 0.D0
  ! Aerosols and CCN
  dM0_aeros_m2 = 0.D0
  dM3_aeros_m2 = 0.D0
  !-------------------------------------!


  !--------TOO MUCH EVAPORATION---------!
  ! We just modify the tendencies. They will be add to moments.
  IF ((M3_m1(1) + dM3_m1_SA) .le. ZERO) THEN
     dM3_m1_SA = (-1.D0) * M3_m1(1) !-- vol acid
  ENDIF
  IF ((M3_m1(2) + dM3_m1_WV) .le. ZERO) THEN
     dM3_m1_WV = (-1.D0) * M3_m1(2) !-- vol water
  ENDIF
  IF ((M0_m1(1) + dM0_m1_drop) .le. ZERO) THEN
     dM0_m1_drop = (-1.D0)* M0_m1(1) !-- nb droplets
     dM0_m1_ccn = (-1.D0) * M0_m1(2) !-- nb ccn
     dM3_m1_SA  = (-1.D0) * M3_m1(1) !-- acid
     dM3_m1_WV  = (-1.D0) * M3_m1(2) !-- water
     dM3_m1_ccn = (-1.D0) * M3_m1(3) !-- ccn
     dM0_aeros = M0_m1(2)        !++ nb ccn
     dM3_aeros = M3_m1(3)        !++ vol ccn
  ENDIF
  ! Second mode of particles !
  IF ((M3_m2(1) + dM3_m2_SA) .le. ZERO) THEN
     dM3_m2_SA = (-1.D0) * M3_m2(1) !-- vol acid
  ENDIF
  IF ((M3_m2(2) + dM3_m2_WV) .le. ZERO) THEN
     dM3_m2_WV = (-1.D0) * M3_m2(2) !-- vol water
  ENDIF
  IF ((M0_m2(1) + dM0_m2_drop) .le. ZERO) THEN
     dM0_m2_drop = (-1.D0)* M0_m2(1) !-- nb droplets
     dM0_m2_ccn = (-1.D0) * M0_m2(2) !-- nb ccn
     dM3_m2_SA  = (-1.D0) * M3_m2(1) !-- acid
     dM3_m2_WV  = (-1.D0) * M3_m2(2) !-- water
     dM3_m2_ccn = (-1.D0) * M3_m2(3) !-- ccn
     dM0_aeros_m2 = M0_m2(2)         !++ nb ccn
     dM3_aeros_m2 = M3_m2(3)         !++ vol ccn
  ENDIF
  !-------------------------------------!


  !--------TOO MUCH CONDENSATION--------!
  ! SUM with new tendencies
  dM3_m1 = dM3_m1_SA + dM3_m1_WV
  dM3_m2 = dM3_m2_SA + dM3_m2_WV
  ! Vapor mode 1
  dM3 = dM3_m1 / (exp(4.5D0*log(sigma1)))
  dM3 = RHOSA * 4.D0/3.D0 * PI * dM3_m1 ! mass
  dM3_SA1 = (dM3_m1*WSA)/(MSA*E)                     ! Total consumed SA vapor 
  dM3_WV1 = (dM3_m1*(1.D0-WSA))/(MWV*E)              ! Total consumed water vapor 
  ! Vapor mode 2
  dM3 = dM3_m2/(exp(4.5D0*log(sigma2)))
  dM3 = RHOSA * 4.D0/3.D0 * PI * dM3_m2 ! mass
  dM3_SA2 = (dM3_m2*WSA)/(MSA*E)           ! Total consumed SA vapor 
  dM3_WV2 = (dM3_m2*(1.D0-WSA))/(MWV*E)    ! Total consumed water vapor 
  ! Sum with distinguished mode: not same process for mode 1 and mode 2
  dM3_SA = dM3_SA1 + dM3_SA2 ! h2so4
  dM3_WV = dM3_WV1 + dM3_WV2 ! h2o
  ! We just modify the tendencies. They will be add to moments.
  ! we define a sort of purcentage of participation for each moment.
  IF (((MRSA_loc + dM3_SA) .le. ZERO).and.(dM3_SA.ne.0.D0).and.(WSA.ne.0.D0)) THEN ! Comsuption of all H2SO4 vapor
     dM3_m1_SA = dM3_SA1/dM3_SA * MRSA_loc * MSA*E/WSA * 1.D0/(RHOSA*4.D0/3.D0*PI)*alpha_k(3,sigma1)
     dM3_m2_SA = dM3_SA2/dM3_SA * MRSA_loc * MSA*E/WSA * 1.D0/(RHOSA*4.D0/3.D0*PI)*alpha_k(3,sigma2)
     dM3_SA = dM3_SA1 + dM3_SA2
  ENDIF
  IF (((MRWV_loc + dM3_WV) .le. ZERO).and.(dM3_WV.ne.0.D0)) THEN  ! Comsuption of all H2O vapor
     dM3_m1_WV = dM3_WV1/dM3_WV * MRWV_loc * MWV*E/(1.D0-WSA) * 1.D0/(RHOSA*4.D0/3.D0*PI)*alpha_k(3,sigma1)
     dM3_m2_WV = dM3_WV2/dM3_WV * MRWV_loc * MWV*E/(1.D0-WSA) * 1.D0/(RHOSA*4.D0/3.D0*PI)*alpha_k(3,sigma2)
     dM3_WV = dM3_WV1 + dM3_WV2
  ENDIF
  !------------------------------------!


  !-------------UPDATES----------------!
  ! Updates Vapors !
  MRSA_loc = MRSA_loc + dM3_SA ! acid
  MRWV_loc = MRWV_loc + dM3_WV ! water
  ! Updtaes of CCN and Aerosols
  M0_aeros = M0_aeros + dM0_aeros +  dM0_aeros_m2     ! number of aerosols
  M3_aeros = M3_aeros + dM3_aeros +  dM3_aeros_m2     ! vol/mass of aerosols  
  ! Updates moments - Mode 1
  M0_m1(1) = M0_m1(1) + dM0_m1_drop ! number of droplets
  M0_m1(2) = M0_m1(2) + dM0_m1_ccn  ! number of CCN in droplets
  M3_m1(1) = M3_m1(1) + dM3_m1_SA   ! vol/mass of Héso4 in droplets
  M3_m1(2) = M3_m1(2) + dM3_m1_WV   ! vol/mass of h2o in droplets
  M3_m1(3) = M3_m1(3) + dM3_m1_ccn  ! vol/mass of CCN in droplets
  ! Updates moments - Mode 2
  M0_m2(1) = M0_m2(1) + dM0_m2_drop
  M0_m2(2) = M0_m2(2) + dM0_m2_ccn
  M3_m2(1) = M3_m2(1) + dM3_m2_SA
  M3_m2(2) = M3_m2(2) + dM3_m2_WV
  M3_m2(3) = M3_m2(3) + dM3_m2_ccn
  !------------------------------------!
!  write(*,*)'l.653: after up', M3_m1(1),M3_m1(2),M3_m1(3),dM3_m1_SA,dM3_m1_WV,dM3_m1_ccn
!  write(*,*)'l.654: after up', M0_m1(1),M0_m1(2),dM0_m1_drop
!-----------------------------------------------------------------------
  ! To controle/check the mass conservation !
  MRTOT = (M3_m1(1)+M3_m1(2))/(exp(4.5D0*log(sigma1)))         ! mode 1
  MRTOT = MRTOT + (M3_m2(1)+M3_m2(2))/(exp(4.5D0*log(sigma2))) ! mode 2
  MRTOT = MRTOT * (4.0D0/3.0D0)*PI*RHOSA                       ! mass 
  MRTOT = MRTOT*WSA/(MSA*E) + MRTOT*(1.D0-WSA)/(MWV*E)         ! binary solution !     
  MRTOT = MRTOT + MRSA_loc + MRWV_loc                          ! with vapor gas
  diff = check - MRTOT                                         ! difference ?

  if (diff .gt. ZERO) then
     write(*,*) 'MASS CHANGEs !!!', diff
     stop
  endif

  ! Updates
  call change_wsa(TAIR,PAIR,MRSA_loc)
  call change_vapor(TAIR,PAIR,MRSA_loc,MRWV_loc)


!-----------------------------------------------------------------------
!                           COAGULATION
!----------------------------------------------------------------------- 
!!$  IF (COAG) THEN
!!$     CALL drop_coagul(TAIR,PAIR,dt,M0_m1,M3_m1,M0_m2,M3_m2)
!!$  END IF


!-----------------------------------------------------------------------
!                        R_mean AND N_tot
!----------------------------------------------------------------------- 
  
  IF (M0_m1(1) .le. ZERO) THEN ! Mean radii MODE 1
     r1 = 0.D0
     N1 = 0.D0
     M0_m1(1) = 0.D0
     M0_aeros = M0_aeros + M0_m1(2)
     M0_m1(2) = 0.D0
     M3_m1(1) = 0.D0
     M3_m1(2) = 0.D0
     M3_aeros = M3_aeros + M3_m1(3)
     M3_m1(3) = 0.D0
  ELSE
     r1 = (1.D0/alpha_k(3,sigma1) * (M3_m1(1)+M3_m1(2)+M3_m1(3))/M0_m1(1))**(1.D0/3.D0)
     N1 = M0_m1(1)
     IF ((r1 .lt. ZERO) .or. (N1.lt.ZERO)) THEN
        r1 = 0.D0
        N1 = 0.D0
     ENDIF
!     write(*,*)'rmean not ZERO',r1,N1,M0_m1(1),M3_m1(1),WSA
  ENDIF
  
  IF (M0_m2(1) .le. ZERO) THEN ! Mean radii MODE 2
     r2 = 0.D0
     N2 = 0.D0
     M0_m2(1) = 0.D0
     M0_aeros = M0_aeros + M0_m2(2)
     M0_m2(2) = 0.D0
     M3_m2(1) = 0.D0
     M3_m2(2) = 0.D0
     M3_aeros = M3_aeros + M3_m2(3)
     M3_m2(3) = 0.D0
  ELSE
     r2 = (1.D0/alpha_k(3,sigma2) * (M3_m2(1)+M3_m2(2)+M3_m2(3))/M0_m2(1))**(1.D0/3.D0)
     N2 = M0_m2(1)
     IF ((r2 .lt. ZERO) .or. (N2.lt.ZERO)) THEN
        r2 = 0.D0
        N2 = 0.D0
     ENDIF
  ENDIF

  IF (M0_aeros .le. ZERO) THEN ! Mean radii MODE AER
     ra = 0.D0
  ELSE
     ra = (1.D0/alpha_k(3,sig_aer) * M3_aeros/M0_aeros)**(1.D0/3.D0)     
  ENDIF

  NTOT = N1 + N2
  ! CCN fraction in distribution => based on mode 1, like WSA ratio !
  frac_ccn_m1 = M0_m1(2)/(M0_m1(1)+M0_m1(2))
  frac_ccn_m2 = M0_m2(2)/(M0_m2(1)+M0_m2(2))

!  write(*,*)'rmean',r1,'NTOT',N1, 'SH2SO4', SH2SO4


!-----------------------------------------------------------------------
!                           MODE MERGING
!----------------------------------------------------------------------- 
  IF (MERGE) THEN
     redge = (r1*r2)**0.5D0 ! Edge radius of mode 1
!     write(*,*)'The intersection size or virtual bonduary', redge

     IF (r1 .EQ. redge) THEN
!        write(*,*)'Mode 1 merge to mode 2', redge, r1
        CALL MERGING(r1,N1,sigma1,sigma2,dM0_merge_m1,dM3_merge_m1, &
             &   dM0_merge_m2,dM3_merge_m2)
        ! mode 1
        M0_m1(1) = dM0_merge_m1
        M0_m1(2) = dM0_merge_m1*frac_ccn_m1
        M3_m1(1) = dM3_merge_m1*(1.D0-frac_ccn_m1)*WSA
        M3_m1(2) = dM3_merge_m1*(1.D0-frac_ccn_m1)*(1.D0-WSA)
        M3_m1(3) = dM3_merge_m1*frac_ccn_m1
        ! mode 2
        M0_m2(1) = M0_m2(1) + dM0_merge_m2
        M0_m2(2) = M0_m2(2) + dM0_merge_m2*(frac_ccn_m2+frac_ccn_m1)
        M3_m2(1) = M3_m2(1) + dM3_merge_m2*(1.D0-(frac_ccn_m2+frac_ccn_m1))*WSA
        M3_m2(2) = M3_m2(2) + dM3_merge_m2*(1.D0-(frac_ccn_m2+frac_ccn_m1))*(1.D0-WSA)
        M3_m2(3) = M3_m2(3) + dM3_merge_m2*(frac_ccn_m2+frac_ccn_m1)
     ENDIF

     IF (r2 .EQ. redge) THEN
!        write(*,*)'Mode 2 merge to mode 1',redge,r2
        CALL MERGING(r2,N2,sigma2,sigma1,dM0_merge_m2,dM3_merge_m2, &
             &   dM0_merge_m1,dM3_merge_m1)
        ! mode 2
        M0_m2(1) = dM0_merge_m2
        M0_m2(2) = dM0_merge_m2*frac_ccn_m2
        M3_m2(1) = dM3_merge_m2*(1.D0-frac_ccn_m2)*WSA
        M3_m2(2) = dM3_merge_m2*(1.D0-frac_ccn_m2)*(1.D0-WSA)
        M3_m2(3) = dM3_merge_m2*frac_ccn_m2
        ! mode 1
        M0_m1(1) = M0_m1(1) + dM0_merge_m1
        M0_m1(2) = M0_m1(2) + dM0_merge_m1*(frac_ccn_m2+frac_ccn_m1)
        M3_m1(1) = M3_m1(1) + dM3_merge_m1*(1.D0-(frac_ccn_m2+frac_ccn_m1))*WSA
        M3_m1(2) = M3_m1(2) + dM3_merge_m1*(1.D0-(frac_ccn_m2+frac_ccn_m1))*(1.D0-WSA)
        M3_m1(3) = M3_m1(3) + dM3_merge_m1*(frac_ccn_m2+frac_ccn_m1)
     ENDIF
  ENDIF

!-----------------------------------------------------------------------
!                        UPDATE OF OUTPUT VALUES
!----------------------------------------------------------------------- 

        ! mode 1
        M0_m1drop = M0_m1(1)
        M0_m1ccn  = M0_m1(2)
        M3_m1sa   = M3_m1(1)
        M3_m1w    = M3_m1(2)
        M3_m1ccn  = M3_m1(3)
        ! mode 2
        M0_m2drop = M0_m2(1)
        M0_m2ccn  = M0_m2(2)
        M3_m2sa   = M3_m2(1)
        M3_m2w    = M3_m2(2)
        M3_m2ccn  = M3_m2(3)

!-----------------------------------------------------------------------
!                          OUTPUTS => GRAPHES ! :D
!-----------------------------------------------------------------------
  if (dt .eq. 1) then
!     write(*,*) 'Ntot, rc mode 1: ', N1,r1
     ! Log-normal distribution function
     CALL build_radius_grid(nbin,rmin,rmax,vratio)
     CALL logdist(sigma1,N1,r1,rad_cld,n_g1)
!     DO i = 1, nbin
!        open(888,file='ND_t3600_nuc-het-flux_dt1_1step_sssat')
!        write(888,'(I10,2X,3(ES15.7,1X))') i, rad_cld(i), n_g1(i)
!     ENDDO
  endif


  RETURN
END SUBROUTINE MAD_MUPHY


!*****************************************************************************
!                              *************                                 !
!                               SUBROUTINES                                  !
!                              *************                                 !
!*****************************************************************************
SUBROUTINE change_layer(TAIR,PAIR)

  use free_param
  use donnees

  IMPLICIT NONE

  real, intent(in) :: TAIR, PAIR
  real :: FPLAIR, DFWVA, MFPAIR, CDTAIR

  MFPAIR = FPLAIR(TAIR,PAIR) ! Mean free path of air molecules
!  write(*,*) 'MFPAIR',MFPAIR, TAIR, PAIR
  Kn = MFPAIR/ri1             ! Knudsen number 

  D = DFWVA(TAIR,PAIR)      ! Molecular diffusivity of the air (m2/s)
  Akn = Kn*(1.333D0+0.71D0*(1.0D0/Kn)/(1.0D0+(1.0D0/Kn)))
  D = D/(1.0D0 + Akn)

  E = PAIR / (KBZ*TAIR)  ! Constant for microphys. processes

  KAIR = CDTAIR(TAIR)        ! Thermal conductivity of air (J/m sec K)

END SUBROUTINE change_layer


!*****************************************************************************
SUBROUTINE change_wsa(TAIR,PAIR,MRSA_loc)

  use free_param
  use donnees

  IMPLICIT NONE

  real, intent(in) :: TAIR, PAIR, MRSA_loc
  real :: ROSAS, PSASAS_ZELE,  STSAS
!  write(*,*)'checkpoint change_wsa', ST, RHOSA, RHOsasat, PSASAS_ZELE(TAIR,WSA), TAIR, WSA
  ST = STSAS(TAIR,WSA)      ! Surface tension of sulfuric acid solution/vapor (N/m)

  RHOSA = ROSAS(TAIR,WSA)          ! Density of liquid sulfuric acid solution (kg/m3)
  RHOsasat = PSASAS_ZELE(TAIR,WSA) ! Saturated sulfuric acid vapor pressure (Pa)

  SH2SO4 = (PAIR*MRSA_loc)/RHOsasat    ! Saturation Ratio

!  write(*,*)'checkpoint change_wsa', PSASAS_ZELE(TAIR,WSA)

END SUBROUTINE change_wsa


!*****************************************************************************
SUBROUTINE change_vapor(TAIR,PAIR,MRSA_loc,MRWV_loc)

  use free_param
  use donnees

  IMPLICIT NONE

  real, intent(in) :: TAIR, PAIR, MRSA_loc, MRWV_loc
  real :: waterps

  H2SO4_m3 = MRSA_loc*PAIR/(KBZ*TAIR) ! Number concentration of H2SO4 (m-3)

  PPWV = PAIR*MRWV_loc                   ! Partial pressure of water vapor (Pa or N/m2)
  PPSA = PAIR*MRSA_loc                   ! Partial Pressure of sulfuric acid (Pa or N/m2)

  RH = PPWV/waterps(TAIR) * 1.0D2    ! Relative humidity (%)

  CHECKWV = PPWV/PAIR
  CHECKSA = PPSA/PAIR

END SUBROUTINE change_vapor


!*****************************************************************************
FUNCTION logn(x,rm,sig)

  !*     Computes the log-normal distribution at x with parameters rm and sig
  !*
  !*     INPUTS
  !*          x    location
  !*          rm   mean radius
  !*          sig  standard deviation
  !*
  !*     OUTPUTS
  !*          The value of the log-normal distribution at location x

  use free_param
  use donnees
  IMPLICIT NONE

  real, intent(in) :: sig, rm, x
  real ::  cste, logn

  cste = 1.D0 / (sqrt(2.D0*PI)*log(sig))
  logn = cste / x * exp(-0.5D0*((log(x)-log(rm))/log(sig))**2)

  RETURN
END FUNCTION logn


!*****************************************************************************
FUNCTION moment(k, ntot, rbar, sigma) 

  !*     Computes the kth order moment of the log-normale with parameters : ntot,rbar and sigma
  !*
  !*     INPUTS
  !*            k     : the order of the moment to compute
  !*            ntot  : the total number of particles of the distribution (ie: M0)
  !*            rbar  : the mean raidus of the ditribution
  !*            sigma : the standard deviation of the distribution.

  IMPLICIT NONE

  real, intent(in) :: ntot, rbar, sigma, k 
  real :: moment

  moment = ntot * rbar ** k * exp(0.5D0 *(k*log(sigma))**2)

  RETURN
END FUNCTION moment


!*****************************************************************************
FUNCTION alpha_k(k,sig)

  IMPLICIT NONE

  real :: k, alpha_k
  real :: sig

  alpha_k = EXP(0.5 * k*k *(LOG(sig))**2)

  return
end FUNCTION alpha_k
