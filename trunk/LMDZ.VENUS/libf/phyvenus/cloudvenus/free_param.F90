MODULE free_param
  
  implicit none
  
  ! TIME (s) 
!  real, parameter :: TSTOP = 3 !Stop
  real, parameter :: TSTOP = 3600 !Stop
!  real, parameter :: TSTOP = 86400.0D0 ! Stop

  integer, parameter :: nbin = 100


  ! FLAGS
  logical, parameter :: MERGE     = .FALSE. ! Mode merging
  logical, parameter :: NUCLEA    = .TRUE.  ! Homogeneous nucleation
  logical, parameter :: MASSFLUX  = .TRUE.  ! Thermodyn. equi. and cond./evap.
  logical, parameter :: HETNUCLEA = .FALSE.  ! Heterogeneous nucleation
  logical, parameter :: SEDIM     = .FALSE.  ! Sedimentation
  logical, parameter :: COAG      = .FALSE. ! Coagulation
  integer, parameter :: PROD      = 0       ! No aer production
!!  integer, parameter :: PROD      = 1       ! Aer production with normal distri
!  integer, parameter :: PROD      = 2       ! Production only in 1 layer
  

  ! MODE 1 & 2 - James 1997 + Knollenberg and Hunten 1980
  ! Be careful, you need to change sigma in sed_and_prod.f90 too
  real, parameter :: sigma1 = 1.56D0  ! Geometric standard deviation
  real, parameter :: ri1 = 3.0D-7     ! Geometric radius (m)
  real, parameter :: sigma2 = 1.29D0  ! Geometric standard deviation
  real, parameter :: ri2 = 1.0D-6     ! Geometric radius (m)  

  
  ! RADIUS GRID
  real, parameter ::  rmin = 0.0001D-6 !m
  real, parameter ::  rmax = 100.0D-6  !m

  
  ! PRODUCTION - AEROSOL     NB:  The source should be at more or less 55km !
  real, parameter :: r_aer = 0.125D-6      !m
  real, parameter :: p_aer = 310000.D0     !Pa, 40km on Venus
!  real, parameter :: p_aer = 260000.D0     !Pa, 40km on Venus
  real, parameter :: sig_aer = 1.56D0      !Of the distribution, same as mode 1
  real, save      :: ntot_aer              !#/m-3
  real, parameter :: ntot_aer_flag2 = 4.D5 !#/m-3
  real, parameter :: tx_prod = 9.D-7       !kg.m-2.s-1, JB:tx_prod = 3.5D-13
  real, parameter :: rho_aer = 2000.D0     !kg/m-3
  REAL, PARAMETER :: sigz  = 20.0D3        ! currently set to 20km


  ! HETEROGENEOUS NUCLEATION
  real, parameter :: mtetalocal = 0.952D0 ! moullabilite


END MODULE free_param
