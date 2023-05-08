MODULE comchem_h

! -------------------------------------------------------------------------------------------------------------
! Purpose : + Stores data relative to :  * 1. Chemistry in the GCM
! -------                                * 2. Upper chemistry pressure grid
!                                        * 3. Coupling with C photochem. module ( cf calchim.F90)
!           
!           + Also contains a routine of initialization for chemistry in the GCM.
!
!           + NB : For newstart there is a specific comchem_newstart_h module.
!
! Author : Jan Vatant d'Ollone (2017-18)
! ------
!
! NB : A given order is assumed for the 44 chemistry tracers :
!      H, H2, CH, CH2s, CH2, CH3, CH4, C2, C2H, C2H2, C2H3, C2H4, C2H5, 
!      C2H6, C3H3, C3H5, C3H6, C3H7, C4H, C4H3, C4H4, C4H2s, CH2CCH2,
!      CH3CCH, C3H8, C4H2, C4H6, C4H10, AC6H6, C3H2, C4H5, AC6H5, N2,
!      N4S, CN, HCN, H2CN, CHCN, CH2CN, CH3CN, C3N, HC3N, NCCN, C4N2
!
!        
! IMPORTANT : Advected chem. tracers are in MASS fraction but upper fields ykim_up are in MOLAR fraction !
! 
! -------------------------------------------------------------------------------------------------------------

IMPLICIT NONE  

   ! ~~~~~~~~~~~~~~~~~~~~~~~~
   ! 1. Chemistry in the GCM
   ! ~~~~~~~~~~~~~~~~~~~~~~~~

   !! Hard-coded number of chemical species for Titan chemistry
   INTEGER, PARAMETER :: nkim = 44

   !! Hard-coded chemical species for Titan chemistry
   CHARACTER(len=10), DIMENSION(nkim), PARAMETER  :: cnames = &
     (/"H         ", "H2        ", "CH        ", "CH2s      ", "CH2       ", "CH3       ", &
       "CH4       ", "C2        ", "C2H       ", "C2H2      ", "C2H3      ", "C2H4      ", &
       "C2H5      ", "C2H6      ", "C3H3      ", "C3H5      ", "C3H6      ", "C3H7      ", & 
       "C4H       ", "C4H3      ", "C4H4      ", "C4H2s     ", "CH2CCH2   ", "CH3CCH    ", &
       "C3H8      ", "C4H2      ", "C4H6      ", "C4H10     ", "AC6H6     ", "C3H2      ", &
       "C4H5      ", "AC6H5     ", "N2        ", "N4S       ", "CN        ", "HCN       ", &
       "H2CN      ", "CHCN      ", "CH2CN     ", "CH3CN     ", "C3N       ", "HC3N      ", &
       "NCCN      ", "C4N2      "/)
   !! Hard-coded chemical species for Titan chemistry + "HV" specie for the photochem module.
   CHARACTER(len=10), DIMENSION(nkim+1)  :: nomqy_c ! Initialized in calchim with null terminator
   !! Hard-coded chemical species molar mass (g.mol-1), shares the same indexing than cnames.
   REAL, DIMENSION(nkim), PARAMETER               :: cmmol = (/ &
       1.01   , 2.0158, 13.02, 14.03, 14.03, 15.03, 16.04  , 24.02, 25.03, 26.04  , 27.05  , &
       28.05  , 29.06 , 30.07, 39.06, 41.07, 42.08, 43.09  , 49.05, 51.07, 52.08  , 50.06  , &
       40.07  , 40.07 , 44.11, 50.06, 54.09, 58.13, 78.1136, 38.05, 53.07, 77.1136, 28.0134, &
       14.01  , 26.02 , 27.04, 28.05, 39.05, 40.04, 41.05  , 50.04, 51.05, 52.04  , 76.1   /)

   !! Hard-coded molar fraction of surface methane
   REAL, PARAMETER :: botCH4 = 0.0565 ! From Niemann et al. 2010 - Huygens GCMS measurements
   


   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
   !  2. Upper chemistry grid 
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~

   INTEGER, SAVE :: nlaykim_up   ! Number of upper atm. layers for chemistry from GCM top to 4.5E-5 Pa (1300km)
   INTEGER, SAVE :: nlaykim_tot  ! Number of total layers for chemistry from surface to 4.5E-5 Pa (1300km)
!$OMP_THREADPRIVATE(nlaykim_up,nlay_kim_tot)

   ! NB : For the startfile we use nlaykim_up grid (upper atm) and for outputs we use nlaykim_tot grid (all layers)
  
   REAL*8, PARAMETER :: grkim_dz = 10.0 ! Vertical discretization of the upper chemistry grid (km)

   REAL,SAVE,ALLOCATABLE,DIMENSION(:)     :: preskim  ! Pressure (Pa) of upper chemistry (mid)-layers
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: zlaykim  ! Pseudo-altitude (km) of upper chemistry (mid)-layers
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: ykim_up  ! Upper chemistry fields (mol/mol)
   
   ! These "_tot" fields are for output only
   REAL,SAVE,ALLOCATABLE,DIMENSION(:)     :: preskim_tot  ! Pressure (Pa) of total chemistry (mid)-layers
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: zlaykim_tot  ! Pseudo-altitude (km) of total chemistry (mid)-layers
   REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: ykim_tot     ! Total chemistry fields (mol/mol)
   
!$OMP_THREADPRIVATE(preskim,zlaykim,ykim_up)
!$OMP_THREADPRIVATE(preskim_tot,zlaykim_tot,ykim_tot)

   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !  3. Interface with photochemical module (cf calchim.F90)
   ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   ! These 3 parameters as well as nkim above, MUST match titan.h in chimtitan !!
   INTEGER, PARAMETER :: nd_kim   = 54     ! Number of photodissociations
   INTEGER, PARAMETER :: nr_kim   = 377    ! Number of reactions in chemistry scheme
   INTEGER, PARAMETER :: nlrt_kim = 650    ! For the UV rad. transf., 650 levels of 2 km

   !! Hardcoded latitude discretisation for actinic fluxes - MUST be coherent with disso.c input files !!
   INTEGER, PARAMETER                         :: nlat_actfluxes = 49
   REAL, DIMENSION(nlat_actfluxes), PARAMETER :: lat_actfluxes  = (/ &
      90.0 ,  86.25,  82.5 ,  78.75,  75.0 ,  71.25,  67.5 ,  63.75,  60.0 ,  56.25,  52.5 ,  48.75, & 
      45.0 ,  41.25,  37.5 ,  33.75,  30.0 ,  26.25,  22.5 ,  18.75,  15.0 ,  11.25,   7.5 ,   3.75, &
       0.0 ,  -3.75,  -7.5 , -11.25, -15.0 , -18.75, -22.5 , -26.25, -30.0 , -33.75, -37.5 , -41.25, &
     -45.0 , -48.75, -52.5 , -56.25, -60.0 , -63.75, -67.5 , -71.25, -75.0 , -78.75, -82.5 , -86.25, &
     -90.0 /)
    

CONTAINS

  SUBROUTINE ini_comchem_h(ngrid)
  
  IMPLICIT NONE
  
    include "dimensions.h"
  
    INTEGER,INTENT(IN) :: ngrid ! number of atmospheric columns
  
    nlaykim_tot = nlaykim_up + llm
  
    IF (.NOT.allocated(preskim)) ALLOCATE(preskim(nlaykim_up))
    IF (.NOT.allocated(zlaykim)) ALLOCATE(zlaykim(ngrid,nlaykim_up))
    IF (.NOT.allocated(ykim_up)) ALLOCATE(ykim_up(nkim,ngrid,nlaykim_up))
    
    IF (.NOT.allocated(preskim_tot)) ALLOCATE(preskim_tot(nlaykim_tot))
    IF (.NOT.allocated(zlaykim_tot)) ALLOCATE(zlaykim_tot(ngrid,nlaykim_tot))
    IF (.NOT.allocated(ykim_tot)) ALLOCATE(ykim_tot(nkim,ngrid,nlaykim_tot))
  
  END SUBROUTINE ini_comchem_h


END MODULE comchem_h
