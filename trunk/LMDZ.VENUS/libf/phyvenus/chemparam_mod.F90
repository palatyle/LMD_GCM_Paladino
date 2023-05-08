MODULE chemparam_mod

!MODULE qui definit les indices des traceurs et leurs masses molaires.
! utilise aussi pour variables communes nuages/photochimie

IMPLICIT NONE

!----------------------------------------------------------------------------  
INTEGER, SAVE :: i_co2, i_co, i_h2, i_h2o, i_o1d,        &
                 i_o, i_o2, i_o2dg, i_o3, i_h,           &
                 i_oh, i_ho2, i_h2o2, i_cl, i_clo,       &
                 i_cl2, i_hcl, i_hocl, i_clco, i_clco3,  &
                 i_cocl2, i_s, i_so, i_so2, i_so3,       &
                 i_s2o2, i_ocs, i_hso3, i_h2so4, i_s2,   &
                 i_clso2, i_oscl, i_n2 
                 
INTEGER, SAVE :: i_h2oliq, i_h2so4liq

INTEGER, SAVE :: i_m0_aer, i_m3_aer,                       &
                 i_m0_mode1drop, i_m0_mode1ccn,            &
		 i_m3_mode1sa, i_m3_mode1w, i_m3_mode1ccn, &
		 i_m0_mode2drop, i_m0_mode2ccn,            &
		 i_m3_mode2sa, i_m3_mode2w, i_m3_mode2ccn

INTEGER, SAVE :: nmicro  ! number of microphysical tracers 

REAL, DIMENSION(:), SAVE, ALLOCATABLE :: M_tr
 

!----------------------------------------------------------------------------
! DEF FOR CL_SCHEME = 1 (AURELIEN)

!     number of clouds mode modelized
      INTEGER, PARAMETER :: nbr_mode = 3
      INTEGER :: i_cloud
      INTEGER, SAVE :: cloudmax
      INTEGER, SAVE :: cloudmin
      REAL, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: R_MEDIAN
      REAL, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: STDDEV
      
!     K_MASS coefficient correspondant à la partie condensee de chaque mode
      REAL, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: K_MASS

      REAL, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: NBRTOT
      REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: WH2SO4
      REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: rho_droplet
!----------------------------------------------------------------------------
! DEF FOR CL_SCHEME = 2 (FULL MICROPHYS)

!----------------------------------------------------------------------------

!---------------------------------------------------------------------------- 
  CONTAINS
!----------------------------------------------------------------------------

  SUBROUTINE cloud_ini(nbr_lon,nbr_lev)

!=============================================================
!	cloud_ini definit le champ 3D des caracteristiques du nuage
  
      INTEGER :: nbr_lon,nbr_lev,i_lev
    
      ALLOCATE(NBRTOT(nbr_lon,nbr_lev,nbr_mode))
      ALLOCATE(R_MEDIAN(nbr_lon,nbr_lev,nbr_mode))
      ALLOCATE(K_MASS(nbr_lon,nbr_lev,nbr_mode))
      ALLOCATE(STDDEV(nbr_lon,nbr_lev,nbr_mode))
      ALLOCATE(WH2SO4(nbr_lon,nbr_lev))
      ALLOCATE(rho_droplet(nbr_lon,nbr_lev))
            
      PRINT*,'=========================='
      PRINT*,'Initialisation cloud layer'
      PRINT*,'=========================='
      PRINT*,'nbr_lon',nbr_lon
      PRINT*,'nbr_lev',nbr_lev
      PRINT*,'nbr_mode',nbr_mode
	
      NBRTOT(:,:,:)    = 0.0E+0
      WH2SO4(:,:)      = 0.0E+0
      rho_droplet(:,:) = 0.0E+0
            
!=============================================================
!		      Initialisation cloud layer 1
!=============================================================
!     cloudmin et cloudmax niveaux du GCM
      cloudmin= 18
      cloudmax= 50

!     radius R_MEDIAN en m (donc *e-6 pour microns)
	
	R_MEDIAN(:,:,:)=0.0E+0		 ! Geometric Average Radius
	STDDEV(:,:,:)=0.0E+0             ! Geometric Std Deviation
	K_MASS(:,:,:)=0.0E+0             ! Coeff multimodal

!	===============================================
!	Knollenberg & Hunten, 1980 and James et al 1997
!	===============================================

!	===============================================
!	Initialisation UNIMODALE
!	===============================================

!     Lower Haze: mode 1
!      DO i_lev=cloudmin,20
!      R_MEDIAN(:,i_lev,1)=0.2e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Lower Cloud: mode 3
!      DO i_lev=21,23
!      R_MEDIAN(:,i_lev,1)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Middle Cloud: mode 2 prime
!      DO i_lev=24,28
!      R_MEDIAN(:,i_lev,1)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.23
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Cloud: mode 2
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,1)=1.0e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Haze: mode 1
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,1)=0.2e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=2.16
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!	===============================================
!	Initialisation TRIMODALE
!	===============================================

!     Lower Haze: mode 1
!      DO i_lev=cloudmin,20
!      R_MEDIAN(:,i_lev,1)=0.3e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
                              
!     Lower Haze: mode 2
!      DO i_lev=cloudmin,20
!      R_MEDIAN(:,i_lev,2)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.23
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO

!     Lower Haze: mode 3
!      DO i_lev=cloudmin,20
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO  

!     Lower Cloud: mode 1
!      DO i_lev=21,23
!      R_MEDIAN(:,i_lev,1)=0.3e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.1
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
                              
!     Lower Cloud: mode 2 prime
!      DO i_lev=21,23
!      R_MEDIAN(:,i_lev,2)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.23
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.4
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO

!     Lower Cloud: mode 3
!      DO i_lev=21,23
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.5
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO 

!     Middle Cloud: mode 1
!      DO i_lev=24,28
!      R_MEDIAN(:,i_lev,1)=0.3e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
         
!     Middle Cloud: mode 2 prime
!      DO i_lev=24,28
!      R_MEDIAN(:,i_lev,2)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.23
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.8
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
   
!     Middle Cloud: mode 3
!      DO i_lev=24,28
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.2
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO


!     Upper Cloud: mode 1
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,1)=0.3e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.15
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Cloud: mode 2
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,2)=1.0e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.85
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
        
!     Upper Cloud: mode 3
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO

!     Upper Haze: mode 1
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,1)=0.3e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.56
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Haze: mode 2
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,2)=1.e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
      
!     Upper Haze: mode 3
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=2.16
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO       
!=============================================================

!	===============================================
!	Initialisation TRIMODALE Knollenberg 
!	===============================================

!     Lower Haze: mode 1
      DO i_lev=cloudmin,22
      R_MEDIAN(:,i_lev,1)=0.1e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=1.57
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=1.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO
                              
!     Lower Haze: mode 2
      DO i_lev=cloudmin,22
      R_MEDIAN(:,i_lev,2)=1.4e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.23
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO

!     Lower Haze: mode 3
      DO i_lev=cloudmin,22
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=1.28
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO  

!     Pre Cloud: mode 1
      DO i_lev=23,23
      R_MEDIAN(:,i_lev,1)=0.15e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=1.8
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=0.04
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO
                              
!     Pre Cloud: mode 2
      DO i_lev=23,23
      R_MEDIAN(:,i_lev,2)=1.0e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.29
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.96
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO

!     Pre Cloud: mode 3
      DO i_lev=23,23
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=1.28
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO 
                  
!      Lower Cloud: mode 1
      DO i_lev=24,24
      R_MEDIAN(:,i_lev,1)=0.2e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=1.8
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=0.014
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO
         	
!     Lower Cloud: mode 2
      DO i_lev=24,24
      R_MEDIAN(:,i_lev,2)=1.0e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.29
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.02
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO

!     Lower Cloud: mode 3
      DO i_lev=24,24
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=1.28
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.966
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO

!     Middle Cloud: mode 1
      DO i_lev=25,28
      R_MEDIAN(:,i_lev,1)=0.15e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=1.9
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=0.0084
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO
         
!     Middle Cloud: mode 2 prime
      DO i_lev=25,28
      R_MEDIAN(:,i_lev,2)=1.4e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.23
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.21
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO
   
!     Middle Cloud: mode 3
      DO i_lev=25,28
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=1.28
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.7816
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO

!      option: upper haze remplacee par extension upper cloud
!         => 35 remplace par cloudmax et upper haze commentee
!	===============================================

!     Upper Cloud: mode 1
      DO i_lev=29,35 !cloudmax
      R_MEDIAN(:,i_lev,1)=0.2e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=2.16
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=0.72
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO

!     Upper Cloud: mode 2
      DO i_lev=29,35 !cloudmax
      R_MEDIAN(:,i_lev,2)=1.0e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.29
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.28
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO
        
!     Upper Cloud: mode 3
      DO i_lev=29,35 !cloudmax
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=1.28
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO

!     Upper Haze: mode 1
      DO i_lev=36, cloudmax
      R_MEDIAN(:,i_lev,1)=0.2e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
      STDDEV(:,i_lev,1)=2.16
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
      K_MASS(:,i_lev,1)=1.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
      END DO

!     Upper Haze: mode 2
      DO i_lev=36, cloudmax
      R_MEDIAN(:,i_lev,2)=1.e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
      STDDEV(:,i_lev,2)=1.29
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
      K_MASS(:,i_lev,2)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
      END DO
      
!     Upper Haze: mode 3
      DO i_lev=36, cloudmax
      R_MEDIAN(:,i_lev,3)=3.65e-6
      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
      STDDEV(:,i_lev,3)=2.16
      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
      K_MASS(:,i_lev,3)=0.0
      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
      END DO       

!=============================================================

!	===============================================================
!	Initialisation TRIMODALE "Knollenberg" sans Mode3, Mode2 etendu
!	===============================================================

!     Lower Haze: mode 1
!      DO i_lev=cloudmin,22
!      R_MEDIAN(:,i_lev,1)=0.1e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.57
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
                              
!     Lower Haze: mode 2
!      DO i_lev=cloudmin,22
!      R_MEDIAN(:,i_lev,2)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.23
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO

!     Lower Haze: mode 3
!      DO i_lev=cloudmin,22
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO  

!     Pre Cloud: mode 1
!      DO i_lev=23,23
!      R_MEDIAN(:,i_lev,1)=0.15e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.8
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.04
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
                              
!     Pre Cloud: mode 2
!      DO i_lev=23,23
!      R_MEDIAN(:,i_lev,2)=1.0e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.96
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO

!     Pre Cloud: mode 3
!      DO i_lev=23,23
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO 
                  
!      Lower Cloud: mode 1
!      DO i_lev=24,24
!      R_MEDIAN(:,i_lev,1)=0.2e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.8
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.014
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
         	
!     Lower Cloud: mode 2
!      DO i_lev=24,24
!      R_MEDIAN(:,i_lev,2)=1.0e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.6
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.986
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO

!     Lower Cloud: mode 3
!      DO i_lev=24,24
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO

!     Middle Cloud: mode 1
!      DO i_lev=25,28
!      R_MEDIAN(:,i_lev,1)=0.15e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=1.9
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.0084
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO
         
!     Middle Cloud: mode 2 prime
!      DO i_lev=25,28
!      R_MEDIAN(:,i_lev,2)=1.4e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.6
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.9916 
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
   
!     Middle Cloud: mode 3
!      DO i_lev=25,28
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO


!     Upper Cloud: mode 1
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,1)=0.2e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=2.16
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=0.72
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Cloud: mode 2
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,2)=1.0e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.28
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
        
!     Upper Cloud: mode 3
!      DO i_lev=29,35
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=1.28
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO

!     Upper Haze: mode 1
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,1)=0.2e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,1)
!      STDDEV(:,i_lev,1)=2.16
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,1)
!      K_MASS(:,i_lev,1)=1.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,1)
!      END DO

!     Upper Haze: mode 2
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,2)=1.e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,2)
!      STDDEV(:,i_lev,2)=1.29
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,2)
!      K_MASS(:,i_lev,2)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,2)
!      END DO
      
!     Upper Haze: mode 3
!      DO i_lev=36, cloudmax
!      R_MEDIAN(:,i_lev,3)=3.65e-6
!      PRINT*,'level',i_lev,'R R_MEDIAN',R_MEDIAN(1,i_lev,3)
!      STDDEV(:,i_lev,3)=2.16
!      PRINT*,'level',i_lev,'Dev Std',STDDEV(1,i_lev,3)
!      K_MASS(:,i_lev,3)=0.0
!      PRINT*,'level',i_lev,'Coeff Mass: k_mass',K_MASS(1,i_lev,3)
!      END DO       

!=============================================================
      PRINT*,'==============================='
      PRINT*,'FIN Initialisation cloud layer'
      PRINT*,'==============================='
      
  END SUBROUTINE cloud_ini
  
! ===========================================================

  SUBROUTINE chemparam_ini
  USE infotrac_phy, ONLY: nqtot, tname
  IMPLICIT NONE
  INTEGER :: i
  
  ALLOCATE(M_tr(nqtot))
  
  DO i=1, nqtot                                 
	
	PRINT*,'i',i
	PRINT*,'tname(i)',tname(i)
	
  	SELECT CASE(tname(i))
  		CASE('co2')
  		i_co2=i
  		PRINT*,'co2',i_co2
  		M_tr(i_co2) = 44.0095  
  		CASE('co')
  		i_co=i
  		PRINT*,'co',i_co
  		M_tr(i_co)=28.0101 
  		CASE('h2')
  		i_h2=i
  		PRINT*,'h2',i_h2
  		M_tr(i_h2)= 2.01588
  		CASE('h2o')
  		i_h2o=i
  		PRINT*,'h2o',i_h2o
  		M_tr(i_h2o)=18.0153
  		CASE('o1d')
  		i_o1d=i
  		PRINT*,'o1d',i_o1d
  		M_tr(i_o1d)=15.994 
  		CASE('o')
  		i_o=i
  		PRINT*,'o',i_o
  		M_tr(i_o)=15.994 
  		CASE('o2')
  		i_o2=i
  		PRINT*,'o2',i_o2
  		M_tr(i_o2)=31.9988
  		CASE('o2dg')
  		i_o2dg=i
  		PRINT*,'o2dg',i_o2dg
  		M_tr(i_o2dg)=31.9988
  		CASE('o3')
  		i_o3=i
  		PRINT*,'o3',i_o3
  		M_tr(i_o3)= 47.9982 
  		CASE('h')
  		i_h=i
  		PRINT*,'h',i_h
  		M_tr(i_h)= 1.00794 
  		CASE('oh')
  		i_oh=i
  		PRINT*,'oh',i_oh
  		M_tr(i_oh)=17.0073
  		CASE('ho2')
  		i_ho2=i
  		PRINT*,'ho2',i_ho2
  		M_tr(i_ho2)=33.0067
  		CASE('h2o2')
  		i_h2o2=i
  		PRINT*,'h2o2',i_h2o2
  		M_tr(i_h2o2)=34.0147 
  		CASE('cl')
  		i_cl=i
  		PRINT*,'cl',i_cl
  		M_tr(i_cl)=35.453
  		CASE('clo')
  		i_clo=i
  		PRINT*,'clo',i_clo
  		M_tr(i_clo)=51.452
  		CASE('cl2')
  		i_cl2=i
  		PRINT*,'cl2',i_cl2
  		M_tr(i_cl2)=70.906
  		CASE('hcl')
  		i_hcl=i
  		PRINT*,'hcl',i_hcl
  		M_tr(i_hcl)=36.461
  		CASE('hocl')
  		i_hocl=i
  		PRINT*,'hocl',i_hocl
  		M_tr(i_hocl)=52.46
  		CASE('clco')
  		i_clco=i
  		PRINT*,'clco',i_clco
  		M_tr(i_clco)=63.463
  		CASE('clco3')
  		i_clco3=i
  		PRINT*,'clco3',i_clco3
  		M_tr(i_clco3)=95.462 
  		CASE('cocl2')
  		i_cocl2=i
  		PRINT*,'cocl2',i_cocl2
  		M_tr(i_cocl2)=98.916
  		CASE('s')
  		i_s=i
  		PRINT*,'s',i_s
  		M_tr(i_s)=32.065
  		CASE('so')
  		i_so=i
  		PRINT*,'so',i_so
  		M_tr(i_so)=48.0644
  		CASE('so2')
  		i_so2=i
  		PRINT*,'so2',i_so2
  		M_tr(i_so2)=64.064 
  		CASE('so3')
  		i_so3=i
  		PRINT*,'so3',i_so3
  		M_tr(i_so3)=80.063
  		CASE('s2o2')
  		i_s2o2=i
  		PRINT*,'s2o2',i_s2o2
  		M_tr(i_s2o2)= 96.1288 
  		CASE('ocs')
  		i_ocs=i
  		PRINT*,'ocs',i_ocs
  		M_tr(i_ocs)=60.0751
  		CASE('hso3')
  		i_hso3=i
  		PRINT*,'hso3',i_hso3
  		M_tr(i_hso3)=81.071
  		CASE('h2so4')
  		i_h2so4=i
  		PRINT*,'h2so4',i_h2so4
  		M_tr(i_h2so4)=98.078
  		CASE('s2')
  		i_s2=i
  		PRINT*,'s2',i_s2
  		M_tr(i_s2)=64.13
  		CASE('clso2')
  		i_clso2=i
  		PRINT*,'clso2',i_clso2
  		M_tr(i_clso2)=99.517
  		CASE('oscl')
  		i_oscl=i
  		PRINT*,'oscl',i_oscl
  		M_tr(i_oscl)=83.517
		CASE('n2')
		i_n2=i
		M_tr(i_n2)=28.013
! MICROPHYSICAL TRACERS FOR CL_SCHEME=1
  		CASE('h2oliq')
  		i_h2oliq=i
  		PRINT*,'h2oliq',i_h2oliq
  		M_tr(i_h2oliq)=18.0153
  		CASE('h2so4liq')
  		i_h2so4liq=i
  		PRINT*,'h2so4liq',i_h2so4liq
  		M_tr(i_h2so4liq)=98.078
! MICROPHYSICAL TRACERS FOR CL_SCHEME=2
                CASE('M0_aer')
		i_m0_aer=i
		PRINT*,'M0_aer',i_m0_aer
		CASE('M3_aer')
		i_m3_aer=i
		PRINT*,'M3_aer',i_m3_aer
                CASE('M0_m1drop')
		i_m0_mode1drop=i
		PRINT*,'M0_m1drop',i_m0_mode1drop
                CASE('M0_m1ccn')
		i_m0_mode1ccn=i
		PRINT*,'M0_m1ccn',i_m0_mode1ccn
		CASE('M3_m1sa')
		i_m3_mode1sa=i
		PRINT*,'M3_m1sa',i_m3_mode1sa
		CASE('M3_m1w')
		i_m3_mode1w=i
		PRINT*,'M3_m1w',i_m3_mode1w
		CASE('M3_m1ccn')
		i_m3_mode1ccn=i
		PRINT*,'M3_m1ccn',i_m3_mode1ccn
                CASE('M0_m2drop')
		i_m0_mode2drop=i
		PRINT*,'M0_m2drop',i_m0_mode2drop
                CASE('M0_m2ccn')
		i_m0_mode2ccn=i
		PRINT*,'M0_m2ccn',i_m0_mode2ccn
		CASE('M3_m2sa')
		i_m3_mode2sa=i
		PRINT*,'M3_m2sa',i_m3_mode2sa
		CASE('M3_m2w')
		i_m3_mode2w=i
		PRINT*,'M3_m2w',i_m3_mode2w
		CASE('M3_m2ccn')
		i_m3_mode2ccn=i
		PRINT*,'M3_m2ccn',i_m3_mode2ccn
  	END SELECT
  	
!  	PRINT*,'M_tr(i)',M_tr(i)
  END DO
  
  END SUBROUTINE chemparam_ini

! ===========================================================

  SUBROUTINE vapors4muphy_ini(nlon,nlev,trac)
  USE infotrac_phy, ONLY: nqtot, tname
  IMPLICIT NONE

  integer :: nlon, nlev
  real    :: trac(nlon,nlev,nqtot) ! traceur ( en vmr)

!  integer :: i
!  real    :: trac1d(nlev,2) ! traceur lu ( en vmr)
 
! lecture d'un fichier texte contenant les profils de trac1d(:1) = H2O et trac1d(:,2) = H2SO4
!  DO i=1,nlon
!     trac(i,:,i_h2o) = trac1d(:,1)
!     trac(i,:,i_h2so4) = trac1d(:,2)
!  ENDDO

!  intitialisation profils altitude H2O et H2SO4
!  profil H2O initial vap+liq == que vap
   trac(:,1:24,i_h2o) = 30.E-6 !
   trac(:,25:50,i_h2o) = 1.E-6 !

   trac(:,:,i_h2so4) = 3.E-9 ! Limite sup Sandor 2012
   trac(:,23:50,i_h2so4) = 2.E-6 ! Profil H2SO4 initial => vap+liq

  END SUBROUTINE vapors4muphy_ini

END MODULE chemparam_mod

