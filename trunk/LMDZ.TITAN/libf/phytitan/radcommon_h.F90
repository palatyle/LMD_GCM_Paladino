module radcommon_h
      use radinc_h, only: L_NSPECTI, L_NSPECTV, NTstar, NTstop
      implicit none

!----------------------------------------------------------------------C
!
!                             radcommon.h
!
!----------------------------------------------------------------------C
!
!  "Include" grid.h and radinc.h before this file in code that uses
!  some or all of this common data set
!
!     WNOI       - Array of wavenumbers at the spectral interval
!                  centers for the infrared.  Array is NSPECTI
!                  elements long.
!     DWNI       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each IR spectral
!                  interval.  NSPECTI elements long.
!     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!                  (in microns) at the center of each IR spectral
!                  interval.
!     WNOV       - Array of wavenumbers at the spectral interval
!                  center for the VISUAL.  Array is NSPECTV
!                  elements long.
!     DWNV       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each VISUAL spectral
!                  interval.  NSPECTV elements long.
!     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!                  (in microns) at the center of each VISUAL spectral
!                  interval.
!     STELLARF   - Array (NSPECTV elements) of stellar flux (W/M^2) in
!                  each spectral interval.  Values are for 1 AU,
!                  scaled to the planetary distance elsewhere.
!     TAURAY     - Array (NSPECTV elements) of the pressure-independent
!                  part of Rayleigh scattering optical depth.
!     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!     FZEROV     - Fraction of zeros in the VISUAL CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!

      REAL*8 BWNI(L_NSPECTI+1), WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI) !BWNI read by master in setspi
      REAL*8 BWNV(L_NSPECTV+1), WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV) !BWNV read by master in setspv
      REAL*8 STELLARF(L_NSPECTV), TAURAY(L_NSPECTV)
!$OMP THREADPRIVATE(WNOI,DWNI,WAVEI,&
	!$OMP WNOV,DWNV,WAVEV,&
	!$OMP STELLARF,TAURAY)

      !! AS: introduced to avoid doing same computations again for continuum
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: indi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: indv
!$OMP THREADPRIVATE(indi,indv)

      !!! ALLOCATABLE STUFF SO THAT DIMENSIONS ARE READ in *.dat FILES -- AS 12/2011  
      REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: gasi, gasv
      REAL*8, DIMENSION(:), ALLOCATABLE :: PGASREF, TGASREF, PFGASREF, GWEIGHT
 
      ! For corrk recombining - JVO 18
      REAL*8,  DIMENSION(:,:,:,:), ALLOCATABLE :: gasi_recomb, gasv_recomb
      REAL*8,  DIMENSION(:,:),     ALLOCATABLE :: pqrold
      REAL*8,  DIMENSION(:),       ALLOCATABLE :: w_cum
      LOGICAL, DIMENSION(:,:),     ALLOCATABLE :: useptold
      INTEGER, DIMENSION(:),       ALLOCATABLE :: permut_idx
!$OMP THREADPRIVATE(gasi_recomb,gasv_recomb,pqrold,w_cum,useptold,permut_idx)
      
      LOGICAL, DIMENSION(:), ALLOCATABLE :: RADVAR_MASK ! for corrk recombin -> read by master in sugas_corrk
      INTEGER, DIMENSION(:), ALLOCATABLE :: RADVAR_INDX ! for corrk recombin -> read by master in sugas_corrk
      
      real*8 FZEROI(L_NSPECTI)
      real*8 FZEROV(L_NSPECTV)
      real*8 pgasmin, pgasmax
      real*8 tgasmin, tgasmax
!$OMP THREADPRIVATE(gasi,gasv,&  !pgasref,tgasref,pfgasref, gweight read by master in sugas_corrk
	!$OMP FZEROI,FZEROV)     !pgasmin,pgasmax,tgasmin,tgasmax read by master in sugas_corrk


      REAL,SAVE :: tstellar ! Stellar brightness temperature (SW)

      real*8,save :: planckir(L_NSPECTI,NTstop-NTstar+1)

      real*8,save :: PTOP
!$OMP THREADPRIVATE(tstellar,planckir,PTOP)

      real*8,parameter :: UBARI = 0.5D0

!     If the gas optical depth (top to the surface) is less than
!     this value, we place that Gauss-point into the "zeros"
!     channel.
      real*8, parameter :: TLIMIT =  1.0D-30

!     Factor to convert pressures from millibars to Pascals
      real*8, parameter :: SCALEP = 1.00D+2

      real*8, parameter :: sigma = 5.67032D-8
      real*8, parameter :: grav = 6.672E-11

      ! extinction of incoming sunlight (Saturn's rings, eclipses, etc...)
      REAL, DIMENSION(:), ALLOCATABLE ,SAVE :: eclipse
!$OMP THREADPRIVATE(eclipse)

      ! Altitude-Latitude-dependent gravity
      REAL, DIMENSION(:,:), ALLOCATABLE , SAVE :: gzlat ! This should be stored elsewhere ...
      REAL, DIMENSION(:), ALLOCATABLE , SAVE :: gzlat_ig
      REAL, DIMENSION(:), ALLOCATABLE , SAVE :: Cmk 
!$OMP THREADPRIVATE(gzlat,gzlat_ig,Cmk)

end module radcommon_h
