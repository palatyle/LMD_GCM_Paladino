module radcommon_h
      use radinc_h, only: L_NSPECTI, L_NSPECTV, NTstart, NTstop, &
                          naerkind, nsizemax
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
!     TAURAYVAR  - Array (NSPECTV elements) of the pressure-independent
!                  part of Rayleigh scattering optical depth for the variable gas.
!     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!     FZEROV     - Fraction of zeros in the VISUAL CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
! 
!   Shortwave
!   ~~~~~~~~~
! 
! For the "naerkind" kind of aerosol radiative properties : 
! QVISsQREF  :  Qext / Qext("longrefvis")
! omegavis   :  single scattering albedo
! gvis       :  assymetry factor
! 
!   Longwave
!   ~~~~~~~~
! 
! For the "naerkind" kind of aerosol radiative properties : 
! QIRsQREF :  Qext / Qext("longrefvis")
! omegaIR  :  mean single scattering albedo
! gIR      :  mean assymetry factor

      REAL*8 BWNI(L_NSPECTI+1), WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI) !BWNI read by master in setspi
      REAL*8 BWNV(L_NSPECTV+1), WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV) !BWNV read by master in setspv
      REAL*8 STELLARF(L_NSPECTV), TAURAY(L_NSPECTV), TAURAYVAR(L_NSPECTV)
!$OMP THREADPRIVATE(WNOI,DWNI,WAVEI,&
	!$OMP WNOV,DWNV,WAVEV,&
	!$OMP STELLARF,TAURAY,TAURAYVAR)

      REAL*8 blami(L_NSPECTI+1)
      REAL*8 blamv(L_NSPECTV+1) ! these are needed by suaer.F90
!$OMP THREADPRIVATE(blami,blamv)

      !! AS: introduced to avoid doing same computations again for continuum
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: indi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: indv
!$OMP THREADPRIVATE(indi,indv)

      !!! ALLOCATABLE STUFF SO THAT DIMENSIONS ARE READ in *.dat FILES -- AS 12/2011  
      REAL*8, DIMENSION(:,:,:,:,:), ALLOCATABLE :: gasi, gasv
      REAL*8, DIMENSION(:), ALLOCATABLE :: PGASREF, TGASREF, WREFVAR, PFGASREF, GWEIGHT
      real*8 FZEROI(L_NSPECTI)
      real*8 FZEROV(L_NSPECTV)
      real*8 pgasmin, pgasmax
      real*8 tgasmin, tgasmax
!$OMP THREADPRIVATE(gasi,gasv,&  !wrefvar,pgasref,tgasref,pfgasref read by master in sugas_corrk
	!$OMP FZEROI,FZEROV)     !pgasmin,pgasmax,tgasmin,tgasmax read by master in sugas_corrk

      real QVISsQREF(L_NSPECTV,naerkind,nsizemax)
      real omegavis(L_NSPECTV,naerkind,nsizemax)
      real gvis(L_NSPECTV,naerkind,nsizemax)
      real QIRsQREF(L_NSPECTI,naerkind,nsizemax)
      real omegair(L_NSPECTI,naerkind,nsizemax)
      real gir(L_NSPECTI,naerkind,nsizemax)
!$OMP THREADPRIVATE(QVISsQREF,omegavis,gvis,QIRsQREF,omegair,gir)


! Reference wavelengths used to compute reference optical depth (m)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL lamrefir(naerkind),lamrefvis(naerkind)

! Actual number of grain size classes in each domain for a
!   given aerosol:

      INTEGER          :: nsize(naerkind,2)

! Particle size axis (depend on the kind of aerosol and the
!   radiation domain)

      DOUBLE PRECISION :: radiustab(naerkind,2,nsizemax)
!$OMP THREADPRIVATE(lamrefir,lamrefvis,radiustab) !nsize read by suaer_corrk

! Extinction coefficient at reference wavelengths;
!   These wavelengths are defined in aeroptproperties, and called
!   longrefvis and longrefir.

      REAL :: QREFvis(naerkind,nsizemax)
      REAL :: QREFir(naerkind,nsizemax)
!      REAL :: omegaREFvis(naerkind,nsizemax)
      REAL :: omegaREFir(naerkind,nsizemax)

      REAL,SAVE :: tstellar ! Stellar brightness temperature (SW)

      REAL*8, DIMENSION(:,:), ALLOCATABLE, SAVE :: planckir

      real*8,save :: PTOP

      real*8,parameter :: UBARI = 0.5D0

!$OMP THREADPRIVATE(QREFvis,QREFir,omegaREFir,& 	! gweight read by master in sugas_corrk
		!$OMP tstellar,planckir,PTOP)

!     If the gas optical depth (top to the surface) is less than
!     this value, we place that Gauss-point into the "zeros"
!     channel.
      real*8, parameter :: TLIMIT =  1.0D-30

!     Factor to convert pressures from millibars to Pascals
      real*8, parameter :: SCALEP = 1.00D+2

      real*8, parameter :: sigma = 5.67032D-8
      real*8, parameter :: grav = 6.672E-11

      real*8,save :: Cmk
      real*8,save :: glat_ig
!$OMP THREADPRIVATE(Cmk,glat_ig)

      ! extinction of incoming sunlight (Saturn's rings, eclipses, etc...)
      REAL, DIMENSION(:), ALLOCATABLE ,SAVE :: eclipse

      !Latitude-dependent gravity
      REAL, DIMENSION(:), ALLOCATABLE , SAVE :: glat
!$OMP THREADPRIVATE(glat,eclipse)

end module radcommon_h
