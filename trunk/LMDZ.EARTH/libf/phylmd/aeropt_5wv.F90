!
! $Id: aeropt_5wv.F90 1470 2010-12-21 15:51:32Z idelkadi $
!

SUBROUTINE AEROPT_5WV(&
   pdel, m_allaer, delt, &
   RHcl, ai, flag_aerosol, &
   pplay, t_seri, &
   tausum, tau, presnivs)

  USE DIMPHY
  USE aero_mod
  USE phys_local_var_mod, only: od550aer,od865aer,ec550aer,od550lt1aer

  !
  !    Yves Balkanski le 12 avril 2006
  !    Celine Deandreis
  !    Anne Cozic  Avril 2009
  !    a partir d'une sous-routine de Johannes Quaas pour les sulfates
  !
  !
  ! Refractive indices for seasalt come from Shettle and Fenn (1979)
  !
  ! Refractive indices from water come from Hale and Querry (1973)
  !
  ! Refractive indices from Ammonium Sulfate Toon and Pollack (1976)
  !
  ! Refractive indices for Dust, internal mixture of minerals coated with 1.5% hematite 
  ! by Volume (Balkanski et al., 2006)
  !
  ! Refractive indices for POM: Kinne (pers. Communication 
  !
  ! Refractive index for BC from Shettle and Fenn (1979)
  !
  ! Shettle, E. P., & Fenn, R. W. (1979), Models for the aerosols of the lower atmosphere and 
  ! the effects of humidity variations on their optical properties, U.S. Air Force Geophysics 
  ! Laboratory Rept. AFGL-TR-79-0214, Hanscomb Air Force Base, MA.
  !
  ! Hale, G. M. and M. R. Querry, Optical constants of water in the 200-nm to 200-m 
  ! wavelength region, Appl. Opt., 12, 555-563, 1973.
  !
  ! Toon, O. B. and J. B. Pollack, The optical constants of several atmospheric aerosol species:
  ! Ammonium sulfate, aluminum oxide, and sodium chloride, J. Geohys. Res., 81, 5733-5748,
  ! 1976.
  !
  ! Balkanski, Y., M. Schulz, T. Claquin And O. Boucher, Reevaluation of mineral aerosol 
  ! radiative forcings suggests a better agreement with satellite and AERONET data, Atmospheric 
  ! Chemistry and Physics Discussions., 6, pp 8383-8419, 2006.
  !
  IMPLICIT NONE
  INCLUDE "YOMCST.h"
  !
  ! Input arguments:
  !
  REAL, DIMENSION(klon,klev), INTENT(in)   :: pdel
  REAL, INTENT(in)                         :: delt
  REAL, DIMENSION(klon,klev,naero_spc), INTENT(in) :: m_allaer
  REAL, DIMENSION(klon,klev), INTENT(in)   :: RHcl     ! humidite relative ciel clair
  INTEGER,INTENT(in)                       :: flag_aerosol
  REAL, DIMENSION(klon,klev), INTENT(in)   :: pplay
  REAL, DIMENSION(klon,klev), INTENT(in)   :: t_seri
  REAL, DIMENSION(klev),      INTENT(in)   :: presnivs
  !
  ! Output arguments:
  !
  REAL, DIMENSION(klon), INTENT(out)          :: ai      ! POLDER aerosol index 
  REAL, DIMENSION(klon,nwave,naero_spc), INTENT(out)      :: tausum
  REAL, DIMENSION(klon,klev,nwave,naero_spc), INTENT(out) :: tau


  !
  ! Local
  !
  INTEGER, PARAMETER :: las = nwave
  LOGICAL :: soluble
  
  INTEGER :: i, k, ierr, m
  INTEGER :: spsol, spinsol, spss, la
  INTEGER :: RH_num(klon,klev)
  INTEGER, PARAMETER :: la443 = 1
  INTEGER, PARAMETER :: la550 = 2
  INTEGER, PARAMETER :: la670 = 3
  INTEGER, PARAMETER :: la765 = 4
  INTEGER, PARAMETER :: la865 = 5
  INTEGER, PARAMETER :: nbre_RH=12
  INTEGER, PARAMETER :: naero_soluble=7   !  1- BC soluble; 2- POM soluble; 3- SO4 acc.
                                          !  4- SO4 coarse; 5 seasalt super-C; 6 seasalt coarse; 7 seasalt acc.
  INTEGER, PARAMETER :: naero_insoluble=3 !  1- Dust; 2- BC insoluble; 3- POM insoluble
  INTEGER, PARAMETER :: nb_level = 19     ! number of vertical levels
  LOGICAL, SAVE :: firstcall=.TRUE.
!$OMP THREADPRIVATE(firstcall)

  REAL :: zrho

  ! Coefficient optiques sur 19 niveaux
  REAL, SAVE, DIMENSION(nb_level) :: presnivs_19  ! Pression milieux couche pour 19 niveaux (nb_level)
!$OMP THREADPRIVATE(presnivs_19)

  REAL, SAVE, DIMENSION(nb_level) :: A1_ASSSM_19, A2_ASSSM_19, A3_ASSSM_19,&
          B1_ASSSM_19, B2_ASSSM_19, C1_ASSSM_19, C2_ASSSM_19,&
          A1_CSSSM_19, A2_CSSSM_19, A3_CSSSM_19,&
          B1_CSSSM_19, B2_CSSSM_19, C1_CSSSM_19, C2_CSSSM_19, &
          A1_SSSSM_19, A2_SSSSM_19, A3_SSSSM_19,&
          B1_SSSSM_19, B2_SSSSM_19, C1_SSSSM_19, C2_SSSSM_19
!$OMP THREADPRIVATE(A1_ASSSM_19, A2_ASSSM_19, A3_ASSSM_19)
!$OMP THREADPRIVATE(B1_ASSSM_19, B2_ASSSM_19, C1_ASSSM_19, C2_ASSSM_19)
!$OMP THREADPRIVATE(A1_CSSSM_19, A2_CSSSM_19, A3_CSSSM_19)
!$OMP THREADPRIVATE(B1_CSSSM_19, B2_CSSSM_19, C1_CSSSM_19, C2_CSSSM_19)
!$OMP THREADPRIVATE(A1_SSSSM_19, A2_SSSSM_19, A3_SSSSM_19)
!$OMP THREADPRIVATE(B1_SSSSM_19, B2_SSSSM_19, C1_SSSSM_19, C2_SSSSM_19)

  ! Coefficient optiques interpole sur le nombre de niveau du modele
  REAL, ALLOCATABLE,  DIMENSION(:), SAVE :: &
          A1_ASSSM, A2_ASSSM, A3_ASSSM,&
          B1_ASSSM, B2_ASSSM, C1_ASSSM, C2_ASSSM,&
          A1_CSSSM, A2_CSSSM, A3_CSSSM,&
          B1_CSSSM, B2_CSSSM, C1_CSSSM, C2_CSSSM, &
          A1_SSSSM, A2_SSSSM, A3_SSSSM,&
          B1_SSSSM, B2_SSSSM, C1_SSSSM, C2_SSSSM
!$OMP THREADPRIVATE(A1_ASSSM, A2_ASSSM, A3_ASSSM)
!$OMP THREADPRIVATE(B1_ASSSM, B2_ASSSM, C1_ASSSM, C2_ASSSM)
!$OMP THREADPRIVATE(A1_CSSSM, A2_CSSSM, A3_CSSSM)
!$OMP THREADPRIVATE(B1_CSSSM, B2_CSSSM, C1_CSSSM, C2_CSSSM)
!$OMP THREADPRIVATE(A1_SSSSM, A2_SSSSM, A3_SSSSM)
!$OMP THREADPRIVATE(B1_SSSSM, B2_SSSSM, C1_SSSSM, C2_SSSSM)


  REAL,PARAMETER :: RH_tab(nbre_RH)=(/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./)
  REAL :: DELTA(klon,klev), rh(klon,klev), H
  REAL :: tau_ae5wv_int ! Intermediate computation of epaisseur optique aerosol
  REAL :: piz_ae5wv_int ! Intermediate single scattering albedo aerosol
  REAL :: cg_ae5wv_int  ! Intermediate asymmetry parameter aerosol
  REAL, PARAMETER :: RH_MAX=95.
  REAL :: taue670(KLON)       ! epaisseur optique aerosol absorption 550 nm
  REAL :: taue865(KLON)       ! epaisseur optique aerosol extinction 865 nm
  REAL :: fac
  REAL :: zdp1(klon,klev) 
  REAL, PARAMETER ::  gravit = 9.80616    ! m2/s
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: aerosol_name
  INTEGER :: nb_aer
  
  REAL :: tau3d(KLON,KLEV), piz3d(KLON,KLEV), cg3d(KLON,KLEV)
  REAL :: abs3d(KLON,KLEV)     ! epaisseur optique d'absorption
  REAL :: dh(KLON,KLEV)
  
  REAL :: alpha_aers_5wv(nbre_RH,las,naero_soluble)   ! ext. coeff. Soluble comp. units *** m2/g 
   !  1- BC soluble; 2- POM soluble; 3- SO4 acc.; 4- SO4 coarse; 5 seasalt super-C; 6 seasalt coarse; 7 seasalt acc.
  REAL :: alpha_aeri_5wv(las,naero_insoluble)         ! ext. coeff. Insoluble comp. 1- Dust: 2- BC; 3- POM
  REAL :: cg_aers_5wv(nbre_RH,las,naero_soluble)      ! Asym. param. soluble comp. 
   !  1- BC soluble; 2- POM soluble; 3- SO4 acc.; 4- SO4 coarse; 5 seasalt super-C; 6 seasalt coarse; 7 seasalt acc.
  REAL :: cg_aeri_5wv(las,naero_insoluble)            ! Asym. param. insoluble comp. 1- Dust: 2- BC; 3- POM
  REAL :: piz_aers_5wv(nbre_RH,las,naero_soluble)   
   !  1- BC soluble; 2- POM soluble; 3- SO4 acc.; 4- SO4 coarse; 5 seasalt super-C; 6 seasalt coarse; 7 seasalt acc.
  REAL :: piz_aeri_5wv(las,naero_insoluble)           ! Insoluble comp. 1- Dust: 2- BC; 3- POM

  REAL, DIMENSION(klon,klev,naero_spc) :: mass_temp
  
  !
  ! Proprietes optiques
  !
  REAL :: radry = 287.054
  REAL :: tau_tmp                     ! dry air mass constant
  REAL :: fact_RH(nbre_RH)
  LOGICAL :: used_tau(naero_spc)
  INTEGER :: n
  
  DATA presnivs_19/&
       100426.5,  98327.6, 95346.5, 90966.8, 84776.9, &
       76536.5,   66292.2, 54559.3, 42501.8, 31806, &
       23787.5,   18252.7, 13996,   10320.8, 7191.1, &
       4661.7,    2732.9,  1345.6,  388.2/

!!ACCUMULATION MODE
  DATA A1_ASSSM_19/ 4.373E+00,  4.361E+00,  4.331E+00, &
                 4.278E+00,  4.223E+00,  4.162E+00, &
                 4.103E+00,  4.035E+00,  3.962E+00, &
                 3.904E+00,  3.871E+00,  3.847E+00, &
                 3.824E+00,  3.780E+00,  3.646E+00, &
                 3.448E+00,  3.179E+00,  2.855E+00,  2.630E+00/
  DATA A2_ASSSM_19/ 2.496E+00,  2.489E+00,  2.472E+00, &
                 2.442E+00,  2.411E+00,  2.376E+00, &
                 2.342E+00,  2.303E+00,  2.261E+00, &
                 2.228E+00,  2.210E+00,  2.196E+00, &
                 2.183E+00,  2.158E+00,  2.081E+00, &
                 1.968E+00,  1.814E+00,  1.630E+00,  1.501E+00/
  DATA A3_ASSSM_19/-4.688E-02, -4.676E-02, -4.644E-02, &
                -4.587E-02, -4.528E-02, -4.463E-02, &
                -4.399E-02, -4.326E-02, -4.248E-02, &
                -4.186E-02, -4.151E-02, -4.125E-02, &
                -4.100E-02, -4.053E-02, -3.910E-02, &
                -3.697E-02, -3.408E-02, -3.061E-02, -2.819E-02/
  DATA B1_ASSSM_19/ 1.165E-08,  1.145E-08,  1.097E-08, &
                 1.012E-08,  9.233E-09,  8.261E-09, &
                 7.297E-09,  6.201E-09,  5.026E-09, &
                 4.098E-09,  3.567E-09,  3.187E-09, &
                 2.807E-09,  2.291E-09,  2.075E-09, &
                 1.756E-09,  1.322E-09,  8.011E-10, 4.379E-10/
  DATA B2_ASSSM_19/ 2.193E-08,  2.192E-08,  2.187E-08, &
                 2.179E-08,  2.171E-08,  2.162E-08, &
                 2.153E-08,  2.143E-08,  2.132E-08, &
                 2.124E-08,  2.119E-08,  2.115E-08, &
                 2.112E-08,  2.106E-08,  2.100E-08, &
                 2.090E-08,  2.077E-08,  2.061E-08,  2.049E-08/
  DATA C1_ASSSM_19/ 7.365E-01,  7.365E-01,  7.365E-01, &
                 7.364E-01,  7.363E-01,  7.362E-01, &
                 7.361E-01,  7.359E-01,  7.358E-01, &
                 7.357E-01,  7.356E-01,  7.356E-01, &
                 7.356E-01,  7.355E-01,  7.354E-01, &
                 7.352E-01,  7.350E-01,  7.347E-01,  7.345E-01/
  DATA C2_ASSSM_19/ 5.833E-02,  5.835E-02,  5.841E-02, &
                 5.850E-02,  5.859E-02,  5.870E-02, &
                 5.880E-02,  5.891E-02,  5.904E-02, &
                 5.914E-02,  5.920E-02,  5.924E-02, &
                 5.928E-02,  5.934E-02,  5.944E-02, &
                 5.959E-02,  5.979E-02,  6.003E-02,  6.020E-02/
!COARSE MODE
  DATA A1_CSSSM_19/ 7.403E-01,  7.422E-01,  7.626E-01, &
                 8.019E-01,  8.270E-01,  8.527E-01, &
                 8.702E-01,  8.806E-01,  8.937E-01, &
                 9.489E-01,  1.030E+00,  1.105E+00, &
                 1.199E+00,  1.357E+00,  1.660E+00, &
                 2.540E+00,  4.421E+00,  2.151E+00,  9.518E-01/
  DATA A2_CSSSM_19/ 4.522E-01,  4.532E-01,  4.644E-01, &
                 4.859E-01,  4.996E-01,  5.137E-01, &
                 5.233E-01,  5.290E-01,  5.361E-01, &
                 5.655E-01,  6.085E-01,  6.483E-01, &
                 6.979E-01,  7.819E-01,  9.488E-01, &
                 1.450E+00,  2.523E+00,  1.228E+00,  5.433E-01/
  DATA A3_CSSSM_19/-8.516E-03, -8.535E-03, -8.744E-03, &
                -9.148E-03, -9.406E-03, -9.668E-03, &
                -9.848E-03, -9.955E-03, -1.009E-02, &
                -1.064E-02, -1.145E-02, -1.219E-02, &
                -1.312E-02, -1.470E-02, -1.783E-02, &
                -2.724E-02, -4.740E-02, -2.306E-02, -1.021E-02/
  DATA B1_CSSSM_19/ 2.535E-07,  2.530E-07,  2.479E-07, &
                 2.380E-07,  2.317E-07,  2.252E-07, &
                 2.208E-07,  2.182E-07,  2.149E-07, &
                 2.051E-07,  1.912E-07,  1.784E-07, &
                 1.624E-07,  1.353E-07,  1.012E-07, &
                 6.016E-08,  2.102E-08,  0.000E+00,  0.000E+00/
  DATA B2_CSSSM_19/ 1.221E-07,  1.217E-07,  1.179E-07, &
                 1.104E-07,  1.056E-07,  1.008E-07, &
                 9.744E-08,  9.546E-08,  9.299E-08, &
                 8.807E-08,  8.150E-08,  7.544E-08, &
                 6.786E-08,  5.504E-08,  4.080E-08, &
                 2.960E-08,  2.300E-08,  2.030E-08,  1.997E-08/
  DATA C1_CSSSM_19/ 7.659E-01,  7.658E-01,  7.652E-01, &
                 7.639E-01,  7.631E-01,  7.623E-01, &
                 7.618E-01,  7.614E-01,  7.610E-01, &
                 7.598E-01,  7.581E-01,  7.566E-01, &
                 7.546E-01,  7.513E-01,  7.472E-01, &
                 7.423E-01,  7.376E-01,  7.342E-01,  7.334E-01/
  DATA C2_CSSSM_19/ 3.691E-02,  3.694E-02,  3.729E-02, &
                 3.796E-02,  3.839E-02,  3.883E-02, &
                 3.913E-02,  3.931E-02,  3.953E-02, &
                 4.035E-02,  4.153E-02,  4.263E-02, &
                 4.400E-02,  4.631E-02,  4.933E-02, &
                 5.331E-02,  5.734E-02,  6.053E-02,  6.128E-02/
!SUPER COARSE MODE
  DATA A1_SSSSM_19/ 2.836E-01,  2.876E-01,  2.563E-01, &
                 2.414E-01,  2.541E-01,  2.546E-01, &
                 2.572E-01,  2.638E-01,  2.781E-01, &
                 3.167E-01,  4.209E-01,  5.286E-01, &
                 6.959E-01,  9.233E-01,  1.282E+00, &
                 1.836E+00,  2.981E+00,  4.355E+00,  4.059E+00/
  DATA A2_SSSSM_19/ 1.608E-01,  1.651E-01,  1.577E-01, &
                 1.587E-01,  1.686E-01,  1.690E-01, &
                 1.711E-01,  1.762E-01,  1.874E-01, &
                 2.138E-01,  2.751E-01,  3.363E-01, &
                 4.279E-01,  5.519E-01,  7.421E-01, &
                 1.048E+00,  1.702E+00,  2.485E+00,  2.317E+00/
  DATA A3_SSSSM_19/-3.025E-03, -3.111E-03, -2.981E-03, &
                -3.005E-03, -3.193E-03, -3.200E-03, &
                -3.239E-03, -3.336E-03, -3.548E-03, &
                -4.047E-03, -5.196E-03, -6.345E-03, &
                -8.061E-03, -1.038E-02, -1.395E-02, &
                -1.970E-02, -3.197E-02, -4.669E-02, -4.352E-02/
  DATA B1_SSSSM_19/ 6.759E-07,  6.246E-07,  5.542E-07, &
                 4.953E-07,  4.746E-07,  4.738E-07, &
                 4.695E-07,  4.588E-07,  4.354E-07, &
                 3.947E-07,  3.461E-07,  3.067E-07, &
                 2.646E-07,  2.095E-07,  1.481E-07, &
                 9.024E-08,  5.747E-08,  2.384E-08,  6.599E-09/
  DATA B2_SSSSM_19/ 5.977E-07,  5.390E-07,  4.468E-07, &
                 3.696E-07,  3.443E-07,  3.433E-07, &
                 3.380E-07,  3.249E-07,  2.962E-07, &
                 2.483E-07,  1.989E-07,  1.623E-07, &
                 1.305E-07,  9.015E-08,  6.111E-08, &
                 3.761E-08,  2.903E-08,  2.337E-08,  2.147E-08/
  DATA C1_SSSSM_19/ 8.120E-01,  8.084E-01,  8.016E-01, &
                 7.953E-01,  7.929E-01,  7.928E-01, &
                 7.923E-01,  7.910E-01,  7.882E-01, &
                 7.834E-01,  7.774E-01,  7.725E-01, &
                 7.673E-01,  7.604E-01,  7.529E-01, &
                 7.458E-01,  7.419E-01,  7.379E-01,  7.360E-01/
  DATA C2_SSSSM_19/ 2.388E-02,  2.392E-02,  2.457E-02,  2.552E-02, &
                 2.615E-02,  2.618E-02,  2.631E-02,  2.663E-02, &
                 2.735E-02,  2.875E-02,  3.113E-02,  3.330E-02, &
                 3.615E-02,  3.997E-02,  4.521E-02,  5.038E-02, &
                 5.358E-02,  5.705E-02,  5.887E-02/
!*********************************************************************
!
!
! 
! 
!  
! 
! From here on we look at the optical parameters at 5 wavelengths:  
! 443nm, 550, 670, 765 and 865 nm 
!                                   le 12 AVRIL 2006 
!  
 DATA alpha_aers_5wv/ & 
                                ! bc soluble 
       7.930,7.930,7.930,7.930,7.930,7.930,     & 
       7.930,7.930,10.893,12.618,14.550,16.613, & 
       7.658,7.658,7.658,7.658,7.658,7.658,     & 
       7.658,7.658,10.351,11.879,13.642,15.510, & 
       7.195,7.195,7.195,7.195,7.195,7.195,     & 
       7.195,7.195,9.551,10.847,12.381,13.994,  & 
       6.736,6.736,6.736,6.736,6.736,6.736,     & 
       6.736,6.736,8.818,9.938,11.283,12.687,   & 
       6.277,6.277,6.277,6.277,6.277,6.277,     & 
       6.277,6.277,8.123,9.094,10.275,11.501,   & 
                                ! pom soluble 
       6.676,6.676,6.676,6.676,6.710,6.934,   & 
       7.141,7.569,8.034,8.529,9.456,10.511,  & 
       5.109,5.109,5.109,5.109,5.189,5.535,   & 
       5.960,6.852,8.008,9.712,12.897,19.676, & 
       3.718,3.718,3.718,3.718,3.779,4.042,   & 
       4.364,5.052,5.956,7.314,9.896,15.688,  & 
       2.849,2.849,2.849,2.849,2.897,3.107,   & 
       3.365,3.916,4.649,5.760,7.900,12.863,  & 
       2.229,2.229,2.229,2.229,2.268,2.437,   & 
       2.645,3.095,3.692,4.608,6.391,10.633,  & 
                                ! Sulfate (Accumulation) 
       5.751,6.215,6.690,7.024,7.599,8.195,      & 
       9.156,10.355,12.660,14.823,18.908,24.508, & 
       4.320,4.675,5.052,5.375,5.787,6.274,      & 
       7.066,8.083,10.088,12.003,15.697,21.133,  & 
       3.079,3.351,3.639,3.886,4.205,4.584,      & 
       5.206,6.019,7.648,9.234,12.391,17.220,    & 
       2.336,2.552,2.781,2.979,3.236,3.540,      & 
       4.046,4.711,6.056,7.388,10.093,14.313,    & 
       1.777,1.949,2.134,2.292,2.503,2.751,      & 
       3.166,3.712,4.828,5.949,8.264,11.922,     & 
                                ! Sulfate (Coarse) 
       5.751,6.215,6.690,7.024,7.599,8.195,      & 
       9.156,10.355,12.660,14.823,18.908,24.508, & 
       4.320,4.675,5.052,5.375,5.787,6.274,      & 
       7.066,8.083,10.088,12.003,15.697,21.133,  & 
       3.079,3.351,3.639,3.886,4.205,4.584,      & 
       5.206,6.019,7.648,9.234,12.391,17.220,    & 
       2.336,2.552,2.781,2.979,3.236,3.540,      & 
       4.046,4.711,6.056,7.388,10.093,14.313,    & 
       1.777,1.949,2.134,2.292,2.503,2.751,      & 
       3.166,3.712,4.828,5.949,8.264,11.922,     & 
                                ! Seasalt soluble super_coarse (computed below for 550nm) 
       0.50,0.90,1.05,1.21,1.40,2.41, &  
       2.66,3.11,3.88,4.52,5.69,8.84, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
     0.52,0.93,1.08,1.24,1.43,2.47, &  
     2.73,3.20,3.99,4.64,5.84,9.04, &  
     0.52,0.93,1.09,1.25,1.44,2.50, &  
     2.76,3.23,4.03,4.68,5.89,9.14, &  
     0.52,0.94,1.09,1.26,1.45,2.51, &  
     2.78,3.25,4.06,4.72,5.94,9.22, &  
                                ! seasalt soluble coarse (computed below for 550nm) 
       0.50,0.90,1.05,1.21,1.40,2.41, &  
       2.66,3.11,3.88,4.52,5.69,8.84, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
     0.52,0.93,1.08,1.24,1.43,2.47, &  
     2.73,3.20,3.99,4.64,5.84,9.04, &  
     0.52,0.93,1.09,1.25,1.44,2.50, &  
     2.76,3.23,4.03,4.68,5.89,9.14, &  
     0.52,0.94,1.09,1.26,1.45,2.51, &  
     2.78,3.25,4.06,4.72,5.94,9.22, &  
                                ! seasalt soluble accumulation (computed below for 550nm) 
     4.28, 7.17, 8.44, 9.85,11.60,22.44,  &  
     25.34,30.54,39.38,46.52,59.33,91.77, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
       0.000,0.000,0.000,0.000,0.000,0.000, &  
     2.48, 4.22, 5.02, 5.94, 7.11,15.29,  &  
     17.70,22.31,30.73,38.06,52.15,90.59, &  
     1.90, 3.29, 3.94, 4.69, 5.65, 12.58, &  
     14.68,18.77,26.41,33.25,46.77,85.50, &  
     1.47, 2.59, 3.12, 3.74, 4.54, 10.42, &  
     12.24,15.82,22.66,28.91,41.54,79.33/ 

  DATA alpha_aeri_5wv/ &
                                 ! dust insoluble 
        0.759, 0.770, 0.775, 0.775, 0.772, & 
                                 !!jb bc insoluble 
        11.536,10.033, 8.422, 7.234, 6.270, & 
                                 ! pom insoluble 
        5.042, 3.101, 1.890, 1.294, 0.934/ 
   ! 
  DATA cg_aers_5wv/ &  
                                 ! bc soluble 
      .651, .651, .651, .651, .651, .651, & 
      .651, .651, .738, .764, .785, .800, & 
      .597, .597, .597, .597, .597, .597, & 
      .597, .597, .695, .725, .751, .770, & 
      .543, .543, .543, .543, .543, .543, & 
      .543, .543, .650, .684, .714, .736, &  
      .504, .504, .504, .504, .504, .504, & 
      .504, .504, .614, .651, .683, .708, &  
      .469, .469, .469, .469, .469, .469, & 
      .469, .469, .582, .620, .655, .681, & 
                                 ! pom soluble 
      .679, .679, .679, .679, .683, .691, & 
      .703, .720, .736, .751, .766, .784, & 
      .656, .656, .656, .656, .659, .669, & 
      .681, .699, .717, .735, .750, .779, &  
      .623, .623, .623, .623, .627, .637, & 
      .649, .668, .688, .709, .734, .762, & 
      .592, .592, .592, .592, .595, .605, & 
      .618, .639, .660, .682, .711, .743, & 
      .561, .561, .561, .561, .565, .575, & 
      .588, .609, .632, .656, .688, .724, & 
                                 ! Accumulation sulfate 
      .671, .684, .697, .704, .714, .723, & 
      .734, .746, .762, .771, .781, .789, & 
      .653, .666, .678, .687, .697, .707, & 
      .719, .732, .751, .762, .775, .789, & 
      .622, .635, .648, .657, .667, .678, & 
      .691, .705, .728, .741, .758, .777, & 
      .591, .604, .617, .627, .638, .650, & 
      .664, .679, .704, .719, .739, .761, & 
      .560, .574, .587, .597, .609, .621, &  
      .637, .653, .680, .697, .719, .745, & 
                                 ! Coarse sulfate 
      .671, .684, .697, .704, .714, .723, & 
      .734, .746, .762, .771, .781, .789, & 
      .653, .666, .678, .687, .697, .707, & 
      .719, .732, .751, .762, .775, .789, & 
      .622, .635, .648, .657, .667, .678, & 
      .691, .705, .728, .741, .758, .777, & 
      .591, .604, .617, .627, .638, .650, & 
      .664, .679, .704, .719, .739, .761, & 
      .560, .574, .587, .597, .609, .621, &  
      .637, .653, .680, .697, .719, .745, & 
                                 ! For super coarse seasalt (computed below for 550nm!) 
      0.730,0.753,0.760,0.766,0.772,0.793, &  
      0.797,0.802,0.809,0.813,0.820,0.830, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.721,0.744,0.750,0.756,0.762,0.784, &  
      0.787,0.793,0.800,0.804,0.811,0.822, &  
      0.717,0.741,0.747,0.753,0.759,0.780, &  
      0.784,0.789,0.795,0.800,0.806,0.817, &  
      0.715,0.739,0.745,0.751,0.757,0.777, &   
      0.781,0.786,0.793,0.797,0.803,0.814, &  
                                 ! For coarse-soluble seasalt (computed below for 550nm!) 
      0.730,0.753,0.760,0.766,0.772,0.793, &  
      0.797,0.802,0.809,0.813,0.820,0.830, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.721,0.744,0.750,0.756,0.762,0.784, &  
      0.787,0.793,0.800,0.804,0.811,0.822, &  
      0.717,0.741,0.747,0.753,0.759,0.780, &  
      0.784,0.789,0.795,0.800,0.806,0.817, &  
      0.715,0.739,0.745,0.751,0.757,0.777, &   
      0.781,0.786,0.793,0.797,0.803,0.814, &  
                                 ! accumulation-seasalt soluble (computed below for 550nm!)  
      0.698,0.722,0.729,0.736,0.743,0.765, &  
      0.768,0.773,0.777,0.779,0.781,0.779, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.000,0.000,0.000,0.000,0.000,0.000, &  
      0.658,0.691,0.701,0.710,0.720,0.756, &  
      0.763,0.771,0.782,0.788,0.795,0.801, &  
      0.632,0.668,0.679,0.690,0.701,0.743, &  
      0.750,0.762,0.775,0.783,0.792,0.804, &  
      0.605,0.644,0.656,0.669,0.681,0.729, &  
      0.737,0.750,0.765,0.775,0.787,0.803/
 !

  DATA cg_aeri_5wv/&
     ! dust insoluble
     0.714, 0.697, 0.688, 0.683, 0.679, &
     ! bc insoluble
     0.511, 0.445, 0.384, 0.342, 0.307, &
     !c pom insoluble
     0.596, 0.536, 0.466, 0.409, 0.359/
  !
  DATA piz_aers_5wv/&
                           ! bc soluble 
  .445, .445, .445, .445, .445, .445, & 
  .445, .445, .470, .487, .508, .531, & 
  .442, .442, .442, .442, .442, .442, & 
  .442, .442, .462, .481, .506, .533, & 
  .427, .427, .427, .427, .427, .427, & 
  .427, .427, .449, .470, .497, .526, & 
  .413, .413, .413, .413, .413, .413, & 
  .413, .413, .437, .458, .486, .516, & 
  .399, .399, .399, .399, .399, .399, & 
  .399, .399, .423, .445, .473, .506, & 
                           ! pom soluble 
  .975, .975, .975, .975, .975, .977, & 
  .979, .982, .984, .987, .990, .994, & 
  .972, .972, .972, .972, .973, .974, & 
  .977, .980, .983, .986, .989, .993, & 
  .963, .963, .963, .963, .964, .966, & 
  .969, .974, .977, .982, .986, .991, & 
  .955, .955, .955, .955, .955, .958, & 
  .962, .967, .972, .977, .983, .989, & 
  .944, .944, .944, .944, .944, .948, & 
  .952, .959, .962, .972, .979, .987, & 
                           ! sulfate soluble accumulation 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
  1.000,1.000,1.000,1.000,1.000,1.000, & 
                           ! sulfate soluble coarse 
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
                           ! seasalt super coarse (computed below for 550nm) 
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, & 
                           ! seasalt coarse (computed below for 550nm) 
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
                           ! seasalt soluble accumulation (computed below for 550nm) 
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000, &  
  1.000,1.000,1.000,1.000,1.000,1.000/ 

 !
  DATA piz_aeri_5wv/&
     ! dust insoluble
     0.944, 0.970, 0.977, 0.982, 0.987, &
     ! bc insoluble
     0.415, 0.387, 0.355, 0.328, 0.301, &
     ! pom insoluble
     0.972, 0.963, 0.943, 0.923, 0.897/

! Interpolation des coefficients optiques de 19 niveaux vers le nombre des niveaux du model
  IF (firstcall) THEN
     firstcall=.FALSE.
! Allocation
    IF (.NOT. ALLOCATED(A1_ASSSM)) THEN
        ALLOCATE(A1_ASSSM(klev),A2_ASSSM(klev), A3_ASSSM(klev),&
          B1_ASSSM(klev), B2_ASSSM(klev), C1_ASSSM(klev), C2_ASSSM(klev),&
          A1_CSSSM(klev), A2_CSSSM(klev), A3_CSSSM(klev),&
          B1_CSSSM(klev), B2_CSSSM(klev), C1_CSSSM(klev), C2_CSSSM(klev),&
          A1_SSSSM(klev), A2_SSSSM(klev), A3_SSSSM(klev),&
          B1_SSSSM(klev), B2_SSSSM(klev), C1_SSSSM(klev), C2_SSSSM(klev), stat=ierr)
        IF (ierr /= 0) CALL abort_gcm('aeropt_5mw', 'pb in allocation 1',1)
     END IF

!Accumulation mode
     CALL pres2lev(A1_ASSSM_19, A1_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_ASSSM_19, A2_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_ASSSM_19, A3_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_ASSSM_19, B1_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_ASSSM_19, B2_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_ASSSM_19, C1_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_ASSSM_19, C2_ASSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
!Coarse mode
     CALL pres2lev(A1_CSSSM_19, A1_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_CSSSM_19, A2_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_CSSSM_19, A3_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_CSSSM_19, B1_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_CSSSM_19, B2_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_CSSSM_19, C1_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_CSSSM_19, C2_CSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
!Super coarse mode
     CALL pres2lev(A1_SSSSM_19, A1_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_SSSSM_19, A2_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_SSSSM_19, A3_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_SSSSM_19, B1_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_SSSSM_19, B2_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_SSSSM_19, C1_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_SSSSM_19, C2_SSSSM, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

  END IF ! firstcall


  ! Initialisations
  ai(:) = 0.
  tausum(:,:,:) = 0.


  DO k=1, klev
    DO i=1, klon
!      IF (t_seri(i,k).EQ.0) stop 'stop aeropt_5wv T '
!      IF (pplay(i,k).EQ.0) stop  'stop aeropt_5wv p '
      zrho=pplay(i,k)/t_seri(i,k)/RD                  ! kg/m3
      dh(i,k)=pdel(i,k)/(gravit*zrho)
!CDIR UNROLL=naero_spc
      mass_temp(i,k,:) = m_allaer(i,k,:) / zrho / 1.e+9
      zdp1(i,k)=pdel(i,k)/(gravit*delt)     ! air mass auxiliary  variable --> zdp1 [kg/(m^2 *s)]

    ENDDO
  ENDDO


  IF (flag_aerosol .EQ. 1) THEN 
     nb_aer = 2
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_ASSO4M
     aerosol_name(2) = id_CSSO4M
  ELSEIF (flag_aerosol .EQ. 2) THEN
     nb_aer = 2
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_ASBCM
     aerosol_name(2) = id_AIBCM
  ELSEIF (flag_aerosol .EQ. 3) THEN 
     nb_aer = 2
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_ASPOMM
     aerosol_name(2) = id_AIPOMM
  ELSEIF (flag_aerosol .EQ. 4) THEN 
     nb_aer = 3
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_CSSSM
     aerosol_name(2) = id_SSSSM
     aerosol_name(3) = id_ASSSM
  ELSEIF (flag_aerosol .EQ. 5) THEN 
     nb_aer = 1
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_CIDUSTM
  ELSEIF (flag_aerosol .EQ. 6) THEN 
     nb_aer = 10
     ALLOCATE (aerosol_name(nb_aer)) 
     aerosol_name(1) = id_ASSO4M      
     aerosol_name(2) = id_ASBCM
     aerosol_name(3) = id_AIBCM
     aerosol_name(4) = id_ASPOMM
     aerosol_name(5) = id_AIPOMM
     aerosol_name(6) = id_CSSSM
     aerosol_name(7) = id_SSSSM
     aerosol_name(8) = id_ASSSM
     aerosol_name(9) = id_CIDUSTM
     aerosol_name(10) = id_CSSO4M
  ENDIF

  ! 
  ! loop over modes, use of precalculated nmd and corresponding sigma
  !    loop over wavelengths
  !    for each mass species in mode
  !      interpolate from Sext to retrieve Sext_at_gridpoint_per_species
  !      compute optical_thickness_at_gridpoint_per_species
  

  !
  ! Calculations that need to be done since we are not in the subroutines INCA
  !      

!CDIR ON_ADB(RH_tab)
!CDIR ON_ADB(fact_RH)
!CDIR NOVECTOR
  DO n=1,nbre_RH-1
    fact_RH(n)=1./(RH_tab(n+1)-RH_tab(n))
  ENDDO
   
  DO k=1, KLEV
!CDIR ON_ADB(RH_tab)
!CDIR ON_ADB(fact_RH)
    DO i=1, KLON
      rh(i,k)=MIN(RHcl(i,k)*100.,RH_MAX)
      RH_num(i,k) = INT( rh(i,k)/10. + 1.)
      IF (rh(i,k).GT.85.) RH_num(i,k)=10
      IF (rh(i,k).GT.90.) RH_num(i,k)=11
      DELTA(i,k)=(rh(i,k)-RH_tab(RH_num(i,k)))*fact_RH(RH_num(i,k))
    ENDDO
  ENDDO

!CDIR SHORTLOOP  
  used_tau(:)=.FALSE.
    
  DO m=1,nb_aer   ! tau is only computed for each mass    
    fac=1.0
    IF (aerosol_name(m).EQ.id_ASBCM) THEN
        soluble=.TRUE.
        spsol=1
        spss=0
    ELSEIF (aerosol_name(m).EQ.id_ASPOMM) THEN 
        soluble=.TRUE.
        spsol=2 
        spss=0
    ELSEIF (aerosol_name(m).EQ.id_ASSO4M) THEN
        soluble=.TRUE.
        spsol=3
        spss=0
        fac=1.375    ! (NH4)2-SO4/SO4 132/96 mass conversion factor for OD
    ELSEIF (aerosol_name(m).EQ.id_CSSO4M) THEN
        soluble=.TRUE.
        spsol=4
        spss=0
        fac=1.375    ! (NH4)2-SO4/SO4 132/96 mass conversion factor for OD
    ELSEIF (aerosol_name(m).EQ.id_SSSSM) THEN 
        soluble=.TRUE.
        spsol=5
        spss=3
    ELSEIF (aerosol_name(m).EQ.id_CSSSM) THEN 
        soluble=.TRUE.
        spsol=6
        spss=2
    ELSEIF (aerosol_name(m).EQ.id_ASSSM) THEN
        soluble=.TRUE.
        spsol=7
        spss=1
    ELSEIF (aerosol_name(m).EQ.id_CIDUSTM) THEN 
        soluble=.FALSE.
        spinsol=1
        spss=0
    ELSEIF  (aerosol_name(m).EQ.id_AIBCM) THEN 
        soluble=.FALSE.
        spinsol=2
        spss=0
    ELSEIF (aerosol_name(m).EQ.id_AIPOMM) THEN 
        soluble=.FALSE.
        spinsol=3
        spss=0
    ELSE 
        CYCLE
    ENDIF

!Bug 21 12 10 AI
!    used_tau(spsol)=.TRUE.
    IF (soluble) then
      used_tau(spsol)=.TRUE.
       ELSE
      used_tau(naero_soluble+spinsol)=.TRUE.
    ENDIF

    DO la=1,las

      IF (soluble) THEN

        IF((la.EQ.2).AND.(spss.NE.0)) THEN !la=2 corresponds to 550 nm
          IF (spss.EQ.1) THEN !accumulation mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_ASSSM)
!CDIR ON_ADB(A2_ASSSM)
!CDIR ON_ADB(A3_ASSSM)
              DO i=1, KLON
                H=rh(i,k)/100
                tau_ae5wv_int=A1_ASSSM(k)+A2_ASSSM(k)*H+A3_ASSSM(k)/(H-1.05)
                tau(i,k,la,spsol) = mass_temp(i,k,spsol)*1000.*zdp1(i,k)   &
                                   *tau_ae5wv_int*delt*fac
                tausum(i,la,spsol)=tausum(i,la,spsol)+tau(i,k,la,spsol)
              ENDDO
            ENDDO
          ENDIF
  
          IF (spss.EQ.2) THEN !coarse mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_CSSSM)
!CDIR ON_ADB(A2_CSSSM)
!CDIR ON_ADB(A3_CSSSM)
              DO i=1, KLON
                H=rh(i,k)/100
                tau_ae5wv_int=A1_CSSSM(k)+A2_CSSSM(k)*H+A3_CSSSM(k)/(H-1.05)
                tau(i,k,la,spsol) = mass_temp(i,k,spsol)*1000.*zdp1(i,k)  &
                                   *tau_ae5wv_int*delt*fac
                tausum(i,la,spsol) = tausum(i,la,spsol)+tau(i,k,la,spsol)
              ENDDO
            ENDDO
          ENDIF

          IF (spss.EQ.3) THEN !super coarse mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_SSSSM)
!CDIR ON_ADB(A2_SSSSM)
!CDIR ON_ADB(A3_SSSSM)
              DO i=1, KLON
                H=rh(i,k)/100
                tau_ae5wv_int=A1_SSSSM(k)+A2_SSSSM(k)*H+A3_SSSSM(k)/(H-1.05)
                tau(i,k,la,spsol) = mass_temp(i,k,spsol)*1000.*zdp1(i,k)  &
                                   *tau_ae5wv_int*delt*fac
                tausum(i,la,spsol)=tausum(i,la,spsol)+tau(i,k,la,spsol)
              ENDDO
            ENDDO
          ENDIF

        ELSE
          DO k=1, KLEV
!CDIR ON_ADB(alpha_aers_5wv)
            DO i=1, KLON
              tau_ae5wv_int = alpha_aers_5wv(RH_num(i,k),la,spsol)+DELTA(i,k)* &
                             (alpha_aers_5wv(RH_num(i,k)+1,la,spsol) - & 
                              alpha_aers_5wv(RH_num(i,k),la,spsol))

              tau(i,k,la,spsol) = mass_temp(i,k,spsol)*1000.*zdp1(i,k)   &
                                 *tau_ae5wv_int*delt*fac
              tausum(i,la,spsol)=tausum(i,la,spsol)+tau(i,k,la,spsol)
            ENDDO
          ENDDO
        ENDIF

      ELSE                                                  ! For insoluble aerosol
        DO k=1, KLEV
!CDIR ON_ADB(alpha_aeri_5wv)
          DO i=1, KLON
            tau_ae5wv_int = alpha_aeri_5wv(la,spinsol)
            tau(i,k,la,naero_soluble+spinsol) = mass_temp(i,k,naero_soluble+spinsol)*1000.*zdp1(i,k)* &
                                                tau_ae5wv_int*delt*fac
            tausum(i,la,naero_soluble+spinsol)= tausum(i,la,naero_soluble+spinsol)  &
                                               +tau(i,k,la,naero_soluble+spinsol)
          ENDDO
        ENDDO
      ENDIF
    ENDDO   ! boucle sur les longueurs d'onde
  ENDDO     ! Boucle  sur les masses de traceurs

  DO m=1,naero_spc
    IF (.NOT.used_tau(m)) tau(:,:,:,m)=0.
  ENDDO  
!
!
!  taue670(:) = SUM(tausum(:,la670,:),dim=2) 
!  taue865(:) = SUM(tausum(:,la865,:),dim=2) 
!
!  DO i=1, klon
!    ai(i)=-LOG(MAX(taue670(i),0.0001)/ &
!       MAX(taue865(i),0.0001))/LOG(670./865.)
!  ENDDO

  DO i=1, klon
     od550aer(i)=0.
     DO m=1,naero_spc
        od550aer(i)=od550aer(i)+tausum(i,2,m)
     END DO
  END DO
  DO i=1, klon
     od865aer(i)=0.
     DO m=1,naero_spc
        od865aer(i)=od865aer(i)+tausum(i,5,m)
     END DO
  END DO
  DO i=1, klon
     DO k=1, KLEV
        ec550aer(i,k)=0.
        DO m=1,naero_spc
           ec550aer(i,k)=ec550aer(i,k)+tau(i,k,2,m)/dh(i,k)
        END DO
     END DO
  END DO
  
   od550lt1aer(:)=tausum(:,2,id_ASSO4M)+tausum(:,2,id_ASBCM)+tausum(:,2,id_AIBCM)+ &
	tausum(:,2,id_ASPOMM)+tausum(:,2,id_AIPOMM)+tausum(:,2,id_ASSSM)+ &
	0.03*tausum(:,2,id_CSSSM)+0.4*tausum(:,2,id_CIDUSTM)



  DEALLOCATE(aerosol_name) 
  
END SUBROUTINE AEROPT_5WV
