!
! $Id: aeropt_2bands.F90 1337 2010-04-02 11:31:05Z fairhead $
!
SUBROUTINE AEROPT_2BANDS( &
     pdel, m_allaer, delt, RHcl, &
     tau_allaer, piz_allaer, &
     cg_allaer, m_allaer_pi, &
     flag_aerosol, pplay, t_seri, presnivs)

  USE dimphy
  USE aero_mod
  USE phys_local_var_mod, only: absvisaer

  !    Yves Balkanski le 12 avril 2006
  !    Celine Deandreis
  !    Anne Cozic Avril 2009
  !    a partir d'une sous-routine de Johannes Quaas pour les sulfates
  !
  IMPLICIT NONE

  INCLUDE "YOMCST.h"
  INCLUDE "iniprint.h"

  !
  ! Input arguments:
  !
  REAL, DIMENSION(klon,klev),     INTENT(in)  :: pdel
  REAL,                           INTENT(in)  :: delt
  REAL, DIMENSION(klon,klev,naero_spc),   INTENT(in)  :: m_allaer
!RAF
  REAL, DIMENSION(klon,klev,naero_spc),   INTENT(in)  :: m_allaer_pi
  REAL, DIMENSION(klon,klev),     INTENT(in)  :: RHcl       ! humidite relative ciel clair
!RAF  REAL, DIMENSION(klon,naero_tot),INTENT(in)  :: fractnat_allaer
  INTEGER,                        INTENT(in)  :: flag_aerosol
  REAL, DIMENSION(klon,klev),     INTENT(in)  :: pplay
  REAL, DIMENSION(klon,klev),     INTENT(in)  :: t_seri
  REAL, DIMENSION(klev),          INTENT(in)  :: presnivs
  !
  ! Output arguments:
  !
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(out) :: tau_allaer ! epaisseur optique aerosol
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(out) :: piz_allaer ! single scattering albedo aerosol
  REAL, DIMENSION(klon,klev,naero_grp,nbands), INTENT(out) :: cg_allaer  ! asymmetry parameter aerosol

  !
  ! Local
  !
  REAL, DIMENSION(klon,klev,naero_tot,nbands) ::  tau_ae
!RAF
  REAL, DIMENSION(klon,klev,naero_tot,nbands) ::  tau_ae_pi
  REAL, DIMENSION(klon,klev,naero_tot,nbands) ::  piz_ae
  REAL, DIMENSION(klon,klev,naero_tot,nbands) ::  cg_ae
  LOGICAL ::  soluble
  INTEGER :: i, k,n, ierr, inu, m, mrfspecies
  INTEGER :: spsol, spinsol, spss
  INTEGER :: RH_num(klon,klev)
  INTEGER, PARAMETER :: nb_level=19 ! number of vertical levels in DATA

  INTEGER, PARAMETER :: nbre_RH=12
  INTEGER, PARAMETER :: naero_soluble=7    ! 1- BC soluble; 2- POM soluble; 3- SO4. acc. 4- SO4 coarse
                                           ! 5- seasalt super coarse  6- seasalt coarse   7- seasalt acc.
  INTEGER, PARAMETER :: naero_insoluble=3  ! 1- Dust; 2- BC insoluble; 3- POM insoluble
  LOGICAL, SAVE :: firstcall=.TRUE. 
!$OMP THREADPRIVATE(firstcall)

! Coefficient optiques sur 19 niveaux
  REAL, SAVE, DIMENSION(nb_level) :: presnivs_19  ! Pression milieux couche pour 19 niveaux (nb_level)
!$OMP THREADPRIVATE(presnivs_19)

  REAL, SAVE, DIMENSION(nb_level) :: A1_ASSSM_b1_19, A2_ASSSM_b1_19, A3_ASSSM_b1_19,&
          B1_ASSSM_b1_19, B2_ASSSM_b1_19, C1_ASSSM_b1_19, C2_ASSSM_b1_19,&
          A1_CSSSM_b1_19, A2_CSSSM_b1_19, A3_CSSSM_b1_19,&
          B1_CSSSM_b1_19, B2_CSSSM_b1_19, C1_CSSSM_b1_19, C2_CSSSM_b1_19,&
          A1_SSSSM_b1_19, A2_SSSSM_b1_19, A3_SSSSM_b1_19,&
          B1_SSSSM_b1_19, B2_SSSSM_b1_19, C1_SSSSM_b1_19, C2_SSSSM_b1_19,&
          A1_ASSSM_b2_19, A2_ASSSM_b2_19, A3_ASSSM_b2_19,&
          B1_ASSSM_b2_19, B2_ASSSM_b2_19, C1_ASSSM_b2_19, C2_ASSSM_b2_19,&
          A1_CSSSM_b2_19, A2_CSSSM_b2_19, A3_CSSSM_b2_19,&
          B1_CSSSM_b2_19, B2_CSSSM_b2_19, C1_CSSSM_b2_19, C2_CSSSM_b2_19,&
          A1_SSSSM_b2_19, A2_SSSSM_b2_19, A3_SSSSM_b2_19,&
          B1_SSSSM_b2_19, B2_SSSSM_b2_19, C1_SSSSM_b2_19, C2_SSSSM_b2_19
!$OMP THREADPRIVATE(A1_ASSSM_b1_19, A2_ASSSM_b1_19, A3_ASSSM_b1_19)
!$OMP THREADPRIVATE(B1_ASSSM_b1_19, B2_ASSSM_b1_19, C1_ASSSM_b1_19, C2_ASSSM_b1_19)
!$OMP THREADPRIVATE(A1_CSSSM_b1_19, A2_CSSSM_b1_19, A3_CSSSM_b1_19)
!$OMP THREADPRIVATE(B1_CSSSM_b1_19, B2_CSSSM_b1_19, C1_CSSSM_b1_19, C2_CSSSM_b1_19)
!$OMP THREADPRIVATE(A1_SSSSM_b1_19, A2_SSSSM_b1_19, A3_SSSSM_b1_19)
!$OMP THREADPRIVATE(B1_SSSSM_b1_19, B2_SSSSM_b1_19, C1_SSSSM_b1_19, C2_SSSSM_b1_19)
!$OMP THREADPRIVATE(A1_ASSSM_b2_19, A2_ASSSM_b2_19, A3_ASSSM_b2_19)
!$OMP THREADPRIVATE(B1_ASSSM_b2_19, B2_ASSSM_b2_19, C1_ASSSM_b2_19, C2_ASSSM_b2_19)
!$OMP THREADPRIVATE(A1_CSSSM_b2_19, A2_CSSSM_b2_19, A3_CSSSM_b2_19)
!$OMP THREADPRIVATE(B1_CSSSM_b2_19, B2_CSSSM_b2_19, C1_CSSSM_b2_19, C2_CSSSM_b2_19)
!$OMP THREADPRIVATE(A1_SSSSM_b2_19, A2_SSSSM_b2_19, A3_SSSSM_b2_19)
!$OMP THREADPRIVATE(B1_SSSSM_b2_19, B2_SSSSM_b2_19, C1_SSSSM_b2_19, C2_SSSSM_b2_19)


! Coefficient optiques interpole sur le nombre de niveau du modele
  REAL, ALLOCATABLE, DIMENSION(:), SAVE :: &
          A1_ASSSM_b1, A2_ASSSM_b1, A3_ASSSM_b1,&
          B1_ASSSM_b1, B2_ASSSM_b1, C1_ASSSM_b1, C2_ASSSM_b1,&
          A1_CSSSM_b1, A2_CSSSM_b1, A3_CSSSM_b1,&
          B1_CSSSM_b1, B2_CSSSM_b1, C1_CSSSM_b1, C2_CSSSM_b1,&
          A1_SSSSM_b1, A2_SSSSM_b1, A3_SSSSM_b1,&
          B1_SSSSM_b1, B2_SSSSM_b1, C1_SSSSM_b1, C2_SSSSM_b1,&
          A1_ASSSM_b2, A2_ASSSM_b2, A3_ASSSM_b2,&
          B1_ASSSM_b2, B2_ASSSM_b2, C1_ASSSM_b2, C2_ASSSM_b2,&
          A1_CSSSM_b2, A2_CSSSM_b2, A3_CSSSM_b2,&
          B1_CSSSM_b2, B2_CSSSM_b2, C1_CSSSM_b2, C2_CSSSM_b2,&
          A1_SSSSM_b2, A2_SSSSM_b2, A3_SSSSM_b2,&
          B1_SSSSM_b2, B2_SSSSM_b2, C1_SSSSM_b2, C2_SSSSM_b2
!$OMP THREADPRIVATE(A1_ASSSM_b1, A2_ASSSM_b1, A3_ASSSM_b1)
!$OMP THREADPRIVATE(B1_ASSSM_b1, B2_ASSSM_b1, C1_ASSSM_b1, C2_ASSSM_b1)
!$OMP THREADPRIVATE(A1_CSSSM_b1, A2_CSSSM_b1, A3_CSSSM_b1)
!$OMP THREADPRIVATE(B1_CSSSM_b1, B2_CSSSM_b1, C1_CSSSM_b1, C2_CSSSM_b1)
!$OMP THREADPRIVATE(A1_SSSSM_b1, A2_SSSSM_b1, A3_SSSSM_b1)
!$OMP THREADPRIVATE(B1_SSSSM_b1, B2_SSSSM_b1, C1_SSSSM_b1, C2_SSSSM_b1)
!$OMP THREADPRIVATE(A1_ASSSM_b2, A2_ASSSM_b2, A3_ASSSM_b2)
!$OMP THREADPRIVATE(B1_ASSSM_b2, B2_ASSSM_b2, C1_ASSSM_b2, C2_ASSSM_b2)
!$OMP THREADPRIVATE(A1_CSSSM_b2, A2_CSSSM_b2, A3_CSSSM_b2)
!$OMP THREADPRIVATE(B1_CSSSM_b2, B2_CSSSM_b2, C1_CSSSM_b2, C2_CSSSM_b2)
!$OMP THREADPRIVATE(A1_SSSSM_b2, A2_SSSSM_b2, A3_SSSSM_b2)
!$OMP THREADPRIVATE(B1_SSSSM_b2, B2_SSSSM_b2, C1_SSSSM_b2, C2_SSSSM_b2)
  
  REAL,PARAMETER :: RH_tab(nbre_RH)=(/0.,10.,20.,30.,40.,50.,60.,70.,80.,85.,90.,95./)
  REAL, PARAMETER :: RH_MAX=95.
  REAL:: DELTA(klon,klev), rh(klon,klev), H
  REAL:: tau_ae2b_int   ! Intermediate computation of epaisseur optique aerosol
  REAL:: piz_ae2b_int   ! Intermediate computation of Single scattering albedo
  REAL:: cg_ae2b_int    ! Intermediate computation of Assymetry parameter
  REAL :: Fact_RH(nbre_RH)
  REAL :: zrho
  REAL :: fac
  REAL :: zdp1(klon,klev) 
  REAL, PARAMETER ::  gravit = 9.80616    ! m2/s
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: aerosol_name
  INTEGER :: nb_aer
  REAL, DIMENSION(klon,klev,naero_spc) :: mass_temp
!RAF
  REAL, DIMENSION(klon,klev,naero_spc) :: mass_temp_pi

  !
  ! Proprietes optiques
  !
  REAL:: alpha_aers_2bands(nbre_RH,nbands,naero_soluble)   !--unit m2/g SO4
  REAL:: alpha_aeri_2bands(nbands,naero_insoluble)
  REAL:: cg_aers_2bands(nbre_RH,nbands,naero_soluble)      !--unit 
  REAL:: cg_aeri_2bands(nbands,naero_insoluble)
  REAL:: piz_aers_2bands(nbre_RH,nbands,naero_soluble)     !-- unit
  REAL:: piz_aeri_2bands(nbands,naero_insoluble)           !-- unit

  INTEGER :: id
  LOGICAL :: used_aer(naero_tot)
  REAL :: tmp_var, tmp_var_pi

  DATA presnivs_19/&
       100426.5,  98327.6, 95346.5, 90966.8, 84776.9, &
       76536.5,   66292.2, 54559.3, 42501.8, 31806, &
       23787.5,   18252.7, 13996,   10320.8, 7191.1, &
       4661.7,    2732.9,  1345.6,  388.2/  


!***********************BAND 1***********************************
!ACCUMULATION MODE
  DATA A1_ASSSM_b1_19/ 4.373E+00,  4.361E+00,  4.331E+00, &
                    4.278E+00,  4.223E+00,  4.162E+00, &
                    4.103E+00,  4.035E+00,  3.962E+00, &
                    3.904E+00,  3.871E+00,  3.847E+00, &
                    3.824E+00,  3.780E+00,  3.646E+00, &
                    3.448E+00,  3.179E+00,  2.855E+00,  2.630E+00/
  DATA A2_ASSSM_b1_19/ 2.496E+00,  2.489E+00,  2.472E+00, &
                    2.442E+00,  2.411E+00,  2.376E+00, &
                    2.342E+00,  2.303E+00,  2.261E+00, &
                    2.228E+00,  2.210E+00,  2.196E+00, &
                    2.183E+00,  2.158E+00,  2.081E+00, &
                    1.968E+00,  1.814E+00,  1.630E+00,  1.501E+00/
  DATA A3_ASSSM_b1_19/-4.688E-02, -4.676E-02, -4.644E-02, &
                   -4.587E-02, -4.528E-02, -4.463E-02, &
                   -4.399E-02, -4.326E-02, -4.248E-02, &
                   -4.186E-02, -4.151E-02, -4.125E-02, &
                   -4.100E-02, -4.053E-02, -3.910E-02, &
                   -3.697E-02, -3.408E-02, -3.061E-02, -2.819E-02/
  DATA B1_ASSSM_b1_19/ 1.165E-08,  1.145E-08,  1.097E-08, &
                    1.012E-08,  9.233E-09,  8.261E-09, &
                    7.297E-09,  6.201E-09,  5.026E-09, &
                    4.098E-09,  3.567E-09,  3.187E-09, &
                    2.807E-09,  2.291E-09,  2.075E-09, &
                    1.756E-09,  1.322E-09,  8.011E-10, 4.379E-10/
  DATA B2_ASSSM_b1_19/ 2.193E-08,  2.192E-08,  2.187E-08, &
                    2.179E-08,  2.171E-08,  2.162E-08, &
                    2.153E-08,  2.143E-08,  2.132E-08, &
                    2.124E-08,  2.119E-08,  2.115E-08, &
                    2.112E-08,  2.106E-08,  2.100E-08, &
                    2.090E-08,  2.077E-08,  2.061E-08,  2.049E-08/
  DATA C1_ASSSM_b1_19/ 7.365E-01,  7.365E-01,  7.365E-01, &
                    7.364E-01,  7.363E-01,  7.362E-01, &
                    7.361E-01,  7.359E-01,  7.358E-01, &
                    7.357E-01,  7.356E-01,  7.356E-01, &
                    7.356E-01,  7.355E-01,  7.354E-01, &
                    7.352E-01,  7.350E-01,  7.347E-01,  7.345E-01/
  DATA C2_ASSSM_b1_19/ 5.833E-02,  5.835E-02,  5.841E-02, &
                    5.850E-02,  5.859E-02,  5.870E-02, &
                    5.880E-02,  5.891E-02,  5.904E-02, &
                    5.914E-02,  5.920E-02,  5.924E-02, &
                    5.928E-02,  5.934E-02,  5.944E-02, &
                    5.959E-02,  5.979E-02,  6.003E-02,  6.020E-02/
!COARSE MODE
  DATA A1_CSSSM_b1_19/ 7.403E-01,  7.422E-01,  7.626E-01, &
                    8.019E-01,  8.270E-01,  8.527E-01, &
                    8.702E-01,  8.806E-01,  8.937E-01, &
                    9.489E-01,  1.030E+00,  1.105E+00, &
                    1.199E+00,  1.357E+00,  1.660E+00, &
                    2.540E+00,  4.421E+00,  2.151E+00,  9.518E-01/
  DATA A2_CSSSM_b1_19/ 4.522E-01,  4.532E-01,  4.644E-01, &
                    4.859E-01,  4.996E-01,  5.137E-01, &
                    5.233E-01,  5.290E-01,  5.361E-01, &
                    5.655E-01,  6.085E-01,  6.483E-01, &
                    6.979E-01,  7.819E-01,  9.488E-01, &
                    1.450E+00,  2.523E+00,  1.228E+00,  5.433E-01/
  DATA A3_CSSSM_b1_19/-8.516E-03, -8.535E-03, -8.744E-03, &
                   -9.148E-03, -9.406E-03, -9.668E-03, &
                   -9.848E-03, -9.955E-03, -1.009E-02, &
                   -1.064E-02, -1.145E-02, -1.219E-02, &
                   -1.312E-02, -1.470E-02, -1.783E-02, &
                   -2.724E-02, -4.740E-02, -2.306E-02, -1.021E-02/
  DATA B1_CSSSM_b1_19/ 2.535E-07,  2.530E-07,  2.479E-07, &
                    2.380E-07,  2.317E-07,  2.252E-07, &
                    2.208E-07,  2.182E-07,  2.149E-07, &
                    2.051E-07,  1.912E-07,  1.784E-07, &
                    1.624E-07,  1.353E-07,  1.012E-07, &
                    6.016E-08,  2.102E-08,  0.000E+00,  0.000E+00/
  DATA B2_CSSSM_b1_19/ 1.221E-07,  1.217E-07,  1.179E-07, &
                    1.104E-07,  1.056E-07,  1.008E-07, &
                    9.744E-08,  9.546E-08,  9.299E-08, &
                    8.807E-08,  8.150E-08,  7.544E-08, &
                    6.786E-08,  5.504E-08,  4.080E-08, &
                    2.960E-08,  2.300E-08,  2.030E-08,  1.997E-08/
  DATA C1_CSSSM_b1_19/ 7.659E-01,  7.658E-01,  7.652E-01, &
                    7.639E-01,  7.631E-01,  7.623E-01, &
                    7.618E-01,  7.614E-01,  7.610E-01, &
                    7.598E-01,  7.581E-01,  7.566E-01, &
                    7.546E-01,  7.513E-01,  7.472E-01, &
                    7.423E-01,  7.376E-01,  7.342E-01,  7.334E-01/
  DATA C2_CSSSM_b1_19/ 3.691E-02,  3.694E-02,  3.729E-02, &
                    3.796E-02,  3.839E-02,  3.883E-02, &
                    3.913E-02,  3.931E-02,  3.953E-02, &
                    4.035E-02,  4.153E-02,  4.263E-02, &
                    4.400E-02,  4.631E-02,  4.933E-02, &
                    5.331E-02,  5.734E-02,  6.053E-02,  6.128E-02/
!SUPER COARSE MODE
  DATA A1_SSSSM_b1_19/ 2.836E-01,  2.876E-01,  2.563E-01, &
                    2.414E-01,  2.541E-01,  2.546E-01, &
                    2.572E-01,  2.638E-01,  2.781E-01, &
                    3.167E-01,  4.209E-01,  5.286E-01, &
                    6.959E-01,  9.233E-01,  1.282E+00, &
                    1.836E+00,  2.981E+00,  4.355E+00,  4.059E+00/
  DATA A2_SSSSM_b1_19/ 1.608E-01,  1.651E-01,  1.577E-01, &
                    1.587E-01,  1.686E-01,  1.690E-01, &
                    1.711E-01,  1.762E-01,  1.874E-01, &
                    2.138E-01,  2.751E-01,  3.363E-01, &
                    4.279E-01,  5.519E-01,  7.421E-01, &
                    1.048E+00,  1.702E+00,  2.485E+00,  2.317E+00/
  DATA A3_SSSSM_b1_19/-3.025E-03, -3.111E-03, -2.981E-03, &
                   -3.005E-03, -3.193E-03, -3.200E-03, &
                   -3.239E-03, -3.336E-03, -3.548E-03, &
                   -4.047E-03, -5.196E-03, -6.345E-03, &
                   -8.061E-03, -1.038E-02, -1.395E-02, &
                   -1.970E-02, -3.197E-02, -4.669E-02, -4.352E-02/
  DATA B1_SSSSM_b1_19/ 6.759E-07,  6.246E-07,  5.542E-07, &
                    4.953E-07,  4.746E-07,  4.738E-07, &
                    4.695E-07,  4.588E-07,  4.354E-07, &
                    3.947E-07,  3.461E-07,  3.067E-07, &
                    2.646E-07,  2.095E-07,  1.481E-07, &
                    9.024E-08,  5.747E-08,  2.384E-08,  6.599E-09/
  DATA B2_SSSSM_b1_19/ 5.977E-07,  5.390E-07,  4.468E-07, &
                    3.696E-07,  3.443E-07,  3.433E-07, &
                    3.380E-07,  3.249E-07,  2.962E-07, &
                    2.483E-07,  1.989E-07,  1.623E-07, &
                    1.305E-07,  9.015E-08,  6.111E-08, &
                    3.761E-08,  2.903E-08,  2.337E-08,  2.147E-08/
  DATA C1_SSSSM_b1_19/ 8.120E-01,  8.084E-01,  8.016E-01, &
                    7.953E-01,  7.929E-01,  7.928E-01, &
                    7.923E-01,  7.910E-01,  7.882E-01, &
                    7.834E-01,  7.774E-01,  7.725E-01, &
                    7.673E-01,  7.604E-01,  7.529E-01, &
                    7.458E-01,  7.419E-01,  7.379E-01,  7.360E-01/
  DATA C2_SSSSM_b1_19/ 2.388E-02,  2.392E-02,  2.457E-02,  2.552E-02, &
                    2.615E-02,  2.618E-02,  2.631E-02,  2.663E-02, &
                    2.735E-02,  2.875E-02,  3.113E-02,  3.330E-02, &
                    3.615E-02,  3.997E-02,  4.521E-02,  5.038E-02, &
                    5.358E-02,  5.705E-02,  5.887E-02/
!*********************BAND 2************************************************
!ACCUMULATION MODE
  DATA A1_ASSSM_b2_19/1.256E+00, 1.246E+00, 1.226E+00, 1.187E+00, 1.148E+00, &
                   1.105E+00, 1.062E+00, 1.014E+00, 9.616E-01, 9.205E-01, &
                   8.970E-01, 8.800E-01, 8.632E-01, 8.371E-01, 7.943E-01, &
                   7.308E-01, 6.448E-01, 5.414E-01, 4.693E-01/
  DATA A2_ASSSM_b2_19/5.321E-01, 5.284E-01, 5.196E-01, 5.036E-01, 4.872E-01, &
                   4.691E-01, 4.512E-01, 4.308E-01, 4.089E-01, 3.917E-01, &
                   3.818E-01, 3.747E-01, 3.676E-01, 3.567E-01, 3.385E-01, &
                   3.116E-01, 2.751E-01, 2.312E-01, 2.006E-01/
  DATA A3_ASSSM_b2_19/-1.053E-02, -1.046E-02, -1.028E-02, -9.964E-03, -9.637E-03, &
                   -9.279E-03, -8.923E-03, -8.518E-03, -8.084E-03, -7.741E-03, &
                   -7.545E-03, -7.405E-03, -7.265E-03, -7.048E-03, -6.687E-03, &
                   -6.156E-03, -5.433E-03, -4.565E-03, -3.961E-03/
  DATA B1_ASSSM_b2_19/1.560E-02, 1.560E-02, 1.561E-02, 1.565E-02, 1.568E-02, &
                   1.572E-02, 1.576E-02, 1.580E-02, 1.584E-02, 1.588E-02, &
                   1.590E-02, 1.592E-02, 1.593E-02, 1.595E-02, 1.599E-02, &
                   1.605E-02, 1.612E-02, 1.621E-02, 1.627E-02/
  DATA B2_ASSSM_b2_19/1.073E-02, 1.074E-02, 1.076E-02, 1.079E-02, 1.082E-02, &
                   1.085E-02, 1.089E-02, 1.093E-02, 1.097E-02, 1.100E-02, &
                   1.102E-02, 1.103E-02, 1.105E-02, 1.107E-02, 1.110E-02, &
                   1.115E-02, 1.122E-02, 1.130E-02, 1.136E-02/
  DATA C1_ASSSM_b2_19/7.429E-01, 7.429E-01, 7.429E-01, 7.427E-01, 7.427E-01, &
                   7.424E-01, 7.423E-01, 7.422E-01, 7.421E-01, 7.420E-01, &
                   7.419E-01, 7.419E-01, 7.418E-01, 7.417E-01, 7.416E-01, &
                   7.415E-01, 7.413E-01, 7.409E-01, 7.408E-01/
  DATA C2_ASSSM_b2_19/3.031E-02, 3.028E-02, 3.022E-02, 3.011E-02, 2.999E-02, &
                   2.986E-02, 2.973E-02, 2.959E-02, 2.943E-02, 2.931E-02, &
                   2.924E-02, 2.919E-02, 2.913E-02, 2.905E-02, 2.893E-02, &
                   2.874E-02, 2.847E-02, 2.817E-02, 2.795E-02/
!COARSE MODE
  DATA A1_CSSSM_b2_19/7.061E-01, 7.074E-01, 7.211E-01, 7.476E-01, 7.647E-01, &
                   7.817E-01, 7.937E-01, 8.007E-01, 8.095E-01, 8.436E-01, &
                   8.932E-01, 9.390E-01, 9.963E-01, 1.093E+00, 1.256E+00, &
                   1.668E+00, 1.581E+00, 3.457E-01, 1.331E-01/
  DATA A2_CSSSM_b2_19/3.617E-01, 3.621E-01, 3.662E-01, 3.739E-01, 3.789E-01, &
                   3.840E-01, 3.874E-01, 3.895E-01, 3.921E-01, 4.001E-01, &
                   4.117E-01, 4.223E-01, 4.356E-01, 4.581E-01, 5.099E-01, &
                   6.831E-01, 6.663E-01, 1.481E-01, 5.703E-02/
  DATA A3_CSSSM_b2_19/-6.953E-03, -6.961E-03, -7.048E-03, -7.216E-03, -7.322E-03, &
                   -7.431E-03, -7.506E-03, -7.551E-03, -7.606E-03, -7.791E-03, &
                   -8.059E-03, -8.305E-03, -8.613E-03, -9.134E-03, -1.023E-02, &
                   -1.365E-02, -1.320E-02, -2.922E-03, -1.125E-03/
  DATA B1_CSSSM_b2_19/1.007E-02, 1.008E-02, 1.012E-02, 1.019E-02, 1.024E-02, &
                   1.029E-02, 1.033E-02, 1.035E-02, 1.038E-02, 1.056E-02, &
                   1.083E-02, 1.109E-02, 1.140E-02, 1.194E-02, 1.270E-02, &
                   1.390E-02, 1.524E-02, 1.639E-02, 1.667E-02/
  DATA B2_CSSSM_b2_19/4.675E-03, 4.682E-03, 4.760E-03, 4.908E-03, 5.004E-03, &
                   5.102E-03, 5.168E-03, 5.207E-03, 5.256E-03, 5.474E-03, &
                   5.793E-03, 6.089E-03, 6.457E-03, 7.081E-03, 7.923E-03, &
                   9.127E-03, 1.041E-02, 1.147E-02, 1.173E-02/
  DATA C1_CSSSM_b2_19/7.571E-01, 7.571E-01, 7.570E-01, 7.568E-01, 7.565E-01, &
                   7.564E-01, 7.563E-01, 7.562E-01, 7.562E-01, 7.557E-01, &
                   7.552E-01, 7.545E-01, 7.539E-01, 7.527E-01, 7.509E-01, &
                   7.478E-01, 7.440E-01, 7.404E-01, 7.394E-01/
  DATA C2_CSSSM_b2_19/4.464E-02, 4.465E-02, 4.468E-02, 4.474E-02, 4.477E-02, &
                   4.480E-02, 4.482E-02, 4.484E-02, 4.486E-02, 4.448E-02, &
                   4.389E-02, 4.334E-02, 4.264E-02, 4.148E-02, 3.957E-02, &
                   3.588E-02, 3.149E-02, 2.751E-02, 2.650E-02/
!SUPER COARSE MODE
  DATA A1_SSSSM_b2_19/2.357E-01, 2.490E-01, 2.666E-01, 2.920E-01, 3.120E-01, &
                   3.128E-01, 3.169E-01, 3.272E-01, 3.498E-01, 3.960E-01, &
                   4.822E-01, 5.634E-01, 6.763E-01, 8.278E-01, 1.047E+00, &
                   1.340E+00, 1.927E+00, 1.648E+00, 1.031E+00/
  DATA A2_SSSSM_b2_19/1.219E-01, 1.337E-01, 1.633E-01, 1.929E-01, 2.057E-01, &
                   2.062E-01, 2.089E-01, 2.155E-01, 2.300E-01, 2.560E-01, &
                   2.908E-01, 3.199E-01, 3.530E-01, 3.965E-01, 4.475E-01, &
                   5.443E-01, 7.943E-01, 6.928E-01, 4.381E-01/
  DATA A3_SSSSM_b2_19/-2.387E-03, -2.599E-03, -3.092E-03, -3.599E-03, -3.832E-03, &
                   -3.842E-03, -3.890E-03, -4.012E-03, -4.276E-03, -4.763E-03, &
                   -5.455E-03, -6.051E-03, -6.763E-03, -7.708E-03, -8.887E-03, &
                   -1.091E-02, -1.585E-02, -1.373E-02, -8.665E-03/
  DATA B1_SSSSM_b2_19/1.260E-02, 1.211E-02, 1.126E-02, 1.056E-02, 1.038E-02, &
                   1.037E-02, 1.033E-02, 1.023E-02, 1.002E-02, 9.717E-03, &
                   9.613E-03, 9.652E-03, 9.983E-03, 1.047E-02, 1.168E-02, &
                   1.301E-02, 1.399E-02, 1.514E-02, 1.578E-02/
  DATA B2_SSSSM_b2_19/2.336E-03, 2.419E-03, 2.506E-03, 2.610E-03, 2.690E-03, &
                   2.694E-03, 2.711E-03, 2.752E-03, 2.844E-03, 3.043E-03, &
                   3.455E-03, 3.871E-03, 4.507E-03, 5.373E-03, 6.786E-03, &
                   8.238E-03, 9.208E-03, 1.032E-02, 1.091E-02/
  DATA C1_SSSSM_b2_19/7.832E-01, 7.787E-01, 7.721E-01, 7.670E-01, 7.657E-01, &
                   7.657E-01, 7.654E-01, 7.648E-01, 7.634E-01, 7.613E-01, &
                   7.596E-01, 7.585E-01, 7.574E-01, 7.560E-01, 7.533E-01, &
                   7.502E-01, 7.476E-01, 7.443E-01, 7.423E-01/
  DATA C2_SSSSM_b2_19/3.144E-02, 3.268E-02, 3.515E-02, 3.748E-02, 3.837E-02, &
                   3.840E-02, 3.860E-02, 3.906E-02, 4.006E-02, 4.173E-02, &
                   4.338E-02, 4.435E-02, 4.459E-02, 4.467E-02, 4.202E-02, &
                   3.864E-02, 3.559E-02, 3.183E-02, 2.964E-02/
!***************************************************************************

  spsol = 0
  spinsol = 0 
  spss = 0 

  DATA alpha_aers_2bands/  & 
       ! bc soluble
       7.675,7.675,7.675,7.675,7.675,7.675,    &
       7.675,7.675,10.433,11.984,13.767,15.567,& 
       4.720,4.720,4.720,4.720,4.720,4.720,    & 
       4.720,4.720,6.081,6.793,7.567,9.344,    & 
       ! pom soluble
       5.503,5.503,5.503,5.503,5.588,5.957,    & 
       6.404,7.340,8.545,10.319,13.595,20.398, & 
       1.402,1.402,1.402,1.402,1.431,1.562,    & 
       1.715,2.032,2.425,2.991,4.193,7.133,    & 
       ! sulfate    
       4.681,5.062,5.460,5.798,6.224,6.733,    & 
       7.556,8.613,10.687,12.265,16.32,21.692, & 
       1.107,1.239,1.381,1.490,1.635,1.8030,   &
       2.071,2.407,3.126,3.940,5.539,7.921,    &
                                ! sulfate coarse
       4.681,5.062,5.460,5.798,6.224,6.733,    & 
       7.556,8.613,10.687,12.265,16.32,21.692, & 
       1.107,1.239,1.381,1.490,1.635,1.8030,   &
       2.071,2.407,3.126,3.940,5.539,7.921,    &
                                ! seasalt Super Coarse Soluble (SS)
       0.5090,0.6554,0.7129,0.7767,0.8529,1.2728, &
       1.3820,1.5792,1.9173,2.2002,2.7173,4.1487, &
       0.5167,0.6613,0.7221,0.7868,0.8622,1.3027, &
       1.4227,1.6317,1.9887,2.2883,2.8356,4.3453, &
                                ! seasalt  Coarse Soluble (CS)
       0.5090,0.6554,0.7129,0.7767,0.8529,1.2728, &
       1.3820,1.5792,1.9173,2.2002,2.7173,4.1487, &
       0.5167,0.6613,0.7221,0.7868,0.8622,1.3027, &
       1.4227,1.6317,1.9887,2.2883,2.8356,4.3453, &
                                ! seasalt  Accumulation Soluble (AS)
       4.125, 4.674, 5.005, 5.434, 5.985, 10.006, &
       11.175,13.376,17.264,20.540,26.604, 42.349,&
       4.187, 3.939, 3.919, 3.937, 3.995,  5.078, &
       5.511, 6.434, 8.317,10.152,14.024, 26.537/

  DATA alpha_aeri_2bands/  & 
       ! dust insoluble
       0.7661,0.7123,&
       ! bc insoluble
       10.360,4.437, &
       ! pom insoluble
       3.741,0.606/

  DATA cg_aers_2bands/ &
       ! bc soluble
       .612, .612, .612, .612, .612, .612, &
       .612, .612, .702, .734, .760, .796, &
       .433, .433, .433, .433, .433, .433, &
       .433, .433, .534, .575, .613, .669, &
       ! pom soluble
       .663, .663, .663, .663, .666, .674, &
       .685, .702, .718, .737, .757, .777, &
       .544, .544, .544, .544, .547, .554, &
       .565, .583, .604, .631, .661, .698, &
       ! sulfate    
       .658, .669, .680, .688, .698, .707, &
       .719, .733, .752, .760, .773, .786, &
       .544, .555, .565, .573, .583, .593, &
       .610, .628, .655, .666, .692, .719, &
                                ! sulfate coarse
       .658, .669, .680, .688, .698, .707, &
       .719, .733, .752, .760, .773, .786, &
       .544, .555, .565, .573, .583, .593, &
       .610, .628, .655, .666, .692, .719, &
                                ! seasalt Super Coarse soluble (SS)
       .727, .747, .755, .761, .770, .788, &
       .792, .799, .805, .809, .815, .826, &
       .717, .738, .745, .752, .761, .779, &
       .781, .786, .793, .797, .803, .813, &
                                ! seasalt Coarse soluble (CS)
       .727, .747, .755, .761, .770, .788, &
       .792, .799, .805, .809, .815, .826, &
       .717, .738, .745, .752, .761, .779, &
       .781, .786, .793, .797, .803, .813, &
                                ! Sesalt Accumulation Soluble (AS)
       .727, .741, .748, .754, .761, .782, &
       .787, .792, .797, .799, .801, .799, &
       .606, .645, .658, .669, .681, .726, &
       .734, .746, .761, .770, .782, .798/

  DATA cg_aeri_2bands/ &
       ! dust insoluble
       .701, .670, &
       ! bc insoluble
       .471, .297, &
       ! pom insoluble
       .568, .365/

  DATA piz_aers_2bands/&
       ! bc soluble
       .445, .445, .445, .445, .445, .445, &
       .445, .445, .461, .480, .505, .528, &
       .362, .362, .362, .362, .362, .362, &
       .362, .362, .381, .405, .437, .483, &
       ! pom soluble
       .972, .972, .972, .972, .972, .974, &
       .976, .979, .982, .986, .989, .992, &
       .924, .924, .924, .924, .925, .927, &
       .932, .938, .945, .952, .961, .970, &
       ! sulfate
       1.000,1.000,1.000,1.000,1.000,1.000, &
       1.000,1.000,1.000,1.000,1.000,1.000, &
       .992, .988, .988, .987, .986, .985,  &
       .985, .985, .984, .984, .984, .984,  &
                                ! sulfate coarse
       1.000,1.000,1.000,1.000,1.000,1.000, &
       1.000,1.000,1.000,1.000,1.000,1.000, &
       .992, .988, .988, .987, .986, .985,  &
       .985, .985, .984, .984, .984, .984,  &
                                ! seasalt Super Coarse Soluble (SS)
       1.000,1.000,1.000,1.000,1.000,1.000, &
       1.000,1.000,1.000,1.000,1.000,1.000, &
       0.992,0.989,0.987,0.986,0.986,0.980, &
       0.980,0.978,0.976,0.976,0.974,0.971, &
                                ! seasalt Coarse soluble (CS)
       1.000,1.000,1.000,1.000,1.000,1.000, &
       1.000,1.000,1.000,1.000,1.000,1.000, &
       0.992,0.989,0.987,0.986,0.986,0.980, &
       0.980,0.978,0.976,0.976,0.974,0.971, &
                                ! seasalt Accumulation Soluble (AS)
       1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
       1.000, 1.000, 1.000, 1.000, 1.000, 1.000, &
       0.970, 0.975, 0.976, 0.977, 0.978, 0.982, &
       0.982, 0.983, 0.984, 0.984, 0.985, 0.985/

  DATA piz_aeri_2bands/ &
       ! dust insoluble
       .963, .987, &
       ! bc insoluble
       .395, .264, &
       ! pom insoluble
       .966, .859/

! Interpolation des coefficients optiques de 19 niveaux vers le nombre des niveaux du model
  IF (firstcall) THEN
     firstcall=.FALSE.
     
     IF (.NOT. ALLOCATED(A1_ASSSM_b1)) THEN
        ALLOCATE(A1_ASSSM_b1(klev),A2_ASSSM_b1(klev), A3_ASSSM_b1(klev),&
          B1_ASSSM_b1(klev), B2_ASSSM_b1(klev), C1_ASSSM_b1(klev), C2_ASSSM_b1(klev),&
          A1_CSSSM_b1(klev), A2_CSSSM_b1(klev), A3_CSSSM_b1(klev),&
          B1_CSSSM_b1(klev), B2_CSSSM_b1(klev), C1_CSSSM_b1(klev), C2_CSSSM_b1(klev),&
          A1_SSSSM_b1(klev), A2_SSSSM_b1(klev), A3_SSSSM_b1(klev),&
          B1_SSSSM_b1(klev), B2_SSSSM_b1(klev), C1_SSSSM_b1(klev), C2_SSSSM_b1(klev),&
          A1_ASSSM_b2(klev), A2_ASSSM_b2(klev), A3_ASSSM_b2(klev),&
          B1_ASSSM_b2(klev), B2_ASSSM_b2(klev), C1_ASSSM_b2(klev), C2_ASSSM_b2(klev),&
          A1_CSSSM_b2(klev), A2_CSSSM_b2(klev), A3_CSSSM_b2(klev),&
          B1_CSSSM_b2(klev), B2_CSSSM_b2(klev), C1_CSSSM_b2(klev), C2_CSSSM_b2(klev),&
          A1_SSSSM_b2(klev), A2_SSSSM_b2(klev), A3_SSSSM_b2(klev),&
          B1_SSSSM_b2(klev), B2_SSSSM_b2(klev), C1_SSSSM_b2(klev), C2_SSSSM_b2(klev), stat=ierr)
        IF (ierr /= 0) CALL abort_gcm('aeropt_2bands', 'pb in allocation 1',1)
     END IF
     
! bande 1
     CALL pres2lev(A1_ASSSM_b1_19, A1_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_ASSSM_b1_19, A2_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_ASSSM_b1_19, A3_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_ASSSM_b1_19, B1_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_ASSSM_b1_19, B2_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_ASSSM_b1_19, C1_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_ASSSM_b1_19, C2_ASSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

     CALL pres2lev(A1_CSSSM_b1_19, A1_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_CSSSM_b1_19, A2_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_CSSSM_b1_19, A3_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_CSSSM_b1_19, B1_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_CSSSM_b1_19, B2_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_CSSSM_b1_19, C1_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_CSSSM_b1_19, C2_CSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

     CALL pres2lev(A1_SSSSM_b1_19, A1_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_SSSSM_b1_19, A2_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_SSSSM_b1_19, A3_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_SSSSM_b1_19, B1_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_SSSSM_b1_19, B2_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_SSSSM_b1_19, C1_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_SSSSM_b1_19, C2_SSSSM_b1, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

! bande 2
     CALL pres2lev(A1_ASSSM_b2_19, A1_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_ASSSM_b2_19, A2_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_ASSSM_b2_19, A3_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_ASSSM_b2_19, B1_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_ASSSM_b2_19, B2_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_ASSSM_b2_19, C1_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_ASSSM_b2_19, C2_ASSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

     CALL pres2lev(A1_CSSSM_b2_19, A1_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_CSSSM_b2_19, A2_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_CSSSM_b2_19, A3_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_CSSSM_b2_19, B1_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_CSSSM_b2_19, B2_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_CSSSM_b2_19, C1_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_CSSSM_b2_19, C2_CSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

     CALL pres2lev(A1_SSSSM_b2_19, A1_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A2_SSSSM_b2_19, A2_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(A3_SSSSM_b2_19, A3_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B1_SSSSM_b2_19, B1_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(B2_SSSSM_b2_19, B2_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C1_SSSSM_b2_19, C1_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)
     CALL pres2lev(C2_SSSSM_b2_19, C2_SSSSM_b2, nb_level, klev, presnivs_19, presnivs, 1, 1, .FALSE.)

  END IF ! firstcall


  DO k=1, klev
    DO i=1, klon
      zrho=pplay(i,k)/t_seri(i,k)/RD                  ! kg/m3
!CDIR UNROLL=naero_spc
      mass_temp(i,k,:) = m_allaer(i,k,:) / zrho / 1.e+9
!RAF zrho
!CDIR UNROLL=naero_spc
      mass_temp_pi(i,k,:) = m_allaer_pi(i,k,:) / zrho / 1.e+9
      zdp1(i,k)=pdel(i,k)/(gravit*delt)      ! air mass auxiliary  variable --> zdp1 [kg/(m^2 *s)]
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
     aerosol_name(10)= id_CSSO4M
  ENDIF


  !
  ! loop over modes, use of precalculated nmd and corresponding sigma
  !    loop over wavelengths
  !    for each mass species in mode
  !      interpolate from Sext to retrieve Sext_at_gridpoint_per_species
  !      compute optical_thickness_at_gridpoint_per_species



!!CDIR ON_ADB(RH_tab)
!CDIR ON_ADB(fact_RH)
!CDIR SHORTLOOP
  DO n=1,nbre_RH-1
    fact_RH(n)=1./(RH_tab(n+1)-RH_tab(n))
  ENDDO
   
  DO k=1, KLEV
!!CDIR ON_ADB(RH_tab)
!CDIR ON_ADB(fact_RH)
    DO i=1, KLON
      rh(i,k)=MIN(RHcl(i,k)*100.,RH_MAX)
      RH_num(i,k) = INT( rh(i,k)/10. + 1.)
      IF (rh(i,k).GT.85.) RH_num(i,k)=10
      IF (rh(i,k).GT.90.) RH_num(i,k)=11
      
      DELTA(i,k)=(rh(i,k)-RH_tab(RH_num(i,k)))*fact_RH(RH_num(i,k))
    ENDDO
  ENDDO

  used_aer(:)=.FALSE.
    
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
     ELSEIF  (aerosol_name(m).EQ.id_CSSO4M) THEN
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

    id=aerosol_name(m)
    used_aer(id)=.TRUE.

     
    IF (soluble) THEN

      IF (spss.NE.0) THEN

         IF (spss.EQ.1) THEN !accumulation mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_ASSSM_b1)
!CDIR ON_ADB(A2_ASSSM_b1)
!CDIR ON_ADB(A3_ASSSM_b1)
!CDIR ON_ADB(B1_ASSSM_b1)
!CDIR ON_ADB(B2_ASSSM_b1)
!CDIR ON_ADB(C1_ASSSM_b1)
!CDIR ON_ADB(C2_ASSSM_b2)
!CDIR ON_ADB(A1_ASSSM_b2)
!CDIR ON_ADB(A2_ASSSM_b2)
!CDIR ON_ADB(A3_ASSSM_b2)
!CDIR ON_ADB(B1_ASSSM_b2)
!CDIR ON_ADB(B2_ASSSM_b2)
!CDIR ON_ADB(C1_ASSSM_b2)
!CDIR ON_ADB(C2_ASSSM_b2)
              DO i=1, KLON
                H=rh(i,k)/100
                tmp_var=mass_temp(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
                tmp_var_pi=mass_temp_pi(i,k,spsol)*1000.*zdp1(i,k)*delt*fac

                ! band 1
                tau_ae2b_int=A1_ASSSM_b1(k)+A2_ASSSM_b1(k)*H+A3_ASSSM_b1(k)/(H-1.05)
                piz_ae2b_int=1-B1_ASSSM_b1(k)-B2_ASSSM_b1(k)*H
                cg_ae2b_int=C1_ASSSM_b1(k)+C2_ASSSM_b1(k)*H

                tau_ae(i,k,id,1) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,1) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,1) = piz_ae2b_int
                cg_ae(i,k,id,1)= cg_ae2b_int
                
                !band 2
                tau_ae2b_int=A1_ASSSM_b2(k)+A2_ASSSM_b2(k)*H+A3_ASSSM_b2(k)/(H-1.05)
                piz_ae2b_int=1-B1_ASSSM_b2(k)-B2_ASSSM_b2(k)*H
                cg_ae2b_int=C1_ASSSM_b2(k)+C2_ASSSM_b2(k)*H

                tau_ae(i,k,id,2) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,2) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,2) = piz_ae2b_int
                cg_ae(i,k,id,2)= cg_ae2b_int

              ENDDO
            ENDDO
          ENDIF

          IF (spss.EQ.2) THEN !coarse mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_CSSSM_b1)
!CDIR ON_ADB(A2_CSSSM_b1)
!CDIR ON_ADB(A3_CSSSM_b1)
!CDIR ON_ADB(B1_CSSSM_b1)
!CDIR ON_ADB(B2_CSSSM_b1)
!CDIR ON_ADB(C1_CSSSM_b1)
!CDIR ON_ADB(C2_CSSSM_b2)
!CDIR ON_ADB(A1_CSSSM_b2)
!CDIR ON_ADB(A2_CSSSM_b2)
!CDIR ON_ADB(A3_CSSSM_b2)
!CDIR ON_ADB(B1_CSSSM_b2)
!CDIR ON_ADB(B2_CSSSM_b2)
!CDIR ON_ADB(C1_CSSSM_b2)
!CDIR ON_ADB(C2_CSSSM_b2)
              DO i=1, KLON
                H=rh(i,k)/100
                tmp_var=mass_temp(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
                tmp_var_pi=mass_temp_pi(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
                ! band 1
                tau_ae2b_int=A1_CSSSM_b1(k)+A2_CSSSM_b1(k)*H+A3_CSSSM_b1(k)/(H-1.05)
                piz_ae2b_int=1-B1_CSSSM_b1(k)-B2_CSSSM_b1(k)*H
                cg_ae2b_int=C1_CSSSM_b1(k)+C2_CSSSM_b1(k)*H

                tau_ae(i,k,id,1) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,1) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,1) = piz_ae2b_int
                cg_ae(i,k,id,1)= cg_ae2b_int

                ! band 2
                tau_ae2b_int=A1_CSSSM_b2(k)+A2_CSSSM_b2(k)*H+A3_CSSSM_b2(k)/(H-1.05)
                piz_ae2b_int=1-B1_CSSSM_b2(k)-B2_CSSSM_b2(k)*H
                cg_ae2b_int=C1_CSSSM_b2(k)+C2_CSSSM_b2(k)*H

                tau_ae(i,k,id,2) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,2) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,2) = piz_ae2b_int
                cg_ae(i,k,id,2)= cg_ae2b_int

             ENDDO
           ENDDO
         ENDIF

         IF (spss.EQ.3) THEN !super coarse mode
            DO k=1, KLEV
!CDIR ON_ADB(A1_SSSSM_b1)
!CDIR ON_ADB(A2_SSSSM_b1)
!CDIR ON_ADB(A3_SSSSM_b1)
!CDIR ON_ADB(B1_SSSSM_b1)
!CDIR ON_ADB(B2_SSSSM_b1)
!CDIR ON_ADB(C1_SSSSM_b1)
!CDIR ON_ADB(C2_SSSSM_b2)
!CDIR ON_ADB(A1_SSSSM_b2)
!CDIR ON_ADB(A2_SSSSM_b2)
!CDIR ON_ADB(A3_SSSSM_b2)
!CDIR ON_ADB(B1_SSSSM_b2)
!CDIR ON_ADB(B2_SSSSM_b2)
!CDIR ON_ADB(C1_SSSSM_b2)
!CDIR ON_ADB(C2_SSSSM_b2)
              DO i=1, KLON
                H=rh(i,k)/100
                tmp_var=mass_temp(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
                tmp_var_pi=mass_temp_pi(i,k,spsol)*1000.*zdp1(i,k)*delt*fac

                ! band 1 
                tau_ae2b_int=A1_SSSSM_b1(k)+A2_SSSSM_b1(k)*H+A3_SSSSM_b1(k)/(H-1.05)
                piz_ae2b_int=1-B1_SSSSM_b1(k)-B2_SSSSM_b1(k)*H
                cg_ae2b_int=C1_SSSSM_b1(k)+C2_SSSSM_b1(k)*H

                tau_ae(i,k,id,1) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,1) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,1) = piz_ae2b_int
                cg_ae(i,k,id,1)= cg_ae2b_int

                ! band 2
                tau_ae2b_int=A1_SSSSM_b2(k)+A2_SSSSM_b2(k)*H+A3_SSSSM_b2(k)/(H-1.05)
                piz_ae2b_int=1-B1_SSSSM_b2(k)-B2_SSSSM_b2(k)*H
                cg_ae2b_int=C1_SSSSM_b2(k)+C2_SSSSM_b2(k)*H

                tau_ae(i,k,id,2) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,2) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,2) = piz_ae2b_int
                cg_ae(i,k,id,2)= cg_ae2b_int

              ENDDO
            ENDDO
          ENDIF

        ELSE
                        
!CDIR ON_ADB(alpha_aers_2bands)
!CDIR ON_ADB(piz_aers_2bands)
!CDIR ON_ADB(cg_aers_2bands)
          DO k=1, KLEV
            DO i=1, KLON
              tmp_var=mass_temp(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
              tmp_var_pi=mass_temp_pi(i,k,spsol)*1000.*zdp1(i,k)*delt*fac
!CDIR UNROLL=nbands
              DO inu=1,nbands

                tau_ae2b_int= alpha_aers_2bands(RH_num(i,k),inu,spsol)+ & 
                              DELTA(i,k)* (alpha_aers_2bands(RH_num(i,k)+1,inu,spsol) - & 
                              alpha_aers_2bands(RH_num(i,k),inu,spsol))
                      
                piz_ae2b_int = piz_aers_2bands(RH_num(i,k),inu,spsol) + & 
                               DELTA(i,k)* (piz_aers_2bands(RH_num(i,k)+1,inu,spsol) - & 
                               piz_aers_2bands(RH_num(i,k),inu,spsol))
                      
                cg_ae2b_int = cg_aers_2bands(RH_num(i,k),inu,spsol) + & 
                              DELTA(i,k)* (cg_aers_2bands(RH_num(i,k)+1,inu,spsol) - & 
                              cg_aers_2bands(RH_num(i,k),inu,spsol))

                tau_ae(i,k,id,inu) = tmp_var*tau_ae2b_int
                tau_ae_pi(i,k,id,inu) =  tmp_var_pi* tau_ae2b_int
                piz_ae(i,k,id,inu) = piz_ae2b_int
                cg_ae(i,k,id,inu)= cg_ae2b_int
                         
              ENDDO
            ENDDO
          ENDDO
        
        ENDIF                     

      ELSE                                                    ! For all aerosol insoluble components

!CDIR ON_ADB(alpha_aers_2bands)
!CDIR ON_ADB(piz_aers_2bands)
!CDIR ON_ADB(cg_aers_2bands)
        DO k=1, KLEV
          DO i=1, KLON
            tmp_var=mass_temp(i,k,naero_soluble+ spinsol)*1000.*zdp1(i,k)*delt*fac
            tmp_var_pi=mass_temp_pi(i,k,naero_soluble+spinsol)*1000.*zdp1(i,k)*delt*fac
!CDIR UNROLL=nbands
            DO inu=1,nbands
              tau_ae2b_int = alpha_aeri_2bands(inu,spinsol)
              piz_ae2b_int = piz_aeri_2bands(inu,spinsol)
              cg_ae2b_int = cg_aeri_2bands(inu,spinsol) 

              tau_ae(i,k,id,inu) = tmp_var*tau_ae2b_int
              tau_ae_pi(i,k,id,inu) = tmp_var_pi*tau_ae2b_int
              piz_ae(i,k,id,inu) = piz_ae2b_int
              cg_ae(i,k,id,inu)= cg_ae2b_int
            ENDDO
          ENDDO
        ENDDO

      ENDIF ! soluble

    ENDDO  ! nb_aer  

  DO m=1,nb_aer   
    IF (.NOT. used_aer(m)) THEN
      tau_ae(:,:,:,:)=0.
      tau_ae_pi(:,:,:,:)=0.
      piz_ae(:,:,:,:)=0.
      cg_ae(:,:,:,:)=0.
    ENDIF
  ENDDO

  DO inu=1, nbands
    DO mrfspecies=1,naero_grp
      IF (mrfspecies .EQ. 2) THEN             ! = total aerosol AER	 
        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_ASSO4M,inu)+tau_ae(i,k,id_CSSO4M,inu)+ &
                                           tau_ae(i,k,id_ASBCM,inu)+tau_ae(i,k,id_AIBCM,inu)+   &						     
                                           tau_ae(i,k,id_ASPOMM,inu)+tau_ae(i,k,id_AIPOMM,inu)+ &	
                                           tau_ae(i,k,id_ASSSM,inu)+tau_ae(i,k,id_CSSSM,inu)+   &
                                           tau_ae(i,k,id_SSSSM,inu)+ tau_ae(i,k,id_CIDUSTM,inu)
	     tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
                 
             piz_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu)+ &
                                             tau_ae(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu)+ &
                                             tau_ae(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu)+ &
                                             tau_ae(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu)+ &
                                             tau_ae(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu)+ &
                                             tau_ae(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu)+ &	
                                             tau_ae(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu)+ &
                                             tau_ae(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu)+ &
                                             tau_ae(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu)+ &
                                             tau_ae(i,k,id_CIDUSTM,inu)*piz_ae(i,k,id_CIDUSTM,inu)) &
                                            /tau_allaer(i,k,mrfspecies,inu)
	     piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)

             cg_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu)*cg_ae(i,k,id_ASSO4M,inu)+ &
                      tau_ae(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu)*cg_ae(i,k,id_CSSO4M,inu)+ &
                      tau_ae(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu)*cg_ae(i,k,id_ASBCM,inu)+ &
                      tau_ae(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu)*cg_ae(i,k,id_AIBCM,inu)+ &
                      tau_ae(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu)*cg_ae(i,k,id_ASPOMM,inu)+ &
                      tau_ae(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu)*cg_ae(i,k,id_AIPOMM,inu)+ &	
                      tau_ae(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu)*cg_ae(i,k,id_ASSSM,inu)+ &
                      tau_ae(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu)*cg_ae(i,k,id_CSSSM,inu)+ &
                      tau_ae(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu)*cg_ae(i,k,id_SSSSM,inu)+ &
                      tau_ae(i,k,id_CIDUSTM,inu)*piz_ae(i,k,id_CIDUSTM,inu)*cg_ae(i,k,id_CIDUSTM,inu))/ &
                      (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO    
        ENDDO 

      ELSEIF (mrfspecies .EQ. 3) THEN             ! = natural aerosol NAT

        DO k=1, KLEV
          DO i=1, KLON
!RAF
	 	 tau_allaer(i,k,mrfspecies,inu)=tau_ae_pi(i,k,id_ASSO4M,inu)+ &
                      tau_ae_pi(i,k,id_CSSO4M,inu)+ &
                      tau_ae_pi(i,k,id_ASBCM,inu)+ &
                      tau_ae_pi(i,k,id_AIBCM,inu)+ &
                      tau_ae_pi(i,k,id_ASPOMM,inu)+ &
                      tau_ae_pi(i,k,id_AIPOMM,inu)+ &	
                      tau_ae_pi(i,k,id_ASSSM,inu)+ &
                      tau_ae_pi(i,k,id_CSSSM,inu)+ &
                      tau_ae_pi(i,k,id_SSSSM,inu)+ &
                      tau_ae_pi(i,k,id_CIDUSTM,inu)
	         tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)

	 	 piz_allaer(i,k,mrfspecies,inu)=(tau_ae_pi(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu)+ &
                      tau_ae_pi(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu)+ &
                      tau_ae_pi(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu)+ &
                      tau_ae_pi(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu)+ &
                      tau_ae_pi(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu)+ &
                      tau_ae_pi(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu)+ &	
                      tau_ae_pi(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu)+ &
                      tau_ae_pi(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu)+ &
                      tau_ae_pi(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu)+ &
                      tau_ae_pi(i,k,id_CIDUSTM,inu)*piz_ae(i,k,id_CIDUSTM,inu)) &
                      /tau_allaer(i,k,mrfspecies,inu)
	         piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)

	 	 cg_allaer(i,k,mrfspecies,inu)=(&
                      tau_ae_pi(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu)*cg_ae(i,k,id_ASSO4M,inu)+ &
                      tau_ae_pi(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu)*cg_ae(i,k,id_CSSO4M,inu)+ &
                      tau_ae_pi(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu)*cg_ae(i,k,id_ASBCM,inu)+ &
                      tau_ae_pi(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu)*cg_ae(i,k,id_AIBCM,inu)+ &
                      tau_ae_pi(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu)*cg_ae(i,k,id_ASPOMM,inu)+ &
                      tau_ae_pi(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu)*cg_ae(i,k,id_AIPOMM,inu)+ &
                      tau_ae_pi(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu)*cg_ae(i,k,id_ASSSM,inu)+ &
                      tau_ae_pi(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu)*cg_ae(i,k,id_CSSSM,inu)+ &
                      tau_ae_pi(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu)*cg_ae(i,k,id_SSSSM,inu)+ &
                      tau_ae_pi(i,k,id_CIDUSTM,inu)*piz_ae(i,k,id_CIDUSTM,inu)*&
                      cg_ae(i,k,id_CIDUSTM,inu))/ &
                      (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO
        ENDDO
                   
      ELSEIF (mrfspecies .EQ. 4) THEN             ! = BC
        DO k=1, KLEV
          DO i=1, KLON
	    tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_ASBCM,inu)+tau_ae(i,k,id_AIBCM,inu)
	    tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
	    piz_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu) &
                      +tau_ae(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu))/ &
                      tau_allaer(i,k,mrfspecies,inu)
	    piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)
            cg_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASBCM,inu)*piz_ae(i,k,id_ASBCM,inu) *cg_ae(i,k,id_ASBCM,inu)&
                      +tau_ae(i,k,id_AIBCM,inu)*piz_ae(i,k,id_AIBCM,inu)*cg_ae(i,k,id_AIBCM,inu))/ &
                      (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO
        ENDDO
              
      ELSEIF (mrfspecies .EQ. 5) THEN             ! = SO4

        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_ASSO4M,inu)+tau_ae(i,k,id_CSSO4M,inu)
	    tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
            piz_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu) &
                      +tau_ae(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu))/ &
                      tau_allaer(i,k,mrfspecies,inu)
	    piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)
            cg_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_CSSO4M,inu)*piz_ae(i,k,id_CSSO4M,inu) *cg_ae(i,k,id_CSSO4M,inu)&
                      +tau_ae(i,k,id_ASSO4M,inu)*piz_ae(i,k,id_ASSO4M,inu)*cg_ae(i,k,id_ASSO4M,inu))/ &
                      (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO
        ENDDO

      ELSEIF (mrfspecies .EQ. 6) THEN             ! = POM

        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_ASPOMM,inu)+tau_ae(i,k,id_AIPOMM,inu)
            tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
	    piz_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu) &
                      +tau_ae(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu))/ &
                      tau_allaer(i,k,mrfspecies,inu)
	    piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)
	    cg_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASPOMM,inu)*piz_ae(i,k,id_ASPOMM,inu) *cg_ae(i,k,id_ASPOMM,inu)&
                      +tau_ae(i,k,id_AIPOMM,inu)*piz_ae(i,k,id_AIPOMM,inu)*cg_ae(i,k,id_AIPOMM,inu))/ &
                      (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO
        ENDDO
              
      ELSEIF (mrfspecies .EQ. 7) THEN             ! = DUST

        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_CIDUSTM,inu)
	    tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
            piz_allaer(i,k,mrfspecies,inu)=piz_ae(i,k,id_CIDUSTM,inu)
	    cg_allaer(i,k,mrfspecies,inu)=cg_ae(i,k,id_CIDUSTM,inu)
          ENDDO
        ENDDO

      ELSEIF (mrfspecies .EQ. 8) THEN             ! = SS

        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=tau_ae(i,k,id_ASSSM,inu)+tau_ae(i,k,id_CSSSM,inu)+tau_ae(i,k,id_SSSSM,inu)
            tau_allaer(i,k,mrfspecies,inu)=MAX(tau_allaer(i,k,mrfspecies,inu),1e-5)
            piz_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu) &
                    +tau_ae(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu) &
                    +tau_ae(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu))/ &
                    tau_allaer(i,k,mrfspecies,inu)
            piz_allaer(i,k,mrfspecies,inu)=MAX(piz_allaer(i,k,mrfspecies,inu),0.1)
            cg_allaer(i,k,mrfspecies,inu)=(tau_ae(i,k,id_ASSSM,inu)*piz_ae(i,k,id_ASSSM,inu) *cg_ae(i,k,id_ASSSM,inu)&
                    +tau_ae(i,k,id_CSSSM,inu)*piz_ae(i,k,id_CSSSM,inu)*cg_ae(i,k,id_CSSSM,inu) &
                    +tau_ae(i,k,id_SSSSM,inu)*piz_ae(i,k,id_SSSSM,inu)*cg_ae(i,k,id_SSSSM,inu))/ &
                    (tau_allaer(i,k,mrfspecies,inu)*piz_allaer(i,k,mrfspecies,inu))
          ENDDO
        ENDDO
      
      ELSEIF (mrfspecies .EQ. 9) THEN             ! = NO3
      
        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=0.   ! preliminary
            piz_allaer(i,k,mrfspecies,inu)=0.
            cg_allaer(i,k,mrfspecies,inu)=0.
          ENDDO
        ENDDO
      
      ELSE

        DO k=1, KLEV
          DO i=1, KLON
            tau_allaer(i,k,mrfspecies,inu)=0.  
            piz_allaer(i,k,mrfspecies,inu)=0.
            cg_allaer(i,k,mrfspecies,inu)=0.
          ENDDO
        ENDDO
           
      ENDIF

    ENDDO
  ENDDO
   

  inu=1
  DO i=1, KLON
     absvisaer(i)=SUM((1-piz_allaer(i,:,:,inu))*tau_allaer(i,:,:,inu))
  END DO	

  DEALLOCATE(aerosol_name) 

END SUBROUTINE AEROPT_2BANDS
