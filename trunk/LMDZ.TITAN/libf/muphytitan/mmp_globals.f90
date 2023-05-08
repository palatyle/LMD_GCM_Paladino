MODULE MMP_GLOBALS
  !! Interface to YAMMS for the LMDZ GCM.
  USE MM_LIB
  USE DATASETS
  IMPLICIT NONE

  PUBLIC

  !> Alpha function parameters. 
  !!
  !! It stores the parameters of the inter-moments relation functions.
  !! 
  !! The inter-moments relation function is represented by the sum of exponential
  !! quadratic expressions:
  !!
  !! $$ 
  !! \displaystyle \alpha(k) = \sum_{i=1}^{n} \exp\left( a_{i}\times k^{2} + 
  !! b_{i}\times k^{2} +c_{i}\right) 
  !! $$
  TYPE, PUBLIC :: aprm
    !> Quadratic coefficients of the quadratic expressions.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: a
    !> Linear coefficients of the quadratic expressions.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: b
    !> Free term of the quadratic expressions.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: c
  END TYPE

  !> Size distribution parameters derived type.
  !!
  !! It stores the parameters of the size distribution law for Titan.
  !! 
  !! The size distribution law is represented by the minimization of a sum of
  !! power law functions:
  !!
  !! $$ 
  !! \displaystyle n\left(r\right) = \frac{A_{0}}{C+\sum_{i=1}^{n} A_{i}\times
  !!                                    \left(\frac{r}{r_{c}}\right)^{-b_{i}}}
  !! $$
  TYPE, PUBLIC :: nprm
    !> Scaling factor.
    REAL(kind=mm_wp)                            :: a0
    !> Characterisitic radius.
    REAL(kind=mm_wp)                            :: rc
    !> Additional constant to the sum of power law.
    REAL(kind=mm_wp)                            :: c
    !> Scaling factor of each power law.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: a
    !> Power of each power law.
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: b
  END TYPE 

  !> Inter-moment relation set of parameters for the spherical mode.
  TYPE(aprm), PUBLIC, SAVE :: mmp_asp
  !> Inter-moment relation set of parameters for the fractal mode.
  TYPE(aprm), PUBLIC, SAVE :: mmp_afp

  !> Size-distribution law parameters of the spherical mode.
  TYPE(nprm), PUBLIC, SAVE :: mmp_pns
  !> Size-distribution law parameters of the fractal mode.
  TYPE(nprm), PUBLIC, SAVE :: mmp_pnf

  !> Data set for @f$<Q>_{SF}^{M0}@f$.
  TYPE(dset2d), PUBLIC, SAVE, TARGET             :: mmp_qbsf0
  !> Extended values of [[mmp_gcm(module):mmp_qbsf0(variable)]] dataset.
  REAL(kind=mm_wp), PUBLIC, SAVE, DIMENSION(2,2) :: mmp_qbsf0_e
  !> Data set for @f$<Q>_{SF}^{M3}@f$.
  TYPE(dset2d), PUBLIC, SAVE, TARGET             :: mmp_qbsf3
  !> Extended values of [[mmp_gcm(module):mmp_qbsf3(variable)]] dataset.
  REAL(kind=mm_wp), PUBLIC, SAVE, DIMENSION(2,2) :: mmp_qbsf3_e
  !> Data set for @f$<Q>_{FF}^{M0}@f$.
  TYPE(dset2d), PUBLIC, SAVE, TARGET             :: mmp_qbff0
  !> Extended values of [[mmp_gcm(module):mmp_qbff0(variable)]] dataset.
  REAL(kind=mm_wp), PUBLIC, SAVE, DIMENSION(2,2) :: mmp_qbff0_e
  
  !> Data set for linear interpolation of transfert probability (M0/CO).
  TYPE(dset1d), PUBLIC, SAVE, TARGET :: mmp_pco0p
  !> Data set for linear interpolation of transfert probability (M3/CO).
  TYPE(dset1d), PUBLIC, SAVE, TARGET :: mmp_pco3p
  !> Data set for linear interpolation of transfert probability (M0/FM).
  TYPE(dset1d), PUBLIC, SAVE, TARGET :: mmp_pfm0p
  !> Data set for linear interpolation of transfert probability (M3/FM). 
  TYPE(dset1d), PUBLIC, SAVE, TARGET :: mmp_pfm3p

  !> \(b_{0}^{t}\) coefficients for Free-molecular regime kernel approximation.
  REAL(kind=mm_wp), PUBLIC, SAVE, DIMENSION(5) :: mmp_bt0 = (/1._mm_wp,1._mm_wp,1._mm_wp,1._mm_wp,1._mm_wp/)
  !> \(b_{3}^{t}\) coefficients for Free-molecular regime kernel approximation.
  REAL(kind=mm_wp), PUBLIC, SAVE, DIMENSION(5) :: mmp_bt3 = (/1._mm_wp,1._mm_wp,1._mm_wp,1._mm_wp,1._mm_wp/)

  !> Spherical probability transfert control flag.
  LOGICAL, SAVE :: mmp_w_ps2s = .true.
  !> Aerosol electric charge correction control flag.
  LOGICAL, SAVE :: mmp_w_qe   = .true.
  !> Optic look-up table file path.
  CHARACTER(len=:), ALLOCATABLE, SAVE :: mmp_optic_file

  CONTAINS 
    
  SUBROUTINE abort_program(err)
    !! Dump error message and abort the program.
    TYPE(error), INTENT(in) :: err !! Error object.
    WRITE(stderr,'(a)') "ERROR: "//TRIM(err%msg) 
    CALL EXIT(err%id)
  END SUBROUTINE abort_program

END MODULE MMP_GLOBALS
