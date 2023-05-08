! Copyright 2013-2015,2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.

!! file: mmp_moments.f90
!! summary: YAMMS/MP2M model external methods
!! author: J. Burgalat
!! date: 2013-2015,2017
!!
!! This file contains the definitions of all external methods that should be defined
!! for mp2m library. 
!! 
!! All the methods defined here satisify the interfaces defined in __m_interfaces__ module 
!! of YAMMS library.

PURE FUNCTION mm_alpha_s(k) RESULT (res)
  !! Inter-moment relation for spherical aerosols size distribution law.
  !! 
  !! The method computes the relation between the kth order moment and the 0th 
  !! order moment of the size-distribution law:
  !!
  !! $$ \dfrac{M_{k}}{M_{0}} = r_{C}^{k} \times \alpha(k,a_{1},...a_{n}) $$
  !!
  !! Here, alpha is computed as a sum of expenontial functions.
  USE MMP_GCM, ONLY : mmp_asp, mm_wp
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in) :: k !! k Order of the moment.
  REAL(kind=mm_wp) :: res           !! Alpha value.
  res = SUM(dexp(mmp_asp%a*k**2+mmp_asp%b*k+mmp_asp%c))
  RETURN
END FUNCTION mm_alpha_s 

PURE FUNCTION mm_alpha_f(k) RESULT (res)
  !! Inter-moment relation for fractal aerosols size distribution law.
  !!
  !! [[mm_alpha_f(function)]] performs the same computations as [[mm_alpha_s(function)]]
  !! using another set of parameters for the exponential functions.
  USE MMP_GCM, ONLY : mmp_afp, mm_wp
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in) :: k !! k Order of the moment.
  REAL(kind=mm_wp) :: res           !! Alpha value.
  res = SUM(dexp(mmp_afp%a*k**2+mmp_afp%b*k+mmp_afp%c))
  RETURN
END FUNCTION mm_alpha_f

FUNCTION mm_ps2s(rcs,k,flow,t,p) RESULT(res)
  !! Get the proportion of aerosols that remains in the spherical mode during SS coagulation.
  !!
  !! From __k__ and __flow__ values, the method selects one of the four probability datasets datasets
  !! in [[mmp_globals(module)]] module (for instance [[mmp_globals(module):mmp_pco0p(variable)]])
  !! and interpolates linearly probability for the given value of __rcs__, __t__ and __p__.
  !!
  !! @warning
  !! Here, the method assumes the datasets define the probability for __spherical__ particles to 
  !! be transferred in the __fractal__ mode, but returns the proportion of particles that remains
  !! in the mode (which is expected by mp2m model).
  !!
  !! @attention
  !! If value cannot be interpolated, the method aborts the program. Normally, it cannot happen 
  !! since we extrapolate the probability for characteristic radius value out of range.
  !!
  !! @attention
  !! Consequently, as the probability can only range from 0 to 1, it is wise to ensure that the 
  !! look-up table limits this range: To do so, one can just add two values at the start and end 
  !! of the table with probabilities respectively set to 0 and 1.
  USE LINTDSET
  USE LOCATORS
  USE MMP_GCM, ONLY : mmp_pco0p,mmp_pfm0p,mmp_pco3p,mmp_pfm3p,mmp_w_ps2s,mm_wp
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in) :: rcs
    !! Characteristic radius of the spherical size-distribution (m).
  INTEGER, INTENT(in)          :: k
    !! Order of the moment (0 or 3).
  INTEGER, INTENT(in)          :: flow
    !! Flow regime indicator (0: Continous, 1: Free-molecular).
  REAL(kind=mm_wp), INTENT(in) :: t
    !! Temperature (K).
  REAL(kind=mm_wp), INTENT(in) :: p
    !! Pressure level (Pa).
  REAL(kind=mm_wp) :: res
    !! Proportion of spherical particles that stay in the spherical mode during SS coagulation.
  TYPE(dset1d), POINTER :: pp
  res = 1._mm_wp
  IF (rcs <= 0.0_mm_wp .OR. .NOT.mmp_w_ps2s) RETURN 
  SELECT CASE(k+flow)
    CASE(0)      ; pp => mmp_pco0p ! 0 = 0 + 0 -> M0 / CO
    CASE(1)      ; pp => mmp_pfm0p ! 1 = 0 + 1 -> M0 / FM
    CASE(3)      ; pp => mmp_pco3p ! 3 = 3 + 0 -> M3 / CO
    CASE(4)      ; pp => mmp_pfm3p ! 4 = 3 + 1 -> M3 / FM
    CASE DEFAULT ; RETURN
  END SELECT
  IF (.NOT.hdcd_lint_dset(rcs,pp,locate_reg_ext,res)) THEN 
    WRITE(*,'(a)') "mm_moments:ps2s_sc: Cannot interpolate transfert probability"
    call EXIT(10)
  ELSE
    ! 05102017: do not care anymore for bad extrapolation: 
    ! Bound probability value between 0 and 1
    ! note: The input look-up table still must have strict monotic variation or
    !       awkward results can be produced.
    res = MAX(0.0_mm_wp,MIN(res,1.0_mm_wp))
    ! we have interpolated f = 1-p and we need p !
    res = 1._mm_wp - res
  ENDIF
END FUNCTION mm_ps2s

FUNCTION mm_qmean(rc1,rc2,order,modes,temp,pres) RESULT(res)
  !! Get the electric correction for coagulation kernel.
  !!
  !! The method computes the eletric charging correction to apply to the coagulation
  !! kernel as a function of the temperature, pressure and the characteristic radius of
  !! the mode involved in the coagulation.
  !! 
  !! Modes are referred by a two letters uppercase string with the combination of:
  !!
  !! - S : spherical mode
  !! - F : fractal mode
  !! 
  !! For example, SS means intra-modal coagulation for spherical particles.
  !!
  !! Here the electric charging correction is computed using linear interpolation from
  !! pre-tabulated values.
  USE LINTDSET
  USE LOCATORS
  USE MMP_GCM, ONLY : mmp_w_qe,mmp_qbsf0,mmp_qbsf3,mmp_qbff0, &
                      mmp_qbsf0_e,mmp_qbsf3_e,mmp_qbff0_e,mm_wp
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in)  :: rc1   !! Characteristic radius of the first mode (m).
  REAL(kind=mm_wp), INTENT(in)  :: rc2   !! Characteristic radius of the the second mode (m).
  INTEGER, INTENT(in)           :: order !! Moment's order (0 or 3 expected).
  CHARACTER(len=2), INTENT(in)  :: modes !! Interaction mode (a combination of [S,F]).
  REAL(kind=mm_wp), INTENT(in)  :: temp  !! Temperature (K).
  REAL(kind=mm_wp), INTENT(in)  :: pres  !! Pressure level (Pa). 
  REAL(kind=mm_wp) :: res                !! Electric charging correction.
  INTEGER       :: chx,np
  REAL(kind=mm_wp) :: vmin,vmax
  REAL(kind=mm_wp) :: r_tmp, t_tmp
  chx = 0 
  IF (.NOT.mmp_w_qe) THEN
    res = 1._mm_wp
    RETURN
  ENDIF

  IF (SCAN(modes(1:1),"sS") /= 0) chx = chx + 1
  IF (SCAN(modes(2:2),"sS") /= 0) chx = chx + 1
  IF (SCAN(modes(1:1),"fF") /= 0) chx = chx + 3
  IF (SCAN(modes(2:2),"fF") /= 0) chx = chx + 3
  chx = chx + order
  SELECT CASE(chx)
    CASE(2)      ! M0/SS
      res = 1._mm_wp 
    CASE(4)      ! M0/SF
      ! Fix max values of input parameters
      r_tmp = MAX(MIN(log(rc1),mmp_qbsf0_e(2,2)),mmp_qbsf0_e(2,1))
      t_tmp = MAX(MIN(temp,mmp_qbsf0_e(1,2)),mmp_qbsf0_e(1,1))
      ! Interpolates values
      IF (.NOT.hdcd_lint_dset(t_tmp,r_tmp,mmp_qbsf0,locate_reg,res)) THEN
        WRITE(*,'(a)') "mm_moments:mm_qmean: Cannot interpolate mean Qelec"
        call EXIT(10)
      ENDIF
    CASE(5)      ! M3/SS
      res = 1._mm_wp
    CASE(6)      ! M0/FF
      r_tmp = MAX(MIN(log(rc1),mmp_qbff0_e(2,2)),mmp_qbff0_e(2,1))
      t_tmp = MAX(MIN(temp,mmp_qbff0_e(1,2)),mmp_qbff0_e(1,1))
      IF (.NOT.hdcd_lint_dset(t_tmp,r_tmp,mmp_qbff0,locate_reg,res)) THEN
        WRITE(*,'(a)') "mm_moments:mm_qmean: Cannot interpolate mean Qelec"
        call EXIT(10)
      ENDIF
    CASE(7)      ! M3/SF
      r_tmp = MAX(MIN(log(rc1),mmp_qbsf3_e(2,2)),mmp_qbsf3_e(2,1))
      t_tmp = MAX(MIN(temp,mmp_qbsf3_e(1,2)),mmp_qbsf3_e(1,1))
      IF (.NOT.hdcd_lint_dset(t_tmp,r_tmp,mmp_qbsf3,locate_reg,res)) THEN
        WRITE(*,'(a)') "mm_moments:mm_qmean: Cannot interpolate mean Qelec"
        call EXIT(10)
      ENDIF
    CASE DEFAULT ! anything else :)
      res = 1._mm_wp
  END SELECT
  RETURN
END FUNCTION mm_qmean

PURE FUNCTION mm_get_btk(t,k) RESULT(res)
  !! Get the \(b_{k}^{T}\) coefficient of the Free Molecular regime.
  !! 
  !! The method get the value of the Free-molecular regime coagulation pre-factor \(b_{k}^{T}\).
  !! For more details about this coefficient, please read [Coagulation](page/haze.html#coagulation) 
  !! documentation page.
  !!
  !! @warning
  !! __k__ can only be one of the following value : 0 or 3. __t__ ranges only from 1 to 5.
  USE MMP_GCM, ONLY : mmp_bt0,mmp_bt3,mm_wp
  IMPLICIT NONE
  INTEGER, INTENT(in) :: t !! Type index of the \(b_{k}^{T}\) coefficient to get
  INTEGER, INTENT(in) :: k !! Moment Order of the \(b_{k}^{T}\) coefficient to get
  REAL(kind=mm_wp) :: res  !! \(b_{k}^{T}\) coefficient
  IF (.NOT.(k == 3 .OR. k == 0)) res = 0._mm_wp
  IF (t > 5 .OR. t < 1) res = 0._mm_wp
  IF (k == 0) THEN
    res = mmp_bt0(t)
  ELSE IF (k == 3) THEN
    res = mmp_bt3(t)
  ENDIF
  RETURN
END FUNCTION mm_get_btk

ELEMENTAL FUNCTION mm_eta_g(t) RESULT (res)
  !! Get the air viscosity at a given temperature.
  !!
  !! The function computes the air viscosity at temperature __t__ using Sutherland method.
  USE MMP_GCM, ONLY: mm_wp
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
  REAL(kind=mm_wp) :: res           !! Air viscosity at given temperature (\(Pa.s^{-1}\)).
  REAL (kind=mm_wp), PARAMETER :: eta0 = 1.75e-5_mm_wp, &
                                  tsut = 109._mm_wp,    &
                                  tref = 293._mm_wp
  res = eta0 *dsqrt(t/tref)*(1._mm_wp+tsut/tref)/(1._mm_wp+tsut/t)
  RETURN
END FUNCTION mm_eta_g

ELEMENTAL FUNCTION mm_lambda_g(t,p) RESULT(res)
  !! Get the air mean free path at given temperature and pressure.
  !!
  !! The method computes the air mean free path:
  !!
  !! $$ \lambda_{g} = \dfrac{k_{b}T}{4\sqrt{2}\pi r_{a}^2 P} $$
  !!
  !! Where \(\lambda_{g}\), is the air mean free path, \(k_{b}\) the Boltzmann constant, T the
  !! temperature, P the pressure level and \(r_{a}\) the radius of an _air molecule_.
  USE MMP_GCM, ONLY: mm_wp,mm_pi,mm_air_rad,mm_kboltz
  IMPLICIT NONE
  REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
  REAL(kind=mm_wp), INTENT(in) :: p !! Pressure level (Pa).
  REAL(kind=mm_wp) :: res           !! Air mean free path (m).
  res = mm_kboltz*t/(4._mm_wp*dsqrt(2._mm_wp)*mm_pi*(mm_air_rad**2)*p)
  RETURN
END FUNCTION mm_lambda_g
