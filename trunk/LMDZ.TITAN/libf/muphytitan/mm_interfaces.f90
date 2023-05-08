! Copyright 2013-2015,2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
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
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: mm_interfaces.f90
!! summary: Interfaces module for external functions
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_INTERFACES
  !! Interfaces to external functions.
  !! 
  !! The module contains the definitions of all "external" functions used by moments model which are 
  !! left to the developer's responsibility. 
  !!
  !! # Functions 
  !!
  !! - [[mm_interfaces(module):mm_alpha_s(interface)]] should compute the inter-moments relation coefficient 
  !!   as a function of the moment's order for the spherical mode.
  !! - [[mm_interfaces(module):mm_alpha_f(interface)]] should perform the same computations as 
  !!   [[mm_interfaces(module):mm_alpha_s(interface)]] but for the fractal mode.
  !! - [[mm_interfaces(module):mm_ps2s(interface)]] should compute the probability for particles of the 
  !!   spherical mode to remain in that mode during coagulation process.
  !! - [[mm_interfaces(module):mm_qmean(interface)]] should compute the mean eletric charge correction to be 
  !!   applied on each coagulation sub-kernels computed in mm_haze module.
  !! - [[mm_interfaces(module):mm_get_btk(interface)]] should compute the \(b_{k}^{T}\) coefficient of the
  !!   free-molecular regime.
  USE MM_MPREC
  IMPLICIT NONE

  PUBLIC

  INTERFACE 

    PURE FUNCTION mm_alpha_s(k) RESULT (res)
      !! Inter-moment relation for spherical aerosols size distribution law.
      !!
      !! The method computes the relation between the kth order moment and the 0th order moment:
      !! $$ \dfrac{M_{k}}{M_{0}} = r_{C}^{k} \times \alpha(k,a_{1},...a_{n}) $$
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: k !! Order of the moment.
      REAL(kind=mm_wp) :: res           !! Alpha value.
    END FUNCTION mm_alpha_s 

    PURE FUNCTION mm_alpha_f(k) RESULT (res)
      !! Inter-moment relation for fractal aerosols size distribution law.
      !!
      !! The method computes the relation between the kth order moment and the 0th order moment:
      !! $$ \dfrac{M_{k}}{M_{0}} = r_{C}^{k} \times \alpha(k,a_{1},...a_{n}) $$
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: k !! Order of the moment.
      REAL(kind=mm_wp) :: res           !! Alpha value.
    END FUNCTION mm_alpha_f

    FUNCTION mm_ps2s(rcs,k,flow,t,p) RESULT(res)
      !! Get the proportion of aerosols that remains in the spherical mode during SS coagulation.
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: rcs  !! Characteristic radius of the spherical size-distribution (m).
      REAL(kind=mm_wp), INTENT(in) :: t    !! Temperature (K).
      REAL(kind=mm_wp), INTENT(in) :: p    !! Pressure level (Pa).
      INTEGER, INTENT(in)          :: k    !! Order of the moment (0 or 3).
      INTEGER, INTENT(in)          :: flow !! Flow regime (0: continuous, 1: Free molecular).
      REAL(kind=mm_wp) :: res              !! Proportion of spherical particles that remains in the spherical mode.
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
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: rc1   !! Characteristic radius of the the first distribution (m).
      REAL(kind=mm_wp), INTENT(in) :: rc2   !! Characteristic radius of the the second distribution (m).
      INTEGER, INTENT(in)          :: order !! Moment's order (0 or 3).
      CHARACTER(len=2), INTENT(in) :: modes !! Interaction mode (a combination of [S,F]).
      REAL(kind=mm_wp), INTENT(in) :: temp  !! Temperature (K).
      REAL(kind=mm_wp), INTENT(in) :: pres  !! Pressure level (Pa).
      REAL(kind=mm_wp) :: res               !! Electric charging correction.
    END FUNCTION mm_qmean

    PURE FUNCTION mm_get_btk(t,k) RESULT(res)
      !! Get the \(b_{k}^{T}\) coefficient of the Free Molecular regime.
      !! 
      !! The method computes and returns the value of the pre-factor \(b_{k}^{T}\) used to
      !! approximate free-molecular regime coagulation kernels.
      !! @note 
      !! For more details about \(b_{k}^{T}\) coefficient, please read the 
      !! [scientific documentation](page/haze.html#free-molecular).
      !!
      !! @attention
      !! In its current version, the model only deals with fixed values of __k__ and __T__.
      !! __k__ can take the values (0,3) and, __T__, the values within [1,5].
      IMPORT mm_wp
      INTEGER, INTENT(in) :: t  !! Interaction identifier.
      INTEGER, INTENT(in) :: k  !! Moment order.
      REAL(kind=mm_wp) :: res   !! \(b_{k}^{T}\) value.
    END FUNCTION mm_get_btk


    ELEMENTAL FUNCTION mm_eta_g(t) RESULT (res)
      !! Get the air viscosity at a given temperature.
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
      REAL(kind=mm_wp) :: res           !! Air viscosity at given temperature (\(Pa.s^{-1}\)).
    END FUNCTION mm_eta_g

    ELEMENTAL FUNCTION mm_lambda_g(t,p) RESULT(res)
      !! Get the air mean free path at given temperature and pressure.
      IMPORT mm_wp
      REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
      REAL(kind=mm_wp), INTENT(in) :: p !! Pressure level (Pa).
      REAL(kind=mm_wp) :: res           !! Air mean free path (m).
    END FUNCTION mm_lambda_g

  END INTERFACE

END MODULE MM_INTERFACES

