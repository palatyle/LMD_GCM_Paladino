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

!! file: mm_methods.f90
!! summary: Model miscellaneous methods module.
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_METHODS
  !! Model miscellaneous methods module.
  !!
  !! The module contains miscellaneous methods used either in the haze and clouds parts of the model.
  !!
  !! All thermodynamic functions related to cloud microphysics (i.e. [[mm_methods(module):mm_lHeatX(interface)]],
  !! [[mm_methods(module):mm_sigX(interface)]] and [[mm_methods(module):mm_psatX(interface)]]) compute related equations 
  !! from \cite{reid1986}. A version of the book is freely available [here](http://f3.tiera.ru/3/Chemistry/References/Poling%20B.E.,%20Prausnitz%20J.M.,%20O'Connell%20J.P.%20The%20Properties%20of%20Gases%20and%20Liquids%20(5ed.,%20MGH,%202000)(ISBN%200070116822)(803s).pdf).
  !!
  !! The module defines the following functions/subroutines/interfaces:
  !!
  !! | name        | description 
  !! | :---------: | :-------------------------------------------------------------------------------------
  !! | mm_lheatx   | Compute latent heat released
  !! | mm_sigx     | Compute surface tension
  !! | mm_psatx    | Compute saturation vapor pressure
  !! | mm_qsatx    | Compute saturation mass mixing ratio
  !! | mm_fshape   | Compute shape factor
  !! | mm_lambda_g | Compute air mean free path
  !! | mm_eta_g    | Compute air viscosity
  !! | mm_get_kfm  | Compute the thermodynamic pre-factor of coagulation kernel in free-molecular regime
  !! | mm_get_kco  | Compute the thermodynamic pre-factor of coagulation kernel in continuous regime
  USE MM_MPREC
  USE MM_GLOBALS
  USE MM_INTERFACES
  IMPLICIT NONE

  PRIVATE 

  PUBLIC  :: mm_sigX, mm_LheatX, mm_psatX, mm_qsatx, mm_fshape, &
             mm_get_kco, mm_get_kfm, mm_eta_g, mm_lambda_g

  ! ---- INTERFACES

  !> Interface to surface tension computation functions.
  !!
  !! The method computes the surface tension of a given specie at given temperature(s).
  !!
  !! ```fortran
  !! FUNCTION mm_sigX(temp,xESP)
  !! ```
  !! 
  !! __xESP__ must always be given as a scalar. If __temp__ is given as a vector, then the method 
  !! computes the result for all the temperatures and returns a vector of same size than __temp__.
  INTERFACE mm_sigX
    MODULE PROCEDURE sigx_sc,sigx_ve
  END INTERFACE

  !> Interface to Latent heat computation functions.
  !! 
  !! The method computes the latent heat released of a given specie at given temperature(s).
  !!
  !! ```fortran
  !! FUNCTION mm_lheatX(temp,xESP)
  !! ```
  !!
  !! __xESP__ must always be given as a scalar. If __temp__ is given as a vector, then the method 
  !! computes the result for all the temperatures and returns a vector of same size than __temp__.
  INTERFACE mm_LheatX
    MODULE PROCEDURE lheatx_sc,lheatx_ve
  END INTERFACE

  !> Interface to saturation vapor pressure computation functions.
  !!
  !! ```fortran
  !! FUNCTION mm_psatX(temp,xESP)
  !! ```
  !! 
  !! The method computes the saturation vapor pressure of a given specie at given temperature(s).
  !!
  !! __xESP__ must always be given as a scalar. If __temp__ is given as a vector, then the method 
  !! computes the result for all the temperatures and returns a vector of same size than __temp__.
  INTERFACE mm_psatX
    MODULE PROCEDURE psatx_sc,psatx_ve
  END INTERFACE

  !! Interface to saturation mass mixing ratio computaiton functions.
  !!
  !! The method computes the mass mixing ratio at saturation of a given specie at given temperature(s) 
  !! and pressure level(s).
  !!
  !! ```fortran
  !! FUNCTION mm_qsatX(temp,pres,xESP) 
  !! ```
  !!
  !! __xESP__ must always be given as a scalar. If __temp__ and __pres__  are given as a vector (of same 
  !! size !), then the method computes the result for each couple of (temperature, pressure) and returns
  !! a vector of same size than __temp__.
  INTERFACE mm_qsatx
    MODULE PROCEDURE qsatx_sc,qsatx_ve
  END INTERFACE

  !> Interface to shape factor computation functions.
  !!
  !! The method computes the shape factor for the heterogeneous nucleation.
  !!
  !! ```fortran
  !! FUNCTION mm_fshape(m,x)
  !! ```
  !!
  !! Where __m__ is cosine of the contact angle and __x__ the curvature radius. __m__ must always be 
  !! given as a scalar. If __x__ is given as a vector, then the method compute the result for each 
  !! value of __x__ and and returns a vector of same size than __x__.
  INTERFACE mm_fshape
    MODULE PROCEDURE fshape_sc,fshape_ve
  END INTERFACE

  CONTAINS

  FUNCTION fshape_sc(cost,rap) RESULT(res)
    !! Get the shape factor of a ccn (scalar).
    !!
    !! The method computes the shape factor for the heterogeneous nucleation on a fractal particle. 
    !! Details about the shape factor can be found in \cite{prup1978}.
    REAL(kind=mm_wp), INTENT(in) :: cost !! Cosine of the contact angle.
    REAL(kind=mm_wp), INTENT(in) :: rap  !! Curvature radius (\(r_{particle}/r^{*}\)).
    REAL(kind=mm_wp) :: res              !! Shape factor value.
    REAL(kind=mm_wp) :: phi,a,b,c
    IF (rap > 3000._mm_wp) THEN
      res = ((2._mm_wp+cost)*(1._mm_wp-cost)**2)/4._mm_wp
    ELSE
      phi = dsqrt(1._mm_wp-2._mm_wp*cost*rap+rap**2)
      a = 1._mm_wp + ( (1._mm_wp-cost*rap)/phi )**3
      b = (rap**3) * (2._mm_wp-3._mm_wp*(rap-cost)/phi+((rap-cost)/phi)**3)
      c = 3._mm_wp * cost * (rap**2) * ((rap-cost)/phi-1._mm_wp)
      res = 0.5_mm_wp*(a+b+c)
    ENDIF
    RETURN
  END FUNCTION fshape_sc

  FUNCTION fshape_ve(cost,rap) RESULT(res)
    !! Get the shape factor of a ccn (vector).
    !!
    !! See [[mm_methods(module):fshape_sc(function)]].
    REAL(kind=mm_wp), INTENT(in)               :: cost !! Cosine of the contact angle.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: rap  !! Curvature radii (\(r_{particle}/r^{*}\)).
    REAL(kind=mm_wp), DIMENSION(SIZE(rap)) :: res      !! Shape factor value.
    REAL(kind=mm_wp), DIMENSION(SIZE(rap)) :: phi,a,b,c
    WHERE(rap > 3000._mm_wp)
      res = ((2._mm_wp+cost)*(1._mm_wp-cost)**2)/4._mm_wp
    ELSEWHERE 
      phi = dsqrt(1._mm_wp-2._mm_wp*cost*rap+rap**2)
      a = 1._mm_wp + ((1._mm_wp-cost*rap)/phi )**3
      b = (rap**3)*(2._mm_wp-3._mm_wp*(rap-cost)/phi+((rap-cost)/phi)**3)
      c = 3._mm_wp*cost*(rap**2)*((rap-cost)/phi-1._mm_wp)
      res = 0.5_mm_wp*(a+b+c)
    ENDWHERE
    RETURN
  END FUNCTION fshape_ve

  FUNCTION LHeatX_sc(temp,xESP) RESULT(res)
    !! Compute latent heat of a given specie at given temperature (scalar).
    !!
    !! The method computes the latent heat equation as given in \cite{reid1986} p. 220 (eq. 7-9.4).
    IMPLICIT NONE
    ! - DUMMY 
    REAL(kind=mm_wp), INTENT(in) :: temp !! temperature (K).
    TYPE(mm_esp), INTENT(in)     :: xESP !! Specie properties.
    REAL(kind=mm_wp) :: res              !! Latent heat of given specie at given temperature (\(J.kg^{-1}\)).
    REAL(kind=mm_wp) :: ftm
    ftm=MIN(1._mm_wp-temp/xESP%tc,1.e-3_mm_wp)
    res = mm_rgas*xESP%tc*(7.08_mm_wp*ftm**0.354_mm_wp+10.95_mm_wp*xESP%w*ftm**0.456_mm_wp)/xESP%masmol
  END FUNCTION LHeatX_sc
    
  FUNCTION LHeatX_ve(temp,xESP) RESULT(res)
    !! Compute latent heat of a given specie at given temperature (vector).
    !!
    !! See [[mm_methods(module):lheatx_sc(function)]].
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: temp !! temperatures (K).
    TYPE(mm_esp), INTENT(in)                   :: xESP !! Specie properties.
    REAL(kind=mm_wp), DIMENSION(SIZE(temp))    :: res  !! Latent heat of given specie at given temperatures (\(J.kg^{-1}\)).
    REAL(kind=mm_wp) :: ftm 
    INTEGER      :: i
    DO i=1,SIZE(temp)
      ftm=MIN(1._mm_wp-temp(i)/xESP%tc,1.e-3_mm_wp)
      res(i) = mm_rgas*xESP%tc*(7.08_mm_wp*ftm**0.354_mm_wp+10.95_mm_wp*xESP%w*ftm**0.456_mm_wp) / &
               xESP%masmol 
    ENDDO
  END FUNCTION LHeatX_ve

  FUNCTION sigX_sc(temp,xESP) RESULT(res)
    !! Get the surface tension between a given specie and the air (scalar).
    !! 
    !! The method computes the surface tension equation as given in \cite{reid1986} p. 637 (eq. 12-3.6).
    REAL(kind=mm_wp), INTENT(in) :: temp !! temperature (K).
    TYPE(mm_esp), INTENT(in)     :: xESP !! Specie properties.
    REAL(kind=mm_wp) :: res              !! Surface tension (\(N.m^{-1}\)).
    REAL(kind=mm_wp) :: tr,tbr,sig
    tr=MIN(temp/xESP%tc,0.99_mm_wp)
    tbr=xESP%tb/xESP%tc
    sig = 0.1196_mm_wp*(1._mm_wp+(tbr*dlog(xESP%pc/1.01325_mm_wp))/(1._mm_wp-tbr))-0.279_mm_wp
    sig = xESP%pc**(2._mm_wp/3._mm_wp)*xESP%tc**(1._mm_wp/3._mm_wp)*sig*(1._mm_wp-tr)**(11._mm_wp/9._mm_wp)
    res = sig*1e3_mm_wp ! dyn/cm2 -> N/m
  END FUNCTION sigX_sc
  
  FUNCTION sigX_ve(temp,xESP) RESULT(res)
    !! Get the surface tension between a given specie and the air (vector).
    !!
    !! See [[mm_methods(module):sigx_sc(function)]].
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: temp !! temperatures (K).
    TYPE(mm_esp), INTENT(in)                   :: xESP !! Specie properties.
    REAL(kind=mm_wp), DIMENSION(SIZE(temp)) :: res     !! Surface tensions (\(N.m^{-1}\)).
    INTEGER      :: i
    REAL(kind=mm_wp) :: tr,tbr,sig
    tbr = xESP%tb/xESP%tc
    sig = 0.1196_mm_wp*(1._mm_wp+(tbr*dlog(xESP%pc/1.01325_mm_wp))/(1._mm_wp-tbr))-0.279_mm_wp
    DO i=1,SIZE(temp)
      tr     = MIN(temp(i)/xESP%tc,0.99_mm_wp)
      sig    = xESP%pc**(2._mm_wp/3._mm_wp)*xESP%tc**(1._mm_wp/3._mm_wp)*sig*(1._mm_wp-tr)**(11._mm_wp/9._mm_wp)
      res(i) = sig*1e3_mm_wp ! dyn/cm2 -> N/m
    ENDDO
  END FUNCTION sigX_ve

  FUNCTION psatX_sc(temp,xESP) RESULT(res)
    !! Get saturation vapor pressure for a given specie at given temperature (scalar).
    !! 
    !! The method computes the saturation vapor pressure equation given in \cite{reid1986} p. 657 (eq. 1).
    !!
    !! @warning
    !! This subroutine accounts for a specific Titan feature: 
    !! If __xESP__ corresponds to \(CH_{4}\), the saturation vapor presure is multiplied by 0.85 
    !! to take into account its dissolution in \(N_{2}\).
    REAL(kind=mm_wp), INTENT(in) :: temp !! Temperature (K).
    TYPE(mm_esp), INTENT(in)     :: xESP !! Specie properties.
    REAL(kind=mm_wp) :: res              !! Saturation vapor pressure (Pa).
    REAL(kind=mm_wp) :: x,qsat
    x = 1._mm_wp-temp/xESP%tc
    IF (x > 0._mm_wp) THEN
      qsat = (1._mm_wp-x)**(-1) * &
      (xESP%a_sat*x + xESP%b_sat*x**1.5_mm_wp + xESP%c_sat*x**3 + xESP%d_sat*x**6)
    ELSE
      qsat = XESP%a_sat*x/abs(1._mm_wp-x)     ! approx for  t > tc
    ENDIF
    !  Special case : ch4 : x0.85 (dissolution in N2)
    IF (xESP%name == "ch4") THEN
      res = xESP%pc*dexp(qsat)*0.85_mm_wp
    ELSE
      res = xESP%pc*dexp(qsat)
    ENDIF
    ! now convert bar to Pa
    res = res * 1e5_mm_wp
  END FUNCTION psatX_sc

  FUNCTION psatX_ve(temp,xESP) RESULT(res)
    !! Get saturation vapor pressure for a given specie at given temperature (vector).
    !! 
    !! See [[mm_methods(module):psatX_sc(function)]].
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: temp !! Temperatures (K).
    TYPE(mm_esp), INTENT(in)                   :: xESP !! Specie properties.
    REAL(kind=mm_wp), DIMENSION(SIZE(temp)) :: res     !! Saturation vapor pressures (Pa).
    INTEGER      :: i
    REAL(kind=mm_wp) :: x,qsat
    DO i=1, SIZE(temp)
      x = 1._mm_wp-temp(i)/xESP%tc
      IF (x > 0._mm_wp) THEN
        qsat = (1._mm_wp-x)**(-1) * & 
        (xESP%a_sat*x + xESP%b_sat*x**1.5_mm_wp + xESP%c_sat*x**3 + xESP%d_sat*x**6)
      ELSE
        qsat = XESP%a_sat*x/abs(1._mm_wp-x)     ! approx for  t > tc
      ENDIF
      res(i) = xESP%pc*dexp(qsat)
      !  Peculiar case : ch4 : x0.85 (dissolution in N2)
      IF (xESP%name == "ch4") res(i) = res(i)* 0.85_mm_wp 
    ENDDO
    res = res * 1e5_mm_wp  ! bar -> Pa
  END FUNCTION psatX_ve

  FUNCTION qsatX_sc(temp,pres,xESP) RESULT(res)
    !! Get the mass mixing ratio of a given specie at saturation (scalar).
    REAL(kind=mm_wp), INTENT(in) :: temp !! Temperature (K).
    REAL(kind=mm_wp), INTENT(in) :: pres !! Pressure level (Pa).
    TYPE(mm_esp), INTENT(in)    :: xESP  !! Specie properties.
    REAL(kind=mm_wp) :: res              !! Mass mixing ratio of the specie.
    REAL(kind=mm_wp) :: x,psat
    psat = mm_psatX(temp,xESP)
    res = (psat / pres) * xESP%fmol2fmas
  END FUNCTION qsatX_sc

  FUNCTION qsatX_ve(temp,pres,xESP) RESULT(res)
    !! Get the mass mixing ratio of a given specie at saturation (vector).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: temp !! Temperatures (K).
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: pres !! Pressure levels (Pa).
    TYPE(mm_esp), INTENT(in)                   :: xESP !! Specie properties.
    REAL(kind=mm_wp), DIMENSION(SIZE(temp)) :: res     !! Mass mixing ratios of the specie.
    REAL(kind=mm_wp), DIMENSION(SIZE(temp)) :: psat
    psat = mm_psatX(temp,xESP)
    res = (psat / pres) * xESP%fmol2fmas
  END FUNCTION qsatX_ve

  ELEMENTAL FUNCTION mm_get_kco(t) RESULT(res)
    !! Get the Continuous regime thermodynamics pre-factor of the coagulation kernel.
    REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
    REAL(kind=mm_wp) :: res           !! Continuous regime thermodynamics pre-factor (\(m^{3}.s^{-1}\)).
    res = 2._mm_wp*mm_kboltz*t / (3._mm_wp*mm_eta_g(t))
    RETURN
  END FUNCTION mm_get_kco

  ELEMENTAL FUNCTION mm_get_kfm(t) RESULT(res)
    !! Get the Free Molecular regime thermodynamics pre-factor of the coagulation kernel.
    REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
    REAL(kind=mm_wp) :: res           !! Free Molecular regime thermodynamics pre-factor (\(m^{5/2}.s^{-1}\)).
    res = (6._mm_wp*mm_kboltz*t/mm_rhoaer)**(0.5_mm_wp)
    RETURN
  END FUNCTION mm_get_kfm

!  ELEMENTAL FUNCTION mm_eta_g(t) RESULT (res)
!    !! Get the air viscosity at a given temperature.
!    !!
!    !! The function computes the air viscosity at temperature __t__ using Sutherland method.
!    REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
!    REAL(kind=mm_wp) :: res           !! Air viscosity at given temperature (\(Pa.s^{-1}\)).
!    REAL (kind=mm_wp), PARAMETER :: eta0 = 1.75e-5_mm_wp, &
!                                    tsut = 109._mm_wp,    &
!                                    tref = 293._mm_wp
!    res = eta0 *dsqrt(t/tref)*(1._mm_wp+tsut/tref)/(1._mm_wp+tsut/t)
!    RETURN
!  END FUNCTION mm_eta_g
!
!  ELEMENTAL FUNCTION mm_lambda_g(t,p) RESULT(res)
!    !! Get the air mean free path at given temperature and pressure.
!    !!
!    !! The method computes the air mean free path:
!    !!
!    !! $$ \lambda_{g} = \dfrac{k_{b}T}{4\sqrt{2}\pi r_{a}^2 P} $$
!    !!
!    !! Where \(\lambda_{g}\), is the air mean free path, \(k_{b}\) the Boltzmann constant, T the
!    !! temperature, P the pressure level and \(r_{a}\) the radius of an _air molecule_.
!    REAL(kind=mm_wp), INTENT(in) :: t !! Temperature (K).
!    REAL(kind=mm_wp), INTENT(in) :: p !! Pressure level (Pa).
!    REAL(kind=mm_wp) :: res           !! Air mean free path (m).
!    res = mm_kboltz*t/(4._mm_wp*dsqrt(2._mm_wp)*mm_pi*(mm_air_rad**2)*p)
!    RETURN
!  END FUNCTION mm_lambda_g

END MODULE MM_METHODS
