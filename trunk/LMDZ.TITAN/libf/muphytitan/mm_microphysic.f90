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

!! file: mm_microphysic.f90
!! brief: Microphysic processes interface module.
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_MICROPHYSIC
  !! Microphysic processes interface module.
  USE MM_MPREC
  USE MM_GLOBALS
  USE MM_HAZE
  USE MM_CLOUDS
  USE MM_METHODS
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mm_muphys, mm_diagnostics, mm_get_radii

  !! Interface to main microphysics subroutine.
  !! 
  !! The interface computes calls either all the microphysics processes ([[mm_microphysic(module):muphys_all(function)]] 
  !! or only aerosols microphysics ([[mm_microphysic(module):muphys_nocld(function)]]) in a single call.
  INTERFACE mm_muphys
    MODULE PROCEDURE muphys_all, muphys_nocld
  END INTERFACE

  CONTAINS


    
  FUNCTION muphys_all(dm0a_s,dm3a_s,dm0a_f,dm3a_f,dm0n,dm3n,dm3i,dgazs) RESULT(ret)
    !! Compute the evolution of moments tracers through haze and clouds microphysics processes.
    !! 
    !! This method computes the evolution of all the microphysics tracers, given under the form of moments 
    !! (and molar fraction for cloud condensible species) during a time step.
    !! 
    !! The method requires that global variables of the model (i.e. variables declared in [[mm_globals(module)]] 
    !! module) are initialized/updated correctly (see [[mm_globals(module):mm_global_init(interface)]],
    !! [[mm_globals(module):mm_column_init(function)]],[[mm_globals(module):mm_aerosols_init(function)]] and
    !! [[mm_globals(module):mm_clouds_init(function)]]).
    !! 
    !! The tendencies returned by the method are defined on the vertical __layers__ of the model from the __GROUND__ to 
    !! the __TOP__ of the atmosphere. They should be added to the input variables used in the initialization methods 
    !! before the latter are called to initialize a new step.
    !! @note
    !! __dm3i__ and __dgazs__ are 2D-arrays with vertical __layers__ in the 1st dimension and the number of
    !! ice components in the 2nd. They share the same _species_ indexing. 
    !!
    !! It should be a 2D-array with the vertical layers in first dimension and the number of ice components in the second.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm0a_s
      !! Tendency of the 0th order moment of the spherical mode distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm3a_s
      !! Tendency of the 3rd order moment of the spherical mode distribution (\(m^{3}.m^{-2}\)). 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm0a_f
      !! Tendency of the 0th order moment of the fractal mode distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm3a_f
      !! Tendency of the 3rd order moment of the fractal mode distribution (\(m^{3}.m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm0n
      !! Tendency of the 0th order moment of the _CCN_ distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)   :: dm3n
      !! Tendency of the 3rd order moment of the _CCN_ distribution (\(m^{3}.m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:,:) :: dm3i
      !! Tendencies of the 3rd order moments of each ice components (\(m^{3}.m^{-2}\)). 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:,:) :: dgazs
      !! Tendencies of each condensible gaz species (\(mol.mol^{-1}\)). 
    LOGICAL :: ret
      !! .true. on success (i.e. model has been initialized at least once previously), .false. otherwise. 
    REAL(kind=mm_wp), DIMENSION(SIZE(dm0a_s)) :: zdm0a_f,zdm3a_f
    INTEGER                                   :: i
    ! Checks initialization
    ret = (mm_ini_col.AND.mm_ini_aer.AND.(.NOT.mm_w_clouds.OR.mm_ini_cld))
    IF (.NOT.ret) RETURN
    ! Calls haze microphysics (-> m-3)
    call mm_haze_microphysics(dm0a_s,dm3a_s,dm0a_f,dm3a_f)
    IF (mm_w_clouds) THEN
      ! Calls cloud microphysics (-> m-3)
      call mm_cloud_microphysics(zdm0a_f,zdm3a_f,dm0n,dm3n,dm3i,dgazs)
      ! add temporary aerosols tendencies (-> m-3)
      dm0a_f = dm0a_f + zdm0a_f  ; dm3a_f = dm3a_f + zdm3a_f
      ! reverse clouds tendencies (-> m-2)
      dm0n   = dm0n(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
      dm3n   = dm3n(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
      DO i=1,mm_nesp
        dm3i(:,i)  = dm3i(mm_nla:1:-1,i)  * mm_dzlev(mm_nla:1:-1)
        dgazs(:,i) = dgazs(mm_nla:1:-1,i) 
      ENDDO
    ELSE
      dm0n = 0._mm_wp ; dm3n = 0._mm_wp ; dm3i = 0._mm_wp ; dgazs = 0._mm_wp
    ENDIF
    ! multiply by altitude thickness and reverse vectors so they go from ground to top :)
    dm0a_s = dm0a_s(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    dm3a_s = dm3a_s(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    dm0a_f = dm0a_f(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    dm3a_f = dm3a_f(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    
    RETURN
  END FUNCTION muphys_all

  FUNCTION muphys_nocld(dm0a_s,dm3a_s,dm0a_f,dm3a_f) RESULT(ret)
    !! Compute the evolution of moments tracers through haze microphysics only.
    !! 
    !! This method is a __light__ version of [[mm_microphysic(module):muphys_all(function)]] where 
    !! only haze microphysics is computed and its tendencies returned.
    !!
    !! The method has the same requirements and remarks than [[mm_microphysic(module):muphys_all(function)]]. 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0a_s
      !! Tendency of the 0th order moment of the spherical mode distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3a_s
      !! Tendency of the 3rd order moment of the spherical mode distribution (\(m^{3}.m^{-2}\)). 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm0a_f
      !! Tendency of the 0th order moment of the fractal mode distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:) :: dm3a_f
      !! Tendency of the 3rd order moment of the fractal mode distribution (\(m^{3}.m^{-2}\)).
    LOGICAL :: ret
      !! .true. on succes (i.e. model has been initialized at least once previously), .false. otherwise. 
    ret = (mm_ini_col.AND.mm_ini_aer)
    IF (.NOT.ret) RETURN
    IF (mm_w_clouds.AND.mm_debug) THEN
      WRITE(*,'(a)') "WARNING: clouds microphysics enabled but will not be &
                     &computed... (wrong interface)"
    ENDIF
    ! Calls haze microphysics
    call mm_haze_microphysics(dm0a_s,dm3a_s,dm0a_f,dm3a_f)
    ! reverse vectors so they go from ground to top :)
    dm0a_s = dm0a_s(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1) 
    dm3a_s = dm3a_s(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    dm0a_f = dm0a_f(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1) 
    dm3a_f = dm3a_f(mm_nla:1:-1) * mm_dzlev(mm_nla:1:-1)
    RETURN
  END FUNCTION muphys_nocld

  SUBROUTINE mm_diagnostics(aer_prec,aer_s_flux,aer_f_flux,         &
                            ccn_prec,ccn_flux, ice_prec,ice_fluxes, &
                            gazs_sat)
    !! Get various diagnostic fields of the microphysics.
    !!
    !! The current diagnostics saved during microphysic computation are :
    !!
    !! - Mass fluxes (aerosols -both mode-, CCN and ices)
    !! - Precipitations (aerosols -total-, CCN and ices)
    !! - condensible gazes saturation ratio
    !!
    !! @note 
    !! Fluxes values are always negative as they account for sedimentation fluxes. They are set as 
    !! vector (for aerosols and CCN) or 2D-array (with the vertical structure in the first dimension
    !! and number of species in the second, for ice) and are ordered from __GROUND__ to __TOP__.
    !!
    !! @note
    !! Precipitations are always positive and defined in meters. For ice, it is set as a vector with 
    !! the precipitations of each cloud ice components.
    !!
    !! @note
    !! __ccnprec__, __iceprec__, __icefluxes__ and __gazsat__ are always set to 0 if clouds 
    !! microphysics is disabled (see [[mm_globals(module):mm_w_clouds(variable)]] documentation).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL                 :: aer_prec   !! Aerosols precipitations (both modes) (m).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:)   :: aer_s_flux !! Spherical aerosol mass flux (\(kg.m^{-2}.s^{-1}\)).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:)   :: aer_f_flux !! Fractal aerosol mass flux (\(kg.m^{-2}.s^{-1}\)).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL                 :: ccn_prec   !! CCN precipitations (m).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:)   :: ccn_flux   !! CCN mass flux (\(kg.m^{-2}.s^{-1}\)).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:)   :: ice_prec   !! Ice precipitations (m).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:,:) :: ice_fluxes !! Ice sedimentation fluxes (\(kg.m^{-2}.s^{-1}\)).
    REAL(kind=mm_wp), INTENT(out), OPTIONAL, DIMENSION(:,:) :: gazs_sat   !! Condensible gaz saturation ratios (--).

    IF (PRESENT(aer_prec))   aer_prec   = ABS(mm_aer_prec) 
    IF (PRESENT(aer_s_flux)) aer_s_flux = -mm_aer_s_flux(mm_nla:1:-1)
    IF (PRESENT(aer_f_flux)) aer_f_flux = -mm_aer_f_flux(mm_nla:1:-1)

    IF (mm_w_clouds) THEN
      IF (PRESENT(ccn_prec))   ccn_prec   = ABS(mm_ccn_prec)
      IF (PRESENT(ice_prec))   ice_prec   = ABS(mm_ice_prec)
      IF (PRESENT(ccn_flux))   ccn_flux   = -mm_ccn_flux(mm_nla:1:-1)
      IF (PRESENT(ice_fluxes)) ice_fluxes = -mm_ice_fluxes(mm_nla:1:-1,:)
      IF (PRESENT(gazs_sat))   gazs_sat   =  mm_gazs_sat(mm_nla:1:-1,:)
    ELSE 
      IF (PRESENT(ccn_prec))   ccn_prec   = 0._mm_wp
      IF (PRESENT(ice_prec))   ice_prec   = 0._mm_wp
      IF (PRESENT(ccn_flux))   ccn_flux   = 0._mm_wp
      IF (PRESENT(ice_fluxes)) ice_fluxes = 0._mm_wp
      IF (PRESENT(gazs_sat))   gazs_sat   = 0._mm_wp
    ENDIF
  END SUBROUTINE mm_diagnostics

  SUBROUTINE mm_get_radii(rcsph,rcfra,rccld)
    !! Get characteristic radii of microphysical tracers on the vertical grid.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:), OPTIONAL :: rcsph !! Spherical mode characteristic radius 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:), OPTIONAL :: rcfra !! Fractal mode characteristic radius 
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:), OPTIONAL :: rccld !! Cloud drops mean radius
    IF (mm_ini_aer) THEN
      IF (PRESENT(rcsph)) rcsph = mm_rcs(mm_nla:1:-1) 
      IF (PRESENT(rcfra)) rcfra = mm_rcf(mm_nla:1:-1) 
    ELSE
      IF (PRESENT(rcsph)) rcsph = 0._mm_wp 
      IF (PRESENT(rcfra)) rcfra = 0._mm_wp 
    ENDIF
    IF (PRESENT(rccld)) THEN
      IF (mm_w_clouds.AND.mm_ini_cld) THEN
        rccld = mm_drad(mm_nla:1:-1)
      ELSE
        rccld = 0._mm_wp
      ENDIF
    ENDIF
  END SUBROUTINE mm_get_radii

END MODULE MM_MICROPHYSIC

