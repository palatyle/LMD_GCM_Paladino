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

!! file: mm_globals.f90
!! summary: Parameters and global variables module.
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_GLOBALS
  !! Parameters and global variables module.
  !!  
  !! # Module overview
  !!
  !! The module defines all the parameters and global variables that are common
  !! to all other modules of the library.
  !! 
  !! It is separated in two parts :
  !!
  !! - Main parameters and global saved variables. Most of these variables should
  !!   be initialized once and should hold the same value during run-time. These
  !!   variables are completly public and initialized by [[mm_globals(module):mm_global_init(interface)]]
  !!   method.
  !! - The second part defines a set of vectors that defines the vertical structure of the atmosphere.
  !!   Each time a new atmospheric column has to be computed (either on a new timestep or on a new couple 
  !!   of longitude/latitude), these vectors should be intialized with new values by calling 
  !!   [[mm_globals(module):mm_column_init(function)]] method. 
  !!   This part is separated in two sets : 
  !!
  !!   - The atmospheric structure with temperature, pressure levels and altitude definitions. 
  !!   - The vertical profiles of tracers with the moments of the two aerosols modes (both \(M_{0}\) 
  !!     and \(M_{3}\) for a total of 4 vectors), the _clouds_ microphysics moments tracers (i.e. 
  !!     \(M_{0}\) and \(M_{3}\) for the ccn and \(M_{3}\) for the ice components).
  !!     Additionally, the module also stores intermediates variables of interest such as the 
  !!     characteristic radii of the aerosols modes, the mean drop radius and the drop density, 
  !!     the molar fraction of each condensible species (related to ice components) and some
  !!     scalar variables that holds arrays sizes.
  !!
  !! @note
  !! All the vectors that represent the vertical structure of the atmosphere (altitude, pressure and 
  !! temperature...) are oriented from the __TOP__ of the atmosphere to the __GROUND__.
  !!
  !! @note 
  !! The module also imports errors module from __FCCP__ library to get definitions of the error object 
  !! everywhere in the library ([[mm_globals(module)]] is always imported, except in [[mm_mprec(module)]]).
  !!
  !! # Global variables 
  !!
  !! [[mm_globals(module)]] module contains the declaration of all global/common variable that are shared
  !! by all other modules of the model. Except for few physical constant which are declared as parameters,
  !! these variables are onlu SAVEd. They are initialized by [[mm_globals(module):mm_global_init(interface)]]
  !! methods.
  !! the following sections list all the global variables by category.
  !!
  !! ## Control flags 
  !! 
  !! | Name               | Description
  !! | :----------------- | :-----------------
  !! | mm_log             | Enable log mode (verbose)
  !! | mm_w_haze_prod     | Enable/Disable haze production
  !! | mm_w_haze_sed      | Enable/Disable haze sedimentation
  !! | mm_w_haze_coag     | Enable/Disable haze coagulation
  !! | mm_w_clouds        | Enable/Disable clouds microphysics
  !! | mm_w_clouds_sed    | Enable/Disable clouds microphysics sedimentation
  !! | mm_w_clouds_nucond | Enable/Disable clouds microphysics nucleation/condensation
  !! | mm_wsed_m0         | Force all aerosols moments to fall at M0 settling velocity 
  !! | mm_wsed_m3         | Force all aerosols moments to fall at M3 settling velocity
  !! | mm_no_fiadero_w    | Enable/Disable __Fiadero__ correction
  !!
  !! ### Related free parameters:
  !!
  !! | Name            | Description
  !! | :-------------- | :-----------------
  !! | mm_fiadero_min  | Minimum ratio for __Fiadero__'s correction 
  !! | mm_fiadero_max  | Maximum ratio for __Fiadero__'s correction
  !! | mm_coag_choice  | Coagulation interaction activation flag. It should be a combination of [[mm_globals(module):mm_coag_no(variable)]], [[mm_globals(module):mm_coag_ss(variable)]], [[mm_globals(module):mm_coag_sf(variable)]] and [[mm_globals(module):mm_coag_ff(variable)]]. 
  !!
  !! ## Physical constants 
  !!
  !! | Name      | Description
  !! | :-------- | :-----------------
  !! | mm_pi     | Pi number
  !! | mm_navo   | Avogadro number
  !! | mm_kboltz | Boltzmann constant (\(J.K^{-1}\))
  !! | mm_rgas   | Perfect gas constant (\(J.mol^{-1}.K^{-1}\))
  !! | mm_fdes   | Desorption energy (\(J\)) (nucleation)
  !! | mm_fdif   | Surface diffusion energy (\(J\)) (nucleation)
  !! | mm_fnus   | Jump frequency (\(s^{-1}\)) (nucleation)
  !! | mm_akn    | Approximated slip-flow correction coefficient (
  !!
  !! ## Free parameters
  !!
  !! | Name        | Description
  !! | :---------- | :-----------------
  !! | mm_rhoaer   | Aerosol density (in \(kg.m^{-3}\))
  !! | mm_df       | Fractal dimension
  !! | mm_rm       | Monomer radius (in m)
  !! | mm_p_prod   | Spherical aerosols production pressure level (Pa)
  !! | mm_p_rcprod | Spherical aerosols equivalent radius production (m)
  !! | mm_tx_prod  | Production rate of spherical aerosols (\(kg.m^{-2}.s^{-1}\))
  !! | mm_d_prod   | Time-dependent sine wve pre-factor.
  !! | mm_w_prod   | Angular frequency of the time-dependent production rate.
  !! | mm_ne       | Electric charging of aerosols (\(e^{-}.m^{-1}\)) (unused)
  !! | mm_rb2ra    | Bulk to apparent radius conversion pre-factor (\(m^X\)) 
  !! | mm_rpla     | Planet radius (m)
  !! | mm_g0       | Planet acceleration due to gravity constant (ground) (\(m.s^{-2}\))
  !! | mm_air_rad  | Air molecules mean radius (m)
  !! | mm_air_mmol | Air molecules molar mass (\(kg.mol^{-1}\))
  !! | mm_dt       | Microphysic time step (s)
  USE MM_MPREC
  USE MM_INTERFACES
  ! from swift
  USE CFGPARSE
  USE STRING_OP
  USE ERRORS
  IMPLICIT NONE

  PUBLIC

  PRIVATE :: cldprop_sc,cldprop_ve,read_esp,check_r1,check_i1,check_l1,check_s1

  ! Protected variables
  ! the following variables are read-only outside this module.
  ! One must call the afferent subroutine to update them.
    
  ! initialization control flags (cannot be updated)
  PROTECTED :: mm_ini,mm_ini_col,mm_ini_aer,mm_ini_cld
  ! model parameters (mm_global_init)
  PROTECTED :: mm_dt,mm_rhoaer,mm_df,mm_rm,mm_p_prod,mm_rc_prod,mm_tx_prod,mm_rpla,mm_g0,mm_rb2ra
  ! atmospheric vertical structure (mm_column_init)
  PROTECTED :: mm_nla,mm_nle,mm_zlay,mm_zlev,mm_play,mm_plev,mm_temp,mm_rhoair,mm_btemp,mm_dzlev,mm_dzlay
  ! Condensible species parameters (mm_global_init)
  PROTECTED :: mm_nesp,mm_spcname,mm_xESPS
  ! Moments parameters (mm_aerosols_init / mm_clouds_init)
  PROTECTED :: mm_m0aer_s, mm_m3aer_s, mm_m0aer_f, mm_m3aer_f, mm_m0ccn, mm_m3ccn, mm_m3ice
  ! Moments parameters (derived, are updated with moments parameters)
  PROTECTED :: mm_rcs, mm_rcf, mm_drad, mm_drho

  LOGICAL, SAVE :: mm_debug = .true.  !! Enable QnD debug mode (can be used for devel).
  LOGICAL, SAVE :: mm_log = .false.   !! Enable log mode (for configuration only).

  LOGICAL, SAVE :: mm_w_haze_prod = .true. !! Enable/Disable haze production.
  LOGICAL, SAVE :: mm_w_haze_sed = .true.  !! Enable/Disable haze sedimentation.
  LOGICAL, SAVE :: mm_w_haze_coag = .true. !! Activate haze coagulation.

  LOGICAL, SAVE :: mm_wsed_m0 = .false. !! Force all aerosols moments to fall at M0 settling velocity.
  LOGICAL, SAVE :: mm_wsed_m3 = .false. !! Force all aerosols moments to fall at M3 settling velocity.

  LOGICAL, SAVE :: mm_var_prod = .false. !! Time variation of production rate control flag.

  LOGICAL, SAVE :: mm_use_effg = .true. !! Enable/Disable effective G for computations.

  !> Enable/Disable __Fiadero__'s correction.
  !!
  !! This flag enables/disables the __Fiadero__ correction alogrithm for fractal mode settling velocity 
  !! computation. 
  !!
  !! @bug
  !! Currently, the Fiadero correction creates instatibilities on the vertical structure. It seems to be 
  !! related to the coupling between the two moments. In order to reduce the instabilities, settling
  !! velocity of moments are forced to be the same, see [[mm_globals(module):mm_wsed_m0(variable)]] and
  !! [[mm_globals(module):mm_wsed_m3(variable)]]).
  LOGICAL, SAVE          :: mm_no_fiadero_w = .false. 

  !> Minimum ratio for __Fiadero__ correction.
  !!
  !! When [[mm_globals(module):mm_no_fiadero_w(variable)]] is disabled, this variable defines the minimum 
  !! value of the moment's ratio between two adjacents vertical cells to be used within the correction.
  REAL(kind=mm_wp), SAVE :: mm_fiadero_min  = 0.1_mm_wp

  !> Maximum ratio for __Fiadero__ correction.
  !!
  !! When [[mm_globals(module):mm_no_fiadero_w(variable)]] is disabled, this variable defines the maximum 
  !! value of the moment's ratio between two adjacents vertical cells to be used within the correction.
  REAL(kind=mm_wp), SAVE :: mm_fiadero_max  = 10._mm_wp

  LOGICAL, SAVE :: mm_w_clouds = .true.       !! Enable/Disable clouds microphysics.
  LOGICAL, SAVE :: mm_w_cloud_sed = .true.    !! Enable/Disable cloud sedimentation.
  LOGICAL, SAVE :: mm_w_cloud_nucond = .true. !! Activate cloud nucleation/condensation.

  INTEGER, PARAMETER :: mm_coag_no = 0 !! no mode interaction for coagulation (i.e. no coagulation at all).
  INTEGER, PARAMETER :: mm_coag_ss = 1 !! SS mode interaction for coagulation.
  INTEGER, PARAMETER :: mm_coag_sf = 2 !! SF mode interaction for coagulation.
  INTEGER, PARAMETER :: mm_coag_ff = 4 !! FF mode interaction for coagulation.
  !> Default interactions to activate (all by default).
  INTEGER, SAVE      :: mm_coag_choice = mm_coag_ss+mm_coag_sf+mm_coag_ff 

  !> Pi number.
  REAL(kind=mm_wp), PARAMETER :: mm_pi = 4._mm_wp*atan(1._mm_wp)
  !> Avogadro number.
  REAL(kind=mm_wp), PARAMETER :: mm_navo = 6.0221367e23_mm_wp 
  !> Boltzmann constant (\(J.K^{-1}\)).
  REAL(kind=mm_wp), PARAMETER :: mm_kboltz = 1.3806488e-23_mm_wp
  !> Perfect gas constant (\(J.mol^{-1}.K^{-1}\)).
  REAL(kind=mm_wp), PARAMETER :: mm_rgas = mm_kboltz * mm_navo
  !> Desorption energy (\(J\)) (nucleation).
  REAL(kind=mm_wp), PARAMETER :: mm_fdes = 0.288e-19_mm_wp
  !> Surface diffusion energy (\(J\)) (nucleation).
  REAL(kind=mm_wp), PARAMETER :: mm_fdif = 0.288e-20_mm_wp
  !> Jump frequency (\(s^{-1}\)) (nucleation).
  REAL(kind=mm_wp), PARAMETER :: mm_nus = 1.e+13_mm_wp
  !> Approximated slip-flow correction coefficient.
  REAL(kind=mm_wp), PARAMETER :: mm_akn = 1.591_mm_wp

  !> Aerosols density (\(kg.m^{-3}\)).
  REAL(kind=mm_wp), SAVE :: mm_rhoaer = 1.e3_mm_wp

  !> Fractal dimension of fractal aerosols.
  REAL(kind=mm_wp), SAVE :: mm_df = 3._mm_wp

  !> Monomer radius (m).
  REAL(kind=mm_wp), SAVE :: mm_rm = 6.66e-8_mm_wp

  !> Spherical aerosols production pressure level (Pa).
  REAL(kind=mm_wp), SAVE :: mm_p_prod = 1._mm_wp

  !> Spherical aerosols equivalent radius production (m)
  REAL(kind=mm_wp), SAVE :: mm_rc_prod = 1.3101721857598102e-9_mm_wp

  !> Production rate of spherical aerosols (\(kg.m^{-2}.s^{-1}\)).
  REAL(kind=mm_wp), SAVE :: mm_tx_prod = 3.5e-13_mm_wp

  !> Aerosol production delta if time variations is enabled (fraction).
  REAL(kind=mm_wp), SAVE :: mm_d_prod  = 0.25_mm_wp

  !> Aerosol production variations angular frequency if time variations is enabled (\(rad.s^{-1}\)).
  REAL(kind=mm_wp), SAVE :: mm_w_prod  = 2.*mm_pi / (86400.*16.)


  !> Electric charging of aerosols (\(e^{-}.m^{-1}\)).
  REAL(kind=mm_wp), SAVE :: mm_ne = -15.e6_mm_wp

  !> Bulk to apparent radius conversion pre-factor (\(m^{X}\)).
  !! 
  !! It is initialized using [[mm_globals(module):mm_rm(variable)]] in 
  !! [[mm_globals(module):mm_global_init(interface)]] from the following equation:
  !!
  !! $$ r_{a} = r_{b}^{3/D_{f}}\times r_{m}^{\frac{D_{f}-3}{D_{f}}} $$
  !!
  !! Where \(r_{a}\) is the apparent radius, \(r_{b}\) the bulk radius and 
  !! \(rb2ra = r_{m}^{\frac{D_{f}-3}{D_{f}}}\) is the returned pre-factor
  REAL(kind=mm_wp), SAVE :: mm_rb2ra = 1._mm_wp 

  !> Characteristic radius threshold.
  REAL(kind=mm_wp), SAVE :: mm_rc_min = 1.e-200_mm_wp

  !> Name of condensible species.
  CHARACTER(len=30), DIMENSION(:), ALLOCATABLE, SAVE :: mm_spcname

  TYPE, PUBLIC :: mm_esp
    !! Cloud related chemical specie properties.
    !!
    !! This derived type is used in thermodynamic methods related to cloud microphysics.
    !! Most of its fields represent parameters of equations from \cite{reid1986}.
    CHARACTER(LEN=10) :: name      !! Specie name.
    REAL(kind=mm_wp)  :: mas       !! Molecular weight (kg).
    REAL(kind=mm_wp)  :: vol       !! Molecular volume (\(m^{3}\)).
    REAL(kind=mm_wp)  :: ray       !! Molecular radius (m).
    REAL(kind=mm_wp)  :: masmol    !! Molar mass (\(kg.mol^{-1}\)).
    REAL(kind=mm_wp)  :: rho       !! density (liquid) (\(kg.m^{-3}\)).
    REAL(kind=mm_wp)  :: tc        !! Critical temperature (K).
    REAL(kind=mm_wp)  :: pc        !! Critical pressure (Bar).
    REAL(kind=mm_wp)  :: tb        !! Boiling point temperature (K).
    REAL(kind=mm_wp)  :: w         !! Acentric factor (--).
    REAL(kind=mm_wp)  :: a_sat     !! Saturation equation A coefficient.
    REAL(kind=mm_wp)  :: b_sat     !! Saturation equation B coefficient.
    REAL(kind=mm_wp)  :: c_sat     !! saturation equation C coefficient.
    REAL(kind=mm_wp)  :: d_sat     !! Saturation equation D coefficient.
    REAL(kind=mm_wp)  :: mteta     !! Wettability.
    REAL(kind=mm_wp)  :: tx_prod   !! Production rate.
    REAL(kind=mm_wp)  :: fmol2fmas !! molar fraction to mass fraction coefficient.
    ! = masmol(X)/masmol(AIR)
  END TYPE mm_esp

  !> Planet radius (m).
  REAL(kind=mm_wp), SAVE                        :: mm_rpla     = 2575000._mm_wp
  !> Planet acceleration due to gravity constant (ground) (\(m.s^{-2}\)).
  REAL(kind=mm_wp), SAVE                        :: mm_g0       = 1.35_mm_wp
  !> Air molecules mean radius (m).
  REAL(kind=mm_wp), SAVE                        :: mm_air_rad  = 1.75e-10_mm_wp
  !> Air molecules molar mass (\(kg.mol^{-1}\)).
  REAL(kind=mm_wp), SAVE                        :: mm_air_mmol = 28e-3_mm_wp
  !> Microphysic time step (s).
  REAL(kind=mm_wp), SAVE                        :: mm_dt       = 5529.6_mm_wp
  !> Model current time tracer (s).
  REAL(kind=mm_wp), SAVE                        :: mm_ct       = 0.0
  !> Total number of clouds condensible species.
  INTEGER, SAVE                                 :: mm_nesp     = -1
  !> Clouds chemical species properties.
  TYPE(mm_esp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_xESPS

  !------------------------
  ! Vertical structure part
  !------------------------

  !> Number of vertical layers.
  INTEGER, SAVE :: mm_nla = -1
  !> Number of vertical levels.
  INTEGER, SAVE :: mm_nle = -1

  !> Altitude layers (m).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_zlay
  !> Altitude levels (m).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_zlev
  !> Pressure layers (Pa).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_play
  !> Pressure levels (Pa).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_plev
  !>  Temperature vertical profile (K).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_temp
  !>  Air density vertical profile (\(kg.m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_rhoair 
  !> Temperature vertical profil at interfaces (K).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_btemp

  !> Atmospheric levels thickness (m).
  !! 
  !! Atmospheric thickness between two adjacent levels (\(m\)) from the 
  !! __TOP__ to the __GROUND__.
  !! @note __mm_dzlev__ is defined on the total number of layers and actually
  !! corresponds to the thickness of a given layer.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_dzlev

  !> Atmospheric layers "thickness" (m).
  !! 
  !! Atmospheric thickness between the center of two adjacent layers (\(m\))
  !! from the __TOP__ to the __GROUND__.
  !! @note 
  !! __mm_dzlay__ is defined on the total number of layers. The last 
  !! value of __mm_dzlay__ is set to twice the altitude of the ground layer.
  !! @note This value corresponds to the thickness between the center of the 
  !! __GROUND__ layer and below the surface. It is arbitrary and not used.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_dzlay

  !> Spherical mode \(0^{th}\) order moment (\(m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m0aer_s
  !> Spherical mode \(3^{rd}\) order moment (\(m^{3}.m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m3aer_s
  !> Fractal mode \(0^{th}\) order moment (\(m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m0aer_f
  !> Fractal mode \(3^{rd}\) order moment (\(m^{3}.m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m3aer_f
  !> CCN \(0^{th}\) order moment (\(m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m0ccn
  !> CCN \(3^{rd}\) order moment (\(m^{3}.m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m3ccn

  !> Ice components 3rd order moments (\(m^{3}.m^{-3}\)).
  !!
  !! It is a 2D array with the vertical layers in first dimension, and the number of ice 
  !! components in the second.
  !! @note 
  !! Both [[mm_globals(module):mm_m3ice(variable)]] and [[mm_globals(module):mm_gazs(variable)]]
  !! share the same indexing (related to species order).
  REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: mm_m3ice

  !> Condensible species molar fraction (\(mol.mol^{-1}\)).
  !!
  !! It is a 2D array with the vertical layers in first dimension, and
  !! the number of condensible species in the second. 
  !! @note 
  !! Both [[mm_globals(module):mm_m3ice(variable)]] and [[mm_globals(module):mm_gazs(variable)]]
  !! share the same indexing (related to species order).
  REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: mm_gazs

  !> Spherical mode characteristic radius (m).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_rcs
  !> Fractal mode characteristic radius (m).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_rcf
  !> Mean Drop radius (m).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_drad
  !> Mean Drop density (\(kg.m^{-3}\)).
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_drho

  !> Aerosols precipitations (m).
  !!
  !! Aerosols precipitations take into account both spherical and fractal modes.
  !! It is updated in [[mm_haze(module):mm_haze_microphysics(subroutine)]].
  REAL(kind=mm_wp), SAVE :: mm_aer_prec = 0._mm_wp

  !> Spherical mode \(M_{0}\) settling velocity (\(m.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the settling velocity for 
  !! the \(0^{th}\) order moment of the spherical mode.
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m0as_vsed

  !> Spherical mode \(M_{3}\) settling velocity (\(m.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the settling velocity for the 
  !! \(3^{rd}\) order moment of the spherical mode.
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m3as_vsed

  !> Fractal mode \(M_{0}\) settling velocity (\(m.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the settling velocity for the 
  !! \(0^{th}\) order moment of the fractal mode.
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m0af_vsed

  !> Fractal mode \(M_{3}\) settling velocity (\(m.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the settling velocity for the 
  !! \(3^{rd}\) order moment of the fractal mode.
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_m3af_vsed

  !> Spherical aerosol mass fluxes (\(kg.m^{-2}.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the mass fluxes for spherical aerosols.
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_aer_s_flux

  !> Fractal aerosol mass fluxes (\(kg.m^{-2}.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the mass fluxes for fractal aerosols
  !! It is updated in [[mm_haze(module):mm_haze_sedimentation(subroutine)]].
  !! @note 
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_aer_f_flux

  !> CCN precipitations (m).
  !! It is updated in [[mm_clouds(module):mm_cloud_microphysics(subroutine)]].
  REAL(kind=mm_wp), SAVE :: mm_ccn_prec = 0._mm_wp

  !> CCN mass fluxes (\(kg.m^{-2}.s^{-1}\)).
  !!
  !! It is a vector with the vertical layers that contains the 
  !! mass fluxes for CCN. 
  !! It is updated in [[mm_clouds(module):mm_cloud_microphysics(subroutine)]].
  !! @note
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_ccn_flux

  !> Ice components precipitations (m).
  !!
  !! It is a vector of [[mm_globals(module):mm_nesp(variable)]] values which share the same indexing 
  !! than [[mm_globals(module):mm_m3ice(variable)]] and [[mm_globals(module):mm_gazs(variable)]].
  !! It is updated in [[mm_clouds(module):mm_cloud_microphysics(subroutine)]].
  !! @note
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE, SAVE :: mm_ice_prec

  !> Ice components sedimentation fluxes (\(kg.m^{-2}.s-1\)).
  !!
  !! It is a 2D-array with the vertical layers in first dimension and the number of ice components 
  !! in the second. It is updated in [[mm_clouds(module):mm_cloud_microphysics(subroutine)]].
  !! @note
  !! This variable is always negative.
  REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: mm_ice_fluxes

  !> Condensible species saturation ratio (--).
  !!
  !! It is a 2D-array with the vertical layers in first dimension and the number of condensible 
  !! species in the second.
  !! It is updated in [[mm_clouds(module):mm_cloud_microphysics(subroutine)]].
  REAL(kind=mm_wp), DIMENSION(:,:), ALLOCATABLE, SAVE :: mm_gazs_sat

  !> [[mm_globals(module):mm_global_init(interface)]] initialization control flag.
  LOGICAL, PUBLIC, SAVE :: mm_ini     = .false.

  !> [[mm_globals(module):mm_column_init(function)]] initialization control flag.
  LOGICAL, PUBLIC, SAVE :: mm_ini_col = .false.

  !> [[mm_globals(module):mm_aerosols_init(function)]] initialization control flag.
  LOGICAL, PUBLIC, SAVE :: mm_ini_aer = .false.

  !> [[mm_globals(module):mm_clouds_init(function)]] initialization control flag.
  LOGICAL, PUBLIC, SAVE :: mm_ini_cld = .false.

  !> Interface to cloud properties methods.
  !!
  !! The method computes clouds properties (mean drop radius and denstity) from their afferent
  !! moments. It is overloaded to compute properties at a single level or over all the vertical
  !! atmospheric structure.
  INTERFACE mm_cloud_properties
    MODULE PROCEDURE cldprop_sc,cldprop_ve
  END INTERFACE

  !> Interface to global initialization.
  !!
  !! The method performs the global initialization of the model.
  !! @warning
  !! If OpenMP is activated, this subroutine must be called in an $OMP SINGLE statement as it 
  !! initializes global variable that are not thread private.
  !!
  !! '''
  !! !$OMP SINGLE
  !! call mm_global_init(...)
  !! !$OMP END SINGLE
  INTERFACE mm_global_init
    MODULE PROCEDURE mm_global_init_0,mm_global_init_1
  END INTERFACE

  !> Check an option from the configuration system.
  !!
  !! The method checks for an option in the configuration system and optionally
  !! set a default value if the option is not found. This is an overloaded method
  !! that can take in input either a floating point, integer, logical or string
  !! option value. 
  INTERFACE mm_check_opt
    MODULE PROCEDURE check_r1,check_i1,check_l1,check_s1
  END INTERFACE

  ! --- OPENMP ---------------
  ! All variable related to column computations should be private to each thread
  !
  !$OMP THREADPRIVATE(mm_ini_col,mm_ini_aer,mm_ini_cld)
  !$OMP THREADPRIVATE(mm_zlay,mm_zlev,mm_play,mm_plev,mm_temp,mm_rhoair,mm_btemp,mm_dzlev,mm_dzlay)
  !$OMP THREADPRIVATE(mm_m0aer_s,mm_m3aer_s,mm_m0aer_f,mm_m3aer_f)
  !$OMP THREADPRIVATE(mm_m0ccn,mm_m3ccn,mm_m3ice,mm_gazs)
  !$OMP THREADPRIVATE(mm_rcs,mm_rcf,mm_drad,mm_drho)
  !$OMP THREADPRIVATE(mm_aer_s_flux,mm_aer_f_flux,mm_ccn_flux,mm_ice_prec,mm_ice_fluxes,mm_gazs_sat)
  !$OMP THREADPRIVATE(mm_m0as_vsed,mm_m3as_vsed,mm_m0af_vsed,mm_m3af_vsed)

  !$OMP THREADPRIVATE(mm_nla,mm_nle)

  ! --------------------------


  CONTAINS 

  FUNCTION mm_global_init_0(dt,df,rm,rho_aer,p_prod,tx_prod,rc_prod,rplanet,g0, &
                            air_rad,air_mmol,coag_interactions,clouds,spcfile,  &
                            w_haze_prod,w_haze_sed,w_haze_coag,w_cloud_nucond,  &
                            w_cloud_sed,force_wsed_to_m0,force_wsed_to_m3,      &
                            no_fiadero,fiadero_min,fiadero_max) RESULT(err)
    !! Initialize global parameters of the model.
    !! 
    !! The function initializes all the global parameters of the model from direct input.
    !! Boolean (and Fiadero) parameters are optional as they are rather testing parameters. Their 
    !! default values are suitable for production runs.  
    !! @note
    !! If the method fails to initialize parameters (i.e. returned error is not 0). Then the model
    !! should probably be aborted as the global variables of the model will not be correctly setup.
    !! @warning
    !! If OpenMP is activated, this subroutine must be called in an $OMP SINGLE statement as it 
    !! initializes global variable that are not thread private.
    !!
    !! '''
    !! !$OMP SINGLE
    !! call mm_global_init_0(...)
    !! !$OMP END SINGLE
    REAL(kind=mm_wp), INTENT(in)           :: dt
      !! Microphysics timestep in seconds.
    REAL(kind=mm_wp), INTENT(in)           :: df
      !! Fractal dimension of fractal aerosol.
    REAL(kind=mm_wp), INTENT(in)           :: rm
      !! Monomer radius in meter.
    REAL(kind=mm_wp), INTENT(in)           :: rho_aer
      !! Aerosol density in \(kg.m^{-3}\).
    REAL(kind=mm_wp), INTENT(in)           :: p_prod
      !!  Aerosol production pressure level in Pa.
    REAL(kind=mm_wp), INTENT(in)           :: tx_prod
      !! Spherical aerosol mode production rate in \(kg.m^{-2}.s^{-1}\).
    REAL(kind=mm_wp), INTENT(in)           :: rc_prod
      !! Spherical mode characteristic radius for production in meter.
    REAL(kind=mm_wp), INTENT(in)           :: rplanet
      !! Planet radius in meter
    REAL(kind=mm_wp), INTENT(in)           :: g0
      !! Planet gravity acceleration at ground level in \(m.s^{-2}\).
    REAL(kind=mm_wp), INTENT(in)           :: air_rad
      !! Air molecules mean radius in meter.
    REAL(kind=mm_wp), INTENT(in)           :: air_mmol
      !! Air molecules mean molar mass in \(kg.mol^{-1}\).
    INTEGER, INTENT(in)                    :: coag_interactions
      !! Coagulation interactions process control flag.
    LOGICAL, INTENT(in)                    :: clouds
      !! Clouds microphysics control flag.
    CHARACTER(len=*), INTENT(in)           :: spcfile
      !! Clouds microphysics condensible species properties file.
    REAL(kind=mm_wp), INTENT(in), OPTIONAL :: fiadero_max
      !! Maximum moment ratio threshold for Fiadero correction (default: 10.) .
    REAL(kind=mm_wp), INTENT(in), OPTIONAL :: fiadero_min
      !! Minimum moment ratio threshold for Fiadero correction (default: 0.1).
    LOGICAL, INTENT(in), OPTIONAL          :: w_haze_prod
      !! Haze microphysics production process control flag (default: T).
    LOGICAL, INTENT(in), OPTIONAL          :: w_haze_sed
      !! Haze microphysics sedimentation process control flag (default: T).
    LOGICAL, INTENT(in), OPTIONAL          :: w_haze_coag
      !! Haze microphysics coagulation process control flag (default: T).
    LOGICAL, INTENT(in), OPTIONAL          :: w_cloud_sed
      !! Cloud microphysics nucleation/conensation process control flag (default: __clouds__ value).
    LOGICAL, INTENT(in), OPTIONAL          :: w_cloud_nucond
      !! Cloud microphysics production process control flag (default: __clouds__ value).
    LOGICAL, INTENT(in), OPTIONAL          :: no_fiadero
      !! Disable Fiadero correction for haze sedimentation process (default: F).
    LOGICAL, INTENT(in), OPTIONAL          :: force_wsed_to_m0
      !! force __all__ aerosols moments to fall at M0 settling velocity (default: T).
    LOGICAL, INTENT(in), OPTIONAL          :: force_wsed_to_m3 
      !! Force __all__ aerosols moments to fall at M3 settling velocity (default: F).
    TYPE(error) :: err
      !! Error status of the function.
    INTEGER                                           :: i
    TYPE(cfgparser)                                   :: cp
    CHARACTER(len=st_slen)                            :: spcpath
    CHARACTER(len=:), ALLOCATABLE                     :: defmsg
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: species
    REAL(kind=mm_wp)                                  :: zfiamin,zfiamax
    LOGICAL                                           :: zwhp,zwhs,zwhc,zwcs,zwcn,znofia, &
                                                         zwstom0,zwstom3

    zwhp = .true. ; zwhs = .true. ; zwhc = .true.
    zwcs = clouds ; zwcn = clouds 
    znofia = .false. ; zfiamin = 0.1_mm_wp ; zfiamax = 10._mm_wp
    zwstom0 = .true. ; zwstom3 = .false.
    err = noerror
    IF (mm_ini) THEN
      err = error("mm_global_init: YAMMS global initialization already performed !",-1)
      RETURN
    ENDIF

    ! Store options values in global variables...
    mm_df          = df 
    mm_rm          = rm 
    mm_rb2ra       = mm_rm**((mm_df-3._mm_wp)/mm_df) ! conversion factor for bulk -> fractal radius
    mm_rhoaer      = rho_aer 
    mm_p_prod      = p_prod
    mm_tx_prod     = tx_prod
    mm_rc_prod     = rc_prod
    mm_rpla        = rplanet
    mm_g0          = g0
    mm_dt          = dt
    mm_air_rad     = mm_air_rad
    mm_air_mmol    = air_mmol
    mm_coag_choice = coag_interactions
    ! check coagulation interactions choice
    IF (mm_coag_choice < 0 .OR. mm_coag_choice > 7) THEN
      err = error("mm_global_init: Invalid choice for coagulation interactions activation",-1)
      RETURN
    ENDIF

    mm_w_clouds = clouds

    ! Check clouds microphysics species file
    ! (only if clouds is activated)
    IF (mm_w_clouds) THEN
      IF (LEN_TRIM(spcfile) == 0) THEN
        err = error("mm_global_init: No species properties file given",-1)
        RETURN
      ENDIF
      ! Reads species properties configuration file  
      err = cfg_read_config(cp,TRIM(spcfile)) ; IF (err /= 0) RETURN
      err = cfg_get_value(cp,"used_species",species) 
      IF (err /= 0) THEN
        err = error("mm_global_init: cannot retrieve 'used_species' values",-1)
        RETURN
      ENDIF
      ! Now attempts to find species properties !!!
      mm_nesp = SIZE(species)
      ALLOCATE(mm_spcname(mm_nesp),mm_xESPS(mm_nesp))
      DO i=1,mm_nesp
        mm_spcname(i) = to_lower(species(i))
        IF(.NOT.cfg_has_section(cp,TRIM(mm_spcname(i)))) THEN
          err = error("mm_global_init: Cannot find "//TRIM(mm_spcname(i))//" properties",-1)
          RETURN
        ELSE
          err = read_esp(cp,TRIM(mm_spcname(i)),mm_xESPS(i))
          ! compute conversion factor: mol.mol-1 => kg.kg-1
          mm_xESPS(i)%fmol2fmas = mm_xESPS(i)%masmol / mm_air_mmol
          IF (err/=0) THEN
            err = error("mm_global_init: "//TRIM(mm_spcname(i))//": "//TRIM(err%msg),-1)
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    ! optional flags
    ! haze control flags
    IF (PRESENT(w_haze_prod)) THEN 
      mm_w_haze_prod = w_haze_prod
    ELSE 
      mm_w_haze_prod = zwhp 
      call printw("mm_haze_production",to_string(mm_w_haze_prod))
    ENDIF
    IF (PRESENT(w_haze_sed)) THEN 
      mm_w_haze_sed = w_haze_sed
    ELSE 
      mm_w_haze_sed = zwhs 
      call printw("mm_haze_sedimentation",to_string(mm_w_haze_sed))
    ENDIF
    IF (PRESENT(w_haze_coag)) THEN 
      mm_w_haze_coag = w_haze_coag
    ELSE 
      mm_w_haze_coag = zwhc
      call printw("mm_haze_coagulation",to_string(mm_w_haze_coag))
    ENDIF
    IF (PRESENT(force_wsed_to_m0)) THEN 
      mm_wsed_m0 = force_wsed_to_m0
    ELSE 
      mm_wsed_m0 = zwstom0
      call printw("mm_wsed_m0",to_string(mm_wsed_m0))
    ENDIF
    IF (PRESENT(force_wsed_to_m3)) THEN 
      mm_wsed_m3 = force_wsed_to_m3
    ELSE 
      mm_wsed_m3 = zwstom3
      call printw("mm_wsed_m3",to_string(mm_wsed_m3))
    ENDIF
    IF (PRESENT(no_fiadero)) THEN 
      mm_no_fiadero_w = no_fiadero
    ELSE 
      mm_no_fiadero_w = znofia 
      call printw("mm_no_fiadero",to_string(mm_no_fiadero_w))
    ENDIF
    IF (PRESENT(fiadero_min)) THEN 
      mm_fiadero_min = fiadero_min
    ELSE 
      mm_fiadero_min = zfiamin
      call printw("mm_fiadero_min",to_string(mm_fiadero_min))
    ENDIF
    IF (PRESENT(fiadero_max)) THEN 
      mm_fiadero_max = fiadero_max
    ELSE 
      mm_fiadero_max = zfiamax
      call printw("mm_fiadero_max",to_string(mm_fiadero_max))
    ENDIF
    ! clouds control flags
    IF (mm_w_clouds) THEN
      IF (PRESENT(w_cloud_sed)) THEN 
        mm_w_cloud_sed = w_cloud_sed
      ELSE 
        mm_w_cloud_sed = zwcs 
        call printw("mm_cloud_sed",to_string(mm_w_cloud_sed)) 
      ENDIF
      IF (PRESENT(w_cloud_nucond)) THEN 
        mm_w_cloud_nucond = w_cloud_nucond
      ELSE 
        mm_w_cloud_nucond = zwcs
        call printw("mm_cloud_nucond",to_string(mm_w_cloud_nucond)) 
      ENDIF
    ENDIF

    ! check w sed flags
    err = noerror
    ! special check for settling velocity
    IF (mm_wsed_m0 .AND. mm_wsed_m3) THEN
      err = error("'wsed_m0' and 'wsed_m3' options are mutually exclusive",-1)
    ENDIF
    mm_ini = err == noerror

    CONTAINS

    SUBROUTINE printw(string,value)
      !! Print a warning message.
      CHARACTER(len=*), INTENT(in) :: string !! Name of the option.
      CHARACTER(len=*), INTENT(in) :: value  !! (string) Value of the option.
      IF (mm_log) &
      WRITE(*,'(a,a,a)') "warning: Parameter "//string//"not given... Using default value: "//value
    END SUBROUTINE printw 
  END FUNCTION mm_global_init_0

  FUNCTION mm_global_init_1(cfg) RESULT(err)
    !! Set global configuration from a configuration file.
    !!
    !! See [[mm_globals(module):mm_global_init_0(function)]].
    TYPE(cfgparser), INTENT(in) :: cfg  !! Configuration file path.
    TYPE(error) :: err                  !! Error status of the function.
    INTEGER                                           :: i
    TYPE(cfgparser)                                   :: spccfg
    CHARACTER(len=st_slen)                            :: spcpath
    CHARACTER(len=:), ALLOCATABLE                     :: defmsg
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: species
    REAL(kind=mm_wp)                                  :: zfiamin,zfiamax
    LOGICAL                                           :: zwhp,zwhs,zwhc,zwcs,zwcn,znofia, &
                                                         zwstom0,zwstom3

    err = noerror

    IF (mm_ini) THEN
      err = error("mm_global_init: YAMMS global initialization already performed !",-1)
      RETURN
    ENDIF

    ! MP2M mandatory parameters
    err = mm_check_opt(cfg_get_value(cfg,"df",mm_df),mm_df,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"rm",mm_rm),mm_rm,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"rho_aer",mm_rhoaer),mm_rhoaer,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"p_prod",mm_p_prod),mm_p_prod,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"tx_prod",mm_tx_prod),mm_tx_prod,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"rc_prod",mm_rc_prod),mm_rc_prod,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"planet_radius",mm_rpla),mm_rpla,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"g0",mm_g0),mm_g0,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"timestep",mm_dt),mm_dt,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"air_radius",mm_air_rad),mm_air_rad,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"air_molarmass",mm_air_mmol),mm_air_mmol,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"haze_coag_interactions",mm_coag_choice),mm_coag_choice,wlog=mm_log)
    IF (err/=0) RETURN
    err = mm_check_opt(cfg_get_value(cfg,"clouds_microphysics",mm_w_clouds),mm_w_clouds,wlog=mm_log)
    IF (err/=0) RETURN

    ! computes the conversion factor for bulk -> fractal radius
    mm_rb2ra = mm_rm**((mm_df-3._mm_wp)/mm_df)

    ! Check coagulation interactions choice
    IF (mm_coag_choice < 0 .OR. mm_coag_choice > 7) THEN
      err = error("mm_global_init: Invalid choice for coagulation interactions activation",-1)
      RETURN
    ENDIF

    ! Check clouds microphysics input
    ! it is read only if clouds is activated. We must to check if it is self-consistent...
    IF (mm_w_clouds) THEN
      ! Gets species property file path
      err = cfg_get_value(cfg,'specie_cfg',spcpath) ; IF (err /= 0) RETURN
      ! Reads species properties configuration file  
      err = cfg_read_config(spccfg,trim(spcpath)) ; IF (err /= 0) RETURN
      err = cfg_get_value(spccfg,"used_species",species) 
      IF (err /= 0) THEN
        err = error("mm_global_init: cannot retrieve 'used_species' values",-1)
        RETURN
      ENDIF
      ! Now attempts to find specides properties !!!
      mm_nesp = SIZE(species)
      ALLOCATE(mm_spcname(mm_nesp),mm_xESPS(mm_nesp))
      !mm_spcname(1:mm_nesp) = species(:)
      DO i=1,mm_nesp
        mm_spcname(i) = to_lower(species(i))
        IF (.NOT.cfg_has_section(spccfg,TRIM(mm_spcname(i)))) THEN
          err = error("mm_global_init: Cannot find "//TRIM(mm_spcname(i))//" properties",-1)
          RETURN
        ELSE
          err = read_esp(spccfg,TRIM(mm_spcname(i)),mm_xESPS(i))
          ! compute conversion factor: mol.mol-1 => kg.kg-1
          mm_xESPS(i)%fmol2fmas = mm_xESPS(i)%masmol / mm_air_mmol
          IF (err/=0) THEN
            err = error(TRIM(mm_spcname(i))//": "//TRIM(err%msg),-2)
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    zwhp = .true. ; zwhs = .true. ; zwhc = .true.
    zwcs = mm_w_clouds ; zwcn = mm_w_clouds
    znofia = .false. ; zfiamin = 0.1_mm_wp ; zfiamax = 10._mm_wp
    zwstom0 = .true. ; zwstom3 = .false.

    ! MP2M Optional parameters
    err = mm_check_opt(cfg_get_value(cfg,"haze_production",mm_w_haze_prod),mm_w_haze_prod,zwhp,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"haze_sedimentation",mm_w_haze_sed),mm_w_haze_sed,zwhs,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"haze_coagulation",mm_w_haze_coag),mm_w_haze_coag,zwhc,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"clouds_sedimentation",mm_w_cloud_sed),mm_w_cloud_sed,zwcs,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"clouds_nucl_cond",mm_w_cloud_nucond),mm_w_cloud_nucond,zwcn,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"wsed_m0",mm_wsed_m0),mm_wsed_m0,zwstom0,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"wsed_m3",mm_wsed_m3),mm_wsed_m3,zwstom3,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"no_fiadero",mm_no_fiadero_w),mm_no_fiadero_w,znofia,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"fiadero_min_ratio",mm_fiadero_min),mm_fiadero_min,zfiamin,wlog=mm_log)
    err = mm_check_opt(cfg_get_value(cfg,"fiadero_max_ratio",mm_fiadero_max),mm_fiadero_max,zfiamax,wlog=mm_log)

    err = noerror
    ! special check for settling velocity
    IF (mm_wsed_m0 .AND. mm_wsed_m3) THEN
      err = error("'wsed_m0' and 'wsed_m3' options are mutually exclusive",-1)
    ENDIF
    mm_ini = err == noerror
  END FUNCTION mm_global_init_1

  FUNCTION mm_column_init(plev,zlev,play,zlay,temp) RESULT(err)
    !! Initialize vertical atmospheric fields.
    !! 
    !! This subroutine initializes vertical fields needed by the microphysics:
    !!
    !! 1. Save reversed input field into "local" array 
    !! 2. Compute thicknesses layers and levels
    !! 3. Interpolate temperature at levels
    !!
    !! The method should be called whenever the vertical structure of the atmosphere changes.
    !!
    !! @attention
    !! All the input vectors should be defined from __GROUND__ to __TOP__ of the atmosphere,
    !! otherwise nasty things will occur in computations. 
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: plev !! Pressure levels (Pa).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: zlev !! Altitude levels (m).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: play !! Pressure layers (Pa).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: zlay !! Altitude at the center of each layer (m).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: temp !! Temperature at the center of each layer (K).
    TYPE(error) :: err                                 !! Error status of the function.
    INTEGER :: i
    mm_ini_col = .false.                          
    err = noerror
    IF (.NOT.mm_ini) THEN
      err = error("mm_column_init: Global initialization not done yet",-1)
      RETURN
    ENDIF
    IF (mm_nla < 0) THEN
      mm_nla = SIZE(play)
    ELSE
      IF (mm_nla /= SIZE(play)) THEN
        err = error("mm_column_init: mm_nla cannot be modified dynamically within the run",-1)
        RETURN
      ENDIF
    ENDIF
    IF (mm_nle < 0) THEN
      mm_nle = SIZE(plev)
    ELSE
      IF (mm_nle /= SIZE(plev)) THEN
        err = error("mm_column_init: mm_nle cannot be modified dynamically within the run",-1)
        RETURN
      ENDIF
    ENDIF
    ! should be trashed soon or later
    IF (mm_nla+1 /= mm_nle) THEN
      err = error("mm_column_init: Inconsistent number of layers/levels",-1)
      RETURN
    ENDIF
    ! Allocates if required
    IF (.NOT.ALLOCATED(mm_plev))   ALLOCATE(mm_plev(mm_nle))
    IF (.NOT.ALLOCATED(mm_zlev))   ALLOCATE(mm_zlev(mm_nle))
    IF (.NOT.ALLOCATED(mm_play))   ALLOCATE(mm_play(mm_nla))
    IF (.NOT.ALLOCATED(mm_zlay))   ALLOCATE(mm_zlay(mm_nla))
    IF (.NOT.ALLOCATED(mm_temp))   ALLOCATE(mm_temp(mm_nla))
    IF (.NOT.ALLOCATED(mm_btemp))  ALLOCATE(mm_btemp(mm_nle))
    IF (.NOT.ALLOCATED(mm_dzlev))  ALLOCATE(mm_dzlev(mm_nla))
    IF (.NOT.ALLOCATED(mm_dzlay))  ALLOCATE(mm_dzlay(mm_nla))
    IF (.NOT.ALLOCATED(mm_rhoair)) ALLOCATE(mm_rhoair(mm_nla))
    ! Saves reversed input vectors
    mm_zlay = zlay(mm_nla:1:-1) ; mm_zlev = zlev(mm_nle:1:-1)
    mm_play = play(mm_nla:1:-1) ; mm_plev = plev(mm_nle:1:-1)
    mm_temp = temp(mm_nla:1:-1)
    ! Computes others vectors
    mm_dzlay(1:mm_nla-1) = mm_zlay(1:mm_nla-1)-mm_zlay(2:mm_nla)
    mm_dzlay(mm_nla)     = mm_dzlay(mm_nla-1) ! actually arbitrary
    mm_dzlev(1:mm_nla)   = mm_zlev(1:mm_nle-1)-mm_zlev(2:mm_nle)
    mm_btemp(2:mm_nla)   = (mm_temp(1:mm_nla-1)+mm_temp(2:mm_nla))/2._mm_wp
    mm_btemp(1)          = mm_temp(1)
    mm_btemp(mm_nle)     = mm_temp(mm_nla)+(mm_temp(mm_nla)-mm_temp(mm_nla-1))/2._mm_wp
    ! Hydrostatic equilibrium
    mm_rhoair(1:mm_nla) = (mm_plev(2:mm_nle)-mm_plev(1:mm_nla)) / &
                          (mm_effg(mm_zlay)*mm_dzlev)
    mm_ini_col = .true.                          
    ! write out profiles (only if BOTH debug and log are enabled).
    IF (mm_log.AND.mm_debug) THEN
      WRITE(*,'(a)') '# TEMP             PLAY             ZLAY             DZLAY            RHOAIR'
      DO i=1,mm_nla
        WRITE(*,'(5(ES15.7,2X))') mm_temp(i),mm_play(i),mm_zlay(i),mm_dzlay(i), mm_rhoair(i)
      ENDDO
      WRITE(*,'(a)') '# TEMP             PLEV             ZLEV             DZLEV'
      DO i=1,mm_nle
        IF (i /= mm_nle) THEN
          WRITE(*,'(4(ES15.7,2X))') mm_btemp(i),mm_plev(i),mm_zlev(i),mm_dzlev(i)
        ELSE
          WRITE(*,'(3(ES15.7,2X))') mm_btemp(i),mm_plev(i),mm_zlev(i)
        ENDIF
      ENDDO
    ENDIF

    RETURN
  END FUNCTION mm_column_init

  FUNCTION mm_aerosols_init(m0aer_s,m3aer_s,m0aer_f,m3aer_f) RESULT(err)
    !! Initialize clouds tracers vertical grid.
    !! 
    !! The subroutine initializes aerosols microphysics tracers columns. It allocates variables if 
    !! required and stores input vectors in reversed order. It also computes the characteristic radii 
    !! of each mode. 
    !! @note
    !! All the input arguments should be defined from ground to top. 
    !!
    !! @attention
    !! [[mm_globals(module):mm_global_init(interface)]] and [[mm_globals(module):mm_column_init(function)]]
    !! must have been called at least once before this method is called. Moreover, this method should be
    !! after each call of [[mm_globals(module):mm_column_init(function)]] to reflect changes in the 
    !! vertical atmospheric structure.
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: m0aer_s !! \(0^{th}\) order moment of the spherical mode (\(m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: m3aer_s !! \(3^{rd}\) order moment of the spherical mode (\(m^{3}.m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: m0aer_f !! \(0^{th}\) order moment of the fractal mode (\(m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in) :: m3aer_f !! \(3^{rd}\) order moment of the fractal mode (\(m^{3}.m^{-2}\)).
    TYPE(error) :: err                                    !! Error status of the function.
    err = noerror
    IF (.NOT.mm_ini) THEN
      err = error("mm_aerosols_init: Global initialization not done yet",-1) ; RETURN
    ENDIF
    IF (.NOT.mm_ini_col) THEN
      err = error("mm_aerosols_init: Column initialization not done yet",-1) ; RETURN
    ENDIF
    ! Check input size ???
    IF (SIZE(m0aer_s) /= mm_nla) THEN
      err = error("mm_aerosols_init: Invalid size for input arrays",-1) ; RETURN
    ENDIF

    ! Allocate variable if required
    IF (.NOT.ALLOCATED(mm_m0aer_s)) ALLOCATE(mm_m0aer_s(mm_nla))
    IF (.NOT.ALLOCATED(mm_m3aer_s)) ALLOCATE(mm_m3aer_s(mm_nla))
    IF (.NOT.ALLOCATED(mm_m0aer_f)) ALLOCATE(mm_m0aer_f(mm_nla))
    IF (.NOT.ALLOCATED(mm_m3aer_f)) ALLOCATE(mm_m3aer_f(mm_nla))
    IF (.NOT.ALLOCATED(mm_rcs))     ALLOCATE(mm_rcs(mm_nla))
    IF (.NOT.ALLOCATED(mm_rcf))     ALLOCATE(mm_rcf(mm_nla))
    ! Allocate memory for diagnostics
    IF (.NOT.ALLOCATED(mm_aer_s_flux)) THEN
      ALLOCATE(mm_aer_s_flux(mm_nla)) ; mm_aer_s_flux(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_aer_f_flux)) THEN
      ALLOCATE(mm_aer_f_flux(mm_nla)) ; mm_aer_f_flux(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_m0as_vsed)) THEN
      ALLOCATE(mm_m0as_vsed(mm_nla)) ; mm_m0as_vsed(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_m3as_vsed)) THEN
      ALLOCATE(mm_m3as_vsed(mm_nla)) ; mm_m3as_vsed(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_m0af_vsed)) THEN
      ALLOCATE(mm_m0af_vsed(mm_nla)) ; mm_m0af_vsed(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_m3af_vsed)) THEN
      ALLOCATE(mm_m3af_vsed(mm_nla)) ; mm_m3af_vsed(:) = 0._mm_wp
    ENDIF
    ! note : mm_dzlev is already from top to ground
    mm_m0aer_s = m0aer_s(mm_nla:1:-1)/mm_dzlev(:)
    mm_m3aer_s = m3aer_s(mm_nla:1:-1)/mm_dzlev(:)
    mm_m0aer_f = m0aer_f(mm_nla:1:-1)/mm_dzlev(:)
    mm_m3aer_f = m3aer_f(mm_nla:1:-1)/mm_dzlev(:)
    ! aerosols characteristic radii
    ! il faudrait peut etre revoir la gestion des zeros ici et la 
    ! remplacer par une valeur seuil des moments.
    !
    !-> JVO 19 : Done. Zero threshold set at espilon from dynamics on the
    ! input moments in calmufi (safer than here). Might still be some unphysical
    ! values after the dynamics near the threshold. Could be a could idea to add
    ! a sanity check filtering too high values of radii.
    !
    ! TBD : Add a sanity check for high radii ????
    WHERE(mm_m3aer_s > 0._mm_wp .AND. mm_m0aer_s > 0._mm_wp)
      mm_rcs = mm_get_rcs(mm_m0aer_s,mm_m3aer_s)
    ELSEWHERE
      mm_rcs = 0._mm_wp
    ENDWHERE
    WHERE(mm_m3aer_f > 0._mm_wp .AND. mm_m0aer_f > 0._mm_wp)
      mm_rcf = mm_get_rcf(mm_m0aer_f,mm_m3aer_f)
    ELSEWHERE
      mm_rcf = 0._mm_wp
    ENDWHERE
    mm_ini_aer = .true.
  END FUNCTION mm_aerosols_init

  FUNCTION mm_clouds_init(m0ccn,m3ccn,m3ice,gazs) RESULT(err)
    !! Initialize clouds tracers vertical grid.
    !! 
    !! The subroutine initializes cloud microphysics tracers columns. It allocates variables if 
    !! required and stores input vectors in reversed order. It also computes the mean drop radius 
    !! and density and allocates diagnostic vectors.
    !! @note
    !! All the input arguments should be defined from ground to top. 
    !!
    !! @attention
    !! [[mm_globals(module):mm_global_init(interface)]] and [[mm_globals(module):mm_column_init(function)]]
    !! must have been called at least once before this method is called. Moreover, this method should be
    !! after each call of [[mm_globals(module):mm_column_init(function)]] to reflect changes in the 
    !! vertical atmospheric structure.
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in)   :: m0ccn !! 0th order moment of the CCN distribution (\(m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:), INTENT(in)   :: m3ccn !! 3rd order moment of the CCN distribution (\(m^{3}.m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(in) :: m3ice !! 3rd order moments of the ice components (\(m^{3}.m^{-2}\)).
    REAL(kind=mm_wp), DIMENSION(:,:), INTENT(in) :: gazs  !! Condensible species gazs molar fraction (\(mol.mol^{-1}\)).
    TYPE(error) :: err                                    !! Error status of the function.
    INTEGER :: i
    err = noerror
    IF (.NOT.mm_ini) THEN
      err = error("Global initialization not done yet",-8)
      RETURN
    ENDIF

    IF (.NOT.mm_w_clouds) THEN
      IF (mm_debug) WRITE(*,'(a)') "WARNING: Cloud microphysic is not enabled..."
      RETURN
    ENDIF

    ! Note:
    !  Here we could check that mm_nla is the same size of gazs(DIM=1)
    !  Actually, mm_nla should always initialized the first time mm_column_init is called, NOT HERE.
    IF (mm_nla < 0)  mm_nla  = SIZE(gazs,DIM=1)
    ! Note: 
    !   here we could check that mm_nesp is the same size of gazs(DIM=2)
    !   Actually, mm_nesp should be always initialized in mm_global_init, NOT HERE.
    IF (mm_nesp < 0) mm_nesp = SIZE(gazs,DIM=2)

    ! Allocate variable if required
    IF (.NOT.ALLOCATED(mm_m0ccn))   ALLOCATE(mm_m0ccn(mm_nla))
    IF (.NOT.ALLOCATED(mm_m3ccn))   ALLOCATE(mm_m3ccn(mm_nla))
    IF (.NOT.ALLOCATED(mm_m3ice))   ALLOCATE(mm_m3ice(mm_nla,mm_nesp))
    IF (.NOT.ALLOCATED(mm_gazs))    ALLOCATE(mm_gazs(mm_nla,mm_nesp))
    IF (.NOT.ALLOCATED(mm_drad))    ALLOCATE(mm_drad(mm_nla))
    IF (.NOT.ALLOCATED(mm_drho))    ALLOCATE(mm_drho(mm_nla))
    ! Allocate memory for diagnostics
    IF (.NOT.ALLOCATED(mm_ccn_flux)) THEN
      ALLOCATE(mm_ccn_flux(mm_nla)) ; mm_ccn_flux(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_ice_prec))   THEN
      ALLOCATE(mm_ice_prec(mm_nesp)) ; mm_ice_prec(:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_ice_fluxes)) THEN
      ALLOCATE(mm_ice_fluxes(mm_nla,mm_nesp)) ; mm_ice_fluxes(:,:) = 0._mm_wp
    ENDIF
    IF (.NOT.ALLOCATED(mm_gazs_sat)) THEN
      ALLOCATE(mm_gazs_sat(mm_nla,mm_nesp)) ; mm_gazs_sat(:,:) = 0._mm_wp
    ENDIF

    ! note mm_dzlev already from top to ground
    mm_m0ccn = m0ccn(mm_nla:1:-1)/mm_dzlev(:)
    mm_m3ccn = m3ccn(mm_nla:1:-1)/mm_dzlev(:)
    DO i=1,mm_nesp
      mm_m3ice(:,i) = m3ice(mm_nla:1:-1,i)/mm_dzlev(:)
      mm_gazs(:,i)  = gazs(mm_nla:1:-1,i)
    ENDDO
    ! drop mean radius
    call mm_cloud_properties(mm_m0ccn,mm_m3ccn,mm_m3ice,mm_drad,mm_drho)
    mm_ini_cld = .true.
  END FUNCTION mm_clouds_init

  SUBROUTINE mm_dump_parameters()
    !! Dump model global parameters on stdout.
    WRITE(*,'(a)')         "========= YAMMS PARAMETERS ============"
    WRITE(*,'(a,L2)')      "mm_w_haze_prod         : ", mm_w_haze_prod
    WRITE(*,'(a,ES14.7)')  "   mm_p_prod           : ", mm_p_prod
    WRITE(*,'(a,ES14.7)')  "   mm_tx_prod          : ", mm_tx_prod
    WRITE(*,'(a,ES14.7)')  "   mm_rc_prod          : ", mm_rc_prod
    WRITE(*,'(a,L2)')      "mm_w_haze_coag         : ", mm_w_haze_coag
    WRITE(*,'(a,I2.2)')    "   mm_coag_interactions: ", mm_coag_choice
    WRITE(*,'(a,L2)')      "mm_w_haze_sed          : ", mm_w_haze_sed 
    WRITE(*,'(a,L2)')      "   mm_wsed_m0          : ", mm_wsed_m0
    WRITE(*,'(a,L2)')      "   mm_wsed_m3          : ", mm_wsed_m3
    WRITE(*,'(a,L2)')      "   mm_no_fiadero_w     : ", mm_no_fiadero_w
    WRITE(*,'(a,ES14.7)')  "   mm_fiadero_min      : ", mm_fiadero_min
    WRITE(*,'(a,ES14.7)')  "   mm_fiadero_max      : ", mm_fiadero_max
    WRITE(*,'(a,L2)')      "mm_w_clouds            : ", mm_w_clouds
    WRITE(*,'(a,L2)')      "   mm_w_cloud_sed      : ", mm_w_cloud_sed
    WRITE(*,'(a,L2)')      "   mm_w_cloud_nucond   : ", mm_w_cloud_nucond
    WRITE(*,'(a)')         "---------------------------------------"
    WRITE(*,'(a,ES14.7)')  "mm_dt                  : ", mm_dt
    IF (mm_nla > -1) THEN
      WRITE(*,'(a,I3.3)')    "mm_nla                 : ", mm_nla
    ELSE
      WRITE(*,'(a)')         "mm_nla                 : not initialized yet"
    ENDIF
    WRITE(*,'(a,ES14.7)')  "mm_df                  : ", mm_df
    WRITE(*,'(a,ES14.7)')  "mm_rm                  : ", mm_rm
    WRITE(*,'(a,ES14.7)')  "mm_rpla                : ", mm_rpla
    WRITE(*,'(a,ES14.7)')  "mm_g0                  : ", mm_g0
    WRITE(*,'(a)')         "======================================="
  END SUBROUTINE mm_dump_parameters

  ELEMENTAL FUNCTION mm_get_rcs(m0,m3) RESULT(res)
    !! Get the characteristic radius for the spherical aerosols size distribution.
    !! 
    !! The method computes the characteristic radius of the size distribution law
    !! of the spherical aerosols mode according to its moments and its inter-moments
    !! relation.
    REAL(kind=mm_wp), INTENT(in) :: m0 !! \(0^{th}\) order moment 
    REAL(kind=mm_wp), INTENT(in) :: m3 !! \(3^{rd}\) order moment
    REAL(kind=mm_wp) :: res            !! Radius
    ! arbitrary: if there is no way to compute radius 
    IF (m3 <= 0._mm_wp .OR. m0 <= 0._mm_wp) res = 1._mm_wp 
    res = (m3/m0/mm_alpha_s(3._mm_wp))**(1._mm_wp/3._mm_wp)
  END FUNCTION mm_get_rcs

  ELEMENTAL FUNCTION mm_get_rcf(m0,m3) RESULT(res)
    !! Get the characteristic radius for the fractal aerosols size distribution.
    !! 
    !! The method computes the characteristic radius of the size distribution law
    !! of the fractal aerosols mode according to its moments and its inter-moments
    !! relation.
    REAL(kind=mm_wp), INTENT(in) :: m0 !! \(0^{th}\) order moment 
    REAL(kind=mm_wp), INTENT(in) :: m3 !! \(3^{rd}\) order moment
    REAL(kind=mm_wp) :: res            !! Radius
    ! arbitrary: if there is no way to compute radius 
    IF (m3 <= 0._mm_wp .OR. m0 <= 0._mm_wp) res = 1._mm_wp 
    res = (m3/m0/mm_alpha_f(3._mm_wp))**(1._mm_wp/3._mm_wp)
  END FUNCTION mm_get_rcf

  ELEMENTAL FUNCTION mm_effg(z) RESULT(effg) 
    !! Compute effective gravitational acceleration.
    REAL(kind=mm_wp), INTENT(in) :: z !! Altitude in meters
    REAL(kind=mm_wp) :: effg          !! Effective gravitational acceleration in \(m.s^{-2}\)
    effg = mm_g0
    IF (mm_use_effg) effg = effg * (mm_rpla/(mm_rpla+z))**2
    RETURN
  END FUNCTION mm_effg 

  !==================================
  ! --- private methods -------------
  !==================================

  SUBROUTINE cldprop_sc(m0ccn,m3ccn,m3ice,drad,drho)
    !! Get cloud drop properties (scalar).
    !!
    !! The method computes the mean radius and mean density of cloud drops.
    !!
    !! @bug 
    !! A possible bug can happen because of threshold snippet. If __drad__ is greater than 
    !! __drmax__ (== 100 microns) it is automatically set to __drmax__, but computation of 
    !! __drho__ remains unmodified. So __drho__ is not correct in that case.
    !!
    !! @todo 
    !! Fix the bug of the subroutine, but it is rather minor, since theoretically we do not 
    !! need the density of the drop.
    !!
    !! @todo 
    !! Think about a better implementation of thresholds, and get sure of their consequences in 
    !! the other parts of the model. 
    REAL(kind=mm_wp), INTENT(in)               :: m0ccn !! \(0^{th}\) order moment of the ccn 
    REAL(kind=mm_wp), INTENT(in)               :: m3ccn !! \(3^{rd}\) order moment of the ccn 
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:) :: m3ice !! \(3^{rd}\) order moments of each ice component
    REAL(kind=mm_wp), INTENT(out)              :: drad  !! Output mean drop radius 
    REAL(kind=mm_wp), INTENT(out), OPTIONAL    :: drho  !! Optional output mean drop density
    REAL(kind=mm_wp)            :: vtot, wtot, ntot 
    REAL(kind=mm_wp), PARAMETER :: threshold = 1.e-25_mm_wp,         &  
                                       drmin = 1.e-10_mm_wp,         &
                                       drmax = 1.e-4_mm_wp,          &
                                      athird = 1._mm_wp/3._mm_wp,    &
                                       pifac = 4._mm_wp/3._mm_wp*mm_pi 
    drad = 0._mm_wp
    ntot = m0ccn
    vtot = pifac*m3ccn+SUM(m3ice)
    wtot = pifac*m3ccn*mm_rhoaer+SUM(m3ice*mm_xESPS(:)%rho)
    IF (ntot <= threshold .OR. vtot <= 0._mm_wp) THEN
      drad  = drmin
      IF (PRESENT(drho)) drho  = mm_rhoaer 
    ELSE
      drad = (vtot/ntot/pifac)**athird
      drad = MAX(MIN(drad,drmax),drmin)
      IF (PRESENT(drho)) drho = wtot/vtot
    ENDIF
    RETURN
  END SUBROUTINE cldprop_sc

  SUBROUTINE cldprop_ve(m0ccn,m3ccn,m3ice,drad,drho)
    !! Get cloud drop properties (vector).
    !!
    !! The method performs the same computations than [[mm_globals(module):cldprop_sc(subroutine)]]
    !! but for the entire vertical atmospheric structure.
    !! Same remarks apply here.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)            :: m0ccn !! 0th order moment of the ccn.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:)            :: m3ccn !! 3rd order moment of the ccn.
    REAL(kind=mm_wp), INTENT(in), DIMENSION(:,:)          :: m3ice !! 3rd order moments of each ice component.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:)           :: drad  !! Output mean drop radius.
    REAL(kind=mm_wp), INTENT(out), DIMENSION(:), OPTIONAL :: drho  !! Optional output mean drop density.
    INTEGER                                     :: i,ns
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE :: vtot,wtot,ntot,rho 
    REAL(kind=mm_wp), PARAMETER                 :: threshold = 1.e-25_mm_wp,         &  
                                                       drmin = 1.e-10_mm_wp,         &
                                                       drmax = 1.e-4_mm_wp,          &
                                                      athird = 1._mm_wp/3._mm_wp,    &
                                                       pifac = 4._mm_wp/3._mm_wp*mm_pi
                                                     
    ns = SIZE(m0ccn) ; ALLOCATE(vtot(ns),wtot(ns),ntot(ns),rho(ns))
    drad = 0._mm_wp
    ntot = m0ccn
    vtot = pifac*m3ccn+SUM(m3ice,DIM=2)
    wtot = pifac*m3ccn*mm_rhoaer
    DO i=1,SIZE(m3ice,DIM=2)
      wtot = wtot+m3ice(:,i)*mm_xESPS(i)%rho
    ENDDO
    WHERE(ntot <= threshold .OR. vtot <= 0._mm_wp) 
      drad  = drmin
      rho = mm_rhoaer
    ELSEWHERE 
      drad = (vtot/ntot/pifac)**athird
      drad = MAX(MIN(drad,drmax),drmin)
      rho = wtot/vtot
    END WHERE
    IF (PRESENT(drho)) drho  = rho
    RETURN
  END SUBROUTINE cldprop_ve

  ! For configuration file (requires fccp library).

  FUNCTION read_esp(parser,sec,pp) RESULT (err)
    !! Read and store [[mm_globals(module):mm_esp(type)]] parameters.
    TYPE(cfgparser), INTENT(in)   :: parser !! Configuration parser.
    CHARACTER(len=*), INTENT(in)  :: sec    !! Name of the specie (should match a section of the configuration.
    TYPE(mm_esp), INTENT(out)     :: pp     !! [[mm_globals(module):mm_esp(type)]] object that stores the parameters.
    TYPE(error)                   :: err    !! Error status of the function.
    err = cfg_get_value(parser,TRIM(sec)//'name',pp%name)       ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'mas',pp%mas)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'vol',pp%vol)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'ray',pp%ray)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'mas',pp%mas)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'vol',pp%vol)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'ray',pp%ray)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'masmol',pp%masmol)   ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'rho',pp%rho)         ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'tc',pp%tc)           ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'pc',pp%pc)           ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'tb',pp%tb)           ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'w',pp%w)             ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'a_sat',pp%a_sat)     ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'b_sat',pp%b_sat)     ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'c_sat',pp%c_sat)     ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'d_sat',pp%d_sat)     ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'mteta',pp%mteta)     ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'tx_prod',pp%tx_prod) ; IF (err /= 0) RETURN
    RETURN
  END FUNCTION read_esp

  ! =========================================================================
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                CONFIGURATION PARSER checking methods 
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! =========================================================================

  FUNCTION check_r1(err,var,def,wlog) RESULT(ret)
    !! Check an option value (float).
    !! 
    !! The method checks an option value and optionally set a default value, __def__ to initialize 
    !! __var__ on error if given.
    TYPE(error), INTENT(in)                :: err  !! Error object from value getter.
    REAL(kind=mm_wp), INTENT(inout)        :: var  !! Input/output option value.
    REAL(kind=mm_wp), INTENT(in), OPTIONAL :: def  !! Default value to set.
    LOGICAL, INTENT(in), OPTIONAL          :: wlog !! .true. to print warning/error message.
    TYPE(error) :: ret                             !! Input error.
    CHARACTER(len=*), PARAMETER :: defmsg = '... Using default value: ' 
    LOGICAL                     :: zlog
    ret = err
    zlog = .false. ; IF (PRESENT(wlog)) zlog = wlog
    IF (err == 0) RETURN
    IF (PRESENT(def)) THEN
      var = def
      IF (zlog) WRITE(*,'(a,a,a)') error_to_string(err,'',.true.),defmsg,to_string(var)
      ret = noerror
    ELSE
      IF (zlog) WRITE(*,'(a)') error_to_string(err,'',.true.)
    ENDIF
  END FUNCTION check_r1

  FUNCTION check_l1(err,var,def,wlog) RESULT(ret)
    !! Check an option value (logical).
    !! 
    !! The method checks an option value and optionally set a default value, __def__ to initialize 
    !! __var__ on error if given.
    TYPE(error), INTENT(in)       :: err  !! Error object from value getter.
    LOGICAL, INTENT(inout)        :: var  !! Input/output option value.
    LOGICAL, INTENT(in), OPTIONAL :: def  !! Default value to set.
    LOGICAL, INTENT(in), OPTIONAL :: wlog !! .true. to print warning/error message.
    TYPE(error) :: ret                    !! Input error.
    CHARACTER(len=*), PARAMETER :: defmsg = '... Using default value: ' 
    LOGICAL                     :: zlog
    ret = err
     zlog = .false. ; IF (PRESENT(wlog)) zlog = wlog
    IF (err == 0) RETURN
    IF (PRESENT(def)) THEN
      var = def
      IF (zlog) WRITE(*,'(a,a,a)') error_to_string(err,'',.true.),defmsg,to_string(var)
      ret = noerror
    ELSE
      IF (zlog) WRITE(*,'(a)') error_to_string(err,'',.true.)
    ENDIF
  END FUNCTION check_l1

  FUNCTION check_i1(err,var,def,wlog) RESULT(ret)
    !! Check an option value (integer).
    !! 
    !! The method checks an option value and optionally set a default value, __def__ to initialize 
    !! __var__ on error if given.
    TYPE(error), INTENT(in)       :: err  !! Error object from value getter.
    INTEGER, INTENT(inout)        :: var  !! Input/output option value.
    INTEGER, INTENT(in), OPTIONAL :: def  !! Default value to set.
    LOGICAL, INTENT(in), OPTIONAL :: wlog !! .true. to print warning/error message.
    TYPE(error) :: ret                    !! Input error.
    CHARACTER(len=*), PARAMETER :: defmsg = '... Using default value: '
    LOGICAL                     :: zlog
    ret = err
     zlog = .false. ; IF (PRESENT(wlog)) zlog = wlog
    IF (err == 0) RETURN
    IF (PRESENT(def)) THEN
      var = def
      IF (zlog) WRITE(*,'(a,a,a)') error_to_string(err,'',.true.),defmsg,to_string(var)
      ret = noerror
    ELSE
      IF (zlog) WRITE(*,'(a)') error_to_string(err,'',.true.)
    ENDIF
  END FUNCTION check_i1

  FUNCTION check_s1(err,var,def,wlog) RESULT(ret)
    !! Check an option value (string).
    !! 
    !! The method checks an option value and optionally set a default value, __def__ to initialize 
    !! __var__ on error if given.
    TYPE(error), INTENT(in)                      :: err  !! Error object from value getter.
    CHARACTER(len=*), INTENT(inout)              :: var  !! Input/output option value.
    CHARACTER(len=*), INTENT(in), OPTIONAL       :: def  !! Default value to set.
    LOGICAL, INTENT(in), OPTIONAL                :: wlog !! .true. to print warning/error message.
    TYPE(error) :: ret                                   !! Input error.
    CHARACTER(len=*), PARAMETER :: defmsg = '... Using default value: ' 
    LOGICAL                     :: zlog
    ret = err 
    zlog = .false. ; IF (PRESENT(wlog)) zlog = wlog
    IF (err == 0) RETURN
    IF (PRESENT(def)) THEN
      var = TRIM(def)
      IF (zlog) WRITE(*,'(a,a,a)') error_to_string(err,'',.true.),defmsg,var
      ret = noerror
    ELSE
      IF (zlog) WRITE(*,'(a)') error_to_string(err,'')
    ENDIF
    RETURN
  END FUNCTION check_s1

END MODULE MM_GLOBALS
