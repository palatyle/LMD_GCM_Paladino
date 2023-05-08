! Copyright 2017 Université de Reims Champagne-Ardenne 
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

!! file: mmp_gcm.f90
!! summary: YAMMS interfaces for the LMDZ GCM.
!! author: J. Burgalat
!! date: 2017
MODULE MMP_GCM
  !! Interface to YAMMS for the LMDZ GCM.
  USE MMP_GLOBALS
  USE MM_LIB
  USE CFGPARSE
  USE DATASETS
  IMPLICIT NONE

  CONTAINS 
    
  SUBROUTINE mmp_initialize(dt,p_prod,tx_prod,rc_prod,rplanet,g0, air_rad,air_mmol,clouds,cfgpath)
    !! Initialize global parameters of the model.
    !! 
    !! The function initializes all the global parameters of the model from direct input.
    !! Boolean (and Fiadero) parameters are optional as they are rather testing parameters. Their 
    !! default values are suitable for production runs.  
    !! @note
    !! If the subroutine fails to initialize parameters, the run is aborted.
    !!
    !! @warning
    !! If OpenMP is activated, this subroutine must be called in an $OMP SINGLE statement as it 
    !! initializes global variable that are not thread private.
    !!
    !! ```
    !! !$OMP SINGLE
    !! call mmp_initialize(...)
    !! !$OMP END SINGLE
    !! ```
    !!
    REAL(kind=mm_wp), INTENT(in)           :: dt
      !! Microphysics timestep in seconds.
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
    LOGICAL, INTENT(in)                    :: clouds
      !! Clouds microphysics control flag.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cfgpath
      !! Internal microphysic configuration file.

    INTEGER                                           :: coag_choice
    REAL(kind=mm_wp)                                  :: fiad_max, fiad_min,df,rm,rho_aer
    LOGICAL                                           :: w_h_prod, w_h_sed, w_h_coag, w_c_sed, w_c_nucond, &
                                                         no_fiadero, fwsed_m0, fwsed_m3
    TYPE(error)                                       :: err
    INTEGER                                           :: i
    TYPE(cfgparser)                                   :: cparser
    CHARACTER(len=st_slen)                            :: spcpath,pssfile,mqfile,opt_file
    CHARACTER(len=st_slen), DIMENSION(:), ALLOCATABLE :: species
    REAL(kind=mm_wp), DIMENSION(:), ALLOCATABLE       :: tmp

    w_h_prod    = .true.
    w_h_sed     = .true.
    w_h_coag    = .true.
    w_c_sed     = clouds
    w_c_nucond  = clouds
    fwsed_m0    = .true.
    fwsed_m3    = .false.
    no_fiadero  = .false.
    fiad_min    = 0.1_mm_wp
    fiad_max    = 10._mm_wp
    coag_choice = 7

    WRITE(*,'(a)') "##### MMP_GCM SPEAKING #####"
    WRITE(*,'(a)') "I will initialize ze microphysics model in moments YAMMS"
    WRITE(*,'(a)') "On error I will simply abort the program. Stay near your computer !"
    WRITE(*,*)
    WRITE(*,'(a)') "Reading muphys configuration file ("//trim(cfgpath)//")..."
    err = cfg_read_config(cparser,TRIM(cfgpath),.true.) 
    IF (err /= 0) THEN
      ! RETURN AN ERROR !!
      call abort_program(err)
    ENDIF
   
    ! YAMMS internal parameters:
    err = mm_check_opt(cfg_get_value(cparser,"rm",rm),rm,50e-9_mm_wp,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"df",df),df,2._mm_wp,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"rho_aer",rho_aer),rho_aer,1000._mm_wp,mm_log)
    ! the following parameters are primarily used to test and debug YAMMS.
    ! They are set in an optional configuration file and default to suitable values for production runs.
    err = mm_check_opt(cfg_get_value(cparser,"haze_production",w_h_prod)    ,w_h_prod   ,.true.   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"haze_sedimentation",w_h_sed)  ,w_h_sed    ,.true.   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"haze_coagulation",w_h_coag)   ,w_h_coag   ,.true.   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"clouds_sedimentation",w_c_sed),w_c_sed    ,clouds   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"clouds_nucl_cond",w_c_nucond) ,w_c_nucond ,clouds   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"wsed_m0",fwsed_m0)            ,fwsed_m0   ,.true.   ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"wsed_m3",fwsed_m3)            ,fwsed_m3   ,.false.  ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"no_fiadero",no_fiadero)       ,no_fiadero ,.false.  ,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"fiadero_min_ratio",fiad_min)  ,fiad_min   ,0.1_mm_wp,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"fiadero_max_ratio",fiad_max)  ,fiad_max   ,10._mm_wp,mm_log)
    err = mm_check_opt(cfg_get_value(cparser,"haze_coag_interactions",coag_choice),coag_choice,7,mm_log)

    ! optic look-up table file path.
    mmp_optic_file = ''
    opt_file = ''
    err = mm_check_opt(cfg_get_value(cparser,"optics/optic_file",opt_file),opt_file,'',mm_log)
    IF (err /= 0) THEN
      WRITE(*,'(a)') "Warning: I was unable to retrieve the path of the optic look-up table file:"
      WRITE(*,'(a)') "  The GCM may abort if it uses YAMMS optical properties calculation module !"
    ELSE
      mmp_optic_file = TRIM(opt_file)
    ENDIF

    ! Retrieve clouds species configuration file
    spcpath = ''
    IF (clouds) THEN
      err = mm_check_opt(cfg_get_value(cparser,"specie_cfg",spcpath), spcpath, wlog=mm_log)
      IF (err/=0) call abort_program(err)
    ENDIF

    ! YAMMS initialization.
    err = mm_global_init_0(dt,df,rm,rho_aer,p_prod,tx_prod,rc_prod,rplanet,g0, &
                           air_rad,air_mmol,coag_choice,clouds,spcpath,  &
                           w_h_prod,w_h_sed,w_h_coag,w_c_nucond,  &
                           w_c_sed,fwsed_m0,fwsed_m3, &
                           no_fiadero,fiad_min,fiad_max)
    IF (err /= 0) call abort_program(err)

    ! Extra initialization (needed for YAMMS method interfaces)
    err = mm_check_opt(cfg_get_value(cparser, "transfert_probability", mmp_w_ps2s), mmp_w_ps2s, wlog=mm_log)
    IF (err/=0) call abort_program(err)
    err = mm_check_opt(cfg_get_value(cparser, "electric_charging"    , mmp_w_qe  ), mmp_w_qe, wlog=mm_log)
    IF (err/=0) call abort_program(err)

    ! initialize transfert probabilities look-up tables
    IF (mm_w_haze_coag .AND. mmp_w_ps2s) THEN
      err = mm_check_opt(cfg_get_value(cparser, "ps2s_file", pssfile), pssfile)
      IF (err /= 0) call abort_program(err)

      IF (.NOT.read_dset(pssfile,'p_m0_co',mmp_pco0p)) THEN
        call abort_program(error("Cannot get 'p_m0_co' from "//pssfile,-1)) 
      ENDIF
      IF (.NOT.read_dset(pssfile,'p_m3_co',mmp_pco3p)) THEN
        call abort_program(error("Cannot get 'p_m3_co' from "//pssfile,-1)) 
      ENDIF
      IF (.NOT.read_dset(pssfile,'p_m0_fm',mmp_pfm0p)) THEN
        call abort_program(error("Cannot get 'p_m0_fm' from "//pssfile,-1))
      ENDIF
      IF (.NOT.read_dset(pssfile,'p_m3_fm',mmp_pfm3p)) THEN
        call abort_program(error("Cannot get 'p_m3_fm' from "//pssfile,-1))
      ENDIF
    ENDIF
    ! initialize mean electric correction look-up tables
    IF (mm_w_haze_coag .AND. mmp_w_qe) THEN
      err = mm_check_opt(cfg_get_value(cparser, "mq_file", mqfile), mqfile)
      IF (err /= 0) call abort_program(err)

      IF (.NOT.read_dset(mqfile,'qbsf0',mmp_qbsf0)) THEN
        call abort_program(error("Cannot get 'qbsf0' from "//mqfile,-1))
      ELSE
        mmp_qbsf0_e(1,1) = MINVAL(mmp_qbsf0%x) 
        mmp_qbsf0_e(1,2) = MAXVAL(mmp_qbsf0%x)
        mmp_qbsf0_e(2,1) = MINVAL(mmp_qbsf0%y)
        mmp_qbsf0_e(2,2) = MAXVAL(mmp_qbsf0%y)
      ENDIF
      IF (.NOT.read_dset(mqfile,'qbsf3',mmp_qbsf3)) THEN
        call abort_program(error("Cannot get 'qbsf3' from "//mqfile,-1))
      ELSE
        mmp_qbsf3_e(1,1) = MINVAL(mmp_qbsf3%x) 
        mmp_qbsf3_e(1,2) = MAXVAL(mmp_qbsf3%x)
        mmp_qbsf3_e(2,1) = MINVAL(mmp_qbsf3%y)
        mmp_qbsf3_e(2,2) = MAXVAL(mmp_qbsf3%y)
      ENDIF
      IF (.NOT.read_dset(mqfile,'qbff0',mmp_qbff0)) THEN
        call abort_program(error("Cannot get 'qbff0' from "//mqfile,-1))
      ELSE
        mmp_qbff0_e(1,1) = MINVAL(mmp_qbff0%x) 
        mmp_qbff0_e(1,2) = MAXVAL(mmp_qbff0%x)
        mmp_qbff0_e(2,1) = MINVAL(mmp_qbff0%y)
        mmp_qbff0_e(2,2) = MAXVAL(mmp_qbff0%y)
      ENDIF
    ENDIF
    ! spherical mode inter-moments function parameters
    IF (.NOT.cfg_has_section(cparser,'alpha_s')) call abort_program(error("Cannot find [alpha_s] section",-1))
    err = read_aprm(cparser,'alpha_s',mmp_asp)
    IF (err /= 0) call abort_program(error("alpha_s: "//TRIM(err%msg),-1))
    ! fractal mode inter-moments function parameters
    IF (.NOT.cfg_has_section(cparser,'alpha_f')) call abort_program(error("Cannot find [alpha_f] section",-1))
    err = read_aprm(cparser,'alpha_f',mmp_afp)
    IF (err /= 0) call abort_program(error("alpha_s: "//TRIM(err%msg),-1))

    ! get size-distribution laws parameters
    IF (.NOT.cfg_has_section(cparser,'dndr_s')) call abort_program(error("Cannot find [dndr_s] section",-2))
    err = read_nprm(cparser,'dndr_s',mmp_pns)
    IF (err /= 0) call abort_program(error("dndr_s: "//TRIM(err%msg),-2))
    IF (.NOT.cfg_has_section(cparser,'dndr_f')) call abort_program(error("Cannot find [dndr_f] section",-2))
    err = read_nprm(cparser,'dndr_f',mmp_pnf)
    IF (err /= 0) call abort_program(error("dndr_f: "//TRIM(err%msg),-2))

    ! btk coefficients
    IF (.NOT.cfg_has_section(cparser,'btks')) call abort_program(error("Cannot find [btks] section",-1))
    err = cfg_get_value(cparser,"btks/bt0",tmp) ; IF (err/=0) call abort_program(error("bt0: "//TRIM(err%msg),-1))
    IF (SIZE(tmp) /= 5)  call abort_program(error("bt0: Inconsistent number of coefficients",-1))
    mmp_bt0 = tmp
    err = cfg_get_value(cparser,"btks/bt3",tmp) ; IF (err/=0) call abort_program(error("bt3: "//TRIM(err%msg),-1))
    IF (SIZE(tmp) /= 5)  call abort_program(error("bt3: Inconsistent number of coefficients",-1))
    mmp_bt3 = tmp

    ! dump parameters ...
    WRITE(*,'(a)')        "========= MUPHYS PARAMETERS ==========="
    WRITE(*,'(a,L2)')     "transfert_probability: ", mmp_w_ps2s
    WRITE(*,'(a,L2)')     "electric_charging    : ", mmp_w_qe
    call mm_dump_parameters()
      
  END SUBROUTINE mmp_initialize

  FUNCTION read_aprm(parser,sec,pp) RESULT(err)
    !! Read and store [[mmp_gcm(module):aprm(type)]] parameters. 
    TYPE(cfgparser), INTENT(in)  :: parser !! Configuration parser 
    CHARACTER(len=*), INTENT(in) :: sec    !! Name of the section that contains the parameters.
    TYPE(aprm), INTENT(out)      :: pp     !! [[mmp_gcm(module):aprm(type)]] object that stores the parameters values.
    TYPE(error) :: err                     !! Error status of the function.
    err = cfg_get_value(parser,TRIM(sec)//'/a',pp%a) ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/b',pp%b) ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/c',pp%c) ; IF (err /= 0) RETURN
    IF (SIZE(pp%a) /= SIZE(pp%b) .OR. SIZE(pp%a) /= SIZE(pp%c)) &
    err = error("Inconsistent number of coefficients (a,b, and c must have the same size)",-1)
    RETURN
  END FUNCTION read_aprm

  FUNCTION read_nprm(parser,sec,pp) RESULT(err)
    !! Read and store [[mmp_gcm(module):nprm(type)]] parameters.
    TYPE(cfgparser), INTENT(in)  :: parser !! Configuration parser 
    CHARACTER(len=*), INTENT(in) :: sec    !! Name of the section that contains the parameters.
    TYPE(nprm), INTENT(out)      :: pp     !! [[mmp_gcm(module):nprm(type)]] object that stores the parameters values.
    TYPE(error) :: err                     !! Error status of the function.
    err = cfg_get_value(parser,TRIM(sec)//'/rc',pp%rc) ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/a0',pp%a0) ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/c',pp%c)   ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/a',pp%a)   ; IF (err /= 0) RETURN
    err = cfg_get_value(parser,TRIM(sec)//'/b',pp%b)   ; IF (err /= 0) RETURN
    IF (SIZE(pp%a) /= SIZE(pp%b)) &
      err = error("Inconsistent number of coefficients (a and b must have the same size)",3)
    RETURN
  END FUNCTION read_nprm

END MODULE MMP_GCM

