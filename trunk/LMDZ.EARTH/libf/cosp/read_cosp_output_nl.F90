!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE READ_COSP_OUTPUT_NL(cosp_nl,cfg)
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  USE mod_phys_lmdz_para
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(out) :: cfg
  ! Local variables
  integer :: i

  logical, save ::   Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,Lcfad_dbze94, &
             Lcfad_lidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp2,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lctpisccp,Ldbze94,Ltauisccp,Ltclisccp, &
             Llongitude,Llatitude,Lparasol_refl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfrac_out,Lbeta_mol532,Ltbrttov

  namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,Lcfad_dbze94, &
             Lcfad_lidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp2, &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lctpisccp,Ldbze94,Ltauisccp, &
             Ltclisccp,Llongitude,Llatitude,Lparasol_refl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfrac_out,Lbeta_mol532,Ltbrttov

  do i=1,N_OUT_LIST
    cfg%out_list(i)=''
  enddo
  
  IF (is_master) THEN
    open(10,file=cosp_nl,status='old')
    read(10,nml=cosp_output)
    close(10)
  ENDIF
  
  CALL bcast(Lradar_sim)
  CALL bcast(Llidar_sim)
  CALL bcast(Lisccp_sim)
  CALL bcast(Lmisr_sim)
  CALL bcast(Lrttov_sim)
  CALL bcast(Lalbisccp)
  CALL bcast(Latb532)
  CALL bcast(Lboxptopisccp)
  CALL bcast(Lboxtauisccp)
  CALL bcast(Lcfad_dbze94)
  CALL bcast(Lcfad_lidarsr532)
  CALL bcast(Lclcalipso2)
  CALL bcast(Lclcalipso)
  CALL bcast(Lclhcalipso)
  CALL bcast(Lclisccp2)
  CALL bcast(Lcllcalipso)
  CALL bcast(Lclmcalipso)
  CALL bcast(Lcltcalipso)
  CALL bcast(Lcltlidarradar)
  CALL bcast(Lctpisccp)
  CALL bcast(Ldbze94)
  CALL bcast(Ltauisccp)
  CALL bcast(Ltclisccp)
  CALL bcast(Llongitude)
  CALL bcast(Llatitude)
  CALL bcast(Lparasol_refl)
  CALL bcast(LclMISR)
  CALL bcast(Lmeantbisccp)
  CALL bcast(Lmeantbclrisccp)
  CALL bcast(Lfrac_out)
  CALL bcast(Lbeta_mol532)
  CALL bcast(Ltbrttov)
!$OMP BARRIER

!  print*,' Cles sorties cosp :'
!  print*,' Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim', &
!           Lradar_sim,Llidar_sim,Lisccp_sim,Lmisr_sim,Lrttov_sim

  ! Deal with dependencies
  if (.not.Lradar_sim) then
    Lcfad_dbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false.
    Ldbze94        = .false.
  endif
  if (.not.Llidar_sim) then
    Latb532 = .false.
    Lcfad_lidarsr532 = .false.
    Lclcalipso2      = .false.
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    Lcltlidarradar   = .false.
    Lparasol_refl    = .false.
    Lbeta_mol532     = .false.
  endif
  if (.not.Lisccp_sim) then
    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp2       = .false.
    Lctpisccp       = .false.
    Ltauisccp       = .false.
    Ltclisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.
  endif
  if (.not.Lmisr_sim) then
    LclMISR = .false.
  endif
  if (.not.Lrttov_sim) then
    Ltbrttov = .false.
  endif
  if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
      (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
    Lfrac_out = .false.
  endif

  ! Diagnostics that use Radar and Lidar
  if (((Lclcalipso2).or.(Lcltlidarradar)).and.((Lradar_sim).or.(Llidar_sim))) then
    Lclcalipso2    = .true.
    Lcltlidarradar = .true.
    Llidar_sim     = .true.
    Lradar_sim     = .true.
  endif

  cfg%Lstats = .false.
  if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.

  ! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim

  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
    cfg%Lwrite_output = .true.
  endif

  ! Output diagnostics
  i = 1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Latb532)          cfg%out_list(i) = 'atb532'
  i = i+1
  if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (Lcfad_dbze94)     cfg%out_list(i) = 'cfad_dbze94'
  i = i+1
  if (Lcfad_lidarsr532) cfg%out_list(i) = 'cfad_lidarsr532'
  i = i+1
  if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lclisccp2)        cfg%out_list(i) = 'clisccp2'
  i = i+1
  if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
  i = i+1
  if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lctpisccp)        cfg%out_list(i) = 'ctpisccp'
  i = i+1
  if (Ldbze94)          cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Ltclisccp)        cfg%out_list(i) = 'tclisccp'
  i = i+1
  if (Llongitude)       cfg%out_list(i) = 'lon'
  i = i+1
  if (Llatitude)        cfg%out_list(i) = 'lat'
  i = i+1
  if (Lparasol_refl)    cfg%out_list(i) = 'parasol_refl'
  i = i+1
  if (LclMISR)          cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (Lfrac_out)        cfg%out_list(i) = 'frac_out'
  i = i+1
  if (Lbeta_mol532)     cfg%out_list(i) = 'beta_mol532'
  i = i+1
  if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'

  if (i /= N_OUT_LIST) then
     print *, 'COSP_IO: wrong number of output diagnostics'
     stop
  endif

  ! Copy diagnostic flags to cfg structure
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lcfad_dbze94 = Lcfad_dbze94
  cfg%Lcfad_lidarsr532 = Lcfad_lidarsr532
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lclisccp2 = Lclisccp2
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lcltlidarradar = Lcltlidarradar
  cfg%Lctpisccp = Lctpisccp
  cfg%Ldbze94 = Ldbze94
  cfg%Ltauisccp = Ltauisccp
  cfg%Ltclisccp = Ltclisccp
  cfg%Llongitude = Llongitude
  cfg%Llatitude = Llatitude
  cfg%Lparasol_refl = Lparasol_refl
  cfg%LclMISR = LclMISR
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lfrac_out = Lfrac_out
  cfg%Lbeta_mol532 = Lbeta_mol532
  cfg%Ltbrttov = Ltbrttov

 END SUBROUTINE READ_COSP_OUTPUT_NL

