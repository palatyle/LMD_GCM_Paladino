MODULE iniphysiq_mod

CONTAINS

subroutine iniphysiq(ngrid,nlayer,nq,piphysiq,&
                     punjours, pdayref, &
                     prad,pg,pr,pcpp,iflag_phys)

!use time_phylmdz_mod, only: iphysiq, day_step, ecritstart,ecritphy,&
!                            dtphys,daysec,day_ini,hour_ini
use comconst_mod, only : cpp
use cpdet_phy_mod, only: init_cpdet_phy
USE control_mod, ONLY: planet_type
!use surfdat_h, only: ini_surfdat_h
!use comsoil_h, only: ini_comsoil_h
!use dimradmars_mod, only: ini_dimradmars_mod
!use time_phylmdz_mod, only : init_time
use dimphy, only: init_dimphy
USE phys_state_var_mod
use module_model_constants, only : nu, TT00
USE variables_mod, only: phour_ini,zdt_split !! zdt_split <> ptimestep
   !real*8,intent(in) :: ptimestep !physics time step (s) [dtphys]

implicit none

REAL,intent(in) :: prad
REAL,intent(in) :: pg
REAL,intent(in) :: pr
REAL,intent(in) :: pcpp
REAL,intent(in) :: punjours
!DOUBLE PRECISION,intent(in) :: ptimestep

!real,intent(in) :: prad ! radius of the planet (m)
!real,intent(in) :: pg ! gravitational acceleration (m/s2)
!real,intent(in) :: pr ! ! reduced gas constant R/mu
!real,intent(in) :: pcpp ! specific heat Cp
!real,intent(in) :: punjours ! length (in s) of a standard day [daysec]
integer,intent(in) :: pdayref ! reference day of for the simulation [day_ini]
!real*8,intent(in) :: ptimestep !physics time step (s) [dtphys]
integer,intent(in) :: iflag_phys ! type of physics to be called

integer,intent(in) :: ngrid ! number of physics columns for this MPI process
integer,intent(in) :: nlayer ! number of atmospheric layers
integer,intent(in) :: nq ! number of tracers
!real,intent(in) :: phour_ini   ! start time (fraction of day) of the run 0=<phour_ini<1
real,intent(in) :: piphysiq   ! call physics every piphysiq dynamical timesteps
real*8 :: CPPP,nuvenus,T0venus

! copy some fundamental parameters to physics 
! and do some initializations 

!! Initialize dimphy module
!call init_dimphy(klon_omp,nlayer)
!...useless in mesoscale

!! initialize physical constants and arrays

      ! Initialize some "temporal and calendar" related variables
      !call init_time(day_ini,hour_ini,pdaysec,ptimestep)
      !call init_time(1,1,1,1,ptimestep)
      !=ptimestep 
      call init_dimphy(ngrid,nlayer)

      Call phys_state_var_init(nq)
      ! allocate "surfdat_h" arrays
      !call ini_surfdat_h(ngrid,nq)

     ! allocate "comgeomfi_h" arrays
      !call ini_comgeomfi_h(ngrid)

      ! allocate "comsoil_h" arrays
      !call ini_comsoil_h(ngrid)

      ! set some variables in "dimradmars_mod"
      !call ini_dimradmars_mod(ngrid,nlayer)
 
      ! set klon and klev for venus physiq
      !call init_dimphy(ngrid,nlayer)

      ! Initialize some physical constants
      call suphec(pcpp)

      cpp=pcpp
      planet_type="venus"
      CPPP=pcpp
      nuvenus=nu
      T0venus=TT00
      call init_cpdet_phy(CPPP,nuvenus,T0venus)

!! a few time constants initialization 
!! not done by init_time in phys_state_var_init
!! and supposed to be done in conf_phys (but not done in mesoscale)
!day_step=punjours/ptimestep
!iphysiq=piphysiq
!ecritstart=0 !! not used in MESOSCALE
!ecritphy=0 !! not used in MESOSCALE

!! check
!PRINT*,'Call to LMD physics:',day_step,' per Martian day'
!PRINT*,'Call to LMD physics: every',iphysiq,' dynamical timesteps'
!PRINT*,'Call to LMD physics: every',dtphys,' seconds'
!PRINT*,'Call to LMD physics: length of day is',daysec,' seconds'
!PRINT*,'Call to LMD physics: initial day is',day_ini
!PRINT*,'Call to LMD physics: initial fractional time is',hour_ini

! Initialize some "temporal and calendar" related variables
!CALL init_time(day_ini,hour_ini,punjours,ptimestep)
!...done in phys_state_var_init

end subroutine iniphysiq

END MODULE iniphysiq_mod
