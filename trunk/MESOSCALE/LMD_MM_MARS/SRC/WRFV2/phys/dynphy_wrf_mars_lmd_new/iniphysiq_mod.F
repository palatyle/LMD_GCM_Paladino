MODULE iniphysiq_mod

CONTAINS

subroutine iniphysiq(ngrid,nlayer,nq,piphysiq,&
                     punjours, pdayref, &
                     prad,pg,pr,pcpp,iflag_phys)

use phys_state_var_init_mod, only: phys_state_var_init
use time_phylmdz_mod, only: iphysiq, day_step, ecritstart,ecritphy,&
                            dtphys,daysec,day_ini,hour_ini
use update_inputs_physiq_mod, only: traceurs
USE variables_mod, only: phour_ini,zdt_split !! zdt_split <> pttimestep
   !real*8,intent(in) :: ptimestep !physics time step (s) [dtphys]

implicit none

real,intent(in) :: prad ! radius of the planet (m)
real,intent(in) :: pg ! gravitational acceleration (m/s2)
real,intent(in) :: pr ! ! reduced gas constant R/mu
real,intent(in) :: pcpp ! specific heat Cp
real,intent(in) :: punjours ! length (in s) of a standard day [daysec]
integer,intent(in) :: pdayref ! reference day of for the simulation [day_ini]
!real*8,intent(in) :: ptimestep !physics time step (s) [dtphys]
integer,intent(in) :: iflag_phys ! type of physics to be called

integer,intent(in) :: ngrid ! number of physics columns for this MPI process
integer,intent(in) :: nlayer ! number of atmospheric layers
integer,intent(in) :: nq ! number of tracers
!real,intent(in) :: phour_ini   ! start time (fraction of day) of the run 0=<phour_ini<1
real,intent(in) :: piphysiq   ! call physics every piphysiq dynamical timesteps

! copy some fundamental parameters to physics 
! and do some initializations 

!! Initialize dimphy module
!call init_dimphy(klon_omp,nlayer)
!...useless in mesoscale

!! initialize physical constants and arrays
call phys_state_var_init(ngrid,nlayer,nq, traceurs, &
                         pdayref,phour_ini,punjours,zdt_split, &
                         prad,pg,pr,pcpp)

!! read callphys.def
call conf_phys(ngrid,nlayer,nq)

!! a few time constants initialization 
!! not done by init_time in phys_state_var_init
!! and supposed to be done in conf_phys (but not done in mesoscale)
day_step=punjours/zdt_split
iphysiq=piphysiq
ecritstart=0 !! not used in MESOSCALE
ecritphy=0 !! not used in MESOSCALE

!! check
PRINT*,'Call to LMD physics:',day_step,' per Martian day'
PRINT*,'Call to LMD physics: every',iphysiq,' dynamical timesteps'
PRINT*,'Call to LMD physics: every',dtphys,' seconds'
PRINT*,'Call to LMD physics: length of day is',daysec,' seconds'
PRINT*,'Call to LMD physics: initial day is',day_ini
PRINT*,'Call to LMD physics: initial fractional time is',hour_ini

! Initialize some "temporal and calendar" related variables
!CALL init_time(day_ini,hour_ini,punjours,ptimestep)
!...done in phys_state_var_init

end subroutine iniphysiq

END MODULE iniphysiq_mod
