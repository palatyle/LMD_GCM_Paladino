










!
! $Id $
!

MODULE control_mod

! LF 01/2010
! Remplacement du fichier et common control

  IMPLICIT NONE

  real,save :: periodav
  real,save :: starttime
  integer,save :: nday
  integer,save :: day_step ! # of dynamical time steps per day
  integer,save :: iperiod ! make a Matsuno step before avery iperiod-1 LF steps
  integer,save :: iapp_tracvl ! apply (cumulated) traceur advection every
                              ! iapp_tracvl dynamical steps
  integer,save :: nsplit_phys ! number of sub-cycle steps in call to physics
  integer,save :: iconser
  integer,save :: iecri
  integer,save :: dissip_period ! apply dissipation every dissip_period
                                ! dynamical step
  integer,save :: iphysiq ! call physics every iphysiq dynamical steps
  integer,save :: iecrimoy
  integer,save :: dayref
  integer,save :: anneeref ! reference year # 
  integer,save :: raz_date
  integer,save :: ip_ebil_dyn
  logical,save :: offline
  logical,save :: cpofT
  logical,save :: force_conserv_tracer ! enforce conservation of tracer mass
  character(len=4),save :: config_inca
  character(len=10),save :: planet_type ! planet type ('earth','mars',...)
  logical,save :: output_grads_dyn ! output dynamics diagnostics in
                                   ! binary grads file 'dyn.dat' (y/n)
  logical,save :: ok_dynzon  ! output zonal transports in dynzon.nc file
  logical,save :: ok_dyn_ins ! output instantaneous values of fields
                             ! in the dynamics in NetCDF files dyn_hist*nc
  logical,save :: ok_dyn_ave ! output averaged values of fields in the dynamics
                             ! in NetCDF files dyn_hist*ave.nc
  logical,save :: resetvarc  ! allows to reset the variables in sortvarc
  logical,save :: less1day   ! allows to run less than 1 day (for Venus)
  real,save :: fractday   ! fraction of the day to run in this case
  
  integer,save :: ndynstep ! Alternative to using less1day&fractday; user may
                           ! specify number of dynamical steps to run

!  integer,save :: ecritphy ! (Mars/generic) output (writediagfi) every
                           ! ecritphy dynamical steps
  integer,save :: ecritstart ! (Mars) output data in "start.nc" every
                             !ecritstart dynamical steps 
  real,save :: timestart ! (Mars) time start for run in "start.nc"

  ! stuff for compatibility with Mars/Generic old dyn cores. To be cleaned!
  integer,save :: idissip ! (Mars/old dyn) dissipation freq.
  real,save :: nday_r ! (Mars/old dyn) number of days to run (possibly including a fraction of day)


END MODULE
