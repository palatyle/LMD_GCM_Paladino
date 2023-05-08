
module control_mod

  implicit none

  integer,save :: nday ! # of days to run
  integer,save :: day_step ! # of dynamical time steps per day
  integer,save :: iperiod  ! make a Matsuno step before avery iperiod-1 LF steps
  integer,save :: iconser !
  integer,save :: dissip_period ! apply dissipation every dissip_period dynamical step
  integer,save :: iphysiq ! call physics every iphysiq dynamical steps
  integer,save :: anneeref ! reference year # ! not used
  real,save :: periodav
!  integer,save :: ecritphy ! output data in "diagfi.nc" every ecritphy dynamical steps 
  character(len=10),save :: planet_type ! planet type ('earth','mars',...)
  character(len=4),save :: config_inca

end module control_mod
