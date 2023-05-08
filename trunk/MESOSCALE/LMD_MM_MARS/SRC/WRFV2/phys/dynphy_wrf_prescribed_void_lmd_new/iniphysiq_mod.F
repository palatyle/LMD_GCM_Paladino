MODULE iniphysiq_mod

CONTAINS

subroutine iniphysiq(ngrid,nlayer,nq,piphysiq,&
                     punjours, pdayref, &
                     prad,pg,pr,pcpp,iflag_phys)

implicit none

real,intent(in) :: prad ! radius of the planet (m)
real,intent(in) :: pg ! gravitational acceleration (m/s2)
real,intent(in) :: pr ! ! reduced gas constant R/mu
real,intent(in) :: pcpp ! specific heat Cp
real,intent(in) :: punjours ! length (in s) of a standard day [daysec]
integer,intent(in) :: pdayref ! reference day of for the simulation [day_ini]
integer,intent(in) :: iflag_phys ! type of physics to be called

integer,intent(in) :: ngrid ! number of physics columns for this MPI process
integer,intent(in) :: nlayer ! number of atmospheric layers
integer,intent(in) :: nq ! number of tracers
real,intent(in) :: piphysiq   ! call physics every piphysiq dynamical timesteps

end subroutine iniphysiq

END MODULE iniphysiq_mod
