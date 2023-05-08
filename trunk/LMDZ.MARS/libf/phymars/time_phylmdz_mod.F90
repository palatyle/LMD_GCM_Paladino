MODULE time_phylmdz_mod

    IMPLICIT NONE
    REAL,SAVE    :: dtphys      ! physics time step (s)
!$OMP THREADPRIVATE(dtphys)
    INTEGER,SAVE :: day_step    ! number of dynamical steps per day
                                ! (set via conf_phys)
!$OMP THREADPRIVATE(day_step)
    REAL,SAVE    :: daysec     ! length of day (s)
!$OMP THREADPRIVATE(daysec)
    INTEGER,SAVE :: day_ini     ! initial day of the run
!$OMP THREADPRIVATE(day_ini)
    REAL,SAVE :: hour_ini       ! start time (fraction of day) of the run
                                ! 0=<hour_ini<1
!$OMP THREADPRIVATE(hour_ini)

    INTEGER,SAVE :: ecritphy    ! for diagfi.nc outputs, write every ecritphy
                                ! dynamical steps (set via conf_phys)
!$OMP THREADPRIVATE(ecritphy)
    INTEGER,SAVE :: iphysiq   ! call physics every iphysiq dynamical step
                              ! (set via conf_phys)
!$OMP THREADPRIVATE(iphysiq)
    INTEGER,SAVE :: ecritstart ! write a restart state every ecritstart
                               ! dynamical steps (set via conf_phys)
!$OMP THREADPRIVATE(ecritstart)
CONTAINS

  SUBROUTINE init_time(day_ini_, hour_ini_, daysec_, dtphys_)
    USE ioipsl_getin_p_mod, ONLY : getin_p
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: day_ini_
    REAL,INTENT(IN) :: hour_ini_
    REAL,INTENT(IN) :: daysec_
    REAL,INTENT(IN) :: dtphys_
    
    day_ini=day_ini_
    hour_ini=hour_ini_
    daysec=daysec_
    dtphys=dtphys_
    
  END SUBROUTINE init_time

END MODULE time_phylmdz_mod      
