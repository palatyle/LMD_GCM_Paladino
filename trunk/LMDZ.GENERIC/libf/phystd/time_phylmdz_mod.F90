MODULE time_phylmdz_mod

    IMPLICIT NONE
    REAL,SAVE    :: dtphys     ! physics time step (s)
!$OMP THREADPRIVATE(dtphys)
    INTEGER,SAVE :: day_step    ! number of dynamical steps per day
                                ! (set via inifis)
!$OMP THREADPRIVATE(day_step)
    INTEGER,SAVE :: nday       ! number of days to run
!$OMP THREADPRIVATE(nday)
    REAL,SAVE    :: daysec     ! length of day (s)
!$OMP THREADPRIVATE(daysec)
    INTEGER,SAVE :: day_ini     ! initial day of the run
!$OMP THREADPRIVATE(day_ini)

    INTEGER,SAVE :: ecritphy  ! for diagfi.nc outputs, write every ecritphy
                              ! dynamical steps (set via inifis)
!$OMP THREADPRIVATE(ecritphy)
    INTEGER,SAVE :: iphysiq   ! call physics every iphysiq dynamical step
                              ! (set via inifis)
!$OMP THREADPRIVATE(iphysiq)

CONTAINS

  SUBROUTINE init_time(day_ini_, daysec_, nday_, dtphys_)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: day_ini_
    REAL,INTENT(IN) :: daysec_
    INTEGER,INTENT(IN) :: nday_
    REAL,INTENT(IN) :: dtphys_
    
    day_ini=day_ini_
    daysec=daysec_
    nday=nday_
    dtphys=dtphys_

  END SUBROUTINE init_time

END MODULE time_phylmdz_mod      
