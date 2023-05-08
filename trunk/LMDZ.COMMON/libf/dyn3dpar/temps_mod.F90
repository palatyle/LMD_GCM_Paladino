MODULE temps_mod

IMPLICIT NONE  

  INTEGER   itaufin ! total number of dynamical steps for the run
  INTEGER   itau_dyn
  INTEGER   itau_phy
  INTEGER   day_ini ! initial day # of simulation sequence
  INTEGER   day_end ! final day # ; i.e. day # when this simulation ends
  INTEGER   annee_ref
  INTEGER   day_ref
  REAL      dt ! (dynamics) time step (changes if doing Matsuno or LF step)
  REAL      jD_ref ! reference julian day date (beginning of experiment)
  REAL      jH_ref ! reference julian "hour" of reference julian date
  REAL      start_time
  CHARACTER (len=10) :: calend ! calendar type

  ! Additionnal Mars stuff:
  REAL hour_ini ! initial fraction of day of simulation sequence (0=<hour_ini<1)

!$OMP THREADPRIVATE(dt,jD_ref,jH_ref,start_time,hour_ini,                        &
!$OMP                day_ini,day_end,annee_ref,day_ref,itau_dyn,itau_phy,itaufin,&
!$OMP                calend)        

!WARNING: when adding a threadprivate variable in this module
!        do not forget to add it to the copyin clause when opening an OpenMP
!        parallel section. e.g. in gcm before call leapfrog_loc and/or
!        possibly in iniphysiq

END MODULE temps_mod
