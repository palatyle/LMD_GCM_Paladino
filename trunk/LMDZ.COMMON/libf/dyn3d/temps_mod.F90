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
      real hour_ini ! initial fraction of day of simulation sequence (0=<hour_ini<1)

END MODULE temps_mod
