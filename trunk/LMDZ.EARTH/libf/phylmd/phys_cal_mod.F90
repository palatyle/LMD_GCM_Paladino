! $Id:$
MODULE phys_cal_mod
! This module contains information on the calendar at the actual time step

  SAVE

  INTEGER :: year_cur      ! current year
  INTEGER :: mth_cur       ! current month
  INTEGER :: day_cur       ! current day
  INTEGER :: days_elapsed  ! number of whole days since start of the simulation 
  INTEGER :: mth_len       ! number of days in the current month
  REAL    :: hour
  REAL    :: jD_1jan
  REAL    :: jH_1jan
  REAL    :: xjour


CONTAINS
  
  SUBROUTINE phys_cal_update(jD_cur, jH_cur)
    ! This subroutine updates the module saved variables.

    USE IOIPSL
    
    REAL, INTENT(IN) :: jD_cur ! jour courant a l'appel de la physique (jour julien)
    REAL, INTENT(IN) :: jH_cur ! heure courante a l'appel de la physique (jour julien)
    
    CALL ju2ymds(jD_cur+jH_cur, year_cur, mth_cur, day_cur, hour)
    CALL ymds2ju(year_cur, 1, 1, 0., jD_1jan)
    
    jH_1jan = jD_1jan - int (jD_1jan)
    jD_1jan = int (jD_1jan) 
    xjour = jD_cur - jD_1jan
    days_elapsed = jD_cur - jD_1jan

    ! Get lenght of acutual month
    mth_len = ioget_mon_len(year_cur,mth_cur)

  END SUBROUTINE phys_cal_update

END MODULE phys_cal_mod
