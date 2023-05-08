MODULE time_phylmdz_mod

    IMPLICIT NONE
    REAL,SAVE    :: pdtphys     ! physics time step (s)
!$OMP THREADPRIVATE(pdtphys)
    INTEGER,SAVE :: annee_ref   ! reference year from the origin
!$OMP THREADPRIVATE(annee_ref)
    INTEGER,SAVE :: day_ref     ! reference day of the origin
!$OMP THREADPRIVATE(day_ref)
    INTEGER,SAVE :: day_ini     ! initial day of the run since first day of annee_ref
!$OMP THREADPRIVATE(day_ini)
    INTEGER,SAVE :: day_end     ! final day of the run since first day of annee_ref
!$OMP THREADPRIVATE(day_end)
    INTEGER,SAVE :: raz_date    ! flag to reset date (0:no, 1:yes)
!$OMP THREADPRIVATE(raz_date)
    INTEGER,SAVE :: itau_phy     ! number of physics iterations
!$OMP THREADPRIVATE(itau_phy)

CONTAINS

  SUBROUTINE init_time(annee_ref_, day_ref_, day_ini_, day_end_, pdtphys_)
    USE ioipsl_getin_p_mod, ONLY : getin_p
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: annee_ref_
    INTEGER,INTENT(IN) :: day_ref_
    INTEGER,INTENT(IN) :: day_ini_
    INTEGER,INTENT(IN) :: day_end_
    REAL,INTENT(IN) :: pdtphys_
    
    annee_ref=annee_ref_
    day_ref=day_ref_
    day_ini=day_ini_
    day_end=day_end_
    pdtphys=pdtphys_

    ! Initialize module variable not inherited from dynamics
    raz_date = 0 ! default value
    CALL getin_p('raz_date', raz_date)
    
  END SUBROUTINE init_time

END MODULE time_phylmdz_mod      
