MODULE timer_filtre
IMPLICIT NONE
  PRIVATE
  REAL :: time
  REAL :: Last_time
  PUBLIC :: Init_timer, start_timer, stop_timer, Print_filtre_timer
CONTAINS

 SUBROUTINE Init_timer
   time=0
   Last_time=0
 END SUBROUTINE Init_timer
 
 SUBROUTINE Start_timer
  
   CALL cpu_time(last_time)

 END SUBROUTINE start_timer
 
 
 SUBROUTINE stop_timer
   REAL :: T 
   
   CALL cpu_time(t)
   Time=Time+t-last_time
 
  END SUBROUTINE stop_timer
  
  SUBROUTINE Print_filtre_timer
  PRINT *,"Temps CPU passe dans le filtre :",Time
  END SUBROUTINE  Print_filtre_timer

END MODULE timer_filtre
