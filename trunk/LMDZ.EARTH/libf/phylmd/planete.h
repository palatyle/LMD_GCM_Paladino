c-----------------------------------------------------------------------
c INCLUDE planet.h

      COMMON/planet/aphelie,periheli,year_day,peri_day,                 &
     &       obliquit,                                                  &
     &       timeperi,e_elips,p_elips,unitastr

      REAL aphelie,periheli,year_day,peri_day,                          &
     &     obliquit,                                                    &
     &       timeperi,e_elips,p_elips,unitastr

c-----------------------------------------------------------------------
!$OMP THREADPRIVATE(/planet/)
