!
! $Id: nuage.h 1286 2009-12-17 13:47:10Z fairhead $
!
      REAL rad_froid, rad_chau1, rad_chau2, t_glace_max, t_glace_min

      common /nuagecom/ rad_froid,rad_chau1, rad_chau2,t_glace_max,     &
     &                  t_glace_min
!$OMP THREADPRIVATE(/nuagecom/)
