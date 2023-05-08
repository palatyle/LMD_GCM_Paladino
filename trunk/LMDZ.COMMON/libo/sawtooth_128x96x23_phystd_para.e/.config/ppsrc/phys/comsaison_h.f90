











       module comsaison_h

       implicit none 

!       integer,save :: isaison
!       logical,save :: callsais
!!$OMP THREADPRIVATE(isaison,callsais)

       real,save :: dist_star,declin,right_ascen
!$OMP THREADPRIVATE(dist_star,declin,right_ascen)

       real, allocatable, dimension(:) :: mu0,fract
!$OMP THREADPRIVATE(mu0,fract)

       end module comsaison_h
