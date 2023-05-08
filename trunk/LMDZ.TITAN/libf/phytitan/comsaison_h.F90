
       module comsaison_h

       implicit none 

       integer isaison
       logical callsais
       real dist_star,declin,right_ascen

       real, allocatable, dimension(:) :: mu0,fract
!$OMP THREADPRIVATE(isaison,callsais,dist_star,declin,mu0,fract)

       end module comsaison_h
