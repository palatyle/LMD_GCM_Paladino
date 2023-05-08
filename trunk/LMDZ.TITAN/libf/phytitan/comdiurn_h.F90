
       module comdiurn_h

       implicit none

       real,allocatable,dimension(:),save :: sinlon, coslon, sinlat, coslat
!$OMP THREADPRIVATE(sinlon,coslon,sinlat,coslat)

       end module comdiurn_h

