











       module comgeomfi_h

       implicit none

       REAL,SAVE :: totarea ! total surface (m2) for this local (MPI/OpenMP) domain
       REAL,SAVE :: totarea_planet ! total planetary surface (m2)
!$OMP THREADPRIVATE(totarea,totarea_planet)

       end module comgeomfi_h

