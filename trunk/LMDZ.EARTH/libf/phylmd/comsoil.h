!
! $Header$
!

      common /comsoil/inertie_sol,inertie_sno,inertie_ice
      real inertie_sol,inertie_sno,inertie_ice
!$OMP THREADPRIVATE(/comsoil/)