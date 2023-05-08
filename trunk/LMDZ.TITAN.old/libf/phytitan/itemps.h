!--------------------------------------------------------------------c
!      itemps.h : variable de temps d'executions.                    c
!--------------------------------------------------------------------c
       real ttphys           ! tps exec de la physique 
       real ttmuphys         ! tps exec de la microphysique
       real tthaze           ! tps exec de la brume
       real ttcclds          ! tps exec CONDENSATION nuages
       real ttsclds          ! tps exec SEDIMENTATION nuages
       real ttdynt           ! tps exec dyn tout confondu (comprend physique)
       real ttphytra         ! tps exec phytrac
       real ttrad            ! tps exec TR
       real ttadvtr          ! tps exec advection des traceurs dans la dynamique

       common /itimes/ttdynt,ttphys,ttphytra,ttrad,                     &
     &                ttmuphys,tthaze,ttcclds,ttsclds,ttadvtr


! ---- variables locales a chaque routines utilisant itemps.h
       real tt0        ! tps de demarage de la routine 
       real ttt0       ! tps de demarage de la routine 
       real tttt0       ! tps de demarage de la routine 
       real tt1        ! tps d'EXECUTION de la routine.
       real ttt1       ! tps d'EXECUTION de la routine.
       real tttt1       ! tps d'EXECUTION de la routine.
   
