!-------------------------------------------------------------------
! INCLUDE microtab.h
!
	 REAL wco,df_GP
         INTEGER nz,nrad,imono,nztop,ntype

         parameter(nz=llm,ntype=5,nrad=10,imono=5,nztop=1)  !VERSION X
!        parameter(nz=llm,ntype=1,nrad=10,imono=5,nztop=1)  !VERSION X

         parameter(wco=177.,df_GP=2.)      !FOR FRACTAL PARTCICLES
!        parameter(wco=1.E+6,df_GP=3.)     !FOR SPHERE PARTICLES

      real rf(nrad),df(nrad),zf,aknc
      common/frac/rf,df,zf,aknc

!********************************************************************
! tcorrect, tx, microfi, cutoff: definis dans physiq.def (clesphys.h)
!------------
         ! WARNING: tx=production rate
         !          tcorrect is readjustment factor: =1 is continiuty
         !                                           =X is q()*X 
!------------
!(*1): si microfi=1, optcv et optci sont appeles a chaque appels de la 
!	physique  pour reactualiser les TAU's. De meme, pg2.F est 
!	active a chaque appel de la physique....
!      si microfi=0., optcv et optci, ainsi que pg2, ne sont appele qu'une
!       fois au debut, comme dans la version originale....
!------------
!       dans optci et optcv:
!      si cutoff=1, brume coupee facon Pascal -> T ok au sol et dans la strato
!                                             -> T tropopause mauvaise
!                                             -> albedo ok
!      si cutoff=2, brume coupee sous 100mbar -> T ok sol/tropopause/strato
!                                             -> mais albedo mauvais

