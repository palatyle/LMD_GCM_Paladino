!
! part of the module_lmd_driver program
! that is (number of domains)-dependent
!


       CALL meso_physiq (ngrid,nlayer,nq,                              &
               firstcall,lastcall,wday_ini,pday,ptime,ptimestep,       &
               pplev,pplay,pphi,pu, pv, pt, pq, pw,                    &
               pdu, pdv, pdt, pdq, pdpsrf, tracerdyn,          &
               wtsurf,wtsoil,wemis,wq2,wqsurf,wco2ice,                 &
               wecri_phys_sec,                &
	       output_tab2d, output_tab3d)




