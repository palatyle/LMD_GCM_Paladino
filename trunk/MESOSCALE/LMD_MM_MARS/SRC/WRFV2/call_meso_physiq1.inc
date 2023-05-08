!
! part of the module_lmd_driver program
! that is (number of domains)-dependent
!


       CALL physiq (ngrid,nlayer,nq,                              &
               firstcall,lastcall,pday,ptime,ptimestep,                &
               pplev,pplay,pphi,pu, pv, pt, pq, pw,                    &
               pdu, pdv, pdt, pdq, pdpsrf, tracerdyn,          &
               wtsurf,wtsoil,wemis,wq2,wqsurf,wco2ice,                 &
#ifdef NEWPHYS
               wtnom, wisoil, wdsoil, &
               wfluxrad,wwstar,               &
#endif
               wday_ini,                      &
	       output_tab2d, output_tab3d,    &
               flag_LES)




