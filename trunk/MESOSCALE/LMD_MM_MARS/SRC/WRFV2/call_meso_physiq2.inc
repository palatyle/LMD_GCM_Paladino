!
! part of the module_lmd_driver program
! that is (number of domains)-dependent
!


SELECT CASE (id)

       CASE(1)
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

       CASE(2)
       CALL d2_physiq (ngrid,nlayer,nq,                           &
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

!       CASE(3:)
!       print *, '** Mars ** Not supported ...'
!       stop
!
END SELECT


