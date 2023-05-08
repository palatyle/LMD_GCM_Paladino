!
! part of the module_lmd_driver program
! that is (number of domains)-dependent
!


SELECT CASE (id)

       CASE(1)
       CALL inifis(ngrid,nlayer, &
#ifdef NEWPHYS
               nq,  &
#endif
               wday_ini,wdaysec,                &
               wappel_phys,                     &
               lat_vec,lon_vec,aire_vec,1/reradius,g,r_d,cp,           &
#ifdef NEWPHYS
               dt,                                                     &
#else
               nq,dt,                                                &
#endif
               womeg,wmugaz,                                           &
               wyear_day,wperiheli,waphelie,wperi_day,wobliquit,       &
               wz0,wemin_turb,wlmixmin,                                &
               wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS,   &
               wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS,        &
#ifdef NEWPHYS
               walbedodat, wphisfi, wvolcapa, wz0tab,                  &
#else
               walbedodat, winertiedat, wphisfi,                       &
#endif
               wzmea,wzstd,wzsig,wzgam,wzthe,                          &
               wtheta, wpsi)

       CASE(2)
       CALL d2_inifis(ngrid,nlayer, &
#ifdef NEWPHYS
               nq,  &
#endif
               wday_ini,wdaysec,                   &
               wappel_phys,                        &
               lat_vec,lon_vec,aire_vec,1/reradius,g,r_d,cp,           &
#ifdef NEWPHYS
               dt,                                                     &
#else
               nq,dt,                                                &
#endif
               womeg,wmugaz,                                           &
               wyear_day,wperiheli,waphelie,wperi_day,wobliquit,       &
               wz0,wemin_turb,wlmixmin,                                &
               wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS,   &
               wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS,        &
#ifdef NEWPHYS
               walbedodat, wphisfi, wvolcapa, wz0tab,                  &
#else
               walbedodat, winertiedat, wphisfi,                       &
#endif
               wzmea,wzstd,wzsig,wzgam,wzthe,                          &
               wtheta, wpsi)

       CASE(3)
       CALL d3_inifis(ngrid,nlayer, &
#ifdef NEWPHYS
               nq,  &
#endif
               wday_ini,wdaysec,                   &
               wappel_phys,                        &
               lat_vec,lon_vec,aire_vec,1/reradius,g,r_d,cp,           &
#ifdef NEWPHYS
               dt,                                                     &
#else
               nq,dt,                                                &
#endif
               womeg,wmugaz,                                           &
               wyear_day,wperiheli,waphelie,wperi_day,wobliquit,       &
               wz0,wemin_turb,wlmixmin,                                &
               wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS,   &
               wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS,        &
#ifdef NEWPHYS
               walbedodat, wphisfi, wvolcapa, wz0tab,                  &
#else
               walbedodat, winertiedat, wphisfi,                       &
#endif
               wzmea,wzstd,wzsig,wzgam,wzthe,                          &
               wtheta, wpsi)


!
!       CASE(3:)
!       print *, '** Mars ** Not supported ...'
!       stop
!
END SELECT
