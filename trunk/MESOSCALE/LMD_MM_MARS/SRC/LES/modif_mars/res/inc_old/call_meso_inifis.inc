!
! part of the module_lmd_driver program
! that is (number of domains)-dependent
!


       CALL meso_inifis(ngridmx,nlayer,nqmx,dt, &
               wday_ini,wdaysec,               &
               wappel_phys,wecri_phys,         &
               lat_vec,lon_vec,aire_vec,1/reradius,g,r_d,cp,           &
               womeg,wmugaz,                                           &
               wyear_day,wperiheli,waphelie,wperi_day,wobliquit,       &
               wz0,wemin_turb,wlmixmin,                                &
               wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS,   &
               wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS,        &
               walbedodat,winertiedat,wphisfi,                         &
               wzmea,wzstd,wzsig,wzgam,wzthe,                          &
               wtheta,wpsi)

