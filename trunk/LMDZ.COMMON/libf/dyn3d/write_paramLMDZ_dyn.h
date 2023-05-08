! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! Attention : il n'y a aucune raison pour ecrire ces constantes
! comme des champs 2D. A corriger un jour ...

c
      ndex2d = 0
      itau_w=itau_dyn+itau
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(prt_level) 
      CALL histwrite(nid_ctesGCM, "prt_level", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(dayref)
      CALL histwrite(nid_ctesGCM, "dayref", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(anneeref)
      CALL histwrite(nid_ctesGCM, "anneeref", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(raz_date)
      CALL histwrite(nid_ctesGCM, "raz_date", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(nday)
      CALL histwrite(nid_ctesGCM, "nday", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(day_step)
      CALL histwrite(nid_ctesGCM, "day_step", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iperiod)
      CALL histwrite(nid_ctesGCM, "iperiod", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iapp_tracvl)
      CALL histwrite(nid_ctesGCM, "iapp_tracvl", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iconser)
      CALL histwrite(nid_ctesGCM, "iconser", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iecri)
      CALL histwrite(nid_ctesGCM, "iecri", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=periodav
      CALL histwrite(nid_ctesGCM, "periodav", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(dissip_period)
      CALL histwrite(nid_ctesGCM, "dissip_period", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      IF(lstardis) THEN
       zx_tmp_2d(1:iip1,1:jjp1)=1.
      ELSE
       zx_tmp_2d(1:iip1,1:jjp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM, "lstardis", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(nitergdiv)
      CALL histwrite(nid_ctesGCM, "nitergdiv", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(nitergrot)
      CALL histwrite(nid_ctesGCM, "nitergrot", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(niterh) 
      CALL histwrite(nid_ctesGCM, "niterh", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=tetagdiv
      CALL histwrite(nid_ctesGCM, "tetagdiv", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=tetagrot
      CALL histwrite(nid_ctesGCM, "tetagrot", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=tetatemp
      CALL histwrite(nid_ctesGCM, "tetatemp", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=coefdis
      CALL histwrite(nid_ctesGCM, "coefdis", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      IF(purmats) THEN
       zx_tmp_2d(1:iip1,1:jjp1)=1.
      ELSE
       zx_tmp_2d(1:iip1,1:jjp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM, "purmats", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      IF(ok_guide) THEN
       zx_tmp_2d(1:iip1,1:jjp1)=1.
      ELSE
       zx_tmp_2d(1:iip1,1:jjp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM, "ok_guide", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      if (calend == 'earth_360d') then
        zx_tmp_2d(1:iip1,1:jjp1)=1.
      else if (calend == 'earth_365d') then
        zx_tmp_2d(1:iip1,1:jjp1)=2.
      else if (calend == 'earth_366d') then
        zx_tmp_2d(1:iip1,1:jjp1)=3.
      endif

      CALL histwrite(nid_ctesGCM, "true_calendar", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iflag_phys)
      CALL histwrite(nid_ctesGCM, "iflag_phys", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=REAL(iphysiq)
      CALL histwrite(nid_ctesGCM, "iphysiq", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/02
! La variable cycle_diurne n'est pas vue par la dynamique
!     IF(cycle_diurne) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "cycle_diurne", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     IF(soil_model) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "soil_model", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     IF(new_oliq) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "new_oliq", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     IF(ok_orodr) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "ok_orodr", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     IF(ok_orolf) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "ok_orolf", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     IF(ok_limitvrai) THEN
!      zx_tmp_2d(1:iip1,1:jjp1)=1.
!     ELSE
!      zx_tmp_2d(1:iip1,1:jjp1)=0.
!     ENDIF
!     CALL histwrite(nid_ctesGCM, "ok_limitvrai", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     zx_tmp_2d(1:iip1,1:jjp1)=nbapp_rad
!     CALL histwrite(nid_ctesGCM, "nbapp_rad", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!
!     zx_tmp_2d(1:iip1,1:jjp1)=iflag_con
!     CALL histwrite(nid_ctesGCM, "iflag_con", itau_w,
!    .               zx_tmp_2d,iip1*jjp1,ndex2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      zx_tmp_2d(1:iip1,1:jjp1)=clon
      CALL histwrite(nid_ctesGCM, "clon", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=clat
      CALL histwrite(nid_ctesGCM, "clat", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=grossismx
      CALL histwrite(nid_ctesGCM, "grossismx", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=grossismy
      CALL histwrite(nid_ctesGCM, "grossismy", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      IF(fxyhypb) THEN
       zx_tmp_2d(1:iip1,1:jjp1)=1.
      ELSE
       zx_tmp_2d(1:iip1,1:jjp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM, "fxyhypb", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=dzoomx
      CALL histwrite(nid_ctesGCM, "dzoomx", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=dzoomy
      CALL histwrite(nid_ctesGCM, "dzoomy", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=taux
      CALL histwrite(nid_ctesGCM, "taux", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=tauy
      CALL histwrite(nid_ctesGCM, "tauy", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      IF(ysinus) THEN
       zx_tmp_2d(1:iip1,1:jjp1)=1.
      ELSE
       zx_tmp_2d(1:iip1,1:jjp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM, "ysinus", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
      zx_tmp_2d(1:iip1,1:jjp1)=ip_ebil_dyn
      CALL histwrite(nid_ctesGCM, "ip_ebil_dyn", itau_w,
     .               zx_tmp_2d,iip1*jjp1,ndex2d)
c
c=================================================================
c
      if (ok_sync) then
        call histsync(nid_ctesGCM)
      endif
c
c=================================================================
