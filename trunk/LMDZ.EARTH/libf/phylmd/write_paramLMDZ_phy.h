c
      IF (is_sequential) THEN
      
      ndex2d = 0
      itau_w = itau_phy + itap
c
c Variables type caractere : plusieurs valeurs possibles
c
      IF(type_ocean.EQ.'force ') THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE IF(type_ocean.EQ.'slab  ') THEN
       zx_tmp_2d(1:iim,1:jjmp1)=2.
      ELSE IF(type_ocean.EQ.'couple') THEN
       zx_tmp_2d(1:iim,1:jjmp1)=3.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ocean",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(type_run.EQ.'CLIM'.OR.type_run.EQ.'ENSP') THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE IF(type_run.EQ.'AMIP'.OR.type_run.EQ.'CFMI') THEN
       zx_tmp_2d(1:iim,1:jjmp1)=2.
      ENDIF
      CALL histwrite(nid_ctesGCM,"type_run",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
c Variables logiques (1=true, 2=false)
c
      IF(ok_veget) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_veget",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_journe) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_journe",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_mensuel) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_mensuel",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_instan) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_instan",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_ade) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_ade",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_aie) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_aie",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c

c
c Champs 2D:
c
      zx_tmp_2d(1:iim,1:jjmp1)=bl95_b0
      CALL histwrite(nid_ctesGCM,"bl95_b0",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=bl95_b1
      CALL histwrite(nid_ctesGCM,"bl95_b1",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ip_ebil_phy
      CALL histwrite(nid_ctesGCM,"ip_ebil_phy",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=R_ecc
      CALL histwrite(nid_ctesGCM,"R_ecc",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=R_peri
      CALL histwrite(nid_ctesGCM,"R_peri",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=R_incl
      CALL histwrite(nid_ctesGCM,"R_incl",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=solaire
      CALL histwrite(nid_ctesGCM,"solaire",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=co2_ppm
      CALL histwrite(nid_ctesGCM,"co2_ppm",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=CH4_ppb
      CALL histwrite(nid_ctesGCM,"CH4_ppb",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=N2O_ppb
      CALL histwrite(nid_ctesGCM,"N2O_ppb",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=CFC11_ppt
      CALL histwrite(nid_ctesGCM,"CFC11_ppt",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=CFC12_ppt
      CALL histwrite(nid_ctesGCM,"CFC12_ppt",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=epmax
      CALL histwrite(nid_ctesGCM,"epmax",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! Mais est-il bien raisonable de stoker ces fichiers comme des
! champs 2D...
! WARNING :
! Il faudrait ici ajoute l'ecriture des champs
!      cycle_diurne = cycle_diurne_omp
!   soil_model = soil_model_omp
!   new_oliq = new_oliq_omp
!   ok_orodr = ok_orodr_omp
!   ok_orolf = ok_orolf_omp
!   ok_limitvrai = ok_limitvrai_omp
!   nbapp_rad = nbapp_rad_omp
!   iflag_con = iflag_con_omp
! qui se trouvaient auparavant dans gcm.def et maintenant dans 
! physiq.def.
! Mais regarder d'abord a quoi ca sert ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c
      IF(ok_adj_ema) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_adj_ema",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=iflag_clw
      CALL histwrite(nid_ctesGCM,"iflag_clw",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cld_lc_lsc
      CALL histwrite(nid_ctesGCM,"cld_lc_lsc",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cld_lc_con
      CALL histwrite(nid_ctesGCM,"cld_lc_con",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cld_tau_lsc
      CALL histwrite(nid_ctesGCM,"cld_tau_lsc",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cld_tau_con
      CALL histwrite(nid_ctesGCM,"cld_tau_con",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ffallv_lsc
      CALL histwrite(nid_ctesGCM,"ffallv_lsc",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ffallv_con
      CALL histwrite(nid_ctesGCM,"ffallv_con",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=coef_eva
      CALL histwrite(nid_ctesGCM,"coef_eva",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(reevap_ice) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"reevap_ice",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=iflag_cldcon
      CALL histwrite(nid_ctesGCM,"iflag_cldcon",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=iflag_pdf
      CALL histwrite(nid_ctesGCM,"iflag_pdf",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=fact_cldcon
      CALL histwrite(nid_ctesGCM,"fact_cldcon",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=facttemps
      CALL histwrite(nid_ctesGCM,"facttemps",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_newmicro) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_newmicro",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ratqsbas
      CALL histwrite(nid_ctesGCM,"ratqsbas",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ratqshaut
      CALL histwrite(nid_ctesGCM,"ratqshaut",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=rad_froid
      CALL histwrite(nid_ctesGCM,"rad_froid",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=rad_chau1
      CALL histwrite(nid_ctesGCM,"rad_chau1",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=rad_chau2
      CALL histwrite(nid_ctesGCM,"rad_chau2",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=top_height
      CALL histwrite(nid_ctesGCM,"top_height",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=overlap
      CALL histwrite(nid_ctesGCM,"overlap",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cdmmax
      CALL histwrite(nid_ctesGCM,"cdmmax",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=cdhmax
      CALL histwrite(nid_ctesGCM,"cdhmax",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ksta
      CALL histwrite(nid_ctesGCM,"ksta",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ksta_ter
      CALL histwrite(nid_ctesGCM,"ksta_ter",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_kzmin) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_kzmin",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=iflag_pbl
      CALL histwrite(nid_ctesGCM,"iflag_pbl",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=lev_histhf
      CALL histwrite(nid_ctesGCM,"lev_histhf",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=lev_histday
      CALL histwrite(nid_ctesGCM,"lev_histday",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=lev_histmth
      CALL histwrite(nid_ctesGCM,"lev_histmth",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      IF(ok_isccp) THEN
       zx_tmp_2d(1:iim,1:jjmp1)=1.
      ELSE
       zx_tmp_2d(1:iim,1:jjmp1)=0.
      ENDIF
      CALL histwrite(nid_ctesGCM,"ok_isccp",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=lonmin_ins
      CALL histwrite(nid_ctesGCM,"lonmin_ins",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=lonmax_ins
      CALL histwrite(nid_ctesGCM,"lonmax_ins",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=latmin_ins
      CALL histwrite(nid_ctesGCM,"latmin_ins",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=latmax_ins
      CALL histwrite(nid_ctesGCM,"latmax_ins",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_ins
      CALL histwrite(nid_ctesGCM,"ecrit_ins",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_hf
      CALL histwrite(nid_ctesGCM,"ecrit_hf",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_day
      CALL histwrite(nid_ctesGCM,"ecrit_day",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_mth
      CALL histwrite(nid_ctesGCM,"ecrit_mth",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_tra
      CALL histwrite(nid_ctesGCM,"ecrit_tra",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_reg
      CALL histwrite(nid_ctesGCM,"ecrit_reg",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=freq_ISCCP
      CALL histwrite(nid_ctesGCM,"freq_ISCCP",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      zx_tmp_2d(1:iim,1:jjmp1)=ecrit_ISCCP
      CALL histwrite(nid_ctesGCM,"ecrit_ISCCP",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
c=================================================================
c=================================================================
c=================================================================
c
      if (ok_sync) then
        call histsync(nid_ctesGCM)
      endif
c
      ENDIF ! mono_cpu
