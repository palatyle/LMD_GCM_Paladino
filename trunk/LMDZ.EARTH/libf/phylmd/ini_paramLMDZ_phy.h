cym    Non implemente en mode parallele

       IF (is_sequential) THEN  
c
       zstophy = dtime
       zout = ecrit_day
c
       idayref = day_ref
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c
       CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlon,zx_lon)
       if (iim.gt.1) then
       DO i = 1, iim
         zx_lon(i,1) = rlon(i+1)
         zx_lon(i,jjmp1) = rlon(i+1)
       ENDDO
       endif
       CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)
c
       CALL histbeg("paramLMDZ_phy.nc", 
     .                 iim,zx_lon(:,1), jjmp1,zx_lat(1,:),
     .                 1,1,1,1,
     .                 itau_phy, zjulian, dtime,
     .                 nhori, nid_ctesGCM)
c
c Variables type caractere : plusieurs valeurs possibles
c
       CALL histdef(nid_ctesGCM, "ocean", 
     .        "Type ocean utilise: 1=force, 2=slab, 3=couple",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "type_run",
     .        "Type run: 1= CLIM ou ENSP, 2= AMIP ou CFMI",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
c Variables logiques (1=true, 0=false)
c
       CALL histdef(nid_ctesGCM, "ok_veget", 
     .        "Type de modele de vegetation: 1=ORCHIDEE, 0=bucket",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_journe", 
     .        "Creation du fichier histday: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_mensuel", 
     .        "Creation du fichier histmth: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_instan", 
     .        "Creation du fichier histins: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_ade", 
     .        "Aerosol direct effect: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_aie", 
     .        "Aerosol indirect effect: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "bl95_b0", 
     .        "Parameter in CDNC-maer link (Boucher&Lohmann 1995)",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "bl95_b1", 
     .        "Parameter in CDNC-maer link (Boucher&Lohmann 1995)",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ip_ebil_phy", 
     .                "Niveau sortie diags bilan energie cote physique",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "R_ecc", 
     .                "Excentricite","-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "R_peri", 
     .                "Equinoxe","-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "R_incl", 
     .                "Inclinaison","deg",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "solaire", 
     .                "Constante solaire","W/m2",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "co2_ppm", 
     .                "Concentration du CO2", "ppm",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "CH4_ppb", 
     .                "Concentration du CH4", "ppb",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "N2O_ppb",
     .                "Concentration du N2O", "ppb",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "CFC11_ppt",
     .                "Concentration du CFC11", "ppt",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "CFC12_ppt",
     .                "Concentration du CFC12", "ppt",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "epmax",
     .                "Efficacite precip", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_adj_ema",
     .                "ok_adj_ema: 1=true, 0=false", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "iflag_clw",
     .                "iflag_clw", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cld_lc_lsc",
     .                "cld_lc_lsc", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cld_lc_con",
     .                "cld_lc_con", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cld_tau_lsc",
     .                "cld_tau_lsc", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cld_tau_con",
     .                "cld_tau_con", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ffallv_lsc",
     .                "ffallv_lsc", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ffallv_con",
     .                "ffallv_con", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "coef_eva",
     .                "coef_eva", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "reevap_ice",
     .                "reevap_ice: 1=true, 0=false", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "iflag_cldcon",
     .                "iflag_cldcon", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "iflag_pdf",
     .                "iflag_pdf", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "fact_cldcon",
     .                "fact_cldcon", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "facttemps",
     .                "facttemps", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_newmicro",
     .                "Nouvelle micro-physique: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ratqsbas",
     .                "ratqsbas", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ratqshaut",
     .                "ratqshaut", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "rad_froid",
     .                "rad_froid", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "rad_chau1",
     .                "rad_chau1", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "rad_chau2",
     .                "rad_chau2", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "top_height",
     .                "top_height", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "overlap",
     .                "overlap", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cdmmax",
     .                "cdmmax", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "cdhmax",
     .                "cdhmax", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ksta",
     .                "ksta", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ksta_ter",
     .                "ksta_ter", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_kzmin",
     .                "ok_kzmin: 1=true, 0=false", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "iflag_pbl",
     .                "iflag_pbl", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "lev_histhf",
     .                "lev_histhf", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "lev_histday",
     .                "lev_histday", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "lev_histmth",
     .                "lev_histmth", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ok_isccp",
     .                "Creation fichier histISCCP: 1=true, 0=false",
     .                "-",iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "lonmin_ins",
     .                "lonmin_ins", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "lonmax_ins",
     .                "lonmax_ins", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "latmin_ins",
     .                "latmin_ins", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "latmax_ins",
     .                "latmax_ins", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_ins",
     .                "ecrit_ins", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_hf",
     .                "ecrit_hf", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_day",
     .                "ecrit_day", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_mth",
     .                "ecrit_mth", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_tra",
     .                "ecrit_tra", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_reg",
     .                "ecrit_reg", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "freq_ISCCP",
     .                "freq_ISCCP", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
       CALL histdef(nid_ctesGCM, "ecrit_ISCCP",
     .                "ecrit_ISCCP", "-",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "once", zstophy,zout)
c
c=================================================================
c
       CALL histend(nid_ctesGCM)
       
       ENDIF ! is_sequential
c
c=================================================================
