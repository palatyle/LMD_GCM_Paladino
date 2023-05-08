!
! $Header$
!
c$OMP MASTER
c
      IF (ok_histNMC(3)) THEN
c
       zstophy = dtime
       zstohf = ecrit_hf
       zstomth = ecrit_mth
c      zout = 6 * 3600.
       zout = freq_outNMC(3)
c
         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c
cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlon,zx_lon)
cym         DO i = 1, iim
cym            zx_lon(i,1) = rlon(i+1)
cym            zx_lon(i,jjmp1) = rlon(i+1)
cym         ENDDO
         DO ll=1,klev
            znivsig(ll)=float(ll)
         ENDDO
cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)
cym         CALL histbeg("histNMC.nc", iim,zx_lon(:,1), jjmp1,zx_lat(1,:),
cym     .                 1,iim,1,jjmp1, itau_phy, zjulian, dtime, 
cym     .                 nhori, nid_hfnmc)

         CALL histbeg_phy("histhfNMC",itau_phy, zjulian, dtime, 
     .                 nhori, nid_hfnmc)
c
         CALL histvert(nid_hfnmc, "plev", "pressure", "Pa",
     .                 nlevSTD3, rlevSTD3, nvert,"down")
c
cIM Astuce MAF: remplacer inst par ave pour les variables NMC pour avoir
cIM             le time_counter et les bounds
c
ccc
ccc Champs 3D interpolles sur des niveaux de pression du NMC
ccc
c
c ATTENTION : pour AMIP2 on interpole t,u,v,wphi,q,rh
c             sur les niveaux du NMC et on somme & moyenne
c             toutes les freq_moyNMC secondes par des routines undefSTD et
c             moy_undefSTD pour eliminer les valeurs "undef"
c             de la moyenne mensuelle
c ======> le "inst(X)" ci-dessous est par consequence factice !
c
          CALL histdef(nid_hfnmc, "tnondef",
     .                 "Valeurs non-definies","-",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "ta",
     .                 "Air temperature","K",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_hfnmc, "zg",
     .                "Geopotential height", "m",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "hus",
     .                 "Specific humidity","1",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_hfnmc, "hur",
     .                 "Relative humidity", "%",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "ua",
     .                 "Eastward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "va",
     .                 "Northward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "wap",
     .                 "Lagrangian tendency of air pressure","Pa s-1",
     .                iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "psbg",
     .         "Pressure sfce below ground","%",
     .         iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "uv",
     .         "uv ",
     .         "m2/s2",iim,jj_nb,nhori, nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "vq",
     .         "vq ",
     .         "m/s * (kg/kg)",iim,jj_nb,nhori, 
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "vT",
     .         "vT ", 
     .         "mK/s",iim,jj_nb,nhori, 
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "wq",
     .         "wq ", 
     .         "(Pa/s)*(kg/kg)",iim,jj_nb,nhori,
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "vphi",
     .         "vphi ", 
     .         "m2/s",iim,jj_nb,nhori, 
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "wT",
     .         "wT ", 
     .         "K*Pa/s",iim,jj_nb,nhori,
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "uxu",
     .         "u2 ", 
     .         "m2/s2",iim,jj_nb,nhori,
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "vxv",
     .         "v2 ", 
     .         "m2/s2",iim,jj_nb,nhori,
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_hfnmc, "TxT",
     .         "T2 ", 
     .         "K2",iim,jj_nb,nhori,
     .          nlevSTD3,1,nlevSTD3, nvert, 32,
     .         "ave(X)", zout,zout)
c
         CALL histend(nid_hfnmc)
c
      ENDIF !(ok_histNMC(2)) THEN
c
c$OMP END MASTER
