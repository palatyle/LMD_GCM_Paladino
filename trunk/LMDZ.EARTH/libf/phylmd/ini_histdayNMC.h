!
! $Header$
!
c$OMP MASTER
      IF (ok_histNMC(2)) THEN
c
       zstophy = dtime
       zstohf = ecrit_hf
       zstomth = ecrit_mth
       zout = freq_outNMC(2)
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
cym     .                 nhori, nid_daynmc)

         CALL histbeg_phy("histdayNMC",itau_phy, zjulian, dtime, 
     .                 nhori, nid_daynmc)
c
        IF(lev_histdayNMC.EQ.nlevSTD) THEN 
         CALL histvert(nid_daynmc, "plev", "pressure", "Pa",
     .                 nlevSTD, rlevSTD, nvert,"down")
        ELSE IF(lev_histdayNMC.EQ.nlevSTD8) THEN 
         CALL histvert(nid_daynmc, "plev", "pressure", "Pa",
     .                 nlevSTD8, rlevSTD8, nvert,"down")
        ENDIF
c
cIM Astuce MAF: remplacer inst par ave pour les variables NMC pour avoir
cIM             le time_counter et les bounds
c
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
        IF(lev_histdayNMC.EQ.nlevSTD) THEN 
          CALL histdef(nid_daynmc, "tnondef",
     .                 "Valeurs non-definies","-",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "ta",
     .                 "Air temperature","K",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_daynmc, "zg",
     .                "Geopotential height", "m",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "hus",
     .                 "Specific humidity","1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_daynmc, "hur",
     .                 "Relative humidity", "%",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "ua",
     .                 "Eastward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "va",
     .                 "Northward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "wap",
     .                 "Lagrangian tendency of air pressure","Pa s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
        ELSE IF(lev_histdayNMC.EQ.nlevSTD8) THEN 
c
          CALL histdef(nid_daynmc, "tnondef",
     .                 "Valeurs non-definies","-",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "ta",
     .                 "Air temperature","K",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_daynmc, "zg",
     .                "Geopotential height", "m",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "hus",
     .                 "Specific humidity","1",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_daynmc, "hur",
     .                 "Relative humidity", "%",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "ua",
     .                 "Eastward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "va",
     .                 "Northward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_daynmc, "wap",
     .                 "Lagrangian tendency of air pressure","Pa s-1",
     .                iim,jj_nb,nhori, nlevSTD8,1,nlevSTD8, nvert, 32,
     .                "ave(X)", zout,zout)
c
        ENDIF
c
         CALL histend(nid_daynmc)
c
      ENDIF !(ok_histNMC(2)) THEN
c$OMP END MASTER
