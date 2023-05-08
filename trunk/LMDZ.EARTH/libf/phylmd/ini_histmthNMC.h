!
! $Id: ini_histmthNMC.h 1403 2010-07-01 09:02:53Z fairhead $
!
c$OMP MASTER
c
      IF (ok_histNMC(1)) THEN
c
       zout = freq_outNMC(1)
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
            znivsig(ll)=REAL(ll)
         ENDDO
cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)
cym         CALL histbeg("histNMC.nc", iim,zx_lon(:,1), jjmp1,zx_lat(1,:),
cym     .                 1,iim,1,jjmp1, itau_phy, zjulian, dtime, 
cym     .                 nhori, nid_mthnmc)

         CALL histbeg_phy("histmthNMC",itau_phy, zjulian, dtime, 
     .                 nhori, nid_mthnmc)
c
         CALL histvert(nid_mthnmc, "plev", "pressure", "Pa",
     .                 nlevSTD, rlevSTD, nvert,"down")
c
cIM Astuce MAF: remplacer inst par ave pour les variables NMC pour avoir
cIM             le time_counter et les bounds
cIM 
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
c
          CALL histdef(nid_mthnmc, "tnondef",
     .                 "Valeurs non-definies","-",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "ta",
     .                 "Air temperature","K",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_mthnmc, "zg",
     .                "Geopotential height", "m",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "hus",
     .                 "Specific humidity","1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
         CALL histdef(nid_mthnmc, "hur",
     .                 "Relative humidity", "%",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "ua",
     .                 "Eastward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "va",
     .                 "Northward wind","m s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "wap",
     .                 "Lagrangian tendency of air pressure","Pa s-1",
     .                iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .                "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "psbg",
     .         "Pressure sfce below ground","%",
     .         iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "uv",
     .         "uv ",
     .         "m2/s2",iim,jj_nb,nhori, nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "vq",
     .         "vq ",
     .         "m/s * (kg/kg)",iim,jj_nb,nhori, 
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "vT",
     .         "vT ", 
     .         "mK/s",iim,jj_nb,nhori, 
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "wq",
     .         "wq ", 
     .         "(Pa/s)*(kg/kg)",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "vphi",
     .         "vphi ", 
     .         "m2/s",iim,jj_nb,nhori, 
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "wT",
     .         "wT ", 
     .         "K*Pa/s",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "uxu",
     .         "u2 ", 
     .         "m2/s2",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "vxv",
     .         "v2 ", 
     .         "m2/s2",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "TxT",
     .         "T2 ", 
     .         "K2",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          CALL histdef(nid_mthnmc, "tro3",
     .         "Ozone mole fraction",
     .         "1e-9",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
c
          if (read_climoz == 2) THEN
           CALL histdef(nid_mthnmc, "tro3_daylight",
     .         "Daylight ozone mole fraction",
     .         "1e-9",iim,jj_nb,nhori,
     .          nlevSTD,1,nlevSTD, nvert, 32,
     .         "ave(X)", zout,zout)
          endif
c
         CALL histend(nid_mthnmc)
c
      ENDIF !(ok_histNMC(1)) THEN
c
c$OMP END MASTER
