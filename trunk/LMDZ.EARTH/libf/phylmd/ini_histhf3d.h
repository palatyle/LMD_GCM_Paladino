c $Header$
c
c sorties hf 3d
c
        zstohf = ecrit_hf
        zout = ecrit_hf
c
c       PRINT*, 'La frequence de sortie hf3d est de ', ecrit_hf
c
        idayref = day_ref
        CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c

cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlon,zx_lon)
cym         DO i = 1, iim
cym            zx_lon(i,1) = rlon(i+1)
cym            zx_lon(i,jjmp1) = rlon(i+1)
cym         ENDDO

cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)

cccIM      CALL histbeg("histhf", iim,zx_lon, jjmp1,zx_lat,
cym         CALL histbeg("histhf3d", iim,zx_lon(:,1), jjmp1,zx_lat(1,:),
cym     .                 1,iim,1,jjmp1, itau_phy, zjulian, dtime, 
cym     .                 nhori, nid_hf3d)
         CALL histbeg_phy("histhf3d", itau_phy, zjulian, dtime, 
     .                 nhori, nid_hf3d)

        CALL histvert(nid_hf3d, "presnivs", "Vertical levels", "mb",
     .                 klev, presnivs/100., nvert)
c
c Champs 3D:
c
        CALL histdef(nid_hf3d, "temp", "Air temperature", "K",
     .                iim,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zstohf,zout)
c
        CALL histdef(nid_hf3d, "ovap", "Specific humidity", "kg/kg",
     .                iim,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zstohf,zout)
c
        CALL histdef(nid_hf3d, "vitu", "Zonal wind", "m/s",
     .                iim,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zstohf,zout)
c
        CALL histdef(nid_hf3d, "vitv", "Meridional wind", "m/s",
     .                iim,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zstohf,zout)
c
        CALL histend(nid_hf3d)
c
