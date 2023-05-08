!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ini_histday.h,v 1.2 2005/01/27 10:06:12 fairhead Exp $
!
      IF (ok_journe) THEN

         zsto = dtime
         zout = dtime * REAL(ecrit_day)

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         call histbeg_phy("histday.nc",itau_phy,
     .                    zjulian,dtime,nhori,nid_day)

!$OMP MASTER
         CALL histvert(nid_day, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------
      IF(lev_histday.GE.1) THEN
c
ccccccccccccc 2D fields, basics
c
         CALL histdef(nid_day, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_day, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_day, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "ve", "Merid energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "cdragh", "Drag coef on T", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "cdragm", "Drag coef on U", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histday.GE.1
c
c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN
c
ccccccccccccc 3D fields, basics
c
         CALL histdef(nid_day, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "dudyn", "Dynamics dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "duvdf", "Boundary-layer dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "mang", "Angular momentum", "kg m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "Kz", "vertical diffusion coef", "m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
c plusieurs traceurs
          if (iflag_trac.eq.1) THEN
            DO iq=1,nqmax
             IF (iq.LE.99) THEN
          WRITE(str2,'(i2.2)') iq
          CALL histdef(nid_day, tname(iq), ttext(iq), "ppm",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
             ELSE
          PRINT*, "Trop de traceurs"
          CALL abort
             ENDIF 
            ENDDO
          endif
c
         CALL histdef(nid_day, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histday.GE.2
c
c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN
c
cccccccccccccccccc  Radiative transfer
c
c 2D
c
         CALL histdef(nid_day, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
c 3D
c
         CALL histdef(nid_day, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "fluxvdf", "PBL net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "fluxdyn", "Dyn. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "fluxajs", "Dry adj. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "fluxec", "Cin. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)
c
      ENDIF !lev_histday.GE.3
c
c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN
c
         CALL histdef(nid_day, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "dtphy", "Physics dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_day, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "dtec", "Cinetic dissip dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "duajs", "Dry convection dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dugwo", "GW oro dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dugwno", "GW non-oro dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "dvvdf", "Boundary-layer dV", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histday.GE.4
c
c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN
c
c        call histdef(nid_day, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)
c
c        call histdef(nid_day, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_day, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

         CALL histend(nid_day)
!$OMP END MASTER

      ENDIF ! fin de test sur ok_journe
