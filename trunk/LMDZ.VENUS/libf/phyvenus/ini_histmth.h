!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ini_histmth.h,v 1.4 2005/03/09 12:30:16 fairhead Exp $
!
      IF (ok_mensuel) THEN

         zsto = dtime
         zout = dtime * ecrit_mth

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         call histbeg_phy("histmth.nc",itau_phy,
     .                    zjulian,dtime,nhori,nid_mth)

!$OMP MASTER
         CALL histvert(nid_mth, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------
      IF(lev_histmth.GE.1) THEN
c
ccccccccccccc 2D fields, basics
c
         CALL histdef(nid_mth, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_mth, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_mth, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "ve", "Merid energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "cdragh", "Drag coef on T", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "cdragm", "Drag coef on U", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histmth.GE.1
c
c-------------------------------------------------------
      IF(lev_histmth.GE.2) THEN
c
ccccccccccccc 3D fields, basics
c
         CALL histdef(nid_mth, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dudyn", "Dynamics dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "duvdf", "Boundary-layer dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "mang", "Angular momentum", "kg m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "Kz", "vertical diffusion coef", "m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "mmean", "Mean molecular mass", "g/mol",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "rho", "Air density [mass/Vol]", "kg/m3",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)


c plusieurs traceurs  !!!outputs in [vmr]
       IF (iflag_trac.eq.1) THEN
            DO iq=1,nqmax
             IF (iq.LE.99) THEN
          WRITE(str2,'(i2.2)') iq
          CALL histdef(nid_mth, tname(iq), ttext(iq), "mol/mol",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
             ELSE
          PRINT*, "Trop de traceurs"
          CALL abort
             ENDIF 
            ENDDO
       ENDIF

       IF (callthermos .and. ok_chem) THEN
          CALL histdef(nid_mth, "d_qmoldifCO2", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
          CALL histdef(nid_mth, "d_qmoldifO3p", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
          CALL histdef(nid_mth, "d_qmoldifN2", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
       ENDIF
c
         CALL histdef(nid_mth, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histmth.GE.2
c
c-------------------------------------------------------
      IF(lev_histmth.GE.3) THEN
c
cccccccccccccccccc  Radiative transfer
c
c 2D
c
         CALL histdef(nid_mth, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)
c
c 3D
c
         CALL histdef(nid_mth, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "fluxvdf", "PBL net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "fluxdyn", "Dyn. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "fluxajs", "Dry adj. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "fluxec", "Cin. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)
c
      ENDIF !lev_histmth.GE.3
c
c-------------------------------------------------------
      IF(lev_histmth.GE.4) THEN
c
         CALL histdef(nid_mth, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "dtphy", "Physics dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth,"dtswrNLTE", "SWNLTE radiation dT", 
     .                "K/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)", zsto,zout)
c
         CALL histdef(nid_mth,"dtswrDCrisp","SWDCrisp radiation dT",
     .                "K/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)", zsto,zout)
c
	 CALL histdef(nid_mth, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth,"dtlwrNLTE", "LWNLTE radiation dT",
c    .                "K/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth,"dtlwrLTE","LW_LTE radiation dT",
c    .                "K/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "dteuv", "UV radiation dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c        CALL histdef(nid_mth, "dtcond", "Therm conduction", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c        CALL histdef(nid_mth, "dumolvis", "molec viscosity (u)"
c    .                ,"m/s2",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c    .                "ave(X)", zsto,zout)
c        CALL histdef(nid_mth, "dvmolvis", "molec viscosity (v)", 
c    .                "m/s2",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c    .                "ave(X)", zsto,zout)

c
c        CALL histdef(nid_mth, "dtec", "Cinetic dissip dT", "K/s",
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

c         CALL histdef(nid_mth, "co2_vmr", "volume mixture ratio co2", 
c     .		"cm-3",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c     .                "ave(X)", zsto,zout)

c         CALL histdef(nid_mth, "o_vmr", "volume mixture ratio o",
c     .                "cm-3",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c     .                "ave(X)", zsto,zout)
c         CALL histdef(nid_mth, "co_vmr", "volume mixture ratio co", 
c     .		"cm-3",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c     .                "ave(X)", zsto,zout)

c         CALL histdef(nid_mth, "n2_vmr", "volume mixture ratio n2",
c     .                "cm-3",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
c     .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "dvvdf", "Boundary-layer dV", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histmth.GE.4
c
c-------------------------------------------------------
      IF(lev_histmth.GE.5) THEN
c
c        call histdef(nid_mth, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)
c
c        call histdef(nid_mth, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
c        CALL histdef(nid_mth, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)
c
      ENDIF !lev_histmth.GE.5
c-------------------------------------------------------

         CALL histend(nid_mth)
!$OMP END MASTER

      ENDIF ! fin de test sur ok_mensuel

