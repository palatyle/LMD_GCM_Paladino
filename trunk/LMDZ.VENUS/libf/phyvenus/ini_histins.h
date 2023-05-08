!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ini_histins.h,v 1.1.1.1 2004/05/19 12:53:08 lmdzadmin Exp $
!
      IF (ok_instan) THEN

         zsto = dtime * ecrit_ins
         zout = dtime * ecrit_ins

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         call histbeg_phy("histins.nc",itau_phy,
     .                    zjulian,dtime,nhori,nid_ins)

!$OMP MASTER
         CALL histvert(nid_ins, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------
      IF(lev_histins.GE.1) THEN
c
ccccccccccccc 2D fields, basics
c
         CALL histdef(nid_ins, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_ins, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)
c
         CALL histdef(nid_ins, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "ve", "Merid energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "cdragh", "Drag coef on T", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "cdragm", "Drag coef on U", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
      ENDIF !lev_histins.GE.1
c
c-------------------------------------------------------
      IF(lev_histins.GE.2) THEN
c
ccccccccccccc 3D fields, basics
c
         CALL histdef(nid_ins, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dudyn", "Dynamics dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "duvdf", "Boundary-layer dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "mang", "Angular momentum", "kg m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "Kz", "vertical diffusion coef", "m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "mmean", "Mean molecular mass", "g/mol",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)

         CALL histdef(nid_ins, "rho", "Air density [mass/Vol]", "kg/m3",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)


c plusieurs traceurs  !!!outputs in [vmr]
       IF (iflag_trac.eq.1) THEN
            DO iq=1,nqmax
             IF (iq.LE.99) THEN
          WRITE(str2,'(i2.2)') iq
          CALL histdef(nid_ins, tname(iq), ttext(iq), "mol/mol",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
             ELSE
          PRINT*, "Trop de traceurs"
          CALL abort
             ENDIF 
            ENDDO
       ENDIF

       IF (callthermos .and. ok_chem) THEN
          CALL histdef(nid_ins, "d_qmoldif CO2", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
          CALL histdef(nid_ins, "d_qmoldif O3p", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
          CALL histdef(nid_ins, "d_qmoldif N2", "Dif molec" , "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
       ENDIF
c
         CALL histdef(nid_ins, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
          if (ok_cloud.and.(cl_scheme.eq.1)) THEN

          if (nb_mode.GE.1) THEN 
           
c
         CALL histdef(nid_ins, "NBRTOTm1", "Nbr total droplet",
     .                "#/cm3",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "R_MEDIANm1", "Median radius
c     .    for log normal distribution" ,
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "STDDEVm1", "Std Deviation
c     .    for log normal distribution",
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

          if (nb_mode.GE.2) THEN

c
         CALL histdef(nid_ins, "NBRTOTm2", "Nbr total droplet",
     .                "#/cm3",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "R_MEDIANm2", "Median radius
c     .    for log normal distribution" ,
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "STDDEVm2", "Std Deviation
c     .    for log normal distribution",
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

          if (nb_mode.GE.3) THEN
          
c
         CALL histdef(nid_ins, "NBRTOTm3", "Nbr total droplet", "#/cm3",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "R_MEDIANm3", "Median radius
c     .   for log normal distribution" ,
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

c
c         CALL histdef(nid_ins, "STDDEVm3", "Std Deviation
c     .    for log normal distribution",
c     .                "fraction",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
c     .                "ins(X)", zsto,zout)
c

         ENDIF
         ENDIF
         ENDIF
         
c
         CALL histdef(nid_ins, "WH2SO4", "Weight fraction H2SO4",
     .                "fraction",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c

c
         CALL histdef(nid_ins, "rho_droplet", "density cloud droplet",
     .                "kg.m-3",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c

		ENDIF
		
          if (ok_sedim.and.(cl_scheme.eq.1)) THEN
c
         CALL histdef(nid_ins, "d_tr_sed_H2SO4", "var mmr from sedim",
     .                "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c

c
         CALL histdef(nid_ins, "d_tr_sed_H2O", "var mmr from sedim",
     .                "kg/kg",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c

c
         CALL histdef(nid_ins, "F_sedim", "tendency from sedim",
     .                "kg.m-2.s-1",
     .                nbp_lon,jj_nb,nhori, klev,1,klev, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
		ENDIF

      ENDIF !lev_histins.GE.2
c
c-------------------------------------------------------
      IF(lev_histins.GE.3) THEN
c
cccccccccccccccccc  Radiative transfer
c
c 2D
c
         CALL histdef(nid_ins, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ins(X)", zsto,zout)
c
c 3D
c
         CALL histdef(nid_ins, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "fluxvdf", "PBL net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "fluxdyn", "Dyn. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "fluxajs", "Dry adj. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "fluxec", "Cin. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ins(X)", zsto,zout)
c
      ENDIF !lev_histins.GE.3
c
c-------------------------------------------------------
      IF(lev_histins.GE.4) THEN
c
         CALL histdef(nid_ins, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "dtphy", "Physics dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "dtec", "Cinetic dissip dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "duajs", "Dry convection dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dugwo", "GW oro dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
         CALL histdef(nid_ins, "dugwno", "GW non-oro dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "dvvdf", "Boundary-layer dV", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ins(X)", zsto,zout)
c
      ENDIF !lev_histins.GE.4
c
c-------------------------------------------------------
      IF(lev_histins.GE.5) THEN
c
c        call histdef(nid_ins, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ins(X)", zsto,zout)
c
c        call histdef(nid_ins, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
c        CALL histdef(nid_ins, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ins(X)", zsto,zout)
c
      ENDIF !lev_histins.GE.5
c-------------------------------------------------------
c
         CALL histend(nid_ins)
!$OMP END MASTER

      ENDIF
