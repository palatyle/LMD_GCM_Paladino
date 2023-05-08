      IF (ok_mensuel) THEN

         zsto = dtime
         zout = dtime * REAL(ecrit_mth)
c zsto1: pour des flux radiatifs calcules tous les radpas appels physiq
         zsto1= dtime * REAL(radpas)

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         CALL histbeg_phy("histmth.nc", itau_phy, zjulian, dtime,
     .                 nhori, nid_mth)

!$OMP MASTER
         CALL histvert(nid_mth, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------
      IF(lev_histmth.GE.1) THEN

ccccccccccccc 2D fields, invariables

         CALL histdef(nid_mth, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

         CALL histdef(nid_mth, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

ccccccc axe Ls
         CALL histdef(nid_mth, "ls", "Solar longitude", "degrees",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

ccccccccccccc 2D fields, variables

         CALL histdef(nid_mth, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "ve", "Merid energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

      ENDIF !lev_histmth.GE.1

c-------------------------------------------------------
      IF(lev_histmth.GE.2) THEN

ccccccccccccc 3D fields, basics

         CALL histdef(nid_mth, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "Kz", "vertical diffusion coef", "m2/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_mth, "duvdf", "Boundary-layer dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dudyn", "Dynamics dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

cccccccccccccccccc  Tracers

         if (iflag_trac.eq.1) THEN
          if (microfi.ge.1) then
c           DO iq=1,nmicro
c             CALL histdef(nid_mth, tname(iq), ttext(iq), "n/m2",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c     .                "ave(X)", zsto,zout)
c           ENDDO
             CALL histdef(nid_mth, "qaer","nb tot aer" , "n/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)

            if (clouds.eq.1) then
             CALL histdef(nid_mth, "qnoy","nb tot noy" , "n/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_mth, "qgl1","V tot gl1" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_mth, "qgl2","V tot gl2" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_mth, "qgl3","V tot gl3" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
c--------------
c ----- SATURATION ESP NUAGES
               CALL histdef(nid_mth,"ch4sat", "saturation CH4", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"c2h6sat", "saturation C2H6", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"c2h2sat", "saturation C2H2", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- RESERVOIR DE SURFACE
               CALL histdef(nid_mth, "reserv", "Reservoir surface","m",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- ECHANGE GAZ SURF/ATM (evaporation)
               CALL histdef(nid_mth, "evapch4", "Evaporation CH4","m",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- PRECIPITATIONS (precipitations moyennes)
               CALL histdef(nid_mth,"prech4","Precip CH4","um/s",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"prec2h6","Precip C2H6",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"prec2h2","Precip C2H2",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"prenoy","Precip NOY",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"preaer","Precip AER",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- FLUX GLACE
               CALL histdef(nid_mth,"flxgl1", "flux gl CH4",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"flxgl2", "flux gl C2H6",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"flxgl3", "flux gl C2H2",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
c --------------
c ----- Source/puits GLACE
               CALL histdef(nid_mth,"solch4", "dQ gl CH4",
     .              "m3/m3",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"solc2h6", "dQ gl C2H6",
     .              "m3/m3",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_mth,"solc2h2", "dQ gl C2H2",
     .              "m3/m3",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .              "ave(X)", zsto,zout)
c --------------
c ----- RAYON DES GOUTTES
               CALL histdef(nid_mth,"rcldbar", "rayon moyen goutte",
     .                "m",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
            endif
	  endif
c --------------
c ----- TRACEURS CHIMIQUES
	  if (nmicro.lt.nqmax) then
           DO iq=nmicro+1,nqmax
         CALL histdef(nid_mth, tname(iq), ttext(iq), "ppm",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
           ENDDO
c Condensation:
c          DO iq=nmicro+1,nqmax
c        CALL histdef(nid_mth, "c_"//tname(iq), "c_"//ttext(iq),
c    .        "ppm/s",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)
c          ENDDO
	  endif
         endif

      ENDIF !lev_histmth.GE.2

c-------------------------------------------------------
      IF(lev_histmth.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

         CALL histdef(nid_mth, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_mth, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_mth, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

c 3D

         CALL histdef(nid_mth, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto1,zout)

c        CALL histdef(nid_mth, "SWup", "upward SW flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto1,zout)

c        CALL histdef(nid_mth, "SWdn", "downward SW flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto1,zout)

         CALL histdef(nid_mth, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto1,zout)

c        CALL histdef(nid_mth, "LWup", "upward LW flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto1,zout)

c        CALL histdef(nid_mth, "LWdn", "downward LW flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto1,zout)

         CALL histdef(nid_mth, "fluxvdf", "PBL net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "fluxdyn", "Dyn. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "fluxajs", "Dry adj. net flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "fluxec", "Cin. net flux","W/m2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
c    .                32, "ave(X)", zsto,zout)

c --------------
c ----- OPACITE BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"thv"//str2,"Haze Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"thi"//str2,"Haze Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- EXTINCTION BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"khv"//str2,"Haze ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"khi"//str2,"Haze ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- OPACITE GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"tgv"//str2,"Gas Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"tgi"//str2,"Haze Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- EXTINCTION GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"kgv"//str2,"Gas ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_mth,"kgi"//str2,"Gas ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- OPACITE NUAGES
         if (clouds.eq.1) then
           CALL histdef(nid_mth,"tcld","Cld Opa proxy",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto,zout)

c --------------
c ----- EXTINCTION NUAGES
           CALL histdef(nid_mth,"kcld","Cld Ext proxy",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto,zout)
         endif

c --------------
c ----- OCCURENCE NUAGES
           do k=1,12
             write(str2,'(i2.2)') k
             CALL histdef(nid_mth,"occcld"//str2,"occ cld",
     .       "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .       "ave(X)",zsto,zout)
           enddo

      ENDIF !lev_histmth.GE.3

c-------------------------------------------------------
      IF(lev_histmth.GE.4) THEN

         CALL histdef(nid_mth, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dtphy", "Physics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_mth, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "dtec", "Cinetic dissip dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "dvvdf", "Boundary-layer dV", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)

      ENDIF !lev_histmth.GE.4

c-------------------------------------------------------
      IF(lev_histmth.GE.5) THEN


c        call histdef(nid_mth, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)

c        call histdef(nid_mth, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

c        CALL histdef(nid_mth, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

      ENDIF !lev_histmth.GE.5
c-------------------------------------------------------

         CALL histend(nid_mth)

      ENDIF ! fin de test sur ok_journe
