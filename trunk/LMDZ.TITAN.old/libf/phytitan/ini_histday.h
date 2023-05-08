      IF (ok_journe) THEN

         zsto = dtime
         zout = dtime * REAL(ecrit_day)
c zsto1: pour des flux radiatifs calcules tous les radpas appels physiq
         zsto1= dtime * REAL(radpas)

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         CALL histbeg_phy("histday.nc", itau_phy, zjulian, dtime,
     .                 nhori, nid_day)

!$OMP MASTER
         CALL histvert(nid_day, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------
      IF(lev_histday.GE.1) THEN

ccccccccccccc 2D fields, invariables

         CALL histdef(nid_day, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

         CALL histdef(nid_day, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

ccccccc axe Ls
         CALL histdef(nid_day, "ls", "Solar longitude", "degrees",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

ccccccccccccc 2D fields, variables

         CALL histdef(nid_day, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_day, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

c        CALL histdef(nid_day, "ve", "Merid energy transport", "-",
c     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c     .                "ave(X)", zsto,zout)

      ENDIF !lev_histday.GE.1

c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN

ccccccccccccc 3D fields, basics

         CALL histdef(nid_day, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_day, "duvdf", "Boundary-layer dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dudyn", "Dynamics dU", "m/s2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

cccccccccccccccccc  Tracers

         if (iflag_trac.eq.1) THEN
          if (microfi.ge.1) then
c           DO iq=1,nmicro
c             CALL histdef(nid_day, tname(iq), ttext(iq), "n/m2",
c     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c     .                "ave(X)", zsto,zout)
c           ENDDO
             CALL histdef(nid_day, "qaer","nb tot aer" , "n/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)

            if (clouds.eq.1) then
             CALL histdef(nid_day, "qnoy","nb tot noy" , "n/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_day, "qgl1","V tot gl1" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_day, "qgl2","V tot gl2" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
             CALL histdef(nid_day, "qgl3","V tot gl3" , "m3/m2",
     .                    nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                    "ave(X)", zsto,zout)
c--------------
c ----- SATURATION ESP NUAGES
               CALL histdef(nid_day,"ch4sat", "saturation CH4", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"c2h6sat", "saturation C2H6", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"c2h2sat", "saturation C2H2", "--",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- RESERVOIR DE SURFACE
               CALL histdef(nid_day, "reserv", "Reservoir surface","m",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- ECHANGE GAZ SURF/ATM (evaporation)
               CALL histdef(nid_day, "evapch4", "Evaporation CH4","m",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- PRECIPITATIONS (precipitations cumulatives)
               CALL histdef(nid_day,"prech4","Precip CH4","m",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"prec2h6","Precip C2H6",
     .                "m",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"prec2h2","Precip C2H2",
     .                "m",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"prenoy","Precip NOY",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
               CALL histdef(nid_day,"preaer","Precip AER",
     .                "um/s",nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
     .                "ave(X)", zsto,zout)
c --------------
c ----- FLUX GLACE
               CALL histdef(nid_day,"flxgl1", "flux gl CH4",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_day,"flxgl2", "flux gl C2H6",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
               CALL histdef(nid_day,"flxgl3", "flux gl C2H2",
     .              "kg/m2/s",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .              "ave(X)", zsto,zout)
c --------------
c ----- RAYON DES GOUTTES
               CALL histdef(nid_day,"rcldbar", "rayon moyen goutte",
     .                "m",nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
            endif
	  endif
c --------------
c ----- TRACEURS CHIMIQUES
	  if (nmicro.lt.nqmax) then
           DO iq=nmicro+1,nqmax
         CALL histdef(nid_day, tname(iq), ttext(iq), "ppm",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)
           ENDDO
	  endif
         endif

      ENDIF !lev_histday.GE.2

c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

         CALL histdef(nid_day, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_day, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

         CALL histdef(nid_day, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "ave(X)", zsto1,zout)

c 3D

         CALL histdef(nid_day, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto1,zout)

         CALL histdef(nid_day, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "ave(X)", zsto1,zout)

c --------------
c ----- OPACITE BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"thv"//str2,"Haze Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"thi"//str2,"Haze Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- EXTINCTION BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"khv"//str2,"Haze ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"khi"//str2,"Haze ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- OPACITE GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"tgv"//str2,"Gas Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"tgi"//str2,"Gas Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- EXTINCTION GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"kgv"//str2,"Gas ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_day,"kgi"//str2,"Gas ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto1,zout)
         ENDDO

c --------------
c ----- OPACITE NUAGES
         if (clouds.eq.1) then
           CALL histdef(nid_day,"tcld","Cld Opa proxy",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto,zout)

c --------------
c ----- EXTINCTION NUAGES
           CALL histdef(nid_day,"kcld","Cld Ext proxy",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ave(X)",zsto,zout)
         endif

      ENDIF !lev_histday.GE.3

c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN

         CALL histdef(nid_day, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dtphy", "Physics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

         CALL histdef(nid_day, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "ave(X)", zsto,zout)

c        CALL histdef(nid_day, "dtec", "Cinetic dissip dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "ave(X)", zsto,zout)

      ENDIF !lev_histday.GE.4

c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN


c        call histdef(nid_day, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)

c        call histdef(nid_day, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "ave(X)", zsto,zout)

c        CALL histdef(nid_day, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

c        CALL histdef(nid_day, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "ave(X)", zsto,zout)

      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

         CALL histend(nid_day)

      ENDIF ! fin de test sur ok_journe
