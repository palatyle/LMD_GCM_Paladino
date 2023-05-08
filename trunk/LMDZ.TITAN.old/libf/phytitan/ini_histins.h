      IF (ok_instan) THEN

          zsto = dtime * REAL(ecrit_ins)
          zout = dtime * REAL(ecrit_ins)

         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, zero, zjulian)

         CALL histbeg_phy("histins.nc", itau_phy, zjulian, dtime,
     .                 nhori, nid_ins)

!$OMP MASTER
         CALL histvert(nid_ins, "presnivs", "Vertical levels", "Pa",
     .                 klev, presnivs, nvert)

c-------------------------------------------------------

      IF(lev_histday.GE.1) THEN

ccccccccccccc 2D fields, invariables

         CALL histdef(nid_ins, "phis", "Surface geop. height", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

         CALL histdef(nid_ins, "aire", "Grid area", "-",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "once",  zsto,zout)

ccccccc axe Ls
         CALL histdef(nid_ins, "ls", "Solar longitude", "degrees",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

ccccccccccccc 2D fields, variables

         CALL histdef(nid_ins, "tsol", "Surface Temperature", "K",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "psol", "Surface Pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "ue", "Zonal energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "ve", "Merid energy transport", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "inst(X)", zsto,zout)

      ENDIF !lev_histday.GE.1

c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN

ccccccccccccc 3D fields, basics

         CALL histdef(nid_ins, "temp", "Air temperature", "K",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "pres", "Air pressure", "Pa",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "geop", "Geopotential height", "m",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "vitu", "Zonal wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "vitv", "Meridional wind", "m/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "vitw", "Vertical wind", "Pa/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "tops", "Solar rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "duvdf", "Boundary-layer dU", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "dudyn", "Dynamics dU", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "inst(X)", zsto,zout)

      ENDIF !lev_histday.GE.2

c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN

cccccccccccccccccc  Tracers

         if (iflag_trac.eq.1) THEN
          if (microfi.ge.1) then
           DO iq=1,nmicro
         CALL histdef(nid_ins, tname(iq), ttext(iq), "n/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)
           ENDDO
	  endif
	  if (nmicro.lt.nqmax) then
           DO iq=nmicro+1,nqmax
         CALL histdef(nid_ins, tname(iq), ttext(iq), "ppm",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)
           ENDDO
	  endif
         endif

cccccccccccccccccc  Radiative transfer

c 2D

         CALL histdef(nid_ins, "topl", "IR rad. at TOA", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "sols", "Solar rad. at surf.", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "soll", "IR rad. at surface", "W/m2",
     .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
     .                "inst(X)", zsto,zout)

c 3D

         CALL histdef(nid_ins, "SWnet", "Net SW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "LWnet", "Net LW flux","W/m2",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert,
     .                32, "inst(X)", zsto,zout)

c --------------
c ----- OPACITE BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"thv"//str2,"Haze Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"thi"//str2,"Haze Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

c --------------
c ----- EXTINCTION BRUME
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"khv"//str2,"Haze ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"khi"//str2,"Haze ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

c --------------
c ----- OPACITE GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"tgv"//str2,"Haze Opa Vis",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"tgi"//str2,"Haze Opa IR",
     .                "--",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

c --------------
c ----- EXTINCTION GAZ
         DO k=7,NSPECV,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"kgv"//str2,"Haze ext Vis ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

         DO k=8,NSPECI,10
           write(str2,'(i2.2)') k
         CALL histdef(nid_ins,"kgi"//str2,"Haze ext IR ",
     .                "m-1",nbp_lon,jj_nb,nhori,klev,1,klev,nvert,32,
     .                "ins(X)",zsto,zout)
         ENDDO

      ENDIF !lev_histday.GE.3

c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN

         CALL histdef(nid_ins, "dtdyn", "Dynamics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "dtphy", "Physics dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "dtvdf", "Boundary-layer dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "dtajs", "Dry adjust. dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "dtswr", "SW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

         CALL histdef(nid_ins, "dtlwr", "LW radiation dT", "K/s",
     .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
     .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "dtec", "Cinetic dissip dT", "K/s",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "dvvdf", "Boundary-layer dV", "m/s2",
c    .                nbp_lon,jj_nb,nhori, klev,1,klev,nvert, 32,
c    .                "inst(X)", zsto,zout)

      ENDIF !lev_histday.GE.4

c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN


c        call histdef(nid_ins, "taux", 
c    $         "Zonal wind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "inst(X)", zsto,zout)

c        call histdef(nid_ins, "tauy", 
c    $         "Meridional xind stress", "Pa",  
c    $         nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32,
c    $         "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "cdrm", "Momentum drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "inst(X)", zsto,zout)

c        CALL histdef(nid_ins, "cdrh", "Heat drag coef.", "-",
c    .                nbp_lon,jj_nb,nhori, 1,1,1, nvert, 32, 
c    .                "inst(X)", zsto,zout)

      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

         CALL histend(nid_ins)

      ENDIF
