!
! $Id $
!
  IF (ecrit_tra>0. .AND. config_inca == 'none') THEN
!$OMP MASTER 
     CALL ymds2ju(annee_ref, 1, day_ref, 0.0, zjulian)
     CALL histbeg_phy("histrac", itau_phy, zjulian, pdtphys,nhori, nid_tra)
     CALL histvert(nid_tra, "presnivs", "Vertical levels", "mb",klev, presnivs, nvert)

     zsto = pdtphys
     zout = ecrit_tra
     CALL histdef(nid_tra, "phis", "Surface geop. height", "-",   &
          iim,jj_nb,nhori, 1,1,1, -99, 32,"once",  zsto,zout)
     CALL histdef(nid_tra, "aire", "Grid area", "-",              &
          iim,jj_nb,nhori, 1,1,1, -99, 32,"once",  zsto,zout)
     CALL histdef(nid_tra, "zmasse", "column density of air in cell", &
          "kg m-2", iim, jj_nb, nhori, klev, 1, klev, nvert, 32, "ave(X)", &
          zsto,zout)

!TRACEURS
!----------------
     DO it = 1,nbtr
        iiq = niadv(it+2)

! CONCENTRATIONS
        CALL histdef(nid_tra, tname(iiq), ttext(iiq), "U/kga",    &
             iim,jj_nb,nhori, klev,1,klev,nvert, 32,"ave(X)", zsto,zout)

! TD LESSIVAGE
        IF (lessivage .AND. aerosol(it)) THEN
           CALL histdef(nid_tra, "fl"//tname(iiq),"Flux "//ttext(iiq), &
                "U/m2/s",iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
        END IF

! TD THERMIQUES
        IF (iflag_thermals.gt.0) THEN
           CALL histdef(nid_tra, "d_tr_th_"//tname(iiq),      &
                "tendance thermique"// ttext(iiq), "?",       &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
                "ave(X)", zsto,zout)
        ENDIF

! TD CONVECTION
        IF (iflag_con.GE.2) THEN
           CALL histdef(nid_tra, "d_tr_cv_"//tname(iiq),   &
                "tendance convection"// ttext(iiq), "?",   &
                iim,jj_nb,nhori, klev,1,klev,nvert, 32,    &
                "ave(X)", zsto,zout)
        ENDIF

! TD COUCHE-LIMITE
        CALL histdef(nid_tra, "d_tr_cl_"//tname(iiq),      &
             "tendance couche limite"// ttext(iiq), "?",   &
             iim,jj_nb,nhori, klev,1,klev,nvert, 32,       &
             "ave(X)", zsto,zout)
     ENDDO
!---------------   
!
! VENT (niveau 1)
     CALL histdef(nid_tra, "pyu1", "Vent niv 1", "-",      &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)     
     CALL histdef(nid_tra, "pyv1", "Vent niv 1", "-",      &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)

! TEMPERATURE DU SOL
     CALL histdef(nid_tra, "ftsol1", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "ftsol2", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "ftsol3", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst",  zout,zout)
     CALL histdef(nid_tra, "ftsol4", "temper sol", "-",    &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)

! NATURE DU SOL
     CALL histdef(nid_tra, "psrf1", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf2", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf3", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 &
          "inst(X)",  zout,zout)
     CALL histdef(nid_tra, "psrf4", "nature sol", "-",     &
          iim,jj_nb,nhori, 1,1,1, -99, 32,                 & 
          "inst(X)",  zout,zout)
! DIVERS
     CALL histdef(nid_tra, "pplay", "pressure","-",        &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "T", "temperature","K",         &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "inst(X)", zout,zout)
     CALL histdef(nid_tra, "mfu", "flux u mont","-",       &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "mfd", "flux u decen","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "en_u", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "en_d", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "de_d", "flux u mont","-",      &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "de_u", "flux u decen","-",     &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)
     CALL histdef(nid_tra, "coefh", "turbulent coef","-",  &
          iim,jj_nb,nhori, klev,1,klev,nvert, 32,          &
          "ave(X)", zsto,zout)   
     
     CALL histend(nid_tra)
!$OMP END MASTER
  END IF ! ecrit_tra>0. .AND. config_inca == 'none'
  
