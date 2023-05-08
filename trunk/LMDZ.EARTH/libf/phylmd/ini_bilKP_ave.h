c
c $Id: ini_bilKP_ave.h 1403 2010-07-01 09:02:53Z fairhead $
c
      IF (ok_journe) THEN
c
         zsto = dtime
         zout = ecrit_day
         typeval=tave
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
cym         write(*,*)'zx_lon = ',zx_lon(:,1)
cym         write(*,*)'zx_lat = ',zx_lat(1,:)
cym         CALL histbeg("histbilKP_ave", iim,zx_lon(:,1), jjmp1,
cym     .                zx_lat(1,:),
cym     .                1,iim,1,jjmp1, itau_phy, zjulian, dtime,
cym     .                nhori, nid_bilKPave)
         CALL histbeg_phy("histbilKP_ave", itau_phy, zjulian, dtime,
     .                nhori, nid_bilKPave)

         write(*,*)'Journee ', itau_phy, zjulian
         CALL histvert(nid_bilKPave, "presnivs",
     .                "Vertical levels","mb",
     .                 klev, presnivs/100., nvert)
c
c
c Champs 3D:
c
         CALL histdef(nid_bilKPave,"ue",
     .   "Zonal energy transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"ve",
     .   "Merid energy transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"uq",
     .   "Zonal humidity transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"vq",
     .   "Merid humidity transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
c Champs 3D:
c
         CALL histdef(nid_bilKPave,"temp",
     .   "Air temperature","K",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"ovap",
     .   "Specific humidity","Kg/Kg",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"geop",
     .   "Geopotential height","m",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"vitu",
     .   "Zonal wind","m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"vitv",
     .   "Meridional wind","m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "vitw", 
     .   "Vertical wind", "m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "pres", 
     .   "Inter-Layer Air pressure",
     .                "Pa",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "play", 
     .   "Mean-Layer Air pressure",
     .                "Pa",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "oliq", 
     .   "Liquid water content", 
     .                "kg/kg",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtdyn", 
     .   "Dynamics dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dqdyn", 
     .   "Dynamics dQ", "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtcon", 
     .   "Convection dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "ducon", 
     .   "Convection du", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dvcon", 
     .   "Convection dv", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dqcon",
     .   "Convection dQ","Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtlsc", 
     .   "Condensation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dqlsc",
     .   "Condensation dQ","Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dtvdf",
     .   "Boundary-layer dT","K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dqvdf", 
     .   "Boundary-layer dQ", 
     .               "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dtajs",
     .   "Ajustement sec dT","K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dqajs",
     .   "Ajustement sec dQ", 
     .               "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dteva", 
     .   "Reevaporation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dqeva",
     .   "Reevaporation dQ",
     .                "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)

c
         CALL histdef(nid_bilKPave, "dtswr", 
     .   "SW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtsw0", 
     .   "SW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtlwr", 
     .   "LW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dtlw0", 
     .   "LW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"duvdf",
     .   "Boundary-layer dU","m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave,"dvvdf",
     .   "Boundary-layer dV","m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         IF (ok_orodr) THEN
         IF (ok_orolf) THEN
         CALL histdef(nid_bilKPave, "duoli",
     .   "Orography dU", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPave, "dvoli", 
     .   "Orography dV", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         ENDIF
         ENDIF
C
         CALL histdef(nid_bilKPave, "duphy",
     .   "Physiq dU","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPave, "dvphy",
     .   "Physiq dV","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPave, "dtphy",
     .   "Physiq dT","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPave, "dqphy",
     .   "Physiq dQ","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPave, "dqlphy",
     .   "Physiq dQl","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
C
         CALL histend(nid_bilKPave)
c
         ndex2d = 0
         ndex3d = 0
c
      ENDIF ! fin de test sur ok_journe
