c
c $Id: ini_bilKP_ins.h 1403 2010-07-01 09:02:53Z fairhead $
c
      IF (ok_journe) THEN
c
         zsto = dtime
         zout = dtime
         typeval=tinst
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
c
cIM 280405 BEG
c
cIM cf. AM 081204 BEG region
          imin_ins=1
          imax_ins=iim
          jmin_ins=1
          jmax_ins=jjmp1
cym          do i=1,iim-1
cym             if(zx_lon(i,1).lt.lonmin_ins) imin_ins=i
cym             if(zx_lon(i,1).le.lonmax_ins) imax_ins=i+1
cym          enddo
cym          do j=1,jjmp1
cym             if(zx_lat(1,j).ge.latmin_ins) jmax_ins=j
cym             if(zx_lat(1,j).gt.latmax_ins) jmin_ins=j
cym          enddo
c
          print*,'On stoke le fichier bilKP instantanne sur ',
     s   imin_ins,imax_ins,jmin_ins,jmax_ins
          print*,'On stoke le fichier bilKP instantanne sur ',
     s   zx_lon(imin_ins,1),zx_lon(imax_ins,1),
     s   zx_lat(1,jmin_ins),zx_lat(1,jmax_ins)
cIM cf. AM 081204 END region
c
cIM 280405 END
c
cym         IF(1.EQ.0) THEN
cym         CALL histbeg("histbilKP_ins", iim,zx_lon(:,1), jjmp1,
cym     .                zx_lat(1,:),
cym     .                1,iim,1,jjmp1, itau_phy, zjulian, dtime,
cym     .                nhori, nid_bilKPins)
         ENDIF
c
cIM 280405 BEG
c
cIM cf. AM 081204 BEG region
cym         CALL histbeg("histbilKP_ins", iim,zx_lon(:,1), 
cym     .                 jjmp1,zx_lat(1,:),
cym     .                 imin_ins,imax_ins-imin_ins+1,
cym     .                 jmin_ins,jmax_ins-jmin_ins+1,
cym     .                 itau_phy, zjulian, dtime,
cym     .                 nhori, nid_bilKPins)
         CALL histbeg_phy("histbilKP_ins", itau_phy, zjulian, dtime,
     .                 nhori, nid_bilKPins)
cIM 081204 END
c
cIM 280405 END
c
         write(*,*)'Journee ', itau_phy, zjulian
         CALL histvert(nid_bilKPins, "presnivs",
     .                "Vertical levels","mb",
     .                 klev, presnivs/100., nvert)
c
c Champs 3D:
c
         CALL histdef(nid_bilKPins,"ue",
     .   "Zonal energy transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"ve",
     .   "Merid energy transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"uq",
     .   "Zonal humidity transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"vq",
     .   "Merid humidity transport","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32, 
     .                typeval, zsto,zout)
c
c Champs 3D:
c
         CALL histdef(nid_bilKPins, "temp",
     .   "Air temperature", "K",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"ovap",
     .   "Specific humidity","Kg/Kg",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"geop",
     .   "Geopotential height", "m",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"vitu", 
     .   "Zonal wind", "m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"vitv", 
     .   "Meridional wind", "m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "vitw",
     .   "Vertical wind", "m/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "pres",
     .   "Inter-Layer Air pressure",
     .                "Pa",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "play",
     .   "Mean-Layer Air pressure",
     .                "Pa",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "oliq",
     .   "Liquid water content", 
     .                "kg/kg",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtdyn",
     .   "Dynamics dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dqdyn",
     .   "Dynamics dQ", "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtcon",
     .   "Convection dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "ducon",
     .   "Convection du", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dvcon",
     .   "Convection dv", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dqcon",
     .   "Convection dQ","Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtlsc",
     .   "Condensation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dqlsc",
     .   "Condensation dQ","Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dtvdf",
     .   "Boundary-layer dT","K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dqvdf", 
     .   "Boundary-layer dQ", 
     .               "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dtajs",
     .   "Ajustement sec dT","K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dqajs",
     .   "Ajustement sec dQ", 
     .               "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dteva",
     .   "Reevaporation dT","K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dqeva",
     .   "Reevaporation dQ",
     .                "Kg/Kg/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)

c
         CALL histdef(nid_bilKPins, "dtswr", 
     .   "SW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtsw0", 
     .   "SW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtlwr", 
     .   "LW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dtlw0", 
     .   "LW radiation dT", "K/s",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"duvdf",
     .   "Boundary-layer dU","m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins,"dvvdf",
     .   "Boundary-layer dV","m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         IF (ok_orodr) THEN
         IF (ok_orolf) THEN
         CALL histdef(nid_bilKPins, "duoli", 
     .   "Orography dU", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         CALL histdef(nid_bilKPins, "dvoli", 
     .   "Orography dV", "m/s2",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
c
         ENDIF
         ENDIF
C
         CALL histdef(nid_bilKPins, "duphy",
     .   "Physiq dU","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPins, "dvphy",
     .   "Physiq dV","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPins, "dtphy",
     .   "Physiq dT","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPins, "dqphy",
     .   "Physiq dQ","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
C
         CALL histdef(nid_bilKPins, "dqlphy",
     .   "Physiq dQl","-",
     .                iim,jjphy_nb,nhori, klev,1,klev,nvert, 32,
     .                typeval, zsto,zout)
cIM 280405 BEG
c
c Champs 2D:
c
c u850, v850
c        DO k=1, nlevSTD
         DO k=1, 12
c
         IF(k.GE.2.AND.k.LE.12) bb2=clevSTD(k)
         IF(k.GE.13.AND.k.LE.17) bb3=clevSTD(k)
c
         IF(bb2.EQ."850") THEN 
c
          CALL histdef(nid_bilKPins, "u"//bb2,
     .                 "Zonal wind "//bb2//"mb","m/s",
     .                iim,jjphy_nb,nhori, 1,1,1, -99, 32,
     .                typeval, zsto,zout)
c
          CALL histdef(nid_bilKPins, "v"//bb2,
     .                 "Meridional wind "//bb2//"mb","m/s",
     .                iim,jjphy_nb,nhori, 1,1,1, -99, 32,
     .                typeval, zsto,zout)
c
         ENDIF !(bb2.EQ."850") 
c
         ENDDO !k=1, 12
c
cIM 280405 END
c
         CALL histend(nid_bilKPins)
c
         ndex2d = 0
         ndex3d = 0
c
      ENDIF ! fin de test sur ok_journe
