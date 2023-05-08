c
c $Id: ini_histday_seri.h 1403 2010-07-01 09:02:53Z fairhead $
c
cym Ne fonctionnera pas en mode parallele
      IF (is_sequential) THEN
      
      IF (type_run.EQ."AMIP") THEN
c
       zstophy = dtime
       zout = ecrit_day
c
         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c
         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlon,zx_lon)
         DO i = 1, iim
            zx_lon(i,1) = rlon(i+1)
            zx_lon(i,jjmp1) = rlon(i+1)
         ENDDO
         DO ll=1,klev
            znivsig(ll)=REAL(ll)
         ENDDO
         CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)
c
         imin_debut=1 
         nbpti=1
         jmin_debut=1 
         nbptj=1
c
         CALL histbeg("histday_seri.nc", 
     .                 iim,zx_lon(:,1), jjmp1,zx_lat(1,:),
     .                 imin_debut,nbpti,jmin_debut,nbptj,
     .                 itau_phy, zjulian, dtime,
     .                 nhori, nid_day_seri)
c
         CALL histvert(nid_day_seri, "presnivs", 
     .                "Vertical levels","mb",
     .                 klev, presnivs/100., nvert)
c
         CALL histdef(nid_day_seri, "bilTOA", 
     .                "Net radiation at model top", "W/m2",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "bils", 
     .                "Net downward energy flux at surface","W/m2",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "ecin", 
     .                "Total kinetic energy (per unit area)","J/m2",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
cIM 151004 BEG
         IF(1.EQ.0) THEN
c
         CALL histdef(nid_day_seri, "momang", 
     .               "Total relative angular momentum (per unit area)",
     .               "kg/s",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "frictor", 
     .               "Friction torque (per unit area)", "N/m",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "mountor", 
     .               "Mountain torque (per unit area)", "N/m",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
         ENDIF !(1.EQ.0) THEN
c
         CALL histdef(nid_day_seri, "momang", 
     .               "Axial angular momentum (per unit area)",
     .               "kg/s",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "torsfc", 
     .        "Total surface torque (including mountain torque)", "N/m",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
cIM 151004 END        
c
         CALL histdef(nid_day_seri, "tamv", 
     .                "Temperature (mass-weighted vert. ave)", "K",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32,
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "psol", 
     .                "Surface pressure", "Pa",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
         CALL histdef(nid_day_seri, "evap", 
     .                "Evaporation and sublimation (per unit area)", 
     .                "kg/(m2*s)",
     .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
     .                "ave(X)", zstophy,zout)
c
c          call histdef(nid_day_seri, 
c    .         "SnowFrac", 
c    .         "Snow-covered area ", "%",  
c    .         iim,jjmp1,nhori, 1,1,1, -99, 32,
c    .         "ave(X)", zstophy,zout)
c
c        CALL histdef(nid_day_seri, "snow_depth", 
cIM 080904  .                "Snow Depth (water equivalent)", "m",
cIM 191104  .                "Snow Depth (water equivalent)", "kg/m2",
c    .                "Snow Mass", "kg/m2",
c    .                iim,jjmp1,nhori, 1,1,1, -99, 32, 
c    .               "ave(X)", zstophy,zout)
c
           call histdef(nid_day_seri, 
     .         "tsol_"//clnsurf(is_oce), 
     .         "SST over open (ice-free) ocean ", "K",  
     .         iim,jjmp1,nhori, 1,1,1, -99, 32,
     .         "ave(X)", zstophy,zout)
c
c=================================================================
c
         CALL histend(nid_day_seri)
c
c=================================================================
      ENDIF ! fin de test sur type_run.EQ.AMIP
      
      ENDIF ! is_sequential
