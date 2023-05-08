!
! $Header$
!

      IF (ok_regdyn) THEN
      
        if (is_sequential) then
c
cIM      PRINT*, 'La frequence de sortie REGDYN est de ', ecrit_mth
c
         idayref = day_ref
         CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c
c axe vertical pour les differents niveaux des histogrammes
      DO iw=1, iwmax
        zx_o500(iw)=wmin+(iw-1./2.)*pas_w
      ENDDO

         CALL histbeg("histREGDYN", kmaxm1,zx_tau, lmaxm1,zx_pc,
     .                 1,kmaxm1,1,lmaxm1, itau_phy, zjulian, dtime, 
     .                 nhoriRD, nid_regdyn)

         CALL histvert(nid_regdyn, "omeganivs", "Omega levels", 
     .                 "mb/day",
     .                 iwmax, zx_o500, komega)
c
c   pour les champs instantannes, il faut mettre la meme valeur pour
c   zout et zsto.
c   dtime est passe par ailleurs a histbeg
c
c        zout = dtime * REAL(NINT(86400./dtime*ecrit_regdyn))
c        zsto = zout
c        print*,'zout,zsto=',zout,zsto
c
c stockage a chaque pas de temps de la physique
c
         zstophy = dtime
cIM 020904      zstophy = dtime * nbapp_isccp

c ecriture mensuelle
c
         zout = dtime * ecrit_mth
cIM 020904      
c        zout = dtime * ecrit_day
c        zout = dtime * REAL(NINT(86400./dtime*ecrit_regdyn))

c
c Champs 3D:
c
c TROP
         CALL histdef(nid_regdyn, "hw1", "Tropics Histogram ", "%",
     &                kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega, 32, 
     &                "ave(X)", zstophy,zout)

         CALL histdef(nid_regdyn, "nh1", "Nb of pixels Tropics Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c

         CALL histdef(nid_regdyn, "nht1",
     &                "Total Nb pixels Tropics Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c
c PAN
         CALL histdef(nid_regdyn, "hw2", "North Pacific Histogram", "%",
     &                kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega, 32, 
     &                "ave(X)", zstophy,zout)

         CALL histdef(nid_regdyn, "nh2", "Nb of pixels North Pacific",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c

         CALL histdef(nid_regdyn, "nht2",
     &                "Total Nb pixels North Pacific Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c CAL
         CALL histdef(nid_regdyn, "hw3", "California Histogram", "%",
     &                kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega, 32, 
     &                "ave(X)", zstophy,zout)

         CALL histdef(nid_regdyn, "nh3", 
     &                "Nb of pixels California Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c

         CALL histdef(nid_regdyn, "nht3",
     &                "Total Nb pixels California Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c HAW
         CALL histdef(nid_regdyn, "hw4", "Hawai Histogram", "%",
     &                kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega, 32, 
     &                "ave(X)", zstophy,zout)

         CALL histdef(nid_regdyn, "nh4", "Nb of pixels Hawai Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c

         CALL histdef(nid_regdyn, "nht4",
     &                "Total Nb pixels Hawai Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c WAP
         CALL histdef(nid_regdyn, "hw5", "Warm Pool Histogram", "%",
     &                kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega, 32, 
     &                "ave(X)", zstophy,zout)

         CALL histdef(nid_regdyn, "nh5", "Nb of pixels Warm Pool Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c

         CALL histdef(nid_regdyn, "nht5",
     &                "Total Nb pixels Warm Pool Histo",
     &                "%",kmaxm1,lmaxm1,nhoriRD, iwmax,1,iwmax, komega,
     &                32,"ave(X)", zstophy,zout)
c
         CALL histend(nid_regdyn)
	 
	 endif ! is_sequential

      endif ! ok_regdyn
