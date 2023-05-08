!
! $Id: ini_histISCCP.h 1403 2010-07-01 09:02:53Z fairhead $
!
      IF (ok_isccp) THEN
c
c$OMP MASTER
      ndex2d = 0
      ndex3d = 0
c
c pour les champs instantannes, il faut mettre la meme valeur pour
c zout et zsto.
c dtime est passe par ailleurs a histbeg
c zstophy = frequence de stockage des champs tous les pdt physiques
c zout = frequence d'ecriture des champs
cIM 300505     zstophy = dtime 
c appel du simulateur toutes les 3heures
!IM on lit la frequence d'appel dans physiq.def
!         zcals(1) = dtime *6.  !toutes les 3h (en s)
          zcals(1) = freq_ISCCP !toutes les freq_ISCCP secondes
        DO n=1, napisccp
          zcalh(n) = zcals(n)/3600. !stoutes les Xh (en heures)
        ENDDO !n
c
c ecriture 8 fois par jour
c       zout = dtime * REAL(NINT(86400./dtime*ecrit_isccp))
c ecriture toutes les 2h (12 fois par jour)
c       zout = dtime * 4.
c ecriture toutes les 1/2 h (48 fois par jour)
c       zout = dtime
c
c       IF(freqout_isccp.EQ.1.) THEN
c ecriture jounaliere
!IM on ecrit les resultats du simulateur ISCCP toutes les 
! ecrit_ISCCP secondes      zout_isccp(1) = ecrit_day !(en s)
          zout_isccp(1) = ecrit_ISCCP !(en s)
c ecriture mensuelle
c         zout = dtime * ecrit_mth !(en s)
        DO n=1, napisccp 
          zoutj(n)=zout_isccp(n)/86400. !(en jours)
c
c le nombre de sous-colonnes ncol : ncol=(100.*zcalh)/zoutd
          ncol(n)=NINT((100.*zcalh(n))/zoutj(n))
          IF(ncol(n).GT.ncolmx) THEN
           PRINT*,'Warning: Augmenter le nombre colonnes du simulateur'
           PRINT*,'         ISCCP ncol=', ncol,' ncolmx=',ncolmx
c          PRINT*,'n ncol',n,ncol(n)
           CALL abort
          ENDIF
c
        DO l=1, ncol(n)
          vertlev(l,n)=REAL(l)
        ENDDO !ncol
c
        ENDDO !n

c       PRINT*, 'La frequence de sortie ISCCP est de ', ecrit_isccp
c
        idayref = day_ref
        CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
c       write(*,*)'ISCCP ', itau_phy, zjulian
c
c
c definition coordonnees lon,lat en globale
c
cym        CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlon,zx_lon)
cym        DO i = 1, iim
cym          zx_lon(i,1) = rlon(i+1)
cym          zx_lon(i,jjmp1) = rlon(i+1)
cym        ENDDO

cym        CALL gr_fi_ecrit(1,klon,iim,jjmp1,rlat,zx_lat)
c
cIM BEG region
cym Desole dans un premier temps le mode region ne marchera pas
cym Il faudra voir dans un second temps pour l'implementer
cym Mais cela posera des problemes au niveau de la reconstruction

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
          print*,'On stoke le fichier histISCCP sur ',
     s   imin_ins,imax_ins,jmin_ins,jmax_ins
cym          print*,'On stoke le fichier histISCCP instantanne sur ',
cym     s   zx_lon(imin_ins,1),zx_lon(imax_ins,1),
cym     s   zx_lat(1,jmin_ins),zx_lat(1,jmax_ins)
cIM END region
c
        IF(1.EQ.0) THEN
cym         CALL histbeg("histISCCP.nc", iim,zx_lon(:,1),jjmp1,zx_lat(1,:),
cym     .                 1, iim, 1, jjmp1,
cym     .                 itau_phy, zjulian, dtime,
cym     .                 nhori, nid_isccp)
         CALL histbeg_phy("histISCCP.nc", itau_phy, zjulian, dtime,
     .                 nhori, nid_isccp)
        ENDIF !(1.EQ.0) THEN
c
cym         CALL histbeg("histISCCP.nc", iim,zx_lon(:,1),
cym     .                 jjmp1,zx_lat(1,:),
cym     .                 imin_ins,imax_ins-imin_ins+1,
cym     .                 jmin_ins,jmax_ins-jmin_ins+1,
cym     .                 itau_phy, zjulian, dtime,
cym     .                 nhori, nid_isccp)

         CALL histbeg_phy("histISCCP.nc", itau_phy, zjulian, dtime,
     .                 nhori, nid_isccp)
c
        IF(type_run.EQ."ENSP".OR.type_run.EQ."CLIM") THEN
         CALL histvert(nid_isccp, "cldtopres","Cloud Top Pressure","mb",
     .                 lmaxm1, cldtopres, nvert,'down')
        ELSE IF(type_run.EQ."AMIP".OR.type_run.EQ."CFMI") THEN
         CALL histvert(nid_isccp,"cldtopres3","Cloud Top Pressure","mb",
     .                 lmax3, cldtopres3, nvert3,'down')
        ENDIF
        DO n=1, napisccp
         CALL histvert(nid_isccp, "Nbcol"//verticaxe(n),
     .        "Nb of Column"//verticaxe(n),"1",
     .        ncol(n), vertlev(:,n), nvlev(n),'up')
        ENDDO
c
        IF(type_run.EQ."ENSP".OR.type_run.EQ."CLIM") THEN
c
c variables a ecrire
c 
         DO n=1, napisccp
c
         DO k=1, kmaxm1
          CALL histdef(nid_isccp, "cldISCCP_"//taulev(k)//verticaxe(n),
     .                "LMDZ ISCCP cld", "%",
     .                iim, jj_nb,nhori,lmaxm1,1,lmaxm1,nvert,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
         ENDDO
c
         CALL histdef(nid_isccp, "nsunlit"//verticaxe(n),
     .                "Nb of calls with sunlit ", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
         CALL histdef(nid_isccp, "meantaucld"//verticaxe(n),
     .                "ISCCP mean cloud optical thickness", "1",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
         ENDDO
c
        ELSE IF(type_run.EQ."AMIP".OR.type_run.EQ."CFMI") THEN
c
         DO n=1, napisccp
c
c         print*,'n=',n,' avant histdef(..LMDZ ISCCP cld'
c
          DO k=1, kmaxm1
           DO l=1, lmaxm1
c
           CALL histdef(nid_isccp, pclev(l)//taulev(k)//verticaxe(n),
     .                "LMDZ ISCCP cld "//cnameisccp(l,k), "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
           ENDDO
          ENDDO
c
c         print*,'n=',n,' avant histdef(..Nb of calls sunlit'
          CALL histdef(nid_isccp, "nsunlit"//verticaxe(n),
     .                "Nb of calls with sunlit ", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
         CALL histdef(nid_isccp, "meantaucld"//verticaxe(n),
     .                "ISCCP mean cloud optical thickness", "1",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
c 9types de nuages ISCCP-D2
          CALL histdef(nid_isccp, "cirr",
     .                "Cirrus lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "cist",
     .                "CiSt lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "deep",
     .                "Deep lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "alcu",
     .                "AlCu lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "alst",
     .                "AlSt lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "nist",
     .                "NiSt lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "cumu",
     .                "Cumu lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "stcu",
     .                "StCu lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "stra",
     .                "Stra lk ISCCP-D2", "%",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
c 3_epaisseurs_optiques x3_pressions_au_sommet_des_nuages  types de nuages 
          CALL histdef(nid_isccp, "thin",
     .                "Opt. thin ISCCP-D2 like clouds", "%",
     .                iim, jj_nb,nhori,lmax3,1,lmax3,nvert3,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "mid",
     .                "Opt. intermediate ISCCP-D2 like clouds", "%",
     .                iim, jj_nb,nhori,lmax3,1,lmax3,nvert3,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
          CALL histdef(nid_isccp, "thick",
     .                "Opt. thick ISCCP-D2 like clouds", "%",
     .                iim, jj_nb,nhori,lmax3,1,lmax3,nvert3,32,
     .                "ave(X)", zcals(n),zout_isccp(n))
c
c        IF(1.EQ.0) THEN
c        IF(n.EQ.3) THEN
c        IF(n.EQ.1) THEN
c
cIM 070905 BEG
         IF(1.EQ.0) THEN
          print*,'n=',n,' avant histdef(..boxptop axe'
cIM verif boxptop
          CALL histdef(nid_isccp,"boxptop"//verticaxe(n),
     .                "Boxptop axe"//verticaxe(n), "mb",
     .                iim, jj_nb,nhori,
     .                ncol(n),1,ncol(n),nvlev(n),32,
cIM  .                ncolmx,1,ncolmx,nvlev,32,
cIM  .                "inst(X)",dtime,dtime)
     .                "ave(X)",zcals(n),zout_isccp(n))
         ENDIF !(1.EQ.0) THEN
cIM 070905 END
c        ENDIF !(n.EQ.3) THEN
c       ENDIF !(1.EQ.0) THEN
c
c         print*,'n=',n,' avant histdef(..seed axe'
          CALL histdef(nid_isccp, "seed"//verticaxe(n),
     .                "seed axe"//verticaxe(n), "-",
     .                iim, jj_nb,nhori,1,1,1,-99,32,
cIM  .                "inst(X)", dtime,dtime)
     .                "ave(X)", zcals(n),zout_isccp(n))
c
         ENDDO !n
        ENDIF 
        CALL histend(nid_isccp)
c
c$OMP END MASTER
      ENDIF ! ok_isccp
