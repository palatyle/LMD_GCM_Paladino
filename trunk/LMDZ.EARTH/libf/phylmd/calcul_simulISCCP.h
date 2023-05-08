c
c $Id: calcul_simulISCCP.h 1403 2010-07-01 09:02:53Z fairhead $
c
c on appelle le simulateur ISCCP toutes les 3h
c et on fait des sorties 1 fois par jour 
c
c ATTENTION : le temps de calcul peut augmenter considerablement !
c =============================================================== c
      DO n=1, napisccp
c
      nbapp_isccp=30 !appel toutes les 15h
cIM 170107      isccppas=NINT((itap*dtime)/3600.) !Nb. d'heures de la physique
      freqin_pdt(n)=ifreq_isccp(n)
c
cIM initialisation nbsunlit pour calculs simulateur ISCCP pdt la journee
c
      DO i=1, klon
       sunlit(i)=1 
       IF(rmu0(i).EQ.0.) sunlit(i)=0
       nbsunlit(1,i,n)=REAL(sunlit(i))
      ENDDO
c
cIM calcul tau, emissivite nuages convectifs
c
      convfra(:,:)=rnebcon(:,:)
      convliq(:,:)=rnebcon(:,:)*clwcon(:,:)
c
      CALL newmicro (paprs, pplay,ok_newmicro,
     .            t_seri, convliq, convfra, dtau_c, dem_c,
     .            cldh_c, cldl_c, cldm_c, cldt_c, cldq_c,
     .            flwp_c, fiwp_c, flwc_c, fiwc_c,
     e            ok_aie,
     e            mass_solu_aero, mass_solu_aero_pi,
     e            bl95_b0, bl95_b1,
     s            cldtaupi, re, fl)
c
cIM calcul tau, emissivite nuages startiformes
c
      CALL newmicro (paprs, pplay,ok_newmicro,
     .            t_seri, cldliq, cldfra, dtau_s, dem_s,
     .            cldh_s, cldl_s, cldm_s, cldt_s, cldq_s,
     .            flwp_s, fiwp_s, flwc_s, fiwc_s,
     e            ok_aie,
     e            mass_solu_aero, mass_solu_aero_pi,
     e            bl95_b0, bl95_b1,
     s            cldtaupi, re, fl)
c
      cldtot(:,:)=min(max(cldfra(:,:),rnebcon(:,:)),1.)
c
cIM inversion des niveaux de pression ==> de haut en bas
c
      CALL haut2bas(klon, klev, pplay, pfull)
      CALL haut2bas(klon, klev, q_seri, qv)
      CALL haut2bas(klon, klev, cldtot, cc)
      CALL haut2bas(klon, klev, rnebcon, conv)
      CALL haut2bas(klon, klev, dtau_s, dtau_sH2B)
      CALL haut2bas(klon, klev, dtau_c, dtau_cH2B)
      CALL haut2bas(klon, klev, t_seri, at)
      CALL haut2bas(klon, klev, dem_s, dem_sH2B)
      CALL haut2bas(klon, klev, dem_c, dem_cH2B)
      CALL haut2bas(klon, klevp1, paprs, phalf)
c
cIM: initialisation de seed
c
        DO i=1, klon
c
         aa=ABS(paprs(i,2)-NINT(paprs(i,2)))
         seed_re(i,n)=1000.*aa+1.
         seed(i,n)=NINT(seed_re(i,n))
c
         IF(seed(i,n).LT.50) THEN
c          print*,'seed<50 avant i seed itap paprs',i,
c    .     seed(i,n),itap,paprs(i,2)
           seed(i,n)=50+seed(i,n)+i+itap
           seed_old(i,n)=seed(i,n)
c
           IF(itap.GT.1) then
            IF(seed(i,n).EQ.seed_old(i,n)) THEN
             seed(i,n)=seed(i,n)+10
             seed_old(i,n)=seed(i,n)
            ENDIF
           ENDIF
c
c          print*,'seed<50 apres i seed itap paprs',i,
c    .     seed(i,n),itap,paprs(i,2)
c
          ELSE IF(seed(i,n).EQ.0) THEN
           print*,'seed=0 i paprs aa seed_re',
     .     i,paprs(i,2),aa,seed_re(i,n)
           abort_message = ''
           CALL abort_gcm (modname,abort_message,1)
          ELSE IF(seed(i,n).LT.0) THEN
           print*,'seed < 0, i seed itap paprs',i,
     .     seed(i,n),itap,paprs(i,2)
           abort_message = ''
           CALL abort_gcm (modname,abort_message,1)
          ENDIF
c
        ENDDO
c     
cIM: pas de debug, debugcol
      debug=0
      debugcol=0
c
cIM o500 ==> distribution nuage ftion du regime dynamique (vit. verticale a 500 hPa)
c
        DO k=1, klevm1
        kp1=k+1
c       PRINT*,'k, presnivs',k,presnivs(k), presnivs(kp1)
        if(presnivs(k).GT.50000.AND.presnivs(kp1).LT.50000.) THEN
         DO i=1, klon
          o500(i)=omega(i,k)*RDAY/100.
c         if(i.EQ.1) print*,' 500hPa lev',k,presnivs(k),presnivs(kp1)
         ENDDO
         GOTO 1000
        endif 
1000  continue
      ENDDO
c
cIM recalcule les nuages vus par satellite, via le simulateur ISCCP
c
      CALL ISCCP_CLOUD_TYPES(
     &     debug,
     &     debugcol,
     &     klon,
     &     sunlit,
     &     klev,
     &     ncol(n),
     &     seed(:,n),
     &     pfull,
     &     phalf,
     &     qv, cc, conv, dtau_sH2B, dtau_cH2B,
     &     top_height,
     &     overlap,
     &     tautab,
     &     invtau,
     &     ztsol,
     &     emsfc_lw,
     &     at, dem_sH2B, dem_cH2B,
     &     fq_isccp(:,:,:,n),
     &     totalcldarea(:,n),
     &     meanptop(:,n),
     &     meantaucld(:,n),
     &     boxtau(:,:,n),
     &     boxptop(:,:,n))
c
      ENDDO !n=1, napisccp

