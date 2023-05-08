c
c $Header$
c
c
cIM on initialise les variables 
c
        CALL ini_undefSTD(itap,freq_outNMC)
c
cIM on interpole les champs sur les niveaux STD de pression 
cIM a chaque pas de temps de la physique
c
c-------------------------------------------------------c
c positionnement de l'argument logique a .false.        c
c pour ne pas recalculer deux fois la meme chose !      c
c a cet effet un appel a plevel_new a ete deplace       c
c a la fin de la serie d'appels                         c
c la boucle 'DO k=1, nlevSTD' a ete internalisee        c
c dans plevel_new, d'ou la creation de cette routine... c
c-------------------------------------------------------c
c
        CALL plevel_new(klon,klev,nlevSTD,.true.,pplay,rlevSTD,
     &              t_seri,tlevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             u_seri,ulevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             v_seri,vlevSTD)
c

c
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zphi/RG,philevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             qx(:,:,ivap),qlevSTD)
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_rh*100.,rhlevSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=u_seri(i,l)*v_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,uvSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*q_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,vqSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,vTSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=omega(i,l)*qx(i,l,ivap)
         ENDDO !i 
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,wqSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*zphi(i,l)/RG
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,vphiSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=omega(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,wTSTD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=u_seri(i,l)*u_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,u2STD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=v_seri(i,l)*v_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,v2STD)
c
        DO l=1, klev
         DO i=1, klon
          zx_tmp_fi3d(i,l)=t_seri(i,l)*t_seri(i,l)
         ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,T2STD)

c
      zx_tmp_fi3d(:,:)=wo(:,:,1) * dobson_u * 1e3 / zmasse / rmo3 * rmd
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,O3STD)
c
      if (read_climoz == 2) THEN
      zx_tmp_fi3d(:,:)=wo(:,:,2) * dobson_u * 1e3 / zmasse / rmo3 * rmd
        CALL plevel_new(klon,klev,nlevSTD,.false.,pplay,rlevSTD,
     &             zx_tmp_fi3d,O3daySTD)
      endif
c
        DO l=1, klev
        DO i=1, klon
         zx_tmp_fi3d(i,l)=paprs(i,l)
        ENDDO !i
        ENDDO !l
        CALL plevel_new(klon,klev,nlevSTD,.true.,zx_tmp_fi3d,rlevSTD,
     &             omega,wlevSTD)
c
cIM on somme les valeurs toutes les freq_calNMC secondes
c
       CALL undefSTD(itap,freq_calNMC, read_climoz)
c
cIM on moyenne a la fin du mois ou du jour (toutes les freq_outNMC secondes)
c
       CALL moy_undefSTD(itap,freq_outNMC,freq_moyNMC)
c
       CALL plevel(klon,klev,.true.,pplay,50000.,
     &              zphi/RG,geo500)

cIM on interpole a chaque pas de temps le SWup(clr) et SWdn(clr) a 200 hPa
c
      CALL plevel(klon,klevp1,.true.,paprs,20000.,
     $     swdn0,SWdn200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     swdn,SWdn200)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     swup0,SWup200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     swup,SWup200)
c
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     lwdn0,LWdn200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     lwdn,LWdn200)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     lwup0,LWup200clr)
      CALL plevel(klon,klevp1,.false.,paprs,20000.,
     $     lwup,LWup200)
c
