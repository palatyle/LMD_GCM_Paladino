c
      dt_cum = dtvr*day_step

!      zan = annee_ref
!      dayref = day_ref
!      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
      tau0 = itau_dyn
c
       pi = 4.0 * ATAN(1.0)
       degres = 180./pi
       rlong = rlonu * degres
       rlatg = rlatu * degres
c
      CALL histbeg("paramLMDZ_dyn.nc", 
     .                 iip1,rlong, jjp1,rlatg,
     .                 1,1,1,1,
     .                 tau0, jD_ref+jH_ref , dt_cum,
     .                 thoriid, nid_ctesGCM)
c
         CALL histdef(nid_ctesGCM, "prt_level", 
     .        "Niveau impression debuggage dynamique",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "dayref", 
     .        "Jour de l etat initial ( = 350  si 20 Decembre par ex.)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "anneeref", 
     .        "Annee de l etat initial",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "anneelim", 
     .        "Annee du fichier limitxxxx.nc  si  ok_limitvrai =y",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "raz_date", 
     .   "Remise a zero (raz) date init.: 0 pas de raz;1=date gcm.def",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "nday", 
     .   "Nombre de jours d integration",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "day_step", 
     .   "nombre de pas par jour pour dt = 1 min",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iperiod", 
     .   "periode pour le pas Matsuno (en pas de temps)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iapp_tracvl", 
     .   "frequence du groupement des flux (en pas de temps)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iconser", 
     .  "periode de sortie des variables de controle (en pas de temps)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iecri", 
     .  "periode d ecriture du fichier histoire (en jour)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "periodav", 
     .  "periode de stockage fichier histmoy (en jour)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "dissip_period", 
     .  "periode de la dissipation (en pas) ... a completer",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "lstardis", 
     .  "choix de l operateur de dissipation: 1= star,0=non-star ??",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "nitergdiv", 
     .  "nombre d iterations de l operateur de dissipation gradiv",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "nitergrot", 
     .  "nombre d iterations de l operateur de dissipation nxgradrot",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "niterh", 
     .  "nombre d iterations de l operateur de dissipation divgrad",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "tetagdiv", 
     ."temps dissipation des + petites long. d ondes pour u,v (gradiv)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "tetagrot", 
     ."temps diss. des + petites long. d ondes pour u,v (nxgradrot)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "tetatemp", 
     ."temps diss. des + petites long. d ondes pour h (divgrad)",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "coefdis", 
     ."coefficient pour gamdissip",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "purmats", 
     ."Choix schema integration temporel: 1=Matsuno,0=Matsuno-leapfrog",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "ok_guide", 
     ."Guidage: 1=true ,0=false",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "true_calendar", 
     ."Choix du calendrier",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "guide_calend", 
     ."Guidage calendrier gregorien: 1=oui ,0=non",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iflag_phys", 
     ."Permet de faire tourner le modele sans physique: 1=avec ,0=sans",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "iphysiq", 
     ."Periode de la physique en pas de temps de la dynamique",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "clon", 
     ."longitude en degres du centre du zoom",
     .                "deg",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "clat", 
     ."latitude en degres du centre du zoom",
     .                "deg",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "grossismx", 
     ."facteur de grossissement du zoom, selon la longitude",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "grossismy", 
     ."facteur de grossissement du zoom, selon la latitude",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "fxyhypb", 
     ."Fonction f(y) hyperbolique  si true=1, sinusoidale si false=0",
     .                "-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "dzoomx", 
     ."extension en longitude de la zone du zoom (fraction zone totale)"
     .                ,"-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "dzoomy", 
     ."extension en latitude de la zone du zoom (fraction zone totale)"
     .                ,"-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "taux", 
     ."raideur du zoom en  X"
     .                ,"-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "tauy", 
     ."raideur du zoom en  Y"
     .                ,"-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "ysinus", 
     ."ysinus=1: Ftion f(y) avec y=Sin(latit.)/ ysinus=0: y = latit"
     .                ,"-",iip1,jjp1,thoriid, 1,1,1, -99, 32, 
     .                "once", dt_cum,dt_cum)
c
         CALL histdef(nid_ctesGCM, "ip_ebil_dyn", 
     ."PRINTlevel for energy conservation diag.; 0/1= pas de print,
     . 2= print","-",iip1,jjp1,thoriid, 1,1,1, -99, 32,
     .                "once", dt_cum,dt_cum)
c
c=================================================================
c
         CALL histend(nid_ctesGCM)
c
c=================================================================
