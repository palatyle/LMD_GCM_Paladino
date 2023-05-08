
!
! $Id: clesphys.h 1421 2010-07-28 13:28:42Z jghattas $
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!
!..include cles_phys.h
!
       LOGICAL cycle_diurne,soil_model,new_oliq,ok_orodr,ok_orolf 
       LOGICAL ok_limitvrai
       INTEGER nbapp_rad, iflag_con
       REAL co2_ppm, co2_ppm0, solaire
       REAL(kind=8) RCO2, RCH4, RN2O, RCFC11, RCFC12  
       REAL(kind=8) CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt

!OM ---> correction du bilan d'eau global
!OM Correction sur precip KE
       REAL cvl_corr
!OM Fonte calotte dans bilan eau
       LOGICAL ok_lic_melt

!IM simulateur ISCCP 
       INTEGER top_height, overlap
!IM seuils cdrm, cdrh
       REAL cdmmax, cdhmax
!IM param. stabilite s/ terres et en dehors
       REAL ksta, ksta_ter
!IM ok_kzmin : clef calcul Kzmin dans la CL de surface cf FH
       LOGICAL ok_kzmin
!IM, MAFo fmagic, pmagic : parametres - additionnel et multiplicatif - 
!                          pour regler l albedo sur ocean
       REAL fmagic, pmagic
! Hauteur (imposee) du contenu en eau du sol
           REAL qsol0
! Frottement au sol (Cdrag)
       Real f_cdrag_ter,f_cdrag_oce
! Rugoro
       Real f_rugoro

!IM lev_histhf  : niveau sorties 6h
!IM lev_histday : niveau sorties journalieres
!IM lev_histmth : niveau sorties mensuelles
!IM lev_histdayNMC : on peut sortir soit sur 8 (comme AR5) ou bien 
!                    sur 17 niveaux de pression
       INTEGER lev_histhf, lev_histday, lev_histmth
       INTEGER lev_histdayNMC
       Integer lev_histins, lev_histLES  


!IM ok_histNMC  : sortie fichiers niveaux de pression (histmthNMC, histdayNMC, histhfNMC)
!IM freq_outNMC : frequences de sortie fichiers niveaux de pression (histmthNMC, histdayNMC, histhfNMC)
!IM freq_calNMC : frequences de calcul fis. hist*NMC.nc
!IM pasphys : pas de temps de physique (secondes)
       REAL pasphys
       LOGICAL ok_histNMC(3)
       REAL freq_outNMC(3) , freq_calNMC(3)

       CHARACTER(len=4) type_run
! aer_type: pour utiliser un fichier constant dans readaerosol 
       CHARACTER*8 :: aer_type 
       LOGICAL ok_isccp, ok_regdyn
       REAL lonmin_ins, lonmax_ins, latmin_ins, latmax_ins
       REAL ecrit_ins, ecrit_hf, ecrit_hf2mth, ecrit_day
       REAL ecrit_mth, ecrit_tra, ecrit_reg 
       REAL ecrit_LES
       REAL freq_ISCCP, ecrit_ISCCP
       REAL freq_COSP
       LOGICAL :: ok_cosp,ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP
       INTEGER :: ip_ebil_phy, iflag_rrtm
       LOGICAL :: ok_strato
       LOGICAL :: ok_hines

       COMMON/clesphys/cycle_diurne, soil_model, new_oliq,              &
     &     ok_orodr, ok_orolf, ok_limitvrai, nbapp_rad, iflag_con       &
     &     , co2_ppm, solaire, RCO2, RCH4, RN2O, RCFC11, RCFC12         &
     &     , CH4_ppb, N2O_ppb, CFC11_ppt, CFC12_ppt                     &
     &     , top_height, overlap, cdmmax, cdhmax, ksta, ksta_ter        &
     &     , ok_kzmin, fmagic, pmagic                                   &
     &     , f_cdrag_ter,f_cdrag_oce,f_rugoro                           &
     &     , lev_histhf, lev_histday, lev_histmth                       &

     &     , lev_histins, lev_histLES, lev_histdayNMC                   &
     &     , pasphys, ok_histNMC, freq_outNMC, freq_calNMC              &
     &     , type_run, ok_isccp, ok_regdyn, ok_cosp                     &
     &     , ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP                     &
     &     , lonmin_ins, lonmax_ins, latmin_ins, latmax_ins             &
     &     , ecrit_ins, ecrit_hf, ecrit_hf2mth, ecrit_day               &
     &     , ecrit_mth, ecrit_tra, ecrit_reg                            &
     &     , freq_ISCCP, ecrit_ISCCP, freq_COSP, ip_ebil_phy            &
     &     , ok_lic_melt, cvl_corr, aer_type                            &
     &     , qsol0, iflag_rrtm, ok_strato,ok_hines,ecrit_LES            &
     &     , co2_ppm0
     
!$OMP THREADPRIVATE(/clesphys/)
 
