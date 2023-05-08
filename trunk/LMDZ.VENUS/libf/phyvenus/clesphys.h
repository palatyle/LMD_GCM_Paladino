!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!..include cles_phys.h
!
       LOGICAL cycle_diurne,soil_model
       LOGICAL ok_orodr,ok_orolf,ok_gw_nonoro
       LOGICAL ok_kzmin
       LOGICAL callnlte,callnirco2,callthermos
       LOGICAL ok_cloud, ok_chem, reinit_trac, ok_sedim
       LOGICAL ok_clmain, physideal, startphy_file
       INTEGER nbapp_rad, nbapp_chim, iflag_con, iflag_ajs
       INTEGER lev_histins, lev_histday, lev_histmth
       INTEGER tr_scheme, cl_scheme
       INTEGER nircorr, nltemodel, solvarmod
       INTEGER nb_mode
       REAL    ecriphy
       REAL    solaire
       REAL    z0, lmixmin
       REAL    ksta, inertie
       REAL    euveff, solarcondate

       COMMON/clesphys_l/ cycle_diurne, soil_model,                     &
     &     ok_orodr, ok_orolf, ok_gw_nonoro, ok_kzmin,                  &
     &     callnlte,callnirco2,callthermos,                             &
     &     ok_cloud, ok_chem, reinit_trac, ok_sedim,                    &
     &     ok_clmain, physideal, startphy_file

       COMMON/clesphys_i/ nbapp_rad, nbapp_chim,                        &
     &     iflag_con, iflag_ajs,                                        &
     &     lev_histins, lev_histday, lev_histmth, tr_scheme,            &
     &     cl_scheme, nircorr, nltemodel, solvarmod, nb_mode

       COMMON/clesphys_r/ ecriphy, solaire, z0, lmixmin,                &
     &     ksta, inertie, euveff, solarcondate

