!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!..include cles_phys.h
!
       LOGICAL cycle_diurne,soil_model 
       LOGICAL ok_orodr,ok_orolf,ok_gw_nonoro
       INTEGER nbapp_rad, nbapp_chim, iflag_con, iflag_ajs
       REAL    ecriphy
       INTEGER lev_histmth, lev_histday
       REAL    solaire

! Parametres pour PBL:
       REAL    z0, lmixmin
       REAL    ksta
       LOGICAL ok_kzmin

! Parametres surface:
       REAL    inertie,emis

! Parametres Chimie:
       logical chimi,ylellouch,hcnrad
       integer vchim,aerprod,htoh2
       
! Parametres Microphysique:
       integer microfi,cutoff,clouds
       real    tx,tcorrect,p_prodaer
       real    xnuf 
       REAL    xvis,xir


       COMMON/clesphys_i/                                               &
     &     nbapp_rad, nbapp_chim, iflag_con, iflag_ajs,                 &
     &     lev_histmth, lev_histday, vchim,aerprod,htoh2,               &
     &     microfi,cutoff,clouds

       COMMON/clesphys_r/                                               &
     &     ecriphy, solaire, z0, lmixmin, ksta, inertie, emis,          &
     &     tx,tcorrect,p_prodaer,xnuf,xvis,xir

       COMMON/clesphys_l/cycle_diurne, soil_model,                      &
     &     ok_orodr, ok_orolf, ok_gw_nonoro, ok_kzmin,                  &
     &     chimi,ylellouch,hcnrad

