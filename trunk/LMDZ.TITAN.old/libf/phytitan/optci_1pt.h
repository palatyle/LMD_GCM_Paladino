!----------------------------------------------
! fichier optci_1p.h
!
!  regroupe les sorties du fichier optci_1pt
!  les nuages oblige a quasiment doubler
!  le nombre de variables c'est juste pour
!  la lisibilite.
!----------------------------------------------

! ----DIAGNOSTIQUES
      real   TAUHID_1pt(NLAYER,NSPECI)
      real   TAUCID_1pt(NLAYER,NSPECI)
      real   TAUGID_1pt(NLAYER,NSPECI)
! ----OPACITES TOTALES
      real   TAUHI_1pt(NSPECI)
      real   TAUCI_1pt(NSPECI)
      real   TAUGI_1pt(NSPECI)
! ----COLONNE NUAGEUSE
      real   DTAUI_1pt(NLAYER,NSPECI)
      real   TAUI_1pt(NLEVEL,NSPECI)
      real   WBARI_1pt(NLAYER,NSPECI)
      real   COSBI_1pt(NLAYER,NSPECI)
! ----COLONNE "CLAIRE"
      real   DTAUIP_1pt(NLAYER,NSPECI)
      real   TAUIP_1pt(NLEVEL,NSPECI)
      real   WBARIP_1pt(NLAYER,NSPECI)
      real   COSBIP_1pt(NLAYER,NSPECI)


      common/opti_1pt/                                                  &
     &        TAUHID_1pt                                                &
     &       ,TAUCID_1pt                                                &
     &       ,TAUGID_1pt                                                &
! ----OPACITES TOTALES
     &       ,TAUHI_1pt                                                 &
     &       ,TAUCI_1pt                                                 &
     &       ,TAUGI_1pt                                                 &
! ----COLONNE NUAGEUSE
     &       ,DTAUI_1pt                                                 &
     &       ,TAUI_1pt                                                  &
     &       ,WBARI_1pt                                                 &
     &       ,COSBI_1pt                                                 &
! ----COLONNE "CLAIRE"
     &       ,DTAUIP_1pt                                                &
     &       ,TAUIP_1pt                                                 &
     &       ,WBARIP_1pt                                                &
     &       ,COSBIP_1pt
