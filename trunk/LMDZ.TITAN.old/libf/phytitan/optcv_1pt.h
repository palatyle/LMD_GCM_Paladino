!----------------------------------------------
! fichier optcv_1p.h
!
!  regroupe les sorties du fichier optcv_1pt
!  les nuages oblige a quasiment doubler
!  le nombre de variables c'est juste pour
!  la lisibilite.
!----------------------------------------------

! ----DIAGNOSTIQUES 
      real   TAUHVD_1pt(NLAYER,NSPECV)
      real   TAUCVD_1pt(NLAYER,NSPECV)
      real   TAUGVD_1pt(NLAYER,NSPECV)
! ----OPACITES TOTALES
      real   TAUHV_1pt(NSPECV)
      real   TAUCV_1pt(NSPECV)
      real   TAURV_1pt(NSPECV)
      real   TAUGV_1pt(NSPECV)
! ----COLONNE NUAGEUSE
      real   TAUV_1pt(NLEVEL,NSPECV,4)
      real   DTAUV_1pt(NLAYER,NSPECV,4)
      real   WBARV_1pt(NLAYER,NSPECV,4)
      real   COSBV_1pt(NLAYER,NSPECV,4)
! ----COLONNE "CLAIRE"
      real   TAUVP_1pt(NLEVEL,NSPECV,4)
      real   DTAUVP_1pt(NLAYER,NSPECV,4)
      real   WBARVP_1pt(NLAYER,NSPECV,4)
      real   COSBVP_1pt(NLAYER,NSPECV,4) 

      common/optv_1pt/                                                  &
     &        TAUHVD_1pt                                                &
     &       ,TAUCVD_1pt                                                &
     &       ,TAUGVD_1pt                                                &
! ----OPACITES TOTALES
     &       ,TAUHV_1pt                                                 &
     &       ,TAUCV_1pt                                                 &
     &       ,TAUGV_1pt                                                 &
! ----COLONNE NUAGEUSE
     &       ,DTAUV_1pt                                                 &
     &       ,TAUV_1pt                                                  &
     &       ,WBARV_1pt                                                 &
     &       ,COSBV_1pt                                                 &
! ----COLONNE "CLAIRE"
     &       ,DTAUVP_1pt                                                &
     &       ,TAUVP_1pt                                                 &
     &       ,WBARVP_1pt                                                &
     &       ,COSBVP_1pt

