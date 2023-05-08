!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
      REAL cld_lc_lsc,cld_lc_con
      REAL cld_tau_lsc,cld_tau_con
      REAL ffallv_lsc,ffallv_con
      REAL coef_eva
      LOGICAL reevap_ice
      INTEGER iflag_pdf

      common/comfisrtilp/                                               &
     &     cld_lc_lsc                                                   &
     &     ,cld_lc_con                                                  &
     &     ,cld_tau_lsc                                                 &
     &     ,cld_tau_con                                                 &
     &     ,ffallv_lsc                                                  &
     &     ,ffallv_con                                                  &
     &     ,coef_eva                                                    &
     &     ,reevap_ice                                                  &
     &     ,iflag_pdf        

!$OMP THREADPRIVATE(/comfisrtilp/)
