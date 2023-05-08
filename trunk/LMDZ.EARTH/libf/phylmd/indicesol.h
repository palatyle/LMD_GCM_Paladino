!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
      INTEGER nbsrf
      PARAMETER (nbsrf=4) ! nombre de sous-fractions pour une maille
!
      INTEGER is_oce
      PARAMETER (is_oce=3) ! ocean
      INTEGER is_sic
      PARAMETER (is_sic=4) ! glace de mer
      INTEGER is_ter
      PARAMETER (is_ter=1) ! terre
      INTEGER is_lic
      PARAMETER (is_lic=2) ! glacier continental
!
      REAL epsfra
      PARAMETER (epsfra=1.0E-05)
!
      CHARACTER(len=3) clnsurf(nbsrf)
      DATA clnsurf/'ter', 'lic', 'oce', 'sic'/
      SAVE clnsurf
!$OMP THREADPRIVATE(clnsurf)
