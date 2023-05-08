!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
! A1.0 Fundamental constants
      REAL RPI,RCLUM,RHPLA,RKBOL,RNAVO
! A1.1 Astronomical constants
      REAL RDAY,REA,REPSM,RSIYEA,RSIDAY,ROMEGA
! A1.1.bis Constantes concernant l'orbite de la Terre:
      REAL R_ecc, R_peri, R_incl
! A1.2 Geoide
      REAL RA,RG,R1SA
! A1.3 Radiation
!     REAL RSIGMA,RI0
      REAL RSIGMA
! A1.4 Thermodynamic gas phase
      REAL R,RMD,RMV,RD,RV,RCPD,RCPV,RCVD,RCVV
      REAL RKAPPA,RETV
! ADAPTATION GCM POUR CP(T)
!      real cpdet
!      external cpdet
! A1.5,6 Thermodynamic liquid,solid phases
      REAL RCW,RCS
! A1.7 Thermodynamic transition of phase
      REAL RLVTT,RLSTT,RLMLT,RTT,RATM
! A1.8 Curve of saturation
      REAL RESTT,RALPW,RBETW,RGAMW,RALPS,RBETS,RGAMS
      REAL RALPD,RBETD,RGAMD
!
      COMMON/YOMCST/RPI ,RCLUM, RHPLA, RKBOL, RNAVO ,RDAY  ,REA         &
     & ,REPSM ,RSIYEA,RSIDAY,ROMEGA , R_ecc, R_peri, R_incl             &
     & ,RA    ,RG ,R1SA                                                 &
     & ,RSIGMA,R ,RMD   ,RMV   ,RD    ,RV    ,RCPD ,RCPV,RCVD           &
     & ,RCVV  ,RKAPPA,RETV ,RCW   ,RCS ,RLVTT ,RLSTT ,RLMLT ,RTT ,RATM  &
     & ,RESTT ,RALPW ,RBETW ,RGAMW ,RALPS ,RBETS ,RGAMS ,RALPD ,RBETD   &
     & ,RGAMD
