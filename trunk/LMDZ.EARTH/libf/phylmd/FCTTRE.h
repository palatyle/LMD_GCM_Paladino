!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!     ------------------------------------------------------------------
!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!     ------------------------------------------------------------------
      REAL PTARG, PDELARG, P5ARG, PQSARG, PCOARG
      REAL FOEEW, FOEDE, qsats, qsatl, dqsats, dqsatl
      LOGICAL thermcep
      PARAMETER (thermcep=.TRUE.)
!
      FOEEW ( PTARG,PDELARG ) = EXP (                                   &
     &          (R3LES*(1.-PDELARG)+R3IES*PDELARG) * (PTARG-RTT)        &
     & / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG)) )
!
      FOEDE ( PTARG,PDELARG,P5ARG,PQSARG,PCOARG ) = PQSARG*PCOARG*P5ARG &
     & / (PTARG-(R4LES*(1.-PDELARG)+R4IES*PDELARG))**2
!
      qsats(ptarg) = 100.0 * 0.622 * 10.0                               &
     &           ** (2.07023 - 0.00320991 * ptarg                       &
     &           - 2484.896 / ptarg + 3.56654 * LOG10(ptarg))
      qsatl(ptarg) = 100.0 * 0.622 * 10.0                               &
     &           ** (23.8319 - 2948.964 / ptarg                         &
     &           - 5.028 * LOG10(ptarg)                                 &
     &           - 29810.16 * EXP( - 0.0699382 * ptarg)                 &
     &           + 25.21935 * EXP( - 2999.924 / ptarg))
!
      dqsats(ptarg,pqsarg) = RLVTT/RCPD*pqsarg * (3.56654/ptarg         &
     &                     +2484.896*LOG(10.)/ptarg**2                  &
     &                     -0.00320991*LOG(10.))
      dqsatl(ptarg,pqsarg) = RLVTT/RCPD*pqsarg*LOG(10.)*                &
     &                (2948.964/ptarg**2-5.028/LOG(10.)/ptarg           &
     &                +25.21935*2999.924/ptarg**2*EXP(-2999.924/ptarg)  &
     &                +29810.16*0.0699382*EXP(-0.0699382*ptarg))
