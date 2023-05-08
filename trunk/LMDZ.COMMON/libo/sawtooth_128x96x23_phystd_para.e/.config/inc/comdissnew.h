!
! $Id: comdissnew.h 1319 2010-02-23 21:29:54Z fairhead $
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!
!-----------------------------------------------------------------------
! INCLUDE 'comdissnew.h'

      COMMON/comdissnew/ lstardis,nitergdiv,nitergrot,niterh,tetagdiv,  &
     &                   tetagrot,tetatemp,coefdis, vert_prof_dissip

      LOGICAL lstardis
      INTEGER nitergdiv, nitergrot, niterh

! For the Earth model:
      integer vert_prof_dissip ! vertical profile of horizontal dissipation
!     Allowed values:
!     0: rational fraction, function of pressure
!     1: tanh of altitude

      REAL     tetagdiv, tetagrot,  tetatemp, coefdis

!
! ... Les parametres de ce common comdissnew sont  lues par defrun_new 
!              sur le fichier  run.def    ....
!
!-----------------------------------------------------------------------
