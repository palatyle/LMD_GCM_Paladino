!-------------------------------------------------------------------
! INCLUDE tabcontrol.h
!-------------------------------------------------------------------

      INTEGER radpas,chimpas ! frequences d'appel rayonnement, chimie
      REAL dtime             ! pas temporel de la physique
      REAL lsinit            ! Solar longitude in the startphy file

! tableau de controle
      INTEGER        length
      PARAMETER    ( length = 100 )
      REAL tabcntr0( length       )


      COMMON/ctltab_i/radpas,chimpas
      COMMON/ctltab_r/dtime,lsinit
      COMMON/ctltab/tabcntr0
