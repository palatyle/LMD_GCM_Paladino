!-------------------------------------------------------------------
! INCLUDE tabcontrol.h
!-------------------------------------------------------------------

      INTEGER radpas,chimpas ! frequences d'appel rayonnement, chimie
      REAL dtime             ! pas temporel de la physique

! tableau de controle
      INTEGER        length
      PARAMETER    ( length = 100 )
      REAL tabcntr0( length       )


      COMMON/ctltab_i/radpas,chimpas
      COMMON/ctltab_r/dtime
      COMMON/ctltab/tabcntr0
