!
! $Id: $
!
!
! handle debug and output levels
! lunout:    unit of file where outputs will be sent
!                           (default is 6, standard output)
! prt_level: level of informative output messages (0 = minimum)
!
      INTEGER lunout, prt_level
      COMMON /comprint/ lunout, prt_level
