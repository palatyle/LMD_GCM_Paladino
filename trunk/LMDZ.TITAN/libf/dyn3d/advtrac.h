!
! Tracers in the dynamics (E.M, oct.2008)
!
      COMMON/advtr/iadv,tnom
      INTEGER :: iadv(nqmx)           ! tracer advection scheme number
      CHARACTER (len=20) :: tnom(nqmx) ! tracer name
      ! NB: tnom length should be compatible with noms() in physics ...
