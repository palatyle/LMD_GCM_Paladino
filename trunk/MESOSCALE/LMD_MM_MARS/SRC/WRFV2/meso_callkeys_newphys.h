!
! For Fortran 77/Fortran 90 compliance always use line continuation
! symbols '&' in columns 73 and 6
!
      COMMON/callkeys/callrad,calldifv,calladj,callcond,callsoil,       &
     &   season,diurnal,iradia,lwrite,calllott                          &
     &   ,iaervar,iddist,topdustref,callstats,calleofdump               &
     &   ,callnirco2,callnlte,callthermos,callconduct,                  &
     &        calleuv,callmolvis,callmoldiff,thermochem                 &
     &   , solarcondate, thermoswater                                   &
     &   , semi                                                         &
     &   , callemis                                                     &
     &   , callg2d                                                      &
     &   , linear                                                       &
     &   , ilwd                                                         &
     &   , ilwb                                                         &
     &   , ilwn                                                         &
     &   , ncouche                                                      &
     &   , alphan                                                       &
     &   , rayleigh                                                     &
     &   , tracer                                                       &
     &   , dustbin,active,doubleq,lifting,callddevil,scavenging         &
     &   , sedimentation,activice,water,caps                            &
     &   , photochem,nqchem_min                                         

      LOGICAL callrad,calldifv,calladj,callcond,callsoil,               &
     &   season,diurnal,lwrite,calllott                                 &
     &   ,callstats,calleofdump                                         &
     &   ,callnirco2,callnlte,callthermos,callconduct,                  &
     &    calleuv,callmolvis,callmoldiff,thermochem,thermoswater


      logical callemis
      logical callg2d
      logical linear

      real topdustref
      real semi
      real alphan
      real solarcondate

      integer ecri_phys 
      integer iddist
      integer iaervar
      integer iradia
      integer ilwd
      integer ilwb
      integer ilwn
      integer ncouche

      logical rayleigh
      logical tracer
      integer dustbin
      logical active,doubleq,lifting,callddevil,scavenging
      logical sedimentation,activice,water,caps
      !!! plus besoin de iceparty ??
      logical photochem
      integer nqchem_min

      integer swrtype ! type of short wave (solar wavelength) radiative
      ! transfer to use 1: Fouquart 2: Toon.
      parameter (swrtype=2)
!      parameter (swrtype=2)
