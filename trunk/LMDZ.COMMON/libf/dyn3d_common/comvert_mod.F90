MODULE comvert_mod

IMPLICIT NONE  

PRIVATE
INCLUDE "dimensions.h"

PUBLIC :: ap,bp,presnivs,dpres,sig,ds,pa,preff,nivsigs,nivsig, &
          aps,bps,scaleheight,pseudoalt,disvert_type, pressure_exner

REAL ap(llm+1) ! hybrid pressure contribution at interlayers
REAL bp (llm+1) ! hybrid sigma contribution at interlayer
REAL presnivs(llm) ! (reference) pressure at mid-layers
REAL dpres(llm)
REAL sig(llm+1)
REAL ds(llm)
REAL pa ! reference pressure (Pa) at which hybrid coordinates
        ! become purely pressure (more or less)
REAL preff  ! reference surface pressure (Pa)
REAL nivsigs(llm)
REAL nivsig(llm+1)
REAL aps(llm) ! hybrid pressure contribution at mid-layers
REAL bps(llm) ! hybrid sigma contribution at mid-layers
REAL scaleheight ! atmospheric (reference) scale height (km)
REAL pseudoalt(llm) ! pseudo-altitude of model levels (km), based on presnivs(),
                     ! preff and scaleheight

INTEGER disvert_type ! type of vertical discretization:
                     ! 1: Earth (default for planet_type==earth),
                     !     automatic generation
                     ! 2: Planets (default for planet_type!=earth),
                     !     using 'z2sig.def' (or 'esasig.def) file

LOGICAL pressure_exner
!     compute pressure inside layers using Exner function, else use mean
!     of pressure values at interfaces

END MODULE comvert_mod
